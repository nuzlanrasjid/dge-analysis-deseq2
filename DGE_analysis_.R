# DGE analysis
# setwd("D:/DGE")


# 1. Load libraries
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(apeglm)

# 2. Load data
count_data <- read.csv("count_data.csv", header = TRUE, row.names = 1)
coldata    <- read.csv("sample_data.csv", header = TRUE, row.names = 1)

# Sanity checks: sample names in coldata must match count_data columns,
# in the same order, or DESeq2 will silently misassign samples to conditions.
stopifnot(all(colnames(count_data) %in% rownames(coldata)))
stopifnot(all(colnames(count_data) == rownames(coldata)))

# set factor levels
coldata$Treatment <- factor(coldata$Treatment)
coldata$Sequencing <- factor(coldata$Sequencing)

# 3. Build DESeq2 dataset
# Design includes "Sequencing" (library type: single vs paired-end) as a
# blocking variable so it doesn't confound the Treatment effect we care about.
dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData   = coldata,
  design    = ~ Sequencing + Treatment
)

# Explicitly set "untreated" as the reference level so log2FoldChange is
# interpreted as treated vs. untreated.
dds$Treatment <- factor(dds$Treatment, levels = c("untreated", "treated"))

# Filter out genes with very low total counts across all samples.
# This is a light pre-filter for speed/memory only -- DESeq2's independent
# filtering (applied later via results()) does the statistically-informed
# filtering for the actual test.
keep <- rowSums(counts(dds)) >= 5
dds  <- dds[keep, ]

# NOTE: this dataset has no technical replicates to collapse. If your own
# data does, use collapseReplicates() here before running DESeq().
# 4. Run DESeq2
dds <- DESeq(dds)
res <- results(dds)
summary(res)

# Stricter FDR threshold for comparison with alpha 5%
res0.05 <- results(dds, alpha = 0.05)
summary(res0.05)

res0.05 <- as.data.frame(res0.05)
res     <- as.data.frame(res)

# Order by p-value for quick inspection of top hits
reorderres0.05 <- res0.05[order(res0.05$pvalue), ]
head(reorderres0.05)

# 5. Filter significant genes
# Two-step filter: adjusted p-value (FDR) < 0.05, then |log2FC| > 1
# (i.e. at least a 2-fold change). This threshold is reused later in the
# volcano plot so the two are consistent.
LFC_THRESHOLD <- 1
PADJ_THRESHOLD <- 0.05

filtered <- res %>%
  filter(padj < PADJ_THRESHOLD) %>%
  filter(abs(log2FoldChange) > LFC_THRESHOLD)

write.csv(res, "res.all.csv")
write.csv(filtered, "filtereddata.filter.csv")

normalized <- counts(dds, normalized = TRUE)
write.csv(normalized, "normalized_data.csv")

# 6. Visualization
## 6.1 Dispersion plot -- sanity check that DESeq2's dispersion shrinkage
## behaved as expected (gene-wise estimates shrunk toward the fitted curve)
plotDispEsts(dds)

## 6.2 PCA plot
## Variance-stabilizing transformation, used here purely for visualization
## (not for the differential expression test itself).
vds <- vst(dds, blind = FALSE)
plotPCA(vds, intgroup = c("Sequencing", "Treatment"))

## 6.3 Sample-to-sample distance heatmap
sampledis    <- dist(t(assay(vds)))
sampledismat <- as.matrix(sampledis)
color <- colorRampPalette(rev(brewer.pal(9, "Reds")))(300)

pheatmap(sampledismat,
         clustering_distance_rows = sampledis,
         clustering_distance_cols = sampledis,
         col = color)

## 6.4 Heatmap of top 10 genes by adjusted p-value (rlog-transformed counts)
top10 <- res[order(res$padj), ][1:10, ]
top10 <- rownames(top10)

rld <- rlog(dds, blind = FALSE)

pheatmap(assay(rld)[top10, ],
         cluster_rows = FALSE, show_rownames = TRUE, cluster_cols = FALSE)

# With sample annotation (Sequencing + Treatment)
annot <- as.data.frame(colData(dds)[, c("Sequencing", "Treatment")])
pheatmap(assay(rld)[top10, ],
         cluster_rows = FALSE, show_rownames = TRUE, cluster_cols = FALSE,
         annotation_col = annot)

## 6.5 Z-score heatmap of top 10 genes
cal_z_score   <- function(x) (x - mean(x)) / sd(x)
allzscore     <- t(apply(normalized, 1, cal_z_score))
subset_zscore <- allzscore[top10, ]
pheatmap(subset_zscore)

## 6.6 MA plot (raw, then shrunk with apeglm to remove noisy low-count genes)
plotMA(dds, ylim = c(-2, 2))

resLFC <- lfcShrink(dds, coef = "Treatment_treated_vs_untreated", type = "apeglm")
plotMA(resLFC, ylim = c(-2, 2))
resLFC <- as.data.frame(resLFC)

## 6.7 Volcano plot
## Uses the same significance thresholds defined in section 5, and labels
## the top 10 most significant genes.
resLFC$diffexpressed <- "NO"
resLFC$diffexpressed[resLFC$log2FoldChange >  LFC_THRESHOLD & resLFC$padj < PADJ_THRESHOLD] <- "UP"
resLFC$diffexpressed[resLFC$log2FoldChange < -LFC_THRESHOLD & resLFC$padj < PADJ_THRESHOLD] <- "DOWN"

resLFC$delabel <- NA
top10_by_padj <- rownames(resLFC[order(resLFC$padj), ])[1:10]
resLFC[top10_by_padj, "delabel"] <- top10_by_padj

ggplot(data = resLFC, aes(x = log2FoldChange, y = -log10(pvalue),
                          col = diffexpressed, label = delabel)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values = c("UP" = "red", "DOWN" = "blue", "NO" = "grey60")) +
  theme(text = element_text(size = 20))



# 7. Biological interpretations
## 7.1 Define gene sets
### use all significant genes, not only the top10 genes, will be splitted by direction
sig <- res%>%filter(!is.na(padj), padj < PADJ_THRESHOLD)

up_genes <- rownames(sig[sig$log2FoldChange> 0, ])
down_genes <- rownames(sig[sig$log2FoldChange< 0, ])


### setting universe, make sure that every gene was actually tested (padj not NA), NOT the whole genome.
### This is the correct background for the hypergeometric test.
universe <- rownames(res[!is.na(res$padj), ])
universe
### check how many genes up and down regulated
length(up_genes); length(down_genes)

head(rownames(count_data))

## 7.2 GO Enrichment process
ego_up <- enrichGO(gene = up_genes,
         universe = universe,
         OrgDb = org.Dm.eg.db,
         keyType = "FLYBASE",
         ont = "BP",
         pAdjustMethod = "BH",
         pvalueCutoff = 0.05,
         qvalueCutoff = 0.05,
         readable = TRUE)

ego_down <- enrichGO(gene = down_genes,
                   universe = universe,
                   OrgDb = org.Dm.eg.db,
                   keyType = "FLYBASE",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

### Remove redundant parent/child terms (makes plots much cleaner)
ego_up_s   <- clusterProfiler::simplify(ego_up,   cutoff = 0.7, by = "p.adjust")
ego_down_s <- clusterProfiler::simplify(ego_down, cutoff = 0.7, by = "p.adjust")

### export file
write.csv(as.data.frame(ego_up_s), "GO_BP_up.csv")
write.csv(as.data.frame(ego_down_s), "GO_BP_down.csv")

## 7.3 KEGG enrichment
map_id <- function(id){
  bitr(id, fromType = "FLYBASE", toType = "ENTREZID", OrgDb = org.Dm.eg.db)$ENTREZID
}
map

options(timeout = 300) ### Internet network connection improvement

kegg_up <- enrichKEGG(gene = map_id(up_genes),
           universe = map_id(universe),
           organism = "dme",
           keyType = "ncbi-geneid",
           pvalueCutoff = 0.05)

kegg_down <- enrichKEGG(gene = map_id(down_genes),
                      universe = map_id(universe),
                      organism = "dme",
                      keyType = "ncbi-geneid",
                      pvalueCutoff = 0.05)

### Convert Entrez IDs back to gene symbols in the results
kegg_up <- setReadable(kegg_up, OrgDb = org.Dm.eg.db, keyType = "ENTREZID")
kegg_down <- setReadable(kegg_down, OrgDb = org.Dm.eg.db, keyType = "ENTREZID")
kegg_down

### checking why kegg_down unidentified
kegg_down_all <- enrichKEGG(gene = map_id(down_genes), universe = map_id(universe),
                            organism = "dme", keyType = "ncbi-geneid",
                            pvalueCutoff = 1, qvalueCutoff = 1)
head(as.data.frame(kegg_down_all)[, c("Description","GeneRatio","BgRatio","pvalue","p.adjust")], 10)

### Export files
write.csv(as.data.frame(kegg_up), "kegg_up.csv")
write.csv(as.data.frame(kegg_down), "kegg_down.csv")

## 7.4 visualization
### 7.4.1 GO visualization
dotplot(ego_up_s, showCategory = 15) + ggtitle("GO_BP: Upregulated")
dotplot(ego_down_s, showCategory = 15) + ggtitle("GO_BP: Downregulated")
barplot(ego_up_s,   showCategory = 15)
barplot(ego_down_s,   showCategory = 15)

### 7.4.2 KEGG visualization
dotplot(kegg_up,   showCategory = 15) + ggtitle("KEGG: up-regulated")
dotplot(kegg_down, showCategory = 15) + ggtitle("KEGG: down-regulated")

###  7.4.3 Gene-concept network (genes colored by log2FC)
FC <- setNames(res$log2FoldChange, rownames(res))
cnetplot(ego_up_s, showCategory = 15, foldChange = FC)

# Enrichment map: needs pairwise term similarity first
ego_up_sim <- pairwise_termsim(ego_up_s)
emapplot(ego_up_sim, showCategory = 20)

# Side-by-side up vs down in a single plot
cc <- compareCluster(geneClusters = list(Up = up_genes, Down = down_genes),
               fun = "enrichGO", universe = universe,
               OrgDb = org.Dm.eg.db, keyType = "FLYBASE",
               ont = "BP", pvalueCutoff = 0.05)
# remove term redundance
cc_s <- clusterProfiler::simplify(cc, cutoff = 0.7, by = "p.adjust")

dotplot(cc_s, showCategory = 6, label_format = 50, font.size = 11) +
  theme(axis.text.y = element_text(size = 10, lineheight = 0.9))

dotplot(cc, showCategory = 10) + theme(axis.text.y = element_text(size = 10),
                                       axis.text.x = element_text(size = 12),
                                       text = element_text(size = 10))


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



# Differential Gene Expression Analysis with DESeq2

A complete RNA-seq differential expression workflow in R using DESeq2, 
from raw count matrix to an annotated volcano plot, which was applied to the
public **pasilla** dataset (*Drosophila melanogaster*).

## Dataset

[GEO accession GSE18508](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18508)
(Brooks et al., 2011) — RNA-seq of *Drosophila melanogaster* S2-DRSC
cells with and without knockdown of the splicing factor *pasilla*.

| | |
|---|---|
| Samples | 7 (4 untreated, 3 treated) |
| Design | `~ Sequencing + Treatment` |
| Covariate | `Sequencing` (single-end vs. paired-end library) |
| Reference level | `untreated` |

`Sequencing` is included as a blocking covariate so differences in
library type don't get mistaken for a treatment effect.

## Requirements

- R 
- Packages: `tidyverse`, `DESeq2`, `pheatmap`, `ggplot2`, `ggrepel`,
  `RColorBrewer`, `apeglm`

\```r
install.packages(c("tidyverse", "ggplot2", "ggrepel", "RColorBrewer"))
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "pheatmap", "apeglm"))
\```

## Usage

1. Place `count_data.csv` and `sample_data.csv` in the project folder.
2. Set your working directory to that folder.
3. Run `DGE_analysis.R` top to bottom.

## Workflow

1. **Data preparation** - load counts and metadata, verify sample
   names/order match between the two files.
2. **Model design** - build a `DESeqDataSet` with design
   `~ Sequencing + Treatment`, set `untreated` as the reference level.
3. **Filtering & testing** - drop very low-count genes, run `DESeq()`,
   extract results at FDR (padj) < 0.05.
4. **Significant gene filtering** - padj < 0.05 **and**
   |log2FoldChange| > 1 (same thresholds reused in the volcano plot).
5. **Effect size shrinkage** - `apeglm` shrinkage of log2 fold changes
   for a cleaner MA plot.
6. **Visualization** - dispersion plot, PCA, sample-distance heatmap,
   top-10-gene heatmaps (plain, annotated, z-score), MA plot
   (before/after shrinkage), and an annotated volcano plot.

## Output files

| File | Description |
|---|---|
| `res.all.csv` | DESeq2 results for all tested genes |
| `filtereddata.filter.csv` | Significant genes only (padj < 0.05, \|log2FC\| > 1) |
| `normalized_data.csv` | DESeq2 size-factor-normalized counts |

## Results

*(Fill in after running the script - e.g. number of significant genes
found, direction of change, and what the PCA/heatmap patterns suggest.)*

## Notes

This is a practice/learning project using a public benchmark dataset,
built to demonstrate the standard DESeq2 differential expression
workflow rather than as an original research study.

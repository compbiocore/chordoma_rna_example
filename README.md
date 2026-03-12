# chordoma_rna_example

RNAseq differential expression analysis pipeline using the [airway](https://bioconductor.org/packages/release/data/experiment/html/airway.html) dataset as a worked example.

## Overview

The analysis (`rnaseq_analysis.qmd`) covers the following steps:

| Step | Tool | Output |
|------|------|--------|
| Differential expression | DESeq2 | `deseq2_results.csv` |
| Heatmap (top 50 variable genes) | ComplexHeatmap | inline figure |
| Volcano plot | EnhancedVolcano | inline figure |
| Gene-set enrichment (GSEA) | clusterProfiler `GSEA()` | `gsea_results.csv`, inline figure |
| Over-representation analysis (ORA) | clusterProfiler `enricher()` | `ora_results.csv`, inline figure |

Gene sets for both GSEA and ORA are the **CP:KEGG_MEDICUS** collection sourced from
[MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/) via the
[msigdbr](https://cran.r-project.org/package=msigdbr) package.

## Requirements

### R packages

Install the required packages from Bioconductor and CRAN:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "airway",
  "DESeq2",
  "ComplexHeatmap",
  "EnhancedVolcano",
  "clusterProfiler",
  "org.Hs.eg.db",
  "AnnotationDbi"
))

install.packages(c("msigdbr", "dplyr", "tibble", "circlize"))
```

> **ashr** (used for LFC shrinkage) is also required:
> ```r
> install.packages("ashr")
> ```

Alternatively, use the provided Dockerfile to build a container, changing the information between the <> characters as needed:
```
docker buildx build --platform linux/amd64,linux/arm64 --push -t <dockerhub_username>/<container_name>:<container_tag> .
```

Or Singularity pull to use on Oscar with Open On Demand RStudio on Singularity (https://docs.ccv.brown.edu/oscar/connecting-to-oscar/open-ondemand):

```
singularity pull <singularity_image_name>.sif docker://<dockerhub_username>/<container_name>:<container_tag>
```

## Usage

Open `rnaseq_analysis.qmd` in RStudio (or another Quarto-aware editor) to run the
analysis interactively chunk by chunk, or render it into a self-contained HTML report:

```r
quarto::quarto_render("rnaseq_analysis.qmd")
```

or from a terminal:

```bash
quarto render rnaseq_analysis.qmd
```

All CSV output files are written to the working directory. Figures are embedded directly
in the rendered HTML document.

## Dataset

The [airway](https://bioconductor.org/packages/release/data/experiment/html/airway.html)
package provides RNA-Seq read counts from an experiment on airway smooth muscle cells
treated with dexamethasone (Himes *et al.*, 2014, PLoS ONE). The contrast of interest is
**treated (`trt`) vs untreated (`untrt`)** while accounting for cell-line (`cell`) as a
covariate.

## Analysis details

### Differential expression (DESeq2)
- Design formula: `~ cell + dex`
- Pre-filtering: genes with fewer than 10 total counts removed
- LFC shrinkage with the `ashr` method

### Heatmap (ComplexHeatmap)
- Top 50 most variable genes selected after variance-stabilising transformation (VST)
- Rows Z-score normalised; columns annotated by treatment and cell line

### Volcano plot (EnhancedVolcano)
- X-axis: shrunken log2 fold-change
- Y-axis: adjusted p-value
- Thresholds: |LFC| > 1, padj < 0.05

### Pathway analysis (clusterProfiler + msigdbr)
- **Gene sets**: CP:KEGG_MEDICUS (MSigDB C2 category, *Homo sapiens*)
- **GSEA**: genes ranked by DESeq2 Wald statistic; `GSEA()` with `eps = 0`
- **ORA**: significant DE genes (padj < 0.05, |LFC| > 1) tested against all genes
  with a valid symbol as the background universe; `enricher()` used

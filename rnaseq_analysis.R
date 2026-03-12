## RNAseq Differential Expression Analysis
## Dataset: airway
## Tools: DESeq2, ComplexHeatmap, EnhancedVolcano, clusterProfiler, msigdbr

# ── 1. Load libraries ─────────────────────────────────────────────────────────

library(airway)
library(DESeq2)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(tibble)
library(circlize)
library(matrixStats)

# ── Configuration ─────────────────────────────────────────────────────────────

SEED <- 42

# Colour palettes used throughout the analysis
COL_DEX  <- c("trt" = "#E64B35", "untrt" = "#4DBBD5")
COL_CELL <- c("N061011" = "#00A087", "N052611" = "#3C5488",
              "N080611" = "#F39B7F", "N61311"  = "#8491B4")

# ── 2. Load the airway dataset ────────────────────────────────────────────────

data(airway)
se <- airway

# ── 3. DESeq2 differential expression analysis ───────────────────────────────

dds <- DESeqDataSet(se, design = ~ cell + dex)

# Pre-filter: keep only rows with at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]

dds <- DESeq(dds)

# Extract results: treated vs untreated
# Keep unshrunken results to retain the Wald statistic used for GSEA ranking
res_unshrunk <- results(dds, contrast = c("dex", "trt", "untrt"))
res_shrunk   <- lfcShrink(dds, contrast = c("dex", "trt", "untrt"),
                           res = res_unshrunk, type = "ashr")

res_df <- as.data.frame(res_shrunk) %>%
  rownames_to_column("ensembl_id") %>%
  mutate(stat = as.data.frame(res_unshrunk)[["stat"]])

# Map Ensembl IDs to gene symbols
res_df$symbol <- mapIds(
  org.Hs.eg.db,
  keys    = res_df$ensembl_id,
  keytype = "ENSEMBL",
  column  = "SYMBOL",
  multiVals = "first"
)

# Save full results table
write.csv(res_df, "deseq2_results.csv", row.names = FALSE)

cat("DESeq2 analysis complete.\n")
cat(sprintf("  Total genes tested: %d\n", nrow(res_df)))
cat(sprintf("  Significant (padj < 0.05, |LFC| > 1): %d\n",
            sum(!is.na(res_df$padj) & res_df$padj < 0.05 &
                  abs(res_df$log2FoldChange) > 1, na.rm = TRUE)))

# ── 4. ComplexHeatmap ─────────────────────────────────────────────────────────

# Variance-stabilising transformation for visualisation
vsd <- vst(dds, blind = FALSE)

# Select the 50 most variable genes
gene_vars  <- rowVars(assay(vsd))
top50_idx  <- order(gene_vars, decreasing = TRUE)[1:50]
mat        <- assay(vsd)[top50_idx, ]

# Z-score normalise rows
mat <- t(scale(t(mat)))

# Add gene symbols as row labels (fall back to Ensembl ID if symbol is NA)
row_symbols <- mapIds(
  org.Hs.eg.db,
  keys    = rownames(mat),
  keytype = "ENSEMBL",
  column  = "SYMBOL",
  multiVals = "first"
)
rownames(mat) <- ifelse(is.na(row_symbols), rownames(mat), row_symbols)

# Column annotations
col_df  <- as.data.frame(colData(dds)[, c("dex", "cell")])

# Define colour palette for the heatmap
heatmap_colors <- colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100)

heatmap_plot <- pheatmap(
  mat,
  color                    = heatmap_colors,
  annotation_col           = col_df,
  annotation_colors        = list(dex = COL_DEX, cell = COL_CELL),
  show_colnames            = TRUE,
  show_rownames            = TRUE,
  fontsize_row             = 8,
  clustering_distance_rows = "euclidean",
  clustering_method        = "complete",
  main                     = "Top 50 Most Variable Genes (airway)"
)

pdf("heatmap_top50_variable_genes.pdf", width = 10, height = 12)
draw(heatmap_plot)
dev.off()

cat("Heatmap saved to heatmap_top50_variable_genes.pdf\n")

# ── 5. EnhancedVolcano plot ───────────────────────────────────────────────────

res_volcano <- res_df %>%
  filter(!is.na(symbol)) %>%
  distinct(symbol, .keep_all = TRUE) %>%
  column_to_rownames("symbol")

sig_genes <- res_volcano %>%
  filter(!is.na(padj) & padj < 0.05 & !is.na(log2FoldChange))

top_up   <- sig_genes %>% arrange(desc(log2FoldChange)) %>% head(10)
top_down <- sig_genes %>% arrange(log2FoldChange)       %>% head(10)

genes_to_label <- unique(c(rownames(top_up), rownames(top_down)))

pdf("volcano_plot.pdf", width = 10, height = 10)
EnhancedVolcano(
  res_volcano,
  lab        = rownames(res_volcano),
  x          = "log2FoldChange",
  y          = "padj",
  title      = "Treated vs Untreated (airway)",
  subtitle   = "DESeq2 with LFC shrinkage (ashr)",
  pCutoff    = 0.05,
  FCcutoff   = 1,
  pointSize  = 2,
  labSize    = 3,
  col        = c("grey30", "#4DBBD5", "#00A087", "#E64B35"),
  legendLabels = c("NS", "Log2FC", "p-value", "p-value & Log2FC"),
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  selectLab  = genes_to_label
)
dev.off()

cat("Volcano plot saved to volcano_plot.pdf\n")

# ── 6. Build TERM2GENE and TERM2NAME from msigdbr (CP:KEGG_MEDICUS) ──────────

kegg_medicus <- msigdbr(
  species    = "Homo sapiens",
  category   = "C2",
  subcategory = "CP:KEGG_MEDICUS"
)

TERM2GENE <- kegg_medicus %>%
  select(gs_id, gene_symbol) %>%
  as.data.frame()

TERM2NAME <- kegg_medicus %>%
  select(gs_id, gs_name) %>%
  distinct() %>%
  as.data.frame()

cat(sprintf("CP:KEGG_MEDICUS gene sets loaded: %d pathways, %d gene-set associations\n",
            nrow(TERM2NAME), nrow(TERM2GENE)))

# ── 7. GSEA analysis (clusterProfiler::GSEA) ──────────────────────────────────

# Build a ranked gene list using the DESeq2 Wald statistic.
# Use genes with valid symbols only; keep one entry per symbol (highest |stat|).
ranked_df <- res_df %>%
  filter(!is.na(symbol), !is.na(stat)) %>%
  arrange(desc(abs(stat))) %>%
  distinct(symbol, .keep_all = TRUE)

gene_list        <- ranked_df$stat
names(gene_list) <- ranked_df$symbol
gene_list        <- sort(gene_list, decreasing = TRUE)

set.seed(SEED)
gsea_result <- GSEA(
  geneList     = gene_list,
  TERM2GENE    = TERM2GENE,
  TERM2NAME    = TERM2NAME,
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  eps          = 0,
  seed         = SEED,
  verbose      = FALSE
)

gsea_df <- as.data.frame(gsea_result)
write.csv(gsea_df, "gsea_results.csv", row.names = FALSE)

cat(sprintf("GSEA complete. Significant pathways (padj < 0.05): %d\n",
            nrow(gsea_df)))

if (nrow(gsea_df) > 0) {
  pdf("gsea_dotplot.pdf", width = 10, height = 8)
  print(dotplot(gsea_result, showCategory = 20, title = "GSEA – CP:KEGG_MEDICUS"))
  dev.off()
  cat("GSEA dotplot saved to gsea_dotplot.pdf\n")
}

# ── 8. ORA analysis (clusterProfiler::enricher) ───────────────────────────────

# Significant DE genes: padj < 0.05 and |LFC| > 1
sig_genes <- res_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1,
         !is.na(symbol)) %>%
  distinct(symbol) %>%
  pull(symbol)

cat(sprintf("ORA input: %d significant DE genes\n", length(sig_genes)))

# Background = all tested genes with a valid symbol
universe <- res_df %>%
  filter(!is.na(symbol)) %>%
  distinct(symbol) %>%
  pull(symbol)

ora_result <- enricher(
  gene         = sig_genes,
  universe     = universe,
  TERM2GENE    = TERM2GENE,
  TERM2NAME    = TERM2NAME,
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

ora_df <- as.data.frame(ora_result)
write.csv(ora_df, "ora_results.csv", row.names = FALSE)

cat(sprintf("ORA complete. Significant pathways (padj < 0.05): %d\n",
            nrow(ora_df)))

if (nrow(ora_df) > 0) {
  pdf("ora_dotplot.pdf", width = 10, height = 8)
  print(dotplot(ora_result, showCategory = 20, title = "ORA – CP:KEGG_MEDICUS"))
  dev.off()
  cat("ORA dotplot saved to ora_dotplot.pdf\n")
}

cat("\nAll analyses complete.\n")

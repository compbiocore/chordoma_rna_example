FROM rocker/tidyverse:4.5.2

RUN Rscript -e "\
  install.packages(c('circlize', 'matrixStats', 'msigdbr')); \
  if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); \
  BiocManager::install(c( \
    'airway', \
    'DESeq2', \
    'ComplexHeatmap', \
    'EnhancedVolcano', \
    'clusterProfiler', \
    'org.Hs.eg.db', \
    'AnnotationDbi' \
  ), ask = FALSE) \
"

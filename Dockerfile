FROM rocker/tidyverse:4.5.2

RUN apt-get update && apt-get install -y --no-install-recommends \
    libglpk-dev \
  && rm -rf /var/lib/apt/lists/*

RUN Rscript -e "\
  install.packages(c('ashr', 'circlize', 'matrixStats', 'msigdbr')); \
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

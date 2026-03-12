FROM rocker/tidyverse:4.5.2

ARG QUARTO_VERSION=1.6.42

RUN apt-get update && apt-get install -y --no-install-recommends \
    libglpk-dev \
  && ARCH=$(dpkg --print-architecture) \
  && curl -fsSL "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-${ARCH}.deb" -o /tmp/quarto.deb \
  && dpkg -i /tmp/quarto.deb \
  && rm /tmp/quarto.deb \
  && rm -rf /var/lib/apt/lists/*

RUN Rscript -e "\
  install.packages(c('ashr', 'circlize', 'here', 'matrixStats', 'msigdbr')); \
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

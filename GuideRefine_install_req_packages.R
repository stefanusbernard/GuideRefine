# Install and update all required packages for the GuideRefine pipeline

# BiocManager is required to install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Note: BSgenome.Hsapiens.NCBI.T2TCHM13v2.0 (needed only for the optional
# T2T-CHM13 mode) is not on Bioconductor — it must be forged locally.
# See genome_annotation/README.md.
list_packages_bioconductor <- c(
  "crisprBase",
  "crisprBowtie",
  "Rbowtie",
  "BSgenome.Hsapiens.UCSC.hg38",
  "GenomicRanges",
  "GenomeInfoDb"
)

for (pkg in list_packages_bioconductor) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = TRUE)
  }
}

list_packages_cran <- c(
  "rmarkdown",
  "stringr",
  "tidyverse",
  "knitr",
  "openxlsx",
  "pandoc",
  "data.table"  # only needed for the optional T2T-CHM13 annotation script
)

new_packages <- list_packages_cran[!(list_packages_cran %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)

# Ensure pandoc (system tool required by rmarkdown) is available
if (!rmarkdown::pandoc_available()) {
  message("Pandoc not found. Installing via the pandoc R package...")
  pandoc::pandoc_install()
  pandoc::pandoc_activate()
}

message("Pandoc version: ", rmarkdown::pandoc_version())

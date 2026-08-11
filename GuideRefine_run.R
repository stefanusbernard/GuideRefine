library(rmarkdown)

if (!rmarkdown::pandoc_available()) {
  stop(
    "Pandoc >= 1.12.3 is required but was not found.\n",
    "Run GuideRefine_install_req_packages.R to install it automatically."
  )
}

my_params <- list(
  output_filename = "my_library",
  library_dir     = "./public_crispr_library/processed/",
  sgrna_library   = "my_library",
  terms = c("CONTROL", "Control", "control", "INTRON", "NO_SITE"),

  ref_bsgenome             = "BSgenome.Hsapiens.UCSC.hg38",
  ccds_filename_dir        = "genome_annotation/CCDS.20221027.txt",
  hgnc_link                = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
  guide_aln_granges_genome = "hg38",
  aln_file_dir             = "./object_intermediate/hg38/",
  check_index_file         = "hg38.1.ebwt",
  bowtie_index_dir         = "bowtie_index/hg38/",
  name_prefix_for_index    = "hg38",
  fasta_file               = "",

  remove_all_single_mismatch = TRUE,
  output_dir = "./output_cleaning/"
)

# T2T-CHM13: build the BSgenome package and CCDS annotation first (see
# genome_annotation/README.md), then swap in this block instead.
#
# my_params <- list(
#   output_filename = "my_library",
#   library_dir     = "./public_crispr_library/processed/",
#   sgrna_library   = "my_library",
#   terms = c("CONTROL", "Control", "control", "INTRON", "NO_SITE"),
#   ref_bsgenome             = "BSgenome.Hsapiens.NCBI.T2TCHM13v2.0",
#   ccds_filename_dir        = "genome_annotation/T2T-CHM13v2.0_gene_annot_granges.rds",
#   hgnc_link                = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
#   guide_aln_granges_genome = "T2T-CHM13v2.0",
#   aln_file_dir             = "./object_intermediate/T2T-CHM13_bsgenome/",
#   check_index_file         = "T2T-CHM13_BSgenome.1.ebwt",
#   bowtie_index_dir         = "bowtie_index/T2T-CHM13_bsgenome/",
#   name_prefix_for_index    = "T2T-CHM13_BSgenome",
#   fasta_file               = "",
#   remove_all_single_mismatch = TRUE,
#   output_dir = "./output_cleaning/"
# )


if (!dir.exists(my_params$output_dir)) {
  dir.create(my_params$output_dir, recursive = TRUE)
}

date_tag <- format(Sys.Date(), "%Y%m%d")
output_file <- paste0(my_params$output_dir, my_params$output_filename, "_", date_tag, ".html")

rmarkdown::render(
  input = "GuideRefine.Rmd",
  output_file = output_file,
  params = my_params,
  envir = new.env()
)

message("Pipeline rendered successfully: ", output_file)

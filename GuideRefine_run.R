library(rmarkdown)

if (!rmarkdown::pandoc_available()) {
  stop(
    "Pandoc >= 1.12.3 is required but was not found.\n",
    "Run GuideRefine_install_req_packages.R to install it automatically."
  )
}

# Define parameters
my_params <- list(
  output_filename = "avana_library",
  
  library_dir = "./public_crispr_library/",
  sgrna_library = "avana_library",
  terms = c("CONTROL", "Control", "control", "INTRON", "Intron", "intron", "NO_SITE", "ONE_SITE_INTERGENIC", "LacZ", "luciferase", "EGFR"),
  
  # reference genome and hgnc annotation
  ref_bsgenome = "BSgenome.Hsapiens.NCBI.T2TCHM13v2.0",
  ccds_filename_dir = "annotation_file/T2T-CHM13v2.0_gene_annot_granges.rds",
  hgnc_link = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
   
  # alignment file location and criteria
  guide_aln_granges_genome = "T2T-CHM13v2.0",
  aln_file_dir = "./object_intermediate/May2026_T2T-CHM13_bsgenome/",
  
  # turn this on for T2T-CHM13
  check_index_file = "T2T-CHM13_BSgenome.1.ebwt",
  bowtie_index_dir = "bowtie_index/T2T-CHM13_bsgenome/",
  name_prefix_for_index = "T2T-CHM13_BSgenome",
  
  # Default parameter is to remove all sgRNA aligned with single mismatch at any 20 nucleotides spacer
  remove_all_single_mismatch = TRUE,
  
  # output dir
  output_dir = "./output_cleaning/May2026_T2T-CHM13/"
)


if (!dir.exists(my_params$output_dir)) {
  dir.create(my_params$output_dir, recursive = TRUE)
}


date_tag <- format(Sys.Date(), "%Y%m%d")
output_file <- paste0(my_params$output_dir, my_params$output_filename, "_", date_tag, ".html")


# Render pipeline
rmarkdown::render(
  input = "GuideRefine.Rmd",
  output_file = output_file,
  params = my_params,
  envir = new.env()
)

message("Pipeline rendered successfully: ", output_file)

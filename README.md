# GuideRefine

A pipeline to remove problematic sgRNAs from genome-wide CRISPR-Cas9 knockout libraries before running a screen analysis (e.g. MAGeCK).

GuideRefine aligns every guide to the reference genome and removes guides that:
- hit more than one genomic locus (**multi-target**)
- align with a single mismatch anywhere in the 20 nt spacer (**off-target**, optional)
- align with a single or double mismatch at the PAM-distal region (**off-target**)

Genes left with fewer than 3 valid guides after filtering are also removed.

Adapted from [DeKegel & Ryan, 2019](https://pubmed.ncbi.nlm.nih.gov/31652272/).

---

## Requirements

- R ≥ 4.2 and Bioconductor ≥ 3.16
- [Bowtie](https://bowtie-bio.sourceforge.net/) ≥ 1.3 installed and on your `PATH`

Install R packages:

```r
install.packages(c("tidyverse", "stringr", "directlabels", "knitr", "xlsx"))

if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("crisprBase", "crisprBowtie", "GenomicRanges",
                       "GenomeInfoDb", "BSgenome.Hsapiens.UCSC.hg38"))
```

---

## How to use

### Step 1 — Prepare your input library

Your sgRNA library must be a **TSV file with no header** and three columns:

```
sgRNA_name    spacer_sequence    gene_symbol
```

Place it in the `public_crispr_library/` folder. Several public libraries (Brunello, Avana, GeckoV2, TKOv3) are already included there.

### Step 2 — Download the CCDS annotation (one-time setup)

Download the latest human CCDS file and save it to `annotation_file/`:

```
https://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS.current.txt
```

A copy from 2022 (`CCDS.20221027.txt`) is already included and works out of the box.

### Step 3 — Configure `GuideRefine_run.R`

Open [GuideRefine_run.R](GuideRefine_run.R) and update the highlighted fields:

```r
my_params <- list(

  output_filename = "my_library",          # <-- name for your output files
  library_dir     = "./public_crispr_library/",
  sgrna_library   = "my_library",          # <-- TSV filename (without .tsv)

  # Control/non-targeting guide identifiers in your gene column
  terms = c("CONTROL", "Control", "control", "INTRON", "NO_SITE"),

  # Reference genome (hg38 — default)
  ref_bsgenome             = "BSgenome.Hsapiens.UCSC.hg38",
  ccds_filename_dir        = "annotation_file/CCDS.20221027.txt",
  guide_aln_granges_genome = "hg38",
  aln_file_dir             = "./object_intermediate/hg38/",
  check_index_file         = "hg38.1.ebwt",
  bowtie_index_dir         = "bowtie_index/hg38/",
  name_prefix_for_index    = "hg38",
  fasta_file               = "",  # only needed if rebuilding the Bowtie index from scratch

  # Filtering strictness:
  # TRUE  = remove guides with a single mismatch anywhere in the spacer (recommended)
  # FALSE = remove only guides with mismatches at the PAM-distal region
  remove_all_single_mismatch = TRUE,

  output_dir = "./output_cleaning/"
)
```

> **T2T-CHM13 genome:** If you want to align against T2T-CHM13 instead of hg38, run the two notebooks in `annotation_file/` first to build the BSgenome package and CCDS annotation, then swap in the T2T parameter block shown in the comments inside `GuideRefine_run.R`.

### Step 4 — Run the pipeline

```r
source("GuideRefine_run.R")
```

Or from the terminal:

```bash
Rscript GuideRefine_run.R
```

The first run aligns all guides against the genome and caches the result in `object_intermediate/`. Subsequent runs on the same library skip alignment automatically. To force re-alignment, delete the `*_aln.csv` file for that library.

---

## Output files

Results are written to `output_cleaning/`:

| File | Description |
|------|-------------|
| `<name>_<date>.html` | HTML report with filtering statistics and QC plots |
| `<name>_<N>K_refined.tsv` | Cleaned sgRNA library ready for MAGeCK |
| `<name>_disposed_sgRNAs.tsv` | Removed guides with the reason for each removal |
| `<name>_full_report.xlsx` | Per-gene guide summary in Excel |

---

## Repository structure

```
GuideRefine/
├── GuideRefine_run.R          # Start here — configure and run the pipeline
├── GuideRefine.Rmd            # Pipeline logic (called automatically by _run.R)
├── GuideRefine_functions.R    # Helper functions
│
├── annotation_file/           # CCDS and genome annotation files
├── bowtie_index/              # Pre-built Bowtie index (hg38 and T2T-CHM13)
├── public_crispr_library/     # Input sgRNA library TSV files
├── object_intermediate/       # Cached alignment files (auto-generated)
└── output_cleaning/           # Pipeline outputs (auto-generated)
```

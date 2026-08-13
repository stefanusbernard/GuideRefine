# GuideRefine

A pipeline to computationally remove problematic sgRNAs from genome-wide CRISPR-Cas9 knockout libraries before running a screen analysis (e.g. MAGeCK).

GuideRefine aligns every guide to the reference genome and removes guides that:
- hit more than one genomic locus (**multi-target**)
- align with a single mismatch anywhere in the spacer (**off-target**, on by default, can be switched off)
- align with a single or double mismatch at the PAM-distal region (**off-target**)

Genes left with fewer than 3 valid guides after filtering are also excluded.

Guide-filtering criteria adapted from [DeKegel & Ryan, 2019](https://pubmed.ncbi.nlm.nih.gov/31652272/)

---

## Dependencies

- R ≥ 4.2 and Bioconductor ≥ 3.16 (developed and tested on R v4.5.2)
- [Bowtie](https://bowtie-bio.sourceforge.net/) ≥ 1.3 on your `PATH`

## Install

```r
source("GuideRefine_install_req_packages.R")
```

One-time step. Installs all required CRAN and Bioconductor packages, plus Pandoc if missing.

## Data to prepare

- **sgRNA library** — TSV, no header, 3 columns (`sgRNA_name  spacer_sequence  gene_symbol`). See [example_guiderefine_run/README.md](example_guiderefine_run/README.md) for a worked example (Brunello) and how to reformat your own library into this format.
- **CCDS annotation** — `CCDS.20221027.txt` is already in `genome_annotation/` and works out of the box. For T2T-CHM13 instead of hg38, see [genome_annotation/README.md](genome_annotation/README.md).
- **Bowtie index** — built automatically on first run: set `fasta_file` in `GuideRefine_run.R` to a genome FASTA, and the pipeline builds the index into `bowtie_index_dir` if it isn't already there. hg38 FASTA: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz. For T2T-CHM13, see [genome_annotation/README.md](genome_annotation/README.md).

## How to use

1. Open [GuideRefine_run.R](GuideRefine_run.R) and set `sgrna_library` to your library's filename (without `.tsv`).
2. Run it:
   ```r
   source("GuideRefine_run.R")
   ```
   or `Rscript GuideRefine_run.R` from a terminal.

The first run aligns all guides and caches the result in `object_intermediate/`. Later runs on the same library skip alignment automatically — delete the `*_aln.csv` file to force re-alignment.

## Output

Written to `output_cleaning/`:

| File | Description |
|------|-------------|
| `<name>_<date>.html` | Report with filtering statistics and QC plots |
| `<name>_<N>K_refined.tsv` | Cleaned sgRNA library ready for MAGeCK |
| `<name>_disposed_sgRNAs.tsv` | Removed guides with the reason for each removal |
| `<name>_full_report.xlsx` | Per-gene guide summary |

---

## Repository structure

```
GuideRefine/
├── GuideRefine_install_req_packages.R  # run once to install dependencies
├── GuideRefine_run.R                   # configure and run the pipeline
├── GuideRefine.Rmd                     # pipeline logic
├── GuideRefine_functions.R             # helper functions
│
├── genome_annotation/         # CCDS / genome annotation
├── example_guiderefine_run/ # worked example: raw + processed Brunello library, reformat script, sample output
│   ├── raw/                 # original download
│   ├── processed/           # GuideRefine-ready TSV
│   └── output/              # sample GuideRefine output
│
├── object_intermediate/     # cached alignments (auto-generated)
└── output_cleaning/         # pipeline outputs (auto-generated)
```

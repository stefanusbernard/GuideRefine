# Gene annotation data

Gene/exon annotation used to map sgRNA alignments onto genes. Most files here are large NCBI downloads or locally-built packages, so they're git-ignored — this explains how to get them.

## hg38 (default)

Download the CCDS file into this folder:

```
https://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS.current.txt
```

That's it — `GuideRefine.Rmd` reads `.txt` CCDS files directly.

## T2T-CHM13v2.0 (optional)

Run the two notebooks in this folder, in order.

**1. `forge_T2T-CHM13_bsgenome.Rmd`** — open and run in RStudio (it uses `rstudioapi` to find its own path, so it can't be `Rscript`'d). It downloads the T2T-CHM13v2.0 genome from NCBI (accession `GCF_009914755.1`, ~3 GB) and forges + installs the `BSgenome.Hsapiens.NCBI.T2TCHM13v2.0` package. This package isn't on Bioconductor, so this is the only way to get it.

**2. `T2T-CHM13_gff_to_ccds_conversion.Rmd`** — produces the T2T equivalent of the CCDS file (`T2T-CHM13v2.0_gene_annot_granges.rds`). First download, for accession `GCF_009914755.1`, from the [NCBI Genome page](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/):
- the RefSeq GFF3 annotation → save as `GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz`
- the assembly report → strip the leading `#`-comment lines and save the remaining tab-separated table as `T2T-CHM13v2_assembly_report.txt`

Run the notebook with this folder as the working directory. It needs `zgrep`/`grep` on `PATH` (Git Bash/WSL on Windows).

The resulting `.rds` is what `ccds_filename_dir` should point to for T2T runs (see the commented T2T block in `GuideRefine_run.R`).

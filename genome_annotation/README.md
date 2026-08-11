# Gene annotation data

Gene/exon annotation used to map sgRNA alignments onto genes. Most files here are large NCBI downloads or locally-built packages, so they're git-ignored — this explains how to get them.

## hg38 (default)

Download the CCDS file into this folder:

```
https://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS.current.txt
```

That's it — `GuideRefine.Rmd` reads `.txt` CCDS files directly.

## T2T-CHM13v2.0 (optional)

Two independent notebooks — run in either order, both produce different outputs GuideRefine needs.

**`forge_T2T-CHM13_bsgenome.Rmd`** → installs the `BSgenome.Hsapiens.NCBI.T2TCHM13v2.0` package (`ref_bsgenome`). Run in RStudio (uses `rstudioapi`, so it can't be `Rscript`'d). Downloads the genome from NCBI (`GCF_009914755.1`, ~3 GB) and forges + installs the package — not on Bioconductor otherwise.

**`T2T-CHM13_gff_to_ccds_conversion.Rmd`** → produces `T2T-CHM13v2.0_gene_annot_granges.rds` (`ccds_filename_dir`). First download, for accession `GCF_009914755.1`, from the [NCBI Genome page](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/):
- the RefSeq GFF3 annotation → save as `GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz`
- the assembly report → strip the leading `#`-comment lines and save as `T2T-CHM13v2_assembly_report.txt`

Run with this folder as the working directory. Needs `zgrep`/`grep` on `PATH` (Git Bash/WSL on Windows). Both outputs are required for T2T runs (see the commented T2T block in `GuideRefine_run.R`).

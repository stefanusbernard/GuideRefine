# GuideRefine: a tool to capture and remove multi-target-/off-target sgRNAs from CRISPR-Cas9 KO libraries
a tool (R Markdown) to clean the CRISPR sgRNA library from:
- suspected multi-target sgRNAs
- suspected off-target sgRNAs:
  - sgRNA align with single mismatch at any location in 20-spacer nucleotides (optional)
  - sgRNA align with single mismatch at PAM-distal site (position 1, 2 in + strand; and position 19, 20 in - strand)
  - sgRNA align with double mismatch at PAM-distal site (position 1 & 2 in + strand; and position 19 & 20 in - strand)

This tool is an adapted version of script developed by Barbara & Ryan, 2019
https://pubmed.ncbi.nlm.nih.gov/31652272/

Input file:
1.  sgRNA library TSV file input for MAGeCK (Format: sgRNA, spacer, gene; without header)
2.  Human genome GRCh38 indexed in bowtie
3.  CCDS annotation; obtained from https://ftp.ncbi.nlm.nih.gov/pub/CCDS/ (pick current_human/CCDS.current.txt, depending on the updated version, the one used for this analysis was released in 12th January 2025)
4.  HGNC dataset; obtained from https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt

*for latest HGNC dataset you can change the link in line 711

`hgnc <- read_tsv(url('https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt'))`

Output file:
1. Refined sgRNA library formatted according to MAGeCK input (sgRNA, spacer, gene; without header)
2. HTML report from R Markdown
3. Excel report containing sgRNAs targeting each genes

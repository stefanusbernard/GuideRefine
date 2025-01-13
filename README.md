# CRISPR sgRNA Library Refining tool
a tool (R Markdown) developed to clean the CRISPR sgRNA library from:
- suspected multi-target sgRNAs
- off-target sgRNAs (align with single mismatch or double mismatch at the PAM distal site)

This tool is an adapted version of script developed by Barbara & Ryan, 2019
https://pubmed.ncbi.nlm.nih.gov/31652272/

Input file:
1.  sgRNA library TSV file input for MAGeCK (Format: sgRNA, spacer, gene; without header)
2.  Human genome GRCh38 indexed in bowtie
3.  CCDS annotation; obtained from https://ftp.ncbi.nlm.nih.gov/pub/CCDS/ (pick current_human/CCDS.current.txt, depending on the updated version, the one used for this analysis was released in 12th January 2025)
4.  HGNC dataset; obtained from https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt

Output file:
1. Refined sgRNA library
2. HTML report from R Markdown
3. Excel report containing sgRNAs targeting each genes

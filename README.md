# CRISPR sgRNA Library Refining script
a tool (R Markdown) developed to clean the CRISPR sgRNA library from:
- suspected multi-target sgRNAs
- off-target sgRNAs (align with single mismatch or double mismatch at the PAM distal site)

Input file:
1.  sgRNA library TSV file input for MAGeCK (Format: sgRNA, spacer, gene; without header)
2.  Human genome GRCh38 indexed in bowtie
3.  CCDS annotation
4.  HGNC dataset

Output file:
1. Refined sgRNA library
2. HTML report from R Markdown
3. Excel report containing sgRNAs targeting each genes

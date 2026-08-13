# Example GuideRefine run

A worked example of a complete GuideRefine run, from raw public library download to
final cleaned output, using the [Brunello](https://pmc.ncbi.nlm.nih.gov/articles/PMC4744125/)
sgRNA library. Use this as a template for preparing your own library and for seeing
what GuideRefine's output looks like before you run the pipeline yourself.

```
example_guiderefine_run/
├── raw/         # broadgpp-brunello-library-contents.txt — original download, unmodified
├── processed/   # broadgpp-brunello-library-contents.tsv — reformatted to GuideRefine's required input
├── reformat_sgRNA_library.rmd   # script used to turn raw/ into processed/
└── output/      # GuideRefine's output from running on processed/broadgpp-brunello-library-contents.tsv
```

## Required input format

GuideRefine needs a **single processed library file**: a TSV with **no header** and
exactly 3 columns, in this order:

| Column | Description |
|---|---|
| `sgRNA` | Guide name/ID |
| `spacer` | The 19-20 nt guide sequence |
| `gene` | Target gene symbol |

`head()` of [`processed/broadgpp-brunello-library-contents.tsv`](processed/broadgpp-brunello-library-contents.tsv):

```
sgA1BG_1	CATCTTCTTTCACCTGAACG	A1BG
sgA1BG_2	CTCCGGGGAGAACTCCGGCG	A1BG
sgA1BG_3	TCTCCATGGTGCATCAGCAC	A1BG
sgA1BG_4	TGGAAGTCCACTCCACTCAG	A1BG
sgA1CF_1	ACAGGAAGAATTCAGTTATG	A1CF
```

## Getting your own library into this format

- **Use `reformat_sgRNA_library.rmd`**: drop your raw library file (`.txt`/`.csv`/`.xlsx`,
  one gene per row, with a gene-symbol column and a 20 nt sequence column) into a `raw/`
  folder, set `filename` in the "manual, reformat single selected library" chunk, and run
  it. It auto-detects the gene-symbol and sequence columns via HGNC matching and column
  length, and writes a 3-column TSV to `processed/`. Rows with a missing gene symbol are
  dropped automatically.
- **Or write your own script/reformat by hand** — as long as the result is a headerless,
  3-column TSV (`sgRNA`, `spacer`, `gene`)

Once you have your processed TSV, point `GuideRefine_run.R` at it:
`library_dir = "./example_guiderefine_run/processed/"`, `sgrna_library = "<filename without .tsv>"`.

## Output

The `output/` folder shows what GuideRefine produced for this Brunello run:

| File | Description |
|------|-------------|
| `broadgpp-brunello-library-contents_20260730.html` | Report with filtering statistics and QC plots |
| `broadgpp-brunello-library-contents_62K_refined.tsv` | Cleaned sgRNA library ready for MAGeCK |
| `broadgpp-brunello-library-contents_disposed_sgRNAs.tsv` | Removed guides with the reason for each removal |
| `broadgpp-brunello-library-contents_full_report.xlsx` | Per-gene guide summary |

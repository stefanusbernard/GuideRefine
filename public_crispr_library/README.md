# Public CRISPR sgRNA libraries

- `raw/` — original downloads, unmodified. See [raw/README.md](raw/README.md) for sources.
- `processed/` — GuideRefine-ready libraries: 3-column TSV, no header (`sgRNA`, `spacer`, `gene`).

| Library | Raw file | Processed file | Script |
|---|---|---|---|
| Avana | `raw/avana_library.txt` | `processed/avana_library.tsv` | `reformat_sgRNA_library.rmd` |
| Brunello | `raw/broadgpp-brunello-library-contents.txt` | `processed/broadgpp-brunello-library-contents.tsv` | `reformat_sgRNA_library.rmd` |
| TKOv3 | `raw/tkov3_guide_sequence.xlsx` | `processed/tkov3_guide_sequence.tsv` | `reformat_sgRNA_library.rmd` |
| Yusa / Project Score | `raw/yusa_hcrispr_ko_grnas.xlsx` | `processed/yusa_hcrispr_ko_grnas.tsv` | `reformat_sgRNA_library.rmd` |
| Jacquère | `raw/Jacquere_PerGuideAnnotations_Quota4.csv` | `processed/Jacquere_PerGuideAnnotations_Quota4.tsv` | `data_processing_jacquere.Rmd` |

## Reprocessing / adding a library

1. Drop the raw file into `raw/`.
2. Reformat it into `processed/` as a 3-column TSV (`sgRNA`, `spacer`, `gene`, no header):
   - **Standard format** (a gene-symbol column + a 20 nt sequence column, in `.txt`/`.csv`/`.xlsx`, one gene per row): use `reformat_sgRNA_library.rmd`. Either set `filename` in the "manual, reformat single selected library" chunk, or just drop the file in `raw/` and run the "automatic pipeline" chunk, which reformats every file in `raw/` it doesn't already recognize. It drops rows with a missing gene symbol automatically.
   - **Non-standard format** (e.g. multi-gene-per-row fields, merging multiple sources): write a dedicated script — see `data_processing_jacquere.Rmd` for an example (resolves `GENE1|GENE2`-style fields).
3. Point `GuideRefine_run.R` at it: `library_dir = "./public_crispr_library/processed/"`, `sgrna_library = "<filename without .tsv>"`.

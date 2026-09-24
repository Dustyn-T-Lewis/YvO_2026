# 01 · Normalization

The raw intensity workbook to a filtered, outlier-screened, cyclic-loess normalized DAList.

```
00_input/YvO_raw.xlsx + YvO_meta.xlsx + HPA_skeletal_muscle_annotations.tsv + YvO_pheno_calc.xlsx
  01_normalize.R             HPA filter -> blood removal -> dedup by uniprot_id ->
                             missingness -> 4-method outlier consensus -> missingness again ->
                             normalize_data(cycloess, span 0.7)
                             -> c_data/01_normalization.xlsx, 02_normalized.csv,
                                03_DAList_normalized.rds, 00_report_intermediates.rds
  02_generate_reports.R      00_report_intermediates.rds -> b_reports/04_diagnostics.pdf
  03_fraction_composition.R  -> c_data/03_fraction_composition.csv, 03_fraction_top_carriers.csv
```

```sh
Rscript 01_normalization/a_script/01_normalize.R
Rscript 01_normalization/a_script/02_generate_reports.R
Rscript 01_normalization/a_script/03_fraction_composition.R
```

## What comes out

`01_normalization.xlsx`: 9 sheets, ships as S8 Table.

`02_normalized.csv`: the text serialisation that stage 02 and `03_DEP/01_run_dep.R` both
read in place of the RDS, so missForest splits and the limma fit do not shift with a
binary-double round trip.

`03_DAList_normalized.rds`: read by stages 02 and 03, `03_DEP/a_script/supp/05` and 04_Figures.

`00_report_intermediates.rds`: read by `02_generate_reports.R`, by
`03_DEP/a_script/02_generate_reports.R` for the 64-against-62 sensitivity table, and by
`F00/a_script/01_supp_panels.R`.

`03_fraction_*.csv` are terminal; no script reads them.

## Orderings that matter

`filter_proteins_by_group(min_reps = 10)` runs twice, once on all samples and again after
`filter_samples()` drops the consensus outliers. Proteins kept on the strength of
observations in a sample later excluded would otherwise survive into the fit
(`01_normalize.R:274-285`).

`align_data_and_metadata()` is called with `prefer_group_blocks = FALSE` and its output
asserted identical to the input column order. The default reorders samples into group
blocks, which would change the matrix every later stage sees.

`02_qc_pre.pdf` is written before `normalize_data()`, `03_qc_post.pdf` after.

## Cost

14.2 s, 2.4 s, 4.8 s.

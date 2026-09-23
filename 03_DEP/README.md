# 03 · Differential abundance

The normalized matrix to four limma contrasts, and the supplement layer built on top.

```
01_normalization/c_data/02_normalized.csv + 03_DAList_normalized.rds
  a_script/01_run_dep.R          ~ 0 + group, duplicateCorrelation on subject, 4 contrasts,
                                 Pi-score -> 01_limma_DAList.rds, 03_combined_results.csv,
                                 03_DEP_results.xlsx (6 sheets)
  a_script/02_generate_reports.R + outlier_sensitivity
  supp/01_effect_size_robustness.R  + blunting, bootstrap_ci, power_analysis,
                                      imputation_sensitivity
  supp/02_supplement_covariate.R    + supplement_{by_age,design,forced_fit,forced_shift,nore},
                                      Overview
  supp/03_reversal_aging_fdr.R      -> 03_reversal_aging_fdr.csv, 03_reversal_null_draws.csv
  supp/04_discordant_ora.R          -> 04_discordant_ora.csv
  supp/05_supplement_sensitivity.R  -> 05_supplement_permanova.csv
```

```sh
YVO_PIPELINE_RUN=1 Rscript 03_DEP/a_script/01_run_dep.R
Rscript 03_DEP/a_script/02_generate_reports.R
Rscript 03_DEP/supp/01_effect_size_robustness.R
Rscript 03_DEP/supp/02_supplement_covariate.R
Rscript 03_DEP/supp/03_reversal_aging_fdr.R
Rscript 03_DEP/supp/04_discordant_ora.R
Rscript 03_DEP/supp/05_supplement_sensitivity.R
```

## What comes out

`03_DEP_results.xlsx`: 17 sheets, ships as S10 Table, built in three passes: 6 by
`01_run_dep.R`, 1 by `02_generate_reports.R`, 4 by `supp/01`, 5 plus Overview by `supp/02`.

`03_combined_results.csv`: one row per protein with `logFC_`, `P.Value_`, `adj.P.Val_`,
`pi_score_` and `sig_pi_` per contrast. Read by `supp/03`, `supp/04` and most of 04_Figures.

The matrix of record, 2,106 proteins across 62 samples. Contrast order is
Aging / Training_Young / Training_Old / Interaction.

| Quantity | Value | Written by |
|---|---|---|
| DEPs at FDR < 0.05 | 278 / 135 / 0 / 1 | `a_script/01_run_dep.R` |
| DEPs at Pi < 0.05 | 195 / 99 / 18 / 33 | `a_script/01_run_dep.R` |
| DEPs at p < 0.05 | 570 / 448 / 196 / 158 | `a_script/01_run_dep.R` |
| Smallest adjusted p, Training_Old | 0.156, tied HPRT1 and FSCN1 | `a_script/01_run_dep.R` |
| Within-subject correlation | 0.278 | `duplicateCorrelation` |
| Pi union across contrasts | 294, of which Aging 174 exclusive | F02 panel D |
| Aging and Training_Young at FDR | 19 shared, 6 same-sign | F04 panel D |

Check a figure caption or a manuscript number against this table, not against
another paragraph.

## Orderings that matter

`01_run_dep.R` rebuilds the workbook from scratch and stops if any of the eleven
later-appended sheets is present, unless `YVO_PIPELINE_RUN` is set. `run_all.R` sets it
because the scripts that rewrite those sheets are about to run.

`02_generate_reports.R`, `supp/01` and `supp/02` each `loadWorkbook()` and `saveWorkbook()`
in turn, so they run in that order. `supp/02` writes Overview last, so the index matches
the sheets.

`supp/03` and `supp/05` require `02_imputation/c_data/01_DAList_imputed.rds`. `supp/01`
skips its `imputation_sensitivity` sheet with a message when that file is absent.

## Cost

`01_run_dep.R` 12.5 s. `supp/01` 160.0 s, which is 10,000 BCa bootstrap replicates per
contrast plus a second full limma fit on the imputed matrix. Everything else is under 20 s.

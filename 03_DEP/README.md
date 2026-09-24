# 03_DEP

Fits limma to the normalized matrix for the four contrasts, then builds the supplementary analyses on that fit.

## Reads

- `01_normalization/c_data/02_normalized.csv` and `03_DAList_normalized.rds`
- `02_imputation/c_data/01_DAList_imputed.rds`, for `supp/01`, `03` and `05`
- `00_input/YvO_meta.xlsx`, for `supp/02`, `03` and `05`

## Writes

- `c_data/03_DEP_results.xlsx`: S10 Table, 17 sheets. `01_run_dep.R` writes 6, `02_generate_reports.R` adds `outlier_sensitivity`, `supp/01` adds 4, and `supp/02` adds 5 and the Overview.
- `c_data/03_combined_results.csv`: one row per protein with `logFC_`, `P.Value_`, `adj.P.Val_`, `pi_score_` and `sig_pi_` for each contrast. Most figures read it.
- `c_data/01_limma_DAList.rds`: the fitted DAList
- `c_data/02_supplement_*.csv`: whether supplement arm can enter the design, and what forcing it in does
- `c_data/03_reversal_aging_fdr.csv`, `03_reversal_null_draws.csv`: whether training in older adults moves the aging proteins back toward young, with its permutation null
- `c_data/04_discordant_ora.csv`: pathway enrichment of the proteins whose training response differs in sign between age groups
- `c_data/05_supplement_permanova.csv`: proteome variance on the supplement stratum in older adults
- `b_reports/S9.pdf`, `S9.png`: S9 Figure, from `supp/02`
- `b_reports/`: proteoDA reports and contrast summaries, not tracked

The fit uses `~ 0 + group` with `duplicateCorrelation` on subject (0.278). The contrasts are Aging (Old_Pre minus Young_Pre), Training_Young, Training_Old and Interaction (Training_Old minus Training_Young). At FDR < 0.05 they give 278, 135, 0 and 1 proteins; at Π < 0.05, 195, 99, 18 and 33. The smallest Training_Old adjusted p is 0.155657, tied between FSCN1 and HPRT1.

## Run

```sh
YVO_PIPELINE_RUN=1 Rscript 03_DEP/a_script/01_run_dep.R
Rscript 03_DEP/a_script/02_generate_reports.R
Rscript 03_DEP/a_script/supp/01_effect_size_robustness.R
Rscript 03_DEP/a_script/supp/02_supplement_covariate.R
Rscript 03_DEP/a_script/supp/03_reversal_aging_fdr.R
Rscript 03_DEP/a_script/supp/04_discordant_ora.R
Rscript 03_DEP/a_script/supp/05_supplement_sensitivity.R
```

## Order

Stages 01 and 02 run first. `01_run_dep.R` rebuilds the workbook from scratch and stops if a later-added sheet is already there, unless `YVO_PIPELINE_RUN` is set, because that sheet would be lost. `02_generate_reports.R`, `supp/01` and `supp/02` each add sheets to the workbook, so they run in that order. `supp/03` and `supp/04` read `03_combined_results.csv`, so they run after `01_run_dep.R`.

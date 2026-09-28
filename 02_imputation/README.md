# 02_imputation

Imputes the normalized matrix with missForest, the benchmark's top-ranked method.

## Reads

- `01_normalization/c_data/02_normalized.csv` and `03_DAList_normalized.rds`
- `c_data/benchmark/04_composite_ranking.csv`, written by the benchmark

## Writes

- `c_data/02_imputation.xlsx`: S9 Table, 7 sheets
- `c_data/01_imputed.csv`: the imputed matrix, read by F05 and F06
- `c_data/01_DAList_imputed.rds`: the imputed DAList with each protein's missingness class, read by `03_DEP/a_script/supp/01`, `03`, `05`, F02, F04 and F05
- `c_data/02_mar_mnar_classification.csv`: MAR, MNAR or complete for each protein, read by the benchmark and F04
- `c_data/00_report_intermediates.rds`: read by F00
- `c_data/benchmark/01_reconstruction.csv` to `04_composite_ranking.csv` and `04_full_report.txt`: the benchmark results
- `b_reports/`: missingness and imputation reports, not tracked

## Run

```sh
Rscript 02_imputation/a_script/01_impute.R
Rscript 02_imputation/a_script/02_generate_reports.R
```

The benchmark of 17 methods is opt-in and not in `run_all.R`; its ranking is tracked, so a normal run does not need it.

```sh
Rscript 02_imputation/a_script/benchmark/_run_all.R
SKIP_IMPUTE=1 Rscript 02_imputation/a_script/benchmark/_run_all.R   # reuse imputed_matrices.rds
```

GSimp is not on CRAN or Bioconductor, so `a_script/benchmark/methods/gsimp_source/` vendors it unmodified from https://github.com/WandeRum/GSimp (Wei et al. 2018, PLOS Computational Biology 14:e1005973) under GPL-3; `methods/09_gsimp.R` sources it.

## Order

The first benchmark run needs `02_mar_mnar_classification.csv` from `01_impute.R`. Both stage scripts need the benchmark ranking; `run_all.R` stops if it is missing.

# 02 · Imputation

The normalized matrix with its gaps filled by missForest, plus the benchmark that chose it.

```
01_normalization/c_data/02_normalized.csv + 03_DAList_normalized.rds
  benchmark/_run_all.R   17 methods -> c_data/benchmark/*.csv, imputed_matrices.rds  (opt-in)
  01_impute.R            k-means on (mean intensity, % missing) -> MAR/MNAR/Complete;
                         Fisher test per MNAR protein; missForest(maxiter 10, ntree 100)
                         -> c_data/02_imputation.xlsx, 01_imputed.csv,
                            02_mar_mnar_classification.csv, 01_DAList_imputed.rds
  02_generate_reports.R  -> b_reports/01_missingness_report.pdf, 02_imputation_report.pdf
```

```sh
Rscript 02_imputation/a_script/benchmark/_run_all.R   # opt-in; SKIP_IMPUTE=1 reuses the cache
Rscript 02_imputation/a_script/01_impute.R
Rscript 02_imputation/a_script/02_generate_reports.R
```

## What comes out

`02_imputation.xlsx`: 7 sheets, ships as S9 Table. Five core sheets, `benchmark_ranking`
when `benchmark/04_composite_ranking.csv` exists, and the Overview index.

`01_DAList_imputed.rds`: the DAList with `miss_classification` and `imputation_reliable`
merged into `$annotation`. Read by `03_DEP/a_script/supp/01,03,05` and by F02, F04, F05, F06.

`02_mar_mnar_classification.csv`: the `km` classifier that `benchmark/_common.R` reads
back as `CLASSIFIERS$km`.

## Orderings that matter

`01_impute.R` must have run once before the benchmark: `benchmark/_common.R` loads
`02_mar_mnar_classification.csv` as one of its two classifiers.

The benchmark must then exist before a full run. `02_generate_reports.R:30` stops on a
missing `04_composite_ranking.csv`, and `run_all.R` preflights the same file because the
benchmark is opt-in and two stages read it.

Rows are sorted with `gene_order` before `missForest` for determinism, so `$annotation` is
re-matched to `rownames($data)` afterwards and asserted identical. `merge()` keeps
left-frame order, which is not the matrix order (`01_impute.R:232-241`).

## Cost

`01_impute.R` 37.9 s, `02_generate_reports.R` 2.4 s. The benchmark is not one of
`run_all.R`'s 22 steps and is run by hand.

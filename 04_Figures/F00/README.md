# F00 · Pipeline QC, supplementary S1a and S1b

Reads the stage 01-03 intermediates and produces two supplementary composites plus the
S1 Table workbook. F00 has no numbered main figure.

```
01_normalization/c_data/00_report_intermediates.rds
02_imputation/c_data/00_report_intermediates.rds + benchmark/04_composite_ranking.csv
03_DEP/c_data/03_DEP_results.xlsx (DA_summary)
  01_supp_panels.R  panels A-N, both composites, the workbook
  90_stitch_F00.R   sources the above
  -> b_reports/supp/{pdf,png}/SUPP_F00_normalization.*   panels A-G
  -> b_reports/supp/{pdf,png}/SUPP_F00_imputation.*      panels H-N
  -> c_data/F00_supplementary.xlsx
```

```sh
Rscript 04_Figures/F00/a_script/90_stitch_F00.R
```

## What comes out

Two 178 x 245 mm composites. A-G covers the filter cascade, per-protein missingness, PCA
before and after cyclic loess, eta-squared retention, outlier consensus and per-sample
missing counts. H-N covers MAR/MNAR classification, the imputation benchmark,
observed-against-imputed density, the MNAR shift audit, sample integrity and DEP counts.

`F00_supplementary.xlsx` is 15 sheets: Overview plus `panel_A` to `panel_N`.

## Orderings that matter

`run_all.R` runs F00 last of the seven figures. It summarises stages 01-03 and reads their
finished intermediates.

`01_supp_panels.R` asserts `02_imputation/c_data/benchmark/04_composite_ranking.csv` before
drawing anything. That file comes from an opt-in script, so panel J would otherwise fail
deep into the run.

`set.seed(42)` is set immediately before panel F, whose `geom_jitter()` draws from the RNG.

## Cost

4.9 s. Nothing is fitted here; every number is read from an upstream artifact.

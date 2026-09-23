# F06 · Figure 6, module discrimination and phenotype coupling

Reads F05's network objects and workbook; writes the three-panel ROC composite, the
supplementary sweep plate and a 13-sheet workbook. Was F07 until 2026-09-22.

```
04_Figures/F05/c_data/F05_supplementary.xlsx + datExpr.rds + module_colors.rds
                                             + me_pre.rds + me_post.rds
02_imputation/c_data/01_imputed.csv, 03_DEP/c_data/03_combined_results.csv
  01_main_panels.R  _supp_module_grid -> _panel_A_auc_bars -> _panel_B_hero_grid
                    -> _supp_panel_B_grid -> _supp_prepare_roc -> _supp_roc_panel
                    -> _supp_multivariate -> _supp_loso_sensitivity -> _supp_loso_wgcna_refit
  02_supp_panels.R  composites the pre-rendered supp PNGs
  90_stitch_F06.R   sources both, sweeps CSVs, copies to Box
  -> c_data/F06_supplementary.xlsx
```

```sh
Rscript 04_Figures/F06/a_script/90_stitch_F06.R
```

## What comes out

`MAIN_F06_composite` with A per-module age ROCs, B per-module training ROCs, C the
module-phenotype hero grid. `SUPP_F06_composite`. A 13-sheet workbook: module grid summary
and curves, multivariate classifier AUC, stability, permutation and curves, classifier
pilot summary and curves, the 180-test panel B screen, and three LOSO sheets.

## Orderings that matter

F05 must run first. Seven scripts stop outright if the F05 workbook is absent —
`_panel_A_auc_bars.R`, `_panel_B_hero_grid.R`, `_supp_prepare_roc.R`,
`_supp_panel_B_grid.R`, `_supp_module_grid.R`, `_supp_multivariate.R` and
`_supp_loso_wgcna_refit.R`. Three of them also read F05's `.rds` objects directly.

`01_main_panels.R` fixes the internal order, each step consuming the previous step's CSVs:
module grid before panel A, hero grid before the panel B sweep, `_supp_prepare_roc` before
`_supp_roc_panel`. `02_supp_panels.R` runs after, reading only rendered PNGs.

`90_stitch_F06.R` ends by removing every remaining `.csv` under `c_data`. The workbook is
the only surviving copy, and it is what `abstract_panels` reads. Verified 2026-09-22: the
workbook rebuilds byte-identically from an empty `c_data`, so a fresh clone is sufficient.

## Cost

255.2 s. `_supp_loso_wgcna_refit.R` refits the network on 30 leave-one-subject-out folds
with no cache.

# F02 · Figure 2, proteome overview and differential expression

Reads the imputed and normalized DALists, the DEP results and the shared fGSEA cache;
produces the six-panel Figure 2, the S3 Figure CV composite and the S3 Table workbook.

```
02_imputation/c_data/01_DAList_imputed.rds + 02_imputation.xlsx
01_normalization/c_data/03_DAList_normalized.rds
03_DEP/c_data/03_combined_results.csv + 03_DEP_results.xlsx
04_Figures/shared/fgsea_tstat_all_v2.csv
  02_supp_panels.R  _supp_A_cv_scatter, _supp_B_cv_violin, _supp_C_imputed, then SUPP composite
  01_main_panels.R  A PCA, B logFC density, C DEP counts, then _panel_D_upset,
                    _panel_E_fgsea, _panel_F_barcode
  90_stitch_F02.R   workbook, cleanup
  -> b_reports/main/{pdf,png}/MAIN_F02_composite.*
  -> c_data/F02_supplementary.xlsx
```

```sh
Rscript 04_Figures/F02/a_script/90_stitch_F02.R
```

## What comes out

A 178 x 115 mm 3 x 2 composite: PCA with bootstrap CIs and a PERMANOVA label, log2FC
distributions with the KS and blunting statistics, DEP percentages at FDR and Pi, the UpSet
overlap from `ComplexHeatmap::make_comb_mat()`, stacked fGSEA counts by database, and the
t-statistic barcode. `F02_supplementary.xlsx` is 14 sheets.

## Orderings that matter

`90_stitch_F02.R` sources `02_supp_panels.R` before `01_main_panels.R`, because the
workbook globs `c_data/` for `audit_|panel_|SUPP_panel_` CSVs and needs the supplementary
ones already written.

That glob makes sheet order follow the filesystem, so `f02_descriptions` is keyed by sheet
name and the stitcher stops on any sheet with no description rather than letting it inherit
a neighbour's.

`_panel_E_fgsea.R` asserts `shared/fgsea_tstat_all_v2.csv` exists. `01_main_panels.R`
sources `shared/build_fgsea_cache.R` first, which rebuilds it when the DEP results are
newer.

## Cost

62.1 s, mostly the 1000-iteration `prcomp()` bootstrap, `adonis2()` and `permutest()` at
999 permutations each, and the 2000-iteration blunting resample. fGSEA is read from cache.

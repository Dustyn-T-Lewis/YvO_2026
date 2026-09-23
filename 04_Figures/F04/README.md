# F04 · Figure 4, training-response concordance between age groups

Reads the combined limma results, the reversal draws and the imputed matrix; writes the
concordance composite, the diagnostics plate, two standalone heatmaps and a 26-sheet
workbook.

```
03_DEP/c_data/03_combined_results.csv + 03_reversal_aging_fdr.csv + 03_reversal_null_draws.csv
02_imputation/c_data/01_DAList_imputed.rds + 02_mar_mnar_classification.csv
00_input/YvO_meta.xlsx
  02_supp_panels.R  13 _supp_* scripts        -> c_data/panel_supp/*.csv
  01_main_panels.R  _panel_A_ORA.R plus shared comparison_panels C/D/E/F
  90_stitch_F04.R   workbook, cleanup, Box copy
  -> c_data/F04_supplementary.xlsx
```

```sh
Rscript 04_Figures/F04/a_script/90_stitch_F04.R
```

## What comes out

`MAIN_F04_composite` with panels A, C, D, E and F, `SUPP_F04_diagnostics` and its captioned
variant, the standalone enrichment and younger-adult DEP heatmaps, and a 26-sheet workbook.

## Orderings that matter

`02_supp_panels.R` runs before `01_main_panels.R`: its CSVs feed both the workbook and the
main composite.

Inside `02_supp_panels.R`, `shared/style.R` is re-sourced after the RRHO2 panel, because
that panel pulls in `print_scale_apply.R`, which multiplies the size globals in place and
does not restore them.

Three supp scripts (`_supp_ora_dedup`, `_supp_threshold_sens`, `_supp_fry_leading`) still
run although their plots left the composite; they write the CSVs behind workbook sheets.

The stitcher calls `cleanup_after_workbook()` and then removes every remaining `.csv` under
`c_data` recursively. Only the workbook survives.

## Cost

96.8 s.

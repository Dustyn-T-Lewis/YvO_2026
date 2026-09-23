# F03 · Figure 3, volcano rings for the four contrasts

Reads the per-contrast DEP tables and the shared fGSEA cache; writes the 2 x 2 volcano
composite, the S4a diagnostics plate, four significance heatmaps and a 17-sheet workbook.

```
03_DEP/c_data/03_DEP_results.xlsx + 03_combined_results.csv
04_Figures/shared/fgsea_tstat_all_v2.csv
  02_supp_panels.R      p / Pi / FDR distributions, MA plots -> c_data/supp/*.csv
  01_main_panels.R      four rings via _build_volcano_panel.R -> c_data/panel_*/
  _supp_sig_heatmaps.R  FDR and Pi clustered heatmaps         -> c_data/sig_heatmaps/
  90_stitch_F03.R       folds every CSV into the workbook
  -> c_data/F03_supplementary.xlsx
```

```sh
Rscript 04_Figures/F03/a_script/90_stitch_F03.R
```

## What comes out

`MAIN_F03_composite` plus its four panels, `SUPP_F03_composite`, the four
`SUPP_F03_*_heatmap` plates, and a 17-sheet workbook: Overview, four ring-term sheets,
eight heatmap grouping sheets and four S4a source sheets. The per-contrast limma tables
are deliberately not copied in; S10 Table is the record for those.

## Orderings that matter

The stitcher sources `02_supp_panels.R`, then `01_main_panels.R`, then
`_supp_sig_heatmaps.R`, and builds the workbook last. Sheet specs are assembled by globbing
`c_data/panel_*/`, `c_data/sig_heatmaps/` and `c_data/supp/`, so every CSV must exist first.
An unlisted sheet stops the run.

`_`-prefixed scripts are sourced by `01_main_panels.R` or the stitcher, never run alone.

`cleanup_after_workbook()` deletes the folded CSVs and removes the `panel_A`-`panel_D`,
`sig_heatmaps` and `supp` subdirectories. The workbook is the durable copy.

## Cost

43.7 s. The fGSEA cache is read, not recomputed.

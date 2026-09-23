# F05 · Figure 5, WGCNA network and module-trait associations

Fits the co-expression network from the imputed matrix, then draws the module-trait
heatmap, the module NES scatters and the network QC plate. Was F06 until 2026-09-22.

```
02_imputation/c_data/01_imputed.csv + 01_DAList_imputed.rds
00_input/wgcna_reference_modules.csv + HPA_skeletal_muscle_annotations.tsv
03_DEP/c_data/03_combined_results.csv
  YvO_WGCNA_run.R   fit modules, eigengenes, trait stats -> c_data/wgcna/*.csv, c_data/*.rds
  01_main_panels.R  _panel_A_module_heatmap, _panel_B_nes_scatters, _module_labels,
                    then triptych / hub / preservation
  02_supp_panels.R  soft threshold, dendrogram, compartment, bicor
  90_stitch_F05.R   workbook, cleanup, Box copy
  -> c_data/F05_supplementary.xlsx
```

```sh
Rscript 04_Figures/F05/a_script/YvO_WGCNA_run.R
Rscript 04_Figures/F05/a_script/90_stitch_F05.R
```

## What comes out

`MAIN_F05_composite`, `SUPP_F05_composite`, `SUPP_F05_modules` with the per-module triptych
pages, and a 32-sheet workbook: module assignments, bio labels, hub proteins and edges,
kME, eigengenes, module-trait associations, preservation, SFT summary and four QC sheets.

`c_data` also keeps the objects F06 reads: `datExpr.rds`, `module_colors.rds`,
`me_pre.rds`, `me_post.rds`, `delta_me.rds`.

## Orderings that matter

`YvO_WGCNA_run.R` runs before `90_stitch_F05.R`; `run_all.R` lists both, the network fit as
its own step immediately before the stitcher. It is single-threaded via
`disableWGCNAThreads()` and seeded, so the TOM is reproducible.

Inside the stitcher, `01_main_panels.R` runs before `02_supp_panels.R`. The triptych, hub
and preservation scripts sit in `01_main_panels.R` so their CSVs exist when the workbook
specs are built.

`_supp_qc_compartment.R` asserts the F05 workbook exists and reads
`WGCNA_module_assignments` out of it, but the stitcher writes that workbook afterwards. On
any given run it therefore reads the previous run's copy. The values have converged, so
this is latent rather than active.

`f06_cleanup()`: the name predates the renumbering, preserves the `.rds` files,
`c_data/wgcna/` and a named list of CSVs, so F06's inputs survive the sweep.

## Cost

`YvO_WGCNA_run.R` 18.8 s; the stitcher 212.6 s. `modulePreservation` at 200 permutations is
cached in `c_data/modulePreservation_cache.rds` and re-run only when an input is newer.

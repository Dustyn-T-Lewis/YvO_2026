# F04

Draws Figure 4, the training-response concordance between age groups, S5a and S5b Figures, and the S5 Table workbook.

## Reads

- `03_DEP/c_data/03_combined_results.csv`
- `03_DEP/c_data/03_reversal_aging_fdr.csv` and `03_reversal_null_draws.csv`
- `02_imputation/c_data/01_DAList_imputed.rds` and `02_mar_mnar_classification.csv`
- `00_input/YvO_meta.xlsx`
- `04_Figures/shared/fgsea_tstat_all_v2.csv`, kept current by F02 and F03

## Writes

- `b_reports/F04.pdf`, `F04.png`: Figure 4
- `b_reports/S5a.pdf`, `S5a.png`: S5a Figure
- `b_reports/S5b.pdf`, `S5b.png`: S5b Figure
- `b_reports/panels/`: each panel on its own, named after its script, and the plots of S5 Table sheets that no figure shows
- `c_data/F04_data.xlsx`: S5 Table

## Run

```sh
Rscript 04_Figures/F04/a_script/F04.R
Rscript 04_Figures/F04/a_script/S5a.R
Rscript 04_Figures/F04/a_script/S5b.R
Rscript 04_Figures/F04/a_script/F04_data.R
```

Each panel script in `a_script/panels/` also runs on its own. Files starting with `_` hold code shared by panels and are only sourced. `S5T_` scripts plot S5 Table sheets that no figure shows. `F04_data.R` runs them.

## Order

Stage 03, F02 and F03 run first. `F04_data.R` runs last: it runs the six `S5T_` scripts, folds the CSVs they and the three composites leave in `c_data/` into the workbook, and deletes them. It stops without writing the workbook if any of those CSVs is missing.

A panel script run on its own leaves its CSVs in `c_data/`. The next `F04_data.R` removes them.

Panels D and E share one engine that writes fixed file names, and each renames its outputs after the engine returns. Do not run the two at the same time.

## Outputs and manuscript items

| File | Manuscript item |
|---|---|
| `b_reports/F04.pdf` | Figure 4 |
| `b_reports/panels/A_quadrant_ora.pdf` | Figure 4A |
| `b_reports/panels/B_rrho2.pdf` | Figure 4B |
| `b_reports/panels/B_legend.png` | Figure 4B colour key |
| `b_reports/panels/C_trajectory.pdf` | Figure 4C |
| `b_reports/panels/D_nes_concordance.pdf` | Figure 4D |
| `b_reports/panels/D_legend.png` | Figure 4D and 4E shape and size key |
| `b_reports/panels/E_nes_reversal.pdf` | Figure 4E |
| `b_reports/panels/E_legend.png` | Figure 4E shape and size key |
| `b_reports/panels/F_fry_barcode.pdf` | Figure 4F |
| `b_reports/S5a.pdf` | S5a Figure |
| `b_reports/panels/S5a_A_rho_bootstrap.pdf` | S5a Figure A |
| `b_reports/panels/S5a_B_goslim_bars.pdf` | S5a Figure B |
| `b_reports/panels/S5a_C_concordance_magnitude.pdf` | S5a Figure C |
| `b_reports/panels/S5a_D_coupling_null.pdf` | S5a Figure D |
| `b_reports/S5b.pdf` | S5b Figure |
| `b_reports/panels/S5b_young_dep_heatmap.pdf` | S5b Figure |
| `b_reports/panels/S5T_enrichment_heatmap.pdf` | S5 Table, sheet `SUPP_enrichment_blunting`, drawn |
| `b_reports/panels/S5T_rrho2_aging.pdf` | S5 Table, sheets `SUPP_rrho2_aging_*`, drawn |
| `b_reports/panels/S5T_ora_dedup.pdf` | S5 Table, sheet `SUPP_ora_dedup`, drawn |
| `b_reports/panels/S5T_threshold_sens.pdf` | S5 Table, sheet `SUPP_threshold_sens`, drawn |
| `b_reports/panels/S5T_fry_leading.pdf` | S5 Table, sheet `SUPP_fry_leading`, drawn |
| `b_reports/panels/S5T_cat_depth.pdf` | S5 Table, sheet `SUPP_cat_depth`, drawn |
| `c_data/F04_data.xlsx`, sheets `panel_A_*` to `panel_F_*` | S5 Table, Figure 4 source data |
| `c_data/F04_data.xlsx`, sheets `SUPP_rho_bootstrap`, `SUPP_goslim_bars`, `SUPP_concordance_magnitude`, `SUPP_coupling_null` | S5 Table, S5a source data |
| `c_data/F04_data.xlsx`, sheets `SUPP_young_dep_*` | S5 Table, S5b source data |
| `c_data/F04_data.xlsx`, sheets `SUPP_enrichment_blunting`, `SUPP_ora_dedup`, `SUPP_threshold_sens`, `SUPP_fry_leading`, `SUPP_rrho2_aging_*`, `SUPP_cat_depth` | S5 Table, diagnostics no figure draws |

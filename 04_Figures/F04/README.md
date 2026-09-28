# F04

Draws Figure 4 (training-response concordance between age groups) and S5a and S5b Figures, and writes S5 Table.

## Reads

- `03_DEP/c_data/03_combined_results.csv`
- `03_DEP/c_data/03_reversal_aging_fdr.csv` and `03_reversal_null_draws.csv`
- `02_imputation/c_data/01_DAList_imputed.rds` and `02_mar_mnar_classification.csv`
- `00_input/YvO_meta.xlsx`
- `04_Figures/shared/fgsea_tstat_all_v2.csv`, kept current by F02 and F03

## Writes

- `b_reports/main/` and `b_reports/supp/`: the composites and their `panels/`, listed in the table below, plus the `S5T_` plots of S5 Table sheets that no figure shows
- `c_data/F04_data.xlsx`: S5 Table

## Run

```sh
Rscript 04_Figures/F04/a_script/main/F04.R
Rscript 04_Figures/F04/a_script/supp/S5a.R
Rscript 04_Figures/F04/a_script/supp/S5b.R
Rscript 04_Figures/F04/a_script/F04_data.R
```

## Order

Stage 03, F02 and F03 run first. `F04_data.R` runs last. It runs the six `S5T_` scripts, folds their CSVs and the composites' into the workbook and deletes them, and stops without writing if any CSV is missing. A panel run on its own leaves CSVs that the next `F04_data.R` removes.

Panels D and E share an engine that writes fixed file names, which each renames afterwards, so never run the two at once.

## Outputs and manuscript items

| File | Manuscript item |
|---|---|
| `b_reports/main/F04.pdf` | Figure 4 |
| `b_reports/main/panels/A_quadrant_ora.pdf` | Figure 4A |
| `b_reports/main/panels/B_rrho2.pdf` | Figure 4B |
| `b_reports/main/panels/B_legend.png` | Figure 4B colour key |
| `b_reports/main/panels/C_trajectory.pdf` | Figure 4C |
| `b_reports/main/panels/D_nes_concordance.pdf` | Figure 4D |
| `b_reports/main/panels/D_legend.png` | Figure 4D and 4E shape and size key |
| `b_reports/main/panels/E_nes_reversal.pdf` | Figure 4E |
| `b_reports/main/panels/E_legend.png` | Figure 4E shape and size key |
| `b_reports/main/panels/F_fry_barcode.pdf` | Figure 4F |
| `b_reports/supp/S5a.pdf` | S5a Figure |
| `b_reports/supp/S5a_Figure.pdf` | S5a Figure with its legend, the official file |
| `b_reports/supp/panels/S5a_A_rho_bootstrap.pdf` | S5a Figure A |
| `b_reports/supp/panels/S5a_B_goslim_bars.pdf` | S5a Figure B |
| `b_reports/supp/panels/S5a_C_concordance_magnitude.pdf` | S5a Figure C |
| `b_reports/supp/panels/S5a_D_coupling_null.pdf` | S5a Figure D |
| `b_reports/supp/S5b.pdf` | S5b Figure |
| `b_reports/supp/S5b_Figure.pdf` | S5b Figure with its legend, the official file |
| `b_reports/supp/panels/S5b_young_dep_heatmap.pdf` | S5b Figure |
| `b_reports/supp/panels/S5T_enrichment_heatmap.pdf` | S5 Table, sheet `SUPP_enrichment_blunting`, drawn |
| `b_reports/supp/panels/S5T_rrho2_aging.pdf` | S5 Table, sheets `SUPP_rrho2_aging_*`, drawn |
| `b_reports/supp/panels/S5T_ora_dedup.pdf` | S5 Table, sheet `SUPP_ora_dedup`, drawn |
| `b_reports/supp/panels/S5T_threshold_sens.pdf` | S5 Table, sheet `SUPP_threshold_sens`, drawn |
| `b_reports/supp/panels/S5T_fry_leading.pdf` | S5 Table, sheet `SUPP_fry_leading`, drawn |
| `b_reports/supp/panels/S5T_cat_depth.pdf` | S5 Table, sheet `SUPP_cat_depth`, drawn |
| `c_data/F04_data.xlsx`, sheets `panel_A_*` to `panel_F_*` | S5 Table, Figure 4 source data |
| `c_data/F04_data.xlsx`, sheets `SUPP_rho_bootstrap`, `SUPP_goslim_bars`, `SUPP_concordance_magnitude`, `SUPP_coupling_null` | S5 Table, S5a source data |
| `c_data/F04_data.xlsx`, sheets `SUPP_young_dep_*` | S5 Table, S5b source data |
| `c_data/F04_data.xlsx`, sheets `SUPP_enrichment_blunting`, `SUPP_ora_dedup`, `SUPP_threshold_sens`, `SUPP_fry_leading`, `SUPP_rrho2_aging_*`, `SUPP_cat_depth` | S5 Table, diagnostics no figure draws |

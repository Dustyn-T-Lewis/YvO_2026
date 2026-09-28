# F01

Draws Figure 1 (training volume and hypertrophy) and S2 Figure, and writes Table 1 and S2 Table.

## Reads

- `00_input/YvO_meta.xlsx`

## Writes

- `b_reports/main/` and `b_reports/supp/`: the composites and their `panels/`, listed in the table below, plus `F01_single_col.pdf` and `.png`, Figure 1 at single-column width
- `c_data/F01_table_1*.csv`: Table 1
- `c_data/F01_data.xlsx`: S2 Table

## Run

```sh
Rscript 04_Figures/F01/a_script/supp/S2.R
Rscript 04_Figures/F01/a_script/main/F01.R
Rscript 04_Figures/F01/a_script/F01_data.R
```

`_prepost.R` holds the pre/post panel and its 2 x 2 mixed ANOVA, used by every panel except A. `04_phenotype_table.R` builds Table 1 and is sourced by `F01_data.R`.

## Order

`F01_data.R` runs last. It builds Table 1B from the summary CSVs the panels leave in `c_data/`, so each table cell and figure p-value come from one calculation, and stops if any is missing. It then folds tables and panel CSVs into the workbook and deletes the panel CSVs.

Panels A, S2 A and S2 D seed their jitter; the others continue the previous panel's random stream, so they match the composites only when run through `F01.R` and `S2.R`.

## Outputs and manuscript items

| File | Manuscript item |
|---|---|
| `b_reports/main/F01.pdf` | Figure 1 |
| `b_reports/main/F01_single_col.pdf` | Figure 1, single-column layout |
| `b_reports/main/panels/A_training_volume.pdf` | Figure 1A |
| `b_reports/main/panels/B_dxa_lbm.pdf` | Figure 1B |
| `b_reports/main/panels/C_vl_thickness.pdf` | Figure 1C |
| `b_reports/supp/S2.pdf` | S2 Figure |
| `b_reports/supp/S2_Figure.pdf` | S2 Figure with its legend, the official file |
| `b_reports/supp/panels/S2_A_deadlift_1rm.pdf` | S2 Figure A |
| `b_reports/supp/panels/S2_B_type_II_fcsa.pdf` | S2 Figure B |
| `b_reports/supp/panels/S2_C_type_I_fcsa.pdf` | S2 Figure C |
| `b_reports/supp/panels/S2_D_dxa_fat_mass.pdf` | S2 Figure D |
| `b_reports/supp/panels/S2_E_fat_to_lean.pdf` | S2 Figure E |
| `c_data/F01_table_1a_characteristics.csv` | Table 1, baseline characteristics |
| `c_data/F01_table_1b_pre_post.csv` | Table 1, pre and post training |
| `c_data/F01_table_1c_composition.csv` | Table 1, cohort composition |
| `c_data/F01_data.xlsx`, sheets `Table_1A_characteristics` to `Table_1C_composition` | S2 Table, Table 1 as printed |
| `c_data/F01_data.xlsx`, sheets `Per_participant` and `Notes` | S2 Table, per-participant values and measure notes |
| `c_data/F01_data.xlsx`, sheets `panel_A_train_volume` to `panel_C_vl_thickness` | S2 Table, Figure 1 source data |
| `c_data/F01_data.xlsx`, sheets `supp_*` | S2 Table, S2 Figure source data |

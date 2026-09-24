# F03

Draws Figure 3, the volcano rings for the four contrasts, and S4a and S4b Figures.

## Reads

- `03_DEP/c_data/03_combined_results.csv`
- `03_DEP/c_data/03_DEP_results.xlsx`, sheets `Aging`, `Training_Young`, `Training_Old`, `Interaction` and `outlier_sensitivity`
- `04_Figures/shared/fgsea_tstat_all_v2.csv`, rebuilt by `shared/build_fgsea_cache.R` when older than the DEP results

## Writes

- `b_reports/main/F03.pdf`, `F03.png`: Figure 3
- `b_reports/supp/S4a.pdf`, `S4a.png`: S4a Figure
- `b_reports/supp/S4a_Figure.pdf`: S4a Figure with its legend, the official supplementary file
- `b_reports/supp/S4b.pdf`: S4b Figure, four pages; `S4b_A.png` to `S4b_D.png`, one per page
- `b_reports/supp/S4b_Figure.pdf`: S4b Figure with its legend, the official supplementary file
- `b_reports/main/panels/`, `b_reports/supp/panels/`: each panel on its own, named after its script
- `c_data/F03_data.xlsx`: S4 Table

## Run

```sh
Rscript 04_Figures/F03/a_script/main/F03.R
Rscript 04_Figures/F03/a_script/supp/S4a.R
Rscript 04_Figures/F03/a_script/supp/S4b.R
Rscript 04_Figures/F03/a_script/F03_data.R
```

Each panel script in `a_script/main/panels/` and `a_script/supp/panels/` also runs on its own. Files starting with `_` hold code shared by several panels and are only sourced.

## Order

Stage 03 runs first. `F03_data.R` runs last: it folds the CSVs the three composites leave in `c_data/` into the workbook and deletes them.

## Outputs and manuscript items

| File | Manuscript item |
|---|---|
| `b_reports/main/F03.pdf` | Figure 3 |
| `b_reports/main/panels/A_volcano_aging.pdf` | Figure 3A |
| `b_reports/main/panels/B_volcano_young.pdf` | Figure 3B |
| `b_reports/main/panels/C_volcano_old.pdf` | Figure 3C |
| `b_reports/main/panels/D_volcano_interaction.pdf` | Figure 3D |
| `b_reports/supp/S4a.pdf` | S4a Figure |
| `b_reports/supp/S4a_Figure.pdf` | S4a Figure with its legend, the official file |
| `b_reports/supp/panels/S4a_A_p_value.pdf` | S4a Figure A |
| `b_reports/supp/panels/S4a_B_pi_score.pdf` | S4a Figure B |
| `b_reports/supp/panels/S4a_C_fdr.pdf` | S4a Figure C |
| `b_reports/supp/panels/S4a_D_ma.pdf` | S4a Figure D |
| `b_reports/supp/panels/S4a_E_outlier_retention.pdf` | S4a Figure E |
| `b_reports/supp/S4b.pdf` | S4b Figure |
| `b_reports/supp/S4b_Figure.pdf` | S4b Figure with its legend, the official file |
| `b_reports/supp/panels/S4b_A_aging_fdr.pdf` | S4b Figure A |
| `b_reports/supp/panels/S4b_B_aging_pi.pdf` | S4b Figure B |
| `b_reports/supp/panels/S4b_C_training_fdr.pdf` | S4b Figure C |
| `b_reports/supp/panels/S4b_D_training_pi.pdf` | S4b Figure D |
| `c_data/F03_data.xlsx`, sheets `RING_A` to `RING_D` | S4 Table, terms each Figure 3 ring labels |
| `c_data/F03_data.xlsx`, sheets `AgeFDR_*`, `AgePi_*`, `TrFDR_*`, `TrPi_*` | S4 Table, S4b groupings |
| `c_data/F03_data.xlsx`, sheets `SUPP_panel_*` | S4 Table, S4a source data |

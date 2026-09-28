# F00

Draws the pipeline QC figures S1a and S1b and writes S1 Table.

## Reads

- `01_normalization/c_data/00_report_intermediates.rds`
- `02_imputation/c_data/00_report_intermediates.rds`
- `02_imputation/c_data/benchmark/04_composite_ranking.csv`, written by the opt-in `02_imputation/a_script/benchmark/_run_all.R`
- `03_DEP/c_data/03_DEP_results.xlsx`, sheet `DA_summary`

## Writes

- `b_reports/supp/`: the composites and their `panels/`, listed in the table below
- `c_data/F00_data.xlsx`: S1 Table

## Run

```sh
Rscript 04_Figures/F00/a_script/supp/S1a.R
Rscript 04_Figures/F00/a_script/supp/S1b.R
Rscript 04_Figures/F00/a_script/F00_data.R
```

`_inputs.R` loads the stage 01 and 02 intermediates for every panel.

## Order

Stages 01 to 03 run first. `F00_data.R` runs last, folding the data frames the panels leave in `c_data/sheets/` into the workbook and deleting them. Panel S1b J stops if the opt-in benchmark ranking is missing.

Panel S1a F jitters its points, drawing a new jitter at every render. It calls `set.seed(42)` before the plot and before its PDF, so the composite always starts from the same state, but the S1a PDF and PNG differ in jitter.

## Outputs and manuscript items

| File | Manuscript item |
|---|---|
| `b_reports/supp/S1a.pdf` | S1a Figure |
| `b_reports/supp/S1a_Figure.pdf` | S1a Figure with its legend, the official file |
| `b_reports/supp/panels/S1a_A_filter_cascade.pdf` | S1a Figure A |
| `b_reports/supp/panels/S1a_B_protein_missingness.pdf` | S1a Figure B |
| `b_reports/supp/panels/S1a_C_pca_pre.pdf` | S1a Figure C |
| `b_reports/supp/panels/S1a_D_pca_post.pdf` | S1a Figure D |
| `b_reports/supp/panels/S1a_E_eta_squared.pdf` | S1a Figure E |
| `b_reports/supp/panels/S1a_F_outlier_consensus.pdf` | S1a Figure F |
| `b_reports/supp/panels/S1a_G_sample_missingness.pdf` | S1a Figure G |
| `b_reports/supp/S1b.pdf` | S1b Figure |
| `b_reports/supp/S1b_Figure.pdf` | S1b Figure with its legend, the official file |
| `b_reports/supp/panels/S1b_H_miss_class_scatter.pdf` | S1b Figure H |
| `b_reports/supp/panels/S1b_I_miss_class_bar.pdf` | S1b Figure I |
| `b_reports/supp/panels/S1b_J_benchmark.pdf` | S1b Figure J |
| `b_reports/supp/panels/S1b_K_imputation_density.pdf` | S1b Figure K |
| `b_reports/supp/panels/S1b_L_mnar_shift.pdf` | S1b Figure L |
| `b_reports/supp/panels/S1b_M_sample_integrity.pdf` | S1b Figure M |
| `b_reports/supp/panels/S1b_N_dep_heatmap.pdf` | S1b Figure N |
| `c_data/F00_data.xlsx`, sheets `panel_A` to `panel_G` | S1 Table, S1a source data |
| `c_data/F00_data.xlsx`, sheets `panel_H` to `panel_N` | S1 Table, S1b source data |

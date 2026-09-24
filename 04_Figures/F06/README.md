# F06

Draws Figure 6, the module ROCs and the module-phenotype coupling grid, and S7 Figure, and writes S7 Table.

## Reads

- `04_Figures/F05/c_data/F05_data.xlsx`, sheets `MEs`, `me_pre`, `me_post`, `delta_me`, `metadata_subj_age`, `metadata_pheno_wide`, `common_subj`, `WGCNA_mod_bio_labels` and `WGCNA_module_assignments`
- `04_Figures/F05/c_data/datExpr.rds`, `module_colors.rds`, `me_pre.rds` and `me_post.rds`
- `02_imputation/c_data/01_imputed.csv`
- `03_DEP/c_data/03_combined_results.csv`

## Writes

- `b_reports/main/F06.pdf`, `F06.png`: Figure 6
- `b_reports/supp/S7.pdf`, `S7.png`: S7 Figure
- `b_reports/supp/S7_Figure.pdf`: S7 Figure with its legend, the official supplementary file
- `b_reports/main/panels/`, `b_reports/supp/panels/`: each panel on its own, named after its script
- `c_data/F06_data.xlsx`: S7 Table

## Run

```sh
Rscript 04_Figures/F06/a_script/main/F06.R
Rscript 04_Figures/F06/a_script/supp/S7.R
Rscript 04_Figures/F06/a_script/F06_data.R
```

Each panel script in `a_script/main/panels/` and `a_script/supp/panels/` also runs on its own. Files starting with `_` hold code shared by several scripts and are only sourced.

## Order

F05 runs first. `A_age_roc.R` and `B_training_roc.R` draw from the module grid that `S7_A_module_grid.R` writes, and run it first when the grid is missing. `S7_B_full_sweep.R` draws from the 180-test screen that `C_hero_grid.R` writes, and runs it first when the screen is missing.

`F06_data.R` runs last. It runs the four analysis steps in `a_script/` that compute sheets but draw nothing: `_supp_prepare_roc.R`, `_supp_multivariate.R`, `_supp_loso_sensitivity.R` and `_supp_loso_wgcna_refit.R`. It then folds their CSVs and the ones the two composites leave in `c_data/` into the workbook and deletes them. `_loocv.R` holds the nested LOOCV classifier the first two share. `_panel_selection.R` names the cells Figure 6A and 6B draw, for the ROC panels and the two LOSO steps.

`_supp_loso_wgcna_refit.R` refits the network on 30 leave-one-subject-out folds with no cache and takes about 255 s.

## Outputs and manuscript items

| File | Manuscript item |
|---|---|
| `b_reports/main/F06.pdf` | Figure 6 |
| `b_reports/main/panels/A_age_roc.pdf` | Figure 6A |
| `b_reports/main/panels/B_training_roc.pdf` | Figure 6B |
| `b_reports/main/panels/C_hero_grid.pdf` | Figure 6C |
| `b_reports/supp/S7.pdf` | S7 Figure |
| `b_reports/supp/S7_Figure.pdf` | S7 Figure with its legend, the official file |
| `b_reports/supp/panels/S7_A_module_grid.pdf` | S7 Figure A |
| `b_reports/supp/panels/S7_B_full_sweep.pdf` | S7 Figure B |
| `c_data/F06_data.xlsx`, sheets `module_grid_summary`, `module_grid_curves` | S7 Table, S7 Figure A and the cells Figure 6A and 6B draw |
| `c_data/F06_data.xlsx`, sheet `panel_B_full_screen` | S7 Table, the 180-test screen behind Figure 6C and S7 Figure B |
| `c_data/F06_data.xlsx`, sheets `panel_A_classifier_auc`, `panel_A_feature_stability`, `panel_A_permutation`, `panel_A_roc_curves` | S7 Table, multivariate age classifiers |
| `c_data/F06_data.xlsx`, sheets `classifier_pilot_summary`, `classifier_pilot_curves` | S7 Table, classifier pilot |
| `c_data/F06_data.xlsx`, sheets `loso_auc_summary`, `loso_wgcna_refit_summary`, `loso_wgcna_refit_mod_stability` | S7 Table, leave-one-subject-out sensitivity of the Figure 6A and 6B cells |

`abstract_panels` reads the sheet `module_grid_summary`.

# F02

Draws Figure 2, the proteome overview and differential expression landscape, and S3 Figure, the coefficient-of-variation diagnostics.

## Reads

- `01_normalization/c_data/03_DAList_normalized.rds`
- `02_imputation/c_data/01_DAList_imputed.rds`
- `02_imputation/c_data/02_imputation.xlsx`, sheets `imputation_mask` and `mar_mnar_classification`
- `03_DEP/c_data/03_combined_results.csv`
- `03_DEP/c_data/03_DEP_results.xlsx`, sheet `blunting`
- `04_Figures/shared/fgsea_tstat_all_v2.csv`, rebuilt by `shared/build_fgsea_cache.R` when older than the DEP results

## Writes

- `b_reports/F02.pdf`, `F02.png`: Figure 2
- `b_reports/S3.pdf`, `S3.png`: S3 Figure
- `b_reports/panels/`: each panel on its own, named after its script
- `c_data/F02_data.xlsx`: S3 Table

## Run

```sh
Rscript 04_Figures/F02/a_script/F02.R
Rscript 04_Figures/F02/a_script/S3.R
Rscript 04_Figures/F02/a_script/F02_data.R
```

Each panel script in `a_script/panels/` also runs on its own. Files starting with `_` hold code shared by several panels and are only sourced.

## Order

Stage 03 runs first. `F02_data.R` runs last: it folds the CSVs the two composites leave in `c_data/` into the workbook and deletes them. The workbook takes its sheets from a glob of `c_data/`, so sheet order follows the file names.

## Outputs and manuscript items

| File | Manuscript item |
|---|---|
| `b_reports/F02.pdf` | Figure 2 |
| `b_reports/panels/A_pca.pdf` | Figure 2A |
| `b_reports/panels/B_logfc_density.pdf` | Figure 2B |
| `b_reports/panels/C_dep_counts.pdf` | Figure 2C |
| `b_reports/panels/D_upset.pdf` | Figure 2D |
| `b_reports/panels/E_fgsea.pdf` | Figure 2E |
| `b_reports/panels/F_barcode.pdf` | Figure 2F |
| `b_reports/S3.pdf` | S3 Figure |
| `b_reports/panels/S3_A_cv_scatter.pdf` | S3 Figure A |
| `b_reports/panels/S3_B_cv_violin.pdf` | S3 Figure B |
| `b_reports/panels/S3_C_imputed.pdf` | S3 Figure C |
| `c_data/F02_data.xlsx`, sheets `panel_A_*` to `panel_F_*` | S3 Table, Figure 2 source data |
| `c_data/F02_data.xlsx`, sheets `SUPP_panel_A_*` to `SUPP_panel_C_*` | S3 Table, S3 Figure source data |

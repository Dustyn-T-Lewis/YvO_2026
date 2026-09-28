# F02

Draws Figure 2 (proteome overview and differential expression) and S3 Figure (coefficient-of-variation diagnostics).

## Reads

- `01_normalization/c_data/03_DAList_normalized.rds`
- `02_imputation/c_data/01_DAList_imputed.rds`
- `02_imputation/c_data/02_imputation.xlsx`, sheets `imputation_mask` and `mar_mnar_classification`
- `03_DEP/c_data/03_combined_results.csv`
- `03_DEP/c_data/03_DEP_results.xlsx`, sheet `blunting`
- `04_Figures/shared/fgsea_tstat_all_v2.csv`, rebuilt by `shared/build_fgsea_cache.R` when older than the DEP results

## Writes

- `b_reports/main/` and `b_reports/supp/`: the composites and their `panels/`, listed in the table below
- `c_data/F02_data.xlsx`: S3 Table

## Run

```sh
Rscript 04_Figures/F02/a_script/main/F02.R
Rscript 04_Figures/F02/a_script/supp/S3.R
Rscript 04_Figures/F02/a_script/F02_data.R
```

## Order

Stage 03 runs first. `F02_data.R` runs last, folding the CSVs the composites leave in `c_data/` into the workbook, in file-name order, and deleting them.

## Outputs and manuscript items

| File | Manuscript item |
|---|---|
| `b_reports/main/F02.pdf` | Figure 2 |
| `b_reports/main/panels/A_pca.pdf` | Figure 2A |
| `b_reports/main/panels/B_logfc_density.pdf` | Figure 2B |
| `b_reports/main/panels/C_dep_counts.pdf` | Figure 2C |
| `b_reports/main/panels/D_upset.pdf` | Figure 2D |
| `b_reports/main/panels/E_fgsea.pdf` | Figure 2E |
| `b_reports/main/panels/F_barcode.pdf` | Figure 2F |
| `b_reports/supp/S3.pdf` | S3 Figure |
| `b_reports/supp/S3_Figure.pdf` | S3 Figure with its legend, the official file |
| `b_reports/supp/panels/S3_A_cv_scatter.pdf` | S3 Figure A |
| `b_reports/supp/panels/S3_B_cv_violin.pdf` | S3 Figure B |
| `b_reports/supp/panels/S3_C_imputed.pdf` | S3 Figure C |
| `c_data/F02_data.xlsx`, sheets `panel_A_*` to `panel_F_*` | S3 Table, Figure 2 source data |
| `c_data/F02_data.xlsx`, sheets `SUPP_panel_A_*` to `SUPP_panel_C_*` | S3 Table, S3 Figure source data |

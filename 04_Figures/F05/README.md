# F05

Fits the WGCNA co-expression network and draws Figure 5 (module-trait heatmap and module NES scatters) and S6 and S8 Figures.

## Reads

- `02_imputation/c_data/01_imputed.csv` and `01_DAList_imputed.rds`, in `YvO_WGCNA_run.R`; `01_imputed.csv` also in `S6_D_bicor.R`
- `00_input/wgcna_reference_modules.csv`, in `YvO_WGCNA_run.R`
- `00_input/HPA_skeletal_muscle_annotations.tsv`, in `S6_C_compartment.R`
- `03_DEP/c_data/03_combined_results.csv`, in `B_nes_scatters.R`
- `c_data/wgcna/` and the other network outputs `YvO_WGCNA_run.R` writes to `c_data/`
- `c_data/F05_data.xlsx`, sheets `WGCNA_module_assignments` and `WGCNA_mod_bio_labels`, in `S6_C_compartment.R`

## Writes

- `b_reports/main/` and `b_reports/supp/`: the composites and their `panels/`, listed in the table below, plus a triptych for each of the five modules S8 Figure leaves out
- `c_data/F05_data.xlsx`: S6 Table
- `c_data/wgcna/`, `c_data/*.rds`, `mod_bio_labels.csv`, `wgcna_kme_all.csv` and the other network outputs of `YvO_WGCNA_run.R`
- `c_data/wgcna_string_clusters.csv`: the STRING cluster cache, from `shared/build_string_cluster_cache.R`, sourced by `F05_data.R`
- `c_data/modulePreservation_cache.rds`: `modulePreservation` at 200 permutations, reused while newer than `datExpr.rds` and `module_colors.rds`
- `c_data/03_panel_B_eigengene_data.csv`, `03_panel_B_triptych_enrichment.csv`: from `panels/_triptych.R`

F06 reads `datExpr.rds`, `me_pre.rds`, `me_post.rds`, `module_colors.rds` and `F05_data.xlsx`. `shared/build_string_cluster_cache.R` reads `mod_bio_labels.csv` and `wgcna/wgcna_module_assignments.csv`.

## Run

```sh
Rscript 04_Figures/F05/a_script/YvO_WGCNA_run.R
Rscript 04_Figures/F05/a_script/main/F05.R
Rscript 04_Figures/F05/a_script/supp/S6.R
Rscript 04_Figures/F05/a_script/supp/S8.R
Rscript 04_Figures/F05/a_script/F05_data.R
```

`_module_labels.R`, `_supp_mod_hub.R` and `_supp_preservation.R` in `a_script/` compute S6 Table sheets for `F05_data.R` and draw nothing.

## Order

Stages 02 and 03 run first, then `YvO_WGCNA_run.R` (single-threaded and seeded, so reproducible), before the composites and F06. `F05_data.R` runs last, folding the composites' CSVs into the workbook and deleting them.

`S6_C_compartment.R` reads module assignments and labels from the `F05_data.xlsx` left by the previous run. Both sheets copy `c_data/wgcna/wgcna_module_assignments.csv` and `c_data/mod_bio_labels.csv`, so S6 Figure C lags one run only when `YvO_WGCNA_run.R` has changed them.

## Outputs and manuscript items

| File | Manuscript item |
|---|---|
| `b_reports/main/F05.pdf` | Figure 5 |
| `b_reports/main/panels/A_module_trait_heatmap.pdf` | Figure 5A |
| `b_reports/main/panels/B_nes_scatters.pdf` | Figure 5B |
| `b_reports/main/panels/B_nes_scatters_legend.png` | Figure 5B, dot size key |
| `b_reports/supp/S6.pdf` | S6 Figure |
| `b_reports/supp/S6_Figure.pdf` | S6 Figure with its legend, the official file |
| `b_reports/supp/panels/S6_A_soft_threshold.pdf` | S6 Figure A |
| `b_reports/supp/panels/S6_B_dendrogram.pdf` | S6 Figure B |
| `b_reports/supp/panels/S6_C_compartment.pdf` | S6 Figure C |
| `b_reports/supp/panels/S6_D_bicor.pdf` | S6 Figure D |
| `b_reports/supp/S8.pdf` | S8 Figure |
| `b_reports/supp/S8_Figure.pdf` | S8 Figure with its legend, the official file |
| `b_reports/supp/panels/S8_A_turquoise.pdf` | S8 Figure A |
| `b_reports/supp/panels/S8_B_black.pdf` | S8 Figure B |
| `b_reports/supp/panels/S8_C_yellow.pdf` | S8 Figure C |
| `b_reports/supp/panels/S8_D_blue.pdf` | S8 Figure D |
| `b_reports/main/panels/other_triptychs_*.pdf` | None |
| `c_data/F05_data.xlsx`, sheet `panel_A_heatmap` | S6 Table, Figure 5A source data |
| `c_data/F05_data.xlsx`, sheet `panel_B_module_fgsea` | S6 Table, Figure 5B source data |
| `c_data/F05_data.xlsx`, sheets `SUPP_panel_A_sft_fit` to `SUPP_panel_D_bicor_sensitivity` | S6 Table, S6 Figure A to D source data |
| `c_data/F05_data.xlsx`, sheets `WGCNA_protein_zscores_by_group` and `WGCNA_module_enrichment` | S6 Table, S8 Figure source data |
| `c_data/F05_data.xlsx`, sheets `WGCNA_module_core`, `WGCNA_module_ora`, `WGCNA_string_clusters` and `WGCNA_ora_terms` | S6 Table, evidence behind each module name |
| `c_data/F05_data.xlsx`, all other sheets | S6 Table, module assignments, hub proteins, kME, eigengenes, module-trait associations, preservation and metadata |

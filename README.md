# YvO_2026

Analysis code for the skeletal muscle proteome response to resistance training in younger and older adults: vastus lateralis biopsies from 32 participants before and after training, 64 samples, quantified by DIA mass spectrometry. Without the two consensus outliers, `Y_S05_Pre` and `Y_S07_Post`, the matrix is 2,106 proteins by 62 samples. limma's four contrasts, Aging, Training_Young, Training_Old and Interaction, give 278, 135, 0 and 1 proteins at FDR < 0.05.

```sh
git clone --depth 1 https://github.com/Dustyn-T-Lewis/YvO_2026.git   # 74 MB; a full clone carries 1.5 GB of history
Rscript setup.R      # once after cloning
Rscript run_all.R    # every step, each in its own R session
```

## Layout

```
00_input/          source data
01_normalization/  filter, outlier removal, cyclic loess     S8 Table
02_imputation/     missForest and the benchmark behind it     S9 Table
03_DEP/            limma and the supplementary analyses       S10 Table, S9 Figure
04_Figures/
  F00              S1a and S1b Figures, S1 Table
  F01 to F06       Figures 1 to 6; the directory number is the figure number
  abstract_panels  the graphical abstract
  shared/          code the figure scripts source
```

Every stage and figure directory has a `README.md`, `a_script/` for code, `b_reports/` for renders and `c_data/` for tables. Figure directories split code and renders into `main/` and `supp/`, each with a `panels/` folder of one script or render per panel; every panel script also runs on its own. Composites are named after the manuscript item they draw. Each supplementary figure also gets `<item>_Figure.pdf`, the official file, with its legend from `04_Figures/shared/supp_legends.txt`. Files starting with `_` are only sourced.

Git tracks the cited renders and workbooks and each stage's tables; a run rewrites the ignored panel renders, reports and R objects.

## Reproducing

`setup.R` runs `renv::restore()`, installing each package at its `renv.lock` version, including proteoDA, RRHO2 and DreamAI from GitHub.

`run_all.R` runs each step as a separate `Rscript` and logs to `.runlogs/`. A full run takes about 22 minutes on an Apple silicon Mac, mostly `F05_data.R`, `F06_data.R` and the effect-size bootstrap in `03_DEP/a_script/supp/01`.

## Manuscript items

| Item | Directory | File |
|---|---|---|
| Figure 1 | `04_Figures/F01` | `b_reports/main/F01.pdf` |
| Figure 2 | `04_Figures/F02` | `b_reports/main/F02.pdf` |
| Figure 3 | `04_Figures/F03` | `b_reports/main/F03.pdf` |
| Figure 4 | `04_Figures/F04` | `b_reports/main/F04.pdf` |
| Figure 5 | `04_Figures/F05` | `b_reports/main/F05.pdf` |
| Figure 6 | `04_Figures/F06` | `b_reports/main/F06.pdf` |
| Table 1 | `04_Figures/F01` | `c_data/F01_table_1a_characteristics.csv`, `_1b_pre_post.csv`, `_1c_composition.csv` |
| Graphical abstract | `04_Figures/abstract_panels` | `b_reports/main/abstract.pdf`, finished by hand |
| S1a Figure | `04_Figures/F00` | `b_reports/supp/S1a_Figure.pdf` |
| S1b Figure | `04_Figures/F00` | `b_reports/supp/S1b_Figure.pdf` |
| S2 Figure | `04_Figures/F01` | `b_reports/supp/S2_Figure.pdf` |
| S3 Figure | `04_Figures/F02` | `b_reports/supp/S3_Figure.pdf` |
| S4a Figure | `04_Figures/F03` | `b_reports/supp/S4a_Figure.pdf` |
| S4b Figure | `04_Figures/F03` | `b_reports/supp/S4b_Figure.pdf` |
| S5a Figure | `04_Figures/F04` | `b_reports/supp/S5a_Figure.pdf` |
| S5b Figure | `04_Figures/F04` | `b_reports/supp/S5b_Figure.pdf` |
| S6 Figure | `04_Figures/F05` | `b_reports/supp/S6_Figure.pdf` |
| S7 Figure | `04_Figures/F06` | `b_reports/supp/S7_Figure.pdf` |
| S8 Figure | `04_Figures/F05` | `b_reports/supp/S8_Figure.pdf` |
| S9 Figure | `03_DEP` | `b_reports/supp/S9_Figure.pdf` |
| S1 Table | `04_Figures/F00` | `c_data/F00_data.xlsx` |
| S2 Table | `04_Figures/F01` | `c_data/F01_data.xlsx` |
| S3 Table | `04_Figures/F02` | `c_data/F02_data.xlsx` |
| S4 Table | `04_Figures/F03` | `c_data/F03_data.xlsx` |
| S5 Table | `04_Figures/F04` | `c_data/F04_data.xlsx` |
| S6 Table | `04_Figures/F05` | `c_data/F05_data.xlsx` |
| S7 Table | `04_Figures/F06` | `c_data/F06_data.xlsx` |
| S8 Table | `01_normalization` | `c_data/01_normalization.xlsx` |
| S9 Table | `02_imputation` | `c_data/02_imputation.xlsx` |
| S10 Table | `03_DEP` | `c_data/03_DEP_results.xlsx` |

Each figure has a PNG beside its PDF, and each supplementary figure a PDF without its legend. S8 and S9 Figures are numbered after S7, not by stage, so citations to S6 Figure and S10 Table kept their numbers.

## Known limitations

- `F05/a_script/YvO_WGCNA_run.R:85` reads `sft$fitIndices$slope` by position, not by power. It is right only because `powers` is `1:20`.
- The scale-free slope at power 12 is -2.43 (R-squared 0.877), outside the commented -1 to -2 range. The check warns only above -1, so it passes silently.
- `F02/a_script/main/panels/A_pca.R:30` blocks PERMANOVA permutations by subject. No permutation changes a subject's age group, so the age p-value is not a permutation test; its R-squared still describes the data.
- `F06/a_script/_supp_prepare_roc.R` uses 200 permutations, so a reported p of 0.005 is the floor.
- `cor <- WGCNA::cor` is set in `YvO_WGCNA_run.R` and `supp/panels/S6_D_bicor.R` and never restored; no later code calls bare `cor()`.
- `shared/print_scale_apply.R` changes `style.R`'s size globals and does not restore them; see `04_Figures/shared/README.md`.

## Citation

Zenodo concept DOI [10.5281/zenodo.19886624](https://doi.org/10.5281/zenodo.19886624); see `CITATION.cff`. MIT licensed (`LICENSE`), except the vendored GPL-3 GSimp source in `02_imputation/a_script/benchmark/methods/gsimp_source/`.

# YvO_2026

Analysis code for the skeletal muscle proteome response to resistance training in younger and older adults. Vastus lateralis biopsies came from 32 participants before and after training, 64 samples, quantified by DIA mass spectrometry. After removing the two consensus outlier samples, `Y_S05_Pre` and `Y_S07_Post`, the analysed matrix is 2,106 proteins across 62 samples. limma tests four contrasts: Aging, Training_Young, Training_Old and Interaction. At FDR < 0.05 they give 278, 135, 0 and 1 proteins.

```sh
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

Each stage and figure directory has a `README.md`, and the same shape: `a_script/` for code, `b_reports/` for renders, `c_data/` for tables. In a figure directory, `a_script/panels/` holds one script per panel, and each composite script is named after the manuscript item it draws. Files starting with `_` hold code that other scripts source.

Git tracks the renders and workbooks the manuscript cites, and the tables each stage writes. Panel renders, stage reports and intermediate R objects are ignored; a run rewrites them.

## Reproducing

`setup.R` calls `renv::restore()`, which installs every package at the version in `renv.lock`: CRAN, Bioconductor, and proteoDA, RRHO2 and DreamAI from GitHub.

`run_all.R` runs each step as a separate `Rscript` and logs to `.runlogs/`. A full run is 39 steps and takes about 19 minutes on an Apple silicon Mac. Three steps take most of it: the leave-one-subject-out network refit in `F06_data.R` (195 s), the module preservation permutations in `F05_data.R` (166 s) and the effect-size bootstrap in `03_DEP/a_script/supp/01` (160 s).

## Manuscript items

| Item | Directory | File |
|---|---|---|
| Figure 1 | `04_Figures/F01` | `b_reports/F01.pdf` |
| Figure 2 | `04_Figures/F02` | `b_reports/F02.pdf` |
| Figure 3 | `04_Figures/F03` | `b_reports/F03.pdf` |
| Figure 4 | `04_Figures/F04` | `b_reports/F04.pdf` |
| Figure 5 | `04_Figures/F05` | `b_reports/F05.pdf` |
| Figure 6 | `04_Figures/F06` | `b_reports/F06.pdf` |
| Table 1 | `04_Figures/F01` | `c_data/F01_table_1a_characteristics.csv`, `_1b_pre_post.csv`, `_1c_composition.csv` |
| Graphical abstract | `04_Figures/abstract_panels` | `b_reports/abstract.pdf`, finished by hand |
| S1a Figure | `04_Figures/F00` | `b_reports/S1a.pdf` |
| S1b Figure | `04_Figures/F00` | `b_reports/S1b.pdf` |
| S2 Figure | `04_Figures/F01` | `b_reports/S2.pdf` |
| S3 Figure | `04_Figures/F02` | `b_reports/S3.pdf` |
| S4a Figure | `04_Figures/F03` | `b_reports/S4a.pdf` |
| S4b Figure | `04_Figures/F03` | `b_reports/S4b.pdf` |
| S5a Figure | `04_Figures/F04` | `b_reports/S5a.pdf` |
| S5b Figure | `04_Figures/F04` | `b_reports/S5b.pdf` |
| S6 Figure | `04_Figures/F05` | `b_reports/S6.pdf` |
| S7 Figure | `04_Figures/F06` | `b_reports/S7.pdf` |
| S8 Figure | `04_Figures/F05` | `b_reports/S8.pdf` |
| S9 Figure | `03_DEP` | `b_reports/S9.pdf` |
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

Every render has a PNG beside the PDF. S8 and S9 Figures follow S7 rather than their stages so that citations to S6 Figure and S10 Table did not move.

## Known limitations

- `F05/a_script/YvO_WGCNA_run.R:91` reads `sft$fitIndices$slope` by position, not by power. It is right only because `powers` is `1:20`.
- The scale-free slope at power 12 is -2.43, outside the -1 to -2 range the code comment gives. The check warns only above -1, so a steeper slope passes without a note. R-squared is 0.877.
- `F02/a_script/panels/A_pca.R:30` runs PERMANOVA with permutations blocked by subject. No permutation changes a subject's age group, so the age p-value is not a permutation test. Its R-squared is still a valid description.
- `F06/a_script/_supp_prepare_roc.R` uses 200 permutations, so its smallest possible p is 0.005. A reported p of 0.005 is that floor.
- `cor <- WGCNA::cor` is set in `YvO_WGCNA_run.R` and `panels/S6_D_bicor.R` and never restored. No code after either calls bare `cor()`, so nothing is affected yet.
- `shared/print_scale_apply.R` changes `style.R`'s size globals and does not restore them; see `04_Figures/shared/README.md`.

## Citation

Archived on Zenodo, concept DOI [10.5281/zenodo.19886624](https://doi.org/10.5281/zenodo.19886624); see `CITATION.cff`. The code is MIT licensed (`LICENSE`), except the vendored GSimp source in `02_imputation/a_script/benchmark/methods/gsimp_source/`, which is GPL-3.

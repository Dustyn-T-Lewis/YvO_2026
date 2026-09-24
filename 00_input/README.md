# 00_input

The source data. No script writes here, so this directory has no `a_script/`, `b_reports/` or `c_data/`.

## Files

- `YvO_raw.xlsx`: protein intensities, 64 samples from 32 participants, Pre and Post. Read by `01_normalize.R`.
- `YvO_meta.xlsx`: one row per sample. `Col_ID` names the intensity column; `01_normalize.R` stops unless the two sets match exactly. Read by `01_normalize.R`, `03_DEP/a_script/supp/02`, `03` and `05`, F01 and F04.
- `YvO_pheno_calc.xlsx`: phenotypes as a formatted report. Read by `01_normalize.R`, which drops the summary rows and stops unless 32 participants remain, and by `03_DEP/a_script/supp/05`.
- `HPA_skeletal_muscle_annotations.tsv`: Human Protein Atlas annotations for the tissue filter. Read by `01_normalize.R` and F05.
- `wgcna_reference_modules.csv`: module assignments from the submitted analysis. `F05/a_script/YvO_WGCNA_run.R` uses them only to keep module colours the same between the two networks.
- `parent_meta/NORE.xlsx`, `EAA.xlsx`, `PPS_older.xlsx`, `PPS_younger.xlsx`: the parent trials' records, for Table 1. Read by F01.

`Subject_ID` is not unique: ten values are shared by two participants. Scripts group samples by the `Col_ID` prefix instead, `sub("_(Pre|Post)$", "", Col_ID)`.

`01_normalize.R` deletes the `supplement_cohort` column before building the DAList, so the model never sees it. No parent trial spans both age groups, so supplement arm is confounded with age. `03_DEP/a_script/supp/02_supplement_covariate.R` reads the column back from this workbook to measure that.

# 00_input

Source data. No script writes here.

## Files

- `YvO_raw.xlsx`: protein intensities, 64 samples from 32 participants, Pre and Post. Read by `01_normalize.R`.
- `YvO_meta.xlsx`: one row per sample. `Col_ID` names the intensity column, and `01_normalize.R` stops unless the two match exactly. Read by `01_normalize.R`, `03_DEP/a_script/supp/02`, `03` and `05`, F01 and F04.
- `YvO_pheno_calc.xlsx`: phenotypes as a formatted report. Read by `01_normalize.R`, which drops summary rows and stops unless 32 participants remain, and `03_DEP/a_script/supp/05`.
- `HPA_skeletal_muscle_annotations.tsv`: Human Protein Atlas annotations for the tissue filter. Read by `01_normalize.R` and `F05/a_script/supp/panels/S6_C_compartment.R`.
- `wgcna_reference_modules.csv`: module assignments from the submitted analysis, used by `F05/a_script/YvO_WGCNA_run.R` only to keep module colours stable.
- `parent_trials.csv`: one row per parent trial and arm, with its age group and number of participants. `YvO_meta.xlsx` links each participant to a trial through `parent_study` and `parent_id`. The parent trials' own records stay with their study teams.

`Subject_ID` is not unique (ten values are shared by two participants), so scripts group samples by the `Col_ID` prefix, `sub("_(Pre|Post)$", "", Col_ID)`.

`01_normalize.R` drops `supplement_cohort` before building the DAList, so the model never sees it. No parent trial spans both age groups, so supplement arm is confounded with age; `03_DEP/a_script/supp/02_supplement_covariate.R` reads the column back to measure that.

# 00 · Input

Hand-entered and hand-exported source files. No script writes anything here.

```
YvO_raw.xlsx                         -> 01_normalize.R
YvO_meta.xlsx                        -> 01_normalize.R, 03_DEP/supp/02,03,05, F01, F04
YvO_pheno_calc.xlsx                  -> 01_normalize.R
HPA_skeletal_muscle_annotations.tsv  -> 01_normalize.R, F05/_supp_qc_compartment.R
wgcna_reference_modules.csv          -> F05/YvO_WGCNA_run.R
parent_meta/{NORE,EAA,PPS_older,PPS_younger}.xlsx -> F01/04_phenotype_table.R
```

## What comes out

`Col_ID` joins `YvO_meta.xlsx` to the intensity columns; `01_normalize.R:61` asserts the
two sets match exactly. `uniprot_id` keys annotation to data rows and stays the row name
through every later stage.

`Subject_ID` is not unique: ten values are shared by two participants. Subjects are
therefore grouped by the `Col_ID` prefix, `sub("_(Pre|Post)$", "", Col_ID)`, in
`01_run_dep.R`, `03_reversal_aging_fdr.R` and `05_supplement_sensitivity.R`.

`YvO_pheno_calc.xlsx` is a formatted report rather than a table. `01_normalize.R:331-336`
drops the `Mean`, `Std Dev`, `n-size`, repeated-header and `T-test` rows, then asserts 32
participants remain.

## Orderings that matter

`supplement_cohort` is deleted at `01_normalize.R:59`, before the DAList is built, so the
fitted model cannot see it. `03_DEP/supp/02_supplement_covariate.R` reads it back out of
this workbook directly.

No parent trial spans both age groups, so supplement arm is confounded with age by
recruitment. The design still fits at full rank, which is the trap: forcing supplement in
takes Aging from 278 proteins to 12 with the Ruple arms pooled, 84 with all six arms, and
6 with control and placebo collapsed. Measured in `02_supplement_forced_fit.csv`.

Adding a column to `YvO_meta.xlsx` makes it available to every stage: `01_normalize.R`
passes the sheet through unchanged apart from that one deletion.

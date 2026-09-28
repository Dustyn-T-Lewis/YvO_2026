# 01_normalization

Filters raw intensities, removes the consensus outliers and normalizes with cyclic loess.

## Reads

- `00_input/YvO_raw.xlsx`, `YvO_meta.xlsx`, `YvO_pheno_calc.xlsx`, `HPA_skeletal_muscle_annotations.tsv`

## Writes

- `c_data/01_normalization.xlsx`: S8 Table, 9 sheets
- `c_data/02_normalized.csv`: the normalized matrix, 2,106 proteins by 62 samples. Stages 02 and 03 read this text copy, not the RDS, so the doubles match between stages.
- `c_data/03_DAList_normalized.rds`: the same matrix as a proteoDA DAList, read by stages 02, 03 and F02
- `c_data/00_report_intermediates.rds`: read by `02_generate_reports.R`, `03_DEP/a_script/02_generate_reports.R` and F00
- `c_data/03_fraction_composition.csv`, `03_fraction_top_carriers.csv`: compartment shares of the extract by protein count and summed intensity
- `b_reports/`: proteoDA QC reports, not tracked

## Run

```sh
Rscript 01_normalization/a_script/01_normalize.R
Rscript 01_normalization/a_script/02_generate_reports.R
Rscript 01_normalization/a_script/03_fraction_composition.R
```

`01_normalize.R` runs first; the other two read its output.

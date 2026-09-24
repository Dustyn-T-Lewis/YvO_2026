# F01 · Figure 1, participants and training response

Reads `00_input/YvO_meta.xlsx` and produces Figure 1, the S2 Figure phenotype stack,
Table 1 and the S2 Table workbook.

```
00_input/YvO_meta.xlsx + parent_meta/*.xlsx
  _prepost_template.R    sourced once per measure with a cfg list
  01_main_panels.R       A training volume, B DXA lean mass, C VL thickness
  02_supp_panels.R       deadlift 1RM, type II fCSA, type I fCSA
  03_body_comp_panels.R  DXA fat mass, fat-to-lean ratio
  05_supp_composite.R    sources 02 + 03, stacks the five rows
  04_phenotype_table.R   Table 1A / 1B / 1C
  90_stitch_F01.R        workbook, cleanup
  -> b_reports/main/{pdf,png}/MAIN_F01_composite{,_single_col}.*
  -> c_data/F01_supplementary.xlsx
```

```sh
Rscript 04_Figures/F01/a_script/90_stitch_F01.R
```

## What comes out

Figure 1 in two widths from the same panels, a five-row 85 x 160 mm supplementary stack,
and the three Table 1 blocks as CSV. `F01_supplementary.xlsx` is 14 sheets.

Each pre/post panel is drawn by `_prepost_template.R`, which fits the 2 x 2 mixed ANOVA
with `rstatix::anova_test()`, writes a source CSV and a summary CSV, then assigns
`<prefix>_left`, `_right`, `_title` and `_subtitle` into the global environment for the
composites to retrieve.

## Orderings that matter

`90_stitch_F01.R` sources `05_supp_composite.R` first, which itself sources
`02_supp_panels.R` and `03_body_comp_panels.R`. The two panel scripts are therefore not
listed in the stitcher.

`04_phenotype_table.R` runs after every panel script, because Table 1B is assembled from
the summary CSVs they write rather than refitted. A table cell and a figure p-value come
from one calculation. It stops on any missing summary.

`cleanup_after_workbook()` deletes the panel and summary CSVs, so it runs after
`build_workbook()`.

## Cost

8.0 s.

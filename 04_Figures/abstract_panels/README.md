# abstract_panels · graphical abstract

Reads finished figure workbooks and the combined DEP table; writes two graphical-abstract
layouts and their source-data workbooks.

```
04_Figures/F01/c_data/F01_supplementary.xlsx  (Per_participant)
04_Figures/F04/c_data/F04_supplementary.xlsx  (panel_C_trajectory)
04_Figures/F06/c_data/F06_supplementary.xlsx  (module_grid_summary)
03_DEP/c_data/03_combined_results.csv
  panel_A_hypertrophy / _B_amplitude / _C_remodelling / _D_direction / _E_modules / _F_aging
  _common.R              shared theme, palette, accessors
  99_assemble.R          2 x 4 grid, panels A B C D E -> GRID_abstract.{pdf,png}
  99_assemble_concise.R  single row, A C B F D E      -> GRID_abstract_concise.{pdf,png}
```

```sh
Rscript 04_Figures/abstract_panels/a_script/99_assemble.R
Rscript 04_Figures/abstract_panels/a_script/99_assemble_concise.R
```

## What comes out

Two plates with `c_data/abstract_panels.xlsx` and `c_data/abstract_panels_concise.xlsx`.
Each workbook carries an Overview naming the script and the counting basis behind each
panel, plus one data sheet per panel. PDFs are drawn on a transparent ground through a
device chosen once at load.

## Orderings that matter

Both assemblers run after the F01, F04 and F06 stitchers, because the panels read those
workbooks. `run_all.R` places them last. F06 deletes its own intermediates, so
`module_grid_summary` is only available from `F06_supplementary.xlsx`.

Each assembler sources its panel scripts into a fresh environment with `sys.source()` and
takes `plots`, `WIDTH` and `panel_data` back out, so no panel inherits another's objects.

`_common.R` reads `shared/style.R` as text rather than sourcing it, because sourcing
executes `devices.R` and defines the size globals `print_scale_apply.R` mutates.
`assert_style()` checks the local palette against `style.R` at the top of each assembler.

The two layouts do not draw the same set: `panel_F_aging` appears only in the concise row.

## Cost

3.8 s and 3.7 s.

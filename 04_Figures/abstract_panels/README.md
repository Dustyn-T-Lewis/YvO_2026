# abstract_panels

Draws the cards behind the graphical abstract, in the 2 x 4 grid the published abstract was built from and in a one-row backup layout.

## Reads

- `03_DEP/c_data/03_combined_results.csv`
- `04_Figures/F01/c_data/F01_data.xlsx`, sheet `Per_participant`
- `04_Figures/F04/c_data/F04_data.xlsx`, sheet `panel_C_trajectory`
- `04_Figures/F06/c_data/F06_data.xlsx`, sheet `module_grid_summary`
- `04_Figures/shared/style.R` and `shared/comparison_panels/panel_E_rrho2.R`, read as text for the palettes and the RRHO colour ramp

## Writes

- `b_reports/abstract.pdf`, `abstract.png`: graphical abstract, cards A to E
- `b_reports/abstract_concise.pdf`, `abstract_concise.png`: the backup one-row layout, cards A C B F D E; not in the manuscript
- `b_reports/panels/`: each card on its own, named after its script
- `c_data/abstract_data.xlsx`: source data for the graphical abstract
- `c_data/abstract_concise_data.xlsx`: source data for the backup layout

## Run

```sh
Rscript 04_Figures/abstract_panels/a_script/abstract.R
Rscript 04_Figures/abstract_panels/a_script/abstract_concise.R
```

Each card script in `a_script/panels/` also runs on its own. `_common.R` holds the theme, palettes, layout helpers and data accessors all cards share, and is only sourced. It reads `shared/style.R` as text and does not source it, because sourcing runs `devices.R` and defines size globals that `print_scale_apply.R` changes in place. `assert_style()` checks the local palette against `style.R` at the top of each composite.

## Order

F01, F04 and F06 run first, because the cards read their workbooks. The two composites can run in either order. Each writes its own workbook, so there is no separate data script.

## Outputs and manuscript items

`Figure_abstract.pdf` in the manuscript package was finished by hand in BioRender from `abstract.pdf`. Card F appears only in the backup layout.

| File | Manuscript item |
|---|---|
| `b_reports/abstract.pdf` | Graphical abstract, Key Results row |
| `b_reports/panels/A_hypertrophy.pdf` | Graphical abstract, VL thickness card |
| `b_reports/panels/B_amplitude.pdf` | Graphical abstract, protein change magnitude card |
| `b_reports/panels/C_remodelling.pdf` | Graphical abstract, proteins altered at FDR and Π cards |
| `b_reports/panels/D_direction.pdf` | Graphical abstract, RRHO map and FDR scatter cards |
| `b_reports/panels/E_modules.pdf` | Graphical abstract, module AUC cards |
| `b_reports/panels/F_aging.pdf` | none; backup layout only |
| `b_reports/abstract_concise.pdf` | none; backup layout |
| `c_data/abstract_data.xlsx`, sheet `Overview` | Graphical abstract, the counting basis of each card |
| `c_data/abstract_data.xlsx`, sheets `panel_A_hypertrophy` to `panel_E_modules` | Graphical abstract, source data for each card |
| `c_data/abstract_concise_data.xlsx` | none; backup layout source data |

The workbooks keep the sheet names and the `script` column they were first written with, so they name the cards `panel_A_hypertrophy` and so on.

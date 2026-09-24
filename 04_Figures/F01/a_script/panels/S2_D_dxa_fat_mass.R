#!/usr/bin/env Rscript
# S2 Figure D: whole-body DXA fat mass before and after training.

setwd(here::here())

set.seed(42)

cfg <- list(
  dv_col = "DXA_FM_kg", y_label = "DXA FM (kg)",
  delta_label = expression(bold("Δ DXA FM (kg)")),
  title = "DXA Fat Mass", name = "S2_D_dxa_fat_mass",
  audit_file = "panel_A_dxa_fat_mass.csv", dat = "04_Figures/F01/c_data/supp",
  filter_complete = TRUE
)
source("04_Figures/F01/a_script/panels/_prepost.R", local = TRUE)

invisible(p_combo)

#!/usr/bin/env Rscript
# Figure 1B: DXA lean body mass before and after training, and the change.

setwd(here::here())

cfg <- list(
  dv_col = "DXA_LBM_kg", y_label = "DXA LBM (kg)",
  delta_label = expression(bold("Δ DXA LBM (kg)")),
  title = "DXA Lean Body Mass", name = "B_dxa_lbm",
  audit_file = "panel_B_dxa_lbm.csv", dat = "04_Figures/F01/c_data"
)
source("04_Figures/F01/a_script/_prepost.R", local = TRUE)

invisible(p_combo)

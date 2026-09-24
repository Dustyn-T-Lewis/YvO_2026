#!/usr/bin/env Rscript
# S2 Figure B: type II fibre cross-sectional area before and after training.

setwd(here::here())

cfg <- list(
  dv_col = "Type_II_fCSA",
  y_label = expression(bold("Type II fCSA (µm"^2 * ")")),
  delta_label = expression(bold("Δ Type II fCSA (µm"^2 * ")")),
  title = "Type II Fiber CSA", name = "S2_B_type_II_fcsa",
  audit_file = "panel_B_type_II_fcsa.csv", dat = "04_Figures/F01/c_data/supp",
  coerce_cols = TRUE, filter_complete = TRUE
)
source("04_Figures/F01/a_script/panels/_prepost.R", local = TRUE)

invisible(p_combo)

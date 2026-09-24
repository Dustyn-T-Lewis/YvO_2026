#!/usr/bin/env Rscript
# S2 Figure C: type I fibre cross-sectional area before and after training.

setwd(here::here())

cfg <- list(
  dv_col = "Type_I_fCSA",
  y_label = expression(bold("Type I fCSA (µm"^2 * ")")),
  delta_label = expression(bold("Δ Type I fCSA (µm"^2 * ")")),
  title = "Type I Fiber CSA", name = "S2_C_type_I_fcsa",
  audit_file = "panel_C_type_I_fcsa.csv", dat = "04_Figures/F01/c_data/supp",
  coerce_cols = TRUE, filter_complete = TRUE
)
source("04_Figures/F01/a_script/panels/_prepost.R", local = TRUE)

invisible(p_combo)

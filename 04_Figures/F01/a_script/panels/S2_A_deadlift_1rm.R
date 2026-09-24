#!/usr/bin/env Rscript
# S2 Figure A: deadlift one-repetition maximum before and after training.

setwd(here::here())

set.seed(42)

cfg <- list(
  dv_col = "deadlift_1rm_kg", y_label = "Deadlift 1RM (kg)",
  delta_label = expression(bold("Δ 1RM (kg)")),
  title = "Deadlift 1RM", name = "S2_A_deadlift_1rm",
  audit_file = "panel_A_deadlift_1rm.csv", dat = "04_Figures/F01/c_data/supp",
  coerce_cols = TRUE, filter_complete = TRUE
)
source("04_Figures/F01/a_script/panels/_prepost.R", local = TRUE)

invisible(p_combo)

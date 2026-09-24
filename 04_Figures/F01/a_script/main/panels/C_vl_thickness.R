#!/usr/bin/env Rscript
# Figure 1C: vastus lateralis thickness before and after training, and the change.

setwd(here::here())

cfg <- list(
  dv_col = "VL_thick_cm", y_label = "VL thickness (cm)",
  delta_label = expression(bold("Δ VL thickness (cm)")),
  title = "VL Thickness", name = "C_vl_thickness",
  audit_file = "panel_C_vl_thickness.csv", dat = "04_Figures/F01/c_data",
  y_breaks = c(0, 0.5, 1.0), y_labels = c("0", ".5", "1")
)
source("04_Figures/F01/a_script/_prepost.R", local = TRUE)

invisible(p_combo)

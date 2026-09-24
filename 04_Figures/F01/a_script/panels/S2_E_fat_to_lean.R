#!/usr/bin/env Rscript
# S2 Figure E: fat-to-lean mass ratio before and after training.
#
# Fat-to-lean is used rather than percent body fat because NORE recorded body
# mass only at baseline, so any percent-of-weight measure at post-training
# rests on a derived denominator for eleven of the 32 participants. The ratio
# needs only the two DXA compartments, both measured at both timepoints.

setwd(here::here())

pacman::p_load(dplyr)

cfg <- list(
  dv_col = "FM_to_LBM", y_label = "FM / LBM",
  delta_label = expression(bold("Δ FM / LBM")),
  title = "Fat-to-Lean Ratio", name = "S2_E_fat_to_lean",
  audit_file = "panel_B_fat_to_lean.csv", dat = "04_Figures/F01/c_data/supp",
  filter_complete = TRUE,
  derive = function(m) mutate(m, FM_to_LBM = DXA_FM_kg / DXA_LBM_kg)
)
source("04_Figures/F01/a_script/panels/_prepost.R", local = TRUE)

invisible(p_combo)

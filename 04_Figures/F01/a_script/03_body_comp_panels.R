#!/usr/bin/env Rscript
# F01 Supp — DXA Fat Mass, Fat-to-Lean Ratio.
# Panels only; 05_supp_composite.R stacks them with the strength/CSA trio.
#
# Fat-to-lean is used rather than percent body fat because NORE recorded body
# mass only at baseline, so any percent-of-weight measure at post-training
# rests on a derived denominator for eleven of the 32 participants. The ratio
# needs only the two DXA compartments, both measured at both timepoints.

withr::local_dir(here::here())

pacman::p_load(withr, patchwork, dplyr)

source("04_Figures/shared/style.R")

set.seed(42)

BASE <- "04_Figures/F01"
RPT_PNG <- file.path(BASE, "b_reports", "supp", "png")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "pdf")
PNL_PNG <- file.path(RPT_PNG, "panels")
PNL_PDF <- file.path(RPT_PDF, "panels")
DAT <- file.path(BASE, "c_data", "supp")
for (d in c(PNL_PNG, PNL_PDF, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

TMPL <- "04_Figures/F01/a_script/_prepost_template.R"

cfg <- list(
  dv_col = "DXA_FM_kg", y_label = "DXA FM (kg)",
  delta_label = expression(bold("Δ DXA FM (kg)")),
  title = "DXA Fat Mass", output_prefix = "pFA",
  file_tag = "panel_A_dxa_fat_mass", audit_file = "panel_A_dxa_fat_mass.csv",
  file_prefix = "SUPP", rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT,
  filter_complete = TRUE
)
source(TMPL)

cfg <- list(
  dv_col = "FM_to_LBM", y_label = "FM / LBM",
  delta_label = expression(bold("Δ FM / LBM")),
  title = "Fat-to-Lean Ratio", output_prefix = "pFB",
  file_tag = "panel_B_fat_to_lean", audit_file = "panel_B_fat_to_lean.csv",
  file_prefix = "SUPP", rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT,
  filter_complete = TRUE,
  derive = function(m) mutate(m, FM_to_LBM = DXA_FM_kg / DXA_LBM_kg)
)
source(TMPL)

message("F01 body-composition panels done")

#!/usr/bin/env Rscript
# F01 Supp — Deadlift 1RM, Type II fCSA, Type I fCSA.
# Panels only; 05_supp_composite.R stacks them with the body-composition pair.

withr::local_dir(here::here())

pacman::p_load(withr, patchwork)

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
  dv_col = "deadlift_1rm_kg", y_label = "Deadlift 1RM (kg)",
  delta_label = expression(bold("Δ 1RM (kg)")),
  title = "Deadlift 1RM", output_prefix = "pSA",
  file_tag = "panel_A_deadlift_1rm", audit_file = "panel_A_deadlift_1rm.csv",
  file_prefix = "SUPP", rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT,
  coerce_cols = TRUE, filter_complete = TRUE
)
source(TMPL)

cfg <- list(
  dv_col = "Type_II_fCSA",
  y_label = expression(bold("Type II fCSA (µm"^2 * ")")),
  delta_label = expression(bold("Δ Type II fCSA (µm"^2 * ")")),
  title = "Type II Fiber CSA", output_prefix = "pSB",
  file_tag = "panel_B_type_II_fcsa", audit_file = "panel_B_type_II_fcsa.csv",
  file_prefix = "SUPP", rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT,
  coerce_cols = TRUE, filter_complete = TRUE
)
source(TMPL)

cfg <- list(
  dv_col = "Type_I_fCSA",
  y_label = expression(bold("Type I fCSA (µm"^2 * ")")),
  delta_label = expression(bold("Δ Type I fCSA (µm"^2 * ")")),
  title = "Type I Fiber CSA", output_prefix = "pSC",
  file_tag = "panel_C_type_I_fcsa", audit_file = "panel_C_type_I_fcsa.csv",
  file_prefix = "SUPP", rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT,
  coerce_cols = TRUE, filter_complete = TRUE
)
source(TMPL)

message("F01 strength/CSA panels done")

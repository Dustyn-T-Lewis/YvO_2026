#!/usr/bin/env Rscript
# F04 Supplementary: Aging x Training(Old) rank-rank overlap, descriptive only.
#
# Aging (Old_Pre - Young_Pre) and Training(Old) (Old_Post - Old_Pre) share
# Old_Pre, so their per-protein statistics are coupled and the corner
# hypergeometric test that panel B reports for Training(Young) x Training(Old)
# is not valid here. The map is shown for its shape; the reversal inference
# lives in F04 panel E at pathway level, where the coupling averages down.
#
# Sourced by 02_supp_panels.R after style.R.

BASE <- "04_Figures/F04"
cfg <- list(
  fig_id = "F04 supp",
  descriptive = TRUE,
  file_stem = "SUPP_rrho2_aging_training_old",
  # Standalone supplement, so the title and subtitle have to fit the canvas.
  panel_w = 150,
  t_col_1 = "t_Aging",
  t_col_2 = "t_Training_Old",
  rrho_labels = c("Aging", "Training (Old)"),
  title = "Aging \u00d7 Training (Old): RRHO2, descriptive",
  subtitle_fmt = paste0(
    "%d shared proteins | the contrasts share Old_Pre, so no overlap test is reported\n",
    "reversal is tested at pathway level in F04 panel E"
  ),
  axis_label_1 = expression("Aging rank" ~ (Down %->% Up)),
  axis_label_2 = expression("Training (Old) rank" ~ (Down %->% Up)),
  quadrant_labels = list(
    UU = "Exacerbated\n(A up, T up)", DD = "Exacerbated\n(A dn, T dn)",
    UD = "Reversed\n(A up, T dn)", DU = "Reversed\n(A dn, T up)"
  ),
  hotspot_export_names = list(UU = "UU", DD = "DD", UD = "UD", DU = "DU"),
  ora_quadrant_names = list(
    UU = "Exacerbated Up", DD = "Exacerbated Down",
    UD = "Reversed (A Up / T Down)", DU = "Reversed (A Down / T Up)"
  ),
  ora_grouped = list(
    file_1_quads = c("ora_UD", "ora_DU"),
    file_2_quads = c("ora_UU", "ora_DD")
  ),
  summary_quadrant_names = list(
    UU = "Exacerbated_Up", UU_slug = "exacerbated_up",
    DD = "Exacerbated_Down", DD_slug = "exacerbated_down",
    UD = "Reversed_AgingUp_TrainingDown", UD_slug = "reversed_a_up",
    DU = "Reversed_AgingDown_TrainingUp", DU_slug = "reversed_a_down"
  ),
  rpt_png = file.path(BASE, "b_reports", "supp", "png", "panels"),
  rpt_pdf = file.path(BASE, "b_reports", "supp", "pdf", "panels"),
  dat = file.path(BASE, "c_data", "panel_supp"),
  supp = NULL
)
source("04_Figures/shared/comparison_panels/panel_E_rrho2.R")

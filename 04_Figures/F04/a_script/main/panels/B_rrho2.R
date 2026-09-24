#!/usr/bin/env Rscript
# Figure 4B: RRHO2 threshold-free rank overlap, Training (Young) x Training
# (Old). Returns the plot with the corner-depth note F04.R prints as the
# subtitle attached as attr(, "sweep_txt").

setwd(here::here())

pacman::p_load(
  dplyr, tidyr, tibble, stringr, readr, ggplot2, patchwork, cowplot, readxl
)

source("04_Figures/shared/style.R")

PNL <- "04_Figures/F04/b_reports/main/panels"
DAT <- "04_Figures/F04/c_data"

cfg <- list(
  fig_id = "F04",
  file_stem = "B_rrho2",
  t_col_1 = "t_Training_Young",
  t_col_2 = "t_Training_Old",
  rrho_labels = c("Training (Young)", "Training (Old)"),
  title = "Threshold-Free Concordance (RRHO2)",
  subtitle_fmt = "%s | %d shared proteins | warm corners = concordant regulation",
  axis_label_1 = expression("Training (Young) rank" ~ (Down %->% Up)),
  axis_label_2 = expression("Training (Old) rank" ~ (Down %->% Up)),
  # Position says which discordance is which, the way panels D and E already
  # label both their off-diagonal corners.
  quadrant_labels = list(
    UU = "Concordant Up", DD = "Concordant Down",
    UD = "Discordant", DU = "Discordant"
  ),
  hotspot_export_names = list(UU = "UU", DD = "DD", UD = "UD", DU = "DU"),
  ora_min_size = 15,
  ora_quadrant_names = list(
    UU = "Concordant Up", DD = "Concordant Down",
    UD = "Discordant (Y Up / O Down)", DU = "Discordant (Y Down / O Up)"
  ),
  ora_grouped = list(
    file_1_quads = c("ora_UU", "ora_DD"),
    file_2_quads = c("ora_UD", "ora_DU"),
    note_if_empty_2 = paste(
      "No pathway reached FDR < 0.05 in either discordant quadrant.",
      "The quadrants are small (18 and 7 proteins), so this is an absence of",
      "power rather than evidence that the discordant proteins share no biology."
    )
  ),
  ora_colors = ORA_QUAD_COLORS_F2,
  summary_quadrant_names = list(
    UU = "Concordant_Up", UU_slug = "concordant_up",
    DD = "Concordant_Down", DD_slug = "concordant_down",
    UD = "Discordant_YoungUp_OldDown", UD_slug = "discordant_y_up",
    DU = "Discordant_YoungDown_OldUp", DU_slug = "discordant_y_down"
  ),
  rpt_png = PNL, rpt_pdf = PNL, dat = DAT,
  # Quadrant ORA needs all four corners populated; the discordant pair is
  # effectively empty here (UD = 0, DU = 1), so the supplementary block is off.
  supp = NULL
)
source("04_Figures/shared/comparison_panels/panel_E_rrho2.R", local = TRUE)

attr(pE_heat, "sweep_txt") <- sweep_txt
invisible(pE_heat)

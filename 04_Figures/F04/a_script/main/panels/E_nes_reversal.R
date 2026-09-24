#!/usr/bin/env Rscript
# Figure 4E: pathway NES scatter, Aging x Training (Old). Returns the plot with
# rho and the reversed fraction attached as attr(, "stats").

setwd(here::here())

pacman::p_load(
  dplyr, tidyr, tibble, stringr, readr, ggplot2, patchwork, cowplot, readxl
)

source("04_Figures/shared/style.R")
source("04_Figures/F04/a_script/main/panels/_nes_scatter.R", local = TRUE)

cfg <- nes_scatter_cfg(
  "Aging", "Training_Old", "Pathway-Level Reversal (fGSEA)",
  "Aging", "Training Old", "reversed",
  "negative ρ = training opposes aging effects", -1,
  SIG_COLORS_F3,
  c("Sig Training only", "Sig Aging only", "Sig Both"),
  list(
    sig_both_label = "Sig Both", sig_x_label = "Sig Aging only",
    sig_y_label = "Sig Training only",
    bg_blue_1 = c(0, Inf, -Inf, 0), bg_blue_2 = c(-Inf, 0, 0, Inf),
    bg_red_1 = c(0, Inf, 0, Inf), bg_red_2 = c(-Inf, 0, -Inf, 0),
    label_tr = "Exacerbated", color_tr = "#D6604D",
    label_tl = "Reversed", color_tl = "#4393C3",
    label_bl = "Exacerbated", color_bl = "#D6604D",
    label_br = "Reversed", color_br = "#4393C3",
    metric_count_fn = function(q1, q2, q3, q4) q2 + q4
  ),
  nudges = NUDGE_E
)
source("04_Figures/shared/comparison_panels/panel_D_nes_scatter.R", local = TRUE)
rename_nes_outputs("E_nes_reversal", "nes_reversal.csv")

attr(pD, "stats") <- list(rho = as.numeric(nes_cor_all$estimate), frac = pw_conc_frac)
invisible(pD)

#!/usr/bin/env Rscript
# Figure 4D: pathway NES scatter, Training (Young) x Training (Old). Returns
# the plot, legend kept, with rho and the concordant fraction attached as
# attr(, "stats").

setwd(here::here())

pacman::p_load(
  dplyr, tidyr, tibble, stringr, readr, ggplot2, patchwork, cowplot, readxl
)

source("04_Figures/shared/style.R")
source("04_Figures/F04/a_script/main/panels/_nes_scatter.R", local = TRUE)

cfg <- nes_scatter_cfg(
  "Training_Young", "Training_Old", "Pathway-Level Concordance (fGSEA)",
  "Training Young", "Training Old", "concordant",
  "positive ρ = shared pathway regulation across age groups", 1,
  SIG_COLORS_F2,
  c("Sig Old only", "Sig Young only", "Sig Both", "Interaction"),
  list(
    sig_both_label = "Sig Both", sig_x_label = "Sig Young only",
    sig_y_label = "Sig Old only",
    bg_red_1 = c(0, Inf, 0, Inf), bg_red_2 = c(-Inf, 0, -Inf, 0),
    bg_blue_1 = c(0, Inf, -Inf, 0), bg_blue_2 = c(-Inf, 0, 0, Inf),
    label_tr = "Concordant Up", color_tr = "#D6604D",
    label_tl = "Discordant", color_tl = "#4393C3",
    label_bl = "Concordant Down", color_bl = "#D6604D",
    label_br = "Discordant", color_br = "#4393C3",
    metric_count_fn = function(q1, q2, q3, q4) q1 + q3
  ),
  nudges = NUDGE_D, seed = 7
)
source("04_Figures/shared/comparison_panels/panel_D_nes_scatter.R", local = TRUE)
rename_nes_outputs("D_nes_concordance", "nes_concordance.csv")

attr(pD, "stats") <- list(rho = as.numeric(nes_cor_all$estimate), frac = pw_conc_frac)
invisible(pD)

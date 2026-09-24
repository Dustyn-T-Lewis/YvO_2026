#!/usr/bin/env Rscript
# S4a Figure B: Pi-score distribution, all four contrasts.

setwd(here::here())
source("04_Figures/F03/a_script/supp/panels/_distribution.R", local = TRUE)

dist_panel("pi_score", "#E05A4E", "\u03A0-score", 0.05,
           "\u03A0 < 0.05: %d", "\u03A0-score distribution", "b", "S4a_B_pi_score")

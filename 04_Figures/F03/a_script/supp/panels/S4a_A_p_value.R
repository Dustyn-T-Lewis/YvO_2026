#!/usr/bin/env Rscript
# S4a Figure A: raw p-value distribution, all four contrasts.

setwd(here::here())
source("04_Figures/F03/a_script/supp/panels/_distribution.R", local = TRUE)

dist_panel("P.Value", "#5DA5DA", "p_value", NULL,
           "p < 0.05: %d", "Raw p-value distribution", "a", "S4a_A_p_value")

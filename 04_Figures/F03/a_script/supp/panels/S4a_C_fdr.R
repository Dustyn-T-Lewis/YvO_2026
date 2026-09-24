#!/usr/bin/env Rscript
# S4a Figure C: BH-adjusted p-value distribution, all four contrasts.

setwd(here::here())
source("04_Figures/F03/a_script/supp/panels/_distribution.R", local = TRUE)

dist_panel("adj.P.Val", "#9B7FBF", "FDR (BH)", 0.05,
           "FDR < 0.05: %d", "FDR distribution", "c", "S4a_C_fdr")

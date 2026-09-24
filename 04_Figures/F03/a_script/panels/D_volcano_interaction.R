#!/usr/bin/env Rscript
# Figure 3D: volcano ring, age by training interaction.

setwd(here::here())
source("04_Figures/F03/a_script/panels/_volcano.R", local = TRUE)

volcano_panel("Interaction", "Age \u00d7 Training Interaction", "Training_Old \u2212 Training_Young", "D", "D_volcano_interaction")

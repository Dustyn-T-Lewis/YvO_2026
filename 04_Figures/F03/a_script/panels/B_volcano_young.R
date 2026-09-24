#!/usr/bin/env Rscript
# Figure 3B: volcano ring, training in younger adults.

setwd(here::here())
source("04_Figures/F03/a_script/panels/_volcano.R", local = TRUE)

volcano_panel("Training_Young", "Training Response (Young)", "Young_Post \u2212 Young_Pre", "B", "B_volcano_young")

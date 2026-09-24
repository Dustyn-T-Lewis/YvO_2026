#!/usr/bin/env Rscript
# Figure 3C: volcano ring, training in older adults.

setwd(here::here())
source("04_Figures/F03/a_script/panels/_volcano.R", local = TRUE)

volcano_panel("Training_Old", "Training Response (Old)", "Old_Post \u2212 Old_Pre", "C", "C_volcano_old")

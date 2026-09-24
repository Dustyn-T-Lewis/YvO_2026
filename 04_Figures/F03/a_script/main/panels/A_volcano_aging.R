#!/usr/bin/env Rscript
# Figure 3A: volcano ring, aging.

setwd(here::here())
source("04_Figures/F03/a_script/main/panels/_volcano.R", local = TRUE)

volcano_panel("Aging", "Aging Effect", "Old_Pre \u2212 Young_Pre", "A", "A_volcano_aging")

#!/usr/bin/env Rscript
# S8 Figure C: yellow module triptych, protein z-scores, eigengene
# trajectories and top ORA terms.

setwd(here::here())
source("04_Figures/F05/a_script/supp/panels/_triptych.R", local = TRUE)

triptych_panel("yellow", "S8_C_yellow")

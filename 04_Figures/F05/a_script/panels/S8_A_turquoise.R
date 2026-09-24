#!/usr/bin/env Rscript
# S8 Figure A: turquoise module triptych, protein z-scores, eigengene
# trajectories and top ORA terms.

setwd(here::here())
source("04_Figures/F05/a_script/panels/_triptych.R", local = TRUE)

triptych_panel("turquoise", "S8_A_turquoise")

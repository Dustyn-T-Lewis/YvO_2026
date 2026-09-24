#!/usr/bin/env Rscript
# Figure 6A: per-module ROCs for the modules that separate younger from older.

setwd(here::here())
source("04_Figures/F06/a_script/main/panels/_roc_grid.R", local = TRUE)

grid_age <- build_grid(age_cells)
save_panel(grid_age, "A_age_roc")

invisible(grid_age)

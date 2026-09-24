#!/usr/bin/env Rscript
# Figure 6B: per-module ROCs for pre against post training, younger over older.

setwd(here::here())
source("04_Figures/F06/a_script/main/panels/_roc_grid.R", local = TRUE)

grid_train <- build_grid(train_cells)
save_panel(grid_train, "B_training_roc")

invisible(grid_train)

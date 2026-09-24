#!/usr/bin/env Rscript
# F06 — Phenotype Prediction: Master Orchestrator
# Sources main panels (all data-generating scripts + composite + xlsx),
# then supp panels (reads pre-rendered PNGs, builds composites).

setwd(here::here())

source("04_Figures/shared/figure_supplement_helpers.R")

DAT <- "04_Figures/F06/c_data"

message("F06: Running main panels + composite + xlsx")
source("04_Figures/F06/a_script/01_main_panels.R")

message("F06: Running supp composite")
source("04_Figures/F06/a_script/02_supp_panels.R")

# Final cleanup: remove any leftover CSVs
remaining <- list.files(DAT, pattern = "\\.csv$",
                        recursive = TRUE, full.names = TRUE)
if (length(remaining)) {
  file.remove(remaining)
  message(sprintf("  final cleanup: removed %d leftover CSV(s)", length(remaining)))
}

message("F06 complete")

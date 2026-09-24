#!/usr/bin/env Rscript
# F00 — Pipeline QC: Orchestrator
# Sources the supp composite script, copies outputs to Box
# F00 has no main figure — produces 2 supp composites:
#   SUPP_F00_normalization (panels A–G)
#   SUPP_F00_imputation    (panels H–N)

setwd(here::here())

source("04_Figures/F00/a_script/01_supp_panels.R")

message("F00 complete")

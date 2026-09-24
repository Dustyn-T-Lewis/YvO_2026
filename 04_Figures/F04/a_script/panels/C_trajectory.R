#!/usr/bin/env Rscript
# Figure 4C: response magnitude retained in older adults, per protein set.

setwd(here::here())

pacman::p_load(
  dplyr, tidyr, tibble, stringr, readr, ggplot2, patchwork, cowplot, readxl
)

source("04_Figures/shared/style.R")
source("04_Figures/F04/a_script/panels/_engine_cfg.R", local = TRUE)

PNL <- "04_Figures/F04/b_reports/panels"
DAT <- "04_Figures/F04/c_data"

cfg <- c(engine_cfg, list(
  rpt_png = PNL, rpt_pdf = PNL, dat = DAT,
  panel_w = 102, panel_h = 140
))
source("04_Figures/shared/comparison_panels/panel_C_trajectory.R", local = TRUE)

# The engine writes a fixed file name.
for (ext in c("png", "pdf")) {
  file.rename(
    file.path(PNL, paste0("MAIN_panel_C_trajectory.", ext)),
    file.path(PNL, paste0("C_trajectory.", ext))
  )
}

invisible(pC_trajectory)

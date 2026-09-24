#!/usr/bin/env Rscript
# Figure 4F: younger-adult FDR and Pi sets against the Training (Old) ranking,
# fry rotation test drawn as a barcode.

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
source("04_Figures/shared/comparison_panels/panel_fry_barcode.R", local = TRUE)

# The engine writes a fixed file name.
for (ext in c("png", "pdf")) {
  file.rename(
    file.path(PNL, paste0("MAIN_panel_F_fry_barcode.", ext)),
    file.path(PNL, paste0("F_fry_barcode.", ext))
  )
}

invisible(pF_barcode)

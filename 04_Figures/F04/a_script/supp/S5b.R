#!/usr/bin/env Rscript
# S5b Figure: protein-level training response across age groups, one panel.

setwd(here::here())

pacman::p_load(grid, gridExtra)

source("04_Figures/shared/style.R")

PANELS <- "04_Figures/F04/a_script/supp/panels"
RPT <- "04_Figures/F04/b_reports/supp"

p <- source_panel(file.path(PANELS, "S5b_young_dep_heatmap.R"))
fig_w_mm <- attr(p, "width_mm")
fig_h_mm <- attr(p, "height_mm")

png(file.path(RPT, "S5b.png"),
  width = fig_w_mm, height = fig_h_mm, units = "mm", res = 300
)
grid.newpage()
grid.draw(p)
dev.off()

pdf_dev <- get_pdf_device()
if (is.character(pdf_dev)) pdf_dev <- match.fun(pdf_dev)
pdf_dev(file.path(RPT, "S5b.pdf"),
  width = fig_w_mm / 25.4, height = fig_h_mm / 25.4
)
grid.newpage()
grid.draw(p)
dev.off()

caption_supp(p, "S5b", fig_w_mm, fig_h_mm, RPT)

message("S5b done")

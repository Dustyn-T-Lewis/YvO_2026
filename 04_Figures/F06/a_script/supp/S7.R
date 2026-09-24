#!/usr/bin/env Rscript
# S7 Figure: A the per-module ROC grid, B the full module-phenotype sweep.
#
# Runs both panel scripts, then reads their PNG renders from b_reports/panels/
# and composites them side by side.
#
# Three panels were dropped: the two leave-one-subject-out sensitivity plots
# and the multivariate classifier decomposition. Each plotted a handful of
# AUCs that ship as loso_auc_summary, loso_wgcna_refit_summary and
# panel_A_classifier_auc in F06_data.xlsx, where they read as a
# table. Their scripts now compute those sheets and draw nothing.

setwd(here::here())

source("04_Figures/shared/style.R")

pacman::p_load(patchwork, cowplot, png, grid)

PANELS <- "04_Figures/F06/a_script/supp/panels"
source_panel(file.path(PANELS, "S7_A_module_grid.R"))
# B runs C_hero_grid.R first when the 180-test screen is missing.
source_panel(file.path(PANELS, "S7_B_full_sweep.R"))

RPT        <- "04_Figures/F06/b_reports/supp"
RPT_PANELS <- file.path(RPT, "panels")

read_panel <- function(file) {
  path <- file.path(RPT_PANELS, file)
  if (!file.exists(path)) stop("Missing: ", path)
  img <- readPNG(path)
  list(grob = rasterGrob(img, interpolate = TRUE),
       aspect = dim(img)[2] / dim(img)[1])   # width/height
}

# Page 1 panels (A + B, aspect-matched)
pA <- read_panel("S7_A_module_grid.png")
pB <- read_panel("S7_B_full_sweep.png")

COMP_H <- 166   # mm tall
wA_mm  <- COMP_H * pA$aspect
wB_mm  <- COMP_H * pB$aspect
COMP_W <- wA_mm + wB_mm

TAG_SZ <- composite_text_sizes(COMP_W, 178)$tag

page1 <- (wrap_elements(full = pA$grob) |
           wrap_elements(full = pB$grob)) +
  plot_layout(widths = c(wA_mm, wB_mm)) &
  theme(plot.margin = margin(0, 0, 0, 0))

TAG_X_A <- 5 / COMP_W
TAG_X_B <- (wA_mm + 5) / COMP_W
TAG_Y   <- 0.992

page1_final <- ggdraw(page1 & theme(plot.margin = margin(2, 2, 2, 2))) +
  draw_label("A", x = TAG_X_A, y = TAG_Y, size = TAG_SZ,
             fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("B", x = TAG_X_B, y = TAG_Y, size = TAG_SZ,
             fontface = "bold", hjust = 0, vjust = 1)

# Write outputs
graphics.off()
pdf_device <- get_raster_pdf_device()

ggsave(file.path(RPT, "S7.pdf"), page1_final,
       width = COMP_W, height = COMP_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
embed_pdf_fonts(file.path(RPT, "S7.pdf"))

ggsave(file.path(RPT, "S7.png"), page1_final,
       width = COMP_W, height = COMP_H, units = "mm",
       dpi = 300, limitsize = FALSE)

caption_supp(page1_final, "S7", COMP_W, COMP_H, RPT)

message(sprintf("Wrote S7 (%.0f x %.0f mm)", COMP_W, COMP_H))

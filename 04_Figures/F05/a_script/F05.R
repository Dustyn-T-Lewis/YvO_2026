#!/usr/bin/env Rscript
# Figure 5: the module-trait heatmap and the module NES scatters, placed as
# rasters on one canvas and cropped.

setwd(here::here())

source("04_Figures/shared/style.R")

pacman::p_load(patchwork, cowplot, png, grid, dplyr, tidyr)

BASE <- "04_Figures/F05"

PANELS <- "04_Figures/F05/a_script/panels"
source_panel(file.path(PANELS, "A_module_trait_heatmap.R"))
source_panel(file.path(PANELS, "B_nes_scatters.R"))

RPT_PDF <- file.path(BASE, "b_reports")
RPT_PNG <- file.path(BASE, "b_reports")

PANEL_DIR <- file.path(BASE, "b_reports", "panels")
read_panel <- function(file, dir = PANEL_DIR) {
  path <- file.path(dir, file)
  if (!file.exists(path)) stop("Missing: ", path)
  rasterGrob(readPNG(path), interpolate = TRUE)
}

pA_grob <- read_panel("A_module_trait_heatmap.png")
pB_grob <- read_panel("B_nes_scatters.png")
pB_leg_grob <- read_panel("B_nes_scatters_legend.png")

layout_cfg <- list(
  comp_w = 470, comp_h = 300,
  a_title_y = 0.835,
  a_legend_y = 0.14,
  # Panel B's box is solved, not chosen, so the scatters' outer panel borders
  # land on the heatmap's top and bottom tile borders. Both panels are rasters
  # fitted to their box with the aspect preserved, so the heatmap borders sit at
  # canvas y 0.751686 and 0.258944, and the scatter borders at 0.056307 and
  # 0.990219 of panel B's own height. Re-measure both if either panel's internal
  # spacing changes -- panel B's axis type sets where its borders fall.
  grid_top = 0.756847,
  grid_bot = 0.229236,
  b_x = 0.43,
  b_w = 0.55,
  b_title_x = 0.636,
  # Panel B's legend is clipped: its three columns are wider than the gap
  # between the heatmap and the crop edge, so the first column loses a few
  # characters. Widening b_leg_w only moves the clipping to the right-hand
  # column; the fix is to re-flow the legend, which is layout surgery.
  b_leg_x = 0.70, b_leg_y = 0.195, b_leg_w = 0.24, b_leg_h = 0.04,
  a_grob_x = 0.01, a_grob_y = 0, a_grob_w = 0.60, a_grob_h = 0.96,
  a_let_x = 15, a_let_y = 248,
  a_ttl_x = 58, a_ttl_y = 248,
  a_sub_x = 58, a_sub_y = 243,
  b_let_x = 285, b_let_y = 248,
  b_ttl_x = 294, b_ttl_y = 248,
  b_sub_x = 294, b_sub_y = 243,
  crop_l = 0.01, crop_r = 0.80,
  crop_b = 0.17, crop_t = 0.85
)

COMP_W <- layout_cfg$comp_w
COMP_H <- layout_cfg$comp_h
A_TITLE_Y <- layout_cfg$a_title_y
A_LEGEND_Y <- layout_cfg$a_legend_y
GRID_TOP <- layout_cfg$grid_top
GRID_BOT <- layout_cfg$grid_bot

B_X <- layout_cfg$b_x
B_W <- layout_cfg$b_w
B_TITLE_X <- layout_cfg$b_title_x
B_GRID_H <- GRID_TOP - GRID_BOT

mm2x <- function(mm) mm / COMP_W
mm2y <- function(mm) mm / COMP_H

A_LET_X <- layout_cfg$a_let_x
A_LET_Y <- layout_cfg$a_let_y
A_TTL_X <- layout_cfg$a_ttl_x
A_TTL_Y <- layout_cfg$a_ttl_y
A_SUB_X <- layout_cfg$a_sub_x
A_SUB_Y <- layout_cfg$a_sub_y

B_LET_X <- layout_cfg$b_let_x
B_LET_Y <- layout_cfg$b_let_y
B_TTL_X <- layout_cfg$b_ttl_x
B_TTL_Y <- layout_cfg$b_ttl_y
B_SUB_X <- layout_cfg$b_sub_x
B_SUB_Y <- layout_cfg$b_sub_y

CROP_L <- layout_cfg$crop_l
CROP_R <- layout_cfg$crop_r
CROP_B <- layout_cfg$crop_b
CROP_T <- layout_cfg$crop_t
SAVE_W <- COMP_W * (CROP_R - CROP_L)
SAVE_H <- COMP_H * (CROP_T - CROP_B)

# Type scales against SAVE_W, not COMP_W: the composite is drawn on a 470 mm
# canvas but cropped to about 371 mm before saving, and it is the saved width
# that shrinks to the print width. Sizing off COMP_W overshoots by ~27% and
# pushes panel B's title past the right edge.
f06_txt <- composite_text_sizes(SAVE_W, 178)
TAG_SZ <- f06_txt$tag
TITLE_SZ <- f06_txt$title
SUBTITLE_SZ <- f06_txt$subtitle

composite_final <- ggdraw(xlim = c(CROP_L, CROP_R), ylim = c(CROP_B, CROP_T)) +
  theme(plot.background = element_rect(fill = "white", color = NA)) +
  # Panel B first (behind) so its white left margin is hidden under Panel A
  draw_grob(pB_grob,
    x = B_X, y = GRID_BOT, width = B_W, height = B_GRID_H,
    hjust = 0, vjust = 0
  ) +
  draw_grob(pB_leg_grob,
    x = layout_cfg$b_leg_x, y = layout_cfg$b_leg_y,
    width = layout_cfg$b_leg_w, height = layout_cfg$b_leg_h,
    hjust = 0.5, vjust = 0.5
  ) +
  draw_grob(pA_grob,
    x = layout_cfg$a_grob_x, y = layout_cfg$a_grob_y,
    width = layout_cfg$a_grob_w, height = layout_cfg$a_grob_h,
    hjust = 0, vjust = 0
  ) +
  draw_label("A",
    x = mm2x(A_LET_X), y = mm2y(A_LET_Y),
    size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1
  ) +
  draw_label("WGCNA Module-Trait Associations",
    x = mm2x(A_TTL_X), y = mm2y(A_TTL_Y),
    size = TITLE_SZ, fontface = "bold",
    hjust = 0, vjust = 1
  ) +
  draw_label(
    paste0(
      "9 modules, labelled enrichment then hub | ",
      "LMM (BH per-trait) | * BH, dashed nominal, bold |r| >= 0.40"
    ),
    x = mm2x(A_SUB_X), y = mm2y(A_SUB_Y),
    size = SUBTITLE_SZ, fontface = "bold.italic", colour = "grey40",
    hjust = 0, vjust = 1
  ) +
  draw_label("B",
    x = mm2x(B_LET_X), y = mm2y(B_LET_Y),
    size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1
  ) +
  draw_label("Module-Level NES Scatters",
    x = mm2x(B_TTL_X), y = mm2y(B_TTL_Y),
    size = TITLE_SZ, fontface = "bold",
    hjust = 0, vjust = 1
  ) +
  draw_label("fGSEA on module-member t-stats",
    x = mm2x(B_SUB_X), y = mm2y(B_SUB_Y),
    size = SUBTITLE_SZ, fontface = "bold.italic", colour = "grey40",
    hjust = 0, vjust = 1
  )

pdf_device <- get_raster_pdf_device()

ggsave(file.path(RPT_PDF, "F05.pdf"), composite_final,
  width = SAVE_W, height = SAVE_H, units = "mm",
  device = pdf_device
)
embed_pdf_fonts(file.path(RPT_PDF, "F05.pdf"))
ggsave(file.path(RPT_PNG, "F05.png"), composite_final,
  width = SAVE_W, height = SAVE_H, units = "mm",
  dpi = 300
)
message("F05 main composite done")

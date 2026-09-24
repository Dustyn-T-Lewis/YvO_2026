#!/usr/bin/env Rscript
# F06 Main — Phenotype Prediction (3-panel composite)
# A : per-module ROCs for what age classifies
# B : per-module ROCs for what training classifies, younger over older
# C : Module-Phenotype coupling, six hero screen cells, 3x2 scatter grid
#
# Sources ALL data-producing scripts in dependency order (SUPP panels first
# to generate CSVs needed by panel_A and the xlsx), then builds the main
# composite + xlsx + cleanup.

setwd(here::here())

pacman::p_load(patchwork, cowplot)

source("04_Figures/shared/style.R")
source("04_Figures/shared/figure_supplement_helpers.R")

BASE <- "04_Figures/F06"
DAT <- file.path(BASE, "c_data")

# Source all data-producing scripts in dependency order

message("F06 Composite: sourcing all panels")

# 1. SUPP module grid (writes c_data/module_grid/*.csv needed by panel_A)
source("04_Figures/F06/a_script/_supp_module_grid.R")

# 2. MAIN panel A (reads module_grid CSVs)
source("04_Figures/F06/a_script/_panel_A_auc_bars.R")

# 3. MAIN panel B (writes panel_B_full_screen_bh.csv)
source("04_Figures/F06/a_script/_panel_B_hero_grid.R")

# 4. SUPP panel B grid (reads panel_B_full_screen_bh.csv)
source("04_Figures/F06/a_script/_supp_panel_B_grid.R")

# 5. SUPP prepare ROC data (writes classifier_pilot_*.csv)
source("04_Figures/F06/a_script/_supp_prepare_roc.R")

# 6. SUPP multivariate classifier (writes panel_A_*.csv)
source("04_Figures/F06/a_script/_supp_multivariate.R")

# 7. SUPP LOSO sensitivity (writes loso_auc/loso_auc_summary.csv)
source("04_Figures/F06/a_script/_supp_loso_sensitivity.R")

# 8. SUPP LOSO WGCNA refit (~13-15 min; writes loso_auc/loso_wgcna_refit_*.csv)
source("04_Figures/F06/a_script/_supp_loso_wgcna_refit.R")

# Build main composite (pA/pB/pC via patchwork + cowplot tags)

# Restore paths clobbered by sourced helper scripts
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")
RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf")
DAT <- file.path(BASE, "c_data")
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)

# Canvas, in millimetres throughout. The previous version derived every title
# position from a chain of fractions and fine adjustments; with a third panel
# to place, stating each edge in mm and converting once is easier to check.
COMP_W <- 210
AB_H <- 63 # the row holding panels A and B, side by side
KEY_H <- 14 # the module colour key, shared by A and B, two rows
C_H <- 88 # the scatter grid, full width
TAG_AB <- 11 # spacer above the A/B row: one subtitle line each
TAG_C <- 13 # spacer above panel C
COMP_H <- TAG_AB + AB_H + KEY_H + TAG_C + C_H

composite <- (
  plot_spacer() /
    ((pA | pB) + plot_layout(widths = c(1, 1))) /
    pKey /
    plot_spacer() /
    pC
) +
  plot_layout(heights = c(TAG_AB, AB_H, KEY_H, TAG_C, C_H))

txt <- composite_text_sizes(COMP_W, 178)

# The figure embeds at 165.1 x 147 mm, so the composite is read at about 79% of
# the size it is drawn. Titles are scaled up to survive that reduction.
TITLE_SCALE <- 1.42
TAG_SZ <- round(txt$tag * TITLE_SCALE, 1)
TTL_SZ <- round(txt$title * TITLE_SCALE, 1)
SUB_SZ <- round(txt$subtitle * TITLE_SCALE, 1)

# y is measured from the bottom, so a panel's title sits at the top of the
# spacer above it. Subtitles drop one title-line below.
mm_y <- function(mm) mm / COMP_H
mm_x <- function(mm) mm / COMP_W
SUB_DROP <- mm_y((TTL_SZ / 72) * 25.4 * 1.35)

TAG_X <- mm_x(4)
X_TTL <- mm_x(8.4)
TAG_Y_AB <- 1 - mm_y(1.5)
TAG_Y_C <- mm_y(C_H + TAG_C) - mm_y(1.5)
# Panel B's title starts at the midline, where its half of the row begins.
TAG_X_B <- 0.5 + mm_x(1)

composite <- composite & theme(plot.margin = margin(0, 4, 0, 5))

panel_head <- function(g, tag, title, subtitle, x, y) {
  g +
    draw_label(tag,
      x = x, y = y + mm_y(0.4), size = TAG_SZ, fontface = "bold",
      hjust = 0, vjust = 1
    ) +
    draw_label(title,
      x = x + X_TTL, y = y,
      size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1
    ) +
    # draw_label defaults to lineheight 0.9, which leaves two subtitle lines
    # overlapping in the descender band. 1.0 clears them.
    draw_label(subtitle,
      x = x + X_TTL, y = y - SUB_DROP,
      size = SUB_SZ, fontface = "bold.italic", colour = "grey40",
      hjust = 0, vjust = 1, lineheight = 1.0
    )
}

composite_final <- ggdraw(composite) |>
  panel_head("A", pA_title, pA_subtitle, TAG_X, TAG_Y_AB) |>
  panel_head("B", pB_title, pB_subtitle, TAG_X_B, TAG_Y_AB) |>
  panel_head("C", pC_title, pC_subtitle, TAG_X, TAG_Y_C)

graphics.off()
pdf_device <- get_pdf_device()

ggsave(file.path(RPT_PDF, "MAIN_F06_composite.pdf"), composite_final,
  width = COMP_W, height = COMP_H, units = "mm",
  device = pdf_device, limitsize = FALSE
)
ggsave(file.path(RPT_PNG, "MAIN_F06_composite.png"), composite_final,
  width = COMP_W, height = COMP_H, units = "mm",
  dpi = 300, limitsize = FALSE
)
message("F06 composite saved: [A: age ROC | B: training ROC] / [C: hero grid]")

# Build supplementary xlsx (12 sheets)

f07 <- function(p) file.path(DAT, p)

message("F06 supplementary workbook")
f07_specs <- list(
  list(name = "module_grid_summary", path = f07("module_grid/module_grid_summary.csv")),
  list(name = "module_grid_curves", path = f07("module_grid/module_grid_curves.csv")),
  list(name = "panel_A_classifier_auc", path = f07("panel_A_multi_classifier_auc.csv")),
  list(name = "panel_A_feature_stability", path = f07("panel_A_feature_stability.csv")),
  list(name = "panel_A_permutation", path = f07("panel_A_permutation.csv")),
  list(name = "panel_A_roc_curves", path = f07("panel_A_roc_curves.csv")),
  list(name = "classifier_pilot_summary", path = f07("classifier_pilot_summary.csv")),
  list(name = "classifier_pilot_curves", path = f07("classifier_pilot_curves.csv")),
  list(name = "panel_B_full_screen", path = f07("panel_B_full_screen_bh.csv")),
  list(name = "loso_auc_summary", path = f07("loso_auc/loso_auc_summary.csv")),
  list(name = "loso_wgcna_refit_summary", path = f07("loso_auc/loso_wgcna_refit_summary.csv")),
  list(name = "loso_wgcna_refit_mod_stability", path = f07("loso_auc/loso_wgcna_refit_module_stability.csv"))
)
build_workbook(
  f07("F06_supplementary.xlsx"),
  title = "S7 Table \u2014 module discrimination and phenotype coupling",
  description = "Phenotype-prediction outputs: univariate module-outcome ROCs, multivariate classifiers, LOSO cross-validation, age-stratified module-phenotype coupling, and per-module Jaccard stability.",
  overview_df = data.frame(
    Sheet = c(
      "module_grid_summary",
      "module_grid_curves",
      "panel_A_classifier_auc",
      "panel_A_feature_stability",
      "panel_A_permutation",
      "panel_A_roc_curves",
      "classifier_pilot_summary",
      "classifier_pilot_curves",
      "panel_B_full_screen",
      "loso_auc_summary",
      "loso_wgcna_refit_summary",
      "loso_wgcna_refit_mod_stability"
    ),
    Description = c(
      "Per-module univariate ROC summary: AUC, permutation p, BH q for each module-outcome pair",
      "Per-module ROC curves: TPR/FPR coordinates for plotting",
      "Multivariate classifier AUC: raw-protein (k=10), phenotype-only, ME-stack, and delta-ME classifiers",
      "Feature stability: top-10 raw-protein classifier features across LOSO folds",
      "Permutation null: observed vs permuted AUC distributions for each classifier",
      "Multivariate ROC curves: TPR/FPR coordinates for each classifier",
      "Classifier pilot summary: module eigengene ROCs for age discrimination (top modules)",
      "Classifier pilot curves: TPR/FPR coordinates for classifier pilot ROCs",
      "Panel B full screen: 180-test module-phenotype correlations (9 modules \u00d7 2 sources \u00d7 5 outcomes \u00d7 2 strata) with BH correction",
      "LOSO fixed-module: leave-one-subject-out AUCs using fixed full-sample module definitions",
      "LOSO + WGCNA refit: full-refit LOSO \u2014 top-12 pair AUCs with network refit on n\u22121 subjects per fold, training modules matched to full-sample by Jaccard",
      "LOSO module stability: per-full-sample-module mean/min Jaccard of training-fold vs full-sample assignments, + count of folds where best Jaccard fell below 0.5"
    ),
    stringsAsFactors = FALSE
  ),
  sheet_specs = f07_specs
)
cleanup_after_workbook(f07_specs,
  extra_subdirs = c(f07("loso_auc"), f07("module_grid"))
)

#!/usr/bin/env Rscript
# S1b Figure: missingness classes, imputation and DEP counts, panels H-N.

setwd(here::here())

pacman::p_load(ggplot2, scales, patchwork)

source("04_Figures/shared/style.R")

PANELS <- "04_Figures/F00/a_script/supp/panels"
pH <- source_panel(file.path(PANELS, "S1b_H_miss_class_scatter.R"))
pI <- source_panel(file.path(PANELS, "S1b_I_miss_class_bar.R"))
pJ <- source_panel(file.path(PANELS, "S1b_J_benchmark.R"))
pK <- source_panel(file.path(PANELS, "S1b_K_imputation_density.R"))
pL <- source_panel(file.path(PANELS, "S1b_L_mnar_shift.R"))
pM <- source_panel(file.path(PANELS, "S1b_M_sample_integrity.R"))
pN <- source_panel(file.path(PANELS, "S1b_N_dep_heatmap.R"))

int_imp <- readRDS("02_imputation/c_data/00_report_intermediates.rds")

RPT <- "04_Figures/F00/b_reports/supp"

COMP_W <- 178
COMP_H <- 245
txt <- composite_text_sizes(COMP_W)

# Composite strip: titles, tags, subtitles, and a key that sits inside its own
# plot box rather than claiming a strip beside or beneath it.
#
# Each panel already carries a subtitle with its own counts; these used to be
# discarded at composite time, which threw away the per-panel numbers a reader
# needs and left the titles floating over unexplained plots.
#
# `at` places the key inside the panel. Every panel picks the corner its data
# leaves empty, so nothing is occluded and no panel loses width to a legend.
sl <- function(p, at = NULL) {
  p +
    theme(
      plot.title = element_text(
        face = "bold", size = FIG_TITLE_SIZE,
        margin = margin(b = 1)
      ),
      plot.subtitle = element_text(
        face = "bold.italic", size = FIG_SUBTITLE_SIZE,
        color = "grey30", margin = margin(t = 0, b = 2)
      ),
      legend.position = if (is.null(at)) "none" else "inside",
      legend.position.inside = at,
      # Anchoring justification to the same point puts that corner of the key
      # on that corner of the panel; without it the key centres on the point
      # and half of it hangs outside.
      legend.justification.inside = at,
      legend.direction = "vertical",
      legend.background = element_rect(fill = alpha("white", 0.75), colour = NA),
      legend.key = element_rect(fill = NA, colour = NA),
      legend.key.size = unit(2, "mm"),
      legend.text = element_text(size = 4.5),
      legend.title = element_blank(),
      legend.margin = margin(1, 2, 1, 2)
    )
}

# I's key repeats its own x-axis labels and N's colourbar restates numbers
# already printed in every tile, so both drop theirs rather than inset it.
page2 <- (sl(pH, at = c(0.99, 0.99)) | sl(pI)) /
  (sl(pJ, at = c(0.99, 0.01)) | sl(pK, at = c(0.99, 0.99))) /
  (sl(pL) | sl(pM, at = c(0.01, 0.99))) /
  sl(pN) +
  plot_layout(heights = c(1, 1.2, 1, 0.75)) +
  plot_annotation(
    title = "Pipeline QC \u2014 Imputation & Differential Expression",
    subtitle = sprintf(
      "missForest (OOB = %.3f) | %d MAR + %d MNAR proteins | limma + duplicateCorrelation",
      int_imp$oob_error, int_imp$n_mar_prots, int_imp$n_mnar_prots
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = txt$title),
      plot.subtitle = element_text(
        face = "italic", size = txt$subtitle - 1,
        color = "grey30"
      )
    )
  ) &
  theme(plot.tag = element_text(face = "bold", size = txt$tag))

pdf_dev <- get_pdf_device()

pdf_dev(file.path(RPT, "S1b.pdf"),
  width = COMP_W / 25.4, height = COMP_H / 25.4
)
print(page2)
dev.off()

ggsave(file.path(RPT, "S1b.png"), page2,
  width = COMP_W, height = COMP_H, units = "mm", dpi = 300
)

caption_supp(page2, "S1b", COMP_W, COMP_H, RPT)

message("S1b done")

#!/usr/bin/env Rscript
# S1a Figure: filtering, normalization, outliers and missingness, panels A-G.

setwd(here::here())

pacman::p_load(ggplot2, scales, patchwork)

source("04_Figures/shared/style.R")

PANELS <- "04_Figures/F00/a_script/supp/panels"
pA <- source_panel(file.path(PANELS, "S1a_A_filter_cascade.R"))
pB <- source_panel(file.path(PANELS, "S1a_B_protein_missingness.R"))
pC <- source_panel(file.path(PANELS, "S1a_C_pca_pre.R"))
pD <- source_panel(file.path(PANELS, "S1a_D_pca_post.R"))
pE <- source_panel(file.path(PANELS, "S1a_E_eta_squared.R"))
pF <- source_panel(file.path(PANELS, "S1a_F_outlier_consensus.R"))
pG <- source_panel(file.path(PANELS, "S1a_G_sample_missingness.R"))

int_norm <- readRDS("01_normalization/c_data/00_report_intermediates.rds")

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
# PCA panels: drop coord_fixed() for uniform alignment in composite
pC_comp <- pC + coord_cartesian()
pD_comp <- pD + coord_cartesian()

page1 <- (sl(pA) | sl(pB)) /
  (sl(pC_comp, at = c(0.01, 0.01)) | sl(pD_comp, at = c(0.01, 0.01))) /
  (sl(pE) | sl(pF, at = c(0.01, 0.99))) /
  sl(pG, at = c(0.99, 0.99)) +
  plot_layout(heights = c(1, 1, 1, 0.75)) +
  plot_annotation(
    title = "Pipeline QC \u2014 Pre-Processing",
    subtitle = sprintf(
      "DIA-NN \u2192 HPA filter \u2192 cyclic loess \u2192 outlier removal | %s \u2192 %s proteins \u00d7 %d samples",
      comma(int_norm$n_raw), comma(int_norm$dal_nrow), int_norm$dal_ncol
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

pdf_dev(file.path(RPT, "S1a.pdf"),
  width = COMP_W / 25.4, height = COMP_H / 25.4
)
print(page1)
dev.off()

ggsave(file.path(RPT, "S1a.png"), page1,
  width = COMP_W, height = COMP_H, units = "mm", dpi = 300
)

caption_supp(page1, "S1a", COMP_W, COMP_H, RPT)

message("S1a done")

#!/usr/bin/env Rscript
# S3 Figure: CV scatter triptych (A) full width, CV violins (B) and
# intra-individual variability (C) side by side.

setwd(here::here())

pacman::p_load(ggplot2, patchwork, cowplot)

source("04_Figures/shared/style.R")

PANELS <- "04_Figures/F02/a_script/panels"
pA <- source_panel(file.path(PANELS, "S3_A_cv_scatter.R"))
pB <- source_panel(file.path(PANELS, "S3_B_cv_violin.R"))
pC <- source_panel(file.path(PANELS, "S3_C_imputed.R"))

pSA_title <- pA[[1]]$labels$title
pSB_title <- pB$labels$title
pSC_title <- pC$labels$title
pA12 <- strip_for_composite(pA[[1]])
pA3 <- strip_for_composite(pA[[2]])
pB <- strip_for_composite(pB)
pC <- strip_for_composite(pC)

RPT <- "04_Figures/F02/b_reports"
pdf_dev <- get_pdf_device()

COMP_W <- 178; COMP_H <- 115
txt <- composite_text_sizes(COMP_W)

# pA12/pA3 from _supp_A (CV scatter), pB from _supp_B (CV violin), pC from _supp_C (imputed)
pSA_full <- (pA12 | pA3) + plot_layout(widths = c(2, 1))
composite <- (wrap_elements(pSA_full) / (pB | pC)) +
  plot_layout(heights = c(1, 0.7))

Y_TOP <- 0.985; Y_BOT <- 0.500
composite <- ggdraw(composite) +
  draw_label("A", x = 0.01, y = Y_TOP, size = txt$tag, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pSA_title, x = 0.04, y = Y_TOP, size = txt$title, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("B", x = 0.01, y = Y_BOT, size = txt$tag, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pSB_title, x = 0.04, y = Y_BOT, size = txt$title, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("C", x = 0.52, y = Y_BOT, size = txt$tag, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pSC_title, x = 0.55, y = Y_BOT, size = txt$title, fontface = "bold", hjust = 0, vjust = 1)

ggsave(file.path(RPT, "S3.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_dev)
ggsave(file.path(RPT, "S3.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

message("S3 done")

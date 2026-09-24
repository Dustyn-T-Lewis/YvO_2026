#!/usr/bin/env Rscript
# Figure 4: training plasticity across age (6-panel composite)
#
# A  protein quadrant scatter + flanking quadrant ORA, Tr.(Y) x Tr.(O)
# B  RRHO2 threshold-free rank overlap, same pair
# C  response magnitude and coherence per subject
# D  pathway NES scatter, Tr.(Y) x Tr.(O)
# E  pathway NES scatter, Aging x Tr.(O)
# F  this study's own FDR and Pi sets against the Training(Old) ranking
#
# A, B, C and D use contrast pairs sharing no group mean. E pairs Aging with
# Training(Old), which share Old_Pre; pathway aggregation shrinks that coupling
# because the shared per-protein noise averages down while genuine
# between-pathway differences do not, so E is reported at pathway level only.
# F's Aging rows are coupled at protein level and take a Pre/Post swap null.
#
# F05 was folded into this figure on 2026-08-24; see docs/decisions.md.
# limma tree: 2106 proteins over 62 samples, and Training(Old) has no FDR hits,
# so panel F's Tr.(O) rows are empty by construction rather than by filtering.

setwd(here::here())

pacman::p_load(
  dplyr, tidyr, tibble, stringr, readr, ggplot2, patchwork, cowplot, readxl
)

# ggrepel places labels by a stochastic search, so an unseeded render puts
# them somewhere new each time. run_all.R runs each script in its own
# Rscript child, which starts from a time-seeded RNG.
set.seed(42)

source("04_Figures/shared/style.R")

RPT <- "04_Figures/F04/b_reports/main"
PANELS <- "04_Figures/F04/a_script/main/panels"

pdf_device <- get_pdf_device()

message("F04 composite: sourcing panels")

panel_a <- source_panel(file.path(PANELS, "A_quadrant_ora.R"))
a_stats <- attr(panel_a, "stats")
n_total_A <- a_stats$n_total
r_spear <- a_stats$r_spear
conc_fdr <- a_stats$conc_fdr
conc_pi <- a_stats$conc_pi
panel_a <- panel_a +
  plot_annotation(
    title = NULL, subtitle = NULL,
    theme = theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank()
    )
  ) &
  theme(legend.position = "none")

panel_b <- source_panel(file.path(PANELS, "B_rrho2.R"))
sweep_txt <- attr(panel_b, "sweep_txt")

panel_c <- source_panel(file.path(PANELS, "C_trajectory.R"))

pD <- source_panel(file.path(PANELS, "D_nes_concordance.R"))
nes_legend <- cowplot::get_plot_component(pD, "guide-box-bottom", return_all = FALSE)
panel_d <- pD + theme(legend.position = "none")
rho_D <- attr(pD, "stats")$rho
frac_D <- attr(pD, "stats")$frac

pE <- source_panel(file.path(PANELS, "E_nes_reversal.R"))
panel_e <- pE + theme(legend.position = "none")
rho_E <- attr(pE, "stats")$rho
frac_E <- attr(pE, "stats")$frac

panel_f <- source_panel(file.path(PANELS, "F_fry_barcode.R"))

COMP_W <- 460
COMP_H <- 320
PRINT_SCALE <- COMP_W / 178
# The panels leave print_scale_apply.R's mutated FIG_* globals behind them, and
# composite_text_sizes multiplies by the canvas factor itself, so the constants
# have to be back at their style.R values or the titles scale twice. The panels
# are finished ggplot objects by now and keep the sizes they were built with.
source("04_Figures/shared/style.R")
# Authored at 460 mm and shrunk to the 178 mm print width, so the type scales
# up by the same factor to land at the shared 7 pt on the page.
f04_txt <- composite_text_sizes(COMP_W, 178)
TAG_SZ <- f04_txt$tag
TTL_SZ <- f04_txt$title
SUB_SZ <- f04_txt$subtitle

ttl <- c(
  A = "Training Concordance",
  B = "Threshold-Free Overlap",
  # C sits in a narrow 4-column box. At the 6 pt subtitle the title clips past
  # roughly 18 characters and the subtitle past roughly 25.
  C = "Blunted Response",
  D = "Pathway Concordance",
  E = "Pathway Reversal",
  F = "Signature Rank"
)
sub <- c(
  A = sprintf(
    # draw_label does not wrap and the eight-column box clips past about 80
    # characters. The fitted slopes live on the lines themselves instead.
    "N = %d, ρ %.2f | FDR %d: %.0f%% dir, ρ %.2f | Π %d: %.0f%%, ρ %.2f",
    n_total_A, r_spear,
    conc_fdr$n, 100 * conc_fdr$agree, conc_fdr$rho,
    conc_pi$n, 100 * conc_pi$agree, conc_pi$rho
  ),
  B = sweep_txt,
  C = "% young magnitude lost",
  D = sprintf("ρ = %.2f | %.0f%% concordant | no shared group mean", rho_D, frac_D * 100),
  E = sprintf("ρ = %.2f | %.0f%% reversed | pathway level only", rho_E, frac_E * 100),
  # draw_label does not wrap, and the six-column box runs out at about 54
  # characters, so the age wording is the short form here.
  F = "fry rotation test: younger sets in the older ranking"
)

layout <- paste(
  "##################",
  strrep("AAAAAAAABBBBBBCCCC\n", 6),
  "##################",
  strrep("DDDDDDEEEEEEFFFFFF\n", 6),
  "GGGGGGGGGGGGGGGGGG",
  sep = "\n"
)
layout <- paste(Filter(nzchar, strsplit(layout, "\n")[[1]]), collapse = "\n")

fig <- wrap_elements(full = panel_a) +
  # B's tile field sat 3.6 pt above panel A's border. The nudge is applied
  # here rather than in panel_E_rrho2.R because that engine also renders the
  # supplementary aging RRHO2, which must not move. B fills its cell height,
  # so a top margin lowers the top edge and shortens the square by the same
  # amount rather than translating it.
  wrap_elements(full = panel_b + theme(plot.margin = margin(1.27, 0, 0, 0, "mm"))) +
  wrap_elements(full = panel_c) +
  wrap_elements(full = panel_d) +
  wrap_elements(full = panel_e) +
  wrap_elements(full = panel_f) +
  wrap_elements(full = nes_legend) +
  plot_layout(
    # Panel C is given 18 pt more, taken from B's columns rather than A's. B is
    # a centred coord_fixed square in a cell 88 pt wider than the square needs,
    # so the width spent here was rendering as nothing and B stays
    # height-limited at the same size. Taking it from A instead would reflow
    # its ggrepel labels and push its bar labels further past their bars.
    design = layout, widths = c(rep(1, 8), rep(0.958, 6), rep(1.063, 4)),
    # Bands 1 and 8 carry the title and subtitle for the row below them. They
    # are fractions of the total, so trimming the middle band shrinks the top
    # one too, and the top band had to grow from 7 to hold its ground. The
    # middle band was 12, which left 26 pt of dead canvas between the D/E/F
    # subtitles and the plots they name; at 9 that gap is about 5 pt.
    heights = c(8, rep(13, 6), 9, rep(13, 6), 6)
  )

# Hand-tuned to the column boundaries, with no link back to the layout, so a
# change to widths has to be carried here by hand. C's cell left edge moved
# 18 pt left with the widths change above; B's did not move.
X_TAG <- c(A = 0.004, B = 0.448, C = 0.768, D = 0.004, E = 0.334, F = 0.659)
# The top row's tags started 0.6 pt from the page edge, against the bottom
# row's 2.1 pt below its band. Lowered to match that relationship and to put
# a readable gap between the subtitles and the panel tops at 52.4.
Y_TAG <- c(A = 0.987, B = 0.987, C = 0.987, D = 0.518, E = 0.518, F = 0.518)
X_TTL <- 0.030
SUB_OFFSET <- 0.019

composite_final <- Reduce(
  function(p, k) {
    p +
      draw_label(k, x = X_TAG[[k]], y = Y_TAG[[k]], size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
      draw_label(ttl[[k]],
        x = X_TAG[[k]] + X_TTL, y = Y_TAG[[k]], size = TTL_SZ,
        fontface = "bold", hjust = 0, vjust = 1
      ) +
      draw_label(sub[[k]],
        x = X_TAG[[k]] + X_TTL, y = Y_TAG[[k]] - SUB_OFFSET, size = SUB_SZ,
        fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40"
      )
  },
  names(ttl),
  init = ggdraw(fig)
)

ggsave(file.path(RPT, "F04.pdf"), composite_final,
  width = COMP_W, height = COMP_H, units = "mm", device = pdf_device
)
ggsave(file.path(RPT, "F04.png"), composite_final,
  width = COMP_W, height = COMP_H, units = "mm", dpi = 300
)

message("F04 composite (6-panel, A/B/C over D/E/F) saved")

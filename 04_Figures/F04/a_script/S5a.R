#!/usr/bin/env Rscript
# S5a Figure: 4-panel diagnostic composite (2x2 grid)
# A: Spearman rho bootstrap        B: GO Slim distribution
# C: Concordance by effect size    D: Shared-baseline coupling diagnostic

setwd(here::here())

pacman::p_load(dplyr, tidyr, tibble, stringr, readr, ggplot2, patchwork, cowplot)

source("04_Figures/shared/style.R")

pdf_device <- get_pdf_device()

PANELS <- "04_Figures/F04/a_script/panels"
RPT <- "04_Figures/F04/b_reports"

message("S5a composite: sourcing panels")
pS_rho_boot <- strip_for_composite(source_panel(file.path(PANELS, "S5a_A_rho_bootstrap.R")))
pS_goslim <- strip_for_composite(source_panel(file.path(PANELS, "S5a_B_goslim_bars.R")))
pS_conc_mag <- strip_for_composite(source_panel(file.path(PANELS, "S5a_C_concordance_magnitude.R")))
pS_coupling <- strip_for_composite(source_panel(file.path(PANELS, "S5a_D_coupling_null.R")))

COMP_W <- 280
COMP_H <- 210

cts <- composite_text_sizes(COMP_W, 178)
TAG_SZ <- cts$tag
TTL_SZ <- cts$title + 1
SUB_SZ <- cts$subtitle + 0.5

axis_fix <- theme(
  legend.position = "none",
  axis.title.y = element_text(margin = margin(0, 2, 0, 0)),
  axis.title.x = element_text(margin = margin(2, 0, 0, 0))
)
# Measured off the render rather than guessed: at 18 pt the top row's panel
# border sat at 8.3 mm and the subtitle ran to 8.7, so its descenders crossed
# the frame. 22 pt puts the border at 9.7. The bottom row was further out --
# border 106.2 mm against a subtitle reaching 108.0 -- so it takes 36.
pS_rho_boot <- pS_rho_boot + axis_fix +
  theme(plot.margin = margin(22, 10, 2, 5))
pS_goslim <- pS_goslim + axis_fix +
  theme(plot.margin = margin(22, 10, 2, 5))
# The bottom row's titles are drawn over its top margin, so it needs more of
# one than the top row, whose titles sit above the canvas edge.
pS_conc_mag <- pS_conc_mag + axis_fix +
  theme(plot.margin = margin(36, 14, 2, 5))
pS_coupling <- pS_coupling + axis_fix +
  theme(plot.margin = margin(36, 10, 2, 5), legend.position = "bottom")

fig <- (pS_rho_boot | pS_goslim) / (pS_conc_mag | pS_coupling) +
  plot_layout(heights = c(1, 1.05), guides = "keep")

X_A <- 0.020
X_B <- 0.510
X_C <- 0.020
X_D <- 0.510
X_TTL <- 0.022
SUB_OFFSET <- 0.017
Y_TOP <- 0.991
# The bottom row's header is drawn into the panels' 30 pt top margin. At
# 0.505 the subtitle's descenders landed on the panel border; 0.518 centres
# the band in the margin instead.
Y_BOT <- 0.518

composite_final <- ggdraw(fig) +
  draw_label("A", x = X_A, y = Y_TOP, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("Spearman \u03c1 Bootstrap", x = X_A + X_TTL, y = Y_TOP, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("1000 replicates, 95% CI", x = X_A + X_TTL, y = Y_TOP - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("B", x = X_B, y = Y_TOP, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("GO Slim Distribution", x = X_B + X_TTL, y = Y_TOP, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("Proteins per GO Slim category by quadrant", x = X_B + X_TTL, y = Y_TOP - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("C", x = X_C, y = Y_BOT, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("Concordance by Effect Magnitude", x = X_C + X_TTL, y = Y_BOT, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("Sign agreement across quintiles of min(|log2FC|)", x = X_C + X_TTL, y = Y_BOT - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("D", x = X_D, y = Y_BOT, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("Shared-Baseline Coupling", x = X_D + X_TTL, y = Y_BOT, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("Observed r against a null that preserves the coupling", x = X_D + X_TTL, y = Y_BOT - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40")

legend_df <- data.frame(
  x = 1:3,
  fill_lab = factor(c("Concordant Up", "Concordant Down", "Discordant"),
    levels = c("Concordant Up", "Concordant Down", "Discordant")
  )
)
legend_plot <- ggplot(legend_df, aes(x = x, y = 1, fill = fill_lab)) +
  geom_col() +
  scale_fill_manual(
    values = c(
      "Concordant Up" = "#E57373",
      "Concordant Down" = "#64B5F6",
      "Discordant" = "#FFB74D"
    ),
    name = "Quadrant"
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold", size = 8),
    legend.text = element_text(size = 7),
    legend.key.size = unit(4, "mm"),
    legend.spacing.x = unit(3, "mm")
  )

legend_grob <- cowplot::get_legend(legend_plot)
composite_final <- composite_final +
  draw_plot(legend_grob, x = 0.08, y = -0.006, width = 0.45, height = 0.060)

ggsave(file.path(RPT, "S5a.pdf"), composite_final,
  width = COMP_W, height = COMP_H, units = "mm", device = pdf_device
)
ggsave(file.path(RPT, "S5a.png"), composite_final,
  width = COMP_W, height = COMP_H, units = "mm", dpi = 300
)

message(sprintf("S5a composite (4-panel) saved: %s x %s mm", COMP_W, COMP_H))

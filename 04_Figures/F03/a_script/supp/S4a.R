#!/usr/bin/env Rscript
# S4a Figure: p, Pi and FDR distributions, MA plots, and DEP retention after
# outlier removal.

setwd(here::here())

pacman::p_load(ggplot2, patchwork, cowplot)

source("04_Figures/shared/style.R")

PANELS <- "04_Figures/F03/a_script/supp/panels"
panels <- lapply(
  c("S4a_A_p_value.R", "S4a_B_pi_score.R", "S4a_C_fdr.R", "S4a_D_ma.R",
    "S4a_E_outlier_retention.R"),
  \(f) source_panel(file.path(PANELS, f))
)
all_titles <- lapply(panels, \(p) p$labels$title)
panels <- lapply(panels, strip_for_composite)

RPT <- "04_Figures/F03/b_reports/supp"

COMP_W <- 178; COMP_H <- 225
txt <- composite_text_sizes(COMP_W)
# 15 pt of head, not 9: the titles below are drawn over the panels, and each
# panel writes its own contrast name at the top of its plotting area. At 9 pt
# the two shared a line and "Raw p-value distribution" printed through "Aging".
grid <- (panels[[1]] | panels[[2]]) / (panels[[3]] | panels[[4]]) / panels[[5]] &
  theme(plot.margin = margin(15, 4, 4, 4))

X_L <- 0.012; X_R <- 0.512; X_TTL <- 0.029
Y_R1 <- 0.974; Y_R2 <- 0.651; Y_R3 <- 0.321

all_tags <- LETTERS[1:5]
all_xs <- c(X_L, X_R, X_L, X_R, X_L)
all_ys <- c(Y_R1, Y_R1, Y_R2, Y_R2, Y_R3)

composite <- ggdraw(grid)
for (i in seq_along(all_tags)) {
  composite <- composite +
    draw_label(all_tags[i], x = all_xs[i], y = all_ys[i],
               fontface = "bold", size = txt$tag, hjust = 0, vjust = 1) +
    draw_label(all_titles[[i]], x = all_xs[i] + X_TTL, y = all_ys[i],
               fontface = "bold", size = txt$title, hjust = 0, vjust = 1)
}

ggsave(file.path(RPT, "S4a.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = get_pdf_device())
ggsave(file.path(RPT, "S4a.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

caption_supp(composite, "S4a", COMP_W, COMP_H, RPT)

message("S4a done")

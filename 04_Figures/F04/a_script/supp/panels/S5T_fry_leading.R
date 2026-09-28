#!/usr/bin/env Rscript
# S5 Table, sheet SUPP_fry_leading: fry leading-edge proteins.
#
# Diagnostic for main panel F, the barcode: which Pi-selected Training(Young)
# proteins sit furthest out in the Training(Old) ranking that panel F tests
# them against, computed straight from the DEP table.

setwd(here::here())

pacman::p_load(dplyr, tidyr, tibble, stringr, readr, ggplot2, patchwork, cowplot)

source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

DAT <- "04_Figures/F04/c_data/panel_supp"
RPT <- "04_Figures/F04/b_reports/supp/panels"
for (d in c(DAT, RPT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

driving_df <- readr::read_csv(DEP_RESULTS, show_col_types = FALSE) |>
  dplyr::filter(
    !is.na(pi_score_Training_Young), pi_score_Training_Young < 0.05,
    !is.na(t_Training_Old)
  ) |>
  dplyr::transmute(
    gene,
    t_test = t_Training_Old,
    logFC_TY = logFC_Training_Young,
    logFC_TO = logFC_Training_Old,
    direction = dplyr::if_else(logFC_Training_Young > 0, "Up", "Down")
  ) |>
  dplyr::arrange(dplyr::desc(abs(t_test))) |>
  as.data.frame()

top_df <- head(driving_df, 20)
top_df$gene <- factor(top_df$gene, levels = rev(top_df$gene))
top_df$dir_label <- ifelse(top_df$direction == "Up", "Concordant Up", "Concordant Down")

write_csv(top_df, file.path(DAT, "SUPP_fry_leading_edge.csv"))
message(sprintf("Top 20 fry driving proteins (of %d total)", nrow(driving_df)))

DIR_COLS <- c("Concordant Up" = "#D6604D", "Concordant Down" = "#4393C3")

pS_fry_lead <- ggplot(top_df, aes(x = t_test, y = gene, color = dir_label)) +
  geom_segment(aes(x = 0, xend = t_test, y = gene, yend = gene),
    linewidth = 0.6
  ) +
  geom_point(size = 2.5) +
  geom_vline(xintercept = 0, color = "grey50", linewidth = 0.3) +
  scale_color_manual(values = DIR_COLS, name = "Direction") +
  labs(
    title = "fry Leading-Edge Proteins",
    subtitle = "Top 20 driving proteins ranked by |t-stat| in Training Old",
    x = "t-statistic (Training Old)", y = NULL
  ) +
  FIG_THEME

PW <- 89
PH <- 70
ggsave(file.path(RPT, "S5T_fry_leading.png"), pS_fry_lead,
  width = PW, height = PH, units = "mm", dpi = 300
)
ggsave(file.path(RPT, "S5T_fry_leading.pdf"), pS_fry_lead,
  width = PW, height = PH, units = "mm", device = get_pdf_device()
)

message("SUPP Panel E (fry leading edge) done")

invisible(pS_fry_lead)

#!/usr/bin/env Rscript
# SUPP panel: fry leading-edge proteins.
#
# Diagnostic for main panel F, the barcode: which Pi-selected Training(Young)
# proteins sit furthest out in the Training(Old) ranking that panel F tests
# them against. It used to read a "panel_C_fry_driving" sheet written by
# panel_C_fry.R, which was retired with F05 on 2026-08-26, so it now computes
# the same quantity straight from the DEP table.
#
# Sourced by 02_supp_panels.R — expects style.R already loaded.
# Exports: pS_fry_lead (ggplot)

BASE <- "04_Figures/F04"
DAT <- file.path(BASE, "c_data", "panel_supp")
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

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

RPT_PNG <- file.path(BASE, "b_reports", "supp", "png", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "pdf", "panels")
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)

PW <- 89
PH <- 70
ggsave(file.path(RPT_PNG, "SUPP_fry_leading.png"), pS_fry_lead,
  width = PW, height = PH, units = "mm", dpi = 300
)
ggsave(file.path(RPT_PDF, "SUPP_fry_leading.pdf"), pS_fry_lead,
  width = PW, height = PH, units = "mm", device = pdf_device
)

pS_fry_lead <- strip_for_composite(pS_fry_lead)

message("SUPP Panel E (fry leading edge) done")

# S4a histogram for one statistic, shared by panels A-C: one facet per
# contrast, 20 bins, the count below 0.05 printed in each.

pacman::p_load(dplyr, readr, readxl, ggplot2)

source("04_Figures/shared/style.R")

CTRS <- c("Aging", "Training_Young", "Training_Old", "Interaction")

read_contrasts <- function() {
  per_contrast <- lapply(CTRS, \(ctr) {
    as.data.frame(read_excel("03_DEP/c_data/03_DEP_results.xlsx", sheet = ctr))
  })
  names(per_contrast) <- CTRS
  per_contrast
}

dist_panel <- function(col, fill, xlab, vline = NULL, stat_fmt, title, tag, name) {
  RPT <- "04_Figures/F03/b_reports/panels"
  DAT <- "04_Figures/F03/c_data/supp"
  for (d in c(RPT, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  per_contrast <- read_contrasts()

  hist_df <- bind_rows(lapply(CTRS, \(ctr) {
    tibble(contrast = ctr, value = per_contrast[[ctr]][[col]])
  })) |> filter(!is.na(value)) |>
    mutate(contrast = factor(contrast, levels = CTRS))

  n_sig <- hist_df |>
    summarise(n_sig = sum(value < 0.05), .by = contrast)

  n_bins <- 20
  p <- ggplot(hist_df, aes(value)) +
    geom_histogram(breaks = seq(0, 1, length.out = n_bins + 1),
                   fill = fill, color = "white", linewidth = 0.3)

  if (!is.null(vline))
    p <- p + geom_vline(xintercept = vline, linetype = "dashed",
                        color = "grey30", linewidth = 0.4)
  if (col == "P.Value") {
    uniform_ref <- hist_df |> summarise(h = n() / n_bins, .by = contrast)
    p <- p + geom_hline(data = uniform_ref, aes(yintercept = h),
                        linetype = "dashed", color = "grey40", linewidth = 0.4)
  }

  p <- p +
    geom_text(data = n_sig,
              aes(x = 0.5, y = Inf, label = CTR_SHORT[as.character(contrast)]),
              inherit.aes = FALSE, hjust = 0.5, vjust = 1.2,
              size = 2.0, fontface = "bold", color = "grey20") +
    geom_text(data = n_sig,
              aes(x = 0.5, y = Inf, label = sprintf(stat_fmt, n_sig)),
              inherit.aes = FALSE, hjust = 0.5, vjust = 3.8,
              size = 2.8, fontface = "bold", color = "grey40") +
    facet_wrap(~contrast, ncol = 2, scales = "free_y",
               labeller = labeller(contrast = CTR_SHORT)) +
    labs(title = title,
         subtitle = sprintf("%s proteins | 20 bins",
                            format(round(nrow(hist_df) / length(CTRS)), big.mark = ",")),
         x = xlab, y = "Proteins", tag = tag) +
    FIG_THEME + theme(strip.background = element_blank(), strip.text = element_blank())

  write_csv(hist_df, file.path(DAT, sprintf("panel_%s.csv", gsub("[. ]", "_", col))))
  ggsave(file.path(RPT, paste0(name, ".png")), p,
         width = 89, height = 75, units = "mm", dpi = 300)
  ggsave(file.path(RPT, paste0(name, ".pdf")), p,
         width = 89, height = 75, units = "mm", device = get_pdf_device())
  invisible(p)
}

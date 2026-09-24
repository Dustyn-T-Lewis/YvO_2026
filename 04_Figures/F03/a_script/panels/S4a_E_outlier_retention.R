#!/usr/bin/env Rscript
# S4a Figure E: DEPs kept after removing the two outlier samples, by contrast.

setwd(here::here())

pacman::p_load(dplyr, tidyr, readxl, ggplot2)

source("04_Figures/shared/style.R")

RPT <- "04_Figures/F03/b_reports/panels"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
CTRS <- c("Aging", "Training_Young", "Training_Old", "Interaction")

out_sens <- tryCatch(
  as.data.frame(read_excel("03_DEP/c_data/03_DEP_results.xlsx", sheet = "outlier_sensitivity")),
  error = \(e) NULL)
if (is.null(out_sens) || nrow(out_sens) == 0) {
  stop("outlier_sensitivity sheet missing from 03_DEP_results.xlsx — ",
       "run 03_DEP/a_script/02_generate_reports.R first")
}

long_df <- out_sens |>
  mutate(Contrast = factor(Contrast, levels = CTRS)) |>
  pivot_longer(c(FDR_full, FDR_reduced, Pi_full, Pi_reduced),
               names_to = "metric_cohort", values_to = "n") |>
  mutate(metric = factor(sub("_.*", "", metric_cohort),
                         levels = c("FDR", "Pi"), labels = c("FDR < 0.05", "Π < 0.05")),
         cohort = factor(sub(".*_", "", metric_cohort), levels = c("full", "reduced")))

p <- ggplot(long_df, aes(cohort, n, fill = cohort)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = n), vjust = -0.3, size = 2.8, fontface = "bold") +
  facet_grid(metric ~ Contrast, scales = "free_y", switch = "y",
             labeller = labeller(Contrast = CTR_SHORT)) +
  scale_fill_manual(values = c(full = "#2166AC", reduced = "#B2182B"),
                    labels = c("Full", "Outlier-removed"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(title = "DEP retention (outlier removal)", x = NULL, y = "DEPs", tag = "e") +
  FIG_THEME + theme(strip.text.x = element_text(face = "bold",
                                                size = FIG_STRIP_SIZE),
                    strip.text.y = element_text(face = "bold", size = FIG_STRIP_SIZE - 1),
                    legend.position = "top", legend.key.size = unit(3, "mm"))
ggsave(file.path(RPT, "S4a_E_outlier_retention.png"), p, width = 178, height = 75,
       units = "mm", dpi = 300)
ggsave(file.path(RPT, "S4a_E_outlier_retention.pdf"), p, width = 178, height = 75,
       units = "mm", device = get_pdf_device())

invisible(p)

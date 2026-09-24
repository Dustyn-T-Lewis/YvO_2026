#!/usr/bin/env Rscript
# S5a Figure D: shared-baseline coupling diagnostic.
#
# Aging is Old_Pre - Young_Pre and Training_Old is Old_Post - Old_Pre, so the
# two share Old_Pre with opposite signs and correlate negatively before any
# biology enters. The panel shows the observed correlation against a null that
# keeps that structure and removes only the training effect, by swapping Pre
# and Post labels within each older participant and re-selecting proteins. The
# observed value falling inside the null is the whole argument, and it is
# easier to see than to read off an interval.

setwd(here::here())

pacman::p_load(dplyr, tidyr, tibble, stringr, readr, ggplot2, patchwork, cowplot)

source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

pdf_device <- get_pdf_device()

BASE <- "04_Figures/F04"
DAT <- file.path(BASE, "c_data", "panel_supp")
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

draws <- read_csv("03_DEP/c_data/03_reversal_null_draws.csv",
  show_col_types = FALSE
)
est <- read_csv("03_DEP/c_data/03_reversal_aging_fdr.csv",
  show_col_types = FALSE
)

pick <- function(pattern) est$r[grepl(pattern, est$estimator)][1]
obs_r <- pick("^Group-mean axes")
null_p <- est$p_value[grepl("^Coupling-preserving", est$estimator)][1]
split_r <- pick("^Split-half")

null_draws <- draws$r[draws$estimator == "coupling_null"]
split_draws <- draws$r[draws$estimator == "split_half"]

plot_df <- tibble(
  r = c(null_draws, split_draws),
  source = factor(
    rep(
      c("Coupling-preserving null", "Split-half over disjoint subjects"),
      c(length(null_draws), length(split_draws))
    ),
    levels = c("Coupling-preserving null", "Split-half over disjoint subjects")
  )
)

write_csv(
  tibble(
    quantity = c("observed", "null_mean", "null_p", "split_half_mean"),
    value = c(obs_r, mean(null_draws), null_p, mean(split_draws))
  ),
  file.path(DAT, "SUPP_coupling_null.csv")
)

pS_coupling <- ggplot(plot_df, aes(r, fill = source)) +
  geom_histogram(bins = 40, colour = NA, alpha = 0.75, position = "identity") +
  geom_vline(
    xintercept = obs_r, colour = "#D6604D", linewidth = 0.8,
    linetype = "dashed"
  ) +
  annotate("text",
    x = obs_r, y = Inf, vjust = 3.4, hjust = -0.08,
    label = sprintf("observed r = %.3f", obs_r),
    size = BASE_STAT, fontface = "bold", colour = "#D6604D"
  ) +
  scale_fill_manual(values = c("grey45", "#4393C3"), name = NULL) +
  labs(
    title = "Shared-Baseline Coupling",
    subtitle = sprintf(
      "observed sits inside the coupling-preserving null (p = %.2f)", null_p
    ),
    x = "Pearson r, Aging against Training in Old", y = "Draws"
  ) +
  FIG_THEME +
  theme(legend.position = "bottom")

RPT_PNG <- file.path(BASE, "b_reports", "supp", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "panels")
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)

PW <- 89
PH <- 70
ggsave(file.path(RPT_PNG, "S5a_D_coupling_null.png"), pS_coupling,
  width = PW, height = PH, units = "mm", dpi = 300
)
ggsave(file.path(RPT_PDF, "S5a_D_coupling_null.pdf"), pS_coupling,
  width = PW, height = PH, units = "mm", device = pdf_device
)

message(sprintf(
  "SUPP Panel D (coupling null) done: observed %.3f, null p %.2f, split-half %.3f",
  obs_r, null_p, split_r
))

invisible(pS_coupling)

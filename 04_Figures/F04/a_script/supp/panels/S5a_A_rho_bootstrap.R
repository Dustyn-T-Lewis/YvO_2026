#!/usr/bin/env Rscript
# S5a Figure A: diagnostic for main panels A and D, showing Spearman rho between
# Training_Young and Training_Old logFC is robust to resampling (1000 replicates).

setwd(here::here())

pacman::p_load(dplyr, tidyr, tibble, stringr, readr, ggplot2, patchwork, cowplot)

source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

pacman::p_load(boot)

DAT <- "04_Figures/F04/c_data/panel_supp"
RPT <- "04_Figures/F04/b_reports/supp/panels"
for (d in c(DAT, RPT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)

boot_df <- dep_df |>
  transmute(logFC_TY = logFC_Training_Young,
            logFC_TO = logFC_Training_Old) |>
  filter(!is.na(logFC_TY), !is.na(logFC_TO))

n_prot <- nrow(boot_df)
obs_rho <- cor(boot_df$logFC_TY, boot_df$logFC_TO, method = "spearman")

rho_stat <- function(data, idx) {
  d <- data[idx, ]
  cor(d$logFC_TY, d$logFC_TO, method = "spearman")
}

set.seed(42)
boot_out <- boot(boot_df, statistic = rho_stat, R = 1000)

rho_vals <- boot_out$t[, 1]
ci <- quantile(rho_vals, c(0.025, 0.975))

export_df <- tibble(replicate = seq_along(rho_vals), rho = rho_vals)
write_csv(export_df, file.path(DAT, "SUPP_rho_bootstrap.csv"))

message(sprintf("Bootstrap rho: %.3f [%.3f, %.3f], n = %d, R = 1000",
                obs_rho, ci[1], ci[2], n_prot))

ci_df <- tibble(xmin = ci[1], xmax = ci[2])

pS_rho_boot <- ggplot(tibble(rho = rho_vals), aes(x = rho)) +
  geom_rect(data = ci_df, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = alpha("#4393C3", 0.25)) +
  geom_histogram(bins = 40, fill = alpha("#4393C3", 0.15), color = "white", linewidth = 0.3) +
  geom_vline(xintercept = obs_rho, color = "#D6604D", linewidth = 0.8,
             linetype = "dashed") +
  annotate("text", x = obs_rho, y = Inf, vjust = 1.5, hjust = -0.1,
           label = sprintf("rho = %.3f", obs_rho),
           size = BASE_STAT, fontface = "bold", color = "#D6604D") +
  labs(title = "Concordance Spearman rho Bootstrap",
       subtitle = sprintf("rho = %.2f [%.3f, %.3f], n = %d, R = 1000",
                          obs_rho, ci[1], ci[2], n_prot),
       x = "Spearman rho", y = "Count") +
  FIG_THEME

PW <- 89; PH <- 70
ggsave(file.path(RPT, "S5a_A_rho_bootstrap.png"), pS_rho_boot,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT, "S5a_A_rho_bootstrap.pdf"), pS_rho_boot,
       width = PW, height = PH, units = "mm", device = get_pdf_device())

message("SUPP Panel B (rho bootstrap) done")

invisible(pS_rho_boot)

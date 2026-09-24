#!/usr/bin/env Rscript
# S1a Figure E: eta-squared, biological signal retained.

setwd(here::here())
source("04_Figures/F00/a_script/supp/panels/_inputs.R", local = TRUE)

eta_df <- tibble(eta2 = as.numeric(int_norm$eta2_vals))
eta_med <- median(eta_df$eta2, na.rm = TRUE)
eta_90 <- quantile(eta_df$eta2, 0.90, na.rm = TRUE)

pE <- ggplot(eta_df, aes(eta2)) +
  geom_density(fill = "#4CAF50", alpha = 0.55, linewidth = 0.4) +
  geom_vline(xintercept = eta_med, color = "grey30", linewidth = 0.4) +
  geom_vline(
    xintercept = eta_90, linetype = "dashed",
    color = "#E05A4E", linewidth = 0.4
  ) +
  annotate("text",
    x = eta_med, y = Inf, label = sprintf("median = %.2f", eta_med),
    vjust = 1.5, hjust = -0.1, size = 2.2
  ) +
  annotate("text",
    x = eta_90, y = Inf, label = sprintf("90th = %.2f", eta_90),
    vjust = 3.0, hjust = -0.1, size = 2.2, color = "#E05A4E"
  ) +
  labs(
    x = expression(eta^2), y = "Density", tag = "E",
    title = expression("Biological signal retention (" * eta^2 * ")"),
    subtitle = sprintf("Group effect | %s proteins", comma(nrow(eta_df)))
  ) +
  FIG_THEME

save_panel(pE, "S1a_E_eta_squared", "panel_E", eta_df)

invisible(pE)

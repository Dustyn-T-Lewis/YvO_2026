#!/usr/bin/env Rscript
# S1b Figure L: shift of MNAR imputed values.

setwd(here::here())
source("04_Figures/F00/a_script/supp/panels/_inputs.R", local = TRUE)

audit <- int_imp$mnar_audit
shift_mean <- mean(audit$shift, na.rm = TRUE)

pL <- ggplot(audit, aes(shift)) +
  geom_histogram(
    binwidth = 0.025, fill = "#E41A1C", alpha = 0.70,
    color = "white", linewidth = 0.2
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.4) +
  geom_vline(xintercept = shift_mean, color = "#B2182B", linewidth = 0.4) +
  annotate("text",
    x = shift_mean, y = Inf,
    label = sprintf("mean = %+.3f", shift_mean),
    vjust = 1.4, hjust = -0.1, size = 2.2, color = "#B2182B"
  ) +
  labs(
    x = "Shift (log2)", y = "MNAR proteins", tag = "L",
    title = "MNAR imputed value shift",
    subtitle = sprintf("%d MNAR proteins | missForest", nrow(audit))
  ) +
  FIG_THEME

save_panel(pL, "S1b_L_mnar_shift", "panel_L", audit)

invisible(pL)

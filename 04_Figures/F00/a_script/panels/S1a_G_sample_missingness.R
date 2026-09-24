#!/usr/bin/env Rscript
# S1a Figure G: missing proteins per sample.

setwd(here::here())
source("04_Figures/F00/a_script/panels/_inputs.R", local = TRUE)

miss_bar <- int_norm$miss_bar_data |>
  filter(status == "Missing") |>
  mutate(Group_Time = factor(Group_Time,
    levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post")
  ))

pG <- ggplot(miss_bar, aes(reorder(Col_ID, -n), n, fill = Group_Time)) +
  geom_col(aes(alpha = is_outlier), width = 0.8) +
  scale_fill_manual(
    values = GROUP_FILL, name = NULL,
    labels = c("Young Pre", "Young Post", "Old Pre", "Old Post")
  ) +
  scale_alpha_manual(values = c("FALSE" = 1, "TRUE" = 0.3), guide = "none") +
  labs(
    x = NULL, y = "Missing proteins", tag = "G",
    title = "Per-sample missing protein counts",
    subtitle = sprintf(
      "%d samples | %d outliers (faded) | %.1f%% overall",
      length(unique(miss_bar$Col_ID)),
      int_norm$n_outliers,
      int_imp$pct_miss
    )
  ) +
  FIG_THEME +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    legend.position = "top", legend.key.size = unit(2.5, "mm"),
    legend.text = element_text(size = 5.5)
  )

save_panel(pG, "S1a_G_sample_missingness", "panel_G", miss_bar,
  width = PW * 2, height = PH * 0.75
)

invisible(pG)

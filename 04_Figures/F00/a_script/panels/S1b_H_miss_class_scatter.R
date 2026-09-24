#!/usr/bin/env Rscript
# S1b Figure H: MAR and MNAR classification against intensity.

setwd(here::here())
source("04_Figures/F00/a_script/panels/_inputs.R", local = TRUE)

mc <- int_imp$miss_class |>
  mutate(classification = factor(classification, levels = c("Complete", "MAR", "MNAR")))
n_total <- nrow(mc)
class_lab <- mc |>
  count(classification) |>
  mutate(lab = sprintf("%s: %d (%.0f%%)", classification, n, 100 * n / n_total))

pH <- ggplot(mc, aes(mean_intensity, pct_miss, color = classification)) +
  geom_point(alpha = 0.45, size = 0.7) +
  scale_color_manual(values = PAL_CLASS, name = NULL, labels = class_lab$lab) +
  scale_y_continuous(labels = \(x) paste0(x, "%")) +
  labs(
    x = "Mean log2 intensity", y = "% missing", tag = "H",
    title = "MAR vs. MNAR classification",
    subtitle = sprintf(
      "%s | %s proteins", int_imp$classification_method,
      comma(n_total)
    )
  ) +
  FIG_THEME +
  theme(
    legend.position = "top", legend.key.size = unit(2.5, "mm"),
    legend.text = element_text(size = 5.5)
  )

save_panel(pH, "S1b_H_miss_class_scatter", "panel_H", mc)

invisible(pH)

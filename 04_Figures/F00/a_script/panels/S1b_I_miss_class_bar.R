#!/usr/bin/env Rscript
# S1b Figure I: proteins per missingness class.

setwd(here::here())
source("04_Figures/F00/a_script/panels/_inputs.R", local = TRUE)

mc <- int_imp$miss_class |>
  mutate(classification = factor(classification, levels = c("Complete", "MAR", "MNAR")))

class_counts <- mc |>
  count(classification) |>
  mutate(pct = 100 * n / sum(n))

pI <- ggplot(class_counts, aes(classification, n, fill = classification)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%d\n(%.0f%%)", n, pct)),
    vjust = -0.2, size = 2.2, fontface = "bold"
  ) +
  scale_fill_manual(values = PAL_CLASS, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(
    x = NULL, y = "Proteins", tag = "I",
    title = "Missingness class counts",
    subtitle = sprintf(
      "MAR: %s vals | MNAR: %s vals | %.1f%% total",
      comma(int_imp$mar_miss_vals), comma(int_imp$mnar_miss_vals),
      int_imp$pct_miss
    )
  ) +
  FIG_THEME

save_panel(pI, "S1b_I_miss_class_bar", "panel_I", class_counts)

invisible(pI)

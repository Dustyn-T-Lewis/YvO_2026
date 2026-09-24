#!/usr/bin/env Rscript
# S1a Figure A: protein filter cascade.

setwd(here::here())
source("04_Figures/F00/a_script/panels/_inputs.R", local = TRUE)

fcasc <- int_norm$filter_log |>
  mutate(
    step = factor(step, levels = step),
    removed = ifelse(is.na(n_removed), 0L, n_removed)
  )

pA <- ggplot(fcasc, aes(step, n_after)) +
  geom_col(fill = "#2166AC", width = 0.7) +
  geom_text(aes(label = comma(n_after)), vjust = -1.8, size = 2.2, fontface = "bold") +
  geom_text(aes(label = ifelse(removed > 0, sprintf("-%s", comma(removed)), "")),
    vjust = -0.4, size = 2.0, color = "#B2182B"
  ) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, 16)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18)), labels = comma) +
  labs(
    x = NULL, y = "Proteins retained", tag = "A",
    title = "Protein filter cascade",
    subtitle = sprintf(
      "%s to %s (%.1f%%)", comma(fcasc$n_after[1]),
      comma(tail(fcasc$n_after, 1)), tail(fcasc$pct_of_raw, 1)
    )
  ) +
  FIG_THEME +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5, lineheight = 0.85),
    # Without this the y title runs into the thousands-separated tick labels.
    axis.title.y = element_text(margin = margin(r = 2))
  )

save_panel(pA, "S1a_A_filter_cascade", "panel_A", fcasc)

invisible(pA)

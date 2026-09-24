#!/usr/bin/env Rscript
# Graphical abstract card C: proteins altered by training at FDR and Π.

setwd(here::here())
source("04_Figures/abstract_panels/a_script/panels/_common.R", local = TRUE)

WIDTH <- 1.65
N_SUBJ <- c(Younger = 17L, Older = 15L)

res <- dep_results()

count_sig <- function(group, age) {
  fdr <- res[[paste0("adj.P.Val_", group)]]
  lfc <- res[[paste0("logFC_", group)]]
  pi <- res[[paste0("sig_pi_", group)]]
  tibble(
    age = age,
    criterion = rep(c("FDR < 0.05", "Π < 0.05"), each = 2),
    direction = rep(c("Up", "Down"), 2),
    n = c(
      sum(fdr < 0.05 & lfc > 0, na.rm = TRUE),
      sum(fdr < 0.05 & lfc < 0, na.rm = TRUE),
      sum(pi == 1, na.rm = TRUE),
      sum(pi == -1, na.rm = TRUE)
    )
  )
}

counts <- bind_rows(
  count_sig("Training_Young", "Younger"),
  count_sig("Training_Old", "Older")
) |>
  mutate(
    criterion = factor(criterion, levels = c("FDR < 0.05", "Π < 0.05")),
    direction = factor(direction, levels = c("Up", "Down")),
    facet_age = factor(
      sprintf("%s  n = %d", age, N_SUBJ[age]),
      levels = sprintf("%s  n = %d", names(N_SUBJ), N_SUBJ)
    )
  )

C_TITLES <- c(
  "FDR < 0.05" = "At FDR < 0.05: 135 in\nyounger, none in older",
  "\u03a0 < 0.05" = "At \u03a0 < 0.05: 99 in\nyounger, 18 in older"
)

count_plot <- function(which_criterion) {
  ggplot(
    filter(counts, criterion == which_criterion),
    aes(direction, n, fill = direction)
  ) +
    geom_col(width = 0.55, colour = OUTLINE, linewidth = OUTLINE_W) +
    geom_text(
      aes(label = n), vjust = -0.35,
      size = pt(GLYPH_PT), fontface = "bold", colour = INK
    ) +
    facet_wrap(vars(facet_age)) +
    scale_fill_manual(values = DIR_PAL) +
    scale_x_discrete(expand = expansion(add = 0.55)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
    labs(
      title = C_TITLES[[which_criterion]],
      x = NULL, y = "Proteins altered"
    ) +
    theme_panel() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.text.x = element_text(margin = margin(b = 0.5)),
      panel.spacing.x = unit(7, "pt")
    )
}

plots <- map(levels(counts$criterion), count_plot)
panel_data <- select(counts, age, criterion, direction, n)

save_card(plots, "C_remodelling", WIDTH)
invisible(plots)

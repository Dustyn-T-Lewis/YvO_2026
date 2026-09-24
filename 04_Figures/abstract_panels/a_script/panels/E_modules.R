#!/usr/bin/env Rscript
# Graphical abstract card E: module AUCs for age and for training.

setwd(here::here())
source("04_Figures/abstract_panels/a_script/panels/_common.R", local = TRUE)

WIDTH <- 3.00
XLIM <- c(0, 1.28)
CONTRASTS <- c(
  Age = "Age", Younger = "Pre vs Post (Young)", Older = "Pre vs Post (Old)"
)

grid_res <- module_grid()
age_sig <- grid_res |>
  filter(row == "Age", sig != "ns") |>
  arrange(auc) |>
  pull(module)
train_sig <- grid_res |>
  filter(row == "Pre vs Post (Young)", module %in% age_sig, sig != "ns") |>
  pull(module)

auc <- grid_res |>
  filter(row %in% CONTRASTS, module %in% age_sig) |>
  mutate(
    contrast = factor(
      names(CONTRASTS)[match(row, CONTRASTS)], levels = names(CONTRASTS)
    ),
    panel = factor(
      if_else(contrast == "Age", "Age", "Training"),
      levels = c("Age", "Training")
    ),
    tag = factor(tag, levels = MODULE_TAG[age_sig])
  ) |>
  filter(panel == "Age" | module %in% train_sig)

# Both blocks draw one bar per row at the same width, so the Training block
# stacks its two contrasts as separate rows rather than dodging them inside
# one. Dodging halves the bar width, which left the two blocks different sizes.
auc_base <- function(dat, title) {
  ggplot(dat, aes(auc, y = .data$row_label, fill = module)) +
    geom_col(width = 0.74, colour = OUTLINE, linewidth = 0.16) +
    geom_text(
      aes(x = auc + 0.02, label = trimws(sprintf("%.2f %s", auc, mark))),
      hjust = 0, size = pt(GLYPH_PT - 1), fontface = "bold", colour = INK
    ) +
    scale_fill_manual(values = MODULE_FILL) +
    scale_x_continuous(limits = XLIM, breaks = c(0, 1), expand = expansion(0)) +
    labs(x = "AUC", y = NULL, title = title) +
    theme_panel() +
    theme(
      axis.text.y = element_text(size = GLYPH_PT, face = "bold", colour = INK)
    )
}

age_rows <- auc |>
  filter(panel == "Age") |>
  mutate(row_label = droplevels(tag))

# One row per module and contrast, ordered so each module's Younger bar sits
# directly above its Older one. Levels must be unique per row or the two bars
# of a pair land on the same row, so the key carries the module and the
# contrast and the axis labels are mapped back onto it.
train_rows <- auc |>
  filter(panel == "Training") |>
  mutate(
    row_key = paste(tag, contrast, sep = "_"),
    row_label = factor(
      row_key,
      levels = paste(
        rep(levels(droplevels(tag)), each = 2),
        c("Older", "Younger"),
        sep = "_"
      ) |>
        (\(x) x[x %in% paste(tag, contrast, sep = "_")])()
    )
  )

age_plot <- auc_base(
  age_rows, "Module scores separate\nyounger from older"
)

train_plot <- auc_base(
  train_rows, "Modules track training\nonly in the younger"
) +
  aes(alpha = contrast) +
  scale_alpha_manual(values = c(Younger = 0.95, Older = 0.4)) +
  scale_y_discrete(
    labels = \(k) sub("_Younger$", "", sub(".*_Older$", "vs older", k))
  )

plots <- list(age_plot, train_plot)
panel_data <- select(auc, panel, contrast, module, tag, auc, perm_p, q_bh, sig)

save_card(plots, "E_modules", WIDTH)
invisible(plots)

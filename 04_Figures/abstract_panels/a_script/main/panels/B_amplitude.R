#!/usr/bin/env Rscript
# Graphical abstract card B: protein change magnitude, younger against older.

setwd(here::here())
source("04_Figures/abstract_panels/a_script/main/panels/_common.R", local = TRUE)

WIDTH <- 2.20
XLIM <- c(0.45, 5.72)
YLIM <- c(0, 0.56)
CAP <- 0.035
PCT_PT <- 5.5

# The sets and their magnitudes come from F04 panel C, which cross-fits the
# selection: half the young subjects choose the members, the other half supply
# the young magnitude. Choosing and measuring on the same subjects inflates the
# young side, so a naive recomputation here would overstate every gap.
sets <- read_excel(F04_XLSX, "panel_C_trajectory") |>
  mutate(
    set_id = row_number(),
    set = c("All", "Top 500", "Top 250", "All FDR", "All \u03a0"),
    Young = young_magnitude,
    Old = old_magnitude,
    xb = set_id,
    pct = percent(1 - retained, accuracy = 1)
  )

long <- sets |>
  pivot_longer(c(Young, Old), names_to = "age", values_to = "med") |>
  mutate(age = factor(age, levels = c("Young", "Old")))

glyph <- ggplot(long, aes(set_id, med, fill = age)) +
  geom_col(
    position = position_dodge(width = 0.82), width = 0.5,
    colour = OUTLINE, linewidth = OUTLINE_W
  ) +
  geom_segment(
    data = sets, inherit.aes = FALSE,
    aes(x = xb, xend = xb, y = Old, yend = Young),
    colour = INK, linewidth = 0.3
  ) +
  geom_segment(
    data = sets, inherit.aes = FALSE,
    aes(x = xb - CAP, xend = xb + CAP, y = Young, yend = Young),
    colour = INK, linewidth = 0.3
  ) +
  geom_segment(
    data = sets, inherit.aes = FALSE,
    aes(x = xb - CAP, xend = xb + CAP, y = Old, yend = Old),
    colour = INK, linewidth = 0.3
  ) +
  geom_text(
    data = sets, inherit.aes = FALSE,
    aes(xb + CAP + 0.10, (Young + Old) / 2, label = pct),
    hjust = 0, size = pt(PCT_PT), fontface = "bold", colour = INK
  ) +
  panel_key(c("Younger", "Older"), unname(AGE_PAL), XLIM, YLIM, x = 0.02) +
  scale_fill_manual(values = AGE_PAL) +
  scale_x_continuous(
    limits = XLIM, breaks = sets$set_id, labels = sets$set,
    expand = expansion(0)
  ) +
  scale_y_continuous(
    limits = YLIM, breaks = c(0, 0.5), expand = expansion(0),
    labels = \(x) drop0(formatC(x, format = "fg"))
  ) +
  labs(
    title = "Proteins change 35\u201364%\nless in older adults",
    x = NULL, y = expression(bold(group("|", log[2] * "FC", "|")))
  ) +
  theme_panel() +
  theme(
    axis.text.y = element_text(size = GLYPH_PT - 1.5),
    axis.text.x = element_text(size = GLYPH_PT - 1),
    axis.title.y = element_text(margin = margin(r = 0))
  )

plots <- list(glyph)
panel_data <- select(
  sets, set, n, median_young = Young, median_old = Old,
  retained, retained_lo, retained_hi, pct
)

save_card(plots, "B_amplitude", WIDTH)
invisible(plots)

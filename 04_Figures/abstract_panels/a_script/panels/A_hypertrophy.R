#!/usr/bin/env Rscript
# Graphical abstract card A: change in VL thickness, younger against older.

setwd(here::here())
source("04_Figures/abstract_panels/a_script/panels/_common.R", local = TRUE)

WIDTH <- 2.20
# The bars, their values and the bracket all finish by CHART_MAX; the rest of
# the panel is left empty for the VL illustration dropped in by hand.
CHART_MAX <- 0.97
ART_FRAC <- 0.40
XLIM <- c(0, CHART_MAX / (1 - ART_FRAC))
BR_X <- 0.84
TIP <- 0.03
BAR_W <- 0.55
AX_BREAKS <- c(0, 0.25, 0.5)

vl <- phenotype() |>
  filter(!is.na(VL_thick_cm_Delta)) |>
  summarise(n = n(), delta = mean(VL_thick_cm_Delta), .by = Age_group) |>
  mutate(
    age = factor(
      if_else(Age_group == "Young", "Younger", "Older"),
      levels = c("Younger", "Older")
    ),
    ypos = if_else(age == "Younger", 2, 1),
    label = sprintf("%s  n = %d", age, n)
  ) |>
  arrange(age)

glyph <- ggplot(vl, aes(delta, age, fill = age)) +
  geom_col(width = BAR_W, colour = OUTLINE, linewidth = OUTLINE_W) +
  geom_text(
    inherit.aes = FALSE,
    aes(0.015, ypos + BAR_W / 2 + 0.06, label = label),
    hjust = 0, vjust = 0, size = pt(GLYPH_PT), fontface = "bold", colour = INK
  ) +
  geom_text(
    aes(x = delta + 0.025, label = drop0(sprintf("+%.2f", delta))),
    hjust = 0, size = pt(GLYPH_PT), fontface = "bold", colour = INK
  ) +
  annotate(
    "segment", x = BR_X, xend = BR_X, y = 1, yend = 2,
    colour = INK, linewidth = 0.3
  ) +
  annotate(
    "segment", x = BR_X, xend = BR_X - TIP, y = c(1, 2), yend = c(1, 2),
    colour = INK, linewidth = 0.3
  ) +
  annotate(
    "text", x = BR_X + 0.02, y = 1.5, label = "**", hjust = 0, vjust = 0.5,
    size = pt(GLYPH_PT + 1), fontface = "bold", colour = INK
  ) +
  annotate(
    "segment", x = 0, xend = max(AX_BREAKS), y = 0.5, yend = 0.5,
    colour = MUTED, linewidth = 0.35
  ) +
  scale_fill_manual(
    values = c(Younger = AGE_PAL[["Young"]], Older = AGE_PAL[["Old"]])
  ) +
  scale_x_continuous(
    limits = XLIM, breaks = AX_BREAKS, expand = expansion(0),
    labels = \(x) drop0(formatC(x, format = "fg"))
  ) +
  scale_y_discrete(limits = rev, expand = expansion(add = 0.5)) +
  labs(
    title = "Younger gain 2.4x the\nmuscle older gain",
    x = expression(bold(Delta * " VL (cm)")), y = NULL
  ) +
  theme_panel() +
  theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.line.x = element_blank(),
    # centred on the axis rather than on the panel, so the title stays out
    # of the space reserved to the right of it
    axis.title.x = element_text(
      hjust = mean(range(AX_BREAKS)) / XLIM[2]
    )
  )

plots <- list(glyph)
panel_data <- select(vl, age, n, delta)

save_card(plots, "A_hypertrophy", WIDTH)
invisible(plots)

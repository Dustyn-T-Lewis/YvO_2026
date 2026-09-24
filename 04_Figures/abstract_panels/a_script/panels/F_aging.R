#!/usr/bin/env Rscript
# Graphical abstract card F: the aging volcano. Drawn only in the concise row.

setwd(here::here())
source("04_Figures/abstract_panels/a_script/panels/_common.R", local = TRUE)

WIDTH <- 1.45

v <- dep_results() |>
  select(gene, lfc = logFC_Aging, p = P.Value_Aging, fdr = adj.P.Val_Aging) |>
  filter(!is.na(lfc), !is.na(p)) |>
  mutate(
    nlp = -log10(p),
    # The archive coloured this card on AGE_PAL. DIR_PAL instead, because red
    # already means Up on the remodelling card and red meaning Younger on one
    # card and Up on another is the collision _common.R:60 warns about.
    call = case_when(
      fdr < 0.05 & lfc > 0 ~ "Up",
      fdr < 0.05 & lfc < 0 ~ "Down",
      .default = "ns"
    )
  )

counts <- c(Up = sum(v$call == "Up"), Down = sum(v$call == "Down"))
stopifnot(sum(counts) == sum(v$fdr < 0.05, na.rm = TRUE))

# The dashed line is where the BH threshold lands in raw p, so it moves with
# the data rather than sitting at a number typed in here.
p_cut <- max(v$p[v$fdr < 0.05], na.rm = TRUE)

XLIM <- c(-1, 1) * max(abs(v$lfc)) * 1.04
YLIM <- c(0, max(v$nlp) * 1.04)

glyph <- ggplot(v, aes(lfc, nlp)) +
  geom_vline(xintercept = 0, colour = MUTED, linewidth = 0.3) +
  geom_hline(
    yintercept = -log10(p_cut), colour = MUTED,
    linetype = "22", linewidth = 0.4
  ) +
  geom_point(
    data = filter(v, call == "ns"),
    colour = "grey65", size = 0.5, alpha = 0.5, shape = 16
  ) +
  geom_point(
    data = filter(v, call != "ns"), aes(fill = call),
    shape = 21, size = 0.9, stroke = 0.12, colour = OUTLINE
  ) +
  panel_key(
    sprintf("%d %s with age", counts, c("up", "down")),
    unname(DIR_PAL[names(counts)]), XLIM, YLIM
  ) +
  scale_fill_manual(values = DIR_PAL) +
  coord_cartesian(xlim = XLIM, ylim = YLIM, expand = FALSE) +
  scale_x_continuous(breaks = 0) +
  labs(
    title = sprintf(
      "At FDR < 0.05: %d differ\nbetween younger and older", sum(counts)
    ),
    x = expression(bold("aging log"[2] * "FC")),
    y = expression(bold("-log"[10] * "P"))
  ) +
  theme_panel()

plots <- list(glyph)
panel_data <- v |>
  filter(call != "ns") |>
  select(gene, logFC_aging = lfc, p = p, fdr = fdr, direction = call)

save_card(plots, "F_aging", WIDTH)
invisible(plots)

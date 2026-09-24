#!/usr/bin/env Rscript
# S5 Table, sheet SUPP_cat_depth: corner enrichment against list depth (CAT curve)
# Diagnostic for main Panel B: the corner hypergeometric reads one point off
# this curve, at TOP_FRAC = 0.10 in shared/comparison_panels/panel_E_rrho2.R.
# Drawing every depth shows what that one number cannot: the enrichment halves
# between the top 10% and the top 25%, and is gone by half the proteome, where
# the corner test degenerates into the sign concordance main Panel A reports.
# Restricted to the two training contrasts. The Aging x Training(Old) pair
# shares Old_Pre, so its chance line would be wrong, and a CAT plot draws that
# line explicitly where the RRHO map leaves it implicit.

setwd(here::here())

pacman::p_load(dplyr, tidyr, tibble, stringr, readr, ggplot2, patchwork, cowplot)

source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

pdf_device <- get_pdf_device()

BASE <- "04_Figures/F04"
DAT <- file.path(BASE, "c_data", "panel_supp")
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

TOP_FRAC <- 0.10

cat_dep <- read_csv("03_DEP/c_data/03_combined_results.csv",
  show_col_types = FALSE
) |>
  select(ty = t_Training_Young, to = t_Training_Old) |>
  filter(!is.na(ty), !is.na(to))

n_cat <- nrow(cat_dep)
r_y_up <- rank(-cat_dep$ty, ties.method = "first")
r_o_up <- rank(-cat_dep$to, ties.method = "first")
r_y_dn <- rank(cat_dep$ty, ties.method = "first")
r_o_dn <- rank(cat_dep$to, ties.method = "first")

# Running intersection size at every depth: sorting once and accumulating is
# O(n), against O(n^2) for a loop that re-intersects at each k.
running_overlap <- function(r1, r2) {
  cumsum(tabulate(pmax(r1, r2), nbins = n_cat))
}

cat_curves <- tibble(
  k = seq_len(n_cat),
  `Concordant up` = running_overlap(r_y_up, r_o_up),
  `Concordant down` = running_overlap(r_y_dn, r_o_dn),
  Discordant = running_overlap(r_y_up, r_o_dn) +
    running_overlap(r_y_dn, r_o_up)
) |>
  pivot_longer(-k, names_to = "corner", values_to = "overlap") |>
  mutate(
    corner = factor(corner, levels = c(
      "Concordant up", "Concordant down", "Discordant"
    )),
    # Two disjoint corners are pooled into Discordant, so it gets twice the
    # chance overlap a single corner does.
    expected = if_else(corner == "Discordant", 2, 1) * k^2 / n_cat,
    fold = overlap / expected,
    depth_pct = 100 * k / n_cat
  )

K_CAT <- round(TOP_FRAC * n_cat)
cat_mark <- filter(cat_curves, k == K_CAT)
cat_p <- cat_mark |>
  mutate(
    p = if_else(
      overlap > expected,
      phyper(overlap - 1, K_CAT, n_cat - K_CAT, K_CAT, lower.tail = FALSE),
      phyper(overlap, 2 * K_CAT, n_cat - 2 * K_CAT, K_CAT, lower.tail = TRUE)
    )
  )

CAT_PAL <- c(
  `Concordant up` = DIR_COLORS[["Up"]],
  `Concordant down` = DIR_COLORS[["Down"]],
  Discordant = "#8C7AA8"
)
# A log axis cannot draw an empty corner, and the discordant one is empty at
# the shallowest depths. The curve breaks there rather than the panel showing
# a finite ratio it does not have.
cat_draw <- filter(cat_curves, k >= 50, fold > 0)

pS_cat <- ggplot(cat_draw, aes(depth_pct, fold, colour = corner)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55") +
  geom_vline(
    xintercept = 100 * TOP_FRAC, linetype = "dotted", colour = "grey35",
    linewidth = 0.4
  ) +
  geom_line(linewidth = 0.45) +
  geom_point(data = cat_mark, size = 0.9) +
  annotate("text",
    x = 97, y = c(5.2, 0.19), hjust = 1, size = BASE_STAT,
    colour = "grey40", label = c("enriched", "depleted")
  ) +
  annotate("text",
    x = 100 * TOP_FRAC + 1.5, y = 7.2, hjust = 0, size = BASE_STAT,
    colour = "grey25",
    label = sprintf("corner test\n(top %.0f%%)", 100 * TOP_FRAC)
  ) +
  scale_colour_manual(values = CAT_PAL, name = NULL) +
  scale_x_continuous(
    limits = c(0, 100), breaks = c(0, 10, 25, 50, 75, 100),
    labels = \(x) paste0(x, "%"), expand = expansion(0)
  ) +
  scale_y_continuous(
    transform = "log2", breaks = 2^(-3:3),
    labels = c("1/8", "1/4", "1/2", "1", "2", "4", "8"),
    limits = c(2^-3.4, 2^3.2), oob = scales::squish
  ) +
  labs(
    title = "Corner enrichment against list depth",
    subtitle = sprintf(
      paste(
        "CAT curve; at the top %.0f%% (k = %d) the corner test reads",
        "%.1fx up, %.1fx down"
      ),
      100 * TOP_FRAC, K_CAT, cat_mark$fold[1], cat_mark$fold[2]
    ),
    x = "List depth (% of the ranked proteome, both contrasts)",
    y = "Overlap relative to chance"
  ) +
  FIG_THEME +
  theme(
    legend.position = "bottom",
    legend.key.height = unit(5, "pt"),
    panel.grid.minor = element_blank()
  )

# A 1% grid redraws the curve without carrying 6,282 rows into the workbook.
write_csv(
  cat_curves |>
    filter(k %in% c(K_CAT, round(seq(0.01, 1, 0.01) * n_cat))) |>
    select(k, depth_pct, corner, overlap, expected, fold) |>
    arrange(k, corner),
  file.path(DAT, "SUPP_cat_depth.csv")
)

RPT_PNG <- file.path(BASE, "b_reports", "supp", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "panels")
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)

PW_CAT <- 89
PH_CAT <- 70
ggsave(file.path(RPT_PNG, "S5T_cat_depth.png"), pS_cat,
  width = PW_CAT, height = PH_CAT, units = "mm", dpi = 300
)
ggsave(file.path(RPT_PDF, "S5T_cat_depth.pdf"), pS_cat,
  width = PW_CAT, height = PH_CAT, units = "mm", device = pdf_device
)

message(sprintf(
  paste(
    "SUPP CAT depth: k = %d (top %.0f%%) | %s |",
    "half-enrichment near %.0f%% depth"
  ),
  K_CAT, 100 * TOP_FRAC,
  paste(sprintf("%s %.1fx p=%.2g", cat_p$corner, cat_p$fold, cat_p$p),
    collapse = " | "
  ),
  cat_curves |>
    filter(
      corner == "Concordant up", k > K_CAT,
      fold <= 1 + (cat_mark$fold[1] - 1) / 2
    ) |>
    slice_min(k, n = 1) |>
    pull(depth_pct)
))

invisible(pS_cat)

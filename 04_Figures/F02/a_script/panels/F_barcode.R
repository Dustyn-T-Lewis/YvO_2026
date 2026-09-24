#!/usr/bin/env Rscript
# Figure 2F: DEP rank location (barcode plot).
# Shows where DEPs sit in the t-statistic-ranked proteome. Density (dark,
# filled) and its peak labels are Pi < 0.05 -- the richer set, so density is
# always drawable even where FDR is sparse or empty. Two independent tick
# tracks below: Pi (upper, near the density) and FDR (lower, separate band).

setwd(here::here())
source("04_Figures/F02/a_script/panels/_main.R", local = TRUE)
DEP_FILE <- "03_DEP/c_data/03_combined_results.csv"
RPT_PNG <- "04_Figures/F02/b_reports/panels"
RPT_PDF <- "04_Figures/F02/b_reports/panels"
DAT <- "04_Figures/F02/c_data"
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

CONTRASTS <- c("Aging", "Training_Young", "Training_Old", "Interaction")
dep_df <- read_csv(DEP_FILE, show_col_types = FALSE)
pdf_device <- get_pdf_device()

PD_W <- 67 # J Physiol: col 3 of 3×2 at 178mm
PD_H <- 55

# Build long-form data: rank position + DEP status per contrast
rank_list <- lapply(CONTRASTS, function(ctr) {
  t_col <- paste0("t_", ctr)
  pi_col <- paste0("pi_score_", ctr)
  lfc_col <- paste0("logFC_", ctr)

  fdr_col <- paste0("adj.P.Val_", ctr)

  dep_df |>
    filter(!is.na(.data[[t_col]])) |>
    arrange(.data[[t_col]]) |>
    mutate(
      rank_frac = seq_len(n()) / n(),
      is_fdr = !is.na(.data[[fdr_col]]) & .data[[fdr_col]] < 0.05,
      is_pi = !is.na(.data[[pi_col]]) & .data[[pi_col]] < 0.05,
      is_dep = is_fdr | is_pi,
      tier = if_else(is_fdr, "FDR", "Pi only"),
      direction = case_when(
        !is_dep ~ NA_character_,
        .data[[lfc_col]] > 0 ~ "Up",
        TRUE ~ "Down"
      ),
      contrast = ctr
    ) |>
    select(gene, contrast, rank_frac, is_dep, is_fdr, is_pi, tier, direction)
})
rank_df <- bind_rows(rank_list)
rank_df$contrast <- factor(rank_df$contrast, levels = CONTRASTS)

dep_only <- rank_df |> filter(is_dep)
dep_only$direction <- factor(dep_only$direction, levels = c("Up", "Down"))
pi_all <- dep_only |> filter(is_pi)
fdr_only <- dep_only |> filter(is_fdr)

# Two independent, non-overlapping tick tracks -- Pi (the richer, primary
# set) nearer the density, FDR a separate band below. Each protein ticks in
# every track it qualifies for; the bands don't share an origin, so a dense
# run of Pi ticks can't visually swallow the FDR track the way nested depths
# from a shared origin did.
PI_TOP <- 0
PI_BOT <- -0.14
FDR_TOP <- -0.18
FDR_BOT <- -0.29

tick_df <- bind_rows(
  pi_all |> mutate(band = "Π", y_start = PI_TOP, y_end = PI_BOT),
  fdr_only |> mutate(band = "FDR", y_start = FDR_TOP, y_end = FDR_BOT)
)

dep_counts <- dep_only |>
  group_by(contrast, tier) |>
  summarise(
    n_up    = sum(direction == "Up"),
    n_down  = sum(direction == "Down"),
    n_total = n(),
    .groups = "drop"
  )

write.csv(dep_counts, file.path(DAT, "panel_F_barcode_enrichment.csv"), row.names = FALSE)

# Pre-compute density curves so we can normalize and control y-range. Keyed
# on Pi rather than FDR -- FDR is sparse enough in some contrasts (Tr.(O),
# Interaction) that a density estimate isn't drawable at all.
DENS_PAD <- 0.06
dens_list <- lapply(split(pi_all, pi_all$contrast), function(ctr_df) {
  lapply(split(ctr_df, ctr_df$direction, drop = TRUE), function(dir_df) {
    if (nrow(dir_df) < 2) {
      return(NULL)
    }
    d <- density(dir_df$rank_frac,
      adjust = 1.8,
      from = -DENS_PAD, to = 1 + DENS_PAD, n = 512
    )
    tibble(
      x = d$x, y = d$y, direction = dir_df$direction[1],
      contrast = dir_df$contrast[1]
    )
  }) |> bind_rows()
}) |> bind_rows()

# Normalize density per contrast: peak = 1.0 (makes panels comparable).
# A Gaussian kernel never reaches exactly zero, so the raw curve traces a
# thin, visible line across the entire x-axis (0-100%) even where virtually
# no proteins sit -- e.g. a Down curve peaked near 0% still draws a hairline
# out past 90%. Blanking the near-zero tail (NA, not filtering rows, so
# geom_line/geom_ribbon actually break instead of connecting across the gap)
# stops the curve where it's no longer meaningfully above baseline.
dens_list <- dens_list |>
  group_by(contrast) |>
  mutate(
    y_norm = y / max(y),
    y_norm = if_else(y_norm < 0.02, NA_real_, y_norm)
  ) |>
  ungroup()
dens_list$direction <- factor(dens_list$direction, levels = c("Up", "Down"))
dens_list$contrast <- factor(dens_list$contrast, levels = CONTRASTS)

ANNOT_SZ <- FIG_AXIS_TEXT / .pt
LABEL_NUDGE <- 0.06

# The two tick bands are told apart only by height and alpha, so name them.
# No contrast column, so the pair repeats in every facet. Together the bands
# are barely 7 pt tall on the composite, and that is what caps the type size:
# any larger and the two labels' boxes touch. The nudge lifts both off their
# arithmetic centre because a label box hangs a little below its anchor --
# without it Pi and FDR close to within a quarter point of each other.
# 3.25 and no lower: the figure prints at 165.1 mm of a 178 mm canvas, so
# 3.0 here lands at 2.78 pt on the page, under the 3 pt floor.
BAND_LBL_SZ <- 3.25 / .pt
BAND_LBL_LIFT <- 0.017
band_lbl_df <- tibble(
  x = 1 + DENS_PAD - 0.002,
  y = c(mean(c(PI_TOP, PI_BOT)), mean(c(FDR_TOP, FDR_BOT))) + BAND_LBL_LIFT,
  label = c("Π", "FDR")
)

# Compute peak positions for label placement
peak_pos <- dens_list |>
  group_by(contrast, direction) |>
  slice_max(y_norm, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(contrast, direction, peak_x = x, peak_y = y_norm)

# Pi DEP counts per contrast x direction -- the density's own set. Contrasts
# with too few Pi hits for a density (essentially never, in practice) state
# their count in place of density labels.
n_down <- pi_all |>
  filter(direction == "Down") |>
  count(contrast) |>
  tibble::deframe()
n_up <- pi_all |>
  filter(direction == "Up") |>
  count(contrast) |>
  tibble::deframe()

sparse_labels <- rank_df |>
  group_by(contrast) |>
  summarise(n_fdr = sum(is_fdr), n_pi = sum(is_pi), .groups = "drop") |>
  filter(n_pi < 2) |>
  mutate(label = sprintf("%d at Π < 0.05  (FDR: %d)", n_pi, n_fdr))
DESC_DOWN <- c(
  Aging = "proteins lower in older vs young",
  Training_Young = "proteins dec. with training",
  Training_Old = "proteins dec. with training",
  Interaction = "proteins with greater Young response"
)
DESC_UP <- c(
  Aging = "proteins higher in older vs young",
  Training_Young = "proteins inc. with training",
  Training_Old = "proteins inc. with training",
  Interaction = "proteins with greater Old response"
)

bg_wash <- tibble(
  contrast = factor(CONTRASTS, levels = CONTRASTS),
  fill = unname(CONTRAST_COLORS[CONTRASTS]),
  xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
)

ad_all <- peak_pos |>
  filter(direction == "Down") |>
  mutate(
    ctr      = as.character(contrast),
    label_x  = peak_x + LABEL_NUDGE,
    # Down runs right from the left peak and Up runs left from the right peak,
    # so long labels meet mid-panel in every contrast. Vertical separation is
    # what keeps them apart, and it has to be absolute: as a fraction of
    # peak_y the gap collapses in contrasts with a short density (Tr.(O)).
    label_y  = 0.38,
    label    = paste(n_down[ctr], DESC_DOWN[ctr])
  ) |>
  filter(!is.na(label))

au_all <- peak_pos |>
  filter(direction == "Up") |>
  mutate(
    ctr      = as.character(contrast),
    label_x  = peak_x - LABEL_NUDGE,
    label_y  = 0.90,
    label    = paste(n_up[ctr], DESC_UP[ctr])
  ) |>
  filter(!is.na(label))

cd_all <- ad_all |>
  mutate(x_start = peak_x, y_start = peak_y, x_end = label_x, y_end = label_y)
cu_all <- au_all |>
  mutate(x_start = peak_x, y_start = peak_y, x_end = label_x, y_end = label_y)

pF <- ggplot() +
  # Background contrast wash (per facet) — darkened to match C/D/E's
  # canonical 0.20 alpha (0.18 is a slight pull-back to avoid over-darkening)
  geom_rect(
    data = bg_wash,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = bg_wash$fill, alpha = 0.18, inherit.aes = FALSE
  ) +
  # Density ribbons -- Pi-based, filled dark and solid to read as the
  # primary signal (ticks below carry the FDR/Pi tier split)
  geom_ribbon(
    data = dens_list,
    aes(x = x, ymin = 0, ymax = y_norm, fill = direction),
    alpha = 0.60, outline.type = "upper"
  ) +
  geom_line(
    data = dens_list,
    aes(x = x, y = y_norm, color = direction),
    linewidth = 0.6
  ) +
  # Barcode ticks: two independent tracks, Pi (all Pi-sig genes) nearer the
  # density, FDR (all FDR-sig genes) in its own band below. FDR drawn dark
  # (the stricter, headline tier), Pi light; narrow linewidth keeps a dense
  # run of ticks from reading as a solid block.
  geom_segment(
    data = tick_df,
    aes(
      x = rank_frac, xend = rank_frac,
      y = y_start, yend = y_end,
      color = direction, alpha = band
    ),
    linewidth = 0.22
  ) +
  scale_alpha_manual(values = c("Π" = 0.45, FDR = 0.9), guide = "none") +
  # Band names, right-anchored inside the panel. Plain text rather than a
  # plated label: the ticks thin out well before the right edge, so the plate
  # was covering nothing. The Up callouts are right-anchored too but sit at
  # y = 0.90, a long way above these.
  geom_text(
    data = band_lbl_df,
    aes(x = x, y = y, label = label),
    hjust = 1, vjust = 0.5, size = BAND_LBL_SZ,
    color = "black", fontface = "bold",
    inherit.aes = FALSE
  ) +
  geom_label(
    data = sparse_labels,
    aes(x = 0.5, y = 0.55, label = label),
    hjust = 0.5, vjust = 0.5, size = ANNOT_SZ,
    fill = alpha("white", 0.85), color = "grey25",
    fontface = "bold.italic", linewidth = 0,
    label.padding = unit(0.10, "lines"),
    inherit.aes = FALSE
  ) +
  # Zero line
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey50") +
  # Down connector segments
  {
    if (nrow(cd_all) > 0) {
      geom_segment(
        data = cd_all,
        aes(x = x_start, xend = x_end, y = y_start, yend = y_end),
        linewidth = 0.3, color = unname(DIR_COLORS["Down"]),
        alpha = 0.4, inherit.aes = FALSE
      )
    }
  } +
  # Down labels (white text in blue box)
  {
    if (nrow(ad_all) > 0) {
      geom_label(
        data = ad_all,
        aes(x = label_x, y = label_y, label = label),
        hjust = 0, vjust = 0.5, size = ANNOT_SZ,
        fill = unname(DIR_COLORS["Down"]), color = "white",
        fontface = "bold", linewidth = 0,
        label.padding = unit(0.08, "lines"),
        inherit.aes = FALSE
      )
    }
  } +
  # Up connector segments
  {
    if (nrow(cu_all) > 0) {
      geom_segment(
        data = cu_all,
        aes(x = x_start, xend = x_end, y = y_start, yend = y_end),
        linewidth = 0.3, color = unname(DIR_COLORS["Up"]),
        alpha = 0.4, inherit.aes = FALSE
      )
    }
  } +
  # Up labels (white text in red box)
  {
    if (nrow(au_all) > 0) {
      geom_label(
        data = au_all,
        aes(x = label_x, y = label_y, label = label),
        hjust = 1, vjust = 0.5, size = ANNOT_SZ,
        fill = unname(DIR_COLORS["Up"]), color = "white",
        fontface = "bold", linewidth = 0,
        label.padding = unit(0.08, "lines"),
        inherit.aes = FALSE
      )
    }
  } +
  # Scales
  scale_fill_manual(values = c(
    Up = unname(DIR_COLORS["Up"]),
    Down = unname(DIR_COLORS["Down"])
  )) +
  scale_color_manual(values = c(
    Up = unname(DIR_COLORS["Up"]),
    Down = unname(DIR_COLORS["Down"])
  )) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(
    xlim = c(-DENS_PAD, 1 + DENS_PAD),
    ylim = c(FDR_BOT, 1.02)
  ) +
  facet_grid(contrast ~ .,
    switch = "y",
    labeller = labeller(contrast = CTR_FACET)
  ) +
  labs(
    title = "DEP Rank Location",
    subtitle = sprintf(
      "%s genes, t-ranked | density = Π, ticks = Π + FDR",
      format(length(unique(rank_df$gene)), big.mark = ",")
    ),
    x = "Rank position (by t-statistic)", y = NULL
  ) +
  FIG_THEME +
  theme(
    plot.subtitle = element_text(
      size = FIG_SUBTITLE_SIZE,
      face = "bold.italic", color = "grey40"
    ),
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text.y.left = element_text(
      face = "bold", size = FIG_AXIS_TEXT - 0.5,
      angle = 0, hjust = 1
    ),
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.y = unit(2, "pt")
  )

ggsave(file.path(RPT_PNG, "F_barcode.png"), pF,
  width = PD_W, height = PD_H, units = "mm", dpi = 300
)
ggsave(file.path(RPT_PDF, "F_barcode.pdf"), pF,
  width = PD_W, height = PD_H, units = "mm", device = pdf_device
)

invisible(pF)

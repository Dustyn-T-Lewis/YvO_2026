# Shared NES Scatter Panel: Panel D in F04 (concordance) and F05 (reversal)

source("04_Figures/shared/style.R")
source("04_Figures/shared/print_scale_apply.R")
source("04_Figures/shared/pathway_utils.R")

pacman::p_load(tidyverse, ggrepel)

PG_W <- cfg$panel_w %||% 146
RPT_PNG <- cfg$rpt_png
RPT_PDF <- cfg$rpt_pdf
DAT <- cfg$dat
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_D"), recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

fgsea_cache <- FGSEA_CACHE
stopifnot("fGSEA cache missing" = file.exists(fgsea_cache))
fgsea_all <- read_csv(fgsea_cache, show_col_types = FALSE)

cx <- cfg$contrast_x
cy <- cfg$contrast_y
nes_x <- paste0("NES_", cx)
nes_y <- paste0("NES_", cy)
padj_x <- paste0("padj_", cx)
padj_y <- paste0("padj_", cy)
size_x <- paste0("size_", cx)
size_y <- paste0("size_", cy)

fgsea_hg <- fgsea_all |>
  filter(
    database %in% c("Hallmark", "GO Slim"),
    contrast %in% c(cx, cy)
  )

fgsea_wide <- fgsea_hg |>
  dplyr::select(pathway, contrast, NES, padj, size, database) |>
  pivot_wider(
    id_cols = c(pathway, database), names_from = contrast,
    values_from = c(NES, padj, size)
  ) |>
  filter(!is.na(.data[[nes_x]]), !is.na(.data[[nes_y]])) |>
  mutate(set_size = coalesce(.data[[size_x]], .data[[size_y]]))

fgsea_wide <- fgsea_wide |>
  mutate(
    sig_1 = !is.na(.data[[padj_x]]) & .data[[padj_x]] < 0.05,
    sig_2 = !is.na(.data[[padj_y]]) & .data[[padj_y]] < 0.05,
    significance = case_when(
      sig_1 & sig_2 ~ cfg$quadrant_defs$sig_both_label,
      sig_1 ~ cfg$quadrant_defs$sig_x_label,
      sig_2 ~ cfg$quadrant_defs$sig_y_label,
      TRUE ~ "NS"
    ) |> factor(levels = names(cfg$sig_colors)),
    pathway_label = clean_pathway_name(pathway),
    db_shape = ifelse(database == "Hallmark", 24, 21)
  )

fgsea_sig <- fgsea_wide |> filter(significance != "NS")

# Redundancy must be recomputed here, not read from the cache's dedup_status:
# that flag was assigned against the full five-database pool, so a GO Slim
# term can be "redundant, merged into a Reactome term" that this panel
# doesn't even show (it's restricted to Hallmark + GO Slim above). Reusing
# the global flag would silently drop terms that are the only representative
# of their kind within what's actually plotted here.
pw_full <- build_pathway_collection(min_size = 15, max_size = 500, include_goslim = TRUE)
pw_hg <- pw_full[classify_database(names(pw_full)) %in% c("Hallmark", "GO Slim")]

dedup_input <- as.data.frame(fgsea_sig)
dedup_input$padj <- pmin(dedup_input[[padj_x]], dedup_input[[padj_y]], na.rm = TRUE)
dedup_flags <- em_dedup_report(dedup_input, pw_hg, cutoff = 0.5)[
  , c("pathway", "dedup_status", "merged_into", "overlap_jaccard")
]
fgsea_sig <- fgsea_sig |> left_join(dedup_flags, by = "pathway")
fgsea_kept <- fgsea_sig |> filter(dedup_status == "kept")

message(sprintf(
  "  %d total pathways (Hallmark: %d, GO Slim: %d) | %d significant -> %d representatives",
  nrow(fgsea_wide),
  sum(fgsea_wide$database == "Hallmark"),
  sum(fgsea_wide$database == "GO Slim"),
  nrow(fgsea_sig), nrow(fgsea_kept)
))

nes_cor_all <- cor.test(fgsea_wide[[nes_x]], fgsea_wide[[nes_y]], method = "spearman")
nes_ci_all <- fisher_z_ci(nes_cor_all$estimate, nrow(fgsea_wide))
nes_cor_sig <- if (nrow(fgsea_kept) >= 3) {
  cor.test(fgsea_kept[[nes_x]], fgsea_kept[[nes_y]], method = "spearman")
} else {
  NULL
}

# The frame runs half again past the furthest point: the labels live in that
# margin, so widening it is what buys room rather than shrinking the type.
nes_lim <- max(abs(c(fgsea_wide[[nes_x]], fgsea_wide[[nes_y]]))) * 1.50

qd <- cfg$quadrant_defs
n_q1 <- sum(fgsea_kept[[nes_x]] > 0 & fgsea_kept[[nes_y]] > 0)
n_q2 <- sum(fgsea_kept[[nes_x]] < 0 & fgsea_kept[[nes_y]] > 0)
n_q3 <- sum(fgsea_kept[[nes_x]] < 0 & fgsea_kept[[nes_y]] < 0)
n_q4 <- sum(fgsea_kept[[nes_x]] > 0 & fgsea_kept[[nes_y]] < 0)
n_metric <- qd$metric_count_fn(n_q1, n_q2, n_q3, n_q4)
n_total_sig <- nrow(fgsea_sig)
n_representative <- nrow(fgsea_kept)
metric_frac <- if (n_representative > 0) n_metric / n_representative else 0

message(sprintf(
  "  NES Spearman (all): rho = %.3f [%.3f, %.3f]",
  nes_cor_all$estimate, nes_ci_all[1], nes_ci_all[2]
))

txt_pw <- scale_text(BASE_PATHWAY, PG_W) # no multiplier — prevents label overlap in composite
txt_quad <- scale_text(BASE_QUADRANT, PG_W) * 1.15 + 1 / .pt

# 24 pathway names in one square panel only fit in a condensed face; fall back
# to the figure font where it is not installed, and say so rather than letting
# the layout silently change.
label_family <- if ("Arial Narrow" %in% systemfonts::system_fonts()$family) {
  "Arial Narrow"
} else {
  warning(
    "Arial Narrow not installed; pathway labels will crowd",
    call. = FALSE
  )
  "Helvetica"
}

label_pw <- fgsea_kept |>
  mutate(
    # Names are set as haloed text rather than filled boxes: 24 rectangles
    # cover far more of the panel than 24 dots do, so the annotation ended up
    # outweighing the data it annotates. The ten most significant carry the
    # contrast colour at full weight; the rest stay legible but recede.
    best_padj = pmin(.data[[padj_x]], .data[[padj_y]], na.rm = TRUE),
    tier = if_else(rank(best_padj, ties.method = "first") <= 10, 1L, 2L),
    label_size = txt_pw * if_else(tier == 1L, 1.20, 1.00) + 0.45,
    # Every name is black; colour belongs to the dots and to the leader that
    # joins each name to its own dot. Rank shows in weight and size alone.
    seg_col = unname(cfg$sig_colors[as.character(significance)]),
    label_face = if_else(tier == 1L, "bold", "plain"),
    pathway_label = pathway_label |>
      str_replace("Amino Acid Metabolic.*", "Amino Acid Metabolism") |>
      str_replace("Muscle System.*", "Muscle System") |>
      str_replace("Ketone Metabolic.*", "Ketone Metabolism") |>
      str_replace("^Trna Metabolic.*", "tRNA Metabolism") |>
      str_replace("^Establishment Or Maintenance Of Cell Polarity$", "Cell Polarity") |>
      str_replace(
        "^Generation Of Precursor Metabolites And Energy$",
        "Precursor Metabolites & Energy"
      ) |>
      str_replace(
        "^Protein Localization To Plasma Membrane$",
        "Plasma Membrane Protein Loc."
      ) |>
      str_replace("^Extracellular Matrix Organization$", "ECM Organization") |>
      str_replace("^Epithelial Mesenchymal Transition$", "EMT") |>
      str_replace("^Microtubule-Based Movement$", "Microtubule Movement") |>
      str_replace("^Mitochondrion Organization$", "Mitochondrial Organization") |>
      str_replace("^UV Response Dn$", "UV Response (Down)") |>
      dplyr::recode(!!!cfg$display_overrides) |>
      str_wrap(width = 18)
  )

# Labels repel cannot place cleanly get an explicit nudge from the caller,
# keyed on the raw pathway ID. A key that no longer matches warns rather than
# silently doing nothing, which is the failure mode of a hard-coded nudge.
nudges <- cfg$label_nudges %||%
  tibble(pathway = character(), nudge_x = numeric(), nudge_y = numeric())
stale <- setdiff(nudges$pathway, label_pw$pathway)
if (length(stale) > 0) {
  warning(
    "label_nudges keys not among labelled pathways: ",
    paste(stale, collapse = ", "),
    call. = FALSE
  )
}
label_pw <- label_pw |>
  left_join(nudges, by = "pathway") |>
  mutate(across(c(nudge_x, nudge_y), \(v) coalesce(v, 0)))

ns_df <- fgsea_wide |> filter(significance == "NS")
sig_df <- fgsea_wide |>
  filter(significance != "NS") |>
  mutate(draw_order = factor(significance, levels = cfg$sig_draw_order)) |>
  arrange(draw_order)

rho_kept_str <- if (!is.null(nes_cor_sig)) {
  sprintf(" | \u03c1(representatives) = %.2f", nes_cor_sig$estimate)
} else {
  ""
}
subtitle_lines <- c(
  sprintf(
    "GO Slim + Hallmark | %d pathways | %d sig., %d representative%s (EnrichmentMap-collapsed)",
    nrow(fgsea_wide), n_total_sig, n_representative,
    ifelse(n_representative == 1, "", "s")
  ),
  sprintf(
    "\u03c1 = %.2f [%.2f, %.2f], %s%s",
    nes_cor_all$estimate, nes_ci_all[1], nes_ci_all[2],
    ifelse(nes_cor_all$p.value < 0.001, "p < 0.001", sprintf("p = %.3f", nes_cor_all$p.value)),
    rho_kept_str
  ),
  sprintf(
    "%.0f%% %s of representatives | %s",
    metric_frac * 100, cfg$subtitle_metric, cfg$subtitle_interpretation
  )
)
# Panel width, not char-count guessing, sets the wrap point -- these lines
# overflowed the canvas before any of today's edits touched them.
subtitle_str <- paste(
  vapply(subtitle_lines, str_wrap, character(1), width = 62),
  collapse = "\n"
)

# Panel A's quadrant tints, literals and alpha both, so that a reader crossing
# from A to D or E reads the same pale red as concordant and the same pale blue
# as discordant. E's corners mean exacerbated/reversed rather than
# concordant/discordant, and keep the same red/blue assignment.
pD <- ggplot(mapping = aes(x = .data[[nes_x]], y = .data[[nes_y]])) +
  annotate("rect",
    xmin = qd$bg_red_1[1], xmax = qd$bg_red_1[2],
    ymin = qd$bg_red_1[3], ymax = qd$bg_red_1[4],
    fill = "#FFE0E0", alpha = 0.55, color = "grey70", linewidth = 0.2
  ) +
  annotate("rect",
    xmin = qd$bg_red_2[1], xmax = qd$bg_red_2[2],
    ymin = qd$bg_red_2[3], ymax = qd$bg_red_2[4],
    fill = "#FFE0E0", alpha = 0.55, color = "grey70", linewidth = 0.2
  ) +
  annotate("rect",
    xmin = qd$bg_blue_1[1], xmax = qd$bg_blue_1[2],
    ymin = qd$bg_blue_1[3], ymax = qd$bg_blue_1[4],
    fill = "#DCEEFF", alpha = 0.55, color = "grey70", linewidth = 0.2
  ) +
  annotate("rect",
    xmin = qd$bg_blue_2[1], xmax = qd$bg_blue_2[2],
    ymin = qd$bg_blue_2[3], ymax = qd$bg_blue_2[4],
    fill = "#DCEEFF", alpha = 0.55, color = "grey70", linewidth = 0.2
  ) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_abline(
    slope = cfg$ref_slope, intercept = 0, linetype = "dashed",
    color = "black", linewidth = 0.3
  ) +
  geom_point(
    data = ns_df, aes(shape = database),
    size = 1.0, fill = "grey70", color = "grey55", alpha = 0.40, stroke = 0.2
  ) +
  geom_point(
    data = sig_df, aes(fill = significance, size = set_size, shape = database),
    color = "black", alpha = 0.95, stroke = 0.5
  ) +
  scale_fill_manual(values = cfg$sig_colors, name = "Significance") +
  scale_shape_manual(values = c("Hallmark" = 24, "GO Slim" = 21), name = "Database") +
  scale_size_continuous(
    range = c(2, 10), name = "Set size",
    breaks = c(20, 50, 100, 200)
  ) +
  # One layer, not two: ggrepel only repels within a layer, so splitting the
  # tiers would place each blind to the other.
  geom_text_repel(
    data = label_pw, aes(label = pathway_label, segment.colour = seg_col),
    color = "black", size = label_pw$label_size,
    family = label_family, fontface = label_pw$label_face,
    bg.colour = "white", bg.r = 0.12,
    lineheight = 0.82,
    max.overlaps = Inf,
    segment.size = 0.45,
    min.segment.length = 0, show.legend = FALSE,
    xlim = c(-nes_lim * 0.93, nes_lim * 0.93),
    ylim = c(-nes_lim * 0.93, nes_lim * 0.93),
    # point.size keeps names off the dots; point.padding is the gap at the
    # end of the leader, so it stays 0 and every line runs into its dot.
    box.padding = 0.85, point.padding = 0, point.size = 8,
    force = 30, force_pull = 0.5,
    nudge_x = label_pw$nudge_x, nudge_y = label_pw$nudge_y,
    max.iter = 400000, max.time = 12,
    # Repel settles into a different local minimum per seed and the spread
    # between them is wider than any parameter here; these two were picked
    # off a sweep as the cleanest for their own panel.
    seed = cfg$label_seed %||% 42
  ) +
  annotate("label",
    x = nes_lim, y = nes_lim,
    label = sprintf("%s  n = %d", qd$label_tr, n_q1),
    hjust = 1, vjust = 1, size = txt_quad, fontface = "bold",
    color = qd$color_tr, fill = alpha("white", 0.92),
    label.padding = unit(1.5, "pt"), lineheight = 0.9
  ) +
  annotate("label",
    x = -nes_lim, y = nes_lim,
    label = sprintf("%s  n = %d", qd$label_tl, n_q2),
    hjust = 0, vjust = 1, size = txt_quad, fontface = "bold",
    color = qd$color_tl, fill = alpha("white", 0.92),
    label.padding = unit(1.5, "pt"), lineheight = 0.9
  ) +
  annotate("label",
    x = -nes_lim, y = -nes_lim,
    label = sprintf("%s  n = %d", qd$label_bl, n_q3),
    hjust = 0, vjust = 0, size = txt_quad, fontface = "bold",
    color = qd$color_bl, fill = alpha("white", 0.92),
    label.padding = unit(1.5, "pt"), lineheight = 0.9
  ) +
  annotate("label",
    x = nes_lim, y = -nes_lim,
    label = sprintf("%s  n = %d", qd$label_br, n_q4),
    hjust = 1, vjust = 0, size = txt_quad, fontface = "bold",
    color = qd$color_br, fill = alpha("white", 0.92),
    label.padding = unit(1.5, "pt"), lineheight = 0.9
  ) +
  scale_x_continuous(expand = expansion(0, 0)) +
  scale_y_continuous(expand = expansion(0, 0)) +
  coord_fixed(ratio = 1, xlim = c(-nes_lim, nes_lim), ylim = c(-nes_lim, nes_lim)) +
  labs(
    title = cfg$title,
    subtitle = subtitle_str,
    x = sprintf("NES (%s)", cfg$axis_x_label),
    y = sprintf("NES (%s)", cfg$axis_y_label)
  ) +
  FIG_THEME +
  theme(
    axis.text = element_text(size = FIG_AXIS_TEXT, face = "bold", color = "grey30"),
    axis.title = element_text(size = FIG_AXIS_TEXT, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = FIG_LEGEND_TITLE, face = "bold", color = "grey25"),
    legend.text = element_text(size = FIG_LEGEND_TEXT, color = "grey20"),
    legend.key.size = unit(2 * PRINT_SCALE, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.spacing.x = unit(3 * PRINT_SCALE, "mm"),
    legend.box.margin = margin(t = -2),
    plot.title.position = "plot",
    plot.margin = margin(0, 0, 0, 0)
  ) +
  guides(
    fill = "none",
    shape = guide_legend(
      nrow = 1, order = 1,
      keyheight = unit(4 * PRINT_SCALE, "mm"),
      keywidth = unit(4 * PRINT_SCALE, "mm"),
      override.aes = list(size = 3 * PRINT_SCALE, fill = "grey50")
    ),
    size = guide_legend(
      nrow = 1, order = 2,
      keyheight = unit(4 * PRINT_SCALE, "mm"),
      keywidth = unit(4 * PRINT_SCALE, "mm")
    )
  )

# The standalone carries title, subtitle and the bottom legend the composite
# strips, and the legend alone is wider than the square panel.
ggsave(file.path(RPT_PNG, "MAIN_panel_D_nes_scatter.png"), pD,
  width = PG_W * 1.3, height = PG_W * 1.15, units = "mm", dpi = 300
)
ggsave(file.path(RPT_PDF, "MAIN_panel_D_nes_scatter.pdf"), pD,
  width = PG_W * 1.3, height = PG_W * 1.15, units = "mm", device = pdf_device
)

export_df <- fgsea_wide |>
  transmute(
    pathway, pathway_label, database,
    !!paste0("NES_", cx) := round(.data[[nes_x]], 3),
    !!paste0("NES_", cy) := round(.data[[nes_y]], 3),
    !!paste0("padj_", cx) := signif(.data[[padj_x]], 4),
    !!paste0("padj_", cy) := signif(.data[[padj_y]], 4),
    significance = as.character(significance),
    set_size
  ) |>
  arrange(significance, desc(abs(.data[[nes_x]]) + abs(.data[[nes_y]])))
write_csv(export_df, file.path(DAT, "panel_D", "nes_scatter.csv"))

pD_legend_grob <- cowplot::get_plot_component(pD, "guide-box-bottom", return_all = FALSE)
if (!is.null(pD_legend_grob)) {
  pD_legend_plot <- cowplot::ggdraw(pD_legend_grob)
  ggsave(file.path(RPT_PNG, "MAIN_panel_D_legend.png"), pD_legend_plot,
    width = 120, height = 14, units = "mm", dpi = 300
  )
}

pD_title <- cfg$title
pD_subtitle <- subtitle_str
pD_legend <- NULL
# Strip titles but KEEP legend (legend provides shape/size key for composite)
pD <- pD + labs(title = NULL, subtitle = NULL, tag = NULL)

pw_conc_frac <- metric_frac # F04
pw_rev_frac <- metric_frac # F05

cat(sprintf("%s Panel D done\n", cfg$fig_id))

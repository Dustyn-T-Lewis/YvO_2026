#!/usr/bin/env Rscript
# Figure 4A: Training Concordance -- FDR-anchored quadrant scatter with
# flanking per-quadrant ORA bars (top 3 terms each). Primary significance
# classification is FDR (R1.3); a Pi < 0.05 (Training_Young) marker rings
# points as secondary context, it never reclassifies a point. Concordance
# rho is reported on three denominators -- whole proteome, the FDR-134 set,
# the Pi-significant set -- because R2.2 asked which one the published
# 61.6% figure came from.
#
# Returns the plot with the numbers F04.R prints in its subtitle attached as
# attr(, "stats").

setwd(here::here())

pacman::p_load(
  dplyr, tidyr, tibble, stringr, readr, ggplot2, patchwork, cowplot, readxl
)

source("04_Figures/shared/style.R")
source("04_Figures/shared/print_scale_apply.R")
source("04_Figures/shared/pathway_utils.R")
pacman::p_load(fgsea, ggrepel, qvalue, purrr)

BASE <- "04_Figures/F04"
RPT_PNG <- file.path(BASE, "b_reports", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "panels")
DAT <- file.path(BASE, "c_data")
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_A"), recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

COMP_RED <- unname(DIR_COLORS["Up"])
COMP_BLUE <- unname(DIR_COLORS["Down"])
N_SHOW <- 5

# F04-specific pathway label shortenings (used in the flanking bar panels)
DISPLAY_LABELS_F04 <- c(
  "Cargo Recognition For Clathrin Mediated Endocytosis" = "Clathrin Endocytosis",
  "The Role Of Gtse1 In G2 M Progression After G2 Checkpoint" = "GTSE1 G2/M Progression",
  "Mitochondrial Complex Ucp1 In Thermogenesis" = "Mito Complex UCP1",
  "Electron Transfer In Complex I" = "e- Transfer in Complex I",
  "Formation Of The Dystrophin Glycoprotein Complex Dgc" = "Dystrophin Complex",
  "Metabolism Of Vitamins And Cofactors" = "Vitamins/Cofactors Metab.",
  "Asparagine N Linked Glycosylation" = "Asparagine N-Linked Glycosylation",
  "Regulation Of Expression Of Slits And Robos" = "Slit/Robo Expression Regulation",
  "Ribosome Associated Quality Control" = "Ribosome QC",
  "Cellular Response To Starvation" = "Starvation Response",
  "ATP Synthesis Coupled Electron Transport" = "ATP Synthesis e- Transport",
  "Respiratory Electron Transport" = "Respiratory ETC",
  "Non Integrin Membrane Ecm Interactions" = "Non-Integrin ECM",
  "Polyol Metabolic Process" = "Polyol Metabolism",
  "Rac1 Gtpase Cycle" = "RAC1 GTPase Cycle",
  "Rac3 Gtpase Cycle" = "RAC3 GTPase Cycle"
)

# Data
dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv",
  show_col_types = FALSE
)
imp_path <- "02_imputation/c_data/02_mar_mnar_classification.csv"
imputation_df <- if (file.exists(imp_path)) {
  read_csv(imp_path, show_col_types = FALSE) |>
    transmute(gene, imputed = classification != "Complete")
} else {
  message("  Note: imputation classification not found, all marked non-imputed")
  tibble(gene = character(), imputed = logical())
}

scatter_df <- dep_df |>
  transmute(gene,
    logFC_TY = logFC_Training_Young, logFC_TO = logFC_Training_Old,
    adj_TY = adj.P.Val_Training_Young, adj_TO = adj.P.Val_Training_Old,
    adj_Int = adj.P.Val_Interaction,
    pi_TY = pi_score_Training_Young, pi_TO = pi_score_Training_Old,
    pi_Int = pi_score_Interaction
  ) |>
  filter(!is.na(logFC_TY), !is.na(logFC_TO)) |>
  left_join(imputation_df, by = "gene") |>
  mutate(
    imputed = replace_na(imputed, FALSE),
    # Scored on FDR or Pi, not FDR alone. Training(Old) returns no FDR hit at
    # all, so an FDR-only rule leaves "Sig Both" and "Sig Old only" empty by
    # construction and their swatches label nothing. Pi populates both.
    sig_y = replace_na(adj_TY, 1) < 0.05 | replace_na(pi_TY, 1) < 0.05,
    sig_o = replace_na(adj_TO, 1) < 0.05 | replace_na(pi_TO, 1) < 0.05,
    sig_i = replace_na(adj_Int, 1) < 0.05 | replace_na(pi_Int, 1) < 0.05,
    sig_class = factor(
      case_when(
        sig_i ~ "Interaction",
        sig_y & sig_o ~ "Sig Both",
        sig_y ~ "Sig Young only",
        sig_o ~ "Sig Old only",
        .default = "NS"
      ),
      levels = names(SIG_COLORS_F2)
    ),
    is_sig = sig_class != "NS",
    # Which rule found it, carried on the point shape rather than the ring the
    # panel used to overlay: the ring said Pi and nothing else, and a shape
    # says FDR, Pi or both in the same mark.
    sig_criterion = factor(
      case_when(
        !is_sig ~ NA_character_,
        (replace_na(adj_TY, 1) < 0.05 | replace_na(adj_TO, 1) < 0.05 |
          replace_na(adj_Int, 1) < 0.05) &
          (replace_na(pi_TY, 1) < 0.05 | replace_na(pi_TO, 1) < 0.05 |
            replace_na(pi_Int, 1) < 0.05) ~ "FDR + \u03a0",
        replace_na(adj_TY, 1) < 0.05 | replace_na(adj_TO, 1) < 0.05 |
          replace_na(adj_Int, 1) < 0.05 ~ "FDR",
        .default = "\u03a0"
      ),
      levels = c("FDR", "\u03a0", "FDR + \u03a0")
    ),
    pi_sig = replace_na(pi_TY, 1) < 0.05,
    quadrant = case_when(
      logFC_TY > 0 & logFC_TO > 0 ~ "Concordant Up",
      logFC_TY < 0 & logFC_TO < 0 ~ "Concordant Down",
      logFC_TY > 0 & logFC_TO < 0 ~ "Discordant (Y Up / O Down)",
      TRUE ~ "Discordant (Y Down / O Up)"
    )
  )

universe <- scatter_df$gene
message(sprintf(
  "  Total proteins: %d | Significant: %d",
  nrow(scatter_df), sum(scatter_df$is_sig)
))

pw_collection <- build_pathway_collection(
  min_size = 15, max_size = 500,
  include_goslim = FALSE,
  exclude_variants = TRUE
)

run_set_ora <- function(genes, set_name) {
  if (length(genes) < 5) {
    return(tibble())
  }
  res <- tryCatch(
    run_ora_deduplicated(
      genes = genes, universe = universe,
      pathways = pw_collection, em_cutoff = 0.5,
      min_size = 15, max_size = 500, padj_cutoff = 1
    ),
    error = function(e) {
      message("  ORA error: ", e$message)
      tibble()
    }
  )
  if (nrow(res) == 0) {
    return(tibble())
  }
  res |>
    filter(is.na(dedup_status) | dedup_status == "kept") |>
    mutate(
      set = set_name, pathway_label = clean_pathway_name(pathway),
      neg_log10_padj = -log10(padj), significant = padj < 0.05
    ) |>
    arrange(desc(neg_log10_padj)) |>
    slice_head(n = N_SHOW)
}

message("\n--- Quadrant ORA (threshold-free) ---")
ora_q1 <- run_set_ora(
  scatter_df$gene[scatter_df$quadrant == "Concordant Up"],
  "Concordant Up"
)
ora_q2 <- run_set_ora(
  scatter_df$gene[scatter_df$quadrant == "Discordant (Y Down / O Up)"],
  "Discordant (Y Down / O Up)"
)
ora_q3 <- run_set_ora(
  scatter_df$gene[scatter_df$quadrant == "Concordant Down"],
  "Concordant Down"
)
ora_q4 <- run_set_ora(
  scatter_df$gene[scatter_df$quadrant == "Discordant (Y Up / O Down)"],
  "Discordant (Y Up / O Down)"
)

all_quad_ora <- bind_rows(ora_q1, ora_q2, ora_q3, ora_q4)
if (nrow(all_quad_ora) > 0) {
  write_csv(all_quad_ora, file.path(DAT, "panel_A", "ora_quadrant.csv"))
}

xlim_range <- c(-3.1, 3.1)
ylim_range <- c(-2.8, 2.8)

# Concordance depends on how proteins are selected, so report the same contrast
# under all three rules. The whole-proteome figure is diluted by proteins with no
# resolvable change in either group; FDR selects on evidence and replicates
# better than the fold-change-weighted Pi-score. R2.6 asked which set the
# published 61.6% came from, so all three belong on the panel that raises it.
rho_on <- function(keep) {
  d <- dep_df |>
    filter(keep, !is.na(logFC_Training_Young), !is.na(logFC_Training_Old))
  list(
    n = nrow(d),
    rho = cor(d$logFC_Training_Young, d$logFC_Training_Old,
      method = "spearman"
    ),
    agree = mean(sign(d$logFC_Training_Young) == sign(d$logFC_Training_Old))
  )
}
conc_fdr <- rho_on(dep_df$adj.P.Val_Training_Young < 0.05)
conc_pi <- rho_on(dep_df$pi_score_Training_Young < 0.05)

# Fillable shapes only, so one scale can carry the criterion while the fill
# carries the contrast.
SIG_SHAPES_F2 <- c("FDR" = 21, "\u03a0" = 24, "FDR + \u03a0" = 22)

ns_df <- filter(scatter_df, sig_class == "NS")
sig_df <- filter(scatter_df, sig_class != "NS")

q_df <- scatter_df |>
  mutate(q = case_when(
    logFC_TY > 0 & logFC_TO > 0 ~ "Q1",
    logFC_TY < 0 & logFC_TO < 0 ~ "Q3",
    logFC_TY > 0 & logFC_TO < 0 ~ "Q4",
    TRUE ~ "Q2"
  ))
q_counts <- q_df |>
  count(q) |>
  deframe()
q_sig <- q_df |>
  filter(sig_class != "NS") |>
  count(q) |>
  deframe()
for (qq in c("Q1", "Q2", "Q3", "Q4")) if (is.na(q_sig[qq])) q_sig[qq] <- 0

label_df <- sig_df |>
  group_by(sig_class) |>
  arrange(desc(abs(logFC_TY) + abs(logFC_TO))) |>
  slice_head(n = 5) |>
  ungroup() |>
  mutate(
    label_fill = SIG_LABEL_FILL_F2[as.character(sig_class)],
    label_text_col = SIG_LABEL_TEXT_F2[as.character(sig_class)]
  )

txt_gene <- scale_text(BASE_GENE, 190) * 0.70 + 1 # +1pt for print legibility
txt_quad <- scale_text(BASE_QUADRANT, 190) * 0.88

# Half-bar builder: quadrant bars flank the scatter left/right, top/bottom
# halves matching each quadrant's screen corner. Left-side bars mirror
# (scale_x_reverse) so they grow outward from the scatter's edge, same as
# the original design. Every bar is solid; significance is the star alone.
# Top margin on the two upper bar plots.
UPPER_DROP_MM <- 1

# Negative top margin pulling the key up under the ORA axis values.
KEY_LIFT_MM <- -46

make_half_bars <- function(df, fill_color, side, ylim, display_labels = character(0)) {
  bar_h <- 0.42
  n_bars <- if (is.null(df) || nrow(df) == 0) 0L else min(nrow(df), 5L)

  if (n_bars == 0) {
    return(ggplot() +
      theme_void() +
      scale_y_continuous(limits = ylim, expand = c(0, 0)))
  }

  y_pos <- if (ylim[1] >= 0) {
    rev(seq(0.3, 2.3, length.out = 5))[seq_len(n_bars)]
  } else {
    seq(-0.3, -2.5, length.out = 5)[seq_len(n_bars)]
  }

  # A bar holds about FULL_BAR_CHARS characters at full length, so each label
  # wraps to its own bar's capacity and only shrinks when even two lines will
  # not fit. Labels never leave the bar; the star sits past the tip.
  FULL_BAR_CHARS <- 22
  wrap_two_lines <- function(nm, w) {
    out <- str_wrap(nm, width = w)
    while (str_count(out, "\n") > 1) {
      w <- w + 1
      out <- str_wrap(nm, width = w)
    }
    out
  }
  bars <- df |>
    arrange(desc(neg_log10_padj)) |>
    slice_head(n = 5) |>
    mutate(
      y = y_pos,
      display_name = str_remove(pathway_label, "^Reference "),
      display_name = coalesce(display_labels[display_name], display_name),
      capacity = pmax(8, floor(neg_log10_padj / max(neg_log10_padj) * FULL_BAR_CHARS)),
      display_name = map2_chr(display_name, capacity, wrap_two_lines),
      longest_line = map_int(
        strsplit(display_name, "\n", fixed = TRUE), \(ln) max(nchar(ln))
      ),
      text_size = scale_text(BASE_PATHWAY, 190) * 0.80 * pmin(1, capacity / longest_line),
      star = sig_stars(padj)
    )

  x_max <- max(bars$neg_log10_padj)
  x_display_max <- x_max * 1.3
  is_upper <- ylim[1] >= 0
  # The lower-right block drops its innermost tick so its axis title can sit
  # in that slot. Everywhere else the tick stays: only that block runs its
  # ticks from 18% of the width out to the edge, leaving the title nowhere.
  brk_fn <- function(limits) {
    b <- scales::pretty_breaks(n = 3)(limits)
    b <- b[b != 0]
    if (!is_upper && side == "right" && length(b)) b[-which.min(b)] else b
  }
  # hjust is applied in device space, which scale_x_reverse does not mirror, so
  # a left-side star left-aligned at its anchor grows back over the bar it
  # annotates. Anchoring each side by its outer edge instead makes one offset
  # clear the bar tip by the same gap whether the star is one asterisk or three.
  star_x <- bars$neg_log10_padj + x_max * 0.035
  star_hjust <- if (side == "left") 1 else 0

  p <- ggplot(bars, aes(y = y)) +
    geom_rect(
      aes(xmin = 0, xmax = neg_log10_padj, ymin = y - bar_h / 2, ymax = y + bar_h / 2),
      fill = scales::alpha(fill_color, 0.85), color = "black", linewidth = 0.3
    ) +
    geom_text(aes(x = neg_log10_padj / 2, y = y, label = display_name),
      hjust = 0.5, size = bars$text_size,
      fontface = "bold", color = "white", lineheight = 0.85
    ) +
    geom_text(aes(x = star_x, y = y, label = star),
      hjust = star_hjust, vjust = 0.5, size = 2.55 * PRINT_SCALE,
      fontface = "bold", color = "black"
    ) +
    labs(x = if (!is_upper) expression(-log[10] ~ q) else NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = FIG_AXIS_TEXT - 1, face = "bold", margin = margin(t = 0, unit = "mm")),
      # Inline with the values rather than on its own line below them: the
      # title sits at the scatter-facing end of the number row, which frees
      # that line for the key.
      axis.title.x = if (!is_upper) {
        element_text(
          size = FIG_AXIS_TEXT - 2.5, face = "bold",
          # Scatter-facing end of each number row. The right block has only
          # about 70 px before its first tick, so its title is pulled left of
          # the block edge into the gap under the scatter.
          hjust = if (side == "left") 1 else 0,
          margin = if (side == "left") {
            margin(t = -4.5, r = 0, unit = "mm")
          } else {
            margin(t = -4.5, l = 0, unit = "mm")
          }
        )
      } else {
        element_blank()
      },
      axis.line.x = element_line(color = "grey50", linewidth = 0.3),
      axis.ticks.x = element_line(color = "grey50", linewidth = 0.3),
      # patchwork aligns every panel in the row, and the scatter spans rows 1-2
      # with a zero margin, so this top margin is what holds the whole A block
      # down away from its subtitle -- not just the two upper bar plots.
      plot.margin = if (is_upper && side == "left") {
        margin(UPPER_DROP_MM, 0, 0, 3, "mm")
      } else if (is_upper) {
        margin(UPPER_DROP_MM, 3, 0, 0, "mm")
      } else if (side == "left") {
        margin(2, 0, 0, 3, "mm")
      } else {
        margin(2, 3, 0, 0, "mm")
      }
    )

  if (side == "left") {
    p + scale_x_reverse(
      limits = c(x_display_max, 0), breaks = brk_fn, expand = expansion(mult = c(0, 0))
    ) +
      scale_y_continuous(limits = ylim, expand = c(0, 0)) +
      coord_cartesian(clip = "off")
  } else {
    p + scale_x_continuous(
      limits = c(0, x_display_max), breaks = brk_fn, expand = expansion(mult = c(0, 0))
    ) +
      scale_y_continuous(limits = ylim, expand = c(0, 0)) +
      coord_cartesian(clip = "off")
  }
}

p_ul <- make_half_bars(ora_q2, COMP_BLUE, "left", c(0, 2.8),
  display_labels = DISPLAY_LABELS_F04
)
p_ll <- make_half_bars(ora_q3, COMP_RED, "left", c(-2.8, 0),
  display_labels = DISPLAY_LABELS_F04
)
p_ur <- make_half_bars(ora_q1, COMP_RED, "right", c(0, 2.8),
  display_labels = DISPLAY_LABELS_F04
)
p_lr <- make_half_bars(ora_q4, COMP_BLUE, "right", c(-2.8, 0),
  display_labels = DISPLAY_LABELS_F04
)

# Center-axis tick labels (behind points, at x=0 / y=0)
x_breaks <- seq(-3, 3, 1)
y_breaks <- seq(-2, 2, 1)
x_tick_df <- tibble(
  x = x_breaks[x_breaks != 0], y = 0,
  label = as.character(x_breaks[x_breaks != 0])
)
y_tick_df <- tibble(
  x = 0, y = y_breaks[y_breaks != 0],
  label = as.character(y_breaks[y_breaks != 0])
)

p_scatter <- ggplot(mapping = aes(x = logFC_TY, y = logFC_TO)) +
  annotate("rect",
    xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
    fill = "#FFE0E0", alpha = 0.55, color = "grey70", linewidth = 0.2
  ) +
  annotate("rect",
    xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
    fill = "#FFE0E0", alpha = 0.55, color = "grey70", linewidth = 0.2
  ) +
  annotate("rect",
    xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0,
    fill = "#DCEEFF", alpha = 0.55, color = "grey70", linewidth = 0.2
  ) +
  annotate("rect",
    xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
    fill = "#DCEEFF", alpha = 0.55, color = "grey70", linewidth = 0.2
  ) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey50", linewidth = 0.3) +
  geom_abline(
    slope = 1, intercept = 0, linetype = "dashed",
    color = "black", linewidth = 0.3
  ) +
  geom_text(
    data = x_tick_df, aes(x = x, y = y, label = label),
    vjust = 1.5, size = 1.3 * PRINT_SCALE, color = "grey40", fontface = "bold"
  ) +
  geom_text(
    data = y_tick_df, aes(x = x, y = y, label = label),
    hjust = -0.5, size = 1.3 * PRINT_SCALE, color = "grey40", fontface = "bold"
  ) +
  geom_point(
    data = ns_df, color = "grey80", fill = "grey85", shape = 21,
    size = 0.35, alpha = 0.3, stroke = 0.10
  ) +
  geom_point(
    # Fill and shape carry the whole encoding, so the marks take no outline and
    # one alpha. The border used to flag imputed proteins; that flag now lives
    # only in the supplementary workbook.
    data = sig_df, aes(fill = sig_class, shape = sig_criterion),
    # "transparent", not NA: an NA colour makes ggplot drop the rows.
    size = 2, colour = "transparent", alpha = 0.85
  ) +
  scale_fill_manual(values = SIG_COLORS_F2, name = "Significance") +
  scale_shape_manual(values = SIG_SHAPES_F2, name = "Criterion") +
  # Gene labels retired: the panel's claim is distributional -- corner counts,
  # concordance, quadrant ORA -- and named extremes support none of it. The
  # flanking bars carry the biology at pathway level, where it has power.
  # Quadrant labels: title over counts, aligned to the corner
  annotate("label",
    x = xlim_range[2], y = ylim_range[2],
    label = sprintf("Concordant Up\n%s/%s", q_sig["Q1"], q_counts["Q1"]),
    hjust = 1, vjust = 1, size = txt_quad, fontface = "bold",
    color = COMP_RED, fill = alpha("white", 0.92),
    label.padding = unit(2.5, "pt"), lineheight = 0.9
  ) +
  annotate("label",
    x = xlim_range[1], y = ylim_range[1],
    label = sprintf("%s/%s\nConcordant Down", q_sig["Q3"], q_counts["Q3"]),
    hjust = 0, vjust = 0, size = txt_quad, fontface = "bold",
    color = COMP_RED, fill = alpha("white", 0.92),
    label.padding = unit(2.5, "pt"), lineheight = 0.9
  ) +
  annotate("label",
    x = xlim_range[1], y = ylim_range[2],
    label = sprintf("Discordant\n(Y down, O up)  %s/%s", q_sig["Q2"], q_counts["Q2"]),
    hjust = 0, vjust = 1, size = txt_quad, fontface = "bold",
    color = COMP_BLUE, fill = alpha("white", 0.92),
    label.padding = unit(2.5, "pt"), lineheight = 0.9
  ) +
  annotate("label",
    x = xlim_range[2], y = ylim_range[1],
    label = sprintf("%s/%s  (Y up, O down)\nDiscordant", q_sig["Q4"], q_counts["Q4"]),
    hjust = 1, vjust = 0, size = txt_quad, fontface = "bold",
    color = COMP_BLUE, fill = alpha("white", 0.92),
    label.padding = unit(2.5, "pt"), lineheight = 0.9
  ) +
  annotate("text",
    x = xlim_range[2] - 0.05, y = -0.22,
    label = expression(log[2] * FC ~ "(Training Young)"),
    hjust = 1, vjust = 1, size = 1.3 * PRINT_SCALE, color = "grey30", fontface = "bold"
  ) +
  annotate("text",
    # Hangs from the top of the panel, which is where the top-left quadrant
    # label now sits on one row. Dropped clear of that label band.
    x = 0, y = ylim_range[2] - 0.35,
    label = expression(log[2] * FC ~ "(Training Old)"),
    hjust = 1, vjust = -0.4, size = 1.3 * PRINT_SCALE, color = "grey30", fontface = "bold",
    angle = 90
  ) +
  coord_cartesian(xlim = xlim_range, ylim = ylim_range, expand = FALSE) +
  labs(x = NULL, y = NULL) +
  FIG_THEME +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(0, 0, 0, 0, "mm"),
    legend.position = "none"
  )

# Custom Significance key
# Two channels, so two runs of glyphs: contrast by fill, then criterion by
# shape drawn in a neutral grey so the shape is the only thing that varies.
key_lvls <- c("Sig Both", "Interaction", "Sig Young only", "Sig Old only")
key_display <- c("Sig Both", "Interaction", "Sig Young", "Sig Old")
key_df <- tibble(
  category = factor(key_lvls, levels = key_lvls),
  display  = key_display,
  fill_col = unname(SIG_COLORS_F2[key_lvls]),
  x        = c(0.22, 1.02, 1.86, 2.60),
  y        = 0
)
shape_df <- tibble(
  display = names(SIG_SHAPES_F2),
  shp     = unname(SIG_SHAPES_F2),
  x       = c(3.42, 3.94, 4.42),
  y       = 0
)
p_key <- ggplot(key_df, aes(x = x, y = y)) +
  geom_point(aes(fill = category),
    shape = 21, size = 1.7 * PRINT_SCALE,
    color = "grey50", stroke = 0.6, alpha = 0.85,
    show.legend = FALSE
  ) +
  geom_text(aes(label = display),
    nudge_x = 0.06, hjust = 0,
    size = 1.55 * PRINT_SCALE, fontface = "bold", color = "grey25"
  ) +
  geom_point(
    data = shape_df, aes(x = x, y = y, shape = display),
    size = 1.7 * PRINT_SCALE, fill = "grey75",
    color = "grey40", stroke = 0.6, show.legend = FALSE
  ) +
  geom_text(
    data = shape_df, aes(x = x, y = y, label = display),
    nudge_x = 0.06, hjust = 0,
    size = 1.55 * PRINT_SCALE, fontface = "bold", color = "grey25"
  ) +
  scale_shape_manual(values = SIG_SHAPES_F2) +
  scale_fill_manual(values = setNames(key_df$fill_col, key_df$category)) +
  scale_x_continuous(
    limits = c(0.12, 5.15),
    expand = c(0, 0)
  ) +
  scale_y_continuous(limits = c(-0.12, 0.12), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = margin(KEY_LIFT_MM, 0, 0, 0, "mm"))

n_total <- nrow(scatter_df)
n_sig <- sum(scatter_df$is_sig)
n_enrich <- if (nrow(all_quad_ora) > 0) sum(all_quad_ora$significant) else 0L
r_spear <- cor(scatter_df$logFC_TY, scatter_df$logFC_TO,
  use = "complete.obs",
  method = "spearman"
)


# Storey's pi1 on the older-adult p-values of the younger-adult responders
# (Storey & Tibshirani 2003, PNAS 100:9440). A count of older-adult DEPs
# understates the shared response badly here: this tree returns zero at
# FDR 0.05, which reads as no response at all. pi1 estimates how many of these
# proteins are genuinely non-null in the older contrast without needing any of
# them to clear a threshold.
p_old_in_young_sig <- dep_df$P.Value_Training_Old[
  dep_df$adj.P.Val_Training_Young < 0.05
]
p_old_in_young_sig <- p_old_in_young_sig[!is.na(p_old_in_young_sig)]
pi1_replication <- 1 - qvalue::pi0est(
  p_old_in_young_sig,
  lambda = seq(0.05, 0.5, 0.05)
)$pi0

# Left/right flanking-bar layout, same as the original design: scatter
# centered, ORA bars flank it left (top=Discordant Y↓O↑, bottom=Concordant
# Down) and right (top=Concordant Up, bottom=Discordant Y↑O↓), key strip
# spans the full width below.
design <- c(
  area(1, 1), # p_ul (top-left ORA bars)
  area(1, 2, 2, 2), # p_scatter (rows 1-2, center)
  area(1, 3), # p_ur (top-right ORA bars)
  area(2, 1), # p_ll (bottom-left ORA bars)
  area(2, 3), # p_lr (bottom-right ORA bars)
  area(3, 1, 3, 3) # key spans full width below scatter, centered
)

composite <- p_ul + p_scatter + p_ur + p_ll + p_lr + p_key +
  plot_layout(
    design = design,
    widths = c(62, 116, 62) / 240,
    # Not an even split. The scatter spans both rows and carries its own
    # margins, so its y = 0 falls below the row boundary; measured on the
    # 460 mm canvas. Row 1 takes those units from row 2 to bring the boundary
    # -- and with it the upper bars' axis -- down onto it. The key row grew to
    # carry two stacked lines, and rows 1 and 2 give up units in proportion so
    # the alignment survives.
    heights = c(89.5, 77.5, 11) / 178
  ) +
  plot_annotation(
    title = "Training Concordance: Quadrant ORA",
    subtitle = sprintf(
      paste0(
        "Threshold-free ORA (hypergeometric) | N = %d | %d DEPs (FDR < 0.05) | %d enriched (FDR < 0.05) | ○ = Π < 0.05 (Training Young)\n",
        "Concordance ρ: %.2f all %d · %.2f the %d FDR-significant ",
        "· %.2f the %d Π-significant (Training-in-Young) | ",
        "pi1 = %.2f of them non-null in Training (Old)"
      ),
      n_total, n_sig, n_enrich,
      r_spear, n_total,
      conc_fdr$rho, conc_fdr$n,
      conc_pi$rho, conc_pi$n,
      pi1_replication
    ),
    theme = theme(
      plot.title = element_text(size = FIG_TITLE_SIZE, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE, hjust = 0, color = "grey30"),
      plot.title.position = "panel"
    )
  )

COMP_W <- 200
COMP_H <- 130
ggsave(file.path(RPT_PNG, "A_quadrant_ora.png"), composite,
  width = COMP_W, height = COMP_H, units = "mm", dpi = 300
)
ggsave(file.path(RPT_PDF, "A_quadrant_ora.pdf"), composite,
  width = COMP_W, height = COMP_H, units = "mm", device = pdf_device
)

message("\nF04 Panel A composite done")

attr(composite, "stats") <- list(
  n_total = n_total, r_spear = r_spear, conc_fdr = conc_fdr, conc_pi = conc_pi
)
invisible(composite)

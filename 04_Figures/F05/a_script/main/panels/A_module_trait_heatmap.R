#!/usr/bin/env Rscript
# Figure 5A: module-trait heatmap beside the module size bars.

setwd(here::here())

source("04_Figures/shared/style.R")

pacman::p_load(readr, dplyr, tidyr, tibble, stringr, patchwork, ggnewscale)

BASE <- "04_Figures/F05"

RPT_PNG <- file.path(BASE, "b_reports", "main", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "main", "panels")
DAT <- file.path(BASE, "c_data")
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

lmm_audit <- read_csv(file.path(DAT, "wgcna/wgcna_lmm_contrast_check.csv"),
  show_col_types = FALSE
)
lmm_wide <- function(value_col) {
  lmm_audit |>
    dplyr::select(module, contrast, all_of(value_col)) |>
    pivot_wider(names_from = contrast, values_from = all_of(value_col)) |>
    column_to_rownames("module") |>
    as.matrix()
}
lmm_r <- lmm_wide("r_equiv")
lmm_p <- lmm_wide("p_bh")
lmm_praw <- lmm_wide("p_raw")

strat_audit <- read_csv(file.path(DAT, "wgcna/wgcna_lmm_stratified_check.csv"),
  show_col_types = FALSE
)

read_matrix <- function(rel) {
  read_csv(file.path(DAT, rel), show_col_types = FALSE) |>
    column_to_rownames("module") |>
    as.matrix()
}
bl_cor_young_full <- read_matrix("wgcna/wgcna_baseline_trait_correlations_young.csv")
bl_pval_young_full <- read_matrix("wgcna/wgcna_baseline_trait_pvalues_bh_young.csv")
bl_cor_old_full <- read_matrix("wgcna/wgcna_baseline_trait_correlations_old.csv")
bl_pval_old_full <- read_matrix("wgcna/wgcna_baseline_trait_pvalues_bh_old.csv")
ch_cor_young_full <- read_matrix("wgcna/wgcna_change_trait_correlations_young.csv")
ch_pval_young_full <- read_matrix("wgcna/wgcna_change_trait_pvalues_bh_young.csv")
ch_cor_old_full <- read_matrix("wgcna/wgcna_change_trait_correlations_old.csv")
ch_pval_old_full <- read_matrix("wgcna/wgcna_change_trait_pvalues_bh_old.csv")

# Uncorrected p alongside the corrected one, for every block. Five cells survive
# BH and ten more are nominal; the panel used to render those ten as
# unremarkable, which is where the Ox. Phos. | ETC module's two near-misses
# (Training-in-Young and delta-VL, both BH = 0.058) were being lost.
bl_praw_young_full <- read_matrix("wgcna/wgcna_baseline_trait_pvalues_raw_young.csv")
bl_praw_old_full <- read_matrix("wgcna/wgcna_baseline_trait_pvalues_raw_old.csv")
ch_praw_young_full <- read_matrix("wgcna/wgcna_change_trait_pvalues_raw_young.csv")
ch_praw_old_full <- read_matrix("wgcna/wgcna_change_trait_pvalues_raw_old.csv")

module_df <- read_csv(file.path(DAT, "wgcna/wgcna_module_assignments.csv"))
MEs <- readRDS(file.path(DAT, "MEs.rds"))
meta <- read_csv(file.path(DAT, "meta.csv"))
mod_bio_labels <- read_csv(file.path(DAT, "mod_bio_labels.csv"))

non_grey <- rownames(lmm_r)[rownames(lmm_r) != "MEgrey"]

mod_size <- module_df |>
  filter(module_color != "grey") |>
  count(module_color, name = "n_proteins") |>
  mutate(module = paste0("ME", module_color)) |>
  filter(module %in% non_grey)
mod_order <- mod_size$module[order(mod_size$n_proteins)]

# Labels sit inside the count bars, whose width is the protein count, so the
# smallest modules have the least room. Nothing is shortened here any more:
# mod_bio_labels.csv carries one label per module, already written to fit, and
# it is the same string the manuscript uses.
pathway_label_map <- setNames(
  mod_bio_labels$bio_label[match(gsub("^ME", "", mod_order), mod_bio_labels$module_color)],
  mod_order
)

# Every label is "<core> | <process>", so it breaks at the rule rather than
# wherever line-balancing lands. str_wrap on its own puts the rule at the head
# of the second line in the four narrowest bars -- "ATP Synth." / "| PMF" --
# where it reads as a stray character.
wrap_at_rule <- function(x, width) {
  vapply(seq_along(x), function(i) {
    if (nchar(x[i]) <= width[i]) {
      return(x[i])
    }
    halves <- str_split_fixed(x[i], " \\| ", 2)
    if (!nzchar(halves[2])) {
      return(str_wrap(x[i], width = width[i]))
    }
    paste0(halves[1], " |\n", halves[2])
  }, character(1))
}

key_mods_me <- paste0("ME", readLines(file.path(DAT, "key_modules.txt")))
key_mods_me <- key_mods_me[nzchar(trimws(key_mods_me))]

mod_counts <- module_df |>
  filter(module_color != "grey") |>
  count(module_color, name = "n_proteins") |>
  mutate(module = paste0("ME", module_color)) |>
  filter(module %in% mod_order) |>
  mutate(
    module = factor(module, levels = mod_order),
    x_sqrt = sqrt(n_proteins)
  )

mod_counts$pathway_label <- pathway_label_map[as.character(mod_counts$module)]
mod_counts$pathway_label[is.na(mod_counts$pathway_label)] <- "N/A"
# The label sits inside the count bar and the bar's length is sqrt(n_proteins),
# so the room available varies 2.2-fold across modules. Wrap to the width each
# bar actually has, or the widest labels overflow the narrowest bars.
mod_counts$pathway_label <- wrap_at_rule(
  mod_counts$pathway_label,
  pmax(10, round(mod_counts$x_sqrt * 1.25))
)
mod_counts$bar_fill <- module_fill(mod_counts$module_color)
mod_counts$bar_text_col <- on_fill(mod_counts$bar_fill)

meta <- meta |>
  mutate(
    time = if (!is.factor(time)) factor(time, levels = c("Pre", "Post")) else time,
    age = if (!is.factor(age)) factor(age, levels = c("Young", "Old")) else age
  )

compute_aging_ttest <- function(me_data, meta_df, timepoint_label) {
  sub_meta <- meta_df |> filter(time == timepoint_label)
  rows <- list()
  for (mod in non_grey) {
    me_vals <- me_data[sub_meta$sample_id, mod]
    young_vals <- me_vals[sub_meta$age == "Young"]
    old_vals <- me_vals[sub_meta$age == "Old"]
    tt <- t.test(old_vals, young_vals) # Old - Young direction
    df_val <- tt$parameter
    t_val <- tt$statistic
    r_eq <- sign(t_val) * sqrt(t_val^2 / (t_val^2 + df_val))
    rows[[mod]] <- tibble(
      module  = mod,
      r_equiv = round(as.numeric(r_eq), 4),
      p_raw   = tt$p.value
    )
  }
  bind_rows(rows)
}

aging_pre <- compute_aging_ttest(MEs, meta, "Pre")
aging_post <- compute_aging_ttest(MEs, meta, "Post")
aging_pre$p_bh <- p.adjust(aging_pre$p_raw, method = "BH")
aging_post$p_bh <- p.adjust(aging_post$p_raw, method = "BH")

message(sprintf(
  "  Aging_Pre:  %d modules, %d with BH<0.05",
  nrow(aging_pre), sum(aging_pre$p_bh < 0.05)
))
message(sprintf(
  "  Aging_Post: %d modules, %d with BH<0.05",
  nrow(aging_post), sum(aging_post$p_bh < 0.05)
))

build_count_bars <- function(txt_count, txt_cell, txt_bar_label) {
  color_label_df <- mod_counts |>
    mutate(color_name = stringr::str_to_title(module_color))
  marker_df <- mod_counts |> filter(module %in% key_mods_me)

  ggplot(mod_counts, aes(x = x_sqrt, y = module)) +
    geom_col(
      fill = mod_counts$bar_fill, color = "black",
      linewidth = 0.3, width = 0.65
    ) +
    geom_text(aes(label = pathway_label, x = x_sqrt * 0.95),
      size = txt_bar_label, fontface = "bold",
      color = mod_counts$bar_text_col, hjust = 0, lineheight = 0.85
    ) +
    geom_text(
      data = color_label_df,
      aes(label = color_name, x = x_sqrt / 2, y = module),
      nudge_y = 0.42, size = txt_count * 0.7, fontface = "bold",
      color = "black", hjust = 0.5
    ) +
    geom_point(
      data = marker_df, aes(y = module),
      x = max(mod_counts$x_sqrt) * 1.0,
      shape = 17, size = 2.5, color = "black"
    ) +
    scale_x_reverse(
      expand = c(0, 0),
      limits = c(sqrt(350), 0),
      breaks = sqrt(c(50, 100, 200, 350)),
      labels = c(50, 100, 200, 350)
    ) +
    scale_y_discrete(labels = NULL) +
    labs(y = NULL, x = "Protein Counts") +
    FIG_THEME +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = txt_cell * 1.7 + 1, face = "bold"),
      axis.ticks.x = element_line(color = "black", linewidth = 0.6),
      axis.ticks.length.x = unit(1.5, "mm"),
      # One printed point is ~3.19 pt of type here: this 400 mm panel is
      # drawn 282 mm wide on the 470 mm canvas, which is cropped to
      # 371.3 mm and printed at 165.1 mm, so 400 mm prints as 125.4 mm.
      axis.title.x = element_text(
        size = txt_cell * 2.3 - 1.19, face = "bold",
        margin = margin(t = -9, unit = "mm")
      ),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.3),
      legend.position = "none",
      plot.margin = margin(2, 0, -6, 2, "mm")
    )
}

build_brackets <- function(brackets_df, xmin_all, xmax_all, txt_brack) {
  brackets_df <- brackets_df |> mutate(mid = (start + end) / 2)

  ggplot() +
    geom_segment(
      data = brackets_df,
      aes(x = start - 0.4, xend = end + 0.4, y = 0.20, yend = 0.20),
      linewidth = 0.4, color = "grey30"
    ) +
    geom_segment(
      data = brackets_df,
      aes(x = start - 0.4, xend = start - 0.4, y = 0.20, yend = 0.10),
      linewidth = 0.4, color = "grey30"
    ) +
    geom_segment(
      data = brackets_df,
      aes(x = end + 0.4, xend = end + 0.4, y = 0.20, yend = 0.10),
      linewidth = 0.4, color = "grey30"
    ) +
    geom_text(
      data = brackets_df,
      aes(x = mid, y = 0.60, label = label),
      size = txt_brack, fontface = "bold", color = "grey25",
      lineheight = 0.85
    ) +
    scale_x_continuous(limits = c(xmin_all, xmax_all), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(2, 2, -12, 0))
}

build_heatmap <- function(heat_df, col_labels, trait_order,
                          xmin_all, xmax_all, txt_cell) {
  contrast_traits <- c(
    "Aging", "Training_Young", "Training_Old", "Interaction",
    "Strat_Trn_Y", "Strat_Trn_O", "Aging_Pre", "Aging_Post"
  )

  # Two channels, two meanings. The outline says how strong the evidence is --
  # solid for BH, dashed for nominal-only -- and weight says how large the
  # correlation is. There used to be four channels, three of which encoded the
  # same BH flag, and a fifth that sized large-but-null cells like significant
  # ones; exactly one LMM cell qualified, and it is nominally significant
  # anyway. The dagger a nominal cell might carry instead of the dashed outline
  # would not survive the composite's device, which renders it as a full stop.
  heat_df <- heat_df |>
    mutate(
      evidence = case_when(
        !is.na(pval) & pval < 0.05 ~ "bh",
        !is.na(praw) & praw < 0.05 ~ "nominal",
        TRUE ~ "none"
      ),
      face = if_else(!is.na(cor) & abs(cor) >= 0.40, "bold", "plain"),
      txt_col = if_else(evidence == "bh", "white", "black"),
      metric = ifelse(as.character(trait) %in% contrast_traits, "r_equiv", "r")
    )

  txt_val <- (txt_cell - 0.5) * 0.81

  ggplot() +
    aes(x = xpos, y = module) +
    geom_tile(
      data = heat_df |> filter(metric == "r_equiv"),
      aes(fill = cor), color = "grey80", linewidth = 0.2
    ) +
    scale_fill_gradient2(
      low = "#1B7837", mid = "white", high = "#762A83",
      midpoint = 0, limits = c(-0.8, 0.8),
      oob = scales::squish,
      name = "r-equiv. (contrasts)",
      na.value = "white",
      guide = "none"
    ) +
    new_scale_fill() +
    geom_tile(
      data = heat_df |> filter(metric == "r"),
      aes(fill = cor), color = "grey80", linewidth = 0.2
    ) +
    scale_fill_gradient2(
      low = "#4393C3", mid = "white", high = "#D6604D",
      midpoint = 0, limits = c(-0.8, 0.8),
      oob = scales::squish,
      name = "r (phenotype)",
      na.value = "white",
      guide = "none"
    ) +
    geom_tile(
      data = heat_df |> filter(evidence == "bh"),
      color = "black", linewidth = 1.0, fill = NA
    ) +
    geom_tile(
      data = heat_df |> filter(evidence == "nominal"),
      color = "black", linewidth = 0.7, linetype = "22", fill = NA
    ) +
    geom_text(
      data = heat_df, aes(label = label), size = txt_val,
      fontface = heat_df$face, color = heat_df$txt_col
    ) +
    geom_text(
      data = heat_df |> filter(stars != ""),
      aes(label = stars), size = txt_val * 1.1, fontface = "bold",
      color = (heat_df |> filter(stars != ""))$txt_col,
      nudge_y = 0.28, vjust = 0.5
    ) +
    scale_x_continuous(
      breaks = heat_df$xpos[!duplicated(heat_df$trait)],
      labels = col_labels[trait_order],
      limits = c(xmin_all, xmax_all),
      expand = c(0, 0)
    ) +
    scale_y_discrete(labels = NULL) +
    labs(x = NULL, y = NULL) +
    FIG_THEME +
    theme( # Trimmed by the same one printed point as the count-bar axis title.
      axis.text.x = element_text(
        angle = 45, hjust = 1,
        size = txt_cell * 2.3 - 1.19, face = "bold"
      ),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      legend.position = "none",
      plot.margin = margin(0, 2, 0, -3)
    )
}

build_legend_strip <- function(low, high, title, txt_cell) {
  df <- data.frame(x = c(0, 0), y = c(0, 0), v = c(-0.8, 0.8))
  ggplot(df, aes(x, y, fill = v)) +
    geom_tile(alpha = 0) +
    scale_fill_gradient2(
      low = low, mid = "white", high = high,
      midpoint = 0, limits = c(-0.8, 0.8),
      oob = scales::squish,
      name = title,
      guide = guide_colorbar(
        barwidth       = unit(52, "mm"),
        barheight      = unit(4.5, "mm"),
        title.position = "left",
        title.vjust    = 0.8
      )
    ) +
    theme_void() +
    theme(
      legend.position = "top",
      legend.title = element_text(size = txt_cell * 2.2, face = "bold.italic"),
      legend.text = element_text(size = txt_cell * 1.8),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(0, "mm"),
      plot.margin = margin(-48, 0, 0, 0, unit = "mm")
    )
}

assemble_figure <- function(p_brackets, p_counts, p_heat,
                            p_leg_contrast, p_leg_phenotype,
                            n_contrast_cols, n_cols, txt_cell) {
  heat_cols <- n_cols
  split_at <- 4 + n_contrast_cols

  design <- c(
    area(1, 5, 1, 4 + heat_cols), # brackets
    area(2, 1, 10, 4), # protein counts
    area(2, 5, 10, 4 + heat_cols), # heatmap
    area(11, 5, 11, split_at), # r-equiv legend
    area(11, split_at + 1, 11, 4 + heat_cols) # r legend + evidence key
  )
  p <- wrap_elements(p_brackets) + p_counts + p_heat +
    p_leg_contrast + p_leg_phenotype +
    plot_layout(design = design, heights = c(0.45, rep(1, 8), 0.10))

  p + plot_annotation(
    title = NULL,
    subtitle = NULL,
    theme = theme(
      # Top margin for the title F05.R draws. The left one holds the "350"
      # tick: it is the outermost break on a scale whose limit is sqrt(350) with
      # no expansion, so the label straddles the panel edge and half of it needs
      # room outside. The count bars are the leftmost element, so that room has
      # to come from the figure margin or the glyph runs off the page.
      plot.margin = margin(12, 2, 2, 14)
    )
  )
}

message("Panel A: combined module-trait heatmap...")

strat_mat_r <- matrix(NA,
  nrow = length(non_grey), ncol = 2,
  dimnames = list(non_grey, c("Strat_Trn_Y", "Strat_Trn_O"))
)
strat_mat_p <- strat_mat_r
strat_mat_praw <- strat_mat_r

for (i in seq_len(nrow(strat_audit))) {
  mod <- strat_audit$module[i]
  if (!mod %in% non_grey) next
  col <- if (strat_audit$age_group[i] == "Young") "Strat_Trn_Y" else "Strat_Trn_O"
  strat_mat_r[mod, col] <- strat_audit$r_equiv[i]
  strat_mat_p[mod, col] <- strat_audit$p_bh[i]
  strat_mat_praw[mod, col] <- strat_audit$p_raw[i]
}

# Aging by timepoint (computed above)
aging_mat_r <- matrix(NA,
  nrow = length(non_grey), ncol = 2,
  dimnames = list(non_grey, c("Aging_Pre", "Aging_Post"))
)
aging_mat_p <- aging_mat_r
aging_mat_praw <- aging_mat_r

fill_aging <- function(df, col) {
  keep <- df$module %in% non_grey
  aging_mat_r[df$module[keep], col] <<- df$r_equiv[keep]
  aging_mat_p[df$module[keep], col] <<- df$p_bh[keep]
  aging_mat_praw[df$module[keep], col] <<- df$p_raw[keep]
}
fill_aging(aging_pre, "Aging_Pre")
fill_aging(aging_post, "Aging_Post")

subset_suffix <- function(m, suffix) {
  m <- m[non_grey, , drop = FALSE]
  colnames(m) <- paste0(colnames(m), suffix)
  m
}
bl_cor_y <- subset_suffix(bl_cor_young_full, "_Y")
bl_pval_y <- subset_suffix(bl_pval_young_full, "_Y")
bl_cor_o <- subset_suffix(bl_cor_old_full, "_O")
bl_pval_o <- subset_suffix(bl_pval_old_full, "_O")
ch_cor_y <- subset_suffix(ch_cor_young_full, "_Y")
ch_pval_y <- subset_suffix(ch_pval_young_full, "_Y")
ch_cor_o <- subset_suffix(ch_cor_old_full, "_O")
ch_pval_o <- subset_suffix(ch_pval_old_full, "_O")
bl_praw_y <- subset_suffix(bl_praw_young_full, "_Y")
bl_praw_o <- subset_suffix(bl_praw_old_full, "_O")
ch_praw_y <- subset_suffix(ch_praw_young_full, "_Y")
ch_praw_o <- subset_suffix(ch_praw_old_full, "_O")

cor_mat <- cbind(
  lmm_r[non_grey, , drop = FALSE], strat_mat_r, aging_mat_r,
  bl_cor_y, bl_cor_o, ch_cor_y, ch_cor_o
)
pval_mat <- cbind(
  lmm_p[non_grey, , drop = FALSE], strat_mat_p, aging_mat_p,
  bl_pval_y, bl_pval_o, ch_pval_y, ch_pval_o
)
praw_mat <- cbind(
  lmm_praw[non_grey, , drop = FALSE], strat_mat_praw, aging_mat_praw,
  bl_praw_y, bl_praw_o, ch_praw_y, ch_praw_o
)

trait_order <- colnames(cor_mat)
n_cols <- ncol(cor_mat)

col_labels <- c(
  Aging = "Aging",
  Training_Young = "Tr. (Y)",
  Training_Old = "Tr. (O)",
  Interaction = "Age\u00d7Trn",
  Strat_Trn_Y = "Young",
  Strat_Trn_O = "Old",
  Aging_Pre = "Pre",
  Aging_Post = "Post",
  BMI_Pre_Y = "BMI",
  VL_Pre_Y = "VL",
  LBM_Pre_Y = "LBM",
  BMI_Pre_O = "BMI",
  VL_Pre_O = "VL",
  LBM_Pre_O = "LBM",
  delta_VL_Y = "\u0394VL",
  delta_LBM_Y = "\u0394LBM",
  delta_VL_O = "\u0394VL",
  delta_LBM_O = "\u0394LBM"
)

gap <- 0.2
x <- 0
xpos <- numeric(n_cols)
section_sizes <- c(4, 2, 2, 3, 3, 2, 2) # LMM, Strat, Age@T, BaseY, BaseO, dY, dO
idx <- 1
for (s in seq_along(section_sizes)) {
  if (s > 1) x <- x + gap
  for (j in seq_len(section_sizes[s])) {
    x <- x + 1
    xpos[idx] <- x
    idx <- idx + 1
  }
}
trait_xpos <- setNames(xpos, trait_order)

# Build heat_df
heat_df <- expand.grid(
  module = non_grey, trait = trait_order,
  stringsAsFactors = FALSE
) |>
  mutate(
    cor = as.vector(cor_mat[non_grey, ]),
    pval = as.vector(pval_mat[non_grey, ]),
    praw = as.vector(praw_mat[non_grey, ]),
    stars = sig_stars(pval),
    label = sprintf("%.2f", cor),
    xpos = trait_xpos[trait]
  )
heat_df$module <- factor(heat_df$module, levels = mod_order)

PA_W <- 400
PA_H <- 280

txt_cell <- scale_text(BASE_GENE, PA_W) * 0.95 * 1.25 + 0.7
txt_count <- scale_text(BASE_COUNT, PA_W) * 0.85 * 1.25 + 2.0
txt_brack <- scale_text(BASE_GENE, PA_W) * 0.95 * 1.25 * 0.85

# The module name is wrapped to its own bar by wrap_at_rule() above, which counts
# characters against a bar width and so cannot see a change in type size. It
# keeps the size that wrap was fitted to; everything else on the panel follows
# scale_text().
txt_bar_label <- (BASE_COUNT * sqrt(PA_W / PANEL_MD) * 0.85 * 1.25 + 2.0) * 0.80

xmin_all <- min(xpos) - 0.55
xmax_all <- max(xpos) + 0.55

section_breaks <- tibble(
  x = c(
    (xpos[4] + xpos[5]) / 2,
    (xpos[6] + xpos[7]) / 2,
    (xpos[8] + xpos[9]) / 2,
    (xpos[11] + xpos[12]) / 2,
    (xpos[14] + xpos[15]) / 2,
    (xpos[16] + xpos[17]) / 2
  )
)

# Every n in the bracket headers is read off meta, so a re-run on a different
# sample set cannot leave a stale count on the figure. The LMM is fitted on all
# observations, not on pairs: Y_S05 has no Pre sample and Y_S07 no Post, so the
# design is unbalanced; the header states obs and subjects, not "paired".
n_obs <- nrow(meta)
n_subj <- dplyr::n_distinct(meta$subject)

n_baseline <- meta |>
  filter(time == "Pre") |>
  count(age)
n_base_y <- n_baseline$n[n_baseline$age == "Young"]
n_base_o <- n_baseline$n[n_baseline$age == "Old"]

# The delta columns correlate change in eigengene against change in the trait,
# so a subject counts only where both timepoints were measured for that trait.
n_delta <- meta |>
  dplyr::select(subject, age, time, VL_thick_cm, DXA_LBM_kg) |>
  pivot_wider(
    names_from = time,
    values_from = c(VL_thick_cm, DXA_LBM_kg)
  ) |>
  summarise(
    n_vl = sum(!is.na(VL_thick_cm_Post - VL_thick_cm_Pre)),
    n_lbm = sum(!is.na(DXA_LBM_kg_Post - DXA_LBM_kg_Pre)),
    .by = age
  )

# One number while both traits have the same subjects; a run where they diverge
# prints the range rather than quietly picking one of the two.
delta_n_label <- function(age_group) {
  row <- n_delta[n_delta$age == age_group, ]
  if (row$n_vl == row$n_lbm) {
    sprintf("n=%d", row$n_vl)
  } else {
    sprintf("n=%d\u2013%d", min(row$n_vl, row$n_lbm), max(row$n_vl, row$n_lbm))
  }
}

brackets <- tribble(
  ~label, ~start, ~end,
  sprintf("Contrasts (LMM)\n(%d obs, %d subj.)", n_obs, n_subj),
  xpos[1], xpos[4],
  "Stratified\nTraining", xpos[5], xpos[6],
  "Aging by\nTimepoint", xpos[7], xpos[8],
  sprintf("Baseline Young\n(n=%d Pre)", n_base_y), xpos[9], xpos[11],
  sprintf("Baseline Old\n(n=%d Pre)", n_base_o), xpos[12], xpos[14],
  sprintf("\u0394 Young\n(%s)", delta_n_label("Young")), xpos[15], xpos[16],
  sprintf("\u0394 Old\n(%s)", delta_n_label("Old")), xpos[17], xpos[18]
)

p_brackets <- build_brackets(brackets, xmin_all, xmax_all, txt_brack)
p_counts <- build_count_bars(txt_count, txt_cell, txt_bar_label)
p_heat <- build_heatmap(
  heat_df, col_labels, trait_order,
  xmin_all, xmax_all, txt_cell
)
p_heat <- p_heat +
  geom_vline(
    data = section_breaks, aes(xintercept = x),
    color = "grey55", linewidth = 0.35, linetype = "solid"
  )

p_leg_contrast <- build_legend_strip(
  low = "#1B7837", high = "#762A83",
  title = "r-equiv. (contrasts)", txt_cell = txt_cell
)
p_leg_phenotype <- build_legend_strip(
  low = "#4393C3", high = "#D6604D",
  title = "r (phenotype)", txt_cell = txt_cell
)

# The title and subtitle are drawn by F05.R, over the margin reserved at
# the top of this panel; the strings this call used to pass were discarded.
fig_A <- assemble_figure(
  p_brackets, p_counts, p_heat,
  p_leg_contrast, p_leg_phenotype,
  n_contrast_cols = 8, n_cols = n_cols, txt_cell = txt_cell
)

write_csv(heat_df, file.path(DAT, "01_panel_A_heatmap_data.csv"))

ggsave(file.path(RPT_PNG, "A_module_trait_heatmap.png"), fig_A,
  width = PA_W, height = PA_H, units = "mm",
  dpi = 300, limitsize = FALSE
)
ggsave(file.path(RPT_PDF, "A_module_trait_heatmap.pdf"), fig_A,
  width = PA_W, height = PA_H, units = "mm",
  device = pdf_device, limitsize = FALSE
)

message(sprintf("  Panel A saved: %d x %d mm (%d columns)", PA_W, PA_H, n_cols))

invisible(fig_A)

#!/usr/bin/env Rscript
# Figure 5B: module NES scatters, training in the young against training in the
# old, and aging against training in the old. Also saves the dot-size legend,
# which F05.R places on its own.

setwd(here::here())

source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

pacman::p_load(readr, dplyr, tidyr, tibble, stringr, patchwork, cowplot, ggrepel)

BASE <- "04_Figures/F05"

RPT_PNG <- file.path(BASE, "b_reports", "main", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "main", "panels")
DAT <- file.path(BASE, "c_data")
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)

stopifnot(
  "DEP results missing: 03_DEP/c_data/03_combined_results.csv" =
    file.exists("03_DEP/c_data/03_combined_results.csv"),
  "WGCNA module assignments missing — run YvO_WGCNA_run.R first" =
    file.exists(file.path(DAT, "wgcna/wgcna_module_assignments.csv")),
  "Module bio-labels missing — run YvO_WGCNA_run.R first" =
    file.exists(file.path(DAT, "mod_bio_labels.csv"))
)

module_df <- read_csv(file.path(DAT, "wgcna/wgcna_module_assignments.csv"), show_col_types = FALSE)
mod_bio <- read_csv(file.path(DAT, "mod_bio_labels.csv"), show_col_types = FALSE)
combined <- read_csv("03_DEP/c_data/03_combined_results.csv")

module_df_filt <- module_df |>
  filter(module_color != "grey", !is.na(gene), gene != "")
module_sets <- split(module_df_filt, module_df_filt$module_color) |>
  lapply(function(x) x$gene)
mod_sizes <- sapply(module_sets, length)

build_ranks <- function(df, col) {
  vals <- df[[col]]
  names(vals) <- df$gene
  vals <- vals[!is.na(vals)]
  sort(vals, decreasing = TRUE)
}
ranks_TY <- build_ranks(combined, "t_Training_Young")
ranks_TO <- build_ranks(combined, "t_Training_Old")
ranks_Aging <- build_ranks(combined, "t_Aging")

run_module_fgsea <- function(ranks, module_sets) {
  set.seed(42) # fgseaMultilevel is permutation-based; seed per call for reproducible NES/p
  res <- fgsea::fgseaMultilevel(
    pathways = module_sets, stats = ranks,
    minSize = 15, maxSize = 500, nPermSimple = 10000, eps = 0
  )
  as.data.frame(res)
}
fgsea_TY <- run_module_fgsea(ranks_TY, module_sets)
fgsea_TO <- run_module_fgsea(ranks_TO, module_sets)
fgsea_Aging <- run_module_fgsea(ranks_Aging, module_sets)

merge_fgsea <- function(res, suffix) {
  res |>
    as_tibble() |>
    dplyr::select(pathway, NES, padj, size) |>
    dplyr::rename_with(~ paste0(., "_", suffix), c(NES, padj, size))
}
fgsea_wide <- merge_fgsea(fgsea_TY, "TY") |>
  left_join(merge_fgsea(fgsea_TO, "TO"), by = "pathway") |>
  left_join(merge_fgsea(fgsea_Aging, "Aging"), by = "pathway") |>
  mutate(module_color = pathway, n_proteins = mod_sizes[pathway])
fgsea_wide <- fgsea_wide |>
  left_join(mod_bio |> dplyr::select(module_color, module_id, bio_label), by = "module_color") |>
  mutate(
    bio_label = ifelse(is.na(bio_label), stringr::str_to_title(module_color), bio_label),
    module_id = ifelse(is.na(module_id), toupper(substr(module_color, 1, 1)), module_id),
    # The dots are labelled by colour, not pathway. Nine wrapped pathway
    # labels collided across this scatter once the labels grew; panel A names
    # them in full a few centimetres to the left, and F06 labels its cells the
    # same way. bio_label keeps the pathway, which is what the sheet exports.
    plot_label = stringr::str_to_title(module_color),
    sig_conc = case_when(
      !is.na(padj_TY) & padj_TY < 0.05 & !is.na(padj_TO) & padj_TO < 0.05 ~ "Both",
      !is.na(padj_TY) & padj_TY < 0.05 ~ "Young only",
      !is.na(padj_TO) & padj_TO < 0.05 ~ "Old only", TRUE ~ "NS"
    ),
    sig_rev = case_when(
      !is.na(padj_Aging) & padj_Aging < 0.05 & !is.na(padj_TO) & padj_TO < 0.05 ~ "Both",
      !is.na(padj_Aging) & padj_Aging < 0.05 ~ "Aging only",
      !is.na(padj_TO) & padj_TO < 0.05 ~ "Training only", TRUE ~ "NS"
    )
  )

# Panels A and B are authored on different canvases and land in the composite at
# different scales, so equal printed type needs unequal canvas type. Panel A is
# 400 mm wide in a 0.60 x 470 = 282 mm box and is width-limited, so it is drawn
# at 282/400 = 0.705; panel B is PB_W x PB_H in a 0.55 x 470 by
# (grid_top - grid_bot) x 300 = 258.5 x 158.46 mm box and is height-limited, so
# it is drawn at 158.46/270 = 0.5869. Both then shrink by the same 165.1/371.3
# at print, which cancels, leaving panel B's type a factor 1.2013 larger on its
# own canvas. Panel A's count-bar ticks and axis title are 13.531 and 15.764 pt,
# printing at 4.242 and 4.942 pt; the three sizes below match that. The box
# height feeds back: changing them moves panel B's own panel borders, so
# grid_top and grid_bot in F05.R have to be re-solved afterwards.
PB_AXIS_TEXT <- 16.25
PB_AXIS_TITLE <- 18.94
# The quadrant names are furniture and the module names are the data, so the
# furniture is set smaller. At the old 6.66 it outweighed every module name and
# its four boxes owned all four corners, which is where repel most needs room:
# dropping it to 5.0 took measured label collisions from 281 px to 5.
PB_QUAD_MM <- 5.0
PB_LABEL_MM <- 5.2

# F04 panels A, D and E tint their quadrants with these two at this alpha.
# Matching them puts every NES scatter in the manuscript on one background,
# and the paler ground gives the coloured module names room to read. Green
# is the third member because reversal is its own claim, not discordance.
QUAD_TINT <- c(pink = "#FFE0E0", blue = "#DCEEFF", green = "#E0FFE0")
QUAD_ALPHA <- 0.55

# WGCNA module colours name the modules everywhere -- panel A, the workbook,
# F06 -- so the dots keep them untouched. As type on white, five of the nine
# fail: yellow and green clear 1.4:1. Darkening to a fixed 4.5:1 while pushing
# saturation back up keeps each name recognisably its module's colour and
# legible at print size.
fgsea_wide <- fgsea_wide |>
  mutate(dot_hex = module_fill(module_color), text_hex = text_safe(dot_hex)) |>
  # Largest first so a 340-protein module cannot bury a 70-protein one. Black
  # M6 sat entirely inside yellow M4 in the aging scatter: 0.198 NES apart,
  # which left 0.23 mm of black showing at print.
  arrange(desc(n_proteins))

# Both scatters share the Training Old axis, so one limit across all three
# contrasts puts them on the same scale and makes the y = x line mean the same
# thing in each.
nes_lim <- max(
  abs(c(fgsea_wide$NES_TY, fgsea_wide$NES_TO, fgsea_wide$NES_Aging)),
  na.rm = TRUE
) * 1.45

build_scatter <- function(df, x_col, y_col, x_lab, y_lab, quad_labels) {
  x_vals <- df[[x_col]]
  y_vals <- df[[y_col]]

  q_tr <- sum(x_vals > 0 & y_vals > 0, na.rm = TRUE)
  q_bl <- sum(x_vals < 0 & y_vals < 0, na.rm = TRUE)
  q_tl <- sum(x_vals < 0 & y_vals > 0, na.rm = TRUE)
  q_br <- sum(x_vals > 0 & y_vals < 0, na.rm = TRUE)

  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    annotate("rect",
      xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
      fill = quad_labels$fill[1], alpha = QUAD_ALPHA
    ) +
    annotate("rect",
      xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
      fill = quad_labels$fill[2], alpha = QUAD_ALPHA
    ) +
    annotate("rect",
      xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0,
      fill = quad_labels$fill[3], alpha = QUAD_ALPHA
    ) +
    annotate("rect",
      xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
      fill = quad_labels$fill[4], alpha = QUAD_ALPHA
    ) +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
    geom_abline(
      slope = 1, intercept = 0, linetype = "dashed",
      color = "black", linewidth = 0.3
    ) +
    geom_point(aes(size = n_proteins),
      fill = df$dot_hex, color = "black",
      shape = 21, alpha = 0.85, stroke = 0.5
    ) +
    # Panel A keys each module by this same string over a bar in this same
    # colour, so the two panels agree without the reader decoding anything.
    # point.padding stays 0 so every leader runs into its own dot; point.size
    # is what keeps the names off them. The vertical bound stops labels
    # climbing into the corners the quadrant boxes occupy, and the seed was
    # taken from a sweep that scored rendered label overlap in pixels -- four
    # seeds reach zero at these settings, so this is a basin, not a fluke.
    ggrepel::geom_text_repel(
      aes(label = plot_label, segment.colour = text_hex),
      color = df$text_hex, size = PB_LABEL_MM, fontface = "bold",
      bg.colour = "white", bg.r = 0.14, lineheight = 0.82,
      segment.size = 0.7, min.segment.length = 0, max.overlaps = Inf,
      show.legend = FALSE,
      box.padding = 0.80, point.padding = 0, point.size = 9,
      force = 25, force_pull = 0.55,
      max.iter = 200000, max.time = 5, seed = 42,
      xlim = c(-nes_lim * 0.95, nes_lim * 0.95),
      ylim = c(-nes_lim * 0.88, nes_lim * 0.88)
    ) +
    annotate("label",
      x = nes_lim, y = nes_lim,
      label = sprintf("%s\nn=%d", quad_labels$label[1], q_tr),
      hjust = 1, vjust = 1, size = PB_QUAD_MM, fontface = "bold", lineheight = 0.9,
      color = quad_labels$text_col[1], fill = scales::alpha("white", 0.92),
      label.padding = unit(1.5, "pt")
    ) +
    annotate("label",
      x = -nes_lim, y = -nes_lim,
      label = sprintf("%s\nn=%d", quad_labels$label[2], q_bl),
      hjust = 0, vjust = 0, size = PB_QUAD_MM, fontface = "bold", lineheight = 0.9,
      color = quad_labels$text_col[2], fill = scales::alpha("white", 0.92),
      label.padding = unit(1.5, "pt")
    ) +
    annotate("label",
      x = nes_lim, y = -nes_lim,
      label = sprintf("%s\nn=%d", quad_labels$label[3], q_br),
      hjust = 1, vjust = 0, size = PB_QUAD_MM, fontface = "bold", lineheight = 0.9,
      color = quad_labels$text_col[3], fill = scales::alpha("white", 0.92),
      label.padding = unit(1.5, "pt")
    ) +
    annotate("label",
      x = -nes_lim, y = nes_lim,
      label = sprintf("%s\nn=%d", quad_labels$label[4], q_tl),
      hjust = 0, vjust = 1, size = PB_QUAD_MM, fontface = "bold", lineheight = 0.9,
      color = quad_labels$text_col[4], fill = scales::alpha("white", 0.92),
      label.padding = unit(1.5, "pt")
    ) +
    scale_size_continuous(
      range = c(4.5, 12), name = "Proteins",
      breaks = c(50, 100, 200, 300)
    ) +
    scale_x_continuous(
      breaks = seq(-6, 6, 3), expand = expansion(0, 0)
    ) +
    scale_y_continuous(
      breaks = seq(-6, 6, 3), expand = expansion(0, 0)
    ) +
    coord_fixed(
      ratio = 1, xlim = c(-nes_lim, nes_lim), ylim = c(-nes_lim, nes_lim)
    ) +
    labs(title = NULL, subtitle = NULL, x = x_lab, y = y_lab) +
    FIG_THEME +
    theme(
      axis.text = element_text(size = PB_AXIS_TEXT, face = "bold", color = "grey30"),
      axis.title.x = element_text(size = PB_AXIS_TITLE, face = "bold"),
      axis.title.y = element_text(
        size = PB_AXIS_TITLE, face = "bold", margin = margin(r = -4)
      ),
      legend.position = "none",
      plot.margin = margin(2, 4, 2, -8)
    )
}

quad_conc <- list(
  label = c("Concordant Up", "Concordant Down", "Discordant", "Discordant"),
  fill = QUAD_TINT[c("pink", "pink", "blue", "blue")],
  text_col = c("#D6604D", "#D6604D", "#4393C3", "#4393C3")
)

p_top <- build_scatter(
  fgsea_wide, "NES_TY", "NES_TO",
  "NES (Training Young)", "NES (Training Old)", quad_conc
)

quad_rev <- list(
  label = c("Exacerbated", "Exacerbated", "Reversed", "Reversed"),
  fill = QUAD_TINT[c("pink", "pink", "green", "green")],
  text_col = c("#D6604D", "#D6604D", "#2E7D32", "#2E7D32")
)

p_bottom <- build_scatter(
  fgsea_wide, "NES_Aging", "NES_TO",
  "NES (Aging)", "NES (Training Old)", quad_rev
)

scatters_panel <- (p_top / p_bottom) +
  plot_layout(heights = c(1, 1))

PB_W <- 220
PB_H <- 270

ggsave(file.path(RPT_PNG, "B_nes_scatters.png"), scatters_panel,
  width = PB_W, height = PB_H, units = "mm", dpi = 300
)
ggsave(file.path(RPT_PDF, "B_nes_scatters.pdf"), scatters_panel,
  width = PB_W, height = PB_H, units = "mm", device = get_pdf_device()
)

p_legend_src <- p_top +
  scale_size_continuous(
    range = c(4.5, 12), name = "Proteins",
    breaks = c(50, 100, 200, 300),
    guide = guide_legend(
      nrow = 1,
      override.aes = list(alpha = 0.7, fill = "grey60", stroke = 0)
    )
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key = element_rect(fill = NA, color = NA),
    legend.key.size = unit(5, "mm"),
    legend.background = element_rect(fill = NA, color = NA)
  )
legend_grob <- cowplot::get_plot_component(p_legend_src, "guide-box-bottom", return_all = FALSE)
p_legend <- cowplot::ggdraw(legend_grob)

ggsave(file.path(RPT_PNG, "B_nes_scatters_legend.png"), p_legend,
  width = 90, height = 16, units = "mm", dpi = 300
)

# dot_hex, text_hex and plot_label are drawing internals; the sheet keeps
# module_color, which is how every other sheet names a module.
write_csv(
  fgsea_wide |> dplyr::select(-dot_hex, -text_hex, -plot_label),
  file.path(DAT, "panel_B_module_fgsea.csv")
)

message("Panel B: scatters + legend saved")

invisible(scatters_panel)

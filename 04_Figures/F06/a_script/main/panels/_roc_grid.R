# Figure 6A and 6B: per-module ROC cells, split by what they classify, and the
# module key both panels share. Reads the module grid that
# S7_A_module_grid.R writes to c_data/module_grid/, and runs that panel first
# when the grid is missing.

pacman::p_load(tidyverse, patchwork)

source("04_Figures/shared/style.R")
source("04_Figures/shared/figure_supplement_helpers.R")
source("04_Figures/F06/a_script/_panel_selection.R", local = TRUE)

BASE <- "04_Figures/F06"
DAT <- file.path(BASE, "c_data", "module_grid")
RPT <- file.path(BASE, "b_reports", "main", "panels")
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(file.path(DAT, "module_grid_summary.csv"))) {
  source_panel("04_Figures/F06/a_script/supp/panels/S7_A_module_grid.R")
}

summ_all <- read_csv(file.path(DAT, "module_grid_summary.csv"), show_col_types = FALSE)
curves_all <- read_csv(file.path(DAT, "module_grid_curves.csv"), show_col_types = FALSE)

F05_SUPP <- "04_Figures/F05/c_data/F05_data.xlsx"
stopifnot(
  "run 04_Figures/F05/a_script/F05_data.R first: missing F05_data.xlsx" =
    file.exists(F05_SUPP)
)
mod_bio_df <- read_sheet_df(F05_SUPP, "WGCNA_mod_bio_labels")
# display_label is the canonical "Turquoise: Beta Ox. | TCA" string every
# other figure prints. Re-deriving colour + separator + bio_label here is how
# F05 and F06 drift apart on the same module.
mod_display <- setNames(mod_bio_df$display_label, mod_bio_df$module_color)

# Two panels, two questions. The grid used to rank all 48 (comparison x
# module) cells together and draw the top twelve, which put six different
# comparisons in one panel and made every cell carry a caption saying which it
# was. Split by comparison, each panel makes one statement.
cells_df <- main_panel_cells(summ_all)

age_cells <- cells_df |> filter(panel == "age")
train_cells <- cells_df |> filter(panel == "training")

GRID_NCOL <- 3L

build_cell <- function(rr, i, n) {
  dd <- curves_all |> filter(row == rr$row, module == rr$module)
  col <- module_fill(rr$module)
  show_x <- i > n - GRID_NCOL # bottom row, whatever the row count is
  show_y <- ((i - 1L) %% GRID_NCOL) == 0L # left column

  auc_line <- sprintf("AUC %.2f", rr$auc)
  ci_line <- sprintf("[%.2f, %.2f]", rr$ci_lo, rr$ci_hi)
  p_line <- if (rr$perm_p < 0.001) "p < .001" else sprintf("p = %.3f", rr$perm_p)
  q_line <- if (is.na(rr$q_bh)) "q = NA" else if (rr$q_bh < 0.001) "q < .001" else sprintf("q = %.3f", rr$q_bh)
  ctx_line <- rr$cell_label

  g <- ggplot(dd, aes(fpr, tpr)) +
    geom_ribbon(aes(ymin = 0, ymax = tpr), fill = col, alpha = 0.32) +
    geom_abline(
      slope = 1, intercept = 0, linetype = "dashed",
      color = "grey70", linewidth = 0.3
    ) +
    geom_line(color = col, linewidth = 1.1) +
    annotate("label",
      x = 0.69, y = 0.45,
      label = paste(auc_line, ci_line, p_line, q_line, sep = "\n"),
      hjust = 0.5, vjust = 0.5, size = 2.0, fontface = "bold",
      color = "grey5", fill = "white",
      linewidth = 0.2, label.padding = unit(3, "pt"),
      label.r = unit(3, "pt")
    ) +
    # Context title centered at bottom of plot (inside, with white background)
    annotate("label",
      x = 0.50, y = 0.03, label = ctx_line,
      hjust = 0.5, vjust = 0, size = 3.1, fontface = "bold.italic",
      color = "grey15", fill = scales::alpha("white", 0.82),
      linewidth = 0, label.padding = unit(1.5, "pt"),
      label.r = unit(3, "pt")
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = c(0, 1), labels = c("0", "1"),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = c(0, 1), labels = c("0", "1"),
      expand = c(0, 0)
    ) +
    coord_cartesian(clip = "off") +
    labs(
      x = if (show_x) "1 \u2212 Specificity" else NULL,
      y = if (show_y) "Sensitivity" else NULL
    ) +
    theme_classic(base_size = 15) +
    theme(
      panel.border = element_rect(
        fill = NA,
        color = rr$border_color,
        linewidth = rr$border_lw,
        linetype = rr$border_lty
      ),
      axis.title.x = element_text(
        size = 5, face = "bold",
        color = "grey25",
        margin = margin(t = -2)
      ),
      axis.title.y = element_text(
        size = 5, face = "bold",
        color = "grey25",
        margin = margin(r = -2)
      ),
      axis.text.x = if (show_x) {
        element_text(size = FIG_AXIS_TEXT, color = "grey30")
      } else {
        element_blank()
      },
      axis.text.y = if (show_y) {
        element_text(size = FIG_AXIS_TEXT, color = "grey30")
      } else {
        element_blank()
      },
      axis.ticks = element_line(color = "grey40", linewidth = 0.3),
      axis.line = element_blank(),
      plot.margin = margin(0, 0, 8, 0)
    ) # bottom margin for context title
  g
}

build_grid <- function(df) {
  n <- nrow(df)
  cells <- lapply(seq_len(n), function(i) build_cell(df[i, ], i, n))
  wrap_plots(cells, ncol = GRID_NCOL) &
    theme(plot.margin = margin(0, 0, 8, 0))
}


# One key for both panels, listing every module either of them draws, in
# canonical order. Data-driven, so it can neither show a module that is absent
# nor omit one that is present.
MODULES_KEY <- intersect(names(mod_display), unique(cells_df$module))
key_df <- tibble(
  module = MODULES_KEY,
  label = unname(mod_display[MODULES_KEY])
)

build_key_cell <- function(mod, lab) {
  ggplot() +
    annotate("rect",
      xmin = 0.02, xmax = 0.09, ymin = 0.22, ymax = 0.78,
      fill = module_fill(mod), color = "grey20", linewidth = 0.4
    ) +
    annotate("text",
      x = 0.11, y = 0.5, label = lab,
      hjust = 0, vjust = 0.5, size = 2.6, color = "grey10"
    ) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(plot.margin = margin(5, 1, 5, 1))
}

key_cells <- purrr::map2(key_df$module, key_df$label, build_key_cell)
# Two rows of three across the full canvas. Six across left each label 35 mm
# and cut "Turquoise: Beta Ox. | TCA" short; three across gives it 70 mm.
key_composite <- wrap_plots(key_cells, nrow = 2, ncol = 3)

# Both grids are half the 210 mm canvas and the same height, so the size is
# the function's, not the caller's.
save_panel <- function(plot, stem) {
  ggsave(file.path(RPT, paste0(stem, ".png")), plot,
    width = 105, height = 62, units = "mm", dpi = 300
  )
  ggsave(file.path(RPT, paste0(stem, ".pdf")), plot,
    width = 105, height = 62, units = "mm", device = get_pdf_device()
  )
}

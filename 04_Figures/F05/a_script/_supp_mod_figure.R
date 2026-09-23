# The four modules that carry an association, one triptych row each. Nine
# triptychs are rendered by _supp_mod_triptych.R and until now shipped nowhere.
#
# Triptychs and not the hub networks, because the sentence that cites this
# figure promises three things per module -- a z-score heatmap, an eigengene
# trajectory and ORA enrichment -- and those are the triptych's three panels.

pacman::p_load(cowplot, png, grid, readr)

MOD_SRC <- file.path(BASE, "b_reports", "supp", "png", "modules")
MOD_PDF <- file.path(BASE, "b_reports", "supp", "pdf")
MOD_PNG <- file.path(BASE, "b_reports", "supp", "png")

# The three age-associated modules first, then the training-associated one,
# which is the order the Discussion takes them in. Every statistic is read from
# the LMM contrast table, which is the test Figure 5 panel A draws and the test
# the Discussion quotes. Typed in by hand, the age rows carried the
# cross-sectional correlations instead (+0.70, -0.76, -0.42) while the text
# carried the LMM (+0.58, -0.61, -0.44), and the figure disagreed with the
# sentence citing it.
mod_rows <- tibble::tribble(
  ~module, ~tag, ~contrast, ~prefix,
  "turquoise", "A", "Aging", "age",
  "black", "B", "Aging", "age",
  "yellow", "C", "Aging", "age",
  "blue", "D", "Training_Young", "training in the young"
)

.lmm <- read_csv(
  file.path(BASE, "c_data", "wgcna", "wgcna_lmm_contrast_check.csv"),
  show_col_types = FALSE
)

mod_rows$assoc <- vapply(seq_len(nrow(mod_rows)), function(i) {
  hit <- .lmm[
    .lmm$module == paste0("ME", mod_rows$module[i]) &
      .lmm$contrast == mod_rows$contrast[i],
  ]
  if (nrow(hit) != 1L) {
    stop("No single LMM row for ", mod_rows$module[i], " / ",
         mod_rows$contrast[i])
  }
  sprintf(
    "%s r = %+.2f, q = %s", mod_rows$prefix[i], hit$r_equiv,
    formatC(hit$p_bh, format = "e", digits = 1)
  )
}, character(1))

# Panel files are named from bio_label, which changes whenever the labels do, so
# the names are rebuilt from the same table rather than written out.
.mod_lab <- read_csv(
  file.path(BASE, "c_data", "mod_bio_labels.csv"),
  show_col_types = FALSE
)
mod_slug <- setNames(
  gsub("[|/ ]+", "_", gsub("\\.", "", tolower(.mod_lab$bio_label))),
  .mod_lab$module_color
)

read_mod_panel <- function(kind, module) {
  path <- file.path(
    MOD_SRC, sprintf("SUPP_%s_%s_%s.png", kind, module, mod_slug[[module]])
  )
  if (!file.exists(path)) stop("Missing: ", path)
  rasterGrob(readPNG(path), interpolate = TRUE)
}

MOD_W <- 200
MOD_H <- 310
mod_txt <- composite_text_sizes(MOD_W, 178)

# A triptych is 2.55:1, so a row 0.95 of this canvas wide is 0.24 of it tall.
# The header takes the rest of the row.
ROW_H <- 0.2419
TITLE_H <- 0.0323
HEAD_H <- 0.026

modules_fig <- ggdraw() +
  theme(plot.background = element_rect(fill = "white", color = NA)) +
  draw_label(
    "WGCNA modules carrying an age or training association",
    x = 0.025, y = 1 - TITLE_H / 2, size = mod_txt$title,
    fontface = "bold", hjust = 0, vjust = 0.5
  )

for (i in seq_len(nrow(mod_rows))) {
  m <- mod_rows$module[i]
  top <- 1 - TITLE_H - (i - 1) * ROW_H
  modules_fig <- modules_fig +
    draw_label(
      mod_rows$tag[i],
      x = 0.02, y = top - HEAD_H / 2, size = mod_txt$tag,
      fontface = "bold", hjust = 0, vjust = 0.5
    ) +
    draw_label(
      # Parenthesised, not rule-separated: the label already carries a vertical
      # rule and a second one would read as a third field.
      sprintf(
        "%s  (%s)",
        .mod_lab$display_label[match(m, .mod_lab$module_color)],
        mod_rows$assoc[i]
      ),
      x = 0.055, y = top - HEAD_H / 2, size = mod_txt$subtitle,
      fontface = "bold", hjust = 0, vjust = 0.5
    ) +
    draw_grob(read_mod_panel("triptych", m),
      x = 0.025, y = top - ROW_H, width = 0.95, height = ROW_H - HEAD_H,
      hjust = 0, vjust = 0
    )
}

mod_device <- get_raster_pdf_device()

ggsave(file.path(MOD_PDF, "SUPP_F05_modules.pdf"), modules_fig,
  width = MOD_W, height = MOD_H, units = "mm",
  device = mod_device, limitsize = FALSE
)
embed_pdf_fonts(file.path(MOD_PDF, "SUPP_F05_modules.pdf"))
ggsave(file.path(MOD_PNG, "SUPP_F05_modules.png"), modules_fig,
  width = MOD_W, height = MOD_H, units = "mm",
  dpi = 300, limitsize = FALSE
)

message("F05 supp modules saved: SUPP_F05_modules.{pdf,png}")

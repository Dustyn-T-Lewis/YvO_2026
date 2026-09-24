pacman::p_load(
  dplyr, tidyr, readr, readxl, stringr, purrr, tibble,
  ggplot2, ggtext, gridtext, grid, gtable, scales, ragg, openxlsx
)

set.seed(42)

PANEL_DIR <- here::here("04_Figures", "abstract_panels")
SCRIPT_DIR <- file.path(PANEL_DIR, "a_script")
REPORT_DIR <- file.path(PANEL_DIR, "b_reports")
DATA_DIR <- file.path(PANEL_DIR, "c_data")
F01_XLSX <- here::here("04_Figures", "F01", "c_data", "F01_data.xlsx")
F06_XLSX <- here::here("04_Figures", "F06", "c_data", "F06_data.xlsx")
F04_XLSX <- here::here("04_Figures", "F04", "c_data", "F04_data.xlsx")

source(here::here("04_Figures", "shared", "devices.R"))
# MODULE_FILL and MODULE_TAG come from the shared palette rather than a copy
# here: a module is named by its colour, so two definitions is two chances for
# the same module to be two colours. The palette file carries no side effects,
# which is why it can be sourced where style.R cannot.
source(here::here("04_Figures", "shared", "module_palette.R"))

# The base pdf device drops Greek, Delta and rho -- devices.R records the
# measurement -- and these cards carry all three. get_pdf_device()'s quartz
# fallback drops its dots, so bg cannot travel through it, and the cards are
# drawn on a transparent ground. Resolved once at load rather than per save.
panel_pdf <- local({
  if (probe_device(function() cairo_pdf(tempfile()))) {
    function(path, width, height) {
      cairo_pdf(path, width = width, height = height, bg = "transparent")
    }
  } else if (probe_device(function() {
    quartz(type = "pdf", file = tempfile(fileext = ".pdf"))
  })) {
    function(path, width, height) {
      quartz(
        file = path, type = "pdf", width = width, height = height,
        bg = "transparent"
      )
    }
  } else {
    function(path, width, height) {
      pdf(path, width = width, height = height, bg = "transparent")
    }
  }
})

INK <- "black"
MUTED <- "grey30"
# Taken from shared/style.R rather than invented here. Every card compares
# training contrasts, not age groups, which is the distinction
# panel_C_trajectory.R:229 draws, so the age split takes CONTRAST_COLORS and
# not AGE_COLORS. Red therefore means Younger on cards A and B and Up on card
# C; F04 already carries both palettes in one figure, panel A on DIR_COLORS
# and panel C on CONTRAST_COLORS. style.R is read, not sourced: it executes
# devices.R and defines size globals that print_scale_apply.R mutates in
# place, and these cards must not inherit that. abstract.R asserts the two
# stay in step.
AGE_PAL <- c(Young = "#E05A4E", Old = "#5DA5DA")
# The direction pair sits at the dark end of the same RdBu ramp style.R draws
# the age pair from. Sharing the mid tones made red mean Younger on one card
# and Up three cards along, a distinction 8% of lightness could not carry.
DIR_PAL <- c(Up = "#B2182B", Down = "#2166AC")
AGREE_PAL <- c(`TRUE` = "#3D7A5A", `FALSE` = "#E0A458")
OUTLINE <- "black"
OUTLINE_W <- 0.2

TITLE_PT <- 9
GLYPH_PT <- 6.5
# One gap per column boundary. The first two are breathing room between
# panels; the last is larger because column 3 carries the widest titles in
# the grid and they overrun their own column by 0.42 in, which at 0.10 left
# the module pathway labels sitting underneath them.
COL_GAP <- c(0.25, 0.25, 0.55)
# The concise row is six cards of one plot each, so one gap serves them all,
# and its panels run taller than the two-row grid's because it has one row of
# them to spend the height on.
CONCISE_GAP <- 0.12
CONCISE_PANEL_H <- 1.10
# A canvas margin. An axis label centred on the outermost break hangs half its
# width past the panel, which the two-row grid absorbs in the padding
# fit_width adds around stacked plots. A row of single plots has no such
# padding, so the first card's "0" ran off the page without this.
ROW_PAD <- 0.06
# Clear space kept to the right of a card title, inside the card.
TITLE_PAD <- 0.12
ROW_GAP <- 0.04
PANEL_H <- 0.86
DPI <- 661

pt <- function(x) unname(x) / ggplot2::.pt

drop0 <- function(x) sub("^([+-]?)0\\.", "\\1.", x)

# The map's ramp is read out of F04's RRHO rather than copied, so the abstract
# and the main figure cannot drift onto two colour languages.
rrho_jet <- function() {
  txt <- paste(
    readLines(
      here::here(
        "04_Figures", "shared", "comparison_panels", "panel_E_rrho2.R"
      ),
      warn = FALSE
    ),
    collapse = "\n"
  )
  block <- sub(
    "(?s)^.*?JET_COLORS\\s*<-\\s*c\\((.*?)\\).*$", "\\1", txt,
    perl = TRUE
  )
  hits <- regmatches(block, gregexpr('"[^"]+"', block))[[1]]
  stopifnot(length(hits) == 9)
  gsub('"', "", hits)
}

npc_to_data <- function(npc, lim) lim[1] + npc * diff(lim)

theme_panel <- function() {
  theme_bw(base_size = GLYPH_PT, base_family = "Helvetica") +
    theme(
      plot.title = element_text(
        size = TITLE_PT, face = "bold", colour = INK, hjust = 0,
        margin = margin(b = 0.5)
      ),
      plot.subtitle = element_blank(),
      axis.title = element_text(
        size = GLYPH_PT - 0.5, face = "bold", colour = INK
      ),
      axis.text = element_text(face = "bold", size = GLYPH_PT, colour = INK),
      axis.ticks = element_line(colour = MUTED, linewidth = 0.25),
      axis.ticks.length = unit(1.2, "pt"),
      axis.line = element_line(colour = MUTED, linewidth = 0.35),
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(
        face = "bold", size = GLYPH_PT - 0.5, colour = INK
      ),
      panel.background = element_rect(fill = NA, colour = NA),
      plot.background = element_rect(fill = NA, colour = NA),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(1, 1, 1, 1)
    )
}

panel_key <- function(labels, colours, xlim, ylim, x = 0.03, top = 0.96,
                      row = 0.13, hjust = 0) {
  geom_text(
    data = tibble(
      label = labels, colour = colours,
      kx = npc_to_data(x, xlim),
      ky = npc_to_data(top - (seq_along(labels) - 1) * row, ylim)
    ),
    inherit.aes = FALSE,
    aes(kx, ky, label = label, colour = I(colour)),
    hjust = hjust, vjust = 1, size = pt(GLYPH_PT), fontface = "bold",
    show.legend = FALSE
  )
}

panel_span <- function(g) {
  ix <- grepl("^panel", g$layout$name)
  list(
    top = min(g$layout$t[ix]), left = min(g$layout$l[ix]),
    right = max(g$layout$r[ix])
  )
}

in_h <- function(u) convertHeight(u, "in", valueOnly = TRUE)
in_w <- function(u) convertWidth(u, "in", valueOnly = TRUE)

# How far right a title reaches from its card's left edge. ggplot positions a
# title against the panel, not the plot, so a card with a wide axis gutter
# starts its title well inside the card and the gutter counts toward the
# reach. Measuring both parts is what lets the layout widen the card to fit
# rather than clip the title, which is what a plot-width-only pack does.
title_w <- function(g) {
  ix <- which(g$layout$name == "title")
  if (!length(ix)) return(0)
  l <- g$layout$l[ix]
  offset <- if (l > 1) in_w(sum(g$widths[seq_len(l - 1)])) else 0
  offset + in_w(grobWidth(g$grobs[[ix]]))
}

above_in <- function(g) {
  top <- panel_span(g)$top
  if (top > 1) in_h(sum(g$heights[seq_len(top - 1)])) else 0
}

set_panel_h <- function(g, h) {
  ix <- grepl("^panel", g$layout$name)
  g$heights[unique(g$layout$t[ix])] <- unit(h, "in")
  g
}

# Pad the furniture flanking the panels so stacked plots share one panel
# column. All arithmetic stays in grid units so it resolves on the drawing
# device, not on whatever device happened to be open when the cell was built.
fit_width <- function(gs) {
  side <- function(g, cols) {
    if (length(cols) == 0) unit(0, "pt") else sum(g$widths[cols])
  }
  if (length(gs) < 2) return(gs)
  lw <- map2(gs, lapply(gs, panel_span), \(g, s) side(g, seq_len(s$left - 1)))
  lmax <- do.call(max, lw)
  gs <- map2(gs, lw, \(g, l) gtable_add_cols(g, lmax - l, pos = 0))
  rw <- map2(gs, lapply(gs, panel_span), function(g, s) {
    side(g, seq_len(ncol(g))[-seq_len(s$right)])
  })
  rmax <- do.call(max, rw)
  map2(gs, rw, \(g, r) gtable_add_cols(g, rmax - r, pos = -1))
}

# Every plot is one cell of a 2 x 4 grid, and the grid only reads as one
# figure if the titles share a line and the panels share two. ggplot stacks a
# title, a strip and an axis above each panel and every one of them is a
# different height, so the gap directly under the title is padded per cell
# until each panel top lands on the same line. Padding under the title rather
# than above it is what keeps the titles themselves in step.
align_above <- function(gs) {
  above <- map_dbl(gs, above_in)
  map2(gs, max(above) - above, function(g, pad) {
    if (pad < 1e-9) return(g)
    at <- g$layout$t[g$layout$name == "title"]
    gtable_add_rows(g, unit(pad, "in"), pos = if (length(at)) max(at) else 0)
  })
}

cell_vp <- function(g, top) {
  viewport(
    x = unit(0, "npc"), y = unit(1, "npc") - unit(top, "in"),
    just = c("left", "top"),
    width = unit(1, "npc"), height = sum(g$heights)
  )
}

# The grid is packed from measured parts, so a plot whose labels grow would
# otherwise run into the row beneath it. Fail instead.
check_grid <- function(gs, col, top) {
  ab <- map_dbl(gs, above_in)
  if (diff(range(ab)) > 1e-3) {
    stop("panel tops differ by ", sprintf("%.3f in", diff(range(ab))),
         call. = FALSE)
  }
  bottom <- top + map_dbl(gs, \(g) in_h(sum(g$heights)))
  walk(unique(col), function(k) {
    ix <- order(top[col == k])
    tp <- top[col == k][ix]
    bt <- bottom[col == k][ix]
    if (length(tp) > 1 && any(tp[-1] - bt[-length(bt)] < -1e-3)) {
      stop("column ", k, ": rows overlap by ",
           sprintf("%.3f in", max(bt[-length(bt)] - tp[-1])), call. = FALSE)
    }
  })
  invisible(TRUE)
}

# The cards borrow two palettes from shared/style.R but cannot source it: it
# runs devices.R and defines size globals that print_scale_apply.R multiplies
# in place. Read the hexes out of the text instead, and fail if either side
# has drifted.
assert_style <- function() {
  style_txt <- paste(
    readLines(here::here("04_Figures", "shared", "style.R"), warn = FALSE),
    collapse = "\n"
  )
  style_hex <- function(const, key) {
    block <- sub(
      paste0("(?s)^.*?", const, "\\s*<-\\s*c\\((.*?)\\).*$"), "\\1",
      style_txt, perl = TRUE
    )
    sub(
      paste0("(?s)^.*?\\b", key, "\\s*=\\s*\"(#[0-9A-Fa-f]{6})\".*$"),
      "\\1", block, perl = TRUE
    )
  }
  style_all <- unlist(regmatches(
    style_txt, gregexpr("#[0-9A-Fa-f]{6}", style_txt)
  ))
  stopifnot(
    unname(AGE_PAL) == c(
      style_hex("CONTRAST_COLORS", "Training_Young"),
      style_hex("CONTRAST_COLORS", "Training_Old")
    ),
    # DIR_PAL deepens the direction pair past DIR_COLORS so it does not read as
    # the age pair three cards along. Both hexes are still style.R's, so the
    # claim to assert is membership, not a name that would not match.
    unname(DIR_PAL) %in% style_all,
    length(rrho_jet()) == 9
  )
  invisible(TRUE)
}

# Both grids draw the same page twice, once per device, and the pdf device is
# resolved at load. Taking the drawing as a function keeps that pairing in one
# place rather than in each assembler.
save_grid <- function(draw, stem, width, height) {
  agg_png(
    file.path(REPORT_DIR, paste0(stem, ".png")),
    width = width, height = height, units = "in", res = DPI,
    background = "transparent"
  )
  draw()
  invisible(dev.off())

  panel_pdf(file.path(REPORT_DIR, paste0(stem, ".pdf")), width, height)
  draw()
  invisible(dev.off())
}

dep_results <- function() {
  read_csv(
    here::here("03_DEP", "c_data", "03_combined_results.csv"),
    show_col_types = FALSE
  )
}

training_pairs <- function() {
  dep_results() |>
    select(
      gene,
      lfc_y = logFC_Training_Young, p_y = P.Value_Training_Young,
      fdr_y = adj.P.Val_Training_Young, pi_y = sig_pi_Training_Young,
      lfc_o = logFC_Training_Old
    ) |>
    filter(!is.na(lfc_y), !is.na(lfc_o))
}

phenotype <- function() read_excel(F01_XLSX, "Per_participant")

module_grid <- function() {
  read_excel(F06_XLSX, "module_grid_summary") |>
    select(row, module, auc, perm_p, q_bh, sig) |>
    mutate(
      tag = MODULE_TAG[module],
      mark = case_when(
        sig == "q<.05" ~ "**", sig == "p<.05" ~ "*", .default = ""
      )
    )
}

# Each card on its own, at the grid's panel height, with a two-plot card
# stacked as it sits in its column. Widened to its title where that is wider.
save_card <- function(plots, stem, width) {
  dir.create(file.path(REPORT_DIR, "panels"), showWarnings = FALSE)
  pdf(NULL)
  gs <- map(plots, ggplotGrob) |>
    map(set_panel_h, h = PANEL_H) |>
    fit_width()
  hs <- map_dbl(gs, \(g) in_h(sum(g$heights)))
  width <- max(width, map_dbl(gs, title_w) + TITLE_PAD)
  invisible(dev.off())
  tops <- cumsum(c(0, head(hs, -1) + ROW_GAP))
  draw <- function() {
    grid.newpage()
    walk2(gs, tops, \(g, top) grid.draw(grobTree(g, vp = cell_vp(g, top))))
  }
  save_grid(draw, file.path("panels", stem), width, max(tops + hs))
}

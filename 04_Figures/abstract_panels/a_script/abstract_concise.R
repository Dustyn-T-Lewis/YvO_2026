#!/usr/bin/env Rscript
# The graphical abstract as one row of six cards, A C B F D E. A backup
# layout; the published abstract was built from abstract.R.

setwd(here::here())
source("04_Figures/abstract_panels/a_script/panels/_common.R")

assert_style()

# A backup layout: the same results as abstract.R on one row. Each card
# takes one plot from a script the two-row grid already sources, so nothing
# here can drift from the main figure except the widths, which are set below
# rather than read from each script -- a script's WIDTH is tuned for its
# column in the 2 x 4 grid and means nothing in a row of six.
CARDS <- tribble(
  ~stem, ~which, ~width, ~basis,
  "panel_A_hypertrophy", 1, 2.20,
  "17 younger, 15 older participants with VL thickness",
  "panel_C_remodelling", 1, 1.55,
  "limma contrasts; 2,106 proteins tested",
  "panel_B_amplitude", 1, 2.05,
  "15 younger, 15 older paired biopsies; 2,106 proteins",
  "panel_F_aging", 1, 1.45,
  "limma aging contrast; 2,106 proteins tested",
  "panel_D_direction", 1, 1.45,
  "2,094 proteins ranked by limma t in both arms",
  "panel_E_modules", 1, 2.05,
  "30 subject eigengenes per module; AUC in-sample"
) |>
  mutate(script = paste0(stem, ".R"))

built <- map(CARDS$script, function(script) {
  e <- new.env(parent = globalenv())
  # The workbook keeps the panel_ names it was published with; the scripts
  # dropped the prefix when they moved to panels/.
  sys.source(
    file.path(SCRIPT_DIR, "panels", sub("^panel_", "", script)),
    envir = e
  )
  list(plots = e$plots, data = e$panel_data)
})

gs <- map2(built, CARDS$which, \(b, i) b$plots[[i]]) |>
  map(ggplotGrob) |>
  map(set_panel_h, h = CONCISE_PANEL_H)

# fit_width pads stacked plots so they share a panel column; with one plot per
# column there is nothing to stack, so only the title line needs aligning.
gs <- align_above(gs)

# Widths above are what each plot needs; a card whose title is wider than
# that takes the title's width instead. Computing it rather than tuning it by
# hand is what keeps a reworded title from silently overlapping its neighbour.
# Text metrics need a device, and none is open until save_grid.
pdf(NULL)
card_w <- pmax(CARDS$width, map_dbl(gs, title_w) + TITLE_PAD)
invisible(dev.off())

cells <- CARDS |>
  mutate(
    col = row_number(),
    row = 1,
    top = 0,
    width = card_w,
    x = cumsum(c(0, head(card_w, -1) + CONCISE_GAP))
  )

check_grid(gs, cells$col, cells$top)

total_w <- max(cells$x + cells$width) + 2 * ROW_PAD
total_h <- max(map_dbl(gs, \(g) in_h(sum(g$heights))))

draw_row <- function() {
  grid.newpage()
  pwalk(list(gs, cells$x, cells$width), function(g, x, w) {
    pushViewport(viewport(
      x = unit(x + ROW_PAD, "in"), y = unit(total_h, "in"),
      width = unit(w, "in"), height = unit(total_h, "in"),
      just = c("left", "top")
    ))
    grid.draw(grobTree(g, vp = cell_vp(g, 0)))
    popViewport()
  })
}

save_grid(draw_row, "abstract_concise", total_w, total_h)

sheets <- c(
  list(Overview = select(CARDS, card = stem, script, plot = which, basis)),
  set_names(map(built, "data"), CARDS$stem)
)

write.xlsx(
  sheets, file.path(DATA_DIR, "abstract_concise_data.xlsx"), overwrite = TRUE
)

message(sprintf(
  "concise row %.2f x %.2f in, %d cards, titles share one line",
  total_w, total_h, nrow(cells)
))

source(here::here("04_Figures", "abstract_panels", "a_script", "_common.R"))

assert_style()

# basis records what each script counted. The scripts do not all rest on the
# same set: B works from paired subject deltas off the imputed matrix, so it
# keeps every protein, while D needs a limma logFC in both arms and loses the
# 12 the model could not fit.
PANELS <- tribble(
  ~stem, ~col, ~basis,
  "panel_A_hypertrophy", 1,
  "17 younger, 15 older participants with VL thickness",
  "panel_B_amplitude", 1,
  "15 younger, 15 older paired biopsies; 2,106 proteins",
  "panel_C_remodelling", 2,
  "limma contrasts; 2,106 proteins tested",
  "panel_D_direction", 3,
  "2,094 proteins with a limma logFC in both arms",
  "panel_E_modules", 4,
  "30 subject eigengenes per module; AUC in-sample"
) |>
  mutate(script = paste0(stem, ".R"))

built <- map(PANELS$script, function(script) {
  e <- new.env(parent = globalenv())
  sys.source(file.path(SCRIPT_DIR, script), envir = e)
  list(plots = e$plots, width = e$WIDTH, data = e$panel_data)
})

# Column 1 holds two scripts of one plot each; every other column holds one
# script whose two plots stack. Flattening here is what turns five scripts
# into the eight cells the grid actually draws.
cells <- tibble(
  stem = rep(PANELS$stem, lengths(map(built, "plots"))),
  col = rep(PANELS$col, lengths(map(built, "plots"))),
  width = rep(map_dbl(built, "width"), lengths(map(built, "plots"))),
  plot = list_flatten(map(built, "plots"))
) |>
  mutate(row = row_number(), .by = col) |>
  mutate(width = max(width), .by = col)

gs <- map(cells$plot, ggplotGrob) |>
  map(set_panel_h, h = PANEL_H)

# Panels line up left to right within a column as well as top to bottom, so
# the furniture either side is padded per column, then the gap under every
# title is padded across all eight so the panel tops share one line.
for (k in unique(cells$col)) {
  ix <- which(cells$col == k)
  gs[ix] <- fit_width(gs[ix])
}
gs <- align_above(gs)

cell_h <- map_dbl(gs, \(g) in_h(sum(g$heights)))
row_h <- map_dbl(1:2, \(r) max(cell_h[cells$row == r]))

cells <- cells |>
  mutate(top = if_else(row == 1, 0, row_h[1] + ROW_GAP))

col_w <- cells |> summarise(w = max(width), .by = col) |> arrange(col)
col_w <- mutate(col_w, x = cumsum(c(0, head(w, -1) + COL_GAP)))
cells <- left_join(cells, select(col_w, col, x), by = "col")

check_grid(gs, cells$col, cells$top)

total_w <- max(cells$x + cells$width)
total_h <- row_h[1] + ROW_GAP + row_h[2]

draw_grid <- function() {
  grid.newpage()
  pwalk(list(gs, cells$x, cells$width, cells$top), function(g, x, w, top) {
    pushViewport(viewport(
      x = unit(x, "in"), y = unit(total_h, "in"),
      width = unit(w, "in"), height = unit(total_h, "in"),
      just = c("left", "top")
    ))
    grid.draw(grobTree(g, vp = cell_vp(g, top)))
    popViewport()
  })
}

save_grid(draw_grid, "GRID_abstract", total_w, total_h)

sheets <- c(
  list(Overview = select(PANELS, panel = stem, script, basis)),
  set_names(map(built, "data"), PANELS$stem)
)

write.xlsx(
  sheets, file.path(DATA_DIR, "abstract_panels.xlsx"), overwrite = TRUE
)

message(sprintf(
  "grid %.2f x %.2f in, %d cells, panel tops share one line",
  total_w, total_h, nrow(cells)
))

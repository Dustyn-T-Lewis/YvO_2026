#!/usr/bin/env Rscript
# Figure 2E: fGSEA pathway counts, as dodged Up/Down bars per contrast stacked
# by database, largest (GO:BP, darkest) at the bottom to smallest (GO Slim).

setwd(here::here())
source("04_Figures/F02/a_script/main/panels/_main.R", local = TRUE)
source("04_Figures/shared/build_fgsea_cache.R")

# Shared fGSEA cache, rebuilt by shared/build_fgsea_cache.R when the stage 03
# results are newer.
fgsea_cache <- "04_Figures/shared/fgsea_tstat_all_v2.csv"
stopifnot(
  "fGSEA cache missing — source shared/build_fgsea_cache.R first" =
    file.exists(fgsea_cache)
)
fgsea_raw <- read_csv(fgsea_cache, show_col_types = FALSE)

DB_ORDER <- c("GO:BP", "Reactome", "Hallmark", "KEGG", "GO Slim")
DISPLAY_CONTRASTS <- c("Aging", "Training_Young", "Training_Old", "Interaction")

shades <- function(cols) {
  setNames(colorRampPalette(cols)(length(DB_ORDER)), DB_ORDER)
}
red_shades <- shades(c("#B2182B", "#D6604D", "#F4A582"))
blue_shades <- shades(c("#2166AC", "#4393C3", "#92C5DE"))

sig_pathways <- fgsea_raw |>
  filter(
    !is.na(padj), padj < 0.05,
    database %in% DB_ORDER, contrast %in% DISPLAY_CONTRASTS
  ) |>
  dplyr::select(contrast, database, pathway, pval, padj, ES, NES, size) |>
  arrange(contrast, database, padj)

write_csv(sig_pathways, file.path(DAT, "panel_E_fgsea_sig.csv"))

count_df <- sig_pathways |>
  mutate(direction = ifelse(NES > 0, "Up", "Down")) |>
  group_by(contrast, direction, database) |>
  summarise(count = n(), .groups = "drop")

sig_counts_wide <- count_df |>
  pivot_wider(names_from = direction, values_from = count, values_fill = 0L) |>
  arrange(contrast, database)
write_csv(sig_counts_wide, file.path(DAT, "panel_E_fgsea_counts.csv"))

full_grid <- expand_grid(
  contrast  = DISPLAY_CONTRASTS,
  direction = c("Up", "Down"),
  database  = DB_ORDER
)
count_df <- full_grid |>
  left_join(count_df, by = c("contrast", "direction", "database")) |>
  mutate(count = replace_na(count, 0))

count_df$database <- factor(count_df$database, levels = DB_ORDER)

BAR_W <- 0.30
DODGE_GAP <- 0.08
ctr_centers <- setNames(seq_along(DISPLAY_CONTRASTS), DISPLAY_CONTRASTS)

count_df <- count_df |>
  mutate(
    x_center = ctr_centers[contrast] +
      ifelse(direction == "Up", -(BAR_W / 2 + DODGE_GAP / 2),
        BAR_W / 2 + DODGE_GAP / 2
      )
  )

count_df <- count_df |>
  arrange(contrast, direction, factor(database, levels = DB_ORDER)) |>
  group_by(contrast, direction) |>
  mutate(
    ymax = cumsum(count),
    ymin = ymax - count
  ) |>
  ungroup()

count_df <- count_df |>
  mutate(fill = ifelse(direction == "Up",
    red_shades[as.character(database)],
    blue_shades[as.character(database)]
  ))

bar_tops <- count_df |>
  group_by(contrast, direction, x_center) |>
  summarise(total = sum(count), .groups = "drop")

# Headroom for the count labels fit to the data: a fixed ceiling of 250 once
# clipped Training_Young's Up bar (261).
y_top <- max(bar_tops$total) * 1.08

PC_W <- 44 # J Physiol: col 2 of 3×2 at 178mm
PC_H <- 55
lbl_sz <- FIG_AXIS_TEXT / .pt

bg_rects <- tibble(
  xmin = seq_along(DISPLAY_CONTRASTS) - 0.5,
  xmax = seq_along(DISPLAY_CONTRASTS) + 0.5,
  fill = CONTRAST_COLORS[DISPLAY_CONTRASTS]
)

p <- ggplot() +
  geom_rect(
    data = bg_rects,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    fill = bg_rects$fill, alpha = 0.20,
    color = "grey70", linewidth = 0.2
  ) +
  geom_rect(
    data = count_df,
    aes(
      xmin = x_center - BAR_W / 2, xmax = x_center + BAR_W / 2,
      ymin = ymin, ymax = ymax
    ),
    fill = count_df$fill, color = "white", linewidth = 0.25
  ) +
  geom_text(
    data = bar_tops |> filter(total > 0),
    aes(x = x_center, y = total, label = total),
    vjust = -0.3, size = lbl_sz, fontface = "bold",
    color = "black"
  ) +
  scale_x_continuous(
    breaks = seq_along(DISPLAY_CONTRASTS),
    labels = CTR_SHORT[DISPLAY_CONTRASTS],
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0)),
    breaks = scales::pretty_breaks(n = 6)
  ) +
  coord_cartesian(clip = "off", ylim = c(0, y_top)) +
  labs(x = NULL, y = "Significant pathways") +
  FIG_THEME +
  theme(
    plot.subtitle = element_text(
      size = FIG_SUBTITLE_SIZE,
      face = "bold.italic", color = "grey40"
    ),
    # Margin respecified, not just hjust: a bare element_text() inherits
    # FIG_THEME's margin(r = -1), which pulls the title into the tick numbers.
    axis.title.y = element_text(hjust = 0.54, margin = margin(r = 1.5)),
    axis.text.x = element_text(
      angle = 35, hjust = 1,
      size = FIG_AXIS_TEXT - 0.5
    ),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    # Minimal margins: E fills its column cell through wrap_elements()
    plot.margin = margin(0, 0, 0, 0)
  )

key_sq_sz <- 1.8 # match panel D key square size
key_txt <- 1.5 # match panel D key font size

# Two keys, database and direction, split at the Tr.(O)/Interaction border.
grey_shades <- shades(c("grey30", "grey75"))
DB_LEGEND_ORDER <- rev(DB_ORDER)
db_df <- tibble(
  label = DB_LEGEND_ORDER, y = -cumsum(c(0, 0.004, 0.004, 0.004, 0.004)),
  fill = grey_shades[DB_LEGEND_ORDER]
)
dir_df <- tibble(
  label = c("Up", "Down"), y = c(0, -0.004),
  fill = unname(DIR_COLORS[c("Up", "Down")])
)

make_key_plot <- function(kdf, ylim) {
  ggplot(kdf) +
    geom_point(aes(x = 0, y = y),
      shape = 22, size = key_sq_sz, fill = kdf$fill,
      color = "grey30", stroke = 0.3
    ) +
    geom_text(aes(x = 0.35, y = y, label = label),
      size = key_txt, color = "grey20", hjust = 0, fontface = "bold"
    ) +
    scale_x_continuous(limits = c(-0.2, 1.5)) +
    coord_cartesian(ylim = ylim, clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

p_key_db <- make_key_plot(db_df, c(min(db_df$y) - 0.005, 0.005))
# The direction key uses panel D's ylim so the two keys space alike.
p_key_dir <- make_key_plot(dir_df, c(-0.010, 0.003))

pe_subtitle <- sprintf(
  "fGSEA, 5 databases, per-db BH | %d / %d sig",
  sum(count_df$count),
  nrow(fgsea_raw |> filter(
    database %in% DB_ORDER,
    contrast %in% DISPLAY_CONTRASTS
  ))
)

# Titles go on the base plot so they anchor to its plot region, as in A and B.
p <- p + labs(title = "Pathway Enrichment (Up/Down)", subtitle = pe_subtitle)

pE <- (p +
  inset_element(p_key_db,
    left = 0.52, right = 0.68,
    top = 1.00, bottom = 0.75
  ) +
  inset_element(p_key_dir,
    left = 0.84, right = 1.00,
    top = 0.98, bottom = 0.85
  )) +
  plot_annotation(theme = theme(plot.margin = margin(t = 6, r = 3, b = 4, l = 3)))

ggsave(file.path(PNL, "E_fgsea.png"), pE,
  width = PC_W, height = PC_H, units = "mm", dpi = 300
)
ggsave(file.path(PNL, "E_fgsea.pdf"), pE,
  width = PC_W, height = PC_H, units = "mm", device = pdf_dev
)
message("F02 Panel E (stacked fGSEA) saved")

invisible(pE)

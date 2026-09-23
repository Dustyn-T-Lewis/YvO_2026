source(here::here("04_Figures", "abstract_panels", "a_script", "_common.R"))
pacman::p_load(RRHO2)

WIDTH <- 1.45
LIM <- c(-2.1, 2.1)

d <- training_pairs() |>
  mutate(agree = sign(lfc_y) == sign(lfc_o))
fdr <- filter(d, fdr_y < 0.05)

stats <- bind_rows(FDR = fdr, All = d, .id = "set") |>
  summarise(
    n = n(), n_agree = sum(agree), pct = mean(agree),
    rho = cor(lfc_y, lfc_o, method = "spearman"),
    .by = set
  ) |>
  mutate(
    block = sprintf(
      "**%s** (%s)<br>%s agree<br><span style='color:%s'>rho %.2f</span>",
      set, comma(n), percent(pct, accuracy = 0.1), MUTED, rho
    )
  )

# Same call F04 panel B makes, so the abstract and the main figure draw one map
# from one set of ranks: shared/comparison_panels/panel_E_rrho2.R.
ranks <- dep_results() |>
  select(gene, ty = t_Training_Young, to = t_Training_Old) |>
  filter(!is.na(ty), !is.na(to))

rrho <- RRHO2_initialize(
  data.frame(gene = ranks$gene, score = ranks$ty),
  data.frame(gene = ranks$gene, score = ranks$to),
  labels = c("Younger", "Older"), log10.ind = TRUE,
  multipleTesting = "none", boundary = 0.02, method = "hyper", stepsize = 20
)

hypermat <- as.data.frame(as.table(rrho$hypermat)) |>
  setNames(c("i", "j", "z")) |>
  mutate(i = as.integer(i), j = as.integer(j), z = pmax(z, 0))

quad_label <- function(dat, x, y, hjust, vjust, fill = NA) {
  geom_richtext(
    data = dat, inherit.aes = FALSE, aes(x, y, label = block),
    hjust = hjust, vjust = vjust, size = pt(GLYPH_PT - 0.5),
    lineheight = 1.25, colour = INK, fill = fill, label.colour = NA,
    label.padding = unit(if (is.na(fill)) rep(0, 4) else rep(1, 4), "pt")
  )
}

# F04's jet ramp, so the two RRHOs speak one colour language. Jet needs a key
# a 0.96 in square has no room for, so the printed percentage calibrates it
# instead: the map's heat tracks n, not effect size, and at n = 2,094 the
# corners burn red on concordance only 1.23x chance.
map_stat <- stats |>
  filter(set == "All") |>
  mutate(x = 3, y = max(hypermat$j) - 2)

glyph_rrho <- ggplot(hypermat, aes(i, j, fill = z)) +
  geom_raster() +
  quad_label(
    map_stat, map_stat$x, map_stat$y,
    hjust = 0, vjust = 1, fill = scales::alpha("white", 0.85)
  ) +
  scale_fill_gradientn(
    colours = rrho_jet(), na.value = rrho_jet()[1], guide = "none"
  ) +
  coord_fixed(expand = FALSE) +
  labs(
    title = "Training moves the same\nproteins in young and old",
    x = expression(bold("younger rank  down" %->% "up")),
    y = expression(bold("older rank"))
  ) +
  theme_panel() +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

quad <- tibble(
  xmin = c(0, LIM[1]), xmax = c(LIM[2], 0),
  ymin = c(0, LIM[1]), ymax = c(LIM[2], 0)
)

glyph_fdr <- ggplot(fdr, aes(lfc_y, lfc_o)) +
  geom_rect(
    data = quad, inherit.aes = FALSE,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "grey80", alpha = 0.45
  ) +
  geom_hline(yintercept = 0, colour = MUTED, linewidth = 0.2) +
  geom_vline(xintercept = 0, colour = MUTED, linewidth = 0.2) +
  geom_abline(
    slope = 1, intercept = 0, linetype = "22", colour = MUTED, linewidth = 0.25
  ) +
  geom_point(
    aes(fill = agree), shape = 21, size = 0.7, stroke = 0.1, colour = OUTLINE
  ) +
  quad_label(
    filter(stats, set == "FDR"), LIM[1] * 0.97, LIM[2] * 0.97,
    hjust = 0, vjust = 1
  ) +
  scale_fill_manual(values = AGREE_PAL) +
  coord_cartesian(xlim = LIM, ylim = LIM, expand = FALSE) +
  scale_x_continuous(breaks = 0) +
  scale_y_continuous(breaks = 0) +
  labs(
    title = paste0(
      "81% of younger\u2019s 135 FDR\n",
      "proteins move alike in older"
    ),
    x = expression(bold("younger log"[2] * "FC")),
    y = expression(bold("older log"[2] * "FC"))
  ) +
  theme_panel()

plots <- list(glyph_rrho, glyph_fdr)
panel_data <- stats |>
  select(set, n, n_agree, pct, rho) |>
  mutate(max_neg_log10_p = max(rrho$hypermat, na.rm = TRUE))

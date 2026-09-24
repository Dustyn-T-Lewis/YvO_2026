#!/usr/bin/env Rscript
# Figure 1A: total training volume over the 12 weeks, younger against older.

setwd(here::here())

pacman::p_load(withr, readxl, dplyr, ggplot2, ggsignif)

source("04_Figures/shared/style.R")

set.seed(42)

RPT <- "04_Figures/F01/b_reports/panels"
DAT <- "04_Figures/F01/c_data"
for (d in c(RPT, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

meta <- read_excel("00_input/YvO_meta.xlsx")
for (col in c(
  "BMI", "Type_I_fCSA", "Type_II_fCSA",
  "deadlift_1rm_kg", "Total_Training_Volume_kg"
)) {
  if (col %in% names(meta) && is.character(meta[[col]])) {
    meta[[col]] <- suppressWarnings(as.numeric(meta[[col]]))
  }
}

meta <- meta |>
  mutate(
    subject_key = sub("_(Pre|Post)$", "", Col_ID),
    Group = factor(Group, levels = c("Young", "Old")),
    Timepoint = factor(Timepoint, levels = c("Pre", "Post"))
  )

tv_df <- meta |>
  filter(Timepoint == "Post") |>
  select(subject_key, Group, tv = Total_Training_Volume_kg) |>
  filter(!is.na(tv))

stats_A <- t.test(tv ~ Group, data = tv_df)
norm_sub <- sprintf("Welch t %s", fmt_p(stats_A$p.value))

audit_A <- tv_df |>
  summarise(
    n = n(), mean = mean(tv), sd = sd(tv), sem = sd / sqrt(n),
    shapiro_p = shapiro.test(tv)$p.value, .by = Group
  ) |>
  mutate(t_test_p = stats_A$p.value)
write.csv(audit_A, file.path(DAT, "panel_A_training_volume.csv"), row.names = FALSE)

tv_df$tv_scaled <- tv_df$tv / 1e5
bar_colors <- c(
  Young = unname(GROUP_FILL["Young_Post"]),
  Old = unname(GROUP_FILL["Old_Post"])
)

pA <- ggplot(tv_df, aes(Group, tv_scaled, fill = Group)) +
  annotate("rect",
    xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
    fill = AGE_COLORS["Young"], alpha = 0.20, color = "grey85", linewidth = 0.15
  ) +
  annotate("rect",
    xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
    fill = AGE_COLORS["Old"], alpha = 0.20, color = "grey85", linewidth = 0.15
  ) +
  geom_bar(stat = "summary", fun = mean, width = 0.6, color = "grey30", linewidth = 0.3) +
  geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2, linewidth = 0.4) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.35, shape = 16, color = "grey30") +
  geom_signif(
    comparisons = list(c("Young", "Old")),
    annotations = fmt_p_plot(stats_A$p.value),
    parse = TRUE, textsize = 1.3, size = 0.3, tip_length = 0.02,
    y_position = max(tv_df$tv_scaled) * 1.15
  ) +
  scale_fill_manual(values = bar_colors) +
  scale_x_discrete(expand = expansion(add = 0.3)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "Training Volume", subtitle = norm_sub,
    y = expression(bold("Volume (×10"^5 * " kg)")), x = NULL
  ) +
  FIG_THEME +
  theme(legend.position = "none")

ggsave(file.path(RPT, "A_training_volume.png"), pA,
  width = 90, height = 150, units = "mm", dpi = 300
)
# The jitter is drawn at render time. Figure 1 used to render only the PNG of
# this panel, so the PDF restores the seed afterwards and panels B and C draw
# the same points they always did.
with_preserve_seed(ggsave(file.path(RPT, "A_training_volume.pdf"), pA,
  width = 90, height = 150, units = "mm", device = get_pdf_device()
))

invisible(pA)

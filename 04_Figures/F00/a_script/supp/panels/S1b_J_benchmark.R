#!/usr/bin/env Rscript
# S1b Figure J: imputation method benchmark.

setwd(here::here())
source("04_Figures/F00/a_script/supp/panels/_inputs.R", local = TRUE)

BENCH_CSV <- "02_imputation/c_data/benchmark/04_composite_ranking.csv"
stopifnot(
  "Benchmark ranking missing - run 02_imputation/a_script/benchmark/_run_all.R" =
    file.exists(BENCH_CSV)
)
bench <- read_csv(BENCH_CSV, show_col_types = FALSE)

bench_plot <- bench |>
  filter(method != "Non_imputed") |>
  arrange(rank) |>
  mutate(
    method = factor(method, levels = rev(method)),
    bar_col = case_when(
      method == "missForest" ~ "Selected",
      rank <= 5 ~ "Top 5",
      TRUE ~ "Other"
    )
  )

mf_rank <- bench$rank[bench$method == "missForest"]
mf_composite <- bench$composite[bench$method == "missForest"]

pJ <- ggplot(bench_plot, aes(composite, method, fill = bar_col)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", composite)),
    hjust = -0.08, size = 1.9, fontface = "bold"
  ) +
  scale_fill_manual(values = c(
    Selected = "#E41A1C", `Top 5` = "#377EB8",
    Other = "grey70"
  ), name = NULL) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.18)),
    limits = c(0, NA)
  ) +
  labs(
    x = "Composite score", y = NULL, tag = "J",
    title = sprintf("Imputation method benchmark (%d)", nrow(bench_plot)),
    subtitle = sprintf(
      "missForest selected (rank #%d, composite = %.3f)",
      mf_rank, mf_composite
    )
  ) +
  FIG_THEME +
  theme(
    legend.position = "top", legend.key.size = unit(2.5, "mm"),
    axis.text.y = element_text(size = 5)
  )

save_panel(pJ, "S1b_J_benchmark", "panel_J", bench,
  height = PH * 1.3
)

invisible(pJ)

#!/usr/bin/env Rscript
# S1b Figure N: DEP counts by contrast and threshold.

setwd(here::here())
source("04_Figures/F00/a_script/panels/_inputs.R", local = TRUE)

DEP_XLSX <- "03_DEP/c_data/03_DEP_results.xlsx"
da_summ <- as.data.frame(read_excel(DEP_XLSX, sheet = "DA_summary"))

dep_counts <- da_summ |>
  filter(type %in% c("up", "down")) |>
  summarise(
    `p<0.05` = sum(sig.PVal),
    `FDR<0.10` = sum(sig.FDR),
    `FDR<0.05` = sum(sig.FDR.05),
    `Pi<0.05` = sum(sig.Pi),
    .by = contrast
  ) |>
  pivot_longer(-contrast, names_to = "threshold", values_to = "n") |>
  mutate(
    contrast = factor(contrast,
      levels = intersect(CONTRAST_ORDER, unique(contrast))
    ),
    threshold = factor(threshold,
      levels = c("p<0.05", "FDR<0.10", "FDR<0.05", "Pi<0.05")
    )
  )

pN <- ggplot(dep_counts, aes(threshold, contrast, fill = n)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = comma(n), color = n > max(n) * 0.5),
    size = 2.5, fontface = "bold"
  ) +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = "grey20"), guide = "none") +
  scale_fill_gradient(low = "#DEEBF7", high = "#08519C", name = "DEPs") +
  labs(
    x = NULL, y = NULL, tag = "N",
    title = "DEP counts by contrast and threshold",
    subtitle = sprintf(
      "limma + duplicateCorrelation | %s proteins",
      comma(int_norm$dal_nrow)
    )
  ) +
  FIG_THEME +
  theme(
    axis.text.x = element_text(angle = 15, hjust = 1, size = 6),
    legend.position = "right", legend.key.size = unit(3, "mm")
  )

save_panel(pN, "S1b_N_dep_heatmap", "panel_N", dep_counts,
  width = PW * 2, height = PH * 0.75
)

invisible(pN)

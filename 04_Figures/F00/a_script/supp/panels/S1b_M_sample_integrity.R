#!/usr/bin/env Rscript
# S1b Figure M: sample means before and after imputation.

setwd(here::here())
source("04_Figures/F00/a_script/supp/panels/_inputs.R", local = TRUE)

samp_integrity <- tibble(
  Col_ID = colnames(int_imp$mat),
  pre = colMeans(int_imp$mat, na.rm = TRUE),
  post = colMeans(int_imp$mat_imp)
) |>
  left_join(int_imp$meta |> select(Col_ID, Group_Time) |> distinct(),
    by = join_by(Col_ID)
  )
r2 <- cor(samp_integrity$pre, samp_integrity$post)^2

pM <- ggplot(samp_integrity, aes(pre, post, color = Group_Time)) +
  geom_abline(
    slope = 1, intercept = 0, linetype = "dashed",
    color = "grey60", linewidth = 0.4
  ) +
  geom_point(size = 1.8, alpha = 0.85) +
  scale_color_manual(
    values = GROUP_FILL, name = NULL,
    labels = c("Young Pre", "Young Post", "Old Pre", "Old Post")
  ) +
  labs(
    x = "Mean log2 (observed)", y = "Mean log2 (imputed)", tag = "M",
    title = "Pre- vs. post-imputation means",
    subtitle = sprintf("R\u00b2 = %.4f | OOB = %.3f", r2, int_imp$oob_error)
  ) +
  FIG_THEME +
  theme(
    legend.position = "top", legend.key.size = unit(2.5, "mm"),
    legend.text = element_text(size = 5.5)
  )

save_panel(pM, "S1b_M_sample_integrity", "panel_M", samp_integrity)

invisible(pM)

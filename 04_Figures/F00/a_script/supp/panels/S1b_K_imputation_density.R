#!/usr/bin/env Rscript
# S1b Figure K: observed and imputed intensity.

setwd(here::here())
source("04_Figures/F00/a_script/supp/panels/_inputs.R", local = TRUE)

obs_vals <- as.numeric(int_imp$mat[!int_imp$was_na])
imp_vals <- as.numeric(int_imp$mat_imp[int_imp$was_na])
dens_df <- bind_rows(
  tibble(value = obs_vals, type = "Observed"),
  tibble(value = imp_vals, type = "Imputed")
) |>
  mutate(type = factor(type, levels = c("Observed", "Imputed")))

k_pal <- c(Observed = "#377EB8", Imputed = "#E41A1C")
pK <- ggplot(dens_df, aes(value, fill = type, color = type)) +
  geom_density(alpha = 0.4, linewidth = 0.4) +
  scale_fill_manual(values = k_pal, name = NULL) +
  scale_color_manual(values = k_pal, name = NULL) +
  labs(
    x = "log2 intensity", y = "Density", tag = "K",
    title = "Observed vs. imputed intensity",
    # OOB sits here because the inset key occupies the top-right corner.
    subtitle = sprintf(
      "%s observed | %s imputed | OOB = %.3f",
      comma(length(obs_vals)), comma(length(imp_vals)), int_imp$oob_error
    )
  ) +
  FIG_THEME +
  theme(legend.position = "top", legend.key.size = unit(2.5, "mm"))

dens_summary <- dens_df |>
  summarise(
    n = n(), mean = mean(value), median = median(value),
    sd = sd(value), .by = type
  ) |>
  mutate(oob_error = ifelse(type == "Imputed", int_imp$oob_error, NA_real_))

save_panel(pK, "S1b_K_imputation_density", "panel_K", dens_summary)

invisible(pK)

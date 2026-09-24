#!/usr/bin/env Rscript
# Triptychs for the five modules S8 Figure leaves out: brown, green, red, pink
# and magenta. No manuscript item draws them; S8.R renders them so each run
# still leaves a triptych for every module.

setwd(here::here())
source("04_Figures/F05/a_script/supp/panels/_triptych.R", local = TRUE)

for (mod in setdiff(KEY_MODULES, c("turquoise", "black", "yellow", "blue"))) {
  triptych_panel(mod, paste0("other_triptychs_", mod))
}

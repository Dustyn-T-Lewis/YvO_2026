# Stage 01 and 02 intermediates, panel size and output folders, shared by
# panels A-N. Each panel saves its render to b_reports/panels and the data
# behind it to c_data/sheets, which F00_data.R folds into S1 Table.

pacman::p_load(dplyr, tibble, tidyr, readr, readxl, ggplot2, scales, stringr)

source("04_Figures/shared/style.R")

BASE <- "04_Figures/F00"
PNL <- file.path(BASE, "b_reports", "supp", "panels")
SHEETS <- file.path(BASE, "c_data", "sheets")
for (d in c(PNL, SHEETS)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

PW <- 89
PH <- 65

int_norm <- readRDS("01_normalization/c_data/00_report_intermediates.rds")
int_imp <- readRDS("02_imputation/c_data/00_report_intermediates.rds")

# Saves one panel as PNG and PDF under its script's name, and the data frame
# its S1 Table sheet holds.
save_panel <- function(p, name, sheet, df, width = PW, height = PH) {
  ggsave(file.path(PNL, paste0(name, ".png")), p,
    width = width, height = height, units = "mm", dpi = 300
  )
  ggsave(file.path(PNL, paste0(name, ".pdf")), p,
    width = width, height = height, units = "mm", device = get_pdf_device()
  )
  saveRDS(df, file.path(SHEETS, paste0(sheet, ".rds")))
}

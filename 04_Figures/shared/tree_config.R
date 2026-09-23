# Repo-wide constants the shared figure library needs: where the DEP results
# are, where the fGSEA cache is, and which paths a figure stitcher must never
# delete. Kept in one file so no figure script hardcodes them.

DEP_RESULTS <- "03_DEP/c_data/03_combined_results.csv"
FGSEA_CACHE <- "04_Figures/shared/fgsea_tstat_all_v2.csv"

# Paths a figure stitcher must never delete when clearing its intermediates.
UPSTREAM_PREFIXES <- c(
  "^00_input/",
  "^01_normalization/",
  "^02_imputation/",
  "^03_DEP/",
  "^04_Figures/shared/"
)

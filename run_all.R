# Runs the pipeline end to end. Each script runs in its own R session so a
# stage cannot inherit objects, masked functions or a stale working directory
# from the one before it. Every script is also runnable on its own — each sets
# its own working directory and loads its own packages — so this file is a
# convenience, not a dependency.
#
# Not --vanilla: that implies --no-init-file, which skips .Rprofile and so
# never activates renv, leaving every step on whatever the system library
# holds. 00_environment.log in each run directory records what was resolved.
#
#   Rscript run_all.R
#
# YvO_WGCNA_run.R builds the network F05 and F06 both read, so it runs before
# either. F06 reads F05's modules as classifier features, so F05 runs before
# F06. F00 summarises stages 01-03, so it runs last.

setwd(here::here())

# Marks a child as part of a full run. Stage 03 rebuilds 03_DEP_results.xlsx
# from scratch, and the five sheets appended after it are only safe to drop
# when the scripts that write them are about to run.
Sys.setenv(YVO_PIPELINE_RUN = "1")

steps <- c(
  "01_normalization/a_script/01_normalize.R",
  "01_normalization/a_script/02_generate_reports.R",
  "01_normalization/a_script/03_fraction_composition.R",
  "02_imputation/a_script/01_impute.R",
  "02_imputation/a_script/02_generate_reports.R",
  "03_DEP/a_script/01_run_dep.R",
  "03_DEP/a_script/02_generate_reports.R",
  "03_DEP/a_script/supp/01_effect_size_robustness.R",
  "03_DEP/a_script/supp/02_supplement_covariate.R",
  "03_DEP/a_script/supp/03_reversal_aging_fdr.R",
  "03_DEP/a_script/supp/04_discordant_ora.R",
  "03_DEP/a_script/supp/05_supplement_sensitivity.R",
  file.path("04_Figures/F05/a_script", c("YvO_WGCNA_run.R", "F05.R", "S6.R", "S8.R", "F05_data.R")),
  file.path("04_Figures/F01/a_script", c("S2.R", "F01.R", "F01_data.R")),
  file.path("04_Figures/F02/a_script", c("F02.R", "S3.R", "F02_data.R")),
  file.path("04_Figures/F03/a_script", c("F03.R", "S4a.R", "S4b.R", "F03_data.R")),
  file.path("04_Figures/F04/a_script", c("F04.R", "S5a.R", "S5b.R", "F04_data.R")),
  file.path("04_Figures/F06/a_script", c("F06.R", "S7.R", "F06_data.R")),
  file.path("04_Figures/F00/a_script", c("S1a.R", "S1b.R", "F00_data.R")),
  "04_Figures/abstract_panels/a_script/abstract.R",
  "04_Figures/abstract_panels/a_script/abstract_concise.R"
)

missing <- steps[!file.exists(steps)]
if (length(missing)) {
  stop("missing pipeline scripts:\n  ", paste(missing, collapse = "\n  "))
}

# The workbooks in 00_input fail loudly in step 1, so they need no preflight.
# The benchmark ranking does: it is written by an opt-in script, read by two
# stages that cannot degrade without it, and its absence is not obvious.
ranking <- "02_imputation/c_data/benchmark/04_composite_ranking.csv"
if (!file.exists(ranking)) {
  stop(
    "missing ", ranking, "\nIt comes from the opt-in benchmark:\n  ",
    "Rscript 02_imputation/a_script/benchmark/_run_all.R"
  )
}

log_dir <- file.path(".runlogs", format(Sys.time(), "run_%Y%m%d_%H%M%S"))
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# One child launched the way the steps are, so the log records the library the
# pipeline actually resolves rather than the one this session happens to see.
system2("Rscript", c("--no-save", "--no-restore", "-e", shQuote(
  'writeLines(c("libPaths:", .libPaths(),
                paste("renv project:", renv::project())))
   renv::status()
   print(sessionInfo())'
)), stdout = file.path(log_dir, "00_environment.log"),
    stderr = file.path(log_dir, "00_environment.log"))

started <- Sys.time()
timings <- data.frame(step = steps, seconds = NA_real_, status = NA_character_)

for (i in seq_along(steps)) {
  step <- steps[i]
  message(sprintf("[%d/%d] %s", i, length(steps), step))
  log_name <- paste0(sprintf("%02d_", i), basename(step), ".log")
  log_file <- file.path(log_dir, log_name)
  t0 <- Sys.time()
  status <- system2("Rscript", c("--no-save", "--no-restore", shQuote(step)),
    stdout = log_file, stderr = log_file
  )
  timings$seconds[i] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  timings$status[i] <- if (status == 0) "ok" else paste0("exit ", status)
  if (status != 0) {
    message("  failed; last 20 lines of ", log_file, ":")
    writeLines(tail(readLines(log_file, warn = FALSE), 20))
    stop("pipeline stopped at ", step)
  }
  message(sprintf("  ok (%.1f s)", timings$seconds[i]))
}

write.csv(timings, file.path(log_dir, "timings.csv"), row.names = FALSE)
message(sprintf(
  "\nAll %d steps completed in %.1f min. Logs in %s",
  nrow(timings),
  as.numeric(difftime(Sys.time(), started, units = "mins")),
  log_dir
))

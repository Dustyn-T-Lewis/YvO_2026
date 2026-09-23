# YvO 2026 — install every package the pipeline needs.
# Run once after cloning:  Rscript setup.R
#
# The stage scripts load what they need with pacman::p_load(), which reaches
# CRAN and nothing else. Bioconductor and GitHub packages have to be installed
# here first, so this is not optional on a fresh machine: stages 01 to 03 do
# not run without proteoDA.

# A fresh R has no mirror set, and install.packages() fails on the "@CRAN@"
# placeholder rather than picking one.
if (identical(unname(getOption("repos")["CRAN"]), "@CRAN@")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

cran <- c(
  "pacman", "here", "withr",
  "readxl", "readr", "openxlsx",
  "dplyr", "tidyr", "tibble", "stringr", "purrr", "forcats", "magrittr",
  "ggplot2", "ggrepel", "ggtext", "ggsignif", "ggbeeswarm",
  "ggforce", "ggnewscale", "patchwork", "cowplot", "gridExtra", "scales", "png",
  "igraph", "ggraph", "tidygraph", "graphlayouts", "concaveman", "circlize",
  "WGCNA", "lme4", "emmeans", "vegan", "boot", "pwr", "deming", "rstatix",
  "pROC", "msigdbr", "missForest"
)

bioc <- c(
  "limma", "fgsea", "GO.db", "org.Hs.eg.db", "AnnotationDbi",
  "impute", "qvalue", "ComplexHeatmap"
)

github <- c(
  proteoDA      = "ByrumLab/proteoDA",
  RRHO2         = "RRHO2/RRHO2",
  # No pipeline script loads enrichVolcano and it is not in renv.lock; F03
  # renders through 04_Figures/shared/volcano_ring.R. Installed for
  # interactive use, so the renv sandbox hides it during a run.
  enrichVolcano = "Dustyn-T-Lewis/enrichVolcano"
)

# Only the opt-in imputation benchmark (RUN_BENCHMARK=1) needs these. Its
# methods skip themselves when a package is absent, so the pipeline runs
# without them.
benchmark_cran <- c(
  "mice", "missMDA", "imputeLCMD", "sn", "abind",
  "glmnet", "foreach", "doParallel", "MASS", "FNN", "randomForest", "rpart"
)
benchmark_bioc <- c("imp4p", "MsCoreUtils", "msImpute", "pcaMethods")
# DreamAI's package root is the repo's Code/ subdirectory, not its top level.
benchmark_github <- c(DreamAI = "WangLab-MSSM/DreamAI/Code")

# Every GitHub package, keyed by package name, so from_github() can look up
# the repo it has to clone.
repos <- c(github, benchmark_github)

install_missing <- function(pkgs, install) {
  missing <- setdiff(pkgs, rownames(installed.packages()))
  if (length(missing)) install(missing)
  invisible(missing)
}

from_cran <- function(pkgs) install.packages(pkgs)
from_bioc <- function(pkgs) {
  install_missing("BiocManager", from_cran)
  BiocManager::install(pkgs, ask = FALSE, update = FALSE)
}
from_github <- function(pkgs) {
  install_missing("remotes", from_cran)
  remotes::install_github(unname(repos[pkgs]))
}

install_missing(cran, from_cran)
install_missing(bioc, from_bioc)
install_missing(names(github), from_github)

# The benchmark is opt-in and every one of its methods skips itself when its
# package is absent, so a failure here must not stop the required installs.
try(install_missing(benchmark_cran, from_cran), silent = TRUE)
try(install_missing(benchmark_bioc, from_bioc), silent = TRUE)
try(install_missing(names(benchmark_github), from_github), silent = TRUE)

still_missing <- setdiff(
  c(cran, bioc, names(github)),
  rownames(installed.packages())
)
if (length(still_missing)) {
  stop("failed to install: ", paste(still_missing, collapse = ", "))
}
message("all required packages present")

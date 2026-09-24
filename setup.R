# Installs the package versions in renv.lock into renv/library, including the
# Bioconductor packages and the three from GitHub (proteoDA, RRHO2, DreamAI).
# Run once after cloning:  Rscript setup.R
#
# .Rprofile activates renv, and renv installs itself first if it is missing.

renv::restore(prompt = FALSE)

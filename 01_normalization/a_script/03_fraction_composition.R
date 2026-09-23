# What compartments did the 500 x g supernatant actually capture? Counts
# understate the contractile apparatus because a handful of myofibrillar
# proteins carry most of the signal, so every compartment is reported by
# protein count and by share of summed linear intensity.

withr::local_dir(here::here())
pacman::p_load(withr, readr, dplyr, tibble, purrr, msigdbr)

dal <- readRDS("01_normalization/c_data/03_DAList_normalized.rds")

abundance <- tibble(
  gene = dal$annotation$gene,
  intensity = rowMeans(2^dal$data, na.rm = TRUE)
) |>
  filter(!is.na(gene), is.finite(intensity))

cc <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:CC")
hallmark <- msigdbr(species = "Homo sapiens", collection = "H")

cc_genes <- function(...) unique(cc$gene_symbol[cc$gs_name %in% c(...)])

sets <- list(
  mitochondrion = cc_genes("GOCC_MITOCHONDRION"),
  contractile = cc_genes(
    "GOCC_CONTRACTILE_MUSCLE_FIBER", "GOCC_MYOSIN_COMPLEX",
    "GOCC_A_BAND", "GOCC_I_BAND", "GOCC_M_BAND", "GOCC_ACTOMYOSIN"
  ),
  sarcoplasmic_reticulum = cc_genes(
    "GOCC_SARCOPLASMIC_RETICULUM_LUMEN",
    "GOCC_JUNCTIONAL_SARCOPLASMIC_RETICULUM_MEMBRANE"
  ),
  cytosolic_ribosome = cc_genes("GOCC_CYTOSOLIC_RIBOSOME"),
  glycolysis = unique(hallmark$gene_symbol[hallmark$gs_name == "HALLMARK_GLYCOLYSIS"])
)

total_intensity <- sum(abundance$intensity)

composition <- imap_dfr(sets, function(members, set_name) {
  hit <- abundance$gene %in% members
  tibble(
    compartment = set_name,
    n_in_set = length(members),
    n_detected = sum(hit),
    pct_of_proteome = 100 * sum(hit) / nrow(abundance),
    pct_of_intensity = 100 * sum(abundance$intensity[hit]) / total_intensity
  )
}) |>
  arrange(desc(pct_of_intensity))

top_carriers <- abundance |>
  slice_max(intensity, n = 15) |>
  mutate(pct_of_intensity = 100 * intensity / total_intensity) |>
  select(gene, pct_of_intensity)

write_csv(composition, "01_normalization/c_data/03_fraction_composition.csv")
write_csv(top_carriers, "01_normalization/c_data/03_fraction_top_carriers.csv")

cat("Proteins quantified:", nrow(abundance), "\n\n")
print(as.data.frame(composition), digits = 3)
cat(
  "\nTop 15 proteins by share of total signal (cumulative",
  round(sum(top_carriers$pct_of_intensity), 1), "%):\n"
)
print(as.data.frame(top_carriers), digits = 3)

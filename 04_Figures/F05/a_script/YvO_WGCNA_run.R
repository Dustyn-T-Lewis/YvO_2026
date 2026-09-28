# Fit WGCNA modules and compute eigengene / trait associations for F05 + F06.

setwd(here::here())

pacman::p_load(WGCNA, tidyverse, lme4, emmeans)
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

disableWGCNAThreads() # single-threaded so correlations/TOM are byte-reproducible across re-runs
set.seed(42)

# WGCNA parameters (Langfelder & Horvath 2008; relaxed R^2 for small-n proteomics).
WGCNA_R2_CUTOFF <- 0.87 # signed R^2 threshold for scale-free topology
WGCNA_NETWORK_TYPE <- "signed"
WGCNA_TOM_TYPE <- "signed"
WGCNA_MIN_MOD_SIZE <- 30L
WGCNA_MERGE_CUT_H <- 0.25

DATA_FILE <- "02_imputation/c_data/01_imputed.csv"
DALIST_RDS <- "02_imputation/c_data/01_DAList_imputed.rds"
# Module assignment of the submitted analysis (YvO_revise @ d34a7dd), used only to keep
# colour labels comparable between the two networks.
REFERENCE_MODULES <- "00_input/wgcna_reference_modules.csv"
DATA_DIR <- "04_Figures/F05/c_data/wgcna"
PANEL_DIR <- "04_Figures/F05/c_data"
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(DATA_FILE), file.exists(DALIST_RDS), file.exists(REFERENCE_MODULES))

df <- read_csv(DATA_FILE)
ann_cols <- c("uniprot_id", "protein", "gene", "description")
ann <- df[, ann_cols]
samp_names <- setdiff(names(df), ann_cols)
mat <- as.matrix(df[, samp_names])
rownames(mat) <- ann$uniprot_id

datExpr <- t(mat)

dal <- readRDS(DALIST_RDS)
dal_meta <- as.data.frame(dal$metadata)

meta <- tibble(
  sample_id = dal_meta$Col_ID,
  subject   = sub("_(Pre|Post)$", "", dal_meta$Col_ID),
  age       = dal_meta$Group,
  time      = dal_meta$Timepoint,
  group     = dal_meta$Group_Time
)

gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  ann <- ann |> filter(uniprot_id %in% colnames(datExpr))
  message(sprintf(
    "After goodSamplesGenes: %d samples x %d proteins",
    nrow(datExpr), ncol(datExpr)
  ))
}

# WGCNA needs its own cor()
cor <- WGCNA::cor

powers <- 1:20
sft <- pickSoftThreshold(datExpr,
  powerVector = powers,
  networkType = WGCNA_NETWORK_TYPE, verbose = 2
)
saveRDS(sft$fitIndices, file.path(DATA_DIR, "sft_fitIndices.rds"))

# Conventional threshold is R^2 > 0.90 (Zhang & Horvath 2005), relaxed to
# WGCNA_R2_CUTOFF for small-n proteomics per Langfelder & Horvath (2008).
#
# On this cohort the two criteria agree, which is worth stating in Methods
# rather than leaning on the relaxed cutoff alone: 0.87 selects power 12, and
# 12 is what the WGCNA FAQ recommends for a signed network with more than 40
# samples (this is 62). Requiring 0.90 would pick 13, one step further out,
# halving mean connectivity from 5.37 to 3.91.
r2_values <- -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq
power_idx <- which(r2_values > WGCNA_R2_CUTOFF)[1]
# The fallback matters only if no power clears the cutoff. 12 is the FAQ's
# signed-network recommendation for n > 40; 6 is its unsigned/signed-hybrid
# figure and would be wrong for this network type.
soft_power <- if (!is.na(power_idx)) powers[power_idx] else 12
# Log-log slope diagnostic: scale-free topology expects slope ~ -1 to -2
sft_slope <- sft$fitIndices$slope[soft_power]
# WGCNA's FAQ recommends power 12 for a signed network with n > 40. Saying so
# when the cutoff happens to agree is a stronger justification than the relaxed
# R^2 alone, and it is the line Methods should quote.
faq_power <- 12L
r2_note <- paste0(
  ifelse(sft_slope > -1, " [NOTE: slope > -1, weak scale-free fit]", ""),
  ifelse(r2_values[soft_power] < 0.90,
    " [below the 0.90 convention; acceptable for small-n, Langfelder & Horvath 2008]", ""
  ),
  if (soft_power == faq_power && nrow(datExpr) > 40) {
    sprintf(
      " [agrees with the WGCNA FAQ: power %d for a signed network, n = %d > 40]",
      faq_power, nrow(datExpr)
    )
  } else {
    ""
  }
)
message(sprintf(
  "  Soft power: %d (R^2 = %.3f, slope = %.2f%s)",
  soft_power, r2_values[soft_power], sft_slope, r2_note
))

# Pearson chosen over bicor. Bicor sensitivity in panels/S6_D_bicor.R confirms concordance.
net <- blockwiseModules(
  datExpr,
  power             = soft_power,
  networkType       = WGCNA_NETWORK_TYPE,
  TOMType           = WGCNA_TOM_TYPE,
  minModuleSize     = WGCNA_MIN_MOD_SIZE,
  mergeCutHeight    = WGCNA_MERGE_CUT_H,
  numericLabels     = TRUE,
  pamRespectsDendro = FALSE,
  verbose           = 3
)

module_colors <- labels2colors(net$colors)
n_modules <- length(unique(net$colors)) - (0 %in% net$colors)

message(sprintf("  Modules detected: %d (+ grey/unassigned)", n_modules))

# WGCNA names modules by size rank at every refit, so a colour is a label and not an
# identity: two modules that swap rank swap colours even when their proteins do not
# move. Relabel against the submitted network so colour names in the manuscript keep
# pointing at the same proteins.
ref <- read_csv(REFERENCE_MODULES, show_col_types = FALSE)
shared <- intersect(colnames(datExpr), ref$protein)
src <- module_colors[match(shared, colnames(datExpr))]
matched <- as.character(matchLabels(src, ref$module_color[match(shared, ref$protein)]))
relabel <- vapply(split(matched, src), \(x) names(sort(table(x), decreasing = TRUE))[1], character(1))
stopifnot(
  setequal(names(relabel), unique(module_colors)),
  !anyDuplicated(relabel)
)
swapped <- names(relabel)[names(relabel) != relabel]
message(sprintf(
  "  Matched to submitted network on %d shared proteins: %d of %d colours held",
  length(shared), length(relabel) - length(swapped), length(relabel)
))
if (length(swapped)) {
  message(sprintf(
    "    %s",
    paste(sprintf("%s -> %s", swapped, relabel[swapped]), collapse = ", ")
  ))
}
module_colors <- unname(relabel[module_colors])

# Dendrogram plot is produced in panels/S6_B_dendrogram.R.

traits <- meta |>
  mutate(
    age_num     = if_else(age == "Old", 1, 0),
    time_num    = if_else(time == "Post", 1, 0),
    interaction = age_num * time_num
  )

# The phenotypes go on the trait matrix and on meta, which the panels read.
pheno_cols <- c(
  "VL_thick_cm", "DXA_LBM_kg", "BMI",
  "Type_I_fCSA", "Type_II_fCSA", "deadlift_1rm_kg"
)
for (pc in intersect(pheno_cols, names(dal_meta))) {
  vals <- dal_meta[[pc]]
  if (!is.numeric(vals)) vals <- as.numeric(as.character(vals))
  traits[[pc]] <- meta[[pc]] <- vals[match(meta$sample_id, dal_meta$Col_ID)]
}

trait_cols <- intersect(c("age_num", "time_num", "interaction", pheno_cols), names(traits))
traits_mat <- as.data.frame(traits[, trait_cols])
rownames(traits_mat) <- meta$sample_id
traits_mat <- traits_mat[rownames(datExpr), ]

MEs <- moduleEigengenes(datExpr, colors = module_colors)$eigengenes
MEs <- orderMEs(MEs)

module_trait_cor <- cor(MEs, traits_mat, use = "pairwise.complete.obs")
n_per_trait <- colSums(!is.na(traits_mat))
module_trait_pval <- module_trait_cor
for (j in seq_len(ncol(module_trait_cor))) {
  module_trait_pval[, j] <- corPvalueStudent(module_trait_cor[, j], n_per_trait[j])
}
module_trait_pval_bh <- matrix(p.adjust(as.vector(module_trait_pval), method = "BH"),
  nrow = nrow(module_trait_pval), dimnames = dimnames(module_trait_pval)
)

write_csv(
  as.data.frame(module_trait_cor) |> rownames_to_column("module"),
  file.path(DATA_DIR, "wgcna_module_trait_correlations.csv")
)
write_csv(
  as.data.frame(module_trait_pval_bh) |> rownames_to_column("module"),
  file.path(DATA_DIR, "wgcna_module_trait_pvalues_bh.csv")
)

kME <- signedKME(datExpr, MEs)

module_df <- tibble(
  uniprot_id   = colnames(datExpr),
  module_color = module_colors,
  module_num   = net$colors
) |> left_join(ann |> dplyr::select(uniprot_id, gene), by = "uniprot_id")

hub_rows <- list()
unique_modules <- setdiff(unique(module_colors), "grey")

for (mod in unique_modules) {
  mod_proteins <- module_df$uniprot_id[module_df$module_color == mod]
  kme_col <- paste0("kME", mod)
  if (!(kme_col %in% colnames(kME))) next

  mod_kme <- tibble(
    uniprot_id = rownames(kME),
    kME = kME[, kme_col]
  ) |>
    filter(uniprot_id %in% mod_proteins) |>
    arrange(desc(abs(kME))) |>
    head(10) |>
    mutate(module = mod) |>
    left_join(ann |> dplyr::select(uniprot_id, gene), by = "uniprot_id")

  hub_rows <- c(hub_rows, list(mod_kme))
}

hub_df <- bind_rows(hub_rows)
cor <- stats::cor # restore after WGCNA computations

bg_genes <- ann$gene[ann$uniprot_id %in% colnames(datExpr)]
bg_genes <- unique(bg_genes[!is.na(bg_genes) & bg_genes != ""])

pw_collection <- build_pathway_collection(min_size = 15, include_goslim = FALSE)

ora_results_list <- list()

for (mod in unique_modules) {
  mod_genes <- module_df$gene[module_df$module_color == mod]
  mod_genes <- unique(mod_genes[!is.na(mod_genes) & mod_genes != ""])

  if (length(mod_genes) < 5) next

  ora_res <- tryCatch(
    run_ora_deduplicated(
      genes = mod_genes, universe = bg_genes, pathways = pw_collection,
      min_size = 15, padj_cutoff = 0.10
    ),
    error = function(e) {
      warning(sprintf("ORA failed for '%s': %s", mod, e$message))
      NULL
    }
  )

  if (!is.null(ora_res) && nrow(ora_res) > 0) {
    ora_res$module <- mod
    ora_res$Description <- clean_pathway_name(ora_res$pathway)
    ora_res$geneID <- vapply(ora_res$overlapGenes, paste, character(1), collapse = "/")
    ora_res$Count <- ora_res$overlap
    ora_res$p.adjust <- ora_res$padj
    ora_res$ID <- ora_res$pathway
    ora_results_list <- c(ora_results_list, list(ora_res))
  }
}

enrich_df <- bind_rows(ora_results_list)
enrich_df$overlapGenes <- NULL

write_csv(module_df, file.path(DATA_DIR, "wgcna_module_assignments.csv"))
write_csv(hub_df, file.path(DATA_DIR, "wgcna_hub_proteins.csv"))
saveRDS(net, file.path(DATA_DIR, "wgcna_network.rds"))
write_csv(enrich_df, file.path(DATA_DIR, "wgcna_module_enrichment.csv"))

key_mod_counts <- enrich_df |>
  count(module, sort = TRUE) |>
  head(5) |>
  pull(module)
writeLines(key_mod_counts, file.path(DATA_DIR, "key_modules.txt"))
writeLines(key_mod_counts[nzchar(trimws(key_mod_counts))], file.path(PANEL_DIR, "key_modules.txt"))

sft_summary <- tibble(
  selected_power    = soft_power,
  R_squared         = r2_values[soft_power],
  mean_connectivity = sft$fitIndices$mean.k.[soft_power],
  n_proteins        = ncol(datExpr),
  n_samples         = nrow(datExpr)
)
write_csv(sft_summary, file.path(DATA_DIR, "wgcna_sft_summary.csv"))

meta$group <- factor(meta$group,
  levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post")
)

mod_sizes <- sort(table(module_colors[module_colors != "grey"]), decreasing = TRUE)

# Every module label is two halves from different evidence on purpose,
# "<enrichment term> | <hub family>", derived in _module_labels.R and shipped as
# sheets WGCNA_module_core and WGCNA_module_ora (where abbreviations such as
# "PMF" or "Gluc." resolve). One label per module, used unchanged on the figure
# and in the manuscript; it is authored here because no trimming rule knows
# that "60S" is the half of "60S Ribosome" worth keeping.
#
# The hub half is the largest protein family among the proteins with
# kME >= 0.6, the module's own definition of a core. It ranges from a
# 48-of-74 respiratory chain in green to a 6-of-53 20S proteasome in brown, and
# the percentage is on the sheet so a thin one reads as thin. Families are
# broad, not gene symbols: "eEF1A" named two of black's ten hubs and reads as
# eIF1A, a different protein in another module.
#
# The enrichment half is a GO:BP term, but not the module's top one. The top
# term puts annotation artefacts on the figure: red's is "Negative Regulation
# Of Syncytium Formation" at padj 9.6e-41, whose whole overlap is ribosomal
# proteins, and blue's is "The Role Of GTSE1 In G2/M Progression", whose
# overlap is tubulins and proteasome subunits; neither happens in post-mitotic
# muscle. So each module's significant terms are first rolled up to a shared
# GO:BP ancestor, which picks the branch the module sits on, and the label
# takes the most significant term beneath it. The ancestor itself was tried as
# the label and climbed too far: blue reached "Macromolecule Metabolic
# Process", brown stopped one step from the GO root.
#
# Three labels are not what the most-significant-term rule alone returns, and
# the sheet's ancestor_or_descendant column says so for each. Green skips its
# own branch: that branch is ATP synthesis, which is pink's half of the same
# process, and green is the proton gradient that drives it. Turquoise is set by
# hand, because its branch is carboxylic-acid catabolism and the fatty-acid half
# of that is what "Beta-Oxidation" already says. Brown comes from the tie-break:
# 0.0003 in adjusted p separated gluconeogenesis from the proteasomal term, a
# gap too small to name a module on.
#
# Black keeps a label with nothing significant behind it. Its best GO:BP term is
# padj 0.083 and the sheet records significant = FALSE. An 11-protein core that
# fails to reform in 24 of 30 leave-one-subject-out refits has little to
# annotate, so the weak result belongs in the legend rather than in a search for
# a threshold that would hide it.
#
# ora_label, the top over-representation term, stays in the table below. It is
# off the figure, but it is what the roll-up is preferred to, and a reader
# checking the choice needs both. If the evidence column stops matching
# bio_label, the modules have moved (colours follow size rank; see above).
#
# The order is load-bearing: with interchangeable halves black led with
# "Transl. Fact.", a biology its enrichment never shows, resting on three
# proteins, over the glycolysis term that ranked first. The convention is stated
# in the panel subtitle rather than prefixed onto all nine labels.
#
# Panel A's wrap_at_rule() splits at the pipe and wraps each half to
# max(10, round(1.25 * sqrt(n_proteins))) characters, 23 on the widest bar and
# 10 on the narrowest; magenta binds.
#
# Blue's hub half is the weakest claim here. Its membership is 25 chaperones,
# but five of its top ten hubs are cytoskeletal (MSN, DYNC1H1, TUBB, TUBA1B,
# MAPT) against two chaperones. "Chaperones" describes the module; a strict
# reading of the hubs would say "Cytoskeleton".
mod_label_lookup <- c(
  turquoise = "Lipid Catabolism | TCA Cycle",
  blue      = "Protein Folding | Chaperones",
  brown     = "Protein Degrad. | Proteasome",
  yellow    = "Muscle Contraction | Sarcomere",
  green     = "Oxid. Phos. | Complex I",
  red       = "Translation | 60S Ribo.",
  black     = "Glycolysis | Transl. Fact.",
  pink      = "Chemiosm. | ATP Synth.",
  magenta   = "Initiation | 40S Ribo."
)

# Trim a pathway name to something that fits inside a bar. Mechanical, so a refit
# that changes the top term still produces a readable label.
short_ora <- function(x) {
  x |>
    str_remove("^(Reference|The Role Of) ") |>
    str_replace("Negative Regulation Of", "Neg. Reg.") |>
    str_replace("Formation Of ATP By Chemiosmotic Coupling", "ATP by Chemiosmosis") |>
    str_replace("ATP Synthesis Coupled Electron Transport", "ATP Synth. Coupled ETC") |>
    str_replace("Activation Of The Mrna.*", "mRNA Activation") |>
    str_replace("Gtse1 In G2 M Progression.*", "GTSE1 in G2/M") |>
    str_replace("Syncytium Formation.*", "Syncytium Form.") |>
    str_replace("Catabolic Process", "Catab.") |>
    str_replace("Striated Muscle Contraction", "Striated Muscle Contr.") |>
    str_squish()
}

mod_top_ora <- enrich_df |>
  filter(dedup_status == "kept", module %in% names(mod_sizes)) |>
  slice_min(padj, n = 1, by = module, with_ties = FALSE) |>
  transmute(module, ora_label = short_ora(Description), ora_padj = padj)

mod_evidence <- enrich_df |>
  filter(dedup_status == "kept", module %in% names(mod_sizes)) |>
  arrange(module, padj, pathway) |>
  summarise(evidence = paste(head(Description, 5), collapse = "; "), .by = module)

mod_bio_labels <- tibble(
  module_color = names(mod_sizes),
  module_id = paste0("M", seq_along(mod_sizes)),
  bio_label = mod_label_lookup[names(mod_sizes)],
  n_proteins = as.integer(mod_sizes)
) |>
  left_join(mod_top_ora, by = c(module_color = "module")) |>
  left_join(mod_evidence, by = c(module_color = "module")) |>
  mutate(display_label = paste0(str_to_title(module_color), ": ", bio_label))

stopifnot(
  "a module colour has no curated label - re-read the ORA profile" =
    !anyNA(mod_bio_labels$bio_label),
  "a module returned no enriched term for the ORA slot" =
    !anyNA(mod_bio_labels$ora_label)
)

# Ten stored hubs is a thin basis for naming a 340-protein module, and kME is
# the module's own internal definition rather than an external annotation, so
# every assigned protein gets its correlation with its own eigengene written
# out. Purely derived from MEs and datExpr above: no network call, no seed.
assigned_df <- module_df |>
  filter(module_color != "grey") |>
  left_join(mod_bio_labels[, c("module_color", "module_id")], by = "module_color")

own_kme <- kME[cbind(
  match(assigned_df$uniprot_id, rownames(kME)),
  match(paste0("kME", assigned_df$module_color), colnames(kME))
)]

stopifnot(
  "a protein has no kME against its own module" = !anyNA(own_kme)
)

kme_all <- assigned_df |>
  mutate(kME = own_kme) |>
  arrange(match(module_color, mod_bio_labels$module_color), desc(kME)) |>
  mutate(kME_rank = row_number(), .by = module_color) |>
  dplyr::select(uniprot_id, gene, module_color, module_id, kME, kME_rank)

write_csv(kme_all, file.path(PANEL_DIR, "wgcna_kme_all.csv"))

saveRDS(MEs, file.path(PANEL_DIR, "MEs.rds"))
saveRDS(kME, file.path(PANEL_DIR, "kME_all.rds"))
saveRDS(datExpr, file.path(PANEL_DIR, "datExpr.rds"))
saveRDS(module_colors, file.path(PANEL_DIR, "module_colors.rds"))
write_csv(meta, file.path(PANEL_DIR, "meta.csv"))
write_csv(mod_bio_labels, file.path(PANEL_DIR, "mod_bio_labels.csv"))
write_csv(ann, file.path(PANEL_DIR, "imp_annotations.csv"))

# Per-gene z-scores averaged within group_time (rows = gene, cols = group level).
# Consumed by panels/_triptych.R for per-module gene heatmaps.
expr_g <- t(datExpr)
rownames(expr_g) <- ann$gene[match(rownames(expr_g), ann$uniprot_id)]
expr_g <- expr_g[!is.na(rownames(expr_g)) & rownames(expr_g) != "", ]
expr_g <- expr_g[!duplicated(rownames(expr_g)), ]
z_g <- t(scale(t(expr_g)))
group_z <- vapply(
  levels(meta$group),
  function(g) rowMeans(z_g[, meta$sample_id[meta$group == g], drop = FALSE], na.rm = TRUE),
  numeric(nrow(z_g))
)
saveRDS(group_z, file.path(PANEL_DIR, "group_z.rds"))

pre_meta <- meta |>
  filter(time == "Pre") |>
  mutate(subject_key = sub("_(Pre|Post)$", "", sample_id))
post_meta <- meta |>
  filter(time == "Post") |>
  mutate(subject_key = sub("_(Pre|Post)$", "", sample_id))

me_pre_raw <- MEs[pre_meta$sample_id, , drop = FALSE]
me_post_raw <- MEs[post_meta$sample_id, , drop = FALSE]
pre_subjects <- pre_meta$subject_key
rownames(me_pre_raw) <- pre_subjects
rownames(me_post_raw) <- post_meta$subject_key

common_subj <- intersect(rownames(me_pre_raw), rownames(me_post_raw))
me_pre <- me_pre_raw[common_subj, , drop = FALSE]
me_post <- me_post_raw[common_subj, , drop = FALSE]
delta_me <- me_post - me_pre

pheno_pre <- pre_meta |>
  dplyr::select(
    subject_key, VL_thick_cm, DXA_LBM_kg, BMI,
    deadlift_1rm_kg, Type_I_fCSA, Type_II_fCSA
  ) |>
  rename(
    VL_Pre = VL_thick_cm, LBM_Pre = DXA_LBM_kg, BMI_Pre = BMI,
    DL_Pre = deadlift_1rm_kg, T1_Pre = Type_I_fCSA, T2_Pre = Type_II_fCSA
  )

pheno_post <- post_meta |>
  dplyr::select(
    subject_key, VL_thick_cm, DXA_LBM_kg,
    deadlift_1rm_kg, Type_I_fCSA, Type_II_fCSA
  ) |>
  rename(
    VL_Post = VL_thick_cm, LBM_Post = DXA_LBM_kg,
    DL_Post = deadlift_1rm_kg, T1_Post = Type_I_fCSA, T2_Post = Type_II_fCSA
  )

pheno_wide <- inner_join(pheno_pre, pheno_post, by = "subject_key") |>
  mutate(
    delta_VL = VL_Post - VL_Pre,
    delta_LBM = LBM_Post - LBM_Pre,
    delta_DL = DL_Post - DL_Pre,
    delta_T1 = T1_Post - T1_Pre,
    delta_T2 = T2_Post - T2_Pre
  ) |>
  filter(subject_key %in% common_subj)

subj_age <- pre_meta |>
  dplyr::select(subject_key, age) |>
  distinct()

delta_vl_vec <- pheno_wide$delta_VL[match(common_subj, pheno_wide$subject_key)]
delta_lbm_vec <- pheno_wide$delta_LBM[match(common_subj, pheno_wide$subject_key)]
delta_dl_vec <- pheno_wide$delta_DL[match(common_subj, pheno_wide$subject_key)]
delta_t1_vec <- pheno_wide$delta_T1[match(common_subj, pheno_wide$subject_key)]
delta_t2_vec <- pheno_wide$delta_T2[match(common_subj, pheno_wide$subject_key)]

pred_cor <- tibble(module = colnames(me_pre)) |>
  rowwise() |>
  mutate(
    r_vl = cor(me_pre[common_subj, module], delta_vl_vec, use = "complete.obs"),
    r_lbm = cor(me_pre[common_subj, module], delta_lbm_vec, use = "complete.obs"),
    max_r = max(abs(r_vl), abs(r_lbm), na.rm = TRUE)
  ) |>
  ungroup()

top3 <- pred_cor |>
  arrange(desc(max_r)) |>
  head(3) |>
  pull(module)

all_mods <- colnames(me_pre_raw)

baseline_traits <- data.frame(
  BMI_Pre = pre_meta$BMI[match(pre_subjects, pre_meta$subject_key)],
  VL_Pre = pre_meta$VL_thick_cm[match(pre_subjects, pre_meta$subject_key)],
  LBM_Pre = pre_meta$DXA_LBM_kg[match(pre_subjects, pre_meta$subject_key)],
  row.names = pre_subjects
)
change_traits <- data.frame(
  delta_VL  = delta_vl_vec,
  delta_LBM = delta_lbm_vec,
  row.names = common_subj
)

# Eigengene-trait Pearson r and Student p over the given subjects. stats::cor
# (restored above), pairwise-complete, with each trait's own non-missing n.
cor_p <- function(me, traits) {
  r <- cor(me, traits, use = "pairwise.complete.obs")
  p <- matrix(NA_real_,
    nrow = ncol(me), ncol = ncol(traits),
    dimnames = list(colnames(me), colnames(traits))
  )
  for (trait in colnames(traits)) {
    p[, trait] <- corPvalueStudent(r[, trait], sum(!is.na(traits[[trait]])))
  }
  list(r = r, p = p)
}

bl_all <- cor_p(me_pre_raw[pre_subjects, all_mods, drop = FALSE], baseline_traits)
bl_cor_mat <- bl_all$r
bl_pval_mat <- bl_all$p
ch_all <- cor_p(delta_me[common_subj, all_mods, drop = FALSE], change_traits)
ch_cor_mat <- ch_all$r
ch_pval_mat <- ch_all$p

young_pre_subj <- subj_age$subject_key[subj_age$age == "Young"]
old_pre_subj <- subj_age$subject_key[subj_age$age == "Old"]

bl_subj_y <- intersect(pre_subjects, young_pre_subj)
bl_y <- cor_p(me_pre_raw[bl_subj_y, all_mods, drop = FALSE], baseline_traits[bl_subj_y, , drop = FALSE])
bl_cor_young <- bl_y$r
bl_pval_young <- bl_y$p

bl_subj_o <- intersect(pre_subjects, old_pre_subj)
bl_o <- cor_p(me_pre_raw[bl_subj_o, all_mods, drop = FALSE], baseline_traits[bl_subj_o, , drop = FALSE])
bl_cor_old <- bl_o$r
bl_pval_old <- bl_o$p

common_young <- intersect(common_subj, young_pre_subj)
ch_y <- cor_p(delta_me[common_young, all_mods, drop = FALSE], change_traits[common_young, , drop = FALSE])
ch_cor_young <- ch_y$r
ch_pval_young <- ch_y$p

common_old <- intersect(common_subj, old_pre_subj)
ch_o <- cor_p(delta_me[common_old, all_mods, drop = FALSE], change_traits[common_old, , drop = FALSE])
ch_cor_old <- ch_o$r
ch_pval_old <- ch_o$p

message(sprintf(
  "  Stratified: Young baseline n=%d, Old baseline n=%d, Young change n=%d, Old change n=%d",
  length(bl_subj_y), length(bl_subj_o), length(common_young), length(common_old)
))

# LMM contrasts model the repeated measures: eigengene ~ group + (1|subject).
# meta$group is already a factor in design order.
lmm_data <- meta

lmm_contrast_list <- list(
  Aging          = c(-1, 0, 1, 0),
  Training_Young = c(-1, 1, 0, 0),
  Training_Old   = c(0, 0, -1, 1),
  Interaction    = c(1, -1, -1, 1)
)

lmm_rows <- list()
for (mod in all_mods) {
  lmm_data[[mod]] <- MEs[lmm_data$sample_id, mod]

  fit <- tryCatch(
    suppressWarnings(
      lmer(as.formula(paste0("`", mod, "` ~ group + (1 | subject)")), data = lmm_data)
    ),
    error = function(e) {
      warning(sprintf("LMM failed for %s: %s", mod, e$message))
      NULL
    }
  )
  if (is.null(fit)) next

  singular <- isSingular(fit)
  if (singular) message(sprintf("  Note: singular fit for %s (near-zero random variance)", mod))

  emm <- emmeans(fit, ~group)

  for (cname in names(lmm_contrast_list)) {
    s <- summary(contrast(emm, list(ctr = lmm_contrast_list[[cname]])), ddf = "Kenward-Roger")
    lmm_rows <- c(lmm_rows, list(tibble(
      module   = mod,
      contrast = cname,
      estimate = round(s$estimate, 5),
      SE       = round(s$SE, 5),
      df       = round(s$df, 2),
      t_ratio  = round(s$t.ratio, 4),
      p_raw    = s$p.value,
      r_equiv  = round(sign(s$estimate) * sqrt(s$t.ratio^2 / (s$t.ratio^2 + s$df)), 4),
      singular = singular
    )))
  }
}

lmm_df <- bind_rows(lmm_rows)
message(sprintf(
  "  LMM contrasts: %d tests (%d modules x %d contrasts)",
  nrow(lmm_df), length(all_mods), length(lmm_contrast_list)
))

# BH per-section for LMM (confirmatory); per-column for stratified (exploratory)
lmm_df$p_bh <- p.adjust(lmm_df$p_raw, method = "BH")

per_col_bh <- function(pmat) {
  for (j in seq_len(ncol(pmat))) pmat[, j] <- p.adjust(pmat[, j], method = "BH")
  pmat
}
bl_pval_bh_young <- per_col_bh(bl_pval_young)
bl_pval_bh_old <- per_col_bh(bl_pval_old)
ch_pval_bh_young <- per_col_bh(ch_pval_young)
ch_pval_bh_old <- per_col_bh(ch_pval_old)

n_l <- nrow(lmm_df)
compat_raw <- c(lmm_df$p_raw, as.vector(bl_pval_mat), as.vector(ch_pval_mat))
compat_bh <- p.adjust(compat_raw, method = "BH")
n_b <- length(as.vector(bl_pval_mat))
n_c <- length(as.vector(ch_pval_mat))
bl_pval_bh <- matrix(compat_bh[(n_l + 1):(n_l + n_b)],
  nrow = nrow(bl_pval_mat),
  dimnames = dimnames(bl_pval_mat)
)
ch_pval_bh <- matrix(compat_bh[(n_l + n_b + 1):(n_l + n_b + n_c)],
  nrow = nrow(ch_pval_mat),
  dimnames = dimnames(ch_pval_mat)
)

write_csv(lmm_df, file.path(DATA_DIR, "wgcna_lmm_contrast_check.csv"))

# Stratified LMM: ME ~ time + (1|subject) within Young and Old separately
strat_rows <- list()
for (age_grp in c("Young", "Old")) {
  strat_data <- lmm_data |>
    filter(grepl(age_grp, group)) |>
    mutate(time = factor(ifelse(grepl("Post", group), "Post", "Pre"),
      levels = c("Pre", "Post")
    ))

  for (mod in all_mods) {
    strat_data[["me_val"]] <- MEs[strat_data$sample_id, mod]

    fit_s <- tryCatch(
      suppressWarnings(
        lmer(me_val ~ time + (1 | subject), data = strat_data)
      ),
      error = function(e) NULL
    )
    if (is.null(fit_s)) next

    s_s <- summary(contrast(emmeans(fit_s, ~time), list(training = c(-1, 1))), ddf = "Kenward-Roger")
    t_s <- s_s$t.ratio
    df_s <- s_s$df
    strat_rows <- c(strat_rows, list(tibble(
      age_group = age_grp,
      module    = mod,
      estimate  = round(s_s$estimate, 5),
      SE        = round(s_s$SE, 5),
      df        = round(df_s, 2),
      t_ratio   = round(t_s, 4),
      p_raw     = s_s$p.value,
      r_equiv   = round(sign(s_s$estimate) * sqrt(t_s^2 / (t_s^2 + df_s)), 4),
      singular  = isSingular(fit_s)
    )))
  }
}
strat_df <- bind_rows(strat_rows)
strat_df$p_bh <- p.adjust(strat_df$p_raw, method = "BH")
write_csv(strat_df, file.path(DATA_DIR, "wgcna_lmm_stratified_check.csv"))
message(sprintf(
  "  Stratified LMM: %d tests (%d modules x %d age groups)",
  nrow(strat_df), length(all_mods), 2
))

for (out in list(
  list(bl_cor_mat, "wgcna_baseline_trait_correlations.csv"),
  list(bl_pval_bh, "wgcna_baseline_trait_pvalues_bh.csv"),
  list(ch_cor_mat, "wgcna_change_trait_correlations.csv"),
  list(ch_pval_bh, "wgcna_change_trait_pvalues_bh.csv")
)) {
  write_csv(as.data.frame(out[[1]]) |> rownames_to_column("module"), file.path(DATA_DIR, out[[2]]))
}

# Per-module pick of the strongest BH-significant trait association across:
# simple correlation (all samples), delta correlation (paired all/young/old),
# and the LMM Aging contrast. Falls back to age_num if nothing reaches BH<0.05.
# Consumed by _supp_mod_hub.R to label hub-network gene-significance scaling.
to_long <- function(mat, value_name) {
  as.data.frame(mat) |>
    rownames_to_column("module_me") |>
    pivot_longer(-module_me, names_to = "trait", values_to = value_name) |>
    mutate(module = sub("^ME", "", module_me)) |>
    dplyr::select(-module_me)
}
join_cor_p <- function(cor_mat, p_mat, label) {
  inner_join(to_long(cor_mat, "gs_r"), to_long(p_mat, "gs_pval_bh"),
    by = c("module", "trait")
  ) |>
    mutate(source_section = label)
}

gs_modules <- mod_bio_labels$module_color

gs_lmm <- lmm_df |>
  filter(contrast == "Aging") |>
  transmute(
    module = sub("^ME", "", module),
    trait = "age_num",
    gs_r = r_equiv, gs_pval_bh = p_bh,
    source_section = "LMM Aging"
  )

gs_sig <- bind_rows(
  join_cor_p(module_trait_cor, module_trait_pval_bh, "Simple correlation"),
  join_cor_p(ch_cor_mat, ch_pval_bh, "Delta correlation"),
  join_cor_p(ch_cor_young, ch_pval_bh_young, "Delta Young"),
  join_cor_p(ch_cor_old, ch_pval_bh_old, "Delta Old"),
  gs_lmm
) |>
  filter(module %in% gs_modules, !is.na(gs_pval_bh), gs_pval_bh < 0.05) |>
  group_by(module) |>
  slice_min(gs_pval_bh, n = 1, with_ties = FALSE) |>
  ungroup() |>
  rename(gs_phenotype = trait)

gs_default <- tibble(
  module = setdiff(gs_modules, gs_sig$module),
  gs_phenotype = "age_num",
  gs_r = NA_real_, gs_pval_bh = NA_real_,
  source_section = "Default"
)

gs_phenotype_choices <- bind_rows(gs_sig, gs_default) |>
  mutate(rationale = ifelse(
    source_section == "Default",
    "no significant traits; aging study default",
    sprintf("%s p_bh=%.3g", source_section, gs_pval_bh)
  )) |>
  arrange(match(module, gs_modules)) |>
  dplyr::select(module, gs_phenotype, gs_r, gs_pval_bh, source_section, rationale)

write_csv(gs_phenotype_choices, file.path(DATA_DIR, "gs_phenotype_choices.csv"))

strat_outputs <- list(
  list(prefix = "baseline_trait_correlations", young = bl_cor_young, old = bl_cor_old),
  list(prefix = "baseline_trait_pvalues_bh", young = bl_pval_bh_young, old = bl_pval_bh_old),
  list(prefix = "baseline_trait_pvalues_raw", young = bl_pval_young, old = bl_pval_old),
  list(prefix = "change_trait_correlations", young = ch_cor_young, old = ch_cor_old),
  list(prefix = "change_trait_pvalues_bh", young = ch_pval_bh_young, old = ch_pval_bh_old),
  list(prefix = "change_trait_pvalues_raw", young = ch_pval_young, old = ch_pval_old)
)
for (out in strat_outputs) {
  for (age in c("young", "old")) {
    write_csv(
      as.data.frame(out[[age]]) |> rownames_to_column("module"),
      file.path(DATA_DIR, paste0("wgcna_", out$prefix, "_", age, ".csv"))
    )
  }
}

saveRDS(me_pre, file.path(PANEL_DIR, "me_pre.rds"))
saveRDS(me_post, file.path(PANEL_DIR, "me_post.rds"))
saveRDS(delta_me, file.path(PANEL_DIR, "delta_me.rds"))
write_csv(pheno_wide, file.path(PANEL_DIR, "pheno_wide.csv"))
write_csv(subj_age, file.path(PANEL_DIR, "subj_age.csv"))

saveRDS(list(
  top3           = top3,
  common_subj    = common_subj,
  pre_subjects   = pre_subjects,
  mod_bio_labels = setNames(mod_bio_labels$bio_label, mod_bio_labels$module_color),
  outcome_labels = c(delta_VL = "Delta VL (cm)", delta_LBM = "Delta LBM (kg)"),
  delta_vl_vec   = delta_vl_vec,
  delta_lbm_vec  = delta_lbm_vec,
  delta_dl_vec   = delta_dl_vec,
  delta_t1_vec   = delta_t1_vec,
  delta_t2_vec   = delta_t2_vec
), file.path(PANEL_DIR, "shared_objects.rds"))

message(sprintf(
  "Done: %d modules, %d hub proteins, %d enriched pathways, %d paired subjects",
  n_modules, nrow(hub_df), nrow(enrich_df), length(common_subj)
))

#!/usr/bin/env Rscript
# F01 — participant characteristics tables for the manuscript.
#
# Three tables rather than one, because the variables carry different column
# shapes and forcing them into a single grid leaves half the cells empty.
#
#   1A  measured once per participant: age, body mass, BMI, training volume.
#       Two age groups, one Welch t-test.
#   1B  measured before and after training. Four cells and the 2x2 mixed
#       ANOVA the F01 panels report, so a table cell and a figure p-value
#       come from one calculation.
#   1C  which parent trial and supplement arm each group is made of. No arm
#       spans both age groups, so there is nothing to test.
#
# Group sizes live in the column headers, and 1B carries per-variable n in its
# own column, so no cell is ever blank.
#
# Fat-to-lean is reported instead of percent body fat: NORE recorded body mass
# only at baseline, so a percent-of-weight measure at post-training would rest
# on a derived denominator for eleven of the 32 participants.

withr::local_dir(here::here())

pacman::p_load(withr, readxl, dplyr, tidyr, stringr, purrr)

DAT <- "04_Figures/F01/c_data"
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

meta <- read_excel("00_input/YvO_meta.xlsx") |>
  mutate(
    subject_key = sub("_(Pre|Post)$", "", Col_ID),
    Group = factor(Group, levels = c("Young", "Old")),
    Timepoint = factor(Timepoint, levels = c("Pre", "Post")),
    FM_to_LBM = DXA_FM_kg / DXA_LBM_kg
  )

n_young <- sum(meta$Timepoint == "Pre" & meta$Group == "Young")
n_old <- sum(meta$Timepoint == "Pre" & meta$Group == "Old")

fmt_p <- function(p) if (p < 0.001) "<0.001" else sprintf("%.3f", p)
dec1 <- function(x) sprintf("%.1f ± %.1f", mean(x), sd(x))

# 1A — one measurement per participant.
single <- function(var, label, fmt = dec1, scale = 1, tp = "Pre") {
  d <- meta |> filter(Timepoint == tp, !is.na(.data[[var]]))
  y <- d[[var]][d$Group == "Young"] / scale
  o <- d[[var]][d$Group == "Old"] / scale
  tibble(
    Variable = label,
    Young = fmt(y), Old = fmt(o),
    p = fmt_p(t.test(y, o)$p.value)
  )
}

table_1a <- bind_rows(
  single("Age", "Age (y)"),
  single("Weight_kg", "Body mass (kg)"),
  single("BMI", "BMI (kg/m²)"),
  single("Total_Training_Volume_kg", "Training volume, 12 wk (×10⁵ kg)",
    scale = 1e5, tp = "Post"
  )
)

# 1B — assembled from the *_summary.csv each panel script writes, so the table
# and the figure report one calculation rather than two.
SUMMARIES <- c(
  "panel_B_dxa_lbm_summary.csv", "panel_C_vl_thickness_summary.csv",
  "supp/panel_A_deadlift_1rm_summary.csv", "supp/panel_B_type_II_fcsa_summary.csv",
  "supp/panel_C_type_I_fcsa_summary.csv", "supp/panel_A_dxa_fat_mass_summary.csv",
  "supp/panel_B_fat_to_lean_summary.csv"
)

missing <- SUMMARIES[!file.exists(file.path(DAT, SUMMARIES))]
if (length(missing)) {
  stop("Run the panel scripts first; missing: ", paste(missing, collapse = ", "))
}

DIGITS <- c(
  DXA_LBM_kg = 1, DXA_FM_kg = 1, FM_to_LBM = 2, VL_thick_cm = 2,
  Type_I_fCSA = 0, Type_II_fCSA = 0, deadlift_1rm_kg = 1
)
LABELS <- c(
  DXA_LBM_kg = "DXA lean body mass (kg)", DXA_FM_kg = "DXA fat mass (kg)",
  FM_to_LBM = "Fat-to-lean mass ratio", VL_thick_cm = "VL thickness (cm)",
  Type_I_fCSA = "Type I fCSA (\u00b5m\u00b2)", Type_II_fCSA = "Type II fCSA (\u00b5m\u00b2)",
  deadlift_1rm_kg = "Deadlift 1RM (kg)"
)
ORDER <- names(LABELS)

cell_fmt <- function(m, s, d) sprintf("%.*f \u00b1 %.*f", d, m, d, s)

table_1b <- map_dfr(file.path(DAT, SUMMARIES), read.csv) |>
  mutate(d = DIGITS[variable]) |>
  transmute(
    variable,
    Variable = LABELS[variable],
    n = sprintf("%d / %d", n_young, n_old),
    Young_Pre = cell_fmt(young_pre_mean, young_pre_sd, d),
    Young_Post = cell_fmt(young_post_mean, young_post_sd, d),
    Old_Pre = cell_fmt(old_pre_mean, old_pre_sd, d),
    Old_Post = cell_fmt(old_post_mean, old_post_sd, d),
    p_Age = vapply(p_age, fmt_p, ""),
    p_Time = vapply(p_time, fmt_p, ""),
    p_Age_x_Time = vapply(p_age_x_time, fmt_p, "")
  ) |>
  arrange(match(variable, ORDER)) |>
  select(-variable)

table_1c <- meta |>
  filter(Timepoint == "Pre") |>
  count(Group, parent_study, supplement, name = "n") |>
  arrange(Group, parent_study, supplement) |>
  rename(`Age group` = Group, `Parent trial` = parent_study, Supplement = supplement)

hdr <- function(d) {
  names(d) <- names(d) |>
    str_replace("^Young$", sprintf("Young (n=%d)", n_young)) |>
    str_replace("^Old$", sprintf("Old (n=%d)", n_old)) |>
    str_replace_all("_", " ")
  d
}
table_1a <- hdr(table_1a)
table_1b <- hdr(table_1b)

REPEATED <- c(
  "DXA_LBM_kg", "DXA_FM_kg", "FM_to_LBM", "VL_thick_cm",
  "Type_I_fCSA", "Type_II_fCSA", "deadlift_1rm_kg"
)

static <- meta |>
  filter(Timepoint == "Pre") |>
  select(
    Participant = subject_key, Age_group = Group, Parent_trial = parent_study,
    Parent_ID = parent_id, Supplement = supplement, Sex, Age,
    Body_mass_kg = Weight_kg, BMI
  )

volume <- meta |>
  filter(Timepoint == "Post") |>
  select(Participant = subject_key, Training_volume_kg = Total_Training_Volume_kg)

repeated_wide <- meta |>
  select(subject_key, Timepoint, all_of(REPEATED)) |>
  pivot_wider(
    names_from = Timepoint, values_from = all_of(REPEATED),
    names_glue = "{.value}_{Timepoint}"
  ) |>
  rename(Participant = subject_key)

for (v in REPEATED) {
  repeated_wide[[paste0(v, "_Delta")]] <-
    repeated_wide[[paste0(v, "_Post")]] - repeated_wide[[paste0(v, "_Pre")]]
}

per_subject <- static |>
  left_join(volume, by = "Participant") |>
  left_join(repeated_wide, by = "Participant") |>
  select(
    Participant, Age_group, Parent_trial, Parent_ID, Supplement, Sex, Age,
    Body_mass_kg, BMI, Training_volume_kg,
    all_of(as.vector(t(outer(REPEATED, c("Pre", "Post", "Delta"), paste, sep = "_"))))
  ) |>
  arrange(Age_group, Participant)

notes <- tibble::tribble(
  ~Item, ~Detail,
  "Cohort", sprintf("%d younger and %d older men. All participants are male, so sex is not tabulated.", n_young, n_old),
  "Design", "Twelve weeks of supervised full-body resistance training, twice weekly, measured before (Pre) and after (Post).",
  "Table 1A", "Variables measured once per participant. Young against Old by Welch two-sample t-test.",
  "Table 1B", "Variables measured at both timepoints. Group x Timepoint mixed ANOVA on participants with a complete pair; n is given per variable as Young / Old.",
  "Table 1C", "Parent trial and supplement arm making up each age group. No supplement arm spans both age groups, so no test is reported.",
  "Values", "Mean +/- SD throughout.",
  "Fat-to-lean", "DXA fat mass divided by DXA lean body mass. Used in place of percent body fat because the NORE trial recorded body mass only at baseline, so a percent-of-weight measure after training would rest on a derived denominator for 11 of the 32 participants.",
  "Missing data", "Fibre cross-sectional areas and deadlift 1RM are incomplete in the parent trials. The per-variable n in Table 1B counts participants with both timepoints.",
  "Participant ids", "Participant is the proteomics sample label. Parent_ID is the participant number inside the parent trial workbook, at 00_input/parent_meta/<Parent_trial>.xlsx.",
  "Source", "00_input/YvO_meta.xlsx, via 04_Figures/F01/a_script/04_phenotype_table.R.",
  "Abbreviations", "DXA, dual-energy x-ray absorptiometry; LBM, lean body mass; VL, vastus lateralis; fCSA, fibre cross-sectional area; 1RM, one-repetition maximum; BMI, body mass index."
)

stopifnot(
  "no blank cells in 1A" = !any(is.na(table_1a) | table_1a == ""),
  "no blank cells in 1B" = !any(is.na(table_1b) | table_1b == ""),
  "no blank cells in 1C" = !any(is.na(table_1c))
)

# The workbook is assembled by 90_stitch_F01.R, which sources this script and
# combines these tables with the panel source data into one supplementary file.
write.csv(table_1a, file.path(DAT, "F01_table_1a_characteristics.csv"), row.names = FALSE)
write.csv(table_1b, file.path(DAT, "F01_table_1b_pre_post.csv"), row.names = FALSE)
write.csv(table_1c, file.path(DAT, "F01_table_1c_composition.csv"), row.names = FALSE)

cat("\n## 1A Participant characteristics\n")
print(as.data.frame(table_1a), row.names = FALSE)
cat("\n## 1B Pre/post training\n")
print(as.data.frame(table_1b), row.names = FALSE)
cat("\n## 1C Cohort composition\n")
print(as.data.frame(table_1c), row.names = FALSE)

# shared

Code and two caches that the figure scripts source; no figure of its own.

## Reads

- `03_DEP/c_data/03_combined_results.csv`, for `build_fgsea_cache.R`
- `04_Figures/F05/c_data/wgcna/wgcna_module_assignments.csv` and `mod_bio_labels.csv`, for `build_string_cluster_cache.R`

## Writes

- `fgsea_tstat_all_v2.csv`: fGSEA on the limma t-statistics for all four contrasts, rebuilt by `build_fgsea_cache.R` when older than `03_combined_results.csv`. Read by F02, F03 and F04.
- `04_Figures/F05/c_data/wgcna_string_clusters.csv`: STRING clusters per module. `build_string_cluster_cache.R` calls the STRING API only when module membership changed, so a normal run needs no network.

## Files

- `style.R`: palettes, `FIG_THEME`, text sizes, `strip_for_composite()`, `source_panel()`, `fmt_p()`. It sources `tree_config.R`, `module_palette.R` and `devices.R`.
- `devices.R`: `get_pdf_device()` and `open_pdf()`, which keep Greek letters in PDFs, and `caption_supp()`, which writes `<item>_Figure.pdf`; also sourced by stages 01 to 03.
- `figure_supplement_helpers.R`: `build_workbook()` and `cleanup_after_workbook()`, used by every `<DIR>_data.R`.
- `supplement_overview.R`: the Overview sheet of the stage workbooks.
- `pathway_utils.R`, `enrichment_dedup.R`: gene-set collections, fGSEA and ORA with redundancy removal.
- `volcano_ring.R`: F03. `go_slim_categories.R`, `print_scale_apply.R`: F04.
- `comparison_panels/`: plot engines configured by a `cfg` list; F03 uses `panel_heatmap_classified.R`, F04 and the abstract the rest.
- `goslim_generic.obo`: the GO Slim definitions.

## Order

Only the two cache builders run, sourced by the figure scripts. `print_scale_apply.R` changes `style.R`'s size globals without restoring them; F04's panel A sources it, so `F04.R` re-sources `style.R` after its panels. `module_palette.R` is separate so the abstract gets the module colours without those globals.

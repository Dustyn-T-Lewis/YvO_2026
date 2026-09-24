# shared · the figure library

Palettes, helpers and two data caches that the F00-F06 scripts source. Not a pipeline
stage; nothing here runs on its own.

```
style.R                      palettes, FIG_THEME, composite_text_sizes(), fmt_p(),
                             strip_for_composite(), boot_median_ci()
  +- tree_config.R           DEP_RESULTS, FGSEA_CACHE, UPSTREAM_PREFIXES
  +- module_palette.R        the nine WGCNA module colours
  +- devices.R               get_pdf_device(), embed_pdf_fonts(), open_pdf()
figure_supplement_helpers.R  build_workbook(), cleanup_after_workbook(), read_sheet_df()
pathway_utils.R              build_pathway_collection(), run_fgsea_deduplicated()
  +- enrichment_dedup.R      run_ora_deduplicated(), classify_database()
supplement_overview.R        add_overview(), write_sheet()
print_scale_apply.R  volcano_ring.R  go_slim_categories.R
build_fgsea_cache.R          -> fgsea_tstat_all_v2.csv (1.4 MB)
comparison_panels/           panel_C_trajectory, panel_D_nes_scatter, panel_E_rrho2,
                             panel_fry_barcode, panel_heatmap_classified
goslim_generic.obo
```

## What comes out

`style.R` is sourced by every panel script in F00-F06 and pulls in `tree_config.R`,
`module_palette.R` and `devices.R` with it. `figure_supplement_helpers.R` is sourced by
every `90_stitch_F0x.R`. `devices.R` and `supplement_overview.R` are also sourced outside
04_Figures, by stages 01, 02 and 03.

`volcano_ring.R` is F03 only; `go_slim_categories.R` and `print_scale_apply.R` F04 only.
Of `comparison_panels/`, F04 sources four and F03 sources `panel_heatmap_classified.R`.
Each takes a `cfg` list defined by the caller.

## Orderings that matter

`print_scale_apply.R` mutates `style.R`'s exported size globals in place and never restores
them, so anything sourced after it at a different canvas width must re-source `style.R`
first. `F04/a_script/02_supp_panels.R` does.

`module_palette.R` is a separate file because `abstract_panels/a_script/panels/_common.R` needs
the module colours but cannot source `style.R`, which runs `devices.R` and defines the
globals `print_scale_apply.R` rewrites.

`build_fgsea_cache.R` rebuilds `fgsea_tstat_all_v2.csv` only when
`03_DEP/c_data/03_combined_results.csv` is newer. Readers such as `F02/_panel_E_fgsea.R`
assert the file exists rather than computing it.

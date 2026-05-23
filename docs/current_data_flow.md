# Current Data Flow Audit

This audit was generated before targeted refactoring of the clusterProfiler to compareGO contract. The repository contains R scripts across preprocessing, ID mapping, QC, differential-expression/enrichment, EWCE, WGCNA, spatial networks, behavior coupling, PRIDE helpers, testing and deprecated folders.

## Global Findings

- Many active scripts still contain machine-specific paths such as `S:/...`, `/Users/...`, and `C:/...`; these are documented by `rg` inventory and should be migrated module by module.
- Active scripts write a mix of canonical analysis tables, exploratory workbooks, figure-only exports, RDS caches, network files and logs directly under legacy `Results/` and `Datasets/` style folders.
- `setwd()` is present in active ID mapping and GCT/pathview workflows; path helpers should be used instead, and unavoidable working-directory changes should restore `oldwd`.
- Package auto-installation is present in several scripts. The refactored priority scripts now expose `AUTO_INSTALL_MISSING_PACKAGES` and default to fail-fast behavior.
- Per-script input/output/dependency/refactor status is maintained in `docs/active_script_io_audit.tsv`.

## Module Refactor Status

- Refactored first: `04_differential_expression_enrichment/01_clusterProfiler.r` and `04_differential_expression_enrichment/02_compareGO.r`, because this handoff had the highest stale-file and duplicate-file risk.
- Safe next candidates: `05_celltype_enrichment_EWCE/01_EWCE_E9.r`, `06_modules_WGCNA/02_module_spatial_networks.r`, `06_modules_WGCNA/03_overlap_modules.r`, `07_spatial_networks/03_bootstrap_network_stability.r`, `07_spatial_networks/04_bootstrap_differential_network_stability.r`, and `07_spatial_networks/05_bootstrap_differential_network_figures.r`; these can mostly be refactored through parameter/config path changes.
- Higher-risk candidates: `06_modules_WGCNA/01_WGCNA v.2.0.0.r`, `07_spatial_networks/01_network_spatial_relations.r`, and `08_behavior_physio_coupling/01_correlate_proteomics_with_behavior.r`; these are central producers or depend on external behavior data and should be changed with input data available.
- Keep documented unless revived: `05_celltype_enrichment_EWCE/90_EWCE_legacy.r`, `90_testing/`, and `99_deprecated/`.

## clusterProfiler to compareGO Before Refactor

`04_differential_expression_enrichment/01_clusterProfiler.r` scanned `cfg$paths$mapped_dir` for mapped contrast CSVs. Each filename was split on `_` to infer a two-token comparison. The script classified each comparison into:

- `phenotype_within_unit`
- `phenotype_between_unit`
- `within_unit_same_phenotype`
- `between_unit_and_phenotype`
- `unclassified`

The route unit was inferred from the comparison tokens and used in nested output folders. This biological routing is preserved.

For each mapped contrast file, clusterProfiler read a gene/protein identifier column and `log2fc` or `logFC`, ran GO GSEA, optional GO simplification, ORA, KEGG GSEA, optional predefined KEGG/pathview, and optional custom NK3R GSEA. It wrote:

- full GO GSEA tables under per-comparison `GO/<ontology>/`
- simplified GO GSEA tables when enabled
- ORA full/simplified tables under `ORA/`
- KEGG tables under `KEGG/`
- plot-used GO GSEA tables under `Datasets/core_enrichment/<ontology>/...`
- plot-used KEGG tables under both the current ontology core-enrichment folder and `Datasets/core_enrichment/KEGG/...`
- plots under legacy routed `Plots/`
- QC summaries and completion flags under routed result folders

`04_differential_expression_enrichment/02_compareGO.r` then selected one `ont`, `ensemble_profiling` route category and `condition` route unit using hard-coded variables. It recursively listed all CSV files in `Datasets/core_enrichment/<ont>/<ensemble_profiling>/<condition>`, named comparisons from filenames, and combined all rows. Separately, it listed mapped/log2FC CSVs from `Datasets/mapped/<ensemble_profiling>/<condition>` for volcano, significant protein, and gene-centric analyses.

Risks in the old flow:

- stale core-enrichment CSVs could be reused silently
- duplicate full/simplified/KEGG/GO files could be mixed by directory pattern alone
- mapped/log2FC files were not formally linked to the enrichment tables
- ontology and route filters were hard-coded and not recorded as a run contract
- empty enrichment tables and missing required columns were not validated before downstream plotting

## clusterProfiler to compareGO After Refactor

clusterProfiler writes canonical source data and a manifest:

`data/processed/04_differential_expression_enrichment/clusterProfiler/clusterProfiler_manifest.csv`

compareGO reads that manifest, filters by ontology/result type/route, validates required enrichment columns, verifies mapped/log2FC input files, writes the consumed subset to:

`data/processed/04_differential_expression_enrichment/compareGO/compareGO_input_manifest.csv`

This makes the data flow explicit and reproducible while preserving the biological route classification. By default, compareGO consumes only `result_type == GSEA_GO` and `used_for_plot == TRUE`. `GSEA_KEGG` rows are recorded by clusterProfiler for provenance and can be selected intentionally by `config/compareGO_config.yml`, but ORA and custom/NK3R outputs are not compareGO inputs unless future analysis logic explicitly supports them.

Checkpoint behavior: if clusterProfiler skips a completed comparison, it reconstructs manifest rows from existing canonical GSEA_GO/GSEA_KEGG source-data tables when present and marks `checkpoint_status = reconstructed_from_checkpoint`. If a checkpoint predates canonical source-data outputs, rerun with `force_rerun: true` to refresh the manifest.

Dry-run commands:

```bash
Rscript 04_differential_expression_enrichment/01_clusterProfiler.r --dry-run
Rscript 04_differential_expression_enrichment/02_compareGO.r --dry-run
```

## PRIDE and Journal Layers

Repository outputs should be separated as:

- canonical processed quantitative matrices: processed proteomics matrices suitable for PRIDE supplementary processed-data upload
- differential-expression result tables: canonical contrast-level statistics suitable for supplements
- enrichment-derived secondary analyses: GO/KEGG/EWCE/WGCNA/network derivatives generally suited for journal supplements, not raw PRIDE redundancy
- figure-only outputs: SVG/PNG/PDF files excluded from PRIDE unless specifically requested by the journal
- source-data tables: exact data behind figures under `results/source_data/`
- exploratory analyses: testing/deprecated outputs excluded from PRIDE and usually excluded from supplements

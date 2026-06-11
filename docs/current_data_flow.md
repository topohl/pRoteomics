# Current Data Flow Audit

This audit was generated before targeted refactoring of the clusterProfiler to compareGO contract. The repository contains R scripts across preprocessing, ID mapping, QC, differential abundance/enrichment, EWCE, WGCNA, spatial networks, behavior coupling, PRIDE helpers, testing and deprecated folders.

## Global Findings

- Many active scripts still contain machine-specific paths such as `S:/...`, `/Users/...`, and `C:/...`; these are documented by `rg` inventory and should be migrated module by module.
- Active scripts write a mix of canonical analysis tables, exploratory workbooks, figure-only exports, RDS caches, network files and logs directly under legacy `Results/` and `Datasets/` style folders.
- `setwd()` is present in active ID mapping and GCT/pathview workflows; path helpers should be used instead, and unavoidable working-directory changes should restore `oldwd`.
- Package auto-installation is present in several scripts. The refactored priority scripts now expose `AUTO_INSTALL_MISSING_PACKAGES` and default to fail-fast behavior.
- Per-script input/output/dependency/refactor status is maintained in `docs/active_script_io_audit.tsv`.
- Follow-up scientific guardrails and refactor candidates are tracked in `docs/refactor_backlog.md`.

## Module Refactor Status

- Refactored first: `04_differential_expression_enrichment/01_clusterProfiler.r` and `04_differential_expression_enrichment/02_compareGO.r`, because this handoff had the highest stale-file and duplicate-file risk.
- Safe next candidates: `05_celltype_enrichment_EWCE/01_EWCE_E9.r`, `06_modules_WGCNA/06_module_spatial_networks.r`, `07_spatial_networks/03_bootstrap_network_stability.r`, `07_spatial_networks/04_bootstrap_differential_network_stability.r`, and `07_spatial_networks/05_bootstrap_differential_network_figures.r`; these can mostly be refactored through parameter/config path changes.
- Higher-risk candidates at the time of the audit were the WGCNA model producer, `07_spatial_networks/01_network_spatial_relations.r`, and `08_behavior_physio_coupling/01_correlate_proteomics_with_behavior.r`; these are central producers or depend on external behavior data and should be changed with input data available.
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

`data/processed/04_differential_expression_enrichment/clusterProfiler/<dataset>/clusterProfiler_manifest.csv`

The current dataset-aware flow is one biological dataset family to one mapped folder to one clusterProfiler manifest/run to one compareGO run. For example, `neuron_neuropil` maps to `data/processed/02_id_mapping/mapped/neuron_neuropil/forward/per_file/`, then clusterProfiler writes the `neuron_neuropil` manifest and source data, and compareGO filters that manifest by `dataset == "neuron_neuropil"`.

compareGO reads that manifest, requires the `dataset` column unless `legacy_mode: true` is explicitly configured, filters by dataset/ontology/result type/route, validates required enrichment columns, verifies mapped/log2FC input files, writes the consumed subset to:

`data/processed/04_differential_expression_enrichment/compareGO/<dataset>/compareGO_input_manifest.csv`

This makes the data flow explicit and reproducible while preserving the biological route classification. By default, compareGO consumes only `result_type == GSEA_GO` and `used_for_plot == TRUE`. `GSEA_KEGG` rows are recorded by clusterProfiler for provenance and can be selected intentionally by `config/compareGO_config.yml`, but ORA and custom/NK3R outputs are not compareGO inputs unless future analysis logic explicitly supports them. Legacy mixed mapped folders such as `data/processed/02_id_mapping/mapped/forward/per_file/` are no longer accepted silently; migrate them into `mapped/<dataset>/<direction>/per_file/`.

`04_differential_expression_enrichment/03_biological_program_summary.r` is an additive interpretation layer over these manifest-selected outputs. It writes dataset-scoped program summaries under `results/tables`, source evidence under `results/source_data`, and a heatmap under `results/figures` for `biological_program_summary/<dataset>/`.

## Canonical Output Roots

Active scripts should write reusable downstream data to `data/processed/<module>/<substep>/<dataset_or_global>/`, result tables to `results/tables/<module>/<substep>/<dataset_or_global>/`, figures to `results/figures/<module>/<substep>/<dataset_or_global>/`, source-data tables to `results/source_data/<module>/<substep>/<dataset_or_global>/`, reports to `results/reports/<module>/<substep>/<dataset_or_global>/`, and logs/manifests/session information to `results/logs/<module>/<substep>/<dataset_or_global>/`.

Deprecated technical roots are not valid new output locations for active scripts: `Results/`, `Output/`, `Final/`, `Final2/`, `publication_ready/`, `Plots/`, `Datasets/`, and `results/module_scores/`. These strings may remain only in documentation, migration maps, or explicit read-only fallback warnings.

## Module-Score Metadata Contract

`01_preprocessing/06_merged_metadata_module_score.r` now writes the dataset-scoped module-score metadata workbook to:

`data/processed/01_preprocessing/06_merged_metadata_module_score/<dataset>/sample_metadata_merged_clean_for_module_scores.xlsx`

Its QC tables, summary report, run manifest, and session info are split under:

- `results/tables/01_preprocessing/06_merged_metadata_module_score/<dataset>/`
- `results/reports/01_preprocessing/06_merged_metadata_module_score/<dataset>/`
- `results/logs/01_preprocessing/06_merged_metadata_module_score/<dataset>/`

`R/dataset_inputs.R` resolves module-score metadata in this order: canonical processed dataset path, explicit `PROTEOMICS_MODULE_SCORE_METADATA_FILE`, legacy dataset-scoped fallbacks with a warning, and legacy global fallback only when `PROTEOMICS_ALLOW_GLOBAL_MODULE_SCORE_METADATA=true`.

Checkpoint behavior: if clusterProfiler skips a completed comparison, it reconstructs manifest rows from existing canonical GSEA_GO/GSEA_KEGG source-data tables when present and marks `checkpoint_status = reconstructed_from_checkpoint`. If a checkpoint predates canonical source-data outputs, rerun with `force_rerun: true` to refresh the manifest.

Dry-run commands:

```bash
Rscript 04_differential_expression_enrichment/01_clusterProfiler.r --dry-run
Rscript 04_differential_expression_enrichment/02_compareGO.r --dry-run
```

## PRIDE and Journal Layers

Repository outputs should be separated as:

- canonical processed quantitative matrices: processed proteomics matrices suitable for PRIDE supplementary processed-data upload
- differential abundance result tables: canonical contrast-level statistics suitable for supplements
- enrichment-derived secondary analyses: GO/KEGG/EWCE/WGCNA/network derivatives generally suited for journal supplements, not raw PRIDE redundancy
- figure-only outputs: SVG/PNG/PDF files excluded from PRIDE unless specifically requested by the journal
- source-data tables: exact data behind figures under `results/source_data/`
- exploratory analyses: testing/deprecated outputs excluded from PRIDE and usually excluded from supplements

## Phase 3 Module I/O Harmonization

Phase 3 applied targeted path-only refactors where script behavior could be preserved without data-dependent interpretation. The following active scripts now source `R/paths.R`, write canonical outputs first, include script I/O headers, add dry-run diagnostics, and avoid committed machine-specific paths:

- `01_preprocessing/04_format_metadata.r`
- `01_preprocessing/05_metadata_create.r`
- `05_celltype_enrichment_EWCE/01_EWCE_E9.r`
- `06_modules_WGCNA/06_module_spatial_networks.r`
- `06_modules_WGCNA/02_curated_overlap_programs.r`
- `06_modules_WGCNA/03_score_module_activity.R`
- `06_modules_WGCNA/04_wgcna_de_gsea_overlap.r`
- `07_spatial_networks/02_differential_networks.r`
- `07_spatial_networks/03_bootstrap_network_stability.r`
- `07_spatial_networks/04_bootstrap_differential_network_stability.r`
- `07_spatial_networks/05_bootstrap_differential_network_figures.r`
- `07_spatial_networks/06_chord_diagram.r`
- `08_behavior_physio_coupling/02_network_behavior_coupling.r`
- `01_preprocessing/06_merged_metadata_module_score.r`

The old technical folder roots map as follows:

| Old root/type | Canonical role | Notes |
| --- | --- | --- |
| `Datasets/raw` | `data/raw/` or `data/processed/01_preprocessing/` | Raw vendor/search exports belong in raw; derived matrices belong in processed. |
| `Datasets/mapped`, `mapped`, `unmapped` | `data/processed/02_id_mapping/` | Mapping reports and summaries should also be copied to `results/reports/02_id_mapping/` or `results/tables/02_id_mapping/`. |
| `Results`, `Output`, `Final`, `Final2`, `publication_ready` | `results/{tables,figures,source_data,reports}/<module>/<substep>/` | These are deprecated technical containers unless preserved as compatibility copies. |
| `Plots`, `figures` | `results/figures/<module>/<substep>/` | Figure source tables should be written separately under `results/source_data/`. |
| `core_enrichment` | `results/source_data/04_differential_expression_enrichment/clusterProfiler/` plus manifest rows | compareGO should consume the manifest, not recursively discover this folder. |
| `module_scores` | `data/processed/01_preprocessing/06_merged_metadata_module_score/<dataset>/` for metadata handoff; `data/processed/06_modules_WGCNA/` and `results/tables/06_modules_WGCNA/` for activity scores | `results/module_scores/` is a read-only legacy fallback only; fallback use emits warnings. |
| `WGCNA_modules_long` | `results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/` plus `data/processed/06_modules_WGCNA/01_WGCNA/<dataset>/wgcna_final_model_state.rds` | Color-stable WGCNA module definitions are reusable inputs for dataset-scoped `03_score_module_activity.R` when `PROTEOMICS_MODULE_DEFINITION_SOURCE=WGCNA`; GO labels are display metadata, not eigengene column names. |
| `biological_claims_table` | `results/tables/biological_claims_table.csv` and optional `.xlsx` | Conservative cross-analysis evidence table for manuscript figure planning; claims remain limited by each source analysis. |
| `EWCE_E9_Results` | `data/processed/05_celltype_enrichment_EWCE/EWCE_E9/` and `results/*/05_celltype_enrichment_EWCE/EWCE_E9/` | EWCE caches stay in processed/cache. |
| `network_spatial_relations` | `data/processed/07_spatial_networks/network_spatial_relations/<dataset>/<spatial_unit>/` and `results/*/07_spatial_networks/network_spatial_relations/<dataset>/<spatial_unit>/` | Dataset-aware spatial network producer. `neuron_neuropil` uses region/layer units; `neuron_soma` and `microglia` use region units. |

Scripts left unchanged in Phase 3 are documented because they are central producers, exploratory notebooks-as-scripts, superseded versions, or require data-aware confirmation of file contracts before changing paths:

- `01_preprocessing/01_impute.r`, `01_preprocessing/02_excel_convert.r`, `01_preprocessing/03_gct_extractR.r`
- `02_id_mapping/01_MapThatProt_batch.r`
- `03_qc_exploration/*.r`
- `06_modules_WGCNA/01_WGCNA.r`, `06_modules_WGCNA/06_module_spatial_networks.r`
- `07_spatial_networks/01_network_spatial_relations.r`
- `08_behavior_physio_coupling/01_correlate_proteomics_with_behavior.r`

## Phase 4 Producer Refactors

Phase 4 fixed the main remaining spatial-network dependency gap. `07_spatial_networks/01_network_spatial_relations.r` now writes the canonical downstream object:

`data/processed/07_spatial_networks/network_spatial_relations/<dataset>/<spatial_unit>/network_spatial_relations_objects.rds`

The object records `dataset`, `spatial_unit`, `spatial_col`, `spatial_label_col`, `spatial_unit_matrix`, and preserves the legacy downstream structure (`expression_matrix`, `sample_metadata`, `region_layer_matrix`, `overall_spearman`, `overall_jaccard`, `nodes`, `group_specific`). Tables, figures, source data, network files, logs, and reports are split across dataset- and spatial-unit-scoped canonical folders.

Phase 4 also canonicalized the preprocessing to ID-mapping contract:

- `01_preprocessing/03_gct_extractR.r` reads a GCT from `data/processed/01_preprocessing/protigy_output/<comparison>/` and writes split contrast CSVs to `data/processed/01_preprocessing/gct_extractR/<comparison>/{forward,reverse}/`.
- `02_id_mapping/01_MapThatProt_batch.r` consumes those split CSVs plus `data/external/MOUSE_10090_idmapping.dat` and writes mapped contrast CSVs to `data/processed/02_id_mapping/mapped/<comparison>/<forward|reverse>/per_file/`.
- The default mapping direction is now `forward`, matching the `clusterProfiler` contract. Reverse mapping remains available with `PROTEOMICS_MAP_DIRECTION=reverse`.
- The default comparison family is `neuron_neuropil`, overridable with `PROTEOMICS_COMPARISON`. The recognized default families are `neuron_neuropil`, `neuron_soma`, and `microglia`; extend `PROTEOMICS_ALLOWED_COMPARISONS` for new families.

No `setwd()` remains in these two scripts, and UniProt idmapping download is opt-in via `AUTO_DOWNLOAD_UNIPROT_MAPPING=true`.

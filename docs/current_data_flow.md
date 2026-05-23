# Current Data Flow Audit

This audit was generated before targeted refactoring of the clusterProfiler to compareGO contract. The repository contains R scripts across preprocessing, ID mapping, QC, differential-expression/enrichment, EWCE, WGCNA, spatial networks, behavior coupling, PRIDE helpers, testing and deprecated folders.

## Global Findings

- Many active scripts still contain machine-specific paths such as `S:/...`, `/Users/...`, and `C:/...`; these are documented by `rg` inventory and should be migrated module by module.
- Active scripts write a mix of canonical analysis tables, exploratory workbooks, figure-only exports, RDS caches, network files and logs directly under legacy `Results/` and `Datasets/` style folders.
- `setwd()` is present in active ID mapping and GCT/pathview workflows; path helpers should be used instead, and unavoidable working-directory changes should restore `oldwd`.
- Package auto-installation is present in several scripts. The refactored priority scripts now expose `AUTO_INSTALL_MISSING_PACKAGES` and default to fail-fast behavior.

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

This makes the data flow explicit and reproducible while preserving the biological route classification.

## PRIDE and Journal Layers

Repository outputs should be separated as:

- canonical processed quantitative matrices: processed proteomics matrices suitable for PRIDE supplementary processed-data upload
- differential-expression result tables: canonical contrast-level statistics suitable for supplements
- enrichment-derived secondary analyses: GO/KEGG/EWCE/WGCNA/network derivatives generally suited for journal supplements, not raw PRIDE redundancy
- figure-only outputs: SVG/PNG/PDF files excluded from PRIDE unless specifically requested by the journal
- source-data tables: exact data behind figures under `results/source_data/`
- exploratory analyses: testing/deprecated outputs excluded from PRIDE and usually excluded from supplements

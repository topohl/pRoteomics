# PRIDE Submission Workflow

## Required Raw Files

PRIDE/ProteomeXchange deposition should include original raw mass-spectrometry files, search result exports, search parameters, FASTA/database files, and sample metadata. These files belong in `pride_submission/` locally and should not be committed to GitHub.

## Processed Matrices

Canonical processed quantitative matrices should come from `data/processed/01_preprocessing/` and include normalization, filtering, missingness and imputation metadata. These are PRIDE-linked processed outputs when they represent the stable quantitative data used for downstream analyses.

## Supplementary Tables

Journal supplementary tables should be generated from canonical outputs:

- mapped protein/contrast tables from `data/processed/02_id_mapping/`
- differential-expression results from `data/processed/04_differential_expression_enrichment/`
- enrichment summaries from `results/tables/04_differential_expression_enrichment/`
- EWCE, WGCNA and network summary tables when used in the manuscript

Avoid uploading duplicated routed enrichment exports to PRIDE. Preserve provenance through manifests and include selected, interpretable supplementary tables for the journal.

## Do Not Upload To PRIDE By Default

- duplicated GO/KEGG plot exports
- figure-only SVG/PNG/PDF artifacts
- testing and deprecated outputs
- transient caches/checkpoints unless required to reproduce a published result
- exploratory tables not used in the manuscript

## Provenance Chain

Raw/vendor files -> preprocessing matrices -> mapped protein IDs -> differential-expression contrasts -> clusterProfiler manifest-backed enrichment tables -> compareGO source data/tables/figures -> PRIDE/journal manifests.

Each exported table should include or be paired with metadata for originating script, input files, analysis date, organism database version, filtering thresholds, imputation method, statistical model, FDR method, software/package versions, ontology/database version and simplification status.

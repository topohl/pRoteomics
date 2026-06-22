# PRIDE Submission Workflow

## Required Raw Files

PRIDE/ProteomeXchange deposition should include original raw mass-spectrometry files, search result exports, search parameters, FASTA/database files, and sample metadata. These files belong in `pride_submission/` locally and should not be committed to GitHub.

## Processed Matrices

Canonical processed quantitative matrices should come from `data/processed/01_preprocessing/` and include normalization, filtering, missingness and imputation metadata. These are PRIDE-linked processed outputs when they represent the stable quantitative data used for downstream analyses.

## Dry-Run Validation

Run these before staging a deposition:

```bash
Rscript 09_export_pride_journal/05_make_pride_manifest.R --dry-run
Rscript 09_export_pride_journal/02_make_sample_metadata.R --dry-run
Rscript 09_export_pride_journal/04_make_supplementary_tables.R --dry-run
Rscript 09_export_pride_journal/10_validate_pride_submission.R --dry-run
Rscript 09_export_pride_journal/06_make_methods_summary.R --dry-run
```

Dry-run mode prints pass/fail diagnostics and target output paths without copying supplementary files or writing final export products.

## Supplementary Tables

Journal supplementary tables should be generated from canonical outputs:

- mapped protein/contrast tables from `data/processed/02_id_mapping/`
- differential abundance results from `data/processed/04_differential_expression_enrichment/`
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

Raw/vendor files -> preprocessing matrices -> mapped protein IDs -> differential abundance contrasts -> clusterProfiler manifest-backed enrichment tables -> compareGO source data/tables/figures -> PRIDE/journal manifests.

Each exported table should include or be paired with metadata for originating script, input files, analysis date, organism database version, filtering thresholds, imputation method, statistical model, FDR method, software/package versions, ontology/database version and simplification status.

# PRIDE and journal export (`09_export_pride_journal`)

Active export module for **processed-data PRIDE deposition** and **journal reproducibility** from the **pg_matrix** stage onward.

## Reproducibility boundary

| Included in this module | Not part of repo entry point |
|---|---|
| pg_matrix TSV reference/copy | Raw/vendor MS files (`.raw`, `.wiff`, …) |
| Sample metadata from pg_matrix era | Original search-engine outputs |
| Imputed/split matrices, GCT handoff, ID mapping | FASTA / search parameters |
| Processed DE/enrichment tables (canonical paths) | Full acquisition reproducibility |

Raw/vendor MS, search results, FASTA, and parameter files may be archived **externally** and linked in the manuscript data availability statement. Their absence here is **expected** for `export_level: pg_matrix_onward` and must not be treated as export failure.

## Export levels

- **`pg_matrix_onward`** (default): PRIDE/journal processed quantification + metadata + provenance from `data/raw/pg_matrix/` through `data/processed/` canonical outputs.

Configuration: [`config/export_config.yml`](config/export_config.yml).

## Commands

Full export (recommended):

```bash
Rscript 09_export_pride_journal/RUN_EXPORT.R --dataset all --export-level pg_matrix_onward
Rscript 09_export_pride_journal/RUN_EXPORT.R --dataset microglia --export-level pg_matrix_onward
```

Dry-run (no writes; validation still reports contract status):

```bash
Rscript 09_export_pride_journal/RUN_EXPORT.R --dataset microglia --export-level pg_matrix_onward --dry-run
```

Individual steps:

```bash
Rscript 09_export_pride_journal/02_make_sample_metadata.R --dataset microglia
Rscript 09_export_pride_journal/03_export_processed_pg_matrix_package.R --dataset microglia
Rscript 09_export_pride_journal/04_make_supplementary_tables.R --dataset microglia
Rscript 09_export_pride_journal/01_make_pride_manifest.R --dataset microglia
Rscript 09_export_pride_journal/05_validate_pride_submission.R --export-level pg_matrix_onward
```

Optional broader manifest (non-default):

```bash
Rscript 09_export_pride_journal/01_make_pride_manifest.R --dataset microglia --include-derived-results
Rscript 09_export_pride_journal/01_make_pride_manifest.R --dataset microglia --recursive
```

## Script map

| Script | Role |
|---|---|
| `RUN_EXPORT.R` | Orchestrates the full export pipeline |
| `01_make_pride_manifest.R` | Canonical file manifest (SHA256, non-recursive by default) |
| `02_make_sample_metadata.R` | SDRF-like metadata from pg_matrix-era sample tables |
| `03_export_processed_pg_matrix_package.R` | **PRIDE** processed matrix package + data dictionary |
| `04_make_supplementary_tables.R` | Journal supplementary tables (config globs only) |
| `05_validate_pride_submission.R` | Contract-aware validation report |
| `06_make_methods_summary.R` | Methods, software versions, pipeline steps, limitations |
| `07_make_biological_claims_table.R` | **Manuscript** claims index (not PRIDE-required) |
| `08_export_manuscript_figures.R` | Figure collection for manuscript |
| `09_export_source_data.R` | Source-data table collection for manuscript |

Legacy numeric names (`03_make_supplementary_tables.R`, `04_validate_…`, etc.) forward to the scripts above for backward compatibility.

## Output folders

Generated payload (gitignored): `pride_submission/`

```text
pride_submission/metadata/              sample_metadata.tsv, sdrf_like_metadata.tsv
pride_submission/processed_data/        per-dataset pg_matrix-onward package
pride_submission/supplementary_tables/  journal supplementary staging
pride_submission/methods/               methods_summary.md, known_limitations.md, …
pride_submission/manifests/             pride_file_manifest.tsv
pride_submission/validation/            validation_report.tsv, validation_summary.md
```

Manuscript artifacts also land under `results/manuscript/` and `results/tables/biological_claims_table.*`.

## PRIDE vs journal

- **PRIDE partial/processed deposition**: staged under `pride_submission/` from pg_matrix onward.
- **Journal reproducibility**: supplementary/source-data/claims/figures under `results/manuscript/` and `results/tables/`.
- Do not imply that reviewers can re-run MS acquisition or search from this repository alone.

See also [`docs/PRIDE_EXPORT.md`](../docs/PRIDE_EXPORT.md) and [`docs/REVIEWER_REPRODUCIBILITY.md`](../docs/REVIEWER_REPRODUCIBILITY.md).

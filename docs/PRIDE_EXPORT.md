# PRIDE Export

`09_export_pride_journal/` is the active code module for **processed-data PRIDE deposition** and **journal reproducibility** from the **pg_matrix** stage onward. `09_pride_submission/` is legacy helper code only.

## Scope

- **In scope**: pg_matrix reference, sample metadata, processed matrices/tables, manifests, methods/limitations, validation.
- **Out of scope for repo reproducibility**: raw/vendor MS files, search-engine outputs, FASTA, search parameters (may exist externally).

See [`09_export_pride_journal/README.md`](../09_export_pride_journal/README.md) for commands and output layout.

## Generated payload (gitignored)

```text
pride_submission/metadata/
pride_submission/processed_data/
pride_submission/supplementary_tables/
pride_submission/methods/
pride_submission/manifests/
pride_submission/validation/
```

## Commands

```bash
Rscript 09_export_pride_journal/RUN_EXPORT.R --dataset all --export-level pg_matrix_onward
Rscript 09_export_pride_journal/RUN_EXPORT.R --dataset microglia --export-level pg_matrix_onward --dry-run
Rscript 09_export_pride_journal/05_validate_pride_submission.R --export-level pg_matrix_onward
```

Large raw/vendor files should be uploaded to PRIDE separately and must not be committed to GitHub.

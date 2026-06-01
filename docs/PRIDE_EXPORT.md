# PRIDE Export

`09_export_pride_journal/` is the active code module for PRIDE, journal, and manuscript export. The folder `09_pride_submission/` is retained as legacy helper code only.

Generated deposition payloads belong in the gitignored `pride_submission/` directory:

```text
pride_submission/metadata/
pride_submission/processed_data/
pride_submission/supplementary_tables/
pride_submission/methods/
pride_submission/manifests/
pride_submission/validation/
```

Active export commands:

```bash
Rscript 09_export_pride_journal/01_make_pride_manifest.R --dry-run
Rscript 09_export_pride_journal/02_make_sample_metadata.R --dry-run
Rscript 09_export_pride_journal/03_make_supplementary_tables.R --dry-run
Rscript 09_export_pride_journal/04_validate_pride_submission.R --dry-run
Rscript 09_export_pride_journal/05_make_methods_summary.R --dry-run
Rscript 09_export_pride_journal/06_make_biological_claims_table.R --dry-run
Rscript 09_export_pride_journal/07_export_manuscript_figures.R --dry-run
Rscript 09_export_pride_journal/08_export_source_data.R --dry-run
```

Large raw/vendor files should be uploaded to PRIDE and never committed to GitHub.

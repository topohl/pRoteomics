# Active Script Audit

Run the reproducibility/path audit from the repository root:

```bash
Rscript tools/audit_active_scripts.R
```

The audit scans active R scripts and excludes `90_testing/`, `99_deprecated/`,
and `05_celltype_enrichment_EWCE/90_EWCE_legacy.r`. It fails on committed
machine-specific paths, `setwd()`, unsafe `sink()` use, package installation
calls, `pacman::p_load()`, and duplicate common helper definitions.

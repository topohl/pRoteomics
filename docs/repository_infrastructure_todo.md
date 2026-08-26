# Repository Infrastructure TODO

This repository currently has no committed `LICENSE` or `CITATION.cff`. A
minimal `renv.lock` is present; refresh it from the data-containing analysis
machine after confirming the intended R/Bioconductor package set.

- `LICENSE`: choose an explicit license before adding a file. Do not infer one from the current repository contents.
- `CITATION.cff`: add once the preferred citation, author list, ORCIDs, and release/version metadata are finalized.
- `renv.lock`: refresh from the data-containing analysis machine after confirming the intended R/Bioconductor package set. The current three-package bootstrap and safe refresh workflow are documented in `docs/RENV_LOCK_STATUS.md`; `tools/audit_renv_lock.R` fails closed while scientific sentinel packages are absent.

Lightweight GitHub Actions workflows are present at `.github/workflows/r-smoke-tests.yml` and `.github/workflows/dry-run.yml`. They parse active R scripts and run data-independent tests and dry-runs without requiring private data.

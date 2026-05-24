# Repository Infrastructure TODO

This repository currently has no committed `LICENSE`, `CITATION.cff`, or `renv.lock`.

- `LICENSE`: choose an explicit license before adding a file. Do not infer one from the current repository contents.
- `CITATION.cff`: add once the preferred citation, author list, ORCIDs, and release/version metadata are finalized.
- `renv.lock`: initialize from the data-containing analysis machine after confirming the intended R/Bioconductor package set.

A lightweight GitHub Actions smoke test is present at `.github/workflows/r-smoke-test.yml`. It parses active R scripts and runs data-independent dry-runs without requiring private data.

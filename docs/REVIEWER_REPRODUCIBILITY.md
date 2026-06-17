# Reviewer Reproducibility

Reviewers can inspect the active pipeline, table contracts, and dry-run behavior without private raw data.

## Minimal Setup

Install R and the lightweight packages used by registry parsing and tests:

```r
install.packages(c("yaml", "testthat"), repos = "https://cloud.r-project.org")
```

For the full private-data analysis environment, use:

```r
install.packages("renv", repos = "https://cloud.r-project.org")
renv::restore()
```

## Dry-Run Commands

```bash
Rscript run_dataset_pipeline.R --list-stages
Rscript run_dataset_pipeline.R --dataset neuron_neuropil --stage qc --dry-run
Rscript run_dataset_pipeline.R --dataset microglia --stage enrichment --dry-run
Rscript tests/testthat.R
```

Layer-resolved stages can also be checked for capability behavior:

```bash
Rscript run_dataset_pipeline.R --dataset neuron_neuropil --stage networks --dry-run
Rscript run_dataset_pipeline.R --dataset microglia --stage networks --dry-run
```

The second command should not run microglia network scripts because `pipeline.yml` does not register `microglia` for the layer-resolved network stage. Direct script execution with `PROTEOMICS_DATASET=microglia` should fail clearly through `assert_dataset_capability()`.

## What These Checks Validate

- `pipeline.yml` is parseable and is the active script registry.
- Active scripts listed in `pipeline.yml` exist.
- `RUN_ORDER.md` does not present unregistered scripts as active.
- Dataset capability helpers reject layer-level microglia analyses.
- Table schemas reject invalid datasets, invalid claim grades, invalid biological claim gate statuses, and p/FDR values outside `[0, 1]`.
- The biological claims table treats `claim_grade` as descriptive only. Manuscript eligibility is controlled by `claim_allowed` and `claim_gate_status`; `missing_evidence` explicitly blocks a claim until the named upstream gate evidence is available.

## Expected Limitations Without Private Data

Dry-runs cannot recompute abundance matrices, differential abundance models, enrichment results, WGCNA modules, spatial networks, behavior coupling, manuscript figures, or PRIDE deposition payloads. Full scientific runs require the private source matrices, metadata, raw/vendor files, and local files described in the data availability statement.

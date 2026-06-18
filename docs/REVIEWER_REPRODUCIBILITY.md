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
Rscript run_dataset_pipeline.R --dataset all --stage export --dry-run --strict-inputs
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
- The biological claims table treats `claim_grade` as descriptive only. Manuscript eligibility is controlled by `claim_type`, `claim_allowed`, and `claim_gate_status`; `missing_required` blocks a claim, while `missing_optional`, `not_applicable`, and `diagnostic_only` document claim-type-specific evidence availability.
- Microglia-neuropil independence separates predeclared adjustment families from the exploratory best-Spearman diagnostic. Only predeclared primary/secondary rows in `microglia_neuropil_independence_claim_gate.csv` can support stronger microglia-specific wording; exploratory best-match rows are diagnostic and never claim-enabling on their own.
- Microglia ROI/local-microenvironment wording remains conservative. `neuropil_sensitive` rows block or downgrade microglia-specific/cell-intrinsic interpretation, while inconclusive rows remain diagnostic/contextual unless a predeclared adjustment passes.
- GO-derived enrichment labels preserve `raw_top_GO_term` and `representative_GO_terms`; reviewers should use `safe_program_label`, `term_label_risk`, `label_confidence`, and `claim_use_class` for manuscript wording. Tissue-mismatched GO labels are flagged and conservatively relabeled, not silently discarded.
- Strict input mode is available through `--strict-inputs` or `PROTEOMICS_STRICT_INPUTS=true`. In strict mode, claim-critical and manuscript-facing scripts fail on missing required canonical inputs instead of selecting the newest matching file or a legacy fallback.
- Reviewer audit CSVs under `results/reviewer_audit/` summarize gate evidence availability, allowed/downgraded/disallowed claim counts, GO-label risk, claim-use classes, microglia-neuropil independence/covariate selection, and input resolution provenance by dataset and claim type.

## Input Provenance Audit

Manuscript-facing reruns append rows to:

```text
results/reviewer_audit/input_resolution_audit.csv
```

This table records expected and resolved paths, resolution mode, strict-mode policy, SHA-256 hash, file modification time, and fallback warnings. Non-strict exploratory fallbacks remain auditable; strict reviewer runs make latest-file fallback impossible for claim-critical inputs.

## Expected Limitations Without Private Data

Dry-runs cannot recompute abundance matrices, differential abundance models, enrichment results, WGCNA modules, spatial networks, behavior coupling, manuscript figures, or PRIDE deposition payloads. Full scientific runs require the private source matrices, metadata, raw/vendor files, and local files described in the data availability statement.

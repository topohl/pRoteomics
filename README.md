# pRoteomics

pRoteomics is the **scientific analysis** repository for the Exp9 spatial
proteomics study: hippocampal neuronal neuropil, neuronal soma, and
microglia/PVM-enriched ROI datasets.

It owns scientific inference. Manuscript prose and journal figure rendering are
intentionally maintained in a separate repository — see
[Publication boundary](#publication-boundary).

`pipeline.yml` is the only active machine-readable source of truth for script
order, dataset support, inputs, outputs, and safe-rerun status. Start with
[WORKFLOW.md](WORKFLOW.md) for the plain-language scientific workflow, then use
[RUN_ORDER.md](RUN_ORDER.md) as the detailed command reference.

## Quick Start

```bash
Rscript run_dataset_pipeline.R --list-stages
Rscript run_dataset_pipeline.R --dataset microglia --stage modules_downstream --dry-run
Rscript run_dataset_pipeline.R --dataset all --stage integration --dry-run
```

Run named stages from `pipeline.yml` after private inputs are available:

```bash
Rscript run_dataset_pipeline.R --dataset <dataset> --stage core
Rscript run_dataset_pipeline.R --dataset <dataset> --stage qc
Rscript run_dataset_pipeline.R --dataset <dataset> --stage enrichment
Rscript run_dataset_pipeline.R --dataset <dataset> --stage modules_wgcna
Rscript run_dataset_pipeline.R --dataset <dataset> --stage modules_downstream
Rscript run_dataset_pipeline.R --dataset all --stage integration
```

Valid dataset families are `neuron_neuropil`, `neuron_soma`, and `microglia`.
Microglia data are region-only microglia/PVM-enriched ROI/local
microenvironment proteomics, not purified microglia.

The authoritative control compartment-marker validation is
`analysis/qc/render_compartment_abundance_figures.R`.
It uses reconstructed observed, non-imputed abundance and animal-level
descriptive summaries; it does not estimate purity, cell fractions,
deconvolution, copy number, or total hippocampal abundance.

## Scientific scope

This repository owns, end to end:

- preprocessing and protein/gene identifier mapping;
- quality control, missingness, marker fidelity and confounding checks;
- spatial systems validation, bilateral aggregation and CA2-SLM robustness;
- differential abundance and enrichment;
- WGCNA module and supermodule construction and interpretation;
- cell-type enrichment;
- spatial and differential networks;
- biological integration and behaviour/physiology coupling;
- canonical result tables and publication source-data exports;
- the scientific regression suite.

## Repository Map

```text
pipeline.yml                         active pipeline registry
WORKFLOW.md                          plain-language workflow guide
RUN_ORDER.md                         detailed command reference
run_dataset_pipeline.R               registry-driven launcher

R/                                   reusable scientific function libraries
  paths.R                              repository-root and path bootstrap
  data_contracts/                      dataset, identifier and module contracts
  qc/                                  QC and compartment-marker helpers
  statistics/                          module statistics, WGCNA, evidence
  spatial/                             spatial identity, atlas and robustness
  enrichment/                          enrichment IO, GO themes, EWCE
  networks/                            network construction helpers
  utilities/                           paths, registry, validation, export, plotting

analysis/                            runnable analysis entrypoints
  01_preprocessing/                    preprocessing and identifier mapping
  02_qc/                               QC, marker and confounding checks
  03_spatial_validation/               spatial systems, bilateral, CA2-SLM
  04_differential_abundance/           differential abundance and enrichment
  05_wgcna/                            WGCNA construction and downstream
  06_gsea/                             cell-type enrichment
  07_spatial_networks/                 spatial and differential networks
  08_integration/                      integration and behaviour coupling
  09_publication_exports/              PRIDE and publication source data

config/                              frozen scientific configuration contracts
data/                                raw, metadata and reference inputs
results/                             canonical analysis products
  publication_source_data/             the only manuscript-facing interface
tests/                               private-data-independent tests
audits/                              provenance and robustness audits
tools/                               maintenance, export and reference utilities
archive/                             superseded code kept for provenance
docs/                                reviewer and maintenance documentation
pride_submission/                    generated, gitignored deposition payload
```

Analysis stage identities describe function rather than history. The former
`04_differential_expression_enrichment/` is now
`analysis/differential_abundance/`; manuscript-facing text should describe
these outputs as differential abundance and enrichment results.

## Publication boundary

This repository's responsibility ends at **canonical publication source data
plus a scientific provenance manifest**:

```
results/publication_source_data/<publication_id>/
results/publication_source_data/manifest.csv
```

**This repository has no publication responsibility beyond that point.** The
manuscript repository (`Exp9_manuscript`, local/private) owns panel
composition, typography, panel dimensions, legends, SVG assembly, PDF/PNG/TIFF
export, journal figure naming and layout, and the submission bundle. It renders from a frozen copy of the
bundle above and never reads a live path in this repository, so a future
restructure here cannot break manuscript rendering.

Behaviour data for Figure 1 and Extended Data 5/9 are owned upstream by
[`topohl/MMMSociability`](https://github.com/topohl/MMMSociability) and enter the
manuscript repository through their own frozen import bundle. They are not
re-analysed here.

The test suite enforces this boundary: no manuscript renderer, no manuscript
prose and no journal assembly code may reappear in this repository.

## Documentation

- [Workflow](WORKFLOW.md)
- [Command reference](RUN_ORDER.md)
- [Documentation map](docs/README.md)
- [Analysis entrypoints](docs/ANALYSIS_ENTRYPOINTS.md)
- [Results ownership](docs/RESULTS_OWNERSHIP.md)
- [Repository architecture](docs/REPOSITORY_ARCHITECTURE.md)
- [Restructure plan and equivalence oracle](docs/RESTRUCTURE_PLAN.md)
- [WGCNA workflow](docs/WGCNA_WORKFLOW.md)
- [Datasets](docs/DATASETS.md)
- [Input contracts](docs/INPUT_CONTRACTS.md)
- [Output contracts](docs/OUTPUT_CONTRACTS.md)
- [Microglia ROI interpretation](docs/MICROGLIA_ROI_INTERPRETATION.md)
- [Reviewer reproducibility](docs/REVIEWER_REPRODUCIBILITY.md)
- [PRIDE export](docs/PRIDE_EXPORT.md)

## Active vs Legacy

Active scripts are listed only in `pipeline.yml`. Scripts excluded from the
canonical automated run are tracked in its `legacy` section with an explicit
replacement and status, and are documented in `docs/NAMING_MIGRATION.md`.

`analysis/publication_source_data/` is the active export module. It builds the
canonical publication source-data bundle, its manifest and the hashes of the
scientific source tables, and audits figure outputs for publication readiness.
It does not assemble journal figures, name or lay them out for a journal, or
build a submission package: those belong to
`Exp9_manuscript/tools/package_journal_figures.R`. Superseded
generations live under `archive/` and must not be treated as canonical;
diagnostic and validation-only scripts live under `audits/`.

## Reproducibility

CI and reviewer dry-runs validate the registry, script availability, dataset
capabilities, and table contracts without private raw data. Full scientific runs
require the private source matrices, metadata, raw/vendor files, and local
deposition payload described in the manuscript data availability statement.

## Author

Tobias Pohl

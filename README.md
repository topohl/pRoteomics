# pRoteomics

pRoteomics is the publication-facing spatial proteomics workflow for hippocampal
neuron neuropil, neuron soma, and microglia-enriched ROI datasets.

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
Microglia data are region-only microglia-enriched ROI/local microenvironment
proteomics, not purified microglia.

## Repository Map

```text
pipeline.yml                         active pipeline registry
WORKFLOW.md                          plain-language workflow guide
RUN_ORDER.md                         detailed command reference
run_dataset_pipeline.R               registry-driven launcher
docs/                                reviewer and maintenance documentation
R/                                   shared path, registry, validation, and dataset helpers
01_preprocessing/                    preprocessing handoff
02_id_mapping/                       protein/gene identifier mapping
03_qc_exploration/                   QC, marker, and confounding checks
04_differential_expression_enrichment/ differential abundance and enrichment analyses
05_celltype_enrichment_EWCE/         EWCE analyses
06_modules_WGCNA/                    WGCNA construction and downstream module interpretation
07_spatial_networks/                 spatial network analyses
08_behavior_physio_coupling/         behavior/physiology coupling
10_biological_integration/           manuscript-level evidence synthesis
09_export_pride_journal/             active PRIDE, manuscript, and source-data export module
09_pride_submission/                 legacy helper code only
pride_submission/                    generated, gitignored local deposition payload
tests/                               private-data-independent tests
```

The folder name `04_differential_expression_enrichment/` is retained for
compatibility, but manuscript-facing text should describe these outputs as
differential abundance and enrichment results.

## Documentation

- [Workflow](WORKFLOW.md)
- [Command reference](RUN_ORDER.md)
- [Documentation map](docs/README.md)
- [WGCNA workflow](docs/WGCNA_WORKFLOW.md)
- [Datasets](docs/DATASETS.md)
- [Input contracts](docs/INPUT_CONTRACTS.md)
- [Output contracts](docs/OUTPUT_CONTRACTS.md)
- [Microglia ROI interpretation](docs/MICROGLIA_ROI_INTERPRETATION.md)
- [Reviewer reproducibility](docs/REVIEWER_REPRODUCIBILITY.md)
- [PRIDE export](docs/PRIDE_EXPORT.md)

## Active vs Legacy

Active scripts are listed only in `pipeline.yml`. Legacy filenames are tracked in
the `legacy` section of `pipeline.yml` and documented in
`docs/NAMING_MIGRATION.md`.

`09_export_pride_journal/` is the active export module. `09_pride_submission/`
is retained as legacy helper code and should not be used as the active PRIDE
workflow.

## Reproducibility

CI and reviewer dry-runs validate the registry, script availability, dataset
capabilities, and table contracts without private raw data. Full scientific runs
require the private source matrices, metadata, raw/vendor files, and local
deposition payload described in the manuscript data availability statement.

## Author

Tobias Pohl

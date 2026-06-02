# pRoteomics

pRoteomics is the publication-facing spatial proteomics workflow for hippocampal neuron neuropil, neuron soma, and microglia-enriched ROI datasets. Numbered folders preserve analysis provenance; `pipeline.yml` is the canonical active pipeline registry.

## How To Run

`run_dataset_pipeline.R` executes registry-defined stages from `pipeline.yml`:

```bash
Rscript run_dataset_pipeline.R --list-stages
Rscript run_dataset_pipeline.R --dataset neuron_neuropil --stage qc --dry-run
Rscript run_dataset_pipeline.R --dataset microglia --stage enrichment --dry-run
```

Use `--stage all` or a named stage (`core`, `qc`, `enrichment`, `modules`, `networks`, `behavior`, `export`) after private inputs are available.

## Repository Map

```text
pipeline.yml                         active pipeline registry
run_dataset_pipeline.R               registry-driven launcher
docs/                                reviewer-facing documentation
inst/schemas/                        machine-readable table contracts
R/                                   shared path, registry, validation, and dataset helpers
01_preprocessing/                    preprocessing handoff
02_id_mapping/                       protein/gene identifier mapping
03_qc_exploration/                   QC and confounding checks
04_differential_expression_enrichment/ differential abundance and enrichment analyses
05_celltype_enrichment_EWCE/         EWCE analyses
06_modules_WGCNA/                    WGCNA, curated overlap programs, and source-scoped module activity analyses
07_spatial_networks/                 layer-resolved spatial network analyses
08_behavior_physio_coupling/         behavior/physiology coupling
09_export_pride_journal/             active PRIDE, manuscript, and source-data export module
09_pride_submission/                 legacy helper code only
pride_submission/                    generated, gitignored local deposition payload
tests/                               private-data-independent tests
```

The folder name `04_differential_expression_enrichment/` is retained for compatibility, but manuscript-facing text should describe these outputs as differential abundance and enrichment results.

## Dataset Guardrails

Dataset capabilities are defined in `R/dataset_config.R` and mirrored by `pipeline.yml`.

| dataset | supported interpretation |
|---|---|
| `neuron_neuropil` | region- and layer-resolved neuron neuropil proteomics |
| `neuron_soma` | region- and layer-resolved neuronal soma-enriched proteomics where metadata support it |
| `microglia` | region-only microglia-enriched ROI/local microenvironment proteomics; not purified microglia |

Layer-level network and strata analyses assert the `layer` capability. Microglia remains region-only unless future metadata explicitly extend its contract.

## Active vs Legacy

Active scripts are listed only in `pipeline.yml`. `09_export_pride_journal/` is the active export module. `09_pride_submission/` is retained as legacy helper code and should not be used as the active PRIDE workflow.

Legacy filenames are tracked in the `legacy` section of `pipeline.yml` and documented in `docs/NAMING_MIGRATION.md`; use the replacement names listed in the active registry.

Module activity scoring is source-scoped: `06_modules_WGCNA/03_score_module_activity.R` writes under `results/tables/06_modules_WGCNA/module_score/<dataset>/<module_definition_source>/`. `module_definition_source` can be `wgcna`, `overlap`, or `custom`; curated overlap programs are biological program sets and are distinct from WGCNA modules.

## Reviewer Docs

- [Run order](RUN_ORDER.md)
- [Datasets](docs/DATASETS.md)
- [Input contracts](docs/INPUT_CONTRACTS.md)
- [Output contracts](docs/OUTPUT_CONTRACTS.md)
- [Microglia ROI interpretation](docs/MICROGLIA_ROI_INTERPRETATION.md)
- [Reviewer reproducibility](docs/REVIEWER_REPRODUCIBILITY.md)
- [PRIDE export](docs/PRIDE_EXPORT.md)
- [Naming migration](docs/NAMING_MIGRATION.md)

## Reproducibility

CI and reviewer dry-runs validate the registry, script availability, dataset capabilities, and table contracts without private raw data. Full scientific runs require the private source matrices, metadata, raw/vendor files, and local deposition payload described in the manuscript data availability statement.

## Author

Tobias Pohl

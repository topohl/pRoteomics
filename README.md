# pRoteomics

pRoteomics contains the publication-facing spatial proteomics analysis pipeline for hippocampal neuron neuropil, neuron soma, and microglia-enriched ROI datasets. The repository keeps numbered module folders for provenance, while `pipeline.yml` is the single machine-readable source of truth for active stages and scripts.

## Datasets

- `neuron_neuropil`: region and layer resolved; neuronal neuropil proteomics.
- `neuron_soma`: region and layer resolved where metadata support it; neuronal soma-enriched proteomics.
- `microglia`: region-only microglia-enriched ROI/local microenvironment proteomics; not purified microglia.

## Quick Start

```bash
Rscript run_dataset_pipeline.R --list-stages
Rscript run_dataset_pipeline.R --dataset neuron_neuropil --stage qc --dry-run
Rscript run_dataset_pipeline.R --dataset microglia --stage enrichment --dry-run
```

For reproducible package setup, use:

```r
install.packages("renv")
renv::restore()
```

## Canonical Run

```bash
for ds in neuron_neuropil neuron_soma microglia; do
  Rscript run_dataset_pipeline.R --dataset "$ds" --stage all --dry-run
  Rscript run_dataset_pipeline.R --dataset "$ds" --stage core
  Rscript run_dataset_pipeline.R --dataset "$ds" --stage qc
  Rscript run_dataset_pipeline.R --dataset "$ds" --stage enrichment
  Rscript run_dataset_pipeline.R --dataset "$ds" --stage modules
done
Rscript 09_export_pride_journal/06_make_biological_claims_table.R
Rscript 09_export_pride_journal/07_export_manuscript_figures.R
Rscript 09_export_pride_journal/08_export_source_data.R
```

## Repo Map

```text
01_preprocessing/
02_id_mapping/
03_qc_exploration/
04_differential_expression_enrichment/
05_celltype_enrichment_EWCE/
06_modules_WGCNA/
07_spatial_networks/
08_behavior_physio_coupling/
09_export_pride_journal/      active export module
09_pride_submission/          legacy helper code only
docs/
inst/schemas/
R/
tests/
```

`pride_submission/` is the generated, gitignored local deposition payload. Large raw/vendor files and generated outputs should not be committed.

## Output Map

```text
data/processed/<module>/<dataset>/
results/reports/<module>/<dataset>/
results/tables/<module>/<dataset>/
results/figures/<module>/<dataset>/
results/source_data/<module>/<dataset>/
results/logs/<module>/<dataset>/
results/manuscript/
pride_submission/
```

Active scripts should write `run_manifest.yml`, `sessionInfo.txt`, and input hashes where applicable, using helpers in `R/paths.R`.

## Reproducibility

- `pipeline.yml` defines active stages, required/optional status, supported datasets, inputs, outputs, and manifests.
- `R/dataset_config.R` defines dataset capabilities so layer-level analyses are not accidentally run on region-only microglia ROI data.
- `inst/schemas/` and `R/schema_validation.R` define machine-readable table contracts.
- CI runs dry-run entrypoints and lightweight tests that do not require private raw data.

## Docs

- [Run order](RUN_ORDER.md)
- [Datasets](docs/DATASETS.md)
- [Input contracts](docs/INPUT_CONTRACTS.md)
- [Output contracts](docs/OUTPUT_CONTRACTS.md)
- [Biological interpretation](docs/BIOLOGICAL_INTERPRETATION.md)
- [Microglia ROI interpretation](docs/MICROGLIA_ROI_INTERPRETATION.md)
- [PRIDE export](docs/PRIDE_EXPORT.md)
- [Reviewer reproducibility](docs/REVIEWER_REPRODUCIBILITY.md)
- [Naming migration](docs/NAMING_MIGRATION.md)

## Author

Tobias Pohl

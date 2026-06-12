# Workflow

This guide explains the scientific flow in plain language. The active
machine-readable registry remains `pipeline.yml`; when this document and the
registry disagree, update `pipeline.yml` first and then refresh this guide.

Use the launcher for registry-defined stages:

```bash
Rscript run_dataset_pipeline.R --list-stages
Rscript run_dataset_pipeline.R --dataset <dataset> --stage <stage> --dry-run
Rscript run_dataset_pipeline.R --dataset <dataset> --stage <stage>
```

Valid datasets are `neuron_neuropil`, `neuron_soma`, and `microglia`. Use
`--dataset all` for global stages such as integration.

## 1. Input And Preprocessing

Purpose: convert upstream ProTigy/limma handoff files and metadata into
dataset-scoped contrast tables and module-score metadata.

Run command:

```bash
Rscript run_dataset_pipeline.R --dataset <dataset> --stage core
```

Key outputs to inspect:

- `data/processed/01_preprocessing/gct_extractR/<dataset>/`
- `data/processed/01_preprocessing/06_merged_metadata_module_score/<dataset>/`
- `results/reports/01_preprocessing/06_merged_metadata_module_score/<dataset>/`

Safe to rerun: partly. Metadata merge is safe; GCT extraction recomputes core
handoff state and should be rerun intentionally.

Before moving on: confirm the expected comparison CSVs exist for the dataset and
metadata rows match the samples used downstream.

## 2. ID Mapping

Purpose: map protein identifiers to UniProt/gene identifiers and create
clusterProfiler-ready contrast files.

Run command:

```bash
Rscript run_dataset_pipeline.R --dataset <dataset> --stage core
```

Key outputs to inspect:

- `data/processed/02_id_mapping/mapped/<dataset>/forward/per_file/`
- `results/tables/02_id_mapping/<dataset>/mapped/forward/summaries/`

Safe to rerun: no for casual downstream reruns. It is core state and should be
rerun only when the mapped inputs or mapping rules intentionally change.

Before moving on: check mapping summaries for unexpected unmapped proteins and
confirm mapped forward contrast files exist.

## 3. QC And Reference Markers

Purpose: check sample quality, missingness, replicate consistency, marker
abundance, PCA/confounding structure, variance partitioning, and reference
marker availability.

Run commands:

```bash
Rscript run_dataset_pipeline.R --dataset all --stage qc_global
Rscript run_dataset_pipeline.R --dataset <dataset> --stage qc
```

Key outputs to inspect:

- `config/marker_panels/wgcna_reference_marker_sets.csv`
- `results/tables/03_qc_exploration/05_empirical_roi_marker_discovery/`
- `results/reports/03_qc_exploration/00_dataset_qc_report/<dataset>/`
- `results/reports/03_qc_exploration/07_qc_biology_confounding_report/<dataset>/`

Safe to rerun: yes. These scripts are interpretation/QC layers and do not
change core matrices.

Before moving on: resolve major QC FAIL/WARN findings, confirm marker registry
availability, and verify microglia is treated as ROI/local microenvironment
signal rather than purified microglia.

## 4. Enrichment And Signatures

Purpose: run differential abundance enrichment, compare GO/GSEA outputs, add
neuropil reference annotation, targeted microglia signatures, program summaries,
external signature overlap, and EWCE where inputs exist.

Run command:

```bash
Rscript run_dataset_pipeline.R --dataset <dataset> --stage enrichment
```

Key outputs to inspect:

- `data/processed/04_differential_expression_enrichment/clusterProfiler/<dataset>/clusterProfiler_manifest.csv`
- `data/processed/04_differential_expression_enrichment/compareGO/<dataset>/compareGO_input_manifest.csv`
- `results/tables/04_differential_expression_enrichment/biological_program_summary/<dataset>/program_summary.csv`
- `results/tables/04_differential_expression_enrichment/microglia_targeted_signature_enrichment/microglia/`

Safe to rerun: yes. This stage consumes mapped contrast files and writes
downstream evidence tables.

Before moving on: confirm the clusterProfiler and compareGO manifests are
dataset-scoped and the biological program summary is present.

## 5. WGCNA Construction

Purpose: build dataset-specific WGCNA networks, modules, module definitions, and
supermodule clustering state.

Run command:

```bash
Rscript run_dataset_pipeline.R --dataset <dataset> --stage modules_wgcna
```

Key outputs to inspect:

- `results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/`
- `results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/supermodules/`
- `results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/WGCNA_module_definitions_for_downstream.csv`

Safe to rerun: no for routine downstream work. This stage recomputes WGCNA core
state and can change module assignments.

Before moving on: review module counts, module definitions, preservation/QC
tables, and supermodule annotation audit outputs.

## 6. Module Scoring And Effects

Purpose: score existing module definitions and model module/supermodule group
effects without rebuilding WGCNA.

Run command:

```bash
Rscript run_dataset_pipeline.R --dataset <dataset> --stage modules_downstream
```

Key outputs to inspect:

- `results/tables/06_modules_WGCNA/module_score/<dataset>/`
- `results/tables/06_modules_WGCNA/group_effects/<dataset>/module_group_effects.csv`
- `results/tables/06_modules_WGCNA/group_effects/<dataset>/supermodule_group_effects.csv`

Safe to rerun: yes. These scripts consume existing WGCNA state and do not change
network/module construction.

Before moving on: check model warnings, `evidence_status`, FDR columns, spatial
unit fields, and whether module scores have adequate coverage.

## 7. Biological Annotation

Purpose: annotate WGCNA modules and supermodules with marker, microenvironment,
neuropil-reference, DE/GSEA-overlap, complex/organelle, and robustness evidence.

Run command:

```bash
Rscript run_dataset_pipeline.R --dataset <dataset> --stage modules_downstream
```

Key outputs to inspect:

- `results/tables/06_modules_WGCNA/module_annotation/<dataset>/WGCNA_module_biological_annotation.csv`
- `results/tables/06_modules_WGCNA/module_annotation/<dataset>/WGCNA_supermodule_biological_annotation.csv`
- `results/tables/06_modules_WGCNA/module_complex_architecture/<dataset>/`
- `results/tables/06_modules_WGCNA/module_robustness_sensitivity/<dataset>/`

Safe to rerun: yes. Annotation is downstream of fixed module definitions.

Before moving on: for microglia, confirm ROI/local microenvironment classes are
used conservatively and not interpreted as purified immune activation.

## 8. Interpretable Summaries

Purpose: join WGCNA effects and biological annotation into manuscript-facing
tables, source data, and figures.

Run command:

```bash
Rscript run_dataset_pipeline.R --dataset <dataset> --stage modules_downstream
Rscript 06_modules_WGCNA/07_wgcna_interpretable_summary.r --dataset all
```

Key outputs to inspect:

- `results/tables/06_modules_WGCNA/interpretable_summary/<dataset>/WGCNA_interpretable_summary.xlsx`
- `results/tables/06_modules_WGCNA/interpretable_summary/<dataset>/WGCNA_supermodule_group_effects_interpretable.csv`
- `results/source_data/06_modules_WGCNA/interpretable_summary/<dataset>/`

Safe to rerun: yes. This is a reporting layer over existing effects and
annotations.

Before moving on: compare plot labels to full annotation labels, inspect label
audit tables, and confirm source-data mirrors were written.

## 9. Integration And Export

Purpose: combine enrichment, WGCNA, microenvironment, robustness, spatial,
behavior, and QC evidence into manuscript-level summaries and export packages.

Run commands:

```bash
Rscript run_dataset_pipeline.R --dataset all --stage integration
Rscript run_dataset_pipeline.R --dataset all --stage export
```

Key outputs to inspect:

- `results/tables/10_biological_integration/cross_compartment_program_atlas/global/`
- `results/tables/10_biological_integration/manuscript_program_summary/global/`
- `results/tables/10_biological_integration/evidence_priority_matrix/global/`
- `pride_submission/`

Safe to rerun: yes. These stages synthesize existing outputs and stage export
files.

Before moving on: review integration input-status rows, resolve export
validation FAIL entries, and confirm supplementary/source-data tables point to
the intended upstream evidence.

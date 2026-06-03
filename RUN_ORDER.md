# Recommended run order

This file gives a conservative execution order for the pRoteomics workflows. Individual projects may skip modules that are not relevant.

Canonical machine-readable outputs now belong under `data/processed/<module>/`. Publication tables, figures, source data, logs and reports belong under `results/tables`, `results/figures`, `results/source_data`, `results/logs` and `results/reports`.

## 1. Preprocessing and metadata harmonization

Dataset-scoped runs can be launched by setting `PROTEOMICS_DATASET` once:

```powershell
$env:PROTEOMICS_DATASET = "microglia"
Rscript 01_preprocessing/03_gct_extractR.r --dry-run
Rscript 02_id_mapping/01_MapThatProt_batch.r --dry-run
Rscript 04_differential_expression_enrichment/01_clusterProfiler.r --dry-run
Rscript 04_differential_expression_enrichment/02_compareGO.r --dry-run
```

Or with the launcher:

```powershell
Rscript run_dataset_pipeline.R --dataset microglia --dry-run
Rscript run_dataset_pipeline.R --dataset microglia
Rscript run_dataset_pipeline.R --list-stages
Rscript run_dataset_pipeline.R --dataset neuron_neuropil --stage modules --dry-run
```

Valid dataset families are `neuron_neuropil`, `neuron_soma`, and `microglia`. The shared `R/dataset_config.R` helper resolves `PROTEOMICS_DATASET` first, with backward-compatible fallback to `PROTEOMICS_COMPARISON` and `PROTEOMICS_GCT_COMPARISON`.

The launcher supports staged execution with `--stage core`, `--stage qc`, `--stage enrichment`, `--stage modules`, `--stage networks`, `--stage behavior`, `--stage export`, or `--stage all`. Each run writes a manifest under `results/logs/pipeline/` with dataset, selected stage, script statuses, timestamps, and exit codes. Dry-runs continue across steps so missing private inputs are visible in one report.

```text
01_preprocessing/
```

Expected stable outputs:

```text
data/metadata/sample_metadata_clean.tsv
data/processed/01_preprocessing/protein_matrix_raw.tsv
data/processed/01_preprocessing/protein_matrix_normalized.tsv
data/processed/01_preprocessing/protein_matrix_imputed.tsv
```

The exact matrix names may differ by project. For PRIDE preparation, the most important file is `sample_metadata_clean.tsv`, including a `raw_file_name` column.

## 2. UniProt and identifier mapping

```text
02_id_mapping/
```

Expected outputs:

```text
data/processed/02_id_mapping/protein_id_mapping.tsv
data/processed/02_id_mapping/mapped_protein_matrix.tsv
data/processed/02_id_mapping/mapped/<dataset>/forward/per_file/*.csv
```

The current canonical contrast handoff is:

```powershell
$env:PROTEOMICS_DATASET = "microglia"
Rscript 01_preprocessing/03_gct_extractR.r --dry-run
Rscript 02_id_mapping/01_MapThatProt_batch.r --dry-run
```

`03_gct_extractR.r` writes split contrast CSVs to:

```text
data/processed/01_preprocessing/gct_extractR/<comparison>/forward/
data/processed/01_preprocessing/gct_extractR/<comparison>/reverse/
```

Run the GCT extraction and ID mapping once per biological dataset family, for example `neuron_neuropil`, `neuron_soma`, and `microglia`. `01_MapThatProt_batch.r` defaults to `PROTEOMICS_DATASET=neuron_neuropil` and `PROTEOMICS_MAP_DIRECTION=forward`; set `PROTEOMICS_DATASET` for each family. It writes clusterProfiler-ready files to:

```text
data/processed/02_id_mapping/mapped/<dataset>/forward/per_file/
```

`03_gct_extractR.r` and `01_MapThatProt_batch.r` resume existing table outputs by default. Existing split/mapped CSV tables are skipped, and missing tables are still processed. To intentionally recompute table outputs, use either:

```powershell
$env:PROTEOMICS_RECOMPUTE = "true"
Rscript 01_preprocessing/03_gct_extractR.r
Rscript 02_id_mapping/01_MapThatProt_batch.r
```

or pass `--recompute` / `--force-rerun` to the script. Step-specific env vars are also available: `PROTEOMICS_GCT_RECOMPUTE=true` and `PROTEOMICS_MAPTHATPROT_RECOMPUTE=true`.

Set `PROTEOMICS_MAP_DIRECTION=reverse` only when intentionally producing reverse contrasts.

## 3. QC and exploratory analysis

```text
03_qc_exploration/
```

Canonical dataset-aware QC scripts are:

```bash
Rscript 03_qc_exploration/00_dataset_qc_report.r --dataset neuron_neuropil --dry-run
Rscript 03_qc_exploration/01_sample_qc_quicksearch.r --dataset neuron_neuropil --dry-run
Rscript 03_qc_exploration/02_missingness_diagnostics.r --dataset neuron_neuropil --dry-run
Rscript 03_qc_exploration/03_replicate_consistency.r --dataset neuron_neuropil --dry-run
Rscript 03_qc_exploration/04_marker_rank_abundance_qc.r --dataset neuron_neuropil --dry-run
Rscript 03_qc_exploration/05_pca_confounding_qc.r --dataset neuron_neuropil --dry-run
Rscript 03_qc_exploration/06_variance_partitioning.r --dataset neuron_neuropil --dry-run
Rscript 03_qc_exploration/07_qc_biology_confounding_report.r --dataset neuron_neuropil --dry-run
```

Repeat with `--dataset neuron_soma` and `--dataset microglia`, or set
`PROTEOMICS_DATASET`. Typical outputs include sample-level quicksearch QC,
missingness diagnostics, replicate consistency, marker abundance sanity checks,
PCA/PC metadata confounding summaries, variance partitioning, and a compact
QC-to-biology PASS/WARN/FAIL report. Optional PCA UMAP/t-SNE/clustering outputs
are disabled by default and should be treated as exploratory only.

`00_dataset_qc_report.r` is the recommended canonical report for manuscript
triage. It summarizes missingness, the PCA-only imputation footprint,
sample/protein counts, PCA, batch/sex/group/region/layer structure where
metadata allows, abundance distributions, and robust outlier flags. It writes
tidy tables, an XLSX bundle, SVGs, a Markdown summary, and `run_manifest.yml`
under `results/*/03_qc_exploration/00_dataset_qc_report/<dataset>/`.

See `03_qc_exploration/README.md` for input overrides, output paths, and legacy
script status.

## 4. Differential expression and enrichment

```text
04_differential_expression_enrichment/
```

Typical outputs include differential abundance tables, clusterProfiler/GSEA outputs, compareGO outputs, and publication-style enrichment figures.

The manifest-backed order is:

```r
source("04_differential_expression_enrichment/01_clusterProfiler.r")
source("04_differential_expression_enrichment/02_compareGO.r")
source("04_differential_expression_enrichment/03_biological_program_summary.r")
source("04_differential_expression_enrichment/04_neuropil_contamination_annotation.r")
source("04_differential_expression_enrichment/05_microglia_targeted_signature_enrichment.r")
```

Dry-run validation without launching analyses:

```bash
Rscript 04_differential_expression_enrichment/01_clusterProfiler.r --dry-run
Rscript 04_differential_expression_enrichment/02_compareGO.r --dry-run
Rscript 04_differential_expression_enrichment/03_biological_program_summary.r --dry-run
Rscript 04_differential_expression_enrichment/04_neuropil_contamination_annotation.r --dry-run
Rscript 04_differential_expression_enrichment/05_microglia_targeted_signature_enrichment.r --dataset microglia --dry-run
```

Set `analysis.dataset` in `config/clusterProfiler_config.yml` and `dataset` in `config/compareGO_config.yml` before each dataset family run. Leave mapped paths empty in the configs to use the dataset-aware defaults. `01_clusterProfiler.r` reads mapped contrast CSVs and writes:

```text
data/processed/04_differential_expression_enrichment/clusterProfiler/<dataset>/clusterProfiler_manifest.csv
results/source_data/04_differential_expression_enrichment/clusterProfiler/<dataset>/
results/figures/04_differential_expression_enrichment/clusterProfiler/<dataset>/
results/tables/04_differential_expression_enrichment/clusterProfiler/<dataset>/
results/logs/04_differential_expression_enrichment/clusterProfiler/<dataset>/
results/reports/04_differential_expression_enrichment/clusterProfiler/<dataset>/
```

`02_compareGO.r` requires the manifest `dataset` column, filters manifest rows by `dataset`, and consumes only manifest-selected clusterProfiler outputs plus the mapped/log2FC files listed in that manifest. It writes:

```text
data/processed/04_differential_expression_enrichment/compareGO/<dataset>/compareGO_input_manifest.csv
results/tables/04_differential_expression_enrichment/compareGO/<dataset>/
results/figures/04_differential_expression_enrichment/compareGO/<dataset>/
results/source_data/04_differential_expression_enrichment/compareGO/<dataset>/
results/logs/04_differential_expression_enrichment/compareGO/<dataset>/
results/reports/04_differential_expression_enrichment/compareGO/<dataset>/
```

`03_biological_program_summary.r` maps manifest-selected GO/GSEA terms to the
seven conservative manuscript programs: RNA/RNP processing,
ribosome/translation, mitochondria/OXPHOS/metabolism,
proteostasis/ubiquitin/protein folding, synapse/vesicle organization,
cytoskeleton/motility, and development/patterning. It writes
`program_summary.csv`, `program_summary_wide.csv`,
`program_summary_heatmap_ready.csv`, `program_term_gene_evidence.csv`, and
`program_atlas_heatmap.svg` under
`results/*/04_differential_expression_enrichment/biological_program_summary/<dataset>/`.

### Microglia-enriched ROI interpretation

Microglia ROIs are interpreted as microglia-enriched local microenvironment
samples, not purified microglia. QC indicates that these ROIs can resemble
neuropil samples for neuropil/synaptic proteins while still showing higher
microglia-marker abundance.

The enrichment stage therefore uses neuropil as a reference annotation layer
plus microglia-targeted signature enrichment. It does not subtract neuropil
intensities or logFC values, because dataset families may be separately
normalized and imputed. Microglia data are region-only; when compared with the
region + layer neuropil reference, neuropil layer units are collapsed to parent
regions rather than treated as one-to-one layer equivalents.

For microglia runs, `04_neuropil_contamination_annotation.r` writes optional
GO/GSEA term annotations and `05_microglia_targeted_signature_enrichment.r`
writes method-backed microglia program enrichments under
`results/tables/04_differential_expression_enrichment/`.

## 5. Cell-type enrichment

```text
05_celltype_enrichment_EWCE/
```

Typical outputs include EWCE result tables, cell-type enrichment plots, and measured-proteome-aware enrichment summaries.

Phase 3 canonicalized `01_EWCE_E9.r` so it reads the processed proteomics matrix and sample metadata through `R/paths.R` and writes to:

```text
data/processed/05_celltype_enrichment_EWCE/EWCE_E9/
results/tables/05_celltype_enrichment_EWCE/EWCE_E9/
results/figures/05_celltype_enrichment_EWCE/EWCE_E9/
results/source_data/05_celltype_enrichment_EWCE/EWCE_E9/
results/logs/05_celltype_enrichment_EWCE/EWCE_E9/
results/reports/05_celltype_enrichment_EWCE/EWCE_E9/
```

Dry-run:

```bash
Rscript 05_celltype_enrichment_EWCE/01_EWCE_E9.r --dry-run
```

## 6. Module and WGCNA analyses

```text
06_modules_WGCNA/
```

Typical outputs include module assignments, module scores, module preservation statistics, and trait correlations.

The intended WGCNA downstream order for each dataset is:

```bash
Rscript 06_modules_WGCNA/01_WGCNA.r --dataset <dataset>
Rscript 03_qc_exploration/06_wgcna_marker_trait_export.r --dataset <dataset>
Rscript 06_modules_WGCNA/05_module_supermodule_group_effects.r --dataset <dataset>
Rscript 06_modules_WGCNA/06_annotate_module_microenvironment.r --dataset <dataset>
Rscript 06_modules_WGCNA/07_wgcna_interpretable_summary.r --dataset <dataset>
```

For microglia-enriched ROI interpretation, also run the neuropil reference
annotation before module annotation when the DE/GSEA manifests are available:

```bash
Rscript 04_differential_expression_enrichment/04_neuropil_contamination_annotation.r --dataset microglia
```

These WGCNA downstream scripts answer group effects and interpretation without
changing primary network construction. Microglia ROI WGCNA uses all proteins
detected after the existing filtering/imputation; neuropil/microglia evidence is
used for annotation, reporting, plotting, and sensitivity interpretation only.

Phase 3 canonicalized the safer downstream/helper scripts, and the newer
downstream interpretation layer adds dry-run checks:

```bash
Rscript 06_modules_WGCNA/01_WGCNA.r --dry-run
Rscript 03_qc_exploration/06_wgcna_marker_trait_export.r --dataset microglia --dry-run
Rscript 06_modules_WGCNA/05_module_supermodule_group_effects.r --dataset microglia --dry-run
Rscript 06_modules_WGCNA/06_annotate_module_microenvironment.r --dataset microglia --dry-run
Rscript 06_modules_WGCNA/07_wgcna_interpretable_summary.r --dataset all --dry-run
Rscript 01_preprocessing/06_merged_metadata_module_score.r --dataset microglia --dry-run
Rscript 01_preprocessing/06_merged_metadata_module_score.r --dataset neuron_soma --dry-run
Rscript 01_preprocessing/06_merged_metadata_module_score.r --dataset neuron_neuropil --dry-run
Rscript 06_modules_WGCNA/03_score_module_activity.R --dataset microglia --dry-run
Rscript 06_modules_WGCNA/04_wgcna_de_gsea_overlap.r --dataset microglia --dry-run
Rscript 06_modules_WGCNA/06_module_spatial_networks.r --dataset microglia --module-definition-source wgcna --dry-run
```

`06_modules_WGCNA/03_score_module_activity.R` is the canonical active module-score script. By default it consumes dataset-aware module definitions (source-scoped by dataset and `PROTEOMICS_MODULE_DEFINITION_SOURCE`) and dataset-scoped merged metadata from:

```text
results/module_scores/<dataset>/sample_metadata_merged_clean_for_module_scores.xlsx
```

Legacy global compatibility metadata at `results/module_scores/sample_metadata_merged_clean_for_module_scores.xlsx` should only be used for older scripts.

To regenerate module-score merged metadata per dataset:

```bash
Rscript 01_preprocessing/06_merged_metadata_module_score.r --dataset microglia
Rscript 01_preprocessing/06_merged_metadata_module_score.r --dataset neuron_soma
Rscript 01_preprocessing/06_merged_metadata_module_score.r --dataset neuron_neuropil
```

`01_WGCNA.r` is dataset-aware and can be launched directly or through `run_dataset_pipeline.R`:

```bash
Rscript 06_modules_WGCNA/01_WGCNA.r --dataset neuron_neuropil --dry-run
Rscript 06_modules_WGCNA/01_WGCNA.r --dataset neuron_neuropil
PROTEOMICS_MODULE_DEFINITION_SOURCE=wgcna Rscript 06_modules_WGCNA/91_module_score.r --dataset neuron_neuropil --dry-run
PROTEOMICS_MODULE_DEFINITION_SOURCE=wgcna Rscript 06_modules_WGCNA/91_module_score.r --dataset neuron_neuropil
```

Preferred module-score invocation examples:

```bash
Rscript 06_modules_WGCNA/03_score_module_activity.R --dataset microglia --dry-run
Rscript 06_modules_WGCNA/03_score_module_activity.R --dataset neuron_soma --dry-run
Rscript 06_modules_WGCNA/03_score_module_activity.R --dataset neuron_neuropil --dry-run
```

It stages dataset-scoped WGCNA input workbooks under `data/processed/06_modules_WGCNA/01_WGCNA/<dataset>/inputs/` when upstream imputed matrices and metadata are available. Set `PROTEOMICS_WGCNA_EXPR_XLSX` and `PROTEOMICS_WGCNA_META_XLSX` only when intentionally using custom inputs. It exports stable color-based WGCNA module definitions plus a downstream contract and ranked biological module summary under:

```text
results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/WGCNA_modules_long.xlsx
results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/WGCNA_module_definitions_for_downstream.csv
results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/WGCNA_module_priority_summary.csv
results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/WGCNA_module_contracts.xlsx
```

Additional module evidence exports:

```text
results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/WGCNA_module_evidence_rank.csv
results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/WGCNA_module_preservation_summary.csv
results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/WGCNA_module_GSEA_coregene_overlap.csv
results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/supermodules/wgcna_module_supermodule_annotation.csv
results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/supermodules/wgcna_supermodule_summary.csv
```

The new group-effect outputs are:

```text
results/tables/06_modules_WGCNA/group_effects/<dataset>/module_group_effects.csv
results/tables/06_modules_WGCNA/group_effects/<dataset>/supermodule_group_effects.csv
results/tables/06_modules_WGCNA/group_effects/<dataset>/supermodule_composition.csv
```

The first table answers which modules differ between CON/RES/SUS; the second
answers which supermodules differ between CON/RES/SUS, including spatial unit
and spatial-adjusted/global contrasts where estimable.

The biological interpretation outputs are:

```text
results/tables/06_modules_WGCNA/module_annotation/<dataset>/WGCNA_module_biological_annotation.csv
results/tables/06_modules_WGCNA/module_annotation/<dataset>/WGCNA_supermodule_biological_annotation.csv
results/tables/06_modules_WGCNA/interpretable_summary/<dataset>/WGCNA_supermodule_group_effects_interpretable.csv
results/tables/06_modules_WGCNA/interpretable_summary/<dataset>/WGCNA_module_group_effects_interpretable.csv
```

For microglia, the annotation tables answer whether changed ROI
modules/supermodules are microglia-supported, shared local microenvironment,
neuropil-sensitive, or ambiguous. The wording is intentionally conservative:
microglia ROIs are not treated as purified microglia, and neuropil evidence is
not subtracted or removed from the primary WGCNA.

Optional existing downstream scripts remain useful:

```bash
Rscript 06_modules_WGCNA/03_score_module_activity.R --dataset <dataset>
Rscript 06_modules_WGCNA/04_wgcna_de_gsea_overlap.r --dataset <dataset>
Rscript 06_modules_WGCNA/06_module_spatial_networks.r --dataset <dataset> --module-definition-source wgcna
```

There is a numbering conflict: the existing `06_module_spatial_networks.r` now
shares the `06_` prefix with the new annotation script. It was not renamed here
to avoid breaking existing calls; a later cleanup can move spatial networks to a
non-conflicting number with a compatibility wrapper.

`04_wgcna_de_gsea_overlap.r` is an optional bridge from WGCNA modules to DE/GSEA manifests. It writes `WGCNA_vs_DE_GSEA_overlap.csv/xlsx` under `results/tables/06_modules_WGCNA/04_wgcna_de_gsea_overlap/<dataset>/` and, when overlap inputs are available, appends strongest overlap columns to `WGCNA_module_priority_summary.csv`. Missing DE/GSEA inputs are recorded as skipped status rather than failing the WGCNA run.

## 7. Spatial network analyses

```text
07_spatial_networks/
```

Typical outputs include spatial-unit network edge tables, differential networks, and bootstrap stability summaries.
`neuron_neuropil` uses region/layer units (`spatial_unit=region_layer`), while `neuron_soma` and `microglia` use region-level units (`spatial_unit=region`).

Phase 4 canonicalized the spatial-network producer and Phase 3 canonicalized downstream network scripts:

```bash
Rscript 07_spatial_networks/01_network_spatial_relations.r --dry-run
Rscript 07_spatial_networks/02_differential_networks.r --dry-run
Rscript 07_spatial_networks/03_bootstrap_network_stability.r --dry-run
Rscript 07_spatial_networks/04_bootstrap_differential_network_stability.r --dry-run
Rscript 07_spatial_networks/05_bootstrap_differential_network_figures.r --dry-run
Rscript 07_spatial_networks/06_chord_diagram.r --dry-run
```

These scripts expect the canonical spatial object:

```text
data/processed/07_spatial_networks/network_spatial_relations/<dataset>/<spatial_unit>/network_spatial_relations_objects.rds
```

`01_network_spatial_relations.r` writes dataset- and spatial-unit-scoped objects while preserving `region_layer_matrix` as a compatibility alias. For `neuron_neuropil`, it also writes the old non-scoped RDS path as a legacy compatibility copy.

## 8. Behavior and physiology coupling

```text
08_behavior_physio_coupling/
```

Typical outputs include correlation tables and figures linking proteomic modules or networks with behavioral/physiological phenotypes.

Phase 3 canonicalized `02_network_behavior_coupling.r`. It expects behavior inputs under:

```text
data/external/behavior/
```

Dry-run:

```bash
Rscript 08_behavior_physio_coupling/02_network_behavior_coupling.r --dry-run
```

## 9. PRIDE / ProteomeXchange and journal export

```r
source("09_export_pride_journal/01_make_pride_manifest.R")
source("09_export_pride_journal/02_make_sample_metadata.R")
source("09_export_pride_journal/03_export_processed_pg_matrix_package.R")
source("09_export_pride_journal/04_make_supplementary_tables.R")
source("09_export_pride_journal/05_validate_pride_submission.R")
source("09_export_pride_journal/06_make_methods_summary.R")
```

Dry-run validation:

```bash
Rscript 09_export_pride_journal/01_make_pride_manifest.R --dry-run
Rscript 09_export_pride_journal/02_make_sample_metadata.R --dry-run
Rscript 09_export_pride_journal/03_export_processed_pg_matrix_package.R --dry-run
Rscript 09_export_pride_journal/04_make_supplementary_tables.R --dry-run
Rscript 09_export_pride_journal/05_validate_pride_submission.R --dry-run
Rscript 09_export_pride_journal/06_make_methods_summary.R --dry-run
Rscript 09_export_pride_journal/07_make_biological_claims_table.R --dry-run
```

`07_make_biological_claims_table.R` generates
`results/tables/biological_claims_table.csv` and, when XLSX support is
installed, `results/tables/biological_claims_table.xlsx`. Columns include
dataset, region, layer/cell compartment, contrast, biological program,
direction, key proteins/genes, evidence type, effect size/NES, raw p/FDR,
robustness/stability metric, source file, figure/table target, and an
interpretation note. Missing statistics are left as `NA`.

## Recommended manuscript workflow

```bash
for ds in neuron_neuropil neuron_soma microglia; do
  Rscript run_dataset_pipeline.R --dataset "$ds" --stage all --dry-run
  Rscript run_dataset_pipeline.R --dataset "$ds" --stage core
  Rscript run_dataset_pipeline.R --dataset "$ds" --stage qc
  Rscript run_dataset_pipeline.R --dataset "$ds" --stage enrichment
  Rscript run_dataset_pipeline.R --dataset "$ds" --stage modules
done
Rscript 09_export_pride_journal/07_make_biological_claims_table.R
```

Use `--stage networks`, `--stage behavior`, and `--stage export` after the core
manuscript evidence tables have passed QC and manifest checks.

Before running this step, copy local large files into the gitignored `pride_submission/` folder as needed:

```text
pride_submission/metadata/
pride_submission/processed_data/
pride_submission/supplementary_tables/
pride_submission/methods/
pride_submission/manifests/
pride_submission/validation/
```

This step generates:

```text
pride_submission/manifests/pride_file_manifest.tsv
pride_submission/metadata/sample_metadata.tsv
pride_submission/metadata/sdrf_like_metadata.tsv
pride_submission/validation/validation_report.tsv
pride_submission/methods/methods_summary.md
```

Review `validation_report.tsv` and resolve all `FAIL` entries before upload.

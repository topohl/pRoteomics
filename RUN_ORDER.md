# Recommended Run Order

`pipeline.yml` is the canonical machine-readable registry for active stages, scripts, supported datasets, expected inputs, expected outputs, and generated manifests. This file is the human-readable run guide.

## 1. Inspect The Active Pipeline

```bash
Rscript run_dataset_pipeline.R --list-stages
```

Valid dataset families are:

```text
neuron_neuropil
neuron_soma
microglia
```

## 2. Reviewer Dry-Runs

```bash
Rscript run_dataset_pipeline.R --dataset neuron_neuropil --stage qc --dry-run
Rscript run_dataset_pipeline.R --dataset microglia --stage enrichment --dry-run
```

Dry-runs validate registry resolution, paths, and lightweight contracts without recomputing scientific analyses or requiring private raw data.

## 3. Canonical Manuscript Workflow

```bash
for ds in neuron_neuropil neuron_soma microglia; do
  Rscript run_dataset_pipeline.R --dataset "$ds" --stage all --dry-run
  Rscript run_dataset_pipeline.R --dataset "$ds" --stage core
  Rscript run_dataset_pipeline.R --dataset "$ds" --stage qc
  Rscript run_dataset_pipeline.R --dataset "$ds" --stage enrichment
  Rscript run_dataset_pipeline.R --dataset "$ds" --stage modules
done

Rscript run_dataset_pipeline.R --dataset neuron_neuropil --stage networks
Rscript run_dataset_pipeline.R --dataset neuron_soma --stage networks
Rscript run_dataset_pipeline.R --dataset neuron_neuropil --stage behavior
Rscript run_dataset_pipeline.R --dataset neuron_soma --stage behavior
Rscript run_dataset_pipeline.R --dataset neuron_neuropil --stage export
```

Network and network-behavior stages are layer-resolved and therefore exclude `microglia` under the current dataset contract.

## 4. Stage Summary

| stage | purpose |
|---|---|
| `core` | GCT/export handoff and UniProt/identifier mapping. |
| `qc` | Dataset QC, missingness, PCA/confounding, and QC-to-biology summaries. |
| `enrichment` | Differential abundance handoff, clusterProfiler/GSEA, compareGO, biological program summaries, neuropil reference annotations, targeted microglia ROI signatures, and EWCE. |
| `modules` | WGCNA, WGCNA-to-differential abundance/GSEA overlap, and module activity scoring. |
| `networks` | Region/layer spatial networks, differential networks, stability, and chord figures for layer-capable datasets. |
| `behavior` | Network/proteomics coupling to behavior and physiology where inputs support the requested analysis. |
| `export` | Active PRIDE, manuscript, biological claims, figure, and source-data export. |

## 5. Active Module Entrypoints

The active entrypoints are the scripts listed in `pipeline.yml`. Publication-facing examples include:

```text
04_differential_expression_enrichment/04_neuropil_reference_annotation.r
06_modules_WGCNA/03_score_module_activity.R
09_export_pride_journal/06_make_biological_claims_table.R
09_export_pride_journal/07_export_manuscript_figures.R
09_export_pride_journal/08_export_source_data.R
```

Legacy filenames and historical helper folders are listed under the `legacy` block of `pipeline.yml` and in `docs/NAMING_MIGRATION.md`; they should not be presented as active scripts.

## 6. PRIDE And Manuscript Export

`09_export_pride_journal/` is the active code module. The generated local deposition payload belongs in gitignored `pride_submission/`:

```text
pride_submission/metadata/
pride_submission/processed_data/
pride_submission/supplementary_tables/
pride_submission/methods/
pride_submission/manifests/
pride_submission/validation/
```

The legacy `09_pride_submission/` folder is retained only for historical helper code.

## 7. Microglia Interpretation Guardrails

The `microglia` dataset is region-only microglia-enriched ROI/local microenvironment proteomics. It is not purified microglia. Neuropil is used as a reference annotation layer, not as a subtraction or decontamination step. See `docs/MICROGLIA_ROI_INTERPRETATION.md`.

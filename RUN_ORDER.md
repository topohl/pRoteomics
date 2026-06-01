# Recommended Run Order

`pipeline.yml` is the canonical machine-readable registry for active stages, scripts, supported datasets, expected inputs, expected outputs, and generated manifests. This file is the human-readable run guide.

Valid dataset families:

```text
neuron_neuropil
neuron_soma
microglia
```

## 1. Inspect The Active Pipeline

```bash
Rscript run_dataset_pipeline.R --list-stages
```

The launcher preserves:

```text
--dataset
--stage
--dry-run
--list-stages
```

## 2. Reviewer Dry-Runs

```bash
Rscript run_dataset_pipeline.R --dataset neuron_neuropil --stage qc --dry-run
Rscript run_dataset_pipeline.R --dataset microglia --stage enrichment --dry-run
```

Dry-runs validate paths and contracts without recomputing scientific analyses.

## 3. Canonical Manuscript Workflow

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

Use `--stage networks`, `--stage behavior`, and `--stage export` when those outputs are needed.

## 4. Stage Summary

| stage | purpose |
|---|---|
| `core` | GCT/export handoff and UniProt/identifier mapping. |
| `qc` | Dataset QC, missingness, PCA/confounding, and QC-to-biology summaries. |
| `enrichment` | Differential abundance handoff, clusterProfiler/GSEA, compareGO, biological program summaries, microglia ROI annotations, and EWCE. |
| `modules` | WGCNA, WGCNA-to-DA/GSEA overlap, and module activity scoring. |
| `networks` | Region/layer spatial networks, differential networks, stability, and chord figures. |
| `behavior` | Network/proteomics coupling to behavior and physiology. |
| `export` | PRIDE, manuscript, biological claims, figure, and source-data export. |

## 5. Active Module Entrypoints

The active entrypoints are listed in `pipeline.yml`. Notable publication-facing names:

```text
04_differential_expression_enrichment/04_neuropil_reference_annotation.r
06_modules_WGCNA/03_score_module_activity.R
09_export_pride_journal/07_export_manuscript_figures.R
09_export_pride_journal/08_export_source_data.R
```

Backward-compatible retained names:

```text
04_differential_expression_enrichment/04_neuropil_contamination_annotation.r
06_modules_WGCNA/91_module_score.r
```

See `docs/NAMING_MIGRATION.md` for the full migration table.

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

The legacy `09_pride_submission/` folder is retained only for historical helper code and should not be treated as the active pipeline module.

## 7. Microglia Interpretation Guardrails

The `microglia` dataset is region-only microglia-enriched ROI/local microenvironment proteomics. It is not purified microglia. Neuropil is used as a reference annotation layer, not as a subtraction or decontamination step. See `docs/MICROGLIA_ROI_INTERPRETATION.md`.

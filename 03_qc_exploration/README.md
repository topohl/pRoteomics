# QC and Exploration

This folder contains dataset-aware QC/exploration scripts for the three active
spatial proteomics dataset families:

- `neuron_neuropil`
- `neuron_soma`
- `microglia`

All canonical scripts accept:

```bash
Rscript 03_qc_exploration/<script>.r --dataset neuron_neuropil --dry-run
Rscript 03_qc_exploration/<script>.r --dataset neuron_soma --dry-run
Rscript 03_qc_exploration/<script>.r --dataset microglia --dry-run
```

`PROTEOMICS_DATASET` can be used instead of `--dataset`. Script-specific input
environment variables are still honored, and the shared defaults come from
`R/paths.R`, `R/dataset_config.R`, `R/dataset_inputs.R`, and
`R/qc_exploration_utils.R`.

## Canonical Scripts

Recommended run order:

0. `00_dataset_qc_report.r`
   - Input: processed expression matrix plus sample metadata.
   - Overrides: `PROTEOMICS_DATASET_QC_MATRIX_FILE`,
     `PROTEOMICS_DATASET_QC_METADATA_FILE`.
   - Output: canonical dataset-level QC tables, XLSX bundle, SVG figures, and
     Markdown summary for missingness, imputation footprint, sample/protein
     counts, PCA, metadata structure, abundance distributions, and outlier
     flags under `results/*/03_qc_exploration/00_dataset_qc_report/<dataset>/`.

1. `01_sample_qc_quicksearch.r`
   - Input: annotated quicksearch stats workbook, default
     `data/raw/pg_matrix/quicksearch.stats.annotated.xlsx`.
   - Override: `PROTEOMICS_QC_STATS_FILE`.
   - Output: sample-level QC figures, robust outlier tables, and QC summary
     tables under `results/figures|tables|logs/03_qc_exploration/01_sample_qc_quicksearch/<dataset>/`.

2. `02_missingness_diagnostics.r`
   - Input: processed expression matrix plus sample metadata.
   - Overrides: `PROTEOMICS_MISSINGNESS_MATRIX_FILE`,
     `PROTEOMICS_MISSINGNESS_METADATA_FILE`.
   - Output: missing fraction per sample/protein, missingness by metadata term,
     association tests, SVG plots, XLSX/CSV tables, and a Markdown PASS/WARN
     summary.
   - Prefer a raw or non-imputed matrix. If the resolved file looks imputed, the
     report states that limitation.

3. `03_replicate_consistency.r`
   - Input: processed expression matrix plus metadata with `AnimalID` and, where
     available, `ReplicateGroup`, region/layer/group/plate fields.
   - Overrides: `PROTEOMICS_REPLICATE_MATRIX_FILE`,
     `PROTEOMICS_REPLICATE_METADATA_FILE`.
   - Output: pairwise sample correlations, within-vs-across animal summaries,
     optional animal-aggregated matrix, SVG plot, and PASS/WARN report.

4. `04_marker_rank_abundance_qc.r`
   - Input: processed expression matrix plus metadata.
   - Overrides: `PROTEOMICS_RANK_ABUNDANCE_MATRIX_FILE`,
     `PROTEOMICS_RANK_ABUNDANCE_METADATA_FILE`.
   - Output: rank-abundance tables and SVG plots plus marker abundance score
     tables for neuronal/synaptic/neuropil, nuclear/soma, microglia, astrocyte,
     oligodendrocyte/myelin, vascular, mitochondrial/OXPHOS, ribosomal, and
     RNP/RNA-processing panels.
   - Marker outputs are compartment sanity checks, not cell-type purity claims.

5. `05_pca_confounding_qc.r`
   - Input: processed expression matrix or strict GCT v1.3 plus metadata.
   - Overrides: `PROTEOMICS_PCA_MATRIX_FILE`,
     `PROTEOMICS_PCA_METADATA_FILE`.
   - Output: PCA scores/loadings, scree plots, metadata-colored PCA plots,
     PC-protein correlations, PC1-PC10 metadata association tables, and
     `PCA_confounding_summary.csv`.
   - UMAP/t-SNE/clustering are default off. Enable only for exploratory checks
     with `--run-embeddings` and/or `--run-clustering`.

6. `06_variance_partitioning.r`
   - Input: processed expression matrix plus metadata.
   - Overrides: `PROTEOMICS_VARPART_MATRIX_FILE`,
     `PROTEOMICS_VARPART_METADATA_FILE`.
   - Output: variance fractions per protein, median variance by term, top
     term-driven proteins, metadata canonical correlation/confounding tables,
     and SVG violin/box summaries.
   - The formula adapts to available metadata and includes `(1|AnimalID)` only
     when repeated samples per animal exist.

7. `07_qc_biology_confounding_report.r`
   - Input: processed expression matrix plus metadata; reuses marker score
     outputs when present.
   - Overrides: `PROTEOMICS_CONFOUNDING_MATRIX_FILE`,
     `PROTEOMICS_CONFOUNDING_METADATA_FILE`.
   - Output: compact CSV/XLSX association report, manuscript-friendly SVG
     heatmap, and Markdown PASS/WARN/FAIL summary connecting QC metrics,
     missingness, marker scores, PCs, and metadata.

## Running Per Dataset

```bash
for ds in neuron_neuropil neuron_soma microglia; do
  Rscript 03_qc_exploration/00_dataset_qc_report.r --dataset "$ds" --dry-run
  Rscript 03_qc_exploration/01_sample_qc_quicksearch.r --dataset "$ds" --dry-run
  Rscript 03_qc_exploration/02_missingness_diagnostics.r --dataset "$ds" --dry-run
  Rscript 03_qc_exploration/03_replicate_consistency.r --dataset "$ds" --dry-run
  Rscript 03_qc_exploration/04_marker_rank_abundance_qc.r --dataset "$ds" --dry-run
  Rscript 03_qc_exploration/05_pca_confounding_qc.r --dataset "$ds" --dry-run
  Rscript 03_qc_exploration/06_variance_partitioning.r --dataset "$ds" --dry-run
  Rscript 03_qc_exploration/07_qc_biology_confounding_report.r --dataset "$ds" --dry-run
done
```

Remove `--dry-run` after resolving missing private inputs.

## Legacy Scripts

Legacy collaborator-specific or unsafe scripts live in `legacy/` and are not
part of the canonical run order.

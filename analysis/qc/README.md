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
`R/paths.R`, `R/data_contracts/dataset_config.R`, `R/data_contracts/dataset_inputs.R`, and
`R/qc/qc_exploration_utils.R`.

## Global joint-compartment QC

Run the raw-derived preprocessing product before the global QC consumer:

```powershell
Rscript analysis/preprocessing/build_joint_protigy_input.R --dataset all --dry-run
Rscript analysis/qc/assess_joint_compartment_quality.R --dataset all --dry-run
```

`build_joint_protigy_input.R` uses the unified raw protein-group matrix,
canonical metadata, and mouse mapping to build a balanced shared core (default:
observed in at least 70% of every dataset-by-technical-block), a complete-case
sensitivity matrix, and a broad detected/not-detected union. It applies log2
and one joint sample-wise median normalization, then only the primary core gets
label-blind protein-wise median imputation. Its strict GCT v1.3 export is one
joint ProTIGY input; ProTIGY must not log-transform, normalize, or impute it.

`assess_joint_compartment_quality.R` consumes that prepared bundle for uncorrected PCA,
metadata associations, correlations/clustering, and exploratory fixed-seed
UMAP/t-SNE. It does not run combined DE/WGCNA, batch-correct the primary PCA,
or interpret microglia-enriched ROI observations as purified microglia.

`render_joint_compartment_qc_figures.R` is a rendering-only consumer
of the completed bundle and `00b` tables. It does not call PCA, UMAP or t-SNE,
and it does not overwrite source-data tables. It creates title-free editable
panels plus assembled 183-mm main and Extended Data figures under
`results/figures/03_qc_exploration/00b_joint_compartment_qc/publication_style/global/`.

## Canonical Scripts

Recommended run order:

0. `assess_joint_compartment_quality.R` (after `analysis/preprocessing/build_joint_protigy_input.R`)
   - Input: raw-derived global joint QC bundle.
   - Override: `PROTEOMICS_JOINT_QC_PROCESSED_DIR`.
   - Output: global PCA/UMAP/t-SNE, associations, correlations, sensitivity
     concordance, figures, and Markdown summary under
     `results/*/03_qc_exploration/00b_joint_compartment_qc/global/`.

0a. `render_joint_compartment_qc_figures.R` (after a completed `00b` run)
   - Inputs: the joint preprocessing RDS and existing `00b` PCA, UMAP, t-SNE,
     correlation, clustering, concordance, normalization, and imputation outputs.
   - Override: `PROTEOMICS_JOINT_QC_PUBLICATION_FIGURE_DIR`.
   - Output: a 183-mm five-panel main figure, a 183-mm Extended Data QC figure,
     and editable 89-/183-mm SVG source panels in the publication-style output
     directory. Coordinates, statistics and sample inclusion are unchanged.
   - The missingness panel uses the 200 broad-union proteins with greatest
     max-minus-min detection-rate variation across datasets; samples are ordered
     by dataset, plate, region and layer, and proteins are clustered by binary
     detection pattern.

1. `assess_dataset_quality.R`
   - Input: processed expression matrix plus sample metadata.
   - Overrides: `PROTEOMICS_DATASET_QC_MATRIX_FILE`,
     `PROTEOMICS_DATASET_QC_METADATA_FILE`.
   - Output: canonical dataset-level QC tables, XLSX bundle, SVG figures, and
     Markdown summary for missingness, imputation footprint, sample/protein
     counts, PCA, metadata structure, abundance distributions, and outlier
     flags under `results/*/03_qc_exploration/00_dataset_qc_report/<dataset>/`.

1. `assess_sample_quality.R`
   - Input: annotated quicksearch stats workbook, default
     `data/raw/pg_matrix/quicksearch.stats.annotated.xlsx`.
   - Override: `PROTEOMICS_QC_STATS_FILE`.
   - Output: sample-level QC figures, robust outlier tables, and QC summary
     tables under `results/figures|tables|logs/03_qc_exploration/01_sample_qc_quicksearch/<dataset>/`.

2. `summarize_missingness.R`
   - Input: processed expression matrix plus sample metadata.
   - Overrides: `PROTEOMICS_MISSINGNESS_MATRIX_FILE`,
     `PROTEOMICS_MISSINGNESS_METADATA_FILE`.
   - Output: missing fraction per sample/protein, missingness by metadata term,
     association tests, SVG plots, XLSX/CSV tables, and a Markdown PASS/WARN
     summary.
   - Prefer a raw or non-imputed matrix. If the resolved file looks imputed, the
     report states that limitation.

3. `assess_replicate_consistency.R`
   - Input: processed expression matrix plus metadata with `AnimalID` and, where
     available, `ReplicateGroup`, region/layer/group/plate fields.
   - Overrides: `PROTEOMICS_REPLICATE_MATRIX_FILE`,
     `PROTEOMICS_REPLICATE_METADATA_FILE`.
   - Output: pairwise sample correlations, within-vs-across animal summaries,
     optional animal-aggregated matrix, SVG plot, and PASS/WARN report.

4. `assess_marker_rank_abundance.R`
   - Input: processed expression matrix plus metadata.
   - Overrides: `PROTEOMICS_RANK_ABUNDANCE_MATRIX_FILE`,
     `PROTEOMICS_RANK_ABUNDANCE_METADATA_FILE`.
   - Output: rank-abundance tables and SVG plots plus marker abundance score
     tables for neuronal/synaptic/neuropil, nuclear/soma, microglia, astrocyte,
     oligodendrocyte/myelin, vascular, mitochondrial/OXPHOS, ribosomal, and
     RNP/RNA-processing panels.
   - Marker outputs are compartment sanity checks, not cell-type purity claims.
   - The default remains all groups in the existing output namespace. Use
     `--group CON` (or `PROTEOMICS_RANK_ABUNDANCE_GROUP_FILTER=CON`) to filter
     the matrix and metadata before ranking/scoring and write beneath
     `group_CON` without replacing the all-group QC outputs.

4e. `render_compartment_abundance_figures.R`
   - Authoritative cross-compartment marker-abundance and detection workflow.
     It reconstructs raw-positive log2 values from the validated joint-QC
     bundle, applies the sample offsets estimated on the joint shared core, and
     retains observed missingness. It never uses the imputed primary matrix as
     quantitative marker abundance.
   - Primary analysis uses CON animals only. Technical replicates and layers are
     aggregated within hemisphere and region, regions within hemisphere, and
     valid hemispheres with equal weight within animal. One valid hemisphere is
     permitted in the primary estimate; a separately named sensitivity requires
     both. A hemisphere requires at least three regions, and intended primary
     and strict marker eligibility require 2/3 and 3/3 CON animals,
     respectively.
   - The primary external-marker universe does not require joint-shared-core
     membership. Shared-core-only and robust-z results are explicitly named
     legacy sensitivities. Primary class evidence uses within-protein centered
     log2 abundance, class medians, and subpanel-balanced medians; no
     conventional inferential p-values are reported.
   - Exports v2 mapping/provenance, hemisphere- and animal-level values,
     detection and direction audits, class/subpanel summaries, all 27 ordered
     descriptive animal-bootstrap draws, three leave-one-animal-out cases,
     deterministic display selection, rank-abundance source data, and compact
     editable SVG/PDF figures. Existing v1 files are not overwritten.
   - Results describe neuronal soma, neuronal neuropil, and microglia/PVM-
     enriched ROI marker fidelity and relative abundance—not cell fractions,
     purity, absolute protein quantity, total hippocampal abundance, or
     deconvolution.
   - `04c` remains a processed-matrix availability/WGCNA bridge. Its nonmissing
     values are post-filter/post-imputation availability, not raw detection.
     `04d` cross-compartment sample-level inference is deprecated and disabled
     by default; it is not an active biological-claim workflow.

5. `assess_pca_confounding.R`
   - Input: processed expression matrix or strict GCT v1.3 plus metadata.
   - Overrides: `PROTEOMICS_PCA_MATRIX_FILE`,
     `PROTEOMICS_PCA_METADATA_FILE`.
   - Output: PCA scores/loadings, scree plots, metadata-colored PCA plots,
     PC-protein correlations, PC1-PC10 metadata association tables, and
     `PCA_confounding_summary.csv`.
   - UMAP/t-SNE/clustering are default off. Enable only for exploratory checks
     with `--run-embeddings` and/or `--run-clustering`.

6. `partition_variance.R`
   - Input: processed expression matrix plus metadata.
   - Overrides: `PROTEOMICS_VARPART_MATRIX_FILE`,
     `PROTEOMICS_VARPART_METADATA_FILE`.
   - Output: variance fractions per protein, median variance by term, top
     term-driven proteins, metadata canonical correlation/confounding tables,
     and SVG violin/box summaries.
   - The formula adapts to available metadata and includes `(1|AnimalID)` only
     when repeated samples per animal exist.

7. `summarize_qc_confounding.R`
   - Input: processed expression matrix plus metadata; reuses marker score
     outputs when present.
   - Overrides: `PROTEOMICS_CONFOUNDING_MATRIX_FILE`,
     `PROTEOMICS_CONFOUNDING_METADATA_FILE`.
   - Output: compact CSV/XLSX association report, manuscript-friendly SVG
     heatmap, and Markdown PASS/WARN/FAIL summary connecting QC metrics,
     missingness, marker scores, PCs, and metadata.

## Running Per Dataset

```powershell
foreach ($dataset in @("neuron_neuropil", "neuron_soma", "microglia")) {
  Rscript analysis/qc/assess_dataset_quality.R --dataset $dataset --dry-run
  Rscript analysis/qc/assess_sample_quality.R --dataset $dataset --dry-run
  Rscript analysis/qc/summarize_missingness.R --dataset $dataset --dry-run
  Rscript analysis/qc/assess_replicate_consistency.R --dataset $dataset --dry-run
  Rscript analysis/qc/assess_marker_rank_abundance.R --dataset $dataset --dry-run
  Rscript analysis/qc/summarize_marker_detectability.R --dataset $dataset --dry-run
  Rscript analysis/qc/assess_pca_confounding.R --dataset $dataset --dry-run
  Rscript analysis/qc/partition_variance.R --dataset $dataset --dry-run
  Rscript analysis/qc/summarize_qc_confounding.R --dataset $dataset --dry-run
}
Rscript analysis/qc/render_compartment_abundance_figures.R --dataset global --dry-run
```

Remove `--dry-run` after resolving missing private inputs.

## Legacy Scripts

Legacy collaborator-specific or unsafe scripts live in `legacy/` and are not
part of the canonical run order.

# Canonical analysis entrypoints

For each result that reaches the manuscript: the one script that produces it.
Where more than one script could plausibly be the source, the authoritative one
is named and the others are labelled.

Derived from `pipeline.yml` (105 registered scripts across 12 stages) and the
frozen v9 figure contract. Machine-readable form:
`results/tables/publication_hardening/manuscript_reachability_inventory.csv`.

## How to read this

- **Entrypoint** — run this script to regenerate the result.
- The pipeline runner is `run_dataset_pipeline.R` at the repository root. It is
  the only root-level entrypoint.
- `RUN_ORDER.md` is authoritative for order. This file is authoritative for
  *which script produces a named manuscript result*.
- For ownership of a whole result family — the canonical writer of a stage
  output directory, plus the other scripts that legitimately contribute files
  to it — the authority is `config/results_ownership.csv`, rendered as
  `docs/RESULTS_OWNERSHIP.md`. Both are generated from `pipeline.yml`.
- `docs/MAINTENANCE.md` lists which document answers which question.

## Analyses

| Result | Entrypoint | Stage |
|---|---|---|
| Joint Protigy input, animal-level matrices | `analysis/01_preprocessing/build_joint_protigy_input.R` | joint_qc_preprocessing |
| GCT extraction, merged metadata | `analysis/01_preprocessing/extract_protigy_contrasts.R`, `analysis/01_preprocessing/build_module_score_metadata.R` | core |
| Protein ID mapping | `analysis/01_preprocessing/map_protein_identifiers.R` | core |
| Per-dataset QC report | `analysis/02_qc/assess_dataset_quality.R` | qc |
| Cross-compartment QC | `analysis/02_qc/assess_joint_compartment_quality.R` | qc_global |
| Marker detectability / WGCNA bridge | `analysis/02_qc/summarize_marker_detectability.R` | qc |
| **Canonical GSEA (GO-BP)** | `analysis/04_differential_abundance/run_clusterprofiler_enrichment.R` | enrichment |
| Cross-contrast GO comparison | `analysis/04_differential_abundance/compare_go_enrichment.R` | enrichment |
| **Spatial GO-program atlas** | `analysis/04_differential_abundance/build_go_program_atlas.R` | enrichment |
| SUS/RES spatial DAP atlas | `analysis/04_differential_abundance/build_sus_res_dap_atlas.R` | enrichment |
| **External signature validation** | `analysis/04_differential_abundance/validate_control_spatial_identity.R` | enrichment |
| Cell-type enrichment (EWCE) | `analysis/06_gsea/run_ewce_celltype_enrichment.R` | enrichment |
| **WGCNA modules** | `analysis/05_wgcna/build_wgcna_modules.R` | modules_wgcna |
| WGCNA identity contract | `analysis/05_wgcna/build_module_identity_contract.R` | modules_downstream |
| WGCNA label adjudication | `analysis/05_wgcna/adjudicate_module_labels.R` | networks |
| Spatial networks | `analysis/07_spatial_networks/build_spatial_networks.R` | networks |
| **Animal-level spatial networks** | `analysis/03_spatial_validation/build_animal_spatial_networks.R` | networks |
| Bilateral spatial identity | `analysis/03_spatial_validation/quantify_bilateral_spatial_identity.R` | networks |
| CA2-SLM robustness | `analysis/03_spatial_validation/audit_ca2_slm_robustness.R` | networks |
| Behaviour coupling | `analysis/08_integration/test_behaviour_proteomics_associations.R` | coupling |
| Edge-behaviour coupling | `analysis/08_integration/test_network_behaviour_coupling.R` | coupling |
| **GSEA × WGCNA concordance** | `analysis/08_integration/test_enrichment_module_concordance.R` | integration |
| Cross-compartment program atlas | `analysis/08_integration/build_cross_compartment_atlas.R` | integration |
| Evidence priority matrix | `analysis/08_integration/build_evidence_priority_matrix.R` | integration |
| PRIDE / journal export | `analysis/09_publication_exports/RUN_EXPORT.R` | export (wrapper) |
| Supplementary tables | `analysis/09_publication_exports/build_supplementary_tables.R` | export |
| Source data release | `analysis/09_publication_exports/09_export_source_data.R` | export |

## Publication figures

**These entry points are in the `Exp9_manuscript` repository, not this one.**
The paths below are relative to that repository's root. They are listed here so
that a reader of this repository can find the consumer of each exported bundle;
nothing in this repository runs them.

The current generation is **final_truth_v9**. The files below are its
*renderers*: they produce the panels. What you run to assemble a figure is that
repository's entry point (`figures/figure_02.R`, `figures/figure_03.R`,
`figures/figure_01.R`), which selects and validates the promoted panels. That
repository's `figures/FIGURE_INDEX.md` lists both, per figure.

| Artefact | Canonical renderer (in `Exp9_manuscript`) |
|---|---|
| Figure 2 | `figures/final_truth_v9_figure_02.R` |
| Figure 3 | `figures/final_truth_v9_figure_03.R` |
| Extended Data | `figures/final_truth_v9_extended_data.R` |
| Legends | `figures/final_truth_v9_legends.R` |
| Supplementary tables | `figures/final_truth_v9_supplementary_tables.R` |
| READMEs | `figures/final_truth_v9_readmes.R` |
| Vector audit + story | `figures/final_truth_v9_vector_audit.R` |
| Claim audit | `figures/final_truth_v9_claim_audit.R` |
| Semantic rules + claim chain | `figures/final_truth_v9_semantics.R` |
| Heatmap scale audit | `figures/final_truth_v9_heatmap_scale_audit.R` |

Panel-level ownership is in that repository's
`figures/figure_final_truth_v9_contract.yml`, which names a renderer per panel.
Five of those renderers live in superseded layer files — see §3 of
[REPOSITORY_ARCHITECTURE.md](REPOSITORY_ARCHITECTURE.md).

## Where more than one script could look authoritative

| Result | Authoritative | Not authoritative |
|---|---|---|
| Figures 2 and 3 | not in this repository — figure rendering lives in `Exp9_manuscript`, which reads only the frozen bundles under `results/publication_source_data/` | this repository exports the source data for those figures; it no longer renders them |
| Compartment marker fidelity | `archive/02_qc/04d_compartment_marker_fidelity.r` | `archive/08_integration/01_compartment_fidelity_summary.R` — a cross-dataset poster summary that emits three identically-named files under a different stage directory |
| PCA panels | `figures/final_truth_v9_*` via `nf_pca_compact` | `archive/03_qc_exploration/legacy/06_pcaPlot_Neha.r`, `archive/deprecated/04_pcaPlot_v2.r`, `archive/deprecated/05_pcaPlot_v3.r` |
| clusterProfiler enrichment | `analysis/04_differential_abundance/run_clusterprofiler_enrichment.R` | the eleven `90_testing/clusterProfiler*` scratch scripts and four `99_deprecated/clusterProfiler*` files |

**Figure rendering is no longer part of this repository.** The
`manuscript_candidates` stage and the `figures/` tree moved to `Exp9_manuscript`
during the repository split, together with all seven figure generations. The
superseded layers remain executable there, so the frozen figures stay
reproducible, but nothing in this repository renders a manuscript figure. What
this repository owns is the source data: `analysis/09_publication_exports/` and
the other registered exporters write `results/publication_source_data/`, which
is the only interface the manuscript repository reads.

## Not entrypoints

- `99_audits/` — one-off verification passes, excluded from the registry by
  design. Never a dependency of a publication artefact.
- `90_testing/`, `99_deprecated/`, `*/legacy/` — retained as history. No
  producer layer sources them, and a test enforces that.
- `audits/wgcna/proteomics_wgcna_downstream_audit.R` — a one-off audit that happens to sit at
  the repository root. It is in no registry step.

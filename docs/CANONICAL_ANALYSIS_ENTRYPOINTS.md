# Canonical analysis entrypoints

For each result that reaches the manuscript: the one script that produces it.
Where more than one script could plausibly be the source, the authoritative one
is named and the others are labelled.

Derived from `pipeline.yml` (155 registered steps across 14 stages) and the
frozen v9 figure contract. Machine-readable form:
`results/tables/publication_hardening/manuscript_reachability_inventory.csv`.

## How to read this

- **Entrypoint** — run this script to regenerate the result.
- The pipeline runner is `run_dataset_pipeline.R` at the repository root. It is
  the only root-level entrypoint.
- `RUN_ORDER.md` is authoritative for order; this file is authoritative for
  *which script owns which result*.

## Analyses

| Result | Entrypoint | Stage |
|---|---|---|
| Joint Protigy input, animal-level matrices | `analysis/01_preprocessing/01_prepare_joint_protigy_input.r` | joint_qc_preprocessing |
| GCT extraction, merged metadata | `analysis/01_preprocessing/03_gct_extractR.r`, `analysis/01_preprocessing/06_merged_metadata_module_score.r` | core |
| Protein ID mapping | `analysis/01_preprocessing/01_MapThatProt_batch.r` | core |
| Per-dataset QC report | `analysis/02_qc/00_dataset_qc_report.r` | qc |
| Cross-compartment QC | `analysis/02_qc/00b_joint_compartment_qc.r` | qc_global |
| Marker detectability / WGCNA bridge | `analysis/02_qc/04c_marker_detectability_and_wgcna_bridge.r` | qc |
| **Canonical GSEA (GO-BP)** | `analysis/04_differential_abundance/01_clusterProfiler.r` | enrichment |
| Cross-contrast GO comparison | `analysis/04_differential_abundance/02_compareGO.r` | enrichment |
| **Spatial GO-program atlas** | `analysis/04_differential_abundance/07_compareGO_spatial_program_atlas.r` | enrichment |
| SUS/RES spatial DAP atlas | `analysis/04_differential_abundance/10_sus_res_spatial_dap_atlas.r` | enrichment |
| **External signature validation** | `analysis/04_differential_abundance/09_control_spatial_identity_validation.r` | enrichment |
| Cell-type enrichment (EWCE) | `analysis/06_gsea/01_EWCE_E9.r` | enrichment |
| **WGCNA modules** | `analysis/05_wgcna/01_WGCNA.r` | modules_wgcna |
| WGCNA identity contract | `analysis/05_wgcna/00_wgcna_identity_contract.R` | modules_downstream |
| WGCNA label adjudication | `analysis/05_wgcna/15_wgcna_label_adjudication.R` | networks |
| Spatial networks | `analysis/07_spatial_networks/01_network_spatial_relations.r` | networks |
| **Animal-level spatial networks** | `analysis/03_spatial_validation/13_animal_spatial_networks.R` | networks |
| Bilateral spatial identity | `analysis/03_spatial_validation/02_bilateral_spatial_identity.R` | networks |
| CA2-SLM robustness | `analysis/03_spatial_validation/16_ca2_slm_robustness_audit.R` | networks |
| Behaviour coupling | `analysis/08_integration/01_correlate_proteomics_with_behavior.r` | coupling |
| Edge-behaviour coupling | `analysis/08_integration/02_network_behavior_coupling.r` | coupling |
| **GSEA × WGCNA concordance** | `analysis/08_integration/05_gsea_wgcna_concordance.R` | integration |
| Cross-compartment program atlas | `analysis/08_integration/01_cross_compartment_program_atlas.r` | integration |
| Evidence priority matrix | `analysis/08_integration/03_evidence_priority_matrix.r` | integration |
| PRIDE / journal export | `analysis/09_publication_exports/RUN_EXPORT.R` | export (wrapper) |
| Supplementary tables | `analysis/09_publication_exports/04_make_supplementary_tables.R` | export |
| Source data release | `analysis/09_publication_exports/09_export_source_data.R` | export |

## Publication figures

The current generation is **final_truth_v9**. These are its entrypoints:

| Artefact | Entrypoint |
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

Panel-level ownership is in `figures/figure_final_truth_v9_contract.yml`, which
names a renderer per panel. Five of those renderers live in superseded layer
files — see §3 of [REPOSITORY_ARCHITECTURE.md](REPOSITORY_ARCHITECTURE.md).

## Where more than one script could look authoritative

| Result | Authoritative | Not authoritative |
|---|---|---|
| Figures 2 and 3 | `figures/final_truth_v9_figure_0[23].R` | `figure_0[23].R`, `candidate_*`, `nature_v2_*`, `story_v3/4/5_*`, `spatial_v6_*`, `nature_final_v7_*`, `editorial_v8_*` — all seven earlier generations are still registered steps of the `manuscript_candidates` stage |
| Compartment marker fidelity | `archive/02_qc/04d_compartment_marker_fidelity.r` | `archive/08_integration/01_compartment_fidelity_summary.R` — a cross-dataset poster summary that emits three identically-named files under a different stage directory |
| PCA panels | `figures/final_truth_v9_*` via `nf_pca_compact` | `archive/03_qc_exploration/legacy/06_pcaPlot_Neha.r`, `archive/deprecated/04_pcaPlot_v2.r`, `archive/deprecated/05_pcaPlot_v3.r` |
| clusterProfiler enrichment | `analysis/04_differential_abundance/01_clusterProfiler.r` | the eleven `90_testing/clusterProfiler*` scratch scripts and four `99_deprecated/clusterProfiler*` files |

**The `manuscript_candidates` stage registers 52 scripts covering all seven
figure generations.** Running that stage regenerates every superseded
generation, not just v9. This is deliberate — the superseded layers remain
executable so the frozen figures stay reproducible — but it means "run the
manuscript stage" is not the same as "rebuild the manuscript figures". For the
current figures, run the ten v9 entrypoints above.

## Not entrypoints

- `99_audits/` — one-off verification passes, excluded from the registry by
  design. Never a dependency of a publication artefact.
- `90_testing/`, `99_deprecated/`, `*/legacy/` — retained as history. No
  producer layer sources them, and a test enforces that.
- `audits/wgcna/proteomics_wgcna_downstream_audit.R` — a one-off audit that happens to sit at
  the repository root. It is in no registry step.

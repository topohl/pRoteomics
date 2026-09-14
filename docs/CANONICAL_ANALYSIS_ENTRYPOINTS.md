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
| Joint Protigy input, animal-level matrices | `01_preprocessing/01_prepare_joint_protigy_input.r` | joint_qc_preprocessing |
| GCT extraction, merged metadata | `01_preprocessing/03_gct_extractR.r`, `01_preprocessing/06_merged_metadata_module_score.r` | core |
| Protein ID mapping | `02_id_mapping/01_MapThatProt_batch.r` | core |
| Per-dataset QC report | `03_qc_exploration/00_dataset_qc_report.r` | qc |
| Cross-compartment QC | `03_qc_exploration/00b_joint_compartment_qc.r` | qc_global |
| Marker detectability / WGCNA bridge | `03_qc_exploration/04c_marker_detectability_and_wgcna_bridge.r` | qc |
| **Canonical GSEA (GO-BP)** | `04_differential_expression_enrichment/01_clusterProfiler.r` | enrichment |
| Cross-contrast GO comparison | `04_differential_expression_enrichment/02_compareGO.r` | enrichment |
| **Spatial GO-program atlas** | `04_differential_expression_enrichment/07_compareGO_spatial_program_atlas.r` | enrichment |
| SUS/RES spatial DAP atlas | `04_differential_expression_enrichment/10_sus_res_spatial_dap_atlas.r` | enrichment |
| **External signature validation** | `04_differential_expression_enrichment/09_control_spatial_identity_validation.r` | enrichment |
| Cell-type enrichment (EWCE) | `05_celltype_enrichment_EWCE/01_EWCE_E9.r` | enrichment |
| **WGCNA modules** | `06_modules_WGCNA/01_WGCNA.r` | modules_wgcna |
| WGCNA identity contract | `06_modules_WGCNA/00_wgcna_identity_contract.R` | modules_downstream |
| WGCNA label adjudication | `06_modules_WGCNA/15_wgcna_label_adjudication.R` | networks |
| Spatial networks | `07_spatial_networks/01_network_spatial_relations.r` | networks |
| **Animal-level spatial networks** | `11_spatial_systems/13_animal_spatial_networks.R` | networks |
| Bilateral spatial identity | `11_spatial_systems/02_bilateral_spatial_identity.R` | networks |
| CA2-SLM robustness | `11_spatial_systems/16_ca2_slm_robustness_audit.R` | networks |
| Behaviour coupling | `08_behavior_physio_coupling/01_correlate_proteomics_with_behavior.r` | coupling |
| Edge-behaviour coupling | `08_behavior_physio_coupling/02_network_behavior_coupling.r` | coupling |
| **GSEA × WGCNA concordance** | `10_biological_integration/05_gsea_wgcna_concordance.R` | integration |
| Cross-compartment program atlas | `10_biological_integration/01_cross_compartment_program_atlas.r` | integration |
| Evidence priority matrix | `10_biological_integration/03_evidence_priority_matrix.r` | integration |
| PRIDE / journal export | `09_export_pride_journal/RUN_EXPORT.R` | export (wrapper) |
| Supplementary tables | `09_export_pride_journal/04_make_supplementary_tables.R` | export |
| Source data release | `09_export_pride_journal/09_export_source_data.R` | export |

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
| Compartment marker fidelity | `03_qc_exploration/04d_compartment_marker_fidelity.r` | `08_biological_interpretation/01_compartment_fidelity_summary.R` — a cross-dataset poster summary that emits three identically-named files under a different stage directory |
| PCA panels | `figures/final_truth_v9_*` via `nf_pca_compact` | `03_qc_exploration/legacy/06_pcaPlot_Neha.r`, `99_deprecated/04_pcaPlot_v2.r`, `99_deprecated/05_pcaPlot_v3.r` |
| clusterProfiler enrichment | `04_differential_expression_enrichment/01_clusterProfiler.r` | the eleven `90_testing/clusterProfiler*` scratch scripts and four `99_deprecated/clusterProfiler*` files |

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
- `proteomics_wgcna_downstream_audit.R` — a one-off audit that happens to sit at
  the repository root. It is in no registry step.

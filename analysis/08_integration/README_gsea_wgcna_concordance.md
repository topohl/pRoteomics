# GSEA/GO to WGCNA concordance

`05_gsea_wgcna_concordance.R` is an additive downstream integration. It reads the canonical spatial ranked-GSEA term table, canonical Stage 07 WGCNA inferential handoffs, canonical WGCNA `ProteinGroupID` feature universes, module-preservation summaries, and available microglia Stage 12/13 robustness/readiness summaries. It does not refit WGCNA, recompute modules, change source p-values/FDRs, loosen thresholds, or create a combined p-value.

## Matching

The checked-in `config/gsea_wgcna_theme_module_mapping.csv` maps ontology-aware manuscript themes to explicit same-dataset Stage 07 WGCNA entities. `mapping_role` records whether a relationship is a primary identity or a secondary multifunctional facet. Multiple theme mappings to one module are biological facets of one network entity, never independent WGCNA confirmations. Module annotations establish biological identity only; they are not group-effect evidence.

Local comparisons require the same normalized spatial unit and the existing exploratory local WGCNA endpoint. Recurrent-cross-spatial comparisons require the existing spatial-adjusted global WGCNA endpoint plus a GSEA theme supported in the same direction in at least two spatial units. This is a descriptive recurrent-cross-spatial GSEA representation, not a formal global GSEA test. Mixed-direction themes are retained but classified as unresolved.

## Overlap families

Each Fisher test uses the dataset-specific canonical WGCNA feature universe. The leading edge is restricted to the deterministic representative GO term for that program/evidence row and is mapped to canonical `ProteinGroupID` through the WGCNA universe's member-accession bridge. BH correction is applied separately within `dataset x formal contrast x GSEA spatial unit` across all prespecified matched program-by-module hypotheses.

Genuine multimodule supermodules remain independent WGCNA endpoints in the concordance table. They are not Fisher-tested as unions of member-module proteins because those tests would duplicate correlated module hypotheses. Singleton compatibility aliases never enter the integration. The legacy `04_wgcna_de_gsea_overlap` output is context only.

## Low-n interpretation

The near-zero boundary is the lower quartile of absolute WGCNA estimates within `dataset x entity_level x analysis_tier x effect_scope`. `concordant_imprecise` requires a direction-compatible point estimate above that boundary, a CI extending to the median observed absolute effect in the supported direction, a valid/stable model, no strong available animal-instability flag, and no tier-specific FDR support. The canonical WGCNA FDR threshold remains 0.05.

Adaptive/resilience patterns jointly use all three formal contrasts. They require recurrent-cross-spatial GSEA direction plus at least one stable direction-concordant global WGCNA module endpoint for every contrast used. `SUS - RES` remains the direct group-difference endpoint. `RES_supported_shift_candidate` and `SUS_supported_shift_candidate` do not imply absence or equivalence in the other group. Significance versus non-significance asymmetry is never treated as evidence of a difference.

## Outputs

Tables and mirrored source-data versions are written under:

`results/{tables,source_data}/10_biological_integration/gsea_wgcna_concordance/global/`

The report, input status, protected-source hash audit, and integration audit are under:

`results/reports/10_biological_integration/gsea_wgcna_concordance/global/`

No manuscript figure is produced.

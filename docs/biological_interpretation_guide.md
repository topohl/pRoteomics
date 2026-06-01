# Biological interpretation guide

This guide separates statistical output from biological claims. The project has low-n strata in several places, so manuscript language should emphasize convergence across analyses rather than single-table significance.

## Analysis Layers

**Differential abundance**

Differential abundance contrasts identify proteins with condition-associated abundance changes within a dataset and route. DA supports statements about measured protein abundance in the sampled compartment. It does not by itself establish pathway activity, cell-type origin, or causality.

**GSEA and clusterProfiler**

GSEA tests whether ranked protein-level changes concentrate in annotated gene sets. NES direction is useful for comparing routes, but GO annotations can be redundant and broad. Treat isolated GO terms as exploratory unless they recur across region/layer routes or align with WGCNA/module-score evidence.

**compareGO**

compareGO harmonizes manifest-selected enrichment outputs across comparisons and routes. It supports relative comparison of enriched terms across anatomical strata. It should not be used as an unfiltered recursive search over old enrichment files.

**Biological program summaries**

`03_biological_program_summary.r` maps GO/GSEA terms into broad programs such as mitochondria, RNA/RNP processing, synapse/vesicle, proteostasis, immune/microglia, and cytoskeleton. This is an interpretation layer, not a new statistical test. Use it for figure planning and thematic synthesis.

**WGCNA**

WGCNA modules are protein co-abundance groups. The strongest evidence comes from modules that combine low condition-model FDR, robust adjusted deltas, high hub strength, interpretable GO labels, preservation, and DA/GSEA core-gene overlap. Module labels are display aids; module color/ID remains the stable contract.

**Module scores**

Module scores summarize module abundance patterns per sample or animal. They can support claims that a pre-defined module changes by region/layer/group. They do not prove all member proteins change in the same direction unless directional robustness outputs also support that.

**Spatial networks**

Spatial network edges are molecular similarity relationships among region/layer profiles. They are not anatomical connectivity or electrophysiology. Single edges should be treated as exploratory unless bootstrap/permutation validation supports them.

**Behavior and physiology coupling**

Behavior coupling links animal-level proteomic summaries to behavior/physiology. Correlations with `n < 6` are exploratory. Prefer FDR-adjusted and model-adjusted tables, and describe associations rather than mechanisms.

## Stress/Social Instability Biology Checks

These axes are biologically motivated lenses for social instability stress interpretation. They should not be treated as stress-specific proof on their own. Use them when effects converge across DA/GSEA, biological program summaries, WGCNA/module scores, spatial network changes, behavior/physiology coupling, and QC/confounding checks.

**HPA/glucocorticoid signaling**

Look for measured or indirectly supported changes in glucocorticoid-response and stress-regulatory proteins, including chaperone/co-chaperone systems and receptor-associated targets where detected. Interpret cautiously because many canonical HPA markers may not be robustly quantified in spatial proteomics.

**Neuroimmune and microglial state**

Check complement, lysosome/phagosome, antigen presentation, interferon/inflammatory, and microglial homeostatic/activation-associated programs. Keep marker abundance, microglial compartment enrichment, and immune activation claims separate; marker shifts alone do not prove activation state.

**Synaptic remodeling and plasticity**

Evaluate vesicle cycling, postsynaptic density, cytoskeleton, adhesion/scaffold, and neurite/synaptic maintenance proteins. Social stress effects may appear as region/layer-specific remodeling rather than a simple global synapse increase or decrease.

**Mitochondria and energy metabolism**

Track OXPHOS, TCA cycle, mitochondrial ribosome/import, mitochondrial dynamics, and energy-buffering proteins. These programs are stress-relevant but can also reflect compartment composition, missingness, or technical PCs.

**Proteostasis and cellular stress response**

Consider heat shock/chaperone, ubiquitin-proteasome, autophagy, lysosome, ER stress, and protein-folding programs. These are useful bridges between cellular stress and protein abundance changes, especially when supported by module-level evidence.

**RNA/RNP processing and translation**

Check ribosomal, RNA-binding, splicing, translation initiation/elongation, and RNP/granule-related proteins. Treat these as regulatory-state hypotheses unless directional robustness and module support are strong.

**Myelin and oligodendrocyte-associated programs**

Evaluate myelin, lipid metabolism, and oligodendrocyte-associated proteins, especially when effects appear in neuropil or spatial network analyses. Distinguish biological remodeling from compartment contamination using marker QC and metadata confounding reports.

**ECM, vascular, and barrier-associated programs**

Track extracellular matrix, basement membrane, endothelial, pericyte, and vascular-associated proteins. These can be stress-relevant but have high sensitivity to sampling, region/layer balance, and vascular contamination.

**Oxidative stress and redox balance**

Consider peroxiredoxin, thioredoxin, glutathione, antioxidant, and reactive oxygen handling proteins. Interpret redox shifts alongside mitochondrial, proteostasis, and QC metrics because handling and abundance-depth effects can mimic biology.

**Behavior and physiology anchors**

When available, connect protein/module effects to animal-level measures such as susceptibility/resilience labels, body weight, adrenal/thymus measures, corticosterone, and behavioral readouts. These should be animal-level validations, not sample-level pseudoreplication.

Before converting any stress-relevant axis into a manuscript claim, check whether the same signal is associated with missingness, plate/batch, region/layer imbalance, animal identity, marker contamination, or early PCs.

## Confirmatory vs Exploratory

Confirmatory-leaning evidence:

- FDR-controlled DA/GSEA or WGCNA effects.
- WGCNA modules with preservation and DA/GSEA overlap support.
- Spatial edges with bootstrap stability and FDR/permutation support.
- Behavior models with adequate n, FDR support, and sensitivity consistency.

Exploratory evidence:

- Single GO terms without recurrence.
- Regex-based biological program membership.
- Unpreserved or low-n WGCNA modules.
- Full spatial edge tables before bootstrap filtering.
- Behavior correlations with `n < 6`.

## Manuscript Figure Roadmap

1. Dataset and contract overview: sample counts, dataset filters, input manifests.
2. DA/GSEA atlas: compareGO heatmaps plus biological program summary.
3. WGCNA module panel: module evidence rank, module labels, top hubs, preservation.
4. Module score validation: per-sample/per-animal module score shifts and directional robustness.
5. Spatial network panel: only bootstrap-supported edges in main figures, full edge tables in supplement.
6. Behavior coupling panel: FDR-adjusted/model-supported associations, with exploratory low-n findings clearly labeled.

Use `results/tables/biological_claims_table.csv` as the cross-analysis checklist before drafting claims.

# Evidence-dependence registry for the spatial systems layer.
#
# WHY THIS EXISTS
#   Most of the evidence streams in this study are computed from the SAME nine
#   animals and often the same proteomic measurements. Treating them as
#   independent corroboration would double-count: a WGCNA module, an empirical
#   compartment contrast and a SUS-RES contrast can all "agree" simply because
#   they are three views of one matrix.
#
#   This registry states, for every stream, what its biological unit is, whether
#   it uses phenotype, how hemispheres are handled, and - crucially - what may
#   and may not be concluded from it. It is prospective: future streams are
#   listed before they are built so their independence class is fixed in advance
#   rather than argued afterwards.
#
# INDEPENDENCE CLASSES
#   primary_phenotype_evidence   the phenotype contrast itself
#   phenotype_independent_context   same animals, no phenotype used
#   internal_reproducibility     paired sides within an animal; NOT replication
#   contextual_same_data         same proteomics, different view
#   external_annotation          an outside reference, independent of this cohort

sps_independence_classes <- function() {
  c("primary_phenotype_evidence", "phenotype_independent_context",
    "internal_reproducibility", "contextual_same_data", "external_annotation")
}

sps_evidence_dependence_registry <- function() {
  rows <- list(
    c("con_spatial_identity",
      "analysis/04_differential_abundance/validate_control_spatial_identity.R",
      "neuron_soma; neuron_neuropil", "AnimalID", "CON only",
      "modelled as a fixed nuisance covariate; sides not separated", "yes", "none",
      "phenotype_independent_context",
      "Defines where a protein normally sits anatomically. Phenotype-blind, so it may be used as the baseline against which a SUS-RES effect is later positioned.",
      "Must not be described as replicated across animals beyond n=3 CON, and must not be used to support any group difference."),

    c("bilateral_spatial_validation",
      "analysis/03_spatial_validation/quantify_bilateral_spatial_identity.R",
      "neuron_soma; neuron_neuropil; microglia", "paired sides within AnimalID", "CON only",
      "left-only and right-only fits compared against the bilateral fit", "yes", "none",
      "internal_reproducibility",
      "Shows whether an anatomical effect reproduces in size and direction on the opposite side of the same brain. Supports measurement reliability.",
      "NOT independent biological replication. Two hemispheres of one animal are one animal; n does not double, and agreement between sides is not evidence that an effect generalises across animals."),

    c("empirical_compartment_identity",
      "analysis/02_qc/discover_empirical_roi_markers.R",
      "cross-dataset", "AnimalID", "none (dataset contrast, group-adjusted)",
      "averaged within animal before modelling", "yes", "none",
      "phenotype_independent_context",
      "Experiment-derived compartment affinity: whether a protein is enriched in the microglia-enriched ROI relative to a neuronal compartment.",
      "Not a cell-proportion estimate and not purified-cell evidence. A microglia-enriched ROI is not purified microglia."),

    c("cross_hemisphere_compartment_validation",
      "analysis/03_spatial_validation/quantify_empirical_compartments.R",
      "cross-dataset", "paired sides within AnimalID", "adjusted for StressGroup",
      "left discovery evaluated on right and vice versa", "yes", "none",
      "internal_reproducibility",
      "Shows a compartment marker keeps its direction on the opposite side of the same animals.",
      "Not independent replication and not external validation. The same animals contribute both sides."),

    c("wgcna_module_identity",
      "analysis/05_wgcna/build_wgcna_modules.R (frozen membership)",
      "neuron_soma; neuron_neuropil; microglia", "protein (network topology)", "none",
      "samples enter the network before any hemisphere aggregation", "yes", "none",
      "contextual_same_data",
      "Groups co-varying proteins into modules. Structure only.",
      "Module membership is not independent evidence for a claim about the proteins in it; it is the same measurements re-expressed."),

    c("wgcna_bilateral_validation",
      "analysis/03_spatial_validation/quantify_module_bilateral_identity.R",
      "neuron_soma; neuron_neuropil; microglia", "paired sides within AnimalID", "none",
      "consumes Stage-05 hemisphere values; sides kept apart", "yes", "none",
      "internal_reproducibility",
      "Separates two distinct claims: that a module's absolute level reproduces across sides, and that its spatial pattern does. A module may show a constant side offset yet a preserved pattern.",
      "NOT independent biological replication: both sides come from the same animal, so n does not increase. It also does not validate the module definition, and poor agreement is not automatically technical failure - it may be real hemispheric asymmetry."),

    c("ewce_differential",
      "analysis/06_gsea/run_ewce_celltype_enrichment.R (Differential arm)",
      "neuron_soma; neuron_neuropil; microglia", "AnimalID", "SUS/RES/CON contrasts",
      "averaged within animal before modelling", "yes",
      "EWCE specificity reference (ewceData CTD)",
      "primary_phenotype_evidence",
      "Cell-type enrichment of phenotype-derived signatures.",
      "Its FDR family must never be shared with phenotype-blind annotation; doing so lets the phenotype arm determine a baseline result's significance."),

    c("ewce_module_annotation",
      "analysis/03_spatial_validation/annotate_module_celltypes.R",
      "neuron_soma; neuron_neuropil; microglia", "gene set (module membership)", "none",
      "not applicable; membership carries no hemisphere", "yes",
      "EWCE specificity reference (ewceData CTD)",
      "external_annotation",
      "External cell-type expression enrichment of a module's gene set against an outside reference, with the measured proteome as background. Phenotype-blind by construction.",
      "External expression enrichment is NOT empirical compartment affinity and not a cell-proportion claim. It says a module's genes are preferentially expressed in a reference cell type, not that the tissue contains that cell type."),

    c("sus_res_da",
      "04_differential_expression_enrichment (ProTIGY/limma DA)",
      "neuron_soma; neuron_neuropil; microglia", "AnimalID", "SUS vs RES vs CON",
      "equal-weight L/R mean before DA", "yes", "none",
      "primary_phenotype_evidence",
      "The primary phenotype contrast.",
      "With 3 animals per group, absence of significance is weak evidence of absence; and an exact SUS-vs-RES label permutation has only 20 assignments, so no two-sided exact p below 0.10 is attainable."),

    c("sus_res_gsea",
      "04_differential_expression_enrichment (preranked GSEA)",
      "neuron_soma; neuron_neuropil; microglia", "AnimalID", "SUS vs RES vs CON",
      "inherited from the DA input", "yes",
      "GO / MSigDB gene sets",
      "primary_phenotype_evidence",
      "Pathway-level summary of the same phenotype contrast.",
      "Not independent of the DA it is computed from; the two must not be counted as two lines of evidence for one claim."),

    c("animal_level_spatial_network",
      "11_spatial_systems (future pass)",
      "neuron_soma; neuron_neuropil; microglia", "AnimalID", "none for construction",
      "one network per side plus an equal-weight bilateral network per animal",
      "yes", "none",
      "contextual_same_data",
      "Similarity between the relative proteomic profiles of two anatomical units within one animal.",
      "Not neural connectivity, not anatomical connectivity, not protein coexpression, and not interregional communication. Edge inference is limited by n=3 per group.")
  )
  out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(out) <- c("evidence_id", "evidence_source", "dataset_scope",
                  "biological_unit", "phenotype_used", "hemisphere_handling",
                  "same_animals_as_primary", "external_reference",
                  "independence_class", "allowed_interpretation",
                  "prohibited_interpretation")
  stopifnot(all(out$independence_class %in% sps_independence_classes()))
  stopifnot(!anyDuplicated(out$evidence_id))
  out
}

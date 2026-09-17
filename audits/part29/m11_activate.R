#!/usr/bin/env Rscript

# Finalization section 18: activate the adjudicated neuropil m11 annotation.
#
# The label is a COMPOSITIONAL claim about which proteins the module contains,
# which is what a co-abundance module can support. It is NOT a cell-type
# identity: these are enriched-ROI proteomics, not sorted cells, so the external
# oligodendrocyte affinity is recorded as supporting CONTEXT in its own field
# and never inside the module name.

setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")
AUD <- file.path("results", "tables", "publication_audits",
                 "upstream_enrichment_v10")
REGDIR <- file.path("config", "wgcna_labels")
OUT <- file.path(REGDIR, "neuron_neuropil.csv")

tmpl <- utils::read.csv(file.path(REGDIR, "microglia.csv"),
                        stringsAsFactors = FALSE)
row <- tmpl[1, , drop = FALSE]
row$dataset <- "neuron_neuropil"
row$level <- "module"
row$entity_id <- "WGCNA_m11"
row$reviewed_biological_label <- "Enriched for myelin-associated proteins"
row$reviewed_short_label <- "Myelin-associated proteins"
row$subcellular_context <- "myelin sheath; axon ensheathment machinery"
row$roi_context <- paste0(
  "Neuropil ROI co-abundance module. External oligodendrocyte affinity (EWCE ",
  "z = 32.2-34.1 across three protein-set scopes, FDR 8.08e-04; reference ",
  "oligodendrocyte panel FDR 7.86e-35) is SUPPORTING CELL CONTEXT only and ",
  "does not establish a cell-intrinsic oligodendrocyte origin. Do not call ",
  "this the oligodendrocyte module.")
row$label_confidence <- "high"
row$manual_review_required <- FALSE
row$aggregation_evidence_class <- "not_applicable_module"
row$structural_status <- "module"
row$rationale <- paste0(
  "Convergent evidence across three internal layers. Enrichment: CC myelin ",
  "sheath GO:0043209 FDR 6.34e-11 (22/86), BP ensheathment of neurons / axon ",
  "ensheathment FDR 1.55e-09, BP myelination FDR 1.24e-08; 13 of 21 ",
  "significant terms carry myelin/ensheathment names and the remaining 8 are ",
  "CC cell-polarity terms built from the same proteins, so there is no ",
  "competing second block. Hubs: 12 of the top 13 are canonical myelin ",
  "proteins (CNP 0.978, MAG 0.974, SIRT2 0.967, ERMN 0.967, ENPP6 0.967, ",
  "MOG 0.963, BCAS1 0.963, PLP1 0.962, CLDN11 0.960, NDRG1 0.959, ",
  "OPALIN 0.950, MBP 0.945); 8 of the top 10 sit inside the leading term gene ",
  "sets. Spatial: peak ca3_so, tau 0.724. The historical automatic label ",
  "'synaptic/cytoskeletal trafficking' is contradicted by all of this and is ",
  "retired. The activated wording is compositional and keeps the module ID ",
  "visible: 'm11, enriched for myelin-associated proteins'.")
row$reviewer <- "Part-29 registry adjudication"
row$proposal_prepared_by <- "Part-28 evidence matrix; Part-29 candidate-form adjudication"
row$review_date <- "2026-09-14"
row$adjudication_status <- "reviewed"

utils::write.csv(row, OUT, row.names = FALSE)
cat("[ok] wrote", OUT, "\n")

prev <- utils::read.csv(file.path(AUD, "wgcna_registry_final_review.csv"),
                        stringsAsFactors = FALSE)
prev$promoted_publication_form <- NA_character_
prev$promotion_status <- "NOT_ACTIVATED"
k <- prev$dataset == "neuron_neuropil" & prev$module_id == "m11"
prev$promoted_publication_form[k] <- "m11, enriched for myelin-associated proteins"
prev$promotion_status[k] <- "ACTIVATED_in_config/wgcna_labels/neuron_neuropil.csv"
prev$prohibited_form <- "the oligodendrocyte module; any cell-intrinsic oligodendrocyte claim"
utils::write.csv(prev, file.path(AUD,
  "wgcna_registry_final_review_promoted.csv"), row.names = FALSE)
cat("[ok] promotion recorded for", sum(k), "module(s);",
    sum(prev$promotion_status == "NOT_ACTIVATED"), "remain unactivated\n")

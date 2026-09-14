#!/usr/bin/env Rscript

# Part-28: what each WGCNA module's label is actually supported by.
#
# ANNOTATION AUDIT ONLY. No module is recomputed, no label registry is written,
# no proposed label is activated. This script reads the four evidence layers the
# canonical workflow already stored - GO enrichment, hub (high-kME) composition,
# spatial profile and external cell-type affinity - and states, per module, what
# the evidence supports and whether the label currently permitted for the
# manuscript agrees with it.
#
# THE NAMING PRINCIPLE (Part-28 section 17). A module is defined by co-abundance
# topology, not by its annotation, so the permitted manuscript form is
# "m02, enriched for mitochondrial respiration proteins" and never "the
# mitochondrial module". The module ID always stays visible.
#
# EXTERNAL AFFINITY IS NOT IDENTITY (section 18). An EWCE or reference-panel
# result alone may describe cell-type CONTEXT; it may never define what the
# module IS. A module whose only evidence layer is external is therefore
# classified CELL_CONTEXT_ONLY, not functional.

source(file.path("R", "paths.R"))
source(repo_path("R", "null_coalescing.R"))
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "nature_v2_figure_utils.R"))
Sys.setenv(PROTEOMICS_SCRIPT_ID = "figures/final_truth_v9_wgcna_annotation_audit.R")

args <- commandArgs(trailingOnly = TRUE)
if ("--dry-run" %in% args || is_dry_run()) {
  message("[DRY-RUN] final_truth_v9 WGCNA annotation evidence audit")
  quit(save = "no", status = 0L)
}

OUT <- path_results("tables", "manuscript_candidates", "final_truth_v9", "audit")
REP <- path_results("reports", "manuscript_candidates", "final_truth_v9")
dir_create(OUT); dir_create(REP)

FDR <- 0.05
DATASETS <- c("neuron_neuropil", "neuron_soma", "microglia")

aff <- nv_read_csv(repo_path("results", "tables", "11_spatial_systems", "atlas",
                             "WGCNA_module_spatial_cell_affinity.csv"))
hubs <- nv_read_csv(repo_path("results", "reviewer_audit",
                              "wgcna_label_adjudication",
                              "WGCNA_module_top25_hubs.csv"))
nam <- nv_read_csv(file.path(OUT, "wgcna_module_naming_audit.csv"))

read_ds <- function(ds, file) {
  p <- repo_path("results", "tables", "06_modules_WGCNA", "01_WGCNA", ds,
                 "modules", file)
  if (!file.exists(p)) return(NULL)
  x <- nv_read_csv(p); x$dataset <- ds; x
}
go <- do.call(rbind, lapply(DATASETS, read_ds,
                            file = "WGCNA_module_GO_enrichment_long.csv"))
summ <- do.call(rbind, lapply(DATASETS, read_ds,
                              file = "WGCNA_module_summary.csv"))
if (is.null(go) || is.null(summ))
  stop("canonical WGCNA module tables not found", call. = FALSE)

mid <- function(x) sub("^WGCNA_", "", as.character(x))
aff$module_id <- mid(aff$ModuleID)
hubs$module_id <- mid(hubs$ModuleID)
go$module_id <- mid(go$ModuleID)
summ$module_id <- mid(summ$ModuleID)

# the widest protein-set scope, so a module is never called unsupported merely
# because a hub-weighted subset was too small to reach significance
SCOPE <- "all"
if ("ModuleProteinSetType" %in% names(go)) {
  sc <- unique(go$ModuleProteinSetType)
  SCOPE <- if ("all" %in% sc) "all" else sc[1]
  go <- go[go$ModuleProteinSetType == SCOPE, , drop = FALSE]
}

top_terms <- function(ds, m, n = 3L) {
  z <- go[go$dataset == ds & go$module_id == m, , drop = FALSE]
  z <- z[is.finite(z$p.adjust), , drop = FALSE]
  if (!nrow(z)) return(list(txt = "", n_sig = 0L, best = NA_real_,
                            best_desc = "", best_ont = ""))
  z <- z[order(z$p.adjust, -z$Count, z$ID), , drop = FALSE]
  sig <- z[z$p.adjust < FDR, , drop = FALSE]
  h <- utils::head(z, n)
  list(txt = paste(sprintf("%s %s (p.adj %.3g, %d/%d)", h$Ontology, h$Description,
                           h$p.adjust, h$Count, h$MappedModuleSize),
                   collapse = "; "),
       n_sig = nrow(sig), best = z$p.adjust[1], best_desc = z$Description[1],
       best_ont = z$Ontology[1])
}

rows <- do.call(rbind, lapply(seq_len(nrow(aff)), function(i) {
  ds <- aff$dataset[i]; m <- aff$module_id[i]
  tt <- top_terms(ds, m)
  hz <- hubs[hubs$dataset == ds & hubs$module_id == m, , drop = FALSE]
  hz <- hz[order(hz$rank), , drop = FALSE]
  nm <- nam[nam$dataset == ds & nam$module_id == m, , drop = FALSE]
  sm <- summ[summ$dataset == ds & summ$module_id == m, , drop = FALSE]

  has_enrich <- tt$n_sig > 0L
  has_hub <- nrow(hz) > 0L
  has_spatial <- is.finite(aff$spatial_tau[i])
  has_external <- isTRUE(aff$external_any_significant[i] %in%
                           c(TRUE, "TRUE", "yes"))
  scope_ok <- identical(as.character(aff$external_celltype_scope_concordance[i]),
                        "all_scopes_agree")

  data.frame(
    dataset = ds, module_id = m,
    module_size = aff$module_size[i],
    median_abs_kME = if (nrow(sm)) sm$median_abs_kME[1] else NA_real_,
    n_enriched_terms_FDR05 = tt$n_sig,
    top_enrichment_terms = tt$txt,
    best_term = tt$best_desc, best_term_ontology = tt$best_ont,
    best_term_FDR = tt$best,
    hub_proteins = paste(utils::head(sprintf("%s (%.3f)", hz$GeneSymbol,
                                             hz$abs_kME), 10), collapse = ", "),
    spatial_peak_unit = aff$peak_unit[i],
    spatial_tau = aff$spatial_tau[i],
    spatial_peak_minus_second = aff$peak_minus_second[i],
    external_cell_type = aff$external_celltype_all[i],
    external_FDR = aff$external_FDR_all[i],
    external_scope_concordance = aff$external_celltype_scope_concordance[i],
    reference_panel = aff$strongest_reference_panel[i],
    reference_panel_FDR = aff$reference_panel_FDR[i],
    active_reviewed_label = if (nrow(nm)) nm$active_reviewed_label[1] else "",
    proposed_label = if (nrow(nm)) nm$proposed_label[1] else "",
    proposal_status = if (nrow(nm)) nm$proposal_status[1] else "",
    evidence_layers_present = paste(c("enrichment", "hubs", "spatial",
                                      "external")[c(has_enrich, has_hub,
                                                    has_spatial, has_external)],
                                    collapse = "+"),
    n_evidence_layers = sum(has_enrich, has_hub, has_spatial, has_external),
    has_enrichment = has_enrich, has_external = has_external,
    external_scopes_agree = scope_ok,
    stringsAsFactors = FALSE)
}))

# ---------------------------------------------------- S18 label confidence
#
# Enrichment coherence carries the classification because it is the only layer
# that speaks to FUNCTION. Hub composition corroborates it, spatial pattern
# describes where the module lives, and external affinity may only add context.
rows$annotation_confidence <- with(rows, ifelse(
  !has_enrichment & has_external, "CELL_CONTEXT_ONLY",
  ifelse(!has_enrichment, "UNRESOLVED",
  ifelse(best_term_FDR < 1e-10 & n_enriched_terms_FDR05 >= 20,
         "HIGH_CONFIDENCE_FUNCTIONAL",
  ifelse(best_term_FDR < 1e-4 & n_enriched_terms_FDR05 >= 5,
         "MODERATE_CONFIDENCE_FUNCTIONAL", "MIXED_FUNCTIONAL")))))

# does the label a reader is currently allowed to use agree with the module's
# own strongest enrichment? compared on content words, so wording differences
# do not count as disagreement
STOP <- c("and", "or", "the", "of", "to", "in", "a", "process", "processes",
          "regulation", "positive", "negative", "cellular", "complex",
          "activity", "protein", "via", "by", "involved", "with")
words <- function(s) {
  w <- unlist(strsplit(tolower(gsub("[^a-z ]", " ", tolower(s))), " +"))
  w <- w[nchar(w) > 3 & !w %in% STOP]
  # a 4-character stem, so synaptic/synapse and axonal/axon match
  unique(substr(w, 1, 4))
}
rows$label_matches_own_enrichment <- mapply(function(lab, best, nsig) {
  if (!nzchar(lab) || !nzchar(best) || nsig == 0L) return(NA)
  length(intersect(words(lab), words(best))) > 0L
}, rows$active_reviewed_label, rows$best_term, rows$n_enriched_terms_FDR05)

# A lexical test cannot adjudicate whether a label agrees with a GO term:
# "mitochondrial / energy metabolism" and "respiratory chain complex" share no
# word but say the same thing. This flag is therefore a REVIEW PROMPT, not a
# verdict - it marks a module whose label shares no stem with its strongest
# term AND whose enrichment is strong enough that the disagreement would
# matter. Every flagged module is listed for human adjudication; none is
# relabelled here.
rows$label_shares_no_stem_with_best_term <-
  !is.na(rows$label_matches_own_enrichment) & !rows$label_matches_own_enrichment
rows$label_review_required <-
  rows$label_shares_no_stem_with_best_term &
  rows$annotation_confidence == "HIGH_CONFIDENCE_FUNCTIONAL"

# ------------------------------------------- S17 permitted manuscript wording
rows$allowed_manuscript_form <- with(rows, ifelse(
  annotation_confidence %in% c("UNRESOLVED", "MIXED_FUNCTIONAL"),
  sprintf("%s %s - refer to by module ID; no functional label is supported",
          dataset, module_id),
  ifelse(annotation_confidence == "CELL_CONTEXT_ONLY",
    sprintf("%s %s - refer to by module ID; external %s affinity may be cited as cell-type CONTEXT only",
            dataset, module_id, external_cell_type),
  ifelse(label_review_required,
    sprintf("%s %s - refer to by module ID; the active label disagrees with its own strongest enrichment (%s)",
            dataset, module_id, best_term),
    sprintf("%s %s, enriched for %s proteins", dataset, module_id,
            best_term)))))
rows$prohibited_manuscript_form <- with(rows, ifelse(
  nzchar(proposed_label),
  sprintf("the %s module; %s (PROPOSED ONLY, not activated in the registry)",
          ifelse(nzchar(active_reviewed_label), active_reviewed_label, module_id),
          proposed_label),
  sprintf("the %s module", ifelse(nzchar(active_reviewed_label),
                                  active_reviewed_label, module_id))))
rows <- rows[order(rows$dataset, rows$module_id), , drop = FALSE]
write_csv_safe(rows, file.path(OUT, "wgcna_annotation_evidence_matrix.csv"))

cat("\n===== PART-28 WGCNA ANNOTATION EVIDENCE =====\n")
cat("modules:", nrow(rows), "| GO scope used:", SCOPE, "\n\n")
print(table(rows$annotation_confidence))
cat("\nevidence layers present:\n"); print(table(rows$n_evidence_layers))
cat("\nactive label disagrees with the module's own strongest enrichment:",
    sum(rows$label_review_required), "\n")
if (any(rows$label_review_required))
  print(rows[rows$label_review_required,
             c("dataset", "module_id", "active_reviewed_label", "best_term",
               "best_term_FDR", "annotation_confidence")], row.names = FALSE)
cat("\nneuropil m11:\n")
print(t(rows[rows$dataset == "neuron_neuropil" & rows$module_id == "m11",
             c("module_size", "n_enriched_terms_FDR05", "best_term",
               "best_term_FDR", "external_cell_type", "external_FDR",
               "active_reviewed_label", "proposed_label",
               "annotation_confidence", "label_review_required")]))
cat("\nwritten to:", OUT, "\n")

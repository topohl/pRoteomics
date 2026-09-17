#!/usr/bin/env Rscript
#
#
# THE QUESTION
#   Is each module's name actually supported by the proteins at the network
#   centre of that module, and by broader annotation evidence (BP / CC / MF,
#   hubs, core membership, marker context)?
#
# STRICT PHENOTYPE BLINDNESS
#   This script must never read a group contrast, a differential-abundance
#   statistic, a SUS/RES/CON effect, or any candidate tier. If naming were
#   informed by outcome, every downstream biological interpretation would be
#   circular with the phenotype analysis it is used to explain. Every emitted
#   table is checked by wcl_assert_phenotype_blind().
#
#   The companion phenotype-AWARE analysis lives in
#   analysis/integration/quantify_candidate_network_position.R and is kept
#   deliberately separate. Neither script reads the other's outputs.
#
# WHAT THIS SCRIPT DOES NOT DO
#   * No WGCNA, no ORA, no label is recomputed. The frozen ORA result is left
#     exactly as it is.
#   * Nothing canonical is overwritten: not ModuleLabel_Final, not the Stage-07
#     final labels, not Figure 3 labels, not config/wgcna_labels/.
#   * Proposed labels are written to a review namespace for MANUAL adjudication
#     before any promotion.
#
# KEY DIAGNOSTIC
#   A GO term can be genuinely overrepresented while being driven mainly by
#   PERIPHERAL module members. The enrichment stays correct, but such a term is
#   a poor module NAME. contributor_centrality_auc measures exactly that.
#
# USAGE
#   Rscript analysis/wgcna/audit_module_label_coherence.R
#   Rscript analysis/wgcna/audit_module_label_coherence.R --dataset neuron_neuropil
#   Rscript analysis/wgcna/audit_module_label_coherence.R --dry-run
# Script: analysis/wgcna/audit_module_label_coherence.R
# Stage: networks
# Scope: per_dataset
# Consumes: required results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/WGCNA_modules_long.csv; results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/WGCNA_module_GO_enrichment_long.csv; results/tables/06_modules_WGCNA/identity_contract/<dataset>/WGCNA_module_supermodule_membership_contract.csv; optional results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/supermodules/wgcna_supermodule_biological_coherence.csv; results/tables/06_modules_WGCNA/interpretable_summary/<dataset>/WGCNA_final_label_lookup.csv; results/tables/06_modules_WGCNA/module_annotation/<dataset>/WGCNA_module_biological_annotation.csv
# Produces: results/reviewer_audit/wgcna_label_review/<dataset>/WGCNA_module_label_coherence_audit.csv; results/reviewer_audit/wgcna_label_review/<dataset>/WGCNA_supermodule_label_coherence_audit.csv; results/reviewer_audit/wgcna_label_review/<dataset>/WGCNA_label_supporting_proteins.csv; +2 more
# Dataset behavior: runs for neuron_neuropil,neuron_soma,microglia according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: PHENOTYPE-BLIND audit of WGCNA module and supermodule biological NAMES.

source("R/paths.R")
source("R/data_contracts/dataset_config.R")
source("R/statistics/integration_utils.R")
source("R/statistics/wgcna_downstream_utils.R")
source("R/statistics/wgcna_candidate_protein_utils.R")
source("R/statistics/wgcna_label_coherence_utils.R")
source("R/utilities/xlsx_package_utils.R")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

SCRIPT_ID <- "analysis/wgcna/audit_module_label_coherence.R"
SUBSTEP <- "label_coherence_audit"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)

cli <- integration_cli(default_dataset = "all")
datasets <- integration_datasets(cli$dataset)

membership_path <- function(dataset) {
  path_results("tables", "06_modules_WGCNA", "01_WGCNA", dataset, "modules",
               "WGCNA_modules_long.csv")
}
go_path <- function(dataset) {
  path_results("tables", "06_modules_WGCNA", "01_WGCNA", dataset, "modules",
               "WGCNA_module_GO_enrichment_long.csv")
}
# docs/OUTPUT_CONTRACTS.md: "The identity contract is the sole
# supermodule-membership authority."
membership_contract_path <- function(dataset) {
  path_results("tables", "06_modules_WGCNA", "identity_contract", dataset,
               "WGCNA_module_supermodule_membership_contract.csv")
}
structural_path <- function(dataset) {
  path_results("tables", "06_modules_WGCNA", "01_WGCNA", dataset, "supermodules",
               "wgcna_supermodule_biological_coherence.csv")
}
final_lookup_path <- function(dataset) {
  path_results("tables", "06_modules_WGCNA", "interpretable_summary", dataset,
               "WGCNA_final_label_lookup.csv")
}
annotation_path <- function(dataset) {
  path_results("tables", "06_modules_WGCNA", "module_annotation", dataset,
               "WGCNA_module_biological_annotation.csv")
}

if (isTRUE(cli$dry_run)) {
  inputs <- list()
  for (dataset in datasets) {
    inputs[[paste0("module_membership__", dataset)]] <- membership_path(dataset)
    inputs[[paste0("module_go_enrichment__", dataset)]] <- go_path(dataset)
    inputs[[paste0("supermodule_membership_contract__", dataset)]] <- membership_contract_path(dataset)
    inputs[[paste0("supermodule_structural_coherence__", dataset)]] <- structural_path(dataset)
    inputs[[paste0("stage07_final_label_lookup__", dataset)]] <- final_lookup_path(dataset)
    inputs[[paste0("stage06_module_annotation__", dataset)]] <- annotation_path(dataset)
  }
  cat("[DRY-RUN] Phenotype-blind WGCNA label coherence audit; no outputs written.\n")
  dry_run_inputs(SCRIPT_ID, inputs)
  cat("[DRY-RUN] Would write to results/reviewer_audit/wgcna_label_review/<dataset>/.\n")
  quit(save = "no", status = 0L)
}

read_required <- function(path, label) {
  if (!file.exists(path)) {
    stop("missing_required_input: ", label, ": ", path, call. = FALSE)
  }
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE, guess_max = Inf)
}
relpath <- function(path) relative_to(normalizePath(path, winslash = "/", mustWork = FALSE))
review_dir <- function(dataset) {
  d <- path_results("reviewer_audit", "wgcna_label_review", dataset)
  dir_create(d)
  d
}

# Top-n GO Descriptions per ontology, lowest p.adjust first, from the whole
# module ("all" protein set) - the same universe Stage 01 labels from.
top_go_terms <- function(go, module_id, ontology, n = 5L) {
  sub <- go[go$ModuleID == module_id & go$Ontology == ontology &
              go$ModuleProteinSetType == "all", , drop = FALSE]
  if (!nrow(sub)) return(NA_character_)
  sub <- sub[order(suppressWarnings(as.numeric(sub$p.adjust)), method = "radix"), ,
             drop = FALSE]
  paste(utils::head(as.character(sub$Description), n), collapse = "; ")
}

# ------------------------------------------------------------------ build

build_dataset <- function(dataset) {
  members <- read_required(membership_path(dataset), "frozen WGCNA membership")
  go <- read_required(go_path(dataset), "module GO enrichment")

  members <- wcp_rank_module_members(members)
  members$rank_fraction <- (members$abs_kME_rank_in_module - 0.5) /
    members$n_module_members

  label_terms <- wcl_label_defining_terms(go, ontology = "BP",
                                          protein_set_type = "all")
  evidence <- wcl_module_label_evidence(members, label_terms)

  # ---- current labels from every naming stage, for side-by-side review
  stage01 <- members %>%
    distinct(.data$ModuleID, .keep_all = TRUE) %>%
    select(any_of(c("ModuleID", "ModuleLabel_Final", "ModuleLabel_Source",
                    "primary_label", "final_label")))
  names(stage01)[names(stage01) == "ModuleLabel_Final"] <- "current_stage01_label"
  names(stage01)[names(stage01) == "ModuleLabel_Source"] <- "current_stage01_label_source"
  evidence <- left_join(evidence, stage01, by = "ModuleID")

  ann_path <- annotation_path(dataset)
  if (file.exists(ann_path)) {
    ann <- read_required(ann_path, "Stage-06 module annotation") %>%
      select(any_of(c("ModuleID", "cleaned_biological_label",
                      "cleaned_biological_label_confidence",
                      "microenvironment_label", "microenvironment_confidence",
                      "annotation_confidence", "label_warning",
                      "GO_label_relevance_flag"))) %>%
      distinct(.data$ModuleID, .keep_all = TRUE)
    names(ann)[names(ann) == "cleaned_biological_label"] <- "current_stage06_cleaned_label"
    evidence <- left_join(evidence, ann, by = "ModuleID")
  }

  lookup_path <- final_lookup_path(dataset)
  if (file.exists(lookup_path)) {
    lk <- read_required(lookup_path, "Stage-07 final label lookup")
    mod_col <- intersect(c("entity_id", "ModuleID", "module_id"), names(lk))[1]
    lab_col <- intersect(c("final_plot_label", "canonical_biological_label",
                           "best_data_driven_label", "final_label"), names(lk))[1]
    if (!is.na(mod_col) && !is.na(lab_col)) {
      lk2 <- lk[, c(mod_col, lab_col)]
      names(lk2) <- c("ModuleID", "current_stage07_final_label")
      lk2 <- lk2[lk2$ModuleID %in% evidence$ModuleID, , drop = FALSE]
      lk2 <- lk2[!duplicated(lk2$ModuleID), , drop = FALSE]
      evidence <- left_join(evidence, lk2, by = "ModuleID")
    }
  }

  # ---- GO context and the two INDEPENDENT semantic assignments
  evidence$top5_GO_BP <- vapply(evidence$ModuleID, top_go_terms, character(1),
                                go = go, ontology = "BP", n = 5L)
  evidence$top5_GO_CC <- vapply(evidence$ModuleID, top_go_terms, character(1),
                                go = go, ontology = "CC", n = 5L)
  evidence$top5_GO_MF <- vapply(evidence$ModuleID, top_go_terms, character(1),
                                go = go, ontology = "MF", n = 5L)

  # GO-ONLY theme: no hub tokens, so it is independent of the network centre.
  # HUB-ONLY theme: only the top-10 hub symbols.
  # Their agreement is a genuine, non-circular hub-concordance measure. The
  # existing Stage-06 scorer cannot do this: its hub_support term is a pure
  # non-emptiness test (R/statistics/wgcna_labeling_utils.R:66) that never compares hubs
  # with the candidate label, and it is TRUE for every row in every dataset.
  go_theme <- character(nrow(evidence))
  hub_theme <- character(nrow(evidence))
  go_conf <- character(nrow(evidence))
  for (i in seq_len(nrow(evidence))) {
    g <- wgcna_assign_semantic_program(
      raw_label = evidence$label_go_description[i],
      bp_terms = evidence$top5_GO_BP[i], cc_terms = evidence$top5_GO_CC[i],
      mf_terms = evidence$top5_GO_MF[i], hubs = character()
    )
    h <- wgcna_assign_semantic_program(
      raw_label = NA_character_, bp_terms = NULL, cc_terms = NULL, mf_terms = NULL,
      hubs = trimws(strsplit(evidence$top10_hub_symbols[i], ";", fixed = TRUE)[[1]])
    )
    go_theme[i] <- g$module_program_primary %||% NA_character_
    go_conf[i] <- g$label_confidence %||% NA_character_
    hub_theme[i] <- h$module_program_primary %||% NA_character_
  }
  evidence$proposed_module_theme <- go_theme
  evidence$proposed_module_theme_confidence <- go_conf
  evidence$hub_derived_theme <- hub_theme

  # The semantic-theme agreement is only meaningful where BOTH calls resolved to
  # a specific theme. The curated vocabulary is built around GO phrases and
  # carries gene-level tokens for only a few programs (myelin genes, hnRNPs), so
  # a hub-only call usually returns "mixed / low-specificity". That is a
  # vocabulary limitation, not evidence of discordance, and it is recorded
  # rather than scored.
  vague <- function(x) is.na(x) | x %in% c("mixed / low-specificity", "")
  evidence$hub_theme_evaluable <- !vague(go_theme) & !vague(hub_theme)
  evidence$hub_theme_concordant <- ifelse(
    evidence$hub_theme_evaluable, go_theme == hub_theme, NA
  )

  # PRIMARY hub-concordance measure, used for classification: the identity of
  # the proteins that actually generated the label-defining GO term, intersected
  # with the module's top-10 hubs. This uses exact GO contributor identity, not
  # keyword matching, so it is available for every module.
  evidence <- wcl_classify_module_labels(evidence, neural_dataset = TRUE)

  # Proposed label: a readable process name from the GO-only theme when the
  # current wording is generic or peripherally driven; otherwise the current
  # label is retained. Deliberately conservative - this is a review aid.
  evidence$proposed_label <- ifelse(
    evidence$review_action %in% c("KEEP"),
    evidence$current_stage01_label,
    ifelse(!is.na(evidence$proposed_module_theme),
           evidence$proposed_module_theme,
           evidence$current_stage01_label)
  )
  evidence$proposed_label_basis <- ifelse(
    evidence$review_action %in% c("KEEP"), "retain_current_label",
    "go_only_semantic_theme_pending_manual_review"
  )
  evidence$naming_evidence_is_phenotype_blind <- TRUE
  evidence$contract_version <- wcl_contract_version()
  evidence$Source <- SCRIPT_ID

  # ---- supporting proteins, long
  supporting <- members %>%
    transmute(
      dataset = .data$dataset, ModuleID = .data$ModuleID,
      ProteinGroupID = .data$ProteinGroupID, GeneSymbol = .data$GeneSymbol,
      EntrezID = .data$EntrezID, abs_kME = .data$abs_kME,
      abs_kME_rank_in_module = .data$abs_kME_rank_in_module,
      n_module_members = .data$n_module_members,
      rank_fraction = .data$rank_fraction,
      is_top5_module_representative = .data$is_top5_module_representative,
      is_top10_module_hub = .data$is_top10_module_hub,
      is_top_hub_25 = .data$is_top_hub_25
    )
  contrib_flag <- rep(FALSE, nrow(supporting))
  for (i in seq_len(nrow(label_terms))) {
    mod <- label_terms$ModuleID[i]
    rows <- which(supporting$ModuleID == mod)
    mapped <- wcl_map_label_contributors(label_terms$label_go_gene_ids[i],
                                         supporting[rows, , drop = FALSE])
    if (length(mapped$rows)) contrib_flag[rows[mapped$rows]] <- TRUE
  }
  supporting$contributes_to_label_term <- contrib_flag
  supporting <- left_join(
    supporting,
    label_terms[, c("ModuleID", "label_go_id", "label_go_description")],
    by = "ModuleID"
  )
  supporting$contract_version <- wcl_contract_version()

  # ---- supermodules
  contract <- read_required(membership_contract_path(dataset),
                            "supermodule membership contract")
  member_map <- contract %>%
    transmute(dataset = .data$dataset, SupermoduleID = .data$supermodule_id,
              ModuleID = .data$module_id)
  structural <- if (file.exists(structural_path(dataset))) {
    read_required(structural_path(dataset), "supermodule structural coherence")
  } else {
    data.frame(dataset = character(), SupermoduleID = character(),
               stringsAsFactors = FALSE)
  }
  sm <- wcl_supermodule_evidence(member_map, evidence, structural)
  sm <- wcl_apply_singleton_inheritance(sm, evidence)
  sm$naming_evidence_is_phenotype_blind <- TRUE
  sm$contract_version <- wcl_contract_version()
  sm$Source <- SCRIPT_ID

  list(module = evidence, supermodule = sm, supporting = supporting,
       label_terms = label_terms)
}

# ------------------------------------------------------- proposed registry

# A PROPOSED reviewed registry, shaped like config/wgcna_labels/microglia.csv so
# it can be adjudicated and then promoted by hand. Written to the review
# namespace only; nothing under config/ is touched.
proposed_registry <- function(dataset, module_evidence, supermodule_evidence) {
  mod <- data.frame(
    dataset = dataset,
    level = "module",
    entity_id = module_evidence$ModuleID,
    reviewed_biological_label = module_evidence$proposed_label,
    reviewed_short_label = module_evidence$proposed_module_theme,
    confidence = module_evidence$proposed_module_theme_confidence,
    manual_review_required = module_evidence$review_action != "KEEP",
    rationale = module_evidence$coherence_rationale,
    proposed_action = module_evidence$review_action,
    stringsAsFactors = FALSE
  )
  sm <- data.frame(
    dataset = dataset,
    level = "supermodule",
    entity_id = supermodule_evidence$SupermoduleID,
    reviewed_biological_label = supermodule_evidence$proposed_supermodule_label,
    reviewed_short_label = supermodule_evidence$dominant_member_theme,
    confidence = ifelse(
      supermodule_evidence$biological_coherence_class %in%
        c("coherent_single_program"), "high",
      ifelse(supermodule_evidence$biological_coherence_class %in%
               c("related_program_family", "singleton"), "moderate", "low")),
    manual_review_required = supermodule_evidence$supermodule_review_action != "KEEP",
    rationale = paste0(
      "structural=", supermodule_evidence$structural_coherence_class,
      "; biological=", supermodule_evidence$biological_coherence_class,
      "; dominant theme in ", supermodule_evidence$n_modules_supporting_dominant_theme,
      "/", supermodule_evidence$n_member_modules, " member modules"),
    proposed_action = supermodule_evidence$supermodule_review_action,
    stringsAsFactors = FALSE
  )
  out <- rbind(mod, sm)
  out$review_status <- "PROPOSED_NOT_PROMOTED"
  out$contract_version <- wcl_contract_version()
  out
}

# ------------------------------------------------------------------ workbook

write_review_workbook <- function(path, dataset, module_evidence,
                                  supermodule_evidence, supporting, registry) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) return(NA_character_)
  wb <- openxlsx::createWorkbook(creator = SCRIPT_ID)
  title_style <- openxlsx::createStyle(fontName = "Arial", fontSize = 14, fontColour = "#1F2933", textDecoration = "bold")
  note_style <- openxlsx::createStyle(fontName = "Arial", fontSize = 9, fontColour = "#5B6770", textDecoration = "italic", wrapText = TRUE, valign = "top")
  header_style <- openxlsx::createStyle(fontName = "Arial", fontSize = 9, fontColour = "#1F2933", fgFill = "#E9EDF0", textDecoration = "bold", wrapText = TRUE, halign = "center", valign = "center", border = "Bottom", borderColour = "#5B6770", borderStyle = "thin")
  body <- openxlsx::createStyle(fontName = "Arial", fontSize = 9, valign = "top")
  wrap <- openxlsx::createStyle(fontName = "Arial", fontSize = 9, wrapText = TRUE, valign = "top")

  add_sheet <- function(sheet, title, note, data) {
    openxlsx::addWorksheet(wb, sheet, gridLines = FALSE, tabColour = "#5B6770")
    nc <- max(1L, ncol(data))
    openxlsx::mergeCells(wb, sheet, cols = seq_len(nc), rows = 1)
    openxlsx::writeData(wb, sheet, title, startRow = 1, startCol = 1)
    openxlsx::addStyle(wb, sheet, title_style, rows = 1, cols = seq_len(nc), gridExpand = TRUE, stack = TRUE)
    openxlsx::mergeCells(wb, sheet, cols = seq_len(nc), rows = 2)
    openxlsx::writeData(wb, sheet, if (nrow(data)) note else paste0(note, "  [NONE IN THIS SCOPE]"), startRow = 2, startCol = 1)
    openxlsx::addStyle(wb, sheet, note_style, rows = 2, cols = seq_len(nc), gridExpand = TRUE, stack = TRUE)
    openxlsx::setRowHeights(wb, sheet, rows = 2, heights = 42)
    if (nrow(data)) {
      openxlsx::writeDataTable(wb, sheet, data, startRow = 4, startCol = 1,
                               tableName = gsub("[^A-Za-z0-9]", "", paste0(sheet, dataset)),
                               tableStyle = "TableStyleLight9", withFilter = TRUE)
      rows <- 5:(nrow(data) + 4L)
      openxlsx::addStyle(wb, sheet, body, rows = rows, cols = seq_len(ncol(data)), gridExpand = TRUE, stack = TRUE)
      for (nm in intersect(c("coherence_rationale", "top10_hub_symbols", "top25_hub_symbols",
                             "label_contributor_symbols", "top5_GO_BP", "top5_GO_CC",
                             "top5_GO_MF", "balanced_hub_panel", "rationale",
                             "member_module_labels", "member_module_proposed_themes"), names(data))) {
        openxlsx::addStyle(wb, sheet, wrap, rows = rows, cols = match(nm, names(data)), gridExpand = TRUE, stack = TRUE)
      }
    } else {
      openxlsx::writeData(wb, sheet, data[0, , drop = FALSE], startRow = 4, startCol = 1, colNames = TRUE)
    }
    openxlsx::addStyle(wb, sheet, header_style, rows = 4, cols = seq_len(ncol(data)), gridExpand = TRUE, stack = TRUE)
    openxlsx::setRowHeights(wb, sheet, rows = 4, heights = 44)
    widths <- pmin(52, pmax(12, nchar(names(data)) + 2))
    openxlsx::setColWidths(wb, sheet, cols = seq_len(ncol(data)), widths = widths)
    openxlsx::freezePane(wb, sheet, firstActiveRow = 5, firstActiveCol = 2)
  }

  # README
  openxlsx::addWorksheet(wb, "README", gridLines = FALSE, tabColour = "#23384D")
  openxlsx::writeData(wb, "README", paste0("WGCNA label coherence review - ", dataset), startRow = 1, startCol = 1)
  openxlsx::addStyle(wb, "README", title_style, rows = 1, cols = 1:3, gridExpand = TRUE, stack = TRUE)
  notes <- data.frame(Note = c(
    "PHENOTYPE-BLIND. No SUS/RES/CON contrast, no differential-abundance statistic and no candidate tier was used to choose or score any name. Every table here is checked by wcl_assert_phenotype_blind().",
    "Nothing canonical is overwritten. ModuleLabel_Final, the Stage-07 final labels, Figure 3 labels and config/wgcna_labels/ are untouched. Everything here is PROPOSED and awaits manual adjudication.",
    "The frozen ORA result is unchanged. A GO term can be genuinely overrepresented yet be driven by PERIPHERAL module members; that makes it a poor NAME, not a wrong enrichment.",
    "contributor_centrality_auc: probability that a label-contributing protein is more central than a non-contributor. 0.5 = no difference; >0.5 = the label describes the network centre; <0.5 = the label is driven by peripheral proteins.",
    "hub_theme_concordant: the GO-only semantic theme and the hub-only semantic theme agree. These are two INDEPENDENT calls of wgcna_assign_semantic_program() - one given only GO terms, one given only the top-10 hub symbols - so their agreement is a genuine hub-concordance measure.",
    "Stage-01 selects the label term by lowest p.adjust within the 'all' protein set with NO significance threshold (enrichGO runs with pvalueCutoff=1). label_go_is_fdr_significant flags terms that do not reach FDR 0.05.",
    "label_go_tie_broken_arbitrarily flags modules where several terms tied at the minimum p.adjust, so the winning name was decided by table row order rather than by evidence.",
    "Supermodule biology is assessed MODULE-BALANCED: member modules are weighted equally so a large module cannot dominate the name, and the illustrative hub panel takes the top hubs from each member module.",
    "A singleton supermodule inherits its single member module's proposed identity exactly rather than being independently named.",
    "ModuleID and SupermoduleID are immutable. Names are metadata."
  ), stringsAsFactors = FALSE)
  openxlsx::writeDataTable(wb, "README", notes, startRow = 3, startCol = 1, tableName = paste0("Notes", gsub("[^A-Za-z0-9]", "", dataset)), tableStyle = "TableStyleLight9")
  classes <- data.frame(
    Class = c(wcl_module_coherence_classes(), "", wcl_supermodule_biological_classes()),
    Meaning = c(
      "FDR-significant label term whose contributors are central, and >=20% of top-10 hubs contribute.",
      "FDR-significant and contributors are central, but the wording is not the most readable description.",
      "Enrichment is real but is driven by peripheral members; poor as a NAME.",
      "The term describes a very broad process and does not discriminate this module.",
      "Several terms tied at the minimum p.adjust and the term is not FDR-significant.",
      "Biology may be real but needs explicit context (e.g. a non-neural signature in a brain dataset).",
      "The label term is not FDR-significant and no clear alternative emerged.",
      "", "All member modules share one theme.", "Most member modules share a theme family.",
      "Themes differ but the eigengenes are structurally coherent.",
      "Themes differ and the structural support is weak.",
      "One member module; inherits its identity.", "No clear resolution."
    ), stringsAsFactors = FALSE)
  openxlsx::writeDataTable(wb, "README", classes, startRow = 3 + nrow(notes) + 3, startCol = 1, tableName = paste0("Classes", gsub("[^A-Za-z0-9]", "", dataset)), tableStyle = "TableStyleLight9")
  openxlsx::setColWidths(wb, "README", cols = 1:3, widths = c(48, 105, 20))
  openxlsx::addStyle(wb, "README", wrap, rows = 3:(3 + nrow(notes) + nrow(classes) + 5), cols = 1:3, gridExpand = TRUE, stack = TRUE)

  review_cols <- intersect(c(
    "dataset", "ModuleID", "module_size", "current_stage01_label",
    "current_stage01_label_source", "current_stage06_cleaned_label",
    "current_stage07_final_label", "label_go_id", "label_go_description",
    "label_go_p_adjust", "label_go_is_fdr_significant", "label_go_count",
    "label_go_gene_ratio", "label_go_n_terms_tied_at_min",
    "label_go_tie_broken_arbitrarily", "n_label_contributors_mapped",
    "label_mapping_fraction", "n_contributors_in_top10",
    "fraction_top10_hubs_contributing", "median_contributor_rank_fraction",
    "contributor_centrality_auc", "top10_hub_symbols", "label_contributor_symbols",
    "top5_GO_BP", "top5_GO_CC", "top5_GO_MF", "proposed_module_theme",
    "hub_derived_theme", "hub_theme_evaluable", "hub_theme_concordant", "coherence_class",
    "review_action", "proposed_label", "coherence_rationale"
  ), names(module_evidence))
  add_sheet("Module_review", paste0("Module label review - ", dataset),
    "One row per module: what it is called, why, its top hubs, which proteins generated the label-defining GO term, and whether that term describes the network centre.",
    module_evidence[, review_cols, drop = FALSE])

  add_sheet("Supermodule_review", paste0("Supermodule label review - ", dataset),
    "Is this supermodule structurally coherent, biologically coherent, both, or neither? Member modules are weighted equally; the hub panel is balanced across members.",
    supermodule_evidence)

  hubs <- supporting[supporting$is_top10_module_hub %in% TRUE, , drop = FALSE]
  hubs <- hubs[order(hubs$ModuleID, hubs$abs_kME_rank_in_module, method = "radix"), , drop = FALSE]
  add_sheet("Module_top10_hubs", paste0("Top-10 hubs per module - ", dataset),
    "The ten most central proteins of each module by |kME|, with whether each contributes to the label-defining GO term.", hubs)

  contribs <- supporting[supporting$contributes_to_label_term %in% TRUE, , drop = FALSE]
  contribs <- contribs[order(contribs$ModuleID, contribs$abs_kME_rank_in_module, method = "radix"), , drop = FALSE]
  add_sheet("Label_GO_contributors", paste0("Label-defining GO term contributors - ", dataset),
    "Every protein that contributed to the GO term that named its module, mapped from the frozen ORA geneID field back to canonical ProteinGroupIDs.", contribs)

  mism <- module_evidence[module_evidence$review_action != "KEEP", , drop = FALSE]
  add_sheet("Potential_mismatches", paste0("Labels needing review - ", dataset),
    "Modules whose name is not clearly supported: not FDR-significant, peripherally driven, generic wording, or decided by an arbitrary tie-break.",
    mism[, review_cols, drop = FALSE])

  add_sheet("Proposed_labels", paste0("Proposed reviewed registry - ", dataset),
    "PROPOSED ONLY. Shaped like config/wgcna_labels/microglia.csv so it can be adjudicated and promoted by hand. Nothing under config/ has been modified.",
    registry)

  xlsx_save_valid_workbook(wb, path)
  path
}

# --------------------------------------------------------------------- main

built <- list()
for (dataset in datasets) {
  message("[", SCRIPT_ID, "] auditing ", dataset)
  built[[dataset]] <- build_dataset(dataset)
}

written_all <- list()
for (dataset in datasets) {
  b <- built[[dataset]]
  out_dir <- review_dir(dataset)

  # Phenotype-blindness gate: refuse to emit anything outcome-derived.
  wcl_assert_phenotype_blind(b$module, "Module label coherence audit",
                             allow = c("naming_evidence_is_phenotype_blind"))
  wcl_assert_phenotype_blind(b$supermodule, "Supermodule label coherence audit",
                             allow = c("naming_evidence_is_phenotype_blind"))
  wcl_assert_phenotype_blind(b$supporting, "Label supporting proteins")

  registry <- proposed_registry(dataset, b$module, b$supermodule)
  wcl_assert_phenotype_blind(registry, "Proposed reviewed registry")

  paths <- c(
    module = file.path(out_dir, "WGCNA_module_label_coherence_audit.csv"),
    supermodule = file.path(out_dir, "WGCNA_supermodule_label_coherence_audit.csv"),
    supporting = file.path(out_dir, "WGCNA_label_supporting_proteins.csv"),
    registry = file.path(out_dir, paste0("PROPOSED_wgcna_reviewed_labels_", dataset, ".csv"))
  )
  write_csv_safe(b$module, paths[["module"]])
  write_csv_safe(b$supermodule, paths[["supermodule"]])
  write_csv_safe(b$supporting, paths[["supporting"]])
  write_csv_safe(registry, paths[["registry"]])
  wb_path <- write_review_workbook(
    file.path(out_dir, "WGCNA_label_review.xlsx"), dataset,
    b$module, b$supermodule, b$supporting, registry
  )
  written_all[[dataset]] <- c(paths, workbook = wb_path)
}

# ------------------------------------------------------------ console report

cat("\n===== WGCNA label coherence audit (phenotype-blind) =====\n")
cat("Contract : ", wcl_contract_version(), "\n", sep = "")
for (dataset in datasets) {
  e <- built[[dataset]]$module
  s <- built[[dataset]]$supermodule
  cat("\n--- ", dataset, " ---\n", sep = "")
  cat("  modules: ", nrow(e), "   supermodules: ", nrow(s),
      " (singletons: ", sum(s$n_member_modules == 1L), ")\n", sep = "")
  cat("  label term NOT FDR-significant : ",
      sum(!(e$label_go_is_fdr_significant %in% TRUE)), "\n", sep = "")
  cat("  label term decided by a tie    : ",
      sum(e$label_go_tie_broken_arbitrarily %in% TRUE), "\n", sep = "")
  cat("  median fraction of top-10 hubs contributing to label term: ",
      sprintf("%.2f", stats::median(e$fraction_top10_hubs_contributing, na.rm = TRUE)),
      "\n", sep = "")
  cat("  hub/GO semantic theme agreement: ",
      sum(e$hub_theme_concordant %in% TRUE), " of ",
      sum(e$hub_theme_evaluable %in% TRUE),
      " evaluable (curated vocabulary is GO-phrase based, so most hub-only calls are unresolved)\n",
      sep = "")
  cat("  coherence classes: ",
      paste(names(table(e$coherence_class)), table(e$coherence_class),
            sep = "=", collapse = ", "), "\n", sep = "")
  cat("  review actions   : ",
      paste(names(table(e$review_action)), table(e$review_action),
            sep = "=", collapse = ", "), "\n", sep = "")
  for (i in seq_len(nrow(e))) {
    cat(sprintf("    %-10s %-40s GOpadj=%8.2g%s AUC=%s hub10=%s  %-38s %s\n",
        e$ModuleID[i], substr(e$current_stage01_label[i], 1, 40),
        e$label_go_p_adjust[i],
        ifelse(e$label_go_is_fdr_significant[i] %in% TRUE, " ", "*"),
        ifelse(is.finite(e$contributor_centrality_auc[i]),
               sprintf("%.2f", e$contributor_centrality_auc[i]), "  NA"),
        ifelse(is.finite(e$fraction_top10_hubs_contributing[i]),
               sprintf("%.2f", e$fraction_top10_hubs_contributing[i]), "  NA"),
        substr(e$coherence_class[i], 1, 38), e$review_action[i]))
  }
}
cat("\n* = label term does not reach FDR 0.05\n")
cat("\nOutputs written under:\n")
for (dataset in datasets) cat("  ", relpath(review_dir(dataset)), "\n", sep = "")
cat("\nAll proposed labels are PROPOSED_NOT_PROMOTED. Nothing under config/ or any\n")
cat("canonical label artifact has been modified.\n")

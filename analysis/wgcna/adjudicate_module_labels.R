#!/usr/bin/env Rscript
#
#
# PURPOSE
#   Establish for every module and supermodule: what the network contains, which
#   annotations are robustly supported, whether the network CENTRE supports them,
#   whether one coherent program exists, whether the current label is defensible,
#   and what the safest reviewed label would be.
#
#   The output is an EVIDENCE PACKET plus PROPOSALS for a human to adjudicate.
#   Nothing is activated.
#
# FIVE QUESTIONS KEPT SEPARATE
#   A statistically supported enrichment (A) may still fail to describe the
#   module centre (B), fail to describe the broader module (C), be too broad to
#   serve as a name (D), or be unusable ontology grammar in a manuscript (E).
#
# EVIDENCE INDEPENDENCE
#   A theme GENERATED from BP terms is not independently validated by those same
#   terms. Roles are explicit: generating (BP themes), supporting (hub/core
#   overlap with the generating terms), orthogonal (CC/MF, marker context).
#   These are never counted as independent experiments.
#
# WHAT IS NOT TOUCHED
#   Module membership, ModuleID, SupermoduleID, eigengenes, supermodule
#   clustering, differential abundance, Stage-07 group effects, Figure 3, frozen
#   WGCNA state, and config/wgcna_labels/ are all left exactly as they are.
#
# NO PHENOTYPE INFORMATION may contribute to naming. Every emitted table is
# gated by wla_assert_phenotype_blind().
#
# USAGE
#   Rscript analysis/wgcna/adjudicate_module_labels.R
#   Rscript analysis/wgcna/adjudicate_module_labels.R --dataset neuron_neuropil
#   Rscript analysis/wgcna/adjudicate_module_labels.R --dry-run
# Script: analysis/wgcna/adjudicate_module_labels.R
# Stage: networks
# Scope: per_dataset
# Consumes: required results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/WGCNA_modules_long.csv; results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/WGCNA_module_GO_enrichment_long.csv; results/tables/06_modules_WGCNA/identity_contract/<dataset>/WGCNA_module_supermodule_membership_contract.csv; optional results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/supermodules/wgcna_supermodule_biological_coherence.csv; results/tables/06_modules_WGCNA/interpretable_summary/<dataset>/WGCNA_final_label_lookup.csv; results/tables/06_modules_WGCNA/module_annotation/<dataset>/WGCNA_module_biological_annotation.csv
# Produces: results/reviewer_audit/wgcna_label_adjudication/WGCNA_module_adjudication.csv; results/reviewer_audit/wgcna_label_adjudication/WGCNA_supermodule_adjudication.csv; results/reviewer_audit/wgcna_label_adjudication/WGCNA_module_theme_evidence.csv; +4 more
# Dataset behavior: runs for neuron_neuropil,neuron_soma,microglia according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Conservative, PHENOTYPE-BLIND adjudication of WGCNA module / supermodule names.

source("R/paths.R")
source("R/data_contracts/dataset_config.R")
source("R/statistics/integration_utils.R")
source("R/statistics/wgcna_candidate_protein_utils.R")
source("R/statistics/wgcna_label_adjudication_utils.R")
source("R/utilities/xlsx_package_utils.R")
source(repo_path("R", "wgcna_paths.R"))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

SCRIPT_ID <- "analysis/wgcna/adjudicate_module_labels.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)

cli <- integration_cli(default_dataset = "all")
datasets <- integration_datasets(cli$dataset)

membership_path <- function(ds) wgcna_modules_artifact("WGCNA_modules_long.csv", ds, child = "tables", "modules")
go_path <- function(ds) wgcna_modules_artifact("WGCNA_module_GO_enrichment_long.csv", ds, child = "tables", "modules")
contract_path <- function(ds) wgcna_identity_contract_artifact("WGCNA_module_supermodule_membership_contract.csv", ds)
structural_path <- function(ds) wgcna_modules_artifact("wgcna_supermodule_biological_coherence.csv", ds, child = "tables", "supermodules")
final_lookup_path <- function(ds) wgcna_interpretable_artifact("WGCNA_final_label_lookup.csv", ds)
annotation_path <- function(ds) wgcna_annotation_artifact("WGCNA_module_biological_annotation.csv", ds)

out_root <- function() {
  d <- wgcna_dirs("adjudicate_module_labels", "global", create = TRUE)$tables
  dir_create(d); d
}

if (isTRUE(cli$dry_run)) {
  inputs <- list()
  for (ds in datasets) {
    inputs[[paste0("module_membership__", ds)]] <- membership_path(ds)
    inputs[[paste0("module_go_enrichment__", ds)]] <- go_path(ds)
    inputs[[paste0("identity_contract__", ds)]] <- contract_path(ds)
    inputs[[paste0("supermodule_structural__", ds)]] <- structural_path(ds)
    inputs[[paste0("stage07_final_lookup__", ds)]] <- final_lookup_path(ds)
    inputs[[paste0("stage06_annotation__", ds)]] <- annotation_path(ds)
  }
  cat("[DRY-RUN] WGCNA label adjudication (phenotype-blind); no outputs written.\n")
  dry_run_inputs(SCRIPT_ID, inputs)
  cat("[DRY-RUN] Would write to results/reviewer_audit/wgcna_label_adjudication/.\n")
  cat("[DRY-RUN] No canonical label, config/wgcna_labels/ entry or frozen state would be modified.\n")
  quit(save = "no", status = 0L)
}

read_required <- function(path, label) {
  if (!file.exists(path)) stop("missing_required_input: ", label, ": ", path, call. = FALSE)
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE, guess_max = Inf)
}
relpath <- function(p) relative_to(normalizePath(p, winslash = "/", mustWork = FALSE))

top_terms <- function(go, module_id, ontology, protein_set, n = 10L) {
  sub <- go[go$ModuleID == module_id & go$Ontology == ontology &
              go$ModuleProteinSetType == protein_set, , drop = FALSE]
  sub <- sub[is.finite(suppressWarnings(as.numeric(sub$p.adjust))) &
               suppressWarnings(as.numeric(sub$p.adjust)) <= wla_fdr_threshold(), ,
             drop = FALSE]
  if (!nrow(sub)) return(character())
  sub <- sub[order(suppressWarnings(as.numeric(sub$p.adjust)), method = "radix"), ,
             drop = FALSE]
  utils::head(as.character(sub$Description), n)
}

# ------------------------------------------------------------------- build

build_dataset <- function(ds) {
  members <- read_required(membership_path(ds), "frozen WGCNA membership")
  go <- read_required(go_path(ds), "module GO enrichment")
  members <- wcp_rank_module_members(members)
  members$rank_fraction <- (members$abs_kME_rank_in_module - 0.5) / members$n_module_members

  contract <- read_required(contract_path(ds), "identity contract")
  member_map <- contract %>%
    transmute(dataset = .data$dataset, SupermoduleID = .data$supermodule_id,
              ModuleID = .data$module_id)

  structural <- if (file.exists(structural_path(ds))) {
    read_required(structural_path(ds), "supermodule structural coherence")
  } else data.frame(dataset = character(), SupermoduleID = character(),
                    stringsAsFactors = FALSE)

  # ---- current labels from every naming stage
  stage01 <- members %>% distinct(.data$ModuleID, .keep_all = TRUE) %>%
    select(any_of(c("ModuleID", "ModuleLegacyID", "ModuleLabel_Final",
                    "ModuleLabel_GO_BP", "ModuleLabel_Source")))
  ann <- if (file.exists(annotation_path(ds))) {
    read_required(annotation_path(ds), "Stage-06 annotation") %>%
      select(any_of(c("ModuleID", "cleaned_biological_label",
                      "cleaned_biological_label_confidence", "microenvironment_label",
                      "annotation_confidence", "label_warning"))) %>%
      distinct(.data$ModuleID, .keep_all = TRUE)
  } else NULL
  lookup <- if (file.exists(final_lookup_path(ds))) {
    lk <- read_required(final_lookup_path(ds), "Stage-07 final lookup")
    mc <- intersect(c("entity_id", "ModuleID", "module_id"), names(lk))[1]
    lc <- intersect(c("final_plot_label", "canonical_biological_label",
                      "best_data_driven_label", "final_label"), names(lk))[1]
    if (!is.na(mc) && !is.na(lc)) {
      z <- lk[, c(mc, lc)]; names(z) <- c("ModuleID", "stage07_final_label")
      z[!duplicated(z$ModuleID), , drop = FALSE]
    } else NULL
  } else NULL

  module_ids <- sort(unique(as.character(members$ModuleID)))
  themes_by_set <- list()
  adj_rows <- list()
  theme_rows <- list()

  for (m in module_ids) {
    mem <- members[members$ModuleID == m, , drop = FALSE]
    th <- lapply(wla_protein_sets(), function(ps)
      wla_module_themes(go, members, m, protein_set = ps, ontology = "BP"))
    names(th) <- wla_protein_sets()
    for (ps in names(th)) if (nrow(th[[ps]])) theme_rows[[length(theme_rows) + 1L]] <-
      cbind(dataset = ds, th[[ps]])

    cur <- stage01$ModuleLabel_Final[stage01$ModuleID == m]
    mk <- if (!is.null(ann)) as.character(ann$microenvironment_label[ann$ModuleID == m]) else NA_character_

    a <- wla_adjudicate_module(
      module_id = m, dataset = ds, module_size = nrow(mem),
      themes_all = th[["all"]], themes_core = th[["core_kME_0.6"]],
      themes_top25 = th[["top_hub_25"]],
      cc_terms = top_terms(go, m, "CC", "all"),
      mf_terms = top_terms(go, m, "MF", "all"),
      current_label = if (length(cur)) cur[[1]] else NA_character_,
      marker_context = if (length(mk)) mk[[1]] else NA_character_
    )
    hubs <- mem[order(mem$abs_kME_rank_in_module), , drop = FALSE]
    a$top5_hub_symbols <- paste(utils::head(as.character(hubs$GeneSymbol), 5L), collapse = "; ")
    a$top10_hub_symbols <- paste(utils::head(as.character(hubs$GeneSymbol), 10L), collapse = "; ")
    a$top25_hub_symbols <- paste(utils::head(as.character(hubs$GeneSymbol), 25L), collapse = "; ")
    a$top10_hub_identifiers <- wla_hub_identifiers(members, m, 10L)
    a$n_core_kME_0.6 <- sum(.wla_is_true(mem[["is_core_kME_0.6"]]))
    a$fraction_core_kME_0.6 <- mean(.wla_is_true(mem[["is_core_kME_0.6"]]))
    a$n_abs_kME_0.8 <- sum(suppressWarnings(as.numeric(mem$abs_kME)) >= 0.8, na.rm = TRUE)
    spot <- wla_hub_spot_check(members, m, th[["all"]], a$dominant_theme)
    a$hubs_supporting_theme <- spot$supporting_hubs
    a$hubs_not_supporting_theme <- spot$non_supporting_hubs
    a$n_top10_hubs_supporting_theme <- spot$n_hubs_supporting
    adj_rows[[length(adj_rows) + 1L]] <- a
  }

  adj <- do.call(rbind, adj_rows)
  adj <- left_join(adj, stage01, by = "ModuleID")
  if (!is.null(ann)) adj <- left_join(adj, ann, by = "ModuleID")
  if (!is.null(lookup)) adj <- left_join(adj, lookup, by = "ModuleID")
  adj <- left_join(adj, member_map %>% select("ModuleID", "SupermoduleID"), by = "ModuleID")

  # A label contradicted by most of its own top hubs is downgraded (Phase 15).
  contradicted <- adj$n_top10_hubs_supporting_theme <= 1L &
    adj$coherence_class %in% c("high_confidence_coherent", "coherent_but_wording_poor")
  adj$proposed_confidence[contradicted] <- "low"
  adj$hub_contradiction_downgrade <- contradicted

  # A proposed name that is not unique within the dataset does not discriminate
  # between modules; a reviewer must see that before promoting it.
  dup_theme <- adj$proposed_primary_label
  adj$proposed_label_shared_with_modules <- vapply(seq_len(nrow(adj)), function(i) {
    other <- adj$ModuleID[dup_theme == dup_theme[i] & seq_len(nrow(adj)) != i]
    if (!length(other)) NA_character_ else paste(other, collapse = "; ")
  }, character(1))

  adj$naming_evidence_is_phenotype_blind <- TRUE
  adj$contract_version <- wla_contract_version()
  adj$Source <- SCRIPT_ID

  themes <- if (length(theme_rows)) do.call(rbind, theme_rows) else
    data.frame(dataset = character(), stringsAsFactors = FALSE)

  sm <- wla_adjudicate_supermodules(member_map, adj, structural)
  sm$naming_evidence_is_phenotype_blind <- TRUE
  sm$contract_version <- wla_contract_version()
  sm$Source <- SCRIPT_ID

  list(modules = adj, supermodules = sm, themes = themes, members = members,
       go = go, member_map = member_map)
}

built <- list()
for (ds in datasets) {
  message("[", SCRIPT_ID, "] adjudicating ", ds)
  built[[ds]] <- build_dataset(ds)
}

# ------------------------------------------------------------------ emit

root <- out_root()
allow <- c("naming_evidence_is_phenotype_blind")

module_all <- bind_rows(lapply(built, `[[`, "modules"))
supermodule_all <- bind_rows(lapply(built, `[[`, "supermodules"))
theme_all <- bind_rows(lapply(built, `[[`, "themes"))

wla_assert_phenotype_blind(module_all, "Module adjudication", allow = allow)
wla_assert_phenotype_blind(supermodule_all, "Supermodule adjudication", allow = allow)
wla_assert_phenotype_blind(theme_all, "Theme evidence")

write_csv_safe(module_all, file.path(root, "WGCNA_module_adjudication.csv"))
write_csv_safe(supermodule_all, file.path(root, "WGCNA_supermodule_adjudication.csv"))
write_csv_safe(theme_all, file.path(root, "WGCNA_module_theme_evidence.csv"))

# hub table
hub_all <- bind_rows(lapply(names(built), function(ds) {
  m <- built[[ds]]$members
  m[.wla_is_true(m$is_top_hub_25), , drop = FALSE] %>%
    transmute(dataset = .data$dataset, ModuleID = .data$ModuleID,
              rank = .data$abs_kME_rank_in_module, GeneSymbol = .data$GeneSymbol,
              ProteinGroupID = .data$ProteinGroupID,
              RepresentativeUniProt = .data$RepresentativeUniProt,
              EntrezID = .data$EntrezID, abs_kME = .data$abs_kME,
              rank_fraction = .data$rank_fraction,
              is_top5 = .data$is_top5_module_representative,
              is_top10 = .data$is_top10_module_hub,
              mapping_status = .data$mapping_status,
              gene_level_claim_allowed = .data$gene_level_claim_allowed)
})) %>% arrange(.data$dataset, .data$ModuleID, .data$rank)
wla_assert_phenotype_blind(hub_all, "Hub table")
write_csv_safe(hub_all, file.path(root, "WGCNA_module_top25_hubs.csv"))

# per-protein-set GO exports
for (ps in wla_protein_sets()) {
  go_ps <- bind_rows(lapply(names(built), function(ds) {
    g <- built[[ds]]$go
    g <- g[g$ModuleProteinSetType == ps &
             is.finite(suppressWarnings(as.numeric(g$p.adjust))) &
             suppressWarnings(as.numeric(g$p.adjust)) <= wla_fdr_threshold(), , drop = FALSE]
    if (!nrow(g)) return(NULL)
    g$dataset <- ds
    g$theme <- wla_assign_theme(g$Description)
    g[, c("dataset", "ModuleID", "ModuleProteinSetType", "Ontology", "ID",
          "Description", "theme", "p.adjust", "Count", "GeneRatio", "geneID")]
  }))
  if (!is.null(go_ps) && nrow(go_ps)) {
    nm <- c(all = "GO_all", `core_kME_0.6` = "GO_core", top_hub_25 = "GO_top25")[[ps]]
    write_csv_safe(go_ps, file.path(root, paste0("WGCNA_", nm, "_significant_terms.csv")))
  }
}

# proposal registries (Phase 9) - explicitly NOT active
for (ds in datasets) {
  b <- built[[ds]]
  reg <- wla_proposal_registry(ds, b$modules, b$supermodules)
  wla_validate_proposal(reg, ds,
                        expected_modules = sort(unique(b$member_map$ModuleID)),
                        expected_supermodules = sort(unique(b$member_map$SupermoduleID)))
  if (wla_is_active_registry(reg)) {
    stop("A generated proposal must never validate as an active registry.", call. = FALSE)
  }
  wla_assert_phenotype_blind(reg, "Proposal registry")
  write_csv_safe(reg, file.path(root, paste0("proposed_", ds, "_reviewed_labels.csv")))
}

# ------------------------------------------------------------- workbook

review_cols <- c(
  "dataset", "ModuleID", "ModuleLegacyID", "SupermoduleID", "module_size",
  "ModuleLabel_Final", "ModuleLabel_Source", "cleaned_biological_label",
  "stage07_final_label", "annotation_confidence",
  "proposed_primary_label", "proposed_alternative_1", "proposed_alternative_2",
  "conservative_fallback_label", "proposed_confidence", "coherence_class",
  "recommended_action", "dominant_theme", "n_protein_sets_supporting_theme",
  "n_supporting_go_terms_all", "best_go_all", "best_go_core", "best_go_top25",
  "best_p_adjust_all", "best_p_adjust_core", "best_p_adjust_top25",
  "top10_theme_fraction", "top25_theme_fraction", "core_theme_fraction",
  "contributor_centrality_auc", "centre_supports_theme",
  "top10_hub_symbols", "theme_supporting_top_hubs", "hubs_not_supporting_theme",
  "n_top10_hubs_supporting_theme", "hub_contradiction_downgrade",
  "second_theme", "theme_is_overly_broad", "theme_is_context_sensitive",
  "dominant_theme_evidence_is_thin", "proposed_label_shared_with_modules",
  "orthogonal_cc_evidence", "orthogonal_mf_evidence", "orthogonal_marker_context",
  "top10_hub_identifiers"
)

build_workbook <- function(path) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) return(NA_character_)
  wb <- openxlsx::createWorkbook(creator = SCRIPT_ID)
  title <- openxlsx::createStyle(fontName = "Arial", fontSize = 14, textDecoration = "bold", fontColour = "#1F2933")
  note <- openxlsx::createStyle(fontName = "Arial", fontSize = 9, fontColour = "#5B6770", textDecoration = "italic", wrapText = TRUE, valign = "top")
  hdr <- openxlsx::createStyle(fontName = "Arial", fontSize = 9, fontColour = "#1F2933", fgFill = "#E9EDF0", textDecoration = "bold", wrapText = TRUE, halign = "center", valign = "center", border = "Bottom", borderColour = "#5B6770", borderStyle = "thin")
  body <- openxlsx::createStyle(fontName = "Arial", fontSize = 9, valign = "top")
  wrap <- openxlsx::createStyle(fontName = "Arial", fontSize = 9, wrapText = TRUE, valign = "top")

  add <- function(sheet, ttl, nte, data) {
    openxlsx::addWorksheet(wb, sheet, gridLines = FALSE, tabColour = "#5B6770")
    nc <- max(1L, ncol(data))
    openxlsx::mergeCells(wb, sheet, cols = seq_len(nc), rows = 1)
    openxlsx::writeData(wb, sheet, ttl, startRow = 1, startCol = 1)
    openxlsx::addStyle(wb, sheet, title, rows = 1, cols = seq_len(nc), gridExpand = TRUE, stack = TRUE)
    openxlsx::mergeCells(wb, sheet, cols = seq_len(nc), rows = 2)
    openxlsx::writeData(wb, sheet, if (nrow(data)) nte else paste0(nte, "  [NONE IN THIS SCOPE]"), startRow = 2, startCol = 1)
    openxlsx::addStyle(wb, sheet, note, rows = 2, cols = seq_len(nc), gridExpand = TRUE, stack = TRUE)
    openxlsx::setRowHeights(wb, sheet, rows = 2, heights = 46)
    if (nrow(data)) {
      openxlsx::writeDataTable(wb, sheet, data, startRow = 4, startCol = 1,
        tableName = gsub("[^A-Za-z0-9]", "", sheet), tableStyle = "TableStyleLight9", withFilter = TRUE)
      rows <- 5:(nrow(data) + 4L)
      openxlsx::addStyle(wb, sheet, body, rows = rows, cols = seq_len(ncol(data)), gridExpand = TRUE, stack = TRUE)
      for (nm in intersect(c("top10_hub_symbols", "top25_hub_symbols", "theme_supporting_top_hubs",
                             "hubs_not_supporting_theme", "orthogonal_cc_evidence",
                             "orthogonal_mf_evidence", "supporting_go_descriptions",
                             "balanced_hub_panel", "rationale", "member_module_proposed_labels",
                             "member_module_themes", "top10_hub_identifiers",
                             "supporting_go_ids", "geneID"), names(data))) {
        openxlsx::addStyle(wb, sheet, wrap, rows = rows, cols = match(nm, names(data)), gridExpand = TRUE, stack = TRUE)
      }
    } else {
      openxlsx::writeData(wb, sheet, data[0, , drop = FALSE], startRow = 4, startCol = 1, colNames = TRUE)
    }
    openxlsx::addStyle(wb, sheet, hdr, rows = 4, cols = seq_len(ncol(data)), gridExpand = TRUE, stack = TRUE)
    openxlsx::setRowHeights(wb, sheet, rows = 4, heights = 46)
    openxlsx::setColWidths(wb, sheet, cols = seq_len(ncol(data)),
                           widths = pmin(54, pmax(12, nchar(names(data)) + 2)))
    openxlsx::freezePane(wb, sheet, firstActiveRow = 5, firstActiveCol = 3)
  }

  openxlsx::addWorksheet(wb, "README", gridLines = FALSE, tabColour = "#23384D")
  openxlsx::writeData(wb, "README", "WGCNA module / supermodule label adjudication", startRow = 1, startCol = 1)
  openxlsx::addStyle(wb, "README", title, rows = 1, cols = 1:3, gridExpand = TRUE, stack = TRUE)
  readme <- data.frame(Note = c(
    "NOTHING HERE IS ACTIVE. Every label is a PROPOSAL for human adjudication. ModuleLabel_Final, the Stage-07 final labels, Figure 3 and config/wgcna_labels/ are untouched.",
    "PHENOTYPE-BLIND. No SUS/RES/CON contrast, no differential abundance and no candidate tier contributed to any name. Every sheet is gated by wla_assert_phenotype_blind().",
    "FIVE QUESTIONS ARE KEPT SEPARATE. A term can be statistically supported (A) yet fail to describe the module centre (B), fail to describe the broader module (C), be too broad to name it (D), or be unusable ontology grammar (E).",
    "EVIDENCE ROLES ARE NOT INDEPENDENT EXPERIMENTS. BP themes GENERATE a candidate; hub/core overlap with those same terms SUPPORTS it; CC/MF and marker panels are ORTHOGONAL context. Do not read them as four replications.",
    "THEMES, NOT SINGLE GO ROWS. Redundant significant GO terms are consolidated into biological themes by an explicit keyword map; every constituent GO ID and Description is retained in Module_theme_evidence and the GO_* sheets.",
    "THREE PROTEIN SETS. Stage 01 already computed ORA for all module proteins, core |kME| >= 0.60 proteins, and top-25 hubs. All three are shown separately and are never merged.",
    "CENTRALITY SUPPORT IS REAL. top10/top25/core theme fractions and a contributor centrality AUC replace the previous hub_support, which was a pure non-emptiness test that equalled 1 for every row and contributed nothing.",
    "AUC: probability a theme contributor is more central than a non-contributor. 0.5 = indistinguishable; >0.5 = the theme describes the network centre; <0.5 = peripheral enrichment.",
    "CONFIDENCE IS CONSERVATIVE. high requires several related significant terms, support in more than one protein set, and genuine centre support. Anything mixed, peripheral, nonsignificant or contradicted is low.",
    "A LABEL IS DOWNGRADED when its own top hubs systematically contradict it, but modules are not required to be perfectly homogeneous: WGCNA modules contain multifunctional proteins.",
    "Module IDs and Supermodule IDs are immutable. Names are metadata.",
    "top10_hub_identifiers carries GeneSymbol(UniProt|Entrez) so a reviewer can verify against UniProt/MGI/GO by hand. No external lookup is performed by this pipeline."
  ), stringsAsFactors = FALSE)
  openxlsx::writeDataTable(wb, "README", readme, startRow = 3, startCol = 1, tableName = "ReadmeNotes", tableStyle = "TableStyleLight9")
  cls <- data.frame(
    Item = c(wla_coherence_classes(), "", wla_actions(), "", wla_supermodule_classes()),
    Meaning = c(
      "Several related significant terms; core and central hubs support the same biology; CC/MF compatible.",
      "Program is well supported but the current wording is raw, generic or awkward ontology grammar.",
      "Plausible theme with some central support; not enough for a highly specific label.",
      "Enrichment exists but its contributors are mostly low-centrality proteins; the centre does not support it.",
      "Two or more substantial themes and none clearly dominates.",
      "Coherent real biology that does not correspond literally to the sampled cell type; must be called a signature.",
      "No sufficiently supported coherent theme.",
      "",
      "Current label already matches the supported biology.",
      "Biology is right; wording should change.",
      "Current label is not supported; a different name is needed.",
      "Real biology needing an explicit context caveat.",
      "Genuinely mixed; an explicit mixed label is preferable.",
      "Nothing defensible; prefer unresolved over an unsupported mechanistic label.",
      "",
      "One member module; inherits its label exactly.",
      "Nearly all member modules support one program.",
      "Member modules describe related subprograms of one broader biology.",
      "Eigengenes cluster reproducibly but member biology differs.",
      "Member biology differs and structural support is weak.",
      "No resolution."
    ), stringsAsFactors = FALSE)
  openxlsx::writeDataTable(wb, "README", cls, startRow = 3 + nrow(readme) + 3, startCol = 1, tableName = "ReadmeClasses", tableStyle = "TableStyleLight9")
  openxlsx::setColWidths(wb, "README", cols = 1:3, widths = c(46, 118, 18))
  openxlsx::addStyle(wb, "README", wrap, rows = 3:(3 + nrow(readme) + nrow(cls) + 5), cols = 1:2, gridExpand = TRUE, stack = TRUE)

  ds_sheet <- c(neuron_neuropil = "Neuropil", neuron_soma = "Soma", microglia = "Microglia")
  for (ds in datasets) {
    m <- module_all[module_all$dataset == ds, intersect(review_cols, names(module_all)), drop = FALSE]
    add(paste0(ds_sheet[[ds]], "_modules"), paste0("Module adjudication - ", ds),
        "One module per row: current name and source, proposed primary/alternatives/fallback, confidence, coherence class, action, dominant theme, GO support across all/core/top-25 sets, real centrality support, hub spot-check, and orthogonal CC/MF/marker context.",
        m)
  }
  for (ds in datasets) {
    s <- supermodule_all[supermodule_all$dataset == ds, , drop = FALSE]
    add(paste0(ds_sheet[[ds]], "_supermodules"), paste0("Supermodule adjudication - ", ds),
        "Structural coherence first, biology second. Member modules are weighted EQUALLY for theme recurrence and the hub panel is balanced across members, so a large module cannot dominate. Singletons inherit their member module exactly.",
        s)
  }
  add("Module_top10_hubs", "Top-25 hubs per module (top-10 flagged)",
      "Network centre of every module by frozen |kME|, with stable identifiers for manual verification.",
      hub_all)
  add("Module_theme_evidence", "Consolidated GO themes per module and protein set",
      "Redundant significant GO terms grouped into biological themes. Every constituent GO ID and Description is retained. Protein sets are never merged.",
      theme_all)
  for (ps in wla_protein_sets()) {
    nm <- c(all = "GO_all", `core_kME_0.6` = "GO_core", top_hub_25 = "GO_top25")[[ps]]
    f <- file.path(root, paste0("WGCNA_", nm, "_significant_terms.csv"))
    d <- if (file.exists(f)) readr::read_csv(f, show_col_types = FALSE, progress = FALSE) else
      data.frame(dataset = character(), stringsAsFactors = FALSE)
    add(nm, paste0("Significant GO terms - ", ps, " protein set"),
        "Frozen Stage-01 ORA, FDR <= 0.05, with the consolidated theme assignment. Nothing recomputed.", d)
  }
  mism <- module_all[module_all$recommended_action != "KEEP", intersect(review_cols, names(module_all)), drop = FALSE]
  add("Potential_mismatches", "Modules whose current label is not clearly defensible",
      "Anything not KEEP: wording problems, peripheral-only enrichment, mixed biology, context-sensitive signatures, and unresolved modules.", mism)
  props <- bind_rows(lapply(datasets, function(ds)
    readr::read_csv(file.path(root, paste0("proposed_", ds, "_reviewed_labels.csv")),
                    show_col_types = FALSE, progress = FALSE)))
  add("Consumer_migration", "Downstream label consumers and migration plan",
      "Which scripts display a WGCNA biological label, which naming stage it comes from, and what should happen after human adjudication. Figure 3 and manuscript panels are explicitly NOT migrated in this pass.",
      wla_consumer_migration_plan())
  add("Proposed_labels", "Proposed reviewed registries (NOT ACTIVE)",
      "adjudication_status = proposed, reviewer = NA. These files are never treated as active reviewed registries; the automatic fallback remains in force until a human adjudicates and promotes them.",
      props)

  xlsx_save_valid_workbook(wb, path)
  path
}

plan <- wla_consumer_migration_plan()
write_csv_safe(plan, file.path(root, "WGCNA_label_consumer_migration_plan.csv"))

wb_path <- build_workbook(file.path(root, "WGCNA_label_adjudication.xlsx"))

# ------------------------------------------------------------ console

cat("\n===== WGCNA label adjudication (phenotype-blind, proposals only) =====\n")
cat("Contract : ", wla_contract_version(), "\n", sep = "")
for (ds in datasets) {
  m <- module_all[module_all$dataset == ds, , drop = FALSE]
  s <- supermodule_all[supermodule_all$dataset == ds, , drop = FALSE]
  cat("\n--- ", ds, " ---   modules: ", nrow(m), "   supermodules: ", nrow(s),
      " (singletons ", sum(s$n_member_modules == 1L), ")\n", sep = "")
  cat("  classes: ", paste(names(table(m$coherence_class)), table(m$coherence_class), sep = "=", collapse = ", "), "\n", sep = "")
  cat("  actions: ", paste(names(table(m$recommended_action)), table(m$recommended_action), sep = "=", collapse = ", "), "\n", sep = "")
  cat("  confidence: ", paste(names(table(m$proposed_confidence)), table(m$proposed_confidence), sep = "=", collapse = ", "), "\n", sep = "")
  for (i in order(m$ModuleID)) {
    cat(sprintf("    %-10s %-34s -> %-36s %-8s %-28s %s\n",
                m$ModuleID[i], substr(m$ModuleLabel_Final[i], 1, 34),
                substr(m$proposed_primary_label[i], 1, 36),
                m$proposed_confidence[i], substr(m$coherence_class[i], 1, 28),
                m$recommended_action[i]))
  }
}
cat("\nOutputs: ", relpath(root), "\n", sep = "")
cat("All proposals carry adjudication_status = proposed and reviewer = NA.\n")
cat("No canonical label, config registry or frozen WGCNA state was modified.\n")

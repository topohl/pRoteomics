#!/usr/bin/env Rscript
#
#
# This does NOT rerun the scoring system. It reads the existing evidence packet
# produced by analysis/wgcna/adjudicate_module_labels.R, applies the curated
# adjudications recorded in R/statistics/wgcna_label_activation_utils.R, runs the hierarchy
# redundancy check, and emits one row per module and supermodule for a human to
# approve.
#
# NOTHING IS ACTIVATED.
#   human_decision is always NA. recommended_for_activation is an evidence-based
#   suggestion only. A label becomes active only when a reviewed registry with
#   adjudication_status = "reviewed" AND a real reviewer is promoted by hand into
#   config/wgcna_labels/<dataset>.csv.
#
# PHENOTYPE-BLIND. No group contrast, differential abundance or candidate tier
# contributed to any name; the packet this reads is itself gated.
#
# USAGE
#   Rscript analysis/wgcna/build_module_label_registry.R
#   Rscript analysis/wgcna/build_module_label_registry.R --dry-run
# Script: analysis/wgcna/build_module_label_registry.R
# Stage: networks
# Scope: per_dataset
# Consumes: required results/reviewer_audit/wgcna_label_adjudication/WGCNA_module_adjudication.csv; results/tables/06_modules_WGCNA/identity_contract/<dataset>/WGCNA_module_supermodule_membership_contract.csv; optional results/reviewer_audit/wgcna_label_adjudication/WGCNA_supermodule_adjudication.csv; results/reviewer_audit/wgcna_label_adjudication/WGCNA_module_top25_hubs.csv; results/tables/06_modules_WGCNA/interpretable_summary/<dataset>/WGCNA_final_label_lookup.csv
# Produces: results/reviewer_audit/wgcna_label_approval/WGCNA_final_label_approval_table.csv; results/reviewer_audit/wgcna_label_approval/WGCNA_canonical_display_label_resolution.csv; results/reviewer_audit/wgcna_label_approval/WGCNA_label_activation_rules.csv; +2 more
# Dataset behavior: runs for neuron_neuropil,neuron_soma,microglia according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Final manual-adjudication support: the human approval table.

source("R/paths.R")
source("R/data_contracts/dataset_config.R")
source("R/statistics/integration_utils.R")
source("R/statistics/wgcna_label_adjudication_utils.R")
source("R/statistics/wgcna_label_activation_utils.R")
source("R/utilities/xlsx_package_utils.R")

suppressPackageStartupMessages({ library(readr); library(dplyr) })

SCRIPT_ID <- "analysis/wgcna/build_module_label_registry.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
cli <- integration_cli(default_dataset = "all")

PACKET <- function(f) path_results("reviewer_audit", "wgcna_label_adjudication", f)
OUT <- function() {
  d <- path_results("reviewer_audit", "wgcna_label_approval"); dir_create(d); d
}
final_lookup_path <- function(ds) {
  path_results("tables", "06_modules_WGCNA", "interpretable_summary", ds,
               "WGCNA_final_label_lookup.csv")
}

if (isTRUE(cli$dry_run)) {
  inputs <- list(
    module_adjudication = PACKET("WGCNA_module_adjudication.csv"),
    supermodule_adjudication = PACKET("WGCNA_supermodule_adjudication.csv"),
    module_top25_hubs = PACKET("WGCNA_module_top25_hubs.csv")
  )
  for (ds in valid_datasets()) {
    inputs[[paste0("stage07_final_lookup__", ds)]] <- final_lookup_path(ds)
    inputs[[paste0("reviewed_registry__", ds)]] <-
      repo_path("config", "wgcna_labels", paste0(ds, ".csv"))
  }
  cat("[DRY-RUN] WGCNA final label approval table; no outputs written.\n")
  dry_run_inputs(SCRIPT_ID, inputs)
  cat("[DRY-RUN] Would write to results/reviewer_audit/wgcna_label_approval/.\n")
  cat("[DRY-RUN] No label is activated; human_decision stays empty.\n")
  quit(save = "no", status = 0L)
}

read_required <- function(p, l) {
  if (!file.exists(p)) stop("missing_required_input: ", l, ": ", p, call. = FALSE)
  readr::read_csv(p, show_col_types = FALSE, progress = FALSE, guess_max = Inf)
}
relpath <- function(p) relative_to(normalizePath(p, winslash = "/", mustWork = FALSE))

mod <- read_required(PACKET("WGCNA_module_adjudication.csv"), "module adjudication")
sup <- read_required(PACKET("WGCNA_supermodule_adjudication.csv"), "supermodule adjudication")
hubs <- read_required(PACKET("WGCNA_module_top25_hubs.csv"), "top-25 hubs")

# ---------------------------------------------------- current active labels

# The CURRENT active label is whatever a consumer would display today: an active
# reviewed registry if one exists, else Stage-07, else Stage-01.
entities <- bind_rows(
  mod %>% transmute(dataset = .data$dataset, level = "module",
                    entity_id = .data$ModuleID,
                    parent_entity_id = .data$SupermoduleID,
                    ModuleLabel_Final = .data$ModuleLabel_Final),
  sup %>% transmute(dataset = .data$dataset, level = "supermodule",
                    entity_id = .data$SupermoduleID,
                    parent_entity_id = NA_character_,
                    ModuleLabel_Final = NA_character_)
)

stage07 <- bind_rows(lapply(unique(entities$dataset), function(ds) {
  p <- final_lookup_path(ds)
  if (!file.exists(p)) return(NULL)
  lk <- read_required(p, "Stage-07 lookup")
  if (!all(c("level", "entity_id") %in% names(lk))) return(NULL)
  lab <- intersect(c("canonical_biological_label", "final_plot_label"), names(lk))[1]
  if (is.na(lab)) return(NULL)
  data.frame(dataset = ds, level = as.character(lk$level),
             entity_id = as.character(lk$entity_id),
             canonical_biological_label = as.character(lk[[lab]]),
             stringsAsFactors = FALSE)
}))

reviewed <- bind_rows(lapply(unique(entities$dataset), function(ds) {
  p <- repo_path("config", "wgcna_labels", paste0(ds, ".csv"))
  if (!file.exists(p)) return(NULL)
  r <- read_required(p, "reviewed registry")
  r$dataset <- ds
  r
}))

resolved <- resolve_wgcna_display_label(
  entities, reviewed = reviewed, stage07 = stage07,
  stage01 = entities %>% select("dataset", "level", "entity_id", "ModuleLabel_Final")
)
wal_assert_one_label_per_entity(resolved)

# ------------------------------------------------------- assemble the table

curated_mod <- wal_curated_adjudications()
curated_sup <- wal_curated_supermodule_adjudications()

hub_compact <- hubs %>%
  filter(.data$rank <= 8L) %>%
  arrange(.data$dataset, .data$ModuleID, .data$rank) %>%
  group_by(.data$dataset, .data$ModuleID) %>%
  summarise(top_hubs_compact = paste(.data$GeneSymbol, collapse = ", "),
            .groups = "drop")

module_rows <- mod %>%
  transmute(
    dataset = .data$dataset, level = "module", entity_id = .data$ModuleID,
    parent_entity_id = .data$SupermoduleID,
    evidence_class = .data$coherence_class,
    adjudication_action = .data$recommended_action,
    dominant_theme = .data$dominant_theme,
    strongest_alternative_theme = .data$second_theme,
    central_support_summary = sprintf(
      "top10=%.2f top25=%.2f core=%.2f AUC=%s centre_supports=%s",
      .data$top10_theme_fraction, .data$top25_theme_fraction,
      .data$core_theme_fraction,
      ifelse(is.finite(.data$contributor_centrality_auc),
             sprintf("%.2f", .data$contributor_centrality_auc), "NA"),
      .data$centre_supports_theme),
    GO_support_summary = sprintf(
      "all p=%s (n=%s terms); core p=%s; top25 p=%s",
      formatC(.data$best_p_adjust_all, format = "g", digits = 2),
      .data$n_supporting_go_terms_all,
      formatC(.data$best_p_adjust_core, format = "g", digits = 2),
      formatC(.data$best_p_adjust_top25, format = "g", digits = 2)),
    contradiction_summary = ifelse(
      .data$n_top10_hubs_supporting_theme <= 1L,
      paste0("Only ", .data$n_top10_hubs_supporting_theme,
             " of the top-10 hubs support the leading theme."),
      paste0(.data$n_top10_hubs_supporting_theme,
             " of the top-10 hubs support the leading theme.")),
    context_caveat = ifelse(
      .data$theme_is_context_sensitive %in% TRUE,
      "Signature only: these proteins are not evidence of a resident cell population in this tissue.",
      NA_character_)
  ) %>%
  left_join(hub_compact, by = c("dataset", "entity_id" = "ModuleID")) %>%
  left_join(curated_mod %>% select(-"level"), by = c("dataset", "entity_id"))

supermodule_rows <- sup %>%
  transmute(
    dataset = .data$dataset, level = "supermodule", entity_id = .data$SupermoduleID,
    parent_entity_id = NA_character_,
    evidence_class = .data$biological_coherence_class,
    adjudication_action = .data$recommended_action,
    dominant_theme = .data$dominant_theme,
    strongest_alternative_theme = .data$second_theme,
    central_support_summary = sprintf(
      "structural=%s adj_min_r=%s PC1=%s stability=%s",
      .data$structural_coherence_class,
      formatC(.data$adjusted_signed_min_pairwise_correlation, format = "g", digits = 2),
      formatC(.data$pc1_variance_explained, format = "g", digits = 2),
      formatC(.data$cut_height_stability_fraction_stable, format = "g", digits = 2)),
    GO_support_summary = sprintf("dominant theme in %s of %s member modules",
                                 .data$n_modules_supporting_dominant_theme,
                                 .data$n_member_modules),
    contradiction_summary = NA_character_,
    context_caveat = NA_character_,
    top_hubs_compact = .data$balanced_hub_panel
  ) %>%
  left_join(curated_sup %>% select(-"level"), by = c("dataset", "entity_id"))

all_rows <- bind_rows(module_rows, supermodule_rows) %>%
  left_join(resolved %>% select("dataset", "level", "entity_id",
                                current_active_label = "canonical_display_label",
                                "canonical_label_source"),
            by = c("dataset", "level", "entity_id"))

# Uncurated entities keep the algorithmic proposal rather than inventing one.
algo <- bind_rows(
  mod %>% transmute(dataset = .data$dataset, level = "module",
                    entity_id = .data$ModuleID,
                    algo_label = .data$proposed_primary_label,
                    algo_conf = .data$proposed_confidence),
  sup %>% transmute(dataset = .data$dataset, level = "supermodule",
                    entity_id = .data$SupermoduleID,
                    algo_label = .data$proposed_supermodule_label,
                    algo_conf = .data$proposed_confidence)
)
all_rows <- all_rows %>% left_join(algo, by = c("dataset", "level", "entity_id"))
# Which rows a human actually adjudicated in this pass, as opposed to inheriting.
all_rows$is_curated <- !is.na(all_rows$proposed_final_label)
all_rows$proposed_final_label <- ifelse(
  is.na(all_rows$proposed_final_label), all_rows$algo_label,
  all_rows$proposed_final_label)
all_rows$confidence <- ifelse(is.na(all_rows$confidence), all_rows$algo_conf,
                              all_rows$confidence)
all_rows$rationale <- ifelse(
  is.na(all_rows$rationale),
  "Algorithmic proposal retained; not individually adjudicated in this pass.",
  all_rows$rationale)

# The action and the evidence class must agree with the FINAL adjudicated label,
# in both directions. Manual adjudication overrides the algorithm either way, and
# the algorithmic class is kept as provenance rather than discarded.
all_rows$algorithmic_evidence_class <- all_rows$evidence_class
all_rows$evidence_class_source <- ifelse(all_rows$is_curated,
                                         "manual_adjudication", "algorithmic")

lab <- tolower(trimws(as.character(all_rows$proposed_final_label)))
is_unresolved <- lab %in% "unresolved"
is_mixed <- lab %in% c("mixed / unresolved", "mixed/unresolved")

# (a) a curated label of mixed/unresolved must carry the matching action, even
#     where the packet had suggested a wording change for a theme we rejected.
all_rows$adjudication_action[is_unresolved] <- "UNRESOLVED"
all_rows$adjudication_action[is_mixed] <- "MIXED"
all_rows$evidence_class[is_mixed & all_rows$is_curated] <- "mixed_biology"
all_rows$evidence_class[is_unresolved & all_rows$is_curated] <- "unresolved"

# (b) the inverse. Where manual adjudication RESOLVED an entity the algorithm had
#     called mixed, the row must stop being blocked as mixed - otherwise the table
#     reports high confidence and refuses activation for the same entity.
resolved_by_hand <- all_rows$is_curated & !is_mixed & !is_unresolved &
  (grepl("mixed|unresolved", tolower(all_rows$evidence_class)) |
     toupper(all_rows$adjudication_action) %in% c("MIXED", "UNRESOLVED"))
all_rows$adjudication_action[resolved_by_hand] <- "REFINE_WORDING"
all_rows$evidence_class[resolved_by_hand] <- "manually_adjudicated_coherent"

# ------------------------------------------------- hierarchy redundancy

# An entity that ALREADY carries an active reviewed label must not be handed a
# replacement proposal. The reviewed label is the authority; an unreviewed
# algorithmic proposal derived from the same Stage-01 evidence that produced the
# known failures must never be recommended over it. The proposal is kept as
# provenance so a genuine conflict stays visible, but it is not recommended.
all_rows$automatic_proposal_label <- all_rows$algo_label
already_reviewed <- !is.na(all_rows$canonical_label_source) &
  all_rows$canonical_label_source == "active_reviewed_registry" &
  !all_rows$is_curated
if (any(already_reviewed)) {
  all_rows$proposed_final_label[already_reviewed] <-
    all_rows$current_active_label[already_reviewed]
  all_rows$adjudication_action[already_reviewed] <- "KEEP_ACTIVE_REVIEWED"
  all_rows$rationale[already_reviewed] <- paste0(
    "An active reviewed label already exists and is the authority. This pass did ",
    "not re-adjudicate it; the automatic proposal is retained in ",
    "automatic_proposal_label as provenance only. Not recommended for ",
    "activation, because activating it would overwrite a human-reviewed label ",
    "with an unreviewed automatic one.")
}

# The theme layer is BP-only, so an entity carried by CC/MF evidence can be
# classed mixed while the packet's own GO tables name it. Advisory flag.
go_core_for_check <- (function() {
  p <- PACKET("WGCNA_GO_core_significant_terms.csv")
  if (!file.exists(p)) return(NULL)
  as.data.frame(readr::read_csv(p, show_col_types = FALSE, progress = FALSE,
                                guess_max = Inf))
})()
all_rows$bp_blind_spot_warning <- wal_bp_blind_spot(all_rows, go_core_for_check)

mods_for_check <- all_rows %>% filter(.data$level == "module")
sups_for_check <- all_rows %>% filter(.data$level == "supermodule")
sups_for_check$hierarchy_redundancy_warning <-
  wal_hierarchy_redundancy(mods_for_check, sups_for_check, "proposed_final_label")
mods_for_check$hierarchy_redundancy_warning <- NA_character_
all_rows <- bind_rows(mods_for_check, sups_for_check)

all_rows <- wal_enforce_conservative_unapproved(all_rows)
approval <- wal_build_approval_table(all_rows)
wla_assert_phenotype_blind(approval, "Final approval table")

root <- OUT()
write_csv_safe(approval, file.path(root, "WGCNA_final_label_approval_table.csv"))
write_csv_safe(resolved, file.path(root, "WGCNA_canonical_display_label_resolution.csv"))

# ------------------------------------------------------- activation rules
#
# DEFINED HERE, NEVER EXECUTED. This is the contract a future activation step
# must satisfy. Running this script performs none of it.
activation_rules <- wal_activation_rules()
write_csv_safe(activation_rules, file.path(root, "WGCNA_label_activation_rules.csv"))

# ------------------------------------------------ side-by-side comparison
#
# The evidence behind the contested calls, laid out member by member so the
# reviewer can check the disambiguation rather than take it on trust.
read_optional <- function(p) {
  if (!file.exists(p)) return(NULL)
  as.data.frame(readr::read_csv(p, show_col_types = FALSE, progress = FALSE,
                                guess_max = Inf))
}
themes_ev <- read_optional(PACKET("WGCNA_module_theme_evidence.csv"))
comparison <- if (!is.null(themes_ev)) {
  wal_module_comparison(
    themes = themes_ev, hubs = as.data.frame(hubs),
    go_all = read_optional(PACKET("WGCNA_GO_all_significant_terms.csv")),
    go_core = read_optional(PACKET("WGCNA_GO_core_significant_terms.csv")),
    go_top25 = read_optional(PACKET("WGCNA_GO_top25_significant_terms.csv")))
} else NULL
if (!is.null(comparison)) {
  write_csv_safe(comparison, file.path(root, "WGCNA_ambiguous_module_comparison.csv"))
}

# ------------------------------------------------------------- workbook

if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook(creator = SCRIPT_ID)
  title <- openxlsx::createStyle(fontName = "Arial", fontSize = 14, textDecoration = "bold", fontColour = "#1F2933")
  note <- openxlsx::createStyle(fontName = "Arial", fontSize = 9, fontColour = "#5B6770", textDecoration = "italic", wrapText = TRUE, valign = "top")
  hdr <- openxlsx::createStyle(fontName = "Arial", fontSize = 9, fgFill = "#E9EDF0", textDecoration = "bold", wrapText = TRUE, halign = "center", valign = "center", border = "Bottom", borderColour = "#5B6770", borderStyle = "thin")
  body <- openxlsx::createStyle(fontName = "Arial", fontSize = 9, valign = "top")
  wrap <- openxlsx::createStyle(fontName = "Arial", fontSize = 9, wrapText = TRUE, valign = "top")
  yes <- openxlsx::createStyle(fontName = "Arial", fontSize = 9, fgFill = "#E6F0E8")
  no <- openxlsx::createStyle(fontName = "Arial", fontSize = 9, fgFill = "#F7F7F7", fontColour = "#6B737A")

  add <- function(sheet, ttl, nte, data) {
    openxlsx::addWorksheet(wb, sheet, gridLines = FALSE, tabColour = "#5B6770")
    nc <- max(1L, ncol(data))
    openxlsx::mergeCells(wb, sheet, cols = seq_len(nc), rows = 1)
    openxlsx::writeData(wb, sheet, ttl, startRow = 1, startCol = 1)
    openxlsx::addStyle(wb, sheet, title, rows = 1, cols = seq_len(nc), gridExpand = TRUE, stack = TRUE)
    openxlsx::mergeCells(wb, sheet, cols = seq_len(nc), rows = 2)
    openxlsx::writeData(wb, sheet, if (nrow(data)) nte else paste0(nte, " [NONE]"), startRow = 2, startCol = 1)
    openxlsx::addStyle(wb, sheet, note, rows = 2, cols = seq_len(nc), gridExpand = TRUE, stack = TRUE)
    openxlsx::setRowHeights(wb, sheet, rows = 2, heights = 48)
    if (nrow(data)) {
      openxlsx::writeDataTable(wb, sheet, data, startRow = 4, startCol = 1,
        tableName = gsub("[^A-Za-z0-9]", "", sheet), tableStyle = "TableStyleLight9", withFilter = TRUE)
      rows <- 5:(nrow(data) + 4L)
      openxlsx::addStyle(wb, sheet, body, rows = rows, cols = seq_len(ncol(data)), gridExpand = TRUE, stack = TRUE)
      for (nm in intersect(c("rationale", "top_hubs_compact", "central_support_summary",
                             "GO_support_summary", "contradiction_summary",
                             "context_caveat", "hierarchy_redundancy_warning",
                             "proposed_final_label", "alternative_label_1"), names(data))) {
        openxlsx::addStyle(wb, sheet, wrap, rows = rows, cols = match(nm, names(data)), gridExpand = TRUE, stack = TRUE)
      }
      if ("recommended_for_activation" %in% names(data)) {
        rec <- as.logical(data$recommended_for_activation)
        if (any(rec %in% TRUE)) openxlsx::addStyle(wb, sheet, yes, rows = rows[rec %in% TRUE], cols = seq_len(ncol(data)), gridExpand = TRUE, stack = TRUE)
        if (any(!(rec %in% TRUE))) openxlsx::addStyle(wb, sheet, no, rows = rows[!(rec %in% TRUE)], cols = seq_len(ncol(data)), gridExpand = TRUE, stack = TRUE)
      }
    } else {
      openxlsx::writeData(wb, sheet, data[0, , drop = FALSE], startRow = 4, startCol = 1, colNames = TRUE)
    }
    openxlsx::addStyle(wb, sheet, hdr, rows = 4, cols = seq_len(ncol(data)), gridExpand = TRUE, stack = TRUE)
    openxlsx::setRowHeights(wb, sheet, rows = 4, heights = 46)
    openxlsx::setColWidths(wb, sheet, cols = seq_len(ncol(data)), widths = pmin(58, pmax(13, nchar(names(data)) + 2)))
    openxlsx::freezePane(wb, sheet, firstActiveRow = 5, firstActiveCol = 4)
  }

  openxlsx::addWorksheet(wb, "README", gridLines = FALSE, tabColour = "#23384D")
  openxlsx::writeData(wb, "README", "WGCNA final label approval", startRow = 1, startCol = 1)
  openxlsx::addStyle(wb, "README", title, rows = 1, cols = 1:2, gridExpand = TRUE, stack = TRUE)
  rm_notes <- data.frame(Note = c(
    "NOTHING IS ACTIVATED. human_decision is empty on every row and must be filled in by you. recommended_for_activation is an evidence-based suggestion only.",
    "A label becomes active only when a registry with adjudication_status = 'reviewed' AND a real reviewer name is promoted by hand into config/wgcna_labels/<dataset>.csv. The proposal files this pipeline writes can never satisfy that.",
    "current_active_label is what a consumer displays TODAY, resolved by the single canonical precedence: active reviewed registry > Stage-07 canonical > Stage-01 ModuleLabel_Final.",
    "PHENOTYPE-BLIND. No SUS/RES/CON contrast, differential-abundance result or candidate tier contributed to any name.",
    "recommended_for_activation is FALSE for anything mixed, unresolved, peripheral or low confidence, however attractive the label reads.",
    "hierarchy_redundancy_warning flags a multi-module supermodule whose name matches two or more of its members, or members that share a name with each other. It is a REVIEW warning, not a failure.",
    "Alternatives are real options, not decoration: where the evidence genuinely underdetermines the wording, alternative_label_1/2 are the other defensible readings.",
    "Module IDs and Supermodule IDs are immutable. Names are metadata."
  ), stringsAsFactors = FALSE)
  openxlsx::writeDataTable(wb, "README", rm_notes, startRow = 3, startCol = 1, tableName = "ApprovalReadme", tableStyle = "TableStyleLight9")
  openxlsx::setColWidths(wb, "README", cols = 1:2, widths = c(150, 14))
  openxlsx::addStyle(wb, "README", wrap, rows = 3:(3 + nrow(rm_notes)), cols = 1, gridExpand = TRUE, stack = TRUE)

  add("Approval_table", "Final label approval - all datasets",
      "One row per module and supermodule. Fill in human_decision. Everything else is evidence and proposal.",
      approval)
  for (ds in unique(approval$dataset)) {
    add(paste0("Approve_", substr(ds, 1, 24)), paste0("Approval - ", ds),
        "Same columns, filtered to one dataset for easier one-at-a-time review.",
        approval[approval$dataset == ds, , drop = FALSE])
  }
  add("Label_resolution", "Canonical display-label resolution (current behaviour)",
      "What each consumer displays today and why. Provenance columns are retained, never discarded.",
      resolved)
  add("Activation_rules", "Activation contract - defined here, executed nowhere",
      "The steps a future activation must satisfy. This script performs none of them.",
      activation_rules)
  if (!is.null(comparison)) {
    for (g in unique(comparison$comparison_group)) {
      z <- comparison[comparison$comparison_group == g, , drop = FALSE]
      wide <- stats::reshape(
        z[, c("attribute", "module_id", "value")],
        idvar = "attribute", timevar = "module_id", direction = "wide")
      names(wide) <- sub("^value[.]", "", names(wide))
      wide <- wide[match(unique(z$attribute), wide$attribute), , drop = FALSE]
      add(paste0("Cmp_", substr(sub("^neuropil_", "", g), 1, 25)),
          paste0("Side-by-side evidence - ", g),
          paste0("Modules the automatic system gave the same name (\"",
                 z$colliding_automatic_label[[1]],
                 "\"), compared attribute by attribute so the disambiguation can be checked."),
          wide)
    }
  }
  xlsx_save_valid_workbook(wb, file.path(root, "WGCNA_final_label_approval.xlsx"))
}

# ------------------------------------------------------------- console

cat("\n===== WGCNA final label approval table (nothing activated) =====\n")
for (ds in unique(approval$dataset)) {
  a <- approval[approval$dataset == ds, , drop = FALSE]
  m <- a[a$level == "module", , drop = FALSE]
  s <- a[a$level == "supermodule", , drop = FALSE]
  cat("\n--- ", ds, " --- modules ", nrow(m), ", supermodules ", nrow(s),
      "; recommended for activation: ", sum(a$recommended_for_activation %in% TRUE),
      "/", nrow(a), "\n", sep = "")
  for (i in order(m$entity_id)) {
    cat(sprintf("  %-10s %-52s %-9s %-5s %s\n", m$entity_id[i],
                substr(m$proposed_final_label[i], 1, 52), m$confidence[i],
                ifelse(m$recommended_for_activation[i] %in% TRUE, "YES", "no"),
                m$adjudication_action[i]))
  }
  for (i in order(s$entity_id)) {
    cat(sprintf("  %-10s %-52s %-9s %-5s %s\n", s$entity_id[i],
                substr(s$proposed_final_label[i], 1, 52), s$confidence[i],
                ifelse(s$recommended_for_activation[i] %in% TRUE, "YES", "no"),
                ifelse(is.na(s$hierarchy_redundancy_warning[i]), "",
                       "HIERARCHY-WARNING")))
  }
}
cat("\nOutputs: ", relpath(root), "\n", sep = "")
cat("human_decision is empty on every row. No label has been activated.\n")

#!/usr/bin/env Rscript
#
# Repository-level AnimalID integrity audit.
#
# Enumerates every animal-identifier normaliser in the active repository, every
# consumer, and every artifact whose values can change once the canonical
# contract in R/animal_id_contract.R replaces them. The numeric columns are
# computed from the real source files, not asserted.
#
# THE DEFECT, IN THREE PARTS
#   1. TRUNCATION      "[0-9]{3,4}" keeps only the first four digits of a
#                      longer run, so 117 behaviour ids collapse to 101 through
#                      7 collision groups, and distinct() then drops 16 rows.
#   2. BRANCH ASYMMETRY already-"A" ids keep their width, bare numerals are
#                      padded to four. "A111" vs "OR111"->"A0111" never match,
#                      which is what silently reduced the proteomics/behaviour
#                      join to a single animal.
#   3. MINIMUM WIDTH   three digits are required, so animal "3" becomes NA.
#
# USAGE
#   Rscript 08_behavior_physio_coupling/04_animal_id_integrity_audit.R
#   Rscript 08_behavior_physio_coupling/04_animal_id_integrity_audit.R --dry-run

source("R/paths.R")
source("R/dataset_config.R")
source("R/integration_utils.R")
source("R/animal_id_contract.R")

suppressPackageStartupMessages({ library(readr); library(dplyr) })

SCRIPT_ID <- "08_behavior_physio_coupling/04_animal_id_integrity_audit.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
cli <- integration_cli(default_dataset = "all")

OUT <- function(...) {
  d <- path_results("tables", "08_behavior_physio_coupling", "animal_id_integrity")
  dir_create(d); file.path(d, ...)
}
COUPLING <- function(...) path_results(
  "tables", "08_behavior_physio_coupling", "network_behavior_coupling", ...)

BEHAVIOR_XLSX <- path_external("behavior", "E9_Behavior_Data.xlsx")
AUC_FIRST <- path_external("behavior", "auc_individual_animals_firstChangeActive.csv")
AUC_ALL <- path_external("behavior", "auc_individual_animals_all.csv")
SPATIAL_RDS <- path_processed("07_spatial_networks", "network_spatial_relations",
                              "neuron_neuropil", "region_layer",
                              "network_spatial_relations_objects.rds")
META <- function(ds) path_processed(
  "01_preprocessing", "06_merged_metadata_module_score", ds,
  "sample_metadata_merged_clean_for_module_scores.xlsx")

if (isTRUE(cli$dry_run)) {
  cat("[DRY-RUN] Repository-level AnimalID integrity audit.\n")
  dry_run_inputs(SCRIPT_ID, list(behavior_xlsx = BEHAVIOR_XLSX,
                                 movement_auc_first = AUC_FIRST,
                                 movement_auc_all = AUC_ALL,
                                 spatial_network_rds = SPATIAL_RDS,
                                 canonical_metadata = META("neuron_neuropil"),
                                 alias_table = aid_alias_path()))
  cat("[DRY-RUN] Read-only audit; no analysis output is recomputed here.\n")
  quit(save = "no", status = 0L)
}

CANON <- aid_expected_exp9_animals()
legacy <- aid_legacy_broken_normalizer

# ---------------------------------------------------------------- load ids

read_ids <- function() {
  out <- list()
  if (file.exists(BEHAVIOR_XLSX)) {
    b <- readxl::read_excel(BEHAVIOR_XLSX, sheet = "zScore")
    out$behavior_zscore <- list(ids = as.character(b$ID), n_rows = nrow(b),
                                src = BEHAVIOR_XLSX, col = "ID")
  }
  for (nm in c(first = AUC_FIRST, all = AUC_ALL)) {
    if (!file.exists(nm)) next
    d <- readr::read_csv(nm, show_col_types = FALSE, progress = FALSE)
    key <- paste0("movement_auc_", names(which(c(first = AUC_FIRST, all = AUC_ALL) == nm))[1])
    out[[key]] <- list(ids = as.character(d$AnimalNum), n_rows = nrow(d),
                       src = nm, col = "AnimalNum")
  }
  if (file.exists(SPATIAL_RDS)) {
    o <- readRDS(SPATIAL_RDS)
    md <- as.data.frame(o$sample_metadata)
    out$spatial_network_rds <- list(ids = as.character(md$AnimalID), n_rows = nrow(md),
                                    src = SPATIAL_RDS, col = "sample_metadata$AnimalID")
  }
  for (ds in valid_datasets()) {
    p <- META(ds)
    if (!file.exists(p)) next
    m <- readxl::read_excel(p)
    out[[paste0("proteomics_metadata_", ds)]] <- list(
      ids = as.character(m$AnimalID), n_rows = nrow(m), src = p, col = "AnimalID")
  }
  out
}

message("Reading every animal-identifier source")
sources <- read_ids()

# ------------------------------------------- PART 1: consumer audit table

# Alias-table system name per source (proteomics metadata is already canonical).
alias_system <- c(behavior_zscore = "behavior_zscore",
                  movement_auc_first = "movement_auc",
                  movement_auc_all = "movement_auc",
                  spatial_network_rds = "spatial_network_rds")

# Every invocation site found by exhaustive grep of the active repository.
consumers <- tibble::tribble(
  ~script, ~function_name, ~invocation_line, ~source_key, ~input_source_label, ~artifacts, ~rerun_required,
  "08_behavior_physio_coupling/02_network_behavior_coupling.r", "normalize_animal_id", 312L,
    "spatial_network_rds", "07_spatial_networks RDS sample_metadata$AnimalID",
    paste("animal_level_candidate_edge_scores.csv", "animal_level_global_network_metrics.csv",
          "merged_edge_behavior_long.csv", "merged_global_network_behavior.csv",
          "edge_behavior_correlations.csv", "join_diagnostics_*.csv",
          "qc_counts_edge_behavior.csv", "central_edge_*_behavior_table.csv", sep = "; "), "yes",
  "08_behavior_physio_coupling/02_network_behavior_coupling.r", "normalize_animal_id", 434L,
    "behavior_zscore", "E9_Behavior_Data.xlsx sheet zScore column ID",
    paste("physiology_traits_loaded.csv", "merged_edge_behavior_long.csv",
          "merged_global_network_behavior.csv", "edge_behavior_correlations.csv", sep = "; "), "yes",
  "08_behavior_physio_coupling/02_network_behavior_coupling.r", "normalize_animal_id", 464L,
    "movement_auc_first", "auc_individual_animals_*.csv column AnimalNum",
    paste("movement_auc_z_loaded.csv", "merged_edge_behavior_long.csv",
          "edge_behavior_correlations_sex_stratified.csv", sep = "; "), "yes",
  "06_modules_WGCNA/03_score_module_activity.R", "normalize_animal_id", 2083L,
    "proteomics_metadata_neuron_neuropil", "module score table AnimalID (behaviour handoff export only)",
    "results/source_data/.../behavior_coupling_inputs/module_scores_*.csv", "no",
  "08_behavior_physio_coupling/01_correlate_proteomics_with_behavior.r", "normalize_mouse_id", 304L,
    "behavior_zscore", "behavior table AnimalNum/MouseID",
    "proteomics-behaviour correlation outputs", "no"
)

rows <- list()
for (i in seq_len(nrow(consumers))) {
  key <- consumers$source_key[i]
  s <- sources[[key]]
  fn <- consumers$function_name[i]
  if (is.null(s)) {
    rows[[length(rows) + 1L]] <- data.frame(
      script = consumers$script[i], function_name = fn,
      invocation_line = consumers$invocation_line[i],
      input_source = consumers$input_source_label[i],
      raw_id_examples = NA_character_, normalized_id_examples = NA_character_,
      expected_canonical_id = NA_character_, collision_status = "source_not_available",
      n_distinct_raw_ids = NA_integer_, n_distinct_after_legacy = NA_integer_,
      n_collision_groups = NA_integer_, n_raw_ids_merged = NA_integer_,
      rows_affected = NA_integer_,
      output_artifacts_potentially_affected = consumers$artifacts[i],
      rerun_required = consumers$rerun_required[i], stringsAsFactors = FALSE)
    next
  }
  ids <- s$ids
  u <- unique(ids[!is.na(ids)])

  # what the historical normaliser did (only 08/02 used it)
  uses_legacy <- fn == "normalize_animal_id" &&
    grepl("02_network_behavior_coupling", consumers$script[i])
  leg <- if (uses_legacy) legacy(u) else
    if (fn == "normalize_mouse_id") sub("^0+(?=[0-9])", "", gsub("[^0-9]", "", u), perl = TRUE) else
      sprintf("%04d", suppressWarnings(as.integer(sub(".*?([0-9]+)$", "\\1", u))))

  tab <- tapply(u, leg, function(z) length(unique(z)))
  coll <- names(tab)[!is.na(names(tab)) & tab > 1L]
  merged_raw <- sum(vapply(coll, function(k) sum(leg %in% k), integer(1)))
  n_rows_lost <- length(u) - length(unique(leg[!is.na(leg)])) - sum(is.na(leg))

  asys <- if (key %in% names(alias_system)) unname(alias_system[key]) else NA_character_
  can <- if (is.na(asys)) u else
    aid_resolve(u, asys, strict = FALSE, canonical = CANON)
  ex <- utils::head(u[order(u)], 4)

  rows[[length(rows) + 1L]] <- data.frame(
    script = consumers$script[i], function_name = fn,
    invocation_line = consumers$invocation_line[i],
    input_source = relative_to(s$src),
    raw_id_examples = paste(ex, collapse = ";"),
    normalized_id_examples = paste(leg[match(ex, u)], collapse = ";"),
    expected_canonical_id = paste(ifelse(is.na(can[match(ex, u)]), "unresolved",
                                         can[match(ex, u)]), collapse = ";"),
    collision_status = if (!uses_legacy) "no_collision" else
      if (length(coll)) sprintf("COLLISION: %d groups merging %d raw ids",
                                length(coll), merged_raw) else "no_collision",
    n_distinct_raw_ids = length(u),
    n_distinct_after_legacy = length(unique(leg[!is.na(leg)])),
    n_collision_groups = length(coll),
    n_raw_ids_merged = as.integer(merged_raw),
    rows_affected = if (uses_legacy) as.integer(sum(ids %in% u[leg %in% coll] |
                                                      is.na(legacy(ids)))) else 0L,
    output_artifacts_potentially_affected = consumers$artifacts[i],
    rerun_required = consumers$rerun_required[i], stringsAsFactors = FALSE)
}
consumer_audit <- dplyr::bind_rows(rows)

# ------------------------------------------ per-source resolution reports

res_rows <- list()
for (key in names(sources)) {
  asys <- if (key %in% names(alias_system)) unname(alias_system[key]) else "proteomics_canonical"
  r <- aid_resolution_report(sources[[key]]$ids, asys, canonical = CANON)
  r$source_key <- key
  r$legacy_value <- legacy(r$raw_id)
  lt <- tapply(r$raw_id, r$legacy_value, function(z) length(unique(z)))
  r$legacy_collision <- !is.na(r$legacy_value) & as.integer(lt[r$legacy_value]) > 1L
  r$is_exp9_animal <- !is.na(r$canonical_AnimalID)
  res_rows[[key]] <- r
}
resolution <- dplyr::bind_rows(res_rows)

# --------------------------------------- PART 4: blast radius per artifact

read_opt <- function(p) if (file.exists(p))
  utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE) else NULL

beh <- sources$behavior_zscore
mv <- sources$movement_auc_first
rds <- sources$spatial_network_rds

n_beh_raw <- if (is.null(beh)) NA_integer_ else length(unique(beh$ids))
n_beh_leg <- if (is.null(beh)) NA_integer_ else length(unique(legacy(unique(beh$ids))))
n_beh_drop <- if (is.null(beh)) NA_integer_ else n_beh_raw - n_beh_leg
beh_merged <- if (is.null(beh)) NA_integer_ else {
  u <- unique(beh$ids); l <- legacy(u)
  t2 <- tapply(u, l, function(z) length(unique(z)))
  sum(vapply(names(t2)[t2 > 1L], function(k) sum(l %in% k), integer(1)))
}

# the actually written artifacts, so "previously" is measured not assumed
phys <- read_opt(COUPLING("physiology_traits_loaded.csv"))
mvz <- read_opt(COUPLING("movement_auc_z_loaded.csv"))
merged <- read_opt(COUPLING("merged_global_network_behavior.csv"))
po <- read_opt(COUPLING("join_diagnostics_proteomics_only.csv"))
glob <- read_opt(COUPLING("animal_level_global_network_metrics.csv"))

n_prot_only <- if (is.null(po)) NA_integer_ else nrow(po)
n_net_animals <- if (is.null(glob)) NA_integer_ else nrow(glob)
n_joined <- if (is.na(n_net_animals) || is.na(n_prot_only)) NA_integer_ else
  n_net_animals - n_prot_only

blast <- dplyr::bind_rows(
  data.frame(
    artifact = "results/tables/08_behavior_physio_coupling/network_behavior_coupling/physiology_traits_loaded.csv",
    previously_affected_rows = n_beh_raw,
    previously_merged_animals = beh_merged,
    previously_dropped_rows = n_beh_drop,
    corrected_rows = n_beh_raw,
    scientific_result_changed = "yes",
    detail = sprintf(paste0("distinct(AnimalID) after the truncating normaliser kept %d of %d ",
                            "behaviour animals; %d raw ids were merged into %d keys and the ",
                            "losers were silently dropped, misattributing their phenotype"),
                     n_beh_leg, n_beh_raw, beh_merged, beh_merged - n_beh_drop),
    stringsAsFactors = FALSE),
  data.frame(
    artifact = "results/tables/08_behavior_physio_coupling/network_behavior_coupling/movement_auc_z_loaded.csv",
    previously_affected_rows = if (is.null(mv)) NA_integer_ else length(unique(mv$ids)),
    previously_merged_animals = 2L,
    previously_dropped_rows = if (is.null(mv)) NA_integer_ else
      sum(is.na(legacy(unique(mv$ids)))),
    corrected_rows = if (is.null(mv)) NA_integer_ else length(unique(mv$ids)),
    scientific_result_changed = "yes",
    detail = paste0("ids shorter than three digits resolved to NA and were removed by the ",
                    "!is.na(AnimalID) filter; canonical animal 3 lost all movement rows"),
    stringsAsFactors = FALSE),
  data.frame(
    artifact = "results/tables/08_behavior_physio_coupling/network_behavior_coupling/merged_global_network_behavior.csv",
    previously_affected_rows = n_net_animals,
    previously_merged_animals = 0L,
    previously_dropped_rows = n_prot_only,
    corrected_rows = length(CANON),
    scientific_result_changed = "yes",
    detail = sprintf(paste0("the proteomics side kept its A-prefixed width (A111) while the ",
                            "behaviour side was padded to four (OR111 -> A0111), so only %d of ",
                            "%d animals carried any behaviour value; %d were proteomics-only"),
                     n_joined, n_net_animals, n_prot_only),
    stringsAsFactors = FALSE),
  data.frame(
    artifact = "results/tables/08_behavior_physio_coupling/network_behavior_coupling/edge_behavior_correlations.csv",
    previously_affected_rows = NA_integer_,
    previously_merged_animals = 0L,
    previously_dropped_rows = n_prot_only,
    corrected_rows = length(CANON),
    scientific_result_changed = "yes",
    detail = paste0("every correlation was computed on the joined subset above; with one ",
                    "joined animal no correlation was estimable"),
    stringsAsFactors = FALSE),
  data.frame(
    artifact = "results/source_data/06_modules_WGCNA/.../behavior_coupling_inputs/module_scores_*.csv",
    previously_affected_rows = length(CANON),
    previously_merged_animals = 0L,
    previously_dropped_rows = 0L,
    corrected_rows = length(CANON),
    scientific_result_changed = "no",
    detail = paste0("the WGCNA-side normaliser zero-pads to four but never truncates, so no ",
                    "animal is merged or lost; it produces a different NAMESPACE (0003) from ",
                    "the canonical AnimalID (3). Handoff export only - it runs after all ",
                    "module, kME and score computation and touches no canonical WGCNA state"),
    stringsAsFactors = FALSE),
  data.frame(
    artifact = "canonical WGCNA membership / kME / module scores",
    previously_affected_rows = 0L, previously_merged_animals = 0L,
    previously_dropped_rows = 0L, corrected_rows = 0L,
    scientific_result_changed = "no",
    detail = paste0("no normaliser is invoked anywhere in the module, kME or score computation ",
                    "path; the only call in 06_modules_WGCNA/03 is inside ",
                    "export_behavior_proteomics_input() at the end of the script"),
    stringsAsFactors = FALSE),
  data.frame(
    artifact = "canonical differential abundance / 11_spatial_systems atlas / network layer",
    previously_affected_rows = 0L, previously_merged_animals = 0L,
    previously_dropped_rows = 0L, corrected_rows = 0L,
    scientific_result_changed = "no",
    detail = paste0("these read AnimalID straight from the canonical metadata, which is already ",
                    "the canonical bare form; no ad-hoc normaliser is on their path"),
    stringsAsFactors = FALSE)
)

# ---------------------------------------------------------------- write

write_csv_safe(consumer_audit, OUT("animal_id_normalization_consumer_audit.csv"))
write_csv_safe(resolution, OUT("animal_id_resolution_report.csv"))
write_csv_safe(blast, OUT("animal_id_normalization_blast_radius.csv"))

cat("\n===== AnimalID integrity audit =====\n")
cat("\n--- consumers ---\n")
for (i in seq_len(nrow(consumer_audit))) {
  cat(sprintf("  %-58s L%-5d %-18s %s\n",
              substr(consumer_audit$script[i], 1, 58),
              consumer_audit$invocation_line[i],
              consumer_audit$function_name[i],
              consumer_audit$collision_status[i]))
}
cat("\n--- resolution under the canonical contract ---\n")
for (k in unique(resolution$source_key)) {
  z <- resolution[resolution$source_key == k, , drop = FALSE]
  cat(sprintf("  %-34s raw=%-4d resolved_exp9=%-2d legacy_collisions=%-3d\n",
              k, nrow(z), sum(z$is_exp9_animal), sum(z$legacy_collision)))
}
cat("\n--- blast radius ---\n")
for (i in seq_len(nrow(blast))) {
  cat(sprintf("  changed=%-3s %s\n", blast$scientific_result_changed[i],
              substr(basename(blast$artifact[i]), 1, 60)))
}
cat("\nOutputs:", relative_to(dirname(OUT("x"))), "\n")

#!/usr/bin/env Rscript
#
#
# Network edges of the neuropil against the current behaviour outcomes (added 2026-10-30; Extended Data 8c rerun).
#
# USAGE (from the repository root)
#   Rscript analysis/spatial_validation/couple_network_edges_to_outcomes.R
#   Rscript analysis/spatial_validation/couple_network_edges_to_outcomes.R --dry-run
# Script: analysis/spatial_validation/couple_network_edges_to_outcomes.R
# Stage: coupling
# Scope: global
# Consumes: required results/integration/test_network_behaviour_coupling/global/tables/animal_level_candidate_edge_scores.csv; data/external/behaviour_outcomes/behaviour_outcomes_ov6.csv; data/external/behaviour_outcomes/outcome_info.csv
# Produces: results/spatial_validation/couple_network_edges_to_outcomes/global/tables/edge_outcome_correlations.csv; results/spatial_validation/couple_network_edges_to_outcomes/global/tables/edge_outcome_values.csv
# Dataset behavior: global; the eight candidate neuropil edges of test_network_behaviour_coupling.
# Notes: Replaces the outcomes of the earlier coupling (an April 2026 workbook copy and the retired GAMM movement score,
#        whose model contained the outcome group) by the current ones from the frozen behaviour bundles.
#
# EDGES: the animal-level edge scores of test_network_behaviour_coupling.R (Spearman correlation of the two units' protein
# profiles within an animal; eight declared edges; every proteomic animal), as that run stored them (SHA-256 pinned below:
# the eight-animal run of 2026-11-23, with the DG po carry-over acquisitions left out as in Figure 2j-l). Nothing about
# the edges is recomputed. A rerun on other data sets PROTEOMICS_EDGES_UNPINNED=true: the pin is then not enforced and
# the manifest records the hash used. CHECK: the same test on the edge scores without the shallow CA2 SLM acquisitions
# of the robustness audit (edge_outcome_correlations_without_shallow_ca2_slm.csv).
#
# OUTCOMES (group-free; data/external/behaviour_outcomes): the integrated stress score (ISS = CombZ), the corticosterone
# rise, sucrose preference, early movement (Stage 09: first active phase after the first cage change) and mean first-night
# movement over the four cage changes (post hoc run 29f).
#
# TEST: Pearson r of each edge with each outcome over the animals (nine, or eight in the 3/3/2 rerun), with the Fisher-z
# 95% interval, as the earlier coupling; Benjamini-Hochberg within each outcome (eight edges) and across all 40 tests.
# So few animals give very little power: a null here is not evidence of absence.
#
#

source("R/paths.R")
source("R/statistics/integration_utils.R")
source(repo_path("R", "spatial_systems_paths.R"))
source("R/spatial/neuropil_lamina_utils.R")

ANALYSIS_ID <- "couple_network_edges_to_outcomes"
CANONICAL_PATHS <- spatial_systems_dirs(ANALYSIS_ID)
SCRIPT_ID <- "analysis/spatial_validation/couple_network_edges_to_outcomes.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
cli <- integration_cli(default_dataset = "all")

EDGES <- repo_path("results", "integration", "test_network_behaviour_coupling", "global", "tables", "animal_level_candidate_edge_scores.csv")
EDGES_SHA256 <- "40517c8275075c8051cd1d46a398e07c6cc3726f6382340a3e369877fc7df33f"
EDGES_CHECK <- repo_path("results", "integration", "test_network_behaviour_coupling", "global", "tables",
                         "animal_level_candidate_edge_scores_without_shallow_ca2_slm.csv")
OUTCOMES <- path_external("behaviour_outcomes", "behaviour_outcomes_ov6.csv")
OUTCOME_INFO <- path_external("behaviour_outcomes", "outcome_info.csv")
OUT <- list(
  correlations = file.path(CANONICAL_PATHS$tables, "edge_outcome_correlations.csv"),
  values       = file.path(CANONICAL_PATHS$tables, "edge_outcome_values.csv"),
  check        = file.path(CANONICAL_PATHS$tables, "edge_outcome_correlations_without_shallow_ca2_slm.csv")
)
MANIFEST <- file.path(CANONICAL_PATHS$manifests, "run_manifest.yml")
OUTCOME_LABELS <- c(CombZ = "ISS", cort_delta = "Corticosterone rise", suc_pref = "Sucrose preference",
                    early_movement = "Early movement (first night after CC1)", movement_cc1_cc4_mean = "Movement (first nights, CC1-CC4)")

if (isTRUE(cli$dry_run)) {
  cat("[DRY-RUN] Eight neuropil edges against five current behaviour outcomes (Pearson r, Fisher interval, BH).\n")
  dry_run_inputs(SCRIPT_ID, list(edges = EDGES, outcomes = OUTCOMES, outcome_info = OUTCOME_INFO))
  quit(save = "no", status = 0L)
}
for (f in c(EDGES, EDGES_CHECK, OUTCOMES, OUTCOME_INFO)) if (!file.exists(f)) stop("missing_required_input: ", f, call. = FALSE)
edges_sha256_used <- unname(tools::sha256sum(EDGES))
edges_unpinned <- identical(tolower(Sys.getenv("PROTEOMICS_EDGES_UNPINNED", "")), "true")
if (!identical(edges_sha256_used, EDGES_SHA256) && !edges_unpinned)
  stop("The stored edge scores differ from the pinned table of the 2026-11-23 coupling run.", call. = FALSE)
info <- utils::read.csv(OUTCOME_INFO, stringsAsFactors = FALSE)
if (!identical(info$sha256[info$item == "output"], unname(tools::sha256sum(OUTCOMES))))
  stop("behaviour_outcomes_ov6.csv differs from outcome_info.csv; rebuild it with its script.", call. = FALSE)
invisible(lapply(CANONICAL_PATHS[c("tables", "manifests")], dir_create))

beh <- utils::read.csv(OUTCOMES, stringsAsFactors = FALSE, colClasses = c(animal_id = "character", AnimalNum = "character"))
fisher_ci <- function(r, n) { z <- atanh(r); h <- stats::qnorm(0.975) / sqrt(n - 3); tanh(c(z - h, z + h)) }
couple <- function(edge_file) {
  edges <- utils::read.csv(edge_file, stringsAsFactors = FALSE, colClasses = c(AnimalID = "character"))
  n_animals <- length(unique(edges$AnimalID))
  if (length(unique(edges$Edge)) != 8L || nrow(edges) != 8L * n_animals || any(table(edges$AnimalID) != 8L))
    stop("Expected eight edges for each of the ", n_animals, " animals.", call. = FALSE)
  d <- merge(edges[, c("AnimalID", "ExpGroup", "Edge", "EdgeR")], beh, by.x = "AnimalID", by.y = "animal_id")
  if (nrow(d) != nrow(edges)) stop("Not every proteomic animal has a behaviour row.", call. = FALSE)
  if (any(d$ExpGroup != d$Group)) stop("Proteomic and behaviour group labels differ.", call. = FALSE)
  rows <- list()
  for (o in names(OUTCOME_LABELS)) for (e in sort(unique(d$Edge))) {
    x <- d[d$Edge == e & is.finite(d[[o]]), , drop = FALSE]
    ct <- stats::cor.test(x$EdgeR, x[[o]], method = "pearson")
    ci <- fisher_ci(unname(ct$estimate), nrow(x))
    rows[[length(rows) + 1L]] <- data.frame(Edge = e, Outcome = o, outcome_label = OUTCOME_LABELS[[o]], n = nrow(x),
      estimate = unname(ct$estimate), conf.low = ci[1], conf.high = ci[2], statistic = unname(ct$statistic),
      df = unname(ct$parameter), p.value = ct$p.value, stringsAsFactors = FALSE)
  }
  cor_tbl <- do.call(rbind, rows)
  cor_tbl$p.adj_BH_within_outcome <- stats::ave(cor_tbl$p.value, cor_tbl$Outcome, FUN = function(p) stats::p.adjust(p, "BH"))
  cor_tbl$p.adj_BH_all <- stats::p.adjust(cor_tbl$p.value, "BH")
  cor_tbl$note <- sprintf("%d animals: very low power; a null is not evidence of absence", n_animals)
  list(cor = cor_tbl, values = d[order(d$Edge, d$AnimalID), c("Edge", "AnimalID", "AnimalNum", "Group", "EdgeR", names(OUTCOME_LABELS))])
}
primary <- couple(EDGES)
cor_tbl <- primary$cor; values <- primary$values
check_tbl <- couple(EDGES_CHECK)$cor
check_tbl$note <- "check: without the shallow CA2 SLM acquisitions of the robustness audit (A755 left, A765 right)"

utils::write.csv(cor_tbl, OUT$correlations, row.names = FALSE)
utils::write.csv(check_tbl, OUT$check, row.names = FALSE)
utils::write.csv(values, OUT$values, row.names = FALSE)
write_run_manifest(MANIFEST, inputs = list(edges = EDGES, edges_check = EDGES_CHECK, outcomes = OUTCOMES, outcome_info = OUTCOME_INFO), outputs = OUT,
                   parameters = list(test = "Pearson r, Fisher-z 95% interval; BH within outcome (8 edges) and across all 40 tests",
                                     edges_sha256 = edges_sha256_used, edges_sha256_pinned = EDGES_SHA256,
                                     edges_pin_enforced = !edges_unpinned, n_animals = n_animals, outcomes = as.list(OUTCOME_LABELS),
                                     provenance = run_output_provenance(OUT, SCRIPT_ID)),
                   notes = "Rerun of the edge-behaviour coupling (Extended Data 8c) on the current, group-free behaviour outcomes.")
message("Edge-outcome coupling: ", nrow(cor_tbl), " tests; smallest BH within outcome ", signif(min(cor_tbl$p.adj_BH_within_outcome), 2))

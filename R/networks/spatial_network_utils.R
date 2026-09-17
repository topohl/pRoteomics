# Reusable spatial-network helpers.

# --- where the canonical network object lives ------------------------------
#
# Phase 6G.2 moved this object from data/processed/07_spatial_networks/... to
# results/spatial_networks/build_spatial_networks/<dataset>/models/. It is a
# result and not a work intermediate because eleven consumers read it, across
# integration, wgcna, spatial_validation and the figure exporter, and the
# lifecycle contract forbids citing work/.
#
# Five scripts had byte-identical copies of this resolution, differing only in
# the name of their dataset variable. One copy means one place to get the
# precedence right.
#
# Precedence is deliberate and ordered:
#   1. an explicit override, for a caller that knows better;
#   2. the normalized canonical location;
#   3. the historical locations, in the order the historical writers used them.
#
# The normalized path comes first so that once a run produces it, it is
# authoritative. A fallback that won ties would leave the historical object
# authoritative forever, which is how a migration quietly fails to take effect.
spatial_network_object_candidates <- function(dataset, spatial_unit = NULL) {
  unit <- if (is.null(spatial_unit) || !length(spatial_unit)) NULL else {
    u <- as.character(spatial_unit)[1]
    if (is.na(u) || !nzchar(u)) NULL else u
  }
  OBJ <- "network_spatial_relations_objects.rds"
  c(
    normalized = canonical_result_path("spatial_networks", "build_spatial_networks",
                                       dataset, "models", unit, OBJ),
    legacy_scoped = path_processed("07_spatial_networks", "network_spatial_relations",
                                   dataset, unit, OBJ),
    legacy_unscoped = path_processed("07_spatial_networks", "network_spatial_relations",
                                     OBJ)
  )
}

# --- where the differential bootstrap tables live --------------------------
#
# The same normalized-first resolution, for the tables that
# test_differential_network_stability writes and
# render_differential_network_figures reads.
#
# This pairing was already broken before Phase 6G.2, and the migration is what
# surfaced it. The writer wrote its tables directly into the substep directory;
# the reader appended "01_Tables" and then hard-stopped when the file was
# missing. That subdirectory does not exist on disk, so the renderer could not
# complete a real run. It is registered required: false, which is why nothing
# reported it. pipeline.yml also declared the phantom 01_Tables path.
#
# 01_Tables is layout chronology, not scientific meaning, so it is not
# reproduced in the normalized namespace. It stays as an explicit compatibility
# candidate, last, in case an older tree has it.
bootstrap_differential_tables_candidates <- function(dataset) {
  legacy <- path_results("tables", "07_spatial_networks",
                         "bootstrap_differential_network_stability")
  c(
    normalized = canonical_result_path("spatial_networks",
                                       "test_differential_network_stability",
                                       dataset, "tables"),
    legacy_flat = legacy,
    legacy_01_tables = file.path(legacy, "01_Tables")
  )
}

# A directory counts only if it actually holds the required table, otherwise an
# empty normalized directory created by a dry run would shadow the real data.
resolve_bootstrap_differential_tables <- function(
    dataset,
    required = "bootstrap_differential_edge_stability_summary.csv",
    env = "PROTEOMICS_BOOTSTRAP_DIFFERENTIAL_TABLES") {
  override <- Sys.getenv(env, unset = "")
  if (nzchar(override)) {
    return(normalizePath(override, winslash = "/", mustWork = FALSE))
  }
  cand <- bootstrap_differential_tables_candidates(dataset)
  ok <- vapply(cand, function(d) file.exists(file.path(d, required)), logical(1))
  if (any(ok)) return(unname(cand[ok][1]))
  unname(cand[["normalized"]])
}

resolve_spatial_network_object <- function(dataset, spatial_unit = NULL,
                                           env = "PROTEOMICS_SPATIAL_NETWORK_OBJECT") {
  override <- Sys.getenv(env, unset = "")
  if (nzchar(override)) {
    return(normalizePath(override, winslash = "/", mustWork = FALSE))
  }
  cand <- spatial_network_object_candidates(dataset, spatial_unit)
  hit <- cand[file.exists(cand)]
  if (length(hit)) return(unname(hit[1]))
  ## nothing on disk yet: name the canonical destination, so a failure message
  ## points at where the object is supposed to be rather than where it used to be
  unname(cand[["normalized"]])
}

make_edge_id <- function(source, target, sep = "--") {
  source <- as.character(source)
  target <- as.character(target)
  paste(pmin(source, target), pmax(source, target), sep = sep)
}

network_interpretation_strength <- function(fdr = NA_real_, stability_score = NA_real_, n_animals = NA_integer_) {
  fdr <- suppressWarnings(as.numeric(fdr))
  stability_score <- suppressWarnings(as.numeric(stability_score))
  n_animals <- suppressWarnings(as.numeric(n_animals))
  if (!is.na(fdr) && fdr <= 0.05 && !is.na(stability_score) && stability_score >= 0.7 && (is.na(n_animals) || n_animals >= 6)) {
    return("strong")
  }
  if ((!is.na(fdr) && fdr <= 0.10) || (!is.na(stability_score) && stability_score >= 0.6)) return("moderate")
  "exploratory"
}

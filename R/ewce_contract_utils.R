# Pure contracts for canonical and sensitivity-only EWCE execution.

ewce_resolve_run_contract <- function(dataset, analysis_unit = "animal", branch = "") {
  dataset <- trimws(as.character(dataset))
  if (length(dataset) != 1L || is.na(dataset) || !nzchar(dataset)) {
    stop("EWCE dataset must be one non-empty value.", call. = FALSE)
  }
  analysis_unit <- tolower(trimws(as.character(analysis_unit)))
  if (length(analysis_unit) != 1L || !analysis_unit %in% c("sample", "animal")) {
    stop("PROTEOMICS_EWCE_ANALYSIS_UNIT must be one of: sample, animal.", call. = FALSE)
  }
  branch <- trimws(as.character(branch))
  if (length(branch) != 1L || is.na(branch)) branch <- ""
  if (nzchar(branch) && !grepl("^[A-Za-z0-9_-]+$", branch)) {
    stop("PROTEOMICS_EWCE_BRANCH must match ^[A-Za-z0-9_-]+$.", call. = FALSE)
  }
  if (identical(analysis_unit, "sample") && !nzchar(branch)) {
    stop(
      "Sample-level EWCE is legacy/sensitivity-only and requires PROTEOMICS_EWCE_BRANCH. ",
      "Canonical EWCE uses animal-level biological units.",
      call. = FALSE
    )
  }

  canonical <- !nzchar(branch)
  list(
    dataset = dataset,
    analysis_unit = analysis_unit,
    branch = if (canonical) "canonical" else branch,
    canonical = canonical,
    substep_id = if (canonical) {
      file.path("EWCE_E9", dataset)
    } else {
      file.path("EWCE_E9_comparison", branch, dataset)
    }
  )
}

ewce_gene_set_sha256 <- function(genes) {
  digest::digest(sort(unique(stats::na.omit(as.character(genes)))), algo = "sha256")
}

ewce_animal_cache_key <- function(contract_version,
                                  input_matrix_sha256,
                                  dataset,
                                  target,
                                  top_n,
                                  annot_level,
                                  reps,
                                  seed,
                                  hits,
                                  background,
                                  ctd_annotation_genes) {
  digest::digest(list(
    analysis_contract_version = contract_version,
    input_matrix_sha256 = input_matrix_sha256,
    dataset = dataset,
    target = target,
    top_n = top_n,
    annot_level = annot_level,
    reps = reps,
    seed = seed,
    hit_gene_set_sha256 = ewce_gene_set_sha256(hits),
    background_gene_set_sha256 = ewce_gene_set_sha256(background),
    ctd_annotation_gene_set_sha256 = ewce_gene_set_sha256(ctd_annotation_genes)
  ), algo = "xxhash64")
}

ewce_is_sample_level <- function(analysis_unit) {
  identical(tolower(trimws(as.character(analysis_unit))), "sample")
}

ewce_legacy_cache_fallback_allowed <- function(analysis_unit) {
  ewce_is_sample_level(analysis_unit)
}

ewce_cache_accounting_table <- function(target_runs, cache_events) {
  target_runs <- as.character(target_runs)
  if (!length(target_runs) || anyNA(target_runs) || any(!nzchar(target_runs)) || anyDuplicated(target_runs)) {
    stop("EWCE cache accounting requires unique, non-empty target-run IDs.", call. = FALSE)
  }
  required <- c("TargetRun", "CacheEvent")
  if (!is.data.frame(cache_events) || !all(required %in% names(cache_events))) {
    stop("EWCE cache accounting requires TargetRun and CacheEvent columns.", call. = FALSE)
  }
  events <- cache_events[, required, drop = FALSE]
  events$TargetRun <- as.character(events$TargetRun)
  events$CacheEvent <- as.character(events$CacheEvent)
  valid_events <- c("cache_hit", "legacy_cache_fallback", "computed")
  if (nrow(events) != length(target_runs) || anyDuplicated(events$TargetRun) ||
      !setequal(events$TargetRun, target_runs) || any(!events$CacheEvent %in% valid_events)) {
    stop("EWCE cache accounting events do not match the expected target runs.", call. = FALSE)
  }
  data.frame(
    expected_target_runs = length(target_runs),
    cache_hits = sum(events$CacheEvent == "cache_hit"),
    cache_misses = sum(events$CacheEvent != "cache_hit"),
    new_bootstrap_computations = sum(events$CacheEvent == "computed"),
    legacy_cache_fallback_count = sum(events$CacheEvent == "legacy_cache_fallback"),
    stringsAsFactors = FALSE
  )
}

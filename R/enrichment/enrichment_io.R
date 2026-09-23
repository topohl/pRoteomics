# Shared enrichment IO and conservative biological-program mapping helpers.

if (!exists("repo_path", mode = "function")) {
  paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
  source(paths_file)
}
if (!exists("input_addressability", mode = "function")) {
  source(repo_path("R", "paths.R"))
}
if (!exists("safe_name", mode = "function")) {
  source(repo_path("R", "validation_utils.R"))
}
## the manifest paths below resolve normalized-first through this helper
if (!exists("resolve_differential_abundance_state", mode = "function")) {
  source(repo_path("R", "differential_abundance_paths.R"))
}

canonical_clusterprofiler_manifest_contract_version <- function() {
  "clusterProfiler_manifest_v3_term_gene_provenance"
}

canonical_comparego_manifest_contract_version <- function() {
  "compareGO_manifest_v3_term_gene_provenance"
}

canonical_comparego_result_types <- function() "GSEA_GO"

clusterprofiler_manifest_columns <- function() c(
  "analysis_id", "dataset", "run_id", "ontology", "result_type", "contrast", "comparison",
  "route_category", "route_unit", "condition", "direction", "simplified", "plot_suffix",
  "used_for_plot", "input_gene_file", "gene_input_file", "input_hash",
  "collapsed_gene_input_file", "collapsed_gene_provenance_file", "term_gene_provenance_file",
  "enrichment_contract_version", "gene_annotation_contract_version", "protein_group_contract_version",
  "gene_mapping_policy", "primary_gene_level_eligibility_rule", "ambiguous_group_policy",
  "duplicate_gene_collapse_rule", "rank_statistic_column", "rank_statistic_type",
  "rank_statistic_fallback_used", "ORA_direction", "universe_definition", "config_file",
  "config_hash", "output_table", "output_plot", "n_genes", "n_terms", "analysis_status",
  "empty_result", "error_message", "checkpoint_status", "created_at"
)

comparego_manifest_columns <- function() c(
  clusterprofiler_manifest_columns(), "input_manifest", "comparego_contract_version",
  "comparego_analysis_status", "term_comparison_file", "term_gene_provenance_output_file",
  "analysis_status_summary_file"
)

clusterprofiler_compatibility_fallback_enabled <- function(strict_mode, fallback_requested) {
  !isTRUE(strict_mode) && isTRUE(fallback_requested)
}

clusterprofiler_compatibility_fallback_audit_paths <- function(audit_root, enabled = FALSE) {
  if (!isTRUE(enabled)) return(character(0))
  file.path(
    audit_root,
    c(
      "compatibility_fallback_accession_annotation_audit.csv",
      "compatibility_fallback_protein_group_annotation_audit.csv"
    )
  )
}

clusterprofiler_output_path_audit <- function(paths, safe_limit = 240L) {
  paths <- unique(as.character(paths))
  paths <- paths[!is.na(paths) & nzchar(paths)]
  normalized <- normalizePath(paths, winslash = "/", mustWork = FALSE)
  out <- data.frame(
    path = normalized,
    path_length = nchar(normalized, type = "chars"),
    safe_limit = as.integer(safe_limit),
    within_safe_limit = nchar(normalized, type = "chars") <= as.integer(safe_limit),
    stringsAsFactors = FALSE
  )
  out[order(out$path_length, decreasing = TRUE, out$path, method = "radix"), , drop = FALSE]
}

validate_clusterprofiler_output_path_lengths <- function(paths, safe_limit = 240L) {
  audit <- clusterprofiler_output_path_audit(paths, safe_limit = safe_limit)
  excessive <- audit[!audit$within_safe_limit, , drop = FALSE]
  if (nrow(excessive)) {
    stop(
      "Expected clusterProfiler output path exceeds the safe Windows/R path limit of ",
      as.integer(safe_limit), " characters (maximum expected length: ", max(excessive$path_length), "). ",
      "Run the repository from a shorter project root such as P:\\ before launching workers. ",
      "Longest path: ", excessive$path[[1]],
      call. = FALSE
    )
  }
  audit
}

clusterprofiler_worker_error <- function(result, fallback) {
  if (!is.list(result) || is.null(result$error) || !length(result$error) ||
      is.na(result$error[[1]]) || !nzchar(as.character(result$error[[1]]))) {
    return(fallback)
  }
  as.character(result$error[[1]])
}

assess_clusterprofiler_worker_result <- function(result, expected_comparison) {
  failed <- function(message, worker_status = "MALFORMED", manifest_has_failure = FALSE) {
    data.frame(
      comparison = as.character(expected_comparison), worker_status = worker_status,
      analysis_status = "failed", error = as.character(message),
      manifest_has_failure = isTRUE(manifest_has_failure), stringsAsFactors = FALSE
    )
  }
  if (!is.list(result)) return(failed("Worker returned a missing or malformed result object."))
  required <- c("status", "comparison", "manifest")
  if (length(setdiff(required, names(result)))) {
    return(failed("Worker result object is missing required fields: status, comparison, manifest."))
  }
  worker_status <- as.character(result$status)[1]
  comparison <- as.character(result$comparison)[1]
  if (is.na(worker_status) || !nzchar(worker_status) || is.na(comparison) || !nzchar(comparison)) {
    return(failed("Worker result contains an empty status or comparison."))
  }
  if (!identical(comparison, as.character(expected_comparison))) {
    return(failed(
      paste0("Worker comparison identity mismatch: expected ", expected_comparison, ", received ", comparison, "."),
      worker_status = worker_status
    ))
  }
  manifest <- result$manifest
  manifest_is_table <- is.data.frame(manifest)
  manifest_has_failure <- manifest_is_table && all(c("result_type", "analysis_status") %in% names(manifest)) &&
    any(manifest$result_type == "GSEA_GO" & manifest$analysis_status == "failed", na.rm = TRUE)
  if (!worker_status %in% c("SUCCESS", "SKIPPED")) {
    return(failed(
      clusterprofiler_worker_error(result, paste0("Worker returned status ", worker_status, ".")),
      worker_status = worker_status, manifest_has_failure = manifest_has_failure
    ))
  }
  if (!manifest_is_table || !all(c("result_type", "analysis_status", "n_terms") %in% names(manifest))) {
    return(failed("Successful worker returned a missing or malformed manifest.", worker_status = worker_status))
  }
  primary <- manifest[manifest$result_type == "GSEA_GO", , drop = FALSE]
  if (nrow(primary) != 1L) {
    return(failed("Successful worker must return exactly one primary GSEA_GO manifest row.", worker_status = worker_status))
  }
  primary_status <- as.character(primary$analysis_status[[1]])
  n_terms <- suppressWarnings(as.integer(primary$n_terms[[1]]))
  if (identical(primary_status, "success_with_terms") && !is.na(n_terms) && n_terms > 0L) {
    category <- "success_with_terms"
  } else if (identical(primary_status, "success_zero_terms") && identical(n_terms, 0L)) {
    category <- "success_zero_terms"
  } else {
    return(failed(
      paste0("Primary GSEA_GO manifest has inconsistent analysis_status/n_terms: ",
        primary_status, "/", ifelse(is.na(n_terms), "NA", n_terms), "."),
      worker_status = worker_status, manifest_has_failure = manifest_has_failure
    ))
  }
  data.frame(
    comparison = comparison, worker_status = worker_status, analysis_status = category,
    error = NA_character_, manifest_has_failure = FALSE, stringsAsFactors = FALSE
  )
}

assess_clusterprofiler_worker_results <- function(results, expected_comparisons) {
  if (!is.list(results)) results <- list(results)
  expected_comparisons <- as.character(expected_comparisons)
  n_rows <- max(length(results), length(expected_comparisons))
  if (!n_rows) {
    return(data.frame(
      comparison = character(), worker_status = character(), analysis_status = character(),
      error = character(), manifest_has_failure = logical(), stringsAsFactors = FALSE
    ))
  }
  rows <- lapply(seq_len(n_rows), function(i) {
    if (i > length(expected_comparisons)) {
      return(data.frame(
        comparison = paste0("unexpected_worker_result_", i), worker_status = "MALFORMED",
        analysis_status = "failed", error = "Received an unexpected extra worker result object.",
        manifest_has_failure = FALSE, stringsAsFactors = FALSE
      ))
    }
    result <- if (i <= length(results)) results[[i]] else NULL
    assess_clusterprofiler_worker_result(result, expected_comparisons[[i]])
  })
  do.call(rbind, rows)
}

clusterprofiler_master_status_counts <- function(assessment) {
  statuses <- c("success_with_terms", "success_zero_terms", "failed")
  counts <- setNames(integer(length(statuses)), statuses)
  observed <- table(factor(assessment$analysis_status, levels = statuses))
  counts[] <- as.integer(observed)
  counts
}

clusterprofiler_master_exit_status <- function(assessment) {
  if (any(assessment$analysis_status == "failed")) 1L else 0L
}

clusterprofiler_master_summary_lines <- function(assessment) {
  counts <- clusterprofiler_master_status_counts(assessment)
  final <- if (counts[["failed"]] > 0L) {
    paste0("RUN FAILED: ", counts[["failed"]], " comparison(s) failed.")
  } else {
    "ALL COMPARISONS COMPLETED SUCCESSFULLY."
  }
  c(
    paste0("success_with_terms: ", counts[["success_with_terms"]]),
    paste0("success_zero_terms: ", counts[["success_zero_terms"]]),
    paste0("failed: ", counts[["failed"]]),
    final
  )
}

write_csv_strict <- function(x, path, label = "CSV output") {
  parent <- dirname(path)
  if (!dir.exists(parent)) {
    stop(label, " parent directory does not exist: ", parent, call. = FALSE)
  }
  tryCatch(
    utils::write.csv(x, path, row.names = FALSE),
    error = function(e) stop("Failed to write ", label, " to ", path, ": ", conditionMessage(e), call. = FALSE)
  )
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size <= 0) {
    stop("Failed to verify written ", label, ": ", path, call. = FALSE)
  }
  written <- tryCatch(
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) stop("Failed to read back written ", label, " at ", path, ": ", conditionMessage(e), call. = FALSE)
  )
  if (!identical(names(written), names(x)) || nrow(written) != nrow(x)) {
    stop("Written ", label, " failed row/column verification: ", path, call. = FALSE)
  }
  invisible(path)
}

term_gene_provenance_columns <- function() c(
  "dataset", "comparison", "result_type", "ontology", "term_id", "term_description",
  "official_gene_symbol", "official_entrez_id", "ProteinGroupID", "member_accessions",
  "protein_group_gene_annotation_status", "gene_level_claim_allowed", "rank_statistic",
  "core_enrichment_member", "enrichment_contract_version", "gene_annotation_contract_version"
)

validate_term_gene_provenance_contract <- function(x, strict = TRUE) {
  missing <- setdiff(term_gene_provenance_columns(), names(x))
  if (length(missing)) {
    stop("Term-gene provenance is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!is.character(x$official_entrez_id)) {
    stop("Term-gene provenance official_entrez_id must remain character.", call. = FALSE)
  }
  if (isTRUE(strict) && any(x$gene_level_claim_allowed %in% FALSE, na.rm = TRUE)) {
    stop("Ineligible protein groups are not permitted in strict term-gene provenance.", call. = FALSE)
  }
  if (isTRUE(strict) && any(as.character(x$protein_group_gene_annotation_status) != "concordant_official_gene")) {
    stop("Ambiguous protein-group annotations are not permitted in strict term-gene provenance.", call. = FALSE)
  }
  identity_columns <- c("dataset", "comparison", "result_type", "ontology", "term_id",
    "official_gene_symbol", "ProteinGroupID")
  if (nrow(x) && any(vapply(x[identity_columns], function(value) {
      any(is.na(value) | !nzchar(trimws(as.character(value))))
    }, logical(1)))) {
    stop("Term-gene provenance contains missing canonical identity values.", call. = FALSE)
  }
  if (nrow(x) && anyDuplicated(x[identity_columns])) {
    stop("Term-gene provenance contains duplicate canonical term/gene/ProteinGroupID rows.", call. = FALSE)
  }
  invisible(TRUE)
}

validate_gsea_result_table_contract <- function(x, context = "GSEA result") {
  required <- c("ID", "Description", "NES", "p.adjust", "setSize", "core_enrichment")
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop(context, " is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

validate_clusterprofiler_manifest_contract <- function(manifest, strict = TRUE, require_files = TRUE,
                                                       file_columns = clusterprofiler_runtime_required_fields()) {
  required <- c(
    "dataset", "comparison", "result_type", "ontology", "analysis_status", "n_terms",
    "output_table", "collapsed_gene_input_file", "collapsed_gene_provenance_file",
    "term_gene_provenance_file", "enrichment_contract_version",
    "gene_annotation_contract_version"
  )
  missing <- setdiff(required, names(manifest))
  if (length(missing)) {
    stop("clusterProfiler manifest is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  expected_version <- canonical_clusterprofiler_manifest_contract_version()
  if (isTRUE(strict) && any(is.na(manifest$enrichment_contract_version) |
      manifest$enrichment_contract_version != expected_version)) {
    stop("Stale clusterProfiler manifest contract; expected ", expected_version, ".", call. = FALSE)
  }
  allowed_status <- c("success_with_terms", "success_zero_terms", "failed")
  invalid_status <- setdiff(unique(as.character(manifest$analysis_status)), allowed_status)
  if (length(invalid_status)) {
    stop("clusterProfiler manifest has unsupported analysis_status: ", paste(invalid_status, collapse = ", "), call. = FALSE)
  }
  supported_results <- c("GSEA_GO", "GSEA_KEGG")
  invalid_results <- setdiff(unique(as.character(manifest$result_type)), supported_results)
  if (length(invalid_results)) {
    stop("clusterProfiler manifest has unsupported result_type: ", paste(invalid_results, collapse = ", "), call. = FALSE)
  }
  success <- manifest$analysis_status %in% c("success_with_terms", "success_zero_terms")
  zero <- manifest$analysis_status == "success_zero_terms"
  with_terms <- manifest$analysis_status == "success_with_terms"
  if (any(zero & (is.na(manifest$n_terms) | manifest$n_terms != 0))) {
    stop("success_zero_terms rows must record n_terms = 0.", call. = FALSE)
  }
  if (any(with_terms & (is.na(manifest$n_terms) | manifest$n_terms <= 0))) {
    stop("success_with_terms rows must record n_terms > 0.", call. = FALSE)
  }
  if (isTRUE(require_files) && any(success)) {
    ## Runtime validity is gated on the fields runtime actually reads. A
    ## provenance-only field stays in the manifest and keeps its declared
    ## value, but nothing opens it, so requiring it to be openable would force
    ## 108 pointless staged copies and would fail for a reason no consumer
    ## cares about. Which fields are which is recorded in
    ## CLUSTERPROFILER_MANIFEST_PATH_FIELDS, derived from a grep for reads.
    ## Semantics are otherwise untouched: still all-or-nothing over the gated
    ## fields, and a single unusable one still invalidates the manifest.
    for (column in intersect(file_columns, names(manifest))) {
      paths <- as.character(manifest[[column]][success])
      ## The contract is unchanged and still all-or-nothing: every declared
      ## file of a successful row must be usable, or the manifest is invalid.
      ## What changed is that the failure is named. file.exists() reports FALSE
      ## for an unmounted declared root, for a path past the 260-character
      ## limit, and for a file that is simply not there, and those demand
      ## different responses from whoever reads this message.
      status <- input_addressability(paths)
      if (any(status != INPUT_STATUS_PRESENT)) {
        stop("Successful clusterProfiler manifest row references unusable ", column, ": ",
          describe_input_status_failures(paths, status), call. = FALSE)
      }
    }
  }
  invisible(TRUE)
}

validate_comparego_manifest_contract <- function(manifest, require_files = TRUE) {
  required <- c(
    "dataset", "comparison", "result_type", "ontology", "analysis_status",
    "comparego_analysis_status", "input_manifest", "term_comparison_file",
    "term_gene_provenance_output_file", "analysis_status_summary_file",
    "enrichment_contract_version", "comparego_contract_version"
  )
  missing <- setdiff(required, names(manifest))
  if (length(missing)) {
    stop("compareGO manifest is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (any(manifest$comparego_contract_version != canonical_comparego_manifest_contract_version())) {
    stop("Stale compareGO manifest contract.", call. = FALSE)
  }
  if (any(!manifest$result_type %in% canonical_comparego_result_types())) {
    stop("compareGO manifest contains unsupported result_type.", call. = FALSE)
  }
  allowed_actions <- c("included", "completed_zero_terms", "recorded_failed")
  if (any(!manifest$comparego_analysis_status %in% allowed_actions)) {
    stop("compareGO manifest contains unsupported comparego_analysis_status.", call. = FALSE)
  }
  if (isTRUE(require_files)) {
    file_columns <- c("input_manifest", "term_comparison_file", "term_gene_provenance_output_file", "analysis_status_summary_file")
    for (column in file_columns) {
      paths <- unique(as.character(manifest[[column]]))
      ## Same all-or-nothing contract as the clusterProfiler manifest above, and
      ## now the same vocabulary. file.exists() returns FALSE for an unmounted
      ## declared root, for a path at or past the character limit, and for a file
      ## that is genuinely not there; only the last is fixed by re-running the
      ## producer, so a message that calls all three "missing" misdirects whoever
      ## reads it. input_addressability() already classifies NA and empty as
      ## absent, so the two guards it replaces are subsumed, not dropped.
      status <- input_addressability(paths)
      if (any(status != INPUT_STATUS_PRESENT)) {
        stop("compareGO manifest references unusable ", column, ": ",
          describe_input_status_failures(paths, status), call. = FALSE)
      }
    }
  }
  invisible(TRUE)
}

read_csv_contract <- function(path, character_columns = character()) {
  header <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, nrows = 0L)
  character_columns <- intersect(character_columns, names(header))
  x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
    colClasses = if (length(character_columns)) c(setNames(rep("character", length(character_columns)), character_columns)) else NA)
  x
}

# Phase 6G.4: these two manifests moved into the normalized namespace as
# models/, because downstream analyses in enrichment, integration and wgcna
# read them as state rather than as provenance.
#
# Every consumer resolves them through these two functions and nowhere else,
# so making the resolution normalized-first here repoints all of them at once.
# Until run_clusterprofiler_enrichment / compare_go_enrichment are actually
# rerun the normalized copy does not exist and the historical one answers,
# which is why current behaviour is unchanged.
canonical_clusterprofiler_manifest_path <- function(dataset, repository_root = repo_path()) {
  resolve_differential_abundance_state(
    "clusterProfiler_manifest.csv",
    owner = "run_clusterprofiler_enrichment", dataset = dataset,
    legacy_substep = "clusterProfiler", repository_root = repository_root)
}

canonical_comparego_manifest_path <- function(dataset, repository_root = repo_path()) {
  resolve_differential_abundance_state(
    "compareGO_input_manifest.csv",
    owner = "compare_go_enrichment", dataset = dataset,
    legacy_substep = "compareGO", repository_root = repository_root)
}

resolve_repository_contract_path <- function(path, repository_root = repo_path()) {
  path <- as.character(path)
  out <- rep(NA_character_, length(path))
  usable <- !is.na(path) & nzchar(trimws(path))
  if (!any(usable)) return(out)
  values <- path[usable]
  absolute <- grepl("^[A-Za-z]:[/\\\\]", values) | grepl("^[/\\\\]{1,2}", values)
  values[!absolute] <- file.path(repository_root, values[!absolute])
  out[usable] <- normalizePath(values, winslash = "/", mustWork = FALSE)
  out
}

# --- runtime manifest resolution --------------------------------------------
#
# The stored manifests are immutable provenance. Every path in them is recorded
# under a substituted P:/ root that describes the machine the run happened on,
# and that root is not mounted here. Re-anchoring the stored suffix on this
# repository root is necessary but not sufficient: the stored paths are short
# (max 207) and this root is 69 characters, so re-anchoring adds ~66 and pushes
# part of the set past the wall R can open.
#
# So the declared path stays exactly as recorded, and a runtime path is derived
# beside it. Only an actively read field whose re-anchored path crosses the
# wall is staged.
#
# Which fields are actually read was derived mechanically from a grep for
# content reads. A field whose value is only carried into another table as a
# string is provenance: nothing opens it, so runtime validity must not depend
# on it being openable. collapsed_gene_provenance_file is the case that
# matters - its re-anchored paths reach 273 characters, and staging 108 cells
# nothing reads would be pure waste.
## input_gene_file is read by analysis/wgcna/compare_module_enrichment_overlap.R,
## not by this file, but it is declared under the same unmounted root and so
## needs the same re-anchoring. Its paths top out at 168 characters, so it is
## runtime-required and never a staging candidate.
## The other three provenance-only fields are listed so the runtime ledger
## accounts for every path-bearing column rather than going quiet about four of
## them. gene_mapping_policy is deliberately absent: it contains prose
## ("SYMBOL/ENTREZ") and is not a path at all.
CLUSTERPROFILER_MANIFEST_PATH_FIELDS <- c(
  output_table                   = "runtime_required",
  collapsed_gene_input_file      = "runtime_required",
  term_gene_provenance_file      = "runtime_required",
  input_gene_file                = "runtime_required",
  collapsed_gene_provenance_file = "provenance_only",
  config_file                    = "provenance_only",
  gene_input_file                = "provenance_only",
  output_plot                    = "provenance_only"
)

clusterprofiler_runtime_required_fields <- function() {
  names(CLUSTERPROFILER_MANIFEST_PATH_FIELDS)[
    CLUSTERPROFILER_MANIFEST_PATH_FIELDS == "runtime_required"]
}

clusterprofiler_provenance_only_fields <- function() {
  names(CLUSTERPROFILER_MANIFEST_PATH_FIELDS)[
    CLUSTERPROFILER_MANIFEST_PATH_FIELDS == "provenance_only"]
}

# Strip whatever root a stored path declares, leaving the repository-relative
# remainder that both environments agree on.
manifest_declared_suffix <- function(path) {
  p <- gsub("\\\\", "/", as.character(path))
  p <- sub("^[A-Za-z]:/+", "", p)
  sub("^/+", "", p)
}

# Resolve one vector of declared paths to usable runtime paths.
#
# Precedence, deliberately: a declared path that already works is used as-is;
# otherwise the re-anchored candidate is used if R can open it; otherwise, if
# the candidate exists physically but is too long, a byte-identical copy is
# staged; otherwise the failure is left truthful. A genuinely absent file is
# never staged to make validation pass.
resolve_runtime_paths <- function(declared, repository_root = repo_path(),
                                  stage_root = NULL, stage = TRUE) {
  declared <- as.character(declared)
  n <- length(declared)
  out <- data.frame(
    declared_path = declared,
    runtime_path = declared,
    runtime_resolution = rep(RUNTIME_RESOLUTION_UNRESOLVED, n),
    declared_status = input_addressability(declared),
    candidate_path = NA_character_,
    candidate_length = NA_integer_,
    source_sha256 = NA_character_,
    staged_sha256 = NA_character_,
    stringsAsFactors = FALSE)
  if (!n) return(out)

  usable <- !is.na(declared) & nzchar(declared)
  ## a declared path that is simply usable here needs nothing done to it
  direct <- usable & out$declared_status == INPUT_STATUS_PRESENT
  out$runtime_resolution[direct] <- RUNTIME_RESOLUTION_DIRECT

  todo <- which(usable & !direct)
  if (!length(todo)) return(out)

  cand <- gsub("\\\\", "/", file.path(repository_root, manifest_declared_suffix(declared[todo])))
  out$candidate_path[todo] <- cand
  out$candidate_length[todo] <- path_length_chars(cand)

  short <- out$candidate_length[todo] < PATH_LENGTH_WALL
  ## Below the wall R can answer for itself, which keeps fixtures fast and free
  ## of any dependency on an external shell.
  if (any(short)) {
    i <- todo[short]
    hit <- file.exists(out$candidate_path[i])
    out$runtime_path[i[hit]] <- out$candidate_path[i[hit]]
    out$runtime_resolution[i[hit]] <- RUNTIME_RESOLUTION_REBASED
  }

  long <- which(!short)
  if (!length(long)) return(out)
  i <- todo[long]
  ## At or beyond the wall only an extended-length-aware runtime can tell us
  ## whether the file is there, so file.exists() must not be consulted.
  facts <- os_path_facts(out$candidate_path[i])
  exists <- !is.na(facts$exists) & facts$exists
  out$source_sha256[i] <- facts$sha256
  if (!isTRUE(stage) || !any(exists)) return(out)

  root <- stage_root %||% stop("resolve_runtime_paths needs a stage_root", call. = FALSE)
  j <- i[exists]
  dest <- staged_destination(out$declared_path[j], root)
  staged <- stage_addressable_copies(out$candidate_path[j], dest,
                                     expected_sha256 = facts$sha256[exists])
  out$staged_sha256[j] <- staged$staged_sha256
  ## Only a verified staged copy becomes the runtime path. A copy that failed,
  ## or a destination holding different bytes, stays unresolved so the
  ## validator can refuse it rather than reading something unexpected.
  good <- j[staged$ok]
  out$runtime_path[good] <- staged$staged_path[staged$ok]
  out$runtime_resolution[good] <- RUNTIME_RESOLUTION_STAGED
  out
}

# Resolve a whole manifest to its runtime view.
#
# Returns the manifest with runtime paths in the runtime-required columns and
# the declared values preserved in <field>_declared, plus a per-cell resolution
# table. The manifest FILE is not touched; this is an in-memory view.
clusterprofiler_manifest_runtime_resolution <- function(manifest, dataset = NULL,
                                                        repository_root = repo_path(),
                                                        stage = TRUE) {
  rows <- list()
  row_dataset <- if ("dataset" %in% names(manifest)) as.character(manifest$dataset) else
    rep(dataset %||% "global", nrow(manifest))
  row_dataset[is.na(row_dataset) | !nzchar(row_dataset)] <- dataset %||% "global"

  for (field in names(CLUSTERPROFILER_MANIFEST_PATH_FIELDS)) {
    if (!field %in% names(manifest)) next
    kind <- CLUSTERPROFILER_MANIFEST_PATH_FIELDS[[field]]
    declared <- as.character(manifest[[field]])

    if (identical(kind, "provenance_only")) {
      ## Recorded, deliberately not resolved and never staged. Its declared
      ## status is still reported so the ledger stays truthful about it.
      rows[[length(rows) + 1L]] <- data.frame(
        dataset = row_dataset, field = field, field_kind = kind,
        declared_path = declared, runtime_path = declared,
        runtime_resolution = RUNTIME_RESOLUTION_NOT_CONSUMED,
        declared_status = input_addressability(declared),
        candidate_path = NA_character_, candidate_length = NA_integer_,
        source_sha256 = NA_character_, staged_sha256 = NA_character_,
        stringsAsFactors = FALSE)
      next
    }

    res <- do.call(rbind, lapply(split(seq_len(nrow(manifest)), row_dataset), function(idx) {
      ds <- row_dataset[[idx[[1]]]]
      r <- resolve_runtime_paths(
        declared[idx], repository_root = repository_root,
        stage_root = path_stage_root("differential_abundance",
                                     "run_clusterprofiler_enrichment", ds),
        stage = stage)
      r$.idx <- idx
      r
    }))
    res <- res[order(res$.idx), , drop = FALSE]
    manifest[[paste0(field, "_declared")]] <- res$declared_path
    manifest[[field]] <- res$runtime_path
    res$.idx <- NULL
    rows[[length(rows) + 1L]] <- data.frame(
      dataset = row_dataset, field = field, field_kind = kind, res,
      stringsAsFactors = FALSE)
  }
  list(manifest = manifest,
       resolution = if (length(rows)) do.call(rbind, rows) else NULL)
}

resolve_manifest_contract_paths <- function(manifest, path_columns, repository_root = repo_path()) {
  for (column in intersect(path_columns, names(manifest))) {
    manifest[[column]] <- resolve_repository_contract_path(manifest[[column]], repository_root)
  }
  manifest
}

read_canonical_clusterprofiler_manifest <- function(path, dataset, strict = TRUE,
                                                    require_files = TRUE,
                                                    repository_root = repo_path(),
                                                    stage_long_paths = TRUE) {
  if (!file.exists(path)) stop("Canonical clusterProfiler manifest not found: ", path, call. = FALSE)
  manifest <- read_csv_contract(path)
  manifest <- resolve_manifest_contract_paths(
    manifest, names(CLUSTERPROFILER_MANIFEST_PATH_FIELDS), repository_root
  )
  ## The stored manifest on disk is never rewritten. Each runtime-required
  ## field keeps its declared value in a <field>_declared column and carries
  ## the usable runtime path in the original column, so every downstream
  ## reader works unchanged while provenance stays recoverable.
  resolution <- clusterprofiler_manifest_runtime_resolution(
    manifest, dataset = dataset, repository_root = repository_root,
    stage = stage_long_paths)
  manifest <- resolution$manifest
  validate_clusterprofiler_manifest_contract(manifest, strict = strict, require_files = require_files)
  attr(manifest, "runtime_resolution") <- resolution$resolution
  manifest <- manifest[as.character(manifest$dataset) == as.character(dataset), , drop = FALSE]
  if (!nrow(manifest)) stop("Canonical clusterProfiler manifest has no rows for dataset ", dataset, ".", call. = FALSE)
  manifest[order(manifest$dataset, manifest$comparison, manifest$result_type, manifest$ontology, method = "radix"), , drop = FALSE]
}

read_canonical_comparego_manifest <- function(path, dataset, require_files = TRUE,
                                              repository_root = repo_path()) {
  if (!file.exists(path)) stop("Canonical compareGO manifest not found: ", path, call. = FALSE)
  manifest <- read_csv_contract(path)
  manifest <- resolve_manifest_contract_paths(
    manifest,
    c("input_manifest", "term_comparison_file", "term_gene_provenance_output_file", "analysis_status_summary_file"),
    repository_root
  )
  validate_comparego_manifest_contract(manifest, require_files = require_files)
  manifest <- manifest[as.character(manifest$dataset) == as.character(dataset), , drop = FALSE]
  if (!nrow(manifest)) stop("Canonical compareGO manifest has no rows for dataset ", dataset, ".", call. = FALSE)
  manifest[order(manifest$dataset, manifest$comparison, manifest$result_type, manifest$ontology, method = "radix"), , drop = FALSE]
}

empty_collapsed_gene_contract <- function() {
  data.frame(
    dataset = character(), comparison = character(), result_type = character(), ontology = character(),
    official_gene_symbol = character(), official_entrez_id = character(), collapsed_statistic = numeric(),
    collapsed_logfc = numeric(), contributing_ProteinGroupIDs = character(), source_file = character(),
    stringsAsFactors = FALSE
  )
}

join_term_provenance_to_collapsed_genes <- function(provenance, collapsed, strict = TRUE) {
  validate_term_gene_provenance_contract(provenance, strict = strict)
  identity <- c("dataset", "comparison", "result_type", "ontology", "official_gene_symbol")
  required_collapsed <- c(identity, "collapsed_statistic", "collapsed_logfc",
    "contributing_ProteinGroupIDs", "source_file")
  missing <- setdiff(required_collapsed, names(collapsed))
  if (length(missing)) {
    stop("Collapsed canonical gene input is missing required downstream columns: ",
      paste(missing, collapse = ", "), call. = FALSE)
  }
  out <- merge(
    provenance,
    collapsed[required_collapsed],
    by = identity,
    all.x = TRUE,
    sort = FALSE
  )
  if (nrow(out) && any(!is.finite(out$collapsed_statistic) | !is.finite(out$collapsed_logfc))) {
    stop("Canonical term-gene provenance could not be joined to finite collapsed rank/log2fc values.", call. = FALSE)
  }
  out <- out[order(out$dataset, out$comparison, out$result_type, out$ontology,
    out$term_id, out$official_gene_symbol, out$ProteinGroupID, method = "radix"), , drop = FALSE]
  rownames(out) <- NULL
  out
}

read_canonical_clusterprofiler_bundle <- function(manifest_path, dataset,
                                                  result_types = "GSEA_GO",
                                                  ontology = NULL,
                                                  comparisons = NULL,
                                                  strict = TRUE,
                                                  repository_root = repo_path()) {
  manifest <- read_canonical_clusterprofiler_manifest(
    manifest_path, dataset, strict = strict, require_files = TRUE,
    repository_root = repository_root
  )
  manifest <- manifest[manifest$result_type %in% result_types, , drop = FALSE]
  if (!is.null(ontology)) manifest <- manifest[toupper(manifest$ontology) %in% toupper(ontology), , drop = FALSE]
  if (!is.null(comparisons)) manifest <- manifest[manifest$comparison %in% comparisons, , drop = FALSE]
  if (!nrow(manifest)) stop("No supported canonical clusterProfiler manifest rows matched the requested identity.", call. = FALSE)

  collected <- collect_canonical_comparego_outputs(manifest, strict = strict, require_files = TRUE)
  collapsed_rows <- list()
  success <- manifest$analysis_status %in% c("success_with_terms", "success_zero_terms")
  for (i in which(success)) {
    row <- manifest[i, , drop = FALSE]
    collapsed <- read_csv_contract(
      as.character(row$collapsed_gene_input_file),
      character_columns = "official_entrez_id"
    )
    required <- c("official_gene_symbol", "official_entrez_id", "collapsed_statistic",
      "collapsed_logfc", "contributing_ProteinGroupIDs")
    missing <- setdiff(required, names(collapsed))
    if (length(missing)) {
      stop("Collapsed canonical gene input is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    if (!is.character(collapsed$official_entrez_id)) {
      stop("Collapsed canonical gene input official_entrez_id must remain character.", call. = FALSE)
    }
    if (anyDuplicated(collapsed$official_gene_symbol)) {
      stop("Collapsed canonical gene input contains duplicate official_gene_symbol values.", call. = FALSE)
    }
    collapsed$dataset <- as.character(row$dataset)
    collapsed$comparison <- as.character(row$comparison)
    collapsed$result_type <- as.character(row$result_type)
    collapsed$ontology <- as.character(row$ontology)
    collapsed$source_file <- as.character(row$collapsed_gene_input_file)
    collapsed_rows[[length(collapsed_rows) + 1L]] <- collapsed
  }
  collapsed <- if (length(collapsed_rows)) do.call(rbind, collapsed_rows) else empty_collapsed_gene_contract()
  if (nrow(collapsed)) {
    collapsed <- collapsed[order(collapsed$dataset, collapsed$comparison, collapsed$ontology,
      collapsed$official_gene_symbol, method = "radix"), , drop = FALSE]
  }
  rownames(manifest) <- rownames(collapsed) <- NULL
  list(
    manifest = manifest,
    terms = collected$terms,
    provenance = collected$provenance,
    status = collected$status,
    collapsed = collapsed,
    manifest_source = normalizePath(manifest_path, winslash = "/", mustWork = FALSE)
  )
}

read_single_declared_contract_table <- function(paths, label, character_columns = character()) {
  paths <- sort(unique(as.character(paths[!is.na(paths) & nzchar(paths)])), method = "radix")
  if (length(paths) != 1L) stop("Canonical compareGO manifest must declare exactly one ", label, ".", call. = FALSE)
  ## "does not exist" is a claim this check is not entitled to make: the same
  ## FALSE is returned for a path past the character limit and for an unmounted
  ## declared root, and those are present-but-unopenable, not absent.
  declared_status <- input_addressability(paths[[1]])
  if (!identical(declared_status, INPUT_STATUS_PRESENT)) {
    stop("Declared ", label, " is not usable: ",
      describe_input_status_failures(paths[[1]], declared_status), call. = FALSE)
  }
  list(path = paths[[1]], data = read_csv_contract(paths[[1]], character_columns = character_columns))
}

read_canonical_comparego_bundle <- function(manifest_path, dataset, repository_root = repo_path()) {
  manifest <- read_canonical_comparego_manifest(
    manifest_path, dataset, require_files = TRUE, repository_root = repository_root
  )
  terms <- read_single_declared_contract_table(manifest$term_comparison_file, "compareGO term comparison")
  provenance <- read_single_declared_contract_table(
    manifest$term_gene_provenance_output_file, "compareGO term-gene provenance", "official_entrez_id"
  )
  status <- read_single_declared_contract_table(manifest$analysis_status_summary_file, "compareGO analysis-status summary")
  required_terms <- c("dataset", "comparison", "result_type", "ontology", "term_id", "term_description", "NES", "p.adjust")
  missing_terms <- setdiff(required_terms, names(terms$data))
  if (length(missing_terms)) stop("compareGO term comparison is missing required columns: ", paste(missing_terms, collapse = ", "), call. = FALSE)
  validate_term_gene_provenance_contract(provenance$data, strict = TRUE)
  required_status <- c("dataset", "comparison", "result_type", "ontology", "analysis_status", "n_terms", "comparego_action")
  missing_status <- setdiff(required_status, names(status$data))
  if (length(missing_status)) stop("compareGO analysis-status summary is missing required columns: ", paste(missing_status, collapse = ", "), call. = FALSE)
  if (nrow(provenance$data)) {
    term_keys <- paste(terms$data$dataset, terms$data$comparison, terms$data$result_type,
      terms$data$ontology, terms$data$term_id, sep = "\r")
    provenance_keys <- paste(provenance$data$dataset, provenance$data$comparison,
      provenance$data$result_type, provenance$data$ontology, provenance$data$term_id, sep = "\r")
    if (any(!provenance_keys %in% term_keys)) {
      stop("compareGO term-gene provenance contains terms absent from the declared term comparison.", call. = FALSE)
    }
  }
  for (x in list(terms$data, provenance$data, status$data)) {
    if (nrow(x) && any(as.character(x$dataset) != as.character(dataset))) {
      stop("Canonical compareGO output contains dataset identity inconsistent with its manifest.", call. = FALSE)
    }
  }
  terms$data <- terms$data[order(terms$data$dataset, terms$data$comparison, terms$data$ontology,
    terms$data$term_id, method = "radix"), , drop = FALSE]
  provenance$data <- provenance$data[order(provenance$data$dataset, provenance$data$comparison,
    provenance$data$ontology, provenance$data$term_id, provenance$data$official_gene_symbol,
    provenance$data$ProteinGroupID, method = "radix"), , drop = FALSE]
  status$data <- status$data[order(status$data$dataset, status$data$comparison,
    status$data$result_type, status$data$ontology, method = "radix"), , drop = FALSE]
  rownames(terms$data) <- rownames(provenance$data) <- rownames(status$data) <- NULL
  list(
    manifest = manifest, terms = terms$data, provenance = provenance$data, status = status$data,
    manifest_source = normalizePath(manifest_path, winslash = "/", mustWork = FALSE),
    term_source = terms$path, provenance_source = provenance$path, status_source = status$path
  )
}

collect_canonical_comparego_outputs <- function(manifest, strict = TRUE, require_files = TRUE) {
  validate_clusterprofiler_manifest_contract(manifest, strict = strict, require_files = require_files)
  supported <- canonical_comparego_result_types()
  unsupported <- manifest[!manifest$result_type %in% supported, , drop = FALSE]
  selected <- manifest[manifest$result_type %in% supported, , drop = FALSE]
  if (!nrow(selected)) stop("No supported canonical compareGO result types were selected.", call. = FALSE)
  selected <- selected[order(selected$dataset, selected$comparison, selected$ontology, selected$result_type, method = "radix"), , drop = FALSE]

  term_rows <- list()
  provenance_rows <- list()
  status_rows <- list()
  for (i in seq_len(nrow(selected))) {
    row <- selected[i, , drop = FALSE]
    status_rows[[i]] <- data.frame(
      dataset = as.character(row$dataset), comparison = as.character(row$comparison),
      result_type = as.character(row$result_type), ontology = as.character(row$ontology),
      analysis_status = as.character(row$analysis_status), n_terms = as.integer(row$n_terms),
      comparego_action = if (row$analysis_status == "failed") "recorded_failed" else if (row$analysis_status == "success_zero_terms") "completed_zero_terms" else "included",
      stringsAsFactors = FALSE
    )
    if (row$analysis_status == "failed") next
    terms <- read_csv_contract(as.character(row$output_table))
    validate_gsea_result_table_contract(terms, paste0(row$dataset, "/", row$comparison))
    if (nrow(terms) != as.integer(row$n_terms)) {
      stop("Manifest n_terms does not match result table for ", row$dataset, "/", row$comparison, ".", call. = FALSE)
    }
    if (nrow(terms)) {
      terms$term_id <- as.character(terms$ID)
      terms$term_description <- as.character(terms$Description)
      terms$dataset <- as.character(row$dataset)
      terms$comparison <- as.character(row$comparison)
      terms$result_type <- as.character(row$result_type)
      terms$ontology <- as.character(row$ontology)
      term_rows[[length(term_rows) + 1L]] <- terms
    }
    provenance <- read_csv_contract(as.character(row$term_gene_provenance_file), character_columns = "official_entrez_id")
    validate_term_gene_provenance_contract(provenance, strict = strict)
    if (!nrow(terms) && nrow(provenance)) {
      stop("Zero-term GSEA result has non-empty term-gene provenance for ",
        row$dataset, "/", row$comparison, ".", call. = FALSE)
    }
    if (nrow(provenance) && any(!provenance$term_id %in% as.character(terms$ID))) {
      stop("Term-gene provenance contains term IDs absent from the declared GSEA result.", call. = FALSE)
    }
    if (nrow(provenance)) {
      expected <- c(dataset = as.character(row$dataset), comparison = as.character(row$comparison),
        result_type = as.character(row$result_type), ontology = as.character(row$ontology))
      for (column in names(expected)) {
        if (any(as.character(provenance[[column]]) != expected[[column]])) {
          stop("Term-gene provenance metadata does not match manifest field ", column, ".", call. = FALSE)
        }
      }
      provenance_rows[[length(provenance_rows) + 1L]] <- provenance
    }
  }
  empty_terms <- data.frame(
    ID = character(), Description = character(), NES = numeric(), p.adjust = numeric(),
    setSize = integer(), core_enrichment = character(), term_id = character(),
    term_description = character(), dataset = character(), comparison = character(),
    result_type = character(), ontology = character(), stringsAsFactors = FALSE
  )
  empty_provenance <- as.data.frame(setNames(lapply(term_gene_provenance_columns(), function(column) {
    if (column %in% c("gene_level_claim_allowed", "core_enrichment_member")) logical()
    else if (column == "rank_statistic") numeric() else character()
  }), term_gene_provenance_columns()), stringsAsFactors = FALSE)
  terms <- if (length(term_rows)) do.call(rbind, term_rows) else empty_terms
  provenance <- if (length(provenance_rows)) do.call(rbind, provenance_rows) else empty_provenance
  status <- do.call(rbind, status_rows)
  if (nrow(unsupported)) {
    unsupported_status <- data.frame(
      dataset = as.character(unsupported$dataset), comparison = as.character(unsupported$comparison),
      result_type = as.character(unsupported$result_type), ontology = as.character(unsupported$ontology),
      analysis_status = as.character(unsupported$analysis_status), n_terms = as.integer(unsupported$n_terms),
      comparego_action = "skipped_unsupported_result_type", stringsAsFactors = FALSE
    )
    status <- rbind(status, unsupported_status)
  }
  if (nrow(terms)) terms <- terms[order(terms$dataset, terms$comparison, terms$ontology, terms$ID, method = "radix"), , drop = FALSE]
  if (nrow(provenance)) provenance <- provenance[order(provenance$dataset, provenance$comparison, provenance$ontology,
    provenance$term_id, provenance$official_gene_symbol, provenance$ProteinGroupID, method = "radix"), , drop = FALSE]
  status <- status[order(status$dataset, status$comparison, status$result_type, status$ontology, method = "radix"), , drop = FALSE]
  rownames(terms) <- rownames(provenance) <- rownames(status) <- NULL
  list(input_manifest = selected, terms = terms, provenance = provenance, status = status)
}

# Broad heuristic text classes retained for backward-compatible technical
# outputs. These regexes are not authoritative manuscript GO-theme mapping;
# use map_go_terms_to_manuscript_themes() for valid GO biological-process IDs.
biological_program_patterns <- function() {
  data.frame(
    biological_program = c(
      "RNA_RNP_processing",
      "Ribosome_Translation",
      "Mitochondria_OXPHOS_Metabolism",
      "Proteostasis_Ubiquitin_Folding",
      "Synapse_Vesicle_Organization",
      "Cytoskeleton_Motility",
      "Development_Patterning",
      "HPA_Glucocorticoid_Response",
      "Neuroimmune_Complement_Phagosome",
      "ECM_Vascular_Barrier",
      "Oxidative_Redox_Stress",
      "Lipid_Myelin_Membrane",
      "Autophagy_Lysosome"
    ),
    pattern = c(
      "rna|ribonucleoprotein|rnp|splice|splicing|mrna|ncrna|rrna|trna|nucleolus|ribonucle|rna processing",
      "translation|ribosom|peptide biosynthetic|cytoplasmic translation|translational initiation|elongation factor|initiation factor",
      "mitochond|oxidative phosphorylation|oxphos|electron transport|respiratory chain|atp synthesis|oxidoreduct|metabol|glycolys|tricarboxylic|acetyl.coa|energy",
      "proteas|ubiquitin|folding|chaperone|heat shock|proteostasis|protein quality",
      "synap|vesicle|neurotransmitter|axon|dendrit|postsynap|presynap|exocytosis|endocytosis|synaptic organization",
      "cytoskeleton|actin|tubulin|microtubule|motility|adhesion|migration|extracellular matrix",
      "develop|pattern|morphogen|differentiation|neurogenesis|gliogenesis|axon guidance|cell fate|regionalization|dorsal.ventral|anterior.posterior",
      "glucocorticoid|corticosterone|cortisol|steroid hormone|hpa axis|stress hormone|nr3c1|nuclear receptor subfamily 3 group c",
      "microglia|immune|inflamm|cytokine|chemokine|complement|phagocyt|phagosome|lysosomal engulfment|antigen presentation|mhc",
      "extracellular matrix|collagen|laminin|basement membrane|vascular|blood vessel|endothelial|pericyte|blood.brain barrier|barrier|integrin",
      "oxidative stress|redox|reactive oxygen|ros|peroxid|glutathione|superoxide|oxidant|antioxidant",
      "lipid|fatty acid|cholesterol|membrane|myelin|oligodendro|sphingolipid|phospholipid",
      "autophag|lysosom|endosom|phagolysosom|vacuolar|proteolysis"
    ),
    stringsAsFactors = FALSE
  )
}

# Legacy first-hit text classifier for generic/non-manuscript consumers.
map_terms_to_programs <- function(df, description_col = "Description") {
  if (!description_col %in% names(df)) {
    df$biological_program <- NA_character_
    return(df)
  }
  patterns <- biological_program_patterns()
  desc <- tolower(as.character(df[[description_col]]))
  hits <- lapply(seq_len(nrow(patterns)), function(i) grepl(patterns$pattern[[i]], desc, ignore.case = TRUE))
  hit_mat <- do.call(cbind, hits)
  program <- rep(NA_character_, length(desc))
  first_hit <- max.col(hit_mat, ties.method = "first")
  any_hit <- rowSums(hit_mat, na.rm = TRUE) > 0
  program[any_hit] <- patterns$biological_program[first_hit[any_hit]]
  df$biological_program <- program
  df
}

read_csv_if_exists <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  readr::read_csv(path, show_col_types = FALSE)
}

first_existing_path <- function(paths) {
  paths <- paths[!is.na(paths) & nzchar(paths)]
  if (!length(paths)) return(NA_character_)
  paths <- unique(normalizePath(paths, winslash = "/", mustWork = FALSE))
  hit <- paths[file.exists(paths)]
  if (!length(hit)) return(NA_character_)
  hit[[1]]
}

latest_file <- function(root, pattern) {
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  if (!dir.exists(root)) return(NA_character_)
  files <- list.files(root, pattern = pattern, full.names = TRUE, recursive = TRUE)
  files <- files[file.exists(files)]
  if (!length(files)) return(NA_character_)
  info <- file.info(files)
  normalizePath(rownames(info)[order(info$mtime, decreasing = TRUE)[1]], winslash = "/", mustWork = FALSE)
}

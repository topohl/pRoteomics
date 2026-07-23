# Canonical WGCNA identity-contract helpers.
#
# This layer is read-only with respect to WGCNA network state and existing
# Stage 01-13 outputs. It derives current identity only from frozen state plus
# an explicitly selected Stage 01 membership artifact.

wgcna_identity_contract_version <- function() {
  "wgcna_identity_contract_v1"
}

wgcna_identity_datasets <- function() {
  c("neuron_soma", "neuron_neuropil", "microglia")
}

wgcna_identity_validate_dataset <- function(dataset) {
  dataset <- as.character(dataset)
  if (length(dataset) != 1L || is.na(dataset) || !dataset %in% wgcna_identity_datasets()) {
    stop(
      "Unsupported WGCNA identity-contract dataset. Expected one of: ",
      paste(wgcna_identity_datasets(), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  dataset
}

wgcna_identity_clean_id <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | !nzchar(x)] <- NA_character_
  x
}

wgcna_identity_normalize_hex_module_id <- function(x) {
  raw <- wgcna_identity_clean_id(x)
  upper <- toupper(raw)
  ok <- !is.na(upper) & grepl("^(ME|WGCNA_)?#[0-9A-F]{6}$", upper, perl = TRUE)
  out <- rep(NA_character_, length(raw))
  out[ok] <- paste0(
    "WGCNA_#",
    sub("^(ME|WGCNA_)?#", "", upper[ok], perl = TRUE)
  )
  out
}

wgcna_identity_bridge_module_ids <- function(raw_ids, state_module_ids) {
  raw_ids <- wgcna_identity_clean_id(raw_ids)
  state_module_ids <- wgcna_identity_clean_id(state_module_ids)
  normalized <- raw_ids
  bridge_type <- rep("exact_current_identity_match", length(raw_ids))

  not_exact <- is.na(raw_ids) | !raw_ids %in% state_module_ids
  normalized[not_exact] <- wgcna_identity_normalize_hex_module_id(raw_ids[not_exact])
  bridge_type[not_exact & !is.na(normalized)] <- "syntactic_module_id_bridge_only"
  bridge_type[is.na(normalized) | !normalized %in% state_module_ids] <-
    "unresolved_no_authoritative_bridge"

  data.frame(
    original_id = raw_ids,
    normalized_id = normalized,
    bridge_type = bridge_type,
    stringsAsFactors = FALSE
  )
}

wgcna_identity_sha256_text <- function(x) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("The digest package is required for deterministic WGCNA identity hashes.", call. = FALSE)
  }
  digest::digest(enc2utf8(as.character(x)), algo = "sha256", serialize = FALSE)
}

wgcna_identity_membership_serialization <- function(membership) {
  required <- c("dataset", "module_id", "supermodule_id")
  missing <- setdiff(required, names(membership))
  if (length(missing)) {
    stop("Membership serialization is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  x <- membership[, required, drop = FALSE]
  for (nm in required) x[[nm]] <- trimws(as.character(x[[nm]]))
  x <- x[order(x$dataset, x$module_id, x$supermodule_id), , drop = FALSE]
  paste0(
    paste(paste(x$dataset, x$module_id, x$supermodule_id, sep = "|"), collapse = "\n"),
    "\n"
  )
}

wgcna_identity_membership_version <- function(membership) {
  wgcna_identity_sha256_text(wgcna_identity_membership_serialization(membership))
}

wgcna_identity_supermodule_membership_key <- function(dataset, member_module_ids) {
  members <- sort(unique(wgcna_identity_clean_id(member_module_ids)))
  members <- members[!is.na(members)]
  wgcna_identity_sha256_text(
    paste0(
      "dataset=", dataset, "\n",
      "members=", paste(members, collapse = ";"), "\n"
    )
  )
}

wgcna_identity_source_paths <- function(dataset) {
  dataset <- wgcna_identity_validate_dataset(dataset)
  stage01_tables <- path_results("tables", "06_modules_WGCNA", "01_WGCNA", dataset)
  stage01_logs <- path_results("logs", "06_modules_WGCNA", "01_WGCNA", dataset)
  list(
    frozen_state = path_processed(
      "06_modules_WGCNA", "01_WGCNA", dataset, "wgcna_final_model_state.rds"
    ),
    module_definitions = file.path(
      stage01_tables, "modules", "WGCNA_module_definitions_for_downstream.csv"
    ),
    current_membership_mapping = file.path(
      stage01_tables, "supermodules", "wgcna_module_supermodule_annotation.csv"
    ),
    current_supermodule_summary = file.path(
      stage01_tables, "supermodules", "wgcna_supermodule_summary.csv"
    ),
    selected_cluster_assignment = file.path(
      stage01_tables, "supermodules", "wgcna_supermodule_eigengene_clusters.csv"
    ),
    clustering_sensitivity = file.path(
      stage01_tables, "supermodules", "supermodule_clustering_sensitivity.csv"
    ),
    run_manifest = file.path(stage01_logs, "run_manifest.yml"),
    wgcna_run_manifest = file.path(stage01_logs, "wgcna_run_manifest.yml"),
    output_manifest = file.path(stage01_logs, "output_manifest.csv"),
    parameter_audit = file.path(
      stage01_tables, "modules", "WGCNA_parameter_audit.csv"
    ),
    analysis_parameters = file.path(stage01_logs, "analysis_parameters.csv")
  )
}

wgcna_identity_relative_path <- function(path) {
  normalize_slashes <- function(x) gsub("\\\\", "/", x)
  path_abs <- normalize_slashes(normalizePath(path, winslash = "/", mustWork = FALSE))
  root_abs <- normalize_slashes(normalizePath(repo_root(), winslash = "/", mustWork = FALSE))
  prefix <- paste0(sub("/+$", "", root_abs), "/")
  if (startsWith(tolower(path_abs), tolower(prefix))) {
    return(substring(path_abs, nchar(prefix) + 1L))
  }
  path_abs
}

wgcna_identity_read_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  utils::read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA", "NaN", "NULL")
  )
}

wgcna_identity_manifest_parameter <- function(path, parameter) {
  if (!file.exists(path)) return(NA_character_)
  if (requireNamespace("yaml", quietly = TRUE)) {
    parsed <- tryCatch(yaml::read_yaml(path), error = function(e) NULL)
    if (is.list(parsed) && is.list(parsed$parameters) && !is.null(parsed$parameters[[parameter]])) {
      return(as.character(parsed$parameters[[parameter]])[[1]])
    }
  }
  lines <- readLines(path, warn = FALSE)
  pattern <- paste0("^\\s*", gsub("([.])", "\\\\\\1", parameter), "\\s*:\\s*")
  hit <- grep(pattern, lines, value = TRUE)
  if (!length(hit)) return(NA_character_)
  trimws(sub(pattern, "", hit[[1]]))
}

wgcna_identity_source_hashes <- function(dataset, paths, selected_membership_path = NA_character_) {
  roles <- names(paths)
  selected_abs <- normalizePath(selected_membership_path, winslash = "/", mustWork = FALSE)
  rows <- lapply(roles, function(role) {
    path <- paths[[role]]
    exists <- file.exists(path) && !dir.exists(path)
    info <- if (exists) file.info(path) else NULL
    path_abs <- normalizePath(path, winslash = "/", mustWork = FALSE)
    validation_status <- if (!exists) {
      "missing"
    } else if (identical(tolower(path_abs), tolower(selected_abs))) {
      "accepted_membership_authority"
    } else if (role %in% c("frozen_state", "module_definitions")) {
      "accepted_module_identity_authority"
    } else if (dataset == "microglia" && role == "current_supermodule_summary") {
      "supporting_current_identity_provenance"
    } else if (dataset != "microglia" &&
               role %in% c("current_membership_mapping", "current_supermodule_summary")) {
      "compatibility_audit_only_conflicting_system"
    } else {
      "supporting_provenance"
    }
    required <- role %in% c(
      "frozen_state", "module_definitions", "current_membership_mapping",
      "current_supermodule_summary", "selected_cluster_assignment",
      "clustering_sensitivity", "run_manifest", "wgcna_run_manifest",
      "output_manifest", "parameter_audit"
    )
    data.frame(
      dataset = dataset,
      source_role = role,
      source_path = wgcna_identity_relative_path(path),
      exists = exists,
      size_bytes = if (exists) as.numeric(info$size) else NA_real_,
      sha256 = if (exists) file_hash_sha256(path) else NA_character_,
      modified_time_utc = if (exists) {
        format(info$mtime, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
      } else {
        NA_character_
      },
      required = required,
      validation_status = validation_status,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

wgcna_identity_state_module_ids <- function(state) {
  candidates <- list(
    if (is.data.frame(state$module_label_table)) state$module_label_table$ModuleID else NULL,
    if (is.data.frame(state$module_summary)) state$module_summary$ModuleID else NULL,
    if (is.data.frame(state$WGCNA_modules_long)) state$WGCNA_modules_long$ModuleID else NULL
  )
  for (candidate in candidates) {
    ids <- sort(unique(wgcna_identity_clean_id(candidate)))
    ids <- ids[!is.na(ids)]
    if (length(ids)) return(ids)
  }
  stop("Frozen WGCNA state does not contain an explicit current ModuleID set.", call. = FALSE)
}

wgcna_identity_state_feature_ids <- function(state) {
  expr_ids <- if (!is.null(state$expression.data)) colnames(state$expression.data) else NULL
  color_ids <- names(state$mergedColors)
  ids <- if (length(expr_ids)) expr_ids else color_ids
  wgcna_identity_clean_id(ids)
}

wgcna_identity_definition_feature_ids <- function(definitions) {
  key <- if ("ProteinGroupID" %in% names(definitions)) {
    "ProteinGroupID"
  } else if ("ProteinID" %in% names(definitions)) {
    "ProteinID"
  } else {
    NA_character_
  }
  if (is.na(key)) return(list(key = NA_character_, ids = character()))
  list(key = key, ids = wgcna_identity_clean_id(definitions[[key]]))
}

wgcna_identity_candidate_membership <- function(candidate, dataset, state_module_ids) {
  empty_membership <- data.frame(
    dataset = character(),
    module_id = character(),
    supermodule_id = character(),
    module_id_original = character(),
    module_id_normalized = character(),
    module_id_bridge_type = character(),
    stringsAsFactors = FALSE
  )
  if (is.null(candidate) || !is.data.frame(candidate) || !nrow(candidate)) {
    return(list(
      membership = empty_membership,
      module_bridge = data.frame(),
      selected_cut_height = NA_real_,
      error = "selected membership candidate is missing or empty"
    ))
  }

  module_col <- if ("ModuleID" %in% names(candidate) &&
                    any(nzchar(trimws(as.character(candidate$ModuleID))), na.rm = TRUE)) {
    "ModuleID"
  } else if ("module_eigengene" %in% names(candidate)) {
    "module_eigengene"
  } else if ("WGCNAInternalColor" %in% names(candidate)) {
    "WGCNAInternalColor"
  } else {
    NA_character_
  }
  super_col <- if ("SupermoduleID" %in% names(candidate) &&
                   any(nzchar(trimws(as.character(candidate$SupermoduleID))), na.rm = TRUE)) {
    "SupermoduleID"
  } else if ("Supermodule_DataDriven" %in% names(candidate)) {
    "Supermodule_DataDriven"
  } else {
    NA_character_
  }
  if (is.na(module_col) || is.na(super_col)) {
    return(list(
      membership = empty_membership,
      module_bridge = data.frame(),
      selected_cut_height = NA_real_,
      error = "selected membership candidate lacks module or supermodule identity columns"
    ))
  }

  raw_module <- wgcna_identity_clean_id(candidate[[module_col]])
  raw_super <- wgcna_identity_clean_id(candidate[[super_col]])
  keep <- !is.na(raw_module) & !is.na(raw_super)
  raw_module <- raw_module[keep]
  raw_super <- raw_super[keep]
  bridge <- wgcna_identity_bridge_module_ids(raw_module, state_module_ids)
  membership <- data.frame(
    dataset = dataset,
    module_id = bridge$normalized_id,
    supermodule_id = raw_super,
    module_id_original = bridge$original_id,
    module_id_normalized = bridge$normalized_id,
    module_id_bridge_type = bridge$bridge_type,
    stringsAsFactors = FALSE
  )

  cut_cols <- intersect(
    c("SupermoduleCutHeight", "supermodule_merge_cut_height", "selected_cut_height"),
    names(candidate)
  )
  cut_values <- numeric()
  for (nm in cut_cols) {
    values <- suppressWarnings(as.numeric(candidate[[nm]][keep]))
    cut_values <- c(cut_values, values[is.finite(values)])
  }
  cut_values <- unique(cut_values)

  list(
    membership = membership,
    module_bridge = bridge,
    selected_cut_height = if (length(cut_values) == 1L) cut_values[[1]] else NA_real_,
    cut_values = cut_values,
    error = if (length(cut_values) > 1L) {
      paste0("selected membership candidate has multiple cut heights: ",
             paste(cut_values, collapse = ", "))
    } else {
      NA_character_
    }
  )
}

wgcna_identity_output_manifest_match <- function(output_manifest_path, candidate_path) {
  if (length(candidate_path) != 1L || is.na(candidate_path) ||
      !file.exists(candidate_path)) {
    return(list(pass = FALSE, expected = NA_character_, observed = NA_character_))
  }
  manifest <- wgcna_identity_read_csv(output_manifest_path)
  if (is.null(manifest) || !all(c("relative_path", "md5") %in% names(manifest))) {
    return(list(pass = FALSE, expected = NA_character_, observed = file_hash(candidate_path)))
  }
  rel <- wgcna_identity_relative_path(candidate_path)
  idx <- which(tolower(gsub("\\\\", "/", manifest$relative_path)) == tolower(rel))
  expected <- if (length(idx) == 1L) as.character(manifest$md5[[idx]]) else NA_character_
  observed <- file_hash(candidate_path)
  list(
    pass = length(idx) == 1L && !is.na(expected) && identical(tolower(expected), tolower(observed)),
    expected = expected,
    observed = observed
  )
}

wgcna_identity_validation_row <- function(dataset, check_name, status,
                                           observed_value, expected_value,
                                           severity = "error", details = "") {
  data.frame(
    dataset = dataset,
    check_name = check_name,
    status = status,
    observed_value = as.character(observed_value),
    expected_value = as.character(expected_value),
    severity = severity,
    details = as.character(details),
    stringsAsFactors = FALSE
  )
}

wgcna_identity_build_validation <- function(dataset, paths, state, definitions,
                                             selected, membership,
                                             entity_contract = NULL) {
  rows <- list()
  add <- function(name, pass, observed, expected, severity = "error",
                  details = "", warning_only = FALSE) {
    status <- if (isTRUE(pass)) "pass" else if (warning_only) "warning" else "fail"
    rows[[length(rows) + 1L]] <<- wgcna_identity_validation_row(
      dataset, name, status, observed, expected,
      if (warning_only) "warning" else severity, details
    )
  }

  required_source_roles <- c(
    "frozen_state", "module_definitions", "current_membership_mapping",
    "current_supermodule_summary", "selected_cluster_assignment",
    "clustering_sensitivity", "run_manifest", "wgcna_run_manifest",
    "output_manifest", "parameter_audit"
  )
  for (role in required_source_roles) {
    add(
      paste0(role, "_exists"),
      file.exists(paths[[role]]),
      file.exists(paths[[role]]),
      TRUE,
      details = wgcna_identity_relative_path(paths[[role]])
    )
  }

  state_modules <- tryCatch(
    wgcna_identity_state_module_ids(state),
    error = function(e) character()
  )
  definition_modules <- if ("ModuleID" %in% names(definitions)) {
    sort(unique(stats::na.omit(wgcna_identity_clean_id(definitions$ModuleID))))
  } else {
    character()
  }
  state_features <- wgcna_identity_state_feature_ids(state)
  definition_features <- wgcna_identity_definition_feature_ids(definitions)
  color_features <- wgcna_identity_clean_id(names(state$mergedColors))

  add(
    "state_ProteinGroupID_values_are_unique",
    length(state_features) > 0L && !anyNA(state_features) && !anyDuplicated(state_features),
    paste0("n=", length(state_features), "; duplicates=", sum(duplicated(state_features))),
    "non-missing unique active feature identifiers",
    details = paste0(
      "Frozen state feature identifiers are read from expression.data column names; ",
      "the aligned definition key is ", definition_features$key, "."
    )
  )
  add(
    "state_module_assignments_align_with_active_feature_universe",
    length(state_features) > 0L &&
      setequal(state_features, color_features) &&
      setequal(state_features, definition_features$ids),
    paste0(
      "state=", length(unique(state_features)),
      "; mergedColors=", length(unique(color_features)),
      "; definitions=", length(unique(definition_features$ids))
    ),
    "identical feature sets"
  )
  add(
    "module_IDs_in_state_equal_module_IDs_in_definitions",
    identical(state_modules, definition_modules),
    paste(state_modules, collapse = ";"),
    paste(definition_modules, collapse = ";")
  )

  mapping_modules <- sort(unique(stats::na.omit(membership$module_id)))
  add(
    "selected_membership_candidate_is_parseable",
    is.na(selected$error),
    if (is.na(selected$error)) "no parse error" else selected$error,
    "no parse error",
    details = wgcna_identity_relative_path(selected$membership_path)
  )
  add(
    "module_IDs_in_state_equal_module_IDs_in_membership_mapping",
    identical(state_modules, mapping_modules),
    paste(mapping_modules, collapse = ";"),
    paste(state_modules, collapse = ";")
  )
  add(
    "every_module_maps_to_exactly_one_supermodule",
    nrow(membership) == length(state_modules) &&
      !anyDuplicated(membership$module_id) &&
      !anyNA(membership$module_id) &&
      !anyNA(membership$supermodule_id),
    paste0(
      "rows=", nrow(membership),
      "; unique_modules=", length(unique(membership$module_id)),
      "; unique_supermodules=", length(unique(membership$supermodule_id))
    ),
    paste0("rows=", length(state_modules), "; one row per state module")
  )
  add(
    "no_missing_identities",
    nrow(membership) > 0L &&
      all(nzchar(membership$module_id)) &&
      all(nzchar(membership$supermodule_id)),
    sum(
      is.na(membership$module_id) | !nzchar(membership$module_id) |
        is.na(membership$supermodule_id) | !nzchar(membership$supermodule_id)
    ),
    0
  )

  bridge_ok <- nrow(membership) == length(state_modules) &&
    all(membership$module_id_bridge_type %in%
          c("exact_current_identity_match", "syntactic_module_id_bridge_only"))
  add(
    "module_identity_bridge_is_exact_and_complete",
    bridge_ok,
    paste(
      names(table(membership$module_id_bridge_type)),
      as.integer(table(membership$module_id_bridge_type)),
      sep = "=",
      collapse = ";"
    ),
    "exact or permitted hexadecimal syntax bridge only"
  )

  source_match <- wgcna_identity_output_manifest_match(
    paths$output_manifest, selected$membership_path
  )
  add(
    "selected_membership_artifact_matches_active_output_manifest",
    source_match$pass,
    source_match$observed,
    source_match$expected,
    details = wgcna_identity_relative_path(selected$membership_path)
  )

  cut <- selected$selected_cut_height
  run_cut <- suppressWarnings(as.numeric(
    wgcna_identity_manifest_parameter(paths$run_manifest, "supermodule_merge_cut_height")
  ))
  wgcna_run_cut <- suppressWarnings(as.numeric(
    wgcna_identity_manifest_parameter(paths$wgcna_run_manifest, "supermodule_merge_cut_height")
  ))
  sensitivity <- wgcna_identity_read_csv(paths$clustering_sensitivity)
  sensitivity_selected <- if (!is.null(sensitivity) &&
                              "selected_supermodule_cut_height" %in% names(sensitivity)) {
    unique(suppressWarnings(as.numeric(sensitivity$selected_supermodule_cut_height)))
  } else {
    numeric()
  }
  sensitivity_selected <- sensitivity_selected[is.finite(sensitivity_selected)]
  sensitivity_primary <- if (!is.null(sensitivity) &&
                             "primary_cut_height" %in% names(sensitivity)) {
    unique(suppressWarnings(as.numeric(sensitivity$primary_cut_height)))
  } else {
    numeric()
  }
  sensitivity_primary <- sensitivity_primary[is.finite(sensitivity_primary)]
  cut_agreement <- is.finite(cut) &&
    length(sensitivity_selected) == 1L &&
    length(sensitivity_primary) == 1L &&
    isTRUE(all.equal(cut, run_cut)) &&
    isTRUE(all.equal(cut, wgcna_run_cut)) &&
    isTRUE(all.equal(cut, sensitivity_selected[[1]])) &&
    isTRUE(all.equal(cut, sensitivity_primary[[1]]))
  add(
    "selected_cut_height_agrees_with_active_provenance",
    cut_agreement,
    paste0(
      "artifact=", cut, "; run_manifest=", run_cut,
      "; wgcna_run_manifest=", wgcna_run_cut,
      "; sensitivity_selected=", paste(sensitivity_selected, collapse = ","),
      "; sensitivity_primary=", paste(sensitivity_primary, collapse = ",")
    ),
    "one identical finite selected cut height"
  )

  primary_rows <- if (!is.null(sensitivity) &&
                      all(c("cut_height", "primary_cut_height") %in% names(sensitivity))) {
    sum(
      abs(
        suppressWarnings(as.numeric(sensitivity$cut_height)) -
          suppressWarnings(as.numeric(sensitivity$primary_cut_height))
      ) < 1e-12,
      na.rm = TRUE
    )
  } else {
    0L
  }
  add(
    "selected_primary_cut_has_explicit_sensitivity_row",
    primary_rows > 0L,
    primary_rows,
    "at least one row",
    details = if (primary_rows > 0L) {
      "The sensitivity grid contains the selected cut."
    } else {
      "The sensitivity table explicitly records selected/primary cut metadata, but its diagnostic grid omits that cut."
    },
    warning_only = primary_rows == 0L
  )

  state_cut <- suppressWarnings(as.numeric(state$parameters$supermodule_merge_cut_height))
  add(
    "frozen_state_cut_height_matches_selected_membership",
    is.finite(state_cut) && isTRUE(all.equal(state_cut, cut)),
    state_cut,
    cut,
    details = "A mismatch is retained as historical frozen-state provenance and does not override active selected-cluster provenance.",
    warning_only = !(is.finite(state_cut) && isTRUE(all.equal(state_cut, cut)))
  )

  summary <- wgcna_identity_read_csv(paths$current_supermodule_summary)
  summary_ids <- if (!is.null(summary) && "SupermoduleID" %in% names(summary)) {
    sort(unique(stats::na.omit(wgcna_identity_clean_id(summary$SupermoduleID))))
  } else {
    character()
  }
  selected_ids <- sort(unique(membership$supermodule_id))
  summary_match <- identical(summary_ids, selected_ids)
  add(
    "current_supermodule_summary_IDs_equal_selected_membership_IDs",
    summary_match,
    paste(summary_ids, collapse = ";"),
    paste(selected_ids, collapse = ";"),
    details = if (dataset == "microglia") {
      "Microglia current Stage 01 mapping and summary are required to agree."
    } else {
      "The neuronal current annotation/summary system is compatibility-only under the conditional authority policy."
    },
    warning_only = dataset != "microglia" && !summary_match
  )

  version_a <- if (nrow(membership)) wgcna_identity_membership_version(membership) else NA_character_
  version_b <- if (nrow(membership)) wgcna_identity_membership_version(
    membership[rev(seq_len(nrow(membership))), , drop = FALSE]
  ) else NA_character_
  add(
    "deterministic_membership_hash",
    !is.na(version_a) && identical(version_a, version_b),
    version_a,
    version_b
  )

  if (!is.null(entity_contract)) {
    add(
      "no_duplicate_entity_contract_keys",
      !anyDuplicated(entity_contract[c("dataset", "level", "entity_id")]),
      sum(duplicated(entity_contract[c("dataset", "level", "entity_id")])),
      0
    )
    add(
      "contract_output_cardinality_equals_observed_state_mapping_cardinality",
      nrow(entity_contract) == length(state_modules) + length(selected_ids),
      nrow(entity_contract),
      length(state_modules) + length(selected_ids)
    )
  }

  do.call(rbind, rows)
}

wgcna_identity_selected_membership <- function(dataset, paths, state_module_ids) {
  membership_path <- if (dataset == "microglia") {
    paths$current_membership_mapping
  } else {
    paths$selected_cluster_assignment
  }
  candidate <- wgcna_identity_read_csv(membership_path)
  parsed <- wgcna_identity_candidate_membership(candidate, dataset, state_module_ids)
  parsed$membership_path <- membership_path
  parsed$identity_source <- if (dataset == "microglia") {
    "frozen_state_plus_current_stage01_membership_mapping"
  } else {
    "frozen_state_plus_selected_stage01_cluster_assignment"
  }
  parsed$selected_cut_height_source <- paste(
    wgcna_identity_relative_path(membership_path),
    "confirmed by run_manifest.yml, wgcna_run_manifest.yml, output_manifest.csv,",
    "and supermodule_clustering_sensitivity.csv"
  )
  parsed
}

wgcna_identity_build_membership_contract <- function(dataset, membership, source_hashes,
                                                       selected) {
  if (!nrow(membership)) {
    return(data.frame(
      dataset = character(),
      module_id = character(),
      supermodule_id = character(),
      module_membership_key = character(),
      supermodule_membership_key = character(),
      membership_version = character(),
      state_sha256 = character(),
      mapping_sha256 = character(),
      membership_status = character(),
      contract_version = character(),
      stringsAsFactors = FALSE
    ))
  }
  membership <- membership[order(membership$dataset, membership$module_id), , drop = FALSE]
  membership_version <- wgcna_identity_membership_version(membership)
  state_sha <- source_hashes$sha256[source_hashes$source_role == "frozen_state"][[1]]
  mapping_sha <- file_hash_sha256(selected$membership_path)
  super_keys <- vapply(
    membership$supermodule_id,
    function(sm) {
      wgcna_identity_supermodule_membership_key(
        dataset, membership$module_id[membership$supermodule_id == sm]
      )
    },
    character(1)
  )
  data.frame(
    dataset = dataset,
    module_id = membership$module_id,
    supermodule_id = membership$supermodule_id,
    module_membership_key = paste(dataset, membership$module_id, sep = "|"),
    supermodule_membership_key = unname(super_keys),
    membership_version = membership_version,
    state_sha256 = state_sha,
    mapping_sha256 = mapping_sha,
    membership_status = "current_authoritative_membership",
    contract_version = wgcna_identity_contract_version(),
    stringsAsFactors = FALSE
  )
}

wgcna_identity_build_entity_contract <- function(dataset, state_module_ids, membership_contract,
                                                   source_hashes, selected) {
  state_sha <- source_hashes$sha256[source_hashes$source_role == "frozen_state"][[1]]
  definitions_sha <- source_hashes$sha256[
    source_hashes$source_role == "module_definitions"
  ][[1]]
  mapping_sha <- if (file.exists(selected$membership_path)) {
    file_hash_sha256(selected$membership_path)
  } else {
    NA_character_
  }
  membership_version <- unique(membership_contract$membership_version)
  if (!length(membership_version)) membership_version <- NA_character_
  super_ids <- sort(unique(membership_contract$supermodule_id))
  super_counts <- table(membership_contract$supermodule_id)

  module_rows <- data.frame(
    dataset = dataset,
    level = "module",
    entity_id = state_module_ids,
    canonical_entity_id = state_module_ids,
    identity_role = "canonical_module",
    n_member_modules = 1L,
    stringsAsFactors = FALSE
  )
  super_rows <- data.frame(
    dataset = dataset,
    level = "supermodule",
    entity_id = super_ids,
    canonical_entity_id = super_ids,
    identity_role = ifelse(
      as.integer(super_counts[super_ids]) == 1L,
      "singleton_supermodule_identity",
      "higher_order_supermodule_identity"
    ),
    n_member_modules = as.integer(super_counts[super_ids]),
    stringsAsFactors = FALSE
  )
  out <- rbind(module_rows, super_rows)
  out$membership_version <- membership_version
  out$state_sha256 <- state_sha
  out$module_definitions_sha256 <- definitions_sha
  out$supermodule_mapping_sha256 <- mapping_sha
  out$identity_source <- selected$identity_source
  out$identity_status <- "current_authoritative_identity"
  out$selected_cut_height <- selected$selected_cut_height
  out$selected_cut_height_source <- selected$selected_cut_height_source
  out$contract_version <- wgcna_identity_contract_version()
  out
}

wgcna_identity_parse_members <- function(value, state_module_ids) {
  value <- wgcna_identity_clean_id(value)
  if (length(value) != 1L || is.na(value)) {
    return(list(original = NA_character_, normalized = NA_character_, complete = FALSE))
  }
  members <- trimws(unlist(strsplit(value, "[;,|]", perl = TRUE), use.names = FALSE))
  members <- members[nzchar(members)]
  bridge <- wgcna_identity_bridge_module_ids(members, state_module_ids)
  complete <- length(members) > 0L &&
    all(bridge$bridge_type %in%
          c("exact_current_identity_match", "syntactic_module_id_bridge_only")) &&
    !anyDuplicated(bridge$normalized_id)
  list(
    original = paste(sort(unique(members)), collapse = ";"),
    normalized = if (complete) {
      paste(sort(unique(bridge$normalized_id)), collapse = ";")
    } else {
      NA_character_
    },
    complete = complete
  )
}

wgcna_identity_classify_compatibility <- function(level, downstream_entity_id,
                                                    downstream_normalized_id,
                                                    downstream_membership_key,
                                                    contract_entity_id,
                                                    contract_membership_key,
                                                    label_conflict = FALSE,
                                                    source_system_status = "") {
  level <- as.character(level)
  exact_id <- !is.na(downstream_entity_id) && !is.na(contract_entity_id) &&
    identical(as.character(downstream_entity_id), as.character(contract_entity_id))
  normalized_id <- !is.na(downstream_normalized_id) && !is.na(contract_entity_id) &&
    identical(as.character(downstream_normalized_id), as.character(contract_entity_id))
  exact_membership <- !is.na(downstream_membership_key) &&
    !is.na(contract_membership_key) &&
    identical(as.character(downstream_membership_key), as.character(contract_membership_key))

  if (is.na(downstream_entity_id) || !nzchar(downstream_entity_id)) {
    return("missing_current_identity")
  }
  if (level == "supermodule" && exact_id && !is.na(downstream_membership_key) &&
      !exact_membership) {
    return("conflicting_membership_same_entity_id")
  }
  if (identical(source_system_status, "historical_membership_system") &&
      !exact_membership) {
    return("historical_membership_system")
  }
  if (level == "supermodule" && exact_membership && isTRUE(label_conflict)) {
    return("label_conflict_same_membership")
  }
  if (level == "supermodule" && exact_membership) {
    return("exact_current_identity_match")
  }
  if (level == "module" && normalized_id && !exact_id) {
    return("syntactic_module_id_bridge_only")
  }
  if (level == "module" && exact_id) {
    return("exact_current_identity_match")
  }
  if (!is.na(contract_entity_id) && nzchar(contract_entity_id)) {
    return("unresolved_no_authoritative_bridge")
  }
  "extra_downstream_identity"
}

wgcna_identity_label_value <- function(df) {
  candidates <- intersect(
    c(
      "canonical_biological_label", "final_plot_label", "final_label",
      "ModuleLabel_Final", "Supermodule_FinalLabel", "Supermodule_DisplayLabel",
      "module_biological_label", "supermodule_biological_label", "biological_label"
    ),
    names(df)
  )
  if (!length(candidates)) return(rep(NA_character_, nrow(df)))
  out <- rep(NA_character_, nrow(df))
  for (nm in candidates) {
    value <- wgcna_identity_clean_id(df[[nm]])
    fill <- is.na(out) & !is.na(value)
    out[fill] <- value[fill]
  }
  out
}

wgcna_identity_downstream_specs <- function(dataset) {
  base <- function(...) path_results("tables", ...)
  source_base <- function(...) path_results("source_data", ...)
  rows <- list(
    data.frame(
      stage = "Stage01_selected_cluster_assignment",
      path = base("06_modules_WGCNA", "01_WGCNA", dataset, "supermodules",
                  "wgcna_supermodule_eigengene_clusters.csv"),
      level = "both",
      source_system_status = "current_membership_system",
      stringsAsFactors = FALSE
    ),
    data.frame(
      stage = "Stage01_current_mapping",
      path = base("06_modules_WGCNA", "01_WGCNA", dataset, "supermodules",
                  "wgcna_module_supermodule_annotation.csv"),
      level = "both",
      source_system_status = if (dataset == "microglia") {
        "current_membership_system"
      } else {
        "historical_membership_system"
      },
      stringsAsFactors = FALSE
    ),
    data.frame(
      stage = "Stage01_current_summary",
      path = base("06_modules_WGCNA", "01_WGCNA", dataset, "supermodules",
                  "wgcna_supermodule_summary.csv"),
      level = "supermodule",
      source_system_status = if (dataset == "microglia") {
        "current_membership_system"
      } else {
        "historical_membership_system"
      },
      stringsAsFactors = FALSE
    ),
    data.frame(
      stage = "Stage05_group_effects",
      path = base("06_modules_WGCNA", "group_effects", dataset,
                  "module_group_effects.csv"),
      level = "module",
      source_system_status = "downstream_identity_usage",
      stringsAsFactors = FALSE
    ),
    data.frame(
      stage = "Stage05_group_effects",
      path = base("06_modules_WGCNA", "group_effects", dataset,
                  "supermodule_group_effects.csv"),
      level = "supermodule",
      source_system_status = "downstream_identity_usage",
      stringsAsFactors = FALSE
    ),
    data.frame(
      stage = "Stage05_group_effects",
      path = base("06_modules_WGCNA", "group_effects", dataset,
                  "supermodule_composition.csv"),
      level = "supermodule",
      source_system_status = "downstream_identity_usage",
      stringsAsFactors = FALSE
    ),
    data.frame(
      stage = "Stage06_annotation",
      path = base("06_modules_WGCNA", "module_annotation", dataset,
                  "WGCNA_module_biological_annotation.csv"),
      level = "module",
      source_system_status = "downstream_identity_usage",
      stringsAsFactors = FALSE
    ),
    data.frame(
      stage = "Stage06_annotation",
      path = base("06_modules_WGCNA", "module_annotation", dataset,
                  "WGCNA_supermodule_biological_annotation.csv"),
      level = "supermodule",
      source_system_status = "downstream_identity_usage",
      stringsAsFactors = FALSE
    ),
    data.frame(
      stage = "Stage07_interpretable_summary",
      path = base("06_modules_WGCNA", "interpretable_summary", dataset,
                  "WGCNA_final_label_lookup.csv"),
      level = "from_column",
      source_system_status = "downstream_identity_usage",
      stringsAsFactors = FALSE
    ),
    data.frame(
      stage = "Stage07_interpretable_summary",
      path = base("06_modules_WGCNA", "interpretable_summary", dataset,
                  "WGCNA_module_group_effects_interpretable.csv"),
      level = "module",
      source_system_status = "downstream_identity_usage",
      stringsAsFactors = FALSE
    ),
    data.frame(
      stage = "Stage07_interpretable_summary",
      path = base("06_modules_WGCNA", "interpretable_summary", dataset,
                  "WGCNA_supermodule_group_effects_interpretable.csv"),
      level = "supermodule",
      source_system_status = "downstream_identity_usage",
      stringsAsFactors = FALSE
    ),
    data.frame(
      stage = "Stage08_score_publication",
      path = base("06_modules_WGCNA", "score_publication_summary", dataset,
                  "WGCNA_score_publication_validation.csv"),
      level = "supermodule",
      source_system_status = "downstream_identity_usage",
      stringsAsFactors = FALSE
    ),
    data.frame(
      stage = "circular_atlas",
      path = source_base("10_biological_integration", "wgcna_circular_atlas",
                         "global", "wgcna_circular_atlas_segments.csv"),
      level = "from_column",
      source_system_status = "downstream_identity_usage",
      stringsAsFactors = FALSE
    ),
    data.frame(
      stage = "Stage13_claim_readiness",
      path = base("06_modules_WGCNA", "claim_readiness", dataset,
                  "WGCNA_entity_claim_readiness.csv"),
      level = "from_column",
      source_system_status = "downstream_identity_usage",
      stringsAsFactors = FALSE
    )
  )
  do.call(rbind, rows)
}

wgcna_identity_first_column <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (!length(hit)) return(rep(NA_character_, nrow(df)))
  out <- rep(NA_character_, nrow(df))
  for (nm in hit) {
    value <- wgcna_identity_clean_id(df[[nm]])
    fill <- is.na(out) & !is.na(value)
    out[fill] <- value[fill]
  }
  out
}

wgcna_identity_build_compatibility_audit <- function(dataset, entity_contract,
                                                       membership_contract,
                                                       state_module_ids) {
  specs <- wgcna_identity_downstream_specs(dataset)
  audit_rows <- list()
  contract_modules <- entity_contract[entity_contract$level == "module", , drop = FALSE]
  contract_supers <- entity_contract[entity_contract$level == "supermodule", , drop = FALSE]
  contract_super_keys <- unique(
    membership_contract[c("supermodule_id", "supermodule_membership_key")]
  )

  for (i in seq_len(nrow(specs))) {
    spec <- specs[i, , drop = FALSE]
    if (!file.exists(spec$path)) next
    df <- tryCatch(wgcna_identity_read_csv(spec$path), error = function(e) NULL)
    if (is.null(df) || !nrow(df)) next
    if ("dataset" %in% names(df)) {
      keep_dataset <- is.na(df$dataset) | as.character(df$dataset) == dataset
      df <- df[keep_dataset, , drop = FALSE]
      if (!nrow(df)) next
    }

    levels <- if (spec$level == "both") {
      c("module", "supermodule")
    } else if (spec$level == "from_column" && "level" %in% names(df)) {
      unique(tolower(as.character(df$level)))
    } else {
      spec$level
    }

    for (level in levels) {
      if (!level %in% c("module", "supermodule")) next
      part <- df
      if ("level" %in% names(part) && spec$level == "from_column") {
        part <- part[tolower(as.character(part$level)) == level, , drop = FALSE]
      }
      if (!nrow(part)) next

      raw_id <- if (level == "module") {
        wgcna_identity_first_column(
          part,
          c(
            "entity_id", "ModuleID", "module_id", "module",
            "module_eigengene"
          )
        )
      } else {
        wgcna_identity_first_column(
          part,
          c(
            "entity_id", "SupermoduleID", "supermodule_id",
            "Supermodule_DataDriven", "supermodule"
          )
        )
      }
      normalized_id <- if (level == "module") {
        bridge <- wgcna_identity_bridge_module_ids(raw_id, state_module_ids)
        bridge$normalized_id
      } else {
        raw_id
      }

      member_text <- wgcna_identity_first_column(
        part,
        c(
          "member_ModuleIDs", "member_module_ids", "member_modules",
          "member_modules.x", "member_modules.y",
          "member_module_internal_colors", "Supermodule_CompositionMembers"
        )
      )
      if (level == "supermodule") {
        row_module <- wgcna_identity_first_column(
          part,
          c(
            "ModuleID", "module_id", "module_eigengene",
            "WGCNAInternalColor"
          )
        )
        for (sm in unique(stats::na.omit(raw_id))) {
          idx <- which(raw_id == sm)
          grouped_modules <- sort(unique(stats::na.omit(row_module[idx])))
          if (length(grouped_modules)) {
            member_text[idx] <- paste(grouped_modules, collapse = ";")
          }
        }
      }
      parsed_members <- lapply(member_text, wgcna_identity_parse_members,
                               state_module_ids = state_module_ids)
      member_original <- vapply(parsed_members, `[[`, character(1), "original")
      member_normalized <- vapply(parsed_members, `[[`, character(1), "normalized")
      member_complete <- vapply(parsed_members, `[[`, logical(1), "complete")
      downstream_key <- rep(NA_character_, nrow(part))
      if (level == "supermodule") {
        for (j in which(member_complete)) {
          downstream_key[[j]] <- wgcna_identity_supermodule_membership_key(
            dataset,
            strsplit(member_normalized[[j]], ";", fixed = TRUE)[[1]]
          )
        }
      }

      contract_id <- rep(NA_character_, nrow(part))
      contract_key <- rep(NA_character_, nrow(part))
      if (level == "module") {
        matched <- match(normalized_id, contract_modules$entity_id)
        contract_id <- contract_modules$entity_id[matched]
      } else {
        by_members <- match(downstream_key, contract_super_keys$supermodule_membership_key)
        contract_id <- contract_super_keys$supermodule_id[by_members]
        contract_key <- contract_super_keys$supermodule_membership_key[by_members]
        by_id <- match(raw_id, contract_supers$entity_id)
        missing_contract <- is.na(contract_id)
        contract_id[missing_contract] <- contract_supers$entity_id[by_id[missing_contract]]
        contract_key[missing_contract] <- contract_super_keys$supermodule_membership_key[
          match(contract_id[missing_contract], contract_super_keys$supermodule_id)
        ]
      }

      labels <- wgcna_identity_label_value(part)
      duplicate_key <- duplicated(data.frame(
        id = raw_id,
        stringsAsFactors = FALSE
      )) | duplicated(data.frame(id = raw_id, stringsAsFactors = FALSE), fromLast = TRUE)

      block <- data.frame(
        dataset = dataset,
        downstream_stage = spec$stage,
        source_file = wgcna_identity_relative_path(spec$path),
        level = level,
        downstream_entity_id = raw_id,
        downstream_entity_id_original = raw_id,
        downstream_entity_id_normalized = normalized_id,
        contract_entity_id = contract_id,
        downstream_membership_key = downstream_key,
        contract_membership_key = contract_key,
        downstream_member_modules_original = member_original,
        downstream_member_modules_normalized = member_normalized,
        exact_id_match = !is.na(raw_id) & !is.na(contract_id) & raw_id == contract_id,
        exact_membership_match = !is.na(downstream_key) &
          !is.na(contract_key) & downstream_key == contract_key,
        identity_compatible = FALSE,
        label_conflict_detected = FALSE,
        duplicate_downstream_key = duplicate_key,
        compatibility_status = NA_character_,
        exclusion_or_warning_reason = NA_character_,
        source_system_status = spec$source_system_status,
        downstream_label = labels,
        stringsAsFactors = FALSE
      )
      audit_rows[[length(audit_rows) + 1L]] <- block

      expected_contract_ids <- if (level == "module") {
        contract_modules$entity_id
      } else {
        contract_supers$entity_id
      }
      missing_ids <- setdiff(expected_contract_ids, unique(stats::na.omit(contract_id)))
      if (length(missing_ids)) {
        missing_keys <- if (level == "supermodule") {
          contract_super_keys$supermodule_membership_key[
            match(missing_ids, contract_super_keys$supermodule_id)
          ]
        } else {
          rep(NA_character_, length(missing_ids))
        }
        audit_rows[[length(audit_rows) + 1L]] <- data.frame(
          dataset = dataset,
          downstream_stage = spec$stage,
          source_file = wgcna_identity_relative_path(spec$path),
          level = level,
          downstream_entity_id = NA_character_,
          downstream_entity_id_original = NA_character_,
          downstream_entity_id_normalized = NA_character_,
          contract_entity_id = missing_ids,
          downstream_membership_key = NA_character_,
          contract_membership_key = missing_keys,
          downstream_member_modules_original = NA_character_,
          downstream_member_modules_normalized = NA_character_,
          exact_id_match = FALSE,
          exact_membership_match = FALSE,
          identity_compatible = FALSE,
          label_conflict_detected = FALSE,
          duplicate_downstream_key = FALSE,
          compatibility_status = "missing_current_identity",
          exclusion_or_warning_reason =
            "Current contract identity is absent from this downstream source.",
          source_system_status = spec$source_system_status,
          downstream_label = NA_character_,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (!length(audit_rows)) {
    return(wgcna_identity_empty_compatibility_audit())
  }

  out <- do.call(rbind, audit_rows)
  label_identity <- ifelse(
    out$level == "module",
    out$contract_entity_id,
    ifelse(out$exact_membership_match, out$contract_membership_key, NA_character_)
  )
  label_group <- ifelse(
    is.na(label_identity),
    NA_character_,
    paste(out$dataset, out$level, label_identity, sep = "|")
  )
  for (group in unique(label_group[!is.na(out$contract_entity_id)])) {
    idx <- which(label_group == group)
    values <- unique(stats::na.omit(wgcna_identity_clean_id(out$downstream_label[idx])))
    if (length(values) > 1L) out$label_conflict_detected[idx] <- TRUE
  }

  out$compatibility_status <- vapply(seq_len(nrow(out)), function(i) {
    wgcna_identity_classify_compatibility(
      out$level[[i]],
      out$downstream_entity_id[[i]],
      out$downstream_entity_id_normalized[[i]],
      out$downstream_membership_key[[i]],
      out$contract_entity_id[[i]],
      out$contract_membership_key[[i]],
      out$label_conflict_detected[[i]],
      out$source_system_status[[i]]
    )
  }, character(1))
  out$identity_compatible <- out$compatibility_status %in% c(
    "exact_current_identity_match",
    "syntactic_module_id_bridge_only",
    "label_conflict_same_membership"
  )
  out$exclusion_or_warning_reason <- ifelse(
    out$identity_compatible,
    ifelse(
      out$compatibility_status == "label_conflict_same_membership",
      "Biological labels differ for an identical membership; no label authority was selected.",
      NA_character_
    ),
    vapply(out$compatibility_status, function(status) {
      switch(
        status,
        historical_membership_system =
          "Source uses a coherent but non-authoritative historical membership system.",
        conflicting_membership_same_entity_id =
          "Equal supermodule ID has a different exact member-module composition.",
        missing_current_identity =
          "Current contract identity is absent from this downstream source.",
        extra_downstream_identity =
          "Downstream identity is not present in the current contract.",
        unresolved_no_authoritative_bridge =
          "Identity cannot be verified using exact IDs, permitted hexadecimal syntax normalization, or exact membership.",
        "Identity is not compatible with the current contract."
      )
    }, character(1))
  )

  required_cols <- c(
    "dataset", "downstream_stage", "source_file", "level",
    "downstream_entity_id", "contract_entity_id",
    "downstream_membership_key", "contract_membership_key",
    "exact_id_match", "exact_membership_match", "identity_compatible",
    "label_conflict_detected", "duplicate_downstream_key",
    "compatibility_status", "exclusion_or_warning_reason"
  )
  extra_cols <- setdiff(names(out), required_cols)
  out[c(required_cols, extra_cols)]
}

wgcna_identity_contract_publishable <- function(validation) {
  !any(validation$status == "fail" & validation$severity == "error")
}

wgcna_identity_empty_compatibility_audit <- function() {
  data.frame(
    dataset = character(), downstream_stage = character(),
    source_file = character(), level = character(),
    downstream_entity_id = character(),
    downstream_entity_id_original = character(),
    downstream_entity_id_normalized = character(),
    contract_entity_id = character(),
    downstream_membership_key = character(),
    contract_membership_key = character(),
    downstream_member_modules_original = character(),
    downstream_member_modules_normalized = character(),
    exact_id_match = logical(), exact_membership_match = logical(),
    identity_compatible = logical(), label_conflict_detected = logical(),
    duplicate_downstream_key = logical(), compatibility_status = character(),
    exclusion_or_warning_reason = character(),
    source_system_status = character(), downstream_label = character(),
    stringsAsFactors = FALSE
  )
}

wgcna_identity_build_contract_bundle <- function(dataset) {
  dataset <- wgcna_identity_validate_dataset(dataset)
  paths <- wgcna_identity_source_paths(dataset)
  selected_path <- if (dataset == "microglia") {
    paths$current_membership_mapping
  } else {
    paths$selected_cluster_assignment
  }
  source_hashes <- wgcna_identity_source_hashes(dataset, paths, selected_path)

  state <- if (file.exists(paths$frozen_state)) {
    tryCatch(readRDS(paths$frozen_state), error = identity)
  } else {
    simpleError("Frozen WGCNA state does not exist.")
  }
  if (inherits(state, "error")) {
    required_roles <- c(
      "frozen_state", "module_definitions", "current_membership_mapping",
      "current_supermodule_summary", "selected_cluster_assignment",
      "clustering_sensitivity", "run_manifest", "wgcna_run_manifest",
      "output_manifest", "parameter_audit"
    )
    validation <- do.call(rbind, lapply(required_roles, function(role) {
      exists <- file.exists(paths[[role]])
      wgcna_identity_validation_row(
        dataset,
        paste0(role, "_exists"),
        if (exists) "pass" else "fail",
        exists,
        TRUE,
        "error",
        wgcna_identity_relative_path(paths[[role]])
      )
    }))
    validation <- rbind(
      validation,
      wgcna_identity_validation_row(
        dataset,
        "frozen_state_is_readable",
        "fail",
        conditionMessage(state),
        "readable RDS",
        "error",
        "Identity cannot be derived without a readable frozen WGCNA state."
      )
    )
    empty_membership <- wgcna_identity_build_membership_contract(
      dataset,
      data.frame(),
      source_hashes,
      list(membership_path = selected_path)
    )
    empty_entities <- data.frame(
      dataset = character(), level = character(), entity_id = character(),
      canonical_entity_id = character(), identity_role = character(),
      n_member_modules = integer(), membership_version = character(),
      state_sha256 = character(), module_definitions_sha256 = character(),
      supermodule_mapping_sha256 = character(), identity_source = character(),
      identity_status = character(), selected_cut_height = numeric(),
      selected_cut_height_source = character(), contract_version = character(),
      stringsAsFactors = FALSE
    )
    return(list(
      dataset = dataset,
      paths = paths,
      selected = list(
        membership_path = selected_path,
        selected_cut_height = NA_real_,
        identity_source = NA_character_,
        selected_cut_height_source = NA_character_,
        error = conditionMessage(state)
      ),
      state_module_ids = character(),
      source_hashes = source_hashes,
      validation = validation,
      membership_contract = empty_membership,
      entity_contract = empty_entities,
      compatibility_audit = wgcna_identity_empty_compatibility_audit(),
      publishable = FALSE
    ))
  }

  definitions <- wgcna_identity_read_csv(paths$module_definitions)
  if (is.null(definitions)) definitions <- data.frame()
  state_module_ids <- wgcna_identity_state_module_ids(state)
  selected <- wgcna_identity_selected_membership(dataset, paths, state_module_ids)
  membership_contract <- wgcna_identity_build_membership_contract(
    dataset, selected$membership, source_hashes, selected
  )
  entity_contract <- wgcna_identity_build_entity_contract(
    dataset, state_module_ids, membership_contract, source_hashes, selected
  )
  validation <- wgcna_identity_build_validation(
    dataset, paths, state, definitions, selected, selected$membership,
    entity_contract
  )

  duplicate_membership_keys <- sum(duplicated(
    membership_contract[c("dataset", "module_membership_key")]
  ))
  validation <- rbind(
    validation,
    wgcna_identity_validation_row(
      dataset,
      "no_duplicate_membership_keys",
      if (duplicate_membership_keys == 0L) "pass" else "fail",
      duplicate_membership_keys,
      0,
      "error",
      "module_membership_key is the unique long-table row key; supermodule_membership_key intentionally repeats for members of the same supermodule."
    )
  )
  compatibility <- wgcna_identity_build_compatibility_audit(
    dataset, entity_contract, membership_contract, state_module_ids
  )

  list(
    dataset = dataset,
    paths = paths,
    selected = selected,
    state_module_ids = state_module_ids,
    source_hashes = source_hashes,
    validation = validation,
    membership_contract = membership_contract,
    entity_contract = entity_contract,
    compatibility_audit = compatibility,
    publishable = wgcna_identity_contract_publishable(validation)
  )
}

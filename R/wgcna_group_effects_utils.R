# Canonical Phase 2 WGCNA group-effect modelling utilities.
#
# This layer is deliberately independent of biological labels and historical
# supermodule maps. The Phase 1 identity contract is the sole supermodule
# membership authority.

if (!exists("wgcna_group_classify_statistical_support", mode = "function")) {
  support_utils <- c(
    if (exists("repo_path", mode = "function")) {
      repo_path("R", "wgcna_support_status_utils.R")
    } else character(),
    file.path("R", "wgcna_support_status_utils.R"),
    file.path("..", "R", "wgcna_support_status_utils.R"),
    file.path("..", "..", "R", "wgcna_support_status_utils.R")
  )
  support_utils <- support_utils[file.exists(support_utils)][[1]]
  source(support_utils)
  rm(support_utils)
}

wgcna_group_effects_contract_version <- function() {
  "wgcna_group_effects_phase2b_v5"
}

wgcna_group_boundary_variance_ratio_tolerance <- function() 1e-4

wgcna_group_lme4_singularity_tolerance <- function() 1e-4

wgcna_group_required_packages <- function() {
  c("lme4", "lmerTest", "emmeans")
}

wgcna_group_require_primary_packages <- function() {
  missing <- wgcna_group_required_packages()[
    !vapply(wgcna_group_required_packages(), requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing)) {
    stop(
      "Canonical WGCNA group-effect modelling requires: ",
      paste(missing, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

wgcna_group_contract_paths <- function(dataset) {
  root <- path_results(
    "tables", "06_modules_WGCNA", "identity_contract", dataset
  )
  list(
    entity = file.path(root, "WGCNA_entity_identity_contract.csv"),
    membership = file.path(
      root, "WGCNA_module_supermodule_membership_contract.csv"
    ),
    status = file.path(root, "WGCNA_identity_contract_status.csv")
  )
}

wgcna_group_read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing required input: ", path, call. = FALSE)
  utils::read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA", "NaN", "NULL")
  )
}

wgcna_group_contract_bundle_hash <- function(file_hashes) {
  serialized <- paste0(
    "entity=", file_hashes[["entity"]], "\n",
    "membership=", file_hashes[["membership"]], "\n",
    "status=", file_hashes[["status"]], "\n"
  )
  wgcna_identity_sha256_text(serialized)
}

wgcna_group_assert_membership_version <- function(values) {
  values <- unique(as.character(values))
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values) != 1L) {
    stop("Membership version is missing or inconsistent.", call. = FALSE)
  }
  values[[1]]
}

wgcna_group_load_identity_contract <- function(dataset, state_path) {
  dataset <- wgcna_identity_validate_dataset(dataset)
  paths <- wgcna_group_contract_paths(dataset)
  missing <- names(paths)[!vapply(paths, file.exists, logical(1))]
  if (length(missing)) {
    stop(
      "Phase 1 identity contract is incomplete for ", dataset, ": ",
      paste(missing, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  entity <- wgcna_group_read_csv(paths$entity)
  membership <- wgcna_group_read_csv(paths$membership)
  status <- wgcna_group_read_csv(paths$status)
  hashes <- vapply(paths, file_hash_sha256, character(1))
  state_hash <- file_hash_sha256(state_path)

  fail <- function(message) {
    stop("Identity-contract precondition failed for ", dataset, ": ", message,
         call. = FALSE)
  }
  if (nrow(status) != 1L || !identical(as.character(status$dataset), dataset)) {
    fail("status must contain exactly one row for the requested dataset")
  }
  if (!as.character(status$status[[1]]) %in%
      c("publishable", "identity_contract_publishable")) {
    fail(paste0("status is not publishable: ", status$status[[1]]))
  }
  if (!nrow(entity) || !nrow(membership)) fail("entity or membership table is empty")
  if (any(as.character(entity$dataset) != dataset) ||
      any(as.character(membership$dataset) != dataset)) {
    fail("dataset values are not exact")
  }
  if (anyDuplicated(entity[c("dataset", "level", "entity_id")])) {
    fail("entity keys are duplicated")
  }
  if (anyDuplicated(membership[c("dataset", "module_id")])) {
    fail("module membership keys are duplicated")
  }
  if (any(is.na(membership$module_id) | !nzchar(trimws(membership$module_id))) ||
      any(is.na(membership$supermodule_id) |
          !nzchar(trimws(membership$supermodule_id)))) {
    fail("membership contains missing identities")
  }

  versions <- unique(c(
    as.character(entity$contract_version),
    as.character(membership$contract_version),
    as.character(status$contract_version)
  ))
  versions <- versions[!is.na(versions) & nzchar(versions)]
  if (length(versions) != 1L ||
      !identical(versions[[1]], wgcna_identity_contract_version())) {
    fail("contract version is missing, mixed, or unsupported")
  }
  membership_versions <- c(
    as.character(entity$membership_version),
    as.character(membership$membership_version),
    as.character(status$membership_version)
  )
  membership_version <- tryCatch(
    wgcna_group_assert_membership_version(membership_versions),
    error = function(e) fail(conditionMessage(e))
  )
  state_hashes <- unique(c(
    as.character(entity$state_sha256),
    as.character(membership$state_sha256)
  ))
  state_hashes <- state_hashes[!is.na(state_hashes) & nzchar(state_hashes)]
  if (length(state_hashes) != 1L || !identical(state_hashes[[1]], state_hash)) {
    fail("source-state SHA256 does not match the frozen state")
  }

  module_entities <- sort(as.character(
    entity$entity_id[entity$level == "module"]
  ))
  super_entities <- sort(as.character(
    entity$entity_id[entity$level == "supermodule"]
  ))
  state_modules <- sort(wgcna_identity_state_module_ids(readRDS(state_path)))
  if (!identical(module_entities, state_modules)) {
    fail("canonical module identities differ from the frozen-state module set")
  }
  if (!identical(module_entities, sort(unique(as.character(membership$module_id))))) {
    fail("module entity and membership sets differ")
  }
  if (!identical(
    super_entities,
    sort(unique(as.character(membership$supermodule_id)))
  )) {
    fail("supermodule entity and membership sets differ")
  }
  if (!identical(as.integer(status$n_modules[[1]]), length(module_entities)) ||
      !identical(
        as.integer(status$n_supermodules[[1]]),
        length(super_entities)
      )) {
    fail("status cardinality differs from current entity cardinality")
  }
  observed_counts <- table(membership$supermodule_id)
  declared_lookup <- stats::setNames(
    as.integer(entity$n_member_modules[entity$level == "supermodule"]),
    as.character(entity$entity_id[entity$level == "supermodule"])
  )
  count_ids <- sort(unique(c(names(observed_counts), names(declared_lookup))))
  if (!identical(
    as.integer(observed_counts[count_ids]),
    as.integer(declared_lookup[count_ids])
  )) {
    fail("declared supermodule member counts differ from membership")
  }

  list(
    dataset = dataset,
    entity = entity,
    membership = membership[order(
      membership$dataset, membership$module_id, membership$supermodule_id
    ), , drop = FALSE],
    status = status,
    paths = paths,
    file_hashes = hashes,
    identity_contract_sha256 = wgcna_group_contract_bundle_hash(hashes),
    membership_version = membership_version,
    contract_version = versions[[1]],
    frozen_state_sha256 = state_hash
  )
}

wgcna_group_module_bridge_error <- function(dataset, detail) {
  stop(
    "Frozen-state module eigengene bridge is incomplete, duplicated, or ",
    "ambiguous for ", dataset, ": ", detail, ".",
    call. = FALSE
  )
}

wgcna_group_metadata_eigengene_columns <- function(labels) {
  out <- rep(NA_character_, nrow(labels))
  if ("module_eigengene" %in% names(labels)) {
    explicit <- wgcna_identity_clean_id(labels$module_eigengene)
    out[!is.na(explicit)] <- explicit[!is.na(explicit)]
  }
  if ("WGCNAInternalColor" %in% names(labels)) {
    internal_color <- wgcna_identity_clean_id(labels$WGCNAInternalColor)
    use_color <- is.na(out) & !is.na(internal_color)
    out[use_color] <- paste0("ME", internal_color[use_color])
  }
  out
}

wgcna_group_classify_ignored_eigengenes <- function(
    dataset, extra_eigengene_cols, labels) {
  empty <- data.frame(
    state_eigengene_col_raw = character(),
    reason_excluded = character(),
    state_metadata_classification = character(),
    stringsAsFactors = FALSE
  )
  if (!length(extra_eigengene_cols)) return(empty)

  metadata_raw <- if (is.data.frame(labels) && nrow(labels)) {
    wgcna_group_metadata_eigengene_columns(labels)
  } else {
    character()
  }
  rows <- lapply(extra_eigengene_cols, function(raw_col) {
    metadata_rows <- which(!is.na(metadata_raw) & metadata_raw == raw_col)
    if (length(metadata_rows) > 1L) {
      wgcna_group_module_bridge_error(
        dataset,
        paste0("extra eigengene ", raw_col, " has ambiguous state metadata")
      )
    }

    metadata_module_id <- NA_character_
    metadata_color <- NA_character_
    metadata_eigengene <- NA_character_
    if (length(metadata_rows) == 1L) {
      row <- metadata_rows[[1]]
      if ("ModuleID" %in% names(labels)) {
        metadata_module_id <- wgcna_identity_clean_id(labels$ModuleID[[row]])
      }
      if (!is.na(metadata_module_id)) {
        wgcna_group_module_bridge_error(
          dataset,
          paste0(
            "extra eigengene ", raw_col,
            " has populated state metadata ModuleID ", metadata_module_id
          )
        )
      }
      if ("WGCNAInternalColor" %in% names(labels)) {
        metadata_color <- wgcna_identity_clean_id(
          labels$WGCNAInternalColor[[row]]
        )
      }
      if ("module_eigengene" %in% names(labels)) {
        metadata_eigengene <- wgcna_identity_clean_id(
          labels$module_eigengene[[row]]
        )
      }
    }

    unassigned_tokens <- c("grey", "gray", "unassigned")
    raw_token <- tolower(sub("^ME", "", raw_col))
    metadata_tokens <- tolower(c(
      metadata_color,
      sub("^ME", "", metadata_eigengene)
    ))
    metadata_tokens <- metadata_tokens[!is.na(metadata_tokens)]
    exact_known_unassigned <- raw_token %in% unassigned_tokens
    metadata_unassigned <- length(metadata_tokens) &&
      any(metadata_tokens %in% unassigned_tokens)
    if (!exact_known_unassigned && !metadata_unassigned) {
      wgcna_group_module_bridge_error(
        dataset,
        paste0("unknown non-contract eigengene ", raw_col)
      )
    }

    classification <- if (length(metadata_rows) == 1L) {
      paste0(
        "ModuleID=<empty>;WGCNAInternalColor=",
        ifelse(is.na(metadata_color), "<empty>", metadata_color),
        ";module_eigengene=",
        ifelse(is.na(metadata_eigengene), "<empty>", metadata_eigengene)
      )
    } else {
      paste0("exact_known_wgcna_unassigned_identity=", raw_col)
    }
    data.frame(
      state_eigengene_col_raw = raw_col,
      reason_excluded = if (metadata_unassigned) {
        "explicit_state_metadata_unassigned_module"
      } else {
        "exact_known_wgcna_unassigned_identity"
      },
      state_metadata_classification = classification,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

wgcna_group_build_module_bridge <- function(dataset, state, contract) {
  module_eigengenes <- extract_module_eigengenes(state)
  eigengene_cols <- names(module_eigengenes)
  eigengene_cols <- eigengene_cols[eigengene_cols != "Sample"]
  if (any(is.na(eigengene_cols) | !nzchar(eigengene_cols)) ||
      anyDuplicated(eigengene_cols)) {
    wgcna_group_module_bridge_error(
      dataset, "raw eigengene columns are missing or duplicated"
    )
  }

  contract_modules <- sort(unique(wgcna_identity_clean_id(
    contract$membership$module_id
  )))
  contract_modules <- contract_modules[!is.na(contract_modules)]
  if (!length(contract_modules)) {
    wgcna_group_module_bridge_error(dataset, "contract ModuleIDs are missing")
  }

  labels <- if (is.data.frame(state$module_label_table)) {
    as.data.frame(
      state$module_label_table,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame()
  }
  metadata_module_ids <- if ("ModuleID" %in% names(labels)) {
    wgcna_identity_clean_id(labels$ModuleID)
  } else {
    rep(NA_character_, nrow(labels))
  }
  has_authoritative_metadata <- any(!is.na(metadata_module_ids))

  if (has_authoritative_metadata) {
    non_contract_metadata_ids <- sort(unique(metadata_module_ids[
      !is.na(metadata_module_ids) &
        !metadata_module_ids %in% contract_modules
    ]))
    if (length(non_contract_metadata_ids)) {
      wgcna_group_module_bridge_error(
        dataset,
        paste0(
          "state module_label_table has populated ModuleID(s) outside the ",
          "contract (", paste(non_contract_metadata_ids, collapse = ", "), ")"
        )
      )
    }
    row_counts <- vapply(
      contract_modules,
      function(module_id) sum(metadata_module_ids == module_id, na.rm = TRUE),
      integer(1)
    )
    if (any(row_counts != 1L)) {
      invalid <- paste0(
        contract_modules[row_counts != 1L], "=", row_counts[row_counts != 1L]
      )
      wgcna_group_module_bridge_error(
        dataset,
        paste0(
          "state module_label_table must contain exactly one row per contract ",
          "ModuleID (", paste(invalid, collapse = ", "), ")"
        )
      )
    }

    label_rows <- labels[
      match(contract_modules, metadata_module_ids),
      ,
      drop = FALSE
    ]
    explicit_eigengene <- if ("module_eigengene" %in% names(label_rows)) {
      wgcna_identity_clean_id(label_rows$module_eigengene)
    } else {
      rep(NA_character_, nrow(label_rows))
    }
    internal_color <- if ("WGCNAInternalColor" %in% names(label_rows)) {
      wgcna_identity_clean_id(label_rows$WGCNAInternalColor)
    } else {
      rep(NA_character_, nrow(label_rows))
    }
    raw_cols <- explicit_eigengene
    bridge_method <- rep(
      "stable_state_metadata_module_eigengene",
      length(contract_modules)
    )
    use_internal_color <- is.na(raw_cols) & !is.na(internal_color)
    raw_cols[use_internal_color] <- paste0(
      "ME", internal_color[use_internal_color]
    )
    bridge_method[use_internal_color] <-
      "stable_state_metadata_internal_color"
    if (any(is.na(raw_cols))) {
      wgcna_group_module_bridge_error(
        dataset,
        "state metadata cannot resolve every contract ModuleID to an eigengene"
      )
    }

    raw_counts <- vapply(
      raw_cols,
      function(raw_col) sum(eigengene_cols == raw_col),
      integer(1)
    )
    if (any(raw_counts != 1L)) {
      invalid <- paste0(
        contract_modules[raw_counts != 1L], "->",
        raw_cols[raw_counts != 1L], "=", raw_counts[raw_counts != 1L]
      )
      wgcna_group_module_bridge_error(
        dataset,
        paste0(
          "state metadata must resolve each contract ModuleID to exactly one ",
          "existing raw eigengene column (", paste(invalid, collapse = ", "), ")"
        )
      )
    }
    out <- data.frame(
      dataset = dataset,
      module_id = contract_modules,
      state_eigengene_col_raw = raw_cols,
      state_eigengene_col_normalized = contract_modules,
      bridge_method = bridge_method,
      stringsAsFactors = FALSE
    )
  } else {
    bridge <- wgcna_identity_bridge_module_ids(
      eigengene_cols, contract_modules
    )
    resolved <- !is.na(bridge$normalized_id) &
      bridge$normalized_id %in% contract_modules
    out <- data.frame(
      dataset = dataset,
      module_id = bridge$normalized_id[resolved],
      state_eigengene_col_raw = bridge$original_id[resolved],
      state_eigengene_col_normalized = bridge$normalized_id[resolved],
      bridge_method = bridge$bridge_type[resolved],
      stringsAsFactors = FALSE
    )
  }

  out <- out[order(out$module_id), , drop = FALSE]
  if (nrow(out) != length(contract_modules) ||
      anyDuplicated(out$module_id) ||
      anyDuplicated(out$state_eigengene_col_raw) ||
      !identical(sort(out$module_id), contract_modules) ||
      any(!out$state_eigengene_col_raw %in% eigengene_cols) ||
      any(!out$bridge_method %in% c(
        "stable_state_metadata_module_eigengene",
        "stable_state_metadata_internal_color",
        "exact_current_identity_match",
        "syntactic_module_id_bridge_only"
      ))) {
    wgcna_group_module_bridge_error(
      dataset, "final biological module mapping is not one-to-one and complete"
    )
  }

  extra_eigengene_cols <- setdiff(
    eigengene_cols, out$state_eigengene_col_raw
  )
  ignored <- wgcna_group_classify_ignored_eigengenes(
    dataset, extra_eigengene_cols, labels
  )
  attr(out, "ignored_non_contract_eigengenes") <- ignored
  out
}

wgcna_group_clean_character <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | !nzchar(x) | toupper(x) %in%
      c("NA", "N/A", "NULL", "NONE", "NAN")] <- NA_character_
  x
}

wgcna_group_parse_animal_from_sample <- function(sample) {
  sample <- as.character(sample)
  parsed <- lapply(sample, function(value) {
    hits <- regmatches(
      toupper(value),
      gregexpr("(?<=_)A0*[0-9]+(?=_[LR]_)", toupper(value), perl = TRUE)
    )[[1]]
    hits <- unique(hits[hits != "-1"])
    if (length(hits) != 1L) {
      return(c(id = NA_character_, status = if (length(hits)) "ambiguous" else "missing"))
    }
    c(id = normalize_wgcna_animal_id(hits), status = "resolved")
  })
  data.frame(
    AnimalID = vapply(parsed, `[[`, character(1), "id"),
    parse_status = vapply(parsed, `[[`, character(1), "status"),
    stringsAsFactors = FALSE
  )
}

wgcna_group_parse_hemisphere_from_sample <- function(sample) {
  parsed <- vapply(as.character(sample), function(value) {
    hits <- regmatches(
      toupper(value),
      gregexpr("(?<=_)[LR](?=_[A-Z0-9]+_)", toupper(value), perl = TRUE)
    )[[1]]
    hits <- unique(hits[hits != "-1"])
    if (length(hits) == 1L) hits else NA_character_
  }, character(1))
  wgcna_group_clean_character(parsed)
}

wgcna_group_animal_provenance <- function(meta, sample) {
  aliases <- c(
    "AnimalID", "AnimalNum", "AnimalNumber", "AnimalNo", "Animal",
    "MouseID", "Mouse", "MouseNum", "MouseNumber", "mouse_id",
    "animal_id", "animal_num", "animal_number", "subject", "SubjectID",
    "donor", "DonorID", "animalid", "animalnum"
  )
  present <- names(meta)[tolower(gsub("[^a-z0-9]", "", names(meta))) %in%
                           tolower(gsub("[^a-z0-9]", "", aliases))]
  dedicated <- lapply(present, function(nm) normalize_wgcna_animal_id(meta[[nm]]))
  parsed <- wgcna_group_parse_animal_from_sample(sample)
  ids <- rep(NA_character_, nrow(meta))
  source <- rep(NA_character_, nrow(meta))
  status <- rep("missing", nrow(meta))
  for (i in seq_len(nrow(meta))) {
    values <- unique(stats::na.omit(vapply(
      dedicated,
      function(x) x[[i]],
      character(1)
    )))
    if (length(values) == 1L) {
      ids[[i]] <- values[[1]]
      source[[i]] <- paste0("dedicated_metadata_column:", paste(present, collapse = ";"))
      status[[i]] <- "resolved"
    } else if (length(values) > 1L) {
      source[[i]] <- paste0("ambiguous_dedicated_metadata_columns:", paste(present, collapse = ";"))
      status[[i]] <- "ambiguous"
    } else if (identical(parsed$parse_status[[i]], "resolved")) {
      ids[[i]] <- parsed$AnimalID[[i]]
      source[[i]] <- "deterministic_sample_name_parsing"
      status[[i]] <- "resolved"
    } else {
      source[[i]] <- if (identical(parsed$parse_status[[i]], "ambiguous")) {
        "ambiguous_sample_name_parsing"
      } else {
        "unavailable"
      }
      status[[i]] <- parsed$parse_status[[i]]
    }
  }
  data.frame(
    AnimalID = ids,
    AnimalID_source = source,
    AnimalID_mapping_status = status,
    stringsAsFactors = FALSE
  )
}

wgcna_group_standardize_stress <- function(x) {
  x <- toupper(wgcna_group_clean_character(x))
  dplyr::case_when(
    grepl("^CON|^CTRL|CONTROL", x, ignore.case = TRUE) ~ "CON",
    grepl("^RES", x, ignore.case = TRUE) ~ "RES",
    grepl("^SUS", x, ignore.case = TRUE) ~ "SUS",
    TRUE ~ x
  )
}

wgcna_group_sample_audit <- function(dataset, state, module_eigengenes) {
  meta <- as.data.frame(state$sample_info, check.names = FALSE, stringsAsFactors = FALSE)
  sample_col <- first_present_col(
    meta, c("Sample", "sample", "SampleID", "sample_id", "row.names")
  )
  metadata_sample <- if (is.na(sample_col)) rownames(meta) else as.character(meta[[sample_col]])
  state_sample <- as.character(module_eigengenes$Sample)
  if (anyNA(state_sample) || any(!nzchar(state_sample)) || anyDuplicated(state_sample)) {
    stop("State eigengene Sample keys must be unique and nonmissing.", call. = FALSE)
  }
  if (anyNA(metadata_sample) || any(!nzchar(metadata_sample)) ||
      anyDuplicated(metadata_sample)) {
    stop("State metadata Sample keys must be unique and nonmissing.", call. = FALSE)
  }
  if (!setequal(state_sample, metadata_sample) ||
      length(state_sample) != length(metadata_sample)) {
    stop(
      "State eigengenes and metadata do not have exact one-to-one Sample coverage.",
      call. = FALSE
    )
  }
  idx <- match(state_sample, metadata_sample)
  meta <- meta[idx, , drop = FALSE]
  animal <- wgcna_group_animal_provenance(meta, state_sample)
  value <- function(candidates) {
    col <- first_present_col(meta, candidates)
    if (is.na(col)) rep(NA_character_, nrow(meta)) else as.character(meta[[col]])
  }
  region <- wgcna_group_clean_character(value(c("Region", "region")))
  layer <- wgcna_group_clean_character(value(c("Layer", "layer")))
  group <- wgcna_group_standardize_stress(value(
    c("StressGroup", "ExpGroup", "condition", "Group", "group", "group2")
  ))
  sex <- wgcna_group_clean_character(value(c("Sex", "sex")))
  batch <- wgcna_group_clean_character(value(
    c("Batch", "batch", "plate", "run", "batch_id")
  ))
  hemisphere <- toupper(wgcna_group_clean_character(value(
    c("Hemisphere", "hemisphere", "Side", "side", "Laterality", "laterality")
  )))
  parsed_hemisphere <- wgcna_group_parse_hemisphere_from_sample(state_sample)
  hemisphere[is.na(hemisphere)] <- parsed_hemisphere[is.na(hemisphere)]
  spatial <- if (dataset == "neuron_neuropil") {
    ifelse(
      !is.na(region) & !is.na(layer),
      paste(tolower(region), tolower(layer), sep = "_"),
      NA_character_
    )
  } else {
    tolower(region)
  }
  reasons <- lapply(seq_along(state_sample), function(i) {
    x <- character()
    if (!identical(state_sample[[i]], metadata_sample[[idx[[i]]]])) {
      x <- c(x, "metadata_sample_mismatch")
    }
    if (!group[[i]] %in% c("CON", "RES", "SUS")) {
      x <- c(x, "missing_or_invalid_StressGroup")
    }
    if (is.na(spatial[[i]]) || !nzchar(spatial[[i]])) {
      x <- c(x, "missing_spatial_unit")
    }
    if (!hemisphere[[i]] %in% c("L", "R")) {
      x <- c(x, "missing_or_invalid_Hemisphere")
    }
    if (!identical(animal$AnimalID_mapping_status[[i]], "resolved") ||
        is.na(animal$AnimalID[[i]])) {
      x <- c(x, "missing_or_ambiguous_AnimalID")
    }
    if (length(x)) paste(x, collapse = ";") else "none"
  })
  audit <- data.frame(
    dataset = dataset,
    Sample = state_sample,
    metadata_Sample = metadata_sample[idx],
    metadata_match_status = "exact_one_to_one_match",
    inclusion_status = ifelse(
      vapply(reasons, identical, logical(1), "none"), "included", "excluded"
    ),
    exclusion_reason = unlist(reasons, use.names = FALSE),
    AnimalID = animal$AnimalID,
    AnimalID_source = animal$AnimalID_source,
    AnimalID_mapping_status = animal$AnimalID_mapping_status,
    StressGroup = group,
    Region = region,
    Layer = layer,
    Hemisphere = hemisphere,
    canonical_spatial_unit = spatial,
    SpatialUnitType = if (dataset == "neuron_neuropil") "region_layer" else "region",
    Sex = sex,
    Batch = batch,
    stringsAsFactors = FALSE
  )
  if (any(audit$inclusion_status != "included")) {
    stop(
      "Sample/biological-replicate preflight excluded ",
      sum(audit$inclusion_status != "included"),
      " state sample(s).",
      call. = FALSE
    )
  }
  animal_groups <- tapply(
    audit$StressGroup, audit$AnimalID,
    function(x) length(unique(stats::na.omit(x)))
  )
  if (any(animal_groups != 1L)) {
    stop("At least one AnimalID maps to multiple StressGroups.", call. = FALSE)
  }
  audit
}

wgcna_group_supermodule_audits <- function(
    dataset, module_eigengenes, super_map, bridge, contract
) {
  ids <- sort(unique(as.character(super_map$SupermoduleID)))
  provenance <- list()
  eigenvalues <- list()
  loadings <- list()
  input_audit <- list()
  for (sid in ids) {
    rows <- super_map[super_map$SupermoduleID == sid, , drop = FALSE]
    rows <- rows[order(rows$module_id), , drop = FALSE]
    members <- as.character(rows$module_eigengene)
    module_ids <- as.character(rows$module_id)
    vals <- module_eigengenes[, members, drop = FALSE]
    variable <- vapply(vals, function(x) {
      value <- stats::var(as.numeric(x), na.rm = TRUE)
      is.finite(value) && value > 0
    }, logical(1))
    if (!all(variable)) {
      stop("Canonical supermodule ", sid, " contains a zero-variance eigengene.",
           call. = FALSE)
    }
    method <- if (length(members) == 1L) {
      "singleton_member_module_eigengene"
    } else {
      "pc1_centered_scaled_oriented_to_positive_member_mean"
    }
    if (length(members) == 1L) {
      pca <- NULL
      variance <- 1
      corr_before <- 1
      corr_after <- 1
      multiplier <- 1
      eigenvalues[[length(eigenvalues) + 1L]] <- data.frame(
        dataset = dataset, supermodule_id = sid, pc = 1L,
        eigenvalue = 1, variance_explained = 1,
        cumulative_variance_explained = 1,
        n_member_modules = 1L, n_samples_used = nrow(vals),
        pca_status = "singleton_identity",
        membership_version = contract$membership_version,
        identity_contract_version = contract$contract_version,
        identity_contract_sha256 = contract$identity_contract_sha256,
        frozen_state_sha256 = contract$frozen_state_sha256,
        stringsAsFactors = FALSE
      )
      loadings[[length(loadings) + 1L]] <- data.frame(
        dataset = dataset, supermodule_id = sid,
        module_id = module_ids, module_eigengene = members,
        pc = 1L, loading = 1, abs_loading = 1, loading_rank = 1L,
        n_member_modules = 1L, n_samples_used = nrow(vals),
        pca_status = "singleton_identity",
        membership_version = contract$membership_version,
        identity_contract_version = contract$contract_version,
        identity_contract_sha256 = contract$identity_contract_sha256,
        frozen_state_sha256 = contract$frozen_state_sha256,
        stringsAsFactors = FALSE
      )
    } else {
      pca <- stats::prcomp(vals, center = TRUE, scale. = TRUE)
      variance_all <- pca$sdev^2 / sum(pca$sdev^2)
      corr_before <- stats::cor(
        pca$x[, 1L], rowMeans(vals), use = "pairwise.complete.obs"
      )
      multiplier <- if (is.finite(corr_before) && corr_before < 0) -1 else 1
      corr_after <- corr_before * multiplier
      variance <- variance_all[[1]]
      eigenvalues[[length(eigenvalues) + 1L]] <- data.frame(
        dataset = dataset, supermodule_id = sid,
        pc = seq_along(variance_all),
        eigenvalue = pca$sdev^2,
        variance_explained = variance_all,
        cumulative_variance_explained = cumsum(variance_all),
        n_member_modules = length(members), n_samples_used = nrow(vals),
        pca_status = "multi_module_pc1",
        membership_version = contract$membership_version,
        identity_contract_version = contract$contract_version,
        identity_contract_sha256 = contract$identity_contract_sha256,
        frozen_state_sha256 = contract$frozen_state_sha256,
        stringsAsFactors = FALSE
      )
      oriented_loading <- pca$rotation[, 1L] * multiplier
      loadings[[length(loadings) + 1L]] <- data.frame(
        dataset = dataset, supermodule_id = sid,
        module_id = module_ids, module_eigengene = members,
        pc = 1L, loading = as.numeric(oriented_loading),
        abs_loading = abs(as.numeric(oriented_loading)),
        loading_rank = rank(-abs(as.numeric(oriented_loading)), ties.method = "first"),
        n_member_modules = length(members), n_samples_used = nrow(vals),
        pca_status = "multi_module_pc1",
        membership_version = contract$membership_version,
        identity_contract_version = contract$contract_version,
        identity_contract_sha256 = contract$identity_contract_sha256,
        frozen_state_sha256 = contract$frozen_state_sha256,
        stringsAsFactors = FALSE
      )
    }
    provenance[[length(provenance) + 1L]] <- data.frame(
      dataset = dataset, level = "supermodule", endpoint_id = sid,
      canonical_entity_id = sid, module_id = NA_character_,
      supermodule_id = sid, state_eigengene_col_raw = NA_character_,
      state_eigengene_col_normalized = NA_character_,
      bridge_method = "canonical_membership_contract",
      member_modules = paste(module_ids, collapse = ";"),
      member_eigengene_columns = paste(members, collapse = ";"),
      n_member_modules = length(members),
      endpoint_col = paste0("SM__", safe_filename(sid)),
      endpoint_construction_method = method,
      pc1_variance_explained = variance,
      pc1_sign_before_orientation = ifelse(corr_before < 0, -1, 1),
      pc1_sign_after_orientation = 1,
      orientation_multiplier = multiplier,
      orientation_correlation_before = corr_before,
      orientation_correlation_after = corr_after,
      membership_version = contract$membership_version,
      identity_contract_version = contract$contract_version,
      identity_contract_sha256 = contract$identity_contract_sha256,
      frozen_state_sha256 = contract$frozen_state_sha256,
      endpoint_provenance_status = "validated_canonical_endpoint",
      stringsAsFactors = FALSE
    )
    input_audit[[length(input_audit) + 1L]] <- data.frame(
      dataset = dataset, supermodule_id = sid,
      n_contract_member_modules = length(module_ids),
      n_matched_member_eigengenes = length(members),
      member_modules = paste(module_ids, collapse = ";"),
      member_eigengene_columns = paste(members, collapse = ";"),
      all_members_present = TRUE,
      membership_version = contract$membership_version,
      identity_contract_version = contract$contract_version,
      frozen_state_sha256 = contract$frozen_state_sha256,
      identity_contract_sha256 = contract$identity_contract_sha256,
      stringsAsFactors = FALSE
    )
  }
  list(
    provenance = dplyr::bind_rows(provenance),
    eigenvalues = dplyr::bind_rows(eigenvalues),
    loadings = dplyr::bind_rows(loadings),
    input_audit = dplyr::bind_rows(input_audit)
  )
}

wgcna_group_construct_endpoints <- function(
    dataset, state, contract, bridge
) {
  module_eigengenes <- extract_module_eigengenes(state)
  super_map <- contract$membership |>
    dplyr::select("dataset", "module_id", "supermodule_id") |>
    dplyr::left_join(
      bridge[, c("module_id", "state_eigengene_col_raw"), drop = FALSE],
      by = "module_id",
      relationship = "many-to-one"
    ) |>
    dplyr::transmute(
      dataset = .data$dataset,
      ModuleID = .data$module_id,
      module_id = .data$module_id,
      SupermoduleID = .data$supermodule_id,
      module_eigengene = .data$state_eigengene_col_raw
    ) |>
    dplyr::arrange(.data$dataset, .data$SupermoduleID, .data$module_id)
  super <- make_supermodule_eigengenes(module_eigengenes, super_map)
  audits <- wgcna_group_supermodule_audits(
    dataset, module_eigengenes, super_map, bridge, contract
  )
  module_provenance <- bridge |>
    dplyr::transmute(
      dataset = .data$dataset, level = "module",
      endpoint_id = .data$module_id,
      canonical_entity_id = .data$module_id,
      module_id = .data$module_id, supermodule_id = NA_character_,
      state_eigengene_col_raw = .data$state_eigengene_col_raw,
      state_eigengene_col_normalized = .data$state_eigengene_col_normalized,
      bridge_method = .data$bridge_method,
      member_modules = .data$module_id,
      member_eigengene_columns = .data$state_eigengene_col_raw,
      n_member_modules = 1L,
      endpoint_col = .data$state_eigengene_col_raw,
      endpoint_construction_method = "frozen_state_module_eigengene",
      pc1_variance_explained = 1,
      pc1_sign_before_orientation = 1,
      pc1_sign_after_orientation = 1,
      orientation_multiplier = 1,
      orientation_correlation_before = 1,
      orientation_correlation_after = 1,
      membership_version = contract$membership_version,
      identity_contract_version = contract$contract_version,
      identity_contract_sha256 = contract$identity_contract_sha256,
      frozen_state_sha256 = contract$frozen_state_sha256,
      endpoint_provenance_status = "validated_canonical_endpoint"
    )
  provenance <- dplyr::bind_rows(module_provenance, audits$provenance)

  composition <- super_map |>
    dplyr::group_by(.data$dataset, .data$SupermoduleID) |>
    dplyr::summarise(
      supermodule_id = dplyr::first(.data$SupermoduleID),
      supermodule_label = dplyr::first(.data$SupermoduleID),
      supermodule_eigengene = paste0(
        "SM__", safe_filename(dplyr::first(.data$SupermoduleID))
      ),
      n_member_modules = dplyr::n(),
      member_modules = paste(sort(.data$module_id), collapse = ";"),
      member_eigengene_columns = paste(
        .data$module_eigengene[order(.data$module_id)], collapse = ";"
      ),
      endpoint_construction_method = dplyr::if_else(
        dplyr::n() == 1L,
        "singleton_member_module_eigengene",
        "pc1_centered_scaled_oriented_to_positive_member_mean"
      ),
      supermodule_merge_rule = "canonical_identity_contract",
      supermodule_merge_cut_height = NA_real_,
      membership_version = contract$membership_version,
      identity_contract_version = contract$contract_version,
      identity_contract_sha256 = contract$identity_contract_sha256,
      frozen_state_sha256 = contract$frozen_state_sha256,
      endpoint_provenance_status = "validated_canonical_endpoint",
      .groups = "drop"
    )
  expected <- contract$membership |>
    dplyr::arrange(.data$dataset, .data$module_id, .data$supermodule_id) |>
    dplyr::select("dataset", "module_id", "supermodule_id")
  observed <- composition |>
    dplyr::select("dataset", "supermodule_id", "member_modules") |>
    tidyr::separate_longer_delim("member_modules", delim = ";") |>
    dplyr::rename(module_id = "member_modules") |>
    dplyr::arrange(.data$dataset, .data$module_id, .data$supermodule_id) |>
    dplyr::select("dataset", "module_id", "supermodule_id")
  if (!identical(as.data.frame(observed), as.data.frame(expected))) {
    stop("Published supermodule composition does not exactly reproduce the identity contract.",
         call. = FALSE)
  }
  list(
    module_eigengenes = module_eigengenes,
    supermodule_eigengenes = super$eigengenes,
    module_map = module_provenance[, c(
      "endpoint_col", "endpoint_id", "endpoint_construction_method",
      "endpoint_provenance_status", "n_member_modules", "member_modules"
    )],
    supermodule_map = audits$provenance[, c(
      "endpoint_col", "endpoint_id", "endpoint_construction_method",
      "endpoint_provenance_status", "n_member_modules", "member_modules"
    )],
    composition = composition,
    provenance = provenance,
    pca_input_audit = audits$input_audit,
    pca_eigenvalues = audits$eigenvalues,
    pca_member_loadings = audits$loadings
  )
}

wgcna_group_count_string <- function(values, groups) {
  tab <- table(groups, useNA = "no")
  if (identical(values, "animals")) {
    tab <- tapply(
      groups$AnimalID, groups$StressGroup,
      function(x) length(unique(stats::na.omit(x)))
    )
  }
  paste(names(tab), as.integer(tab), sep = "=", collapse = ";")
}

wgcna_group_counts <- function(dat) {
  samples_group <- table(dat$StressGroup)
  animals_group <- tapply(
    dat$AnimalID, dat$StressGroup,
    function(x) length(unique(stats::na.omit(x)))
  )
  sample_spatial <- table(dat$SpatialUnit)
  animal_spatial <- tapply(
    dat$AnimalID, dat$SpatialUnit,
    function(x) length(unique(stats::na.omit(x)))
  )
  list(
    n_samples = nrow(dat),
    n_animals = length(unique(stats::na.omit(dat$AnimalID))),
    n_samples_per_group = paste(
      names(samples_group), as.integer(samples_group), sep = "=", collapse = ";"
    ),
    n_animals_per_group = paste(
      names(animals_group), as.integer(animals_group), sep = "=", collapse = ";"
    ),
    n_samples_per_spatial_unit = paste(
      names(sample_spatial), as.integer(sample_spatial), sep = "=", collapse = ";"
    ),
    n_animals_per_spatial_unit = paste(
      names(animal_spatial), as.integer(animal_spatial), sep = "=", collapse = ";"
    ),
    min_animals_per_group = if (length(animals_group)) {
      min(as.integer(animals_group))
    } else {
      0L
    }
  )
}

wgcna_group_varying_covariates <- function(dat, candidates) {
  present <- intersect(candidates, names(dat))
  present[vapply(present, function(nm) {
    length(unique(stats::na.omit(dat[[nm]]))) > 1L
  }, logical(1))]
}

wgcna_group_capture <- function(expr) {
  warnings <- character()
  value <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = identity
  )
  list(value = value, warnings = unique(warnings))
}

wgcna_group_named_contrast_weights <- function(
    observed_levels, numerator, denominator
) {
  observed_levels <- as.character(observed_levels)
  if (!all(c(numerator, denominator) %in% observed_levels)) return(NULL)
  weights <- stats::setNames(rep(0, length(observed_levels)), observed_levels)
  weights[[numerator]] <- 1
  weights[[denominator]] <- -1
  weights[observed_levels]
}

wgcna_group_predeclared_contrasts <- function(observed_levels) {
  pairs <- list(
    "RES - CON" = c("RES", "CON"),
    "SUS - CON" = c("SUS", "CON"),
    "SUS - RES" = c("SUS", "RES")
  )
  methods <- list()
  for (label in names(pairs)) {
    weights <- wgcna_group_named_contrast_weights(
      observed_levels, pairs[[label]][[1]], pairs[[label]][[2]]
    )
    if (!is.null(weights)) methods[[label]] <- weights
  }
  methods
}

wgcna_group_emmeans_contrast_summary <- function(
    emm_grid, methods, ci_level = 0.95
) {
  if (!is.numeric(ci_level) || length(ci_level) != 1L ||
      !is.finite(ci_level) || ci_level <= 0 || ci_level >= 1) {
    stop("ci_level must be one finite number strictly between zero and one.",
         call. = FALSE)
  }
  contrast_grid <- emmeans::contrast(
    emm_grid, method = methods, adjust = "none"
  )
  out <- as.data.frame(summary(
    contrast_grid, infer = c(TRUE, TRUE), level = ci_level,
    adjust = "none"
  ))
  lower_name <- intersect(c("lower.CL", "asymp.LCL"), names(out))
  upper_name <- intersect(c("upper.CL", "asymp.UCL"), names(out))
  if (length(lower_name) != 1L || length(upper_name) != 1L) {
    stop(
      "emmeans contrast summary did not return one confidence-interval pair.",
      call. = FALSE
    )
  }
  out$CI_low <- suppressWarnings(as.numeric(out[[lower_name[[1L]]]]))
  out$CI_high <- suppressWarnings(as.numeric(out[[upper_name[[1L]]]]))
  out$CI_method <- "emmeans_contrast_summary"
  out$CI_level <- as.numeric(ci_level)
  out$CI_df_method <- if ("df" %in% names(out)) {
    ifelse(
      is.finite(suppressWarnings(as.numeric(out$df))),
      "finite_df_from_emmeans", "asymptotic_from_emmeans"
    )
  } else {
    "asymptotic_from_emmeans"
  }
  out
}

wgcna_group_model_type_for_scope <- function(animal_ids) {
  animal_ids <- wgcna_group_clean_character(animal_ids)
  if (anyNA(animal_ids)) {
    stop("AnimalID is missing or ambiguous.", call. = FALSE)
  }
  if (any(table(animal_ids) > 1L)) "lmerTest_lmer" else "lm"
}

wgcna_group_model_validation_schema <- function() {
  data.frame(
    dataset = character(), level = character(), endpoint_id = character(),
    effect_scope = character(), spatial_unit = character(),
    attempt_status = character(), failure_reason = character(),
    requested_formula = character(), actual_formula = character(),
    model_type = character(), fixed_effect_rank = integer(),
    fixed_effect_columns = integer(), rank_deficient_model = logical(),
    singular_model = logical(), convergence_status = character(),
    convergence_warnings = character(), optimizer_messages = character(),
    optimizer_code = character(), model_warnings = character(),
    residual_df = numeric(), n_samples = integer(), n_animals = integer(),
    n_samples_per_group = character(), n_animals_per_group = character(),
    n_samples_per_spatial_unit = character(),
    n_animals_per_spatial_unit = character(),
    n_contrasts_requested = integer(), n_contrasts_eligible = integer(),
    n_contrasts_emitted = integer(), excluded_contrasts = character(),
    emmeans_status = character(), has_repeated_animals = logical(),
    animal_random_effect_used = logical(), AnimalID_source = character(),
    model_valid_for_inference = logical(),
    model_diagnostic_scope = character(),
    model_stability_status = character(),
    primary_model_stable = logical(), claim_allowed_model = logical(),
    enters_primary_fdr = logical(),
    manuscript_claim_ready = character(),
    random_effect_structure = character(),
    random_intercept_variance = numeric(), residual_variance = numeric(),
    variance_ratio = numeric(), ICC = numeric(),
    is_singular_lme4 = logical(), boundary_by_variance_ratio = logical(),
    boundary_variance_ratio_tolerance = numeric(),
    lme4_singularity_tolerance = numeric(),
    singularity_diagnostic_status = character(),
    diagnostic_review_required = logical(),
    singularity_class = character(), boundary_warning = character(),
    reduced_formula = character(), full_formula = character(),
    identical_rows_verified = logical(),
    reduced_row_hash = character(), full_row_hash = character(),
    reduced_fixed_effect_rank = integer(),
    reduced_fixed_effect_columns = integer(),
    reduced_optimizer_code = character(),
    reduced_optimizer_messages = character(),
    reduced_convergence_status = character(),
    reduced_convergence_warnings = character(),
    reduced_random_effect_structure = character(),
    reduced_random_intercept_variance = numeric(),
    reduced_residual_variance = numeric(),
    reduced_variance_ratio = numeric(), reduced_ICC = numeric(),
    reduced_is_singular_lme4 = logical(),
    reduced_boundary_by_variance_ratio = logical(),
    reduced_singularity_class = character(),
    likelihood_ratio_statistic = numeric(),
    likelihood_ratio_df = numeric(), likelihood_ratio_p_value = numeric(),
    reduced_model_stability_status = character(),
    reduced_diagnostic_review_required = logical(),
    full_fixed_effect_rank = integer(),
    full_fixed_effect_columns = integer(),
    full_optimizer_code = character(),
    full_optimizer_messages = character(),
    full_convergence_status = character(),
    full_convergence_warnings = character(),
    full_random_effect_structure = character(),
    full_random_intercept_variance = numeric(),
    full_residual_variance = numeric(),
    full_variance_ratio = numeric(), full_ICC = numeric(),
    full_is_singular_lme4 = logical(),
    full_boundary_by_variance_ratio = logical(),
    full_singularity_class = character(),
    full_model_stability_status = character(),
    full_diagnostic_review_required = logical(),
    membership_version = character(), identity_contract_version = character(),
    identity_contract_sha256 = character(),
    frozen_state_sha256 = character(),
    stringsAsFactors = FALSE
  )
}

wgcna_group_primary_schema <- function() {
  data.frame(
    dataset = character(), level = character(), endpoint_id = character(),
    endpoint_label = character(), module_id = character(),
    supermodule_id = character(), module_label = character(),
    supermodule_label = character(),
    hypothesis_level = character(),
    canonical_claim_entity_id = character(),
    claim_entity_role = character(), support_source_entity_id = character(),
    independent_hypothesis = logical(), enters_primary_fdr = logical(),
    display_allowed = logical(), manuscript_claim_ready = character(),
    analysis_tier = character(), primary_analysis_tier = character(),
    primary_contrast = logical(), test_type = character(),
    spatial_unit = character(),
    effect_scope = character(), SpatialUnitType = character(),
    model_type = character(), has_repeated_animals = logical(),
    n_animals = integer(), n_animals_total = integer(),
    n_animals_per_group = character(), min_animals_per_group = integer(),
    min_unique_animals_compared_group = integer(),
    n_samples_total = integer(), n_samples_per_group = character(),
    n_samples_per_spatial_unit = character(),
    n_animals_per_spatial_unit = character(),
    animal_level_status = character(), pseudoreplication_guard = character(),
    contrast = character(), estimate = numeric(), SE = numeric(),
    CI_low = numeric(), CI_high = numeric(),
    CI_method = character(), CI_level = numeric(),
    CI_df_method = character(),
    statistic = numeric(), df_num = numeric(), df_den = numeric(),
    p_value = numeric(),
    FDR_primary_global = numeric(), FDR_primary_family_id = character(),
    n_tests_FDR_primary = integer(),
    FDR_secondary_global = numeric(), FDR_secondary_family_id = character(),
    n_tests_FDR_secondary_global = integer(),
    FDR_interaction_omnibus = numeric(),
    FDR_interaction_family_id = character(),
    n_tests_FDR_interaction_omnibus = integer(),
    FDR_local_exploratory = numeric(), FDR_local_family_id = character(),
    n_tests_FDR_local_exploratory = integer(),
    FDR_conservative_all_tests = numeric(),
    FDR_conservative_family_id = character(),
    n_tests_FDR_conservative_all_tests = integer(),
    FDR_within_dataset_level = numeric(),
    FDR_dataset_all_levels = numeric(), FDR_global = numeric(),
    FDR_family_within_level_id = character(),
    FDR_family_dataset_id = character(),
    n_tests_FDR_within_dataset_level = integer(),
    n_tests_FDR_dataset_all_levels = integer(), FDR_method = character(),
    evidence_status = character(), statistical_support_status = character(),
    direction = character(),
    n_samples = integer(), formula_requested = character(),
    formula_used = character(), dropped_covariates = character(),
    model_family = character(), model_formula = character(),
    model_valid_for_inference = logical(),
    model_diagnostic_scope = character(),
    model_stability_status = character(),
    primary_model_stable = logical(), claim_allowed_model = logical(),
    model_downgrade_reason = character(), fallback_used = logical(),
    fallback_type = character(), rank_deficient_model = logical(),
    singular_model = logical(), emmeans_success = logical(),
    animal_random_effect_used = logical(),
    biological_replicate_unit = character(), model_warning = character(),
    fixed_effect_rank = integer(), fixed_effect_columns = integer(),
    convergence_status = character(), convergence_warnings = character(),
    optimizer_messages = character(), optimizer_code = character(),
    random_effect_structure = character(),
    random_intercept_variance = numeric(), residual_variance = numeric(),
    variance_ratio = numeric(), ICC = numeric(),
    is_singular_lme4 = logical(), boundary_by_variance_ratio = logical(),
    boundary_variance_ratio_tolerance = numeric(),
    lme4_singularity_tolerance = numeric(),
    singularity_diagnostic_status = character(),
    diagnostic_review_required = logical(),
    singularity_class = character(), boundary_warning = character(),
    residual_df = numeric(), AnimalID_source = character(),
    animal_id_mapping_status = character(),
    reduced_formula = character(), full_formula = character(),
    identical_rows_verified = logical(),
    reduced_row_hash = character(), full_row_hash = character(),
    reduced_fixed_effect_rank = integer(),
    reduced_fixed_effect_columns = integer(),
    reduced_optimizer_code = character(),
    reduced_optimizer_messages = character(),
    reduced_convergence_status = character(),
    reduced_convergence_warnings = character(),
    reduced_random_effect_structure = character(),
    reduced_random_intercept_variance = numeric(),
    reduced_residual_variance = numeric(),
    reduced_variance_ratio = numeric(), reduced_ICC = numeric(),
    reduced_is_singular_lme4 = logical(),
    reduced_boundary_by_variance_ratio = logical(),
    reduced_singularity_class = character(),
    likelihood_ratio_statistic = numeric(),
    likelihood_ratio_df = numeric(), likelihood_ratio_p_value = numeric(),
    reduced_model_stability_status = character(),
    reduced_diagnostic_review_required = logical(),
    full_fixed_effect_rank = integer(),
    full_fixed_effect_columns = integer(),
    full_optimizer_code = character(),
    full_optimizer_messages = character(),
    full_convergence_status = character(),
    full_convergence_warnings = character(),
    full_random_effect_structure = character(),
    full_random_intercept_variance = numeric(),
    full_residual_variance = numeric(),
    full_variance_ratio = numeric(), full_ICC = numeric(),
    full_is_singular_lme4 = logical(),
    full_boundary_by_variance_ratio = logical(),
    full_singularity_class = character(),
    full_model_stability_status = character(),
    full_diagnostic_review_required = logical(),
    membership_version = character(),
    identity_contract_version = character(),
    identity_contract_sha256 = character(),
    frozen_state_sha256 = character(),
    endpoint_construction_method = character(),
    endpoint_provenance_status = character(),
    stringsAsFactors = FALSE
  )
}

wgcna_group_hemisphere_value_schema <- function() {
  data.frame(
    dataset = character(), level = character(), endpoint_id = character(),
    AnimalID = character(), StressGroup = character(),
    SpatialUnit = character(), canonical_spatial_unit = character(),
    Hemisphere = character(), contributing_sample_ids = character(),
    source_sample_count = integer(), source_values = character(),
    hemisphere_value = numeric(), n_hemisphere_source_rows = integer(),
    duplicate_source_row_status = character(),
    aggregation_method = character(), AnimalID_source = character(),
    Region = character(), Layer = character(),
    SpatialUnitType = character(), stringsAsFactors = FALSE
  )
}

wgcna_group_animal_spatial_value_schema <- function() {
  data.frame(
    AnimalID = character(), StressGroup = character(),
    canonical_spatial_unit = character(), eigengene = numeric(),
    n_source_samples = integer(), n_hemispheres_observed = integer(),
    hemispheres_observed = character(), contributing_samples = character(),
    ordered_hemispheres = character(),
    ordered_contributing_sample_ids = character(),
    ordered_hemisphere_source_values = character(),
    ordered_hemisphere_means = character(),
    ordered_hemisphere_source_row_counts = character(),
    ordered_hemisphere_duplicate_statuses = character(),
    hemisphere_values = character(), AnimalID_source = character(),
    Region = character(), Layer = character(),
    SpatialUnitType = character(), dataset = character(),
    level = character(), endpoint_id = character(),
    SpatialUnit = character(), bilateral_pair_status = character(),
    hemisphere_completeness_status = character(),
    aggregation_method = character(),
    duplicate_source_row_status = character(),
    aggregated_row_id = character(), aggregated_row_sha256 = character(),
    aggregated_row_hash_contract = character(),
    stringsAsFactors = FALSE
  )
}

wgcna_group_bind_schema <- function(rows, schema) {
  if (!length(rows)) return(schema)
  out <- dplyr::bind_rows(rows)
  for (nm in setdiff(names(schema), names(out))) {
    prototype <- schema[[nm]]
    value <- if (is.logical(prototype)) {
      NA
    } else if (is.integer(prototype)) {
      NA_integer_
    } else if (is.numeric(prototype)) {
      NA_real_
    } else {
      NA_character_
    }
    out[[nm]] <- rep(value, nrow(out))
  }
  out[, names(schema), drop = FALSE]
}

wgcna_group_fit_attempt_phase2_legacy <- function(
    dat, dataset, level, endpoint_id, endpoint_col, construction_method,
    provenance_status, effect_scope, contract, spatial_unit = NA_character_
) {
  if (!is.na(spatial_unit)) {
    dat <- dat[dat$SpatialUnit == spatial_unit, , drop = FALSE]
  }
  dat$eigengene <- as.numeric(dat[[endpoint_col]])
  dat <- dat[
    is.finite(dat$eigengene) &
      dat$StressGroup %in% c("CON", "RES", "SUS"),
    ,
    drop = FALSE
  ]
  dat$StressGroup <- droplevels(factor(dat$StressGroup))
  dat$SpatialUnit <- droplevels(factor(dat$SpatialUnit))
  counts <- wgcna_group_counts(dat)
  model_type <- wgcna_group_model_type_for_scope(dat$AnimalID)
  repeated <- identical(model_type, "lmerTest_lmer")
  covariates <- wgcna_group_varying_covariates(dat, c("Sex", "Batch"))
  fixed_rhs <- switch(
    effect_scope,
    within_spatial_unit = c("StressGroup", covariates),
    spatial_adjusted_global = c("StressGroup", "SpatialUnit", covariates),
    stress_by_spatial_interaction = c("StressGroup * SpatialUnit", covariates)
  )
  fixed_rhs <- unique(fixed_rhs)
  fixed_formula <- stats::as.formula(
    paste("eigengene ~", paste(fixed_rhs, collapse = " + "))
  )
  requested <- paste(
    deparse(fixed_formula),
    if (repeated) "+ (1 | AnimalID)" else ""
  )
  validation_base <- list(
    dataset = dataset, level = level, endpoint_id = endpoint_id,
    effect_scope = effect_scope,
    spatial_unit = ifelse(is.na(spatial_unit), {
      if (effect_scope == "spatial_adjusted_global") {
        "global_spatial_adjusted"
      } else {
        "all_spatial_units"
      }
    }, spatial_unit),
    requested_formula = requested,
    n_samples = counts$n_samples, n_animals = counts$n_animals,
    n_samples_per_group = counts$n_samples_per_group,
    n_animals_per_group = counts$n_animals_per_group,
    n_samples_per_spatial_unit = counts$n_samples_per_spatial_unit,
    n_animals_per_spatial_unit = counts$n_animals_per_spatial_unit,
    has_repeated_animals = repeated,
    AnimalID_source = paste(sort(unique(dat$AnimalID_source)), collapse = ";"),
    membership_version = contract$membership_version,
    identity_contract_version = contract$contract_version,
    identity_contract_sha256 = contract$identity_contract_sha256,
    frozen_state_sha256 = contract$frozen_state_sha256
  )
  make_validation <- function(...) {
    values <- c(validation_base, list(...))
    row <- as.data.frame(values, stringsAsFactors = FALSE)
    wgcna_group_bind_schema(list(row), wgcna_group_model_validation_schema())
  }

  if (nrow(dat) < 4L || length(levels(dat$StressGroup)) < 2L) {
    return(list(
      primary = wgcna_group_primary_schema(),
      validation = make_validation(
        attempt_status = "failed_preflight",
        failure_reason = "too_few_samples_or_groups",
        actual_formula = NA_character_, model_type = NA_character_,
        fixed_effect_rank = NA_integer_, fixed_effect_columns = NA_integer_,
        rank_deficient_model = NA, singular_model = NA,
        convergence_status = "not_fitted", convergence_warnings = NA_character_,
        optimizer_messages = NA_character_, optimizer_code = NA_character_,
        model_warnings = NA_character_, residual_df = NA_real_,
        n_contrasts_requested = 3L, n_contrasts_eligible = 0L,
        n_contrasts_emitted = 0L,
        excluded_contrasts = "all:too_few_samples_or_groups",
        emmeans_status = "not_run",
        animal_random_effect_used = FALSE
      )
    ))
  }
  if (effect_scope == "stress_by_spatial_interaction") {
    cells <- table(dat$StressGroup, dat$SpatialUnit)
    if (nrow(cells) < 2L || ncol(cells) < 2L || any(cells == 0L)) {
      return(list(
        primary = wgcna_group_primary_schema(),
        validation = make_validation(
          attempt_status = "failed_preflight",
          failure_reason = "incomplete_stress_by_spatial_design",
          actual_formula = NA_character_, model_type = NA_character_,
          fixed_effect_rank = NA_integer_, fixed_effect_columns = NA_integer_,
          rank_deficient_model = NA, singular_model = NA,
          convergence_status = "not_fitted",
          convergence_warnings = NA_character_,
          optimizer_messages = NA_character_, optimizer_code = NA_character_,
          model_warnings = NA_character_, residual_df = NA_real_,
          n_contrasts_requested = 3L * nlevels(dat$SpatialUnit),
          n_contrasts_eligible = 0L, n_contrasts_emitted = 0L,
          excluded_contrasts = "all:incomplete_stress_by_spatial_design",
          emmeans_status = "not_run",
          animal_random_effect_used = FALSE
        )
      ))
    }
  }

  fixed_matrix <- tryCatch(
    stats::model.matrix(fixed_formula, data = dat),
    error = identity
  )
  if (inherits(fixed_matrix, "error")) {
    return(list(
      primary = wgcna_group_primary_schema(),
      validation = make_validation(
        attempt_status = "failed_preflight",
        failure_reason = paste0("fixed_effect_matrix_failed:", conditionMessage(fixed_matrix)),
        actual_formula = NA_character_, model_type = NA_character_,
        fixed_effect_rank = NA_integer_, fixed_effect_columns = NA_integer_,
        rank_deficient_model = TRUE, singular_model = NA,
        convergence_status = "not_fitted", convergence_warnings = NA_character_,
        optimizer_messages = NA_character_, optimizer_code = NA_character_,
        model_warnings = NA_character_, residual_df = NA_real_,
        n_contrasts_requested = 3L, n_contrasts_eligible = 0L,
        n_contrasts_emitted = 0L,
        excluded_contrasts = "all:fixed_effect_matrix_failed",
        emmeans_status = "not_run",
        animal_random_effect_used = FALSE
      )
    ))
  }
  fixed_rank <- qr(fixed_matrix)$rank
  rank_deficient <- fixed_rank < ncol(fixed_matrix)
  actual_formula <- if (repeated) {
    stats::as.formula(
      paste(deparse(fixed_formula), "+ (1 | AnimalID)")
    )
  } else {
    fixed_formula
  }
  fit_capture <- wgcna_group_capture(
    if (repeated) {
      lmerTest::lmer(actual_formula, data = dat, REML = FALSE)
    } else {
      stats::lm(actual_formula, data = dat)
    }
  )
  if (inherits(fit_capture$value, "error")) {
    return(list(
      primary = wgcna_group_primary_schema(),
      validation = make_validation(
        attempt_status = "failed_fit",
        failure_reason = conditionMessage(fit_capture$value),
        actual_formula = paste(deparse(actual_formula), collapse = ""),
        model_type = model_type, fixed_effect_rank = fixed_rank,
        fixed_effect_columns = ncol(fixed_matrix),
        rank_deficient_model = rank_deficient, singular_model = NA,
        convergence_status = "fit_failed",
        convergence_warnings = paste(fit_capture$warnings, collapse = ";"),
        optimizer_messages = NA_character_, optimizer_code = NA_character_,
        model_warnings = paste(fit_capture$warnings, collapse = ";"),
        residual_df = NA_real_, n_contrasts_requested = 3L,
        n_contrasts_eligible = 0L, n_contrasts_emitted = 0L,
        excluded_contrasts = "all:model_fit_failed",
        emmeans_status = "not_run",
        animal_random_effect_used = repeated
      )
    ))
  }
  fit <- fit_capture$value
  singular <- if (repeated) {
    lme4::isSingular(
      fit, tol = wgcna_group_lme4_singularity_tolerance()
    )
  } else {
    FALSE
  }
  convergence_messages <- if (repeated) {
    as.character(fit@optinfo$conv$lme4$messages %||% character())
  } else {
    character()
  }
  optimizer_code <- if (repeated) {
    paste(as.character(fit@optinfo$conv$opt %||% 0), collapse = ";")
  } else {
    "0"
  }
  optimizer_messages <- if (repeated) {
    paste(as.character(fit@optinfo$warnings %||% character()), collapse = ";")
  } else {
    NA_character_
  }
  convergence_failure_messages <- convergence_messages[
    !grepl(
      "^boundary \\(singular\\) fit",
      convergence_messages,
      ignore.case = TRUE
    )
  ]
  convergence_ok <- !length(convergence_failure_messages) &&
    all(suppressWarnings(as.numeric(strsplit(optimizer_code, ";", fixed = TRUE)[[1]])) == 0)
  model_valid <- !rank_deficient && !singular && convergence_ok
  failure <- c(
    if (rank_deficient) "rank_deficient_fixed_effects",
    if (singular) "singular_random_effect",
    if (!convergence_ok) "failed_convergence"
  )
  residual_df <- suppressWarnings(as.numeric(stats::df.residual(fit)))
  observed_levels <- if (repeated) {
    levels(dat$StressGroup)
  } else {
    levels(model.frame(fit)$StressGroup)
  }
  methods <- wgcna_group_predeclared_contrasts(observed_levels)
  by_spatial <- effect_scope == "stress_by_spatial_interaction"
  requested_contrasts <- if (by_spatial) {
    3L * nlevels(dat$SpatialUnit)
  } else {
    3L
  }
  if (!model_valid) {
    return(list(
      primary = wgcna_group_primary_schema(),
      validation = make_validation(
        attempt_status = "failed_model_validity",
        failure_reason = paste(failure, collapse = ";"),
        actual_formula = paste(deparse(stats::formula(fit)), collapse = ""),
        model_type = model_type, fixed_effect_rank = fixed_rank,
        fixed_effect_columns = ncol(fixed_matrix),
        rank_deficient_model = rank_deficient, singular_model = singular,
        convergence_status = if (convergence_ok) "converged" else "failed",
        convergence_warnings = paste(convergence_messages, collapse = ";"),
        optimizer_messages = optimizer_messages, optimizer_code = optimizer_code,
        model_warnings = paste(fit_capture$warnings, collapse = ";"),
        residual_df = residual_df,
        n_contrasts_requested = requested_contrasts,
        n_contrasts_eligible = 0L, n_contrasts_emitted = 0L,
        excluded_contrasts = paste0("all:", paste(failure, collapse = ";")),
        emmeans_status = "not_run_invalid_model",
        animal_random_effect_used = repeated
      )
    ))
  }

  emm_capture <- wgcna_group_capture(
    if (by_spatial) {
      emmeans::emmeans(fit, ~ StressGroup | SpatialUnit)
    } else {
      emmeans::emmeans(fit, ~ StressGroup)
    }
  )
  if (inherits(emm_capture$value, "error") || !length(methods)) {
    reason <- if (inherits(emm_capture$value, "error")) {
      conditionMessage(emm_capture$value)
    } else {
      "no_predeclared_contrast_has_both_groups"
    }
    return(list(
      primary = wgcna_group_primary_schema(),
      validation = make_validation(
        attempt_status = "failed_emmeans",
        failure_reason = reason,
        actual_formula = paste(deparse(stats::formula(fit)), collapse = ""),
        model_type = model_type, fixed_effect_rank = fixed_rank,
        fixed_effect_columns = ncol(fixed_matrix),
        rank_deficient_model = FALSE, singular_model = FALSE,
        convergence_status = "converged",
        convergence_warnings = paste(convergence_messages, collapse = ";"),
        optimizer_messages = optimizer_messages, optimizer_code = optimizer_code,
        model_warnings = paste(c(fit_capture$warnings, emm_capture$warnings), collapse = ";"),
        residual_df = residual_df,
        n_contrasts_requested = requested_contrasts,
        n_contrasts_eligible = 0L, n_contrasts_emitted = 0L,
        excluded_contrasts = paste0("all:", reason),
        emmeans_status = "failed",
        animal_random_effect_used = repeated
      )
    ))
  }
  contrast_capture <- wgcna_group_capture(
    as.data.frame(
      wgcna_group_emmeans_contrast_summary(
        emm_capture$value, methods, ci_level = 0.95
      )
    )
  )
  if (inherits(contrast_capture$value, "error")) {
    return(list(
      primary = wgcna_group_primary_schema(),
      validation = make_validation(
        attempt_status = "failed_emmeans",
        failure_reason = conditionMessage(contrast_capture$value),
        actual_formula = paste(deparse(stats::formula(fit)), collapse = ""),
        model_type = model_type, fixed_effect_rank = fixed_rank,
        fixed_effect_columns = ncol(fixed_matrix),
        rank_deficient_model = FALSE, singular_model = FALSE,
        convergence_status = "converged",
        convergence_warnings = paste(convergence_messages, collapse = ";"),
        optimizer_messages = optimizer_messages, optimizer_code = optimizer_code,
        model_warnings = paste(
          c(fit_capture$warnings, emm_capture$warnings, contrast_capture$warnings),
          collapse = ";"
        ),
        residual_df = residual_df,
        n_contrasts_requested = requested_contrasts,
        n_contrasts_eligible = 0L, n_contrasts_emitted = 0L,
        excluded_contrasts = "all:contrast_failed",
        emmeans_status = "failed",
        animal_random_effect_used = repeated
      )
    ))
  }
  contrast_df <- contrast_capture$value
  stat_col <- intersect(c("t.ratio", "z.ratio"), names(contrast_df))
  stat_col <- if (length(stat_col)) stat_col[[1]] else NA_character_
  spatial_values <- if (by_spatial) {
    as.character(contrast_df$SpatialUnit)
  } else {
    rep(ifelse(
      is.na(spatial_unit),
      if (effect_scope == "spatial_adjusted_global") {
        "global_spatial_adjusted"
      } else {
        "all_spatial_units"
      },
      spatial_unit
    ), nrow(contrast_df))
  }
  pairs <- list(
    "RES - CON" = c("RES", "CON"),
    "SUS - CON" = c("SUS", "CON"),
    "SUS - RES" = c("SUS", "RES")
  )
  primary <- list()
  excluded <- character()
  eligible <- 0L
  for (i in seq_len(nrow(contrast_df))) {
    label <- as.character(contrast_df$contrast[[i]])
    if (!label %in% names(pairs)) next
    unit <- spatial_values[[i]]
    support_dat <- if (by_spatial) {
      dat[as.character(dat$SpatialUnit) == unit, , drop = FALSE]
    } else {
      dat
    }
    pair <- pairs[[label]]
    support <- tapply(
      support_dat$AnimalID[support_dat$StressGroup %in% pair],
      support_dat$StressGroup[support_dat$StressGroup %in% pair],
      function(x) length(unique(stats::na.omit(x)))
    )
    support <- support[pair]
    min_animals <- if (length(support) == 2L && !anyNA(support)) {
      min(as.integer(support))
    } else {
      0L
    }
    finite <- all(is.finite(c(
      contrast_df$estimate[[i]], contrast_df$SE[[i]],
      if (!is.na(stat_col)) contrast_df[[stat_col]][[i]] else NA_real_,
      contrast_df$p.value[[i]]
    )))
    if (min_animals < 3L) {
      excluded <- c(excluded, paste0(unit, "|", label, ":insufficient_unique_animals"))
      next
    }
    eligible <- eligible + 1L
    if (!finite) {
      excluded <- c(excluded, paste0(unit, "|", label, ":nonfinite_or_nonestimable"))
      next
    }
    primary[[length(primary) + 1L]] <- data.frame(
      dataset = dataset, level = level, endpoint_id = endpoint_id,
      endpoint_label = endpoint_id,
      module_id = if (level == "module") endpoint_id else NA_character_,
      supermodule_id = if (level == "supermodule") endpoint_id else NA_character_,
      module_label = if (level == "module") endpoint_id else NA_character_,
      supermodule_label = if (level == "supermodule") endpoint_id else NA_character_,
      spatial_unit = unit, effect_scope = effect_scope,
      SpatialUnitType = if (dataset == "neuron_neuropil") "region_layer" else "region",
      model_type = model_type, has_repeated_animals = repeated,
      n_animals = counts$n_animals, n_animals_total = counts$n_animals,
      n_animals_per_group = counts$n_animals_per_group,
      min_animals_per_group = counts$min_animals_per_group,
      min_unique_animals_compared_group = min_animals,
      n_samples_total = counts$n_samples,
      n_samples_per_group = counts$n_samples_per_group,
      n_samples_per_spatial_unit = counts$n_samples_per_spatial_unit,
      n_animals_per_spatial_unit = counts$n_animals_per_spatial_unit,
      animal_level_status = "biological_replicate_valid",
      pseudoreplication_guard = if (repeated) {
        "animal_random_intercept_used"
      } else {
        "one_observation_per_animal"
      },
      contrast = label, estimate = as.numeric(contrast_df$estimate[[i]]),
      SE = as.numeric(contrast_df$SE[[i]]),
      CI_low = as.numeric(contrast_df$CI_low[[i]]),
      CI_high = as.numeric(contrast_df$CI_high[[i]]),
      CI_method = as.character(contrast_df$CI_method[[i]]),
      CI_level = as.numeric(contrast_df$CI_level[[i]]),
      CI_df_method = as.character(contrast_df$CI_df_method[[i]]),
      statistic = as.numeric(contrast_df[[stat_col]][[i]]),
      p_value = as.numeric(contrast_df$p.value[[i]]),
      evidence_status = NA_character_,
      direction = ifelse(
        contrast_df$estimate[[i]] > 0, "higher",
        ifelse(contrast_df$estimate[[i]] < 0, "lower", "zero")
      ),
      n_samples = counts$n_samples,
      formula_requested = requested,
      formula_used = paste(deparse(stats::formula(fit)), collapse = ""),
      dropped_covariates = paste(
        setdiff(c("Sex", "Batch"), covariates), collapse = ";"
      ),
      model_family = if (repeated) "linear_mixed_model" else "linear_model",
      model_formula = paste(deparse(stats::formula(fit)), collapse = ""),
      primary_model_stable = TRUE, claim_allowed_model = TRUE,
      model_downgrade_reason = "none", fallback_used = FALSE,
      fallback_type = "none", rank_deficient_model = FALSE,
      singular_model = FALSE, emmeans_success = TRUE,
      animal_random_effect_used = repeated,
      biological_replicate_unit = "AnimalID",
      model_warning = paste(
        unique(c(fit_capture$warnings, emm_capture$warnings, contrast_capture$warnings)),
        collapse = ";"
      ),
      fixed_effect_rank = fixed_rank,
      fixed_effect_columns = ncol(fixed_matrix),
      convergence_status = "converged",
      convergence_warnings = paste(convergence_messages, collapse = ";"),
      optimizer_messages = optimizer_messages, optimizer_code = optimizer_code,
      residual_df = if ("df" %in% names(contrast_df)) {
        as.numeric(contrast_df$df[[i]])
      } else {
        residual_df
      },
      AnimalID_source = paste(sort(unique(dat$AnimalID_source)), collapse = ";"),
      animal_id_mapping_status = "resolved",
      membership_version = contract$membership_version,
      identity_contract_version = contract$contract_version,
      identity_contract_sha256 = contract$identity_contract_sha256,
      frozen_state_sha256 = contract$frozen_state_sha256,
      endpoint_construction_method = construction_method,
      endpoint_provenance_status = provenance_status,
      stringsAsFactors = FALSE
    )
  }
  primary <- wgcna_group_bind_schema(primary, wgcna_group_primary_schema())
  attempt_status <- if (nrow(primary)) {
    if (length(excluded)) "success_with_excluded_contrasts" else "success"
  } else {
    "failed_contrasts"
  }
  validation <- make_validation(
    attempt_status = attempt_status,
    failure_reason = if (nrow(primary)) "none" else "no_finite_estimable_contrasts",
    actual_formula = paste(deparse(stats::formula(fit)), collapse = ""),
    model_type = model_type, fixed_effect_rank = fixed_rank,
    fixed_effect_columns = ncol(fixed_matrix),
    rank_deficient_model = FALSE, singular_model = FALSE,
    convergence_status = "converged",
    convergence_warnings = paste(convergence_messages, collapse = ";"),
    optimizer_messages = optimizer_messages, optimizer_code = optimizer_code,
    model_warnings = paste(
      unique(c(fit_capture$warnings, emm_capture$warnings, contrast_capture$warnings)),
      collapse = ";"
    ),
    residual_df = residual_df,
    n_contrasts_requested = requested_contrasts,
    n_contrasts_eligible = eligible,
    n_contrasts_emitted = nrow(primary),
    excluded_contrasts = if (length(excluded)) {
      paste(excluded, collapse = ";")
    } else {
      "none"
    },
    emmeans_status = "success",
    animal_random_effect_used = repeated
  )
  list(primary = primary, validation = validation)
}

wgcna_group_run_level_phase2_legacy <- function(
    dataset, level, eigengenes, endpoint_map, sample_audit, contract
) {
  dat <- sample_audit
  idx <- match(dat$Sample, eigengenes$Sample)
  if (anyNA(idx) || nrow(eigengenes) != nrow(dat)) {
    stop("Endpoint/sample alignment would lose state samples.", call. = FALSE)
  }
  dat$SpatialUnit <- dat$canonical_spatial_unit
  for (nm in setdiff(names(eigengenes), "Sample")) {
    dat[[nm]] <- eigengenes[[nm]][idx]
  }
  primary <- list()
  validation <- list()
  spatial_units <- sort(unique(dat$SpatialUnit))
  endpoint_map <- endpoint_map[order(endpoint_map$endpoint_id), , drop = FALSE]
  for (i in seq_len(nrow(endpoint_map))) {
    endpoint <- endpoint_map[i, , drop = FALSE]
    for (unit in spatial_units) {
      result <- wgcna_group_fit_attempt(
        dat, dataset, level, endpoint$endpoint_id,
        endpoint$endpoint_col, endpoint$endpoint_construction_method,
        endpoint$endpoint_provenance_status, "within_spatial_unit",
        contract, unit
      )
      primary[[length(primary) + 1L]] <- result$primary
      validation[[length(validation) + 1L]] <- result$validation
    }
    for (scope in c(
      "spatial_adjusted_global", "stress_by_spatial_interaction"
    )) {
      result <- wgcna_group_fit_attempt(
        dat, dataset, level, endpoint$endpoint_id,
        endpoint$endpoint_col, endpoint$endpoint_construction_method,
        endpoint$endpoint_provenance_status, scope, contract
      )
      primary[[length(primary) + 1L]] <- result$primary
      validation[[length(validation) + 1L]] <- result$validation
    }
  }
  list(
    primary = wgcna_group_bind_schema(primary, wgcna_group_primary_schema()),
    validation = wgcna_group_bind_schema(
      validation, wgcna_group_model_validation_schema()
    )
  )
}

wgcna_group_apply_fdr_phase2_legacy <- function(module_rows, supermodule_rows) {
  all <- dplyr::bind_rows(module_rows, supermodule_rows)
  if (!nrow(all)) {
    return(list(module = module_rows, supermodule = supermodule_rows))
  }
  if (any(!all$claim_allowed_model) ||
      any(!is.finite(all$p_value)) ||
      any(all$fallback_used)) {
    stop("Only finite, claim-allowed, non-fallback rows may enter FDR families.",
         call. = FALSE)
  }
  all$FDR_within_dataset_level <- NA_real_
  for (level in unique(all$level)) {
    idx <- which(all$level == level)
    all$FDR_within_dataset_level[idx] <- stats::p.adjust(
      all$p_value[idx], method = "BH"
    )
    all$n_tests_FDR_within_dataset_level[idx] <- length(idx)
    all$FDR_family_within_level_id[idx] <- paste0(
      unique(all$dataset), "::", level, "::all_primary_scopes_contrasts"
    )
  }
  all$FDR_dataset_all_levels <- stats::p.adjust(all$p_value, method = "BH")
  all$FDR_global <- all$FDR_dataset_all_levels
  all$n_tests_FDR_dataset_all_levels <- nrow(all)
  all$FDR_family_dataset_id <- paste0(
    unique(all$dataset), "::all_levels::all_primary_scopes_contrasts"
  )
  all$FDR_method <- "BH"
  all$evidence_status <- dplyr::case_when(
    all$FDR_dataset_all_levels <= 0.05 ~ "robust_FDR",
    all$FDR_dataset_all_levels <= 0.10 ~ "suggestive_FDR10",
    all$FDR_within_dataset_level <= 0.05 ~ "within_level_FDR_only",
    all$p_value < 0.05 ~ "nominal_only",
    TRUE ~ "not_supported"
  )
  key <- paste(
    all$dataset, all$level, all$endpoint_id, all$effect_scope,
    all$spatial_unit, all$contrast, sep = "|"
  )
  if (anyDuplicated(key)) stop("Primary endpoint statistical keys are duplicated.",
                               call. = FALSE)
  if (any(!is.finite(all$FDR_within_dataset_level)) ||
      any(!is.finite(all$FDR_dataset_all_levels)) ||
      any(all$FDR_within_dataset_level < 0 |
          all$FDR_within_dataset_level > 1) ||
      any(all$FDR_dataset_all_levels < 0 |
          all$FDR_dataset_all_levels > 1) ||
      !identical(all$FDR_global, all$FDR_dataset_all_levels)) {
    stop("FDR validation failed.", call. = FALSE)
  }
  all <- all[order(
    all$level, all$endpoint_id, all$effect_scope,
    all$spatial_unit, all$contrast
  ), names(wgcna_group_primary_schema()), drop = FALSE]
  list(
    module = all[all$level == "module", , drop = FALSE],
    supermodule = all[all$level == "supermodule", , drop = FALSE]
  )
}

# Phase 2B statistical-production contract. The Phase 2 implementations above
# remain named as legacy helpers solely to make the contract transition
# reviewable; all canonical Stage 05 calls resolve to the definitions below.
wgcna_group_canonical_text_raw <- function(x) {
  x <- enc2utf8(as.character(x))
  if (anyNA(x)) {
    stop("Canonical hash content contains NA.", call. = FALSE)
  }
  charToRaw(paste0(paste(x, collapse = "\n"), "\n"))
}

wgcna_group_raw_sha256 <- function(x) {
  if (!is.raw(x)) {
    stop("SHA-256 payload must be a raw vector.", call. = FALSE)
  }
  path <- tempfile("wgcna_group_text_", fileext = ".txt")
  on.exit(unlink(path, force = TRUE), add = TRUE)
  connection <- file(path, open = "wb")
  tryCatch(
    writeBin(x, connection),
    finally = close(connection)
  )
  file_hash_sha256(path)
}

wgcna_group_text_sha256 <- function(x) {
  wgcna_group_raw_sha256(wgcna_group_canonical_text_raw(x))
}

wgcna_group_model_row_sha256_v4 <- function(x) {
  x <- enc2utf8(as.character(x))
  if (anyNA(x)) {
    stop("Model-row hash content contains NA.", call. = FALSE)
  }
  x <- gsub("\r\n|\r", "\n", x, perl = TRUE)
  payload <- paste0(paste(x, collapse = "\n"), "\n")
  payload <- gsub("\n", "\r\n", payload, fixed = TRUE)
  wgcna_group_raw_sha256(charToRaw(payload))
}

wgcna_group_aggregate_hash_fields <- function(
    dataset, level, endpoint_id, AnimalID, StressGroup, SpatialUnit,
    ordered_hemispheres, ordered_contributing_sample_ids,
    ordered_hemisphere_source_values, ordered_hemisphere_means,
    bilateral_value, aggregation_method
) {
  c(
    paste0("dataset=", dataset),
    paste0("level=", level),
    paste0("endpoint_id=", endpoint_id),
    paste0("AnimalID=", AnimalID),
    paste0("StressGroup=", StressGroup),
    paste0("SpatialUnit=", SpatialUnit),
    paste0("ordered_hemispheres=", ordered_hemispheres),
    paste0(
      "ordered_contributing_sample_ids=",
      ordered_contributing_sample_ids
    ),
    paste0(
      "ordered_hemisphere_source_values=",
      ordered_hemisphere_source_values
    ),
    paste0("hemisphere_means=", ordered_hemisphere_means),
    paste0("bilateral_value=", sprintf("%.17g", bilateral_value)),
    paste0("aggregation_method=", aggregation_method)
  )
}

wgcna_group_aggregate_row_sha256 <- function(row) {
  if (nrow(row) != 1L) {
    stop("Aggregate hash calculation requires exactly one row.", call. = FALSE)
  }
  wgcna_group_text_sha256(wgcna_group_aggregate_hash_fields(
    dataset = row$dataset[[1]],
    level = row$level[[1]],
    endpoint_id = row$endpoint_id[[1]],
    AnimalID = row$AnimalID[[1]],
    StressGroup = row$StressGroup[[1]],
    SpatialUnit = row$SpatialUnit[[1]],
    ordered_hemispheres = row$ordered_hemispheres[[1]],
    ordered_contributing_sample_ids =
      row$ordered_contributing_sample_ids[[1]],
    ordered_hemisphere_source_values =
      row$ordered_hemisphere_source_values[[1]],
    ordered_hemisphere_means = row$ordered_hemisphere_means[[1]],
    bilateral_value = row$eigengene[[1]],
    aggregation_method = row$aggregation_method[[1]]
  ))
}

wgcna_group_aggregate_endpoint <- function(
    dat, dataset, level, endpoint_id, endpoint_col
) {
  dat$eigengene <- suppressWarnings(as.numeric(dat[[endpoint_col]]))
  if (any(!is.na(dat$Hemisphere) & !dat$Hemisphere %in% c("L", "R")) ||
      any(is.na(dat$Hemisphere) & is.finite(dat$eigengene))) {
    stop("Animal-spatial endpoint aggregation has unresolved hemisphere.",
         call. = FALSE)
  }
  dat <- dat[is.finite(dat$eigengene) &
    dat$StressGroup %in% c("CON", "RES", "SUS"), , drop = FALSE]
  if (!nrow(dat) || anyNA(dat$Sample) || any(!nzchar(dat$Sample))) {
    stop("Animal-spatial endpoint aggregation has missing source samples.",
         call. = FALSE)
  }
  invariant <- dat |>
    dplyr::group_by(.data$AnimalID, .data$canonical_spatial_unit) |>
    dplyr::summarise(
      n_stress_groups = dplyr::n_distinct(.data$StressGroup),
      n_spatial_types = dplyr::n_distinct(.data$SpatialUnitType),
      n_regions = dplyr::n_distinct(.data$Region),
      n_layers = dplyr::n_distinct(.data$Layer),
      .groups = "drop"
    )
  if (any(invariant$n_stress_groups != 1L) ||
      any(invariant$n_spatial_types != 1L) ||
      any(invariant$n_regions != 1L) ||
      any(invariant$n_layers != 1L)) {
    stop("Animal-spatial endpoint aggregation has inconsistent metadata.",
         call. = FALSE)
  }

  hemisphere_keys <- c(
    "AnimalID", "StressGroup", "canonical_spatial_unit", "Hemisphere"
  )
  split_key <- do.call(
    interaction,
    c(dat[hemisphere_keys], list(
      drop = TRUE, lex.order = TRUE, sep = "\035"
    ))
  )
  hemisphere_rows <- lapply(split(dat, split_key), function(cell) {
    cell <- cell[order(cell$Sample, cell$eigengene), , drop = FALSE]
    duplicated_sample <- duplicated(cell$Sample) |
      duplicated(cell$Sample, fromLast = TRUE)
    if (any(duplicated_sample)) {
      ambiguous <- vapply(
        split(cell$eigengene[duplicated_sample],
              cell$Sample[duplicated_sample]),
        function(x) length(unique(sprintf("%.17g", x))) > 1L,
        logical(1)
      )
      if (any(ambiguous)) {
        stop(
          "Duplicate source-sample ambiguity cannot be deterministically ",
          "resolved for an animal-spatial-hemisphere cell.",
          call. = FALSE
        )
      }
    }
    source_values <- sprintf("%.17g", cell$eigengene)
    data.frame(
      dataset = dataset,
      level = level,
      endpoint_id = endpoint_id,
      AnimalID = cell$AnimalID[[1]],
      StressGroup = cell$StressGroup[[1]],
      SpatialUnit = cell$canonical_spatial_unit[[1]],
      canonical_spatial_unit = cell$canonical_spatial_unit[[1]],
      Hemisphere = cell$Hemisphere[[1]],
      contributing_sample_ids = paste(cell$Sample, collapse = ";"),
      source_sample_count = nrow(cell),
      source_values = paste(source_values, collapse = ";"),
      hemisphere_value = mean(cell$eigengene),
      n_hemisphere_source_rows = nrow(cell),
      duplicate_source_row_status = if (nrow(cell) > 1L) {
        "technical_duplicates_collapsed_by_mean"
      } else {
        "single_source_row"
      },
      aggregation_method = "arithmetic_mean_within_hemisphere",
      AnimalID_source = paste(
        sort(unique(cell$AnimalID_source)), collapse = ";"
      ),
      Region = cell$Region[[1]],
      Layer = cell$Layer[[1]],
      SpatialUnitType = cell$SpatialUnitType[[1]],
      stringsAsFactors = FALSE
    )
  })
  hemisphere_values <- wgcna_group_bind_schema(
    hemisphere_rows, wgcna_group_hemisphere_value_schema()
  ) |>
    dplyr::arrange(
      .data$endpoint_id, .data$AnimalID, .data$SpatialUnit,
      .data$Hemisphere
    )

  final_keys <- c("AnimalID", "StressGroup", "SpatialUnit")
  final_split <- do.call(
    interaction,
    c(hemisphere_values[final_keys], list(
      drop = TRUE, lex.order = TRUE, sep = "\035"
    ))
  )
  aggregated_rows <- lapply(
    split(hemisphere_values, final_split),
    function(cell) {
      cell <- cell[order(cell$Hemisphere), , drop = FALSE]
      if (anyDuplicated(cell$Hemisphere) || !nrow(cell) ||
          nrow(cell) > 2L) {
        stop(
          "Animal-spatial endpoint aggregation failed hemisphere cardinality.",
          call. = FALSE
        )
      }
      final_value <- mean(cell$hemisphere_value)
      aggregation_method <-
        "equal_weight_mean_available_LR_after_within_hemisphere_mean"
      canonical_content <- wgcna_group_aggregate_hash_fields(
        dataset = dataset,
        level = level,
        endpoint_id = endpoint_id,
        AnimalID = cell$AnimalID[[1]],
        StressGroup = cell$StressGroup[[1]],
        SpatialUnit = cell$SpatialUnit[[1]],
        ordered_hemispheres = paste(cell$Hemisphere, collapse = ";"),
        ordered_contributing_sample_ids = paste(
          cell$contributing_sample_ids, collapse = "|"
        ),
        ordered_hemisphere_source_values = paste(
          cell$source_values, collapse = "|"
        ),
        ordered_hemisphere_means = paste(
          sprintf("%.17g", cell$hemisphere_value), collapse = ";"
        ),
        bilateral_value = final_value,
        aggregation_method = aggregation_method
      )
      data.frame(
        AnimalID = cell$AnimalID[[1]],
        StressGroup = cell$StressGroup[[1]],
        canonical_spatial_unit = cell$SpatialUnit[[1]],
        eigengene = final_value,
        n_source_samples = sum(cell$source_sample_count),
        n_hemispheres_observed = nrow(cell),
        hemispheres_observed = paste(cell$Hemisphere, collapse = ";"),
        contributing_samples = paste(
          cell$contributing_sample_ids, collapse = ";"
        ),
        ordered_hemispheres = paste(cell$Hemisphere, collapse = ";"),
        ordered_contributing_sample_ids = paste(
          cell$contributing_sample_ids, collapse = "|"
        ),
        ordered_hemisphere_source_values = paste(
          cell$source_values, collapse = "|"
        ),
        ordered_hemisphere_means = paste(
          sprintf("%.17g", cell$hemisphere_value), collapse = ";"
        ),
        ordered_hemisphere_source_row_counts = paste(
          cell$n_hemisphere_source_rows, collapse = ";"
        ),
        ordered_hemisphere_duplicate_statuses = paste(
          cell$duplicate_source_row_status, collapse = ";"
        ),
        hemisphere_values = paste(
          sprintf("%.17g", cell$hemisphere_value), collapse = ";"
        ),
        AnimalID_source = paste(
          sort(unique(cell$AnimalID_source)), collapse = ";"
        ),
        Region = cell$Region[[1]],
        Layer = cell$Layer[[1]],
        SpatialUnitType = cell$SpatialUnitType[[1]],
        dataset = dataset,
        level = level,
        endpoint_id = endpoint_id,
        SpatialUnit = cell$SpatialUnit[[1]],
        bilateral_pair_status = if (nrow(cell) == 2L) {
          "bilateral_complete"
        } else {
          "single_hemisphere_observed"
        },
        hemisphere_completeness_status = if (nrow(cell) == 2L) {
          "complete_bilateral_pair"
        } else {
          "one_sided_observed_no_imputation"
        },
        aggregation_method = aggregation_method,
        duplicate_source_row_status = if (
          any(cell$n_hemisphere_source_rows > 1L)
        ) {
          "technical_duplicates_collapsed_within_hemisphere"
        } else {
          "no_technical_duplicates"
        },
        aggregated_row_id = paste(
          dataset, level, endpoint_id, cell$AnimalID[[1]],
          cell$SpatialUnit[[1]], sep = "|"
        ),
        aggregated_row_sha256 =
          wgcna_group_text_sha256(canonical_content),
        aggregated_row_hash_contract = "sha256_utf8_lf_v1",
        stringsAsFactors = FALSE
      )
    }
  )
  aggregated <- wgcna_group_bind_schema(
    aggregated_rows, wgcna_group_animal_spatial_value_schema()
  )
  if (anyDuplicated(aggregated$aggregated_row_id) ||
      anyDuplicated(aggregated$aggregated_row_sha256) ||
      any(!is.finite(aggregated$eigengene)) ||
      any(!aggregated$n_hemispheres_observed %in% 1:2)) {
    stop("Animal-spatial endpoint aggregation failed its cardinality contract.",
         call. = FALSE)
  }
  list(
    hemisphere_values = hemisphere_values,
    animal_spatial_values = aggregated |>
      dplyr::arrange(.data$endpoint_id, .data$AnimalID, .data$SpatialUnit)
  )
}

wgcna_group_random_structure_class <- function(formula) {
  bars <- if (requireNamespace("reformulas", quietly = TRUE)) {
    reformulas::findbars(formula)
  } else {
    suppressWarnings(lme4::findbars(formula))
  }
  if (!length(bars)) return("none")
  if (length(bars) == 1L) {
    lhs <- paste(deparse(bars[[1]][[2]]), collapse = "")
    grouping <- paste(deparse(bars[[1]][[3]]), collapse = "")
    if (identical(lhs, "1") && identical(grouping, "AnimalID")) {
      return("single_animal_intercept")
    }
  }
  "complex_random_effects"
}

wgcna_group_classify_model_diagnostics <- function(
    model_type, fixed_effect_full_rank, optimizer_success,
    nonboundary_convergence_success, finite_model_quantities,
    contrast_estimable = TRUE, finite_contrast = TRUE,
    random_effect_structure = "none", random_intercept_variance = NA_real_,
    residual_variance = NA_real_, is_singular_lme4 = FALSE,
    sample_contract_valid = TRUE,
    boundary_variance_ratio_tolerance =
      wgcna_group_boundary_variance_ratio_tolerance(),
    lme4_singularity_tolerance =
      wgcna_group_lme4_singularity_tolerance()
) {
  variance_ratio <- if (
    is.finite(random_intercept_variance) &&
      is.finite(residual_variance) && residual_variance != 0
  ) {
    random_intercept_variance / residual_variance
  } else {
    NA_real_
  }
  icc <- if (
    is.finite(random_intercept_variance) &&
      is.finite(residual_variance) &&
      (random_intercept_variance + residual_variance) != 0
  ) {
    random_intercept_variance /
      (random_intercept_variance + residual_variance)
  } else {
    NA_real_
  }
  variance_diagnostics_valid <- is.finite(random_intercept_variance) &&
    random_intercept_variance >= 0 &&
    is.finite(residual_variance) && residual_variance > 0 &&
    is.finite(variance_ratio) && variance_ratio >= 0 &&
    is.finite(icc) && icc >= 0 && icc <= 1
  boundary_by_ratio <- if (variance_diagnostics_valid) {
    variance_ratio <= boundary_variance_ratio_tolerance
  } else {
    NA
  }
  diagnostic_status <- if (identical(model_type, "lm")) {
    "not_applicable_lm"
  } else if (!isTRUE(variance_diagnostics_valid) ||
             !is.logical(is_singular_lme4) ||
             length(is_singular_lme4) != 1L ||
             is.na(is_singular_lme4)) {
    NA_character_
  } else if (isTRUE(boundary_by_ratio) && isTRUE(is_singular_lme4)) {
    "concordant_boundary"
  } else if (!isTRUE(boundary_by_ratio) && !isTRUE(is_singular_lme4)) {
    "concordant_nonboundary"
  } else if (isTRUE(boundary_by_ratio)) {
    "variance_ratio_boundary_only"
  } else {
    "lme4_singular_only"
  }
  diagnostic_disagreement <- diagnostic_status %in% c(
    "variance_ratio_boundary_only", "lme4_singular_only"
  )
  invalid_reasons <- c(
    if (!isTRUE(sample_contract_valid)) "failed_sample_contract",
    if (!isTRUE(fixed_effect_full_rank)) "rank_deficient_fixed_effects",
    if (!isTRUE(optimizer_success)) "optimizer_failure",
    if (!isTRUE(nonboundary_convergence_success)) {
      "nonboundary_convergence_failure"
    },
    if (!isTRUE(finite_model_quantities)) "nonfinite_model_quantities",
    if (!isTRUE(contrast_estimable)) "failed_or_nonestimable_contrast",
    if (!isTRUE(finite_contrast)) "nonfinite_contrast",
    if (identical(random_effect_structure, "complex_random_effects") &&
        isTRUE(is_singular_lme4)) {
      "complex_singular_random_effects"
    }
  )
  if (length(invalid_reasons)) {
    return(list(
      model_valid_for_inference = FALSE,
      model_stability_status = "invalid",
      primary_model_stable = FALSE,
      claim_allowed_model = FALSE,
      singularity_class = if (isTRUE(is_singular_lme4)) {
        "invalid_complex_or_unresolved_singularity"
      } else {
        "not_boundary_classified"
      },
      boundary_warning = NA_character_,
      failure_reason = paste(invalid_reasons, collapse = ";"),
      variance_ratio = variance_ratio,
      ICC = icc,
      is_singular_lme4 = is_singular_lme4,
      boundary_by_variance_ratio = boundary_by_ratio,
      boundary_variance_ratio_tolerance =
        boundary_variance_ratio_tolerance,
      lme4_singularity_tolerance = lme4_singularity_tolerance,
      singularity_diagnostic_status = diagnostic_status,
      diagnostic_review_required = isTRUE(diagnostic_disagreement)
    ))
  }
  if (identical(model_type, "lm")) {
    return(list(
      model_valid_for_inference = TRUE,
      model_stability_status = "stable_animal_level_lm",
      primary_model_stable = TRUE,
      claim_allowed_model = TRUE,
      singularity_class = "not_applicable_lm",
      boundary_warning = NA_character_,
      failure_reason = "none",
      variance_ratio = NA_real_,
      ICC = NA_real_,
      is_singular_lme4 = FALSE,
      boundary_by_variance_ratio = NA,
      boundary_variance_ratio_tolerance =
        boundary_variance_ratio_tolerance,
      lme4_singularity_tolerance = lme4_singularity_tolerance,
      singularity_diagnostic_status = "not_applicable_lm",
      diagnostic_review_required = FALSE
    ))
  }
  if (!identical(random_effect_structure, "single_animal_intercept") ||
      !isTRUE(variance_diagnostics_valid)) {
    return(list(
      model_valid_for_inference = FALSE,
      model_stability_status = "invalid",
      primary_model_stable = FALSE,
      claim_allowed_model = FALSE,
      singularity_class = "invalid_random_effect_structure_or_variance",
      boundary_warning = NA_character_,
      failure_reason = "invalid_random_effect_structure_or_variance",
      variance_ratio = variance_ratio,
      ICC = icc,
      is_singular_lme4 = is_singular_lme4,
      boundary_by_variance_ratio = boundary_by_ratio,
      boundary_variance_ratio_tolerance =
        boundary_variance_ratio_tolerance,
      lme4_singularity_tolerance = lme4_singularity_tolerance,
      singularity_diagnostic_status = diagnostic_status,
      diagnostic_review_required = isTRUE(diagnostic_disagreement)
    ))
  }
  if (isTRUE(boundary_by_ratio)) {
    warning_text <- paste0(
      "AnimalID random-intercept variance ratio is at the prespecified ",
      "boundary (variance ratio <= ",
      format(boundary_variance_ratio_tolerance, scientific = TRUE),
      "); fixed-effect inference is retained but later animal-level ",
      "readiness evidence is required."
    )
    if (isTRUE(diagnostic_disagreement)) {
      warning_text <- paste0(
        warning_text, " The variance-ratio classification and lme4::",
        "isSingular() diagnostic disagree; Phase 3 diagnostic review is ",
        "required."
      )
    }
    return(list(
      model_valid_for_inference = TRUE,
      model_stability_status = "boundary_random_intercept_zero",
      primary_model_stable = FALSE,
      claim_allowed_model = TRUE,
      singularity_class = "boundary_single_animal_intercept",
      boundary_warning = warning_text,
      failure_reason = "none",
      variance_ratio = variance_ratio,
      ICC = icc,
      is_singular_lme4 = is_singular_lme4,
      boundary_by_variance_ratio = TRUE,
      boundary_variance_ratio_tolerance =
        boundary_variance_ratio_tolerance,
      lme4_singularity_tolerance = lme4_singularity_tolerance,
      singularity_diagnostic_status = diagnostic_status,
      diagnostic_review_required = isTRUE(diagnostic_disagreement)
    ))
  }
  warning_text <- if (isTRUE(diagnostic_disagreement)) {
    paste0(
      "The variance-ratio classification is nonboundary but lme4::",
      "isSingular() is TRUE; inference is retained and Phase 3 diagnostic ",
      "review is required."
    )
  } else {
    NA_character_
  }
  list(
    model_valid_for_inference = TRUE,
    model_stability_status = "stable_mixed_model",
    primary_model_stable = !isTRUE(diagnostic_disagreement),
    claim_allowed_model = TRUE,
    singularity_class = "stable_single_animal_intercept",
    boundary_warning = warning_text,
    failure_reason = "none",
    variance_ratio = variance_ratio,
    ICC = icc,
    is_singular_lme4 = is_singular_lme4,
    boundary_by_variance_ratio = FALSE,
    boundary_variance_ratio_tolerance =
      boundary_variance_ratio_tolerance,
    lme4_singularity_tolerance = lme4_singularity_tolerance,
    singularity_diagnostic_status = diagnostic_status,
    diagnostic_review_required = isTRUE(diagnostic_disagreement)
  )
}

wgcna_group_fit_diagnostics <- function(
    fit, model_type, fixed_rank, fixed_columns, fit_warnings = character(),
    contrast_estimable = TRUE, finite_contrast = TRUE
) {
  repeated <- identical(model_type, "lmerTest_lmer")
  rank_full <- identical(as.integer(fixed_rank), as.integer(fixed_columns))
  if (!repeated) {
    residual_variance <- suppressWarnings(as.numeric(stats::sigma(fit)^2))
    classification <- wgcna_group_classify_model_diagnostics(
      model_type = "lm",
      fixed_effect_full_rank = rank_full,
      optimizer_success = TRUE,
      nonboundary_convergence_success = TRUE,
      finite_model_quantities = all(is.finite(stats::coef(fit))) &&
        is.finite(residual_variance) && residual_variance > 0,
      contrast_estimable = contrast_estimable,
      finite_contrast = finite_contrast,
      residual_variance = residual_variance
    )
    return(c(classification, list(
      random_effect_structure = "none",
      random_intercept_variance = NA_real_,
      residual_variance = residual_variance,
      singular_model = FALSE,
      convergence_status = "converged",
      convergence_warnings = NA_character_,
      optimizer_messages = NA_character_,
      optimizer_code = "0"
    )))
  }
  convergence_messages <- as.character(
    fit@optinfo$conv$lme4$messages %||% character()
  )
  nonboundary <- convergence_messages[
    !grepl("^boundary \\(singular\\) fit", convergence_messages,
           ignore.case = TRUE)
  ]
  optimizer_values <- suppressWarnings(as.numeric(
    unlist(fit@optinfo$conv$opt %||% 0)
  ))
  optimizer_success <- length(optimizer_values) > 0L &&
    all(is.finite(optimizer_values)) && all(optimizer_values == 0)
  optimizer_code <- paste(optimizer_values, collapse = ";")
  optimizer_messages <- paste(
    as.character(fit@optinfo$warnings %||% character()), collapse = ";"
  )
  varcorr <- as.data.frame(lme4::VarCorr(fit))
  random_row <- varcorr[
    varcorr$grp == "AnimalID" &
      varcorr$var1 == "(Intercept)" & is.na(varcorr$var2),
    ,
    drop = FALSE
  ]
  random_variance <- if (nrow(random_row) == 1L) {
    suppressWarnings(as.numeric(random_row$vcov[[1]]))
  } else {
    NA_real_
  }
  residual_variance <- suppressWarnings(as.numeric(stats::sigma(fit)^2))
  random_structure <- wgcna_group_random_structure_class(
    stats::formula(fit)
  )
  singular <- lme4::isSingular(
    fit, tol = wgcna_group_lme4_singularity_tolerance()
  )
  finite_model <- all(is.finite(lme4::fixef(fit))) &&
    is.finite(residual_variance) && residual_variance > 0 &&
    is.finite(random_variance) && random_variance >= 0
  classification <- wgcna_group_classify_model_diagnostics(
    model_type = model_type,
    fixed_effect_full_rank = rank_full,
    optimizer_success = optimizer_success,
    nonboundary_convergence_success = !length(nonboundary),
    finite_model_quantities = finite_model,
    contrast_estimable = contrast_estimable,
    finite_contrast = finite_contrast,
    random_effect_structure = random_structure,
    random_intercept_variance = random_variance,
    residual_variance = residual_variance,
    is_singular_lme4 = singular
  )
  c(classification, list(
    random_effect_structure = random_structure,
    random_intercept_variance = random_variance,
    residual_variance = residual_variance,
    singular_model = singular,
    convergence_status = if (!length(nonboundary) && optimizer_success) {
      "converged"
    } else {
      "failed"
    },
    convergence_warnings = paste(convergence_messages, collapse = ";"),
    optimizer_messages = optimizer_messages,
    optimizer_code = optimizer_code,
    model_warnings = paste(unique(fit_warnings), collapse = ";")
  ))
}

wgcna_group_hypothesis_level <- function(level) {
  if (identical(level, "module")) "module" else "higher_order_multimodule"
}

wgcna_group_composite_interaction_diagnostics <- function(
    reduced_diagnostics, full_diagnostics,
    identical_rows_verified = TRUE,
    likelihood_ratio_valid = TRUE,
    failure_reason = NULL
) {
  reduced_valid <- isTRUE(
    reduced_diagnostics$model_valid_for_inference
  )
  full_valid <- isTRUE(full_diagnostics$model_valid_for_inference)
  composite_valid <- reduced_valid && full_valid &&
    isTRUE(identical_rows_verified) && isTRUE(likelihood_ratio_valid)
  boundary <- composite_valid && any(c(
    reduced_diagnostics$model_stability_status,
    full_diagnostics$model_stability_status
  ) == "boundary_random_intercept_zero")
  review_required <- composite_valid && (
    isTRUE(reduced_diagnostics$diagnostic_review_required) ||
      isTRUE(full_diagnostics$diagnostic_review_required)
  )
  warnings <- unique(stats::na.omit(c(
    reduced_diagnostics$boundary_warning,
    full_diagnostics$boundary_warning
  )))
  warnings <- warnings[nzchar(warnings)]
  out <- full_diagnostics
  out$model_valid_for_inference <- composite_valid
  out$model_stability_status <- if (!composite_valid) {
    "invalid"
  } else if (boundary) {
    "boundary_random_intercept_zero"
  } else {
    "stable_mixed_model"
  }
  out$primary_model_stable <- composite_valid &&
    !boundary && !review_required
  out$claim_allowed_model <- composite_valid
  out$diagnostic_review_required <- review_required
  out$singularity_class <- if (!composite_valid) {
    "invalid_composite_reduced_full"
  } else if (boundary) {
    "composite_reduced_full_boundary"
  } else if (review_required) {
    "composite_reduced_full_diagnostic_review"
  } else {
    "composite_reduced_full_stable"
  }
  out$boundary_warning <- if (length(warnings)) {
    paste(warnings, collapse = ";")
  } else {
    NA_character_
  }
  out$failure_reason <- if (composite_valid) {
    "none"
  } else if (!is.null(failure_reason) && nzchar(failure_reason)) {
    failure_reason
  } else {
    paste(c(
      if (!reduced_valid) {
        paste0("reduced:", reduced_diagnostics$failure_reason)
      },
      if (!full_valid) {
        paste0("full:", full_diagnostics$failure_reason)
      },
      if (!isTRUE(identical_rows_verified)) {
        "reduced_full_rows_not_identical"
      },
      if (!isTRUE(likelihood_ratio_valid)) {
        "invalid_likelihood_ratio_result"
      }
    ), collapse = ";")
  }
  out$random_intercept_variance <- NA_real_
  out$residual_variance <- NA_real_
  out$variance_ratio <- NA_real_
  out$ICC <- NA_real_
  out$is_singular_lme4 <- NA
  out$boundary_by_variance_ratio <- NA
  out
}

wgcna_group_invalid_composite_diagnostics <- function(
    diagnostics, failure_reason = diagnostics$failure_reason
) {
  invalid <- diagnostics
  invalid$model_valid_for_inference <- FALSE
  invalid$model_stability_status <- "invalid"
  invalid$primary_model_stable <- FALSE
  invalid$claim_allowed_model <- FALSE
  invalid$diagnostic_review_required <- FALSE
  invalid$singularity_class <- "invalid_composite_reduced_full"
  invalid$failure_reason <- failure_reason
  invalid$random_intercept_variance <- NA_real_
  invalid$residual_variance <- NA_real_
  invalid$variance_ratio <- NA_real_
  invalid$ICC <- NA_real_
  invalid$is_singular_lme4 <- NA
  invalid$boundary_by_variance_ratio <- NA
  invalid
}

wgcna_group_make_effect_row <- function(
    dat, dataset, level, endpoint, contract, effect_scope, spatial_unit,
    contrast, estimate, SE, statistic, p_value, model_type, fit,
    diagnostics, fixed_rank, fixed_columns, min_animals,
    CI_low = NA_real_, CI_high = NA_real_,
    CI_method = "not_applicable", CI_level = NA_real_,
    CI_df_method = "not_applicable",
    test_type = "named_contrast", df_num = NA_real_, df_den = NA_real_,
    reduced_formula = NA_character_, full_formula = NA_character_,
    identical_rows_verified = NA, reduced_row_hash = NA_character_,
    full_row_hash = NA_character_, likelihood_ratio_statistic = NA_real_,
    likelihood_ratio_df = NA_real_, likelihood_ratio_p_value = NA_real_,
    reduced_model_stability_status = NA_character_,
    full_model_stability_status = NA_character_,
    reduced_diagnostics = NULL, full_diagnostics = NULL,
    reduced_fixed_rank = NA_integer_, reduced_fixed_columns = NA_integer_,
    full_fixed_rank = NA_integer_, full_fixed_columns = NA_integer_
) {
  counts <- wgcna_group_counts(dat)
  hypothesis_level <- wgcna_group_hypothesis_level(level)
  model_diagnostic_scope <- dplyr::case_when(
    identical(test_type, "interaction_omnibus") ~
      "composite_reduced_full",
    identical(test_type, "conditional_interaction_followup") ~
      "full_interaction_model",
    TRUE ~ "single_fitted_model"
  )
  analysis_tier <- dplyr::case_when(
    identical(test_type, "interaction_omnibus") ~
      "secondary_spatial_heterogeneity",
    identical(effect_scope, "spatial_adjusted_global") &&
      identical(contrast, "SUS - RES") ~ "primary_wgcna_global",
    identical(effect_scope, "spatial_adjusted_global") &&
      contrast %in% c("RES - CON", "SUS - CON") ~
      "secondary_contextual_global",
    identical(effect_scope, "within_spatial_unit") ~
      "exploratory_spatial_localization",
    TRUE ~ "exploratory_interaction_followup"
  )
  formula_used <- if (!is.null(fit)) {
    paste(deparse(stats::formula(fit)), collapse = "")
  } else {
    full_formula
  }
  row <- data.frame(
    dataset = dataset, level = level, endpoint_id = endpoint$endpoint_id,
    endpoint_label = endpoint$endpoint_id,
    module_id = if (level == "module") endpoint$endpoint_id else NA_character_,
    supermodule_id = if (level == "supermodule") {
      endpoint$endpoint_id
    } else {
      NA_character_
    },
    module_label = if (level == "module") endpoint$endpoint_id else NA_character_,
    supermodule_label = if (level == "supermodule") {
      endpoint$endpoint_id
    } else {
      NA_character_
    },
    hypothesis_level = hypothesis_level,
    canonical_claim_entity_id = endpoint$endpoint_id,
    claim_entity_role = if (level == "module") {
      "canonical_module"
    } else {
      "higher_order_block"
    },
    support_source_entity_id = endpoint$endpoint_id,
    independent_hypothesis = TRUE,
    enters_primary_fdr = identical(analysis_tier, "primary_wgcna_global"),
    display_allowed = TRUE,
    manuscript_claim_ready = "not_assessed_stage05",
    analysis_tier = analysis_tier,
    primary_analysis_tier = if (identical(
      analysis_tier, "primary_wgcna_global"
    )) {
      "primary_wgcna_global"
    } else {
      NA_character_
    },
    primary_contrast = identical(analysis_tier, "primary_wgcna_global"),
    test_type = test_type,
    spatial_unit = spatial_unit, effect_scope = effect_scope,
    SpatialUnitType = dplyr::first(dat$SpatialUnitType),
    model_type = model_type,
    has_repeated_animals = any(table(dat$AnimalID) > 1L),
    n_animals = counts$n_animals, n_animals_total = counts$n_animals,
    n_animals_per_group = counts$n_animals_per_group,
    min_animals_per_group = counts$min_animals_per_group,
    min_unique_animals_compared_group = min_animals,
    n_samples_total = counts$n_samples,
    n_samples_per_group = counts$n_samples_per_group,
    n_samples_per_spatial_unit = counts$n_samples_per_spatial_unit,
    n_animals_per_spatial_unit = counts$n_animals_per_spatial_unit,
    animal_level_status = if (identical(model_type, "lmerTest_lmer")) {
      "repeated_animal_spatial_unit_mixed_model"
    } else {
      "animal_level"
    },
    pseudoreplication_guard = if (identical(model_type, "lmerTest_lmer")) {
      "bilateral_mean_then_animal_random_intercept"
    } else {
      "one_bilateral_mean_observation_per_animal"
    },
    contrast = contrast, estimate = as.numeric(estimate),
    SE = as.numeric(SE), CI_low = as.numeric(CI_low),
    CI_high = as.numeric(CI_high),
    CI_method = as.character(CI_method), CI_level = as.numeric(CI_level),
    CI_df_method = as.character(CI_df_method),
    statistic = as.numeric(statistic), df_num = as.numeric(df_num),
    df_den = as.numeric(df_den), p_value = as.numeric(p_value),
    evidence_status = NA_character_,
    statistical_support_status = NA_character_,
    direction = if (!is.finite(estimate)) {
      NA_character_
    } else if (estimate > 0) {
      "higher"
    } else if (estimate < 0) {
      "lower"
    } else {
      "zero"
    },
    n_samples = counts$n_samples,
    formula_requested = formula_used, formula_used = formula_used,
    dropped_covariates = "none_not_requested",
    model_family = if (identical(model_type, "lmerTest_lmer")) {
      "linear_mixed_model"
    } else {
      "linear_model"
    },
    model_formula = formula_used,
    model_valid_for_inference = diagnostics$model_valid_for_inference,
    model_diagnostic_scope = model_diagnostic_scope,
    model_stability_status = diagnostics$model_stability_status,
    primary_model_stable = diagnostics$primary_model_stable,
    claim_allowed_model = diagnostics$claim_allowed_model,
    model_downgrade_reason = ifelse(
      diagnostics$model_stability_status == "boundary_random_intercept_zero",
      "boundary_random_intercept_zero",
      "none"
    ),
    fallback_used = FALSE, fallback_type = "none",
    rank_deficient_model = fixed_rank < fixed_columns,
    singular_model = diagnostics$singular_model %||% FALSE,
    emmeans_success = identical(test_type, "named_contrast"),
    animal_random_effect_used = identical(model_type, "lmerTest_lmer"),
    biological_replicate_unit = "AnimalID",
    model_warning = diagnostics$model_warnings %||% NA_character_,
    fixed_effect_rank = fixed_rank, fixed_effect_columns = fixed_columns,
    convergence_status = diagnostics$convergence_status,
    convergence_warnings = diagnostics$convergence_warnings,
    optimizer_messages = diagnostics$optimizer_messages,
    optimizer_code = diagnostics$optimizer_code,
    random_effect_structure = diagnostics$random_effect_structure,
    random_intercept_variance = diagnostics$random_intercept_variance,
    residual_variance = diagnostics$residual_variance,
    variance_ratio = diagnostics$variance_ratio,
    ICC = diagnostics$ICC,
    is_singular_lme4 = if (!is.null(diagnostics$is_singular_lme4)) {
      diagnostics$is_singular_lme4
    } else {
      diagnostics$singular_model %||% FALSE
    },
    boundary_by_variance_ratio =
      diagnostics$boundary_by_variance_ratio %||% NA,
    boundary_variance_ratio_tolerance =
      wgcna_group_boundary_variance_ratio_tolerance(),
    lme4_singularity_tolerance =
      wgcna_group_lme4_singularity_tolerance(),
    singularity_diagnostic_status =
      diagnostics$singularity_diagnostic_status %||% NA_character_,
    diagnostic_review_required =
      diagnostics$diagnostic_review_required %||% FALSE,
    singularity_class = diagnostics$singularity_class,
    boundary_warning = diagnostics$boundary_warning,
    residual_df = if (is.finite(df_den)) df_den else NA_real_,
    AnimalID_source = paste(sort(unique(dat$AnimalID_source)), collapse = ";"),
    animal_id_mapping_status = "resolved",
    reduced_formula = reduced_formula, full_formula = full_formula,
    identical_rows_verified = identical_rows_verified,
    reduced_row_hash = reduced_row_hash, full_row_hash = full_row_hash,
    reduced_fixed_effect_rank = reduced_fixed_rank,
    reduced_fixed_effect_columns = reduced_fixed_columns,
    reduced_optimizer_code =
      reduced_diagnostics$optimizer_code %||% NA_character_,
    reduced_optimizer_messages =
      reduced_diagnostics$optimizer_messages %||% NA_character_,
    reduced_convergence_status =
      reduced_diagnostics$convergence_status %||% NA_character_,
    reduced_convergence_warnings =
      reduced_diagnostics$convergence_warnings %||% NA_character_,
    reduced_random_effect_structure =
      reduced_diagnostics$random_effect_structure %||% NA_character_,
    reduced_random_intercept_variance =
      reduced_diagnostics$random_intercept_variance %||% NA_real_,
    reduced_residual_variance =
      reduced_diagnostics$residual_variance %||% NA_real_,
    reduced_variance_ratio =
      reduced_diagnostics$variance_ratio %||% NA_real_,
    reduced_ICC = reduced_diagnostics$ICC %||% NA_real_,
    reduced_is_singular_lme4 =
      reduced_diagnostics$is_singular_lme4 %||%
      reduced_diagnostics$singular_model %||% NA,
    reduced_boundary_by_variance_ratio =
      reduced_diagnostics$boundary_by_variance_ratio %||% NA,
    reduced_singularity_class =
      reduced_diagnostics$singularity_class %||% NA_character_,
    likelihood_ratio_statistic = likelihood_ratio_statistic,
    likelihood_ratio_df = likelihood_ratio_df,
    likelihood_ratio_p_value = likelihood_ratio_p_value,
    reduced_model_stability_status = reduced_model_stability_status,
    reduced_diagnostic_review_required =
      reduced_diagnostics$diagnostic_review_required %||% NA,
    full_fixed_effect_rank = full_fixed_rank,
    full_fixed_effect_columns = full_fixed_columns,
    full_optimizer_code =
      full_diagnostics$optimizer_code %||% NA_character_,
    full_optimizer_messages =
      full_diagnostics$optimizer_messages %||% NA_character_,
    full_convergence_status =
      full_diagnostics$convergence_status %||% NA_character_,
    full_convergence_warnings =
      full_diagnostics$convergence_warnings %||% NA_character_,
    full_random_effect_structure =
      full_diagnostics$random_effect_structure %||% NA_character_,
    full_random_intercept_variance =
      full_diagnostics$random_intercept_variance %||% NA_real_,
    full_residual_variance =
      full_diagnostics$residual_variance %||% NA_real_,
    full_variance_ratio = full_diagnostics$variance_ratio %||% NA_real_,
    full_ICC = full_diagnostics$ICC %||% NA_real_,
    full_is_singular_lme4 =
      full_diagnostics$is_singular_lme4 %||%
      full_diagnostics$singular_model %||% NA,
    full_boundary_by_variance_ratio =
      full_diagnostics$boundary_by_variance_ratio %||% NA,
    full_singularity_class =
      full_diagnostics$singularity_class %||% NA_character_,
    full_model_stability_status = full_model_stability_status,
    full_diagnostic_review_required =
      full_diagnostics$diagnostic_review_required %||% NA,
    membership_version = contract$membership_version,
    identity_contract_version = contract$contract_version,
    identity_contract_sha256 = contract$identity_contract_sha256,
    frozen_state_sha256 = contract$frozen_state_sha256,
    endpoint_construction_method = endpoint$endpoint_construction_method,
    endpoint_provenance_status = endpoint$endpoint_provenance_status,
    stringsAsFactors = FALSE
  )
  wgcna_group_bind_schema(list(row), wgcna_group_primary_schema())
}

wgcna_group_validation_row <- function(
    dataset, level, endpoint, contract, effect_scope, spatial_unit,
    attempt_status, failure_reason, requested_formula, actual_formula,
    model_type, fixed_rank, fixed_columns, diagnostics, dat,
    n_contrasts_requested, n_contrasts_eligible, n_contrasts_emitted,
    excluded_contrasts = "none", emmeans_status = "not_run",
    reduced_formula = NA_character_, full_formula = NA_character_,
    identical_rows_verified = NA, reduced_row_hash = NA_character_,
    full_row_hash = NA_character_, likelihood_ratio_statistic = NA_real_,
    likelihood_ratio_df = NA_real_, likelihood_ratio_p_value = NA_real_,
    reduced_model_stability_status = NA_character_,
    full_model_stability_status = NA_character_,
    reduced_diagnostics = NULL, full_diagnostics = NULL,
    reduced_fixed_rank = NA_integer_, reduced_fixed_columns = NA_integer_,
    full_fixed_rank = NA_integer_, full_fixed_columns = NA_integer_
) {
  counts <- wgcna_group_counts(dat)
  model_diagnostic_scope <- if (identical(
    effect_scope, "stress_by_spatial_interaction"
  )) {
    "composite_reduced_full"
  } else {
    "single_fitted_model"
  }
  if (identical(model_diagnostic_scope, "composite_reduced_full") &&
      !isTRUE(diagnostics$model_valid_for_inference) &&
      !identical(
        diagnostics$singularity_class,
        "invalid_composite_reduced_full"
      )) {
    diagnostics <- wgcna_group_invalid_composite_diagnostics(
      diagnostics, failure_reason
    )
  }
  values <- data.frame(
    dataset = dataset, level = level, endpoint_id = endpoint$endpoint_id,
    effect_scope = effect_scope, spatial_unit = spatial_unit,
    attempt_status = attempt_status, failure_reason = failure_reason,
    requested_formula = requested_formula, actual_formula = actual_formula,
    model_type = model_type, fixed_effect_rank = fixed_rank,
    fixed_effect_columns = fixed_columns,
    rank_deficient_model = is.finite(fixed_rank) &&
      is.finite(fixed_columns) && fixed_rank < fixed_columns,
    singular_model = diagnostics$singular_model %||% NA,
    convergence_status = diagnostics$convergence_status %||% "not_fitted",
    convergence_warnings = diagnostics$convergence_warnings %||% NA_character_,
    optimizer_messages = diagnostics$optimizer_messages %||% NA_character_,
    optimizer_code = diagnostics$optimizer_code %||% NA_character_,
    model_warnings = diagnostics$model_warnings %||% NA_character_,
    residual_df = NA_real_, n_samples = counts$n_samples,
    n_animals = counts$n_animals,
    n_samples_per_group = counts$n_samples_per_group,
    n_animals_per_group = counts$n_animals_per_group,
    n_samples_per_spatial_unit = counts$n_samples_per_spatial_unit,
    n_animals_per_spatial_unit = counts$n_animals_per_spatial_unit,
    n_contrasts_requested = n_contrasts_requested,
    n_contrasts_eligible = n_contrasts_eligible,
    n_contrasts_emitted = n_contrasts_emitted,
    excluded_contrasts = excluded_contrasts,
    emmeans_status = emmeans_status,
    has_repeated_animals = any(table(dat$AnimalID) > 1L),
    animal_random_effect_used = identical(model_type, "lmerTest_lmer"),
    AnimalID_source = paste(sort(unique(dat$AnimalID_source)), collapse = ";"),
    model_valid_for_inference =
      diagnostics$model_valid_for_inference %||% FALSE,
    model_diagnostic_scope = model_diagnostic_scope,
    model_stability_status = diagnostics$model_stability_status %||% "invalid",
    primary_model_stable = diagnostics$primary_model_stable %||% FALSE,
    claim_allowed_model = diagnostics$claim_allowed_model %||% FALSE,
    enters_primary_fdr = FALSE,
    manuscript_claim_ready = "not_assessed_stage05",
    random_effect_structure =
      diagnostics$random_effect_structure %||% NA_character_,
    random_intercept_variance =
      diagnostics$random_intercept_variance %||% NA_real_,
    residual_variance = diagnostics$residual_variance %||% NA_real_,
    variance_ratio = diagnostics$variance_ratio %||% NA_real_,
    ICC = diagnostics$ICC %||% NA_real_,
    is_singular_lme4 = if (!is.null(diagnostics$is_singular_lme4)) {
      diagnostics$is_singular_lme4
    } else {
      diagnostics$singular_model %||% NA
    },
    boundary_by_variance_ratio =
      diagnostics$boundary_by_variance_ratio %||% NA,
    boundary_variance_ratio_tolerance =
      wgcna_group_boundary_variance_ratio_tolerance(),
    lme4_singularity_tolerance =
      wgcna_group_lme4_singularity_tolerance(),
    singularity_diagnostic_status =
      diagnostics$singularity_diagnostic_status %||% NA_character_,
    diagnostic_review_required =
      diagnostics$diagnostic_review_required %||% FALSE,
    singularity_class = diagnostics$singularity_class %||% NA_character_,
    boundary_warning = diagnostics$boundary_warning %||% NA_character_,
    reduced_formula = reduced_formula, full_formula = full_formula,
    identical_rows_verified = identical_rows_verified,
    reduced_row_hash = reduced_row_hash, full_row_hash = full_row_hash,
    reduced_fixed_effect_rank = reduced_fixed_rank,
    reduced_fixed_effect_columns = reduced_fixed_columns,
    reduced_optimizer_code =
      reduced_diagnostics$optimizer_code %||% NA_character_,
    reduced_optimizer_messages =
      reduced_diagnostics$optimizer_messages %||% NA_character_,
    reduced_convergence_status =
      reduced_diagnostics$convergence_status %||% NA_character_,
    reduced_convergence_warnings =
      reduced_diagnostics$convergence_warnings %||% NA_character_,
    reduced_random_effect_structure =
      reduced_diagnostics$random_effect_structure %||% NA_character_,
    reduced_random_intercept_variance =
      reduced_diagnostics$random_intercept_variance %||% NA_real_,
    reduced_residual_variance =
      reduced_diagnostics$residual_variance %||% NA_real_,
    reduced_variance_ratio =
      reduced_diagnostics$variance_ratio %||% NA_real_,
    reduced_ICC = reduced_diagnostics$ICC %||% NA_real_,
    reduced_is_singular_lme4 =
      reduced_diagnostics$is_singular_lme4 %||%
      reduced_diagnostics$singular_model %||% NA,
    reduced_boundary_by_variance_ratio =
      reduced_diagnostics$boundary_by_variance_ratio %||% NA,
    reduced_singularity_class =
      reduced_diagnostics$singularity_class %||% NA_character_,
    likelihood_ratio_statistic = likelihood_ratio_statistic,
    likelihood_ratio_df = likelihood_ratio_df,
    likelihood_ratio_p_value = likelihood_ratio_p_value,
    reduced_model_stability_status = reduced_model_stability_status,
    reduced_diagnostic_review_required =
      reduced_diagnostics$diagnostic_review_required %||% NA,
    full_fixed_effect_rank = full_fixed_rank,
    full_fixed_effect_columns = full_fixed_columns,
    full_optimizer_code =
      full_diagnostics$optimizer_code %||% NA_character_,
    full_optimizer_messages =
      full_diagnostics$optimizer_messages %||% NA_character_,
    full_convergence_status =
      full_diagnostics$convergence_status %||% NA_character_,
    full_convergence_warnings =
      full_diagnostics$convergence_warnings %||% NA_character_,
    full_random_effect_structure =
      full_diagnostics$random_effect_structure %||% NA_character_,
    full_random_intercept_variance =
      full_diagnostics$random_intercept_variance %||% NA_real_,
    full_residual_variance =
      full_diagnostics$residual_variance %||% NA_real_,
    full_variance_ratio = full_diagnostics$variance_ratio %||% NA_real_,
    full_ICC = full_diagnostics$ICC %||% NA_real_,
    full_is_singular_lme4 =
      full_diagnostics$is_singular_lme4 %||%
      full_diagnostics$singular_model %||% NA,
    full_boundary_by_variance_ratio =
      full_diagnostics$boundary_by_variance_ratio %||% NA,
    full_singularity_class =
      full_diagnostics$singularity_class %||% NA_character_,
    full_model_stability_status = full_model_stability_status,
    full_diagnostic_review_required =
      full_diagnostics$diagnostic_review_required %||% NA,
    membership_version = contract$membership_version,
    identity_contract_version = contract$contract_version,
    identity_contract_sha256 = contract$identity_contract_sha256,
    frozen_state_sha256 = contract$frozen_state_sha256,
    stringsAsFactors = FALSE
  )
  wgcna_group_bind_schema(list(values), wgcna_group_model_validation_schema())
}

wgcna_group_invalid_diagnostics <- function(reason, model_type = NA_character_) {
  list(
    model_valid_for_inference = FALSE,
    model_stability_status = "invalid",
    primary_model_stable = FALSE,
    claim_allowed_model = FALSE,
    singularity_class = "not_classified_invalid",
    boundary_warning = NA_character_,
    failure_reason = reason,
    variance_ratio = NA_real_, ICC = NA_real_,
    is_singular_lme4 = NA, boundary_by_variance_ratio = NA,
    boundary_variance_ratio_tolerance =
      wgcna_group_boundary_variance_ratio_tolerance(),
    lme4_singularity_tolerance =
      wgcna_group_lme4_singularity_tolerance(),
    singularity_diagnostic_status = if (identical(model_type, "lm")) {
      "not_applicable_lm"
    } else {
      NA_character_
    },
    diagnostic_review_required = FALSE,
    random_effect_structure = if (identical(model_type, "lm")) {
      "none"
    } else {
      NA_character_
    },
    random_intercept_variance = NA_real_, residual_variance = NA_real_,
    singular_model = NA, convergence_status = "not_fitted",
    convergence_warnings = NA_character_, optimizer_messages = NA_character_,
    optimizer_code = NA_character_, model_warnings = NA_character_
  )
}

wgcna_group_fit_attempt <- function(
    dat, dataset, level, endpoint, effect_scope, contract,
    spatial_unit = NA_character_
) {
  if (!is.na(spatial_unit)) {
    dat <- dat[dat$SpatialUnit == spatial_unit, , drop = FALSE]
  }
  dat <- dat[is.finite(dat$eigengene), , drop = FALSE]
  dat$StressGroup <- droplevels(factor(dat$StressGroup))
  dat$SpatialUnit <- droplevels(factor(dat$SpatialUnit))
  output_spatial <- if (!is.na(spatial_unit)) {
    spatial_unit
  } else {
    "global_spatial_adjusted"
  }
  model_type <- tryCatch(
    wgcna_group_model_type_for_scope(dat$AnimalID),
    error = function(e) NA_character_
  )
  repeated <- identical(model_type, "lmerTest_lmer")
  fixed_formula <- if (identical(effect_scope, "spatial_adjusted_global")) {
    eigengene ~ StressGroup + SpatialUnit
  } else {
    eigengene ~ StressGroup
  }
  actual_formula <- if (repeated) {
    stats::as.formula(
      paste(deparse(fixed_formula), "+ (1 | AnimalID)")
    )
  } else {
    fixed_formula
  }
  requested <- paste(deparse(actual_formula), collapse = "")
  empty <- wgcna_group_primary_schema()
  if (nrow(dat) < 4L || nlevels(dat$StressGroup) < 2L ||
      is.na(model_type)) {
    diagnostics <- wgcna_group_invalid_diagnostics(
      "too_few_rows_groups_or_invalid_AnimalID", model_type
    )
    validation <- wgcna_group_validation_row(
      dataset, level, endpoint, contract, effect_scope, output_spatial,
      "failed_preflight", diagnostics$failure_reason,
      requested, NA_character_, model_type,
      NA_integer_, NA_integer_, diagnostics, dat,
      3L, 0L, 0L, "all:failed_preflight"
    )
    return(list(primary = empty, validation = validation))
  }
  fixed_matrix <- tryCatch(
    stats::model.matrix(fixed_formula, data = dat), error = identity
  )
  if (inherits(fixed_matrix, "error")) {
    diagnostics <- wgcna_group_invalid_diagnostics(
      paste0("fixed_effect_matrix_failed:", conditionMessage(fixed_matrix)),
      model_type
    )
    validation <- wgcna_group_validation_row(
      dataset, level, endpoint, contract, effect_scope, output_spatial,
      "failed_preflight", diagnostics$failure_reason,
      requested, NA_character_, model_type,
      NA_integer_, NA_integer_, diagnostics, dat,
      3L, 0L, 0L, "all:fixed_effect_matrix_failed"
    )
    return(list(primary = empty, validation = validation))
  }
  fixed_rank <- qr(fixed_matrix)$rank
  fixed_columns <- ncol(fixed_matrix)
  fit_capture <- wgcna_group_capture(
    if (repeated) {
      lmerTest::lmer(actual_formula, data = dat, REML = FALSE)
    } else {
      stats::lm(actual_formula, data = dat)
    }
  )
  if (inherits(fit_capture$value, "error")) {
    diagnostics <- wgcna_group_invalid_diagnostics(
      paste0("model_fit_failed:", conditionMessage(fit_capture$value)),
      model_type
    )
    validation <- wgcna_group_validation_row(
      dataset, level, endpoint, contract, effect_scope, output_spatial,
      "failed_fit", diagnostics$failure_reason,
      requested, requested, model_type,
      fixed_rank, fixed_columns, diagnostics, dat,
      3L, 0L, 0L, "all:model_fit_failed"
    )
    return(list(primary = empty, validation = validation))
  }
  fit <- fit_capture$value
  diagnostics <- wgcna_group_fit_diagnostics(
    fit, model_type, fixed_rank, fixed_columns, fit_capture$warnings
  )
  if (!isTRUE(diagnostics$model_valid_for_inference)) {
    validation <- wgcna_group_validation_row(
      dataset, level, endpoint, contract, effect_scope, output_spatial,
      "failed_model_validity", diagnostics$failure_reason,
      requested, paste(deparse(stats::formula(fit)), collapse = ""),
      model_type, fixed_rank, fixed_columns, diagnostics, dat,
      3L, 0L, 0L, paste0("all:", diagnostics$failure_reason)
    )
    return(list(primary = empty, validation = validation))
  }
  emm_capture <- wgcna_group_capture(
    emmeans::emmeans(fit, ~ StressGroup)
  )
  methods <- wgcna_group_predeclared_contrasts(levels(dat$StressGroup))
  if (inherits(emm_capture$value, "error") || !length(methods)) {
    reason <- if (inherits(emm_capture$value, "error")) {
      conditionMessage(emm_capture$value)
    } else {
      "no_predeclared_contrast_has_both_groups"
    }
    invalid <- wgcna_group_fit_diagnostics(
      fit, model_type, fixed_rank, fixed_columns,
      c(fit_capture$warnings, emm_capture$warnings),
      contrast_estimable = FALSE, finite_contrast = FALSE
    )
    validation <- wgcna_group_validation_row(
      dataset, level, endpoint, contract, effect_scope, output_spatial,
      "failed_emmeans", reason, requested,
      paste(deparse(stats::formula(fit)), collapse = ""),
      model_type, fixed_rank, fixed_columns, invalid, dat,
      3L, 0L, 0L, paste0("all:", reason), "failed"
    )
    return(list(primary = empty, validation = validation))
  }
  contrast_capture <- wgcna_group_capture(
    wgcna_group_emmeans_contrast_summary(
      emm_capture$value, methods, ci_level = 0.95
    )
  )
  if (inherits(contrast_capture$value, "error")) {
    invalid <- wgcna_group_fit_diagnostics(
      fit, model_type, fixed_rank, fixed_columns,
      c(fit_capture$warnings, emm_capture$warnings, contrast_capture$warnings),
      contrast_estimable = FALSE, finite_contrast = FALSE
    )
    reason <- conditionMessage(contrast_capture$value)
    validation <- wgcna_group_validation_row(
      dataset, level, endpoint, contract, effect_scope, output_spatial,
      "failed_emmeans", reason, requested,
      paste(deparse(stats::formula(fit)), collapse = ""),
      model_type, fixed_rank, fixed_columns, invalid, dat,
      3L, 0L, 0L, paste0("all:", reason), "failed"
    )
    return(list(primary = empty, validation = validation))
  }
  contrast_df <- contrast_capture$value
  stat_col <- intersect(c("t.ratio", "z.ratio"), names(contrast_df))
  stat_col <- if (length(stat_col)) stat_col[[1]] else NA_character_
  pairs <- list(
    "RES - CON" = c("RES", "CON"),
    "SUS - CON" = c("SUS", "CON"),
    "SUS - RES" = c("SUS", "RES")
  )
  rows <- list()
  excluded <- character()
  eligible <- 0L
  for (i in seq_len(nrow(contrast_df))) {
    label <- as.character(contrast_df$contrast[[i]])
    if (!label %in% names(pairs)) next
    pair <- pairs[[label]]
    support <- tapply(
      dat$AnimalID[dat$StressGroup %in% pair],
      dat$StressGroup[dat$StressGroup %in% pair],
      function(x) length(unique(stats::na.omit(x)))
    )
    support <- support[pair]
    min_animals <- if (length(support) == 2L && !anyNA(support)) {
      min(as.integer(support))
    } else {
      0L
    }
    if (min_animals < 3L) {
      excluded <- c(
        excluded, paste0(output_spatial, "|", label,
                         ":insufficient_unique_animals")
      )
      next
    }
    eligible <- eligible + 1L
    finite <- !is.na(stat_col) && all(is.finite(c(
      contrast_df$estimate[[i]], contrast_df$SE[[i]],
      contrast_df[[stat_col]][[i]], contrast_df$p.value[[i]]
    )))
    row_diagnostics <- wgcna_group_fit_diagnostics(
      fit, model_type, fixed_rank, fixed_columns,
      unique(c(
        fit_capture$warnings, emm_capture$warnings,
        contrast_capture$warnings
      )),
      contrast_estimable = finite, finite_contrast = finite
    )
    if (!isTRUE(row_diagnostics$model_valid_for_inference)) {
      excluded <- c(
        excluded, paste0(output_spatial, "|", label, ":",
                         row_diagnostics$failure_reason)
      )
      next
    }
    rows[[length(rows) + 1L]] <- wgcna_group_make_effect_row(
      dat, dataset, level, endpoint, contract, effect_scope, output_spatial,
      label, contrast_df$estimate[[i]], contrast_df$SE[[i]],
      contrast_df[[stat_col]][[i]], contrast_df$p.value[[i]],
      model_type, fit, row_diagnostics, fixed_rank, fixed_columns,
      min_animals,
      CI_low = contrast_df$CI_low[[i]],
      CI_high = contrast_df$CI_high[[i]],
      CI_method = contrast_df$CI_method[[i]],
      CI_level = contrast_df$CI_level[[i]],
      CI_df_method = contrast_df$CI_df_method[[i]],
      df_num = 1,
      df_den = if ("df" %in% names(contrast_df)) {
        suppressWarnings(as.numeric(contrast_df$df[[i]]))
      } else {
        NA_real_
      }
    )
  }
  primary <- wgcna_group_bind_schema(rows, wgcna_group_primary_schema())
  status <- if (nrow(primary)) {
    if (length(excluded)) "success_with_excluded_contrasts" else "success"
  } else {
    "failed_contrasts"
  }
  validation <- wgcna_group_validation_row(
    dataset, level, endpoint, contract, effect_scope, output_spatial,
    status, if (nrow(primary)) "none" else "no_valid_estimable_contrasts",
    requested, paste(deparse(stats::formula(fit)), collapse = ""),
    model_type, fixed_rank, fixed_columns, diagnostics, dat,
    3L, eligible, nrow(primary),
    if (length(excluded)) paste(excluded, collapse = ";") else "none",
    "success"
  )
  list(primary = primary, validation = validation)
}

wgcna_group_fit_interaction_omnibus <- function(
    dat, dataset, level, endpoint, contract
) {
  dat <- dat[is.finite(dat$eigengene), , drop = FALSE]
  dat$StressGroup <- droplevels(factor(dat$StressGroup))
  dat$SpatialUnit <- droplevels(factor(dat$SpatialUnit))
  reduced_formula <- eigengene ~ StressGroup + SpatialUnit + (1 | AnimalID)
  full_formula <- eigengene ~ StressGroup * SpatialUnit + (1 | AnimalID)
  reduced_text <- paste(deparse(reduced_formula), collapse = "")
  full_text <- paste(deparse(full_formula), collapse = "")
  empty <- wgcna_group_primary_schema()
  empty_followup <- wgcna_group_primary_schema()
  model_type <- "lmerTest_lmer"
  cells <- table(dat$StressGroup, dat$SpatialUnit)
  if (nrow(cells) < 2L || ncol(cells) < 2L || any(cells == 0L) ||
      !any(table(dat$AnimalID) > 1L)) {
    diagnostics <- wgcna_group_invalid_diagnostics(
      "incomplete_or_nonrepeated_interaction_design", model_type
    )
    validation <- wgcna_group_validation_row(
      dataset, level, endpoint, contract, "stress_by_spatial_interaction",
      "all_spatial_units", "failed_preflight", diagnostics$failure_reason,
      full_text, NA_character_, model_type,
      NA_integer_, NA_integer_, diagnostics, dat,
      1L, 0L, 0L, "interaction_omnibus:failed_preflight",
      reduced_formula = reduced_text, full_formula = full_text
    )
    return(list(
      primary = empty, followup = empty_followup, validation = validation
    ))
  }
  reduced_fixed <- stats::model.matrix(
    eigengene ~ StressGroup + SpatialUnit, data = dat
  )
  full_fixed <- stats::model.matrix(
    eigengene ~ StressGroup * SpatialUnit, data = dat
  )
  reduced_capture <- wgcna_group_capture(
    lmerTest::lmer(reduced_formula, data = dat, REML = FALSE)
  )
  full_capture <- wgcna_group_capture(
    lmerTest::lmer(full_formula, data = dat, REML = FALSE)
  )
  if (inherits(reduced_capture$value, "error") ||
      inherits(full_capture$value, "error")) {
    reason <- paste(
      if (inherits(reduced_capture$value, "error")) {
        paste0("reduced:", conditionMessage(reduced_capture$value))
      },
      if (inherits(full_capture$value, "error")) {
        paste0("full:", conditionMessage(full_capture$value))
      },
      collapse = ";"
    )
    diagnostics <- wgcna_group_invalid_diagnostics(reason, model_type)
    validation <- wgcna_group_validation_row(
      dataset, level, endpoint, contract, "stress_by_spatial_interaction",
      "all_spatial_units", "failed_fit", reason, full_text, NA_character_,
      model_type, qr(full_fixed)$rank, ncol(full_fixed), diagnostics, dat,
      1L, 0L, 0L, "interaction_omnibus:model_fit_failed",
      reduced_formula = reduced_text, full_formula = full_text
    )
    return(list(
      primary = empty, followup = empty_followup, validation = validation
    ))
  }
  reduced <- reduced_capture$value
  full <- full_capture$value
  reduced_diagnostics <- wgcna_group_fit_diagnostics(
    reduced, model_type, qr(reduced_fixed)$rank, ncol(reduced_fixed),
    reduced_capture$warnings
  )
  full_diagnostics <- wgcna_group_fit_diagnostics(
    full, model_type, qr(full_fixed)$rank, ncol(full_fixed),
    full_capture$warnings
  )
  model_row_text <- function(fit) {
    fitted_rows <- rownames(stats::model.frame(fit))
    fitted_dat <- dat[fitted_rows, , drop = FALSE]
    paste(
      fitted_rows, fitted_dat$AnimalID, fitted_dat$StressGroup,
      fitted_dat$SpatialUnit,
      format(fitted_dat$eigengene, digits = 17, scientific = TRUE),
      sep = "|", collapse = "\n"
    )
  }
  reduced_rows <- rownames(stats::model.frame(reduced))
  full_rows <- rownames(stats::model.frame(full))
  # Reduced/full row hashes are frozen v4 provenance identifiers. Keep their
  # explicit legacy byte contract separate from the v5 LF aggregate hash.
  reduced_hash <- wgcna_group_model_row_sha256_v4(model_row_text(reduced))
  full_hash <- wgcna_group_model_row_sha256_v4(model_row_text(full))
  identical_rows <- identical(reduced_rows, full_rows) &&
    identical(reduced_hash, full_hash)
  if (!isTRUE(reduced_diagnostics$model_valid_for_inference) ||
      !isTRUE(full_diagnostics$model_valid_for_inference) ||
      !identical_rows) {
    reason <- paste(c(
      if (!isTRUE(reduced_diagnostics$model_valid_for_inference)) {
        paste0("reduced:", reduced_diagnostics$failure_reason)
      },
      if (!isTRUE(full_diagnostics$model_valid_for_inference)) {
        paste0("full:", full_diagnostics$failure_reason)
      },
      if (!identical_rows) "reduced_full_rows_not_identical"
    ), collapse = ";")
    invalid <- full_diagnostics
    invalid$model_valid_for_inference <- FALSE
    invalid$model_stability_status <- "invalid"
    invalid$primary_model_stable <- FALSE
    invalid$claim_allowed_model <- FALSE
    invalid$failure_reason <- reason
    validation <- wgcna_group_validation_row(
      dataset, level, endpoint, contract, "stress_by_spatial_interaction",
      "all_spatial_units", "failed_model_validity", reason,
      full_text, full_text, model_type,
      qr(full_fixed)$rank, ncol(full_fixed), invalid, dat,
      1L, 0L, 0L, paste0("interaction_omnibus:", reason),
      reduced_formula = reduced_text, full_formula = full_text,
      identical_rows_verified = identical_rows,
      reduced_row_hash = reduced_hash, full_row_hash = full_hash,
      reduced_model_stability_status =
        reduced_diagnostics$model_stability_status,
      full_model_stability_status = full_diagnostics$model_stability_status,
      reduced_diagnostics = reduced_diagnostics,
      full_diagnostics = full_diagnostics,
      reduced_fixed_rank = qr(reduced_fixed)$rank,
      reduced_fixed_columns = ncol(reduced_fixed),
      full_fixed_rank = qr(full_fixed)$rank,
      full_fixed_columns = ncol(full_fixed)
    )
    return(list(
      primary = empty, followup = empty_followup, validation = validation
    ))
  }
  lrt_capture <- wgcna_group_capture(
    as.data.frame(stats::anova(reduced, full, refit = FALSE))
  )
  if (inherits(lrt_capture$value, "error") ||
      nrow(lrt_capture$value) < 2L) {
    reason <- if (inherits(lrt_capture$value, "error")) {
      conditionMessage(lrt_capture$value)
    } else {
      "nested_likelihood_ratio_result_missing"
    }
    invalid <- full_diagnostics
    invalid$model_valid_for_inference <- FALSE
    invalid$model_stability_status <- "invalid"
    invalid$primary_model_stable <- FALSE
    invalid$claim_allowed_model <- FALSE
    invalid$failure_reason <- reason
    validation <- wgcna_group_validation_row(
      dataset, level, endpoint, contract, "stress_by_spatial_interaction",
      "all_spatial_units", "failed_likelihood_ratio_test", reason,
      full_text, full_text, model_type,
      qr(full_fixed)$rank, ncol(full_fixed), invalid, dat,
      1L, 0L, 0L, paste0("interaction_omnibus:", reason),
      reduced_formula = reduced_text, full_formula = full_text,
      identical_rows_verified = identical_rows,
      reduced_row_hash = reduced_hash, full_row_hash = full_hash,
      reduced_model_stability_status =
        reduced_diagnostics$model_stability_status,
      full_model_stability_status = full_diagnostics$model_stability_status,
      reduced_diagnostics = reduced_diagnostics,
      full_diagnostics = full_diagnostics,
      reduced_fixed_rank = qr(reduced_fixed)$rank,
      reduced_fixed_columns = ncol(reduced_fixed),
      full_fixed_rank = qr(full_fixed)$rank,
      full_fixed_columns = ncol(full_fixed)
    )
    return(list(
      primary = empty, followup = empty_followup, validation = validation
    ))
  }
  lrt <- lrt_capture$value
  lrt_stat <- suppressWarnings(as.numeric(lrt$Chisq[[2]]))
  lrt_df_column <- intersect(c("Df", "Chi Df"), names(lrt))
  lrt_df <- if (length(lrt_df_column) == 1L) {
    suppressWarnings(as.numeric(lrt[[lrt_df_column]][[2]]))
  } else {
    NA_real_
  }
  lrt_p <- suppressWarnings(as.numeric(lrt$`Pr(>Chisq)`[[2]]))
  finite_lrt <- isTRUE(
    all(is.finite(c(lrt_stat, lrt_df, lrt_p))) &&
      lrt_df > 0 && lrt_p >= 0 && lrt_p <= 1
  )
  if (!isTRUE(finite_lrt)) {
    invalid <- full_diagnostics
    invalid$model_valid_for_inference <- FALSE
    invalid$model_stability_status <- "invalid"
    invalid$primary_model_stable <- FALSE
    invalid$claim_allowed_model <- FALSE
    invalid$failure_reason <- "nonfinite_likelihood_ratio_result"
    validation <- wgcna_group_validation_row(
      dataset, level, endpoint, contract, "stress_by_spatial_interaction",
      "all_spatial_units", "failed_likelihood_ratio_test",
      invalid$failure_reason, full_text, full_text, model_type,
      qr(full_fixed)$rank, ncol(full_fixed), invalid, dat,
      1L, 0L, 0L, "interaction_omnibus:nonfinite_result",
      reduced_formula = reduced_text, full_formula = full_text,
      identical_rows_verified = identical_rows,
      reduced_row_hash = reduced_hash, full_row_hash = full_hash,
      reduced_model_stability_status =
        reduced_diagnostics$model_stability_status,
      full_model_stability_status = full_diagnostics$model_stability_status,
      reduced_diagnostics = reduced_diagnostics,
      full_diagnostics = full_diagnostics,
      reduced_fixed_rank = qr(reduced_fixed)$rank,
      reduced_fixed_columns = ncol(reduced_fixed),
      full_fixed_rank = qr(full_fixed)$rank,
      full_fixed_columns = ncol(full_fixed)
    )
    return(list(
      primary = empty, followup = empty_followup, validation = validation
    ))
  }
  combined_diagnostics <- wgcna_group_composite_interaction_diagnostics(
    reduced_diagnostics,
    full_diagnostics,
    identical_rows_verified = identical_rows,
    likelihood_ratio_valid = finite_lrt
  )
  counts <- wgcna_group_counts(dat)
  primary <- wgcna_group_make_effect_row(
    dat, dataset, level, endpoint, contract,
    "stress_by_spatial_interaction", "all_spatial_units",
    "StressGroup x SpatialUnit omnibus",
    NA_real_, NA_real_, lrt_stat, lrt_p,
    model_type, full, combined_diagnostics,
    qr(full_fixed)$rank, ncol(full_fixed),
    counts$min_animals_per_group,
    test_type = "interaction_omnibus",
    df_num = lrt_df, df_den = NA_real_,
    reduced_formula = reduced_text, full_formula = full_text,
    identical_rows_verified = identical_rows,
    reduced_row_hash = reduced_hash, full_row_hash = full_hash,
    likelihood_ratio_statistic = lrt_stat,
    likelihood_ratio_df = lrt_df, likelihood_ratio_p_value = lrt_p,
    reduced_model_stability_status =
      reduced_diagnostics$model_stability_status,
    full_model_stability_status = full_diagnostics$model_stability_status,
    reduced_diagnostics = reduced_diagnostics,
    full_diagnostics = full_diagnostics,
    reduced_fixed_rank = qr(reduced_fixed)$rank,
    reduced_fixed_columns = ncol(reduced_fixed),
    full_fixed_rank = qr(full_fixed)$rank,
    full_fixed_columns = ncol(full_fixed)
  )
  followup <- list()
  emm_capture <- wgcna_group_capture(
    emmeans::emmeans(full, ~ StressGroup | SpatialUnit)
  )
  methods <- wgcna_group_predeclared_contrasts(levels(dat$StressGroup))
  if (!inherits(emm_capture$value, "error") && length(methods)) {
    contrast_capture <- wgcna_group_capture(
      wgcna_group_emmeans_contrast_summary(
        emm_capture$value, methods, ci_level = 0.95
      )
    )
    if (!inherits(contrast_capture$value, "error")) {
      contrast_df <- contrast_capture$value
      stat_col <- intersect(c("t.ratio", "z.ratio"), names(contrast_df))
      stat_col <- if (length(stat_col)) stat_col[[1]] else NA_character_
      pairs <- list(
        "RES - CON" = c("RES", "CON"),
        "SUS - CON" = c("SUS", "CON"),
        "SUS - RES" = c("SUS", "RES")
      )
      for (i in seq_len(nrow(contrast_df))) {
        label <- as.character(contrast_df$contrast[[i]])
        unit <- as.character(contrast_df$SpatialUnit[[i]])
        if (!label %in% names(pairs) || is.na(stat_col)) next
        support_dat <- dat[as.character(dat$SpatialUnit) == unit, ,
                           drop = FALSE]
        pair <- pairs[[label]]
        support <- tapply(
          support_dat$AnimalID[support_dat$StressGroup %in% pair],
          support_dat$StressGroup[support_dat$StressGroup %in% pair],
          function(x) length(unique(stats::na.omit(x)))
        )
        support <- support[pair]
        min_animals <- if (length(support) == 2L && !anyNA(support)) {
          min(as.integer(support))
        } else {
          0L
        }
        finite <- min_animals >= 3L && all(is.finite(c(
          contrast_df$estimate[[i]], contrast_df$SE[[i]],
          contrast_df[[stat_col]][[i]], contrast_df$p.value[[i]]
        )))
        if (!finite) next
        followup[[length(followup) + 1L]] <- wgcna_group_make_effect_row(
          support_dat, dataset, level, endpoint, contract,
          "stress_by_spatial_interaction", unit, label,
          contrast_df$estimate[[i]], contrast_df$SE[[i]],
          contrast_df[[stat_col]][[i]], contrast_df$p.value[[i]],
          model_type, full, full_diagnostics,
          qr(full_fixed)$rank, ncol(full_fixed), min_animals,
          CI_low = contrast_df$CI_low[[i]],
          CI_high = contrast_df$CI_high[[i]],
          CI_method = contrast_df$CI_method[[i]],
          CI_level = contrast_df$CI_level[[i]],
          CI_df_method = contrast_df$CI_df_method[[i]],
          test_type = "conditional_interaction_followup",
          df_num = 1,
          df_den = if ("df" %in% names(contrast_df)) {
            suppressWarnings(as.numeric(contrast_df$df[[i]]))
          } else {
            NA_real_
          },
          reduced_formula = reduced_text, full_formula = full_text,
          identical_rows_verified = identical_rows,
          reduced_row_hash = reduced_hash, full_row_hash = full_hash,
          reduced_model_stability_status =
            reduced_diagnostics$model_stability_status,
          full_model_stability_status = full_diagnostics$model_stability_status,
          reduced_diagnostics = reduced_diagnostics,
          full_diagnostics = full_diagnostics,
          reduced_fixed_rank = qr(reduced_fixed)$rank,
          reduced_fixed_columns = ncol(reduced_fixed),
          full_fixed_rank = qr(full_fixed)$rank,
          full_fixed_columns = ncol(full_fixed)
        )
      }
    }
  }
  followup <- wgcna_group_bind_schema(
    followup, wgcna_group_primary_schema()
  )
  validation <- wgcna_group_validation_row(
    dataset, level, endpoint, contract, "stress_by_spatial_interaction",
    "all_spatial_units", "success", "none", full_text, full_text,
    model_type, qr(full_fixed)$rank, ncol(full_fixed),
    combined_diagnostics, dat, 1L, 1L, 1L, "none",
    if (nrow(followup)) "conditional_followup_success" else
      "conditional_followup_not_available",
    reduced_formula = reduced_text, full_formula = full_text,
    identical_rows_verified = identical_rows,
    reduced_row_hash = reduced_hash, full_row_hash = full_hash,
    likelihood_ratio_statistic = lrt_stat,
    likelihood_ratio_df = lrt_df, likelihood_ratio_p_value = lrt_p,
    reduced_model_stability_status =
      reduced_diagnostics$model_stability_status,
    full_model_stability_status = full_diagnostics$model_stability_status,
    reduced_diagnostics = reduced_diagnostics,
    full_diagnostics = full_diagnostics,
    reduced_fixed_rank = qr(reduced_fixed)$rank,
    reduced_fixed_columns = ncol(reduced_fixed),
    full_fixed_rank = qr(full_fixed)$rank,
    full_fixed_columns = ncol(full_fixed)
  )
  list(primary = primary, followup = followup, validation = validation)
}

wgcna_group_clear_fdr_fields <- function(rows) {
  if (!nrow(rows)) return(rows)
  numeric_cols <- c(
    "FDR_primary_global", "FDR_secondary_global",
    "FDR_interaction_omnibus", "FDR_local_exploratory",
    "FDR_conservative_all_tests", "FDR_global",
    "FDR_within_dataset_level", "FDR_dataset_all_levels"
  )
  character_cols <- c(
    "FDR_primary_family_id", "FDR_secondary_family_id",
    "FDR_interaction_family_id", "FDR_local_family_id",
    "FDR_conservative_family_id", "FDR_family_within_level_id",
    "FDR_family_dataset_id"
  )
  integer_cols <- c(
    "n_tests_FDR_primary", "n_tests_FDR_secondary_global",
    "n_tests_FDR_interaction_omnibus",
    "n_tests_FDR_local_exploratory",
    "n_tests_FDR_conservative_all_tests",
    "n_tests_FDR_within_dataset_level",
    "n_tests_FDR_dataset_all_levels"
  )
  for (nm in numeric_cols) rows[[nm]] <- NA_real_
  for (nm in character_cols) rows[[nm]] <- NA_character_
  for (nm in integer_cols) rows[[nm]] <- NA_integer_
  rows$enters_primary_fdr <- FALSE
  rows$primary_contrast <- FALSE
  rows$primary_analysis_tier <- NA_character_
  rows$statistical_support_status <- "sensitivity_not_fdr_adjusted"
  rows$evidence_status <- "sensitivity_not_fdr_adjusted"
  rows
}

wgcna_group_endpoint_sensitivities <- function(
    raw_dat, aggregated, dataset, level, endpoint, contract
) {
  rows <- list()
  hemisphere_aggregated <- list()
  for (hemisphere in c("L", "R")) {
    hemi_raw <- raw_dat[raw_dat$Hemisphere == hemisphere, , drop = FALSE]
    hemi <- wgcna_group_aggregate_endpoint(
      hemi_raw, dataset, level, endpoint$endpoint_id, endpoint$endpoint_col
    )$animal_spatial_values
    hemisphere_aggregated[[hemisphere]] <- hemi
    result <- wgcna_group_fit_attempt(
      hemi, dataset, level, endpoint, "spatial_adjusted_global", contract
    )$primary
    result <- result[result$contrast == "SUS - RES", , drop = FALSE]
    if (nrow(result)) {
      result$analysis_tier <- paste0(
        "sensitivity_", tolower(hemisphere), "_only"
      )
      rows[[length(rows) + 1L]] <- wgcna_group_clear_fdr_fields(result)
    }
  }
  paired <- aggregated[
    aggregated$n_hemispheres_observed == 2L, , drop = FALSE
  ]
  if (nrow(paired)) {
    result <- wgcna_group_fit_attempt(
      paired, dataset, level, endpoint,
      "spatial_adjusted_global", contract
    )$primary
    result <- result[result$contrast == "SUS - RES", , drop = FALSE]
    if (nrow(result)) {
      result$analysis_tier <- "sensitivity_complete_bilateral_pairs"
      rows[[length(rows) + 1L]] <- wgcna_group_clear_fdr_fields(result)
    }
  }
  animal_wide <- aggregated |>
    dplyr::group_by(.data$AnimalID, .data$StressGroup) |>
    dplyr::summarise(
      eigengene = mean(.data$eigengene),
      SpatialUnit = "animal_wide",
      canonical_spatial_unit = "animal_wide",
      SpatialUnitType = "animal_wide",
      AnimalID_source = paste(
        sort(unique(.data$AnimalID_source)), collapse = ";"
      ),
      Region = NA_character_, Layer = NA_character_,
      .groups = "drop"
    )
  result <- wgcna_group_fit_attempt(
    animal_wide, dataset, level, endpoint,
    "within_spatial_unit", contract, "animal_wide"
  )$primary
  result <- result[result$contrast == "SUS - RES", , drop = FALSE]
  if (nrow(result)) {
    result$analysis_tier <- "sensitivity_animal_wide_mean"
    result$effect_scope <- "animal_wide_mean_sensitivity"
    rows[[length(rows) + 1L]] <- wgcna_group_clear_fdr_fields(result)
  }
  left <- hemisphere_aggregated[["L"]]
  right <- hemisphere_aggregated[["R"]]
  concordance <- dplyr::inner_join(
    left |>
      dplyr::select(
        "AnimalID", "SpatialUnit", left_eigengene = "eigengene"
      ),
    right |>
      dplyr::select(
        "AnimalID", "SpatialUnit", right_eigengene = "eigengene"
      ),
    by = c("AnimalID", "SpatialUnit"),
    relationship = "one-to-one"
  )
  pearson <- if (nrow(concordance) >= 3L) {
    stats::cor(
      concordance$left_eigengene, concordance$right_eigengene,
      method = "pearson"
    )
  } else {
    NA_real_
  }
  spearman <- if (nrow(concordance) >= 3L) {
    stats::cor(
      concordance$left_eigengene, concordance$right_eigengene,
      method = "spearman"
    )
  } else {
    NA_real_
  }
  concordance_summary <- data.frame(
    dataset = dataset, level = level, endpoint_id = endpoint$endpoint_id,
    hypothesis_level = wgcna_group_hypothesis_level(level),
    n_paired_animal_spatial_units = nrow(concordance),
    pearson_r = pearson, spearman_rho = spearman,
    concordance_status = if (nrow(concordance) >= 3L) {
      "available"
    } else {
      "insufficient_pairs"
    },
    manuscript_claim_ready = "not_assessed_stage05",
    stringsAsFactors = FALSE
  )
  list(
    sensitivity = wgcna_group_bind_schema(
      rows, wgcna_group_primary_schema()
    ),
    concordance = concordance_summary
  )
}

wgcna_group_run_level <- function(
    dataset, level, eigengenes, endpoint_map, sample_audit, contract
) {
  dat <- sample_audit
  idx <- match(dat$Sample, eigengenes$Sample)
  if (anyNA(idx) || nrow(eigengenes) != nrow(dat)) {
    stop("Endpoint/sample alignment would lose frozen-state samples.",
         call. = FALSE)
  }
  dat$SpatialUnit <- dat$canonical_spatial_unit
  for (nm in setdiff(names(eigengenes), "Sample")) {
    dat[[nm]] <- eigengenes[[nm]][idx]
  }
  primary <- list()
  validation <- list()
  followup <- list()
  aggregation <- list()
  hemisphere_aggregation <- list()
  sensitivity <- list()
  concordance <- list()
  endpoint_map <- endpoint_map[order(endpoint_map$endpoint_id), , drop = FALSE]
  if (identical(level, "supermodule")) {
    endpoint_map <- endpoint_map[
      endpoint_map$n_member_modules > 1L, , drop = FALSE
    ]
  }
  for (i in seq_len(nrow(endpoint_map))) {
    endpoint <- endpoint_map[i, , drop = FALSE]
    aggregation_result <- wgcna_group_aggregate_endpoint(
      dat, dataset, level, endpoint$endpoint_id, endpoint$endpoint_col
    )
    aggregated <- aggregation_result$animal_spatial_values
    aggregation[[length(aggregation) + 1L]] <- aggregated
    hemisphere_aggregation[[length(hemisphere_aggregation) + 1L]] <-
      aggregation_result$hemisphere_values
    sensitivity_result <- wgcna_group_endpoint_sensitivities(
      dat, aggregated, dataset, level, endpoint, contract
    )
    sensitivity[[length(sensitivity) + 1L]] <-
      sensitivity_result$sensitivity
    concordance[[length(concordance) + 1L]] <-
      sensitivity_result$concordance
    spatial_units <- sort(unique(aggregated$SpatialUnit))
    for (unit in spatial_units) {
      result <- wgcna_group_fit_attempt(
        aggregated, dataset, level, endpoint,
        "within_spatial_unit", contract, unit
      )
      primary[[length(primary) + 1L]] <- result$primary
      validation[[length(validation) + 1L]] <- result$validation
    }
    global <- wgcna_group_fit_attempt(
      aggregated, dataset, level, endpoint,
      "spatial_adjusted_global", contract
    )
    primary[[length(primary) + 1L]] <- global$primary
    validation[[length(validation) + 1L]] <- global$validation
    interaction <- wgcna_group_fit_interaction_omnibus(
      aggregated, dataset, level, endpoint, contract
    )
    primary[[length(primary) + 1L]] <- interaction$primary
    validation[[length(validation) + 1L]] <- interaction$validation
    followup[[length(followup) + 1L]] <- interaction$followup
  }
  list(
    primary = wgcna_group_bind_schema(
      primary, wgcna_group_primary_schema()
    ),
    validation = wgcna_group_bind_schema(
      validation, wgcna_group_model_validation_schema()
    ),
    interaction_followup = wgcna_group_bind_schema(
      followup, wgcna_group_primary_schema()
    ),
    animal_spatial_values = dplyr::bind_rows(aggregation),
    hemisphere_values = dplyr::bind_rows(hemisphere_aggregation),
    sensitivity = wgcna_group_bind_schema(
      sensitivity, wgcna_group_primary_schema()
    ),
    left_right_concordance = dplyr::bind_rows(concordance)
  )
}

wgcna_group_apply_bh_family <- function(
    data, idx, fdr_col, family_col, n_col, family_id
) {
  if (!length(idx)) return(data)
  data[[fdr_col]][idx] <- stats::p.adjust(data$p_value[idx], method = "BH")
  data[[family_col]][idx] <- family_id
  data[[n_col]][idx] <- length(idx)
  data
}

wgcna_group_apply_fdr <- function(module_rows, supermodule_rows) {
  all <- dplyr::bind_rows(module_rows, supermodule_rows)
  if (!nrow(all)) {
    return(list(module = module_rows, supermodule = supermodule_rows))
  }
  eligible <- all$model_valid_for_inference %in% TRUE &
    all$independent_hypothesis %in% TRUE &
    is.finite(all$p_value) &
    !(all$fallback_used %in% TRUE)
  if (any(!eligible)) {
    stop("Canonical effect rows must be valid independent hypotheses before alias expansion.",
         call. = FALSE)
  }
  fdr_numeric <- c(
    "FDR_primary_global", "FDR_secondary_global",
    "FDR_interaction_omnibus", "FDR_local_exploratory",
    "FDR_conservative_all_tests", "FDR_global",
    "FDR_within_dataset_level", "FDR_dataset_all_levels"
  )
  for (nm in fdr_numeric) all[[nm]] <- NA_real_
  family_character <- c(
    "FDR_primary_family_id", "FDR_secondary_family_id",
    "FDR_interaction_family_id", "FDR_local_family_id",
    "FDR_conservative_family_id", "FDR_family_within_level_id",
    "FDR_family_dataset_id"
  )
  for (nm in family_character) all[[nm]] <- NA_character_
  family_counts <- c(
    "n_tests_FDR_primary", "n_tests_FDR_secondary_global",
    "n_tests_FDR_interaction_omnibus",
    "n_tests_FDR_local_exploratory",
    "n_tests_FDR_conservative_all_tests",
    "n_tests_FDR_within_dataset_level",
    "n_tests_FDR_dataset_all_levels"
  )
  for (nm in family_counts) all[[nm]] <- NA_integer_
  all$enters_primary_fdr <- FALSE

  datasets <- unique(all$dataset)
  hypothesis_levels <- unique(all$hypothesis_level)
  for (dataset in datasets) {
    for (hypothesis_level in hypothesis_levels) {
      primary_idx <- which(
        eligible &
          all$dataset == dataset &
          all$hypothesis_level == hypothesis_level &
          all$analysis_tier == "primary_wgcna_global" &
          all$contrast == "SUS - RES" &
          all$effect_scope == "spatial_adjusted_global" &
          all$spatial_unit == "global_spatial_adjusted"
      )
      if (length(primary_idx)) {
        family_id <- paste(
          dataset, hypothesis_level, "SUS-RES",
          "spatial_adjusted_global", sep = "::"
        )
        all <- wgcna_group_apply_bh_family(
          all, primary_idx, "FDR_primary_global",
          "FDR_primary_family_id", "n_tests_FDR_primary", family_id
        )
        all$enters_primary_fdr[primary_idx] <- TRUE
      }
      for (contrast in c("RES - CON", "SUS - CON")) {
        secondary_idx <- which(
          eligible &
            all$dataset == dataset &
            all$hypothesis_level == hypothesis_level &
            all$analysis_tier == "secondary_contextual_global" &
            all$contrast == contrast &
            all$effect_scope == "spatial_adjusted_global" &
            all$spatial_unit == "global_spatial_adjusted"
        )
        if (length(secondary_idx)) {
          family_id <- paste(
            dataset, hypothesis_level,
            gsub(" ", "", contrast, fixed = TRUE),
            "spatial_adjusted_global", sep = "::"
          )
          all <- wgcna_group_apply_bh_family(
            all, secondary_idx, "FDR_secondary_global",
            "FDR_secondary_family_id",
            "n_tests_FDR_secondary_global", family_id
          )
        }
      }
      interaction_idx <- which(
        eligible &
          all$dataset == dataset &
          all$hypothesis_level == hypothesis_level &
          all$analysis_tier == "secondary_spatial_heterogeneity" &
          all$test_type == "interaction_omnibus"
      )
      if (length(interaction_idx)) {
        family_id <- paste(
          dataset, hypothesis_level, "interaction_omnibus", sep = "::"
        )
        all <- wgcna_group_apply_bh_family(
          all, interaction_idx, "FDR_interaction_omnibus",
          "FDR_interaction_family_id",
          "n_tests_FDR_interaction_omnibus", family_id
        )
      }
      for (contrast in c("RES - CON", "SUS - CON", "SUS - RES")) {
        local_idx <- which(
          eligible &
            all$dataset == dataset &
            all$hypothesis_level == hypothesis_level &
            all$analysis_tier == "exploratory_spatial_localization" &
            all$contrast == contrast &
            all$effect_scope == "within_spatial_unit"
        )
        if (length(local_idx)) {
          family_id <- paste(
            dataset, hypothesis_level,
            gsub(" ", "", contrast, fixed = TRUE),
            "within_spatial_unit", sep = "::"
          )
          all <- wgcna_group_apply_bh_family(
            all, local_idx, "FDR_local_exploratory",
            "FDR_local_family_id",
            "n_tests_FDR_local_exploratory", family_id
          )
        }
      }
    }
    conservative_idx <- which(
      eligible &
        all$dataset == dataset &
        all$analysis_tier %in% c(
          "primary_wgcna_global",
          "secondary_contextual_global",
          "secondary_spatial_heterogeneity"
        )
    )
    if (length(conservative_idx)) {
      family_id <- paste(
        dataset, "independent_primary_secondary_all_tests", sep = "::"
      )
      all <- wgcna_group_apply_bh_family(
        all, conservative_idx, "FDR_conservative_all_tests",
        "FDR_conservative_family_id",
        "n_tests_FDR_conservative_all_tests", family_id
      )
      all$FDR_global[conservative_idx] <-
        all$FDR_conservative_all_tests[conservative_idx]
    }
  }
  primary_rows <- all$analysis_tier == "primary_wgcna_global"
  secondary_rows <- all$analysis_tier == "secondary_contextual_global"
  if (any(primary_rows & secondary_rows) ||
      any(primary_rows & !is.na(all$FDR_secondary_global)) ||
      any(secondary_rows & !is.na(all$FDR_primary_global)) ||
      any(primary_rows & all$contrast != "SUS - RES") ||
      any(secondary_rows & !all$contrast %in% c("RES - CON", "SUS - CON"))) {
    stop("Primary and secondary global FDR families overlap or contain a prohibited contrast.",
         call. = FALSE)
  }
  applicable_fdr <- dplyr::case_when(
    all$analysis_tier == "primary_wgcna_global" ~ all$FDR_primary_global,
    all$analysis_tier == "secondary_contextual_global" ~
      all$FDR_secondary_global,
    all$analysis_tier == "secondary_spatial_heterogeneity" ~
      all$FDR_interaction_omnibus,
    all$analysis_tier == "exploratory_spatial_localization" ~
      all$FDR_local_exploratory,
    TRUE ~ NA_real_
  )
  all$statistical_support_status <-
    wgcna_group_classify_statistical_support(
      all$p_value, applicable_fdr
    )
  all$evidence_status <- all$statistical_support_status
  all$FDR_method <- "BH"
  key <- paste(
    all$dataset, all$level, all$endpoint_id, all$effect_scope,
    all$spatial_unit, all$contrast, all$test_type, sep = "|"
  )
  if (anyDuplicated(key)) {
    stop("Canonical endpoint statistical keys are duplicated.", call. = FALSE)
  }
  all <- all[order(
    all$level, all$endpoint_id, all$effect_scope,
    all$spatial_unit, all$contrast
  ), names(wgcna_group_primary_schema()), drop = FALSE]
  list(
    module = all[all$level == "module", , drop = FALSE],
    supermodule = all[all$level == "supermodule", , drop = FALSE]
  )
}

wgcna_group_make_singleton_alias_rows <- function(
    module_rows, supermodule_composition
) {
  singleton <- supermodule_composition[
    supermodule_composition$n_member_modules == 1L, , drop = FALSE
  ]
  if (!nrow(singleton)) return(wgcna_group_primary_schema())
  rows <- list()
  fdr_cols <- c(
    "FDR_primary_global", "FDR_secondary_global",
    "FDR_interaction_omnibus", "FDR_local_exploratory",
    "FDR_conservative_all_tests", "FDR_global",
    "FDR_within_dataset_level", "FDR_dataset_all_levels"
  )
  family_cols <- c(
    "FDR_primary_family_id", "FDR_secondary_family_id",
    "FDR_interaction_family_id", "FDR_local_family_id",
    "FDR_conservative_family_id", "FDR_family_within_level_id",
    "FDR_family_dataset_id"
  )
  count_cols <- c(
    "n_tests_FDR_primary", "n_tests_FDR_secondary_global",
    "n_tests_FDR_interaction_omnibus",
    "n_tests_FDR_local_exploratory",
    "n_tests_FDR_conservative_all_tests",
    "n_tests_FDR_within_dataset_level",
    "n_tests_FDR_dataset_all_levels"
  )
  for (i in seq_len(nrow(singleton))) {
    member <- as.character(singleton$member_modules[[i]])
    alias_id <- as.character(singleton$supermodule_id[[i]])
    source <- module_rows[module_rows$module_id == member, , drop = FALSE]
    if (!nrow(source)) {
      stop("Singleton alias member has no canonical module effects: ", member,
           call. = FALSE)
    }
    alias <- source
    alias$level <- "supermodule"
    alias$endpoint_id <- alias_id
    alias$endpoint_label <- alias_id
    alias$module_id <- NA_character_
    alias$supermodule_id <- alias_id
    alias$module_label <- NA_character_
    alias$supermodule_label <- alias_id
    alias$hypothesis_level <- "compatibility_alias"
    alias$canonical_claim_entity_id <- member
    alias$claim_entity_role <- "compatibility_alias"
    alias$support_source_entity_id <- member
    alias$independent_hypothesis <- FALSE
    alias$enters_primary_fdr <- FALSE
    alias$display_allowed <- TRUE
    alias$manuscript_claim_ready <- "not_assessed_stage05"
    alias$statistical_support_status <- "inherited_from_canonical_entity"
    alias$evidence_status <- "inherited_from_canonical_entity"
    for (nm in fdr_cols) alias[[nm]] <- NA_real_
    for (nm in family_cols) alias[[nm]] <- NA_character_
    for (nm in count_cols) alias[[nm]] <- NA_integer_
    rows[[length(rows) + 1L]] <- alias
  }
  aliases <- wgcna_group_bind_schema(rows, wgcna_group_primary_schema())
  if (any(vapply(fdr_cols, function(nm) {
    any(!is.na(aliases[[nm]]))
  }, logical(1)))) {
    stop("Singleton compatibility aliases inherited an independent FDR value.",
         call. = FALSE)
  }
  aliases
}

wgcna_group_supermodule_values <- function(
    dataset, super_eigengenes, super_map, sample_audit, contract
) {
  idx <- match(sample_audit$Sample, super_eigengenes$Sample)
  long <- super_eigengenes[idx, , drop = FALSE] |>
    tidyr::pivot_longer(
      cols = -"Sample", names_to = "supermodule_eigengene",
      values_to = "eigengene"
    ) |>
    dplyr::left_join(
      super_map |>
        dplyr::transmute(
          supermodule_eigengene = .data$endpoint_col,
          supermodule_id = .data$endpoint_id
        ),
      by = "supermodule_eigengene",
      relationship = "many-to-one"
    ) |>
    dplyr::left_join(
      sample_audit |>
        dplyr::select(
          "Sample", "StressGroup", "AnimalID", "AnimalID_source",
          "Region", "Layer", "canonical_spatial_unit", "Sex", "Batch",
          "SpatialUnitType"
        ),
      by = "Sample",
      relationship = "many-to-one"
    ) |>
    dplyr::mutate(
      dataset = dataset,
      supermodule_label = .data$supermodule_id,
      canonical_supermodule_label = NA_character_,
      canonical_supermodule_plot_label = NA_character_,
      membership_version = contract$membership_version,
      identity_contract_version = contract$contract_version,
      identity_contract_sha256 = contract$identity_contract_sha256,
      frozen_state_sha256 = contract$frozen_state_sha256
    ) |>
    dplyr::group_by(.data$dataset, .data$supermodule_id) |>
    dplyr::mutate(
      eigengene_z = if (stats::sd(.data$eigengene) > 0) {
        as.numeric(scale(.data$eigengene))
      } else {
        NA_real_
      }
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(.data$supermodule_id, .data$Sample)
  long
}

wgcna_group_consumer_scan_rules <- function() {
  list(
    included = paste(
      "roots: R/, inst/schemas/, numbered active stage directories 00_-98_,",
      "tests/testthat/, and repository-root .R/.r/.yml/.yaml files; extensions:",
      ".R, .r, .yml, .yaml"
    ),
    excluded = paste(
      "not included by root selection: .git/, 99_deprecated/, results/, data/,",
      "tmp/, docs/, review bundles, generated outputs, tests/fixtures/, and",
      "tests/testthat/fixtures/; explicitly excluded producers, registry, and",
      "Stage 05 fixtures containing deliberate legacy tokens:",
      "R/wgcna_group_effects_utils.R and",
      "06_modules_WGCNA/05_module_supermodule_group_effects.r;",
      "tests/testthat/test-wgcna-group-effects-phase2b.R,",
      "tests/testthat/test-wgcna-group-effects-contract.R, and",
      "tests/testthat/test-schema-validation.R"
    )
  )
}

wgcna_group_consumer_tokens <- function() {
  data.frame(
    scan_token = c(
      "FDR_global",
      "FDR_within_dataset_level",
      "evidence_status",
      "claim_allowed_model",
      "primary_model_stable",
      "stress_by_spatial_interaction",
      "singleton_member_module_eigengene|is_singleton_supermodule|n_member_modules\\s*(==|<=)\\s*1|compatibility_alias"
    ),
    consumed_column = c(
      "FDR_global",
      "FDR_within_dataset_level",
      "evidence_status",
      "claim_allowed_model",
      "primary_model_stable",
      "<row_contract:interaction_conditional>",
      "<row_contract:singleton_supermodule>"
    ),
    current_semantics = c(
      "deprecated combined dataset-level all-tests FDR gate",
      "legacy BH family spanning all scopes and contrasts within level",
      "legacy support classification derived from broad FDR families",
      "legacy model/contrast eligibility flag",
      "legacy stable-fit flag that excluded all singular fits",
      "conditional StressGroup contrasts emitted from interaction models",
      "singleton supermodule rows treated as fitted supermodule endpoints"
    ),
    required_new_column = c(
      "FDR_primary_global or tier-specific FDR_* after row-scope review",
      "tier-specific FDR_primary_global/FDR_secondary_global/FDR_interaction_omnibus/FDR_local_exploratory",
      "statistical_support_status plus analysis_tier",
      "model_valid_for_inference plus model_stability_status",
      "model_stability_status",
      "test_type=interaction_omnibus; conditional rows are exploratory follow-up only",
      "canonical_claim_entity_id, support_source_entity_id, independent_hypothesis"
    ),
    proposed_phase3_action = c(
      "Select the biologically appropriate tier-specific FDR; never use the deprecated alias as the main claim gate.",
      "Migrate row selection and family provenance atomically across all consumers.",
      "Consume the new controlled status only after filtering the intended analysis tier.",
      "Separate validity from stability and manuscript readiness.",
      "Handle boundary_random_intercept_zero with the later animal-level readiness requirement.",
      "Replace primary-table conditional-row assumptions with the omnibus row and optional exploratory follow-up table.",
      "Resolve aliases to canonical modules and prevent duplicate independent claims."
    ),
    stringsAsFactors = FALSE
  )
}

wgcna_group_scan_downstream_consumers <- function(root = repo_root()) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  top_dirs <- list.dirs(root, recursive = FALSE, full.names = FALSE)
  numbered_dirs <- top_dirs[
    grepl("^(0[0-9]|[1-8][0-9]|9[0-8])_", top_dirs)
  ]
  scan_roots <- c(
    file.path(root, "R"),
    file.path(root, "inst", "schemas"),
    file.path(root, "tests", "testthat"),
    file.path(root, numbered_dirs)
  )
  scan_roots <- scan_roots[dir.exists(scan_roots)]
  files <- unlist(lapply(scan_roots, function(path) {
    list.files(
      path, pattern = "\\.(R|r|ya?ml)$", recursive = TRUE,
      full.names = TRUE, include.dirs = FALSE, all.files = TRUE, no.. = TRUE
    )
  }), use.names = FALSE)
  root_files <- list.files(
    root, pattern = "\\.(R|r|ya?ml)$", recursive = FALSE,
    full.names = TRUE, include.dirs = FALSE, all.files = TRUE, no.. = TRUE
  )
  files <- unique(c(files, root_files))
  files <- normalizePath(files, winslash = "/", mustWork = FALSE)
  relative <- substring(files, nchar(root) + 2L)
  relative <- gsub("\\\\", "/", relative)
  excluded <- relative %in% c(
    "R/wgcna_group_effects_utils.R",
    "06_modules_WGCNA/05_module_supermodule_group_effects.r",
    "tests/testthat/test-wgcna-group-effects-phase2b.R",
    "tests/testthat/test-wgcna-group-effects-contract.R",
    "tests/testthat/test-schema-validation.R"
  ) | grepl("^tests/(fixtures|testthat/fixtures)/", relative)
  keep <- !excluded
  files <- files[keep]
  relative <- relative[keep]
  ordering <- order(relative)
  files <- files[ordering]
  relative <- relative[ordering]
  tokens <- wgcna_group_consumer_tokens()
  rules <- wgcna_group_consumer_scan_rules()
  rows <- list()
  for (i in seq_along(files)) {
    text <- readLines(files[[i]], warn = FALSE, encoding = "UTF-8")
    consumer_kind <- dplyr::case_when(
      grepl("^([0-9]{2})_[^/]+/", relative[[i]]) ~ "runtime_stage",
      grepl("^R/", relative[[i]]) ~ "runtime_helper",
      identical(relative[[i]], "pipeline.yml") ~ "pipeline_registry",
      grepl("^inst/schemas/", relative[[i]]) ~ "schema",
      grepl("^tests/", relative[[i]]) ~ "test",
      TRUE ~ "other_source"
    )
    runtime_consumer <- consumer_kind %in% c(
      "runtime_stage", "runtime_helper"
    )
    audit_role <- dplyr::case_when(
      runtime_consumer ~ "runtime_phase3_migration",
      identical(consumer_kind, "pipeline_registry") ~
        "pipeline_registry_contract_update",
      consumer_kind %in% c("schema", "test") ~
        "schema_test_contract_update",
      TRUE ~ "nonruntime_support_reference"
    )
    for (j in seq_len(nrow(tokens))) {
      hits <- grep(tokens$scan_token[[j]], text, perl = TRUE)
      if (!length(hits)) next
      for (line_number in hits) {
        rows[[length(rows) + 1L]] <- data.frame(
          consumer_script = relative[[i]],
          consumed_file = paste0(
            "results/tables/06_modules_WGCNA/group_effects/<dataset>/",
            "{module_group_effects.csv|supermodule_group_effects.csv}"
          ),
          legacy_token = tokens$consumed_column[[j]],
          line_number = as.integer(line_number),
          consumed_column = tokens$consumed_column[[j]],
          current_semantics = tokens$current_semantics[[j]],
          required_new_column = tokens$required_new_column[[j]],
          migration_required = TRUE,
          blocking_for_execution = runtime_consumer,
          proposed_phase3_action = tokens$proposed_phase3_action[[j]],
          consumer_kind = consumer_kind,
          runtime_execution_consumer = runtime_consumer,
          phase3_runtime_migration_required = runtime_consumer,
          audit_record_role = audit_role,
          should_block_execution = runtime_consumer,
          enforcement_status = "advisory_not_runtime_enforced",
          scan_line_numbers = as.character(line_number),
          scan_included_path_rules = rules$included,
          scan_excluded_path_rules = rules$excluded,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  out <- dplyr::bind_rows(rows)
  if (!nrow(out)) {
    stop("Downstream consumer scan found no legacy contract consumers.",
         call. = FALSE)
  }
  out |>
    dplyr::arrange(
      .data$consumer_script, .data$line_number, .data$consumed_column
    )
}

wgcna_group_contract_status <- function(
    dataset, source_paths, canonical_primary_outputs_complete,
    required_source_paths = source_paths
) {
  wgcna_group_assert_stage05_source_dependencies(
    source_paths, required_source_paths
  )
  source_paths <- sort(unique(normalizePath(
    source_paths, winslash = "/", mustWork = TRUE
  )))
  relative <- vapply(source_paths, relative_to, character(1))
  hashes <- vapply(source_paths, file_hash_sha256, character(1))
  if (anyNA(hashes) || any(!nzchar(hashes))) {
    stop("Stage 05 contract status requires complete source SHA-256 hashes.",
         call. = FALSE)
  }
  data.frame(
    dataset = dataset,
    stage05_contract_version = wgcna_group_effects_contract_version(),
    stage05_output_status = "phase2b_statistical_outputs_complete",
    publication_status = "not_assessed_stage05",
    downstream_compatible = FALSE,
    downstream_contract_status = "phase3_migration_required",
    downstream_migration_required = TRUE,
    should_block_execution = TRUE,
    canonical_primary_outputs_complete =
      isTRUE(canonical_primary_outputs_complete),
    generated_at = format(
      Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"
    ),
    source_hashes = paste(
      paste(relative, hashes, sep = "="), collapse = ";"
    ),
    stringsAsFactors = FALSE
  )
}

wgcna_group_stage05_source_dependencies <- function(
    state_path, contract_paths
) {
  source_files <- c(
    repo_path("06_modules_WGCNA", "05_module_supermodule_group_effects.r"),
    repo_path("R", "paths.R"),
    repo_path("R", "dataset_config.R"),
    repo_path("R", "dataset_inputs.R"),
    repo_path("R", "module_contracts.R"),
    repo_path("R", "wgcna_downstream_utils.R"),
    repo_path("R", "wgcna_identity_contract_utils.R"),
    repo_path("R", "wgcna_group_effects_utils.R")
  )
  schema_files <- c(
    "module_group_effects.yml",
    "supermodule_group_effects.yml",
    "wgcna_group_effect_model_validation.yml",
    "wgcna_group_effect_animal_spatial_unit_values.yml",
    "wgcna_group_effect_hemisphere_values.yml",
    "wgcna_group_effect_downstream_consumer_migration_audit.yml",
    "wgcna_group_effect_contract_status.yml"
  )
  sort(unique(c(
    source_files,
    repo_path("inst", "schemas", schema_files),
    state_path,
    unname(unlist(contract_paths))
  )))
}

wgcna_group_assert_stage05_source_dependencies <- function(
    source_paths, required_source_paths
) {
  normalize <- function(x) {
    sort(unique(tolower(normalizePath(
      x, winslash = "/", mustWork = TRUE
    ))))
  }
  missing <- setdiff(normalize(required_source_paths), normalize(source_paths))
  if (length(missing)) {
    stop(
      "Stage 05 source-hash set omits direct dependencies: ",
      paste(missing, collapse = ";"),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

wgcna_group_allowed_output_names <- function() {
  c(
    "module_group_effects.csv",
    "supermodule_group_effects.csv",
    "supermodule_composition.csv",
    "WGCNA_group_effect_endpoint_provenance.csv",
    "WGCNA_group_effect_model_validation.csv",
    "WGCNA_group_effect_sample_inclusion_audit.csv",
    "WGCNA_group_effect_animal_spatial_unit_values.csv",
    "WGCNA_group_effect_hemisphere_values.csv",
    "WGCNA_group_effect_interaction_conditional_followup.csv",
    "WGCNA_group_effect_sensitivity.csv",
    "WGCNA_group_effect_left_right_concordance.csv",
    "WGCNA_group_effect_downstream_consumer_migration_audit.csv",
    "WGCNA_group_effect_contract_status.csv",
    "supermodule_pca_input_audit.csv",
    "supermodule_pca_eigenvalues.csv",
    "supermodule_pca_member_loadings.csv",
    "all_supermodule_eigengene_group_values.csv",
    "all_supermodule_eigengene_spatial_group_values.csv",
    "WGCNA_group_effect_legacy_output_staleness_audit.csv"
  )
}

wgcna_group_legacy_reason <- function(path) {
  name <- tolower(basename(path))
  if (grepl("\\.svg$", name)) {
    return("legacy Stage 05 figure; Phase 2 does not regenerate figures")
  }
  if (grepl("significance|bracket", name)) {
    return("legacy significance/bracket output; not part of canonical Phase 2 inference")
  }
  if (grepl("selected_sus|interpretation", name)) {
    return("legacy selected-supermodule or interpretation output; identity and claim migration deferred")
  }
  if (grepl("label|contents", name)) {
    return("legacy biological-label/content output; Phase 2 is quantitative only")
  }
  if (grepl("marker_trait", name)) {
    return("legacy annotation-only marker-trait correlation; excluded from Phase 2")
  }
  if (grepl("module_to_supermodule", name)) {
    return("legacy historical membership/annotation map; canonical identity contract supersedes it")
  }
  if (grepl("cohend|region_layer|effect", name)) {
    return("legacy auxiliary effect/heatmap source; publication rendering deferred")
  }
  "legacy Stage 05 auxiliary output; not regenerated by canonical Phase 2"
}

wgcna_group_legacy_staleness_audit <- function(dataset, paths) {
  roots <- c(paths$tables, paths$source_data, paths$figures)
  allowed <- wgcna_group_allowed_output_names()
  rows <- list()
  for (root in roots) {
    if (!dir.exists(root)) next
    files <- list.files(root, recursive = TRUE, full.names = TRUE)
    files <- files[file.exists(files) & !dir.exists(files)]
    files <- files[!basename(files) %in% allowed]
    for (path in files) {
      info <- file.info(path)
      rows[[length(rows) + 1L]] <- data.frame(
        dataset = dataset,
        output_scope = if (identical(root, paths$tables)) {
          "tables"
        } else if (identical(root, paths$source_data)) {
          "source_data"
        } else {
          "figures"
        },
        source_file = wgcna_identity_relative_path(path),
        file_name = basename(path),
        size_bytes = as.numeric(info$size),
        sha256 = file_hash_sha256(path),
        modified_time_utc = format(
          info$mtime, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"
        ),
        staleness_reason = wgcna_group_legacy_reason(path),
        preservation_status = "preserved_byte_identical_not_regenerated",
        stringsAsFactors = FALSE
      )
    }
  }
  dplyr::bind_rows(rows) |>
    dplyr::arrange(.data$output_scope, .data$source_file)
}

wgcna_group_copy_tree <- function(source, destination) {
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(source)) return(invisible(destination))
  entries <- list.files(
    source, recursive = TRUE, full.names = TRUE,
    all.files = TRUE, no.. = TRUE, include.dirs = TRUE
  )
  if (!length(entries)) return(invisible(destination))
  info <- file.info(entries)
  relative <- substring(
    normalizePath(entries, winslash = "/", mustWork = TRUE),
    nchar(normalizePath(source, winslash = "/", mustWork = TRUE)) + 2L
  )
  dirs <- entries[info$isdir %in% TRUE]
  dir_rel <- relative[info$isdir %in% TRUE]
  for (i in seq_along(dirs)) {
    dir.create(
      file.path(destination, dir_rel[[i]]),
      recursive = TRUE, showWarnings = FALSE
    )
  }
  files <- entries[!info$isdir]
  file_rel <- relative[!info$isdir]
  for (i in seq_along(files)) {
    target <- file.path(destination, file_rel[[i]])
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(files[[i]], target, overwrite = TRUE, copy.date = TRUE)) {
      stop("Could not stage existing output: ", files[[i]], call. = FALSE)
    }
  }
  invisible(destination)
}

wgcna_group_prepare_stage <- function(paths, dataset, run_id) {
  parent <- dirname(paths$logs)
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  root <- file.path(parent, paste0(".phase2_stage_", dataset, "_", run_id))
  wgcna_group_safe_remove_dir(root, parent)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  staged <- list(
    tables = file.path(root, "tables"),
    source_data = file.path(root, "source_data"),
    logs = file.path(root, "logs")
  )
  wgcna_group_copy_tree(paths$tables, staged$tables)
  wgcna_group_copy_tree(paths$source_data, staged$source_data)
  wgcna_group_copy_tree(paths$logs, staged$logs)
  list(root = root, directories = staged)
}

wgcna_group_safe_remove_dir <- function(path, expected_parent) {
  path_norm <- normalizePath(path, winslash = "/", mustWork = FALSE)
  parent_norm <- normalizePath(expected_parent, winslash = "/", mustWork = TRUE)
  if (!identical(
    tolower(dirname(path_norm)),
    tolower(parent_norm)
  )) {
    stop("Refusing to remove a transaction directory outside its expected parent.",
         call. = FALSE)
  }
  if (dir.exists(path)) unlink(path, recursive = TRUE, force = TRUE)
}

wgcna_group_atomic_publish <- function(
    staged_directories, target_directories, run_id, fail_after = Inf
) {
  if (!identical(names(staged_directories), names(target_directories))) {
    stop("Atomic publication directory roles differ.", call. = FALSE)
  }
  roles <- names(target_directories)
  backups <- stats::setNames(character(length(roles)), roles)
  target_existed <- stats::setNames(logical(length(roles)), roles)
  swapped <- character()
  error <- NULL
  tryCatch({
    for (i in seq_along(roles)) {
      role <- roles[[i]]
      target <- target_directories[[role]]
      stage <- staged_directories[[role]]
      parent <- dirname(target)
      dir.create(parent, recursive = TRUE, showWarnings = FALSE)
      backup <- file.path(parent, paste0(".phase2_backup_", basename(target), "_", run_id))
      backups[[role]] <- backup
      wgcna_group_safe_remove_dir(backup, parent)
      target_existed[[role]] <- dir.exists(target)
      if (target_existed[[role]] && !file.rename(target, backup)) {
        stop("Could not move canonical directory to transaction backup: ", target)
      }
      if (!file.rename(stage, target)) {
        if (target_existed[[role]]) file.rename(backup, target)
        stop("Could not publish staged canonical directory: ", target)
      }
      swapped <- c(swapped, role)
      if (i >= fail_after) stop("Injected atomic-publication failure.")
    }
  }, error = function(e) {
    error <<- e
  })
  if (!is.null(error)) {
    for (role in rev(swapped)) {
      target <- target_directories[[role]]
      parent <- dirname(target)
      wgcna_group_safe_remove_dir(target, parent)
      if (target_existed[[role]] && dir.exists(backups[[role]])) {
        if (!file.rename(backups[[role]], target)) {
          stop("Atomic rollback failed for ", target, call. = FALSE)
        }
      }
    }
    stop(conditionMessage(error), call. = FALSE)
  }
  for (role in roles) {
    if (dir.exists(backups[[role]])) {
      wgcna_group_safe_remove_dir(backups[[role]], dirname(backups[[role]]))
    }
  }
  invisible(TRUE)
}

wgcna_group_diagnostic_scope_violations <- function(
    table, require_valid = TRUE
) {
  allowed_scopes <- c(
    "single_fitted_model",
    "composite_reduced_full",
    "full_interaction_model"
  )
  rows <- list()
  add_violation <- function(i, reason) {
    rows[[length(rows) + 1L]] <<- data.frame(
      row_index = as.integer(i),
      dataset = as.character(table$dataset[[i]]),
      level = as.character(table$level[[i]]),
      endpoint_id = as.character(table$endpoint_id[[i]]),
      test_type = as.character(table$test_type[[i]]),
      model_diagnostic_scope =
        as.character(table$model_diagnostic_scope[[i]]),
      reason = reason,
      stringsAsFactors = FALSE
    )
  }
  scalar_equal <- function(x, y) {
    isTRUE(all.equal(x, y, check.attributes = FALSE, tolerance = 0))
  }
  raw_fields <- c(
    "random_intercept_variance", "residual_variance",
    "variance_ratio", "ICC", "is_singular_lme4",
    "boundary_by_variance_ratio"
  )
  for (i in seq_len(nrow(table))) {
    scope <- table$model_diagnostic_scope[[i]]
    if (!scope %in% allowed_scopes) {
      add_violation(i, "invalid_model_diagnostic_scope")
      next
    }
    expected_scope <- dplyr::case_when(
      identical(table$test_type[[i]], "interaction_omnibus") ~
        "composite_reduced_full",
      identical(
        table$test_type[[i]], "conditional_interaction_followup"
      ) ~ "full_interaction_model",
      TRUE ~ "single_fitted_model"
    )
    if (!identical(scope, expected_scope)) {
      add_violation(i, "scope_does_not_match_test_type")
    }
    if (isTRUE(require_valid) &&
        !isTRUE(table$model_valid_for_inference[[i]])) {
      add_violation(i, "output_row_not_valid_for_inference")
    }
    if (identical(scope, "composite_reduced_full")) {
      if (any(!is.na(unlist(table[i, raw_fields, drop = FALSE])))) {
        add_violation(i, "composite_top_level_raw_diagnostics_not_na")
      }
      complete_character <- c(
        "reduced_formula", "full_formula",
        "reduced_row_hash", "full_row_hash",
        "reduced_optimizer_code", "full_optimizer_code",
        "reduced_convergence_status", "full_convergence_status",
        "reduced_random_effect_structure",
        "full_random_effect_structure",
        "reduced_singularity_class", "full_singularity_class",
        "reduced_model_stability_status",
        "full_model_stability_status"
      )
      complete_numeric <- c(
        "reduced_fixed_effect_rank", "full_fixed_effect_rank",
        "reduced_fixed_effect_columns", "full_fixed_effect_columns",
        "reduced_random_intercept_variance",
        "full_random_intercept_variance",
        "reduced_residual_variance", "full_residual_variance",
        "reduced_variance_ratio", "full_variance_ratio",
        "reduced_ICC", "full_ICC"
      )
      complete_logical <- c(
        "reduced_is_singular_lme4", "full_is_singular_lme4",
        "reduced_boundary_by_variance_ratio",
        "full_boundary_by_variance_ratio",
        "reduced_diagnostic_review_required",
        "full_diagnostic_review_required"
      )
      if (any(is.na(unlist(
        table[i, complete_character, drop = FALSE]
      ))) || any(!nzchar(as.character(unlist(
        table[i, complete_character, drop = FALSE]
      ))))) {
        add_violation(i, "composite_character_diagnostics_incomplete")
      }
      if (any(!is.finite(as.numeric(unlist(
        table[i, complete_numeric, drop = FALSE]
      ))))) {
        add_violation(i, "composite_numeric_diagnostics_incomplete")
      }
      if (any(is.na(unlist(
        table[i, complete_logical, drop = FALSE]
      ))) || !isTRUE(table$identical_rows_verified[[i]])) {
        add_violation(i, "composite_logical_or_row_contract_incomplete")
      }
      reduced_status <- table$reduced_model_stability_status[[i]]
      full_status <- table$full_model_stability_status[[i]]
      valid_models <- all(c(reduced_status, full_status) %in% c(
        "stable_mixed_model", "boundary_random_intercept_zero"
      ))
      boundary <- valid_models && any(c(
        reduced_status, full_status
      ) == "boundary_random_intercept_zero")
      review <- valid_models && (
        isTRUE(table$reduced_diagnostic_review_required[[i]]) ||
          isTRUE(table$full_diagnostic_review_required[[i]])
      )
      expected_status <- if (!valid_models) {
        "invalid"
      } else if (boundary) {
        "boundary_random_intercept_zero"
      } else {
        "stable_mixed_model"
      }
      expected_class <- if (!valid_models) {
        "invalid_composite_reduced_full"
      } else if (boundary) {
        "composite_reduced_full_boundary"
      } else if (review) {
        "composite_reduced_full_diagnostic_review"
      } else {
        "composite_reduced_full_stable"
      }
      if (!identical(
        table$model_stability_status[[i]], expected_status
      ) || !identical(
        table$singularity_class[[i]], expected_class
      ) || !identical(
        table$primary_model_stable[[i]],
        valid_models && !boundary && !review
      ) || !identical(
        table$diagnostic_review_required[[i]], review
      )) {
        add_violation(i, "composite_status_mismatch")
      }
    } else if (identical(scope, "full_interaction_model")) {
      mapping <- c(
        random_intercept_variance =
          "full_random_intercept_variance",
        residual_variance = "full_residual_variance",
        variance_ratio = "full_variance_ratio",
        ICC = "full_ICC",
        is_singular_lme4 = "full_is_singular_lme4",
        boundary_by_variance_ratio =
          "full_boundary_by_variance_ratio",
        singularity_class = "full_singularity_class",
        model_stability_status = "full_model_stability_status",
        diagnostic_review_required =
          "full_diagnostic_review_required"
      )
      mismatched <- vapply(names(mapping), function(top) {
        !scalar_equal(table[[top]][[i]], table[[mapping[[top]]]][[i]])
      }, logical(1))
      if (any(mismatched)) {
        add_violation(i, "full_interaction_top_level_mismatch")
      }
    } else if (isTRUE(table$model_valid_for_inference[[i]]) &&
               identical(table$model_family[[i]], "linear_mixed_model")) {
      if (any(!is.finite(as.numeric(unlist(
        table[i, raw_fields[1:4], drop = FALSE]
      )))) ||
          is.na(table$is_singular_lme4[[i]]) ||
          is.na(table$boundary_by_variance_ratio[[i]])) {
        add_violation(i, "single_fit_diagnostics_incomplete")
      }
      expected_status <- if (
        isTRUE(table$boundary_by_variance_ratio[[i]])
      ) {
        "boundary_random_intercept_zero"
      } else {
        "stable_mixed_model"
      }
      if (!identical(
        table$model_stability_status[[i]], expected_status
      )) {
        add_violation(i, "single_fit_status_mismatch")
      }
    }
  }
  dplyr::bind_rows(rows)
}

wgcna_group_model_validation_scope_violations <- function(table) {
  rows <- list()
  add_violation <- function(i, reason) {
    rows[[length(rows) + 1L]] <<- data.frame(
      row_index = as.integer(i),
      dataset = as.character(table$dataset[[i]]),
      level = as.character(table$level[[i]]),
      endpoint_id = as.character(table$endpoint_id[[i]]),
      model_diagnostic_scope =
        as.character(table$model_diagnostic_scope[[i]]),
      reason = reason,
      stringsAsFactors = FALSE
    )
  }
  raw_fields <- c(
    "random_intercept_variance", "residual_variance",
    "variance_ratio", "ICC", "is_singular_lme4",
    "boundary_by_variance_ratio"
  )
  for (i in seq_len(nrow(table))) {
    expected_scope <- if (identical(
      table$effect_scope[[i]], "stress_by_spatial_interaction"
    )) {
      "composite_reduced_full"
    } else {
      "single_fitted_model"
    }
    if (!identical(table$model_diagnostic_scope[[i]], expected_scope)) {
      add_violation(i, "validation_scope_mismatch")
      next
    }
    if (identical(expected_scope, "composite_reduced_full")) {
      if (any(!is.na(unlist(table[i, raw_fields, drop = FALSE])))) {
        add_violation(i, "validation_composite_raw_diagnostics_not_na")
      }
      if (!isTRUE(table$model_valid_for_inference[[i]])) {
        if (!identical(
          table$singularity_class[[i]],
          "invalid_composite_reduced_full"
        )) {
          add_violation(i, "validation_invalid_composite_class_mismatch")
        }
        next
      }
      reduced_status <- table$reduced_model_stability_status[[i]]
      full_status <- table$full_model_stability_status[[i]]
      valid_models <- all(c(reduced_status, full_status) %in% c(
        "stable_mixed_model", "boundary_random_intercept_zero"
      ))
      boundary <- valid_models && any(c(
        reduced_status, full_status
      ) == "boundary_random_intercept_zero")
      review <- valid_models && (
        isTRUE(table$reduced_diagnostic_review_required[[i]]) ||
          isTRUE(table$full_diagnostic_review_required[[i]])
      )
      expected_status <- if (boundary) {
        "boundary_random_intercept_zero"
      } else {
        "stable_mixed_model"
      }
      expected_class <- if (boundary) {
        "composite_reduced_full_boundary"
      } else if (review) {
        "composite_reduced_full_diagnostic_review"
      } else {
        "composite_reduced_full_stable"
      }
      if (!valid_models ||
          !isTRUE(table$identical_rows_verified[[i]]) ||
          !identical(table$model_stability_status[[i]], expected_status) ||
          !identical(table$singularity_class[[i]], expected_class) ||
          !identical(
            table$primary_model_stable[[i]], !boundary && !review
          ) ||
          !identical(table$diagnostic_review_required[[i]], review)) {
        add_violation(i, "validation_composite_status_mismatch")
      }
    } else if (isTRUE(table$model_valid_for_inference[[i]]) &&
               identical(table$model_type[[i]], "lmerTest_lmer")) {
      if (any(!is.finite(as.numeric(unlist(
        table[i, raw_fields[1:4], drop = FALSE]
      )))) ||
          is.na(table$is_singular_lme4[[i]]) ||
          is.na(table$boundary_by_variance_ratio[[i]])) {
        add_violation(i, "validation_single_fit_diagnostics_incomplete")
      }
    }
  }
  dplyr::bind_rows(rows)
}

wgcna_group_ci_p_compatibility_violations <- function(
    table, alpha = 0.05, tolerance = 1e-10
) {
  if (!nrow(table)) return(data.frame())
  required <- c(
    "dataset", "level", "endpoint_id", "effect_scope", "spatial_unit",
    "contrast", "test_type", "estimate", "SE", "CI_low", "CI_high",
    "CI_method", "CI_level", "CI_df_method", "p_value", "df_den"
  )
  missing <- setdiff(required, names(table))
  if (length(missing)) {
    stop(
      "CI compatibility check is missing columns: ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  named <- table[
    table$test_type %in% c(
      "named_contrast", "conditional_interaction_followup"
    ) &
      is.finite(table$estimate) & is.finite(table$SE) &
      is.finite(table$CI_low) & is.finite(table$CI_high) &
      is.finite(table$p_value),
    , drop = FALSE
  ]
  if (!nrow(named)) return(data.frame())
  ci_excludes <- named$CI_low > tolerance | named$CI_high < -tolerance
  p_supported <- named$p_value < alpha - tolerance
  boundary <- abs(named$p_value - alpha) <= tolerance |
    abs(named$CI_low) <= tolerance | abs(named$CI_high) <= tolerance
  method_ok <- named$CI_method == "emmeans_contrast_summary" &
    abs(named$CI_level - (1 - alpha)) <= tolerance &
    named$CI_df_method %in% c(
      "finite_df_from_emmeans", "asymptotic_from_emmeans"
    )
  interval_ok <- named$CI_low <= named$estimate + tolerance &
    named$CI_high >= named$estimate - tolerance
  incompatible <- (!boundary & ci_excludes != p_supported) |
    !method_ok | !interval_ok
  out <- named[incompatible, required, drop = FALSE]
  if (nrow(out)) {
    out$CI_excludes_zero <- ci_excludes[incompatible]
    out$p_below_alpha <- p_supported[incompatible]
    out$boundary_edge_case <- boundary[incompatible]
  }
  rownames(out) <- NULL
  out
}

wgcna_group_validate_output_bundle <- function(outputs, contract) {
  required <- wgcna_group_allowed_output_names()
  missing <- setdiff(required, names(outputs))
  if (length(missing)) {
    stop("Canonical output bundle is missing: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  module <- outputs[["module_group_effects.csv"]]
  super <- outputs[["supermodule_group_effects.csv"]]
  followup <-
    outputs[["WGCNA_group_effect_interaction_conditional_followup.csv"]]
  model_validation <-
    outputs[["WGCNA_group_effect_model_validation.csv"]]
  composition <- outputs[["supermodule_composition.csv"]]
  for (table in list(module, super)) {
    named <- table$test_type == "named_contrast"
    omnibus <- table$test_type == "interaction_omnibus"
    if (nrow(table) && (
      any(!is.finite(table$p_value)) ||
      any(!table$model_valid_for_inference) ||
      any(!table$claim_allowed_model) ||
      any(table$fallback_used) ||
      any(named & (
        !is.finite(table$estimate) |
          !is.finite(table$SE) |
          !is.finite(table$statistic)
      )) ||
      any(omnibus & (
        !is.finite(table$likelihood_ratio_statistic) |
          !is.finite(table$likelihood_ratio_df) |
          !is.finite(table$likelihood_ratio_p_value) |
          !(table$identical_rows_verified %in% TRUE) |
          is.na(table$reduced_fixed_effect_rank) |
          is.na(table$full_fixed_effect_rank) |
          is.na(table$reduced_model_stability_status) |
          is.na(table$full_model_stability_status)
      ))
    )) {
      stop("Primary group-effect output contains invalid or diagnostic rows.",
           call. = FALSE)
    }
  }
  diagnostic_violations <- dplyr::bind_rows(
    wgcna_group_diagnostic_scope_violations(module),
    wgcna_group_diagnostic_scope_violations(super),
    wgcna_group_diagnostic_scope_violations(followup)
  )
  if (nrow(diagnostic_violations)) {
    stop(
      "Model diagnostic-scope contract failed: ",
      paste(unique(diagnostic_violations$reason), collapse = ";"),
      call. = FALSE
    )
  }
  ci_violations <- dplyr::bind_rows(
    wgcna_group_ci_p_compatibility_violations(module),
    wgcna_group_ci_p_compatibility_violations(super),
    wgcna_group_ci_p_compatibility_violations(followup)
  )
  if (nrow(ci_violations)) {
    stop(
      "Confidence-interval and p-value inference contract failed.",
      call. = FALSE
    )
  }
  model_validation_violations <-
    wgcna_group_model_validation_scope_violations(model_validation)
  if (nrow(model_validation_violations)) {
    stop(
      "Model-validation diagnostic-scope contract failed: ",
      paste(
        unique(model_validation_violations$reason), collapse = ";"
      ),
      call. = FALSE
    )
  }
  primary_rows <- dplyr::bind_rows(module, super) |>
    dplyr::filter(.data$independent_hypothesis %in% TRUE)
  if (any(
    primary_rows$analysis_tier == "primary_wgcna_global" &
      (!is.na(primary_rows$FDR_secondary_global) |
         primary_rows$contrast != "SUS - RES")
  ) || any(
    primary_rows$analysis_tier == "secondary_contextual_global" &
      (!is.na(primary_rows$FDR_primary_global) |
         !primary_rows$contrast %in% c("RES - CON", "SUS - CON"))
  )) {
    stop("Primary and contextual global FDR fields are not mutually exclusive.",
         call. = FALSE)
  }
  conservative_rows <- !is.na(primary_rows$FDR_conservative_all_tests)
  if (!isTRUE(all.equal(
    primary_rows$FDR_global[conservative_rows],
    primary_rows$FDR_conservative_all_tests[conservative_rows],
    check.attributes = FALSE
  ))) {
    stop("FDR_global is not an exact deprecated conservative-FDR alias.",
         call. = FALSE)
  }
  aliases <- super[super$hypothesis_level == "compatibility_alias", ,
                   drop = FALSE]
  alias_fdr <- c(
    "FDR_primary_global", "FDR_secondary_global",
    "FDR_interaction_omnibus", "FDR_local_exploratory",
    "FDR_conservative_all_tests", "FDR_global",
    "FDR_within_dataset_level", "FDR_dataset_all_levels"
  )
  if (nrow(aliases) && (
    any(aliases$independent_hypothesis) ||
      any(aliases$enters_primary_fdr) ||
      any(aliases$statistical_support_status !=
            "inherited_from_canonical_entity") ||
      any(vapply(alias_fdr, function(nm) {
        any(!is.na(aliases[[nm]]))
      }, logical(1)))
  )) {
    stop("Singleton compatibility-alias inference contract failed.",
         call. = FALSE)
  }
  inherited_cols <- c(
    "estimate", "SE", "CI_low", "CI_high",
    "CI_method", "CI_level", "CI_df_method", "p_value",
    "model_formula", "formula_used", "random_intercept_variance",
    "residual_variance", "variance_ratio", "ICC",
    "is_singular_lme4", "boundary_by_variance_ratio",
    "boundary_variance_ratio_tolerance", "lme4_singularity_tolerance",
    "singularity_diagnostic_status", "diagnostic_review_required",
    "singularity_class", "model_valid_for_inference",
    "model_diagnostic_scope",
    "model_stability_status", "claim_allowed_model",
    "identity_contract_sha256", "frozen_state_sha256"
  )
  if (nrow(aliases)) {
    source_key <- paste(
      module$module_id, module$effect_scope, module$spatial_unit,
      module$contrast, module$test_type, sep = "|"
    )
    alias_source_key <- paste(
      aliases$support_source_entity_id, aliases$effect_scope,
      aliases$spatial_unit, aliases$contrast, aliases$test_type, sep = "|"
    )
    idx <- match(alias_source_key, source_key)
    if (anyNA(idx)) {
      stop("Singleton aliases do not resolve to canonical module rows.",
           call. = FALSE)
    }
    for (nm in inherited_cols) {
      if (!isTRUE(all.equal(
        aliases[[nm]], module[[nm]][idx], check.attributes = FALSE
      ))) {
        stop("Singleton alias failed to inherit canonical field: ", nm,
             call. = FALSE)
      }
    }
  }
  status <- outputs[["WGCNA_group_effect_contract_status.csv"]]
  if (nrow(status) != 1L ||
      status$stage05_output_status[[1]] !=
        "phase2b_statistical_outputs_complete" ||
      status$publication_status[[1]] != "not_assessed_stage05" ||
      status$downstream_compatible[[1]] %in% TRUE ||
      status$downstream_contract_status[[1]] !=
        "phase3_migration_required" ||
      !(status$downstream_migration_required[[1]] %in% TRUE) ||
      !(status$should_block_execution[[1]] %in% TRUE) ||
      !(status$canonical_primary_outputs_complete[[1]] %in% TRUE)) {
    stop("Stage 05 contract-status row is incomplete or semantically invalid.",
         call. = FALSE)
  }
  consumer_audit <-
    outputs[["WGCNA_group_effect_downstream_consumer_migration_audit.csv"]]
  required_audit_cols <- c(
    "consumer_script", "consumed_file", "consumed_column",
    "current_semantics", "required_new_column", "migration_required",
    "blocking_for_execution", "proposed_phase3_action", "legacy_token",
    "line_number", "consumer_kind", "runtime_execution_consumer",
    "phase3_runtime_migration_required", "audit_record_role"
  )
  if (!all(required_audit_cols %in% names(consumer_audit)) ||
      !nrow(consumer_audit) ||
      any(!consumer_audit$migration_required) ||
      any(consumer_audit$blocking_for_execution !=
            consumer_audit$runtime_execution_consumer) ||
      any(consumer_audit$phase3_runtime_migration_required !=
            consumer_audit$runtime_execution_consumer) ||
      any(consumer_audit$enforcement_status !=
            "advisory_not_runtime_enforced")) {
    stop("Downstream-consumer migration audit is incomplete.", call. = FALSE)
  }
  aggregates <-
    outputs[["WGCNA_group_effect_animal_spatial_unit_values.csv"]]
  hemispheres <- outputs[["WGCNA_group_effect_hemisphere_values.csv"]]
  if (!nrow(aggregates) || !nrow(hemispheres) ||
      anyDuplicated(aggregates$aggregated_row_id) ||
      anyDuplicated(aggregates$aggregated_row_sha256) ||
      any(!grepl("^[0-9a-f]{64}$", aggregates$aggregated_row_sha256)) ||
      any(aggregates$aggregated_row_hash_contract != "sha256_utf8_lf_v1") ||
      any(vapply(seq_len(nrow(aggregates)), function(i) {
        !identical(
          aggregates$aggregated_row_sha256[[i]],
          wgcna_group_aggregate_row_sha256(
            aggregates[i, , drop = FALSE]
          )
        )
      }, logical(1))) ||
      any(!aggregates$n_hemispheres_observed %in% 1:2) ||
      any(!is.finite(aggregates$eigengene)) ||
      any(!is.finite(hemispheres$hemisphere_value))) {
    stop("Hemisphere or animal-spatial aggregation contract failed.",
         call. = FALSE)
  }
  observed <- composition |>
    dplyr::select("dataset", "supermodule_id", "member_modules") |>
    tidyr::separate_longer_delim("member_modules", delim = ";") |>
    dplyr::rename(module_id = "member_modules") |>
    dplyr::arrange(.data$dataset, .data$module_id, .data$supermodule_id) |>
    dplyr::select("dataset", "module_id", "supermodule_id")
  expected <- contract$membership |>
    dplyr::arrange(.data$dataset, .data$module_id, .data$supermodule_id) |>
    dplyr::select("dataset", "module_id", "supermodule_id")
  if (!identical(as.data.frame(observed), as.data.frame(expected))) {
    stop("Output composition differs from canonical membership.", call. = FALSE)
  }
  invisible(TRUE)
}

wgcna_group_build_bundle <- function(dataset, state, contract) {
  wgcna_group_require_primary_packages()
  bridge <- wgcna_group_build_module_bridge(dataset, state, contract)
  endpoints <- wgcna_group_construct_endpoints(dataset, state, contract, bridge)
  sample_audit <- wgcna_group_sample_audit(
    dataset, state, endpoints$module_eigengenes
  )
  module <- wgcna_group_run_level(
    dataset, "module", endpoints$module_eigengenes,
    endpoints$module_map, sample_audit, contract
  )
  supermodule <- wgcna_group_run_level(
    dataset, "supermodule", endpoints$supermodule_eigengenes,
    endpoints$supermodule_map, sample_audit, contract
  )
  adjusted <- wgcna_group_apply_fdr(module$primary, supermodule$primary)
  singleton_aliases <- wgcna_group_make_singleton_alias_rows(
    adjusted$module, endpoints$composition
  )
  supermodule_effects <- dplyr::bind_rows(
    adjusted$supermodule, singleton_aliases
  ) |>
    dplyr::arrange(.data$endpoint_id, .data$effect_scope,
                   .data$spatial_unit, .data$contrast)
  values <- wgcna_group_supermodule_values(
    dataset, endpoints$supermodule_eigengenes,
    endpoints$supermodule_map, sample_audit, contract
  )
  list(
    module_group_effects = adjusted$module,
    supermodule_group_effects = supermodule_effects,
    supermodule_composition = endpoints$composition,
    endpoint_provenance = endpoints$provenance,
    model_validation = dplyr::bind_rows(
      module$validation, supermodule$validation
    ),
    sample_inclusion_audit = sample_audit,
    animal_spatial_unit_values = dplyr::bind_rows(
      module$animal_spatial_values,
      supermodule$animal_spatial_values
    ),
    hemisphere_values = dplyr::bind_rows(
      module$hemisphere_values,
      supermodule$hemisphere_values
    ),
    interaction_conditional_followup = dplyr::bind_rows(
      module$interaction_followup,
      supermodule$interaction_followup
    ),
    sensitivity = dplyr::bind_rows(
      module$sensitivity, supermodule$sensitivity
    ),
    left_right_concordance = dplyr::bind_rows(
      module$left_right_concordance,
      supermodule$left_right_concordance
    ),
    supermodule_pca_input_audit = endpoints$pca_input_audit,
    supermodule_pca_eigenvalues = endpoints$pca_eigenvalues,
    supermodule_pca_member_loadings = endpoints$pca_member_loadings,
    all_supermodule_eigengene_group_values = values,
    all_supermodule_eigengene_spatial_group_values = values
  )
}

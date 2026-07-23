# Canonical Phase 2 WGCNA group-effect modelling utilities.
#
# This layer is deliberately independent of biological labels and historical
# supermodule maps. The Phase 1 identity contract is the sole supermodule
# membership authority.

wgcna_group_effects_contract_version <- function() {
  "wgcna_group_effects_v2"
}

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

wgcna_group_build_module_bridge <- function(dataset, state, contract) {
  module_eigengenes <- extract_module_eigengenes(state)
  eigengene_cols <- setdiff(names(module_eigengenes), "Sample")
  contract_modules <- sort(unique(as.character(contract$membership$module_id)))

  if (dataset %in% c("neuron_soma", "neuron_neuropil")) {
    bridge <- wgcna_identity_bridge_module_ids(eigengene_cols, contract_modules)
    out <- data.frame(
      dataset = dataset,
      module_id = bridge$normalized_id,
      state_eigengene_col_raw = bridge$original_id,
      state_eigengene_col_normalized = bridge$normalized_id,
      bridge_method = bridge$bridge_type,
      stringsAsFactors = FALSE
    )
  } else {
    labels <- as.data.frame(
      state$module_label_table,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    required <- c("ModuleID", "WGCNAInternalColor")
    if (!all(required %in% names(labels))) {
      stop(
        "Microglia frozen state lacks ModuleID/WGCNAInternalColor metadata.",
        call. = FALSE
      )
    }
    labels <- labels[labels$ModuleID %in% contract_modules, required, drop = FALSE]
    labels$state_eigengene_col_raw <- paste0("ME", labels$WGCNAInternalColor)
    out <- data.frame(
      dataset = dataset,
      module_id = as.character(labels$ModuleID),
      state_eigengene_col_raw = as.character(labels$state_eigengene_col_raw),
      state_eigengene_col_normalized = as.character(labels$ModuleID),
      bridge_method = "stable_state_metadata_bridge",
      stringsAsFactors = FALSE
    )
  }
  out <- out[order(out$module_id), , drop = FALSE]
  if (nrow(out) != length(contract_modules) ||
      anyDuplicated(out$module_id) ||
      anyDuplicated(out$state_eigengene_col_raw) ||
      !identical(sort(out$module_id), contract_modules) ||
      !identical(sort(out$state_eigengene_col_raw), sort(eigengene_cols)) ||
      any(!out$state_eigengene_col_raw %in% eigengene_cols)) {
    stop(
      "Frozen-state module eigengene bridge is incomplete, duplicated, or ambiguous for ",
      dataset, ".",
      call. = FALSE
    )
  }
  if (dataset != "microglia" &&
      any(!out$bridge_method %in%
          c("exact_current_identity_match", "syntactic_module_id_bridge_only"))) {
    stop("Neuronal module bridge uses a prohibited mapping.", call. = FALSE)
  }
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
      "endpoint_provenance_status"
    )],
    supermodule_map = audits$provenance[, c(
      "endpoint_col", "endpoint_id", "endpoint_construction_method",
      "endpoint_provenance_status"
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
    supermodule_label = character(), spatial_unit = character(),
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
    statistic = numeric(), p_value = numeric(),
    FDR_within_dataset_level = numeric(),
    FDR_dataset_all_levels = numeric(), FDR_global = numeric(),
    FDR_family_within_level_id = character(),
    FDR_family_dataset_id = character(),
    n_tests_FDR_within_dataset_level = integer(),
    n_tests_FDR_dataset_all_levels = integer(), FDR_method = character(),
    evidence_status = character(), direction = character(),
    n_samples = integer(), formula_requested = character(),
    formula_used = character(), dropped_covariates = character(),
    model_family = character(), model_formula = character(),
    primary_model_stable = logical(), claim_allowed_model = logical(),
    model_downgrade_reason = character(), fallback_used = logical(),
    fallback_type = character(), rank_deficient_model = logical(),
    singular_model = logical(), emmeans_success = logical(),
    animal_random_effect_used = logical(),
    biological_replicate_unit = character(), model_warning = character(),
    fixed_effect_rank = integer(), fixed_effect_columns = integer(),
    convergence_status = character(), convergence_warnings = character(),
    optimizer_messages = character(), optimizer_code = character(),
    residual_df = numeric(), AnimalID_source = character(),
    animal_id_mapping_status = character(),
    membership_version = character(),
    identity_contract_version = character(),
    identity_contract_sha256 = character(),
    frozen_state_sha256 = character(),
    endpoint_construction_method = character(),
    endpoint_provenance_status = character(),
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

wgcna_group_fit_attempt <- function(
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
  singular <- if (repeated) lme4::isSingular(fit) else FALSE
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
      summary(
        emmeans::contrast(
          emm_capture$value, method = methods, adjust = "none"
        )
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

wgcna_group_run_level <- function(
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

wgcna_group_apply_fdr <- function(module_rows, supermodule_rows) {
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

wgcna_group_allowed_output_names <- function() {
  c(
    "module_group_effects.csv",
    "supermodule_group_effects.csv",
    "supermodule_composition.csv",
    "WGCNA_group_effect_endpoint_provenance.csv",
    "WGCNA_group_effect_model_validation.csv",
    "WGCNA_group_effect_sample_inclusion_audit.csv",
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

wgcna_group_validate_output_bundle <- function(outputs, contract) {
  required <- wgcna_group_allowed_output_names()
  missing <- setdiff(required, names(outputs))
  if (length(missing)) {
    stop("Canonical output bundle is missing: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  module <- outputs[["module_group_effects.csv"]]
  super <- outputs[["supermodule_group_effects.csv"]]
  composition <- outputs[["supermodule_composition.csv"]]
  for (table in list(module, super)) {
    if (nrow(table) && (
      any(!is.finite(table$estimate)) ||
      any(!is.finite(table$SE)) ||
      any(!is.finite(table$statistic)) ||
      any(!is.finite(table$p_value)) ||
      any(!table$claim_allowed_model) ||
      any(table$fallback_used)
    )) {
      stop("Primary group-effect output contains invalid or diagnostic rows.",
           call. = FALSE)
    }
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
  values <- wgcna_group_supermodule_values(
    dataset, endpoints$supermodule_eigengenes,
    endpoints$supermodule_map, sample_audit, contract
  )
  list(
    module_group_effects = adjusted$module,
    supermodule_group_effects = adjusted$supermodule,
    supermodule_composition = endpoints$composition,
    endpoint_provenance = endpoints$provenance,
    model_validation = dplyr::bind_rows(
      module$validation, supermodule$validation
    ),
    sample_inclusion_audit = sample_audit,
    supermodule_pca_input_audit = endpoints$pca_input_audit,
    supermodule_pca_eigenvalues = endpoints$pca_eigenvalues,
    supermodule_pca_member_loadings = endpoints$pca_member_loadings,
    all_supermodule_eigengene_group_values = values,
    all_supermodule_eigengene_spatial_group_values = values
  )
}

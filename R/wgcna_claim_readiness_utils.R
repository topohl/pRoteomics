# Shared loader and validator for the finalized microglia WGCNA Stage 13
# manuscript-readiness handoff.

if (!exists("repo_path", mode = "function")) {
  paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
  source(paths_file)
}
source(repo_path("R", "schema_validation.R"))

microglia_wgcna_claim_readiness_path <- function() {
  path_results(
    "tables", "06_modules_WGCNA", "claim_readiness", "microglia",
    "WGCNA_entity_claim_readiness.csv"
  )
}

wgcna_claim_readiness_alias_map <- function() {
  data.frame(
    entity_id = c("SM02", "SM04", "SM05", "SM06", "SM07", "SM08"),
    canonical_claim_entity_id = c("WGCNA_m03", "WGCNA_m09", "WGCNA_m13", "WGCNA_m01", "WGCNA_m02", "WGCNA_m06"),
    stringsAsFactors = FALSE
  )
}

assert_unique_wgcna_keys <- function(df, keys, artifact) {
  missing <- setdiff(keys, names(df))
  if (length(missing)) {
    stop(artifact, " missing stable key column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  key <- do.call(paste, c(lapply(df[keys], as.character), sep = "\r"))
  duplicated_keys <- unique(key[duplicated(key) | duplicated(key, fromLast = TRUE)])
  if (length(duplicated_keys)) {
    stop(
      artifact, " contains duplicate stable key(s) for ", paste(keys, collapse = " + "),
      "; examples: ", paste(utils::head(gsub("\r", " | ", duplicated_keys, fixed = TRUE), 8L), collapse = "; "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

load_microglia_wgcna_claim_readiness <- function(path = microglia_wgcna_claim_readiness_path()) {
  if (!file.exists(path)) {
    stop("Required microglia Stage 13 claim-readiness table is missing: ", path, call. = FALSE)
  }
  readiness <- if (requireNamespace("readr", quietly = TRUE)) {
    readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  } else {
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  }

  validate_table_schema(readiness, "wgcna_entity_claim_readiness", strict = TRUE)
  if (nrow(readiness) != 22L) {
    stop("Stage 13 must contain exactly 22 technical identities; observed ", nrow(readiness), ".", call. = FALSE)
  }
  assert_unique_wgcna_keys(readiness, c("dataset", "level", "entity_id"), "Stage 13 claim readiness")
  if (!identical(unique(as.character(readiness$dataset)), "microglia")) {
    stop("Stage 13 downstream contract is restricted to dataset=microglia.", call. = FALSE)
  }

  roles <- table(as.character(readiness$claim_entity_role))
  expected_roles <- c(canonical_module = 13L, higher_order_block = 3L, compatibility_alias = 6L)
  observed_roles <- stats::setNames(rep(0L, length(expected_roles)), names(expected_roles))
  observed_roles[names(roles)] <- as.integer(roles)
  if (!identical(observed_roles, expected_roles)) {
    stop(
      "Stage 13 role counts must be 13 canonical modules, 3 higher-order blocks, and 6 compatibility aliases; observed ",
      paste(names(observed_roles), observed_roles, sep = "=", collapse = ", "), ".",
      call. = FALSE
    )
  }

  modules <- readiness[readiness$claim_entity_role == "canonical_module", , drop = FALSE]
  blocks <- readiness[readiness$claim_entity_role == "higher_order_block", , drop = FALSE]
  aliases <- readiness[readiness$claim_entity_role == "compatibility_alias", , drop = FALSE]
  if (!all(modules$level == "module") || !setequal(modules$entity_id, sprintf("WGCNA_m%02d", 1:13))) {
    stop("Stage 13 canonical-module identity set is invalid.", call. = FALSE)
  }
  if (!all(blocks$level == "supermodule") || !setequal(blocks$entity_id, c("SM01", "SM03", "SM09"))) {
    stop("Stage 13 higher-order blocks must be exactly SM01, SM03, and SM09.", call. = FALSE)
  }

  expected_aliases <- wgcna_claim_readiness_alias_map()
  observed_aliases <- as.data.frame(
    aliases[match(expected_aliases$entity_id, aliases$entity_id), c("entity_id", "canonical_claim_entity_id"), drop = FALSE],
    stringsAsFactors = FALSE
  )
  rownames(observed_aliases) <- NULL
  if (anyNA(observed_aliases$entity_id) || !identical(observed_aliases, expected_aliases)) {
    stop("Stage 13 compatibility aliases do not match the finalized six exact mappings.", call. = FALSE)
  }
  if (any(aliases$separate_manuscript_claim_allowed %in% TRUE) ||
      any(as.character(aliases$compatibility_target_level) != "module") ||
      any(as.character(aliases$compatibility_target_entity_id) != as.character(aliases$canonical_claim_entity_id))) {
    stop("Stage 13 compatibility aliases must be module aliases and never separately claimable.", call. = FALSE)
  }

  block_allowed <- stats::setNames(as.logical(blocks$separate_manuscript_claim_allowed), blocks$entity_id)
  if (!identical(names(block_allowed)[block_allowed %in% TRUE], "SM09") ||
      any(block_allowed[c("SM01", "SM03")] %in% TRUE)) {
    stop("SM09 must be the only separately claimable Stage 13 higher-order block.", call. = FALSE)
  }
  if (sum(as.logical(modules$separate_manuscript_claim_allowed), na.rm = TRUE) != 10L) {
    stop("Stage 13 must contain exactly ten independently architecture-eligible canonical modules.", call. = FALSE)
  }
  if (any(as.character(readiness$group_effect_status) == "FDR_supported")) {
    stop("Finalized Stage 13 contract must not contain an FDR-supported module or supermodule group effect.", call. = FALSE)
  }
  versions <- unique(as.character(readiness$readiness_contract_version))
  if (length(versions) != 1L || is.na(versions) || !nzchar(versions)) {
    stop("Stage 13 must contain exactly one non-empty readiness_contract_version.", call. = FALSE)
  }

  source_path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  attr(readiness, "stage13_source_file") <- source_path
  list(
    all = readiness,
    independently_claimable = readiness[as.logical(readiness$separate_manuscript_claim_allowed), , drop = FALSE],
    canonical_modules = modules,
    higher_order_blocks = blocks,
    compatibility_aliases = aliases,
    source_path = source_path,
    readiness_contract_version = versions[[1]]
  )
}

stage13_join_columns <- function() {
  c(
    "dataset", "level", "entity_id", "canonical_claim_entity_id", "claim_entity_role",
    "separate_manuscript_claim_allowed", "primary_architecture_status",
    "spatial_dependence_class", "animal_stability_status", "group_effect_status",
    "allowed_claim_scope", "prohibited_claim_scope", "manuscript_placement",
    "readiness_contract_version"
  )
}

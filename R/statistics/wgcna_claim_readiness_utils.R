# Shared loader and validator for the finalized microglia WGCNA Stage 13
# manuscript-readiness handoff.

if (!exists("repo_path", mode = "function")) {
  paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
  source(paths_file)
}
source(repo_path("R", "schema_validation.R"))
if (!exists("wgcna_dirs", mode = "function")) {
  if (!exists("repo_path", mode = "function")) {
    paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
    source(paths_file)
  }
  source(repo_path("R", "wgcna_paths.R"))
}

microglia_wgcna_claim_readiness_path <- function() {
  wgcna_claim_readiness_artifact("WGCNA_entity_claim_readiness.csv")
}

# The current module-to-supermodule membership (identity contract). Every structural expectation
# below is derived from it, so the contract follows the network that was actually built.
wgcna_claim_readiness_membership <- function(dataset = "microglia") {
  p <- wgcna_identity_contract_artifact("WGCNA_module_supermodule_membership_contract.csv", dataset)
  if (!file.exists(p)) stop("Membership contract missing for Stage 13 validation: ", p, call. = FALSE)
  m <- utils::read.csv(p, stringsAsFactors = FALSE)
  unique(m[, c("module_id", "supermodule_id")])
}

# Singleton supermodules are compatibility aliases of their one member module.
wgcna_claim_readiness_alias_map <- function(membership = wgcna_claim_readiness_membership()) {
  n <- table(membership$supermodule_id)
  single <- sort(names(n)[n == 1L])
  data.frame(
    entity_id = single,
    canonical_claim_entity_id = membership$module_id[match(single, membership$supermodule_id)],
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
  membership <- wgcna_claim_readiness_membership("microglia")
  n_sm <- table(membership$supermodule_id)
  expected_modules <- sort(unique(membership$module_id))
  expected_blocks <- sort(names(n_sm)[n_sm > 1L])
  n_expected <- length(expected_modules) + length(n_sm)
  if (nrow(readiness) != n_expected) {
    stop("Stage 13 must contain exactly ", n_expected, " technical identities (current membership); observed ", nrow(readiness), ".", call. = FALSE)
  }
  assert_unique_wgcna_keys(readiness, c("dataset", "level", "entity_id"), "Stage 13 claim readiness")
  if (!identical(unique(as.character(readiness$dataset)), "microglia")) {
    stop("Stage 13 downstream contract is restricted to dataset=microglia.", call. = FALSE)
  }

  roles <- table(as.character(readiness$claim_entity_role))
  expected_roles <- c(canonical_module = length(expected_modules), higher_order_block = length(expected_blocks),
                      compatibility_alias = as.integer(sum(n_sm == 1L)))
  observed_roles <- stats::setNames(rep(0L, length(expected_roles)), names(expected_roles))
  observed_roles[names(roles)] <- as.integer(roles)
  if (!identical(observed_roles, expected_roles)) {
    stop(
      "Stage 13 role counts must match the current membership (", paste(names(expected_roles), expected_roles, sep = "=", collapse = ", "), "); observed ",
      paste(names(observed_roles), observed_roles, sep = "=", collapse = ", "), ".",
      call. = FALSE
    )
  }

  modules <- readiness[readiness$claim_entity_role == "canonical_module", , drop = FALSE]
  blocks <- readiness[readiness$claim_entity_role == "higher_order_block", , drop = FALSE]
  aliases <- readiness[readiness$claim_entity_role == "compatibility_alias", , drop = FALSE]
  if (!all(modules$level == "module") || !setequal(modules$entity_id, expected_modules)) {
    stop("Stage 13 canonical-module identity set is invalid.", call. = FALSE)
  }
  if (!all(blocks$level == "supermodule") || !setequal(blocks$entity_id, expected_blocks)) {
    stop("Stage 13 higher-order blocks must be the multi-module supermodules of the current membership: ",
         paste(expected_blocks, collapse = ", "), ".", call. = FALSE)
  }

  expected_aliases <- wgcna_claim_readiness_alias_map(membership)
  observed_aliases <- as.data.frame(
    aliases[match(expected_aliases$entity_id, aliases$entity_id), c("entity_id", "canonical_claim_entity_id"), drop = FALSE],
    stringsAsFactors = FALSE
  )
  rownames(observed_aliases) <- NULL
  if (anyNA(observed_aliases$entity_id) || !identical(observed_aliases, expected_aliases)) {
    stop("Stage 13 compatibility aliases do not match the singleton supermodules of the current membership.", call. = FALSE)
  }
  if (any(aliases$separate_manuscript_claim_allowed %in% TRUE) ||
      any(as.character(aliases$compatibility_target_level) != "module") ||
      any(as.character(aliases$compatibility_target_entity_id) != as.character(aliases$canonical_claim_entity_id))) {
    stop("Stage 13 compatibility aliases must be module aliases and never separately claimable.", call. = FALSE)
  }

  # Which blocks and modules are claimable is a result of the Stage 13 audit, not a fixed contract;
  # it is reported, not pinned (the network frozen before 2026-10-19 had SM09 and ten modules).
  # An FDR-supported group effect is likewise a result to surface for review, not an error.
  if (any(as.character(readiness$group_effect_status) == "FDR_supported")) {
    warning("Stage 13 contains FDR-supported module or supermodule group effect(s): ",
            paste(readiness$entity_id[as.character(readiness$group_effect_status) == "FDR_supported"], collapse = ", "),
            ". Review before any claim is made.", call. = FALSE)
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

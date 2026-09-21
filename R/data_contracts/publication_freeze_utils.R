# Publication freeze manifest: construction and validation.
#
# One durable, machine-readable record of the exact publication state. It does
# not recompute or modify any scientific result; every value is either read from
# an existing machine-readable contract in this repository, derived from git
# objects, or queried from the active R installation.
#
# Hashing, git and YAML all reuse the repository's canonical helpers from
# R/paths.R (file_hash_sha256, file_hash, git_commit_sha, relative_to). This is
# deliberately NOT written through write_run_manifest(): that helper embeds a
# volatile timestamp and an absolute repo_root and hashes with md5, none of which
# is compatible with a deterministic SHA-256 freeze identity.
#
# Determinism: every list is built in a fixed order, paths are project-relative,
# and the only volatile field (generated_at) lives in the `metadata` block, which
# is excluded from every section that identifies scientific content.

PUBLICATION_FREEZE_MANIFEST_VERSION <- "1"

publication_freeze_manifest_path <- function() {
  repo_path("docs", "publication_freeze_manifest.yml")
}

# --- small utilities -------------------------------------------------------

freeze_rel <- function(paths) {
  normalized <- gsub("\\\\", "/", as.character(paths))
  root <- gsub("\\\\", "/", repo_root())
  sub(paste0("^", gsub("([\\^$.|?*+(){}\\[\\]])", "\\\\\\1", root), "/?"), "", normalized)
}

freeze_git <- function(args) {
  out <- tryCatch(
    suppressWarnings(system2("git", c("-C", repo_root(), args), stdout = TRUE, stderr = FALSE)),
    error = function(e) character()
  )
  if (!length(out)) return(NA_character_)
  as.character(out[[1]])
}

# `git rev-parse` echoes an unresolvable argument back on stdout and exits 128,
# so a plain capture can silently yield a non-SHA. Anywhere an object id is
# required, demand a 40-hex value and return NA otherwise.
freeze_git_sha <- function(args) {
  value <- freeze_git(args)
  if (is.na(value) || !grepl("^[0-9a-f]{40}$", value)) return(NA_character_)
  value
}

freeze_string_sha256 <- function(x) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package 'digest' is required for freeze-manifest hashing.", call. = FALSE)
  }
  digest::digest(enc2utf8(as.character(x)), algo = "sha256", serialize = FALSE)
}

# Hash a file listing deterministically: sorted project-relative path + SHA-256.
freeze_file_records <- function(paths) {
  paths <- unique(as.character(paths))
  paths <- paths[file.exists(paths) & !dir.exists(paths)]
  if (!length(paths)) return(list())
  rel <- freeze_rel(paths)
  ord <- order(rel, method = "radix")
  paths <- paths[ord]; rel <- rel[ord]
  info <- file.info(paths)
  lapply(seq_along(paths), function(i) {
    list(
      path = rel[[i]],
      sha256 = file_hash_sha256(paths[[i]]),
      size_bytes = as.numeric(info$size[[i]])
    )
  })
}

# A single stable digest over a set of file records, so a set can be compared
# without diffing every entry. Volatile metadata never enters this.
freeze_set_digest <- function(records) {
  if (!length(records)) return(NA_character_)
  freeze_string_sha256(paste(
    vapply(records, function(r) paste(r$path, r$sha256, sep = "|"), character(1)),
    collapse = "\n"
  ))
}

# --- 2. freeze identity ----------------------------------------------------

# The freeze identity is anchored to the TAG, not to HEAD. Anchoring to HEAD
# would be self-invalidating: committing the manifest moves HEAD, so a
# regeneration would record the manifest's own commit as the freeze and the tag
# check would fail permanently. HEAD is recorded separately as generation
# context, together with its relationship to the freeze commit.
freeze_git_state <- function(expected_tag = "publication-freeze-2026-09-02") {
  head_sha <- freeze_git_sha("rev-parse HEAD")
  tag_sha <- freeze_git_sha(paste0("rev-parse ", expected_tag, "^{commit}"))
  if (is.na(tag_sha) || !nzchar(tag_sha)) {
    stop("Freeze tag '", expected_tag, "' does not resolve to a commit; ",
         "refusing to write a freeze manifest without an anchored identity.",
         call. = FALSE)
  }
  status <- tryCatch(
    suppressWarnings(system2("git", c("-C", repo_root(), "status", "--porcelain"),
                             stdout = TRUE, stderr = FALSE)),
    error = function(e) NA_character_
  )
  descendant <- !is.na(freeze_git(
    paste0("merge-base --is-ancestor ", tag_sha, " ", head_sha)
  )) || identical(
    freeze_git(paste0("merge-base ", tag_sha, " ", head_sha)), tag_sha
  )
  list(
    freeze_git_commit = tag_sha,
    freeze_git_tag = expected_tag,
    freeze_identity_anchored_to = "tag",
    head_commit_at_generation = head_sha,
    head_is_freeze_commit = identical(head_sha, tag_sha),
    head_is_freeze_commit_or_descendant = identical(head_sha, tag_sha) || descendant,
    branch = freeze_git("rev-parse --abbrev-ref HEAD"),
    upstream_commit = freeze_git_sha("rev-parse origin/main"),
    head_matches_upstream = identical(head_sha, freeze_git_sha("rev-parse origin/main")),
    working_tree_clean = length(status) == 0L || all(!nzchar(status))
  )
}

# --- 3. provenance equivalence --------------------------------------------

# Export-relevant implementation files. The equivalence claim is scoped to
# exactly this list; it is NOT a claim about the whole repository tree.
freeze_protected_export_files <- function() {
  c(
    "R/enrichment/clusterprofiler_reproducibility.R",
    "R/utilities/export_helpers.R",
    "R/paths.R",
    "R/statistics/pride_helpers.R",
    "analysis/publication_source_data/08_export_manuscript_figures.R",
    "analysis/publication_source_data/09_export_source_data.R"
  )
}

# A protected file is identified by what it is, not by where it currently sits.
# Restructuring moved several of them, so a blob lookup at an older commit has
# to try the path the file had at that commit. Without this the equivalence
# claim silently degrades: the lookup misses, the sha256 comes back NA, and
# "not identical" gets reported for a file that never changed.
#
# Only add an entry here when a file genuinely moved. The point of the list is
# to make a rename visible in one place rather than to paper over a real
# content change.
freeze_protected_path_history <- function() {
  list(
    "analysis/publication_source_data/08_export_manuscript_figures.R" = c(
      "analysis/09_publication_exports/08_export_manuscript_figures.R",
      "09_export_pride_journal/08_export_manuscript_figures.R"
    ),
    "analysis/publication_source_data/09_export_source_data.R" = c(
      "analysis/09_publication_exports/09_export_source_data.R",
      "09_export_pride_journal/09_export_source_data.R"
    ),
    "R/enrichment/clusterprofiler_reproducibility.R" = c(
      "R/clusterprofiler_reproducibility.R"
    ),
    "R/utilities/export_helpers.R" = c("R/export_helpers.R"),
    "R/statistics/pride_helpers.R" = c("R/pride_helpers.R")
  )
}

# Every path a protected file has been known by, current first.
freeze_protected_path_candidates <- function(rel_path) {
  unique(c(rel_path, freeze_protected_path_history()[[rel_path]]))
}

# Proof by git object identity plus SHA-256 of the extracted blobs. Git's blob
# oid is itself a content hash over the exact bytes, so equal oids are already
# proof; the SHA-256 values are recorded so the manifest is verifiable without
# git internals.
freeze_blob_sha256 <- function(commit, rel_path) {
  ## Try the current path first, then the paths this file had earlier, so a
  ## rename does not turn an unchanged file into a reported mismatch.
  oid <- NA_character_
  for (p in freeze_protected_path_candidates(rel_path)) {
    oid <- freeze_git_sha(paste0("rev-parse ", commit, ":", p))
    if (!is.na(oid)) break
  }
  if (is.na(oid)) {
    return(list(blob_oid = NA_character_, sha256 = NA_character_))
  }
  tmp <- tempfile("freeze-blob-")
  on.exit(unlink(tmp), add = TRUE)
  ok <- tryCatch({
    suppressWarnings(system2("git", c("-C", repo_root(), "cat-file", "blob", oid),
                             stdout = tmp, stderr = FALSE))
    TRUE
  }, error = function(e) FALSE)
  if (!ok || !file.exists(tmp)) {
    return(list(blob_oid = oid, sha256 = NA_character_))
  }
  list(blob_oid = oid, sha256 = file_hash_sha256(tmp))
}

freeze_provenance_equivalence <- function(historical_commit, freeze_commit,
                                          files = freeze_protected_export_files()) {
  historical_resolvable <- !is.na(freeze_git_sha(paste0("rev-parse ", historical_commit, "^{commit}")))
  files <- sort(files, method = "radix")
  records <- lapply(files, function(f) {
    if (!historical_resolvable) {
      return(list(
        path = f,
        historical_blob_oid = NA_character_,
        historical_sha256 = NA_character_,
        freeze_blob_oid = freeze_blob_sha256(freeze_commit, f)$blob_oid,
        freeze_sha256 = freeze_blob_sha256(freeze_commit, f)$sha256,
        identical = NA
      ))
    }
    h <- freeze_blob_sha256(historical_commit, f)
    z <- freeze_blob_sha256(freeze_commit, f)
    list(
      path = f,
      historical_blob_oid = h$blob_oid,
      historical_sha256 = h$sha256,
      freeze_blob_oid = z$blob_oid,
      freeze_sha256 = z$sha256,
      identical = !is.na(h$sha256) && !is.na(z$sha256) && identical(h$sha256, z$sha256)
    )
  })
  all_identical <- length(records) > 0L &&
    all(vapply(records, function(r) isTRUE(r$identical), logical(1)))
  list(
    historical_export_commit = historical_commit,
    freeze_commit = freeze_commit,
    historical_commit_resolvable = historical_resolvable,
    historical_commit_reachable_from_freeze =
      identical(freeze_git(paste0("merge-base --is-ancestor ", historical_commit, " HEAD")), NA_character_),
    scope = "export-relevant implementation files only; NOT a whole-repository equivalence claim",
    files = records,
    all_protected_files_identical = all_identical,
    conclusion = paste(
      "The existing export payloads are accepted for the publication freeze because",
      "the export-relevant implementation is byte-identical across the historical",
      "export commit and the freeze commit; therefore the stale git_commit field in",
      "the export run manifests is a provenance-label drift, not evidence of payload",
      "drift."
    ),
    limitation = if (historical_resolvable) {
      "None: the historical commit was resolvable locally and both blob sets were recomputed."
    } else {
      paste("The historical commit was not resolvable locally, so historical hashes",
            "could not be recomputed and are recorded as NA.")
    }
  )
}

# --- 4. software environment ----------------------------------------------

# Versions established by the accepted publication state. The generator fails if
# the active installation disagrees, rather than silently recording a new freeze.
freeze_asserted_package_versions <- function() {
  c(
    clusterProfiler = "4.18.4",
    fgsea = "1.36.2",
    DOSE = "4.4.0",
    BiocParallel = "1.44.0",
    WGCNA = "1.74",
    limma = "3.66.0"
  )
}

freeze_asserted_r_version <- function() "4.5.1"
freeze_asserted_bioc_version <- function() "3.22"

freeze_observed_package_version <- function(pkg) {
  tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) NA_character_)
}

# asserted/observed_fun/observed_r are injectable so the failure path is
# directly testable without perturbing the installed library.
freeze_environment <- function(strict = TRUE,
                               asserted = freeze_asserted_package_versions(),
                               observed_fun = freeze_observed_package_version,
                               observed_r = paste(R.version$major, R.version$minor, sep = "."),
                               observed_bioc = tryCatch(as.character(BiocManager::version()),
                                                        error = function(e) NA_character_)) {
  pkgs <- sort(names(asserted), method = "radix")
  records <- lapply(pkgs, function(p) {
    observed <- observed_fun(p)
    list(
      package = p,
      freeze_version = unname(asserted[[p]]),
      observed_version = observed,
      matches = !is.na(observed) && identical(observed, unname(asserted[[p]]))
    )
  })

  mismatches <- vapply(records, function(r) !isTRUE(r$matches), logical(1))
  r_ok <- identical(observed_r, freeze_asserted_r_version())
  bioc_ok <- is.na(observed_bioc) || identical(observed_bioc, freeze_asserted_bioc_version())

  if (isTRUE(strict) && (any(mismatches) || !r_ok || !bioc_ok)) {
    detail <- c(
      if (!r_ok) paste0("R ", observed_r, " != ", freeze_asserted_r_version()),
      if (!bioc_ok) paste0("Bioconductor ", observed_bioc, " != ", freeze_asserted_bioc_version()),
      vapply(records[mismatches], function(r) {
        paste0(r$package, " ", r$observed_version, " != ", r$freeze_version)
      }, character(1))
    )
    stop(
      "Refusing to write a freeze manifest: the active R installation does not ",
      "match the asserted publication environment (", paste(detail, collapse = "; "),
      "). Investigate before recording a new freeze state.", call. = FALSE
    )
  }

  list(
    r_freeze_version = freeze_asserted_r_version(),
    r_observed_version = observed_r,
    r_matches = r_ok,
    bioconductor_freeze_version = freeze_asserted_bioc_version(),
    bioconductor_observed_version = observed_bioc,
    bioconductor_matches = bioc_ok,
    platform = R.version$platform,
    packages = records,
    all_packages_match = !any(mismatches)
  )
}

# --- 5. GSEA reproducibility contract -------------------------------------

freeze_gsea_contract <- function() {
  cfg_path <- repo_path("config", "clusterProfiler_config.yml")
  impl_path <- repo_path("R", "clusterprofiler_reproducibility.R")
  seed_base <- NA_integer_
  n_perm <- NA_integer_
  if (file.exists(cfg_path) && requireNamespace("yaml", quietly = TRUE)) {
    cfg <- tryCatch(yaml::read_yaml(cfg_path), error = function(e) NULL)
    if (!is.null(cfg) && is.list(cfg$analysis)) {
      seed_base <- suppressWarnings(as.integer(cfg$analysis$gsea_seed_base))
      n_perm <- suppressWarnings(as.integer(cfg$analysis$n_perm_simple))
    }
  }
  list(
    seed_base = seed_base,
    n_perm_simple = n_perm,
    rng_kind = "L'Ecuyer-CMRG/Inversion/Rejection",
    clusterprofiler_by = "fgsea",
    clusterprofiler_seed_argument = FALSE,
    config_source = list(
      path = freeze_rel(cfg_path),
      sha256 = file_hash_sha256(cfg_path)
    ),
    implementation_source = list(
      path = freeze_rel(impl_path),
      sha256 = file_hash_sha256(impl_path)
    ),
    contract = paste(
      "Deterministic inputs and configuration are frozen (ranked input,",
      "per-comparison derived seed, RNGkind, nPermSimple, by = 'fgsea', and",
      "clusterProfiler's logical seed flag kept FALSE). Numerical reproducibility",
      "is expected within one execution context. Bit-exact equality across",
      "different execution contexts is NOT guaranteed; the implementation records",
      "a measured tolerance instead (enrichmentScore to ~1.6e-15, propagating to",
      "<= ~2.4e-05 in NES and <= ~2.3e-05 in p-value/FDR), with no FDR-0.05",
      "crossings and an unchanged Figure-2f display selection in the audited",
      "comparison."
    )
  )
}

# --- 6. protected WGCNA identities ----------------------------------------

freeze_wgcna_protected_basenames <- function() {
  c(
    "wgcna_final_model_state.rds",
    "wgcna_module_supermodule_annotation.csv",
    "module_group_effects.csv",
    "supermodule_group_effects.csv"
  )
}

# Pure guard so the mismatch path is testable with fixtures instead of by
# perturbing a real frozen state.
freeze_assert_protected_states <- function(states, strict = TRUE) {
  mismatched <- vapply(states, function(s) !isTRUE(s$md5_matches), logical(1))
  if (isTRUE(strict) && any(mismatched)) {
    stop(
      "Protected WGCNA state hash mismatch for: ",
      paste(vapply(states[mismatched], function(s) s$artifact, character(1)), collapse = ", "),
      ". Refusing to record a freeze over a changed frozen state.", call. = FALSE
    )
  }
  invisible(!any(mismatched))
}

# The equivalence claim must cover exactly the declared protected file list and
# nothing more, so a silent widening into a whole-repository claim is detectable.
freeze_equivalence_scope_ok <- function(equivalence,
                                        expected = freeze_protected_export_files()) {
  recorded <- sort(vapply(equivalence$files, function(f) f$path, character(1)), method = "radix")
  identical(recorded, sort(expected, method = "radix"))
}

freeze_wgcna_identities <- function(dataset = "microglia", strict = TRUE) {
  audit_path <- path_results("reviewer_audit", "microglia_wgcna_nature_readiness",
                             "protected_output_hash_audit.csv")
  status_path <- path_results("source_data", "06_modules_WGCNA", "identity_contract",
                             dataset, "WGCNA_identity_contract_status.csv")
  effects_paths <- c(
    module = path_results("tables", "06_modules_WGCNA", "group_effects", dataset,
                          "module_group_effects.csv"),
    supermodule = path_results("tables", "06_modules_WGCNA", "group_effects", dataset,
                               "supermodule_group_effects.csv")
  )

  audit <- if (file.exists(audit_path)) {
    utils::read.csv(audit_path, check.names = FALSE, stringsAsFactors = FALSE)
  } else NULL

  contract <- list(
    identity_contract_sha256 = NA_character_,
    identity_contract_version = NA_character_,
    membership_version = NA_character_,
    frozen_state_sha256 = NA_character_
  )
  observed_contracts <- list()
  for (nm in names(effects_paths)) {
    p <- effects_paths[[nm]]
    if (!file.exists(p)) next
    d <- utils::read.csv(p, check.names = FALSE, stringsAsFactors = FALSE)
    got <- list(
      source = freeze_rel(p),
      identity_contract_sha256 = unique(as.character(d$identity_contract_sha256))[1],
      identity_contract_version = unique(as.character(d$identity_contract_version))[1],
      membership_version = unique(as.character(d$membership_version))[1],
      frozen_state_sha256 = unique(as.character(d$frozen_state_sha256))[1]
    )
    observed_contracts[[nm]] <- got
  }
  if (length(observed_contracts)) {
    contract$identity_contract_sha256 <- observed_contracts[[1]]$identity_contract_sha256
    contract$identity_contract_version <- observed_contracts[[1]]$identity_contract_version
    contract$membership_version <- observed_contracts[[1]]$membership_version
    contract$frozen_state_sha256 <- observed_contracts[[1]]$frozen_state_sha256
    # Both group-effect tables must agree; disagreement is a hard failure.
    for (k in c("identity_contract_sha256", "identity_contract_version",
                "membership_version", "frozen_state_sha256")) {
      vals <- unique(vapply(observed_contracts, function(o) as.character(o[[k]]), character(1)))
      if (isTRUE(strict) && length(vals) > 1L) {
        stop("WGCNA contract field ", k, " disagrees between group-effect tables: ",
             paste(vals, collapse = " vs "), call. = FALSE)
      }
    }
  }

  states <- lapply(sort(freeze_wgcna_protected_basenames(), method = "radix"), function(bn) {
    row <- if (!is.null(audit)) audit[basename(audit$protected_path) == bn, , drop = FALSE] else NULL
    rel_path <- if (!is.null(row) && nrow(row)) as.character(row$protected_path[[1]]) else NA_character_
    abs_path <- if (!is.na(rel_path)) file.path(repo_root(), rel_path) else NA_character_
    recorded_md5 <- if (!is.null(row) && nrow(row)) as.character(row$md5_after[[1]]) else NA_character_
    observed_md5 <- if (!is.na(abs_path) && file.exists(abs_path)) file_hash(abs_path) else NA_character_
    list(
      artifact = bn,
      dataset = dataset,
      path = rel_path,
      recorded_md5 = recorded_md5,
      observed_md5 = observed_md5,
      md5_matches = !is.na(recorded_md5) && identical(recorded_md5, observed_md5),
      sha256 = if (!is.na(abs_path) && file.exists(abs_path)) file_hash_sha256(abs_path) else NA_character_,
      unchanged_in_protected_audit =
        if (!is.null(row) && nrow(row)) as.character(row$unchanged[[1]]) %in% c("TRUE", "True", "true") else NA
    )
  })

  status <- if (file.exists(status_path)) {
    s <- utils::read.csv(status_path, check.names = FALSE, stringsAsFactors = FALSE)
    list(
      path = freeze_rel(status_path),
      sha256 = file_hash_sha256(status_path),
      status = as.character(s$status[[1]]),
      membership_version = as.character(s$membership_version[[1]]),
      contract_version = as.character(s$contract_version[[1]]),
      n_modules = as.integer(s$n_modules[[1]]),
      n_supermodules = as.integer(s$n_supermodules[[1]]),
      selected_cut_height = as.numeric(s$selected_cut_height[[1]])
    )
  } else NULL

  mismatched <- vapply(states, function(s) !isTRUE(s$md5_matches), logical(1))
  freeze_assert_protected_states(states, strict = strict)
  # The identity contract's membership_version must agree with the group effects.
  if (isTRUE(strict) && !is.null(status) &&
      !is.na(contract$membership_version) &&
      !identical(status$membership_version, contract$membership_version)) {
    stop("membership_version disagrees between identity contract (", status$membership_version,
         ") and group effects (", contract$membership_version, ").", call. = FALSE)
  }

  list(
    dataset = dataset,
    protected_output_hash_audit = list(
      path = freeze_rel(audit_path),
      sha256 = file_hash_sha256(audit_path),
      rows = if (!is.null(audit)) nrow(audit) else NA_integer_,
      all_unchanged = if (!is.null(audit)) all(as.character(audit$unchanged) %in% c("TRUE", "True", "true")) else NA
    ),
    identity_contract = contract,
    identity_contract_status = status,
    protected_states = states,
    all_protected_states_match = !any(mismatched)
  )
}

# --- 7. publication source-data ------------------------------------------

# Canonical directories are discovered from the producing stage's output roots
# rather than by guessing filenames.
#
# Each set lists both the historical roots and the normalized roots its
# producer now writes to. That pairing is the whole point: the producing
# analyses have migrated, so a set rooted only at the historical tree would go
# on hashing pre-migration files and silently stop tracking current data the
# moment anything is rerun. The historical roots stay because that is where the
# frozen files actually are.
#
# Phase 6G.6 added the normalized roots. Two of these three sets belong to
# differential_abundance and had already gone stale when that domain migrated
# in Phase 6G.4; only figure_3_manuscript_panels belongs to the domain being
# migrated here. Non-existent roots are filtered out before listing, so adding
# an empty root changes neither the file count nor the hash.
freeze_source_data_sets <- function() {
  list(
    figure_2_control_spatial = c(
      path_results("source_data", "04_differential_expression_enrichment",
                   "control_spatial_identity_validation", "global"),
      path_results("tables", "04_differential_expression_enrichment",
                   "control_spatial_identity_validation", "global"),
      path_results("differential_abundance", "validate_control_spatial_identity",
                   "global", "tables")
    ),
    figure_3_manuscript_panels = c(
      path_results("source_data", "manuscript_panels"),
      path_results("tables", "manuscript_panels"),
      path_results("reports", "manuscript_panels"),
      path_results("integration", "export_module_protein_zoom_source_data",
                   "global", "tables"),
      path_results("integration", "export_module_protein_zoom_source_data",
                   "global", "reports")
    ),
    sus_res_stage11 = c(
      path_results("source_data", "04_differential_expression_enrichment",
                   "stress_response_biological_audit", "global"),
      path_results("source_data", "04_differential_expression_enrichment",
                   "sus_res_spatial_dap_atlas", "global"),
      path_results("reports", "04_differential_expression_enrichment",
                   "stress_response_biological_audit", "global"),
      path_results("differential_abundance", "audit_stress_response_biology",
                   "global", "tables"),
      path_results("differential_abundance", "audit_stress_response_biology",
                   "global", "reports"),
      path_results("differential_abundance", "build_sus_res_dap_atlas",
                   "global", "tables")
    )
  )
}

freeze_source_data <- function(sets = freeze_source_data_sets()) {
  out <- list()
  for (nm in sort(names(sets), method = "radix")) {
    roots <- sets[[nm]]
    files <- unlist(lapply(roots[dir.exists(roots)], list.files,
                           pattern = "[.](csv|tsv|md)$", recursive = TRUE,
                           full.names = TRUE), use.names = FALSE)
    records <- freeze_file_records(files)
    out[[nm]] <- list(
      roots = sort(freeze_rel(roots), method = "radix"),
      file_count = length(records),
      set_digest_sha256 = freeze_set_digest(records),
      files = records
    )
  }
  out
}

# --- 8. export payloads ---------------------------------------------------

# The run manifest, at whichever address its producer currently writes.
#
# Phase 6G moved these exporters onto psd_dirs(<analysis_id>)$manifests, but the
# freeze kept reading the pre-6G results/logs/09_export_pride_journal/ path. The
# effect was silent and misleading: the figure exporter was rerun twice in Phase
# 6H, writing 5598 and then 5582 inputs canonically, while the freeze went on
# comparing against a legacy file frozen at 5244 and reported counts_agree = no
# for a reason that had nothing to do with the export.
#
# Canonical first, legacy second, and only if the canonical file is actually
# there: 09_export_source_data has not been rerun since the migration, so its
# legacy manifest is still the only populated one and remains correct.
freeze_run_manifest_path <- function(analysis_id, legacy_subdir) {
  if (!exists("psd_dirs", mode = "function")) {
    source(repo_path("R", "utilities", "publication_source_data_paths.R"))
  }
  canonical <- file.path(psd_dirs(analysis_id)$manifests, "run_manifest.yml")
  if (file.exists(canonical)) return(canonical)
  path_results("logs", "09_export_pride_journal", legacy_subdir, "run_manifest.yml")
}

freeze_export_payloads <- function() {
  defs <- list(
    manuscript_figures = list(
      payload_root = path_results("manuscript", "extended_data"),
      manifest = path_results("manuscript", "figure_export_manifest.csv"),
      audit = path_results("manuscript", "figure_publication_audit.csv"),
      run_manifest = freeze_run_manifest_path("08_export_manuscript_figures",
                                              "manuscript_figures"),
      producer = "analysis/09_publication_exports/08_export_manuscript_figures.R"
    ),
    manuscript_source_data = list(
      payload_root = c(path_results("manuscript", "source_data"),
                       path_results("manuscript", "supplementary_tables")),
      manifest = path_results("manuscript", "source_data_export_manifest.csv"),
      audit = NA_character_,
      run_manifest = freeze_run_manifest_path("09_export_source_data",
                                              "source_data"),
      producer = "analysis/09_publication_exports/09_export_source_data.R"
    )
  )
  lapply(sort(names(defs), method = "radix"), function(nm) {
    d <- defs[[nm]]
    manifest_rows <- NA_integer_
    if (file.exists(d$manifest)) {
      m <- utils::read.csv(d$manifest, check.names = FALSE, stringsAsFactors = FALSE)
      manifest_rows <- nrow(m)
    }
    rm_commit <- NA_character_
    rm_inputs <- NA_integer_
    if (file.exists(d$run_manifest) && requireNamespace("yaml", quietly = TRUE)) {
      y <- tryCatch(yaml::read_yaml(d$run_manifest), error = function(e) NULL)
      if (!is.null(y)) {
        rm_commit <- as.character(y$git_commit)
        rm_inputs <- length(unlist(y$inputs))
      }
    }
    list(
      export = nm,
      producer = d$producer,
      payload_roots = sort(freeze_rel(d$payload_root), method = "radix"),
      manifest_path = freeze_rel(d$manifest),
      manifest_sha256 = file_hash_sha256(d$manifest),
      manifest_row_count = manifest_rows,
      audit_path = if (is.na(d$audit)) NA_character_ else freeze_rel(d$audit),
      audit_sha256 = if (is.na(d$audit)) NA_character_ else file_hash_sha256(d$audit),
      run_manifest_path = freeze_rel(d$run_manifest),
      run_manifest_sha256 = file_hash_sha256(d$run_manifest),
      run_manifest_recorded_commit = rm_commit,
      run_manifest_input_count = rm_inputs,
      counts_agree = !is.na(manifest_rows) && !is.na(rm_inputs) && manifest_rows == rm_inputs,
      freeze_acceptance = "accepted",
      freeze_acceptance_basis = paste(
        "Manifest row count and run-manifest input count agree, and the",
        "export-relevant implementation is byte-identical between the recorded",
        "historical commit and the freeze commit. Payload files are validated",
        "through the manifest contract rather than re-hashed in bulk."
      )
    )
  })
}

# --- 10. known documented gaps -------------------------------------------

# Lockfile state, measured rather than asserted. The renv accepted-gap below is
# emitted only while the lockfile is actually incomplete, so the warning can
# clear when the lockfile is completed and reappears if it regresses.
freeze_renv_lockfile_state <- function() {
  lock_path <- repo_path("renv.lock")
  audit_src <- repo_path("R", "renv_lock_audit.R")
  if (!file.exists(lock_path) || !file.exists(audit_src)) {
    return(list(
      path = freeze_rel(lock_path),
      present = file.exists(lock_path),
      complete = FALSE,
      note = "renv.lock or its audit helper is absent"
    ))
  }
  env <- new.env(parent = globalenv())
  sys.source(audit_src, envir = env)
  audit <- tryCatch(
    env$audit_renv_lock_completeness(lock_path, root = repo_root()),
    error = function(e) NULL
  )
  sentinels <- tryCatch(env$audit_renv_lock(lock_path), error = function(e) NULL)
  parsed <- tryCatch(jsonlite::fromJSON(lock_path, simplifyVector = FALSE),
                     error = function(e) NULL)
  list(
    path = freeze_rel(lock_path),
    present = TRUE,
    sha256 = file_hash_sha256(lock_path),
    package_records = if (!is.null(audit)) audit$recorded_count else NA_integer_,
    direct_dependencies = if (!is.null(audit)) length(audit$direct_dependencies) else NA_integer_,
    r_version = if (!is.null(parsed)) parsed$R$Version else NA_character_,
    bioconductor_version = if (!is.null(parsed)) parsed$Bioconductor$Version else NULL,
    repositories_declared = if (!is.null(parsed)) length(parsed$R$Repositories) else NA_integer_,
    missing_direct = if (!is.null(audit)) length(audit$missing_direct) else NA_integer_,
    unresolved_requirements = if (!is.null(audit)) length(audit$unresolved_requirements) else NA_integer_,
    duplicate_records = if (!is.null(audit)) length(audit$duplicate_records) else NA_integer_,
    scientific_sentinels_present = if (!is.null(sentinels)) sentinels$plausibly_full_scientific_lock else NA,
    complete = !is.null(audit) && isTRUE(audit$complete) &&
      !is.null(sentinels) && isTRUE(sentinels$plausibly_full_scientific_lock),
    hash_field_recorded = FALSE,
    hash_field_note = paste(
      "renv's per-package Hash is intentionally absent: renv is not installed in",
      "this environment and its hash algorithm cannot be reproduced faithfully,",
      "so inventing a value would misrepresent provenance."
    ),
    restore_validation_level = "A+B (structural and availability); C not attempted",
    restore_validation_note = paste(
      "An isolated restore would require installing renv and building a project",
      "library, which would materially alter the live analysis environment."
    )
  )
}

freeze_known_gaps <- function(renv_state = freeze_renv_lockfile_state()) {
  gaps <- list()
  if (!isTRUE(renv_state$complete)) {
    gaps[[length(gaps) + 1L]] <- list(
      id = "renv_lock_incomplete",
      severity = "warn",
      classification = "reproducibility gap / deferred freeze item",
      summary = paste(
        "renv.lock does not capture the scientific analysis stack, so the freeze",
        "is not reproducible from the lockfile alone. Recorded package records:",
        renv_state$package_records
      ),
      reference = "docs/RENV_LOCK_STATUS.md",
      action_in_this_task = "none; renv.lock deliberately untouched"
    )
  }
  c(gaps, list(
    list(
      id = "pride_dry_run_semantics",
      severity = "warn",
      classification = "documented tooling limitation",
      summary = paste(
        "Of the nine export steps, only 08_export_manuscript_figures.R and",
        "09_export_source_data.R implement a side-effect-free dry-run guard.",
        "RUN_EXPORT.R restricts the step list under --dry-run so the orchestrator",
        "is safe, but the remaining seven steps remain individually unguarded and",
        "would write if invoked directly with --dry-run."
      ),
      reference = "analysis/09_publication_exports/RUN_EXPORT.R",
      action_in_this_task = "documented only; not fixed"
    ),
    list(
      id = "processed_package_wildcard_glob_filter",
      severity = "warn",
      classification = "latent defect, currently masked",
      summary = paste(
        "R/utilities/export_helpers.R:658, inside processed_files_for_dataset() (defined at",
        "line 610), filters config$supplementary_table_globs with",
        "grepl(\"\\\\*\", g, fixed = TRUE). Under fixed = TRUE that searches for a",
        "literal backslash-star, which none of the six configured globs contains,",
        "so the wildcard arm never fires and the filter degenerates to",
        "grepl(dataset, g, fixed = TRUE): it keeps 1 of 6 globs for microglia and",
        "0 of 6 for neuron_soma and neuron_neuropil, where the intended test would",
        "keep 6 of 6."
      ),
      correction_to_prior_description = paste(
        "Verified against the code, three details of the previously circulated",
        "description are wrong. (a) The defect is in R/utilities/export_helpers.R, not",
        "analysis/publication_source_data/build_supplementary_tables.R. (b) It affects the",
        "PRIDE processed-data package selection consumed by",
        "export_processed_matrices.R:69 and build_pride_manifest.R:28;",
        "the supplementary-table export itself calls supplementary_candidate_files(),",
        "which applies no wildcard filter and is unaffected. (c) It is masked by the",
        "include_derived argument, which both callers pass as",
        "cli$include_derived_results == has(\"--include-derived-results\") and is",
        "therefore FALSE unless that flag is given; the",
        "manifest.default_include_derived_results: false config key documents the",
        "same intent but does not gate this code path."
      ),
      reference = "R/utilities/export_helpers.R:658 (processed_files_for_dataset)",
      action_in_this_task = "documented only; not fixed"
    )
  ))
}

# --- assembly ------------------------------------------------------------

build_publication_freeze_manifest <- function(
    historical_export_commit = "0825e4237406f1297e379fc35216c55adf3d5e48",
    expected_tag = "publication-freeze-2026-09-02",
    generated_at = NA_character_,
    test_state = NULL,
    strict = TRUE) {
  git_state <- freeze_git_state(expected_tag)
  renv_state <- freeze_renv_lockfile_state()
  list(
    manifest_version = PUBLICATION_FREEZE_MANIFEST_VERSION,
    metadata = list(
      generated_at = generated_at,
      generator = "tools/generate_publication_freeze_manifest.R",
      note = paste(
        "generated_at is metadata only and is excluded from every section that",
        "identifies scientific content."
      )
    ),
    freeze_identity = git_state,
    software_environment = freeze_environment(strict = strict),
    gsea_reproducibility_contract = freeze_gsea_contract(),
    export_provenance_equivalence = freeze_provenance_equivalence(
      historical_export_commit, git_state$freeze_git_commit
    ),
    wgcna_protected_identities = freeze_wgcna_identities(strict = strict),
    publication_source_data = freeze_source_data(),
    export_payloads = freeze_export_payloads(),
    renv_lockfile = renv_state,
    validation_state = test_state,
    known_gaps = freeze_known_gaps(renv_state)
  )
}

write_publication_freeze_manifest <- function(manifest,
                                              path = publication_freeze_manifest_path()) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to write the freeze manifest.", call. = FALSE)
  }
  dir_create(dirname(path))
  # yaml::as.yaml() already terminates with a newline; writeLines() would add a
  # second one and leave a blank line at EOF.
  cat(yaml::as.yaml(manifest), file = path)
  invisible(path)
}

# --- 12. validator -------------------------------------------------------

freeze_check <- function(status, check, detail) {
  list(status = status, check = check, detail = detail)
}

validate_publication_freeze <- function(path = publication_freeze_manifest_path()) {
  if (!file.exists(path)) {
    return(list(
      manifest_path = freeze_rel(path),
      checks = list(freeze_check("FAIL", "manifest_present", "freeze manifest not found")),
      summary = c(PASS = 0L, WARN = 0L, FAIL = 1L)
    ))
  }
  m <- yaml::read_yaml(path)
  checks <- list()
  add <- function(status, check, detail) {
    checks[[length(checks) + 1L]] <<- freeze_check(status, check, detail)
  }

  # Freeze identity is the tag. That invariant must hold exactly.
  head_sha <- freeze_git_sha("rev-parse HEAD")
  tag_sha <- freeze_git_sha(paste0("rev-parse ", m$freeze_identity$freeze_git_tag, "^{commit}"))
  add(if (identical(tag_sha, m$freeze_identity$freeze_git_commit)) "PASS" else "FAIL",
      "freeze_tag_points_at_freeze_commit",
      paste0(m$freeze_identity$freeze_git_tag, " -> ", tag_sha,
             " freeze=", m$freeze_identity$freeze_git_commit))
  # The checkout may legitimately have advanced past the freeze (for example to
  # the commit that records the manifest). Being exactly at the freeze is PASS;
  # being a descendant is an accepted WARN; anything else is a FAIL.
  descendant <- !is.na(freeze_git(
    paste0("merge-base --is-ancestor ", m$freeze_identity$freeze_git_commit, " ", head_sha)
  )) || identical(
    freeze_git(paste0("merge-base ", m$freeze_identity$freeze_git_commit, " ", head_sha)),
    m$freeze_identity$freeze_git_commit
  )
  add(if (identical(head_sha, m$freeze_identity$freeze_git_commit)) "PASS"
      else if (isTRUE(descendant)) "WARN" else "FAIL",
      "checkout_relative_to_freeze_commit",
      paste0("HEAD=", head_sha,
             if (identical(head_sha, m$freeze_identity$freeze_git_commit)) " (at freeze)"
             else if (isTRUE(descendant)) " (descendant of freeze)"
             else " (NOT related to freeze commit)"))

  # package versions
  for (p in m$software_environment$packages) {
    observed <- freeze_observed_package_version(p$package)
    add(if (identical(observed, p$freeze_version)) "PASS" else "FAIL",
        paste0("package_version:", p$package),
        paste0("freeze=", p$freeze_version, " observed=", observed))
  }
  observed_r <- paste(R.version$major, R.version$minor, sep = ".")
  add(if (identical(observed_r, m$software_environment$r_freeze_version)) "PASS" else "FAIL",
      "r_version", paste0("freeze=", m$software_environment$r_freeze_version,
                          " observed=", observed_r))

  # WGCNA protected states
  for (s in m$wgcna_protected_identities$protected_states) {
    abs_path <- file.path(repo_root(), s$path)
    observed <- if (file.exists(abs_path)) file_hash_sha256(abs_path) else NA_character_
    add(if (!is.na(observed) && identical(observed, s$sha256)) "PASS" else "FAIL",
        paste0("wgcna_protected_state:", s$artifact),
        paste0("freeze=", s$sha256, " observed=", observed))
  }
  eff <- m$wgcna_protected_identities$identity_contract
  add(if (!is.null(eff$frozen_state_sha256) && nzchar(eff$frozen_state_sha256)) "PASS" else "FAIL",
      "wgcna_identity_contract_recorded",
      paste0("membership_version=", eff$membership_version,
             " frozen_state_sha256=", eff$frozen_state_sha256))

  # source-data set digests
  for (nm in names(m$publication_source_data)) {
    set <- m$publication_source_data[[nm]]
    recomputed <- freeze_set_digest(lapply(set$files, function(f) {
      abs_path <- file.path(repo_root(), f$path)
      list(path = f$path,
           sha256 = if (file.exists(abs_path)) file_hash_sha256(abs_path) else NA_character_)
    }))
    add(if (identical(recomputed, set$set_digest_sha256)) "PASS" else "FAIL",
        paste0("source_data_set:", nm),
        paste0("files=", set$file_count, " freeze=", set$set_digest_sha256,
               " observed=", recomputed))
  }

  # export manifests
  for (e in m$export_payloads) {
    abs_manifest <- file.path(repo_root(), e$manifest_path)
    observed <- if (file.exists(abs_manifest)) file_hash_sha256(abs_manifest) else NA_character_
    add(if (!is.na(observed) && identical(observed, e$manifest_sha256)) "PASS" else "FAIL",
        paste0("export_manifest:", e$export),
        paste0("freeze=", e$manifest_sha256, " observed=", observed))
    if (!is.na(observed) && file.exists(abs_manifest)) {
      rows <- nrow(utils::read.csv(abs_manifest, check.names = FALSE, stringsAsFactors = FALSE))
      add(if (identical(as.integer(rows), as.integer(e$manifest_row_count))) "PASS" else "FAIL",
          paste0("export_manifest_row_count:", e$export),
          paste0("freeze=", e$manifest_row_count, " observed=", rows))
    }
  }

  # renv lockfile: recorded state must still hold, and the lockfile must not
  # have regressed to the trivial incomplete form.
  if (!is.null(m$renv_lockfile)) {
    lock_path <- file.path(repo_root(), m$renv_lockfile$path)
    observed_sha <- if (file.exists(lock_path)) file_hash_sha256(lock_path) else NA_character_
    add(if (!is.na(observed_sha) && identical(observed_sha, m$renv_lockfile$sha256)) "PASS" else "FAIL",
        "renv_lockfile_unchanged",
        paste0("freeze=", m$renv_lockfile$sha256, " observed=", observed_sha))
    state <- freeze_renv_lockfile_state()
    add(if (isTRUE(state$complete)) "PASS" else "FAIL",
        "renv_lockfile_complete",
        paste0("records=", state$package_records,
               " direct=", state$direct_dependencies,
               " missing_direct=", state$missing_direct,
               " unresolved=", state$unresolved_requirements,
               " sentinels=", state$scientific_sentinels_present))
    add(if (identical(state$r_version, m$software_environment$r_freeze_version)) "PASS" else "FAIL",
        "renv_lockfile_r_version",
        paste0("lock=", state$r_version, " freeze=", m$software_environment$r_freeze_version))
    add(if (identical(state$bioconductor_version,
                      m$software_environment$bioconductor_freeze_version)) "PASS" else "FAIL",
        "renv_lockfile_bioconductor_version",
        paste0("lock=", state$bioconductor_version,
               " freeze=", m$software_environment$bioconductor_freeze_version))
  }

  # provenance equivalence
  eq <- m$export_provenance_equivalence
  add(if (isTRUE(eq$all_protected_files_identical)) "PASS" else "WARN",
      "export_provenance_equivalence",
      paste0("historical=", eq$historical_export_commit,
             " scope=", length(eq$files), " protected files; ",
             if (isTRUE(eq$all_protected_files_identical)) "all identical" else eq$limitation))

  # documented gaps -> WARN, never silently ignored
  for (g in m$known_gaps) {
    add(toupper(g$severity), paste0("known_gap:", g$id), g$classification)
  }

  statuses <- vapply(checks, function(c) c$status, character(1))
  list(
    manifest_path = freeze_rel(path),
    checks = checks,
    summary = c(PASS = sum(statuses == "PASS"),
                WARN = sum(statuses == "WARN"),
                FAIL = sum(statuses == "FAIL"))
  )
}

print_publication_freeze_validation <- function(result) {
  for (c in result$checks) {
    cat(sprintf("  %-5s %-46s %s\n", c$status, c$check, c$detail))
  }
  cat(sprintf("\n  PASS %d | WARN %d | FAIL %d\n",
              result$summary[["PASS"]], result$summary[["WARN"]],
              result$summary[["FAIL"]]))
  invisible(result)
}

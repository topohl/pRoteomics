# Shared path helpers for pRoteomics.
#
# All committed scripts should default to repository-relative paths. For local
# machines, set PROTEOMICS_PROJECT_ROOT or provide a config override outside Git.

.null_coalescing_source_files <- unlist(lapply(
  sys.frames(),
  function(frame) if (!is.null(frame$ofile)) frame$ofile else character()
), use.names = FALSE)
.null_coalescing_search_dirs <- c(
  dirname(.null_coalescing_source_files[file.exists(
    .null_coalescing_source_files
  )]),
  Sys.getenv("PROTEOMICS_PROJECT_ROOT", unset = ""),
  getwd()
)
.null_coalescing_search_dirs <- unique(
  .null_coalescing_search_dirs[nzchar(.null_coalescing_search_dirs)]
)
.null_coalescing_file <- ""
for (.null_coalescing_start in .null_coalescing_search_dirs) {
  .null_coalescing_search_dir <- normalizePath(
    .null_coalescing_start, winslash = "/", mustWork = FALSE
  )
  repeat {
    .null_coalescing_candidate <- Filter(file.exists, c(
      file.path(
        .null_coalescing_search_dir, "R", "utilities", "null_coalescing.R"
      ),
      file.path(.null_coalescing_search_dir, "R", "null_coalescing.R")
    ))[1]
    if (!is.na(.null_coalescing_candidate)) {
      .null_coalescing_file <- .null_coalescing_candidate
      break
    }
    .null_coalescing_parent <- dirname(.null_coalescing_search_dir)
    if (identical(.null_coalescing_parent, .null_coalescing_search_dir)) break
    .null_coalescing_search_dir <- .null_coalescing_parent
  }
  if (nzchar(.null_coalescing_file)) break
}
if (!nzchar(.null_coalescing_file)) {
  stop("Could not locate R/null_coalescing.R.", call. = FALSE)
}
source(.null_coalescing_file)
rm(
  .null_coalescing_source_files, .null_coalescing_search_dirs,
  .null_coalescing_start, .null_coalescing_search_dir,
  .null_coalescing_candidate, .null_coalescing_parent,
  .null_coalescing_file
)

repo_root <- function() {
  env_root <- Sys.getenv("PROTEOMICS_PROJECT_ROOT", unset = "")
  if (nzchar(env_root)) {
    return(normalizePath(env_root, winslash = "/", mustWork = FALSE))
  }

  if (requireNamespace("rprojroot", quietly = TRUE)) {
    root <- tryCatch(
      rprojroot::find_root(rprojroot::has_file("README.md") | rprojroot::is_git_root),
      error = function(e) NULL
    )
    if (!is.null(root)) return(normalizePath(root, winslash = "/", mustWork = FALSE))
  }

  cur <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  repeat {
    markers <- c(".git", "README.md", "pipeline.yml", "analysis")
    if (any(file.exists(file.path(cur, markers)))) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) {
      stop("Could not locate repository root. Set PROTEOMICS_PROJECT_ROOT.", call. = FALSE)
    }
    cur <- parent
  }
}

# R/ is organised into domain subdirectories (data_contracts, qc, statistics,
# spatial, enrichment, networks, utilities). Call sites address a library by
# bare name -- repo_path("R", "module_stats.R") -- so the domain layout can be
# changed without editing the 800+ source() lines that name these libraries.
# This resolver maps a bare library name onto its domain directory and is the
# minimum path abstraction the layout requires.
.r_library_cache <- new.env(parent = emptyenv())

r_library_path <- function(name, root = repo_root()) {
  flat <- file.path(root, "R", name)
  if (file.exists(flat)) return(flat)

  key <- paste0("index:", root)
  if (!exists(key, envir = .r_library_cache, inherits = FALSE)) {
    rel <- list.files(
      file.path(root, "R"), pattern = "[.]R$", recursive = TRUE
    )
    index <- stats::setNames(file.path(root, "R", rel), basename(rel))
    assign(key, index[!duplicated(names(index))], envir = .r_library_cache)
  }
  index <- get(key, envir = .r_library_cache, inherits = FALSE)

  hit <- unname(index[name])
  if (!is.na(hit)) return(hit)
  flat
}

repo_path <- function(...) {
  parts <- list(...)
  if (length(parts) == 2L &&
      identical(as.character(parts[[1]]), "R") &&
      length(parts[[2]]) == 1L &&
      grepl("[.]R$", as.character(parts[[2]]))) {
    return(r_library_path(as.character(parts[[2]])))
  }
  file.path(repo_root(), ...)
}

path_raw <- function(...) repo_path("data", "raw", ...)
path_metadata <- function(...) repo_path("data", "metadata", ...)
path_external <- function(...) repo_path("data", "external", ...)
path_processed <- function(...) repo_path("data", "processed", ...)
path_results <- function(...) repo_path("results", ...)

# --- output lifecycles (config/output_layout.yml) -------------------------
# work/    regenerable intermediates, never cited, never canonical
# results/ canonical scientific results
# exports/ frozen outward-facing bundles, the only thing the manuscript reads
path_work <- function(...) repo_path("work", ...)
path_export <- function(...) repo_path("exports", ...)

.output_layout_cache <- new.env(parent = emptyenv())

output_layout <- function() {
  if (!is.null(.output_layout_cache$layout)) return(.output_layout_cache$layout)
  f <- repo_path("config", "output_layout.yml")
  if (!file.exists(f) || !requireNamespace("yaml", quietly = TRUE)) return(NULL)
  .output_layout_cache$layout <- yaml::read_yaml(f)
  .output_layout_cache$layout
}

output_layout_domains <- function() {
  l <- output_layout()
  if (is.null(l)) return(character(0))
  as.character(unlist(l$domains, use.names = FALSE))
}

output_layout_children <- function() {
  l <- output_layout()
  if (is.null(l)) return(character(0))
  names(l$canonical_result_layout$children_in_use)
}

# results/<domain>/<analysis_id>/<scope>/<child>/...
#
# analysis_id is the canonical owner's script stem, so an output is addressed
# by the analysis that owns it rather than by a historical stage number. The
# domain and child are checked against the contract, because a typo here would
# silently create a sibling namespace that no reader would ever find.
canonical_result_path <- function(domain, analysis_id, scope = "global",
                                  child = "tables", ...) {
  domains <- output_layout_domains()
  children <- output_layout_children()
  if (length(domains) && !domain %in% domains) {
    stop("unknown output domain '", domain, "'; config/output_layout.yml declares: ",
         paste(domains, collapse = ", "), call. = FALSE)
  }
  if (length(children) && !child %in% children) {
    stop("unknown result child '", child, "'; config/output_layout.yml declares: ",
         paste(children, collapse = ", "), call. = FALSE)
  }
  analysis_id <- sub("[.][Rr]$", "", basename(analysis_id))
  if (!nzchar(scope)) scope <- "global"
  path_results(domain, analysis_id, scope, child, ...)
}

canonical_work_path <- function(domain, analysis_id, scope = "global", ...) {
  analysis_id <- sub("[.][Rr]$", "", basename(analysis_id))
  if (!nzchar(scope)) scope <- "global"
  path_work(domain, analysis_id, scope, ...)
}

# The normalized replacement for create_module_dirs(), which keyed output on
# stage identity. One call gives a writer every destination it needs, so
# repointing a writer is one edit rather than one per write site.
#
# `suffix` is a segment below the child, for a writer whose historical layout
# nested one: build_spatial_networks wrote
# network_spatial_relations/<dataset>/<spatial_unit>, and the spatial unit
# carries meaning, so it is preserved.
#
# figures and logs are aliases, because renaming a key to plots or manifests is
# a pure rename with no judgement in it.
#
# There is deliberately no `processed` alias. Under the historical layout every
# writer put both its persistent objects and its scratch files under
# data/processed, and those two now diverge: a consumed object belongs in
# models/ and a disposable one in work/. Aliasing the old key would silently
# pick one, and that choice is exactly what section 2 of the Phase 6G.2 brief
# requires a human to make per artefact. A writer must therefore name $models
# or $work.
canonical_module_dirs <- function(domain, analysis_id, scope = "global",
                                  suffix = NULL, create = TRUE) {
  aid <- sub("[.][Rr]$", "", basename(analysis_id))
  if (!nzchar(scope)) scope <- "global"
  join <- function(p) if (is.null(suffix) || !length(suffix) || !nzchar(suffix)) p else file.path(p, suffix)

  out <- list(
    tables    = join(canonical_result_path(domain, aid, scope, "tables")),
    plots     = join(canonical_result_path(domain, aid, scope, "plots")),
    models    = join(canonical_result_path(domain, aid, scope, "models")),
    manifests = join(canonical_result_path(domain, aid, scope, "manifests")),
    reports   = join(canonical_result_path(domain, aid, scope, "reports")),
    work      = join(canonical_work_path(domain, aid, scope))
  )
  out$source_data <- file.path(out$tables, "source_data")
  out$figures <- out$plots
  out$logs <- out$manifests

  if (isTRUE(create)) invisible(lapply(out, dir_create))
  out
}

# --- legacy output roots --------------------------------------------------
# Registered in config/legacy_output_registry.csv: real artefacts that no
# registered writer produces any more. Reads are fine, writes are not.
legacy_output_roots <- function() {
  f <- repo_path("config", "legacy_output_registry.csv")
  if (!file.exists(f)) return(character(0))
  d <- utils::read.csv(f, stringsAsFactors = FALSE)
  d$legacy_path[d$policy == "LEGACY_READ_ONLY"]
}

is_legacy_output_path <- function(path) {
  roots <- legacy_output_roots()
  if (!length(roots)) return(rep(FALSE, length(path)))
  rel <- relative_to(path)
  vapply(rel, function(p) any(p == roots | startsWith(p, paste0(roots, "/"))),
         logical(1), USE.NAMES = FALSE)
}

dir_create <- function(...) {
  path <- file.path(...)
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

ensure_dir <- dir_create

safe_filename <- function(x, max_chars = 180) {
  x <- as.character(x)
  x <- gsub("[/\\\\:*?\"<>|]+", "_", x)
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x <- ifelse(nzchar(x), x, "unnamed")
  substr(x, 1, max_chars)
}

file_hash <- function(path) {
  if (is.null(path) || !length(path) || is.na(path) || !file.exists(path)) return(NA_character_)
  unname(tools::md5sum(path))
}

file_hash_sha256 <- function(path) {
  if (
    is.null(path) ||
      !length(path) ||
      is.na(path) ||
      !file.exists(path) ||
      dir.exists(path)
  ) {
    return(NA_character_)
  }

  if ("sha256sum" %in% getNamespaceExports("tools")) {
    return(unname(tools::sha256sum(path)))
  }

  if (requireNamespace("digest", quietly = TRUE)) {
    return(unname(digest::digest(
      file = path,
      algo = "sha256"
    )))
  }

  warning(
    "SHA-256 hashing is unavailable. Install the R package 'digest'.",
    call. = FALSE
  )

  NA_character_
}

strict_inputs_enabled <- function(config = NULL) {
  args <- c(commandArgs(trailingOnly = FALSE), commandArgs(trailingOnly = TRUE))
  from_args <- "--strict-inputs" %in% args
  from_env <- tolower(Sys.getenv("PROTEOMICS_STRICT_INPUTS", unset = "")) %in% c("1", "true", "yes")
  from_config <- FALSE
  if (is.list(config)) {
    from_config <- isTRUE(config$strict_inputs) ||
      isTRUE(config$runtime$strict_inputs)
  }
  isTRUE(from_args || from_env || from_config)
}

input_resolution_audit_path <- function() {
  path_results("reviewer_audit", "input_resolution_audit.csv")
}

input_resolution_audit_columns <- function() {
  c(
    "script", "dataset", "stage", "input_name", "expected_path", "resolved_path",
    "resolution_mode", "strict_mode", "allowed_in_strict_mode", "file_exists",
    "file_hash_sha256", "file_mtime", "producer_script_or_artifact_id", "warning"
  )
}

input_file_mtime <- function(path) {
  if (is.null(path) || !length(path) || is.na(path) || !file.exists(path)) return(NA_character_)
  format(file.info(path)$mtime[[1]], "%Y-%m-%d %H:%M:%S %z")
}

append_input_resolution_audit <- function(rows, path = input_resolution_audit_path()) {
  if (is.null(rows) || !length(rows)) return(invisible(path))
  rows <- as.data.frame(rows, stringsAsFactors = FALSE)
  cols <- input_resolution_audit_columns()
  for (col in setdiff(cols, names(rows))) rows[[col]] <- NA
  rows <- rows[, cols, drop = FALSE]
  dir_create(dirname(path))
  write_header <- !file.exists(path) || file.info(path)$size == 0
  utils::write.table(
    rows,
    file = path,
    sep = ",",
    row.names = FALSE,
    col.names = write_header,
    append = !write_header,
    na = "",
    qmethod = "double"
  )
  invisible(path)
}

record_input_resolution <- function(
  script = Sys.getenv("PROTEOMICS_SCRIPT_ID", unset = NA_character_),
  dataset = Sys.getenv("PROTEOMICS_DATASET", unset = NA_character_),
  stage = NA_character_,
  input_name,
  expected_path = NA_character_,
  resolved_path = NA_character_,
  resolution_mode = "canonical",
  strict_mode = strict_inputs_enabled(),
  allowed_in_strict_mode = TRUE,
  producer_script_or_artifact_id = NA_character_,
  warning = NA_character_
) {
  one_path <- function(path) {
    path <- as.character(path)
    path <- path[!is.na(path) & nzchar(path)]
    if (!length(path)) return(NA_character_)
    normalizePath(path[[1]], winslash = "/", mustWork = FALSE)
  }
  expected_path <- one_path(expected_path)
  resolved_path <- one_path(resolved_path)
  append_input_resolution_audit(data.frame(
    script = script,
    dataset = dataset,
    stage = stage,
    input_name = input_name,
    expected_path = expected_path,
    resolved_path = resolved_path,
    resolution_mode = resolution_mode,
    strict_mode = isTRUE(strict_mode),
    allowed_in_strict_mode = isTRUE(allowed_in_strict_mode),
    file_exists = !is.na(resolved_path) && file.exists(resolved_path),
    file_hash_sha256 = file_hash_sha256(resolved_path),
    file_mtime = input_file_mtime(resolved_path),
    producer_script_or_artifact_id = producer_script_or_artifact_id,
    warning = warning,
    stringsAsFactors = FALSE
  ))
}

latest_input_candidate <- function(roots, pattern, recursive = TRUE) {
  roots <- roots[!is.na(roots) & nzchar(roots)]
  if (!length(roots) || is.na(pattern) || !nzchar(pattern)) return(NA_character_)
  files <- unlist(lapply(roots, function(root) {
    root <- normalizePath(root, winslash = "/", mustWork = FALSE)
    if (!dir.exists(root)) return(character())
    list.files(root, pattern = pattern, full.names = TRUE, recursive = recursive)
  }), use.names = FALSE)
  files <- files[file.exists(files)]
  if (!length(files)) return(NA_character_)
  info <- file.info(files)
  normalizePath(rownames(info)[order(info$mtime, decreasing = TRUE)[1]], winslash = "/", mustWork = FALSE)
}

resolve_input_path <- function(
  input_name,
  expected_path = NA_character_,
  explicit_path = NA_character_,
  fallback_paths = character(),
  latest_roots = character(),
  latest_pattern = NA_character_,
  recursive = TRUE,
  required = TRUE,
  script = Sys.getenv("PROTEOMICS_SCRIPT_ID", unset = NA_character_),
  dataset = Sys.getenv("PROTEOMICS_DATASET", unset = NA_character_),
  stage = NA_character_,
  producer_script_or_artifact_id = NA_character_,
  allow_fallback_in_strict = FALSE,
  allow_latest_in_strict = FALSE,
  record_resolution = TRUE
) {
  strict <- strict_inputs_enabled()
  norm <- function(x) {
    x <- x[!is.na(x) & nzchar(x)]
    if (!length(x)) return(character())
    normalizePath(x, winslash = "/", mustWork = FALSE)
  }
  expected_path <- norm(expected_path)[1]
  if (is.na(expected_path)) expected_path <- NA_character_
  explicit_path <- norm(explicit_path)[1]
  if (is.na(explicit_path)) explicit_path <- NA_character_
  fallback_paths <- norm(fallback_paths)

  finish <- function(resolved, mode, allowed, warn = NA_character_) {
    if (!is.na(warn) && nzchar(warn)) warning(warn, call. = FALSE)
    if (isTRUE(record_resolution)) {
      record_input_resolution(
        script = script,
        dataset = dataset,
        stage = stage,
        input_name = input_name,
        expected_path = expected_path,
        resolved_path = resolved,
        resolution_mode = mode,
        strict_mode = strict,
        allowed_in_strict_mode = allowed,
        producer_script_or_artifact_id = producer_script_or_artifact_id,
        warning = warn
      )
    }
    resolved
  }

  if (!is.na(explicit_path) && file.exists(explicit_path)) {
    return(finish(explicit_path, "explicit_override", TRUE))
  }
  if (!is.na(explicit_path) && !file.exists(explicit_path)) {
    warn <- paste0("Explicit input override does not exist for ", input_name, ": ", explicit_path)
    finish(explicit_path, "explicit_missing", TRUE, warn)
    if (isTRUE(required)) stop(warn, call. = FALSE)
    return(explicit_path)
  }
  if (!is.na(expected_path) && file.exists(expected_path)) {
    return(finish(expected_path, "canonical", TRUE))
  }

  fallback_hit <- fallback_paths[file.exists(fallback_paths)][1]
  if (!is.na(fallback_hit)) {
    warn <- paste0("Using non-canonical fallback for ", input_name, ": ", fallback_hit)
    if (isTRUE(strict) && !isTRUE(allow_fallback_in_strict)) {
      warn <- paste0("Strict input mode forbids fallback for ", input_name, ". Expected canonical input: ", expected_path)
      finish(fallback_hit, "fallback_forbidden_strict", FALSE, warn)
      if (isTRUE(required)) stop(warn, call. = FALSE)
      return(NA_character_)
    }
    return(finish(fallback_hit, "fallback", isTRUE(allow_fallback_in_strict), warn))
  }

  latest_hit <- latest_input_candidate(latest_roots, latest_pattern, recursive = recursive)
  if (!is.na(latest_hit)) {
    warn <- paste0("Using newest matching fallback for ", input_name, ": ", latest_hit)
    if (isTRUE(strict) && !isTRUE(allow_latest_in_strict)) {
      warn <- paste0("Strict input mode forbids newest-file fallback for ", input_name, ". Expected canonical input: ", expected_path)
      finish(latest_hit, "latest_forbidden_strict", FALSE, warn)
      if (isTRUE(required)) stop(warn, call. = FALSE)
      return(NA_character_)
    }
    return(finish(latest_hit, "latest_fallback", isTRUE(allow_latest_in_strict), warn))
  }

  warn <- paste0("Missing input for ", input_name, if (!is.na(expected_path)) paste0(": ", expected_path) else ".")
  finish(NA_character_, "missing", TRUE, if (isTRUE(required)) warn else NA_character_)
  if (isTRUE(required)) stop(warn, call. = FALSE)
  NA_character_
}

path_or_env <- function(env, default, must_exist = FALSE, kind = c("file", "dir", "any")) {
  kind <- match.arg(kind)
  value <- Sys.getenv(env, unset = "")
  path <- if (nzchar(value)) value else default
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)

  if (isTRUE(must_exist)) {
    exists_ok <- switch(
      kind,
      file = file.exists(path),
      dir = dir.exists(path),
      any = file.exists(path) || dir.exists(path)
    )
    if (!exists_ok) {
      stop(
        "Required ", kind, " path does not exist: ", path,
        ". Set ", env, " to override this location.",
        call. = FALSE
      )
    }
  }

  path
}

relative_to <- function(path, root = repo_root()) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  sub(paste0("^", gsub("([\\^$.|?*+(){}])", "\\\\\\1", root), "/?"), "", path)
}

write_session_info <- function(path = path_results("logs", "sessionInfo.txt")) {
  dir_create(dirname(path))
  capture.output(utils::sessionInfo(), file = path)
  invisible(path)
}

git_commit_sha <- function() {
  sha <- tryCatch(
    system2("git", c("-C", repo_root(), "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(e) NA_character_
  )
  sha <- if (length(sha)) as.character(sha[[1]]) else NA_character_
  if (is.na(sha) || !nzchar(sha)) NA_character_ else sha
}

run_context_metadata <- function() {
  env_flags <- c(
    "PROTEOMICS_DATASET",
    "PROTEOMICS_DRY_RUN",
    "PROTEOMICS_STRICT_INPUTS",
    "PROTEOMICS_RECOMPUTE",
    "PROTEOMICS_WGCNA_FORCE_FULL"
  )
  env <- as.list(Sys.getenv(env_flags, unset = NA_character_))
  names(env) <- env_flags
  list(
    git_commit = git_commit_sha(),
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    platform = R.version$platform,
    environment = env
  )
}

write_config_snapshot <- function(config, path) {
  dir_create(dirname(path))
  if (requireNamespace("yaml", quietly = TRUE)) {
    writeLines(yaml::as.yaml(config), path)
  } else {
    capture.output(str(config), file = path)
  }
  invisible(path)
}

# Provenance for one canonical result family, written into its manifests/
# child: config/output_layout.yml section 9.
#
# biological unit and statistical scope are deliberately not restated here.
# They are declared once, for the analyses that reach the manuscript, in
# docs/MANUSCRIPT_STATISTICAL_CONTRACT.md, which is generated from the frozen
# v9 contract. A manifest that paraphrased them would duplicate a frozen
# scientific artefact with no authority to do so, so it points at it instead.
write_result_manifest <- function(domain, analysis_id, scope = "global",
                                  inputs = list(), outputs = list(),
                                  parameters = list(), config_files = character(0),
                                  notes = NULL) {
  aid <- sub("[.][Rr]$", "", basename(analysis_id))
  path <- canonical_result_path(domain, aid, scope, "manifests", "manifest.yml")
  dir_create(dirname(path))

  hash_map <- function(x) {
    x <- as.character(unlist(x, use.names = TRUE))
    x <- x[!is.na(x) & nzchar(x)]
    if (!length(x)) return(list())
    stats::setNames(lapply(x, function(p)
      list(path = relative_to(p),
           sha256 = if (file.exists(p) && !dir.exists(p)) file_hash_sha256(p) else NA_character_)),
      names(x) %||% basename(x))
  }

  owner <- NA_character_
  own_file <- repo_path("config", "results_ownership.csv")
  if (file.exists(own_file)) {
    o <- utils::read.csv(own_file, stringsAsFactors = FALSE)
    hit <- grep(paste0("/", aid, "[.][Rr]$"), o$canonical_owner)
    if (length(hit)) owner <- o$canonical_owner[hit[1]]
  }

  manifest <- list(
    contract_version = "result_manifest_v1",
    analysis_id = aid,
    domain = domain,
    scope = scope,
    canonical_owner = if (is.na(owner)) paste0("analysis/", domain, "/", aid, ".R") else owner,
    source_commit = git_commit_sha(),
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    config = hash_map(config_files),
    upstream_inputs = hash_map(inputs),
    outputs = hash_map(outputs),
    parameters = parameters,
    statistical_contract = "docs/MANUSCRIPT_STATISTICAL_CONTRACT.md",
    notes = notes
  )
  if (requireNamespace("yaml", quietly = TRUE)) {
    writeLines(yaml::as.yaml(manifest), path, useBytes = TRUE)
  }
  invisible(path)
}

write_run_manifest <- function(path, inputs = list(), outputs = list(), parameters = list(), notes = NULL) {
  dir_create(dirname(path))

  flatten_paths <- function(x) {
    vals <- unlist(x, recursive = TRUE, use.names = TRUE)
    vals <- vals[!is.na(vals) & nzchar(as.character(vals))]
    as.character(vals)
  }

  input_paths <- flatten_paths(inputs)
  input_hashes <- lapply(input_paths, file_hash)
  names(input_hashes) <- names(input_paths)

  session_path <- file.path(dirname(path), "sessionInfo.txt")
  capture.output(utils::sessionInfo(), file = session_path)

  manifest <- list(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    repo_root = repo_root(),
    git_commit = git_commit_sha(),
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    platform = R.version$platform,
    environment = run_context_metadata()$environment,
    inputs = inputs,
    input_hashes = input_hashes,
    outputs = outputs,
    parameters = parameters,
    notes = notes,
    session_info = relative_to(session_path)
  )

  if (requireNamespace("yaml", quietly = TRUE)) {
    writeLines(yaml::as.yaml(manifest), path)
  } else {
    capture.output(str(manifest, max.level = 4), file = path)
  }

  invisible(path)
}

is_dry_run <- function(config = NULL) {
  from_args <- "--dry-run" %in% commandArgs(trailingOnly = FALSE) ||
    "--dry-run" %in% commandArgs(trailingOnly = TRUE)
  from_env <- tolower(Sys.getenv("PROTEOMICS_DRY_RUN", unset = "")) %in% c("1", "true", "yes")
  from_config <- FALSE
  if (is.list(config)) {
    from_config <- isTRUE(config$dry_run) ||
      isTRUE(config$runtime$dry_run)
  }
  isTRUE(from_args || from_env || from_config)
}

dry_run_line <- function(label, value = "", status = NULL) {
  prefix <- if (is.null(status)) "[DRY-RUN]" else paste0("[DRY-RUN ", status, "]")
  message(prefix, " ", label, if (nzchar(as.character(value))) paste0(": ", value) else "")
}

module_paths <- function(module, substep = NULL) {
  tail <- c(module, substep)
  tail <- tail[nzchar(tail)]
  list(
    processed = do.call(path_processed, as.list(tail)),
    figures = do.call(path_results, as.list(c("figures", tail))),
    tables = do.call(path_results, as.list(c("tables", tail))),
    source_data = do.call(path_results, as.list(c("source_data", tail))),
    logs = do.call(path_results, as.list(c("logs", tail))),
    reports = do.call(path_results, as.list(c("reports", tail)))
  )
}

create_module_dirs <- function(module, substep = NULL) {
  paths <- module_paths(module, substep)
  invisible(lapply(paths, dir_create))
  paths
}

pride_submission_dir <- function(...) repo_path("pride_submission", ...)

ensure_pride_dirs <- function() {
  dirs <- pride_submission_dir(c(
    "metadata",
    "processed_data",
    "supplementary_tables",
    "methods",
    "manifests",
    "validation"
  ))
  invisible(lapply(dirs, dir_create))
  invisible(dirs)
}

# Backward-compatible alias used by earlier PRIDE helper scripts.
pride_package_dir <- function() pride_submission_dir()

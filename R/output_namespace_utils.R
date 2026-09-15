# Output-namespace contracts and path classification.
#
# These helpers make the current saving structure explicit without moving or
# rewriting historical scientific artifacts.

if (!exists("repo_path", mode = "function")) {
  paths_file <- if (file.exists(file.path("R", "paths.R"))) {
    file.path("R", "paths.R")
  } else {
    file.path("..", "R", "paths.R")
  }
  source(paths_file)
}

output_namespace_contract_path <- function() {
  repo_path("config", "output_namespaces.yml")
}

validate_output_namespace_contract <- function(contract) {
  required <- c(
    "contract_version", "analytical_roots", "manuscript_authoring_roots",
    "manuscript_export_root", "legacy_manuscript_authoring_roots",
    "historical_or_failed_export_prefixes", "diagnostic_roots",
    "comparison_markers", "rules"
  )
  missing <- setdiff(required, names(contract))
  if (length(missing)) {
    stop(
      "Output namespace contract is missing field(s): ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  expected_kinds <- c("figures", "tables", "source_data", "reports", "logs")
  missing_kinds <- setdiff(expected_kinds, names(contract$analytical_roots))
  if (length(missing_kinds)) {
    stop(
      "Output namespace analytical roots are missing: ",
      paste(missing_kinds, collapse = ", "), call. = FALSE
    )
  }
  authoring_kinds <- c("figures", "source_data", "reports", "logs")
  missing_authoring <- setdiff(
    authoring_kinds, names(contract$manuscript_authoring_roots)
  )
  if (length(missing_authoring)) {
    stop(
      "Output namespace manuscript-authoring roots are missing: ",
      paste(missing_authoring, collapse = ", "), call. = FALSE
    )
  }
  scalar_paths <- c(
    unlist(contract$analytical_roots, use.names = FALSE),
    unlist(contract$manuscript_authoring_roots, use.names = FALSE),
    as.character(contract$manuscript_export_root),
    as.character(unlist(contract$legacy_manuscript_authoring_roots, use.names = FALSE)),
    as.character(unlist(contract$historical_or_failed_export_prefixes, use.names = FALSE)),
    as.character(unlist(contract$diagnostic_roots, use.names = FALSE))
  )
  if (anyNA(scalar_paths) || any(!nzchar(scalar_paths))) {
    stop("Output namespace contract contains an empty path.", call. = FALSE)
  }
  if (anyDuplicated(unlist(contract$manuscript_authoring_roots, use.names = FALSE))) {
    stop("Manuscript-authoring roots must be unique.", call. = FALSE)
  }
  rules <- unlist(contract$rules, use.names = FALSE)
  if (!length(rules) || any(!vapply(rules, isTRUE, logical(1)))) {
    stop("Every output namespace safety rule must be enabled.", call. = FALSE)
  }
  invisible(contract)
}

read_output_namespace_contract <- function(
    path = output_namespace_contract_path()) {
  if (!file.exists(path)) {
    stop("Output namespace contract not found: ", path, call. = FALSE)
  }
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to read the output namespace contract.", call. = FALSE)
  }
  contract <- yaml::read_yaml(path)
  validate_output_namespace_contract(contract)
  contract$contract_path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  contract
}

output_namespace_normalize_relative <- function(paths, root = repo_root()) {
  paths <- gsub("\\\\", "/", as.character(paths))
  root <- gsub("\\\\", "/", normalizePath(root, winslash = "/", mustWork = FALSE))
  absolute <- grepl("^(?:[A-Za-z]:/|//)", paths, perl = TRUE)
  if (any(absolute)) {
    normalized <- normalizePath(paths[absolute], winslash = "/", mustWork = FALSE)
    prefix <- paste0(
      "^", gsub("([\\^$.|?*+(){}])", "\\\\\\1", root), "/?"
    )
    paths[absolute] <- sub(prefix, "", normalized, perl = TRUE)
  }
  paths <- sub("^\\./", "", paths)
  sub("^/+", "", paths)
}

output_namespace_has_prefix <- function(paths, prefixes) {
  prefixes <- tolower(gsub("\\\\", "/", as.character(prefixes)))
  paths <- tolower(paths)
  vapply(paths, function(path) {
    any(vapply(prefixes, function(prefix) {
      startsWith(path, prefix) && (
        endsWith(prefix, "/") || nchar(path) == nchar(prefix) ||
          substr(path, nchar(prefix) + 1L, nchar(prefix) + 1L) == "/"
      )
    }, logical(1)))
  }, logical(1))
}

classify_output_namespace <- function(
    paths, contract = read_output_namespace_contract(), root = repo_root()) {
  relative <- output_namespace_normalize_relative(paths, root = root)
  lower <- tolower(relative)
  status <- rep("outside_managed_outputs", length(relative))

  results_path <- startsWith(lower, "results/") | lower == "results"
  status[results_path] <- "unclassified_results"

  historical_prefixes <- tolower(gsub(
    "\\\\", "/", as.character(unlist(
      contract$historical_or_failed_export_prefixes, use.names = FALSE
    ))
  ))
  historical <- vapply(lower, function(path) {
    any(startsWith(path, historical_prefixes))
  }, logical(1))
  status[historical] <- "historical_or_failed_export"

  authoring <- output_namespace_has_prefix(
    lower, unlist(contract$manuscript_authoring_roots, use.names = FALSE)
  )
  status[authoring] <- "manuscript_authoring"

  legacy_authoring <- output_namespace_has_prefix(
    lower, contract$legacy_manuscript_authoring_roots
  )
  status[legacy_authoring] <- "legacy_manuscript_authoring"

  manuscript_export <- output_namespace_has_prefix(
    lower, contract$manuscript_export_root
  ) & !historical
  status[manuscript_export] <- "manuscript_export"

  diagnostic <- output_namespace_has_prefix(lower, contract$diagnostic_roots) |
    grepl("(^|/)(diagnostic|diagnostics|audit)(_|/|$)", lower, perl = TRUE)
  status[diagnostic & results_path & !historical & !authoring &
           !legacy_authoring & !manuscript_export] <- "diagnostic_or_audit"

  comparison_markers <- tolower(as.character(unlist(
    contract$comparison_markers, use.names = FALSE
  )))
  comparison <- vapply(lower, function(path) {
    any(vapply(comparison_markers, function(marker) {
      grepl(marker, path, fixed = TRUE)
    }, logical(1)))
  }, logical(1))
  status[comparison & results_path & !historical & !authoring &
           !legacy_authoring & !manuscript_export] <- "comparison_or_candidate"

  analytical <- output_namespace_has_prefix(
    lower, unlist(contract$analytical_roots, use.names = FALSE)
  )
  status[analytical & status == "unclassified_results"] <-
    "canonical_stage_output"

  data_processed <- output_namespace_has_prefix(lower, "data/processed")
  status[data_processed] <- "processed_data"
  pride <- output_namespace_has_prefix(lower, "pride_submission")
  status[pride] <- "pride_export"
  configuration <- output_namespace_has_prefix(lower, "config")
  status[configuration] <- "source_control_configuration"

  data.frame(
    path = as.character(paths),
    repository_relative_path = relative,
    namespace = status,
    stringsAsFactors = FALSE
  )
}

# The manuscript figure IDs this repository recognises, declared once and shared
# with the export router in R/export_helpers.R. This is an explicit allow-list,
# not a pattern: an unrecognised figure number must still fail closed, because a
# typo silently creating results/figures/manuscript/figure_07 is exactly the
# failure this validation exists to prevent.
#
# Figure 01 is included because pipeline.yml declares results/manuscript/figure_1
# as an export destination. Before it was listed here that slot was unreachable -
# declared by the pipeline and rejected by the validator - so a Figure 1 renderer
# could not have written anywhere legal. No Figure 1 renderer exists yet; this
# makes the namespace valid ahead of one.
MANUSCRIPT_FIGURE_IDS <- c("01", "02", "03")

output_namespace_manuscript_figure_paths <- function(
    output_root, figure_id) {
  figure_id <- suppressWarnings(as.integer(figure_id))
  figure_id <- if (length(figure_id) != 1L || is.na(figure_id)) {
    NA_character_
  } else {
    sprintf("%02d", figure_id)
  }
  if (is.na(figure_id) || !figure_id %in% MANUSCRIPT_FIGURE_IDS) {
    stop("Manuscript figure ID must be one of ",
         paste(MANUSCRIPT_FIGURE_IDS, collapse = ", "), ".", call. = FALSE)
  }
  figure_stub <- paste0("figure_", figure_id)
  list(
    figures = file.path(output_root, "figures", "manuscript", figure_stub),
    panels = file.path(
      output_root, "figures", "manuscript", figure_stub, "panels"
    ),
    assembled = file.path(
      output_root, "figures", "manuscript", figure_stub, "assembled"
    ),
    source_data = file.path(
      output_root, "source_data", "manuscript", figure_stub
    ),
    reports = file.path(
      output_root, "reports", "manuscript_figures", figure_stub
    ),
    logs = file.path(
      output_root, "logs", "manuscript_figures", figure_stub
    )
  )
}

output_namespace_manuscript_export_root <- function(
    results_root = path_results(), contract = read_output_namespace_contract()) {
  configured <- as.character(contract$manuscript_export_root)
  configured_tail <- sub("^results/?", "", configured)
  file.path(results_root, configured_tail)
}

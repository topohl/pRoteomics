#!/usr/bin/env Rscript

# Regenerate registry-derived documentation without changing execution order.

paths_file <- if (file.exists(file.path("R", "paths.R"))) {
  file.path("R", "paths.R")
} else {
  file.path("..", "R", "paths.R")
}
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "pipeline_registry.R"))

registry <- read_pipeline_registry(repo_path("pipeline.yml"))
entries <- pipeline_registry_entries(registry)

run_order_path <- repo_path("RUN_ORDER.md")
run_order <- readLines(run_order_path, warn = FALSE)
begin_marker <- "<!-- BEGIN GENERATED PIPELINE REGISTRY INDEX -->"
end_marker <- "<!-- END GENERATED PIPELINE REGISTRY INDEX -->"
begin <- which(run_order == begin_marker)
end <- which(run_order == end_marker)
if (length(begin) != 1L || length(end) != 1L || begin >= end) {
  stop("RUN_ORDER.md has invalid generated registry index markers.", call. = FALSE)
}
index_lines <- c(
  begin_marker,
  "",
  "## Authoritative script index (generated)",
  "",
  "Do not hand-edit this block; run `Rscript tools/generate_pipeline_docs.R`.",
  "",
  sprintf(
    "%d. `%s` - stage `%s`; scope `%s`",
    seq_len(nrow(entries)), entries$script, entries$stage, entries$scope
  ),
  "",
  end_marker
)
run_order <- c(
  if (begin > 1L) run_order[seq_len(begin - 1L)] else character(),
  index_lines,
  if (end < length(run_order)) run_order[seq.int(end + 1L, length(run_order))] else character()
)
writeLines(run_order, run_order_path, useBytes = TRUE)

steps <- pipeline_steps(
  registry, pipeline_stage_names(registry), dataset = "all",
  include_unsupported = TRUE
)
steps <- steps[!duplicated(steps$script), , drop = FALSE]
active_status <- ifelse(steps$required, "active_required", "active_optional")
active <- data.frame(
  script_path = steps$script,
  status = active_status,
  consumes = paste0(
    "required: ", steps$consumes_required,
    "; optional: ", steps$consumes_optional
  ),
  produces = steps$produces,
  uses_manifest_config_env = paste0(
    "stage=", steps$stage, "; scope=", steps$scope,
    "; datasets=", steps$supported_datasets,
    ifelse(nzchar(steps$env), paste0("; env=", steps$env), "")
  ),
  canonical_output_directory = steps$produces,
  remaining_TODOs = steps$notes,
  stringsAsFactors = FALSE
)

legacy <- do.call(rbind, lapply(registry$legacy, function(item) {
  status_text <- as.character(item$status)
  is_deprecated <- grepl(
    "deprecated|legacy|older", status_text, ignore.case = TRUE
  )
  data.frame(
    script_path = as.character(item$script),
    status = if (is_deprecated) "legacy_deprecated" else "registry_excluded",
    consumes = paste0("replacement: ", as.character(item$replacement)),
    produces = "",
    uses_manifest_config_env = "excluded from active pipeline stages",
    canonical_output_directory = "not canonical",
    remaining_TODOs = paste0(
      status_text, "; replacement: ", as.character(item$replacement)
    ),
    stringsAsFactors = FALSE
  )
}))

note <- data.frame(
  script_path = "AUDIT_GENERATED_SNAPSHOT",
  status = "documentation_note",
  consumes = "Generated from pipeline.yml; pipeline.yml is the active source of truth.",
  produces = "Refresh with Rscript tools/generate_pipeline_docs.R after structural changes.",
  uses_manifest_config_env = "not applicable",
  canonical_output_directory = "not canonical",
  remaining_TODOs = "Generated/audit snapshot; do not treat as active workflow instructions.",
  stringsAsFactors = FALSE
)
utils::write.table(
  rbind(note, active, legacy),
  repo_path("docs", "active_script_io_audit.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

validate_run_order_against_registry(registry, run_order_path)
message("Regenerated RUN_ORDER.md registry index and active script I/O audit.")

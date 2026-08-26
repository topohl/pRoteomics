#!/usr/bin/env Rscript

# Read-only by default. Pass --output <path> to write the row-level impact CSV.

paths_file <- if (file.exists(file.path("R", "paths.R"))) {
  file.path("R", "paths.R")
} else {
  file.path("..", "R", "paths.R")
}
source(paths_file)
source(repo_path("R", "wgcna_support_status_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
output_index <- match("--output", args)
output_path <- if (!is.na(output_index)) {
  if (output_index == length(args)) {
    stop("--output requires a path.", call. = FALSE)
  }
  args[[output_index + 1L]]
} else {
  NA_character_
}

group_effect_files <- Sys.glob(file.path(
  repo_root(), "results", "tables", "06_modules_WGCNA", "group_effects",
  "*", "*group_effects.csv"
))
if (!length(group_effect_files)) {
  stop("No Stage 05 group-effect tables were found.", call. = FALSE)
}

audits <- lapply(group_effect_files, function(path) {
  data <- utils::read.csv(
    path, check.names = FALSE, stringsAsFactors = FALSE,
    na.strings = c("", "NA", "NaN", "NULL")
  )
  wgcna_group_support_status_impact_audit(
    data,
    gsub("\\\\", "/", substring(path, nchar(repo_root()) + 2L))
  )
})
audit <- do.call(rbind, audits)
rownames(audit) <- NULL

if (nrow(audit)) {
  audit$source_key <- paste(
    paste0("dataset=", audit$dataset),
    paste0("level=", audit$level),
    paste0("endpoint_id=", audit$endpoint_id),
    paste0("effect_scope=", audit$effect_scope),
    paste0("spatial_unit=", audit$spatial_unit),
    paste0("contrast=", audit$contrast),
    paste0("test_type=", audit$test_type),
    sep = "|"
  )
}

handoff_files <- Sys.glob(file.path(
  repo_root(), "results", "tables", "06_modules_WGCNA",
  "interpretable_summary", "*", "WGCNA_inferential_handoff.csv"
))
handoff_keys <- unique(unlist(lapply(handoff_files, function(path) {
  data <- utils::read.csv(
    path, check.names = FALSE, stringsAsFactors = FALSE,
    na.strings = c("", "NA", "NaN", "NULL")
  )
  if ("source_key" %in% names(data)) as.character(data$source_key) else character()
}), use.names = FALSE))
audit$stage07_handoff_consumed <- audit$source_key %in% handoff_keys

consumer_files <- c(
  file.path(
    repo_root(), "results", "source_data", "10_biological_integration",
    "wgcna_circular_atlas", "global",
    c(
      "wgcna_circular_atlas_plot_source.csv",
      "wgcna_circular_heatmap_source_module.csv",
      "wgcna_circular_heatmap_source_supermodule.csv",
      "wgcna_circular_global_supermodule_support_source.csv"
    )
  ),
  file.path(
    repo_root(), "results", "tables", "10_biological_integration",
    "gsea_wgcna_concordance", "global", "gsea_wgcna_concordance_long.csv"
  )
)
consumer_files <- consumer_files[file.exists(consumer_files)]
usage <- lapply(consumer_files, function(path) {
  data <- utils::read.csv(
    path, check.names = FALSE, stringsAsFactors = FALSE,
    na.strings = c("", "NA", "NaN", "NULL")
  )
  key_column <- intersect(c("source_key", "wgcna_source_key"), names(data))
  if (!length(key_column)) return(NULL)
  keys <- unique(as.character(data[[key_column[[1]]]]))
  matched <- audit$source_key %in% keys
  if (!any(matched)) return(NULL)
  data.frame(
    source_key = audit$source_key[matched],
    consumer_artifact = gsub(
      "\\\\", "/", substring(path, nchar(repo_root()) + 2L)
    ),
    stringsAsFactors = FALSE
  )
})
usage <- Filter(Negate(is.null), usage)
usage <- if (length(usage)) do.call(rbind, usage) else data.frame(
  source_key = character(), consumer_artifact = character()
)
audit$manuscript_consumer_artifacts <- vapply(audit$source_key, function(key) {
  paste(unique(usage$consumer_artifact[usage$source_key == key]), collapse = ";")
}, character(1))

cat("Rows changing classification: ", nrow(audit), "\n", sep = "")
if (nrow(audit)) {
  print(with(audit, table(dataset, level, old_status, new_status)))
  cat(
    "Rows present in Stage 07 handoffs: ",
    sum(audit$stage07_handoff_consumed), "\n", sep = ""
  )
  cat(
    "Rows present in scanned manuscript consumer artifacts: ",
    sum(nzchar(audit$manuscript_consumer_artifacts)), "\n", sep = ""
  )
  print(audit, row.names = FALSE)
}

if (!is.na(output_path)) {
  output_path <- normalizePath(output_path, winslash = "/", mustWork = FALSE)
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(audit, output_path, row.names = FALSE, na = "")
  message("Wrote impact audit: ", output_path)
}

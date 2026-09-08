#!/usr/bin/env Rscript

# Read-only audit of the repository's declared and top-level output namespaces.

paths_file <- if (file.exists(file.path("R", "paths.R"))) {
  file.path("R", "paths.R")
} else {
  file.path("..", "R", "paths.R")
}
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "pipeline_registry.R"))
source(repo_path("R", "output_namespace_utils.R"))

contract <- read_output_namespace_contract()
registry <- read_pipeline_registry(repo_path("pipeline.yml"))

declared_rows <- list()
for (stage_name in names(registry$stages)) {
  scripts <- registry$stages[[stage_name]]$scripts
  for (step in scripts) {
    produces <- as.character(unlist(step$produces %||% character(), use.names = FALSE))
    if (!length(produces)) next
    declared_rows[[length(declared_rows) + 1L]] <- data.frame(
      stage = stage_name,
      script = as.character(step$script),
      output = produces,
      stringsAsFactors = FALSE
    )
  }
}
declared <- do.call(rbind, declared_rows)
declared <- cbind(
  declared,
  classify_output_namespace(declared$output)[, c(
    "repository_relative_path", "namespace"
  ), drop = FALSE]
)

figure_steps <- startsWith(declared$script, "figures/")
figure_outputs_ok <- all(
  declared$namespace[figure_steps] == "manuscript_authoring"
)

cat("Output namespace contract: ", contract$contract_version, "\n", sep = "")
cat("Contract path: ", relative_to(contract$contract_path), "\n", sep = "")
cat("Declared pipeline outputs: ", nrow(declared), "\n", sep = "")
cat("Declared namespace summary:\n")
print(sort(table(declared$namespace), decreasing = TRUE))
cat(
  "Manuscript entry-point output isolation: ",
  if (figure_outputs_ok) "PASS" else "FAIL", "\n", sep = ""
)

manuscript_root <- output_namespace_manuscript_export_root()
if (dir.exists(manuscript_root)) {
  children <- list.files(
    manuscript_root, full.names = TRUE, recursive = FALSE,
    all.files = TRUE, no.. = TRUE
  )
  child_audit <- classify_output_namespace(children)
  child_audit$kind <- ifelse(dir.exists(children), "directory", "file")
  cat("Existing manuscript-export children:\n")
  print(child_audit[, c(
    "repository_relative_path", "namespace", "kind"
  ), drop = FALSE], row.names = FALSE)
  historical <- child_audit$namespace == "historical_or_failed_export"
  if (any(historical)) {
    cat(
      "INFO: historical/failed export trees are present and remain untouched: ",
      paste(child_audit$repository_relative_path[historical], collapse = "; "),
      "\n", sep = ""
    )
  }
}

if (!figure_outputs_ok) {
  bad <- declared[figure_steps & declared$namespace != "manuscript_authoring", , drop = FALSE]
  print(bad, row.names = FALSE)
  stop(
    "A manuscript entry point declares output outside manuscript-authoring roots.",
    call. = FALSE
  )
}

cat("Output namespace audit completed without writing files.\n")

# Audit active R scripts for path/I-O reproducibility hazards.
#
# Run from repository root:
#   Rscript tools/audit_active_scripts.R

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)

pipeline_file <- repo_path("pipeline.yml")
pipeline_lines <- readLines(pipeline_file, warn = FALSE)

active_script_rows <- function(lines) {
  rows <- data.frame(stage = character(), script = character(), stringsAsFactors = FALSE)
  in_stages <- FALSE
  stage <- NA_character_
  current_script <- NA_character_
  for (line in lines) {
    if (grepl("^stages:", line)) in_stages <- TRUE
    if (grepl("^legacy:", line)) in_stages <- FALSE
    if (!in_stages) next
    stage_match <- regexec("^  ([A-Za-z0-9_]+):\\s*$", line)
    stage_hit <- regmatches(line, stage_match)[[1]]
    if (length(stage_hit)) stage <- stage_hit[[2]]
    script_match <- regexec("^      - script: \"([^\"]+)\"", line)
    script_hit <- regmatches(line, script_match)[[1]]
    if (length(script_hit)) current_script <- script_hit[[2]]
    if (grepl("^        scope: \"", line) && !is.na(current_script)) {
      rows <- rbind(rows, data.frame(stage = stage, script = current_script, stringsAsFactors = FALSE))
      current_script <- NA_character_
    }
  }
  rows
}

active_registry <- active_script_rows(pipeline_lines)
rel_files <- active_registry$script
active_files <- repo_path(rel_files)

patterns <- list(
  "S:/Lab_Member/Tobi" = "S:/Lab_Member/Tobi",
  "C:/" = "C:/",
  "/Users/" = "/Users/",
  "setwd(" = "setwd\\s*\\(",
  "install.packages(" = "install\\.packages\\s*\\(",
  "pacman::p_load(" = "pacman::p_load\\s*\\(",
  "hard-coded absolute Windows path" = "(^|[^A-Za-z])[A-Za-z]:[/\\\\]",
  "hard-coded absolute macOS path" = "/Users/"
)

deprecated_roots <- c(
  "Results/",
  "Output/",
  "Final/",
  "Final2/",
  "publication_ready/",
  "Plots/",
  "Datasets/",
  "results/module_scores/",
  "path_results(\"module_scores\""
)

write_context <- paste(c(
  "write",
  "save",
  "dir.create",
  "ensure_dir",
  "dir_create",
  "file.copy",
  "file.path",
  "path_results",
  "out_dir",
  "out_file",
  "output"
), collapse = "|")

allowed_deprecated_context <- function(line) {
  grepl("#|legacy|fallback|compat|deprecated|migration|warning\\s*\\(", line, ignore.case = TRUE, perl = TRUE)
}

hits <- list()
add_hit <- function(file, line, pattern, context) {
  hits[[length(hits) + 1L]] <<- data.frame(
    file = file,
    line = line,
    pattern = pattern,
    context = trimws(context),
    stringsAsFactors = FALSE
  )
}

for (i in seq_along(active_files)) {
  file <- active_files[[i]]
  rel <- rel_files[[i]]
  if (!file.exists(file)) {
    add_hit(rel, NA_integer_, "active script missing", "Listed in pipeline.yml but file is absent")
    next
  }
  lines <- readLines(file, warn = FALSE)
  header_block <- paste(utils::head(lines, 60), collapse = "\n")
  required_header_fields <- c("Script:", "Stage:", "Scope:", "Consumes:", "Produces:", "Dataset behavior:", "Notes:")
  for (field in required_header_fields) {
    if (!grepl(field, header_block, fixed = TRUE)) {
      add_hit(rel, 1L, paste0("missing I/O header field: ", field), "Add compact script-level I/O header near top of active script")
    }
  }

  for (nm in names(patterns)) {
    idx <- grep(patterns[[nm]], lines, perl = TRUE)
    if (length(idx) > 0) {
      for (j in idx) add_hit(rel, j, nm, lines[[j]])
    }
  }

  for (root in deprecated_roots) {
    idx <- grep(root, lines, fixed = TRUE)
    if (length(idx) > 0) {
      for (j in idx) {
        line <- lines[[j]]
        if (grepl(write_context, line, ignore.case = TRUE, perl = TRUE) && !allowed_deprecated_context(line)) {
          add_hit(rel, j, paste0("deprecated output root write: ", root), line)
        }
      }
    }
  }

  sink_idx <- grep("sink\\s*\\(", lines, perl = TRUE)
  if (length(sink_idx) > 0) {
    for (j in sink_idx) {
      window <- lines[j:min(length(lines), j + 5L)]
      if (!any(grepl("on\\.exit\\s*\\(", window))) {
        add_hit(rel, j, "sink() without nearby on.exit cleanup", lines[[j]])
      }
    }
  }

  helper_names <- c("ensure_dir", "save_plot", "save_table")
  for (helper in helper_names) {
    idx <- grep(paste0("^\\s*", helper, "\\s*<-\\s*function\\b"), lines, perl = TRUE)
    if (length(idx) > 1L) {
      for (j in idx) add_hit(rel, j, paste0("duplicate helper definition: ", helper), lines[[j]])
    }
  }
}

audit <- if (length(hits) > 0) do.call(rbind, hits) else data.frame(
  file = character(), line = integer(), pattern = character(), context = character()
)

if (nrow(audit) > 0) {
  print(audit, row.names = FALSE)
  stop("Active script audit failed with ", nrow(audit), " finding(s).", call. = FALSE)
}

message("Active script audit passed: no blocked patterns found.")

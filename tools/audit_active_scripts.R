# Audit active R scripts for path/I-O reproducibility hazards.
#
# Run from repository root:
#   Rscript tools/audit_active_scripts.R

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)

active_files <- list.files(repo_root(), pattern = "\\.[rR]$", recursive = TRUE, full.names = TRUE)
rel_files <- relative_to(active_files)
excluded <- grepl("(^|/)90_testing/", rel_files) |
  grepl("(^|/)99_deprecated/", rel_files) |
  grepl("^tools/", rel_files) |
  rel_files == "05_celltype_enrichment_EWCE/90_EWCE_legacy.r"
active_files <- active_files[!excluded]
rel_files <- rel_files[!excluded]

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
  lines <- readLines(file, warn = FALSE)

  for (nm in names(patterns)) {
    idx <- grep(patterns[[nm]], lines, perl = TRUE)
    if (length(idx) > 0) {
      for (j in idx) add_hit(rel, j, nm, lines[[j]])
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

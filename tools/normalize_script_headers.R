#!/usr/bin/env Rscript

# Bring every canonical analysis entrypoint up to the script-header contract
# that tools/audit_active_scripts.R enforces:
#
#   Script: Stage: Scope: Consumes: Produces: Dataset behavior: Notes:
#
# Six of the seven fields are derived from pipeline.yml, so they cannot state
# anything the registry does not already state. `Notes:` is taken from prose
# the file already carries; where a file has none, it gets a factual line about
# its registration rather than an invented description of its biology.
#
# Why these field names and not `Inputs:`/`Outputs:`/`Biological unit:`:
# Consumes/Produces already mean inputs/outputs and are machine-audited on 77
# files, so renaming them would be vocabulary churn. Biological unit and
# statistical scope are scientific statements; they are declared once, for the
# 11 analyses that reach the manuscript, in docs/MANUSCRIPT_STATISTICAL_CONTRACT.md,
# which is generated from the frozen v9 contract. Restating them per script
# would duplicate a frozen scientific artefact with no authority to do so, so
# the header points at that contract instead.
#
# Run with DRY=1 to preview.

source(file.path("R", "paths.R"))
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "pipeline_registry.R"))

DRY <- nzchar(Sys.getenv("DRY"))

FIELDS <- c("Script:", "Stage:", "Scope:", "Consumes:", "Produces:",
            "Dataset behavior:", "Notes:")

# These two are byte-compared against historical commits by
# freeze_protected_export_files(); never rewrite them.
NEVER_TOUCH <- c(
  "analysis/publication_source_data/08_export_manuscript_figures.R",
  "analysis/publication_source_data/09_export_source_data.R"
)

registry <- read_pipeline_registry(repo_path("pipeline.yml"))
steps <- pipeline_steps(registry, pipeline_stage_names(registry),
                        dataset = "all", include_unsupported = TRUE)

split_paths <- function(x) {
  p <- trimws(unlist(strsplit(paste(x, collapse = "|"), "[|;]")))
  unique(p[nzchar(p)])
}
summarise_paths <- function(x, n = 3L) {
  x <- split_paths(x)
  if (!length(x)) return("none declared in pipeline.yml")
  out <- paste(utils::head(x, n), collapse = "; ")
  if (length(x) > n) out <- paste0(out, "; +", length(x) - n, " more")
  out
}

# the file's own first descriptive comment, if it has one
own_prose <- function(lines) {
  cand <- character(0)
  for (x in utils::head(lines, 40)) {
    if (grepl("^#!", x)) next                        # shebang is not prose
    if (!grepl("^\\s*#", x)) {
      if (nzchar(trimws(x))) break else next
    }
    t <- trimws(sub("^\\s*#+", "", x))
    if (!nzchar(t)) next
    if (grepl("^[=-]{3,}$", t)) next                 # rule lines
    if (grepl("^[A-Z][A-Za-z ]{2,24}:", t)) next     # an existing field
    if (nchar(t) < 12) next
    cand <- c(cand, t)
  }
  if (!length(cand)) return(NA_character_)
  cand[1]
}

leading_block_end <- function(lines) {
  last <- 0L
  for (i in seq_along(lines)) {
    if (grepl("^\\s*#", lines[i])) { last <- i; next }
    if (!nzchar(trimws(lines[i]))) next
    break
  }
  last
}

scripts <- unique(steps$script)
scripts <- scripts[grepl("^analysis/", scripts) & file.exists(scripts)]
scripts <- setdiff(scripts, NEVER_TOUCH)

changed <- character(0)
added <- 0L
for (f in scripts) {
  lines <- readLines(f, warn = FALSE)
  hb <- paste(utils::head(lines, 60), collapse = "\n")
  missing <- FIELDS[!vapply(FIELDS, function(fd) grepl(fd, hb, fixed = TRUE), logical(1))]
  if (!length(missing)) next

  s <- steps[steps$script == f, , drop = FALSE]
  ds <- unique(s$dataset[s$supported %in% TRUE])
  if (!length(ds)) ds <- unique(s$dataset)
  prose <- own_prose(lines)

  value <- c(
    "Script:" = f,
    "Stage:" = paste(unique(s$stage), collapse = ", "),
    "Scope:" = paste(unique(s$scope), collapse = ", "),
    "Consumes:" = paste0("required ", summarise_paths(s$consumes_required),
                         "; optional ", summarise_paths(s$consumes_optional)),
    "Produces:" = summarise_paths(s$produces, 3L),
    "Dataset behavior:" = paste0("runs for ", paste(ds, collapse = ","),
      " according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported."),
    "Notes:" = if (!is.na(prose)) prose else paste0(
      "Registered in pipeline.yml stage ", paste(unique(s$stage), collapse = "/"),
      "; declares ", length(split_paths(s$produces)), " output path(s).")
  )
  new_lines <- paste0("# ", names(value)[FIELDS %in% missing], " ",
                      value[FIELDS %in% missing])

  at <- leading_block_end(lines)
  if (at == 0L) {
    # no leading comment block: open one at the very top, after any shebang
    top <- if (length(lines) && grepl("^#!", lines[1])) 1L else 0L
    block <- c("# ================================================================",
               new_lines,
               "# ================================================================",
               "")
    out <- append(lines, block, after = top)
  } else {
    # extend the existing leading comment block, which contains only comments
    # and blank lines, so nothing can land inside a string literal
    out <- append(lines, new_lines, after = at)
  }

  changed <- c(changed, f)
  added <- added + length(new_lines)
  if (!DRY) {
    con <- file(f, open = "wb")
    writeLines(out, con, sep = "\n")
    close(con)
  }
}

cat(if (DRY) "WOULD update: " else "updated: ", length(changed), " entrypoints\n", sep = "")
cat("header field lines added:", added, "\n")
cat("freeze-protected files skipped:", length(NEVER_TOUCH), "\n")
if (length(changed)) cat("\n", paste0("  ", utils::head(changed, 10), collapse = "\n"), "\n")

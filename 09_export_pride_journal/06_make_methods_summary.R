#!/usr/bin/env Rscript
# Methods/provenance summary for pg_matrix-onward reproducibility (not raw-MS search).

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "export_helpers.R"))
source(repo_path("R", "pipeline_registry.R"))

ensure_pride_dirs()
config <- load_export_config()
cli <- export_cli_args()
methods_dir <- pride_submission_dir("methods")

collect_session_info_files <- function() {
  logs_root <- path_results("logs")
  if (!dir.exists(logs_root)) return(character(0))
  list.files(logs_root, pattern = "^sessionInfo\\.txt$", recursive = TRUE, full.names = TRUE)
}

collect_software_versions <- function() {
  rows <- data.frame(package = character(), version = character(), source = character(), stringsAsFactors = FALSE)
  lockfile <- repo_path("renv.lock")
  if (file.exists(lockfile) && requireNamespace("jsonlite", quietly = TRUE)) {
    lock <- jsonlite::fromJSON(lockfile)
    pkgs <- lock$Packages
    if (!is.null(pkgs) && length(pkgs)) {
      nm <- names(pkgs)
      rows <- rbind(rows, data.frame(
        package = nm,
        version = vapply(pkgs, function(x) as.character(x$Version %||% NA_character_), character(1)),
        source = "renv.lock",
        stringsAsFactors = FALSE
      ))
    }
  }
  si_files <- collect_session_info_files()
  if (length(si_files)) {
    rows <- rbind(rows, data.frame(package = "R_session", version = basename(dirname(si_files[[1]])), source = si_files[[1]], stringsAsFactors = FALSE))
  }
  rows
}

build_pipeline_steps <- function() {
  registry <- tryCatch(read_pipeline_registry(repo_path("pipeline.yml")), error = function(e) NULL)
  if (is.null(registry)) return(data.frame(stage = character(), script = character(), stringsAsFactors = FALSE))
  stages <- pipeline_stage_names(registry)
  out <- do.call(rbind, lapply(stages, function(st) {
    steps <- pipeline_steps(registry, st)
    if (!nrow(steps)) return(NULL)
    data.frame(stage = st, script = steps$script, required = steps$required, stringsAsFactors = FALSE)
  }))
  if (is.null(out)) data.frame(stage = character(), script = character(), stringsAsFactors = FALSE) else out
}

if (isTRUE(cli$dry_run)) {
  dry_run_line("Script", "09_export_pride_journal/06_make_methods_summary.R")
  dry_run_line("Methods summary target", file.path(methods_dir, "methods_summary.md"))
  dry_run_line("Software versions target", file.path(methods_dir, "software_versions.tsv"))
  dry_run_line("Pipeline steps target", file.path(methods_dir, "pipeline_steps.tsv"))
  dry_run_line("Known limitations target", file.path(methods_dir, "known_limitations.md"))
  quit(status = 0, save = "no")
}

summary_path <- file.path(methods_dir, "methods_summary.md")
software_path <- file.path(methods_dir, "software_versions.tsv")
pipeline_path <- file.path(methods_dir, "pipeline_steps.tsv")
limitations_path <- file.path(methods_dir, "known_limitations.md")

lines <- c(
  "# Proteomics Methods and Provenance Summary",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  "",
  "## Reproducibility entry point",
  "",
  config$scope_note,
  "",
  paste0("Export level: **", cli$export_level %||% config$export_level, "**"),
  "",
  "## Analysis layers (in-repo)",
  "",
  "1. **pg_matrix input** — protein-group quantification matrix (`data/raw/pg_matrix/`).",
  "2. **Preprocessing** — imputation, dataset splits, ProTigy/GCT handoff (`data/processed/01_preprocessing/`).",
  "3. **Identifier mapping** — protein-to-gene mapping for contrasts (`data/processed/02_id_mapping/`).",
  "4. **Differential abundance & enrichment** — clusterProfiler/compareGO and downstream summaries.",
  "5. **Modules, networks, behavior coupling** — results under `results/`.",
  "",
  "## External / optional materials",
  "",
  "- Raw/vendor MS files, search-engine outputs, FASTA, and search parameters may be deposited separately.",
  "- They are **not** required to rerun the committed workflow from pg_matrix onward.",
  "",
  "## Configuration",
  "",
  paste0("- Export config: `", relative_to(export_config_path(), repo_root()), "`"),
  paste0("- Pipeline registry: `pipeline.yml`"),
  ""
)
writeLines(lines, summary_path)

sw <- collect_software_versions()
if (nrow(sw)) write_tsv(sw, software_path) else warning("No renv.lock or sessionInfo found for software_versions.tsv", call. = FALSE)

ps <- build_pipeline_steps()
if (nrow(ps)) utils::write.table(ps, pipeline_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

limit_lines <- c(
  "# Known limitations",
  "",
  "## pg_matrix-onward reproducibility",
  "",
  "This repository does **not** reconstruct raw mass-spectrometry acquisition, instrument methods,",
  "search-engine execution, or FASTA/search-parameter choices. Scientific reproducibility of the",
  "committed analysis begins at the **pg_matrix** quantification stage.",
  "",
  "## PRIDE partial / processed deposition",
  "",
  "A PRIDE submission from this module is intended as **processed-data** and metadata supporting",
  "journal reproducibility, optionally complemented by externally archived raw/search files.",
  "",
  "## Dataset interpretation",
  "",
  paste0("- **microglia**: ", config$interpretation_scope$microglia),
  paste0("- **neuron_neuropil**: ", config$interpretation_scope$neuron_neuropil),
  paste0("- **neuron_soma**: ", config$interpretation_scope$neuron_soma),
  ""
)
writeLines(limit_lines, limitations_path)

message("Wrote methods summary: ", summary_path)
if (file.exists(software_path)) message("Wrote software versions: ", software_path)
if (file.exists(pipeline_path)) message("Wrote pipeline steps: ", pipeline_path)
message("Wrote known limitations: ", limitations_path)

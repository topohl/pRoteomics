# Sample-centring contract: the one setting that separates the two analysis versions (config/normalisation.yml).
#
# sample_centring_offsets() is the only implementation of the centring; archive/01_preprocessing/01_impute.r applies
# it and stamps the version it wrote next to the imputed workbooks. The launcher compares that stamp with the setting
# and with the inputs the downstream readers resolve (imputed workbooks, newest metadata-joined workbook, animal-level
# input, the input of the limma fit); run manifests record setting and stamp. Freeze and export on the pRoteomics side
# refuse to run until a decision for one centring is approved (the Exp9_manuscript import tools are not guarded).
# Outside the check: the legacy derived trees (results/tables, data/processed/02-10, gct_extractR, protigy_output),
# which hold original-centring results and can still be reached by hard-coded or fallback reads, and the frozen export
# bundle (exports/, results/manuscript, pride_submission), which holds the original version.

NORMALISATION_METHODS <- c(
  complete_case_median = "median of the proteins quantified in every sample of the compartment",
  all_observed_median = "median of all observed values of the sample"
)
IMPUTED_WORKBOOK_PATTERN <- "^[0-9]{8}_pgmatrix_imputed_.+_[0-9]+samples_missing70pct[.]xlsx$"
# Environment variables that redirect an analysis input (recorded in manifests: such a run is not described by the stamp)
INPUT_OVERRIDE_PATTERN <- paste0("^PROTEOMICS_(IMPUTE_|EXCEL_CONVERT_|PROTIGY_INPUT_|WGCNA_UPSTREAM_|MODULE_SCORE_|EWCE_MATRIX|GCT_|",
                                 "JOINT_QC_PROCESSED|PROJECT_ROOT|SAMPLE_CENTRING|MODULE_DEFINITION_SOURCE)")
COMPARTMENTS <- c("microglia", "neuron_neuropil", "neuron_soma")

normalisation_setting <- function() {
  cfg_path <- repo_path("config", "normalisation.yml")
  if (!file.exists(cfg_path)) stop("Missing ", cfg_path, call. = FALSE)
  cfg <- yaml::read_yaml(cfg_path)
  env <- Sys.getenv("PROTEOMICS_SAMPLE_CENTRING", unset = "")
  method <- if (nzchar(env)) env else as.character(cfg$sample_centring %||% "")
  if (!method %in% names(NORMALISATION_METHODS)) {
    stop("Unknown sample_centring '", method, "'; use one of: ", paste(names(NORMALISATION_METHODS), collapse = ", "), call. = FALSE)
  }
  list(sample_centring = method, description = unname(NORMALISATION_METHODS[[method]]),
       decision = as.character(cfg$decision %||% "pending"),
       approved_sample_centring = as.character(cfg$approved_sample_centring %||% ""),
       source = if (nzchar(env)) "PROTEOMICS_SAMPLE_CENTRING" else "config/normalisation.yml")
}

# mat: proteins x samples on the log2 scale, NA where a protein was not quantified. Returns one offset per sample.
sample_centring_offsets <- function(mat, method) {
  mat <- as.matrix(mat)
  if (identical(method, "complete_case_median")) {
    complete <- stats::complete.cases(mat)
    if (!any(complete)) stop("No protein is quantified in every sample; complete-case centring is undefined.", call. = FALSE)
    out <- apply(mat[complete, , drop = FALSE], 2L, stats::median)
    attr(out, "n_reference_proteins") <- sum(complete)
  } else if (identical(method, "all_observed_median")) {
    out <- apply(mat, 2L, stats::median, na.rm = TRUE)
    attr(out, "n_reference_proteins") <- NA_integer_
  } else {
    stop("Unknown sample_centring: ", method, call. = FALSE)
  }
  out
}

normalisation_stamp_path <- function(dir = path_processed("01_preprocessing", "impute")) {
  file.path(dir, "normalisation_stamp.yml")
}

git_tree_dirty <- function() {
  out <- tryCatch(suppressWarnings(system2("git", c("-C", shQuote(repo_root()), "status", "--porcelain", "--untracked-files=no"),
                                           stdout = TRUE, stderr = FALSE)), error = function(e) NA_character_)
  if (length(out) == 1L && is.na(out)) NA else length(out) > 0L
}

write_normalisation_stamp <- function(dir, setting, outputs, note = NULL, extra = list()) {
  stamp <- c(list(
    contract_version = "sample_centring_stamp_v1",
    sample_centring = setting$sample_centring,
    description = setting$description,
    setting_source = setting$source,
    written_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    git_commit = git_commit_sha(),
    git_uncommitted_changes = git_tree_dirty(),
    outputs = lapply(outputs, function(p) list(file = basename(p), sha256 = file_hash_sha256(p)))), extra, list(note = note))
  yaml::write_yaml(stamp, normalisation_stamp_path(dir))
  invisible(normalisation_stamp_path(dir))
}

# Imputed workbooks written before the stamp existed carry the centring only in imputation_qc.csv (columns added on
# 2026-10-19/21); workbooks without them were written with the original centring. Such a stamp is marked inferred.
read_normalisation_stamp <- function(dir = path_processed("01_preprocessing", "impute")) {
  p <- normalisation_stamp_path(dir)
  if (file.exists(p)) return(yaml::read_yaml(p))
  qc <- file.path(dir, "imputation_qc.csv")
  if (!file.exists(qc)) return(NULL)
  q <- utils::read.csv(qc, stringsAsFactors = FALSE)
  method <- if (nrow(q) && "sample_centring" %in% names(q) && length(unique(q$sample_centring)) == 1L) q$sample_centring[1] else
    if (!nrow(q)) NA_character_ else if (!"centring" %in% names(q)) "all_observed_median" else
    if (all(grepl("quantified in all samples", q$centring, fixed = TRUE))) "complete_case_median" else NA_character_
  list(contract_version = "inferred_from_imputation_qc", sample_centring = method, written_at = NA_character_, inferred = TRUE)
}

# The workbook a downstream reader would resolve for a compartment: the newest match (that is what the resolvers do).
newest_matching <- function(roots, pattern, recursive) {
  f <- unlist(lapply(roots[dir.exists(roots)], function(r) list.files(r, pattern = pattern, full.names = TRUE, recursive = recursive)))
  if (!length(f)) return(NA_character_)
  f[order(file.info(f)$mtime, decreasing = TRUE)][1]
}

assert_normalisation_consistent <- function(context = "pipeline") {
  s <- normalisation_setting(); st <- read_normalisation_stamp()
  if (is.null(st) || is.na(st$sample_centring %||% NA_character_)) {
    stop(context, ": the imputed workbooks carry no normalisation stamp; rerun archive/01_preprocessing/01_impute.r.", call. = FALSE)
  }
  imp_dir <- path_processed("01_preprocessing", "impute")
  on_disk <- list.files(imp_dir, pattern = IMPUTED_WORKBOOK_PATTERN)
  stamped <- vapply(st$outputs %||% list(), function(o) o$file, character(1))
  fail <- function(...) stop(context, ": ", ..., call. = FALSE)
  if (length(stamped)) {
    # The resolver takes the newest matching workbook, so the directory must hold exactly the stamped workbooks, unmodified.
    foreign <- setdiff(on_disk, stamped); absent <- setdiff(stamped, on_disk)
    changed <- Filter(function(o) o$file %in% on_disk && !identical(file_hash_sha256(file.path(imp_dir, o$file)), o$sha256), st$outputs)
    if (length(foreign) || length(absent) || length(changed)) {
      fail("the imputed workbooks do not match the normalisation stamp (not stamped: ", paste(c(foreign, "none")[seq_len(max(1, length(foreign)))], collapse = ", "),
           "; missing: ", paste(c(absent, "none")[seq_len(max(1, length(absent)))], collapse = ", "),
           "; modified: ", paste(c(vapply(changed, function(o) o$file, ""), "none")[seq_len(max(1, length(changed)))], collapse = ", "), ").")
    }
  } else {
    # An inferred stamp records no files: accept it only when the choice of workbook is unambiguous.
    per_cmp <- vapply(COMPARTMENTS, function(ds) sum(grepl(paste0("_pgmatrix_imputed_", ds, "_"), on_disk, fixed = TRUE)), integer(1))
    if (any(per_cmp != 1L)) fail("the stamp is inferred from imputation_qc.csv and the imputed workbooks are ambiguous (",
                                 paste(sprintf("%s: %d", names(per_cmp), per_cmp), collapse = ", "), "); write a stamp.")
    stamped <- on_disk
  }
  for (ds in COMPARTMENTS) {
    mine <- stamped[grepl(paste0("_pgmatrix_imputed_", ds, "_"), stamped, fixed = TRUE)]
    if (!length(mine)) next
    # Metadata-joined workbook as the module-score resolver finds it (same roots, recursive, newest first).
    joined <- newest_matching(c(path_processed("01_preprocessing"), path_processed("morpheus")),
                              paste0("^[0-9]{8}_pgmatrix_imputed_", ds, "_[0-9]+samples_missing70pct_with_metadata[.]xlsx$"), recursive = TRUE)
    if (!is.na(joined) && !identical(sub("_with_metadata[.]xlsx$", ".xlsx", basename(joined)), mine[1])) {
      fail("the newest metadata-joined workbook for ", ds, " (", joined, ") was not made from the stamped workbook ", mine[1], ".")
    }
    # Animal-level input of the differential model: built from the stamped workbook, and the limma fit from that input?
    gct <- path_processed("01_preprocessing", "protigy_input_animal_level", ds, paste0(ds, "_animal_level.gct"))
    summ <- path_results("tables", "01_preprocessing", "02a_prepare_animal_level_protigy_input", ds, "aggregation_summary.csv")
    want <- if (length(st$outputs)) Filter(function(o) identical(o$file, mine[1]), st$outputs)[[1]]$sha256 else file_hash_sha256(file.path(imp_dir, mine[1]))
    if (file.exists(gct) && !file.exists(summ)) fail("the animal-level input for ", ds, " has no aggregation summary to trace it to an imputed workbook; rerun analysis/preprocessing/build_animal_level_protigy_input.R.")
    if (file.exists(summ)) {
      h <- utils::read.csv(summ, stringsAsFactors = FALSE)$source_expression_hash[1]
      if (!identical(h, want)) fail("the animal-level input for ", ds, " was built from another imputed workbook (hash ", substr(h, 1, 12),
                                    "); rerun analysis/preprocessing/build_animal_level_protigy_input.R.")
    }
    fit_man <- path_results("preprocessing", "fit_animal_level_limma", ds, "manifests", "run_manifest.csv")
    if (file.exists(gct) && file.exists(fit_man)) {
      md5 <- utils::read.csv(fit_man, stringsAsFactors = FALSE)$input_md5[1]
      if (!identical(md5, unname(tools::md5sum(gct)))) fail("the differential model for ", ds, " was fitted on another animal-level input; rerun analysis/preprocessing/fit_animal_level_limma.R.")
    }
  }
  if (!identical(st$sample_centring, s$sample_centring)) {
    stop(sprintf(paste0("%s: normalisation mismatch. The setting is %s (from %s) but the imputed workbooks in the working ",
                        "tree were written with %s. Regenerate the inputs (archive/01_preprocessing/01_impute.r, 02_excel_convert.r, ",
                        "analysis/preprocessing/build_animal_level_protigy_input.R) and rerun the pipeline, or restore the ",
                        "setting. Nothing was run."), context, s$sample_centring, s$source, st$sample_centring), call. = FALSE)
  }
  invisible(s)
}

# TRUE when a directory holds a stamped analysis version (or imputed workbooks), whatever path spelling reaches it.
dir_holds_version <- function(dir) {
  file.exists(normalisation_stamp_path(dir)) || length(list.files(dir, pattern = IMPUTED_WORKBOOK_PATTERN)) > 0L
}

# Regenerating a directory under another setting overwrites the version it holds; allowed only once that version has a
# preserved copy listed in config/normalisation_versions.yml whose imputed-workbook hashes equal the ones stamped here.
assert_working_version_preserved <- function(target_method, dir = path_processed("01_preprocessing", "impute")) {
  st <- read_normalisation_stamp(dir)
  if (is.null(st) || identical(st$sample_centring, target_method)) return(invisible(TRUE))
  versions <- yaml::read_yaml(repo_path("config", "normalisation_versions.yml"))$versions
  held <- vapply(st$outputs %||% list(), function(o) o$sha256, character(1))
  if (!length(held)) held <- vapply(list.files(dir, pattern = IMPUTED_WORKBOOK_PATTERN, full.names = TRUE), file_hash_sha256, character(1))
  ok <- vapply(versions, function(v) {
    identical(v$sample_centring, st$sample_centring) && identical(v$status, "preserved") && dir.exists(v$root %||% "") &&
      file.exists(file.path(v$root, v$manifest %||% "SHA256_MANIFEST.csv")) && length(held) > 0L &&
      all(held %in% unlist(v$imputed_workbook_sha256 %||% character()))
  }, logical(1))
  if (!any(ok)) {
    stop("The directory ", dir, " holds the ", st$sample_centring, " version, and config/normalisation_versions.yml lists no ",
         "preserved copy of exactly these workbooks. Preserve it (tools/normalisation_review/snapshot_working_version.py) and ",
         "record the copy before regenerating with ", target_method, ".", call. = FALSE)
  }
  invisible(TRUE)
}

# Freeze/export/import entry points call this first. Only a script with its own side-effect-free dry run
# (dry_run_capable = TRUE) may continue in dry mode; every other script stops while the decision is pending.
assert_normalisation_decided <- function(context, dry_run_capable = FALSE, dry = NULL) {
  if (is.null(dry)) dry <- tolower(Sys.getenv("PROTEOMICS_DRY_RUN", unset = "")) %in% c("1", "true", "yes") ||
    "--dry-run" %in% commandArgs(trailingOnly = TRUE)
  report_only <- isTRUE(dry_run_capable) && isTRUE(dry)
  s <- tryCatch(normalisation_setting(), error = function(e) if (report_only) e else stop(e))
  if (inherits(s, "error")) { message("[BLOCKED] ", context, ": ", conditionMessage(s)); return(invisible(NULL)) }
  problem <- if (!identical(s$decision, "approved")) {
    paste0("config/normalisation.yml has decision '", s$decision, "'. Freeze, export and manuscript import wait for the normalisation decision.")
  } else if (!identical(s$approved_sample_centring, s$sample_centring) || !identical(s$source, "config/normalisation.yml")) {
    paste0("the approved centring (", s$approved_sample_centring, ") is not the active setting (", s$sample_centring, " from ", s$source, ").")
  } else NULL
  if (!is.null(problem)) {
    if (report_only) { message("[BLOCKED] ", context, " is blocked: ", problem); return(invisible(s)) }
    stop(context, " is blocked: ", problem, call. = FALSE)
  }
  if (report_only) return(invisible(tryCatch(assert_normalisation_consistent(context), error = function(e) { message("[BLOCKED] ", conditionMessage(e)); s })))
  assert_normalisation_consistent(context)
}

# Time the analysis version in the working tree was written (its oldest stamped imputed workbook). An exported or
# frozen file older than this was produced under the previous version.
normalisation_version_time <- function() {
  st <- read_normalisation_stamp(); imp_dir <- path_processed("01_preprocessing", "impute")
  f <- file.path(imp_dir, vapply(st$outputs %||% list(), function(o) o$file, character(1)))
  if (!length(f)) f <- list.files(imp_dir, pattern = IMPUTED_WORKBOOK_PATTERN, full.names = TRUE)
  if (!length(f)) stop("No imputed workbooks to date the analysis version.", call. = FALSE)
  min(file.info(f)$mtime)
}

assert_sources_current <- function(paths, context) {
  paths <- paths[file.exists(paths) & !dir.exists(paths)]
  if (!length(paths)) return(invisible(TRUE))
  old <- paths[file.info(paths)$mtime < normalisation_version_time()]
  if (length(old)) {
    stop(context, ": ", length(old), " source file(s) predate the analysis version in the working tree (e.g. ", old[1],
         "); regenerate them under the chosen centring before exporting.", call. = FALSE)
  }
  invisible(TRUE)
}

# Recorded in run and result manifests. Input overrides are listed because a run that reads an overridden input is not
# described by the working-tree stamp.
normalisation_provenance <- function(dir = Sys.getenv("PROTEOMICS_NORMALISATION_STAMP_DIR", unset = path_processed("01_preprocessing", "impute"))) {
  tryCatch({
    s <- normalisation_setting(); st <- read_normalisation_stamp(dir)
    env <- Sys.getenv()
    ov <- names(env)[grepl(INPUT_OVERRIDE_PATTERN, names(env)) & nzchar(env) & !names(env) %in% c("PROTEOMICS_STRICT_INPUTS", "PROTEOMICS_NORMALISATION_STAMP_DIR")]
    list(sample_centring_setting = s$sample_centring, setting_source = s$source, decision = s$decision,
         stamp_directory = relative_to(dir), stamp_sample_centring = st$sample_centring %||% NA_character_,
         stamp_written_at = st$written_at %||% NA_character_,
         input_overrides = if (length(ov)) paste0(ov, "=", unname(env[ov])) else "none")
  }, error = function(e) list(error = conditionMessage(e)))
}

#!/usr/bin/env Rscript

# Phase 6G section 21: where every registered analysis actually writes.
#
# The registry is a declaration; the code is the fact. This audit reads the
# entrypoints and reports which namespace each one resolves its destinations
# into, so "canonical writers still targeting the historical namespace" is a
# number rather than an impression.
#
# Reads are not writes. A keyword scan cannot tell them apart and gets this
# badly wrong: the migrated EWCE entrypoint reads
# path_processed("01_preprocessing", ...), which is preprocessing's output, and
# a keyword scan calls that a legacy write. So the detection works on the
# parse tree:
#
#   1. collect the variables assigned from a legacy path construction;
#   2. find the calls that actually write or create a directory;
#   3. a legacy write is a write call whose argument subtree reaches a legacy
#      construction, or one of those variables.
#
# Comments are absent from a parse tree by construction, so a migrated writer
# naming its historical namespace in a comment cannot be miscounted.

source(file.path("R", "paths.R"))
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "pipeline_registry.R"))

registry <- read_pipeline_registry(repo_path("pipeline.yml"))
steps <- pipeline_steps(registry, pipeline_stage_names(registry),
                        dataset = "all", include_unsupported = TRUE)
steps <- steps[!duplicated(steps$script), , drop = FALSE]
steps <- steps[grepl("^analysis/", steps$script), , drop = FALSE]

split_paths <- function(x) {
  p <- trimws(unlist(strsplit(paste(x, collapse = "|"), "[|;]")))
  unique(p[nzchar(p)])
}

# canonical_module_dirs() is the normalized replacement for
# create_module_dirs(): it returns every destination a writer needs, built
# from canonical_result_path() and canonical_work_path(). A writer that uses
# it resolves through the output contract even though it never names those
# two functions itself.
# spatial_systems_dirs() is spatial_validation's thin wrapper over
# canonical_module_dirs(): it fixes the domain and the global scope, which
# every analysis in that domain shares, and defaults to creating nothing so
# each writer creates only the lifecycle directories it uses. A writer calling
# it resolves through the output contract without naming the lower-level
# helpers itself.
NORMALIZED_CALLS <- c("canonical_result_path", "canonical_work_path",
                      "canonical_module_dirs", "spatial_systems_dirs")
# these create directories, so calling one is itself a write
LEGACY_FACTORIES <- c("create_module_dirs", "module_paths", "qc_paths")
PATH_BUILDERS <- c("path_results", "path_processed")
# A write call is whatever actually puts bytes on disk, which includes the
# repository's own wrappers. Phase 6G.3 found this list too short: three
# spatial_validation scripts write workbooks through
# xlsx_save_valid_workbook(wb, wb_path), where wb_path is
# path_results("tables", "11_spatial_systems", ...). The path construction was
# detected, the write was not, so those legacy writes were invisible and the
# "active legacy write sites = 0" gate would have passed with them still in
# place. A detector that only knows base R is not AST-aware detection of this
# codebase.
#
# The wrappers below were derived, not guessed: every function defined in R/
# whose body reaches write.csv/write.table/writeLines/saveRDS/ggsave/
# saveWorkbook/write.xlsx/file.copy/file.rename/dir.create, and which analysis/
# actually calls. Pure-removal helpers (unlink) are deliberately excluded: a
# delete is not a write to a destination.
BASE_WRITE_CALLS <- c("dir_create", "write.csv", "write.table", "write_csv_safe",
                      "write_csv_safe2", "writeLines", "saveRDS", "ggsave",
                      "saveWorkbook", "write.xlsx", "file.copy", "file.rename",
                      "write_run_manifest", "write_result_manifest", "write_yaml",
                      "png", "pdf", "svg", "jpeg", "tiff", "cairo_pdf")
WRAPPER_WRITE_CALLS <- c(
  "xlsx_save_valid_workbook", "write_input_status", "save_nature_svg",
  "write_sus_res_biological_audit_workbook", "write_config_snapshot",
  "write_csv_strict", "save_plot_dual", "write_gct_v1.3",
  "joint_qc_write_matrix_tsv", "joint_qc_write_gct_v13",
  "write_gct_extract_contract_manifest", "write_tsv", "copy_export_targets",
  "copy_export_file", "write_validation_summary_md", "qc_write_csv",
  "qc_write_xlsx", "qc_save_square_svg", "joint_pub_save_svg",
  "tokenize_wgcna_mouse_only", "wgcna_group_prepare_stage",
  "wgcna_group_atomic_publish")
WRITE_CALLS <- c(BASE_WRITE_CALLS, WRAPPER_WRITE_CALLS)
STAGE_NS <- "^[0-9]{2}[a-z]?_[A-Za-z]"

call_name <- function(e) {
  if (!is.call(e)) return(NA_character_)
  fn <- e[[1]]
  if (is.name(fn)) return(as.character(fn))
  # namespaced calls such as openxlsx::saveWorkbook
  if (is.call(fn) && length(fn) == 3L && as.character(fn[[1]]) %in% c("::", ":::")) {
    return(as.character(fn[[3]]))
  }
  NA_character_
}

## A parse tree can hold the empty symbol: a formal with no default, or an
## omitted index in x[i, ]. It can be assigned to a variable but any reference
## to it raises "argument is missing", including is.null(). So the test is
## wrapped, covers NULL as well, and treats a throw as "skip".
is_skippable <- function(x) {
  tryCatch(is.null(x) || (is.symbol(x) && !nzchar(as.character(x))), error = function(...) TRUE)
}

walk <- function(e, fn) {
  ## Force the argument here, under a guard. R passes it as a promise, so an
  ## empty symbol would otherwise blow up at an arbitrary later point inside
  ## the recursion rather than where it can be skipped.
  if (!tryCatch({ e; TRUE }, error = function(...) FALSE)) return(invisible(NULL))
  fn(e)
  ## a function definition: walk the body, never the formals
  if (is.call(e) && is.name(e[[1]]) && identical(as.character(e[[1]]), "function")) {
    if (length(e) >= 3L) walk(e[[3]], fn)
    return(invisible(NULL))
  }
  if (is.call(e) || is.expression(e) || is.list(e)) {
    for (i in seq_along(e)) {
      el <- tryCatch(e[[i]], error = function(...) NULL)
      ## is.null() would itself force el, so the whole test is guarded: an
      ## empty symbol can be assigned to a variable but not referenced.
      if (is_skippable(el)) next
      walk(el, fn)
    }
  }
  invisible(NULL)
}

# a path construction that lands in the historical namespace
is_legacy_construction <- function(e) {
  nm <- call_name(e)
  if (is.na(nm)) return(FALSE)
  if (nm %in% LEGACY_FACTORIES) return(TRUE)
  if (nm %in% PATH_BUILDERS) {
    args <- as.list(e)[-1]
    lits <- unlist(lapply(args, function(a) if (is.character(a)) a else NULL))
    return(any(grepl(STAGE_NS, lits)))
  }
  FALSE
}

subtree_has <- function(e, pred) {
  found <- FALSE
  walk(e, function(x) if (!found && isTRUE(pred(x))) found <<- TRUE)
  found
}

subtree_names <- function(e) {
  out <- character(0)
  walk(e, function(x) if (is.name(x)) out <<- c(out, as.character(x)))
  unique(out)
}

analyse <- function(f) {
  if (!file.exists(f)) {
    return(list(api = FALSE, legacy = FALSE, n_legacy = 0L, n_api = 0L, parsed = FALSE))
  }
  exprs <- tryCatch(parse(f), error = function(e) NULL)
  if (is.null(exprs)) {
    return(list(api = FALSE, legacy = FALSE, n_legacy = 0L, n_api = 0L, parsed = FALSE))
  }

  ## 1. variables holding a legacy path
  legacy_vars <- character(0)
  for (e in exprs) {
    walk(e, function(x) {
      if (is.call(x) && length(x) == 3L &&
          as.character(x[[1]])[1] %in% c("<-", "=", "<<-") && is.name(x[[2]])) {
        if (subtree_has(x[[3]], is_legacy_construction)) {
          legacy_vars <<- c(legacy_vars, as.character(x[[2]]))
        }
      }
    })
  }
  legacy_vars <- unique(legacy_vars)

  ## 2. and 3. write calls reaching a legacy path
  n_legacy <- 0L
  n_api <- 0L
  for (e in exprs) {
    walk(e, function(x) {
      nm <- call_name(x)
      if (is.na(nm)) return(invisible(NULL))
      if (nm %in% NORMALIZED_CALLS) n_api <<- n_api + 1L
      if (nm %in% LEGACY_FACTORIES) {
        n_legacy <<- n_legacy + 1L                 # creates directories itself
        return(invisible(NULL))
      }
      if (nm %in% WRITE_CALLS) {
        args <- as.list(x)[-1]
        hit <- any(vapply(args, function(a)
          subtree_has(a, is_legacy_construction), logical(1))) ||
          any(legacy_vars %in% unlist(lapply(args, subtree_names)))
        if (hit) n_legacy <<- n_legacy + 1L
      }
      invisible(NULL)
    })
  }
  list(api = n_api > 0L, legacy = n_legacy > 0L,
       n_legacy = n_legacy, n_api = n_api, parsed = TRUE)
}

rows <- lapply(seq_len(nrow(steps)), function(i) {
  f <- steps$script[i]
  domain <- sub("^analysis/([^/]+)/.*", "\\1", f)
  aid <- sub("[.][Rr]$", "", basename(f))
  a <- analyse(f)

  declared <- split_paths(steps$produces[i])
  declared_normalized <- length(declared) > 0 &&
    all(startsWith(declared, paste0("results/", domain, "/")) |
        startsWith(declared, paste0("work/", domain, "/")))

  status <- if (!a$parsed) {
    "UNPARSED"
  } else if (a$api && !a$legacy && declared_normalized) {
    "MIGRATED"
  } else if ((a$api && !a$legacy) != declared_normalized) {
    "PARTIAL"                        # code and declaration disagree
  } else {
    "PENDING"
  }

  data.frame(
    domain = domain,
    analysis_id = aid,
    entrypoint = f,
    canonical_output = paste(utils::head(declared, 3), collapse = " | "),
    n_declared_outputs = length(declared),
    actual_resolved_write_root = if (a$api && !a$legacy) {
      paste0("results/", domain, "/", aid, "/")
    } else if (a$legacy) {
      "historical stage namespace"
    } else {
      "no write destination constructed in this file"
    },
    expected_write_root = paste0("results/", domain, "/", aid, "/"),
    n_layout_api_calls = a$n_api,
    n_legacy_writes = a$n_legacy,
    legacy_write = a$legacy,
    declared_normalized = declared_normalized,
    migration_status = status,
    stringsAsFactors = FALSE)
})
d <- do.call(rbind, rows)
d <- d[order(d$migration_status, d$domain, d$analysis_id), ]

if (!dir.exists("audits")) dir.create("audits")
utils::write.csv(d, "audits/phase6g_writer_namespace_audit.csv", row.names = FALSE)

cat("registered analysis entrypoints:", nrow(d), "\n\n")
cat("=== migration_status ===\n"); print(table(d$migration_status))
cat("\n=== by domain ===\n"); print(table(d$domain, d$migration_status))
cat("\nwriters resolving through the output-layout API:", sum(d$migration_status == "MIGRATED"), "\n")
cat("writers still writing the historical namespace :", sum(d$legacy_write), "\n")
cat("legacy write sites in total                     :", sum(d$n_legacy_writes), "\n")

split <- d[d$migration_status == "PARTIAL", , drop = FALSE]
if (nrow(split)) {
  cat("\nFAIL: code and registry disagree for:\n")
  print(split[, c("entrypoint", "legacy_write", "declared_normalized")], row.names = FALSE)
  stop("split-brain writer: migrate the code and the declaration together",
       call. = FALSE)
}
cat("split-brain writers: 0\n")

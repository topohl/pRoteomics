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
## Phase 6G.8: the canonical-versus-historical rule lives in one library so it
## can be tested on fixtures without running this audit.
source(repo_path("R", "output_namespace_classification.R"))

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
# differential_abundance_dirs() is the same kind of wrapper for that domain,
# and differential_abundance_relative_path() is the repo-relative form used by
# the one writer that supports redirecting its output root.
# qc_dirs() is the QC domain's wrapper, and it replaced qc_paths(), the legacy
# factory eleven of the fifteen QC writers shared.
NORMALIZED_CALLS <- c("canonical_result_path", "canonical_work_path",
                      "canonical_module_dirs", "spatial_systems_dirs",
                      # Phase 6G.9. The publication boundary resolves through
                      # its own thin layer, which delegates to
                      # canonical_module_dirs for production and keeps the
                      # frozen exports root separate.
                      "psd_dirs", "psd_find", "psd_claims_audit",
                      "psd_claims_table", "psd_export_root",
                      "differential_abundance_dirs",
                      "differential_abundance_relative_path",
                      "qc_dirs", "integration_dirs",
                      # Phase 6G.7. Two preprocessing writers take a
                      # configurable output root, so their destination is a
                      # resolver call rather than a fixed dirs() lookup.
                      "preprocessing_dirs",
                      "preprocessing_gct_extract_dir",
                      "preprocessing_gct_extract_manifest_dir",
                      "preprocessing_gct_extract_manifest_path",
                      "preprocessing_mapping_dir",
                      "preprocessing_mapping_result_dir",
                      # Phase 6G.8. One domain-level directory helper plus a
                      # few named artifact helpers; see R/wgcna/wgcna_paths.R.
                      "wgcna_dirs", "wgcna_find", "wgcna_dir_any")
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
                      "png", "pdf", "svg", "jpeg", "tiff", "cairo_pdf",
                      # Phase 6G.8. The list above knew write.csv and
                      # write.xlsx, the base/openxlsx spellings, but not the
                      # readr and writexl ones the codebase actually uses most:
                      # write_csv has 346 call sites in analysis/ and
                      # write_xlsx 37, and every one of them was invisible to
                      # this audit. A destination spec is added for each below,
                      # because a write call with no DEST_ARG entry has all of
                      # its arguments tested and that is how correctly migrated
                      # writes got reported as legacy in Phase 6G.4.
                      "write_csv", "write_xlsx", "write_json", "write_delim",
                      "write_lines", "write_rds", "svglite")
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

# Phase 6G.8: classify a destination from the contract, not from numbering.
#
# The previous rule called a path construction legacy only when one of its
# literal arguments matched a numbered stage namespace. That is a proxy for the
# thing that matters, and it missed a whole family: results/reviewer_audit/
# carries no stage number, so five WGCNA writers wrote there while their
# pipeline.yml declarations already said results/wgcna/..., and the split-brain
# gate still reported 0. A detector that recognizes historical destinations by
# their spelling will keep missing every historical root that is not numbered.
#
# The contract is the authority instead:
#   config/output_layout.yml declares the canonical shape
#     results/<domain>/<analysis_id>/<scope>/<child>/
#   so in canonical code the FIRST argument of path_results() is a declared
#   domain, while in historical code it is a "kind" (tables, figures, logs,
#   source_data, reports, reviewer_audit).
#   config/legacy_output_registry.csv confirms the historical roots
#   independently, and is used when the first argument is computed.
#   config/output_layout.yml declares exactly three lifecycles - work, results
#   and exports - so data/processed is not a canonical output destination at
#   all, and a write through path_processed() is noncanonical by contract.
#
# STAGE_NS is retained as an additional signal, never as the only one, so the
# numbered roots keep being caught when the first argument is not a literal.
STAGE_NS <- "^[0-9]{2}[a-z]?_[A-Za-z]"

OUTPUT_DOMAINS <- output_layout_domains()
OUTPUT_CHILDREN <- output_layout_children()

# Registered historical roots, as their segment under results/.
LEGACY_REGISTERED_SEGMENTS <- local({
  reg <- repo_path("config", "legacy_output_registry.csv")
  if (!file.exists(reg)) return(character(0))
  d <- utils::read.csv(reg, stringsAsFactors = FALSE)
  s <- vapply(strsplit(d$legacy_path, "/", fixed = TRUE),
              function(p) if (length(p) >= 2L) p[[2]] else NA_character_,
              character(1))
  unique(s[!is.na(s)])
})

# Adjudicated destinations outside results/ that are NOT legacy writes:
#   config/           a generated configuration contract, adjudicated in 6G.5
#   exports/          the frozen outward-facing bundle
#   pride_submission/ gitignored export staging
ALLOWED_NONRESULT_ROOTS <- c("config", "exports", "pride_submission")

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
# Only the destination argument decides where a call writes.
#
# The earlier rule tested every argument of a write call, which cannot tell a
# destination from a payload. Phase 6G.4 showed the cost: two correctly
# migrated writes in test_microglia_targeted_signatures were reported as
# legacy writes because
#
#   writeLines(c("...", EMPIRICAL_ROI_MARKER_PATH, ...), PATHS$reports)
#   write_run_manifest(file.path(PATHS$logs, "..."), inputs = list(contrast_dir(...)))
#
# mention historical paths in their *content*: a report that prints which
# input it read, and a manifest that records its inputs as provenance. Both
# write to the normalized namespace. Recording where you read from is the
# opposite of writing there, and a gate that calls it a legacy write would
# block a correct migration.
#
# Positions follow each function's signature; a named argument wins over the
# position. A call absent from this table keeps the conservative behaviour of
# testing every argument, which is right for dir_create() and friends where
# every argument is part of the destination.
DEST_ARG <- list(
  write.csv = list(2L, "file"), write.table = list(2L, "file"),
  write_csv_safe = list(2L, "path"), write_csv_safe2 = list(2L, "path"),
  writeLines = list(2L, "con"), saveRDS = list(2L, "file"),
  ggsave = list(1L, "filename"), saveWorkbook = list(2L, "file"),
  write.xlsx = list(2L, "file"), file.copy = list(2L, "to"),
  file.rename = list(2L, "to"), write_run_manifest = list(1L, "path"),
  write_result_manifest = list(1L, "path"), write_yaml = list(2L, "file"),
  png = list(1L, "filename"), pdf = list(1L, "file"),
  svg = list(1L, "filename"), jpeg = list(1L, "filename"),
  tiff = list(1L, "filename"), cairo_pdf = list(1L, "filename"),
  xlsx_save_valid_workbook = list(2L, "path"),
  write_csv_strict = list(2L, "path"), save_nature_svg = list(2L, "filename"),
  save_plot_dual = list(2L, "path"), qc_write_csv = list(2L, "path"),
  qc_write_xlsx = list(2L, "path"), write_tsv = list(2L, "path"),
  # Phase 6G.8 additions, positions taken from each function's signature.
  write_csv = list(2L, "file"), write_xlsx = list(2L, "path"),
  write_json = list(2L, "path"), write_delim = list(2L, "file"),
  write_lines = list(2L, "file"), write_rds = list(2L, "file"),
  svglite = list(1L, "filename")
)

# Calls that write only when they are handed a destination. cat() has 552 call
# sites in analysis/ and almost all of them print to the console; treating it
# as an unconditional write would make every legacy path named in a message a
# legacy write, which is the false-positive class Phase 6G.4 removed. So these
# count as a write only when the named argument below is actually present, and
# that argument is the only destination considered.
CONDITIONAL_WRITE_CALLS <- c(cat = "file", capture.output = "file", sink = "file")

# The subtrees that determine this call's destination.
dest_args <- function(nm, call) {
  spec <- DEST_ARG[[nm]]
  args <- as.list(call)[-1]
  if (is.null(spec) || !length(args)) return(args)
  nms <- names(args)
  if (!is.null(nms) && spec[[2]] %in% nms) return(args[nms == spec[[2]]])
  ## positional: count only the unnamed arguments
  unnamed <- if (is.null(nms)) args else args[!nzchar(nms)]
  if (length(unnamed) >= spec[[1]]) return(unnamed[spec[[1]]])
  args
}

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

# The canonical-versus-historical rule is R/utilities/output_namespace_classification.R.
# It is delegated rather than duplicated: two copies of this predicate would
# drift, and the whole reason Phase 6G.8 needed a correction is that the rule
# was a proxy nobody could test in isolation.
is_legacy_construction <- function(e) {
  is_legacy_path_construction(e, domains = OUTPUT_DOMAINS,
                              legacy_segments = LEGACY_REGISTERED_SEGMENTS,
                              allowed = ALLOWED_NONRESULT_ROOTS)
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

# The R libraries an entrypoint sources, one level deep.
#
# Some writers construct no destination themselves:
# analysis/qc/render_compartment_abundance_figures.R has no write call at all,
# because R/qc/control_compartment_abundance_workflow.R does the rendering and
# builds the paths. Judging only the entrypoint reports such a writer as
# having no destination, which then looks like a split brain against a
# perfectly good registry declaration. One level is enough for the cases that
# exist and keeps this a reader rather than a link-time resolver.
sourced_libraries <- function(f) {
  exprs <- tryCatch(parse(f), error = function(e) NULL)
  if (is.null(exprs)) return(character(0))
  out <- character(0)
  for (e in exprs) walk(e, function(x) {
    if (!is.call(x)) return(invisible(NULL))
    nm <- call_name(x)
    if (is.na(nm) || !nm %in% c("source", "sys.source")) return(invisible(NULL))
    p <- tryCatch(eval(x[[2]], envir = globalenv()), error = function(...) NA_character_)
    if (length(p) == 1L && !is.na(p) && file.exists(p) && grepl("[.][Rr]$", p)) {
      out <<- c(out, p)
    }
    invisible(NULL)
  })
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
      if (nm %in% WRITE_CALLS || nm %in% names(CONDITIONAL_WRITE_CALLS)) {
        args <- if (nm %in% names(CONDITIONAL_WRITE_CALLS)) {
          ## a write only when the destination argument is actually supplied
          key <- CONDITIONAL_WRITE_CALLS[[nm]]
          a <- as.list(x)[-1]
          nms <- names(a)
          if (is.null(nms) || !key %in% nms) list() else a[nms == key]
        } else {
          dest_args(nm, x)
        }
        if (length(args)) {
          hit <- any(vapply(args, function(a)
            subtree_has(a, is_legacy_construction), logical(1))) ||
            any(legacy_vars %in% unlist(lapply(args, subtree_names)))
          if (hit) n_legacy <<- n_legacy + 1L
        }
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
  ## A generated configuration contract is not a result and is exempt.
  ##
  ## analysis/qc/build_reference_marker_registry.R declares
  ## config/marker_panels/wgcna_reference_marker_sets.csv. That file is
  ## committed, sits beside three hand-maintained marker panels, and is read as
  ## configuration by twelve scripts across five domains through
  ## consumes_required/consumes_optional. The output-layout contract governs
  ## results, not configuration, so requiring it under results/<domain>/ would
  ## be wrong rather than merely strict.
  declared_results <- declared[!startsWith(declared, "config/")]

  ## Phase 6G.9. results/<domain>/ is not the only canonical destination in the
  ## repository, and treating it as such mis-reported the publication boundary
  ## as unmigrated.
  ##
  ## Three further namespaces are declared, not incidental:
  ##   exports/            output_layout.yml, "frozen outward-facing bundles",
  ##                       may_be_canonical: false - a COPY of canonical results
  ##   pride_submission/   an allowed non-result root in the output-namespace
  ##                       classifier; the proteomics repository deposit bundle
  ##   results/manuscript/ output_namespaces.yml manuscript_export_root, with
  ##                       the explicit rule
  ##                       exporters_write_only_to_manuscript_export_root: true
  ##
  ## A writer targeting one of these is obeying a contract, not evading one, so
  ## it counts as normalized. The roots are read from the contracts rather than
  ## hard-coded, so this cannot drift from them. Everything else still has to
  ## sit under results/<domain>/ or work/<domain>/.
  sanctioned_roots <- local({
    r <- c("exports/", "pride_submission/")
    ns <- tryCatch(read_output_namespace_contract(), error = function(e) NULL)
    mer <- if (!is.null(ns) && !is.null(ns$manuscript_export_root)) {
      paste0(sub("/*$", "", as.character(ns$manuscript_export_root)), "/")
    } else "results/manuscript/"
    c(r, mer)
  })
  declared_normalized <- length(declared_results) > 0 &&
    all(startsWith(declared_results, paste0("results/", domain, "/")) |
        startsWith(declared_results, paste0("work/", domain, "/")) |
        Reduce(`|`, lapply(sanctioned_roots, function(p) startsWith(declared_results, p))))

  ## A writer that constructs no destination at all cannot contradict its
  ## declaration, so it is not a split brain.
  ##
  ## analysis/qc/render_compartment_abundance_figures.R has no write call and
  ## no path construction: R/qc/control_compartment_abundance_workflow.R does
  ## both. Calling that PARTIAL would report a disagreement between a
  ## declaration and code that says nothing. It gets its own status so it stays
  ## visible rather than being quietly counted as migrated, and the domain's
  ## test asserts that the library it delegates to resolves normalized.
  constructs_nothing <- a$parsed && a$n_api == 0L && a$n_legacy == 0L

  status <- if (!a$parsed) {
    "UNPARSED"
  } else if (constructs_nothing) {
    if (declared_normalized) "DELEGATED" else "PENDING"
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

#!/usr/bin/env Rscript

# Mechanically enumerate everything that may depend on one analysis's output
# namespace.
#
#   Rscript tools/enumerate_output_consumers.R --analysis-id <id> [--quiet]
#   -> audits/consumer_inventory/<analysis_id>.csv
#
# Why this exists. The enrichment pilot missed
# tests/testthat/test-ewce-dataset-saving.R through three separate searches,
# because that file never names a path: it names the fragment "EWCE_E9". A
# search for the full historical path cannot find it, and neither can a search
# for an output filename. Missing one consumer of one terminal analysis was
# recoverable; missing one of WGCNA's twenty-four writers would not be.
#
# The vocabulary is derived, never hand-listed. Hand-listing is what failed:
# the tokens you forget are exactly the ones you do not search for.
#
# Structural extraction where it pays. R files are parsed, so a string literal
# is attributed to the call that encloses it, which is what separates a read
# from a write from a mention in an expectation. Comments do not appear in a
# parse tree at all, so they are scanned separately and can never be promoted
# to a runtime consumer. This is not a static interpreter and does not try to
# be: a token reaching a helper function is reported as helper-mediated rather
# than resolved through it.
#
# Known limitation, and it grows as migration proceeds. A consumer that reads
# through a shared resolver contains no path token of its own, so this tool
# cannot see it. Phase 6G.2 produced the first instances: after
# test_network_behaviour_coupling.R and
# render_differential_network_figures.R were repointed onto
# resolve_spatial_network_object() and
# resolve_bootstrap_differential_tables(), both dropped out of the inventory
# while still reading exactly what they read before. The dependency did not
# disappear, it moved to R/networks/spatial_network_utils.R, which the
# inventory does report.
#
# So a shrinking inventory is not evidence that consumers went away, and must
# be reconciled against the pre-migration set rather than read on its own
# (audits/phase6g_spatial_networks_consumer_reconciliation.csv). The
# compensating control is by function name instead of by path: the callers of
# a resolver are enumerable, and the writer-namespace test asserts that every
# caller actually loads the library that defines it.

source(file.path("R", "paths.R"))
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "pipeline_registry.R"))

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = "") {
  i <- which(args == flag)
  if (!length(i) || i[1] == length(args)) return(default)
  args[i[1] + 1L]
}
ANALYSIS_ID <- arg_value("--analysis-id")
QUIET <- "--quiet" %in% args
if (!nzchar(ANALYSIS_ID)) {
  stop("usage: --analysis-id <analysis id, e.g. run_ewce_celltype_enrichment>",
       call. = FALSE)
}

## The sibling repository is scanned only when it is pointed at explicitly.
## Constructing its path by default would make this tool a runtime cross-repo
## dependency, which the boundary audit forbids and which would be a real
## regression rather than a technicality: the two repositories are coupled only
## through the frozen export bundle. Set EXP9_MANUSCRIPT_ROOT to include the
## manuscript surface, as its own importer requires PROTEOMICS_ROOT to be set.
MR <- Sys.getenv("EXP9_MANUSCRIPT_ROOT", unset = "")

say <- function(...) if (!QUIET) cat(...)

# --------------------------------------------------------------- vocabulary

registry <- read_pipeline_registry(repo_path("pipeline.yml"))
steps <- pipeline_steps(registry, pipeline_stage_names(registry),
                        dataset = "all", include_unsupported = TRUE)
steps <- steps[!duplicated(steps$script), , drop = FALSE]

split_paths <- function(x) {
  p <- trimws(unlist(strsplit(paste(x, collapse = "|"), "[|;]")))
  unique(p[nzchar(p)])
}

entry <- steps$script[sub("[.][Rr]$", "", basename(steps$script)) == ANALYSIS_ID]
if (!length(entry)) {
  # an unregistered writer is still enumerable if the file exists
  cand <- list.files(repo_path("analysis"), pattern = paste0("^", ANALYSIS_ID, "[.][Rr]$"),
                     recursive = TRUE, full.names = FALSE)
  if (length(cand)) entry <- file.path("analysis", cand[1])
}
if (!length(entry)) stop("no entrypoint found for analysis id '", ANALYSIS_ID, "'", call. = FALSE)
entry <- entry[1]
DOMAIN <- sub("^analysis/([^/]+)/.*", "\\1", entry)

## declared outputs now, and the historical paths Phase 6F recorded
declared_now <- split_paths(steps$produces[steps$script == entry])
declared_hist <- character(0)
inv6f <- repo_path("audits", "phase6f_output_inventory.csv")
if (file.exists(inv6f)) {
  d6 <- utils::read.csv(inv6f, stringsAsFactors = FALSE)
  ## Match on the basename: Phase 6F recorded these rows before the analysis
  ## directories were renamed, so the full path no longer agrees and matching
  ## on it loses the historical namespace entirely.
  d6 <- d6[basename(d6$producing_analysis) == basename(entry), , drop = FALSE]
  declared_hist <- unique(c(d6$current_path, d6$proposed_path))
  declared_hist <- declared_hist[nzchar(declared_hist)]
}
all_paths <- unique(c(declared_now, declared_hist))
all_paths <- all_paths[nzchar(all_paths)]

ARTIFACT_WORDS <- c("results", "data", "processed", "work", "exports", "tables",
                    "figures", "plots", "models", "manifests", "reports",
                    "source_data", "logs", "reviewer_audit")
SCOPE_WORDS <- c("global", "<dataset>", valid_datasets())

## Which path segments are this analysis's identity, and which are internal
## structure. The distinction is positional, not lexical.
##
##   results/tables/07_spatial_networks/bootstrap_differential_network_stability/01_Tables/x.csv
##                  ^ stage namespace   ^ substep                               ^ structure
##
## The stage namespace and the substep name an analysis; a subdirectory below
## them organises its files. Searching for "01_Tables" attributes any analysis
## that happens to use the same subdivision to this one, which is how
## differential_abundance/compare_go_enrichment.R was wrongly reported as a
## spatial_networks consumer. Note that a lexical test cannot separate these:
## "01_Tables" looks exactly like a stage namespace, and it appears in only one
## declared path, so a distinctiveness test passes it too.
identity_segments <- function(p) {
  q <- strsplit(p, "/", fixed = TRUE)[[1]]
  q <- q[nzchar(q)]
  if (!length(q)) return(character(0))
  ## drop the lifecycle root and the artifact type, then keep the next two
  while (length(q) && q[1] %in% ARTIFACT_WORDS) q <- q[-1]
  utils::head(q, 2)
}
segments <- unique(unlist(lapply(all_paths, identity_segments)))
segments <- segments[nzchar(segments)]
segments <- setdiff(segments, c(ARTIFACT_WORDS, SCOPE_WORDS))
segments <- segments[!grepl("^[*<]", segments)]
segments <- segments[!grepl("[.][A-Za-z]{2,5}$", segments)]
stage_ns <- segments[grepl("^[0-9]{2}[a-z]?_", segments)]
leaf_ns <- setdiff(segments, c(stage_ns, ANALYSIS_ID, DOMAIN))
leaf_ns <- leaf_ns[nchar(leaf_ns) >= 4L]

## The same distinctiveness rule the filenames get, for the same reason. A
## segment such as "01_Tables" is a generic subdirectory that several analyses
## use, so searching for it attributes their code to this analysis. A segment
## is this analysis's identity only if no other analysis's declared paths
## contain it. Measured against the inventory, not asserted.
if (length(leaf_ns) && !is.null(d6_all <- tryCatch(
      utils::read.csv(inv6f, stringsAsFactors = FALSE), error = function(e) NULL))) {
  other_paths <- d6_all$current_path[basename(d6_all$producing_analysis) != basename(entry)]
  other_segments <- unique(unlist(strsplit(other_paths, "/", fixed = TRUE)))
  shared_segments <- intersect(leaf_ns, other_segments)
  leaf_ns <- setdiff(leaf_ns, shared_segments)
} else {
  shared_segments <- character(0)
}

## canonical output filenames, read out of the writer's own write calls
EMPTY_OK <- function(x) tryCatch(is.null(x) || (is.symbol(x) && !nzchar(as.character(x))),
                                 error = function(...) TRUE)
walk <- function(e, fn) {
  if (!tryCatch({ e; TRUE }, error = function(...) FALSE)) return(invisible(NULL))
  fn(e)
  if (is.call(e) && is.name(e[[1]]) && identical(as.character(e[[1]]), "function")) {
    if (length(e) >= 3L) walk(e[[3]], fn)
    return(invisible(NULL))
  }
  if (is.call(e) || is.expression(e) || is.list(e)) {
    for (i in seq_along(e)) {
      el <- tryCatch(e[[i]], error = function(...) NULL)
      if (EMPTY_OK(el)) next
      walk(el, fn)
    }
  }
  invisible(NULL)
}
call_name <- function(e) {
  if (!is.call(e)) return(NA_character_)
  fn <- e[[1]]
  if (is.name(fn)) return(as.character(fn))
  if (is.call(fn) && length(fn) == 3L && as.character(fn[[1]]) %in% c("::", ":::")) {
    return(as.character(fn[[3]]))
  }
  NA_character_
}

DATA_EXT <- "[.](csv|tsv|txt|xlsx|xls|rds|rda|RData|yml|yaml|json|pdf|svg|png|tiff|gct)$"

## Output filenames come only from literals inside a write call. Taking every
## literal with a data extension pulls in this analysis's *inputs* (the imputed
## matrix pattern) and bare suffixes such as ".rds", which then match a third of
## the repository. A filename is an output because something writes it.
WRITE_CALLS_FOR_VOCAB <- c("write.csv", "write.table", "writeLines", "saveRDS",
                           "ggsave", "saveWorkbook", "write.xlsx",
                           "write_csv_safe", "write_csv_safe2", "write_yaml")
output_filenames <- character(0)
if (file.exists(entry)) {
  ex <- tryCatch(parse(entry), error = function(e) NULL)
  if (!is.null(ex)) for (e in ex) walk(e, function(x) {
    nm <- call_name(x)
    if (is.na(nm) || !nm %in% WRITE_CALLS_FOR_VOCAB) return(invisible(NULL))
    walk(x, function(y) {
      if (is.character(y) && length(y) == 1L && grepl(DATA_EXT, y)) {
        output_filenames <<- c(output_filenames, basename(y))
      }
    })
    invisible(NULL)
  })
}
## a stem of at least six characters: ".rds" and "qc.csv" are not identities
output_filenames <- unique(output_filenames[nzchar(output_filenames)])
output_filenames <- output_filenames[nchar(sub(DATA_EXT, "", output_filenames)) >= 6L]

## Distinctiveness, measured against the repository rather than asserted: a
## filename shared with another registered entrypoint is shared infrastructure
## (sessionInfo.txt, run_manifest.yml) and not this analysis's identity.
other_entries <- setdiff(steps$script[grepl("^analysis/", steps$script)], entry)
other_text <- vapply(other_entries, function(p) {
  if (!file.exists(p)) return("")
  paste(readLines(p, warn = FALSE), collapse = "\n")
}, character(1))
if (length(output_filenames) && length(other_text)) {
  shared <- vapply(output_filenames, function(tok)
    sum(vapply(other_text, function(t) grepl(tok, t, fixed = TRUE), logical(1))),
    integer(1))
  generic_filenames <- names(shared)[shared > 0L]
  output_filenames <- setdiff(output_filenames, generic_filenames)
} else {
  generic_filenames <- character(0)
}

vocab <- list(
  FULL_PATH = unique(sub("/$", "", all_paths[grepl("/", all_paths)])),
  NAMESPACE_BASENAME = unique(c(stage_ns, leaf_ns)),
  OUTPUT_FILENAME = output_filenames,
  ANALYSIS_REFERENCE = unique(c(ANALYSIS_ID, basename(entry)))
)
vocab <- lapply(vocab, function(v) v[nzchar(v)])

say("analysis id        :", ANALYSIS_ID, "\n")
say("entrypoint         :", entry, "\n")
say("domain             :", DOMAIN, "\n")
say("declared paths     :", length(all_paths), "\n")
say("stage namespaces   :", paste(stage_ns, collapse = ", "), "\n")
say("leaf namespaces    :", paste(leaf_ns, collapse = ", "), "\n")
say("output filenames   :", length(output_filenames), "\n")

# ------------------------------------------------------------------ surfaces

tracked <- system2("git", c("ls-files"), stdout = TRUE)
SURFACES <- "^(analysis|R|tests|tools|audits|config|docs)/|^pipeline\\.yml$"
files <- tracked[grepl(SURFACES, tracked)]

## This tool's own output must not be part of its input. A consumer inventory
## quotes the matched line for every finding, so once one is committed the next
## run matches every token in it and reports the inventory as a consumer of the
## analysis it describes. That is a feedback loop, not a dependency: committing
## the six spatial_networks inventories turned 0 unclassified rows into 67332
## and blocked the preflight. Excluded by path, so the exclusion cannot be
## defeated by the file merely being untracked at the time.
TOOL_OUTPUT <- "^audits/(consumer_inventory/|phase6g_preflight_)"
files <- files[!grepl(TOOL_OUTPUT, files)]
files <- files[grepl("[.]([Rr]|ya?ml|md|csv|tsv|txt|json)$", files) | files == "pipeline.yml"]
files <- data.frame(path = files, repo = "pRoteomics", stringsAsFactors = FALSE)

## the manuscript repository, where the frozen interface legitimately names
## scientific outputs
if (nzchar(MR) && dir.exists(MR)) {
  mt <- suppressWarnings(system2("git", c("-C", MR, "ls-files"), stdout = TRUE))
  mt <- mt[grepl("[.]([Rr]|ya?ml|md|csv|tsv)$", mt)]
  if (length(mt)) {
    files <- rbind(files, data.frame(path = mt, repo = "Exp9_manuscript",
                                     stringsAsFactors = FALSE))
  }
}
files$abs <- ifelse(files$repo == "pRoteomics", repo_path(files$path),
                    file.path(MR, files$path))
say("surfaces scanned   :", nrow(files), "files across",
    length(unique(files$repo)), "repositories\n\n")

# --------------------------------------------------------------- classifiers

HISTORICAL_SURFACE <- paste0(
  "^(archive/|audits/(phase6|restructure_|test_migration|part29|program_evidence|",
  "publication_hardening|migration)|docs/(NAMING_MIGRATION|PRERESTRUCTURE_FREEZE|",
  "RESTRUCTURE_PLAN|restructure_inventory|publication_freeze_manifest|",
  "MANUSCRIPT_DRAFTING_PROGRESS|PUBLICATION_HARDENING_PROGRESS)|",
  "provenance/|source_data/)")

## Which record files are regenerated by code, and which are maintained by
## hand. A generated record repoints itself the next time its generator runs;
## a hand-maintained contract table does not, so it is a real consumer. The
## distinction is measured by looking for a writer, not assumed.
generated_artifacts <- local({
  cand <- tracked_all <- system2("git", c("ls-files"), stdout = TRUE)
  cand <- cand[grepl("^(config|docs|audits)/.*[.](csv|tsv|ya?ml|md)$", cand)]
  code <- tracked_all[grepl("^(tools|audits|analysis|R|tests)/.*[.][Rr]$", tracked_all)]
  code_txt <- vapply(code, function(p) {
    if (!file.exists(p) || file.info(p)$size > 4e6) return("")
    paste(readLines(p, warn = FALSE), collapse = "\n")
  }, character(1))
  writes <- grepl("(write[._][A-Za-z]+|writeLines|saveWorkbook|write_yaml)", code_txt)
  code_txt <- code_txt[writes]
  if (!length(code_txt)) return(character(0))
  keep <- vapply(cand, function(f) {
    b <- basename(f)
    any(vapply(code_txt, function(t) grepl(b, t, fixed = TRUE), logical(1)))
  }, logical(1))
  cand[keep]
})

## Scripts that address outputs by their *baseline* path.
##
## audits/verify_scientific_contracts.R names frozen carriers such as
## results/tables/11_spatial_systems/atlas/... and then resolves each one
## through audits/restructure_migration_map.csv (baseline_path ->
## destination_path). It is an equivalence oracle keyed on the pre-restructure
## address, not a reader of the current namespace, so migrating a writer must
## not repoint it: the map already carries the translation, and rewriting the
## keys would destroy the very provenance it exists to prove.
##
## Before this rule those 212 hits landed in UNKNOWN, because the classifier's
## fallthrough knew tools/, R/, config/, docs/ and analysis/ but not an
## executable script under audits/. That blocked the Phase 6G.3 preflight.
##
## Derived rather than listed, in keeping with the rest of this tool: the
## property that matters is "resolves through the migration map", which is
## observable. An audits/ script that does NOT do so is classified
## AUDIT_READER and stays runtime-relevant, so a genuine new reader cannot be
## excused by living in the same directory.
baseline_keyed_surfaces <- local({
  code <- system2("git", c("ls-files"), stdout = TRUE)
  code <- code[grepl("^(audits|tools)/.*[.][Rr]$", code)]
  keep <- vapply(code, function(p) {
    if (!file.exists(p) || file.info(p)$size > 4e6) return(FALSE)
    txt <- paste(readLines(p, warn = FALSE), collapse = "\n")
    grepl("restructure_migration_map", txt, fixed = TRUE)
  }, logical(1))
  unname(code[keep])
})

READ_CALLS <- c("read.csv", "read.delim", "read.table", "readRDS", "readLines",
                "read_csv", "read_tsv", "read_yaml", "read.xlsx", "read_excel",
                "fread", "load", "list.files", "dir", "Sys.glob", "file.exists",
                "dir.exists", "file.info", "readxl::read_excel", "loadWorkbook",
                "read_csv_required", "read_csv_optional", "latest_matching_file",
                "resolve_input_path")
WRITE_CALLS <- c("write.csv", "write.table", "writeLines", "saveRDS", "ggsave",
                 "saveWorkbook", "write.xlsx", "write_csv_safe", "write_csv_safe2",
                 "dir_create", "file.copy", "file.rename", "write_yaml",
                 "write_run_manifest", "write_result_manifest")
PATH_BUILDERS <- c("file.path", "repo_path", "path_results", "path_processed",
                   "canonical_result_path", "canonical_work_path", "path_export",
                   "path_work", "here", "normalizePath")
EXPECT_CALLS <- c("expect_true", "expect_false", "expect_match", "expect_equal",
                  "expect_identical", "expect_gt", "expect_lt", "expect_setequal",
                  "grepl", "sub", "gsub", "regexpr", "startsWith", "endsWith")

# A full path and an output filename are distinctive enough to match as plain
# substrings. A NAMESPACE_BASENAME is not: it is a bare directory name such as
# "atlas", "precision" or "data_contract", and a substring match on those
# produces dependencies that do not exist.
#
# Phase 6G.4 found two concrete false positives, and they had already misled a
# phase brief into asserting three cross-domain dependencies that were never
# there: "data_contract" matched inside R/data_contracts/dataset_config.R, an
# unrelated library path, and "precision" matched the word precision inside a
# sentence about adjustment families.
#
# So a directory name only counts when it is addressed as a path segment:
# followed by a separator, or standing alone as a complete quoted element.
# \Q...\E quotes the token for the regex engine so no token needs escaping by
# hand.
match_token <- function(txt, token, kind = NA_character_) {
  if (!identical(kind, "NAMESPACE_BASENAME")) {
    return(grepl(token, txt, fixed = TRUE))
  }
  as_segment <- grepl(paste0("(^|[^A-Za-z0-9_])\\Q", token, "\\E/"), txt, perl = TRUE)
  as_element <- grepl(paste0("[\"']\\Q", token, "\\E[\"']"), txt, perl = TRUE)
  as_segment | as_element
}

classify_match_type <- function(token, kind_hint, file, in_comment, is_glob) {
  if (in_comment) return("DOCUMENTATION_ONLY")
  if (is_glob) return("CONFIG_GLOB")
  if (kind_hint == "FULL_PATH") return(if (startsWith(file, "tests/")) "TEST_FRAGMENT" else "FULL_PATH")
  if (kind_hint == "OUTPUT_FILENAME") return("OUTPUT_FILENAME")
  if (kind_hint == "NAMESPACE_BASENAME") {
    return(if (startsWith(file, "tests/")) "TEST_FRAGMENT" else "NAMESPACE_BASENAME")
  }
  kind_hint
}

rows <- list()
add_row <- function(...) rows[[length(rows) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)

for (i in seq_len(nrow(files))) {
  f <- files$path[i]
  abs <- files$abs[i]
  repo <- files$repo[i]
  if (!file.exists(abs)) next
  if (file.info(abs)$size > 6e6) next
  if (identical(f, "audits/consumer_inventory")) next
  lines <- tryCatch(readLines(abs, warn = FALSE), error = function(e) character(0))
  if (!length(lines)) next
  if (identical(basename(f), paste0(ANALYSIS_ID, ".R"))) next   # the writer itself

  is_r <- grepl("[.][Rr]$", f)
  is_cfg <- grepl("[.]ya?ml$", f)
  is_doc <- grepl("[.]md$", f)
  is_hist <- grepl(HISTORICAL_SURFACE, f)
  is_test <- startsWith(f, "tests/")

  comment_line <- if (is_r) grepl("^\\s*#", lines) else rep(FALSE, length(lines))

  ## ---- structural pass over R code -----------------------------------
  ast_hits <- list()
  if (is_r) {
    ex <- tryCatch(parse(abs), error = function(e) NULL)
    if (!is.null(ex)) {
      for (e in ex) {
        walk(e, function(x) {
          nm <- call_name(x)
          if (is.na(nm)) return(invisible(NULL))
          ctx <- if (nm %in% READ_CALLS) "READ"
                 else if (nm %in% WRITE_CALLS) "WRITE"
                 else if (nm %in% PATH_BUILDERS) "BUILD"
                 else if (nm %in% EXPECT_CALLS) "ASSERT"
                 else NA_character_
          if (is.na(ctx)) return(invisible(NULL))
          lits <- character(0)
          walk(x, function(y) if (is.character(y) && length(y) == 1L) lits <<- c(lits, y))
          if (!length(lits)) return(invisible(NULL))
          joined <- paste(lits, collapse = "/")
          for (kind in names(vocab)) {
            for (tok in vocab[[kind]]) {
              if (match_token(joined, tok, kind) ||
                  any(vapply(lits, function(z) match_token(z, tok, kind), logical(1)))) {
                ast_hits[[length(ast_hits) + 1L]] <<- list(
                  ctx = ctx, call = nm, kind = kind, token = tok,
                  expr = paste(utils::head(lits, 4), collapse = " , "))
              }
            }
          }
          invisible(NULL)
        })
      }
    }
  }

  ## ---- textual pass, for every surface --------------------------------
  for (kind in names(vocab)) {
    for (tok in vocab[[kind]]) {
      hit_lines <- which(match_token(lines, tok, kind))
      if (!length(hit_lines)) next
      for (ln in hit_lines) {
        in_comment <- isTRUE(comment_line[ln])
        is_glob <- is_cfg && grepl("[*]", lines[ln])
        ast <- Filter(function(h) h$kind == kind && h$token == tok, ast_hits)
        ctx <- if (length(ast)) ast[[1]]$ctx else NA_character_
        callnm <- if (length(ast)) ast[[1]]$call else NA_character_

        mt <- classify_match_type(tok, kind, f, in_comment, is_glob)
        if (identical(kind, "ANALYSIS_REFERENCE") && !in_comment) mt <- "ANALYSIS_REFERENCE"

        ## dependency kind
        dk <- if (in_comment || is_doc) {
          "DOCUMENTATION"
        } else if (is_hist || f %in% baseline_keyed_surfaces) {
          "PROVENANCE"
        } else if (identical(f, "pipeline.yml")) {
          "CONFIG"
        } else if (is_cfg) {
          "CONFIG"
        } else if (is_test) {
          "TEST"
        } else if (!is.na(ctx) && ctx == "READ") {
          "DIRECT_READER"
        } else if (!is.na(ctx) && ctx == "WRITE") {
          "WRITER_NOT_CONSUMER"
        } else if (!is.na(ctx) && ctx == "BUILD") {
          "INDIRECT_READER_VIA_HELPER"
        } else if (startsWith(f, "tools/")) {
          "EXPORTER"
        } else if (startsWith(f, "R/")) {
          "INDIRECT_READER_VIA_HELPER"
        } else if (f %in% generated_artifacts) {
          ## regenerated by a tool: it will repoint itself
          "GENERATED_RECORD"
        } else if (grepl("^(config|docs)/.*[.](csv|tsv)$", f)) {
          ## a contract table nothing regenerates, so it must be updated by hand
          "CONTRACT_RECORD"
        } else if (startsWith(f, "analysis/")) {
          "DOWNSTREAM_ANALYSIS"
        } else if (startsWith(f, "audits/") && is_r) {
          ## an executable audit that is not keyed on the migration map really
          ## does read the current namespace, so it must be repointed
          "AUDIT_READER"
        } else {
          "UNKNOWN"
        }

        ## a helper function definition that synthesises the path
        if (is_r && !in_comment && grepl("<- function", lines[ln])) mt <- "HELPER_FUNCTION"
        ## registry dependency, when the hit is in a consumes field
        if (identical(f, "pipeline.yml") && grepl("consumes", lines[ln])) {
          mt <- "REGISTRY_DEPENDENCY"
        }

        ## A mention of the script itself is not a dependency on its output
        ## namespace, so it is recorded but not counted as runtime-relevant.
        runtime <- !(dk %in% c("DOCUMENTATION", "PROVENANCE", "WRITER_NOT_CONSUMER", "GENERATED_RECORD")) &&
          !identical(kind, "ANALYSIS_REFERENCE")

        conf <- if (!is.na(ctx) && ctx %in% c("READ", "BUILD")) "HIGH"
                else if (mt %in% c("CONFIG_GLOB", "REGISTRY_DEPENDENCY", "FULL_PATH")) "HIGH"
                else if (mt %in% c("OUTPUT_FILENAME", "HELPER_FUNCTION")) "MEDIUM"
                else if (mt %in% c("DOCUMENTATION_ONLY")) "N/A"
                else "MEDIUM"

        add_row(
          analysis_id = ANALYSIS_ID,
          repo = repo,
          consumer_file = f,
          line_or_expression = paste0(ln, ": ", trimws(substr(lines[ln], 1, 160))),
          match_type = mt,
          matched_token = tok,
          enclosing_call = if (is.na(callnm)) "" else callnm,
          dependency_kind = dk,
          active_or_historical = if (is_hist) "historical" else "active",
          runtime_relevant = runtime,
          confidence = conf,
          adjudication = "")
      }
    }
  }
}

inv <- if (length(rows)) do.call(rbind, rows) else
  data.frame(analysis_id = character(0), repo = character(0),
             consumer_file = character(0), line_or_expression = character(0),
             match_type = character(0), matched_token = character(0),
             enclosing_call = character(0), dependency_kind = character(0),
             active_or_historical = character(0), runtime_relevant = logical(0),
             confidence = character(0), adjudication = character(0),
             stringsAsFactors = FALSE)

## one row per (file, line, token): the same token found by both passes is one
## finding, not two
if (nrow(inv)) {
  key <- paste(inv$consumer_file, inv$line_or_expression, inv$matched_token)
  inv <- inv[!duplicated(key), , drop = FALSE]
  inv <- inv[order(!inv$runtime_relevant, inv$dependency_kind, inv$consumer_file), ]
}

out_dir <- repo_path("audits", "consumer_inventory")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, paste0(ANALYSIS_ID, ".csv"))
utils::write.csv(inv, out_file, row.names = FALSE)

say("findings           :", nrow(inv), "\n")
say("distinct files     :", length(unique(inv$consumer_file)), "\n\n")
if (nrow(inv)) {
  say("=== match_type ===\n"); if (!QUIET) print(table(inv$match_type))
  say("\n=== dependency_kind ===\n"); if (!QUIET) print(table(inv$dependency_kind))
  say("\nruntime-relevant findings :", sum(inv$runtime_relevant), "\n")
  say("runtime-relevant files    :",
      length(unique(inv$consumer_file[inv$runtime_relevant])), "\n")
  say("UNKNOWN dependency kinds  :", sum(inv$dependency_kind == "UNKNOWN"), "\n")
  rt <- unique(inv$consumer_file[inv$runtime_relevant])
  if (length(rt)) {
    say("\n=== runtime-relevant consumer files ===\n")
    for (x in rt) {
      k <- unique(inv$dependency_kind[inv$runtime_relevant & inv$consumer_file == x])
      say(sprintf("  %-72s %s\n", x, paste(k, collapse = ",")))
    }
  }
}
say("\nwrote ", relative_to(out_file), "\n")

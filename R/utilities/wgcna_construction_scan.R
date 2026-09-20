# Generalized enumerator for historical WGCNA path constructions.
#
# Phase 6G.8 Batch 3A froze a 176-row construction oracle built by scanning for
# the literal call path_results(...). That method was sound but incomplete: the
# same historical path can be built as file.path("results", ...), through a
# local alias variable, through a one-line alias function, or behind Sys.glob().
# A migration that only converts the path_results spelling leaves live readers
# pointing at the historical tree while the family is reported complete.
#
# This file answers one question for any R source file: which expressions in it
# construct a path into the historical WGCNA output tree, and what object does
# each one name? It classifies; it never edits.
#
# Substring search is deliberately NOT the classifier. "06_modules_WGCNA" also
# appears in comments, in frozen manifests and in test fixtures that assert the
# historical contract on purpose. Detection therefore runs on parse structure:
# a construction is a CALL whose resolved segment list contains the stage, where
# segments come from literal arguments plus alias expansion.

if (!exists("repo_path", mode = "function")) {
  paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
  source(paths_file)
}

WCS_STAGE <- "06_modules_WGCNA"

# Historical result families. reviewer_audit carries no stage segment at all,
# so it is matched on its own rather than derived from the stage position.
WCS_RESULT_FAMILIES <- c("tables", "figures", "source_data", "logs", "reports")
WCS_REVIEWER_AUDIT <- "reviewer_audit"

# Path builders and the segments each one implies before its own arguments.
WCS_BUILDER_PREFIX <- list(
  path_results    = c("results"),
  path_processed  = c("data", "processed"),
  path_figures    = c("results", "figures"),
  path_tables     = c("results", "tables"),
  repo_path       = character(0),
  file.path       = character(0),
  Sys.glob        = character(0),
  normalizePath   = character(0),
  test_path       = character(0)
)

# Calls that consume a path and reveal wildcard discovery.
WCS_GLOB_CALLS <- c("Sys.glob", "list.files", "list.dirs", "dir")

# ---------------------------------------------------------------- span text
#
# A parse node's span is (line1, col1) to (line2, col2). The interior-line
# range (line1+1):(line2-1) is only valid when an interior line exists; for a
# 2-line span it DESCENDS and silently splices the wrong lines. That defect
# produced a whole family of false "site not found" results in Batch 3A, so the
# guard lives here with the extractor rather than in each caller.
wcs_span_text <- function(lines, line1, col1, line2, col2) {
  if (is.na(line1) || is.na(line2) || line2 < line1) {
    stop("malformed span: line1=", line1, " line2=", line2, call. = FALSE)
  }
  if (line1 < 1L || line2 > length(lines)) {
    stop("span outside file: ", line1, "-", line2, " of ", length(lines), call. = FALSE)
  }
  if (line1 == line2) return(substr(lines[line1], col1, col2))
  mid <- if (line2 > line1 + 1L) lines[(line1 + 1L):(line2 - 1L)] else character(0)
  paste(c(substr(lines[line1], col1, nchar(lines[line1])),
          mid,
          substr(lines[line2], 1L, col2)), collapse = "\n")
}

# ------------------------------------------------------------ segment logic
#
# A literal character argument contributes its own value; anything else
# contributes NA, which marks a runtime-determined segment. NA segments are
# kept rather than dropped: their POSITION still matters when locating the
# stage and the family.
wcs_arg_segments <- function(a, aliases) {
  ## A literal may carry several segments at once:
  ##   file.path(repo_root, "results/tables/06_modules_WGCNA/group_effects/ds/x.csv")
  ## builds the same path as the segment-per-argument spelling, but no single
  ## argument equals the stage. Splitting on "/" makes the two forms compare
  ## equal. Found by the independent cross-check, which saw ten constructions in
  ## analysis/wgcna/audit_microglia_module_claims.R that this scanner did not.
  if (is.character(a) && length(a) == 1L) {
    if (!grepl("/", a, fixed = TRUE)) return(a)
    s <- strsplit(a, "/", fixed = TRUE)[[1]]
    return(s[nzchar(s)])
  }
  if (is.name(a)) {
    nm <- as.character(a)
    if (!is.null(aliases[[nm]])) return(aliases[[nm]])
    return(NA_character_)
  }
  if (is.call(a)) {
    inner <- wcs_call_segments(a, aliases)
    if (!is.null(inner)) return(inner)
    ## repo_root() and getwd() resolve to the repository root, contributing no
    ## named segment of their own.
    fn <- if (is.name(a[[1]])) as.character(a[[1]]) else ""
    if (fn %in% c("repo_root", "getwd")) return(character(0))
    return(NA_character_)
  }
  NA_character_
}

# Resolve a call to the segment list it builds, or NULL if it is not a path
# construction at all.
# A namespaced call, pkg::fn(...), has a CALL in the function position rather
# than a name, so the bare as.character() test misses it. test fixtures reach
# the historical tree through testthat::test_path("..", "..", "results", ...).
wcs_call_name <- function(x) {
  if (is.name(x)) return(as.character(x))
  if (is.call(x) && is.name(x[[1]]) &&
      as.character(x[[1]]) %in% c("::", ":::") && length(x) == 3L) {
    return(as.character(x[[3]]))
  }
  NA_character_
}

wcs_call_segments <- function(cl, aliases) {
  if (!is.call(cl)) return(NULL)
  fn <- wcs_call_name(cl[[1]])
  if (is.na(fn)) return(NULL)
  prefix <- WCS_BUILDER_PREFIX[[fn]]
  if (is.null(prefix)) {
    ## a local alias FUNCTION, e.g. base <- function(...) path_results("tables", ...)
    af <- aliases[[paste0("fn:", fn)]]
    if (is.null(af)) return(NULL)
    prefix <- af
  }
  args <- as.list(cl)[-1]
  args <- args[!nzchar(names(args) %||% rep("", length(args))) |
                 is.null(names(args))]
  segs <- prefix
  for (a in args) segs <- c(segs, wcs_arg_segments(a, aliases))
  segs
}

# `%||%` comes from the canonical R/null_coalescing.R, loaded via R/paths.R in
# the bootstrap above. It is deliberately NOT redefined here: the repository
# keeps exactly one definition so that its NA-aware semantics cannot fork.

# Does a resolved segment list land inside the historical WGCNA tree?
wcs_is_historical <- function(segs) {
  if (!length(segs)) return(FALSE)
  s <- segs[!is.na(segs)]
  if (WCS_STAGE %in% s) return(TRUE)
  ## reviewer_audit has no stage segment; only the WGCNA-owned subset counts,
  ## and that is decided by the caller from the family segment.
  FALSE
}

# Pull the object identity out of a resolved segment list.
wcs_describe <- function(segs) {
  i <- which(segs == WCS_STAGE)
  if (!length(i)) return(NULL)
  i <- i[[1]]
  before <- segs[seq_len(i - 1L)]
  after <- if (i < length(segs)) segs[(i + 1L):length(segs)] else character(0)
  kind <- {
    b <- before[!is.na(before)]
    hit <- b[b %in% c(WCS_RESULT_FAMILIES, WCS_REVIEWER_AUDIT)]
    if (length(hit)) hit[[length(hit)]] else if ("processed" %in% b) "processed" else NA_character_
  }
  family <- if (length(after) && !is.na(after[[1]])) after[[1]] else NA_character_
  artifact <- {
    lit <- after[!is.na(after)]
    hit <- lit[grepl("[.][A-Za-z0-9]{2,5}$", lit)]
    if (length(hit)) hit[[length(hit)]] else NA_character_
  }
  list(kind = kind, family = family, artifact = artifact,
       n_after = length(after),
       has_wildcard = any(grepl("[*?]", segs[!is.na(segs)])),
       dynamic_segments = sum(is.na(segs)))
}

# ------------------------------------------------------------------- aliases
#
# A variable alias is any assignment whose right-hand side is a path
# construction that stops at a directory; a function alias is a one-line
# wrapper forwarding `...` to a path builder, e.g.
#   base <- function(...) path_results("tables", ...)
# Both hide the stage segment from any search that looks for it literally.
#
# Exported rather than inlined in the scanner because the converter must resolve
# call arguments against exactly the same alias table. Two implementations would
# drift, and a drifting alias table silently changes which object a rewrite
# names.
wcs_file_aliases <- function(file) {
  p <- tryCatch(parse(file, keep.source = TRUE), error = function(e) NULL)
  if (is.null(p)) return(list())
  pd <- utils::getParseData(p)
  if (is.null(pd) || !nrow(pd)) return(list())
  lines <- readLines(file, warn = FALSE)
  aliases <- list()
  assigns <- pd[pd$token %in% c("LEFT_ASSIGN", "RIGHT_ASSIGN", "EQ_ASSIGN"), ]
  for (k in seq_len(nrow(assigns))) {
    par <- pd[pd$id == assigns$parent[k], ]
    if (!nrow(par)) next
    txt <- tryCatch(wcs_span_text(lines, par$line1, par$col1, par$line2, par$col2),
                    error = function(e) NA_character_)
    if (is.na(txt)) next
    ex <- tryCatch(str2lang(txt), error = function(e) NULL)
    if (is.null(ex) || !is.call(ex) || length(ex) < 3L) next
    lhs <- ex[[2]]; rhs <- ex[[3]]
    if (!is.name(lhs)) next
    nm <- as.character(lhs)
    if (is.call(rhs) && is.name(rhs[[1]]) && as.character(rhs[[1]]) == "function") {
      body_ex <- rhs[[3]]
      if (is.call(body_ex) && is.name(body_ex[[1]])) {
        bfn <- as.character(body_ex[[1]])
        pre <- WCS_BUILDER_PREFIX[[bfn]]
        if (!is.null(pre)) {
          ba <- as.list(body_ex)[-1]
          lits <- character(0)
          for (a in ba) {
            if (is.character(a) && length(a) == 1L) lits <- c(lits, a) else break
          }
          if (any(vapply(ba, function(a) is.name(a) && as.character(a) == "...", logical(1)))) {
            aliases[[paste0("fn:", nm)]] <- c(pre, lits)
          }
        }
      }
      next
    }
    ## A plain string constant is an alias too. build_wgcna_modules.R holds
    ##   wgcna_module <- "06_modules_WGCNA"
    ## and then builds path_processed(wgcna_module, "01_WGCNA", ...). No
    ## argument of that call is the stage literal, so a call-only alias pass
    ## cannot see it and four live reads of the Stage-01 inputs stayed
    ## invisible. Found by the independent cross-check.
    if (is.character(rhs) && length(rhs) == 1L && nzchar(rhs)) {
      s <- strsplit(rhs, "/", fixed = TRUE)[[1]]
      s <- s[nzchar(s)]
      if (length(s) && !grepl("[.][A-Za-z0-9]{2,5}$", s[[length(s)]])) aliases[[nm]] <- s
      next
    }
    segs <- tryCatch(wcs_call_segments(rhs, aliases), error = function(e) NULL)
    if (is.null(segs) || !length(segs)) next
    lit <- segs[!is.na(segs)]
    if (!length(lit)) next
    if (grepl("[.][A-Za-z0-9]{2,5}$", lit[[length(lit)]])) next   # a file, not a root
    aliases[[nm]] <- segs
  }
  aliases
}

# ------------------------------------------------------------------- scanner
wcs_scan_file <- function(file) {
  p <- tryCatch(parse(file, keep.source = TRUE), error = function(e) NULL)
  if (is.null(p)) return(NULL)
  pd <- utils::getParseData(p)
  if (is.null(pd) || !nrow(pd)) return(NULL)
  lines <- readLines(file, warn = FALSE)

  calls <- pd[pd$token == "SYMBOL_FUNCTION_CALL", ]
  if (!nrow(calls)) return(NULL)
  node_id <- pd$parent[match(calls$parent, pd$id)]
  nodes <- pd[match(node_id, pd$id), ]
  nodes$fn <- calls$text
  nodes <- nodes[!is.na(nodes$line1), ]
  nodes <- nodes[order(nodes$line1, nodes$col1), ]

  lang_of <- function(i) {
    txt <- tryCatch(wcs_span_text(lines, nodes$line1[i], nodes$col1[i],
                                  nodes$line2[i], nodes$col2[i]),
                    error = function(e) NA_character_)
    if (is.na(txt)) return(NULL)
    tryCatch(str2lang(txt), error = function(e) NULL)
  }

  aliases <- wcs_file_aliases(file)

  ## pass 2: constructions
  out <- list()
  for (i in seq_len(nrow(nodes))) {
    cl <- lang_of(i)
    if (is.null(cl)) next
    segs <- tryCatch(wcs_call_segments(cl, aliases), error = function(e) NULL)
    if (is.null(segs) || !wcs_is_historical(segs)) next
    desc <- wcs_describe(segs)
    if (is.null(desc)) next
    ## skip a node whose parent call is itself a detected construction with the
    ## same span start, to avoid double counting builder nesting
    txt <- wcs_span_text(lines, nodes$line1[i], nodes$col1[i], nodes$line2[i], nodes$col2[i])
    out[[length(out) + 1L]] <- data.frame(
      file = file,
      line1 = nodes$line1[i], col1 = nodes$col1[i],
      line2 = nodes$line2[i], col2 = nodes$col2[i],
      fn = nodes$fn[i],
      form = if (nodes$fn[i] == "file.path" && any(names(aliases) %in%
                   vapply(as.list(cl)[-1], function(a) if (is.name(a)) as.character(a) else "", "")))
               "alias" else nodes$fn[i],
      kind = desc$kind, artifact_family = desc$family, artifact = desc$artifact,
      has_wildcard = desc$has_wildcard, dynamic_segments = desc$dynamic_segments,
      in_provenance_call = wcs_in_provenance_call(pd, nodes$id[i]),
      n_lines = nodes$line2[i] - nodes$line1[i] + 1L,
      text = gsub("[[:space:]]+", " ", txt),
      stringsAsFactors = FALSE)
  }
  if (!length(out)) return(NULL)
  res <- do.call(rbind, out)

  ## Drop a construction fully contained in another detected construction from
  ## the same file: file.path(Sys.glob(...)) style nesting would otherwise be
  ## counted twice for one object.
  keep <- rep(TRUE, nrow(res))
  for (a in seq_len(nrow(res))) {
    for (b in seq_len(nrow(res))) {
      if (a == b || !keep[a]) next
      inside <- (res$line1[b] < res$line1[a] ||
                   (res$line1[b] == res$line1[a] && res$col1[b] <= res$col1[a])) &&
        (res$line2[b] > res$line2[a] ||
           (res$line2[b] == res$line2[a] && res$col2[b] >= res$col2[a]))
      strictly <- inside && !(res$line1[a] == res$line1[b] && res$col1[a] == res$col1[b] &&
                                res$line2[a] == res$line2[b] && res$col2[a] == res$col2[b])
      if (strictly) keep[a] <- FALSE
    }
  }
  res[keep, , drop = FALSE]
}

wcs_scan_paths <- function(paths) {
  files <- unlist(lapply(paths, function(p) {
    if (dir.exists(p)) list.files(p, pattern = "[.][Rr]$", recursive = TRUE, full.names = TRUE)
    else if (file.exists(p)) p else character(0)
  }), use.names = FALSE)
  files <- unique(sub("^[.]/", "", files))
  res <- lapply(files, function(f) tryCatch(wcs_scan_file(f), error = function(e) NULL))
  res <- res[!vapply(res, is.null, logical(1))]
  if (!length(res)) return(NULL)
  out <- do.call(rbind, res)
  out$file <- sub("^.*/proteomics/", "", gsub("\\\\", "/", out$file))
  out[order(out$artifact_family, out$file, out$line1), , drop = FALSE]
}

# --------------------------------------------------------------- classifier
#
# Coverage analysis is not a migration list. A test fixture that asserts the
# historical contract, a freeze oracle that hashes historical identity and a
# manifest that records where an input came from are all correct as they are.
WCS_CLASSES <- c("LIVE_READER", "TEST_FIXTURE", "PROVENANCE_ORACLE",
                 "FREEZE_ORACLE", "PROVENANCE_RECORD", "TOOLING_REFERENCE",
                 "LEGACY_NO_ACTIVE_PRODUCER", "RESOLVER_DEFINITION")

# The resolver itself. wg_legacy_root() exists precisely to build the historical
# root, because that is the fallback every migrated reader depends on. It is the
# migration's mechanism, not a consumer of it, and normalizing it would remove
# the fallback entirely. Invisible to the scanner until string-constant aliases
# were bound, which is why it needs saying out loud.
WCS_RESOLVER_FILES <- c("R/wgcna/wgcna_paths.R")

# Files whose historical paths are CITATIONS, not reads. The manuscript
# provenance tables record, for each stated quantity, which artifact is its
# authoritative source; the path is data in a row, never opened. Rewriting one
# would change what the manuscript claims its evidence was.
WCS_PROVENANCE_FILES <- c(
  "audits/publication_hardening/10_manuscript_provenance.R",
  "audits/publication_hardening/11_manuscript_phase2_contracts.R"
)

# Individual functions whose constructed path is a stored identifier rather than
# a location to open. Recorded at function granularity because the enclosing
# file does perform real reads elsewhere.
#
# wgcna_inferential_handoff_source_artifact builds the repo-RELATIVE
# source_artifact column of WGCNA_inferential_handoff.csv. Existing handoff
# artifacts already carry that exact string and the validator compares it for
# equality, so resolving it to an absolute normalized-first path makes every
# stored row fail validation. It was converted once and reverted.
WCS_PROVENANCE_FUNCTIONS <- c(
  "wgcna_inferential_handoff_source_artifact"
)

# Calls whose arguments are provenance, not reads. write_run_manifest() records
# what a run consumed and produced; a path appearing in its inputs/outputs list
# is a statement ABOUT the run, written into the manifest as text and never
# opened. Repointing one would change what the manifest claims the run read.
WCS_PROVENANCE_CALLS <- c("write_run_manifest")

# Is the construction at `line` lexically inside a provenance-recording call?
# Matching is on the enclosing call in the parse tree, not on proximity.
wcs_in_provenance_call <- function(pd, node_id) {
  seen <- 0L
  id <- node_id
  while (length(id) == 1L && !is.na(id) && id != 0L && seen < 200L) {
    seen <- seen + 1L
    row <- pd[pd$id == id, , drop = FALSE]
    if (!nrow(row)) return(FALSE)
    ## A call expr's direct children are expr wrappers; the
    ## SYMBOL_FUNCTION_CALL naming the function sits one level below that, so
    ## the function name has to be read from the grandchildren.
    kids <- pd$id[pd$parent == id]
    fn <- pd$text[pd$token == "SYMBOL_FUNCTION_CALL" &
                    (pd$parent %in% kids | pd$parent == id)]
    if (length(fn) && any(fn %in% WCS_PROVENANCE_CALLS)) return(TRUE)
    id <- row$parent[[1]]
  }
  FALSE
}

# Is this construction inside the body of a known provenance-value function?
wcs_in_provenance_function <- function(file, line) {
  if (!file.exists(file)) return(rep(FALSE, length(line)))
  src <- readLines(file, warn = FALSE)
  starts <- integer(0)
  for (fn in WCS_PROVENANCE_FUNCTIONS) {
    hit <- grep(paste0("^", fn, "[[:space:]]*<-[[:space:]]*function"), src)
    starts <- c(starts, hit)
  }
  if (!length(starts)) return(rep(FALSE, length(line)))
  ## a top-level function body ends at the first column-1 closing brace
  ends <- vapply(starts, function(s) {
    close <- grep("^\\}", src)
    close <- close[close > s]
    if (length(close)) close[[1]] else length(src)
  }, integer(1))
  vapply(line, function(l) any(l >= starts & l <= ends), logical(1))
}

# vapply-safe elementwise isTRUE
isTRUE_vec <- function(x) !is.na(x) & as.logical(x)

wcs_classify <- function(df) {
  cls <- rep(NA_character_, nrow(df))
  f <- df$file
  cls[is.na(cls) & startsWith(f, "tests/")] <- "TEST_FIXTURE"
  cls[is.na(cls) & f %in% WCS_RESOLVER_FILES] <- "RESOLVER_DEFINITION"
  cls[is.na(cls) & f == "R/data_contracts/publication_freeze_utils.R"] <- "FREEZE_ORACLE"
  cls[is.na(cls) & f %in% WCS_PROVENANCE_FILES] <- "PROVENANCE_RECORD"
  if (!is.null(df$in_provenance_call)) {
    cls[is.na(cls) & isTRUE_vec(df$in_provenance_call)] <- "PROVENANCE_RECORD"
  }
  if (!is.null(df$line1)) {
    inprov <- vapply(seq_len(nrow(df)), function(i) {
      if (!is.na(cls[i])) return(FALSE)
      isTRUE(wcs_in_provenance_function(f[i], df$line1[i]))
    }, logical(1))
    cls[inprov] <- "PROVENANCE_RECORD"
  }
  cls[is.na(cls) & startsWith(f, "tools/")] <- "TOOLING_REFERENCE"
  cls[is.na(cls)] <- "LIVE_READER"
  cls
}

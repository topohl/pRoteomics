# Shared path helpers for pRoteomics.
#
# All committed scripts should default to repository-relative paths. For local
# machines, set PROTEOMICS_PROJECT_ROOT or provide a config override outside Git.

.null_coalescing_source_files <- unlist(lapply(
  sys.frames(),
  function(frame) if (!is.null(frame$ofile)) frame$ofile else character()
), use.names = FALSE)
.null_coalescing_search_dirs <- c(
  dirname(.null_coalescing_source_files[file.exists(
    .null_coalescing_source_files
  )]),
  Sys.getenv("PROTEOMICS_PROJECT_ROOT", unset = ""),
  getwd()
)
.null_coalescing_search_dirs <- unique(
  .null_coalescing_search_dirs[nzchar(.null_coalescing_search_dirs)]
)
.null_coalescing_file <- ""
for (.null_coalescing_start in .null_coalescing_search_dirs) {
  .null_coalescing_search_dir <- normalizePath(
    .null_coalescing_start, winslash = "/", mustWork = FALSE
  )
  repeat {
    .null_coalescing_candidate <- Filter(file.exists, c(
      file.path(
        .null_coalescing_search_dir, "R", "utilities", "null_coalescing.R"
      ),
      file.path(.null_coalescing_search_dir, "R", "null_coalescing.R")
    ))[1]
    if (!is.na(.null_coalescing_candidate)) {
      .null_coalescing_file <- .null_coalescing_candidate
      break
    }
    .null_coalescing_parent <- dirname(.null_coalescing_search_dir)
    if (identical(.null_coalescing_parent, .null_coalescing_search_dir)) break
    .null_coalescing_search_dir <- .null_coalescing_parent
  }
  if (nzchar(.null_coalescing_file)) break
}
if (!nzchar(.null_coalescing_file)) {
  stop("Could not locate R/null_coalescing.R.", call. = FALSE)
}
source(.null_coalescing_file)
rm(
  .null_coalescing_source_files, .null_coalescing_search_dirs,
  .null_coalescing_start, .null_coalescing_search_dir,
  .null_coalescing_candidate, .null_coalescing_parent,
  .null_coalescing_file
)

repo_root <- function() {
  env_root <- Sys.getenv("PROTEOMICS_PROJECT_ROOT", unset = "")
  if (nzchar(env_root)) {
    return(normalizePath(env_root, winslash = "/", mustWork = FALSE))
  }

  if (requireNamespace("rprojroot", quietly = TRUE)) {
    root <- tryCatch(
      rprojroot::find_root(rprojroot::has_file("README.md") | rprojroot::is_git_root),
      error = function(e) NULL
    )
    if (!is.null(root)) return(normalizePath(root, winslash = "/", mustWork = FALSE))
  }

  cur <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  repeat {
    markers <- c(".git", "README.md", "pipeline.yml", "analysis")
    if (any(file.exists(file.path(cur, markers)))) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) {
      stop("Could not locate repository root. Set PROTEOMICS_PROJECT_ROOT.", call. = FALSE)
    }
    cur <- parent
  }
}

# R/ is organised into domain subdirectories (data_contracts, qc, statistics,
# spatial, enrichment, networks, utilities). Call sites address a library by
# bare name -- repo_path("R", "module_stats.R") -- so the domain layout can be
# changed without editing the 800+ source() lines that name these libraries.
# This resolver maps a bare library name onto its domain directory and is the
# minimum path abstraction the layout requires.
.r_library_cache <- new.env(parent = emptyenv())

r_library_path <- function(name, root = repo_root()) {
  flat <- file.path(root, "R", name)
  if (file.exists(flat)) return(flat)

  key <- paste0("index:", root)
  if (!exists(key, envir = .r_library_cache, inherits = FALSE)) {
    rel <- list.files(
      file.path(root, "R"), pattern = "[.]R$", recursive = TRUE
    )
    index <- stats::setNames(file.path(root, "R", rel), basename(rel))
    assign(key, index[!duplicated(names(index))], envir = .r_library_cache)
  }
  index <- get(key, envir = .r_library_cache, inherits = FALSE)

  hit <- unname(index[name])
  if (!is.na(hit)) return(hit)
  flat
}

repo_path <- function(...) {
  parts <- list(...)
  if (length(parts) == 2L &&
      identical(as.character(parts[[1]]), "R") &&
      length(parts[[2]]) == 1L &&
      grepl("[.]R$", as.character(parts[[2]]))) {
    return(r_library_path(as.character(parts[[2]])))
  }
  file.path(repo_root(), ...)
}

path_raw <- function(...) repo_path("data", "raw", ...)
path_metadata <- function(...) repo_path("data", "metadata", ...)
path_external <- function(...) repo_path("data", "external", ...)
path_processed <- function(...) repo_path("data", "processed", ...)
path_results <- function(...) repo_path("results", ...)

# --- output lifecycles (config/output_layout.yml) -------------------------
# work/    regenerable intermediates, never cited, never canonical
# results/ canonical scientific results
# exports/ frozen outward-facing bundles, the only thing the manuscript reads
path_work <- function(...) repo_path("work", ...)
path_export <- function(...) repo_path("exports", ...)

.output_layout_cache <- new.env(parent = emptyenv())

output_layout <- function() {
  if (!is.null(.output_layout_cache$layout)) return(.output_layout_cache$layout)
  f <- repo_path("config", "output_layout.yml")
  if (!file.exists(f) || !requireNamespace("yaml", quietly = TRUE)) return(NULL)
  .output_layout_cache$layout <- yaml::read_yaml(f)
  .output_layout_cache$layout
}

output_layout_domains <- function() {
  l <- output_layout()
  if (is.null(l)) return(character(0))
  as.character(unlist(l$domains, use.names = FALSE))
}

output_layout_children <- function() {
  l <- output_layout()
  if (is.null(l)) return(character(0))
  names(l$canonical_result_layout$children_in_use)
}

# results/<domain>/<analysis_id>/<scope>/<child>/...
#
# analysis_id is the canonical owner's script stem, so an output is addressed
# by the analysis that owns it rather than by a historical stage number. The
# domain and child are checked against the contract, because a typo here would
# silently create a sibling namespace that no reader would ever find.
canonical_result_path <- function(domain, analysis_id, scope = "global",
                                  child = "tables", ...) {
  domains <- output_layout_domains()
  children <- output_layout_children()
  if (length(domains) && !domain %in% domains) {
    stop("unknown output domain '", domain, "'; config/output_layout.yml declares: ",
         paste(domains, collapse = ", "), call. = FALSE)
  }
  if (length(children) && !child %in% children) {
    stop("unknown result child '", child, "'; config/output_layout.yml declares: ",
         paste(children, collapse = ", "), call. = FALSE)
  }
  analysis_id <- sub("[.][Rr]$", "", basename(analysis_id))
  if (!nzchar(scope)) scope <- "global"
  path_results(domain, analysis_id, scope, child, ...)
}

canonical_work_path <- function(domain, analysis_id, scope = "global", ...) {
  analysis_id <- sub("[.][Rr]$", "", basename(analysis_id))
  if (!nzchar(scope)) scope <- "global"
  path_work(domain, analysis_id, scope, ...)
}

# The normalized replacement for create_module_dirs(), which keyed output on
# stage identity. One call gives a writer every destination it needs, so
# repointing a writer is one edit rather than one per write site.
#
# `suffix` is a segment below the child, for a writer whose historical layout
# nested one: build_spatial_networks wrote
# network_spatial_relations/<dataset>/<spatial_unit>, and the spatial unit
# carries meaning, so it is preserved.
#
# figures and logs are aliases, because renaming a key to plots or manifests is
# a pure rename with no judgement in it.
#
# There is deliberately no `processed` alias. Under the historical layout every
# writer put both its persistent objects and its scratch files under
# data/processed, and those two now diverge: a consumed object belongs in
# models/ and a disposable one in work/. Aliasing the old key would silently
# pick one, and that choice is exactly what section 2 of the Phase 6G.2 brief
# requires a human to make per artefact. A writer must therefore name $models
# or $work.
canonical_module_dirs <- function(domain, analysis_id, scope = "global",
                                  suffix = NULL, create = TRUE) {
  aid <- sub("[.][Rr]$", "", basename(analysis_id))
  if (!nzchar(scope)) scope <- "global"
  join <- function(p) if (is.null(suffix) || !length(suffix) || !nzchar(suffix)) p else file.path(p, suffix)

  out <- list(
    tables    = join(canonical_result_path(domain, aid, scope, "tables")),
    plots     = join(canonical_result_path(domain, aid, scope, "plots")),
    models    = join(canonical_result_path(domain, aid, scope, "models")),
    manifests = join(canonical_result_path(domain, aid, scope, "manifests")),
    reports   = join(canonical_result_path(domain, aid, scope, "reports")),
    work      = join(canonical_work_path(domain, aid, scope))
  )
  out$source_data <- file.path(out$tables, "source_data")
  out$figures <- out$plots
  out$logs <- out$manifests

  if (isTRUE(create)) invisible(lapply(out, dir_create))
  out
}

# --- legacy output roots --------------------------------------------------
# Registered in config/legacy_output_registry.csv: real artefacts that no
# registered writer produces any more. Reads are fine, writes are not.
legacy_output_roots <- function() {
  f <- repo_path("config", "legacy_output_registry.csv")
  if (!file.exists(f)) return(character(0))
  d <- utils::read.csv(f, stringsAsFactors = FALSE)
  d$legacy_path[d$policy == "LEGACY_READ_ONLY"]
}

is_legacy_output_path <- function(path) {
  roots <- legacy_output_roots()
  if (!length(roots)) return(rep(FALSE, length(path)))
  rel <- relative_to(path)
  vapply(rel, function(p) any(p == roots | startsWith(p, paste0(roots, "/"))),
         logical(1), USE.NAMES = FALSE)
}

dir_create <- function(...) {
  path <- file.path(...)
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

ensure_dir <- dir_create

safe_filename <- function(x, max_chars = 180) {
  x <- as.character(x)
  x <- gsub("[/\\\\:*?\"<>|]+", "_", x)
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x <- ifelse(nzchar(x), x, "unnamed")
  substr(x, 1, max_chars)
}

# --- input addressability ---------------------------------------------------
#
# Windows MAX_PATH is 260 characters and R cannot open a path at or beyond it.
# Measured in this repository with LongPathsEnabled = 0: every sampled path at
# <= 259 characters opened, none at >= 260 did, and file.copy() returns FALSE
# with a warning rather than raising. list.files() still enumerates them, so an
# over-limit file is visible but unreadable, and a relative path from a deep
# working directory does not escape the limit either - the resolved path is
# what counts.
#
# The consequence is that file.exists() alone cannot separate three materially
# different conditions:
#
#   * the declared root is not mounted in this session. Manifests written while
#     the repository was addressed through a substituted P:/ root are the known
#     case; the artefacts may be perfectly fine, just not reachable by the
#     address recorded for them.
#   * the file is present and enumerable, but its resolved path reaches the
#     wall, so R cannot open it.
#   * the file genuinely is not there.
#
# All three previously reported as "missing", which let the provenance ledger
# assert something untrue about the run. These helpers keep them distinct.
# Classification is diagnosis only: it does not decide whether a contract
# passes, and a caller that required a file still fails when it is unusable.

PATH_LENGTH_WALL <- 260L

# The same measurement the clusterProfiler preflight uses, so the two cannot
# drift: normalizePath() with mustWork = FALSE, counted in characters.
path_length_chars <- function(path) {
  nchar(normalizePath(path, winslash = "/", mustWork = FALSE), type = "chars")
}

# The root a path declares for itself: a drive letter, or a UNC //server/share.
# NA means the path declares no root of its own and is therefore addressed
# relative to the repository.
path_declared_root <- function(path) {
  path <- as.character(path)
  out <- rep(NA_character_, length(path))
  if (!length(path)) return(out)
  usable <- !is.na(path) & nzchar(path)
  drive <- usable & grepl("^[A-Za-z]:[/\\\\]", path)
  out[drive] <- paste0(substr(path[drive], 1, 2), "/")
  unc <- usable & !drive & grepl("^[/\\\\]{2}[^/\\\\]+[/\\\\][^/\\\\]+", path)
  if (any(unc)) {
    out[unc] <- gsub("\\\\", "/",
                     sub("^([/\\\\]{2}[^/\\\\]+[/\\\\][^/\\\\]+).*$", "\\1", path[unc]))
  }
  out
}

# A path with no declared root is addressed against the repository root, which
# is available by construction, so absence of a declared root is not a failure.
path_declared_root_available <- function(path) {
  root <- path_declared_root(path)
  ifelse(is.na(root), TRUE, dir.exists(root))
}

INPUT_STATUS_PRESENT <- "present"
INPUT_STATUS_ROOT_UNMOUNTED <- "declared_root_unmounted"
INPUT_STATUS_OVER_LIMIT <- "path_over_limit"
INPUT_STATUS_ABSENT <- "absent"
INPUT_STATUS_LEVELS <- c(
  INPUT_STATUS_PRESENT, INPUT_STATUS_ROOT_UNMOUNTED,
  INPUT_STATUS_OVER_LIMIT, INPUT_STATUS_ABSENT
)

# Precedence is fixed and total, and the order is forced by what each test can
# honestly answer:
#   1. declared root availability, because nothing below it can be measured
#      when the root the path names is not mounted;
#   2. the path budget, because file.exists() is not trustworthy at or beyond
#      the wall - it returns FALSE for files that demonstrably exist;
#   3. existence, consulted only once the address is known to be usable.
# An empty or NA path declares nothing and cannot be resolved, so it is
# reported as absent rather than given a class it has not earned.
input_addressability <- function(path, wall = PATH_LENGTH_WALL) {
  path <- as.character(path)
  out <- rep(NA_character_, length(path))
  if (!length(path)) return(out)

  blank <- is.na(path) | !nzchar(trimws(path))
  out[blank] <- INPUT_STATUS_ABSENT
  idx <- which(!blank)
  if (!length(idx)) return(out)

  unmounted <- !path_declared_root_available(path[idx])
  out[idx[unmounted]] <- INPUT_STATUS_ROOT_UNMOUNTED
  idx <- idx[!unmounted]
  if (!length(idx)) return(out)

  over <- path_length_chars(path[idx]) >= as.integer(wall)
  out[idx[over]] <- INPUT_STATUS_OVER_LIMIT
  idx <- idx[!over]
  if (!length(idx)) return(out)

  hit <- file.exists(path[idx]) | dir.exists(path[idx])
  out[idx[hit]] <- INPUT_STATUS_PRESENT
  out[idx[!hit]] <- INPUT_STATUS_ABSENT
  out
}

# The one compatibility predicate. Anything that used to ask file.exists() and
# meant "can I use this input" asks this instead, and it is defined as exactly
# status == present so a Boolean can never disagree with the classification.
input_is_present <- function(path, wall = PATH_LENGTH_WALL) {
  input_addressability(path, wall = wall) == INPUT_STATUS_PRESENT
}

input_status_message <- function(status) {
  vapply(as.character(status), function(s) switch(
    s,
    present = "input available",
    declared_root_unmounted =
      "declared root is not mounted in this session; the input was not looked for",
    path_over_limit = paste0(
      "path reaches the ", PATH_LENGTH_WALL,
      "-character limit and cannot be opened by R; presence is undetermined"),
    absent = "input not available",
    "input status unknown"
  ), character(1), USE.NAMES = FALSE)
}

# Groups unusable paths by class so a failure message says which of the three
# failures happened, instead of calling all of them "missing".
describe_input_status_failures <- function(path, status = NULL,
                                           max_shown = 3L) {
  path <- as.character(path)
  status <- status %||% input_addressability(path)
  bad <- which(status != INPUT_STATUS_PRESENT)
  if (!length(bad)) return("")
  label <- c(
    declared_root_unmounted = "declared root not mounted",
    path_over_limit = paste0("path at or beyond the ", PATH_LENGTH_WALL,
                             "-character limit (present but unopenable)"),
    absent = "genuinely absent"
  )
  parts <- character(0)
  for (cls in names(label)) {
    hit <- bad[status[bad] == cls]
    if (!length(hit)) next
    shown <- path[utils::head(hit, max_shown)]
    more <- length(hit) - length(shown)
    parts <- c(parts, paste0(
      label[[cls]], " (", length(hit), "): ",
      paste(shown, collapse = ", "),
      if (more > 0L) paste0(" [+", more, " more]") else ""
    ))
  }
  paste(parts, collapse = "; ")
}

# --- runtime addressability staging -----------------------------------------
#
# Some inputs are recorded under an address that cannot be used here, and the
# corrected address is too long for R to open. The historical enrichment
# manifests are the case in hand: every path is stored under a substituted P:/
# root, and re-anchoring the short stored suffix on this 69-character
# repository root pushes part of it past the wall.
#
# The answer is to keep the declared path as provenance and derive a usable
# runtime path, staging a byte-identical copy into the regenerable work/
# lifecycle only where the re-anchored path actually crosses the wall.
#
# Staging is addressability, never transformation: no conversion, filtering,
# renaming by meaning, recompression or schema change. The staged basename is
# copied verbatim, so an extension cannot be clipped.
#
# Provenance keeps pointing at the declared scientific artifact. Nothing under
# work/ may be cited, and the layout contract already says so.

RUNTIME_RESOLUTION_DIRECT <- "direct"
RUNTIME_RESOLUTION_REBASED <- "rebased"
RUNTIME_RESOLUTION_STAGED <- "staged"
RUNTIME_RESOLUTION_UNRESOLVED <- "unresolved"
RUNTIME_RESOLUTION_NOT_CONSUMED <- "not_consumed"
RUNTIME_RESOLUTION_LEVELS <- c(
  RUNTIME_RESOLUTION_DIRECT, RUNTIME_RESOLUTION_REBASED,
  RUNTIME_RESOLUTION_STAGED, RUNTIME_RESOLUTION_UNRESOLVED,
  RUNTIME_RESOLUTION_NOT_CONSUMED
)

# R cannot stat, open or hash a path at or beyond the wall, so the facts about
# a long source have to come from a runtime that can. .NET Core is
# extended-length aware; this is the same reason the existing robocopy stager
# in audits/part29 works. The helper is written to a temporary file rather than
# passed through -Command because the paths would otherwise have to survive two
# layers of shell quoting.
.os_path_helper <- local({
  cached <- NULL
  function() {
    if (!is.null(cached) && file.exists(cached)) return(cached)
    f <- tempfile("os_path_helper_", fileext = ".ps1")
    writeLines(c(
      "param([string]$Mode, [string]$InFile, [string]$OutFile)",
      "$tab = [string][char]9",
      "$w = [System.IO.StreamWriter]::new($OutFile, $false, [System.Text.UTF8Encoding]::new($false))",
      "foreach ($line in [System.IO.File]::ReadLines($InFile)) {",
      "  if ([string]::IsNullOrWhiteSpace($line)) { continue }",
      "  $parts = $line.Split([char]9)",
      "  $src = $parts[0]",
      "  if ($Mode -eq 'facts') {",
      "    $e = [System.IO.File]::Exists($src)",
      "    $len = -1; $mt = ''; $sha = ''",
      "    if ($e) {",
      "      try {",
      "        $fi = [System.IO.FileInfo]::new($src)",
      "        $len = $fi.Length",
      "        $mt = $fi.LastWriteTimeUtc.ToString('o')",
      "        $sha = (Get-FileHash -LiteralPath $src -Algorithm SHA256).Hash",
      "      } catch { $len = -1 }",
      "    }",
      "    $w.WriteLine($src + $tab + $e.ToString() + $tab + $len + $tab + $mt + $tab + $sha)",
      "  } elseif ($Mode -eq 'copy') {",
      "    $dst = $parts[1]",
      "    $ok = 'FALSE'; $msg = ''",
      "    try {",
      "      $d = [System.IO.Path]::GetDirectoryName($dst)",
      "      if (-not [System.IO.Directory]::Exists($d)) { [void][System.IO.Directory]::CreateDirectory($d) }",
      "      [System.IO.File]::Copy($src, $dst, $false)",
      "      $ok = 'TRUE'",
      "    } catch { $msg = $_.Exception.Message }",
      "    $w.WriteLine($src + $tab + $dst + $tab + $ok + $tab + $msg)",
      "  }",
      "}",
      "$w.Close()"
    ), f)
    cached <<- f
    f
  }
})

.os_path_run <- function(mode, lines) {
  helper <- .os_path_helper()
  infile <- tempfile("os_path_in_", fileext = ".txt")
  outfile <- tempfile("os_path_out_", fileext = ".txt")
  on.exit(unlink(c(infile, outfile)), add = TRUE)
  writeLines(lines, infile, useBytes = TRUE)
  status <- suppressWarnings(system2(
    "pwsh", c("-NoProfile", "-File", helper, "-Mode", mode,
              "-InFile", infile, "-OutFile", outfile),
    stdout = TRUE, stderr = TRUE))
  if (!file.exists(outfile)) {
    stop("OS path helper did not run (mode=", mode, "): ",
         paste(utils::head(status, 5), collapse = " | "), call. = FALSE)
  }
  readLines(outfile, warn = FALSE)
}

# exists / bytes / mtime / sha256 for paths of any length.
os_path_facts <- function(paths) {
  paths <- as.character(paths)
  empty <- data.frame(path = character(), exists = logical(), bytes = numeric(),
                      mtime_utc = character(), sha256 = character(),
                      stringsAsFactors = FALSE)
  keep <- !is.na(paths) & nzchar(paths)
  if (!any(keep)) return(empty)
  raw <- .os_path_run("facts", unique(paths[keep]))
  sp <- strsplit(raw, "\t", fixed = TRUE)
  sp <- sp[vapply(sp, length, 0L) >= 3L]
  out <- data.frame(
    path = vapply(sp, `[`, "", 1L),
    exists = toupper(vapply(sp, `[`, "", 2L)) == "TRUE",
    bytes = suppressWarnings(as.numeric(vapply(sp, `[`, "", 3L))),
    mtime_utc = vapply(sp, function(x) if (length(x) >= 4L) x[[4]] else "", ""),
    sha256 = tolower(vapply(sp, function(x) if (length(x) >= 5L) x[[5]] else "", "")),
    stringsAsFactors = FALSE)
  out[match(paths, out$path), , drop = FALSE]
}

# Copy that refuses to overwrite. Staleness is resolved by the caller, which
# can see the hashes; silently replacing a mismatched staged file would hide
# exactly the problem worth knowing about.
os_copy_no_clobber <- function(src, dst) {
  stopifnot(length(src) == length(dst))
  if (!length(src)) return(logical(0))
  raw <- .os_path_run("copy", paste(src, dst, sep = "\t"))
  sp <- strsplit(raw, "\t", fixed = TRUE)
  ok <- setNames(
    toupper(vapply(sp, function(x) if (length(x) >= 3L) x[[3]] else "FALSE", "")) == "TRUE",
    vapply(sp, `[`, "", 1L))
  unname(ok[src])
}

# The staging root lives inside the declared work/ lifecycle and is owned by
# the analysis that produced the artifacts, so it cannot drift into being a
# second scientific namespace.
path_stage_root <- function(domain, analysis_id, scope = "global", create = FALSE) {
  root <- canonical_work_path(domain, analysis_id, scope, "path_stage")
  if (isTRUE(create)) dir_create(root)
  root
}

# Deterministic destination: one source maps to exactly one staged path.
# The digest is taken over the DECLARED path, the stable provenance identifier,
# and the basename is reused verbatim so the extension survives untouched.
# A per-digest directory avoids both collisions and rebuilding the original
# deep hierarchy.
staged_destination <- function(declared_path, root, hash_chars = 12L) {
  declared_path <- as.character(declared_path)
  digest <- vapply(declared_path, function(p) {
    h <- if (requireNamespace("digest", quietly = TRUE)) {
      digest::digest(p, algo = "sha256", serialize = FALSE)
    } else {
      f <- tempfile("stage_key_")
      on.exit(unlink(f), add = TRUE)
      writeBin(charToRaw(p), f)
      unname(tools::sha256sum(f))
    }
    substr(h, 1L, hash_chars)
  }, character(1), USE.NAMES = FALSE)
  file.path(root, digest, basename(declared_path))
}

# Stage byte-identical copies for addressability, idempotently.
#
# A destination that already holds the expected bytes is reused. A destination
# that holds anything else is reported as a mismatch and is NOT overwritten:
# silently replacing it would hide the one condition worth knowing about, and
# silently accepting it would feed stale content to the science.
#
# The source is treated as read-only and is checked for that: size, mtime and
# hash are compared either side of the copy.
stage_addressable_copies <- function(source, destination, expected_sha256 = NULL) {
  source <- as.character(source)
  destination <- as.character(destination)
  stopifnot(length(source) == length(destination))
  n <- length(source)
  out <- data.frame(
    source = source, staged_path = destination,
    action = rep(NA_character_, n), staged_sha256 = NA_character_,
    source_sha256 = if (is.null(expected_sha256)) NA_character_ else as.character(expected_sha256),
    source_bytes = NA_real_, staged_bytes = NA_real_,
    same_hash = FALSE, same_size = FALSE, source_unchanged = NA,
    ok = FALSE, stringsAsFactors = FALSE)
  if (!n) return(out)

  before <- os_path_facts(source)
  if (is.null(expected_sha256)) out$source_sha256 <- before$sha256
  out$source_bytes <- before$bytes

  present <- file.exists(destination)
  if (any(present)) {
    i <- which(present)
    have <- tolower(vapply(destination[i], function(p) unname(tools::sha256sum(p)),
                           character(1), USE.NAMES = FALSE))
    out$staged_sha256[i] <- have
    match_i <- have == tolower(out$source_sha256[i])
    out$action[i] <- ifelse(match_i, "reused", "mismatch")
  }
  todo <- which(!present)
  if (length(todo)) {
    ok <- os_copy_no_clobber(source[todo], destination[todo])
    out$action[todo] <- ifelse(ok, "copied", "copy_failed")
    got <- todo[file.exists(destination[todo])]
    if (length(got)) {
      out$staged_sha256[got] <- tolower(vapply(destination[got],
        function(p) unname(tools::sha256sum(p)), character(1), USE.NAMES = FALSE))
    }
  }

  staged_ok <- file.exists(out$staged_path)
  out$staged_bytes[staged_ok] <- file.size(out$staged_path[staged_ok])
  out$same_hash <- !is.na(out$staged_sha256) & !is.na(out$source_sha256) &
    tolower(out$staged_sha256) == tolower(out$source_sha256)
  out$same_size <- !is.na(out$staged_bytes) & !is.na(out$source_bytes) &
    out$staged_bytes == out$source_bytes

  after <- os_path_facts(source)
  out$source_unchanged <- !is.na(after$sha256) & !is.na(before$sha256) &
    after$sha256 == before$sha256 & after$bytes == before$bytes &
    after$mtime_utc == before$mtime_utc

  out$ok <- out$action %in% c("copied", "reused") & out$same_hash & out$same_size &
    out$source_unchanged
  out
}

# --- budgeted figure targets ------------------------------------------------
#
# A figure filename has to be budgeted against the ABSOLUTE path, not against
# a fixed basename length, because the space available depends on how deep the
# output directory is. The WGCNA module-score directories measure 190-198
# characters, so a basename limit that looks generous in isolation still
# overruns; conversely the canonical successor root is 127 characters and has
# ample room. One number cannot serve both.
#
# What went wrong before: three ggsave() sites built their target with paste0()
# and no budget at all. The underlying Windows file API then truncated the
# path at MAX_PATH-1 instead of failing, so 354 figures landed on disk with
# ".svg" cut off - and the manuscript exporter selects on \.(svg|pdf|png)$, so
# every one of them silently dropped out of publication discovery. A naming
# bug became a publication-selection bug.
#
# Two budgets, deliberately distinct:
#   FIGURE_WRITE_BUDGET (240) is the conservative budget for NEW writes.
#   PATH_LENGTH_WALL    (260) is the empirical wall R cannot read past.
FIGURE_WRITE_BUDGET <- 240L

# The filename budget left by a directory, under a chosen total budget.
figure_filename_budget <- function(directory, budget = FIGURE_WRITE_BUDGET) {
  as.integer(budget) - path_length_chars(directory) - 1L
}

# Fit `filename` into the space `directory` leaves, keeping the extension whole.
#
# `taken` lets a caller pass names already claimed in that directory so a
# shortened stem cannot collide; disambiguation appends the repository's usual
# short stable digest, and the extension is reattached after it.
budgeted_figure_target <- function(directory, filename,
                                   budget = FIGURE_WRITE_BUDGET,
                                   taken = character(0), digest_chars = 8L) {
  filename <- as.character(filename)
  room <- figure_filename_budget(directory, budget)
  ext <- tools::file_ext(filename)
  ext_part <- if (nzchar(ext)) paste0(".", ext) else ""
  stem <- tools::file_path_sans_ext(filename)

  ## An extension is atomic. If the directory leaves no room for one stem
  ## character plus the whole extension, that is a contract error: emitting a
  ## clipped ".s" or ".sv" is exactly the failure this function exists to stop.
  if (room < nchar(ext_part, type = "chars") + 1L) {
    stop("Figure directory leaves ", room, " characters for a filename, which cannot hold ",
         "a stem character plus the extension '", ext_part, "'. Directory: ", directory,
         call. = FALSE)
  }

  out <- if (nchar(filename, type = "chars") <= room) filename else
    paste0(substr(stem, 1L, room - nchar(ext_part, type = "chars")), ext_part)

  if (out %in% taken) {
    suffix <- paste0("__", substr(
      if (requireNamespace("digest", quietly = TRUE))
        digest::digest(file.path(directory, filename), algo = "sha256", serialize = FALSE)
      else sprintf("%08x", sum(utf8ToInt(file.path(directory, filename)))),
      1L, digest_chars))
    keep <- room - nchar(suffix, type = "chars") - nchar(ext_part, type = "chars")
    if (keep < 1L) {
      stop("Figure directory leaves no room to disambiguate '", filename, "' in ",
           directory, call. = FALSE)
    }
    out <- paste0(substr(tools::file_path_sans_ext(out), 1L, keep), suffix, ext_part)
  }
  out
}

# The full target, for a caller that just wants a path it can write.
budgeted_figure_path <- function(directory, filename,
                                 budget = FIGURE_WRITE_BUDGET, taken = character(0)) {
  file.path(directory, budgeted_figure_target(directory, filename, budget, taken))
}

file_hash <- function(path) {
  if (is.null(path) || !length(path) || is.na(path) || !file.exists(path)) return(NA_character_)
  unname(tools::md5sum(path))
}

file_hash_sha256 <- function(path) {
  if (
    is.null(path) ||
      !length(path) ||
      is.na(path) ||
      !file.exists(path) ||
      dir.exists(path)
  ) {
    return(NA_character_)
  }

  if ("sha256sum" %in% getNamespaceExports("tools")) {
    return(unname(tools::sha256sum(path)))
  }

  if (requireNamespace("digest", quietly = TRUE)) {
    return(unname(digest::digest(
      file = path,
      algo = "sha256"
    )))
  }

  warning(
    "SHA-256 hashing is unavailable. Install the R package 'digest'.",
    call. = FALSE
  )

  NA_character_
}

strict_inputs_enabled <- function(config = NULL) {
  args <- c(commandArgs(trailingOnly = FALSE), commandArgs(trailingOnly = TRUE))
  from_args <- "--strict-inputs" %in% args
  from_env <- tolower(Sys.getenv("PROTEOMICS_STRICT_INPUTS", unset = "")) %in% c("1", "true", "yes")
  from_config <- FALSE
  if (is.list(config)) {
    from_config <- isTRUE(config$strict_inputs) ||
      isTRUE(config$runtime$strict_inputs)
  }
  isTRUE(from_args || from_env || from_config)
}

input_resolution_audit_path <- function() {
  path_results("reviewer_audit", "input_resolution_audit.csv")
}

input_resolution_audit_columns <- function() {
  c(
    "script", "dataset", "stage", "input_name", "expected_path", "resolved_path",
    "resolution_mode", "strict_mode", "allowed_in_strict_mode", "file_exists",
    "file_hash_sha256", "file_mtime", "producer_script_or_artifact_id", "warning"
  )
}

input_file_mtime <- function(path) {
  if (is.null(path) || !length(path) || is.na(path) || !file.exists(path)) return(NA_character_)
  format(file.info(path)$mtime[[1]], "%Y-%m-%d %H:%M:%S %z")
}

append_input_resolution_audit <- function(rows, path = input_resolution_audit_path()) {
  if (is.null(rows) || !length(rows)) return(invisible(path))
  rows <- as.data.frame(rows, stringsAsFactors = FALSE)
  cols <- input_resolution_audit_columns()
  for (col in setdiff(cols, names(rows))) rows[[col]] <- NA
  rows <- rows[, cols, drop = FALSE]
  dir_create(dirname(path))
  write_header <- !file.exists(path) || file.info(path)$size == 0

  # Render first, append once. write.table(append = TRUE) issues several
  # writes per call, and this ledger is appended to by concurrently running
  # analyses. Three records in the file on disk are spliced remnants of that,
  # and the three stray quote characters they leave behind desynchronise CSV
  # quoting for every one of the 50,510 lines that follow: 0.005% of the rows
  # cost 97% of the ledger. Rendering to a buffer keeps write.table's exact
  # formatting and quoting but narrows the write to a single call and makes
  # the block checkable before it is committed.
  buffer <- character(0)
  render <- textConnection("buffer", open = "w", local = TRUE)
  utils::write.table(
    rows,
    file = render,
    sep = ",",
    row.names = FALSE,
    col.names = write_header,
    na = "",
    qmethod = "double"
  )
  close(render)

  # A record carrying an odd number of quote characters is not merely a bad
  # row: it reopens a quoted field and swallows the rest of the file.
  # Dropping it loses one row, admitting it loses the ledger.
  balanced <- vapply(buffer, function(line) {
    hit <- gregexpr("\"", line, fixed = TRUE)[[1]]
    (if (hit[[1]] == -1L) 0L else length(hit)) %% 2L == 0L
  }, logical(1), USE.NAMES = FALSE)
  if (any(!balanced)) {
    warning("Dropped ", sum(!balanced), " unbalanced input-resolution audit ",
      "record(s); admitting them would have made ", basename(path),
      " unparseable from that point on.", call. = FALSE)
    buffer <- buffer[balanced]
  }
  if (!length(buffer)) return(invisible(path))

  # Text mode, matching what write.table(file = path) did, so the ledger keeps
  # the line endings it already has.
  con <- file(path, open = if (write_header) "wt" else "at")
  on.exit(close(con), add = TRUE)
  writeLines(buffer, con)
  invisible(path)
}

record_input_resolution <- function(
  script = Sys.getenv("PROTEOMICS_SCRIPT_ID", unset = NA_character_),
  dataset = Sys.getenv("PROTEOMICS_DATASET", unset = NA_character_),
  stage = NA_character_,
  input_name,
  expected_path = NA_character_,
  resolved_path = NA_character_,
  resolution_mode = "canonical",
  strict_mode = strict_inputs_enabled(),
  allowed_in_strict_mode = TRUE,
  producer_script_or_artifact_id = NA_character_,
  warning = NA_character_
) {
  one_path <- function(path) {
    path <- as.character(path)
    path <- path[!is.na(path) & nzchar(path)]
    if (!length(path)) return(NA_character_)
    normalizePath(path[[1]], winslash = "/", mustWork = FALSE)
  }
  expected_path <- one_path(expected_path)
  resolved_path <- one_path(resolved_path)
  append_input_resolution_audit(data.frame(
    script = script,
    dataset = dataset,
    stage = stage,
    input_name = input_name,
    expected_path = expected_path,
    resolved_path = resolved_path,
    resolution_mode = resolution_mode,
    strict_mode = isTRUE(strict_mode),
    allowed_in_strict_mode = isTRUE(allowed_in_strict_mode),
    file_exists = !is.na(resolved_path) && file.exists(resolved_path),
    file_hash_sha256 = file_hash_sha256(resolved_path),
    file_mtime = input_file_mtime(resolved_path),
    producer_script_or_artifact_id = producer_script_or_artifact_id,
    warning = warning,
    stringsAsFactors = FALSE
  ))
}

latest_input_candidate <- function(roots, pattern, recursive = TRUE) {
  roots <- roots[!is.na(roots) & nzchar(roots)]
  if (!length(roots) || is.na(pattern) || !nzchar(pattern)) return(NA_character_)
  files <- unlist(lapply(roots, function(root) {
    root <- normalizePath(root, winslash = "/", mustWork = FALSE)
    ## An unmounted declared root and an empty directory are not the same
    ## thing, and dir.exists() cannot tell them apart. Enumerating a root
    ## that was never reachable would silently narrow the candidate set.
    if (!path_declared_root_available(root) || !dir.exists(root)) return(character())
    list.files(root, pattern = pattern, full.names = TRUE, recursive = recursive)
  }), use.names = FALSE)
  files <- files[input_is_present(files)]
  if (!length(files)) return(NA_character_)
  info <- file.info(files)
  normalizePath(rownames(info)[order(info$mtime, decreasing = TRUE)[1]], winslash = "/", mustWork = FALSE)
}

resolve_input_path <- function(
  input_name,
  expected_path = NA_character_,
  explicit_path = NA_character_,
  fallback_paths = character(),
  latest_roots = character(),
  latest_pattern = NA_character_,
  recursive = TRUE,
  required = TRUE,
  script = Sys.getenv("PROTEOMICS_SCRIPT_ID", unset = NA_character_),
  dataset = Sys.getenv("PROTEOMICS_DATASET", unset = NA_character_),
  stage = NA_character_,
  producer_script_or_artifact_id = NA_character_,
  allow_fallback_in_strict = FALSE,
  allow_latest_in_strict = FALSE,
  record_resolution = TRUE
) {
  strict <- strict_inputs_enabled()
  norm <- function(x) {
    x <- x[!is.na(x) & nzchar(x)]
    if (!length(x)) return(character())
    normalizePath(x, winslash = "/", mustWork = FALSE)
  }
  expected_path <- norm(expected_path)[1]
  if (is.na(expected_path)) expected_path <- NA_character_
  explicit_path <- norm(explicit_path)[1]
  if (is.na(explicit_path)) explicit_path <- NA_character_
  fallback_paths <- norm(fallback_paths)

  finish <- function(resolved, mode, allowed, warn = NA_character_) {
    if (!is.na(warn) && nzchar(warn)) warning(warn, call. = FALSE)
    if (isTRUE(record_resolution)) {
      record_input_resolution(
        script = script,
        dataset = dataset,
        stage = stage,
        input_name = input_name,
        expected_path = expected_path,
        resolved_path = resolved,
        resolution_mode = mode,
        strict_mode = strict,
        allowed_in_strict_mode = allowed,
        producer_script_or_artifact_id = producer_script_or_artifact_id,
        warning = warn
      )
    }
    resolved
  }

  ## Every presence decision below goes through the four-state contract
  ## declared at the top of this file. It used to ask file.exists() here, and
  ## the resolution_mode tokens it emits are the second vocabulary that the
  ## contract was written to end: they still name the PRECEDENCE that was
  ## used, which is their job, but a failure now also carries the
  ## addressability state, so "missing" can no longer stand for an unmounted
  ## root or a path past the character limit.
  if (!is.na(explicit_path) && input_is_present(explicit_path)) {
    return(finish(explicit_path, "explicit_override", TRUE))
  }
  if (!is.na(explicit_path) && !input_is_present(explicit_path)) {
    explicit_state <- input_addressability(explicit_path)
    warn <- paste0("Explicit input override for ", input_name, " is not usable: ",
      input_status_message(explicit_state), ": ", explicit_path)
    finish(explicit_path, paste0("explicit_missing:", explicit_state), TRUE, warn)
    if (isTRUE(required)) stop(warn, call. = FALSE)
    return(explicit_path)
  }
  if (!is.na(expected_path) && input_is_present(expected_path)) {
    return(finish(expected_path, "canonical", TRUE))
  }

  ## A canonical input that is present-but-unopenable must not be silently
  ## replaced by a non-canonical fallback: that substitutes different data
  ## under the same name. Only genuine absence earns a fallback.
  expected_state <- if (is.na(expected_path)) NA_character_ else input_addressability(expected_path)
  if (!is.na(expected_state) && expected_state %in% c(INPUT_STATUS_ROOT_UNMOUNTED,
                                                      INPUT_STATUS_OVER_LIMIT)) {
    warn <- paste0("Canonical input for ", input_name, " is present but unopenable (",
      input_status_message(expected_state), "), so no fallback may stand in for it: ",
      expected_path)
    finish(NA_character_, paste0("canonical_undetermined:", expected_state), FALSE, warn)
    if (isTRUE(required)) stop(warn, call. = FALSE)
    return(NA_character_)
  }
  fallback_hit <- fallback_paths[input_is_present(fallback_paths)][1]
  if (!is.na(fallback_hit)) {
    warn <- paste0("Using non-canonical fallback for ", input_name, ": ", fallback_hit)
    if (isTRUE(strict) && !isTRUE(allow_fallback_in_strict)) {
      warn <- paste0("Strict input mode forbids fallback for ", input_name, ". Expected canonical input: ", expected_path)
      finish(fallback_hit, "fallback_forbidden_strict", FALSE, warn)
      if (isTRUE(required)) stop(warn, call. = FALSE)
      return(NA_character_)
    }
    return(finish(fallback_hit, "fallback", isTRUE(allow_fallback_in_strict), warn))
  }

  latest_hit <- latest_input_candidate(latest_roots, latest_pattern, recursive = recursive)
  if (!is.na(latest_hit)) {
    warn <- paste0("Using newest matching fallback for ", input_name, ": ", latest_hit)
    if (isTRUE(strict) && !isTRUE(allow_latest_in_strict)) {
      warn <- paste0("Strict input mode forbids newest-file fallback for ", input_name, ". Expected canonical input: ", expected_path)
      finish(latest_hit, "latest_forbidden_strict", FALSE, warn)
      if (isTRUE(required)) stop(warn, call. = FALSE)
      return(NA_character_)
    }
    return(finish(latest_hit, "latest_fallback", isTRUE(allow_latest_in_strict), warn))
  }

  warn <- paste0("Missing input for ", input_name,
    if (!is.na(expected_path)) paste0(": ", expected_path) else ".")
  finish(NA_character_,
    if (is.na(expected_state)) "missing" else paste0("missing:", expected_state),
    TRUE, if (isTRUE(required)) warn else NA_character_)
  if (isTRUE(required)) stop(warn, call. = FALSE)
  NA_character_
}

path_or_env <- function(env, default, must_exist = FALSE, kind = c("file", "dir", "any")) {
  kind <- match.arg(kind)
  value <- Sys.getenv(env, unset = "")
  path <- if (nzchar(value)) value else default
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)

  if (isTRUE(must_exist)) {
    exists_ok <- switch(
      kind,
      file = file.exists(path),
      dir = dir.exists(path),
      any = file.exists(path) || dir.exists(path)
    )
    if (!exists_ok) {
      stop(
        "Required ", kind, " path does not exist: ", path,
        ". Set ", env, " to override this location.",
        call. = FALSE
      )
    }
  }

  path
}

relative_to <- function(path, root = repo_root()) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  sub(paste0("^", gsub("([\\^$.|?*+(){}])", "\\\\\\1", root), "/?"), "", path)
}

write_session_info <- function(path = path_results("logs", "sessionInfo.txt")) {
  dir_create(dirname(path))
  capture.output(utils::sessionInfo(), file = path)
  invisible(path)
}

git_commit_sha <- function() {
  sha <- tryCatch(
    system2("git", c("-C", repo_root(), "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(e) NA_character_
  )
  sha <- if (length(sha)) as.character(sha[[1]]) else NA_character_
  if (is.na(sha) || !nzchar(sha)) NA_character_ else sha
}

run_context_metadata <- function() {
  env_flags <- c(
    "PROTEOMICS_DATASET",
    "PROTEOMICS_DRY_RUN",
    "PROTEOMICS_STRICT_INPUTS",
    "PROTEOMICS_RECOMPUTE",
    "PROTEOMICS_WGCNA_FORCE_FULL"
  )
  env <- as.list(Sys.getenv(env_flags, unset = NA_character_))
  names(env) <- env_flags
  list(
    git_commit = git_commit_sha(),
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    platform = R.version$platform,
    environment = env
  )
}

write_config_snapshot <- function(config, path) {
  dir_create(dirname(path))
  if (requireNamespace("yaml", quietly = TRUE)) {
    writeLines(yaml::as.yaml(config), path)
  } else {
    capture.output(str(config), file = path)
  }
  invisible(path)
}

# Provenance for one canonical result family, written into its manifests/
# child: config/output_layout.yml section 9.
#
# biological unit and statistical scope are deliberately not restated here.
# They are declared once, for the analyses that reach the manuscript, in
# docs/MANUSCRIPT_STATISTICAL_CONTRACT.md, which is generated from the frozen
# v9 contract. A manifest that paraphrased them would duplicate a frozen
# scientific artefact with no authority to do so, so it points at it instead.
write_result_manifest <- function(domain, analysis_id, scope = "global",
                                  inputs = list(), outputs = list(),
                                  parameters = list(), config_files = character(0),
                                  notes = NULL) {
  aid <- sub("[.][Rr]$", "", basename(analysis_id))
  path <- canonical_result_path(domain, aid, scope, "manifests", "manifest.yml")
  dir_create(dirname(path))

  hash_map <- function(x) {
    x <- as.character(unlist(x, use.names = TRUE))
    x <- x[!is.na(x) & nzchar(x)]
    if (!length(x)) return(list())
    stats::setNames(lapply(x, function(p)
      list(path = relative_to(p),
           sha256 = if (file.exists(p) && !dir.exists(p)) file_hash_sha256(p) else NA_character_)),
      names(x) %||% basename(x))
  }

  owner <- NA_character_
  own_file <- repo_path("config", "results_ownership.csv")
  if (file.exists(own_file)) {
    o <- utils::read.csv(own_file, stringsAsFactors = FALSE)
    hit <- grep(paste0("/", aid, "[.][Rr]$"), o$canonical_owner)
    if (length(hit)) owner <- o$canonical_owner[hit[1]]
  }

  manifest <- list(
    contract_version = "result_manifest_v1",
    analysis_id = aid,
    domain = domain,
    scope = scope,
    canonical_owner = if (is.na(owner)) paste0("analysis/", domain, "/", aid, ".R") else owner,
    source_commit = git_commit_sha(),
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    config = hash_map(config_files),
    upstream_inputs = hash_map(inputs),
    outputs = hash_map(outputs),
    parameters = parameters,
    statistical_contract = "docs/MANUSCRIPT_STATISTICAL_CONTRACT.md",
    notes = notes
  )
  if (requireNamespace("yaml", quietly = TRUE)) {
    writeLines(yaml::as.yaml(manifest), path, useBytes = TRUE)
  }
  invisible(path)
}

write_run_manifest <- function(path, inputs = list(), outputs = list(), parameters = list(), notes = NULL) {
  dir_create(dirname(path))

  flatten_paths <- function(x) {
    vals <- unlist(x, recursive = TRUE, use.names = TRUE)
    vals <- vals[!is.na(vals) & nzchar(as.character(vals))]
    as.character(vals)
  }

  input_paths <- flatten_paths(inputs)
  input_hashes <- lapply(input_paths, file_hash)
  names(input_hashes) <- names(input_paths)

  session_path <- file.path(dirname(path), "sessionInfo.txt")
  capture.output(utils::sessionInfo(), file = session_path)

  manifest <- list(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    repo_root = repo_root(),
    git_commit = git_commit_sha(),
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    platform = R.version$platform,
    environment = run_context_metadata()$environment,
    inputs = inputs,
    input_hashes = input_hashes,
    outputs = outputs,
    parameters = parameters,
    notes = notes,
    session_info = relative_to(session_path)
  )

  if (requireNamespace("yaml", quietly = TRUE)) {
    writeLines(yaml::as.yaml(manifest), path)
  } else {
    capture.output(str(manifest, max.level = 4), file = path)
  }

  invisible(path)
}

is_dry_run <- function(config = NULL) {
  from_args <- "--dry-run" %in% commandArgs(trailingOnly = FALSE) ||
    "--dry-run" %in% commandArgs(trailingOnly = TRUE)
  from_env <- tolower(Sys.getenv("PROTEOMICS_DRY_RUN", unset = "")) %in% c("1", "true", "yes")
  from_config <- FALSE
  if (is.list(config)) {
    from_config <- isTRUE(config$dry_run) ||
      isTRUE(config$runtime$dry_run)
  }
  isTRUE(from_args || from_env || from_config)
}

dry_run_line <- function(label, value = "", status = NULL) {
  prefix <- if (is.null(status)) "[DRY-RUN]" else paste0("[DRY-RUN ", status, "]")
  message(prefix, " ", label, if (nzchar(as.character(value))) paste0(": ", value) else "")
}

module_paths <- function(module, substep = NULL) {
  tail <- c(module, substep)
  tail <- tail[nzchar(tail)]
  list(
    processed = do.call(path_processed, as.list(tail)),
    figures = do.call(path_results, as.list(c("figures", tail))),
    tables = do.call(path_results, as.list(c("tables", tail))),
    source_data = do.call(path_results, as.list(c("source_data", tail))),
    logs = do.call(path_results, as.list(c("logs", tail))),
    reports = do.call(path_results, as.list(c("reports", tail)))
  )
}

create_module_dirs <- function(module, substep = NULL) {
  paths <- module_paths(module, substep)
  invisible(lapply(paths, dir_create))
  paths
}

pride_submission_dir <- function(...) repo_path("pride_submission", ...)

ensure_pride_dirs <- function() {
  dirs <- pride_submission_dir(c(
    "metadata",
    "processed_data",
    "supplementary_tables",
    "methods",
    "manifests",
    "validation"
  ))
  invisible(lapply(dirs, dir_create))
  invisible(dirs)
}

# Backward-compatible alias used by earlier PRIDE helper scripts.
pride_package_dir <- function() pride_submission_dir()

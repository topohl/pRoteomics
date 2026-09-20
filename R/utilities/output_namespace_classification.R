# Is this destination canonical, or historical?
#
# Extracted from tools/audit_writer_namespaces.R in Phase 6G.8 so the rule can
# be tested on fixtures without running the whole audit, and so every tool that
# needs the answer gets the same one.
#
# The rule is derived from the contract, not from the spelling of a path. The
# earlier version classified a destination as historical only when one of its
# literal segments matched a numbered stage namespace (^[0-9]{2}[a-z]?_), which
# is a proxy: results/reviewer_audit/ carries no stage number, so five WGCNA
# writers wrote there while their declarations already said results/wgcna/...
# and the split-brain gate still reported zero. Any historical root that is not
# numbered would have slipped through the same way.
#
# The authorities are:
#   config/output_layout.yml        the canonical shape and the three lifecycles
#   config/legacy_output_registry.csv  the historical roots, independently
#   an explicit, adjudicated exception list for destinations outside results/

if (!exists("repo_path", mode = "function")) {
  paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
  source(paths_file)
}

# Path builders whose result is a destination when a write consumes it.
ONC_PATH_BUILDERS <- c("path_results", "path_processed")

# These create directories, so calling one is itself a write.
ONC_LEGACY_FACTORIES <- c("create_module_dirs", "module_paths", "qc_paths",
                          "wgcna_downstream_paths")

# Destinations outside results/ that are adjudicated NOT historical:
#   config/           a generated configuration contract (adjudicated 6G.5)
#   exports/          the frozen outward-facing bundle
#   pride_submission/ gitignored export staging
ONC_ALLOWED_NONRESULT_ROOTS <- c("config", "exports", "pride_submission")

# A numbered stage namespace. Retained as an ADDITIONAL signal for the case
# where the first argument is computed, never as the only rule.
ONC_STAGE_NS <- "^[0-9]{2}[a-z]?_[A-Za-z]"

onc_canonical_domains <- function() {
  d <- tryCatch(output_layout_domains(), error = function(...) character(0))
  if (!length(d)) character(0) else d
}

onc_canonical_children <- function() {
  d <- tryCatch(output_layout_children(), error = function(...) character(0))
  if (!length(d)) character(0) else d
}

# The registered historical roots, as their segment under results/.
#
# Only the LEGACY_READ_ONLY rows count. The registry also carries
# ACTIVE_NOT_LEGACY rows, and after Phase 6G those include the normalized
# domain roots themselves - results/qc, results/integration and the rest,
# each with live writers. Taking every row would have made the fallback branch
# classify a canonical domain root as historical.
onc_legacy_registered_segments <- function() {
  reg <- repo_path("config", "legacy_output_registry.csv")
  if (!file.exists(reg)) return(character(0))
  d <- tryCatch(utils::read.csv(reg, stringsAsFactors = FALSE),
                error = function(...) NULL)
  if (is.null(d) || !"legacy_path" %in% names(d)) return(character(0))
  if ("policy" %in% names(d)) d <- d[d$policy == "LEGACY_READ_ONLY", , drop = FALSE]
  if (!nrow(d)) return(character(0))
  s <- vapply(strsplit(d$legacy_path, "/", fixed = TRUE),
              function(p) if (length(p) >= 2L) p[[2]] else NA_character_,
              character(1))
  s <- unique(s[!is.na(s)])
  ## results/publication_source_data is LEGACY_READ_ONLY - it is where the
  ## bundle lived before Phase 6F moved it to exports/ - but it is also a
  ## declared domain name, so it must not shadow that domain's canonical root.
  setdiff(s, onc_canonical_domains())
}

onc_is_skippable <- function(x) {
  tryCatch(is.null(x) || (is.symbol(x) && !nzchar(as.character(x))),
           error = function(...) TRUE)
}

onc_call_name <- function(e) {
  if (!is.call(e)) return(NA_character_)
  fn <- e[[1]]
  if (is.name(fn)) return(as.character(fn))
  if (is.call(fn) && length(fn) == 3L && is.name(fn[[1]]) &&
      as.character(fn[[1]]) %in% c("::", ":::")) {
    return(as.character(fn[[3]]))
  }
  NA_character_
}

# TRUE when this call constructs a historical or otherwise noncanonical
# destination.
#
# Decision order:
#   1. a legacy directory factory is historical by definition
#   2. path_processed() is noncanonical: config/output_layout.yml declares
#      exactly three lifecycles - work, results, exports - and data/processed
#      is not one of them
#   3. a literal first argument that is a declared output domain is canonical;
#      a literal first argument that is anything else is a historical kind or
#      root, which is what catches results/reviewer_audit/
#   4. when the first argument is computed, fall back to the registered
#      historical segments and the numbered-stage signal
is_legacy_path_construction <- function(e,
                                        domains = onc_canonical_domains(),
                                        legacy_segments = onc_legacy_registered_segments(),
                                        allowed = ONC_ALLOWED_NONRESULT_ROOTS) {
  nm <- onc_call_name(e)
  if (is.na(nm)) return(FALSE)
  if (nm %in% ONC_LEGACY_FACTORIES) return(TRUE)
  if (!nm %in% ONC_PATH_BUILDERS) return(FALSE)

  args <- Filter(function(a) !onc_is_skippable(a), as.list(e)[-1])
  lits <- unlist(lapply(args, function(a) if (is.character(a)) a else NULL))

  if (identical(nm, "path_processed")) return(TRUE)

  first <- if (length(args) && is.character(args[[1]]) && length(args[[1]]) == 1L) {
    args[[1]]
  } else NA_character_
  if (!is.na(first)) {
    if (first %in% domains) return(FALSE)
    if (first %in% allowed) return(FALSE)
    return(TRUE)
  }

  any(lits %in% legacy_segments) || any(grepl(ONC_STAGE_NS, lits))
}

# Convenience for tests and tools: classify a literal repo-relative path.
is_legacy_output_path <- function(path,
                                  domains = onc_canonical_domains(),
                                  legacy_segments = onc_legacy_registered_segments(),
                                  allowed = ONC_ALLOWED_NONRESULT_ROOTS) {
  seg <- strsplit(sub("^[.]/", "", path), "/", fixed = TRUE)[[1]]
  seg <- seg[nzchar(seg)]
  if (!length(seg)) return(FALSE)
  if (seg[[1]] %in% allowed) return(FALSE)
  if (identical(seg[[1]], "work")) return(FALSE)
  if (identical(seg[[1]], "data") && length(seg) >= 2L && identical(seg[[2]], "processed")) {
    return(TRUE)
  }
  if (!identical(seg[[1]], "results")) return(FALSE)
  if (length(seg) < 2L) return(FALSE)
  if (seg[[2]] %in% domains) return(FALSE)
  TRUE
}

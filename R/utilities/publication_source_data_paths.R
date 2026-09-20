# Canonical and historical locations for publication_source_data outputs.
#
# This domain sits on the publication boundary, so its layout answers two
# different questions and they must not be confused.
#
#   PRODUCTION  results/publication_source_data/<analysis_id>/<scope>/<child>/
#               where a writer canonically produces. config/output_layout.yml
#               declares exports with may_be_canonical: false, so a writer can
#               never canonically produce INTO exports/.
#
#   BOUNDARY    exports/publication_source_data/
#               the frozen outward-facing bundle. It is a COPY of canonical
#               results, carrying a sha256 per file, and it is the only thing
#               Exp9_manuscript imports. Produced by packaging, not by an
#               analysis writer.
#
# Path semantics only. No selection logic, no schema logic, no science: this
# file answers "where does this object live" and nothing else.

if (!exists("repo_path", mode = "function")) {
  paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
  source(paths_file)
}

PSD_DOMAIN <- "publication_source_data"

# Historical roots this domain wrote to before Phase 6G.9. reviewer_audit is a
# shared cross-domain tree - WGCNA and the publication hardening audits also
# live there - so only this domain's own artifacts move out of it.
PSD_LEGACY_ROOTS <- c("reviewer_audit", "tables", "logs", "source_data")

# The canonical production directories for one analysis.
psd_dirs <- function(analysis_id, scope = "global", create = FALSE) {
  if (is.null(scope) || !length(scope) || !nzchar(scope)) scope <- "global"
  canonical_module_dirs(PSD_DOMAIN, analysis_id, scope = scope, create = create)
}

# The frozen outward bundle root. Not a production target: packaging copies
# into it and the manuscript imports from it.
psd_export_root <- function(...) repo_path("exports", PSD_DOMAIN, ...)

# Normalized-first resolution for a publication_source_data artifact that
# something else reads.
#
# Structured exactly like the QC and WGCNA resolvers: the caller names the
# artifact and its owner, the normalized location is preferred, and the
# historical location remains a working fallback until a canonical run exists.
# Returns the normalized path when neither is present, so a caller's own
# missing-input handling runs unchanged.
psd_find <- function(filename, owner, legacy_root = "reviewer_audit",
                     scope = "global", child = "tables") {
  norm <- file.path(canonical_result_path(PSD_DOMAIN, owner, scope, child), filename)
  legacy <- if (identical(legacy_root, "tables")) {
    path_results("tables", filename)
  } else {
    path_results(legacy_root, filename)
  }
  for (p in c(norm, legacy)) if (file.exists(p)) return(unname(p))
  unname(norm)
}

# The claims-table reviewer audits, which several validators read back.
psd_claims_audit <- function(filename) {
  psd_find(filename, owner = "build_biological_claims_table",
           legacy_root = "reviewer_audit")
}

# The biological claims table itself. Historically it sat at the bare
# results/tables/ root with no stage segment at all.
psd_claims_table <- function(filename = "biological_claims_table.csv") {
  psd_find(filename, owner = "build_biological_claims_table",
           legacy_root = "tables")
}

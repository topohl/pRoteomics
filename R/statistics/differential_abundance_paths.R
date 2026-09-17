# Where differential_abundance's results live, before and after Phase 6G.4.
#
# Historical layout:
#   results/<kind>/04_differential_expression_enrichment/<substep>/<rest...>
#   data/processed/04_differential_expression_enrichment/<substep>/<dataset>/...
#
# Normalized layout, keyed on analysis identity:
#   results/differential_abundance/<analysis_id>/<scope>/<child>/<suffix...>
#
# The historical substep names are dropped because the analysis_id replaces
# them, and one of them carried chronology that should not survive:
# 01b_gsea_protein_direction_audit becomes audit_gsea_protein_direction. No
# scientific vocabulary is lost; ontology and dataset stay in the path because
# they carry meaning.
#
# --- what the processed tree becomes ---------------------------------------
#
# Two analyses here wrote into data/processed rather than results: the
# clusterProfiler manifest and bundle, and the compareGO input manifest. They
# are not intermediates. Downstream analyses in enrichment, integration and
# wgcna read them as state, through
# canonical_clusterprofiler_manifest_path() and
# canonical_comparego_manifest_path(), so they are models/ and not work/: the
# lifecycle contract forbids a downstream contract naming work/.
#
# --- what deliberately does not move --------------------------------------
#
# run_clusterprofiler_enrichment also supports a branch-comparison mode. With
# PROTEOMICS_ENRICHMENT_BRANCH set it writes to
# 04_differential_expression_enrichment_comparison/<branch>, a sandbox that
# exists so an experimental run cannot be mistaken for canonical output.
# R/utilities/export_helpers.R and test-source-data-journal-scope.R already
# exclude that tree from journal scope by name. It is not a canonical output,
# nothing declares or consumes it, and repointing it would mean editing those
# publication-scope rules, so it stays exactly where it is. Only the canonical
# branch is migrated.

# The canonical output directories for one differential_abundance analysis.
#
# create = FALSE by default: the caller creates what it writes, so an analysis
# that only produces tables does not leave empty plots/ and models/ behind.
# Every member create_module_dirs() offered is present except `processed`,
# which is deliberately absent so that each former user has to choose between
# models/ and work/ rather than inheriting a processed-data destination by
# accident.
differential_abundance_dirs <- function(analysis_id, scope = "global",
                                        suffix = NULL, create = FALSE) {
  canonical_module_dirs("differential_abundance", analysis_id, scope = scope,
                        suffix = suffix, create = create)
}

# The canonical shape as a repo-relative path.
#
# For the one caller that cannot use an absolute repo path:
# validate_control_spatial_identity supports PROTEOMICS_* redirection of its
# output root so the validation harness can write elsewhere, and
# canonical_result_path() always resolves against the repository root.
#
# It is derived from canonical_result_path() rather than written out
# independently, so the shape cannot drift from the contract, and the contract
# still gets to reject an unknown domain or child.
differential_abundance_relative_path <- function(analysis_id, scope = "global",
                                                 child = "tables", ...) {
  ## file.path() returns character(0) if any argument has length zero, so an
  ## absent optional segment has to be dropped rather than passed through as
  ## NULL. This silently produced empty paths the first time round.
  extra <- Filter(function(x) length(x) && nzchar(x), list(...))
  abs <- do.call(canonical_result_path,
                 c(list("differential_abundance", analysis_id, scope, child),
                   extra))
  root <- normalizePath(repo_root(), winslash = "/", mustWork = FALSE)
  sub(paste0("^", root, "/"), "",
      normalizePath(abs, winslash = "/", mustWork = FALSE))
}

# Normalized-first candidates for a processed-state artifact that downstream
# domains read. `owner` is the producing analysis_id.
differential_abundance_state_candidates <- function(filename, owner, dataset,
                                                    legacy_substep,
                                                    repository_root = repo_path()) {
  c(
    normalized = file.path(repository_root, "results", "differential_abundance",
                           owner, as.character(dataset), "models", filename),
    legacy = file.path(repository_root, "data", "processed",
                       "04_differential_expression_enrichment", legacy_substep,
                       as.character(dataset), filename)
  )
}

# Resolve one such artifact, preferring the normalized copy.
#
# Until an analysis is rerun the normalized copy does not exist and the
# historical one answers, which is what keeps every current consumer working
# unchanged. Existence is tested on the file, not on its directory, so an
# empty skeleton left by a dry run cannot shadow real data.
resolve_differential_abundance_state <- function(filename, owner, dataset,
                                                 legacy_substep,
                                                 repository_root = repo_path()) {
  cand <- differential_abundance_state_candidates(filename, owner, dataset,
                                                  legacy_substep, repository_root)
  hit <- cand[file.exists(cand)]
  if (length(hit)) return(unname(hit[1]))
  unname(cand[["normalized"]])
}

# Canonical and historical locations for WGCNA outputs.
#
# Path semantics only. No module construction, statistics or label logic
# belongs here: this file answers "where does this object live" and nothing
# else, so the hundred-odd consumer files across seven analysis domains can
# resolve a WGCNA path without loading the WGCNA stack.
#
# WGCNA is the largest domain in the repository - 24 writers, 123 declared
# outputs, 5,926 historical files - and its historical layout uses more roots
# than any other:
#
#   results/{tables,figures,source_data,logs,reports}/06_modules_WGCNA/<substep>/...
#   data/processed/06_modules_WGCNA/<substep>/...        (network state carriers)
#   results/reviewer_audit/<family>/...                  (reviewer-facing audits,
#                                                         no 06_modules_WGCNA segment)
#
# The reviewer_audit tree is shared with publication_source_data - 23 of its 41
# entries are WGCNA-owned and the rest are not - so it is offered as an explicit
# legacy family rather than derived, and normalization moves only the
# WGCNA-owned artifacts out of it.
#
# Deliberately NOT a path DSL (Phase 6G.8 section 33): one directory helper,
# one candidate builder, and a small number of named artifact helpers for the
# objects that many consumers read.

if (!exists("repo_path", mode = "function")) {
  paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
  source(paths_file)
}

WGCNA_DOMAIN <- "wgcna"
WGCNA_LEGACY_STAGE <- "06_modules_WGCNA"

# Historical families a WGCNA object could have lived in.
WGCNA_LEGACY_FAMILIES <- c("tables", "figures", "source_data", "logs", "reports",
                           "processed", "reviewer_audit")

# Normalized children, as config/output_layout.yml declares them.
WGCNA_CHILDREN <- c("tables", "plots", "models", "manifests", "reports")

# Drop zero-length arguments before file.path() sees them. file.path() returns
# character(0) if any argument has length zero, which silently turns a real
# destination into an empty path.
wg_join <- function(...) {
  parts <- Filter(function(x) length(x) && !is.na(x) && nzchar(x), list(...))
  if (!length(parts)) return(character(0))
  do.call(file.path, parts)
}

# Trailing segments, with the empty ones dropped.
#
# canonical_result_path() forwards its ... straight to path_results(), which
# reaches file.path(), and file.path() returns character(0) if ANY argument has
# length zero. So a caller passing a NULL trailing segment - which is the
# natural way to say "no extra segment here" - silently collapsed the whole
# destination to an empty path. wg_join() already guards the legacy side; this
# guards the normalized side.
wg_extra <- function(...) {
  Filter(function(x) length(x) && !is.na(x) && nzchar(x), list(...))
}

wg_canonical <- function(owner, scope, child, ...) {
  do.call(canonical_result_path,
          c(list(WGCNA_DOMAIN, owner, scope, child), wg_extra(...)))
}

# Does this directory actually hold a contract file?
#
# Regular files only, and defined exactly once. dir.exists() is not evidence
# that data exist: a dry run, or dirs(create = TRUE), leaves a skeleton behind,
# and one of those skeletons contains a tables/source_data subdirectory, so
# counting directory entries made an empty normalized directory look populated
# and shadow real historical data.
wg_has_files <- function(d) {
  if (!length(d) || is.na(d) || !nzchar(d) || !dir.exists(d)) return(FALSE)
  e <- list.files(d, full.names = TRUE, no.. = TRUE)
  if (!length(e)) return(FALSE)
  if (any(!file.info(e)$isdir)) return(TRUE)
  ## Only sub-directories at the top level. That is not empty: a Stage-01
  ## dataset directory holds modules/ and supermodules/ and no loose files, and
  ## a top-level-only test called it empty and fell through to the normalized
  ## default - selecting a directory that does not exist over the populated
  ## historical one. The recursive check runs only in this case, so the common
  ## path stays a single listing.
  ##
  ## The 6G.5 guard is preserved: a normalized tree of freshly created but
  ## empty directories still contains no file anywhere and is still rejected.
  length(list.files(d, recursive = TRUE, include.dirs = FALSE)) > 0L
}

# Where WGCNA results live after Phase 6G.8:
#   results/wgcna/<analysis_id>/<scope>/<child>/
#
# This replaces create_module_dirs("06_modules_WGCNA", "<substep>/<dataset>")
# and module_paths(...) for every writer. create = TRUE mirrors
# create_module_dirs' side effect; the default is FALSE so a dry run can build
# a destination without leaving a skeleton behind.
wgcna_dirs <- function(analysis_id, scope = "global", suffix = NULL,
                       create = FALSE) {
  if (is.null(scope) || !length(scope) || !nzchar(scope)) scope <- "global"
  canonical_module_dirs(WGCNA_DOMAIN, analysis_id, scope = scope,
                        suffix = suffix, create = create)
}

# The historical root for one legacy family and substep.
wg_legacy_root <- function(legacy_family, legacy_substep) {
  legacy_family <- match.arg(legacy_family, WGCNA_LEGACY_FAMILIES)
  if (identical(legacy_family, "processed")) {
    wg_join(path_processed(WGCNA_LEGACY_STAGE), legacy_substep)
  } else if (identical(legacy_family, "reviewer_audit")) {
    ## no 06_modules_WGCNA segment: the reviewer tree is flat by family name
    wg_join(path_results("reviewer_audit"), legacy_substep)
  } else {
    wg_join(path_results(legacy_family, WGCNA_LEGACY_STAGE), legacy_substep)
  }
}

# Normalized-first candidates for one named WGCNA artifact.
#
# `owner` is the WGCNA analysis_id that declares the file and `legacy_substep`
# the directory it used to live in. Both are checked against the registry by
# tests/testthat/test-wgcna-writer-namespace.R, so a wrong owner is a test
# failure rather than a silent mis-resolution.
#
# Trailing `...` segments are preserved in both candidates, because
# config/output_layout.yml keeps meaningful post-scope segments: a module
# family, a spatial unit or a contrast direction is part of the object's
# identity, not chronology.
#
# Existence is tested on the file, never on its directory.
wgcna_artifact_candidates <- function(filename, owner, legacy_substep,
                                      scope = "global", child = "tables",
                                      legacy_family = "tables",
                                      legacy_scoped = TRUE, ...,
                                      source_data = FALSE) {
  child <- match.arg(child, WGCNA_CHILDREN)
  legacy_family <- match.arg(legacy_family, WGCNA_LEGACY_FAMILIES)
  if (is.null(scope) || !length(scope) || !nzchar(scope)) scope <- "global"

  ## source_data inverts the usual child/family relationship. Historically it
  ## was its own root, results/source_data/06_modules_WGCNA/<substep>/<scope>/,
  ## with no extra segment; normalized it is a child OF tables,
  ## results/wgcna/<owner>/<scope>/tables/source_data/. So the extra segment
  ## belongs on the normalized side only, and passing it through ... would
  ## wrongly add it to the historical path as well.
  norm <- if (isTRUE(source_data)) {
    wg_canonical(owner, scope, child, "source_data", ...) |> wg_join(filename)
  } else {
    wg_canonical(owner, scope, child, ...) |> wg_join(filename)
  }

  root <- wg_legacy_root(legacy_family, legacy_substep)
  ## The historical layout put the dataset in different places: most substeps
  ## append it, a few wrote straight into the substep directory, so both shapes
  ## are offered rather than assumed.
  legacy_a <- if (isTRUE(legacy_scoped)) {
    wg_join(root, scope, ..., filename)
  } else {
    wg_join(root, ..., filename)
  }
  legacy_b <- if (isTRUE(legacy_scoped)) {
    wg_join(root, ..., filename)
  } else {
    wg_join(root, scope, ..., filename)
  }

  c(normalized = unname(norm), legacy = unname(legacy_a),
    legacy_alt = unname(legacy_b))
}

# Resolve one named WGCNA artifact, normalized first.
#
#   normalized present + historical present -> normalized
#   normalized absent  + historical present -> historical
#   normalized present + historical absent  -> normalized
#   neither                                 -> the normalized path, so the
#                                              caller's own missing-input
#                                              message names the future home
wgcna_find <- function(filename, owner, legacy_substep, scope = "global",
                       child = "tables", legacy_family = "tables",
                       legacy_scoped = TRUE, ..., source_data = FALSE) {
  cand <- wgcna_artifact_candidates(filename, owner, legacy_substep, scope,
                                    child, legacy_family, legacy_scoped, ...,
                                    source_data = source_data)
  hit <- cand[file.exists(cand)]
  if (length(hit)) return(unname(hit[1]))
  unname(cand[["normalized"]])
}

# One shared internal for every family helper, so the child-to-legacy-family
# mapping and the source_data inversion are expressed once.
wg_family_find <- function(filename, owner, substep, scope, child = "tables",
                           source_data = FALSE, legacy_scoped = TRUE, ...) {
  if (isTRUE(source_data)) {
    return(wgcna_find(filename, owner = owner, legacy_substep = substep,
                      scope = scope, child = "tables",
                      legacy_family = "source_data",
                      legacy_scoped = legacy_scoped, ...,
                      source_data = TRUE))
  }
  wgcna_find(filename, owner = owner, legacy_substep = substep, scope = scope,
             child = child, legacy_family = wg_legacy_family_for(child),
             legacy_scoped = legacy_scoped, ...)
}

# The four writer-own substeps that carry a handful of cross-domain reads.
# They get named helpers too, so every converted call site names an owner.
wgcna_claim_readiness_artifact <- function(filename, scope = "microglia",
                                           child = "tables", source_data = FALSE, ...) {
  wg_family_find(filename, "audit_module_claim_readiness", "claim_readiness",
                 scope, child, source_data, TRUE, ...)
}
wgcna_complex_architecture_artifact <- function(filename, dataset,
                                                child = "tables", source_data = FALSE, ...) {
  wg_family_find(filename, "summarize_module_complex_architecture",
                 "module_complex_architecture", dataset, child, source_data, TRUE, ...)
}
wgcna_robustness_artifact <- function(filename, dataset, child = "tables",
                                      source_data = FALSE, ...) {
  wg_family_find(filename, "audit_module_robustness",
                 "module_robustness_sensitivity", dataset, child, source_data, TRUE, ...)
}
wgcna_score_summary_artifact <- function(filename, dataset, child = "tables",
                                         source_data = FALSE, ...) {
  wg_family_find(filename, "summarize_module_scores",
                 "score_publication_summary", dataset, child, source_data, TRUE, ...)
}

# The identity-contract directory, for a caller that then names the entity and
# membership contracts under it. The scope is always explicit: the per-dataset
# contracts have byte-identical filenames and their module ids collide, so a
# crossed scope would join cleanly and silently.
wgcna_identity_contract_dir <- function(dataset, child = "tables") {
  wgcna_dir_any(owner = "build_module_identity_contract",
                legacy_substep = "identity_contract", scope = dataset,
                child = child, legacy_family = wg_legacy_family_for(child),
                legacy_scoped = TRUE)
}

# The microglia-neuropil independence directory, owned by
# test_microglia_neuropil_independence. Its figure renderer is a different
# analysis and reads this directory rather than writing it.
wgcna_independence_dir <- function(scope = "microglia", child = "tables") {
  wgcna_dir_any(owner = "test_microglia_neuropil_independence",
                legacy_substep = "microglia_neuropil_independence", scope = scope,
                child = child, legacy_family = wg_legacy_family_for(child),
                legacy_scoped = TRUE)
}

# The module-score directory, for a caller that then names several artifacts
# under it. The module-definition source is required for the same reason it is
# required on the artifact helper: neuron_neuropil has both overlap/ and wgcna/
# on disk and they are different objects.
wgcna_module_score_dir <- function(dataset, module_definition_source,
                                   child = "tables") {
  if (!length(module_definition_source) || !nzchar(module_definition_source)) {
    stop("module_definition_source is required: neuron_neuropil has more than ",
         "one module-score source on disk and they are different objects.",
         call. = FALSE)
  }
  wgcna_dir_any(owner = "score_module_activity", legacy_substep = "module_score",
                scope = dataset, child = child,
                legacy_family = wg_legacy_family_for(child),
                legacy_scoped = TRUE, module_definition_source)
}

# A directory, normalized first, for a reader that wants the whole folder
# rather than one named file. Preferred only when it holds regular files, so an
# empty normalized skeleton cannot shadow populated historical data. Anchoring
# on a single guessed probe filename is not used here at all (section 28).
wgcna_dir_any <- function(owner, legacy_substep, scope = "global",
                          child = "tables", legacy_family = "tables",
                          legacy_scoped = TRUE, ...) {
  child <- match.arg(child, WGCNA_CHILDREN)
  legacy_family <- match.arg(legacy_family, WGCNA_LEGACY_FAMILIES)
  if (is.null(scope) || !length(scope) || !nzchar(scope)) scope <- "global"

  norm <- wg_canonical(owner, scope, child, ...)
  if (wg_has_files(norm)) return(unname(norm))

  root <- wg_legacy_root(legacy_family, legacy_substep)
  ordered <- if (isTRUE(legacy_scoped)) {
    list(wg_join(root, scope, ...), wg_join(root, ...))
  } else {
    list(wg_join(root, ...), wg_join(root, scope, ...))
  }
  for (d in ordered) if (wg_has_files(d)) return(unname(d))
  unname(norm)
}

# --------------------------------------------------- payload-aware root choice
#
# A broad root is not an artifact. Choosing between a normalized and a
# historical root by dir.exists(), or by "contains any file", is unsafe for two
# separate reasons this migration has already hit:
#
#   * a normalized root can exist and be empty, or hold only freshly created
#     sub-directories, and would then shadow the populated historical root
#     (Batch 3C, wg_has_files);
#   * a historical WGCNA root contains microglia_failed_20260720_133211, so
#     "contains files" can be satisfied entirely by a failed run.
#
# A root therefore counts only when it contains the payload the CALLER needs.
# The caller states that payload explicitly - filenames or regexes - so the
# contract lives at the call site rather than being guessed here.
#
# `required` is a character vector of patterns matched recursively, or a
# predicate taking the directory.
wg_has_payload <- function(d, required) {
  if (!length(d) || is.na(d) || !nzchar(d) || !dir.exists(d)) return(FALSE)
  if (is.function(required)) return(isTRUE(required(d)))
  if (!length(required)) return(FALSE)
  any(vapply(required, function(p) {
    length(list.files(d, pattern = p, recursive = TRUE, include.dirs = FALSE)) > 0L
  }, logical(1)))
}

# Normalized first, historical second, and neither unless it actually carries
# the payload. When neither does, the normalized root is returned so the
# caller's own missing-input handling runs unchanged - every caller of this
# already guards with dir.exists()/file.exists().
wg_root_by_payload <- function(normalized, legacy, required) {
  if (wg_has_payload(normalized, required)) return(unname(normalized))
  if (wg_has_payload(legacy, required)) return(unname(legacy))
  unname(normalized)
}

# curated_overlap_programs - build_curated_overlap_programs, global scope.
# Used as a last-resort newest-match search root, so the payload is named: only
# a directory actually holding one of the two program workbooks qualifies.
WGCNA_CURATED_OVERLAP_PAYLOAD <- c("^curated_overlap_programs[.]xlsx$",
                                   "^Overlap_based_neuropil_modules_classified[.]xlsx$")

wgcna_curated_overlap_dir <- function(child = "tables",
                                      required = WGCNA_CURATED_OVERLAP_PAYLOAD) {
  wg_root_by_payload(
    wg_canonical("build_curated_overlap_programs", "global", child),
    wg_join(wg_legacy_root(wg_legacy_family_for(child), "curated_overlap_programs"), "global"),
    required)
}

# wgcna_publication_figures - render_module_figures. The scope is explicit so
# that no other scope, and in particular no failed run, can be pulled in by a
# recursive listing of the family root.
wgcna_publication_figures_dir <- function(scope, child = "plots",
                                          required = "[.](svg|pdf|png)$") {
  wg_root_by_payload(
    wg_canonical("render_module_figures", scope, child),
    wg_join(wg_legacy_root(wg_legacy_family_for(child), "wgcna_publication_figures"), scope),
    required)
}

# Every root a publication figure collector must look in.
#
# This is a root INVENTORY, not a fallback choice: the collector has to find
# figures wherever they currently are, and they are still historical. Both
# roots are returned, normalized first. Dataset-scope filtering is the
# caller's responsibility and is NOT done here - see
# drop_noncanonical_wgcna_dataset_scopes(), which is what keeps the failed run
# out of publication selection.
wgcna_figure_discovery_roots <- function() {
  unname(c(repo_path("results", WGCNA_DOMAIN),
           wg_legacy_root("figures", character(0))))
}

# The historical family that corresponds to a normalized child.
#
# Historically an artifact's root encoded its kind: a table went to
# results/tables/, a figure to results/figures/, a run manifest to
# results/logs/, the serialised state to data/processed/. Normalizing renames
# those kinds to children, so a family helper that hard-codes one legacy family
# resolves the wrong historical path as soon as a caller asks for a different
# child - which is exactly what happened for the two figure families, whose
# validation CSV lives under results/tables/ rather than results/figures/.
wg_legacy_family_for <- function(child) {
  switch(child,
         tables = "tables",
         plots = "figures",
         models = "processed",
         manifests = "logs",
         reports = "reports",
         "tables")
}

# ------------------------------------------------ owner-specific family API
#
# One helper per historical artifact family, eleven in all, each binding the
# canonical producer, the historical substep and the legacy family so a caller
# cannot get them wrong. The caller supplies the exact artifact name and the
# scope, which keeps object identity explicit at the call site.
#
# This is deliberately NOT a generic resolver. There is no
# find_any_wgcna_file(), no first_existing_wgcna_object() and no
# scan_all_results_wgcna(), because the family investigation found eight
# concrete ways a generic one returns the wrong object: an ontology token
# (_BP/_CC/_MF) distinguishing three siblings in one flat directory; a
# skip-stub written to the same filename as a real result on five failure
# branches; a corrected-versus-uncorrected figure family owned by a different
# analysis_id; module ids that collide across datasets and join cleanly; a
# module_definition_source segment that neuron_neuropil has in two variants on
# disk; caution labels that must never cross from microglia to the neuron
# datasets; and a retained failed-run directory sitting alongside the three
# real dataset scopes.
#
# Each helper is a thin binding over wgcna_find()/wgcna_dir_any(). They hold
# path semantics, owner identity, scope resolution and the compatibility
# fallback, and no scientific logic.

# 01_WGCNA - build_wgcna_modules. Module construction: membership, supermodule
# assignment, preservation, feature universe, the serialised network state.
# `...` carries the meaningful post-scope segment, normally "modules" or
# "supermodules".
# Stage 01 keeps its objects in sub-directories - modules/, supermodules/,
# inputs/ - so the caller names the sub-directory through ... after the child.
# Two of its objects are figure source data, and for those the normalized and
# historical sides invert (normalized tables/source_data vs a historical
# source_data root of its own), which is what wg_family_find expresses.
# `...` precedes source_data deliberately, as in wgcna_artifact_candidates().
# This family is the only one that routinely passes a sub-directory, and with
# source_data ahead of `...` a positional "modules" binds to source_data
# instead: the sub-directory vanishes from the path and the resolver silently
# looks for the file one level too high.
wgcna_modules_artifact <- function(filename, dataset, child = "tables", ...,
                                   source_data = FALSE) {
  wg_family_find(filename, "build_wgcna_modules", "01_WGCNA", dataset, child,
                 source_data, TRUE, ...)
}

# group_effects - test_module_phenotypes. Stage-05 non-inferential values and
# the phenotype test outputs. Hemisphere-resolved values live here under their
# own filenames; the caller names which one it needs, so bilateral and
# hemisphere forms cannot be substituted for one another.
wgcna_group_effects_artifact <- function(filename, dataset, child = "tables", ...) {
  wgcna_find(filename, owner = "test_module_phenotypes",
             legacy_substep = "group_effects", scope = dataset, child = child,
             legacy_family = wg_legacy_family_for(child), legacy_scoped = TRUE, ...)
}

# interpretable_summary - summarize_module_interpretation. Adjudicated labels
# and the claim-facing inferential handoff.
# Two of its objects are figure source data rather than tables, and the
# normalized and historical sides invert for those: normalized keeps them under
# tables/source_data while the historical tree had a source_data family of its
# own. wg_family_find expresses that inversion once, so this helper routes
# through it exactly as the other source_data-carrying families do.
wgcna_interpretable_artifact <- function(filename, dataset, child = "tables",
                                         source_data = FALSE, ...) {
  wg_family_find(filename, "summarize_module_interpretation",
                 "interpretable_summary", dataset, child, source_data, TRUE, ...)
}

# module_annotation - annotate_module_microenvironment. Biological annotation
# and the microenvironment caution labels, which are dataset-specific by
# contract: a neuron dataset's annotation must never carry microglia caution
# text, so the scope is always explicit here.
wgcna_annotation_artifact <- function(filename, dataset, child = "tables", ...) {
  wgcna_find(filename, owner = "annotate_module_microenvironment",
             legacy_substep = "module_annotation", scope = dataset,
             child = child, legacy_family = wg_legacy_family_for(child),
             legacy_scoped = TRUE, ...)
}

# identity_contract - build_module_identity_contract. The per-dataset entity
# and membership contracts. Their filenames are byte-identical across datasets
# and module ids collide, so a crossed scope would join cleanly and silently.
wgcna_identity_contract_artifact <- function(filename, dataset, child = "tables", ...) {
  wgcna_find(filename, owner = "build_module_identity_contract",
             legacy_substep = "identity_contract", scope = dataset,
             child = child, legacy_family = wg_legacy_family_for(child),
             legacy_scoped = TRUE, ...)
}

# microglia_neuropil_independence - test_microglia_neuropil_independence.
# The dataset-scoped trio, plus three global-scope reviewer audits that lived
# flat under results/reviewer_audit/.
wgcna_independence_artifact <- function(filename, scope = "microglia",
                                        child = "tables", ...) {
  wgcna_find(filename, owner = "test_microglia_neuropil_independence",
             legacy_substep = "microglia_neuropil_independence", scope = scope,
             child = child, legacy_family = wg_legacy_family_for(child),
             legacy_scoped = TRUE, ...)
}
wgcna_independence_audit <- function(filename) {
  ## the flat reviewer-audit form: results/reviewer_audit/<file>.csv
  wgcna_find(filename, owner = "test_microglia_neuropil_independence",
             legacy_substep = NULL, scope = "global", child = "tables",
             legacy_family = "reviewer_audit", legacy_scoped = FALSE)
}

# curated_overlap_programs - build_curated_overlap_programs. Global scope only.
# Every sibling resolver in score_module_activity passes a dataset, and this
# one must not: a dataset scope here names a location no producer writes.
wgcna_curated_overlap_artifact <- function(filename, child = "tables", ...) {
  wgcna_find(filename, owner = "build_curated_overlap_programs",
             legacy_substep = "curated_overlap_programs", scope = "global",
             child = child,
             legacy_family = wg_legacy_family_for(child),
             legacy_scoped = TRUE, ...)
}

# 01b_module_supermodule_GO_heatmaps - render_module_go_heatmaps. The three
# ontology variants (_BP/_CC/_MF) share one flat directory and differ only by
# that token, so the caller always names the full filename.
wgcna_go_heatmap_artifact <- function(filename, dataset, child = "tables", ...) {
  wgcna_find(filename, owner = "render_module_go_heatmaps",
             legacy_substep = "01b_module_supermodule_GO_heatmaps",
             scope = dataset, child = child, legacy_family = wg_legacy_family_for(child),
             legacy_scoped = TRUE, ...)
}

# module_score - score_module_activity. The module-definition source (wgcna,
# overlap, custom) is a required segment, not an optional one: neuron_neuropil
# has both overlap/ and wgcna/ on disk, so dropping it would silently return
# overlap-derived scores where wgcna-derived were needed.
wgcna_module_score_artifact <- function(filename, dataset,
                                        module_definition_source,
                                        child = "tables", ...) {
  if (!length(module_definition_source) || !nzchar(module_definition_source)) {
    stop("module_definition_source is required: neuron_neuropil has more than ",
         "one module-score source on disk and they are different objects.",
         call. = FALSE)
  }
  wgcna_find(filename, owner = "score_module_activity",
             legacy_substep = "module_score", scope = dataset, child = child,
             legacy_family = wg_legacy_family_for(child), legacy_scoped = TRUE,
             module_definition_source, ...)
}

# 04_wgcna_de_gsea_overlap - compare_module_enrichment_overlap.
#
# The producer writes a one-row skip stub to this same filename on five
# failure branches, so a normalized candidate can exist and still be a stub.
# Existence alone is therefore not sufficient evidence here, and the caller is
# expected to validate the contents it requires; the helper's job is only to
# name the right object for the right dataset.
wgcna_gsea_overlap_artifact <- function(filename, dataset, child = "tables", ...) {
  wgcna_find(filename, owner = "compare_module_enrichment_overlap",
             legacy_substep = "04_wgcna_de_gsea_overlap", scope = dataset,
             child = child, legacy_family = wg_legacy_family_for(child),
             legacy_scoped = TRUE, ...)
}

# wgcna_publication_figures - render_module_figures.
#
# Distinct from wgcna_publication_figures_corrected, which is a different
# scientific object owned by render_microglia_module_figures. The two historical
# directories share a prefix, so any startsWith or glob over the shorter name
# would collide; this helper binds the exact substep.
wgcna_publication_figure_artifact <- function(filename, scope,
                                              child = "plots", ...) {
  wgcna_find(filename, owner = "render_module_figures",
             legacy_substep = "wgcna_publication_figures", scope = scope,
             child = child, legacy_family = wg_legacy_family_for(child),
             legacy_scoped = TRUE, ...)
}
wgcna_corrected_figure_artifact <- function(filename, scope,
                                            child = "plots", ...) {
  wgcna_find(filename, owner = "render_microglia_module_figures",
             legacy_substep = "wgcna_publication_figures_corrected",
             scope = scope, child = child,
             legacy_family = wg_legacy_family_for(child),
             legacy_scoped = TRUE, ...)
}

# The nature-readiness reviewer family, owned by audit_microglia_module_claims
# and read by two other WGCNA scripts.
wgcna_readiness_artifact <- function(filename, child = "tables") {
  wgcna_find(filename, owner = "audit_microglia_module_claims",
             legacy_substep = "microglia_wgcna_nature_readiness",
             scope = "global", child = child, legacy_family = "reviewer_audit",
             legacy_scoped = FALSE)
}

# ------------------------------------------------- named objects readers want
#
# A handful of WGCNA objects carry most of the downstream dependency graph.
# Their readers build the path themselves rather than going through one helper,
# so these give them a chokepoint - one line instead of six resolver arguments
# repeated across seven domains. Everything else uses wgcna_find() directly.

# The serialised network state: modules, eigengenes, TOM-derived structure.
# A frozen scientific state carrier, so MODEL/PERSISTENT_OBJECT, never a cache.
wgcna_final_state <- function(dataset) {
  wgcna_find("wgcna_final_model_state.rds", owner = "build_wgcna_modules",
             legacy_substep = "01_WGCNA", scope = dataset, child = "models",
             legacy_family = "processed")
}

# Module membership, long form: the canonical module assignment table.
wgcna_modules_long <- function(dataset) {
  wgcna_find("WGCNA_modules_long.csv", owner = "build_wgcna_modules",
             legacy_substep = "01_WGCNA", scope = dataset, child = "tables",
             legacy_family = "tables", legacy_scoped = TRUE, "modules")
}

# The downstream module definition handoff.
wgcna_module_definitions <- function(dataset) {
  wgcna_find("WGCNA_module_definitions_for_downstream.csv",
             owner = "build_wgcna_modules", legacy_substep = "01_WGCNA",
             scope = dataset, child = "tables", legacy_family = "tables",
             legacy_scoped = TRUE, "modules")
}

# The adjudicated Stage-07 label lookup. Phase 6F adjudicated "final label" as
# scientifically meaningful - the adjudicated label as opposed to the raw
# Stage-01 label - so the filename is preserved verbatim (section 8/43).
wgcna_final_label_lookup <- function(dataset) {
  wgcna_find("WGCNA_final_label_lookup.csv",
             owner = "summarize_module_interpretation",
             legacy_substep = "interpretable_summary", scope = dataset,
             child = "tables", legacy_family = "tables")
}

# The module directory, for readers that glob it rather than naming one file.
wgcna_modules_dir <- function(dataset) {
  wgcna_dir_any(owner = "build_wgcna_modules", legacy_substep = "01_WGCNA",
                scope = dataset, child = "tables", legacy_family = "tables",
                legacy_scoped = TRUE, "modules")
}

# The Stage-01 supermodule directory. Kept distinct from the module directory:
# a module and a supermodule are different entities and their ids collide, so
# a caller must never be able to reach one by asking for the other.
wgcna_supermodules_dir <- function(dataset) {
  wgcna_dir_any(owner = "build_wgcna_modules", legacy_substep = "01_WGCNA",
                scope = dataset, child = "tables", legacy_family = "tables",
                legacy_scoped = TRUE, "supermodules")
}

# The Stage-01 per-dataset directory for one output kind. child names which
# kind: tables, manifests (historically logs) or models (historically
# data/processed).
wgcna_modules_scope_dir <- function(dataset, child = "tables") {
  wgcna_dir_any(owner = "build_wgcna_modules", legacy_substep = "01_WGCNA",
                scope = dataset, child = child,
                legacy_family = wg_legacy_family_for(child), legacy_scoped = TRUE)
}

# The directory whose CHILDREN are dataset scopes, for callers that enumerate
# which datasets have Stage-01 output.
#
# This is not an artifact read and not a broad results root: it is one owner's
# scope level, which sits at a different depth on each side -
# results/wgcna/build_wgcna_modules/<scope>/... normalized, and
# results/tables/06_modules_WGCNA/01_WGCNA/<scope>/... historically.
#
# Callers must still intersect the listing with valid_datasets(). The failed run
# microglia_failed_20260720_133211 is a real directory at this level, and an
# exact set intersection is the only thing that excludes it - "microglia" is a
# prefix of it, so no prefix or substring test is safe.
wgcna_modules_scope_root <- function(child = "tables") {
  child <- match.arg(child, WGCNA_CHILDREN)
  norm <- repo_path("results", WGCNA_DOMAIN, "build_wgcna_modules")
  if (dir.exists(norm) && length(list.dirs(norm, recursive = FALSE, full.names = TRUE))) {
    return(unname(norm))
  }
  unname(wg_legacy_root(wg_legacy_family_for(child), "01_WGCNA"))
}

# The Stage-05 non-inferential hemisphere values directory.
wgcna_group_effects_dir <- function(dataset) {
  wgcna_dir_any(owner = "test_module_phenotypes", legacy_substep = "group_effects",
                scope = dataset, child = "tables", legacy_family = "tables")
}

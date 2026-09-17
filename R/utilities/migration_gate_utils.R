# Gates for a writer migration.
#
# Extracted from tools/preflight_writer_migration.R so the rules can be tested
# against synthetic cases. A gate that has never been shown to fail is not a
# gate: spatial_networks has no genuine collision, so without a fixture the
# collision rule would enter the first coordinating domain unexercised.

# Two canonical writers resolving to the same destination is the failure the
# Phase 6E ownership registry exists to prevent.
#
# It is permitted in exactly one shape: the ownership registry names one of
# them the canonical owner of that result family, and the other writes
# somewhere else. "Somewhere else" is the operative part - a contributor that
# writes the *same* destination is a collision no matter what the registry
# says, because two writers would race for one file.
#
# `proposals` is a named list: analysis_id -> character vector of destinations.
# `ownership` is config/results_ownership.csv, or NULL when unavailable.
migration_destination_collisions <- function(proposals, ownership = NULL) {
  proposals <- Filter(function(p) length(p) > 0L, lapply(proposals, function(p) {
    p <- trimws(as.character(p))
    unique(p[nzchar(p)])
  }))
  if (!length(proposals)) return(list())

  flat <- unlist(proposals, use.names = FALSE)
  shared <- unique(flat[duplicated(flat)])
  if (!length(shared)) return(list())

  family_of <- function(dest) sub("/[^/]+/?$", "", dest)

  out <- list()
  for (d in shared) {
    writers <- names(proposals)[vapply(proposals, function(p) d %in% p, logical(1))]
    if (length(writers) < 2L) next

    declared_owner <- character(0)
    if (!is.null(ownership) && all(c("result_family", "canonical_owner") %in% names(ownership))) {
      fam <- family_of(d)
      hit <- ownership$canonical_owner[ownership$result_family == fam |
                                       ownership$result_family == d]
      declared_owner <- unique(hit[nzchar(hit)])
    }

    ## The ownership registry names the canonical owner of a result *family*.
    ## It never licenses two writers to target one destination: that is a race
    ## for a single file regardless of who owns the family. So a shared
    ## destination is always unresolvable, and the permitted shape from the
    ## brief - a declared owner plus a contributor writing a distinct
    ## intermediate - produces no shared destination and therefore never
    ## reaches this branch at all.
    owner_ids <- sub("[.][Rr]$", "", basename(declared_owner))

    out[[d]] <- list(
      destination = d,
      writers = writers,
      declared_owner = declared_owner,
      owner_is_one_of_the_writers = any(owner_ids %in% writers),
      resolvable = FALSE)
  }
  out
}

# A domain is blocked when a collision has no resolution, when a dependency is
# unclassified, or when an inventory could not be produced.
migration_gate_blockers <- function(collisions = list(),
                                    unknown_dependency_kinds = 0L,
                                    inventories_failed = 0L) {
  blockers <- character(0)
  if (isTRUE(unknown_dependency_kinds > 0L)) {
    blockers <- c(blockers, "unknown dependency kinds in a consumer inventory")
  }
  if (isTRUE(inventories_failed > 0L)) {
    blockers <- c(blockers, "a writer's consumer inventory could not be produced")
  }
  unresolved <- Filter(function(x) !isTRUE(x$resolvable), collisions)
  if (length(unresolved)) {
    blockers <- c(blockers,
                  "two canonical writers share a destination with no declared owner")
  }
  blockers
}

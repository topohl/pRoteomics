testthat::local_edition(3)

source(testthat::test_path(
  "..", "..", "R", "gsea_wgcna_concordance_utils.R"
))

# Minimal canonical WGCNA universe in the shape gww_build_universe_bundle()
# requires: dataset-scoped ProteinGroupIDs plus historical accession fields.
universe_fixture <- function(dataset = "neuron_soma", pgids = NULL) {
  if (is.null(pgids)) {
    pgids <- sprintf(
      "PG:%s:%s", dataset, c("AAAA1111", "BBBB2222", "CCCC3333")
    )
  }
  n <- length(pgids)
  data.frame(
    ProteinGroupID = pgids,
    included_in_wgcna = TRUE,
    ModuleID = c("WGCNA_m01", "WGCNA_m01", "WGCNA_m02")[seq_len(n)],
    member_accessions = c("P11111", "P22222", "P33333")[seq_len(n)],
    RepresentativeUniProt = c("P11111", "P22222", "P33333")[seq_len(n)],
    MemberUniProts = c("P11111", "P22222", "P33333")[seq_len(n)],
    representative_accession = c("P11111", "P22222", "P33333")[seq_len(n)],
    stringsAsFactors = FALSE
  )
}

evidence_fixture <- function(dataset, leading_edge, token_type, id = "EV1") {
  data.frame(
    gsea_evidence_id = id,
    dataset = dataset,
    leading_edge_accessions = leading_edge,
    leading_edge_token_type = token_type,
    stringsAsFactors = FALSE
  )
}

testthat::test_that("source contract routes to the right identity family", {
  testthat::expect_identical(
    gww_leading_edge_token_type("canonical_compareGO_ranked_GSEA_GO_BP"),
    "protein_group_id"
  )
  testthat::expect_identical(
    gww_leading_edge_token_type("ranked_GSEA"), "accession"
  )
  testthat::expect_error(
    gww_leading_edge_token_type("ORA"),
    "Cannot route leading-edge identity"
  )
})

testthat::test_that("the token map is typed and carries explicit provenance", {
  b <- gww_build_universe_bundle(universe_fixture(), "neuron_soma")
  testthat::expect_true(all(
    c("dataset", "ProteinGroupID", "token", "token_type",
      "token_source_column") %in% names(b$token_map)
  ))
  testthat::expect_setequal(
    unique(b$token_map$token_type), c("accession", "protein_group_id")
  )
  pg <- b$token_map[b$token_map$token_type == "protein_group_id", ]
  testthat::expect_identical(pg$token, pg$ProteinGroupID)
  testthat::expect_identical(unique(pg$token_source_column), "ProteinGroupID")
})

testthat::test_that("1. canonical PG tokens map exactly to same-dataset ProteinGroupID", {
  b <- gww_build_universe_bundle(universe_fixture(), "neuron_soma")
  ev <- evidence_fixture(
    "neuron_soma",
    "PG:neuron_soma:AAAA1111;PG:neuron_soma:BBBB2222",
    "protein_group_id"
  )
  out <- gww_map_leading_edges(ev, b$token_map)
  testthat::expect_identical(out$n_leading_edge_source, 2L)
  testthat::expect_identical(out$n_leading_edge_mapped, 2L)
  testthat::expect_identical(out$leading_edge_mapping_fraction, 1)
  testthat::expect_identical(
    out$leading_edge_proteins,
    "PG:neuron_soma:AAAA1111;PG:neuron_soma:BBBB2222"
  )
  testthat::expect_identical(out$unmapped_leading_edge_tokens, "")
})

testthat::test_that("2. a foreign-dataset ProteinGroupID is rejected, never matched", {
  b <- gww_build_universe_bundle(universe_fixture("neuron_soma"), "neuron_soma")
  ev <- evidence_fixture(
    "neuron_soma",
    "PG:neuron_soma:AAAA1111;PG:microglia:AAAA1111",
    "protein_group_id"
  )
  testthat::expect_error(
    gww_map_leading_edges(ev, b$token_map),
    "outside the PG:neuron_soma: namespace"
  )
})

testthat::test_that("2b. the same suffix in two datasets never cross-matches", {
  soma <- gww_build_universe_bundle(
    universe_fixture("neuron_soma"), "neuron_soma"
  )
  micro <- gww_build_universe_bundle(
    universe_fixture("microglia"), "microglia"
  )
  pooled <- rbind(soma$token_map, micro$token_map)
  ev <- evidence_fixture(
    "microglia", "PG:microglia:AAAA1111", "protein_group_id"
  )
  out <- gww_map_leading_edges(ev, pooled)
  testthat::expect_identical(out$n_leading_edge_mapped, 1L)
  testthat::expect_identical(out$leading_edge_proteins, "PG:microglia:AAAA1111")
  # The soma universe holds the identical hex suffix; it must be unreachable.
  testthat::expect_false(grepl("neuron_soma", out$leading_edge_proteins))
})

testthat::test_that("3. an unknown PG ID stays explicitly unmapped, never aliased", {
  b <- gww_build_universe_bundle(universe_fixture(), "neuron_soma")
  ev <- evidence_fixture(
    "neuron_soma",
    "PG:neuron_soma:AAAA1111;PG:neuron_soma:ZZZZ9999",
    "protein_group_id"
  )
  out <- gww_map_leading_edges(ev, b$token_map)
  testthat::expect_identical(out$n_leading_edge_source, 2L)
  testthat::expect_identical(out$n_leading_edge_mapped, 1L)
  testthat::expect_identical(out$leading_edge_mapping_fraction, 0.5)
  testthat::expect_identical(
    out$unmapped_leading_edge_tokens, "PG:neuron_soma:ZZZZ9999"
  )
  testthat::expect_identical(
    out$leading_edge_proteins, "PG:neuron_soma:AAAA1111"
  )
})

testthat::test_that("4. historical accession input still uses the accession path", {
  b <- gww_build_universe_bundle(universe_fixture(), "neuron_soma")
  ev <- evidence_fixture("neuron_soma", "P11111/P22222", "accession")
  out <- gww_map_leading_edges(ev, b$token_map)
  testthat::expect_identical(out$n_leading_edge_mapped, 2L)
  testthat::expect_identical(
    out$leading_edge_proteins,
    "PG:neuron_soma:AAAA1111;PG:neuron_soma:BBBB2222"
  )
})

testthat::test_that("4b. the two identity spaces are mutually unreachable", {
  b <- gww_build_universe_bundle(universe_fixture(), "neuron_soma")
  # An accession under the canonical contract is a vocabulary violation.
  ev <- evidence_fixture("neuron_soma", "P11111", "protein_group_id")
  testthat::expect_error(
    gww_map_leading_edges(ev, b$token_map), "outside the PG:neuron_soma:"
  )
  # A ProteinGroupID is absent from the accession space, so the zero-mapped
  # guard fires rather than any value being silently aliased across spaces.
  ev2 <- evidence_fixture("neuron_soma", "PG:neuron_soma:AAAA1111", "accession")
  testthat::expect_error(
    gww_map_leading_edges(ev2, b$token_map),
    "no proteins in the canonical WGCNA universe"
  )
})

testthat::test_that("5. mixed vocabularies in one canonical leading edge fail closed", {
  b <- gww_build_universe_bundle(universe_fixture(), "neuron_soma")
  ev <- evidence_fixture(
    "neuron_soma", "PG:neuron_soma:AAAA1111;P22222", "protein_group_id"
  )
  testthat::expect_error(
    gww_map_leading_edges(ev, b$token_map), "outside the PG:neuron_soma:"
  )
})

testthat::test_that("5b. evidence without a resolved identity family fails closed", {
  b <- gww_build_universe_bundle(universe_fixture(), "neuron_soma")
  ev <- evidence_fixture(
    "neuron_soma", "PG:neuron_soma:AAAA1111", NA_character_
  )
  testthat::expect_error(
    gww_map_leading_edges(ev, b$token_map),
    "leading edge but no resolved leading-edge identity family"
  )
  # A row with no leading edge at all has nothing to route and is not an
  # identity violation; the mapped-count guard covers it instead.
  empty <- evidence_fixture("neuron_soma", NA_character_, NA_character_)
  testthat::expect_error(
    gww_map_leading_edges(empty, b$token_map),
    "no proteins in the canonical WGCNA universe"
  )
  ev2 <- evidence_fixture("neuron_soma", "PG:neuron_soma:AAAA1111", "accession")
  ev2$leading_edge_token_type <- NULL
  testthat::expect_error(
    gww_map_leading_edges(ev2, b$token_map),
    "missing the leading-edge identity"
  )
})

testthat::test_that("6. a duplicated ProteinGroupID in a dataset universe fails validation", {
  dup <- universe_fixture(
    "neuron_soma",
    pgids = c("PG:neuron_soma:AAAA1111", "PG:neuron_soma:AAAA1111",
              "PG:neuron_soma:CCCC3333")
  )
  testthat::expect_error(
    gww_build_universe_bundle(dup, "neuron_soma"), "duplicated"
  )
})

testthat::test_that("6b. a universe ProteinGroupID outside the namespace fails validation", {
  bad <- universe_fixture(
    "neuron_soma",
    pgids = c("PG:neuron_soma:AAAA1111", "PG:microglia:BBBB2222",
              "PG:neuron_soma:CCCC3333")
  )
  testthat::expect_error(
    gww_build_universe_bundle(bad, "neuron_soma"),
    "outside the PG:neuron_soma: namespace"
  )
})

testthat::test_that("7. canonical mapping does not alter module membership", {
  u <- universe_fixture()
  b <- gww_build_universe_bundle(u, "neuron_soma")
  testthat::expect_identical(b$module_sets, split(u$ProteinGroupID, u$ModuleID))
  before <- b$module_sets
  invisible(gww_map_leading_edges(
    evidence_fixture(
      "neuron_soma", "PG:neuron_soma:AAAA1111", "protein_group_id"
    ),
    b$token_map
  ))
  testthat::expect_identical(b$module_sets, before)
  testthat::expect_setequal(
    b$module_sets[["WGCNA_m01"]],
    c("PG:neuron_soma:AAAA1111", "PG:neuron_soma:BBBB2222")
  )
})

testthat::test_that("8. a previously failing canonical leading edge now maps", {
  b <- gww_build_universe_bundle(universe_fixture(), "neuron_soma")
  ev <- evidence_fixture(
    "neuron_soma",
    "PG:neuron_soma:AAAA1111;PG:neuron_soma:BBBB2222;PG:neuron_soma:CCCC3333",
    "protein_group_id"
  )
  out <- gww_map_leading_edges(ev, b$token_map)
  testthat::expect_gte(out$n_leading_edge_mapped, 1L)
  testthat::expect_identical(out$n_leading_edge_mapped, 3L)

  # Before typed routing the map held accessions only, so this leading edge
  # mapped to zero proteins and tripped the guard. Reproducing that state
  # reproduces the original failure exactly.
  accession_only <- b$token_map[b$token_map$token_type == "accession", ]
  testthat::expect_error(
    gww_map_leading_edges(ev, accession_only),
    "no proteins in the canonical WGCNA universe"
  )
})

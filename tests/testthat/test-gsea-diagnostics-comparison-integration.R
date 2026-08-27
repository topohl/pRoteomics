testthat::local_edition(3)

# Integration-level: exercise the real Stage-06 diagnostics chain
#   assert columns -> resolve comparison identity -> transmute -> validate
#   -> build local/recurrent GSEA evidence
# against each accepted comparison vocabulary, using the real helpers rather
# than re-implementing them.

source(testthat::test_path(
  "..", "..", "R", "gsea_wgcna_concordance_utils.R"
))

CANONICAL_FAMILY <- "canonical_compareGO_ranked_GSEA_GO_BP"

# The ontology-aware theme-assignment export shape Stage 06 reads, restricted
# to the columns the diagnostics transmute dereferences.
theme_export_fixture <- function(comparison_style = c("source", "historical",
                                                     "both", "conflicting",
                                                     "neither")) {
  comparison_style <- match.arg(comparison_style)
  units <- c("ca1", "ca1", "ca2", "ca2", "ca3", "ca3")
  cmp <- paste0(units, "_sus_res")
  x <- data.frame(
    dataset = "microglia",
    phenotype_contrast = "SUS_vs_RES",
    spatial_unit = units,
    GO_ID = c("GO:0002181", "GO:0006364", "GO:0002181", "GO:0006364",
              "GO:0002181", "GO:0006364"),
    GO_description = c("cytoplasmic translation", "rRNA processing",
                       "cytoplasmic translation", "rRNA processing",
                       "cytoplasmic translation", "rRNA processing"),
    NES = c(1.9, 1.7, 1.8, 1.6, 1.85, 1.65),
    raw_p = 0.001,
    GSEA_FDR = 0.01,
    leading_edge_proteins = paste(
      "PG:microglia:AAAA1111", "PG:microglia:BBBB2222", sep = ";"
    ),
    leading_edge_genes = "Rps3;Rpl4",
    evidence_source_family = CANONICAL_FAMILY,
    theme_id = "ribosome_translation",
    manuscript_theme = "Ribosome and translation",
    theme_role = "primary",
    theme_claim_eligible = TRUE,
    anchor_GO_ID = "GO:0002181",
    mapping_type = "go_id_ontology",
    registry_version = "fixture_v1",
    theme_assignment_id = paste0("fixture_", seq_along(units)),
    stringsAsFactors = FALSE
  )
  if (comparison_style == "source") {
    x$source_comparison <- cmp
  } else if (comparison_style == "historical") {
    x$comparison <- cmp
  } else if (comparison_style == "both") {
    x$source_comparison <- cmp
    x$comparison <- cmp
  } else if (comparison_style == "conflicting") {
    x$source_comparison <- cmp
    x$comparison <- c(cmp[1], "ca2_res_con", cmp[3:6])
  }
  x
}

# Mirrors the Stage-06 assertion + resolution + transmute sequence exactly.
diagnostics_theme_terms <- function(gsea_terms) {
  gww_assert_columns(
    gsea_terms,
    c(
      "dataset", "phenotype_contrast", "spatial_unit", "GO_ID",
      "GO_description", "NES", "raw_p", "GSEA_FDR", "theme_id",
      "manuscript_theme", "theme_role", "theme_claim_eligible",
      "theme_assignment_id", "leading_edge_proteins", "leading_edge_genes",
      "evidence_source_family", "anchor_GO_ID", "mapping_type",
      "registry_version"
    ),
    "Ontology-aware all-contrast GSEA theme assignments"
  )
  gsea_terms <- gww_resolve_comparison_field(
    gsea_terms, "Ontology-aware all-contrast GSEA theme assignments"
  )
  out <- gsea_terms |>
    dplyr::transmute(
      dataset, phenotype_contrast, spatial_unit,
      dplyr::across(dplyr::any_of("source_comparison")),
      comparison = .data$comparison,
      ID = .data$GO_ID, Description = .data$GO_description,
      GO_ID, GO_description, NES,
      pvalue = .data$raw_p, p.adjust = .data$GSEA_FDR,
      core_enrichment = .data$leading_edge_proteins,
      core_enrichment_gene = .data$leading_edge_genes,
      evidence_source_family,
      theme_id, manuscript_theme, theme_role, theme_claim_eligible,
      anchor_GO_IDs = .data$anchor_GO_ID,
      mapping_method = .data$mapping_type,
      registry_version, theme_assignment_id,
      source_supplementary_file = NA_character_
    )
  out$leading_edge_token_type <- gww_leading_edge_token_type(
    out$evidence_source_family
  )
  gww_validate_theme_terms(out)
  out
}

build_evidence <- function(terms) {
  loc <- gww_build_local_gsea_evidence(terms)
  rec <- gww_build_recurrent_cross_spatial_gsea_evidence(
    terms, loc, min_recurrent_units = 2L
  )
  dplyr::bind_rows(loc, rec)
}

testthat::test_that("source_comparison-only diagnostics input works end-to-end", {
  terms <- diagnostics_theme_terms(theme_export_fixture("source"))
  testthat::expect_true("comparison" %in% names(terms))
  testthat::expect_true("source_comparison" %in% names(terms))
  testthat::expect_identical(terms$comparison, terms$source_comparison)
  ev <- build_evidence(terms)
  testthat::expect_gt(nrow(ev), 0L)
  testthat::expect_true(all(
    ev$leading_edge_token_type[!is.na(ev$leading_edge_accessions)] ==
      "protein_group_id"
  ))
})

testthat::test_that("historical comparison-only diagnostics input works end-to-end", {
  terms <- diagnostics_theme_terms(theme_export_fixture("historical"))
  testthat::expect_true("comparison" %in% names(terms))
  # nothing fabricates a source_comparison provenance that the export lacked
  testthat::expect_false("source_comparison" %in% names(terms))
  ev <- build_evidence(terms)
  testthat::expect_gt(nrow(ev), 0L)
})

testthat::test_that("both vocabularies present and identical are accepted", {
  terms <- diagnostics_theme_terms(theme_export_fixture("both"))
  testthat::expect_identical(terms$comparison, terms$source_comparison)
  testthat::expect_gt(nrow(build_evidence(terms)), 0L)
})

testthat::test_that("conflicting dual comparison fields fail closed", {
  testthat::expect_error(
    diagnostics_theme_terms(theme_export_fixture("conflicting")),
    "conflicting comparison identity"
  )
})

testthat::test_that("no comparison field at all fails closed", {
  testthat::expect_error(
    diagnostics_theme_terms(theme_export_fixture("neither")),
    "missing a comparison identity column"
  )
})

testthat::test_that("the two vocabularies yield identical downstream evidence", {
  a <- build_evidence(diagnostics_theme_terms(theme_export_fixture("source")))
  b <- build_evidence(
    diagnostics_theme_terms(theme_export_fixture("historical"))
  )
  # Provenance may differ in name; every derived value must not.
  drop <- c("source_comparison")
  common <- setdiff(intersect(names(a), names(b)), drop)
  testthat::expect_setequal(setdiff(names(a), drop), setdiff(names(b), drop))
  testthat::expect_equal(
    a[order(a$gsea_evidence_id), common],
    b[order(b$gsea_evidence_id), common],
    ignore_attr = "row.names"
  )
})

testthat::test_that("the required-column assertion covers every dereferenced field", {
  dereferenced <- c(
    "dataset", "phenotype_contrast", "spatial_unit", "GO_ID",
    "GO_description", "NES", "raw_p", "GSEA_FDR", "leading_edge_proteins",
    "leading_edge_genes", "evidence_source_family", "theme_id",
    "manuscript_theme", "theme_role", "theme_claim_eligible", "anchor_GO_ID",
    "mapping_type", "registry_version", "theme_assignment_id"
  )
  for (col in dereferenced) {
    x <- theme_export_fixture("source")
    x[[col]] <- NULL
    testthat::expect_error(
      diagnostics_theme_terms(x), "missing required column", info = col
    )
  }
})

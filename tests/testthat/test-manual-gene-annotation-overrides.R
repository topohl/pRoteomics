source(testthat::test_path("..", "..", "R", "paths.R"))

# Accession-keyed, versioned gene annotation overrides (format of 2026-11-16). These tests use the installed
# org.Mm.eg.db, because the override symbol must be that release's symbol for the Entrez ID.

skip_without_orgdb <- function() {
  testthat::skip_if_not_installed("org.Mm.eg.db")
  testthat::skip_if_not_installed("AnnotationDbi")
  testthat::skip_if_not_installed("readr")
}
release <- function() paste0("org.Mm.eg.db ", as.character(utils::packageVersion("org.Mm.eg.db")))
# a valid one-row table for Q8R092; columns given in `rows` replace (or add to) the defaults, one row per row of `rows`
override_file <- function(rows) {
  base <- list(member_accession = "Q8R092", official_gene_symbol = "4933434E20Rik", official_entrez_id = "99650",
    current_official_symbol = "Ltap1", mgi_id = "MGI:1914027", override_status = "approved", shared_gene_reason = "",
    annotation_release = release())
  df <- as.data.frame(lapply(base, function(v) rep(v, nrow(rows))), stringsAsFactors = FALSE)
  for (nm in names(rows)) df[[nm]] <- rows[[nm]]
  f <- tempfile(fileext = ".csv")
  utils::write.csv(df, f, row.names = FALSE)
  f
}
maps_for <- function(accessions) {
  build_mouse_gene_annotation_maps(getExportedValue("org.Mm.eg.db", "org.Mm.eg.db"), accessions = accessions)
}

testthat::test_that("an approved accession override resolves an accession the UniProt lookup cannot reach", {
  skip_without_orgdb()
  source(repo_path("R", "protein_mapping_utils.R"))
  ov <- read_manual_gene_annotation_overrides(override_file(data.frame(member_accession = "Q8R092", stringsAsFactors = FALSE)))
  testthat::expect_equal(attr(ov, "status"), "loaded")
  testthat::expect_equal(nrow(ov), 1L)
  # as in production, Q8R092 has no Gene_Name line, so it is not in the accession index of the maps
  res <- resolve_mouse_gene_annotation("Q8R092", "CA043", maps_for("P17182"), ov)
  testthat::expect_equal(res$gene_annotation_status, "resolved")
  testthat::expect_equal(res$gene_annotation_strategy, "manual_annotation_override")
  testthat::expect_equal(res$official_gene_symbol, "4933434E20Rik")
  testthat::expect_equal(res$official_entrez_id, "99650")
  testthat::expect_true(res$gene_annotation_manual_override_used)
  # without the override the same accession stays unresolved
  testthat::expect_equal(resolve_mouse_gene_annotation("Q8R092", "CA043", maps_for("P17182"), NULL)$gene_annotation_status, "unresolved")
})

testthat::test_that("proposed and held rows are validated but not applied", {
  skip_without_orgdb()
  source(repo_path("R", "protein_mapping_utils.R"))
  f <- override_file(data.frame(member_accession = c("Q8R092", "P23198"), official_gene_symbol = c("4933434E20Rik", "Cbx3"),
    official_entrez_id = c("99650", "12417"), override_status = c("proposed", "held"), stringsAsFactors = FALSE))
  ov <- read_manual_gene_annotation_overrides(f)
  testthat::expect_equal(nrow(ov), 0L)
  testthat::expect_equal(attr(ov, "rows_not_applied"), 2L)
  testthat::expect_equal(resolve_mouse_gene_annotation("Q8R092", "CA043", maps_for("P17182"), ov)$gene_annotation_status, "unresolved")
})

testthat::test_that("a current symbol that the installed release does not know is rejected", {
  skip_without_orgdb()
  source(repo_path("R", "protein_mapping_utils.R"))
  f <- override_file(data.frame(member_accession = "Q8R092", official_gene_symbol = "Ltap1", stringsAsFactors = FALSE))
  testthat::expect_error(read_manual_gene_annotation_overrides(f), "unknown to the installed")
})

testthat::test_that("the analysis symbol must belong to the Entrez ID", {
  skip_without_orgdb()
  source(repo_path("R", "protein_mapping_utils.R"))
  f <- override_file(data.frame(member_accession = "Q8R092", official_gene_symbol = "AU022252", stringsAsFactors = FALSE))
  testthat::expect_error(read_manual_gene_annotation_overrides(f), "does not belong to the Entrez ID")
})

testthat::test_that("overrides are keyed by accession only", {
  skip_without_orgdb()
  source(repo_path("R", "protein_mapping_utils.R"))
  testthat::expect_error(read_manual_gene_annotation_overrides(override_file(data.frame(member_accession = "CA043", stringsAsFactors = FALSE))),
    "not a UniProt accession")
  f <- override_file(data.frame(member_accession = "Q8R092", lookup_value = "CA043", stringsAsFactors = FALSE))
  testthat::expect_error(read_manual_gene_annotation_overrides(f), "only member_accession may key")
})

testthat::test_that("missing Entrez IDs, other releases, bad status and duplicated accessions are rejected", {
  skip_without_orgdb()
  source(repo_path("R", "protein_mapping_utils.R"))
  testthat::expect_error(read_manual_gene_annotation_overrides(override_file(data.frame(member_accession = "Q8R092", official_entrez_id = "", stringsAsFactors = FALSE))),
    "numeric Entrez ID")
  testthat::expect_error(read_manual_gene_annotation_overrides(override_file(data.frame(member_accession = "Q8R092", annotation_release = "org.Mm.eg.db 3.19.1", stringsAsFactors = FALSE))),
    "does not match the installed")
  testthat::expect_error(read_manual_gene_annotation_overrides(override_file(data.frame(member_accession = "Q8R092", override_status = "maybe", stringsAsFactors = FALSE))),
    "override_status must be one of")
  testthat::expect_error(read_manual_gene_annotation_overrides(override_file(data.frame(member_accession = c("Q8R092", "Q8R092"), stringsAsFactors = FALSE))),
    "duplicated accession")
  f <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(member_accession = "Q8R092", official_gene_symbol = "4933434E20Rik", official_entrez_id = "99650"), f, row.names = FALSE)
  testthat::expect_error(read_manual_gene_annotation_overrides(f), "versioned format required")
})

testthat::test_that("several accessions may share one Entrez ID only with a stated reason", {
  skip_without_orgdb()
  source(repo_path("R", "protein_mapping_utils.R"))
  rows <- data.frame(member_accession = c("Q8R092", "Q05D39"), stringsAsFactors = FALSE)
  testthat::expect_error(read_manual_gene_annotation_overrides(override_file(rows)), "without shared_gene_reason")
  rows$shared_gene_reason <- "isoform entries of one gene in one protein group"
  ov <- read_manual_gene_annotation_overrides(override_file(rows))
  testthat::expect_equal(nrow(ov), 2L)
})

testthat::test_that("an override never replaces a unique UniProt mapping", {
  skip_without_orgdb()
  source(repo_path("R", "protein_mapping_utils.R"))
  f <- override_file(data.frame(member_accession = "Q8C3W1", official_gene_symbol = "AU022252", official_entrez_id = "230696", stringsAsFactors = FALSE))
  ov <- read_manual_gene_annotation_overrides(f)
  res <- resolve_mouse_gene_annotation("Q8C3W1", "CA198", maps_for("Q8C3W1"), ov)
  testthat::expect_equal(res$gene_annotation_strategy, "unique_uniprot_to_symbol")
  testthat::expect_equal(res$official_gene_symbol, "2310022B05Rik")
  testthat::expect_false(res$gene_annotation_manual_override_used)
})

testthat::test_that("shared genes: pre-existing ones are reported, override-created ones stop unless explained", {
  source(repo_path("R", "protein_mapping_utils.R"))
  wide <- data.frame(ProteinGroupID = c("PG:a", "PG:b", "PG:c", "PG:d"), member_accessions = c("D3Z2H9", "P21107", "Q8R092", "Q05D39"),
    official_gene_symbol = c("Tpm3", "Tpm3", "4933434E20Rik", "4933434E20Rik"), official_entrez_id = c("59069", "59069", "99650", "99650"),
    gene_level_claim_allowed = TRUE, stringsAsFactors = FALSE)
  bridge <- data.frame(ProteinGroupID = wide$ProteinGroupID, gene_annotation_manual_override_used = c(FALSE, FALSE, TRUE, FALSE))
  testthat::expect_error(audit_override_gene_duplicates(wide, bridge, NULL), "without shared_gene_reason")
  ok <- audit_override_gene_duplicates(wide[1:2, ], bridge[1:2, ], NULL)
  testthat::expect_equal(sort(unique(ok$official_entrez_id)), "59069")
  testthat::expect_false(any(ok$override_involved))
  ov <- data.frame(lookup_value = "Q8R092", official_gene_symbol = "4933434E20Rik", official_entrez_id = "99650",
    shared_gene_reason = "isoform entries of one gene", stringsAsFactors = FALSE)
  res <- audit_override_gene_duplicates(wide, bridge, ov)
  testthat::expect_true(all(c("59069", "99650") %in% res$official_entrez_id))
})

testthat::test_that("the review override table itself passes validation", {
  skip_without_orgdb()
  source(repo_path("R", "protein_mapping_utils.R"))
  f <- repo_path("data", "metadata", "manual_gene_annotation_overrides.csv")
  testthat::skip_if_not(file.exists(f))
  ov <- read_manual_gene_annotation_overrides(f)
  raw <- utils::read.csv(f, colClasses = "character")
  testthat::expect_equal(nrow(ov), sum(raw$override_status == "approved"))
  testthat::expect_false(anyDuplicated(raw$member_accession) > 0)
})

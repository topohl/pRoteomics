source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "enrichment_io.R"))

load_direction_comparison_parser <- function() {
  script <- testthat::test_path(
    "..", "..", "04_differential_expression_enrichment",
    "01b_gsea_protein_direction_audit.r"
  )
  expressions <- parse(script)
  definition <- Filter(function(expr) {
    is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      is.symbol(expr[[2]]) && identical(as.character(expr[[2]]), "parse_comparison_sides")
  }, as.list(expressions))
  testthat::expect_length(definition, 1L)
  environment <- new.env(parent = baseenv())
  eval(definition[[1]], envir = environment)
  environment$parse_comparison_sides
}

load_direction_status_helpers <- function() {
  script <- testthat::test_path(
    "..", "..", "04_differential_expression_enrichment",
    "01b_gsea_protein_direction_audit.r"
  )
  targets <- c(
    "parse_comparison_sides",
    "make_status_summary_rows",
    "assert_bind_rows_compatible",
    "write_csv_in_chunks"
  )
  definitions <- Filter(function(expr) {
    is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% targets
  }, as.list(parse(script)))
  testthat::expect_equal(
    sort(vapply(definitions, function(expr) as.character(expr[[2]]), character(1))),
    sort(targets)
  )
  environment <- new.env(parent = baseenv())
  for (definition in definitions) eval(definition, envir = environment)
  mget(targets, envir = environment, inherits = FALSE)
}

evaluate_direction_summary_constructor <- function(term_summary) {
  script <- testthat::test_path(
    "..", "..", "04_differential_expression_enrichment",
    "01b_gsea_protein_direction_audit.r"
  )
  definitions <- Filter(function(expr) {
    is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      is.symbol(expr[[2]]) && identical(as.character(expr[[2]]), "summary_with_terms")
  }, as.list(parse(script)))
  testthat::expect_length(definitions, 1L)
  environment <- new.env(parent = baseenv())
  environment$term_summary <- term_summary
  eval(definitions[[1]], envir = environment)
  environment$summary_with_terms
}

direction_term_summary_fixture <- function() {
  data.frame(
    dataset = "microglia",
    ontology = "BP",
    contrast = "CA1res_CA1con",
    formal_contrast = "CA1res_CA1con",
    positive_side_label = "RES",
    negative_side_label = "CON",
    analysis_status = "success_with_terms",
    classification = "consistent",
    stringsAsFactors = FALSE
  )
}

direction_summary_with_terms_fixture <- function() {
  evaluate_direction_summary_constructor(direction_term_summary_fixture())
}

write_contract_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
  path
}

empty_gsea_fixture <- function() {
  data.frame(
    ID = character(), Description = character(), NES = numeric(), p.adjust = numeric(),
    setSize = integer(), core_enrichment = character(), stringsAsFactors = FALSE
  )
}

empty_term_provenance_fixture <- function() {
  as.data.frame(setNames(lapply(term_gene_provenance_columns(), function(column) {
    if (column %in% c("gene_level_claim_allowed", "core_enrichment_member")) logical()
    else if (column == "rank_statistic") numeric()
    else character()
  }), term_gene_provenance_columns()), stringsAsFactors = FALSE)
}

make_downstream_fixture <- function(root, reverse_rows = FALSE, zero_only = FALSE) {
  dataset <- "microglia"
  rel <- function(...) gsub("\\\\", "/", file.path(...))
  base <- file.path(root, "fixture")
  terms_path <- file.path(base, "misleading_filename.csv")
  zero_terms_path <- file.path(base, "zero.csv")
  collapsed_path <- file.path(base, "collapsed.csv")
  zero_collapsed_path <- file.path(base, "zero_collapsed.csv")
  collapsed_provenance_path <- file.path(base, "collapsed_provenance.csv")
  provenance_path <- file.path(base, "term_provenance.csv")
  zero_provenance_path <- file.path(base, "zero_provenance.csv")

  terms <- data.frame(
    ID = "GO:0001", Description = "synaptic organization", NES = 2,
    p.adjust = 0.01, setSize = 1L, core_enrichment = "Aif1",
    stringsAsFactors = FALSE
  )
  collapsed <- data.frame(
    official_gene_symbol = "Aif1", official_entrez_id = "00123",
    collapsed_statistic = 2, collapsed_logfc = 1.5,
    contributing_ProteinGroupIDs = "PG_A;PG_B", stringsAsFactors = FALSE
  )
  provenance <- data.frame(
    dataset = dataset, comparison = "manifest_identity", result_type = "GSEA_GO",
    ontology = "BP", term_id = "GO:0001", term_description = "synaptic organization",
    official_gene_symbol = "Aif1", official_entrez_id = "00123",
    ProteinGroupID = c("PG_A", "PG_B"), member_accessions = c("P1", "P2"),
    protein_group_gene_annotation_status = "concordant_official_gene",
    gene_level_claim_allowed = TRUE, rank_statistic = 2,
    core_enrichment_member = TRUE,
    enrichment_contract_version = "protein_group_enrichment_v3_term_gene_provenance",
    gene_annotation_contract_version = "mouse_gene_annotation_v1",
    stringsAsFactors = FALSE
  )
  if (reverse_rows) provenance <- provenance[rev(seq_len(nrow(provenance))), , drop = FALSE]
  write_contract_csv(terms, terms_path)
  write_contract_csv(empty_gsea_fixture(), zero_terms_path)
  write_contract_csv(collapsed, collapsed_path)
  write_contract_csv(collapsed, zero_collapsed_path)
  write_contract_csv(collapsed, collapsed_provenance_path)
  write_contract_csv(provenance, provenance_path)
  write_contract_csv(empty_term_provenance_fixture(), zero_provenance_path)

  rows <- data.frame(
    dataset = dataset,
    comparison = c("manifest_identity", "zero_identity", "failed_identity"),
    result_type = "GSEA_GO", ontology = "BP",
    analysis_status = c("success_with_terms", "success_zero_terms", "failed"),
    n_terms = c(1L, 0L, NA_integer_),
    output_table = c(rel("fixture", "misleading_filename.csv"), rel("fixture", "zero.csv"), NA),
    collapsed_gene_input_file = c(rel("fixture", "collapsed.csv"), rel("fixture", "zero_collapsed.csv"), NA),
    collapsed_gene_provenance_file = c(rel("fixture", "collapsed_provenance.csv"), rel("fixture", "collapsed_provenance.csv"), NA),
    term_gene_provenance_file = c(rel("fixture", "term_provenance.csv"), rel("fixture", "zero_provenance.csv"), NA),
    enrichment_contract_version = canonical_clusterprofiler_manifest_contract_version(),
    gene_annotation_contract_version = "mouse_gene_annotation_v1",
    route_category = "phenotype_within_unit", route_unit = "CA1",
    stringsAsFactors = FALSE
  )
  if (zero_only) rows <- rows[rows$comparison == "zero_identity", , drop = FALSE]
  if (reverse_rows) rows <- rows[rev(seq_len(nrow(rows))), , drop = FALSE]
  manifest_path <- canonical_clusterprofiler_manifest_path(dataset, root)
  write_contract_csv(rows, manifest_path)
  list(root = root, dataset = dataset, manifest = manifest_path, rows = rows,
    provenance = provenance, terms = terms)
}

make_comparego_fixture <- function(root, cluster_fixture) {
  base <- file.path(root, "compare")
  term_path <- write_contract_csv(data.frame(
    ID = "GO:0001", Description = "synaptic organization", NES = 2, p.adjust = 0.01,
    setSize = 1L, core_enrichment = "Aif1", term_id = "GO:0001",
    term_description = "synaptic organization", dataset = "microglia",
    comparison = "manifest_identity", result_type = "GSEA_GO", ontology = "BP",
    stringsAsFactors = FALSE
  ), file.path(base, "terms.csv"))
  provenance_path <- write_contract_csv(cluster_fixture$provenance, file.path(base, "provenance.csv"))
  status_path <- write_contract_csv(data.frame(
    dataset = "microglia", comparison = c("manifest_identity", "zero_identity", "failed_identity"),
    result_type = "GSEA_GO", ontology = "BP",
    analysis_status = c("success_with_terms", "success_zero_terms", "failed"),
    n_terms = c(1L, 0L, NA_integer_),
    comparego_action = c("included", "completed_zero_terms", "recorded_failed"),
    stringsAsFactors = FALSE
  ), file.path(base, "status.csv"))
  rel <- function(path) gsub("\\\\", "/", substring(path, nchar(root) + 2L))
  manifest <- data.frame(
    dataset = "microglia", comparison = c("manifest_identity", "zero_identity", "failed_identity"),
    result_type = "GSEA_GO", ontology = "BP",
    analysis_status = c("success_with_terms", "success_zero_terms", "failed"),
    comparego_analysis_status = c("included", "completed_zero_terms", "recorded_failed"),
    input_manifest = rel(cluster_fixture$manifest),
    term_comparison_file = rel(term_path),
    term_gene_provenance_output_file = rel(provenance_path),
    analysis_status_summary_file = rel(status_path),
    enrichment_contract_version = canonical_clusterprofiler_manifest_contract_version(),
    comparego_contract_version = canonical_comparego_manifest_contract_version(),
    route_category = "phenotype_within_unit", route_unit = "CA1",
    stringsAsFactors = FALSE
  )
  path <- canonical_comparego_manifest_path("microglia", root)
  write_contract_csv(manifest, path)
  list(path = path, manifest = manifest)
}

testthat::test_that("canonical manifest path ignores newer decoys", {
  root <- tempfile("b2a-canonical-")
  fixture <- make_downstream_fixture(root)
  decoy <- write_contract_csv(fixture$rows, file.path(root, "decoys", "clusterProfiler_manifest_newer.csv"))
  Sys.setFileTime(decoy, Sys.time() + 3600)
  testthat::expect_equal(
    normalizePath(canonical_clusterprofiler_manifest_path("microglia", root), winslash = "/"),
    normalizePath(fixture$manifest, winslash = "/")
  )
})

testthat::test_that("clusterProfiler bundle preserves manifest identity and all contributors", {
  fixture <- make_downstream_fixture(tempfile("b2a-bundle-"))
  bundle <- read_canonical_clusterprofiler_bundle(
    fixture$manifest, fixture$dataset, repository_root = fixture$root
  )
  testthat::expect_equal(bundle$terms$comparison, "manifest_identity")
  testthat::expect_setequal(bundle$provenance$ProteinGroupID, c("PG_A", "PG_B"))
  testthat::expect_equal(unique(bundle$provenance$official_gene_symbol), "Aif1")
  testthat::expect_type(bundle$provenance$official_entrez_id, "character")
  testthat::expect_equal(unique(bundle$provenance$official_entrez_id), "00123")
  testthat::expect_setequal(bundle$status$analysis_status, c("success_with_terms", "success_zero_terms", "failed"))
  contributors <- join_term_provenance_to_collapsed_genes(bundle$provenance, bundle$collapsed)
  testthat::expect_setequal(contributors$ProteinGroupID, c("PG_A", "PG_B"))
  testthat::expect_equal(unique(contributors$collapsed_logfc), 1.5)
})

testthat::test_that("clusterProfiler bundle is deterministic under reversed input order", {
  a <- make_downstream_fixture(tempfile("b2a-order-a-"), reverse_rows = FALSE)
  b <- make_downstream_fixture(tempfile("b2a-order-b-"), reverse_rows = TRUE)
  out_a <- read_canonical_clusterprofiler_bundle(a$manifest, a$dataset, repository_root = a$root)
  out_b <- read_canonical_clusterprofiler_bundle(b$manifest, b$dataset, repository_root = b$root)
  comparable <- c("dataset", "comparison", "ontology", "term_id", "official_gene_symbol", "ProteinGroupID")
  testthat::expect_equal(out_a$provenance[comparable], out_b$provenance[comparable])
  testthat::expect_equal(out_a$status, out_b$status)
})

testthat::test_that("strict contracts reject stale, missing, malformed, and ineligible provenance", {
  stale <- make_downstream_fixture(tempfile("b2a-stale-"))
  stale_rows <- utils::read.csv(stale$manifest, stringsAsFactors = FALSE)
  stale_rows$enrichment_contract_version <- "stale"
  write_contract_csv(stale_rows, stale$manifest)
  testthat::expect_error(
    read_canonical_clusterprofiler_bundle(stale$manifest, stale$dataset, repository_root = stale$root),
    "Stale clusterProfiler manifest contract"
  )

  missing <- make_downstream_fixture(tempfile("b2a-missing-"))
  unlink(file.path(missing$root, "fixture", "term_provenance.csv"))
  testthat::expect_error(
    read_canonical_clusterprofiler_bundle(missing$manifest, missing$dataset, repository_root = missing$root),
    "missing term_gene_provenance_file"
  )

  malformed <- make_downstream_fixture(tempfile("b2a-malformed-"))
  write_contract_csv(data.frame(dataset = "microglia"), file.path(malformed$root, "fixture", "term_provenance.csv"))
  testthat::expect_error(
    read_canonical_clusterprofiler_bundle(malformed$manifest, malformed$dataset, repository_root = malformed$root),
    "Term-gene provenance is missing required columns"
  )

  ineligible <- make_downstream_fixture(tempfile("b2a-ineligible-"))
  bad <- ineligible$provenance
  bad$gene_level_claim_allowed[[1]] <- FALSE
  bad$protein_group_gene_annotation_status[[1]] <- "ambiguous"
  write_contract_csv(bad, file.path(ineligible$root, "fixture", "term_provenance.csv"))
  testthat::expect_error(
    read_canonical_clusterprofiler_bundle(ineligible$manifest, ineligible$dataset, repository_root = ineligible$root),
    "Ineligible protein groups"
  )

  ambiguous <- make_downstream_fixture(tempfile("b2a-ambiguous-"))
  bad <- ambiguous$provenance
  bad$protein_group_gene_annotation_status[[1]] <- "conflicting_member_annotations"
  write_contract_csv(bad, file.path(ambiguous$root, "fixture", "term_provenance.csv"))
  testthat::expect_error(
    read_canonical_clusterprofiler_bundle(ambiguous$manifest, ambiguous$dataset, repository_root = ambiguous$root),
    "Ambiguous protein-group annotations"
  )
})

testthat::test_that("zero-term clusterProfiler bundle has valid typed empty contracts", {
  fixture <- make_downstream_fixture(tempfile("b2a-zero-"), zero_only = TRUE)
  bundle <- read_canonical_clusterprofiler_bundle(
    fixture$manifest, fixture$dataset, repository_root = fixture$root
  )
  testthat::expect_equal(nrow(bundle$terms), 0L)
  testthat::expect_equal(names(bundle$provenance), term_gene_provenance_columns())
  testthat::expect_type(bundle$provenance$official_entrez_id, "character")
  testthat::expect_equal(bundle$status$analysis_status, "success_zero_terms")
})

testthat::test_that("compareGO bundle reads only manifest-declared canonical outputs", {
  fixture <- make_downstream_fixture(tempfile("b2a-compare-"))
  compare <- make_comparego_fixture(fixture$root, fixture)
  bundle <- read_canonical_comparego_bundle(compare$path, "microglia", repository_root = fixture$root)
  testthat::expect_equal(bundle$terms$comparison, "manifest_identity")
  testthat::expect_setequal(bundle$provenance$ProteinGroupID, c("PG_A", "PG_B"))
  testthat::expect_type(bundle$provenance$official_entrez_id, "character")
  testthat::expect_setequal(bundle$status$analysis_status, c("success_with_terms", "success_zero_terms", "failed"))

  stale <- utils::read.csv(compare$path, stringsAsFactors = FALSE)
  stale$comparego_contract_version <- "stale"
  write_contract_csv(stale, compare$path)
  testthat::expect_error(
    read_canonical_comparego_bundle(compare$path, "microglia", repository_root = fixture$root),
    "Stale compareGO manifest contract"
  )
})

testthat::test_that("active B2a scripts contain no legacy discovery or identifier remapping", {
  paths <- testthat::test_path("..", "..", "04_differential_expression_enrichment",
    c("01b_gsea_protein_direction_audit.r", "04_neuropil_reference_annotation.r", "06_biological_program_summary.r"))
  scripts <- vapply(paths, function(path) paste(readLines(path, warn = FALSE), collapse = "\n"), character(1))
  testthat::expect_false(any(grepl("list\\.files\\(", scripts)))
  testthat::expect_false(any(grepl("latest_manifest\\s*<-", scripts)))
  testthat::expect_false(any(grepl("fromType\\s*=\\s*['\"]UNIPROT|keyType\\s*=\\s*['\"]UNIPROT", scripts)))
  testthat::expect_false(any(grepl("\\bbitr\\s*\\(|\\bmapIds\\s*\\(", scripts)))
  testthat::expect_false(any(grepl("!duplicated\\s*\\(", scripts)))
  testthat::expect_false(any(grepl("strsplit.*core_enrichment", scripts)))
  testthat::expect_match(scripts[[1]], "term_gene_provenance_file")
  testthat::expect_match(scripts[[2]], "prepare_terms\\(microglia_manifest\\$bundle\\)")
  testthat::expect_match(scripts[[3]], "read_canonical_comparego_bundle")
})

testthat::test_that("comparison-side parser is type-stable for empty and scalar inputs", {
  parse_sides <- load_direction_comparison_parser()
  empty <- parse_sides(character(0))
  testthat::expect_s3_class(empty, "data.frame")
  testthat::expect_equal(
    names(empty),
    c("formal_contrast", "positive_side_label", "negative_side_label")
  )
  testthat::expect_equal(nrow(empty), 0L)
  testthat::expect_true(all(vapply(empty, is.character, logical(1))))

  recognized <- parse_sides("CA1res_CA1con")
  testthat::expect_equal(recognized$formal_contrast, "CA1res_CA1con")
  testthat::expect_equal(recognized$positive_side_label, "RES")
  testthat::expect_equal(recognized$negative_side_label, "CON")

  unknown <- parse_sides("unknown_comparison")
  testthat::expect_equal(unknown$positive_side_label, "positive_log2fc_side")
  testthat::expect_equal(unknown$negative_side_label, "negative_log2fc_side")

  missing <- parse_sides(NA_character_)
  testthat::expect_true(is.na(missing$formal_contrast))
  testthat::expect_equal(missing$positive_side_label, "positive_log2fc_side")
  testthat::expect_equal(missing$negative_side_label, "negative_log2fc_side")

  blank <- parse_sides("")
  testthat::expect_identical(blank$formal_contrast, "")
  testthat::expect_equal(blank$positive_side_label, "positive_log2fc_side")
  testthat::expect_equal(blank$negative_side_label, "negative_log2fc_side")
})

testthat::test_that("comparison-side parser maps vectors independently and preserves order", {
  parse_sides <- load_direction_comparison_parser()
  comparisons <- c("CA1res_CA1con", "CA1sus_CA1res", "CA1sus_CA1con")
  expected_positive <- c("RES", "SUS", "SUS")
  expected_negative <- c("CON", "RES", "CON")
  parsed <- parse_sides(comparisons)
  testthat::expect_equal(nrow(parsed), length(comparisons))
  testthat::expect_equal(parsed$formal_contrast, comparisons)
  testthat::expect_equal(parsed$positive_side_label, expected_positive)
  testthat::expect_equal(parsed$negative_side_label, expected_negative)

  reversed <- parse_sides(rev(comparisons))
  expected_reversed <- parsed[rev(seq_len(nrow(parsed))), , drop = FALSE]
  rownames(expected_reversed) <- NULL
  testthat::expect_equal(reversed, expected_reversed)
})

testthat::test_that("canonical direction-audit attachment handles zero and multiple status rows", {
  parse_sides <- load_direction_comparison_parser()
  attach_sides <- function(input) {
    dplyr::bind_cols(input, parse_sides(input$comparison)) |>
      dplyr::transmute(
        comparison, analysis_status, formal_contrast,
        positive_side_label, negative_side_label
      )
  }

  empty_input <- data.frame(comparison = character(), analysis_status = character())
  empty_output <- attach_sides(empty_input)
  testthat::expect_equal(nrow(empty_output), 0L)
  testthat::expect_equal(
    names(empty_output),
    c("comparison", "analysis_status", "formal_contrast", "positive_side_label", "negative_side_label")
  )

  status_input <- data.frame(
    comparison = c("zero_identity", "CA1sus_CA1res"),
    analysis_status = c("success_zero_terms", "success_with_terms"),
    stringsAsFactors = FALSE
  )
  status_output <- attach_sides(status_input)
  testthat::expect_equal(nrow(status_output), 2L)
  testthat::expect_equal(status_output$analysis_status[[1]], "success_zero_terms")
  testthat::expect_equal(status_output$formal_contrast, status_input$comparison)
  testthat::expect_equal(status_output$positive_side_label, c("positive_log2fc_side", "SUS"))

  reversed <- attach_sides(status_input[2:1, , drop = FALSE])
  expected_reversed <- status_output[2:1, , drop = FALSE]
  rownames(expected_reversed) <- NULL
  testthat::expect_equal(reversed, expected_reversed)
})

testthat::test_that("direction summary status rows retain compatible explicit types", {
  helpers <- load_direction_status_helpers()
  summary_with_terms <- direction_summary_with_terms_fixture()
  status <- data.frame(
    dataset = "microglia",
    ontology = "BP",
    comparison = "CA1sus_CA1res",
    analysis_status = "success_zero_terms",
    stringsAsFactors = FALSE
  )
  status_only <- helpers$make_status_summary_rows(status)

  testthat::expect_no_error(
    helpers$assert_bind_rows_compatible(summary_with_terms, status_only, "test summary")
  )
  combined <- dplyr::bind_rows(summary_with_terms, status_only)
  testthat::expect_identical(names(status_only), names(summary_with_terms))
  testthat::expect_type(combined$interpretation, "character")
  testthat::expect_type(combined$frac_consistent, "double")
  testthat::expect_true(all(vapply(
    combined[c("n_terms", "n_consistent", "n_weak_or_mixed_core", "n_inconsistent")],
    is.integer,
    logical(1)
  )))
  testthat::expect_equal(combined$analysis_status, c("success_with_terms", "success_zero_terms"))
  testthat::expect_match(combined$interpretation[[2]], "zero terms", fixed = TRUE)
})

testthat::test_that("zero-term-only and zero-row direction summaries are type-stable", {
  helpers <- load_direction_status_helpers()
  empty_terms <- evaluate_direction_summary_constructor(
    direction_term_summary_fixture()[0, , drop = FALSE]
  )
  zero_term_status <- data.frame(
    dataset = "microglia",
    ontology = "BP",
    comparison = "zero_term_comparison",
    analysis_status = "success_zero_terms",
    stringsAsFactors = FALSE
  )
  status_only <- helpers$make_status_summary_rows(zero_term_status)
  testthat::expect_no_error(
    helpers$assert_bind_rows_compatible(empty_terms, status_only, "zero-term summary")
  )
  zero_term_summary <- dplyr::bind_rows(empty_terms, status_only)
  testthat::expect_equal(nrow(zero_term_summary), 1L)
  testthat::expect_identical(zero_term_summary$analysis_status, "success_zero_terms")
  testthat::expect_type(zero_term_summary$interpretation, "character")
  testthat::expect_type(zero_term_summary$frac_consistent, "double")
  testthat::expect_type(zero_term_summary$n_terms, "integer")

  empty_status <- zero_term_status[0, , drop = FALSE]
  empty_output <- helpers$make_status_summary_rows(empty_status)
  testthat::expect_equal(nrow(empty_output), 0L)
  testthat::expect_identical(names(empty_output), names(empty_terms))
  testthat::expect_type(empty_output$interpretation, "character")
  testthat::expect_type(empty_output$frac_consistent, "double")
  testthat::expect_type(empty_output$n_terms, "integer")
})

testthat::test_that("mixed direction statuses preserve failures and deterministic order", {
  helpers <- load_direction_status_helpers()
  summary_with_terms <- direction_summary_with_terms_fixture()
  statuses <- data.frame(
    dataset = rep("microglia", 3L),
    ontology = rep("BP", 3L),
    comparison = c("CA1res_CA1con", "zero_term_comparison", "failed_comparison"),
    analysis_status = c("success_with_terms", "success_zero_terms", "failed"),
    stringsAsFactors = FALSE
  )
  build_summary <- function(input) {
    status_only <- helpers$make_status_summary_rows(input)
    helpers$assert_bind_rows_compatible(summary_with_terms, status_only, "mixed summary")
    dplyr::bind_rows(summary_with_terms, status_only) |>
      dplyr::arrange(.data$dataset, .data$ontology, .data$contrast)
  }

  result <- build_summary(statuses)
  reversed <- build_summary(statuses[3:1, , drop = FALSE])
  testthat::expect_equal(reversed, result)
  testthat::expect_equal(
    sort(result$analysis_status),
    sort(c("success_with_terms", "success_zero_terms", "failed"))
  )
  failed <- result[result$analysis_status == "failed", , drop = FALSE]
  testthat::expect_equal(nrow(failed), 1L)
  testthat::expect_equal(failed$n_terms, 0L)
  testthat::expect_match(failed$interpretation, "failed", fixed = TRUE)
  testthat::expect_equal(
    nrow(result[result$analysis_status == "success_with_terms", , drop = FALSE]),
    1L
  )
})

testthat::test_that("direction audit chunked CSV writer preserves rows and typed empty schemas", {
  helpers <- load_direction_status_helpers()
  path <- tempfile(fileext = ".csv")
  input <- data.frame(
    ProteinGroupID = paste0("PG", 1:5),
    n_terms = 1:5,
    score = c(1.5, NA_real_, -2, 0, 3.25),
    stringsAsFactors = FALSE
  )
  helpers$write_csv_in_chunks(input, path, chunk_size = 2L)
  observed <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  testthat::expect_equal(observed, input)
  testthat::expect_equal(length(readLines(path, warn = FALSE)), nrow(input) + 1L)

  empty_path <- tempfile(fileext = ".csv")
  helpers$write_csv_in_chunks(input[0, , drop = FALSE], empty_path, chunk_size = 2L)
  empty <- utils::read.csv(empty_path, stringsAsFactors = FALSE, check.names = FALSE)
  testthat::expect_equal(nrow(empty), 0L)
  testthat::expect_identical(names(empty), names(input))
  testthat::expect_equal(length(readLines(empty_path, warn = FALSE)), 1L)
})

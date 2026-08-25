source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "enrichment_io.R"))

suppressPackageStartupMessages(library(dplyr))

load_assignment_helpers <- function(script, targets) {
  expressions <- as.list(parse(script))
  definitions <- Filter(function(expr) {
    is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% targets
  }, expressions)
  observed <- vapply(definitions, function(expr) as.character(expr[[2]]), character(1))
  testthat::expect_setequal(observed, targets)
  environment <- new.env(parent = globalenv())
  for (definition in definitions) eval(definition, envir = environment)
  environment
}

load_neuropil_annotation_helpers <- function() {
  script <- testthat::test_path(
    "..", "..", "04_differential_expression_enrichment",
    "04_neuropil_reference_annotation.r"
  )
  environment <- load_assignment_helpers(script, c(
    "REFERENCE_MATCHING_CONTRACT_VERSION",
    "CLASSIFICATION_CONTRACT_VERSION",
    "marker_sets",
    "normalize_id",
    "parse_reference_identity",
    "marker_fraction",
    "marker_hits",
    "term_overlap_fraction",
    "jaccard",
    "classify_term",
    "empty_neuropil_evidence",
    "safe_min",
    "safe_max",
    "find_best_neuropil_match",
    "annotate_microglia_terms"
  ))
  environment$REFERENCE_DATASET <- "neuron_neuropil"
  environment
}

load_neuropil_program_helpers <- function() {
  script <- testthat::test_path(
    "..", "..", "04_differential_expression_enrichment",
    "06_biological_program_summary.r"
  )
  load_assignment_helpers(script, c(
    "mode_value",
    "NEUROPIL_PROGRAM_ANNOTATION_CONTRACT_VERSION",
    "empty_neuropil_program_annotation",
    "safe_mean",
    "summarise_neuropil_annotation_by_program",
    "attach_neuropil_program_annotation"
  ))
}

make_reference_term <- function(
    dataset,
    comparison,
    route_unit,
    region,
    contrast,
    p_adjust,
    score,
    genes,
    term_key = "GSEA_GO::BP::GO:0001") {
  tibble::tibble(
    dataset = dataset,
    comparison = comparison,
    result_type = "GSEA_GO",
    ontology = "BP",
    term_id = "GO:0001",
    term_description = "test process",
    p_adjust = p_adjust,
    score = score,
    route_category = "phenotype_within_unit",
    route_unit = route_unit,
    source_table = paste0(comparison, ".csv"),
    term_key = term_key,
    anatomical_region = region,
    formal_group_contrast = contrast,
    reference_identity_status = "valid",
    reference_identity_key = paste(region, contrast, sep = "::"),
    neuropil_match_key = paste(term_key, region, contrast, sep = "::"),
    genes = list(genes)
  )
}

annotate_fixture <- function(helpers, microglia_term, neuropil_terms) {
  helpers$annotate_microglia_terms(microglia_term, neuropil_terms)
}

testthat::test_that("canonical reference identity parser uses route region and formal condition sides", {
  helpers <- load_neuropil_annotation_helpers()
  parsed <- helpers$parse_reference_identity(
    c("CA1microgliares_CA1microgliacon", "DGmicrogliasus_DGmicrogliares"),
    "phenotype_within_unit",
    c("CA1_microglia", "DG_microglia")
  )
  testthat::expect_equal(parsed$anatomical_region, c("CA1", "DG"))
  testthat::expect_equal(parsed$formal_group_contrast, c("RES-vs-CON", "SUS-vs-RES"))
  testthat::expect_true(all(parsed$reference_identity_status == "valid"))

  invalid <- helpers$parse_reference_identity(
    "unstructured_comparison", "phenotype_within_unit", "CA1_microglia"
  )
  testthat::expect_equal(invalid$reference_identity_status, "invalid")
})

testthat::test_that("reference matching rejects cross-region and cross-contrast borrowing", {
  helpers <- load_neuropil_annotation_helpers()
  microglia <- make_reference_term(
    "microglia", "CA1microgliasus_CA1microgliares", "CA1_microglia",
    "CA1", "SUS-vs-RES", 0.01, 2, c("GeneA", "GeneB")
  )
  cross_region <- make_reference_term(
    "neuron_neuropil", "CA3srsus_CA3srres", "CA3_sr",
    "CA3", "SUS-vs-RES", 0.001, 2, c("GeneA", "GeneB")
  )
  cross_contrast <- make_reference_term(
    "neuron_neuropil", "CA1srres_CA1srcon", "CA1_sr",
    "CA1", "RES-vs-CON", 0.001, 2, c("GeneA", "GeneB")
  )

  region_result <- helpers$find_best_neuropil_match(microglia, cross_region)
  contrast_result <- helpers$find_best_neuropil_match(microglia, cross_contrast)
  testthat::expect_false(region_result$neuropil_term_present)
  testthat::expect_false(contrast_result$neuropil_term_present)
  testthat::expect_equal(region_result$n_matched_neuropil_comparisons, 0L)
  testthat::expect_equal(contrast_result$n_matched_neuropil_comparisons, 0L)
})

testthat::test_that("same-region same-contrast neuropil layers are retained transparently", {
  helpers <- load_neuropil_annotation_helpers()
  microglia <- make_reference_term(
    "microglia", "CA1microgliasus_CA1microgliares", "CA1_microglia",
    "CA1", "SUS-vs-RES", 0.01, 2, c("GeneA", "GeneB", "GeneC")
  )
  neuropil <- dplyr::bind_rows(
    make_reference_term(
      "neuron_neuropil", "CA1slmsus_CA1slmres", "CA1_slm",
      "CA1", "SUS-vs-RES", 0.03, 1, c("GeneA", "GeneB")
    ),
    make_reference_term(
      "neuron_neuropil", "CA1srsus_CA1srres", "CA1_sr",
      "CA1", "SUS-vs-RES", 0.20, 2, c("GeneA")
    )
  )
  result <- helpers$find_best_neuropil_match(microglia, neuropil)
  testthat::expect_equal(result$n_matched_neuropil_comparisons, 2L)
  testthat::expect_equal(result$n_matched_neuropil_layers, 2L)
  testthat::expect_equal(result$n_significant_neuropil_matches, 1L)
  testthat::expect_true(result$any_significant_neuropil_match)
  testthat::expect_equal(result$matched_neuropil_layers, "CA1_slm;CA1_sr")
})

testthat::test_that("term presence without neuropil FDR support cannot produce a neuropil-supported claim", {
  helpers <- load_neuropil_annotation_helpers()
  genes <- c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE")
  microglia <- make_reference_term(
    "microglia", "CA1microgliares_CA1microgliacon", "CA1_microglia",
    "CA1", "RES-vs-CON", 0.01, 2, genes
  )
  neuropil <- make_reference_term(
    "neuron_neuropil", "CA1srres_CA1srcon", "CA1_sr",
    "CA1", "RES-vs-CON", 0.20, 2, genes
  )
  result <- annotate_fixture(helpers, microglia, neuropil)
  testthat::expect_true(result$neuropil_term_present)
  testthat::expect_false(result$any_significant_neuropil_match)
  testthat::expect_false(result$interpretation_class %in% c("mixed_microenvironment", "neuropil_sensitive"))
})

testthat::test_that("nonsignificant microglia terms always receive the non-claim class", {
  helpers <- load_neuropil_annotation_helpers()
  microglia <- make_reference_term(
    "microglia", "CA1microgliares_CA1microgliacon", "CA1_microglia",
    "CA1", "RES-vs-CON", 0.20, 2, c("Syn1", "Syp", "GeneA", "GeneB")
  )
  neuropil <- make_reference_term(
    "neuron_neuropil", "CA1srres_CA1srcon", "CA1_sr",
    "CA1", "RES-vs-CON", 0.001, 2, c("Syn1", "Syp", "GeneA", "GeneB")
  )
  result <- annotate_fixture(helpers, microglia, neuropil)
  testthat::expect_false(result$microglia_significant)
  testthat::expect_equal(result$interpretation_class, "not_evaluated_non_significant")
})

testthat::test_that("significant matched overlap supports mixed microenvironment classification", {
  helpers <- load_neuropil_annotation_helpers()
  microglia <- make_reference_term(
    "microglia", "CA1microgliares_CA1microgliacon", "CA1_microglia",
    "CA1", "RES-vs-CON", 0.01, 2, c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE")
  )
  neuropil <- make_reference_term(
    "neuron_neuropil", "CA1srres_CA1srcon", "CA1_sr",
    "CA1", "RES-vs-CON", 0.01, 1, c("GeneA", "GeneB", "Other1", "Other2")
  )
  result <- annotate_fixture(helpers, microglia, neuropil)
  testthat::expect_equal(result$max_significant_gene_overlap_fraction, 0.4)
  testthat::expect_equal(result$interpretation_class, "mixed_microenvironment")
})

testthat::test_that("neuropil-sensitive classification retains marker overlap and direction requirements", {
  helpers <- load_neuropil_annotation_helpers()
  microglia <- make_reference_term(
    "microglia", "CA1microgliares_CA1microgliacon", "CA1_microglia",
    "CA1", "RES-vs-CON", 0.01, 2, c("Syn1", "Syp", "GeneA", "GeneB")
  )
  same_direction <- make_reference_term(
    "neuron_neuropil", "CA1srres_CA1srcon", "CA1_sr",
    "CA1", "RES-vs-CON", 0.01, 1, c("Syn1", "Syp", "Other1", "Other2")
  )
  opposite_direction <- same_direction
  opposite_direction$score <- -1

  supported <- annotate_fixture(helpers, microglia, same_direction)
  direction_failed <- annotate_fixture(helpers, microglia, opposite_direction)
  testthat::expect_equal(supported$interpretation_class, "neuropil_sensitive")
  testthat::expect_false(direction_failed$interpretation_class == "neuropil_sensitive")
  testthat::expect_equal(direction_failed$interpretation_class, "mixed_microenvironment")
})

testthat::test_that("neuropil annotation is invariant to reference input row order", {
  helpers <- load_neuropil_annotation_helpers()
  microglia <- make_reference_term(
    "microglia", "CA1microgliares_CA1microgliacon", "CA1_microglia",
    "CA1", "RES-vs-CON", 0.01, 2, c("GeneA", "GeneB", "GeneC")
  )
  neuropil <- dplyr::bind_rows(
    make_reference_term(
      "neuron_neuropil", "CA1srres_CA1srcon", "CA1_sr",
      "CA1", "RES-vs-CON", 0.01, 1, c("GeneA", "GeneB")
    ),
    make_reference_term(
      "neuron_neuropil", "CA1slmres_CA1slmcon", "CA1_slm",
      "CA1", "RES-vs-CON", 0.02, 1, c("GeneA")
    )
  )
  forward <- annotate_fixture(helpers, microglia, neuropil)
  reversed <- annotate_fixture(helpers, microglia, neuropil[2:1, , drop = FALSE])
  testthat::expect_equal(forward, reversed)
})

testthat::test_that("program integration is program-specific and excludes nonsignificant majority votes", {
  helpers <- load_neuropil_program_helpers()
  route <- list(
    dataset = "microglia",
    comparison = "CA1microgliares_CA1microgliacon",
    route_category = "phenotype_within_unit",
    route_unit = "CA1_microglia"
  )
  annotation <- tibble::tibble(
    dataset = route$dataset,
    comparison = route$comparison,
    route_category = route$route_category,
    route_unit = route$route_unit,
    term_description = c(
      "synaptic organization",
      "ribosomal translation",
      rep("ribosome biogenesis", 20)
    ),
    p_adjust = c(0.01, 0.02, rep(0.5, 20)),
    interpretation_class = c(
      "mixed_microenvironment",
      "microglia_robust",
      rep("mixed_microenvironment", 20)
    ),
    classification_contract_version = "test_contract",
    gene_overlap_fraction = 0.25,
    gene_jaccard = 0.15,
    neuropil_marker_fraction = 0.10,
    microglia_marker_fraction = 0.20
  )
  program_summary <- tibble::tibble(
    dataset = route$dataset,
    comparison = route$comparison,
    route_category = route$route_category,
    route_unit = route$route_unit,
    biological_program = c(
      "Synapse_Vesicle_Organization",
      "Ribosome_Translation",
      "Oxidative_Redox_Stress"
    )
  )

  result <- helpers$attach_neuropil_program_annotation(program_summary, annotation)
  synapse <- result %>% dplyr::filter(.data$biological_program == "Synapse_Vesicle_Organization")
  ribosome <- result %>% dplyr::filter(.data$biological_program == "Ribosome_Translation")
  redox <- result %>% dplyr::filter(.data$biological_program == "Oxidative_Redox_Stress")

  testthat::expect_equal(synapse$interpretation_class, "mixed_microenvironment")
  testthat::expect_equal(synapse$n_significant_microglia_terms, 1L)
  testthat::expect_equal(ribosome$interpretation_class, "microglia_robust")
  testthat::expect_equal(ribosome$n_significant_microglia_terms, 1L)
  testthat::expect_equal(ribosome$n_mixed_microenvironment, 0L)
  testthat::expect_equal(redox$interpretation_class, "no_significant_microglia_terms")
  testthat::expect_equal(redox$neuropil_program_annotation_status, "no_significant_microglia_terms")
})

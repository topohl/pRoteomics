testthat::local_edition(3)

testthat::skip_if_not_installed("GO.db")
testthat::skip_if_not_installed("AnnotationDbi")

source(testthat::test_path("..", "..", "R", "manuscript_go_theme_utils.R"))
source(testthat::test_path("..", "..", "R", "gsea_wgcna_concordance_utils.R"))

ontology_fixture <- function(ids, descriptions) {
  data.frame(
    dataset = "neuron_neuropil",
    phenotype_contrast = "SUS_vs_RES",
    spatial_unit = "CA3-SR",
    program_class = "Mitochondria_OXPHOS_Metabolism",
    Comparison = "fixture",
    ID = ids,
    Description = descriptions,
    NES = -1.5,
    pvalue = 0.001,
    p.adjust = 0.01,
    core_enrichment = "P1/P2",
    core_enrichment_gene = "G1/G2",
    evidence_source_family = "ranked_GSEA",
    source_supplementary_file = "fixture.csv",
    stringsAsFactors = FALSE
  )
}

testthat::test_that("valid GO IDs use ontology assignments, never legacy regex", {
  registry <- read_manuscript_go_theme_registry(testthat::test_path(
    "..", "..", "config", "manuscript_go_theme_registry.tsv"
  ))
  x <- ontology_fixture(
    c("GO:0006119", "GO:0006029"),
    c("oxidative phosphorylation", "proteoglycan metabolic process")
  )
  out <- gww_build_ontology_theme_term_table(x, registry)
  oxphos <- out[out$GO_ID == "GO:0006119", , drop = FALSE]
  proteoglycan <- out[out$GO_ID == "GO:0006029", , drop = FALSE]
  testthat::expect_true(any(
    oxphos$theme_id == "mitochondrial_respiration_oxphos"
  ))
  testthat::expect_true(all(
    oxphos$mapping_method == "go_id_ontology", na.rm = TRUE
  ))
  testthat::expect_true(all(is.na(proteoglycan$theme_id)))
  testthat::expect_identical(
    unique(proteoglycan$legacy_program_class),
    "Mitochondria_OXPHOS_Metabolism"
  )
})

testthat::test_that("QC-review ontology branches cannot become claim evidence", {
  registry <- read_manuscript_go_theme_registry(testthat::test_path(
    "..", "..", "config", "manuscript_go_theme_registry.tsv"
  ))
  x <- ontology_fixture("GO:0008544", "epidermis development")
  out <- gww_build_ontology_theme_term_table(x, registry)
  testthat::expect_true(any(out$theme_role == "qc_review"))
  testthat::expect_false(any(out$theme_claim_eligible %in% TRUE))
  testthat::expect_equal(nrow(gww_build_local_gsea_evidence(out)), 0L)
})

testthat::test_that("checked-in theme mappings are composite-key unique", {
  mapping <- gww_read_theme_entity_mapping(testthat::test_path(
    "..", "..", "config", "gsea_wgcna_theme_module_mapping.csv"
  ))
  key <- c("dataset", "theme_id", "entity_level", "entity_id")
  testthat::expect_identical(anyDuplicated(mapping[key]), 0L)
  testthat::expect_true(all(
    mapping$mapping_role %in% c("primary", "secondary_multifunctional")
  ))
  testthat::expect_false(any(
    mapping$theme_id %in% c(
      "cytoskeleton_structure", "epithelial_epidermal_qc"
    ) &
      mapping$approved_for_manuscript_interpretation %in% TRUE
  ))
})

testthat::test_that("resolved multifaceted mappings retain one WGCNA entity", {
  mapping <- gww_read_theme_entity_mapping(testthat::test_path(
    "..", "..", "config", "gsea_wgcna_theme_module_mapping.csv"
  ))
  approved_secondary <- mapping |>
    dplyr::filter(
      .data$mapping_status == "approved",
      .data$mapping_role == "secondary_multifunctional"
    )
  expected <- data.frame(
    dataset = c("neuron_neuropil", "neuron_neuropil", "microglia"),
    theme_id = c(
      "autophagy_lysosome_endosome", "chromatin_organization",
      "chromatin_organization"
    ),
    entity_id = c("WGCNA_m05", "WGCNA_m12", "WGCNA_m11"),
    stringsAsFactors = FALSE
  )
  testthat::expect_equal(
    nrow(dplyr::anti_join(
      expected, approved_secondary,
      by = c("dataset", "theme_id", "entity_id")
    )), 0L
  )
  testthat::expect_false(any(
    mapping$dataset == "neuron_soma" &
      mapping$theme_id == "ribosome_translation" &
      mapping$entity_id == "WGCNA_m04" &
      mapping$approved_for_manuscript_interpretation %in% TRUE
  ))
  testthat::expect_false(any(
    mapping$dataset == "microglia" &
      mapping$theme_id == "autophagy_lysosome_endosome" &
      mapping$entity_id == "WGCNA_m08" &
      mapping$approved_for_manuscript_interpretation %in% TRUE
  ))
})

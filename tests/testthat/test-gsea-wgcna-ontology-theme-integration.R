source(testthat::test_path("..", "..", "R", "paths.R"))

testthat::local_edition(3)

testthat::skip_if_not_installed("GO.db")
testthat::skip_if_not_installed("AnnotationDbi")

source(repo_path("R", "manuscript_go_theme_utils.R"))
source(repo_path("R", "gsea_wgcna_concordance_utils.R"))

# The checked-in mapping, read as the pipeline reads it: against the current
# supermodule membership. It names the primary network's entities, so a tree
# whose Stage 07 manifest declares provisional labels (the sensitivity tree)
# skips these checks (helper-wgcna-current-artifacts.R) and tests the lineage
# proposal it reads instead (below).
checked_in_theme_mapping <- function() {
  decl <- wgcna_declared_label_mode("microglia")
  testthat::skip_if(identical(decl$mode, "provisional"), paste0(
    "Declared provisional label mode (", decl$declared_by,
    "): config/gsea_wgcna_theme_module_mapping.csv maps the primary network and is not read here."))
  contracts <- gww_membership_contract_paths()
  testthat::skip_if_not(all(file.exists(contracts)), "identity contracts are not built in this tree")
  gww_read_theme_entity_mapping(
    testthat::test_path("..", "..", "config", "gsea_wgcna_theme_module_mapping.csv"),
    gww_read_supermodule_membership(contracts)
  )
}

# Fixtures for the membership guard: microglia SM03 = m04 + m08 + m12, SM08 = m07 + m11.
membership_fixture <- data.frame(
  dataset = "microglia",
  module_id = c("WGCNA_m04", "WGCNA_m08", "WGCNA_m12", "WGCNA_m07", "WGCNA_m11"),
  supermodule_id = c("SM03", "SM03", "SM03", "SM08", "SM08"),
  stringsAsFactors = FALSE
)
mapping_fixture <- function(entity_id, members) {
  path <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(
    dataset = "microglia", theme_id = "synaptic_signaling_vesicle",
    entity_level = c("module", "supermodule"), entity_id = c("WGCNA_m04", entity_id),
    mapping_status = "approved", mapping_rationale = "fixture", annotation_source = "fixture",
    confidence = "high", approved_for_manuscript_interpretation = TRUE, mapping_role = "primary",
    member_module_ids = c(NA, members), stringsAsFactors = FALSE
  ), path, row.names = FALSE, na = "")
  path
}

ontology_fixture <- function(ids, descriptions) {
  data.frame(
    dataset = "neuron_neuropil",
    phenotype_contrast = "SUS_vs_RES",
    spatial_unit = "CA3-SR",
    program_class = "Mitochondria_OXPHOS_Metabolism",
    comparison = "fixture",
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
  mapping <- checked_in_theme_mapping()
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
  mapping <- checked_in_theme_mapping()
  approved_secondary <- mapping |>
    dplyr::filter(
      .data$mapping_status == "approved",
      .data$mapping_role == "secondary_multifunctional"
    )
  # eight-animal network (mapping of 2026-11-23): the secondary facets carried by lineage
  expected <- data.frame(
    dataset = c("neuron_neuropil", "neuron_neuropil", "microglia"),
    theme_id = c(
      "autophagy_lysosome_endosome", "chromatin_organization",
      "chromatin_organization"
    ),
    entity_id = c("WGCNA_m04", "WGCNA_m12", "WGCNA_m10"),
    stringsAsFactors = FALSE
  )
  testthat::expect_equal(
    nrow(dplyr::anti_join(
      expected, approved_secondary,
      by = c("dataset", "theme_id", "entity_id")
    )), 0L
  )
  # facets left for review stay out of the manuscript
  testthat::expect_false(any(
    mapping$dataset == "neuron_soma" &
      mapping$theme_id == "chromatin_organization" &
      mapping$entity_id == "WGCNA_m02" &
      mapping$approved_for_manuscript_interpretation %in% TRUE
  ))
  testthat::expect_false(any(
    mapping$dataset == "neuron_neuropil" &
      mapping$theme_id == "autophagy_lysosome_endosome" &
      mapping$entity_id == "WGCNA_m01" &
      mapping$approved_for_manuscript_interpretation %in% TRUE
  ))
  testthat::expect_true(all(
    mapping$approved_for_manuscript_interpretation[mapping$mapping_status == "manual_review_candidate"] %in% FALSE
  ))
})

# ---- supermodule membership guard ------------------------------------------

testthat::test_that("checked-in supermodule rows name the current members", {
  mapping <- checked_in_theme_mapping()
  sm <- mapping[mapping$entity_level == "supermodule", , drop = FALSE]
  testthat::expect_gt(nrow(sm), 0L)
  testthat::expect_true(all(nzchar(sm$member_module_ids)))
  testthat::expect_true(all(is.na(mapping$member_module_ids[mapping$entity_level == "module"])))
})

# A provisional tree reads the lineage proposal (PROTEOMICS_GSEA_WGCNA_THEME_MAPPING, by default the
# output of tools/map_wgcna_module_lineage.R); it must pass the same guard, or concordance stops.
testthat::test_that("the lineage mapping a provisional tree reads passes the membership guard", {
  decl <- wgcna_declared_label_mode("microglia")
  testthat::skip_if_not(identical(decl$mode, "provisional"),
                        "Reviewed label mode: the checked-in mapping is read and tested above.")
  path <- Sys.getenv("PROTEOMICS_GSEA_WGCNA_THEME_MAPPING", unset = repo_path(
    "results", "wgcna", "map_wgcna_module_lineage", "global", "tables", "provisional_theme_module_mapping.csv"))
  testthat::skip_if_not(file.exists(path), "tools/map_wgcna_module_lineage.R has not been run in this tree")
  contracts <- gww_membership_contract_paths()
  testthat::skip_if_not(all(file.exists(contracts)), "identity contracts are not built in this tree")
  mapping <- gww_read_theme_entity_mapping(path, gww_read_supermodule_membership(contracts))
  testthat::expect_gt(nrow(mapping), 0L)
  testthat::expect_true(all(mapping$mapping_status == "provisional_supported"))
  testthat::expect_false(any(mapping$approved_for_manuscript_interpretation %in% TRUE))
})

testthat::test_that("an unchanged member set listed in another order is accepted", {
  x <- gww_read_theme_entity_mapping(mapping_fixture("SM03", "WGCNA_m12; WGCNA_m04;WGCNA_m08"),
                                     membership_fixture)
  testthat::expect_equal(nrow(x), 2L)
})

testthat::test_that("a surviving supermodule ID whose members changed is refused", {
  testthat::expect_error(
    gww_read_theme_entity_mapping(mapping_fixture("SM03", "WGCNA_m04;WGCNA_m08"), membership_fixture),
    "membership changed since review")
  # the renumbering case: an old ID that now names a different block
  testthat::expect_error(
    gww_read_theme_entity_mapping(mapping_fixture("SM08", "WGCNA_m09;WGCNA_m13"), membership_fixture),
    "membership changed since review")
})

testthat::test_that("a supermodule ID absent from the current membership is refused", {
  testthat::expect_error(
    gww_read_theme_entity_mapping(mapping_fixture("SM09", "WGCNA_m07;WGCNA_m11"), membership_fixture),
    "absent from the current membership")
})

testthat::test_that("supermodule rows must state members and the reader needs the membership", {
  testthat::expect_error(
    gww_read_theme_entity_mapping(mapping_fixture("SM03", NA), membership_fixture),
    "need member_module_ids")
  testthat::expect_error(
    gww_read_theme_entity_mapping(mapping_fixture("SM03", "WGCNA_m04;WGCNA_m08;WGCNA_m12")),
    "current supermodule membership")
})

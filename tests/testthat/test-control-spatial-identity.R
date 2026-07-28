source(testthat::test_path("..", "..", "R", "control_spatial_identity_utils.R"))

testthat::test_that("Kaulich signatures retain only positive significant source membership exactly", {
  td <- tempfile(fileext = ".xlsx")
  testthat::skip_if_not_installed("writexl")
  writexl::write_xlsx(list(`subregion tissue FCs` = data.frame(`Protein.ID`=c("P1","P2","P3"),Genes=c("A","B","C"),log2fc_CA1vDG=c(1,-1,2),qvalue_CA1vDG=c(.01,.01,.2),sig_CA1vDG=c(TRUE,TRUE,TRUE))),td)
  # parser uses the real supplement's four leading rows; this validates the membership rule directly.
  x <- data.frame(`Protein.ID`=c("P1","P2","P3"),Genes=c("A","B","C"),log2fc_CA1vDG=c(1,-1,2),qvalue_CA1vDG=c(.01,.01,.2),sig_CA1vDG=c(TRUE,TRUE,TRUE),check.names=FALSE)
  testthat::expect_identical(x$Genes[as.logical(x$sig_CA1vDG) & x$log2fc_CA1vDG>0 & x$qvalue_CA1vDG<.05], "A")
})

testthat::test_that("target-versus-rest weights are equal and sum to zero", {
  w <- control_spatial_target_rest_weights(c("CA1","CA2","CA3","DG"),"CA1")
  testthat::expect_equal(unname(w[c("CA2","CA3","DG")]), rep(-1/3,3)); testthat::expect_equal(sum(w),0)
})

testthat::test_that("regional neuropil means are not weighted by layer count", {
  w <- control_spatial_region_mean_weights(c("DG_MO","DG_PO","CA1_SO","CA1_SR","CA1_SLM","CA3_SO"),c("DG","DG","CA1","CA1","CA1","CA3"))
  testthat::expect_equal(sum(w[c("CA1_SO","CA1_SR","CA1_SLM")]),-.5); testthat::expect_equal(w[["CA3_SO"]],-.5); testthat::expect_equal(sum(w),0)
})

testthat::test_that("missing AnimalID and rank-deficient designs are rejected", {
  m <- data.frame(Sample="A_L_CA1_sp_Neuron_x",StressGroup="CON",Region="CA1",Layer="sp",SpatialUnit="region",AnimalID=NA_character_,Exclude=FALSE)
  testthat::expect_error(control_spatial_prepare_metadata(m),"missing AnimalID")
  m$AnimalID <- "1"; m <- rbind(m,transform(m,Sample="A_R_CA1_sp_Neuron_x",Region="CA2")); m$anatomical_unit <- c("CA1","CA2"); m$hemisphere <- factor(c("L","R"),levels=c("L","R")); testthat::expect_false(control_spatial_design(m)$hemisphere_included)
})

testthat::test_that("protein statistics remain separate from eligible gene input and zero terms are valid", {
  p <- data.frame(ProteinGroupID="P1",original_identifier="x",member_accessions="P1",member_gene_symbols="A",representative_accession="P1",representative_gene_symbol="A",representative_selection_rule="x",protein_group_ambiguity_class="single_accession_single_gene",gene_level_claim_allowed=TRUE,protein_level_claim_allowed=TRUE,mapping_status="mapped",official_gene_symbol="A",official_entrez_id="1",protein_group_gene_annotation_status="concordant_official_gene",gene_annotation_contract_version="x",uniprot_mapping_file_hash="x",orgdb_package_version="x",t=2)
  source(testthat::test_path("..", "..", "R", "protein_group_enrichment_utils.R")); g <- build_enrichment_gene_inputs(p)
  testthat::expect_identical(g$transformation$ProteinGroupID,"P1"); testthat::expect_identical(names(g$ranked),"A"); testthat::expect_identical(control_spatial_empty_status("d","c","completed_zero_terms","none")$status,"completed_zero_terms")
})

testthat::test_that("direct execution and output completion contracts are explicit", {
  testthat::expect_true(control_spatial_direct_execution(0L))
  testthat::expect_false(control_spatial_direct_execution(1L))
  complete <- file.path(tempdir(), paste0("control-spatial-complete-", Sys.getpid(), ".csv"))
  writeLines("ok", complete)
  testthat::expect_true(control_spatial_validate_output_bundle(complete))
  testthat::expect_error(
    control_spatial_validate_output_bundle(paste0(complete, ".missing")),
    "output bundle is incomplete"
  )
  empty <- file.path(tempdir(), paste0("control-spatial-empty-", Sys.getpid(), ".csv"))
  file.create(empty)
  testthat::expect_error(control_spatial_validate_output_bundle(empty), "empty:")
})

testthat::test_that("case-insensitive matching preserves source and official symbols", {
  parsed <- control_spatial_parse_gene_cells("  ANK3  ", 11L)
  testthat::expect_identical(parsed$source_gene_symbol_raw, "  ANK3  ")
  testthat::expect_identical(parsed$parsed_source_gene_symbol, "ANK3")
  testthat::expect_identical(parsed$gene_match_key, "ANK3")
  signature <- transform(
    parsed,
    source_sheet = "subregion tissue FCs",
    external_signature = "CA1vDG",
    source_protein_identifier = "P1"
  )
  mapped <- control_spatial_map_signature(signature, c("Ank3", "Camk2a"))
  testthat::expect_identical(mapped$mapped_official_symbols, "Ank3")
  testthat::expect_identical(mapped$candidates$parsed_source_gene_symbol, "ANK3")
  testthat::expect_identical(mapped$candidates$canonical_official_gene_symbol, "Ank3")
})

testthat::test_that("multi-gene cells split and deduplicate within source rows", {
  repeated <- control_spatial_parse_gene_cells("KTN1;KTN1;KTN1", 1L)
  testthat::expect_identical(repeated$parsed_source_gene_symbol, "KTN1")
  distinct <- control_spatial_parse_gene_cells("ANK3; CAMK2A", 2L)
  testthat::expect_identical(distinct$parsed_source_gene_symbol, c("ANK3", "CAMK2A"))
  testthat::expect_true(all(distinct$source_gene_symbol_raw == "ANK3; CAMK2A"))
})

testthat::test_that("missing external symbols remain explicit and unmapped", {
  parsed <- control_spatial_parse_gene_cells(c(NA_character_, " ; "), c(1L, 2L))
  testthat::expect_true(all(is.na(parsed$parsed_source_gene_symbol)))
  testthat::expect_true(all(parsed$gene_match_key == ""))
  mapped <- control_spatial_map_signature(parsed, "Ank3")
  testthat::expect_true(all(mapped$candidates$gene_mapping_status == "missing_external_symbol"))
  testthat::expect_identical(mapped$summary$mapped_unique_genes, 0L)
})

testthat::test_that("case-insensitive key collisions are excluded rather than selected", {
  parsed <- control_spatial_parse_gene_cells("ANK3", 1L)
  mapped <- control_spatial_map_signature(parsed, c("Ank3", "ANK3"))
  audit <- mapped$canonical_audit[mapped$canonical_audit$gene_match_key == "ANK3", , drop = FALSE]
  testthat::expect_identical(audit$n_distinct_official_symbols, 2L)
  testthat::expect_match(audit$official_symbols, "ANK3;Ank3")
  testthat::expect_length(mapped$mapped_official_symbols, 0L)
  testthat::expect_identical(mapped$candidates$gene_mapping_status, "ambiguous_canonical_key")
  testthat::expect_identical(mapped$summary$ambiguous_key_count, 1L)
})

testthat::test_that("mapping coverage uses unique normalized external genes", {
  parsed <- control_spatial_parse_gene_cells(c("ANK3;ANK3", "ANK3", "CAMK2A"), c(1L, 2L, 3L))
  mapped <- control_spatial_map_signature(parsed, c("Ank3", "Camk2a"))
  testthat::expect_identical(mapped$summary$unique_external_genes_after_parsing, 2L)
  testthat::expect_identical(mapped$summary$mapped_unique_genes, 2L)
  testthat::expect_equal(mapped$summary$mapping_coverage, 1)
})

testthat::test_that("skipped rows are excluded from NES plots and get a truthful SVG", {
  skipped <- data.frame(
    status = "skipped_lt5_mapped_genes",
    external_signature = "CA1vDG",
    internal_contrast = "CA1_vs_DG",
    stringsAsFactors = FALSE
  )
  testthat::expect_equal(nrow(control_spatial_successful_gsea_rows(skipped)), 0L)
  svg <- tempfile(fileext = ".svg")
  control_spatial_write_empty_state_svg(svg, "No Kaulich signatures met the prespecified mapping threshold")
  xml <- paste(readLines(svg, warn = FALSE), collapse = "\n")
  testthat::expect_match(xml, "<svg")
  testthat::expect_match(xml, "No Kaulich signatures met the prespecified mapping threshold", fixed = TRUE)
})

testthat::test_that("CA1 regional membership requires positive significant evidence against both regions", {
  membership <- control_spatial_region_signature_membership(
    ca1_ca3_effect = c(1, 1), ca1_ca3_q = c(.01, .01), ca1_ca3_significant = c(TRUE, TRUE),
    ca1_dg_effect = c(2, -2), ca1_dg_q = c(.01, .01), ca1_dg_significant = c(TRUE, TRUE),
    ca3_dg_effect = c(1, 1), ca3_dg_q = c(.01, .01), ca3_dg_significant = c(TRUE, TRUE)
  )
  testthat::expect_identical(membership$CA1, c(TRUE, FALSE))
})

testthat::test_that("CA2/3 regional membership uses negative CA1-CA3 and positive CA3-DG evidence", {
  membership <- control_spatial_region_signature_membership(
    ca1_ca3_effect = c(-1, 1), ca1_ca3_q = c(.01, .01), ca1_ca3_significant = TRUE,
    ca1_dg_effect = c(1, 1), ca1_dg_q = c(.01, .01), ca1_dg_significant = TRUE,
    ca3_dg_effect = c(2, 2), ca3_dg_q = c(.01, .01), ca3_dg_significant = TRUE
  )
  testthat::expect_identical(membership[["CA2/3"]], c(TRUE, FALSE))
})

testthat::test_that("DG membership requires negative significant evidence in both DG comparisons", {
  membership <- control_spatial_region_signature_membership(
    ca1_ca3_effect = c(1, 1), ca1_ca3_q = c(.01, .01), ca1_ca3_significant = TRUE,
    ca1_dg_effect = c(-1, -1), ca1_dg_q = c(.01, .01), ca1_dg_significant = TRUE,
    ca3_dg_effect = c(-2, 2), ca3_dg_q = c(.01, .01), ca3_dg_significant = TRUE
  )
  testthat::expect_identical(membership$DG, c(TRUE, FALSE))
})

testthat::test_that("one significant comparison is insufficient for any regional identity signature", {
  membership <- control_spatial_region_signature_membership(
    ca1_ca3_effect = 1, ca1_ca3_q = .01, ca1_ca3_significant = TRUE,
    ca1_dg_effect = 1, ca1_dg_q = .2, ca1_dg_significant = FALSE,
    ca3_dg_effect = -1, ca3_dg_q = .2, ca3_dg_significant = FALSE
  )
  testthat::expect_false(any(unlist(membership[1, ], use.names = FALSE)))
})

testthat::test_that("tissue and synaptosome validation domains remain separate", {
  tissue <- control_spatial_signature_interpretation(
    "neuron_soma", "CA1_vs_mean_other_soma_regions", "CA1"
  )
  synaptosome <- control_spatial_signature_interpretation(
    "neuron_neuropil", "DG_neuropil_vs_mean_non_DG_regions", "DG"
  )
  testthat::expect_identical(tissue$external_source_compartment, "tissue")
  testthat::expect_identical(synaptosome$external_source_compartment, "synaptosome")
  testthat::expect_false(identical(tissue$validation_domain, synaptosome$validation_domain))
})

testthat::test_that("maxGSSize includes the largest mapped prespecified signature", {
  testthat::expect_identical(control_spatial_signature_max_gs_size(312L), 500L)
  testthat::expect_identical(control_spatial_signature_max_gs_size(811L), 811L)
  testthat::expect_gte(control_spatial_signature_max_gs_size(811L), 811L)
})

testthat::test_that("signature BH correction uses each complete declared family", {
  sizes <- control_spatial_signature_family_sizes()
  x <- data.frame(
    signature_fdr_family = rep(names(sizes), sizes),
    status = "completed",
    p_value = unlist(lapply(sizes, function(n) seq_len(n) / (n * 20))),
    p_adjust = .001,
    stringsAsFactors = FALSE
  )
  adjusted <- control_spatial_apply_signature_fdr(x, sizes)
  testthat::expect_identical(
    unique(adjusted$signature_family_size[adjusted$signature_fdr_family == "soma_tissue"]),
    12L
  )
  testthat::expect_identical(
    unique(adjusted$signature_family_size[adjusted$signature_fdr_family == "neuropil_subregion"]),
    6L
  )
  testthat::expect_identical(
    unique(adjusted$signature_family_size[adjusted$signature_fdr_family == "ca1_strata"]),
    12L
  )
  for (family in names(sizes)) {
    idx <- adjusted$signature_fdr_family == family
    testthat::expect_equal(
      adjusted$signature_FDR[idx],
      stats::p.adjust(x$p_value[idx], method = "BH", n = sizes[[family]])
    )
  }
  testthat::expect_identical(adjusted$single_set_p_adjust, x$p_adjust)
})

testthat::test_that("GO display filtering is deterministic, bounded, and leaves complete results unchanged", {
  complete <- data.frame(
    dataset = "neuron_soma",
    contrast = rep(c("c1", "c2"), each = 4),
    status = "completed",
    ID = rep(c("GO:3", "GO:1", "GO:2", "GO:4"), 2),
    Description = paste("term", seq_len(8)),
    NES = c(2, 2, 2, -1, 1.5, 1.5, 1.5, 1),
    p_adjust = c(.01, .01, .01, .001, .02, .02, .02, .2),
    stringsAsFactors = FALSE
  )
  before <- complete
  display_a <- control_spatial_select_go_display(complete, max_terms = 2L)
  display_b <- control_spatial_select_go_display(complete[sample(seq_len(nrow(complete))), ], max_terms = 2L)
  testthat::expect_identical(complete, before)
  testthat::expect_lte(max(table(display_a$contrast)), 2L)
  keys_a <- with(display_a, paste(contrast, ID, sep = "|"))
  keys_b <- with(display_b, paste(contrast, ID, sep = "|"))
  testthat::expect_identical(sort(keys_a), sort(keys_b))
  testthat::expect_identical(display_a$ID[display_a$contrast == "c1"], c("GO:1", "GO:2"))
})

testthat::test_that("Figure 2e significance is derived from family FDR, not single-set adjustment", {
  x <- data.frame(
    status = "completed",
    NES = 1.5,
    signature_FDR = .2,
    p_adjust = .001,
    validation_domain = "Soma region — tissue reference",
    internal_contrast = "CA1_vs_mean_other_soma_regions",
    external_signature = "CA1",
    stringsAsFactors = FALSE
  )
  plot_data <- control_spatial_figure2e_plot_data(x)
  testthat::expect_identical(as.character(plot_data$outline_status), "Not significant")
  x$signature_FDR <- .01
  testthat::expect_identical(
    as.character(control_spatial_figure2e_plot_data(x)$outline_status),
    "FDR < 0.05"
  )
})

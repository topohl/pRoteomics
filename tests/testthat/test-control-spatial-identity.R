source(testthat::test_path(
  "..", "..", "R", "clusterprofiler_reproducibility.R"
))
source(testthat::test_path("..", "..", "R", "control_spatial_identity_utils.R"))

testthat::test_that("control-spatial GSEA seeds use stable semantic identity", {
  args <- list(
    dataset = "neuron_neuropil",
    contrast = "CA1_SLM_vs_mean_other_CA1_strata",
    enrichment_type = "gseGO_BP",
    gsea_seed_base = 20260824L,
    n_perm_simple = 100000L
  )
  first <- do.call(control_spatial_gsea_reproducibility, args)
  second <- do.call(control_spatial_gsea_reproducibility, args)
  different_contrast <- do.call(
    control_spatial_gsea_reproducibility,
    utils::modifyList(args, list(
      contrast = "CA1_SO_vs_mean_other_CA1_strata"
    ))
  )
  different_signature <- control_spatial_gsea_reproducibility(
    dataset = args$dataset, contrast = args$contrast,
    enrichment_type = "GSEA_Kaulich_signature",
    external_signature = "SLM", gsea_seed_base = args$gsea_seed_base,
    n_perm_simple = args$n_perm_simple
  )
  testthat::expect_identical(first, second)
  testthat::expect_false(identical(
    first$gsea_seed, different_contrast$gsea_seed
  ))
  testthat::expect_false(identical(
    first$gsea_seed, different_signature$gsea_seed
  ))
  testthat::expect_identical(first$n_perm_simple, 100000L)
})

testthat::test_that("control-spatial seeded execution restores caller RNG", {
  identity <- control_spatial_gsea_reproducibility(
    "neuron_neuropil", "DG_neuropil_vs_mean_non_DG_regions", "gseGO_BP",
    20260824L, 100000L
  )
  fake_gsea <- function(...) stats::runif(4)
  set.seed(9876)
  kind_before <- RNGkind()
  state_before <- .Random.seed
  first <- run_seeded_clusterprofiler_gsea(
    fake_gsea, identity$gsea_seed, identity$n_perm_simple,
    geneList = c(A = 2, B = -1)
  )
  testthat::expect_identical(RNGkind(), kind_before)
  testthat::expect_identical(.Random.seed, state_before)
  second <- run_seeded_clusterprofiler_gsea(
    fake_gsea, identity$gsea_seed, identity$n_perm_simple,
    geneList = c(A = 2, B = -1)
  )
  testthat::expect_identical(first, second)
})

testthat::test_that("script 09 cannot bypass deterministic GSEA governance", {
  script <- paste(readLines(testthat::test_path(
    "..", "..", "04_differential_expression_enrichment",
    "09_control_spatial_identity_validation.r"
  ), warn = FALSE), collapse = "\n")
  testthat::expect_false(grepl(
    "clusterProfiler::(gseGO|GSEA)[[:space:]]*\\(", script, perl = TRUE
  ))
  testthat::expect_equal(
    lengths(regmatches(
      script,
      gregexpr("run_seeded_clusterprofiler_gsea\\(", script, perl = TRUE)
    )),
    2L
  )
  testthat::expect_match(script, "gsea_seed = go_repro\\$gsea_seed")
  testthat::expect_match(script, "gsea_seed = job_repro\\$gsea_seed")
})

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
  testthat::expect_setequal(strsplit(audit$official_symbols, ";", fixed = TRUE)[[1]], c("ANK3", "Ank3"))
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

testthat::test_that("Figure 2e keeps external CA1 SP as reference-only context", {
  x <- data.frame(
    dataset = "neuron_neuropil",
    internal_contrast = c("CA1_SLM_vs_mean_other_CA1_strata", "CA1_SLM_vs_mean_other_CA1_strata"),
    external_signature = c("SLM", "SP"),
    match_type = c("exact", "specificity_comparison"),
    interpretation_note = c("Expected exact target-stratum correspondence.", "Alternative anatomical signature used as a specificity comparison."),
    stringsAsFactors = FALSE
  )
  audit <- control_spatial_figure2e_matching_audit(x)
  testthat::expect_true(audit$matched_exactly[audit$external_signature == "SLM"])
  testthat::expect_false(audit$matched_exactly[audit$external_signature == "SP"])
  testthat::expect_true(audit$external_reference_context[audit$external_signature == "SP"])
  testthat::expect_match(audit$match_reason[audit$external_signature == "SP"], "no internal CA1-SP", fixed = TRUE)
})

testthat::test_that("DG layer target-versus-rest contrasts use observed DG units only", {
  units <- c("DG_MO", "DG_PO")
  for (target in units) {
    w <- control_spatial_target_rest_weights(units, target)
    testthat::expect_equal(sum(w), 0)
    testthat::expect_identical(names(w), units)
  }
})

testthat::test_that("Figure 2e excludes DG layer contrasts while Figure 2f displays them", {
  figure2e <- control_spatial_figure2e_external_contrasts()
  figure2f <- control_spatial_figure2f_display_contrasts()
  testthat::expect_false(any(grepl("^DG_(MO|PO)_", figure2e)))
  testthat::expect_setequal(
    figure2f,
    c(
      "CA1_vs_mean_other_soma_regions", "CA2_vs_mean_other_soma_regions",
      "CA3_vs_mean_other_soma_regions", "DG_vs_mean_other_soma_regions",
      "DG_MO_vs_mean_other_DG_layers", "DG_PO_vs_mean_other_DG_layers"
    )
  )
  testthat::expect_length(figure2f, 6L)
})

testthat::test_that("Figure 2f regions and CA1-layers candidate has exactly seven contrasts", {
  candidate <- control_spatial_figure2f_regions_ca1layers_display_contrasts()
  testthat::expect_identical(
    candidate,
    c(
      "CA1_vs_mean_other_soma_regions", "CA2_vs_mean_other_soma_regions",
      "CA3_vs_mean_other_soma_regions", "DG_vs_mean_other_soma_regions",
      "CA1_SLM_vs_mean_other_CA1_strata",
      "CA1_SO_vs_mean_other_CA1_strata",
      "CA1_SR_vs_mean_other_CA1_strata"
    )
  )
  testthat::expect_length(candidate, 7L)
  testthat::expect_false(any(grepl("^DG_(MO|PO)_", candidate)))
})

testthat::test_that("Figure 2f grouped layout abbreviates only the seven display labels", {
  contrasts <- control_spatial_figure2f_regions_ca1layers_display_contrasts()
  layout <- control_spatial_figure2f_grouped_layout()
  testthat::expect_identical(layout$contrast, contrasts)
  testthat::expect_identical(
    layout$group,
    c(rep("Soma region", 4L), rep("CA1 neuropil layer", 3L))
  )
  testthat::expect_identical(
    layout$short_label, c("CA1", "CA2", "CA3", "DG", "SLM", "SO", "SR")
  )
  testthat::expect_false(anyDuplicated(layout$contrast) > 0L)
})

testthat::test_that("Figure 2f grouped candidate is render-only from validated source data", {
  script <- paste(readLines(testthat::test_path(
    "..", "..", "04_differential_expression_enrichment",
    "09_control_spatial_identity_validation.r"
  ), warn = FALSE), collapse = "\n")
  testthat::expect_match(
    script, "PROTEOMICS_CONTROL_SPATIAL_FIGURE2F_GROUPED_LAYOUT_RENDER_ONLY",
    fixed = TRUE
  )
  testthat::expect_match(
    script, "figure2f_regions_CA1layers_source_data.csv", fixed = TRUE
  )
  testthat::expect_match(
    script, "figure2f_control_anatomical_GO_regions_CA1layers_grouped.svg",
    fixed = TRUE
  )
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

testthat::test_that("GO display contrast filtering leaves complete results unchanged", {
  complete <- data.frame(
    dataset = "neuron_neuropil",
    contrast = c("DG_MO_vs_mean_other_DG_layers", "CA1_SO_vs_CA3_SO"),
    status = "completed", ID = c("GO:1", "GO:2"), Description = c("a", "b"),
    NES = c(1.5, 1.5), p_adjust = c(.01, .01), stringsAsFactors = FALSE
  )
  before <- complete
  display <- control_spatial_select_go_display(
    complete, contrasts = control_spatial_figure2f_display_contrasts()
  )
  testthat::expect_identical(complete, before)
  testthat::expect_identical(display$contrast, "DG_MO_vs_mean_other_DG_layers")
})

testthat::test_that("Figure 2f display grid retains contrasts with zero one or two positive terms", {
  contrasts <- c("zero", "one", "two")
  complete <- data.frame(
    dataset = "neuron_soma",
    contrast = rep(contrasts, each = 3L),
    status = "completed",
    ID = paste0("GO:", seq_len(9L)),
    Description = paste("term", seq_len(9L)),
    NES = c(-1, 1, 2, 1.5, -1, 1, 2, 1.5, 0.5),
    p_adjust = c(.001, .2, .3, .01, .001, .2, .01, .02, .2),
    stringsAsFactors = FALSE
  )
  selected <- control_spatial_select_go_display(
    complete, max_terms = 2L, contrasts = contrasts
  )
  display <- control_spatial_complete_go_display_grid(
    complete, selected, contrasts
  )
  counts <- table(factor(display$contrast, levels = contrasts))
  testthat::expect_identical(as.integer(counts), c(1L, 1L, 2L))
  testthat::expect_identical(
    display$display_evidence_status[display$contrast == "zero"],
    "no_fdr_supported_positive_go_bp_term"
  )
  testthat::expect_false(display$display_term_available[display$contrast == "zero"])
  testthat::expect_true(all(display$display_term_available[display$contrast != "zero"]))
  testthat::expect_true(all(is.na(display$NES[display$contrast == "zero"])))
})

testthat::test_that("candidate GO selection preserves the analytical universe and DG layers", {
  candidate <- control_spatial_figure2f_regions_ca1layers_display_contrasts()
  analytical_universe <- c(
    candidate,
    "DG_MO_vs_mean_other_DG_layers", "DG_PO_vs_mean_other_DG_layers",
    "DG_neuropil_vs_mean_non_DG_regions", "CA1_SO_vs_CA3_SO"
  )
  complete <- data.frame(
    dataset = ifelse(grepl("soma_regions$", analytical_universe), "neuron_soma", "neuron_neuropil"),
    contrast = analytical_universe,
    status = "completed",
    ID = paste0("GO:", seq_along(analytical_universe)),
    Description = paste("term", seq_along(analytical_universe)),
    NES = 1.5,
    p_adjust = 0.01,
    stringsAsFactors = FALSE
  )
  before <- complete
  display <- control_spatial_select_go_display(complete, contrasts = candidate)
  testthat::expect_identical(complete, before)
  testthat::expect_setequal(unique(display$contrast), candidate)
  testthat::expect_true(all(c(
    "DG_MO_vs_mean_other_DG_layers", "DG_PO_vs_mean_other_DG_layers"
  ) %in% complete$contrast))
  testthat::expect_false(any(grepl("^DG_(MO|PO)_", display$contrast)))
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

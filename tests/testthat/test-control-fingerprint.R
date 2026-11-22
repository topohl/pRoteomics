# build_control_fingerprint.R and R/spatial/control_fingerprint_utils.R: the
# control spatial fingerprint of Extended Data Fig. 6c/e (kept, Tobias
# 2026-11-22). Pure fixtures; no data or results needed.

source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "control_fingerprint_utils.R"))

CFP_SCRIPT <- "analysis/differential_abundance/build_control_fingerprint.R"
CFP_ORDER <- cfp_spatial_order(repo_path("config", "manuscript_spatial_order.yml"))

cfp_base_fixture <- function(rows) {
  out <- data.frame(
    dataset = vapply(rows, `[[`, "", "dataset"),
    spatial_unit = vapply(rows, `[[`, "", "unit"),
    ProteinGroupID = vapply(rows, `[[`, "", "pg"),
    GeneSymbol = vapply(rows, `[[`, "", "symbol"),
    con_z = vapply(rows, function(r) as.numeric(r$z %||% 0), 0),
    con_mean_log2 = vapply(rows, function(r) as.numeric(r$mean %||% 1), 0),
    stringsAsFactors = FALSE)
  out$gene <- cfp_gene_key(out$GeneSymbol)
  out
}

# one protein group measured in two units of one dataset
cfp_rows <- function(dataset, units, pg, symbol) {
  lapply(units, function(u) list(dataset = dataset, unit = u, pg = pg, symbol = symbol))
}

cfp_dg_soma_contrast <- function() {
  data.frame(
    dataset = "neuron_soma",
    contrast = "DG_vs_mean_other_soma_regions",
    ProteinGroupID = c("PG:neuron_soma:H13", "PG:neuron_soma:PIT", "PG:neuron_soma:XYZ"),
    official_gene_symbol = c("H1f3", "Pitpnm2", "Xyz"),
    logFC = c(2.6, 1.0, -3.0),
    adj.P.Val = c(1e-19, 1e-20, 1e-30),
    stringsAsFactors = FALSE)
}

testthat::test_that("a selected protein whose symbol differs in the baseline stops the run (H1F3)", {
  pick <- cfp_select_contrast_proteins(cfp_dg_soma_contrast())
  testthat::expect_identical(pick$gene, c("PITPNM2", "H1F3"))
  # the retired producer's baseline spelled histone H1.3 the UniProt way
  base <- cfp_base_fixture(c(
    cfp_rows("neuron_soma", c("CA1_sp", "DG_sg"), "PG:neuron_soma:H13", "H1-3"),
    cfp_rows("neuron_soma", c("CA1_sp", "DG_sg"), "PG:neuron_soma:PIT", "Pitpnm2")))
  prot <- cfp_fingerprint_proteins(base, pick)
  testthat::expect_false("H1F3" %in% prot$gene)
  testthat::expect_error(cfp_check_selection_in_proteins(pick, prot),
                         "H1F3 \\(DG_vs_mean_other_soma_regions\\)")
})

testthat::test_that("with canonical symbols the selected protein is shown in every compartment", {
  pick <- cfp_select_contrast_proteins(cfp_dg_soma_contrast())
  base <- cfp_base_fixture(c(
    cfp_rows("neuron_soma", c("CA1_sp", "DG_sg"), "PG:neuron_soma:H13", "H1f3"),
    cfp_rows("neuron_soma", c("CA1_sp", "DG_sg"), "PG:neuron_soma:PIT", "Pitpnm2"),
    cfp_rows("neuron_neuropil", c("CA1_so", "DG_po"), "PG:neuron_neuropil:H13", "H1f3"),
    cfp_rows("microglia", c("CA1", "DG"), "PG:microglia:XYZ", "Xyz")))
  prot <- cfp_fingerprint_proteins(base, pick)
  testthat::expect_silent(cfp_check_selection_in_proteins(pick, prot))
  h1 <- prot[prot$gene == "H1F3", , drop = FALSE]
  testthat::expect_setequal(unique(h1$dataset), c("neuron_soma", "neuron_neuropil"))
  testthat::expect_true(all(h1$selected_ProteinGroupID == "PG:neuron_soma:H13"))
  testthat::expect_true(all(h1$selected_dataset == "neuron_soma"))
  # a symbol on the negative side of the contrast is never shown
  testthat::expect_false("XYZ" %in% prot$gene)
  # the figure readers' columns are all there
  testthat::expect_true(all(c("gene", "dataset", "spatial_unit", "con_z", "contrast")
                            %in% names(cfp_annotate_units(prot, CFP_ORDER))))
})

testthat::test_that("the selected protein group itself must carry the symbol", {
  pick <- cfp_select_contrast_proteins(cfp_dg_soma_contrast())
  # H1F3 is present, but under another protein group than the one selected
  base <- cfp_base_fixture(c(
    cfp_rows("neuron_soma", "DG_sg", "PG:neuron_soma:H13", "H1-3"),
    cfp_rows("neuron_soma", "DG_sg", "PG:neuron_soma:OTHER", "H1f3"),
    cfp_rows("neuron_soma", "DG_sg", "PG:neuron_soma:PIT", "Pitpnm2")))
  prot <- cfp_fingerprint_proteins(base, pick)
  testthat::expect_error(cfp_check_selection_in_proteins(pick, prot),
                         "H1F3 = PG:neuron_soma:H13 in neuron_soma")
})

testthat::test_that("two protein groups on one cell stop the run", {
  contrasts <- data.frame(
    dataset = "neuron_neuropil", contrast = "DG_PO_vs_mean_other_DG_layers",
    ProteinGroupID = "PG:neuron_neuropil:T1", official_gene_symbol = "Tpm3",
    logFC = 1, adj.P.Val = 1e-5, stringsAsFactors = FALSE)
  pick <- cfp_select_contrast_proteins(contrasts)
  base <- cfp_base_fixture(c(
    cfp_rows("neuron_neuropil", "DG_po", "PG:neuron_neuropil:T1", "Tpm3"),
    cfp_rows("neuron_neuropil", "DG_po", "PG:neuron_neuropil:T2", "Tpm3")))
  prot <- cfp_fingerprint_proteins(base, pick)
  testthat::expect_error(cfp_check_selection_in_proteins(pick, prot),
                         "more than one protein group per cell: TPM3 in neuron_neuropil/DG_po")
})

testthat::test_that("selection: enriched side, BH p then |logFC| then symbol, one row per symbol", {
  contrasts <- data.frame(
    dataset = "neuron_neuropil", contrast = "CA1_SLM_vs_mean_other_CA1_strata",
    ProteinGroupID = paste0("PG:neuron_neuropil:", 1:6),
    official_gene_symbol = c("Bbb", "Aaa", "Ccc", "Ccc", NA, "Ddd"),
    logFC = c(1, 1, 3, 2, 5, -4),
    adj.P.Val = c(1e-5, 1e-5, 1e-5, 1e-6, 1e-9, 1e-12),
    stringsAsFactors = FALSE)
  pick <- cfp_select_contrast_proteins(contrasts)
  # Ccc's best group (1e-6) wins, then the three-way tie at 1e-5 breaks on
  # |logFC|, which Ccc's second group would win but its symbol is taken
  testthat::expect_identical(pick$gene, c("CCC", "AAA"))
  testthat::expect_identical(pick$ProteinGroupID, c("PG:neuron_neuropil:4", "PG:neuron_neuropil:2"))
  testthat::expect_identical(nrow(cfp_select_contrast_proteins(contrasts, top_n = 3L)), 3L)
})

testthat::test_that("the CON baseline uses CON animals only and standardises within protein", {
  meta <- data.frame(
    AnimalID = c("A1", "A2", "A3", "A4", "A1", "A2", "A3", "A4"),
    StressGroup = c("CON", "CON", "CON", "SUS", "CON", "CON", "CON", "SUS"),
    SpatialUnit = c(rep("CA1", 4), rep("CA2", 4)),
    stringsAsFactors = FALSE)
  mat <- rbind(P1 = c(1, 2, 3, 100, 3, 4, 5, 100),
               P2 = c(2, 2, 2, 100, 2, 2, 2, -100))
  out <- cfp_con_profile(mat, meta, "neuron_soma", CFP_ORDER)
  testthat::expect_identical(out$spatial_unit, c("CA1_sp", "CA1_sp", "CA2_sp", "CA2_sp"))
  testthat::expect_equal(out$con_mean_log2, c(2, 2, 4, 2))
  testthat::expect_equal(out$con_z[out$ProteinGroupID == "P1"], c(-1, 1) / sqrt(2))
  testthat::expect_true(all(is.na(out$con_z[out$ProteinGroupID == "P2"])))
  testthat::expect_true(all(out$n_con_animals == 3L))
})

testthat::test_that("a signature with fewer than five measured members scores NA", {
  base <- cfp_base_fixture(c(
    lapply(1:5, function(i) list(dataset = "neuron_neuropil", unit = "CA1_so",
                                 pg = paste0("PG", i), symbol = paste0("G", i), z = i)),
    lapply(1:4, function(i) list(dataset = "neuron_neuropil", unit = "CA1_sr",
                                 pg = paste0("PG", i), symbol = paste0("G", i), z = i))))
  sig <- data.frame(signature = "S", signature_family = "f", gene = paste0("G", 1:6),
                    stringsAsFactors = FALSE)
  out <- cfp_score_signatures(base, sig)
  testthat::expect_equal(out$score[out$spatial_unit == "CA1_so"], 3)
  testthat::expect_identical(out$n_measured_members, c(5L, 4L))
  testthat::expect_true(is.na(out$score[out$spatial_unit == "CA1_sr"]))
})

testthat::test_that("spatial units resolve within their dataset and fail closed", {
  testthat::expect_identical(cfp_resolve_unit(c("ca1", "DG"), "neuron_soma", CFP_ORDER),
                             c("CA1_sp", "DG_sg"))
  testthat::expect_identical(cfp_resolve_unit("ca1", "microglia", CFP_ORDER), "CA1")
  testthat::expect_error(cfp_resolve_unit("CA1_slm", "neuron_soma", CFP_ORDER), "Unrecognised")
  a <- cfp_annotate_units(data.frame(dataset = c("neuron_neuropil", "microglia"),
                                     spatial_unit = c("CA2_slm", "CA1"),
                                     stringsAsFactors = FALSE), CFP_ORDER)
  testthat::expect_identical(a$sg_unit, c("CA2_slm", "CA1"))
  testthat::expect_identical(a$sg_layer_display, c("SLM", NA))
  testthat::expect_identical(a$sg_layer_is_resolution, c(TRUE, FALSE))
  testthat::expect_identical(a$sg_compartment_label, c("Neuropil", "Microglia ROI"))
})

testthat::test_that("the baseline must reproduce the atlas CON baseline", {
  base <- data.frame(dataset = "neuron_soma", ProteinGroupID = "P1",
                     spatial_unit = c("CA1_sp", "CA2_sp"), con_mean_log2 = c(2, 4),
                     stringsAsFactors = FALSE)
  atlas <- data.frame(dataset = "neuron_soma", ProteinGroupID = "P1", n_con_animals = 3L,
                      ca1 = 2, ca2 = 4, ca1_slm = NA_real_)
  testthat::expect_identical(cfp_reconcile_atlas(base, atlas, CFP_ORDER)$shared_cells, 2L)
  atlas$ca2 <- 4.001
  testthat::expect_error(cfp_reconcile_atlas(base, atlas, CFP_ORDER), "does not reproduce")
})

testthat::test_that("the selection table records a phenotype-blind source for every row", {
  pick <- cfp_select_contrast_proteins(cfp_dg_soma_contrast())
  sig <- cfp_signature_sets(
    data.frame(external_signature = c("CA1", "SLM"),
               canonical_official_gene_symbol = c("Wfs1", NA), stringsAsFactors = FALSE),
    data.frame(marker_set = c("reference_astrocyte", "reference_cortical_excitatory_neuron"),
               gene_symbol = c("Aqp4", "Satb2"), source_name = "Allen", use_for = "x",
               stringsAsFactors = FALSE))
  testthat::expect_identical(sort(unique(sig$signature)), c("Astrocyte", "Kaulich CA1"))
  sel <- cfp_selection_table(pick, sig)
  testthat::expect_true(all(sel$phenotype_blind == "YES"))
  testthat::expect_false(any(grepl("SUS *- *RES|SUS-RES|sus_res", sel$selection_rule)))
  testthat::expect_true(all(grepl("protein_baseline_spatial_profile", sel$excluded_source_note, fixed = TRUE)))
  testthat::expect_true(all(grepl("SUS - RES", sel$excluded_source_note, fixed = TRUE)))
  testthat::expect_identical(sel$ProteinGroupID[sel$display_form == "A_individual_proteins"],
                             pick$ProteinGroupID)
})

testthat::test_that("the script checks the selection before it writes anything", {
  l <- readLines(repo_path(CFP_SCRIPT), warn = FALSE)
  code <- paste(l[!grepl("^\\s*#", l)], collapse = "\n")
  testthat::expect_match(code, "control_fingerprint_utils.R", fixed = TRUE)
  check_at <- regexpr("cfp_check_selection_in_proteins(", code, fixed = TRUE)
  write_at <- regexpr("write_csv_safe(", code, fixed = TRUE)
  testthat::expect_gt(check_at, 0L)
  testthat::expect_lt(check_at, write_at)
  # the dry run stops before any input is resolved or any directory is created
  dry_at <- regexpr("if (is_dry_run())", code, fixed = TRUE)
  testthat::expect_lt(dry_at, regexpr("resolve_dataset_inputs(", code, fixed = TRUE))
  testthat::expect_lt(dry_at, regexpr("dir_create(", code, fixed = TRUE))
})

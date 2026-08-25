testthat::local_edition(3)

source(testthat::test_path("..", "..", "R", "empirical_roi_marker_utils.R"))

make_empirical_metadata <- function(dataset, animals = sprintf("A%04d", 1:9),
                                    regions = c("CA1", "CA2", "CA3", "DG"),
                                    layers = NULL, one_sided = NULL,
                                    misleading_names = FALSE) {
  if (is.null(layers)) {
    layers <- if (dataset == "neuron_neuropil") c("slm", "so", "sr") else if (dataset == "neuron_soma") "sp" else "microglia"
  }
  meta <- expand.grid(
    AnimalID = animals, region = regions, layer = layers,
    ReplicateGroup = c("Left", "Right"),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  if (!is.null(one_sided)) {
    drop <- meta$AnimalID == one_sided$AnimalID & meta$region == one_sided$region &
      meta$ReplicateGroup == one_sided$drop_side
    meta <- meta[!drop, , drop = FALSE]
  }
  meta$celltype_layer <- dataset
  meta$ExpGroup <- rep(c("CON", "RES", "SUS"), length.out = nrow(meta))
  animal_group <- setNames(rep(c("CON", "RES", "SUS"), length.out = length(animals)), animals)
  meta$ExpGroup <- unname(animal_group[meta$AnimalID])
  name_layer <- if (misleading_names) "WRONG_LAYER" else meta$layer
  meta$sample_id <- paste(meta$AnimalID, meta$region, name_layer, meta$ReplicateGroup, seq_len(nrow(meta)), sep = "_")
  meta <- meta[, c("sample_id", "AnimalID", "ReplicateGroup", "region", "layer", "celltype_layer", "ExpGroup")]
  rownames(meta) <- meta$sample_id
  meta
}

make_empirical_matrix <- function(meta, values = NULL, protein_ids = c("PG1", "PG2")) {
  if (is.null(values)) {
    values <- t(vapply(seq_along(protein_ids), function(i) {
      i + match(meta$AnimalID, unique(meta$AnimalID)) / 10 +
        match(meta$region, c("CA1", "CA2", "CA3", "DG")) / 100 +
        (meta$ReplicateGroup == "Right") / 1000
    }, numeric(nrow(meta))))
  }
  if (is.vector(values)) values <- matrix(values, nrow = length(protein_ids), byrow = TRUE)
  dimnames(values) <- list(protein_ids, meta$sample_id)
  values
}

make_aggregated_model_fixture <- function() {
  animals <- sprintf("A%04d", 1:9)
  datasets <- empirical_roi_dataset_levels()
  proteins <- paste0("PG", 1:6)
  setNames(lapply(seq_along(datasets), function(dataset_index) {
    dataset <- datasets[[dataset_index]]
    observation_id <- paste(animals, dataset, sep = "::")
    mat <- outer(seq_along(proteins), seq_along(animals), function(feature, animal) {
      10 + feature * 0.37 + animal * 0.19 + dataset_index * feature * 0.23 +
        ((feature * animal + dataset_index) %% 5) * 0.07
    })
    dimnames(mat) <- list(proteins, observation_id)
    meta <- data.frame(
      AnimalID = animals,
      ExpGroup = rep(c("CON", "RES", "SUS"), each = 3),
      dataset = dataset,
      observation_id = observation_id,
      stringsAsFactors = FALSE,
      row.names = observation_id
    )
    list(mat = mat, meta = meta)
  }), datasets)
}

marker_stats_fixture <- function() {
  data.frame(
    ProteinGroupID = "PG1", official_gene_symbol = "Aif1",
    gene_mapping_claim_allowed = TRUE,
    detection_rate_microglia = 0.9, detection_rate_neuropil = 0.9,
    detection_rate_soma = 0.9,
    logFC_microglia_vs_neuropil = 1, logFC_microglia_vs_soma = 1,
    logFC_neuropil_vs_microglia = -1, logFC_soma_vs_microglia = -1,
    p_value_microglia_vs_neuropil = 0.001, p_value_microglia_vs_soma = 0.001,
    p_value_neuropil_vs_microglia = 0.001, p_value_soma_vs_microglia = 0.001,
    FDR_microglia_vs_neuropil = 0.01, FDR_microglia_vs_soma = 0.01,
    FDR_neuropil_vs_microglia = 0.01, FDR_soma_vs_microglia = 0.01,
    residual_df = 16, stringsAsFactors = FALSE
  )
}

testthat::test_that("cross-dataset identity uses exact canonical accession groups, not genes", {
  tables <- setNames(lapply(empirical_roi_dataset_levels(), function(dataset) {
    data.frame(
      ProteinGroupID = paste0("PG:", dataset, ":", c("one", "two", "partial")),
      member_accessions = c("P00001", "P00002", "P00003"),
      mapping_status = c("mapped", "mapped", "partially_mapped"),
      official_gene_symbol = c("SameGene", "SameGene", "OtherGene"),
      stringsAsFactors = FALSE
    )
  }), empirical_roi_dataset_levels())
  index <- build_empirical_roi_protein_index(tables)
  eligible <- unique(index$empirical_protein_group_key[index$cross_dataset_model_eligible])

  testthat::expect_setequal(
    eligible,
    c("EPG:exact_accession_set:P00001", "EPG:exact_accession_set:P00002")
  )
  testthat::expect_equal(length(eligible), 2L)
  testthat::expect_false(any(index$cross_dataset_model_eligible[index$exact_member_accessions == "P00003"]))

  aggregated <- setNames(lapply(empirical_roi_dataset_levels(), function(dataset) {
    ids <- tables[[dataset]]$ProteinGroupID
    mat <- matrix(seq_along(ids), ncol = 1, dimnames = list(ids, paste0("A1::", dataset)))
    meta <- data.frame(
      AnimalID = "A1", ExpGroup = "CON", dataset = dataset,
      observation_id = paste0("A1::", dataset), stringsAsFactors = FALSE,
      row.names = paste0("A1::", dataset)
    )
    list(mat = mat, meta = meta)
  }), empirical_roi_dataset_levels())
  aligned <- align_empirical_roi_aggregated(aggregated, index)
  testthat::expect_equal(rownames(aligned$neuron_neuropil$mat), sort(eligible))
  testthat::expect_equal(nrow(aligned$microglia$mat), 2L)

  global_features <- data.frame(
    ProteinGroupID = c("PG:global:one", "PG:global:two"),
    member_accessions = c("P00001", "P00002"), mapping_status = "mapped",
    stringsAsFactors = FALSE
  )
  detection <- expand.grid(
    ProteinGroupID = global_features$ProteinGroupID,
    dataset = empirical_roi_dataset_levels(),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  detection$observed_detection_rate_raw <- seq_len(nrow(detection)) / nrow(detection)
  bridged <- map_empirical_roi_raw_detection(index, global_features, detection)
  testthat::expect_equal(nrow(bridged), 6L)
  testthat::expect_true(all(is.finite(bridged$observed_detection_rate_raw)))
  testthat::expect_true(all(
    bridged$detectability_source == "raw_unified_pre_imputation_exact_accession_bridge"
  ))
})

testthat::test_that("AnimalID, not hemisphere, defines the final replicate grain", {
  aggregated <- lapply(empirical_roi_dataset_levels(), function(dataset) {
    meta <- make_empirical_metadata(dataset)
    aggregate_empirical_roi_dataset(make_empirical_matrix(meta), meta, dataset)
  })
  names(aggregated) <- empirical_roi_dataset_levels()
  model <- build_empirical_roi_model_input(aggregated)

  testthat::expect_equal(ncol(model$expression), 27L)
  testthat::expect_equal(nrow(model$observation_meta), 27L)
  testthat::expect_equal(length(unique(model$observation_meta$AnimalID)), 9L)
  testthat::expect_true(any(grepl("AnimalID", colnames(model$design), fixed = TRUE)))
  testthat::expect_false(any(grepl("ReplicateGroup", colnames(model$design), fixed = TRUE)))
})

testthat::test_that("hemispheres average within spatial unit and a one-sided unit is retained", {
  meta <- make_empirical_metadata(
    "neuron_soma", animals = c("A111", "A0002"),
    one_sided = list(AnimalID = "A111", region = "CA1", drop_side = "Right")
  )
  values <- matrix(0, nrow = 1, ncol = nrow(meta), dimnames = list("PG1", meta$sample_id))
  values[1, meta$AnimalID == "A0002" & meta$region == "CA1" & meta$ReplicateGroup == "Left"] <- 2
  values[1, meta$AnimalID == "A0002" & meta$region == "CA1" & meta$ReplicateGroup == "Right"] <- 6
  values[1, meta$AnimalID == "A111" & meta$region == "CA1"] <- 7
  aggregated <- aggregate_empirical_roi_dataset(values, meta, "neuron_soma")

  testthat::expect_equal(unname(aggregated$hemisphere_mat["PG1", "A0002\037CA1"]), 4)
  testthat::expect_equal(unname(aggregated$hemisphere_mat["PG1", "A111\037CA1"]), 7)
  testthat::expect_equal(aggregated$aggregation_qc$n_bilateral_spatial_units, 7L)
  testthat::expect_equal(aggregated$aggregation_qc$n_one_sided_spatial_units, 1L)
  testthat::expect_equal(aggregated$aggregation_qc$one_sided_spatial_units, "A111::CA1")
})

testthat::test_that("regions receive equal animal-level weight despite unequal side counts", {
  meta <- make_empirical_metadata(
    "neuron_soma", animals = "A111",
    one_sided = list(AnimalID = "A111", region = "CA1", drop_side = "Right")
  )
  values <- matrix(ifelse(meta$region == "CA1", 100, 0), nrow = 1,
                   dimnames = list("PG1", meta$sample_id))
  aggregated <- aggregate_empirical_roi_dataset(values, meta, "neuron_soma")
  testthat::expect_equal(unname(aggregated$mat["PG1", ]), 25)
})

testthat::test_that("neuropil layers are equally averaged before equal region weighting", {
  meta <- make_empirical_metadata("neuron_neuropil", animals = "A0001")
  layer_value <- c(slm = 0, so = 0, sr = 30)
  values <- matrix(ifelse(meta$region == "CA1", layer_value[meta$layer], 0), nrow = 1,
                   dimnames = list("PG1", meta$sample_id))
  aggregated <- aggregate_empirical_roi_dataset(values, meta, "neuron_neuropil")
  testthat::expect_equal(unname(aggregated$region_mat["PG1", "A0001\037CA1"]), 10)
  testthat::expect_equal(unname(aggregated$mat["PG1", ]), 2.5)
})

testthat::test_that("structured anatomy overrides misleading filename labels", {
  meta <- data.frame(
    sample_id = c(
      "raw_A0003_CA2_sp_Neuron_S129.raw", "raw_A0003_CA2_sp_Neuron_S130.raw",
      "raw_A0003_CA2_sr_Neuron_S131.raw", "raw_A0003_CA2_sr_Neuron_S132.raw"
    ),
    AnimalID = "A0003", ReplicateGroup = c("Left", "Right", "Left", "Right"),
    region = "CA2", layer = c("sr", "sr", "sp", "sp"),
    celltype_layer = c("neuron_neuropil", "neuron_neuropil", "neuron_soma", "neuron_soma"),
    ExpGroup = "CON", stringsAsFactors = FALSE
  )
  audit <- validate_known_empirical_metadata_corrections(meta)
  testthat::expect_true(all(audit$correction_valid))
  testthat::expect_equal(audit$layer, c("sr", "sr", "sp", "sp"))

  neuropil <- meta[meta$celltype_layer == "neuron_neuropil", , drop = FALSE]
  rownames(neuropil) <- neuropil$sample_id
  grouped <- aggregate_empirical_matrix_equal_weight(
    matrix(c(2, 6), nrow = 1, dimnames = list("PG1", neuropil$sample_id)),
    neuropil, c("AnimalID", "region", "layer")
  )
  testthat::expect_equal(colnames(grouped$mat), "A0003\037CA2\037sr")
})

testthat::test_that("complete paired 9 by 3 limma design has residual df 16", {
  model_input <- build_empirical_roi_model_input(make_aggregated_model_fixture())
  fit <- build_empirical_roi_limma_inference(model_input)
  testthat::expect_equal(qr(model_input$design)$rank, 11L)
  testthat::expect_equal(unique(fit$stats$nominal_residual_df), 16)
  testthat::expect_equal(unique(fit$stats$residual_df), 16)
  testthat::expect_true(all(fit$stats$n_animals == 9L))
  testthat::expect_true(all(fit$stats$n_observations == 27L))
  testthat::expect_true(all(fit$stats$finite_p_values_all_contrasts))
  testthat::expect_true(all(fit$stats$finite_FDR_all_contrasts))
})

testthat::test_that("NA FDR fails closed and both microglia contrasts are required", {
  testthat::expect_identical(empirical_roi_fdr_pass(c(NA_real_, 0.05, 0.2), 0.1), c(FALSE, TRUE, FALSE))
  stats <- marker_stats_fixture()
  stats$FDR_microglia_vs_soma <- NA_real_
  proposed <- select_empirical_roi_marker_sets(stats)
  testthat::expect_equal(nrow(proposed), 0L)
  testthat::expect_true(all(c(
    "marker_set", "marker_evidence_type", "inferential_test_valid", "claim_allowed"
  ) %in% names(proposed)))

  stats$FDR_microglia_vs_soma <- 0.01
  proposed <- select_empirical_roi_marker_sets(stats)
  micro <- proposed[proposed$marker_set == "empirical_microglia_roi_enriched", , drop = FALSE]
  testthat::expect_equal(nrow(micro), 1L)
  testthat::expect_true(micro$inferential_test_valid)
  testthat::expect_true(micro$claim_allowed)
})

testthat::test_that("shared co-detection is descriptive and never claim-ready", {
  stats <- marker_stats_fixture()
  stats$logFC_microglia_vs_neuropil <- 0.1
  stats$logFC_microglia_vs_soma <- 0.1
  stats$logFC_neuropil_vs_microglia <- -0.1
  stats$logFC_soma_vs_microglia <- -0.1
  stats$FDR_microglia_vs_neuropil <- NA_real_
  shared <- select_empirical_roi_marker_sets(stats)
  shared <- shared[shared$marker_set == "empirical_microglia_neuropil_shared", , drop = FALSE]
  testthat::expect_equal(nrow(shared), 1L)
  testthat::expect_equal(shared$marker_evidence_type, "descriptive_only")
  testthat::expect_false(shared$inferential_test_valid)
  testthat::expect_false(shared$claim_allowed)
  testthat::expect_true(is.na(shared$GeneSymbol))
  testthat::expect_equal(shared$mapped_gene_symbol, "Aif1")
})

testthat::test_that("aggregation is invariant to sample and hemisphere row ordering", {
  meta <- make_empirical_metadata("neuron_neuropil", animals = c("A0001", "A0002"))
  mat <- make_empirical_matrix(meta)
  forward <- aggregate_empirical_roi_dataset(mat, meta, "neuron_neuropil")
  order <- rev(seq_len(nrow(meta)))
  reversed_meta <- meta[order, , drop = FALSE]
  reversed_mat <- mat[, order, drop = FALSE]
  reverse <- aggregate_empirical_roi_dataset(reversed_mat, reversed_meta, "neuron_neuropil")
  testthat::expect_equal(reverse$mat, forward$mat)
  testthat::expect_equal(reverse$region_mat, forward$region_mat)
  testthat::expect_equal(reverse$hemisphere_mat, forward$hemisphere_mat)
})

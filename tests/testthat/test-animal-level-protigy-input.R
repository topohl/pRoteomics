source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "dataset_config.R"))
source(testthat::test_path("..", "..", "R", "dataset_inputs.R"))
source(testthat::test_path("..", "..", "R", "protigy_input_utils.R"))

make_protigy_fixture <- function(dataset = "neuron_neuropil", specs = NULL) {
  if (is.null(specs)) {
    specs <- data.frame(
      AnimalID = c("A0001", "A0001"),
      ExpGroup = c("1", "1"),
      hemisphere = c("Left", "Right"),
      region = c("CA1", "CA1"),
      layer = c("slm", "slm"),
      value = c(1, 3),
      exclude = FALSE,
      stringsAsFactors = FALSE
    )
  }
  if (!"exclude" %in% names(specs)) specs$exclude <- FALSE
  hemi_token <- ifelse(specs$hemisphere == "Left", "L", "R")
  sample_id <- paste(
    "run",
    specs$AnimalID,
    hemi_token,
    specs$region,
    specs$layer,
    if (dataset == "microglia") "Microglia" else "Neuron",
    paste0("sample", seq_len(nrow(specs)), ".d"),
    sep = "_"
  )
  celltype <- if (dataset == "microglia") "microglia" else "neuron"
  celltype_layer <- dataset
  metadata <- data.frame(
    sample_id = sample_id,
    exclude = specs$exclude,
    AnimalID = specs$AnimalID,
    ExpGroup = specs$ExpGroup,
    ReplicateGroup = specs$hemisphere,
    region = specs$region,
    layer = if (dataset == "microglia") "microglia" else specs$layer,
    celltype = celltype,
    celltype_layer = celltype_layer,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  sample_values <- lapply(specs$value, function(value) c(value, value + 10))
  names(sample_values) <- sample_id
  expression <- as.data.frame(sample_values, check.names = FALSE)
  expression$Protein.Group <- c("PG2", "PG1")
  expression$`T: Protein.Names` <- c("protein_b", "protein_a")
  expression$Genes <- c("B", "A")
  expression$First.Protein.Description <- c("B protein", "A protein")
  list(expression = expression, metadata = metadata)
}

prepare_fixture <- function(dataset = "neuron_neuropil", specs = NULL, fail = TRUE) {
  fixture <- make_protigy_fixture(dataset, specs)
  protigy_prepare_animal_level(
    fixture$expression,
    fixture$metadata,
    dataset,
    validate_e9_design = FALSE,
    fail_on_invalid = fail
  )
}

testthat::test_that("bilateral values use an equal-weight mean on existing log2 values", {
  result <- prepare_fixture()
  testthat::expect_equal(result$animal_expression$A0001_CA1_slm, c(2, 12))
  testthat::expect_identical(
    result$aggregation_audit$aggregation_method,
    "equal_weight_mean_LR_on_existing_imputed_log2_values"
  )
  testthat::expect_identical(
    result$aggregation_audit$hemisphere_status,
    "bilateral_complete"
  )
})

testthat::test_that("audit-driven column aggregation preserves gene identity and bilateral means", {
  fixture <- make_protigy_fixture()
  bundle <- protigy_prepare_animal_level(
    fixture$expression,
    fixture$metadata,
    "neuron_neuropil",
    validate_e9_design = FALSE
  )
  sample_ids <- fixture$metadata$sample_id
  gene_matrix <- cbind(c(1, 11), c(3, 13))
  colnames(gene_matrix) <- sample_ids
  rownames(gene_matrix) <- c("Gene_B", "Gene_A")

  aggregated <- protigy_aggregate_expression_columns(gene_matrix, bundle$aggregation_audit)

  testthat::expect_identical(rownames(aggregated), c("Gene_B", "Gene_A"))
  testthat::expect_identical(colnames(aggregated), "A0001_CA1_slm")
  testthat::expect_equal(aggregated[, 1], c(Gene_B = 2, Gene_A = 12))
})

testthat::test_that("audit-driven one-sided aggregation copies the observed column exactly", {
  specs <- data.frame(
    AnimalID = c("A0001", "A0002"),
    ExpGroup = c("1", "2"),
    hemisphere = c("Left", "Right"),
    region = c("CA1", "CA2"),
    layer = c("sp", "sp"),
    value = c(7, 17),
    stringsAsFactors = FALSE
  )
  fixture <- make_protigy_fixture("neuron_soma", specs)
  bundle <- protigy_prepare_animal_level(
    fixture$expression,
    fixture$metadata,
    "neuron_soma",
    validate_e9_design = FALSE
  )
  gene_matrix <- cbind(c(7, 107), c(17, 117))
  colnames(gene_matrix) <- fixture$metadata$sample_id
  rownames(gene_matrix) <- c("Gene_1", "Gene_2")

  aggregated <- protigy_aggregate_expression_columns(gene_matrix, bundle$aggregation_audit)

  testthat::expect_identical(aggregated[, "A0001_CA1_sp"], gene_matrix[, 1])
  testthat::expect_identical(aggregated[, "A0002_CA2_sp"], gene_matrix[, 2])
})

testthat::test_that("audit-driven aggregation fails closed on row or sample identity defects", {
  fixture <- make_protigy_fixture()
  bundle <- protigy_prepare_animal_level(
    fixture$expression,
    fixture$metadata,
    "neuron_neuropil",
    validate_e9_design = FALSE
  )
  gene_matrix <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    dimnames = list(c("Gene_A", "Gene_B"), fixture$metadata$sample_id)
  )
  missing_sample <- gene_matrix[, 1, drop = FALSE]
  testthat::expect_error(
    protigy_aggregate_expression_columns(missing_sample, bundle$aggregation_audit),
    "absent from the gene matrix",
    fixed = TRUE
  )
  rownames(gene_matrix) <- c("Gene_A", "Gene_A")
  testthat::expect_error(
    protigy_aggregate_expression_columns(gene_matrix, bundle$aggregation_audit),
    "row names must be unique",
    fixed = TRUE
  )
})

testthat::test_that("animals are never averaged together", {
  specs <- data.frame(
    AnimalID = rep(c("A0001", "A0002"), each = 2),
    ExpGroup = rep(c("1", "2"), each = 2),
    hemisphere = rep(c("Left", "Right"), 2),
    region = "CA1",
    layer = "slm",
    value = c(1, 3, 10, 14),
    stringsAsFactors = FALSE
  )
  result <- prepare_fixture(specs = specs)
  testthat::expect_equal(result$animal_expression$A0001_CA1_slm[[1]], 2)
  testthat::expect_equal(result$animal_expression$A0002_CA1_slm[[1]], 12)
})

testthat::test_that("neuropil regions and layers remain separate", {
  specs <- data.frame(
    AnimalID = "A0001",
    ExpGroup = "1",
    hemisphere = rep(c("Left", "Right"), 3),
    region = rep(c("CA1", "CA2", "CA1"), each = 2),
    layer = rep(c("slm", "slm", "sr"), each = 2),
    value = c(1, 3, 10, 14, 20, 24),
    stringsAsFactors = FALSE
  )
  result <- prepare_fixture(specs = specs)
  testthat::expect_identical(
    names(result$animal_expression),
    c("id", "Description", "A0001_CA1_slm", "A0001_CA1_sr", "A0001_CA2_slm")
  )
  testthat::expect_equal(result$animal_expression$A0001_CA1_slm[[1]], 2)
  testthat::expect_equal(result$animal_expression$A0001_CA1_sr[[1]], 22)
  testthat::expect_equal(result$animal_expression$A0001_CA2_slm[[1]], 12)
})

testthat::test_that("soma groups by region and retains its unique descriptive layer", {
  result <- prepare_fixture(dataset = "neuron_soma")
  testthat::expect_identical(result$spatial_unit, "region")
  testthat::expect_identical(result$aggregation_audit$canonical_spatial_unit, "CA1")
  testthat::expect_identical(result$output_metadata$layer, "slm")
  testthat::expect_identical(result$output_metadata$sample_id, "A0001_CA1_slm")

  mismatch <- make_protigy_fixture("neuron_soma")
  mismatch$metadata$layer[[2]] <- "sg"
  testthat::expect_error(
    protigy_prepare_animal_level(
      mismatch$expression,
      mismatch$metadata,
      "neuron_soma",
      validate_e9_design = FALSE
    ),
    "descriptive layer"
  )
})

testthat::test_that("microglia groups by region and uses the microglia token", {
  result <- prepare_fixture(dataset = "microglia")
  testthat::expect_identical(result$spatial_unit, "region")
  testthat::expect_identical(result$aggregation_audit$canonical_spatial_unit, "CA1")
  testthat::expect_identical(result$output_metadata$layer, "microglia")
  testthat::expect_identical(result$output_metadata$sample_id, "A0001_CA1_microglia")
})

testthat::test_that("one-sided observed hemispheres are retained without imputation", {
  left_only <- data.frame(
    AnimalID = "A0001", ExpGroup = "1", hemisphere = "Left",
    region = "CA1", layer = "slm", value = 1,
    stringsAsFactors = FALSE
  )
  right_only <- left_only
  right_only$hemisphere <- "Right"
  left_result <- prepare_fixture(specs = left_only)
  right_result <- prepare_fixture(specs = right_only)
  testthat::expect_equal(left_result$animal_expression$A0001_CA1_slm, c(1, 11))
  testthat::expect_equal(right_result$animal_expression$A0001_CA1_slm, c(1, 11))
  testthat::expect_identical(
    left_result$aggregation_audit$hemisphere_status,
    "left_only_observed"
  )
  testthat::expect_identical(
    right_result$aggregation_audit$hemisphere_status,
    "right_only_observed"
  )
  testthat::expect_identical(
    left_result$aggregation_audit$aggregation_method,
    "single_observed_hemisphere_no_imputation"
  )
  testthat::expect_equal(nrow(left_result$output_metadata), 1L)
  testthat::expect_equal(nrow(right_result$output_metadata), 1L)
})

testthat::test_that("neither valid hemisphere is not fabricated", {
  testthat::expect_identical(protigy_unit_status(0L, 0L), "missing_both")
  excluded <- data.frame(
    AnimalID = "A0001", ExpGroup = "1", hemisphere = "Left",
    region = "CA1", layer = "slm", value = 1, exclude = TRUE,
    stringsAsFactors = FALSE
  )
  testthat::expect_error(
    prepare_fixture(specs = excluded),
    "No samples remain after applying exclusions"
  )
})

testthat::test_that("duplicate Left and duplicate Right observations fail", {
  duplicate_left <- data.frame(
    AnimalID = "A0001", ExpGroup = "1",
    hemisphere = c("Left", "Left", "Right"),
    region = "CA1", layer = "slm", value = c(1, 2, 3),
    stringsAsFactors = FALSE
  )
  duplicate_right <- duplicate_left
  duplicate_right$hemisphere <- c("Left", "Right", "Right")
  testthat::expect_error(prepare_fixture(specs = duplicate_left), "Left=2")
  testthat::expect_error(prepare_fixture(specs = duplicate_right), "Right=2")
})

testthat::test_that("hemisphere metadata must agree with the sample name", {
  fixture <- make_protigy_fixture()
  fixture$metadata$ReplicateGroup[[1]] <- "Right"
  testthat::expect_error(
    protigy_prepare_animal_level(
      fixture$expression,
      fixture$metadata,
      "neuron_neuropil",
      validate_e9_design = FALSE
    ),
    "conflicts"
  )
  fixture <- make_protigy_fixture()
  fixture$metadata$ReplicateGroup[[1]] <- "unknown"
  testthat::expect_error(
    protigy_prepare_animal_level(
      fixture$expression,
      fixture$metadata,
      "neuron_neuropil",
      validate_e9_design = FALSE
    ),
    "Only Left/L and Right/R"
  )
})

testthat::test_that("one AnimalID mapping to multiple ExpGroups fails", {
  fixture <- make_protigy_fixture()
  fixture$metadata$ExpGroup <- c("1", "2")
  testthat::expect_error(
    protigy_prepare_animal_level(
      fixture$expression,
      fixture$metadata,
      "neuron_neuropil",
      validate_e9_design = FALSE
    ),
    "ExpGroup"
  )
})

testthat::test_that("excluded samples are audited but never used", {
  specs <- data.frame(
    AnimalID = "A0001", ExpGroup = "1",
    hemisphere = c("Left", "Right", "Left"),
    region = "CA1", layer = "slm", value = c(1, 3, 100),
    exclude = c(FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  result <- prepare_fixture(specs = specs)
  testthat::expect_equal(result$animal_expression$A0001_CA1_slm[[1]], 2)
  excluded <- result$source_sample_assignment$inclusion_status == "excluded_by_metadata"
  testthat::expect_equal(sum(excluded), 1L)
  testthat::expect_equal(result$source_sample_assignment$assigned_output_count[excluded], 0L)
})

testthat::test_that("an excluded opposite hemisphere is not reintroduced", {
  specs <- data.frame(
    AnimalID = "A0001", ExpGroup = "1",
    hemisphere = c("Left", "Right"),
    region = "CA1", layer = "slm", value = c(7, 99),
    exclude = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  result <- prepare_fixture(specs = specs)
  testthat::expect_equal(result$animal_expression$A0001_CA1_slm, c(7, 17))
  testthat::expect_identical(
    result$aggregation_audit$hemisphere_status,
    "left_only_observed"
  )
  testthat::expect_equal(nrow(result$output_metadata), 1L)
  right <- result$source_sample_assignment$hemisphere == "Right"
  right[is.na(right)] <- FALSE
  testthat::expect_false(any(right))
  excluded_assignment <- result$source_sample_assignment$inclusion_status ==
    "excluded_by_metadata"
  testthat::expect_equal(sum(excluded_assignment), 1L)
  testthat::expect_equal(
    result$source_sample_assignment$assigned_output_count[excluded_assignment],
    0L
  )
})

testthat::test_that("strict bilateral sensitivity differs only by one-sided units", {
  specs <- data.frame(
    AnimalID = c("A0001", "A0001", "A0002"),
    ExpGroup = c("1", "1", "2"),
    hemisphere = c("Left", "Right", "Left"),
    region = "CA1", layer = "slm", value = c(1, 3, 10),
    stringsAsFactors = FALSE
  )
  result <- prepare_fixture(specs = specs)
  testthat::expect_identical(
    names(result$animal_expression),
    c("id", "Description", "A0001_CA1_slm", "A0002_CA1_slm")
  )
  testthat::expect_identical(
    names(result$strict_bilateral_expression),
    c("id", "Description", "A0001_CA1_slm")
  )
  testthat::expect_identical(
    result$strict_bilateral_expression$id,
    result$animal_expression$id
  )
  testthat::expect_equal(
    result$strict_bilateral_expression$A0001_CA1_slm,
    result$animal_expression$A0001_CA1_slm
  )
  testthat::expect_identical(
    setdiff(result$output_metadata$sample_id, result$strict_bilateral_metadata$sample_id),
    "A0002_CA1_slm"
  )
})

testthat::test_that("every included source sample maps once to one output", {
  result <- prepare_fixture()
  included <- result$source_sample_assignment$inclusion_status == "included"
  testthat::expect_true(all(result$source_sample_assignment$assigned_output_count[included] == 1L))
  testthat::expect_false(anyDuplicated(result$source_sample_assignment$sample_id[included]) > 0L)
  testthat::expect_equal(
    nrow(result$output_metadata),
    nrow(unique(result$aggregation_audit[c("AnimalID", "canonical_spatial_unit")]))
  )
})

testthat::test_that("output ordering and protein identity are deterministic", {
  specs <- data.frame(
    AnimalID = rep(c("A0002", "A0001"), each = 2),
    ExpGroup = rep(c("2", "1"), each = 2),
    hemisphere = rep(c("Right", "Left"), 2),
    region = "CA1", layer = "slm", value = c(3, 1, 14, 10),
    stringsAsFactors = FALSE
  )
  first <- prepare_fixture(specs = specs)
  second <- prepare_fixture(specs = specs[nrow(specs):1, ])
  testthat::expect_identical(
    names(first$animal_expression),
    c("id", "Description", "A0001_CA1_slm", "A0002_CA1_slm")
  )
  testthat::expect_identical(names(first$animal_expression), names(second$animal_expression))
  testthat::expect_identical(first$feature_ids, c("protein_b", "protein_a"))
  testthat::expect_identical(first$feature_descriptions, c("B protein", "A protein"))
  testthat::expect_identical(first$description_source, "First.Protein.Description")
  testthat::expect_identical(first$description_fallback_count, 0L)
  testthat::expect_identical(second$feature_ids, first$feature_ids)
})

testthat::test_that("strict GCT v1.3 dimensions match animal-level contents", {
  result <- prepare_fixture()
  path <- tempfile(fileext = ".gct")
  write_gct_v1.3(result$annotated_table, path, result$output_metadata)
  dims <- as.integer(strsplit(readLines(path, n = 2)[[2]], "\t", fixed = TRUE)[[1]])
  testthat::expect_identical(
    dims,
    c(result$n_proteins, nrow(result$output_metadata), 1L, ncol(result$output_metadata))
  )
  testthat::expect_true(validate_gct_v1.3(
    path,
    result$n_proteins,
    nrow(result$output_metadata),
    1L,
    ncol(result$output_metadata)
  ))
})

testthat::test_that("ProTigy writer adds Description equals id fallback", {
  metadata <- data.frame(
    sample_id = c("s1", "s2", "s3"),
    group = c("A", "B", "A"),
    AnimalID = c("A1", "A2", "A3"),
    stringsAsFactors = FALSE
  )
  df <- data.frame(
    id = c("group", "AnimalID", "P1", "P2"),
    s1 = c("A", "A1", 1, 4),
    s2 = c("B", "A2", 2, 5),
    s3 = c("A", "A3", 3, 6),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  path <- tempfile(fileext = ".gct")
  write_gct_v1.3(df, path, metadata)
  lines <- readLines(path, warn = FALSE)
  testthat::expect_identical(lines[[2]], "2\t3\t1\t2")
  testthat::expect_identical(lengths(lapply(lines[3:7], gct_v1.3_fields)), rep(5L, 5L))
  testthat::expect_identical(
    lines[3:7],
    c(
      "id\tDescription\ts1\ts2\ts3",
      "group\tna\tA\tB\tA",
      "AnimalID\tna\tA1\tA2\tA3",
      "P1\tP1\t1\t2\t3",
      "P2\tP2\t4\t5\t6"
    )
  )
  testthat::expect_true(validate_gct_v1.3(path, 2L, 3L, 1L, 2L))
})

testthat::test_that("nrhd zero remains structurally valid outside the ProTigy handoff", {
  lines <- c(
    "#1.3",
    "2\t3\t0\t2",
    "id\ts1\ts2\ts3",
    "group\tA\tB\tA",
    "AnimalID\tA1\tA2\tA3",
    "P1\t1\t2\t3",
    "P2\t4\t5\t6"
  )
  path <- tempfile(fileext = ".gct")
  writeLines(lines, path)
  testthat::expect_true(validate_gct_v1.3(path, 2L, 3L, 0L, 2L))
})

testthat::test_that("genuine row descriptors are declared and padded for column metadata", {
  metadata <- data.frame(
    sample_id = c("s1", "s2", "s3"),
    group = c("A", "B", "A"),
    batch = c("x", "x", "y"),
    stringsAsFactors = FALSE
  )
  df <- data.frame(
    id = c("group", "batch", "P1", "P2"),
    Description = c(NA, NA, "desc1", "desc2"),
    s1 = c("A", "x", 1, 4),
    s2 = c("B", "x", 2, 5),
    s3 = c("A", "y", 3, 6),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  path <- tempfile(fileext = ".gct")
  write_gct_v1.3(df, path, metadata)
  lines <- readLines(path, warn = FALSE)
  testthat::expect_identical(lines[[2]], "2\t3\t1\t2")
  testthat::expect_identical(lengths(lapply(lines[3:7], gct_v1.3_fields)), rep(5L, 5L))
  testthat::expect_identical(
    lines[3:7],
    c(
      "id\tDescription\ts1\ts2\ts3",
      "group\tna\tA\tB\tA",
      "batch\tna\tx\tx\ty",
      "P1\tdesc1\t1\t2\t3",
      "P2\tdesc2\t4\t5\t6"
    )
  )
  testthat::expect_true(validate_gct_v1.3(path, 2L, 3L, 1L, 2L))

  fields <- lapply(lines[6:7], gct_v1.3_fields)
  parsed_rdesc <- data.frame(
    Description = vapply(fields, `[[`, character(1), 2L),
    row.names = vapply(fields, `[[`, character(1), 1L),
    stringsAsFactors = FALSE
  )
  parsed_mat <- do.call(rbind, lapply(fields, function(x) as.numeric(x[3:5])))
  rownames(parsed_mat) <- vapply(fields, `[[`, character(1), 1L)
  testthat::expect_identical(rownames(parsed_rdesc), rownames(parsed_mat))
  testthat::expect_identical(parsed_rdesc$Description, c("desc1", "desc2"))
})

testthat::test_that("GCT validator rejects structural and matrix corruption", {
  valid_lines <- c(
    "#1.3",
    "2\t3\t0\t2",
    "id\ts1\ts2\ts3",
    "group\tA\tB\tA",
    "batch\tx\tx\ty",
    "P1\t1\t2\t3",
    "P2\t4\t5\t6"
  )
  expect_invalid <- function(lines, pattern) {
    path <- tempfile(fileext = ".gct")
    writeLines(lines, path)
    testthat::expect_error(validate_gct_v1.3(path), pattern)
  }
  expect_invalid(c("#1.2", valid_lines[-1]), "marker")
  expect_invalid(c(valid_lines[1], "2\t3\t0", valid_lines[-c(1, 2)]), "dimension")
  expect_invalid(valid_lines[-7], "line-count")
  bad_width <- valid_lines
  bad_width[[4]] <- "group\tA\tB"
  expect_invalid(bad_width, "field-count")
  bad_header <- valid_lines
  bad_header[[3]] <- "Name\ts1\ts2\ts3"
  expect_invalid(bad_header, "reserved field 'id'")
  duplicate_samples <- valid_lines
  duplicate_samples[[3]] <- "id\ts1\ts2\ts2"
  expect_invalid(duplicate_samples, "sample IDs")
  duplicate_rows <- valid_lines
  duplicate_rows[[7]] <- "P1\t4\t5\t6"
  expect_invalid(duplicate_rows, "row IDs")
  nonnumeric <- valid_lines
  nonnumeric[[7]] <- "P2\t4\tnot_numeric\t6"
  expect_invalid(nonnumeric, "finite numeric")
  nonfinite <- valid_lines
  nonfinite[[7]] <- "P2\t4\tInf\t6"
  expect_invalid(nonfinite, "finite numeric")

  wrong_descriptor <- c(
    "#1.3", "1\t1\t1\t0", "id\tOther\ts1", "P1\tdesc1\t1"
  )
  expect_invalid(wrong_descriptor, "must name it 'Description'")
  empty_descriptor <- c(
    "#1.3", "1\t1\t1\t0", "id\tDescription\ts1", "P1\t\t1"
  )
  expect_invalid(empty_descriptor, "row-descriptor values")
})

testthat::test_that("phenotypeWithinUnit preserves legacy region_layer_group naming", {
  result <- prepare_fixture()
  testthat::expect_identical(result$output_metadata$phenotypeWithinUnit, "CA1_slm_1")
  microglia <- prepare_fixture("microglia")
  testthat::expect_identical(
    microglia$output_metadata$phenotypeWithinUnit,
    "CA1_microglia_1"
  )
})

testthat::test_that("02_excel_convert shared writer uses cmapR-compatible id semantics", {
  metadata <- data.frame(
    sample_id = c("sample1", "sample2", "sample3"),
    group = c("A", "B", "A"),
    stringsAsFactors = FALSE
  )
  df <- data.frame(
    id = c("group", "protein_a", "protein_b"),
    sample1 = c("A", 1, 4),
    sample2 = c("B", 2, 5),
    sample3 = c("A", 3, 6),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  path <- tempfile(fileext = ".gct")
  write_gct_v1.3(df, path, metadata)
  testthat::expect_identical(
    readLines(path, warn = FALSE),
    c(
      "#1.3",
      "2\t3\t1\t1",
      "id\tDescription\tsample1\tsample2\tsample3",
      "group\tna\tA\tB\tA",
      "protein_a\tprotein_a\t1\t2\t3",
      "protein_b\tprotein_b\t4\t5\t6"
    )
  )

  legacy <- data.frame(
    sample1 = 1,
    Protein.Group = "PG1",
    `T: Protein.Names` = "protein_a",
    Genes = "A",
    check.names = FALSE
  )
  prepared <- protigy_prepare_legacy_expression(legacy)
  testthat::expect_identical(names(prepared), c("id", "Description", "sample1"))
  testthat::expect_identical(prepared$id, "protein_a")
  testthat::expect_identical(prepared$Description, "protein_a")

  alternate_headers <- data.frame(
    sample1 = 1,
    `T: Protein.Group` = "PG1",
    Protein.Names = "protein_a",
    Genes = "A",
    check.names = FALSE
  )
  alternate_prepared <- protigy_prepare_legacy_expression(alternate_headers)
  testthat::expect_identical(names(alternate_prepared), c("id", "Description", "sample1"))
  testthat::expect_identical(alternate_prepared$id, "protein_a")
  testthat::expect_identical(alternate_prepared$Description, "protein_a")

  excel_convert_text <- paste(
    readLines(testthat::test_path("..", "..", "01_preprocessing", "02_excel_convert.r")),
    collapse = "\n"
  )
  testthat::expect_match(
    excel_convert_text,
    'source(repo_path("R", "protigy_input_utils.R"))',
    fixed = TRUE
  )
})

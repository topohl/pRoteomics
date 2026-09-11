source(testthat::test_path("..", "..", "R", "paths.R"))

# Regression tests for the anatomical-axis parsing in the spatial network
# producer.
#
# Every acquisition sample name in this study embeds the instrument method token
# "80SPDzoom". An UNANCHORED search for the layer vocabulary matches the "SP"
# inside that token before the real layer token, because str_extract returns the
# first match anywhere in the string. That collapsed all 180 neuron_neuropil
# samples to Layer = "sp", reducing the region-layer network from 10 anatomical
# units to 4 ("CA1_sp", "CA2_sp", "CA3_sp", "DG_sp") and leaving the similarity
# matrix empty. These tests pin the behaviour so the axis cannot silently
# collapse again.

.spatial_producer <- function() {
  repo_path("07_spatial_networks", "01_network_spatial_relations.r")
}

# The real acquisition naming convention, verbatim in shape.
.fx_sample_names <- function() {
  c(
    "Bluto_20250703_FCo_Evo2_80SPDzoom_Tobias_A111_L_CA1_slm_Neuron_S025_S1-A4_1_12892.d",
    "Bluto_20250703_FCo_Evo2_80SPDzoom_Tobias_A111_R_CA1_so_Neuron_S027_S1-C4_1_12894.d",
    "Bluto_20250703_FCo_Evo2_80SPDzoom_Tobias_A127_L_CA2_sr_Neuron_S029_S1-E4_1_12856.d",
    "Bluto_20250703_FCo_Evo2_80SPDzoom_Tobias_A127_R_CA3_so_Neuron_S031_S1-G4_1_12881.d",
    "Bluto_20250703_FCo_Evo2_80SPDzoom_Tobias_A765_L_DG_mo_Neuron_S037_S1-E5_1_12928.d",
    "Bluto_20250703_FCo_Evo2_80SPDzoom_Tobias_A765_R_DG_po_Neuron_S039_S1-G5_1_12865.d",
    # soma layers, to prove "sp" is still reachable when it is the real token
    "Bluto_20250703_FCo_Evo2_80SPDzoom_Tobias_A129_L_CA1_sp_Neuron_S041_S1-A6_1_12900.d",
    "Bluto_20250703_FCo_Evo2_80SPDzoom_Tobias_A129_R_DG_sg_Neuron_S043_S1-C6_1_12901.d"
  )
}

# The patterns as they appear in the producer, applied the same way.
.extract_with_producer_patterns <- function(sample_names) {
  txt <- paste(readLines(.spatial_producer(), warn = FALSE), collapse = "\n")
  region_pat <- regmatches(
    txt, regexpr('Region = stringr::str_match\\(SampleKey, regex\\("[^"]+"', txt))
  layer_pat <- regmatches(
    txt, regexpr('Layer = stringr::str_match\\(SampleKey, regex\\("[^"]+"', txt))
  testthat::expect_length(region_pat, 1L)
  testthat::expect_length(layer_pat, 1L)
  grab <- function(decl) sub('^.*regex\\("', "", sub('"$', "", decl))
  key <- basename(sample_names)
  list(
    Region = toupper(stringr::str_match(key, stringr::regex(grab(region_pat), ignore_case = TRUE))[, 2]),
    Layer  = tolower(stringr::str_match(key, stringr::regex(grab(layer_pat), ignore_case = TRUE))[, 2])
  )
}

testthat::test_that("the layer token is never taken from the acquisition method string", {
  testthat::skip_if_not_installed("stringr")
  got <- .extract_with_producer_patterns(.fx_sample_names())

  testthat::expect_identical(
    got$Layer,
    c("slm", "so", "sr", "so", "mo", "po", "sp", "sg"))
  # the specific historical failure: not everything collapses to "sp"
  testthat::expect_false(all(got$Layer == "sp"))
  testthat::expect_gt(length(unique(got$Layer)), 1L)
})

testthat::test_that("the region token is parsed correctly alongside it", {
  testthat::skip_if_not_installed("stringr")
  got <- .extract_with_producer_patterns(.fx_sample_names())
  testthat::expect_identical(
    got$Region,
    c("CA1", "CA1", "CA2", "CA3", "DG", "DG", "CA1", "DG"))
})

testthat::test_that("region-layer units recover the full anatomical axis", {
  testthat::skip_if_not_installed("stringr")
  got <- .extract_with_producer_patterns(.fx_sample_names())
  units <- paste(got$Region, got$Layer, sep = "_")
  testthat::expect_true(all(c("CA1_slm", "CA1_so", "CA2_sr", "CA3_so",
                              "DG_mo", "DG_po") %in% units))
  # 4 collapsed pseudo-units are the signature of the bug
  testthat::expect_false(identical(sort(unique(units)),
                                   c("CA1_sp", "CA2_sp", "CA3_sp", "DG_sp")))
})

testthat::test_that("canonical metadata outranks filename parsing", {
  txt <- paste(readLines(.spatial_producer(), warn = FALSE), collapse = "\n")
  # Parsed values may only FILL a missing metadata value, never override one.
  testthat::expect_true(grepl(
    "Region = ifelse\\(is.na\\(Region\\) \\| Region == \"\", parsed\\$Region, Region\\)", txt))
  testthat::expect_true(grepl(
    "Layer = ifelse\\(is.na\\(Layer\\) \\| Layer == \"\", parsed\\$Layer, Layer\\)", txt))
  # the inverted form is what allowed the collapse to survive the join
  testthat::expect_false(grepl("Layer = ifelse\\(!is.na\\(parsed\\$Layer\\)", txt))
  testthat::expect_false(grepl("Region = ifelse\\(!is.na\\(parsed\\$Region\\)", txt))
})

testthat::test_that("the built neuropil network carries the full anatomical axis", {
  p <- path_results("logs", "07_spatial_networks", "network_spatial_relations",
                    "neuron_neuropil", "region_layer",
                    "network_spatial_relations_objects.rds")
  testthat::skip_if_not(file.exists(p), "neuropil spatial network not built")
  o <- readRDS(p)
  testthat::skip_if(is.null(o$nodes), "network object has no nodes table")

  testthat::expect_gt(nrow(o$nodes), 4L)
  testthat::expect_gt(length(unique(o$nodes$Layer)), 1L)
  testthat::expect_false(all(o$nodes$Layer == "sp"))
  # the similarity matrix must actually exist and be square over those nodes
  cm <- if (is.list(o$overall_spearman)) o$overall_spearman$cor_mat else o$overall_spearman
  testthat::expect_true(is.matrix(cm))
  testthat::expect_identical(nrow(cm), ncol(cm))
  testthat::expect_identical(nrow(cm), nrow(o$nodes))
})

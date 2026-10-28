source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "spatial", "neuropil_lamina_utils.R"))

# Figure 2j-l (map_input_territories, map_afferent_sources, compare_lamina_synapses): the shared
# carry-over rule on synthetic data, and the invariants of the written tables when they exist.

LAMINAE <- c("CA1_so", "CA1_sr", "CA1_slm", "CA2_so", "CA2_sr", "CA2_slm", "CA3_so", "CA3_sr", "DG_mo", "DG_po")
SOMA <- c("CA1_sp", "CA2_sp", "CA3_sp", "DG_sg")

testthat::test_that("the DG po carry-over rule excludes only DG po acquisitions above every other lamina", {
  units <- c(SOMA, LAMINAE)
  info <- expand.grid(unit = units, hemisphere = c("Left", "Right"), AnimalID = c("A1", "A2"), stringsAsFactors = FALSE)
  info$Sample <- paste(info$AnimalID, info$hemisphere, info$unit, sep = "_")
  info$group <- "CON"
  # 20 soma-only proteins (seen in every soma acquisition, nowhere else) and 30 everywhere
  obs <- matrix(FALSE, 50, nrow(info), dimnames = list(paste0("p", 1:50), info$Sample))
  obs[21:50, ] <- TRUE
  obs[1:20, info$unit %in% SOMA] <- TRUE
  # carry-over: A1 Left DG po shows 12 soma-only proteins, A2 Right DG po 8; CA1 SR of A2 Left shows 8
  obs[1:12, "A1_Left_DG_po"] <- TRUE
  obs[1:8, "A2_Right_DG_po"] <- TRUE
  obs[1:8, "A2_Left_CA1_sr"] <- TRUE
  carry <- dg_po_carryover(info[, c("Sample", "AnimalID", "group", "hemisphere", "unit")], obs, LAMINAE, SOMA)
  testthat::expect_setequal(carry$unit, LAMINAE)
  testthat::expect_equal(unique(carry$n_soma_only_proteins), 20)
  testthat::expect_equal(unique(carry$carryover_threshold), 8)            # max over the other nine laminae
  testthat::expect_identical(carry$Sample[carry$excluded], "A1_Left_DG_po")   # strictly above; 8 is kept
  testthat::expect_false(any(carry$excluded[carry$unit != "DG_po"]))
})

testthat::test_that("run provenance records the SHA-256 of every output", {
  f <- tempfile(fileext = ".csv")
  writeLines("a,b\n1,2", f)
  prov <- run_output_provenance(list(x = f), "analysis/spatial_validation/map_input_territories.R")
  testthat::expect_identical(prov$output_sha256[[basename(f)]], unname(tools::sha256sum(f)))
  testthat::expect_true(is.logical(prov$code_uncommitted))
})

tables <- function(owner, file) repo_path("results", "spatial_validation", owner, "global", "tables", file)
read_if <- function(p) if (file.exists(p)) utils::read.csv(p, stringsAsFactors = FALSE) else NULL
check_sets <- function(tab, con_animals) {
  # control rows hold only the control animals; replication rows none of them; BH within each set
  testthat::expect_setequal(unique(tab$animal_set), c("CON", "RES+SUS"))
  for (i in seq_len(nrow(tab))) {
    an <- strsplit(tab$animals[i], "; ", fixed = TRUE)[[1]]
    testthat::expect_equal(length(an), tab$n_animals[i])
    if (tab$animal_set[i] == "CON") testthat::expect_true(all(an %in% con_animals)) else
      testthat::expect_false(any(an %in% con_animals))
  }
  for (s in unique(tab$animal_set)) {
    k <- tab$animal_set == s
    testthat::expect_equal(tab$p_adjust[k], stats::p.adjust(tab$p_value[k], method = "BH"), tolerance = 1e-12)
  }
}

testthat::test_that("input-territory tables keep the animal sets apart and carry what Figure 2j reads", {
  sc <- read_if(tables("map_input_territories", "input_territory_scores_by_animal.csv"))
  ts <- read_if(tables("map_input_territories", "input_territory_tests.csv"))
  sv <- read_if(tables("map_input_territories", "input_territory_sensitivity.csv"))
  testthat::skip_if(is.null(sc) || is.null(ts) || is.null(sv), "map_input_territories has not been run here")
  testthat::expect_true(all(c("AnimalID", "group", "territory", "territory_label", "protein_class", "n_proteins",
                              "same_input_z", "same_region_z", "input_organisation_d", "null_median_d", "null_q025_d",
                              "null_q975_d", "excess_over_random", "null_q025_excess", "null_q975_excess",
                              "p_greater_than_random", "p_less_than_random", "p_presynaptic_minus_postsynaptic") %in% names(sc)))
  testthat::expect_true(all(c("animal_set", "territory", "territory_label", "comparison", "estimate", "std_error", "ci_low",
                              "ci_high", "df", "t_value", "p_value", "p_adjust", "n_animals", "animals_positive", "animals",
                              "method") %in% names(ts)))
  check_sets(ts, unique(sc$AnimalID[sc$group == "CON"]))
  testthat::expect_equal(nrow(ts), 12)                       # 2 territories x 3 comparisons x 2 animal sets
  testthat::expect_true(all(c("animal_set", "territory", "pairs", "analysis", "statistic", "estimate", "ci_low", "ci_high",
                              "p_value", "p_lamina_permutation", "n_animals", "animals_positive", "per_animal", "ci_note",
                              "method", "role") %in% names(sv)))
  perm <- sv$analysis == "lamina-label permutation"
  testthat::expect_true(any(perm) && all(is.na(sv$ci_low[perm])) && all(nzchar(sv$ci_note[perm])))
  testthat::expect_true(all(is.finite(sv$ci_low[!perm])))
})

testthat::test_that("afferent-source contrasts keep the animal sets apart", {
  sc <- read_if(tables("map_afferent_sources", "afferent_source_scores_by_animal.csv"))
  ct <- read_if(tables("map_afferent_sources", "afferent_source_contrasts.csv"))
  testthat::skip_if(is.null(sc) || is.null(ct), "map_afferent_sources has not been run here")
  check_sets(ct, unique(sc$AnimalID[sc$group == "CON"]))
})

testthat::test_that("inhibitory-share contrasts keep the animal sets apart and name each pair's animals", {
  an <- read_if(tables("compare_lamina_synapses", "inhibitory_share_by_animal.csv"))
  ct <- read_if(tables("compare_lamina_synapses", "inhibitory_share_contrasts.csv"))
  testthat::skip_if(is.null(an) || is.null(ct), "compare_lamina_synapses has not been run here")
  testthat::expect_true(all(c("gad2_vgat_share_log2", "vglut2_minus_vglut1_log2", "vglut2_acquisitions_observed") %in% names(an)))
  check_sets(ct, unique(an$AnimalID[an$group == "CON"]))
})

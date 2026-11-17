source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "spatial", "neuropil_lamina_utils.R"))

# Figure 2j-l (map_input_territories, map_afferent_sources, compare_lamina_synapses): the shared
# carry-over rule on synthetic data, and the invariants of the written tables when they exist.

LAMINAE <- c("CA1_so", "CA1_sr", "CA1_slm", "CA2_so", "CA2_sr", "CA2_slm", "CA3_so", "CA3_sr", "DG_mo", "DG_po")
SOMA <- c("CA1_sp", "CA2_sp", "CA3_sp", "DG_sg")

# Two animals, both hemispheres of every unit; 20 soma-only proteins (seen in every soma acquisition,
# nowhere else), 30 everywhere and 10 spare proteins seen nowhere. CA1 SR of A2 Left shows 8 soma-only
# proteins, so the threshold is 8.
synthetic_bundle <- function() {
  info <- expand.grid(unit = c(SOMA, LAMINAE), hemisphere = c("Left", "Right"), AnimalID = c("A1", "A2"),
                      stringsAsFactors = FALSE)
  info$Sample <- paste(info$AnimalID, info$hemisphere, info$unit, sep = "_")
  info$group <- "CON"
  obs <- matrix(FALSE, 60, nrow(info), dimnames = list(paste0("p", 1:60), info$Sample))
  obs[21:50, ] <- TRUE
  obs[1:20, info$unit %in% SOMA] <- TRUE
  obs[1:8, "A2_Left_CA1_sr"] <- TRUE
  list(info = info[, c("Sample", "AnimalID", "group", "hemisphere", "unit")], obs = obs)
}

testthat::test_that("the DG po carry-over rule excludes only DG po acquisitions above every other lamina", {
  b <- synthetic_bundle()
  # carry-over: A1 Left DG po shows 12 soma-only proteins, A2 Right DG po 8
  b$obs[1:12, "A1_Left_DG_po"] <- TRUE
  b$obs[1:8, "A2_Right_DG_po"] <- TRUE
  carry <- dg_po_carryover(b$info, b$obs, LAMINAE, SOMA)
  testthat::expect_setequal(carry$unit, LAMINAE)
  testthat::expect_equal(unique(carry$n_soma_only_proteins), 20)
  testthat::expect_equal(unique(carry$carryover_threshold), 8)            # max over the other nine laminae
  testthat::expect_identical(carry$Sample[carry$excluded], "A1_Left_DG_po")   # strictly above; 8 is kept
  testthat::expect_identical(carry$Sample[carry$flagged], "A1_Left_DG_po")
  testthat::expect_false(any(carry$retained_under_floor))                 # A1 keeps its unflagged Right DG po
  testthat::expect_false(any(carry$excluded[carry$unit != "DG_po"]))
})

testthat::test_that("an animal whose DG po acquisitions are all flagged keeps its least-contaminated one", {
  b <- synthetic_bundle()
  b$obs[1:12, "A1_Left_DG_po"] <- TRUE
  b$obs[1:10, "A1_Right_DG_po"] <- TRUE
  carry <- dg_po_carryover(b$info, b$obs, LAMINAE, SOMA)
  po <- carry[carry$unit == "DG_po", , drop = FALSE]
  testthat::expect_setequal(po$Sample[po$flagged], c("A1_Left_DG_po", "A1_Right_DG_po"))
  testthat::expect_identical(po$Sample[po$retained_under_floor], "A1_Right_DG_po")   # 10 < 12
  testthat::expect_false(po$excluded[po$Sample == "A1_Right_DG_po"])
  testthat::expect_identical(carry$Sample[carry$excluded], "A1_Left_DG_po")
  testthat::expect_false(any(carry$retained_under_floor[carry$AnimalID == "A2"]))
  testthat::expect_true(all(carry$excluded == (carry$flagged & !carry$retained_under_floor)))
})

testthat::test_that("the floor breaks ties deterministically", {
  b <- synthetic_bundle()
  b$obs[1:12, c("A1_Left_DG_po", "A1_Right_DG_po")] <- TRUE
  # equal counts and equal shares: the Sample name decides, whatever the row order
  floor_of <- function(bb) {
    carry <- dg_po_carryover(bb$info, bb$obs, LAMINAE, SOMA)
    carry$Sample[carry$retained_under_floor]
  }
  testthat::expect_identical(floor_of(b), "A1_Left_DG_po")
  rev_b <- list(info = b$info[rev(seq_len(nrow(b$info))), ], obs = b$obs[, rev(seq_len(ncol(b$obs)))])
  testthat::expect_identical(floor_of(rev_b), "A1_Left_DG_po")
  # equal counts: the lower soma-only share of the detected proteins decides (Right detects 10 more proteins)
  b$obs[51:60, "A1_Right_DG_po"] <- TRUE
  testthat::expect_identical(floor_of(b), "A1_Right_DG_po")
})

testthat::test_that("an animal with an unflagged DG po keeps only that one", {
  b <- synthetic_bundle()
  b$obs[1:12, "A1_Left_DG_po"] <- TRUE
  b$obs[1:5, "A1_Right_DG_po"] <- TRUE
  carry <- dg_po_carryover(b$info, b$obs, LAMINAE, SOMA)
  a1 <- carry[carry$AnimalID == "A1" & carry$unit == "DG_po", , drop = FALSE]
  testthat::expect_identical(a1$Sample[!a1$excluded], "A1_Right_DG_po")
  testthat::expect_false(any(a1$retained_under_floor))
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
  # control rows hold only the control animals; replication rows none of them; BH within each set and
  # sensitivity, over the rows that carry an adjusted P (post hoc rows have none)
  testthat::expect_setequal(unique(tab$animal_set), c("CON", "RES+SUS"))
  for (i in seq_len(nrow(tab))) {
    an <- strsplit(tab$animals[i], "; ", fixed = TRUE)[[1]]
    testthat::expect_equal(length(an), tab$n_animals[i])
    if (tab$animal_set[i] == "CON") testthat::expect_true(all(an %in% con_animals)) else
      testthat::expect_false(any(an %in% con_animals))
  }
  sens <- if ("sensitivity" %in% names(tab)) tab$sensitivity else rep("none", nrow(tab))
  testthat::expect_true(any(sens == "none" & !is.na(tab$p_adjust)))
  for (s in unique(paste(tab$animal_set, sens))) {
    k <- paste(tab$animal_set, sens) == s & !is.na(tab$p_adjust)
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
  testthat::expect_false(anyNA(ts$p_adjust))
  testthat::expect_true(all(c("animal_set", "territory", "pairs", "analysis", "statistic", "estimate", "ci_low", "ci_high",
                              "p_value", "p_lamina_permutation", "n_animals", "animals_positive", "per_animal", "ci_note",
                              "method", "role") %in% names(sv)))
  perm <- sv$analysis == "lamina-label permutation"
  testthat::expect_true(any(perm) && all(is.na(sv$ci_low[perm])) && all(nzchar(sv$ci_note[perm])))
  testthat::expect_true(all(is.finite(sv$ci_low[!perm])))
  # the strict carry-over rule leaves out the animals kept only by the floor
  st <- sv[sv$analysis == "strict carry-over rule", , drop = FALSE]
  testthat::expect_true(all(st$n_animals < ts$n_animals[match(st$animal_set, ts$animal_set)]))
})

testthat::test_that("afferent-source contrasts keep the animal sets apart", {
  sc <- read_if(tables("map_afferent_sources", "afferent_source_scores_by_animal.csv"))
  ct <- read_if(tables("map_afferent_sources", "afferent_source_contrasts.csv"))
  testthat::skip_if(is.null(sc) || is.null(ct), "map_afferent_sources has not been run here")
  check_sets(ct, unique(sc$AnimalID[sc$group == "CON"]))
  testthat::expect_true(all(ct$sensitivity %in% c("none", "strict carry-over rule")))
  # every pre-declared row carries an adjusted P, as many per animal set and run as the method names; context rows none
  decl <- !grepl("context", ct$role)
  testthat::expect_false(anyNA(ct$p_adjust[decl]))
  testthat::expect_true(all(is.na(ct$p_adjust[!decl])))
  k <- paste(ct$animal_set, ct$sensitivity)[decl]
  testthat::expect_equal(as.vector(table(k)[k]), as.integer(sub(".*BH over the ([0-9]+) contrasts.*", "\\1", ct$method[decl])))
  # the strict-rule rerun has a subset of each set's animals and leaves at least one animal out
  st <- ct[ct$sensitivity == "strict carry-over rule", , drop = FALSE]
  pr <- ct[ct$sensitivity == "none", , drop = FALSE]
  animals_of <- function(d) unique(unlist(strsplit(d$animals, "; ", fixed = TRUE)))
  for (s in unique(st$animal_set))
    testthat::expect_true(all(animals_of(st[st$animal_set == s, ]) %in% animals_of(pr[pr$animal_set == s, ])))
  if (nrow(st)) testthat::expect_lt(length(animals_of(st)), length(animals_of(pr)))
})

testthat::test_that("inhibitory-share contrasts keep the animal sets apart and name each pair's animals", {
  an <- read_if(tables("compare_lamina_synapses", "inhibitory_share_by_animal.csv"))
  ct <- read_if(tables("compare_lamina_synapses", "inhibitory_share_contrasts.csv"))
  testthat::skip_if(is.null(an) || is.null(ct), "compare_lamina_synapses has not been run here")
  testthat::expect_true(all(c("gad2_vgat_share_log2", "vglut2_minus_vglut1_log2", "vglut2_acquisitions_observed") %in% names(an)))
  check_sets(ct, unique(an$AnimalID[an$group == "CON"]))
  # the carry-over sensitivities concern only the pairs with DG po and are not corrected
  s <- ct$sensitivity != "none"
  testthat::expect_true(all(ct$lamina_a[s] == "DG_po" | ct$lamina_b[s] == "DG_po"))
  testthat::expect_true(all(is.na(ct$p_adjust[s])))
  testthat::expect_equal(sum(!s), 2 * 15)                    # 15 pre-declared pairs x 2 animal sets
  testthat::expect_false(anyNA(ct$p_adjust[!s]))
})

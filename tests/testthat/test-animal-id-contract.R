source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "animal_id_contract.R"))

CANON <- aid_expected_exp9_animals()
rd <- function(p) utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
code_of <- function(...) {
  l <- readLines(repo_path(...), warn = FALSE)
  l <- sub("#.*$", "", l)
  paste(l[nzchar(trimws(l))], collapse = "\n")
}

# =====================================================================
# the specific defects
# =====================================================================

testthat::test_that("13856 and 13857 never collapse to the same AnimalID", {
  a <- aid_format_normalize("13856")
  b <- aid_format_normalize("13857")
  testthat::expect_identical(a, "13856")
  testthat::expect_identical(b, "13857")
  testthat::expect_false(identical(a, b))
  # and the historical normaliser really did merge them, so this test is
  # guarding a defect that actually occurred
  testthat::expect_identical(aid_legacy_broken_normalizer("13856"),
                             aid_legacy_broken_normalizer("13857"))
  testthat::expect_identical(aid_legacy_broken_normalizer("13856"), "A1385")
})

testthat::test_that("no terminal digit is ever deleted to force a match", {
  for (v in c("3", "13", "133", "1333", "13333", "133333", "00690", "00696")) {
    r <- aid_format_normalize(v)
    stripped <- sub("^0+(?=[0-9])", "", v, perl = TRUE)
    testthat::expect_identical(r, stripped, info = v)
  }
  # every 5-digit zero-padded id stays distinct; the legacy rule merged them
  five <- c("00690", "00691", "00692", "00693", "00694", "00695", "00696")
  testthat::expect_length(unique(aid_format_normalize(five)), 7L)
  testthat::expect_length(unique(aid_legacy_broken_normalizer(five)), 1L)
})

testthat::test_that("the padding asymmetry that broke the join is gone", {
  # A111 and OR111 are the SAME animal and must compare equal
  testthat::expect_identical(aid_format_normalize("A111"),
                             aid_format_normalize("OR111"))
  testthat::expect_identical(aid_format_normalize("A0003"),
                             aid_format_normalize("0003"))
  testthat::expect_identical(aid_format_normalize("A0003"), "3")
  # the historical rule disagreed on exactly these pairs
  testthat::expect_false(identical(aid_legacy_broken_normalizer("A111"),
                                   aid_legacy_broken_normalizer("OR111")))
  # normalisation is idempotent
  once <- aid_format_normalize(c("A111", "OR111", "0003", "OQ755"))
  testthat::expect_identical(aid_format_normalize(once), once)
})

testthat::test_that("single- and double-digit ids resolve instead of becoming NA", {
  testthat::expect_identical(aid_format_normalize("3"), "3")
  testthat::expect_true(is.na(aid_legacy_broken_normalizer("3")))
  testthat::expect_identical(
    as.character(aid_resolve("3", "movement_auc", strict = TRUE, canonical = CANON)),
    "3")
})

# =====================================================================
# fail-closed behaviour
# =====================================================================

testthat::test_that("unknown ids fail closed in strict mode and are NA in audit mode", {
  testthat::expect_error(
    aid_resolve("ZZZ", "behavior_zscore", strict = TRUE, canonical = CANON),
    "unresolved AnimalID")
  testthat::expect_true(is.na(
    aid_resolve("ZZZ", "behavior_zscore", strict = FALSE, canonical = CANON)[1]))
  # an id with interior letters is not silently squeezed into digits
  testthat::expect_true(is.na(aid_format_normalize("A1B2")))
  testthat::expect_error(
    aid_resolve("A1B2", "behavior_zscore", strict = TRUE, canonical = CANON),
    "unresolved AnimalID")
})

testthat::test_that("OQ754 resolves only when an authoritative alias supports it", {
  # 754 is a real animal in the behaviour cohort but is NOT an Exp9 animal,
  # so against the Exp9 vocabulary it must fail closed rather than be invented
  testthat::expect_true(is.na(
    aid_resolve("OQ754", "behavior_zscore", strict = FALSE, canonical = CANON)[1]))
  testthat::expect_error(
    aid_resolve("OQ754", "behavior_zscore", strict = TRUE, canonical = CANON),
    "unresolved AnimalID")
  # the historical rule invented "A0754" out of it
  testthat::expect_identical(aid_legacy_broken_normalizer("OQ754"), "A0754")
  # with an explicit alias it resolves, and only then
  al <- data.frame(source_system = "behavior_zscore", raw_id = "OQ754",
                   canonical_AnimalID = "754", evidence = "test", status = "active",
                   stringsAsFactors = FALSE)
  testthat::expect_identical(
    aid_resolve("OQ754", "behavior_zscore", strict = TRUE,
                canonical = c(CANON, "754"), aliases = al)[1], "754")
})

testthat::test_that("collisions are a hard error, never a silent merge", {
  testthat::expect_error(
    aid_assert_no_collision(c("00690", "00691"), c("690", "690"), "unit test"),
    "collision")
  # the guard is per source system: mixing two systems in ONE call is itself
  # the error case, because within a single system two raw ids for one animal
  # means the vocabulary is ambiguous
  testthat::expect_error(
    aid_assert_no_collision(c("A111", "OR111"), c("111", "111"), "mixed"),
    "collision")
  testthat::expect_silent(aid_assert_no_collision(c("A111", "A127"),
                                                  c("111", "127"), "one system"))
})

testthat::test_that("two source ids for the SAME animal are allowed across systems", {
  # A111 (proteomics) and OR111 (behaviour) are one animal; that is a legitimate
  # many-to-one and must not trip the guard when each system is checked alone
  testthat::expect_silent(aid_assert_no_collision("A111", "111", "proteomics"))
  testthat::expect_silent(aid_assert_no_collision("OR111", "111", "behaviour"))
  # and both systems agree on the animal
  testthat::expect_identical(
    as.character(aid_resolve("A111", "spatial_network_rds", TRUE, CANON)),
    as.character(aid_resolve("OR111", "behavior_zscore", TRUE, CANON)))
})

# =====================================================================
# the alias table
# =====================================================================

testthat::test_that("the alias table is well formed and one-to-one per system", {
  p <- aid_alias_path()
  testthat::expect_true(file.exists(p))
  a <- aid_alias_table()
  testthat::expect_true(nrow(a) > 0L)
  testthat::expect_setequal(names(a), c("source_system", "raw_id",
    "canonical_AnimalID", "evidence", "status"))
  # one raw id -> at most one AnimalID within a system
  key <- paste(a$source_system, a$raw_id)
  testthat::expect_identical(anyDuplicated(key), 0L)
  # every alias target is a real canonical animal
  testthat::expect_true(all(a$canonical_AnimalID %in% CANON))
  # every alias carries provenance
  testthat::expect_true(all(nzchar(a$evidence)))
  # within a system, two DIFFERENT raw ids may not claim the same animal
  for (s in unique(a$source_system)) {
    z <- a[a$source_system == s, , drop = FALSE]
    testthat::expect_identical(anyDuplicated(z$canonical_AnimalID), 0L, info = s)
  }
})

testthat::test_that("all nine Exp9 animals resolve uniquely in every source system", {
  a <- aid_alias_table()
  for (s in unique(a$source_system)) {
    z <- a[a$source_system == s, , drop = FALSE]
    r <- aid_resolve(z$raw_id, s, strict = TRUE, canonical = CANON)
    testthat::expect_setequal(as.character(r), z$canonical_AnimalID)
    testthat::expect_silent(aid_assert_no_collision(z$raw_id, r, s))
  }
  # and the union over systems covers all nine
  testthat::expect_setequal(unique(a$canonical_AnimalID), CANON)
})

testthat::test_that("the canonical vocabulary still matches the authoritative metadata", {
  p <- path_processed("01_preprocessing", "06_merged_metadata_module_score",
                      "neuron_neuropil",
                      "sample_metadata_merged_clean_for_module_scores.xlsx")
  testthat::skip_if_not(file.exists(p), "canonical metadata not available")
  m <- readxl::read_excel(p)
  testthat::expect_setequal(unique(as.character(m$AnimalID)), CANON)
  # the contract constant must not drift from the metadata
  g <- aid_expected_exp9_groups()
  obs <- unique(as.data.frame(m[, c("AnimalID", "StressGroup")]))
  testthat::expect_identical(
    unname(g[as.character(obs$AnimalID)]), as.character(obs$StressGroup))
  # normalising the metadata cannot create a duplicate animal
  r <- aid_format_normalize(unique(as.character(m$AnimalID)))
  testthat::expect_identical(anyDuplicated(r), 0L)
  testthat::expect_length(r, 9L)
})

testthat::test_that("no dataset gains a duplicate AnimalID through normalisation", {
  for (ds in c("neuron_neuropil", "neuron_soma", "microglia")) {
    p <- path_processed("01_preprocessing", "06_merged_metadata_module_score", ds,
                        "sample_metadata_merged_clean_for_module_scores.xlsx")
    if (!file.exists(p)) next
    m <- readxl::read_excel(p)
    raw <- unique(as.character(m$AnimalID))
    res <- aid_resolve(raw, "proteomics_canonical", strict = TRUE, canonical = CANON)
    testthat::expect_silent(aid_assert_no_collision(raw, res, ds))
    testthat::expect_identical(anyDuplicated(as.character(res)), 0L, info = ds)
  }
})

# =====================================================================
# consumers
# =====================================================================

testthat::test_that("the coupling consumer uses the contract, not a local regex", {
  src <- code_of("analysis/08_integration", "test_network_behaviour_coupling.R")
  testthat::expect_true(grepl("animal_id_contract.R", src, fixed = TRUE))
  testthat::expect_true(grepl("aid_resolve(", src, fixed = TRUE))
  testthat::expect_true(grepl("aid_assert_no_collision(", src, fixed = TRUE))
  # the truncating regex is gone from the code
  testthat::expect_false(grepl("A[0-9]{3,4}|[0-9]{3,4}", src, fixed = TRUE))
  testthat::expect_false(grepl('str_pad(x, 4, pad = "0")', src, fixed = TRUE))
  # the old name is retained only as a hard error, not as a working function
  testthat::expect_true(grepl("normalize_animal_id <- function", src, fixed = TRUE))
  testthat::expect_true(grepl("is withdrawn", src, fixed = TRUE))
})

testthat::test_that("the corrected coupling recovers every animal", {
  d <- path_results("tables", "08_behavior_physio_coupling", "network_behavior_coupling")
  testthat::skip_if_not(dir.exists(d), "coupling not run")
  s <- rd(file.path(d, "join_diagnostics_summary.csv"))
  val <- stats::setNames(as.numeric(s$value), s$metric)
  testthat::expect_identical(unname(val[["animals_lost_from_proteomics"]]), 0)
  testthat::expect_identical(unname(val[["animals_after_join"]]), 9)
  po <- rd(file.path(d, "join_diagnostics_proteomics_only.csv"))
  testthat::expect_identical(nrow(po), 0L)

  # every downstream table now carries the canonical bare AnimalID vocabulary
  for (f in c("physiology_traits_loaded.csv", "movement_auc_z_loaded.csv",
              "merged_global_network_behavior.csv")) {
    z <- rd(file.path(d, f))
    testthat::expect_true(all(as.character(z$AnimalID) %in% CANON), info = f)
    testthat::expect_false(any(grepl("^A", as.character(z$AnimalID))), info = f)
  }
  # physiology is one row per animal, and no animal was dropped by distinct()
  ph <- rd(file.path(d, "physiology_traits_loaded.csv"))
  testthat::expect_identical(nrow(ph), 9L)
  testthat::expect_identical(anyDuplicated(ph$AnimalID), 0L)
  # the group split comes from proteomics, not from the CON/SIS behaviour coding
  g <- rd(file.path(d, "merged_global_network_behavior.csv"))
  testthat::expect_setequal(unique(g$Group), c("CON", "RES", "SUS"))
})

testthat::test_that("the WGCNA-side normaliser touches no canonical module state", {
  src <- readLines(repo_path("analysis/05_wgcna", "score_module_activity.R"),
                   warn = FALSE)
  def <- grep("^normalize_animal_id <- function", src)
  use <- grep("normalize_animal_id\\(", src)
  use <- setdiff(use, def)
  testthat::expect_length(use, 1L)
  # the single call site sits inside the behaviour handoff export, which runs
  # after every module, kME and score computation
  ex <- grep("^export_behavior_proteomics_input <- function", src)
  testthat::expect_true(use > ex)
  testthat::expect_true(def > ex - 20L)
})

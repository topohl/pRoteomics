# Proteins cutreeDynamic leaves unassigned (label 0) must never become a WGCNA module.

if (!exists("repo_path", mode = "function")) {
  source(testthat::test_path("..", "..", "R", "paths.R"))
}

testthat::test_that("cut label 0 takes the unassigned colour and real labels keep their earlier colours", {
  source(repo_path("R", "module_contracts.R"))
  pal <- wgcna_publication_module_palette()
  earlier_mapping <- function(labels) {
    u <- sort(unique(labels))
    as.character(stats::setNames(unname(pal)[seq_along(u)], as.character(u))[as.character(labels)])
  }
  with_zero <- c(0, 1, 1, 2, 0, 3, 2)
  cols <- wgcna_cut_label_colors(with_zero, pal)
  testthat::expect_identical(cols[with_zero == 0], rep("grey", 2))
  testthat::expect_false(unname(pal)[1] %in% cols)
  testthat::expect_identical(cols[with_zero != 0], earlier_mapping(with_zero)[with_zero != 0])
  without_zero <- c(1, 2, 2, 3, 1)
  testthat::expect_identical(wgcna_cut_label_colors(without_zero, pal), earlier_mapping(without_zero))
  testthat::expect_error(wgcna_cut_label_colors(0:length(pal), pal), "unique publication colours")
})

testthat::test_that("unassigned proteins never get a ModuleID", {
  testthat::skip_if_not_installed("dplyr")
  source(repo_path("R", "module_contracts.R"))
  meta <- wgcna_module_color_metadata(c("#7FA6C1", "grey", "#5D7894"))
  testthat::expect_identical(meta$WGCNAInternalColor, c("#5D7894", "#7FA6C1"))
  testthat::expect_identical(meta$ModuleID, c("WGCNA_m01", "WGCNA_m02"))
  testthat::expect_silent(validate_wgcna_module_color_metadata(meta, c("#7FA6C1", "grey", "#5D7894")))
  bad <- rbind(meta, data.frame(ModuleID = "WGCNA_m03", ModuleLegacyID = "WGCNA_grey", ModuleColor = "#8AA0AF",
                                ModuleColorName = "blue_grey", ModuleColorLabel = "blue grey", WGCNAInternalColor = "grey"))
  testthat::expect_error(validate_wgcna_module_color_metadata(bad), "unassigned")

  ids <- paste0("PG:microglia:P", 1:4)
  merged <- stats::setNames(c("#7FA6C1", "grey", "#5D7894", "#7FA6C1"), ids)
  testthat::expect_error(wgcna_integrate_feature_module_metadata(ids, merged, meta), "ModuleID")
  members <- wgcna_integrate_feature_module_metadata(ids[merged != "grey"], merged, meta)
  testthat::expect_equal(nrow(members), 3L)
  testthat::expect_false(anyNA(members$ModuleID))
})

testthat::test_that("mergeCloseModules keeps unassigned proteins out of modules and the build drops their eigengene", {
  testthat::skip_if_not_installed("WGCNA")
  source(repo_path("R", "module_contracts.R"))
  set.seed(7)
  n <- 30
  a <- stats::rnorm(n); b <- stats::rnorm(n)
  x <- cbind(sapply(1:10, function(i) a + stats::rnorm(n, sd = 0.3)),
             sapply(1:10, function(i) b + stats::rnorm(n, sd = 0.3)),
             sapply(1:3, function(i) stats::rnorm(n)))
  colnames(x) <- paste0("PG:microglia:G", seq_len(ncol(x)))
  labels <- c(rep(1, 10), rep(2, 10), rep(0, 3))
  cols <- wgcna_cut_label_colors(labels)
  merged <- suppressMessages(WGCNA::mergeCloseModules(x, cols, cutHeight = 0.25,
                                                      unassdColor = wgcna_unassigned_module_color(), verbose = 0))
  testthat::expect_identical(as.character(merged$colors[labels == 0]), rep("grey", 3))
  testthat::expect_identical(as.character(merged$colors[labels != 0]), cols[labels != 0])
  mes <- WGCNA::orderMEs(merged$newMEs)
  mes <- mes[, colnames(mes) != paste0("ME", wgcna_unassigned_module_color()), drop = FALSE]
  testthat::expect_false("MEgrey" %in% colnames(mes))
  testthat::expect_equal(ncol(mes), 2L)

  script <- paste(readLines(repo_path("analysis/wgcna", "build_wgcna_modules.R"), warn = FALSE), collapse = "\n")
  testthat::expect_match(script, "ModuleColors <- wgcna_cut_label_colors(Modules, colorSeq)", fixed = TRUE)
  testthat::expect_match(script, "unassdColor = wgcna_unassigned_module_color()", fixed = TRUE)
  testthat::expect_match(script, "colnames(mergedMEs) != unassigned_me", fixed = TRUE)
  testthat::expect_false(grepl("palette_vals[seq_len(nmods)]", script, fixed = TRUE))
})

testthat::test_that("the WGCNA universe keeps unassigned proteins out of every module set", {
  testthat::skip_if_not_installed("dplyr")
  source(repo_path("R", "gsea_wgcna_concordance_utils.R"))
  universe <- data.frame(
    ProteinGroupID = paste0("PG:microglia:P", 1:4),
    included_in_wgcna = c(TRUE, TRUE, TRUE, FALSE),
    ModuleID = c("WGCNA_m01", NA, "WGCNA_m02", NA),
    member_accessions = c("Q00001", "Q00002", "Q00003", "Q00004"),
    stringsAsFactors = FALSE
  )
  bundle <- gww_build_universe_bundle(universe, "microglia")
  testthat::expect_setequal(bundle$universe$ProteinGroupID, paste0("PG:microglia:P", 1:3))
  testthat::expect_setequal(names(bundle$module_sets), c("WGCNA_m01", "WGCNA_m02"))
  testthat::expect_false("PG:microglia:P2" %in% unlist(bundle$module_sets))
  universe$ProteinGroupID[1] <- NA
  testthat::expect_error(gww_build_universe_bundle(universe, "microglia"), "missing identity")
})

testthat::test_that("the theme composition reads an empty ModuleID as no module", {
  testthat::skip_if_not_installed("data.table")
  testthat::skip_if_not_installed("readr")
  # the universe as build_wgcna_modules.R writes it: unassigned network proteins have no ModuleID
  universe <- data.frame(
    ProteinGroupID = paste0("PG:microglia:P", 1:6),
    included_in_wgcna = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
    ModuleID = c("WGCNA_m01", NA, "WGCNA_m02", "WGCNA_m01", NA, NA),
    stringsAsFactors = FALSE
  )
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  readr::write_csv(universe, path, na = "")

  # the script's own lines, from the fread of the universe to the module members
  src <- readLines(repo_path("analysis/integration", "test_theme_module_composition.R"), warn = FALSE)
  first <- grep("^  fu <- data.table::fread\\(inp\\[\\[\"wgcna_universe\"\\]\\]", src)
  last <- grep("^  fu_in_module <- ", src)
  testthat::expect_length(first, 1L)
  testthat::expect_length(last, 1L)
  env <- new.env()
  env$inp <- c(wgcna_universe = path)
  eval(parse(text = src[first:last]), envir = env)
  testthat::expect_identical(sum(is.na(env$fu$ModuleID)), 2L)
  testthat::expect_false(any(env$fu$ModuleID == "", na.rm = TRUE))
  testthat::expect_setequal(env$fu$ProteinGroupID, paste0("PG:microglia:P", 1:5))
  testthat::expect_setequal(unique(env$fu_in_module$ModuleID), c("WGCNA_m01", "WGCNA_m02"))

  # module entities and the partition are built from fu_in_module only
  testthat::expect_true(any(grepl("mods_all <- split(fu_in_module$ProteinGroupID, fu_in_module$ModuleID)",
                                  src, fixed = TRUE)))
  testthat::expect_true(any(grepl("merge(st, fu_in_module[, .(ProteinGroupID, module_id = ModuleID)]",
                                  src, fixed = TRUE)))
  testthat::expect_false(any(grepl("split(fu$ProteinGroupID, fu$ModuleID)", src, fixed = TRUE)))
})

testthat::test_that("a cached state with a module below the minimum module size is not reused", {
  source(repo_path("R", "module_contracts.R"))
  ids <- paste0("PG:microglia:P", 1:9)
  meta <- data.frame(ModuleID = c("WGCNA_m01", "WGCNA_m02"), WGCNAInternalColor = c("#486A8A", "#7FA6C1"),
                     stringsAsFactors = FALSE)
  ok <- list(parameters = list(min_module_size = 3),
             mergedColors = stats::setNames(c(rep("#486A8A", 3), rep("#7FA6C1", 4), "grey", "grey"), ids),
             module_color_metadata = meta)
  testthat::expect_silent(validate_wgcna_cached_module_sizes(ok))   # two unassigned proteins are not a module

  small <- ok
  small$mergedColors <- stats::setNames(c(rep("#486A8A", 2), rep("#7FA6C1", 7)), ids)
  testthat::expect_error(validate_wgcna_cached_module_sizes(small), "WGCNA_m01 #486A8A \\(2 proteins\\)")
  testthat::expect_error(validate_wgcna_cached_module_sizes(small), "PROTEOMICS_WGCNA_FORCE_FULL=true", fixed = TRUE)

  no_size <- ok
  no_size$parameters <- list()
  testthat::expect_error(validate_wgcna_cached_module_sizes(no_size), "min_module_size")

  # the build checks the cached state before it hydrates it, and the dry run reports the same check
  src <- readLines(repo_path("analysis/wgcna", "build_wgcna_modules.R"), warn = FALSE)
  guard <- grep("^  validate_wgcna_cached_module_sizes\\(cached_state\\)", src)
  hydrate <- grep("^  hydrated_state <- wgcna_hydrate_cached_state\\(cached_state\\)", src)
  testthat::expect_length(guard, 1L)
  testthat::expect_length(hydrate, 1L)
  testthat::expect_lt(guard, hydrate)
  testthat::expect_true(any(grepl("validate_wgcna_cached_module_sizes(state_early)", src, fixed = TRUE)))
})

# Current artifacts: skipped where the network is not built.
for (ds in c("neuron_neuropil", "neuron_soma", "microglia")) {
  testthat::test_that(paste("the current", ds, "WGCNA state has no module below the minimum module size"), {
    state_path <- canonical_result_path("wgcna", "build_wgcna_modules", ds, "models", "wgcna_final_model_state.rds")
    testthat::skip_if_not(file.exists(state_path), paste("Current WGCNA state not built:", state_path))
    source(repo_path("R", "module_contracts.R"))
    state <- readRDS(state_path)
    sizes <- table(as.character(state$WGCNA_modules_long$ModuleID))
    testthat::expect_gt(length(sizes), 0L)
    testthat::expect_true(all(sizes >= as.numeric(state$parameters$min_module_size)), info = ds)
    testthat::expect_no_error(validate_wgcna_cached_module_sizes(state))
    testthat::expect_false("grey" %in% state$module_color_metadata$WGCNAInternalColor, info = ds)
    testthat::expect_false("MEgrey" %in% colnames(state$mergedMEs), info = ds)
  })
}

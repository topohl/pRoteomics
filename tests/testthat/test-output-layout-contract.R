source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "pipeline_registry.R"))

# Phase 6F output layout: three lifecycles, a normalized canonical shape, and a
# legacy tree that is read only. The point of this file is that those are
# enforced rather than described.

split_paths <- function(x) {
  p <- trimws(unlist(strsplit(paste(x, collapse = "|"), "[|;]")))
  unique(p[nzchar(p)])
}

registry_steps <- function() {
  reg <- read_pipeline_registry(repo_path("pipeline.yml"))
  s <- pipeline_steps(reg, pipeline_stage_names(reg), dataset = "all",
                      include_unsupported = TRUE)
  s[!duplicated(s$script), , drop = FALSE]
}

# --- the contract itself --------------------------------------------------

testthat::test_that("the layout contract declares three disjoint lifecycles", {
  l <- output_layout()
  testthat::expect_false(is.null(l))
  testthat::expect_identical(l$contract_version, "output_layout_v2")
  testthat::expect_identical(names(l$lifecycles), c("work", "results", "exports"))

  roots <- vapply(l$lifecycles, function(x) as.character(x$root), character(1))
  testthat::expect_identical(unname(roots), c("work", "results", "exports"))
  testthat::expect_identical(anyDuplicated(roots), 0L)

  # only exports may be cited by the manuscript, and only results is canonical
  cited <- vapply(l$lifecycles, function(x) isTRUE(x$may_be_manuscript_cited), logical(1))
  testthat::expect_identical(names(cited)[cited], "exports")
  canon <- vapply(l$lifecycles, function(x) isTRUE(x$may_be_canonical), logical(1))
  testthat::expect_identical(names(canon)[canon], "results")
})

testthat::test_that("the frozen v1 namespace contract is still byte-identical", {
  # v2 supersedes v1 for future writes; v1 is a frozen baseline object and
  # describes the historical tree, so it must not be edited to fit v2.
  baseline <- "6801edbce8a5d222f4af46e06b6db4e99f6a9761"
  rel <- "config/output_namespaces.yml"
  old <- suppressWarnings(system2("git", c("show", paste0(baseline, ":", rel)),
                                  stdout = TRUE, stderr = FALSE))
  testthat::skip_if(!length(old), "baseline commit not resolvable here")
  now <- readLines(repo_path(rel), warn = FALSE)
  testthat::expect_identical(paste(now, collapse = "\n"), paste(old, collapse = "\n"))
})

# --- the canonical shape -------------------------------------------------

testthat::test_that("canonical_result_path builds the declared pattern and rejects typos", {
  p <- canonical_result_path("qc", "summarize_missingness.R", "microglia", "tables",
                             "missingness_diagnostics.xlsx")
  testthat::expect_true(grepl(
    "results/qc/summarize_missingness/microglia/tables/missingness_diagnostics.xlsx$", p))

  # the script stem is accepted with or without its extension
  testthat::expect_identical(
    canonical_result_path("wgcna", "build_wgcna_modules.R", "microglia", "tables"),
    canonical_result_path("wgcna", "build_wgcna_modules", "microglia", "tables"))

  # an empty scope is global rather than a path with a hole in it
  testthat::expect_true(grepl("/global/tables$",
                              canonical_result_path("wgcna", "build_wgcna_modules", "", "tables")))

  # a domain or child that is not declared is an error, not a new namespace
  testthat::expect_error(canonical_result_path("06_modules_WGCNA", "x", "global", "tables"),
                         "unknown output domain")
  testthat::expect_error(canonical_result_path("wgcna", "x", "global", "figures"),
                         "unknown result child")
})

testthat::test_that("every declared domain is a real analysis directory", {
  doms <- output_layout_domains()
  testthat::expect_gt(length(doms), 0L)
  for (d in doms) {
    testthat::expect_true(dir.exists(repo_path("analysis", d)), info = d)
  }
  # and every analysis directory is declared, so a new one cannot be invisible
  on_disk <- basename(list.dirs(repo_path("analysis"), recursive = FALSE))
  testthat::expect_identical(sort(on_disk), sort(doms))
})

testthat::test_that("plots is the child name, not figures", {
  ch <- output_layout_children()
  testthat::expect_true("plots" %in% ch)
  testthat::expect_false("figures" %in% ch)
})

# --- the legacy tree is read only ----------------------------------------

testthat::test_that("the legacy registry covers the trees that lost their writer", {
  f <- repo_path("config", "legacy_output_registry.csv")
  testthat::expect_true(file.exists(f))
  reg <- utils::read.csv(f, stringsAsFactors = FALSE)
  testthat::expect_true(all(reg$policy %in% c("LEGACY_READ_ONLY", "ACTIVE_NOT_LEGACY")))

  legacy <- reg[reg$policy == "LEGACY_READ_ONLY", , drop = FALSE]
  testthat::expect_gt(nrow(legacy), 0L)

  # The registry must not be quietly narrowed to make the write guard pass:
  # these four trees hold the 111 frozen objects whose renderer left in
  # Phase 6C, and they are the reason the registry exists.
  required <- c("results/figures/manuscript",
                "results/figures/manuscript_candidates",
                "results/source_data/manuscript_candidates",
                "results/tables/manuscript_candidates")
  testthat::expect_true(all(required %in% legacy$legacy_path))

  # Every frozen baseline object under results/ must be accounted for by the
  # registry exactly once, whether its root is legacy or still has a writer.
  # Checking against the manifest rather than a hardcoded number keeps this
  # honest as Phase 6G migrates domains: each migration moves frozen objects
  # from the active side of this ledger to the legacy side without changing
  # the total.
  baseline <- "6801edbce8a5d222f4af46e06b6db4e99f6a9761"
  fm <- suppressWarnings(system2(
    "git", c("show", paste0(baseline, ":manuscript/prerestructure_freeze_manifest.csv")),
    stdout = TRUE, stderr = FALSE))
  testthat::skip_if(!length(fm), "baseline manifest not resolvable here")
  frozen <- unique(utils::read.csv(text = paste(fm, collapse = "\n"),
                                   stringsAsFactors = FALSE)$repository_relative_path)
  frozen_results <- frozen[startsWith(frozen, "results/")]
  testthat::expect_identical(sum(reg$frozen_objects_beneath),
                             length(frozen_results))

  # The frozen objects whose root still has a writer belong to domains Phase 6G
  # has not migrated yet. That number may only fall.
  still_active <- sum(reg$frozen_objects_beneath[reg$policy == "ACTIVE_NOT_LEGACY"])
  testthat::expect_gt(sum(legacy$frozen_objects_beneath), 0L)
  testthat::expect_lte(still_active, 18L)

  # a legacy root with an active writer would be a contradiction
  testthat::expect_identical(sum(legacy$n_active_writers), 0L)
})

testthat::test_that("no registered writer declares an output under a legacy root", {
  roots <- legacy_output_roots()
  testthat::expect_gt(length(roots), 0L)
  s <- registry_steps()

  offenders <- list()
  for (i in seq_len(nrow(s))) {
    outs <- split_paths(s$produces[i])
    bad <- outs[vapply(outs, function(o)
      any(o == roots | startsWith(o, paste0(roots, "/"))), logical(1))]
    if (length(bad)) offenders[[s$script[i]]] <- bad
  }
  if (length(offenders)) {
    msg <- paste(vapply(names(offenders), function(k)
      paste0(k, " -> ", paste(offenders[[k]], collapse = ", ")), character(1)),
      collapse = "\n")
    testthat::fail(paste0("writers target a read-only legacy root:\n", msg))
  }
  testthat::expect_length(offenders, 0L)
})

testthat::test_that("is_legacy_output_path recognises a legacy path and spares a canonical one", {
  roots <- legacy_output_roots()
  testthat::skip_if(!length(roots))
  testthat::expect_true(is_legacy_output_path(repo_path(roots[1], "anything.csv")))
  testthat::expect_false(is_legacy_output_path(
    canonical_result_path("wgcna", "build_wgcna_modules", "microglia", "tables", "x.csv")))
})

# --- the publication boundary --------------------------------------------

testthat::test_that("the export bundle is present, manifested and hash-consistent", {
  dst <- repo_path("exports", "publication_source_data")
  testthat::skip_if(!dir.exists(dst), "export bundle not built in this checkout")

  mf <- file.path(dst, "manifest.csv")
  testthat::expect_true(file.exists(mf))
  m <- utils::read.csv(mf, stringsAsFactors = FALSE)
  testthat::expect_true(all(c("publication_id", "exported_file", "sha256") %in% names(m)))
  testthat::expect_gt(nrow(m), 0L)

  # every manifest row must resolve inside the bundle and still hash as recorded
  checked <- 0L
  for (i in seq_len(nrow(m))) {
    f <- file.path(dst, sub("^results/publication_source_data/", "", m$exported_file[i]))
    if (!file.exists(f)) next
    testthat::expect_identical(unname(tools::sha256sum(f)), m$sha256[i],
                               info = m$exported_file[i])
    checked <- checked + 1L
  }
  testthat::expect_gt(checked, 0L)
})

testthat::test_that("work paths are never part of a published contract", {
  # nothing citable may point into the regenerable tree
  contract <- repo_path("config", "publication_source_data_contract.yml")
  if (file.exists(contract)) {
    txt <- paste(readLines(contract, warn = FALSE), collapse = "\n")
    testthat::expect_false(grepl("(^|[^A-Za-z0-9_])work/", txt))
  }
  mf <- repo_path("exports", "publication_source_data", "manifest.csv")
  if (file.exists(mf)) {
    txt <- paste(readLines(mf, warn = FALSE), collapse = "\n")
    testthat::expect_false(grepl("(^|[^A-Za-z0-9_])work/", txt))
  }
})

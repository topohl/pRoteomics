source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "export_helpers.R"))

repo <- normalizePath(testthat::test_path("..", ".."), winslash = "/", mustWork = TRUE)

# ---------------------------------------------------------------------------
# 1. Static regression guard: --dry-run must be side-effect free.
#
# The manuscript export scripts previously called dir_create() *before* their
# dry-run guard, so --dry-run created results/manuscript subdirectories. Assert
# structurally that no filesystem-mutating call precedes the guard, so the
# defect cannot be reintroduced without failing here.
# ---------------------------------------------------------------------------
mutating_call_pattern <- paste(
  "dir_create", "dir\\.create", "file\\.copy", "file\\.rename", "file\\.remove",
  "file\\.create", "unlink", "write\\.csv", "write\\.table", "writeLines",
  "saveRDS", "write_run_manifest", "qc_write_csv",
  sep = "|"
)

dry_run_scripts <- c(
  "09_export_pride_journal/08_export_manuscript_figures.R",
  "09_export_pride_journal/09_export_source_data.R"
)

for (rel in dry_run_scripts) {
  testthat::test_that(paste(basename(rel), "does not mutate before its dry-run guard"), {
    src <- readLines(file.path(repo, rel), warn = FALSE)
    guard <- grep("isTRUE\\(dry_run\\)", src)
    testthat::expect_length(guard, 1L)

    code_lines <- !grepl("^\\s*#", src)
    mutating <- which(grepl(mutating_call_pattern, src) & code_lines)
    before_guard <- mutating[mutating < guard[[1]]]

    testthat::expect_identical(
      before_guard, integer(0),
      info = paste0(
        "mutating call(s) before the dry-run guard in ", rel, ": ",
        paste(sprintf("L%d: %s", before_guard, trimws(src[before_guard])), collapse = " | ")
      )
    )
  })
}

testthat::test_that("the export runner only runs dry-run-capable steps under --dry-run", {
  src <- readLines(file.path(repo, "09_export_pride_journal/RUN_EXPORT.R"), warn = FALSE)
  joined <- paste(src, collapse = "\n")
  # Steps without a side-effect-free dry-run path must be skipped, not executed.
  testthat::expect_match(joined, "dry_run_capable", fixed = TRUE)
  testthat::expect_match(joined, "steps <- intersect(steps, dry_run_capable)", fixed = TRUE)
  for (rel in dry_run_scripts) {
    testthat::expect_match(joined, rel, fixed = TRUE)
  }
})

# ---------------------------------------------------------------------------
# 2. Dynamic proof: running each script with --dry-run in an isolated project
#    root creates no directory and no file.
# ---------------------------------------------------------------------------
testthat::test_that("--dry-run creates no directories or files in an isolated root", {
  testthat::skip_if_not(nzchar(Sys.which("Rscript")), "Rscript not on PATH")

  root <- withr::local_tempdir("export_dry_run_")
  dir.create(file.path(root, "09_export_pride_journal"), recursive = TRUE)
  file.copy(file.path(repo, "R"), root, recursive = TRUE)
  for (rel in dry_run_scripts) {
    file.copy(file.path(repo, rel), file.path(root, "09_export_pride_journal"))
  }

  # Minimal readable inputs so the selection step has something to find.
  fixtures <- c(
    "results/tables/04_differential_expression_enrichment/demo/global/table_a.csv",
    "results/source_data/04_differential_expression_enrichment/demo/global/source_a.csv",
    "results/figures/06_modules_WGCNA/group_effects/microglia/panel.svg"
  )
  for (f in fixtures) {
    dir.create(dirname(file.path(root, f)), recursive = TRUE, showWarnings = FALSE)
    writeLines("x", file.path(root, f))
  }

  snapshot <- function() {
    sort(list.files(root, recursive = TRUE, all.files = TRUE,
                    include.dirs = TRUE, no.. = TRUE))
  }
  before <- snapshot()

  withr::local_envvar(c(PROTEOMICS_PROJECT_ROOT = root))
  withr::local_dir(root)

  for (rel in dry_run_scripts) {
    status <- system2("Rscript", c(rel, "--dry-run"), stdout = FALSE, stderr = FALSE)
    testthat::expect_identical(as.integer(status), 0L,
                               info = paste("dry-run exited non-zero for", rel))
  }

  after <- snapshot()
  testthat::expect_identical(after, before)
  testthat::expect_false(dir.exists(file.path(root, "results", "manuscript")))
})

# ---------------------------------------------------------------------------
# 3. Selection excludes removed legacy WGCNA dataset-suffixed aliases.
# ---------------------------------------------------------------------------
testthat::test_that("dataset-suffixed aliases are dropped only when a canonical sibling exists", {
  root <- withr::local_tempdir("alias_filter_")
  d <- file.path(root, "results", "figures", "06_modules_WGCNA", "group_effects", "microglia")
  dir.create(d, recursive = TRUE)

  canonical <- file.path(d, "top_supermodule_effects_dotplot.svg")
  alias_svg <- file.path(d, "top_supermodule_effects_dotplot_soma.svg")
  alias_pdf <- file.path(d, "top_supermodule_effects_dotplot_microglia.pdf")
  canonical_pdf <- file.path(d, "top_supermodule_effects_dotplot.pdf")
  orphan <- file.path(d, "orphan_plot_neuropil.svg")          # no canonical sibling
  descriptive <- file.path(d, "roi_fingerprint_with_soma.svg") # not an alias
  for (f in c(canonical, alias_svg, alias_pdf, canonical_pdf, orphan, descriptive)) {
    writeLines("x", f)
  }

  kept <- drop_legacy_dataset_suffixed_aliases(
    c(canonical, alias_svg, alias_pdf, canonical_pdf, orphan, descriptive)
  )

  testthat::expect_true(canonical %in% kept)
  testthat::expect_true(canonical_pdf %in% kept)
  testthat::expect_false(alias_svg %in% kept)
  testthat::expect_false(alias_pdf %in% kept)
  # Sole surviving output for a dataset is not an alias: keep it rather than
  # silently dropping the only version of that figure.
  testthat::expect_true(orphan %in% kept)
  # A descriptive name that merely ends in a dataset word is not an alias.
  testthat::expect_true(descriptive %in% kept)
})

testthat::test_that("alias filter matches extension and is empty-input safe", {
  root <- withr::local_tempdir("alias_filter_ext_")
  dir.create(root, showWarnings = FALSE, recursive = TRUE)
  canonical_svg <- file.path(root, "fig.svg")
  alias_pdf <- file.path(root, "fig_soma.pdf")   # canonical .pdf absent
  writeLines("x", canonical_svg)
  writeLines("x", alias_pdf)

  # Extension must match: fig.svg does not make fig_soma.pdf an alias.
  testthat::expect_true(alias_pdf %in% drop_legacy_dataset_suffixed_aliases(alias_pdf))
  testthat::expect_identical(drop_legacy_dataset_suffixed_aliases(character(0)), character(0))
})

testthat::test_that("the stale manuscript payload cannot be re-collected as an export source", {
  # Without this, rebuilding would re-ingest legacy aliases already sitting in
  # results/manuscript/extended_data, recreating the defect it is meant to fix.
  stale <- file.path(
    tempdir(), "results", "manuscript", "extended_data",
    "06_modules_WGCNA_interpretable_summary_neuron_soma_top_supermodule_effects_dotplot_soma.svg"
  )
  live <- file.path(
    tempdir(), "results", "figures", "06_modules_WGCNA", "interpretable_summary",
    "neuron_soma", "top_supermodule_effects_dotplot.svg"
  )
  testthat::expect_false(is_exportable_result_path(stale))
  testthat::expect_true(is_exportable_result_path(live))
  # The renamed-aside superseded payload must stay excluded too.
  superseded <- file.path(
    tempdir(), "results", "manuscript", "_superseded_20260622",
    "extended_data", "top_supermodule_effects_dotplot_soma.svg"
  )
  testthat::expect_false(is_exportable_result_path(superseded))
})

# ---------------------------------------------------------------------------
# 4. Export routing points at current canonical outputs.
# ---------------------------------------------------------------------------
testthat::test_that("figure export routes from canonical result roots only", {
  src <- readLines(
    file.path(repo, "09_export_pride_journal/08_export_manuscript_figures.R"), warn = FALSE
  )
  joined <- paste(src, collapse = "\n")
  # Canonical EWCE root, not the comparison branch.
  testthat::expect_match(joined, "canonical_ewce_figure_root()", fixed = TRUE)
  testthat::expect_false(grepl("EWCE_E9_comparison", joined, fixed = TRUE))
  # Both selection filters applied.
  testthat::expect_match(joined, "is_exportable_result_path(candidates)", fixed = TRUE)
  testthat::expect_match(joined, "drop_legacy_dataset_suffixed_aliases(candidates)", fixed = TRUE)
})

testthat::test_that("current Figure-2 source data carries the accepted CA1 slot-2 term", {
  f <- file.path(
    repo, "results", "source_data", "04_differential_expression_enrichment",
    "control_spatial_identity_validation", "global", "figure2f_source_data.csv"
  )
  testthat::skip_if_not(file.exists(f), "Figure-2f source data not present in this checkout")
  d <- utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  ca1 <- d[d$contrast == "CA1_vs_mean_other_soma_regions", , drop = FALSE]
  testthat::expect_identical(
    trimws(ca1$Description),
    c("RNA splicing", "mRNA metabolic process")
  )
  # The superseded slot-2 selection must not reappear for CA1.
  testthat::expect_false("mRNA processing" %in% trimws(ca1$Description))
})

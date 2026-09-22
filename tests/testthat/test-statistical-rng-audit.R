# Statistical RNG reproducibility, as distinct from rendering RNG.
#
# Phases 6H.6 and 6H.7 closed RENDERING RNG: a seed there can move only
# geometry. This file governs STATISTICAL RNG, where a seed can move a sampled
# observation, a permutation p-value, an FDR status, and therefore a reported
# conclusion. The two must not be conflated, and the constants must not be
# shared - which is why NATURE_JITTER_SEED and NATURE_REPEL_SEED are asserted
# below to be absent from every statistical call site.
#
# The audit itself is Phase 6H.8; these tests pin its conclusions so that a
# newly added unseeded statistical draw fails rather than passing unnoticed.

source(testthat::test_path("..", "..", "R", "paths.R"))

AUDIT_CSV <- repo_path("audits", "phase6h_statistical_rng_audit.csv")

# Phase 6H.8 found exactly one active statistical RNG call without a seed - the
# compareGO bootstrap - and allowlisted it here because it was unreachable.
# Phase 6H.10 archived that region, so the allowlist is now EMPTY: every
# statistical draw in active code carries a deterministic seed source, and any
# new unseeded one fails this file with nowhere to hide.
KNOWN_UNSEEDED <- character(0)

RENDER_FNS <- c("position_jitter", "position_jitterdodge", "geom_jitter",
                "geom_text_repel", "geom_label_repel")
DRAW_FNS <- c("sample", "sample.int", "slice_sample", "sample_n", "sample_frac",
              "runif", "rnorm", "rbinom", "rexp", "rpois", "rgamma", "rbeta",
              "rmultinom", "rhyper", "rnbinom", "rlnorm", "rweibull")
LIB_FNS <- c("kmeans", "Rtsne", "umap", "randomForest", "cv.glmnet", "boot",
             "fgsea", "fgseaSimple", "fgseaMultilevel", "GSEA", "gseGO",
             "gseKEGG", "nmf", "Mclust")

# Discovery is parse-tree based: a comment or a string mentioning sample() is
# not a call, and a line grep counts it.
rng_calls <- function(roots = c("analysis", "R", "tools"),
                      fns = c(DRAW_FNS, LIB_FNS)) {
  paths <- vapply(roots, repo_path, character(1))
  files <- unlist(lapply(paths[dir.exists(paths)], list.files,
                         pattern = "[.][Rr]$", recursive = TRUE, full.names = TRUE))
  out <- list()
  for (f in files) {
    pd <- tryCatch(utils::getParseData(parse(f, keep.source = TRUE)),
                   error = function(e) NULL)
    if (is.null(pd) || !nrow(pd)) next
    hits <- pd[pd$token == "SYMBOL_FUNCTION_CALL" & pd$text %in% fns, , drop = FALSE]
    if (!nrow(hits)) next
    rel <- sub(paste0("^", repo_path(), "/"), "", gsub("\\\\", "/", f))
    for (k in seq_len(nrow(hits)))
      out[[length(out) + 1L]] <- list(file = rel, line = hits$line1[k], fn = hits$text[k])
  }
  out
}

testthat::test_that("the statistical RNG audit table exists with its contract columns", {
  testthat::skip_if_not(file.exists(AUDIT_CSV), "audit table not present")
  a <- utils::read.csv(AUDIT_CSV, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("file", "function", "line_call", "purpose", "active",
                "publication_facing", "rendering_or_statistical", "explicit_seed",
                "seed_value_or_source", "consumes_global_rng",
                "downstream_artifact", "inferential_role", "action_needed")
  testthat::expect_true(all(required %in% names(a)),
    info = paste("missing:", paste(setdiff(required, names(a)), collapse = ", ")))
  testthat::expect_gt(nrow(a), 0L)
})

testthat::test_that("every active statistical RNG call is explained", {
  # The Phase 6H.8 hard gate: a statistical draw with no recorded purpose, no
  # seed provenance or no inferential role is unadjudicated, and an
  # unadjudicated draw is exactly what makes a frozen number unreproducible.
  testthat::skip_if_not(file.exists(AUDIT_CSV), "audit table not present")
  a <- utils::read.csv(AUDIT_CSV, stringsAsFactors = FALSE, check.names = FALSE)
  s <- a[a$rendering_or_statistical == "statistical", , drop = FALSE]
  testthat::expect_gt(nrow(s), 0L)
  blank <- function(x) !nzchar(trimws(as.character(x)))
  testthat::expect_identical(sum(blank(s$purpose)), 0L)
  testthat::expect_identical(sum(blank(s$seed_value_or_source)), 0L)
  testthat::expect_identical(sum(blank(s$inferential_role)), 0L)
  testthat::expect_identical(sum(blank(s$downstream_artifact)), 0L)
  testthat::expect_true(all(s$inferential_role %in%
    c("PRIMARY_INFERENCE", "SENSITIVITY_ANALYSIS", "DESCRIPTIVE_DIAGNOSTIC",
      "COMPUTATIONAL_SHORTCUT", "VISUALIZATION_ONLY", "UNUSED")))
})

testthat::test_that("only the adjudicated call lacks a statistical seed", {
  testthat::skip_if_not(file.exists(AUDIT_CSV), "audit table not present")
  a <- utils::read.csv(AUDIT_CSV, stringsAsFactors = FALSE, check.names = FALSE)
  s <- a[a$rendering_or_statistical == "statistical", , drop = FALSE]
  unseeded <- unique(s$file[!as.logical(s$explicit_seed)])
  testthat::expect_identical(setdiff(unseeded, KNOWN_UNSEEDED), character(0),
    info = paste("new unseeded statistical RNG:",
                 paste(setdiff(unseeded, KNOWN_UNSEEDED), collapse = ", ")))
  # and every call whose role is PRIMARY_INFERENCE is seeded, without exception
  prim <- s[s$inferential_role == "PRIMARY_INFERENCE", , drop = FALSE]
  testthat::expect_gt(nrow(prim), 0L)
  testthat::expect_true(all(as.logical(prim$explicit_seed)),
    info = "a primary-inference draw is unseeded")
})

testthat::test_that("rendering seeds are never used for statistical sampling", {
  # Brief section 13: presentation geometry and scientific sampling must not
  # share a constant, or changing a figure's layout would move a number.
  testthat::skip_if_not(file.exists(AUDIT_CSV), "audit table not present")
  a <- utils::read.csv(AUDIT_CSV, stringsAsFactors = FALSE, check.names = FALSE)
  s <- a[a$rendering_or_statistical == "statistical", , drop = FALSE]
  testthat::expect_false(any(grepl("NATURE_JITTER_SEED|NATURE_REPEL_SEED",
                                   s$seed_value_or_source)))
  # and no statistical call site in the tree passes a rendering constant
  for (c in rng_calls()) {
    src <- readLines(file.path(repo_path(), c$file), warn = FALSE)
    line <- src[c$line]
    testthat::expect_false(grepl("NATURE_JITTER_SEED|NATURE_REPEL_SEED", line),
      info = paste(c$file, c$line, "uses a rendering seed for a statistical draw"))
  }
})

testthat::test_that("statistical and rendering seeds stay independently declared", {
  # gsea_seed_base happens to be the same integer as NATURE_REPEL_SEED and
  # NATURE_JITTER_SEED. That is a coincidence of value, not a dependency, and it
  # must stay that way: deriving one from another - or "tidying" all three into
  # one shared constant - would couple figure layout to GSEA p-values.
  pn <- readLines(repo_path("R", "utilities", "plotting_nature.R"), warn = FALSE)
  render_defs <- grep("^NATURE_(REPEL|JITTER)_SEED <-", pn, value = TRUE)
  testthat::expect_identical(length(render_defs), 2L)
  # each render seed is a bare literal, not a reference to a statistical seed
  testthat::expect_false(any(grepl("gsea|seed_base|config|yaml|GO_", render_defs,
                                   ignore.case = TRUE)))

  cfg <- repo_path("config", "compareGO_config.yml")
  rp <- repo_path("analysis", "differential_abundance", "run_clusterprofiler_enrichment.R")
  testthat::skip_if_not(file.exists(rp), "clusterProfiler script not present")
  rps <- readLines(rp, warn = FALSE)
  # and the statistical seed base is never sourced from the rendering library
  testthat::expect_false(any(grepl("NATURE_REPEL_SEED|NATURE_JITTER_SEED", rps)))
  if (file.exists(cfg))
    testthat::expect_false(any(grepl("NATURE_", readLines(cfg, warn = FALSE))))
})

testthat::test_that("the discovered call set matches the audited call set", {
  # Guards the audit against drift: a statistical draw added to active code
  # that never reaches the table would otherwise be invisible.
  testthat::skip_if_not(file.exists(AUDIT_CSV), "audit table not present")
  a <- utils::read.csv(AUDIT_CSV, stringsAsFactors = FALSE, check.names = FALSE)
  audited <- paste0(a$file, ":", sub(":.*$", "", a$line_call))
  found <- vapply(rng_calls(), function(c) paste0(c$file, ":", c$line), character(1))
  testthat::expect_identical(sort(setdiff(found, audited)), character(0),
    info = paste("statistical RNG calls missing from the audit:",
                 paste(setdiff(found, audited), collapse = ", ")))
})

testthat::test_that("the compareGO bootstrap is gone from the active surface", {
  # Phase 6H.8 recommended R1 because the slice_sample bootstrap sat below an
  # unconditional quit() and could not execute; this test was the tripwire that
  # asserted it was still PRESENT but unreachable. Phase 6H.10 implemented
  # decision C2 and archived the whole unreachable region, so the premise is
  # spent: the call is not unreachable, it is absent.
  f <- repo_path("analysis", "differential_abundance", "compare_go_enrichment.R")
  testthat::skip_if_not(file.exists(f), "compare_go_enrichment.R not present")
  ex <- parse(f, keep.source = TRUE)
  pd <- utils::getParseData(ex)

  # Phase 6H.10 implemented C2 and archived the unreachable region, so the
  # bootstrap is no longer merely unreachable - it is absent. The replacement
  # invariant is stronger and still live: no statistical draw anywhere in the
  # canonical path, and the exit as the final expression rather than a divider
  # with dead code behind it.
  testthat::expect_identical(
    sum(pd$token == "SYMBOL_FUNCTION_CALL" & pd$text == "slice_sample"), 0L)
  testthat::expect_identical(
    sum(pd$token == "SYMBOL_FUNCTION_CALL" &
          pd$text %in% c("sample", "sample.int", "runif", "rnorm", "set.seed")), 0L,
    info = "a statistical draw has appeared in the canonical compareGO path")

  last <- ex[[length(ex)]]
  testthat::expect_true(is.call(last) && identical(as.character(last[[1]])[1], "quit"),
    info = "the canonical exit is no longer the final expression - has a new tail appeared?")
  testthat::expect_true(all(vapply(as.list(last)[-1],
    function(a) !is.call(a) && !is.name(a), logical(1))))

  # the obsolete marker must not come back with a new tail
  testthat::expect_false(any(grepl("LEGACY_COMPAREGO_TAIL_DISABLED_BY_CANONICAL_EXIT",
                                   readLines(f, warn = FALSE), fixed = TRUE)))
})

testthat::test_that("WGCNA::pickSoftThreshold is deterministic, not stochastic", {
  # It was flagged by a library-name heuristic and cleared empirically. Recording
  # that here stops it being re-flagged, and catches an upstream change that
  # made it stochastic.
  testthat::skip_if_not_installed("WGCNA")
  set.seed(7)
  X <- matrix(stats::rnorm(30 * 40), nrow = 30,
              dimnames = list(paste0("s", 1:30), paste0("g", 1:40)))
  run <- function(s) {
    set.seed(s)
    utils::capture.output(
      r <- suppressWarnings(suppressMessages(WGCNA::pickSoftThreshold(
        X, networkType = "signed", powerVector = c(4, 6), verbose = 0))))
    r$fitIndices
  }
  testthat::expect_equal(run(101), run(999))
})

testthat::test_that("the GSEA RNG scope is deterministic and stream-preserving", {
  # This is the repository's existing statistical-seed convention and the model
  # for any future remediation: L'Ecuyer-CMRG plus a per-comparison seed inside
  # withr::with_preserve_seed.
  source(repo_path("R", "enrichment", "clusterprofiler_reproducibility.R"))
  testthat::expect_true(is.function(run_with_stable_gsea_rng))

  draw <- function() runif(3)
  set.seed(11); a <- run_with_stable_gsea_rng(draw, gsea_seed = 4242L)
  set.seed(99); b <- run_with_stable_gsea_rng(draw, gsea_seed = 4242L)
  # same seed, different incoming state -> identical result
  testthat::expect_identical(a, b)
  # different seed -> different result, so the seed is live
  set.seed(11); c1 <- run_with_stable_gsea_rng(draw, gsea_seed = 777L)
  testthat::expect_false(identical(a, c1))

  # and the caller's stream is left exactly as it was
  set.seed(4242)
  before <- .Random.seed
  invisible(run_with_stable_gsea_rng(draw, gsea_seed = 4242L))
  testthat::expect_true(identical(.Random.seed, before))
  # including the RNG kind, which the wrapper switches and restores
  kind_before <- RNGkind()
  invisible(run_with_stable_gsea_rng(draw, gsea_seed = 4242L))
  testthat::expect_identical(RNGkind(), kind_before)
})

testthat::test_that("fgsea really does depend on the RNG, so the scope is needed", {
  # If fgsea were deterministic the wrapper would be pointless; it is not.
  testthat::skip_if_not_installed("fgsea")
  set.seed(3)
  stats_v <- stats::setNames(stats::rnorm(200), paste0("g", 1:200))
  pw <- list(P1 = paste0("g", 1:20), P2 = paste0("g", 40:80))
  run <- function(s) {
    set.seed(s)
    # fgsea emits a progress bar that neither suppressMessages nor
    # capture.output intercepts, so this test adds two bars to the suite log.
    # Cosmetic; the suite already carries similar output from other packages.
    r <- suppressWarnings(suppressMessages(fgsea::fgseaSimple(
      pathways = pw, stats = stats_v, minSize = 3, maxSize = 500,
      nperm = 200, nproc = 1)))
    stats::setNames(r$pval, r$pathway)
  }
  testthat::expect_false(isTRUE(all.equal(run(1), run(999))))
})

testthat::test_that("the compareGO bootstrap statistic is closed-form in expectation", {
  # Brief section 14. The block at compare_go_enrichment.R:2939-2975 resamples
  # enrichment-result ROWS that already carry p.adjust; it recomputes no
  # enrichment, no p-value and no FDR. So "recovery" is just the probability
  # that a pre-existing significant row appears at least once in a size-n
  # resample, whose expectation is fixed by the per-term row multiplicities:
  #     P(recovered) = 1 - prod_c (1 - m_tc/n_c)^n_c
  # The seed moves the Monte Carlo error of that estimate, never the estimand.
  testthat::skip_if_not_installed("dplyr")
  set.seed(4242)
  terms <- paste0("GO_", sprintf("%03d", 1:120))
  df <- do.call(rbind, lapply(paste0("C", 1:4), function(cc) {
    k <- sample(seq_along(terms), 80L)
    data.frame(Comparison = cc, Description = terms[k],
               p.adjust = ifelse(stats::runif(80L) < 0.08,
                                 stats::runif(80L, 1e-4, 0.049),
                                 stats::runif(80L, 0.06, 0.9)),
               stringsAsFactors = FALSE)
  }))
  top <- utils::head(unique(df$Description[df$p.adjust < 0.05]), 15L)

  statistic <- function(n_boot = 40L) {
    res <- integer(n_boot)
    for (b in seq_len(n_boot)) {
      s <- df |>
        dplyr::group_by(.data$Comparison) |>
        dplyr::slice_sample(prop = 1, replace = TRUE) |>
        dplyr::ungroup()
      res[b] <- sum(top %in% unique(s$Description[s$p.adjust < 0.05]))
    }
    mean(res) / length(top)
  }

  n_c <- table(df$Comparison)
  sig <- df[df$p.adjust < 0.05, , drop = FALSE]
  expected <- mean(vapply(top, function(t) {
    m <- table(sig$Comparison[sig$Description == t])
    1 - prod(vapply(names(m), function(cc)
      (1 - as.numeric(m[[cc]]) / as.numeric(n_c[[cc]]))^as.numeric(n_c[[cc]]),
      numeric(1)))
  }, numeric(1)))

  set.seed(1);   r1 <- statistic()
  set.seed(999); r2 <- statistic()
  # the seed does move the reported number ...
  testthat::expect_false(identical(r1, r2))
  # ... but only within Monte Carlo error of a fixed closed-form quantity
  testthat::expect_lt(abs(mean(c(r1, r2)) - expected), 0.08)
  # and the estimand itself does not depend on the RNG at all
  testthat::expect_true(is.finite(expected))
})

testthat::test_that("the compareGO bootstrap recomputes no inferential quantity", {
  # The decisive structural fact: a resample can only re-draw p.adjust values
  # that already exist, so no new p-value or FDR status can be created.
  testthat::skip_if_not_installed("dplyr")
  set.seed(7)
  df <- data.frame(
    Comparison = rep(c("A", "B"), each = 50L),
    Description = paste0("GO_", sprintf("%03d", c(1:50, 1:50))),
    p.adjust = stats::runif(100L), stringsAsFactors = FALSE)
  for (s in c(1L, 42L, 999L)) {
    set.seed(s)
    r <- df |>
      dplyr::group_by(.data$Comparison) |>
      dplyr::slice_sample(prop = 1, replace = TRUE) |>
      dplyr::ungroup()
    testthat::expect_identical(sum(!(r$p.adjust %in% df$p.adjust)), 0L)
    testthat::expect_identical(nrow(r), nrow(df))
  }
})

testthat::test_that("the stochastic call's downstream artifacts are untouched", {
  # Brief section 16: this audit must change no artifact. The only outputs the
  # compareGO bootstrap ever produced are superseded copies.
  base <- path_results("manuscript", "_superseded_20260622", "supplementary_tables")
  testthat::skip_if_not(dir.exists(base), "superseded bundle not present")
  expected <- c(
    "08_Bootstrap_Stability_Summary.xlsx" = 15431,
    "04_differential_expression_enrichment_compareGO_microglia_BP_phenotype_within_unit_08_Bootstrap_Stability_Summary.xlsx" = 5135,
    "04_differential_expression_enrichment_compareGO_neuron_soma_BP_phenotype_within_unit_08_Bootstrap_Stability_Summary.xlsx" = 5136,
    "04_differential_expression_enrichment_compareGO_neuron_neuropil_BP_phenotype_within_unit_08_Bootstrap_Stability_Summary.xlsx" = 5143)
  for (nm in names(expected)) {
    f <- file.path(base, nm)
    testthat::skip_if_not(file.exists(f), paste(nm, "absent"))
    testthat::expect_identical(file.size(f), expected[[nm]], info = nm)
  }
  # and no live copy has appeared in the manuscript package
  live <- path_results("manuscript", "supplementary_tables")
  if (dir.exists(live))
    testthat::expect_identical(
      length(list.files(live, pattern = "Bootstrap_Stability_Summary")), 0L)
})

testthat::test_that("the curated workbook is out of PRIDE and retained internally", {
  # Phase 6H.8 found this workbook sitting in the PRIDE payload, flagged
  # intended_for_PRIDE = TRUE by a path rule rather than by curation, and
  # pinned it here so it could not quietly change. Phase 6H.9 adjudicated P2
  # and Phase 6H.10 implemented it, so the assertion is inverted: the outward
  # copy must be ABSENT and the internal provenance copy must carry the same
  # bytes. The detailed gates live in test-comparego-tail-archival.R.
  f <- repo_path("pride_submission", "supplementary_tables",
                 paste0("results_tables_04_differential_expression_enrichment_",
                        "compareGO_neuron_neuropil_BP_phenotype_within_unit_",
                        "08_Bootstrap_Stability_Summary.xlsx"))
  # Deliberately NOT skip_if_not(file.exists(f)). An earlier version gated the
  # whole block on the file existing, so its disappearance would have produced
  # a silent skip instead of a signal. A tripwire that vanishes with the thing
  # it watches is not a tripwire.
  testthat::expect_false(file.exists(f),
    info = "the workbook is back in pride_submission - P2 has been reverted")

  curated <- path_results("manuscript", "_curated",
                 paste0("results_tables_04_differential_expression_enrichment_",
                        "compareGO_neuron_neuropil_BP_phenotype_within_unit_",
                        "08_Bootstrap_Stability_Summary.xlsx"))
  testthat::expect_true(file.exists(curated),
    info = "the internal provenance copy is missing - P2 retains the bytes")
  if (file.exists(curated)) {
    testthat::expect_identical(file.size(curated), 5141)
    testthat::expect_identical(
      unname(tools::sha256sum(curated)),
      "024d3671f2cd6026e8bfd7feb6d9839c7ff23eeb41c98d8f66235d854185baba")
  }

  man <- repo_path("pride_submission", "manifests", "pride_file_manifest.tsv")
  testthat::skip_if_not(file.exists(man), "PRIDE manifest not present")
  m <- utils::read.delim(man, stringsAsFactors = FALSE, check.names = FALSE)
  row <- m[grepl("Bootstrap_Stability_Summary", m$file_path, fixed = TRUE), , drop = FALSE]

  # The invariant that must hold under EITHER decision, checked without any
  # dependence on the file being present: on-disk presence and manifest
  # endorsement agree. A manifest row for an absent file, or an unlisted file on
  # disk, is a bookkeeping break either way.
  testthat::expect_identical(nrow(row) > 0L, file.exists(f),
    info = "pride_submission presence and manifest endorsement have diverged")

  if (nrow(row)) {
    testthat::expect_identical(nrow(row), 1L)
    testthat::expect_identical(row$export_category[1], "pride_staging")
    testthat::expect_true(as.logical(row$intended_for_PRIDE[1]))
  }
})

testthat::test_that("the accepted package and freeze are unchanged by this audit", {
  mp <- path_results("manuscript", "figure_export_manifest.csv")
  testthat::skip_if_not(file.exists(mp), "figure export manifest not present")
  testthat::expect_identical(
    unname(tools::sha256sum(mp)),
    "0fd0c9ed9febbc05ed6928b7fc6bfffdab845daba6255d8b47c607dee5d3c02c")
  fz <- repo_path("docs", "publication_freeze_manifest.yml")
  testthat::skip_if_not(file.exists(fz), "freeze manifest not present")
  testthat::expect_identical(
    unname(tools::sha256sum(fz)),
    "b4d37250360e2e07136ccdfbca64bd0a629e22946ced1f5ffd1fb730b445c49a")
})

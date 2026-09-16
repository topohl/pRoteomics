# Guards for the imported Figure 1 behavioural evidence bridge.
#
# The behavioural analysis lives in a separate repository and is not reproduced
# here, so the risk this phase carries is not a wrong model but a wrong quote:
# a number drifting away from the frozen contract, a superseded Phase-2B
# assertion surviving, or prohibited wording creeping back into prose. Until
# now nothing under manuscript/ was read by any test at all. These guards close
# that gap.

repo <- function(...) file.path(testthat::test_path("..", ".."), ...)
MS <- repo("manuscript")
BRIDGE <- file.path(MS, "figure1_bridge_mmmsociability")
DRAFT <- file.path(MS, "manuscript_draft.md")
rd <- function(p) utils::read.csv(p, stringsAsFactors = FALSE)

BUNDLE_FILES <- c("figure1_claim_contract.csv", "figure1_methods_contract.csv",
                  "figure1_repo_provenance.csv",
                  "figure1_source_data_manifest.csv",
                  "figure1_timeline_contract.csv")

test_that("the bridge is a byte-exact mirror of exactly the frozen bundle", {
  skip_if_not(dir.exists(BRIDGE), "figure 1 bridge not imported")
  # Exactly these files and nothing else - an extra entry here would be a
  # proteomics-side artefact masquerading as frozen upstream evidence. The
  # source_data/ subdirectory is the one permitted addition: it is the frozen
  # Figure 1 panel source data, itself a byte-exact upstream mirror, and it is
  # pinned to its own expected contents below rather than waved through.
  expect_setequal(list.files(BRIDGE), c(BUNDLE_FILES, "source_data"))
  expect_true(dir.exists(file.path(BRIDGE, "source_data")))
  expect_setequal(
    list.files(file.path(BRIDGE, "source_data")),
    c("figure1a_timeline_source.csv", "figure1b_combz_classification_source.csv",
      "figure1c_movement_combz_source.csv", "figure1d_loao_predictions_source.csv",
      "figure1e_permutation_source.csv", "figure1f_repeated_cv_source.csv",
      "figure1_panel_statistics.csv", "figure1_panel_source_manifest.csv"))
  # and nothing executable may live in the evidence interface
  expect_equal(list.files(BRIDGE, pattern = "[.][Rr]$", recursive = TRUE),
               character(0))

  man <- rd(file.path(MS, "figure1_bridge_import_manifest.csv"))
  expect_equal(nrow(man), length(BUNDLE_FILES))
  skip_if_not(requireNamespace("digest", quietly = TRUE), "digest unavailable")
  for (i in seq_len(nrow(man))) {
    p <- repo(man$imported_file[i])
    expect_true(file.exists(p))
    expect_equal(file.size(p), man$bytes[i])
    expect_equal(digest::digest(p, algo = "sha256", file = TRUE),
                 man$sha256[i])
  }
})

test_that("the bridge records the full upstream provenance chain", {
  skip_if_not(dir.exists(BRIDGE), "figure 1 bridge not imported")
  p <- rd(file.path(MS, "figure1_bridge_provenance.csv"))
  val <- function(k) p$value[match(k, p$element)]
  expect_equal(val("source_repository"), "topohl/MMMSociability")
  # The commits the import is pinned to. These are full 40-character shas and
  # move whenever the upstream bundle is re-frozen, so assert their shape and
  # their mutual consistency rather than literal values that would have to be
  # edited on every resync.
  for (k in c("analysis_commit", "bundle_commit", "verified_source_HEAD"))
    expect_match(val(k), "^[0-9a-f]{40}$")
  # the re-freeze is itself the verified head, so these two agree
  expect_equal(val("bundle_commit"), val("verified_source_HEAD"))
  # the facts that make the quote trustworthy at all
  expect_equal(val("behavioural_code_copied"), "NONE")
  expect_equal(val("behavioural_results_recomputed"), "NONE")
  expect_equal(val("scientific_values_changed_by_this_resync"), "NONE")
  # the bundle must have been exported from a clean source tree; the first
  # freeze recorded "dirty at export time" and that is what FC-10 was about
  expect_equal(val("bundle_internal_worktree_state"), "clean")
  expect_equal(val("source_worktree_state_at_import"), "clean")
})

test_that("the imported contract carries the corrected prediction seeds", {
  skip_if_not(dir.exists(BRIDGE), "figure 1 bridge not imported")
  fm <- rd(file.path(BRIDGE, "figure1_methods_contract.csv"))
  design <- fm$value[fm$methods_id == "FM-07"]
  # 521 assigns the CV folds. 123 is the association-bootstrap seed and must
  # never appear in the prediction-design row: that substitution is exactly the
  # defect this bridge was re-frozen to repair.
  expect_true(grepl("seed 521", design, fixed = TRUE))
  expect_false(grepl("123", design, fixed = TRUE))
  assoc <- fm$value[fm$methods_id == "FM-06"]
  expect_true(grepl("seed 123", assoc, fixed = TRUE))
  expect_false(grepl("521", assoc, fixed = TRUE))
})

test_that("no behavioural analysis code was copied into this repository", {
  # The bridge is a CSV evidence interface. An .R file appearing under it, or a
  # copy of the upstream analysis entrypoint anywhere here, would mean the
  # separation had been broken.
  skip_if_not(dir.exists(BRIDGE), "figure 1 bridge not imported")
  expect_equal(list.files(BRIDGE, pattern = "[.][Rr]$"), character(0))
  upstream <- c("09_early_prediction_model_ladder",
                "build_figure1_export_bundle", "build_figure1_manuscript_bridge",
                "build_later_outcome_combz")
  here <- list.files(repo("."), pattern = "[.][Rr]$", recursive = TRUE)
  here <- here[!grepl("^(\\.git|renv)/", here)]
  for (u in upstream) expect_false(any(grepl(u, basename(here), fixed = TRUE)))
})

test_that("line-ending normalisation is disabled for the bridge", {
  # .gitattributes sets * text=auto repo-wide. Under that rule git would store
  # these CRLF files as LF and check them out differently on another platform,
  # so the recorded SHA-256 would not reproduce and the integrity guard above
  # would fail for a reason that has nothing to do with the evidence.
  ga <- readLines(repo(".gitattributes"), warn = FALSE)
  expect_true(any(grepl("figure1_bridge_mmmsociability/\\*\\*\\s+-text", ga)))
})

# ------------------------------------------------ manuscript prose obligations

test_that("Results 1 is written and is no longer blocked", {
  skip_if_not(file.exists(DRAFT), "draft absent")
  d <- readLines(DRAFT, warn = FALSE)
  h <- grep("^## 1\\. ", d)
  expect_length(h, 1L)
  # the heading must not reintroduce the term the outcome is not named by
  expect_false(grepl("stress burden", d[h], ignore.case = TRUE))
  nxt <- grep("^## 2\\. ", d)
  expect_length(nxt, 1L)
  sec <- d[(h + 1L):(nxt - 1L)]
  expect_false(any(grepl("\\[BLOCKED", sec, fixed = FALSE)))
  expect_false(any(grepl("METHOD DETAIL UNRESOLVED", sec, fixed = TRUE)))
  # a real section, not a stub
  expect_gt(length(unlist(strsplit(paste(sec, collapse = " "), "\\s+"))), 600L)
})

test_that("Results 4 stays closed and claims no integration result", {
  skip_if_not(file.exists(DRAFT), "draft absent")
  d <- readLines(DRAFT, warn = FALSE)
  h <- grep("^## 4\\. ", d)
  expect_length(h, 1L)
  sec <- paste(d[h:(h + 12L)], collapse = " ")
  expect_true(grepl("CLOSED", sec, fixed = TRUE))
  expect_true(grepl("BH-006", sec, fixed = TRUE))
})

test_that("manuscript prose carries no wording the frozen contract prohibits", {
  skip_if_not(file.exists(DRAFT), "draft absent")
  d <- readLines(DRAFT, warn = FALSE)
  # Standing repository convention: a line that instructs against a phrase, or
  # quotes it in order to replace it, necessarily contains it. Exempt those,
  # exactly as the v9 semantics layer and the PH-008 guard do.
  DENIAL <- paste0("never|must not|do not |does not|cannot|prohibited|",
                   "instead of|rather than|avoid|banned|forbidden|not licensed|",
                   "\\bNOT\\b|\\bNEVER\\b|is not |are not |no longer")
  BANNED <- c(
    "female-specific", "sex-specific", "stronger in females",
    "driven by females", "predicts susceptibility", "predicts resilience",
    "externally validated", "external validation", "independent replication",
    "strong prediction", "highly accurate", "biomarker",
    # names a different, group-adjusted model upstream
    "movement-only", "Movement only")
  # A line-level denial regex is too blunt for prose: "a sensitivity analysis and
  # not independent replication" is plainly a denial, but matches none of
  # "is not", "are not" or "does not". Check instead whether the banned phrase is
  # itself negated, by looking at the words immediately before each occurrence.
  # That catches the real failure - a bald assertion - without exempting a whole
  # line merely because the word "not" appears somewhere else on it.
  NEGATED <- function(line, term) {
    pos <- gregexpr(term, line, ignore.case = TRUE)[[1]]
    if (pos[1] == -1L) return(FALSE)
    all(vapply(pos, function(p) {
      lead <- substr(line, max(1L, p - 40L), p - 1L)
      grepl("(not|never|neither|nor|without|rather than|instead of)[^[:alnum:]]*$",
            lead, ignore.case = TRUE, perl = TRUE)
    }, logical(1)))
  }
  for (b in BANNED) {
    hits <- grep(b, d, ignore.case = TRUE, value = TRUE)
    hits <- hits[!grepl(DENIAL, hits, perl = TRUE)]
    hits <- hits[!vapply(hits, NEGATED, logical(1), term = b)]
    expect_equal(hits, character(0),
                 info = paste("prohibited wording in draft:", b))
  }

  # The provenance tables are manuscript-facing too - a claim row is what a
  # reviewer reads when tracing a sentence - and scanning only the draft let
  # "movement-only" survive in C1-8 and M-26 through a whole phase. Scan the
  # claim and statement tables on the wording that names a real upstream model.
  # Scan column-aware, not line-aware: a prohibited_wording cell exists in order
  # to name banned phrasing, so matching it there is correct rather than a
  # defect. Everything that asserts something is fair game.
  NAME_BANNED <- c("movement-only", "movement only")
  ASSERTIVE <- function(nm) !grepl("prohibit|banned|forbidden|phase2b|conflict",
                                   nm, ignore.case = TRUE)
  for (f in c("results_claim_provenance.csv", "methods_statement_provenance.csv",
              "results_statement_provenance.csv")) {
    p <- file.path(MS, f)
    if (!file.exists(p)) next
    tb <- rd(p)
    for (nm in names(tb)[vapply(tb, is.character, logical(1))]) {
      if (!ASSERTIVE(nm)) next
      for (b in NAME_BANNED) {
        hits <- grep(b, tb[[nm]], ignore.case = TRUE, value = TRUE)
        hits <- hits[!grepl(DENIAL, hits, perl = TRUE)]
        expect_equal(hits, character(0),
                     info = paste0("prohibited model name in ", f, " column ",
                                   nm, ": ", b))
      }
    }
  }
})

test_that("no provenance row states the bootstrap seed as the CV seed", {
  skip_if_not(dir.exists(BRIDGE), "figure 1 bridge not imported")
  # 521 assigns folds; 123 drives the association bootstrap. A row describing
  # the repeated grouped CV must never quote 123.
  for (f in c("methods_statement_provenance.csv",
              "behavior_prediction_contract.csv",
              "figure1_red_team_review.csv")) {
    p <- file.path(MS, f)
    if (!file.exists(p)) next
    ln <- readLines(p, warn = FALSE)
    bad <- grep("(five-fold|5-fold|grouped)[^\"]{0,120}seed 123", ln,
                ignore.case = TRUE, value = TRUE)
    # a row that explicitly records the historical defect is legitimate
    bad <- bad[!grepl("previously stated|was the|transcription defect|FC-02|default",
                      bad, ignore.case = TRUE)]
    expect_equal(bad, character(0),
                 info = paste("CV seed misstated as 123 in", f))
  }
})

test_that("CombZ is described in the licensed direction and never renamed", {
  skip_if_not(file.exists(DRAFT), "draft absent")
  d <- paste(readLines(DRAFT, warn = FALSE), collapse = " ")
  # the sign contract: this is the one error that would invert the whole result
  expect_true(grepl("[Hh]igher CombZ indicates a more resilient-like", d))
  expect_false(grepl("stress burden score", d, ignore.case = TRUE))
  expect_false(grepl("[Hh]igher CombZ[^.]{0,40}worse", d))
})

# --------------------------------------------------- numbers match the bundle

test_that("every behavioural number in Results 1 matches the frozen contract", {
  skip_if_not(dir.exists(BRIDGE), "figure 1 bridge not imported")
  skip_if_not(file.exists(DRAFT), "draft absent")
  cc <- rd(file.path(BRIDGE, "figure1_claim_contract.csv"))
  fm <- rd(file.path(BRIDGE, "figure1_methods_contract.csv"))
  # Prose uses a typographic minus (U+2212); the contracts use ASCII hyphen.
  # Normalise both sides or every negative number compares unequal for a purely
  # typographic reason.
  ascii <- function(z) gsub("−", "-", z, fixed = TRUE)
  d <- ascii(readLines(DRAFT, warn = FALSE))
  h <- grep("^## 1\\. ", d); nxt <- grep("^## 2\\. ", d)
  sec <- paste(d[(h + 1L):(nxt - 1L)], collapse = " ")

  # Each value must appear in the drafted section AND be present in the frozen
  # contract text, so the draft cannot drift from the bundle in either
  # direction.
  claims <- ascii(paste(cc$claim_text_candidate, collapse = " "))
  for (v in c("-0.39", "0.159", "0.156", "1/1001")) {
    expect_true(grepl(v, sec, fixed = TRUE), info = paste("draft missing", v))
    expect_true(grepl(v, claims, fixed = TRUE),
                info = paste("contract missing", v))
  }
  # n and composition come from the methods contract row FM-11
  n_row <- fm$value[fm$methods_id == "FM-11"]
  for (v in c("111", "58", "53", "24", "49", "38")) {
    expect_true(grepl(v, sec, fixed = TRUE))
    expect_true(grepl(v, n_row, fixed = TRUE))
  }

  # Thresholds. The frozen bundle carries them ROUNDED to three decimals; the
  # manuscript quotes full precision, which is licensed for reproducibility and
  # corroborated by the upstream canonical-definition test. Assert both that
  # Methods carries the full-precision form and that it is consistent with the
  # rounded value the bundle actually states - so a typo in the long form is
  # still caught.
  meth <- paste(d[nxt:length(d)], collapse = " ")
  full <- c(male = -0.436641698, female = -0.222390844)
  expect_true(grepl("-0.436641698", meth, fixed = TRUE))
  expect_true(grepl("-0.222390844", meth, fixed = TRUE))
  fm05 <- ascii(fm$value[fm$methods_id == "FM-05"])
  for (nm in names(full))
    expect_true(grepl(sprintf("%.3f", full[[nm]]), fm05, fixed = TRUE),
                info = paste("FM-05 does not corroborate the", nm, "threshold"))
})

test_that("the repeated-CV spread is not described as a confidence interval", {
  # Upstream records this as "2.5-97.5% quantiles across repeated CV splits" and
  # its generating script explicitly warns that the legacy ci aliases are not
  # confidence intervals. Calling it a CI would overstate it.
  skip_if_not(file.exists(DRAFT), "draft absent")
  d <- paste(readLines(DRAFT, warn = FALSE), collapse = " ")
  expect_true(grepl("percentile range", d, fixed = TRUE))
  expect_false(grepl("95% (confidence interval|CI)[^.]{0,60}0\\.116", d,
                     perl = TRUE))
})

# ----------------------------------------------------- provenance completeness

test_that("every Figure 1 claim resolves to a real frozen bundle row", {
  skip_if_not(dir.exists(BRIDGE), "figure 1 bridge not imported")
  cp <- rd(file.path(MS, "results_claim_provenance.csv"))
  f1 <- cp[grepl("^C1-", cp$claim_id), , drop = FALSE]
  expect_gt(nrow(f1), 0)
  expect_equal(anyDuplicated(f1$claim_id), 0L)
  cc <- rd(file.path(BRIDGE, "figure1_claim_contract.csv"))
  fm <- rd(file.path(BRIDGE, "figure1_methods_contract.csv"))
  known <- c(cc$claim_id, fm$methods_id)
  for (i in seq_len(nrow(f1))) {
    src <- paste(f1$primary_source[i], f1$secondary_source[i])
    # the row must name the bridge and at least one real F1-xx or FM-xx row
    expect_true(grepl("figure1_bridge_mmmsociability", src, fixed = TRUE))
    expect_true(any(vapply(known, function(k) grepl(k, src, fixed = TRUE),
                           logical(1))),
                info = paste(f1$claim_id[i], "names no known bundle row"))
  }
  # and every claim must carry a verification status naming the bridge
  expect_true(all(grepl("FROZEN_BRIDGE", f1$verification_status)))
})

test_that("superseded Phase-2B assertions are marked, not silently left true", {
  skip_if_not(dir.exists(BRIDGE), "figure 1 bridge not imported")
  m <- rd(file.path(MS, "methods_statement_provenance.csv"))
  for (id in c("M-16", "M-17", "M-18"))
    expect_true(grepl("SUPERSEDED", m$verified[m$methods_id == id]),
                info = paste(id, "still asserts the pre-bridge state"))
  # the contracts that declared Figure 1 unresolved must no longer say so
  bp <- rd(file.path(MS, "behavior_prediction_contract.csv"))
  expect_false(any(grepl("NO PREDICTION ANALYSIS EXISTS", bp$value)))
  os <- rd(file.path(MS, "outcome_score_contract.csv"))
  expect_false(any(grepl("\\[UNRESOLVED\\]", os$value)))
  # every point of disagreement is recorded rather than quietly resolved
  fc <- rd(file.path(MS, "figure1_bridge_conflicts.csv"))
  expect_gt(nrow(fc), 0)
  expect_true(all(nzchar(fc$authority_applied)))
  expect_true(all(nzchar(fc$resolution_in_manuscript)))
})

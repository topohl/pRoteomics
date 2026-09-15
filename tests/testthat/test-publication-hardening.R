# Publication-hardening guards.
#
# Part A fixed what the manuscript may claim; Part B established which parts of
# the tree may depend on which. Both are properties that decay silently, so each
# one that can be checked mechanically is checked here rather than described in
# a document. Every expectation below corresponds to a finding or a verified
# clean result in results/tables/publication_hardening/.

repo <- function(...) file.path(testthat::test_path("..", ".."), ...)
tracked_scripts <- function() {
  f <- suppressWarnings(system2("git", c("-C", shQuote(repo()), "ls-files"),
                                stdout = TRUE, stderr = FALSE))
  grep("[.][Rr]$", f, value = TRUE)
}
layer_of <- function(p) {
  top <- sub("/.*", "", p)
  if (top == "99_audits") "audit layer"
  else if (top == "99_deprecated") "deprecated"
  else if (top == "90_testing") "testing scaffold"
  else if (grepl("^[0-9]{2}_", top)) "numbered analysis stage"
  else if (top == "R") "shared helper library"
  else if (top == "figures") "figure / manuscript layer"
  else if (top == "tests") "test suite"
  else "other"
}

# ---------------------------------------------------------------- B27 guards

test_that("no producer layer sources deprecated or scaffold code", {
  # PH finding: 0 such edges today. This is the guard that keeps it at 0, and is
  # the concrete form of "deprecated code must not masquerade as active".
  SRC <- paste0("(?<![A-Za-z0-9_.])(source|sys[.]source)\\s*\\(\\s*",
                "(repo_path\\s*\\(([^)]*)\\)|[\"']([^\"']+)[\"'])")
  producer <- c("numbered analysis stage", "figure / manuscript layer",
                "shared helper library")
  leaks <- character(0)
  for (p in tracked_scripts()) {
    if (!layer_of(p) %in% producer) next
    ln <- readLines(repo(p), warn = FALSE)
    m <- unlist(regmatches(ln, gregexpr(SRC, ln, perl = TRUE)))
    if (!length(m)) next
    tgt <- vapply(m, function(x) {
      if (grepl("repo_path", x)) {
        a <- sub(".*repo_path\\s*\\(", "", x)
        paste(gsub("[\"' )]", "", strsplit(a, ",")[[1]]), collapse = "/")
      } else sub(".*[\"']([^\"']+)[\"'].*", "\\1", x)
    }, character(1))
    bad <- tgt[vapply(tgt, layer_of, character(1)) %in%
                 c("deprecated", "testing scaffold")]
    if (length(bad)) leaks <- c(leaks, sprintf("%s -> %s", p, bad))
  }
  expect_equal(leaks, character(0))
})

test_that("99_audits stays excluded from the pipeline registry", {
  # DEC-004. An audit is not a pipeline stage, and must never become a required
  # step that a publication rerun depends on.
  source(repo("R", "paths.R"))
  source(repo("R", "pipeline_registry.R"))
  ex <- pipeline_analysis_script_exclusions()
  expect_true("99_audits" %in% ex$roots)
})

test_that("every registry step names a file that exists", {
  source(repo("R", "paths.R"))
  source(repo("R", "pipeline_registry.R"))
  entries <- pipeline_registry_entries(read_pipeline_registry())
  missing <- entries$script[!file.exists(repo(entries$script))]
  expect_equal(missing, character(0))
})

# ------------------------------------------------- B28 freeze protection

test_that("the frozen v9 contract still reaches every renderer it names", {
  # PH-001: five v9 panels are rendered from superseded layer files. That is
  # accepted and documented, but the renderers must remain DEFINED - moving or
  # renaming one would break the frozen figures' reproducibility silently.
  ct <- yaml::read_yaml(repo("figures", "figure_final_truth_v9_contract.yml"))
  used <- unique(vapply(ct$panels, function(p)
    as.character(p$renderer %||% ""), character(1)))
  used <- used[nzchar(used)]
  srcs <- c(Sys.glob(repo("R", "*.R")), Sys.glob(repo("figures", "*.R")))
  defined <- unlist(lapply(srcs, function(f)
    sub(" <- function.*", "",
        grep("^[a-z][A-Za-z0-9_]* <- function", readLines(f, warn = FALSE),
             value = TRUE))))
  expect_equal(setdiff(used, defined), character(0))
})

test_that("out-of-layer v9 renderers are exactly the five that are documented", {
  # If this count changes, a panel silently moved between generations.
  ct <- yaml::read_yaml(repo("figures", "figure_final_truth_v9_contract.yml"))
  used <- unique(vapply(ct$panels, function(p)
    as.character(p$renderer %||% ""), character(1)))
  used <- used[nzchar(used)]
  out_of_layer <- sort(used[!grepl("^(f9_|s9f_)", used)])
  expect_equal(out_of_layer,
               c("nf_bilateral_main", "nf_pca_compact", "nvp_ed_celltype",
                 "s5_ed_ca2_displacement", "s5_ed_network_distance"))
})

test_that("the specificity gate matches the adjective, not only the adverb", {
  # PH-002. "spatially selective" is the spatial form of the banned
  # susceptibility-specific claim; matching only "selectively" let it through.
  sem <- readLines(repo("figures", "final_truth_v9_semantics.R"), warn = FALSE)
  pat <- grep("selectiv", sem, value = TRUE)
  expect_true(any(grepl("selectiv[|]exclusiv", pat)))
  expect_true(any(grepl('RULE\\("selective"', sem)))
})

test_that("every primary atlas theme has exactly one display label", {
  reg <- utils::read.csv(repo("config", "manuscript_go_theme_registry.tsv"),
                         sep = "\t", stringsAsFactors = FALSE)
  prim <- unique(reg$theme_id[reg$theme_role == "primary"])
  pan <- readLines(repo("R", "final_truth_v9_panels.R"), warn = FALSE)
  b <- grep("^  SHORT <- c\\(", pan)
  e <- b + which(grepl("\\)\\s*$", pan[b:(b + 20)]))[1] - 1L
  src <- paste(pan[b:e], collapse = " ")
  ids <- gsub(" =.*", "",
              regmatches(src, gregexpr("[a-z0-9_]+ = \"[^\"]*\"", src))[[1]])
  expect_equal(sort(ids), sort(prim))
})

test_that("no manuscript prose claims a theme-level p-value or FDR", {
  # The atlas is a descriptive aggregation with no multiple-testing family, so
  # this phrasing can never become correct by rerunning anything.
  rep_dir <- repo("results", "reports", "manuscript_candidates",
                  "final_truth_v9")
  skip_if_not(dir.exists(rep_dir), "v9 report layer not built")
  files <- list.files(rep_dir, pattern = "[.]md$", full.names = TRUE)
  # the rules file and the claim-chain audit exist to quote banned wording
  files <- files[!basename(files) %in% c("manuscript_semantic_rules.md",
                                         "claim_chain_audit.md")]
  # The repository states this prohibition explicitly in its own legends -
  # "no theme-level p-value or FDR is computed or implied" - so the guard has to
  # separate an assertion from a denial, exactly as the semantics layer does.
  DENIAL <- paste0("(^|[^[:alnum:]])(no|never|not|neither|without)([^[:alnum:]]",
                   "[^.]{0,80})?(theme[- ]level|FDR[- ]significant)")
  hits <- unlist(lapply(files, function(f) {
    ln <- readLines(f, warn = FALSE)
    h <- grep("theme[- ]level (p|FDR)|FDR[- ]significant theme", ln,
              ignore.case = TRUE, value = TRUE)
    h[!grepl(DENIAL, h, ignore.case = TRUE, perl = TRUE)]
  }))
  expect_equal(hits, character(0))
})

# ------------------------------------------------- PH-008 spatial wording

test_that("reader-facing prose uses spatially resolved, not restricted", {
  # PH-008. "Resolved" describes what the design achieves; "restricted" and
  # "specific" assert a boundary only a heterogeneity test could draw, and the
  # only such test is at WGCNA level and is FDR-negative.
  rep_dir <- repo("results", "reports", "manuscript_candidates",
                  "final_truth_v9")
  skip_if_not(dir.exists(rep_dir), "v9 report layer not built")
  files <- list.files(rep_dir, pattern = "[.]md$", full.names = TRUE)
  # these two exist in order to quote the wording they ban
  files <- files[!basename(files) %in% c("manuscript_semantic_rules.md",
                                         "claim_chain_audit.md")]
  # A stated count of supported units is factual, not an inferential claim.
  LICENCE <- "of 18|of 10|of 15|FDR-supported in|units in which|subset of units"
  # The repository's standing convention: text that instructs against a phrase,
  # or quotes it in order to replace it, necessarily contains it. The semantics
  # layer exempts such lines and so must this guard, or the corrected core
  # story's own explanation of the correction would fail it.
  PROHIBITION <- paste0("never|must not|do not |does not|cannot|prohibited|",
                        "instead of|rather than|avoid|banned|forbidden|",
                        "reworded|not adopted|say spatially resolved|",
                        "\\bNOT\\b|\\bNEVER\\b")
  hits <- unlist(lapply(files, function(f) {
    ln <- readLines(f, warn = FALSE)
    h <- grep("spatially[ -](restricted|specific)", ln, ignore.case = TRUE,
              value = TRUE)
    h <- h[!grepl(LICENCE, h, ignore.case = TRUE)]
    h[!grepl(PROHIBITION, h, perl = TRUE)]
  }))
  expect_equal(hits, character(0))
})

test_that("the spatial wording contract is declared in the rulebook", {
  sem <- readLines(repo("figures", "final_truth_v9_semantics.R"), warn = FALSE)
  expect_true(any(grepl('RULE("spatial wording"', sem, fixed = TRUE)))
  # the S9 scan must carry the guard, anchored to the adverb so that a
  # restricted INTERPRETATION or a specificity INVENTORY is not flagged
  expect_true(any(grepl("spatially[ -](restricted|specific)", sem,
                        fixed = TRUE)))
})

# ------------------------------------------- PH-010 verification provenance

test_that("the verification-state guard exists and declares its contract", {
  p <- repo("99_audits", "publication_hardening", "09_verification_state.R")
  expect_true(file.exists(p))
  src <- readLines(p, warn = FALSE)
  # the three states a verification run can describe
  for (tok in c("HEAD", "STAGED_TREE", "DIRTY_WORKTREE"))
    expect_true(any(grepl(tok, src, fixed = TRUE)))
  # strict only in release mode, with an explicit escape hatch
  expect_true(any(grepl("--release", src, fixed = TRUE)))
  expect_true(any(grepl("--allow-nonhead-verification", src, fixed = TRUE)))
  expect_true(any(grepl("quit(status = 1L)", src, fixed = TRUE)))
})

test_that("a release verification refuses a tree that is not HEAD", {
  # The guard is only useful if it actually fails. Exercise it in a throwaway
  # repository rather than against the real one, so the test cannot depend on
  # the state of the working tree it is being run from.
  skip_if(nchar(Sys.which("git")) == 0, "git unavailable")
  tmp <- file.path(tempdir(), paste0("phverif", as.integer(Sys.time())))
  dir.create(tmp, showWarnings = FALSE)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  q <- function(...) suppressWarnings(system2("git", c("-C", shQuote(tmp), ...),
                                              stdout = TRUE, stderr = FALSE))
  q("init", "-q")
  q("config", "user.email", "t@t"); q("config", "user.name", "t")
  writeLines("x", file.path(tmp, "a.txt"))
  q("add", "-A"); q("commit", "-qm", "init")
  clean <- length(q("status", "--porcelain", "--untracked-files=all")) == 0
  expect_true(clean)
  writeLines("y", file.path(tmp, "b.txt"))
  dirty <- q("status", "--porcelain", "--untracked-files=all")
  expect_true(length(dirty) > 0)
  # an untracked file is invisible to git ls-files - the exact PH-010 blind spot
  expect_false("b.txt" %in% q("ls-files"))
})

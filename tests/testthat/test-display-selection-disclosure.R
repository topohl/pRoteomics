# The display-selection release must not overstate what it supports.
#
# The six inventories are released as the DENOMINATORS behind selections shown
# in Figures 2 and 3. Two things make that release honest, and both are easy to
# lose in a later edit:
#
#   1. it says what it supports - protein membership and universe size;
#   2. it says what it does NOT - the per-protein GSEA rank statistic and the
#      original ranked-list order, neither of which reaches this analysis
#      because upstream supplies leading_edge_genes alphabetically sorted.
#
# Phase 6H.11 found the second point stated on two of the three disclosure
# surfaces and missing from the data dictionary, whose definition of
# leading_edge_size_of_term ended "Figure 3 g/h/i display 7 of these" - true,
# but an invitation to look for WHICH seven in a table that cannot say.
#
# These tests assert the semantic contract rather than exact prose, so
# rewording stays allowed and dropping the caveat does not.

source(testthat::test_path("..", "..", "R", "paths.R"))

TABLES <- path_results("integration", "build_display_selection_inventories",
                       "global", "tables")
BUNDLE <- repo_path("exports", "publication_source_data",
                    "supplementary_selection_inventories")
PRODUCER <- repo_path("analysis", "integration",
                      "build_display_selection_inventories.R")

# Does a piece of text affirm membership / universe-size support?
affirms_membership <- function(x) {
  grepl("membership", x, ignore.case = TRUE) &&
    grepl("universe size|set size|number of leading-edge", x, ignore.case = TRUE)
}
# Does it explicitly deny carrying the rank statistic AND the order?
denies_rank_and_order <- function(x) {
  neg <- "not preserve|does not carry|is not carried|cannot be reproduced|not present"
  grepl(neg, x, ignore.case = TRUE) &&
    grepl("rank statistic", x, ignore.case = TRUE) &&
    grepl("order", x, ignore.case = TRUE)
}

dictionary_row <- function(dir) {
  f <- file.path(dir, "inventory_data_dictionary.csv")
  if (!file.exists(f)) return(NULL)
  d <- utils::read.csv(f, stringsAsFactors = FALSE)
  d[d$table_file == "leading_edge_protein_inventory.csv" &
      d$column == "leading_edge_size_of_term", , drop = FALSE]
}

testthat::test_that("the producer defines the caveat, so regeneration keeps it", {
  # Fixing only the generated CSV would be undone by the next run. The
  # authoritative text lives in the producer.
  testthat::skip_if_not(file.exists(PRODUCER), "producer absent")
  src <- paste(readLines(PRODUCER, warn = FALSE), collapse = "\n")
  hit <- regmatches(src, regexpr("leading_edge_size_of_term = \"[^\"]+\"", src))
  testthat::expect_identical(length(hit), 1L)
  testthat::expect_true(affirms_membership(hit),
    info = "producer definition no longer states membership / universe size")
  testthat::expect_true(denies_rank_and_order(hit),
    info = "producer definition no longer denies the rank statistic and order")
})

testthat::test_that("the data dictionary states both halves of the contract", {
  for (dir in c(TABLES, BUNDLE)) {
    r <- dictionary_row(dir)
    testthat::skip_if_not(!is.null(r), paste("dictionary absent in", dir))
    testthat::expect_identical(nrow(r), 1L)
    testthat::expect_true(affirms_membership(r$definition[1]), info = dir)
    testthat::expect_true(denies_rank_and_order(r$definition[1]), info = dir)
  }
})

testthat::test_that("no disclosure surface claims the inventory carries rank or order", {
  # The panel's RULE may be described as ranking by the stored statistic - that
  # is a true statement about the figure. What must never appear is a claim
  # that these released tables let a reader reproduce it.
  surfaces <- c(file.path(TABLES, "inventory_data_dictionary.csv"),
                file.path(TABLES, "display_selection_disclosure.csv"),
                repo_path("docs", "FIGURE_SELECTION_RULES.md"))
  bad <- c("inventory (?:contains|carries|preserves|provides)[^.]*rank",
           "rank statistic[^.]*(?:released|included|carried) (?:in|by) (?:the )?inventor",
           "reproduce the (?:displayed )?ordering from these tables")
  for (f in surfaces) {
    testthat::skip_if_not(file.exists(f), paste(basename(f), "absent"))
    txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
    for (p in bad)
      testthat::expect_false(grepl(p, txt, ignore.case = TRUE, perl = TRUE),
        info = paste(basename(f), "matches an unsupported rank claim:", p))
  }
})

testthat::test_that("the three surfaces agree on what is and is not supported", {
  dict <- dictionary_row(TABLES)
  disc_f <- file.path(TABLES, "display_selection_disclosure.csv")
  rules_f <- repo_path("docs", "FIGURE_SELECTION_RULES.md")
  testthat::skip_if_not(!is.null(dict) && file.exists(disc_f) && file.exists(rules_f),
                        "a disclosure surface is absent")

  disc <- utils::read.csv(disc_f, stringsAsFactors = FALSE)
  ghi <- disc[grepl("3 g/h/i", disc[[1]]), , drop = FALSE]
  testthat::expect_identical(nrow(ghi), 1L)
  rules <- paste(readLines(rules_f, warn = FALSE), collapse = "\n")

  # all three deny rank/order; the disclosure CSV and dictionary do it in the
  # field a reader lands on, the rules doc in prose
  testthat::expect_true(denies_rank_and_order(dict$definition[1]))
  testthat::expect_true(denies_rank_and_order(
    paste(unlist(ghi), collapse = " ")))
  testthat::expect_true(grepl("alphabetically", rules, ignore.case = TRUE) ||
                        denies_rank_and_order(rules))
})

testthat::test_that("the released scientific tables are unchanged by documentation edits", {
  # A dictionary reword must never move a number. These are the row counts the
  # Phase 6H closure accepted.
  expected <- c(pathway_enrichment_inventory.csv = 203073L,
                leading_edge_protein_inventory.csv = 282296L,
                leading_edge_protein_recurrence.csv = 12227L,
                pathway_enrichment_inventory_fdr_supported.csv = 3559L,
                inventory_data_dictionary.csv = 72L,
                display_selection_disclosure.csv = 7L)
  man <- repo_path("exports", "publication_source_data", "manifest.csv")
  testthat::skip_if_not(file.exists(man), "export manifest absent")
  m <- utils::read.csv(man, stringsAsFactors = FALSE)
  s <- m[m$publication_id == "supplementary_selection_inventories", , drop = FALSE]
  testthat::expect_identical(nrow(s), 6L)
  for (nm in names(expected)) {
    row <- s[basename(s$exported_file) == nm, , drop = FALSE]
    testthat::expect_identical(nrow(row), 1L, info = nm)
    testthat::expect_identical(as.integer(row$rows[1]), expected[[nm]], info = nm)
  }
  # and the bundle stays the accepted shape
  testthat::expect_identical(nrow(m), 60L)
  testthat::expect_identical(length(unique(m$publication_id)), 11L)
})

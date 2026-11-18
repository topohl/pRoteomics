source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "protigy_input_utils.R"))

test_that("protein groups without a mouse protein are left out, as in WGCNA", {
  expression <- data.frame(
    `T: Protein.Names` = c("ALBU_MOUSE", "K1C10_HUMAN", "TRYP_PIG;ALBU_MOUSE", "KRT34_HUMAN;K1C14_HUMAN",
                           NA, "", "tbb5_mouse", "ALBU_BOVIN"),
    S001 = 1:8,
    check.names = FALSE, stringsAsFactors = FALSE
  )
  res <- protigy_drop_non_mouse_groups(expression)
  expect_equal(res$expression[["T: Protein.Names"]], c("ALBU_MOUSE", "TRYP_PIG;ALBU_MOUSE", NA, "", "tbb5_mouse"))
  expect_equal(res$dropped$protein_names, c("K1C10_HUMAN", "KRT34_HUMAN;K1C14_HUMAN", "ALBU_BOVIN"))
  expect_equal(res$dropped$source_row, c(2L, 4L, 8L))
  expect_equal(res$expression$S001, c(1L, 3L, 5L, 6L, 7L))
})

test_that("the plain Protein.Names column works and a missing name column stops", {
  ok <- protigy_drop_non_mouse_groups(data.frame(Protein.Names = c("A_MOUSE", "B_HUMAN"), S1 = 1:2))
  expect_equal(nrow(ok$expression), 1L)
  expect_error(protigy_drop_non_mouse_groups(data.frame(id = "A_MOUSE", S1 = 1)), "protein-name column")
})

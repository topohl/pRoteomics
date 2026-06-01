testthat::test_that("path helpers resolve inside repository", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  root <- repo_root()
  testthat::expect_true(file.exists(file.path(root, "README.md")))
  testthat::expect_match(repo_path("R", "paths.R"), "R[/\\\\]paths[.]R$")
  testthat::expect_equal(relative_to(repo_path("README.md")), "README.md")
  testthat::expect_match(safe_filename("a/b c"), "a_b_c")
})

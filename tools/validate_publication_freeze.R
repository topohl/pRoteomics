#!/usr/bin/env Rscript

# Validate the current checkout against docs/publication_freeze_manifest.yml.
#
# Strictly read-only: it hashes existing files and compares them with the
# recorded freeze identity. It never writes, moves or regenerates anything.
#
# Exit status:
#   0  no FAIL checks (WARN for documented accepted gaps is not fatal)
#   1  one or more FAIL checks
#
# Usage:
#   Rscript tools/validate_publication_freeze.R

paths_file <- if (file.exists(file.path("R", "paths.R"))) {
  file.path("R", "paths.R")
} else {
  file.path("..", "R", "paths.R")
}
source(paths_file)
source(repo_path("R", "publication_freeze_utils.R"))

result <- validate_publication_freeze()

cat("Publication freeze validation:", result$manifest_path, "\n\n")
print_publication_freeze_validation(result)

if (result$summary[["FAIL"]] > 0L) {
  cat("\nFREEZE VALIDATION FAILED\n")
  quit(status = 1, save = "no")
}
cat("\nFREEZE VALIDATION OK",
    if (result$summary[["WARN"]] > 0L) {
      paste0(" (", result$summary[["WARN"]], " documented accepted gap(s))")
    } else "", "\n", sep = "")
quit(status = 0, save = "no")

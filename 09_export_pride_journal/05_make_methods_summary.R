#!/usr/bin/env Rscript
# Deprecated path — use 06_make_methods_summary.R
args <- commandArgs(trailingOnly = TRUE)
script <- if (file.exists(file.path("R", "paths.R"))) {
  file.path("09_export_pride_journal", "06_make_methods_summary.R")
} else {
  file.path("..", "09_export_pride_journal", "06_make_methods_summary.R")
}
status <- system2("Rscript", c(script, args))
quit(status = if (is.null(status)) 0 else status, save = "no")

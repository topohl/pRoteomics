#!/usr/bin/env Rscript
# Deprecated path — use 09_export_source_data.R
args <- commandArgs(trailingOnly = TRUE)
script <- if (file.exists(file.path("R", "paths.R"))) {
  file.path("09_export_pride_journal", "09_export_source_data.R")
} else {
  file.path("..", "09_export_pride_journal", "09_export_source_data.R")
}
status <- system2("Rscript", c(script, args))
quit(status = if (is.null(status)) 0 else status, save = "no")

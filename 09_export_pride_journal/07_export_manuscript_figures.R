#!/usr/bin/env Rscript
# Deprecated path — use 08_export_manuscript_figures.R
args <- commandArgs(trailingOnly = TRUE)
script <- if (file.exists(file.path("R", "paths.R"))) {
  file.path("09_export_pride_journal", "08_export_manuscript_figures.R")
} else {
  file.path("..", "09_export_pride_journal", "08_export_manuscript_figures.R")
}
status <- system2("Rscript", c(script, args))
quit(status = if (is.null(status)) 0 else status, save = "no")

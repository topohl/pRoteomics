#!/usr/bin/env Rscript

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)

message("Legacy entrypoint: 06_modules_WGCNA/91_module_score.r. Use 06_modules_WGCNA/03_score_module_activity.R instead.")
source(repo_path("06_modules_WGCNA", "03_score_module_activity.R"))

#!/usr/bin/env Rscript

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)

message("Deprecated entrypoint: 06_modules_WGCNA/04_overlap_modules.r. Use 06_modules_WGCNA/02_curated_overlap_programs.r instead.")
source(repo_path("06_modules_WGCNA", "02_curated_overlap_programs.r"))

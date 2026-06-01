#!/usr/bin/env Rscript

# Publication-facing alias for the active module activity scoring workflow.
# The implementation remains in 91_module_score.r for backward compatibility.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("06_modules_WGCNA", "91_module_score.r"))

# Compatibility wrapper for the active overlap-module script.
#
# RUN_ORDER.md now points to 02_curated_overlap_programs.r, but this wrapper preserves
# older overlap-program calls without duplicating analysis code.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("06_modules_WGCNA", "02_curated_overlap_programs.r"))

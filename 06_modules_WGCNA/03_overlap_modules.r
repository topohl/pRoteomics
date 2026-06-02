# Deprecated compatibility wrapper.
#
# The overlap-program producer now lives in 06_modules_WGCNA/04_overlap_modules.r
# and writes global-scoped curated_overlap_programs artifacts.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)

message(
  "03_overlap_modules.r is deprecated. ",
  "Forwarding to 06_modules_WGCNA/04_overlap_modules.r, which writes ",
  "curated_overlap_programs under the explicit global output path."
)

source(repo_path("06_modules_WGCNA", "04_overlap_modules.r"))

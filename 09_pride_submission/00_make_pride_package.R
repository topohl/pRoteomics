# Assemble the local PRIDE/ProteomeXchange package skeleton and metadata files.
#
# Run from the repository root:
#   source("09_pride_submission/00_make_pride_package.R")

source("R/paths.R")
source("R/pride_helpers.R")

ensure_pride_dirs()

source("09_pride_submission/01_generate_sdrf.R")
source("09_pride_submission/02_generate_pride_manifest.R")
source("09_pride_submission/03_generate_checksums.R")
source("09_pride_submission/04_validate_pride_package.R")

message("PRIDE package metadata generation complete.")
message("Review: ", repo_path("PRIDE_package", "00_metadata"))

# Generate checksum.md5 for the local PRIDE package.
#
# Checksums are used to verify file integrity after copying or upload.

source("R/paths.R")
source("R/pride_helpers.R")

ensure_pride_dirs()
package_dir <- pride_package_dir()
out_file <- repo_path("PRIDE_package", "00_metadata", "checksum.md5")

files <- list_pride_files(package_dir)
files <- files[basename(files) != "checksum.md5"]

md5 <- make_md5_table(files, package_dir)

if (nrow(md5) == 0) {
  cat("", file = out_file)
} else {
  # Standard md5sum-style two-column format: checksum, relative path.
  lines <- paste(md5$md5, md5$relative_path)
  writeLines(lines, out_file)
}

message("Wrote checksum file: ", out_file)

# Generate a PRIDE package file manifest by scanning PRIDE_package/.
#
# The manifest is a local quality-control object. It is not a substitute for
# PRIDE's submission metadata, but it makes the package traceable.

source("R/paths.R")
source("R/pride_helpers.R")

ensure_pride_dirs()
package_dir <- pride_package_dir()
metadata_dir <- repo_path("PRIDE_package", "00_metadata")
out_file <- file.path(metadata_dir, "pride_file_manifest.tsv")

files <- list_pride_files(package_dir)
files <- files[basename(files) != "pride_file_manifest.tsv"]

if (length(files) == 0) {
  manifest <- data.frame(
    file_name = character(),
    relative_path = character(),
    pride_file_type = character(),
    description = character(),
    generated_by_script = character(),
    sample_id = character(),
    raw_file_name = character(),
    software = character(),
    include_in_pride = logical(),
    md5 = character(),
    stringsAsFactors = FALSE
  )
} else {
  rel <- relative_to(files, package_dir)
  md5 <- make_md5_table(files, package_dir)
  manifest <- data.frame(
    file_name = basename(files),
    relative_path = rel,
    pride_file_type = vapply(rel, classify_pride_file, character(1)),
    description = "",
    generated_by_script = "",
    sample_id = "",
    raw_file_name = "",
    software = "",
    include_in_pride = TRUE,
    stringsAsFactors = FALSE
  )
  manifest <- merge(manifest, md5, by = "relative_path", all.x = TRUE, sort = FALSE)
  manifest <- manifest[, c(
    "file_name", "relative_path", "pride_file_type", "description",
    "generated_by_script", "sample_id", "raw_file_name", "software",
    "include_in_pride", "md5"
  )]
}

write_tsv(manifest, out_file)
message("Wrote PRIDE file manifest: ", out_file)

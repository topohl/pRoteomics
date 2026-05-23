# Validates PRIDE/journal staging folders and manifest presence.
# Consumes pride_submission/; produces pride_submission/validation/validation_report.tsv.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)

ensure_pride_dirs()
checks <- data.frame(check = character(), status = character(), detail = character(), stringsAsFactors = FALSE)
add_check <- function(name, ok, detail) {
  checks <<- rbind(checks, data.frame(check = name, status = if (ok) "PASS" else "FAIL", detail = detail, stringsAsFactors = FALSE))
}

manifest <- pride_submission_dir("manifests", "pride_file_manifest.tsv")
add_check("manifest_exists", file.exists(manifest), manifest)
add_check("sample_metadata_exists", length(list.files(pride_submission_dir("metadata"), pattern = "sample_metadata|sdrf", full.names = TRUE)) > 0, pride_submission_dir("metadata"))
add_check("processed_data_dir_exists", dir.exists(pride_submission_dir("processed_data")), pride_submission_dir("processed_data"))
add_check("methods_dir_exists", dir.exists(pride_submission_dir("methods")), pride_submission_dir("methods"))

if (file.exists(manifest)) {
  mf <- utils::read.delim(manifest, stringsAsFactors = FALSE)
  add_check("manifest_has_md5", "md5" %in% names(mf) && all(nzchar(mf$md5)), "md5 column populated")
  add_check("manifest_has_pride_flags", "intended_for_PRIDE" %in% names(mf), "intended_for_PRIDE column")
}

out_file <- pride_submission_dir("validation", "validation_report.tsv")
utils::write.table(checks, out_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
message("Wrote validation report: ", out_file)

# Creates pride_submission/manifests/pride_file_manifest.tsv.
# Consumes canonical data/results/pride_submission files; produces a PRIDE/journal file manifest.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)

ensure_pride_dirs()

classify_export_file <- function(path) {
  p <- tolower(gsub("\\\\", "/", path))
  ext <- tolower(tools::file_ext(p))
  if (ext %in% c("raw", "wiff", "mzml", "mzxml", "d", "tdf", "baf")) return("raw_ms")
  if (grepl("matrix|quant|imputed|normalized|processed", p)) return("processed_quantitative_matrix")
  if (grepl("differential|contrast|log2fc|logfc", p)) return("differential_expression_result")
  if (grepl("enrichment|clusterprofiler|comparego|ewce|wgcna|network", p)) return("secondary_analysis")
  if (grepl("metadata|sdrf|sample", p)) return("metadata")
  if (grepl("method|sessioninfo|config|manifest", p)) return("provenance")
  if (ext %in% c("svg", "png", "pdf")) return("figure_only")
  "other"
}

candidate_roots <- c(
  pride_submission_dir(),
  path_processed(),
  path_results("tables"),
  path_results("source_data"),
  path_results("reports")
)
candidate_roots <- candidate_roots[dir.exists(candidate_roots)]
files <- unique(unlist(lapply(candidate_roots, list.files, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE)))
files <- files[file.exists(files) & !file.info(files)$isdir]

manifest <- data.frame(
  file_path = relative_to(files, repo_root()),
  file_type = vapply(files, classify_export_file, character(1)),
  intended_for_PRIDE = vapply(files, function(f) classify_export_file(f) %in% c("raw_ms", "processed_quantitative_matrix", "metadata", "provenance"), logical(1)),
  intended_for_supplement = vapply(files, function(f) classify_export_file(f) %in% c("differential_expression_result", "secondary_analysis", "metadata"), logical(1)),
  md5 = unname(tools::md5sum(files)),
  size_bytes = file.info(files)$size,
  stringsAsFactors = FALSE
)

out_file <- pride_submission_dir("manifests", "pride_file_manifest.tsv")
utils::write.table(manifest, out_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
message("Wrote PRIDE/journal manifest: ", out_file)

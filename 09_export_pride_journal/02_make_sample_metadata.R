# Creates sample metadata and SDRF-like tables where feasible.
# Consumes data/metadata/sample_metadata_clean.{tsv,csv}; produces pride_submission/metadata/.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(file.path(repo_root(), "R", "pride_helpers.R"))

ensure_pride_dirs()
meta <- read_sample_metadata()
if (nrow(meta) == 0) stop("No sample metadata found in data/metadata or legacy candidate paths.", call. = FALSE)

out_clean <- pride_submission_dir("metadata", "sample_metadata.tsv")
write_tsv(meta, out_clean)

sdrf <- data.frame(
  source_name = coalesce_column(meta, c("sample_id", "sample", "mouse_id", "raw_file_name")),
  characteristics_organism = coalesce_column(meta, c("organism"), "Mus musculus"),
  characteristics_strain = coalesce_column(meta, c("strain"), NA_character_),
  characteristics_sex = coalesce_column(meta, c("sex"), NA_character_),
  characteristics_phenotype = coalesce_column(meta, c("phenotype", "group", "condition"), NA_character_),
  assay_name = coalesce_column(meta, c("raw_file_name", "file", "filename")),
  comment_data_file = coalesce_column(meta, c("raw_file_name", "file", "filename")),
  stringsAsFactors = FALSE
)
out_sdrf <- pride_submission_dir("metadata", "sdrf_like_metadata.tsv")
write_tsv(sdrf, out_sdrf)
message("Wrote sample metadata: ", out_clean)
message("Wrote SDRF-like metadata: ", out_sdrf)

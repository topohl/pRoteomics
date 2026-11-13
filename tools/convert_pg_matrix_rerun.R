#!/usr/bin/env Rscript
# Turn a facility DIA-NN delivery (protein-group matrix and run statistics) into the raw inputs of this design rerun
# tree, under the file names the pipeline reads.
#
# Usage, from the tree root:
#   Rscript tools/convert_pg_matrix_rerun.R <report.pg_matrix.tsv|.xlsx> <report.stats.tsv|.xlsx> [--check]
#
# --check validates the delivery and reports the run mapping without writing anything.
#
# Rules
# - A TSV must be DIA-NN's own output. Every value cell has to be empty or a plain number with a decimal point and,
#   for large values, a lower-case exponent (1.42367e+06). A file that was opened and saved in a spreadsheet program
#   carries decimal commas, capital-E numbers rounded to three digits or dates, and is refused.
# - Runs are matched to the metadata by sample number plus acquisition run number (S113 + 13026). That key is unique
#   and does not depend on the path prefix or on the layer token, which the facility corrected for A0003 S129-S132.
# - Runs of animals outside the design (DROP_ANIMALS) are dropped and recorded.
# - The metadata workbook is trimmed to the delivered runs, because the joint QC needs a one-to-one match between
#   matrix runs and metadata rows. The full workbook is read from BASE_COMMIT, so the script can be rerun on a later
#   delivery.
#
# Outputs (file names kept because the pipeline and several scripts name them):
#   data/raw/pg_matrix/quicksearch.pg_matrix.tsv            run columns renamed to the metadata sample_id
#   data/raw/pg_matrix/quicksearch.stats.annotated.xlsx     16 metadata columns plus the DIA-NN run statistics
#   data/metadata/TPE9_sample_metadata_males.xlsx           trimmed to the analysed runs
#   data/metadata/pg_matrix_rerun_run_map.csv               one row per metadata row and delivered run, with action
#   data/raw/pg_matrix_rerun/                               the delivery, copied verbatim, plus SHA256SUMS and
#                                                           README_PROVENANCE.txt

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
})

BASE_COMMIT <- "1ebc175"
METADATA_REL <- "data/metadata/TPE9_sample_metadata_males.xlsx"
DROP_ANIMALS <- "A764"
ANNOTATION <- c("Protein.Group", "Protein.Names", "Genes", "First.Protein.Description")
NUMBER_RE <- "^([0-9]+([.][0-9]+)?(e[+-]?[0-9]+)?)?$"
RUN_KEY_RE <- "_(S[0-9]{3})_S[0-9]-[A-H][0-9]{1,2}_[0-9]+_([0-9]+)[.]d$"

args <- commandArgs(trailingOnly = TRUE)
check_only <- "--check" %in% args
args <- setdiff(args, "--check")
if (length(args) != 2L) stop("usage: convert_pg_matrix_rerun.R <pg_matrix> <stats> [--check]", call. = FALSE)
matrix_file <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
stats_file <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
if (!file.exists("pipeline.yml") || !file.exists(METADATA_REL)) stop("run from the tree root", call. = FALSE)

sha256 <- function(path) digest::digest(file = path, algo = "sha256")
problems <- character()

# Reads a delivery file as text so that nothing is reformatted before it is checked.
read_delivery <- function(path) {
  if (grepl("[.]xlsx$", path, ignore.case = TRUE)) {
    x <- as.data.frame(read_excel(path, sheet = 1, guess_max = 1e6), check.names = FALSE, stringsAsFactors = FALSE)
    attr(x, "format") <- "xlsx"
    return(x)
  }
  x <- read_tsv(path, col_types = cols(.default = col_character()), na = character(), trim_ws = FALSE,
                progress = FALSE)
  x <- as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)
  attr(x, "format") <- "tsv"
  x
}

# Checks DIA-NN number formatting in text columns and returns them as numbers.
as_checked_numbers <- function(x, cols, label) {
  out <- x
  n_bad <- 0L
  n_cells <- 0L
  bad_cols <- character()
  examples <- character()
  for (cl in cols) {
    v <- x[[cl]]
    if (is.numeric(v)) next
    if (is.logical(v) && all(is.na(v))) { out[[cl]] <- as.numeric(v); next }
    bad <- !grepl(NUMBER_RE, v)
    n_cells <- n_cells + sum(nzchar(v))
    if (any(bad)) {
      n_bad <- n_bad + sum(bad)
      bad_cols <- c(bad_cols, cl)
      examples <- c(examples, utils::head(unique(v[bad]), 2))
    }
    v[!nzchar(v)] <- NA_character_
    out[[cl]] <- suppressWarnings(as.numeric(v))
  }
  if (n_bad) {
    problems <<- c(problems, sprintf(
      "%s: %d of %d filled cells in %d of %d columns are not DIA-NN numbers (spreadsheet-saved file?), e.g. %s",
      label, n_bad, n_cells, length(bad_cols), length(cols), paste(utils::head(unique(examples), 6), collapse = " | ")))
  }
  out
}

run_key <- function(run) {
  m <- regmatches(run, regexec(RUN_KEY_RE, run))
  vapply(m, function(z) if (length(z) == 3L) paste0(z[[2]], "|", z[[3]]) else NA_character_, character(1))
}

# Full metadata of the base commit.
tmp_meta <- tempfile(fileext = ".xlsx")
status <- system2("git", c("cat-file", "blob", paste0(BASE_COMMIT, ":", METADATA_REL)), stdout = tmp_meta)
if (!identical(as.integer(status), 0L)) stop("could not read the metadata workbook of ", BASE_COMMIT, call. = FALSE)
meta_full <- read_excel(tmp_meta)
if (nrow(meta_full) != 326L || ncol(meta_full) != 16L) stop("unexpected base metadata shape", call. = FALSE)
meta_full$key <- paste0(meta_full$sample_number, "|", meta_full$run_order)
if (anyDuplicated(meta_full$key)) stop("sample number plus run number is not unique in the metadata", call. = FALSE)
if (!identical(run_key(meta_full$sample_id), meta_full$key)) stop("metadata sample_id and key columns disagree", call. = FALSE)

# Matrix.
pg <- read_delivery(matrix_file)
if (!identical(names(pg)[1:4], ANNOTATION)) stop("the first four matrix columns are not ", paste(ANNOTATION, collapse = ", "), call. = FALSE)
runs <- names(pg)[-(1:4)]
pg <- as_checked_numbers(pg, runs, "matrix")
if (anyDuplicated(pg$Protein.Group)) problems <- c(problems, "matrix: duplicated Protein.Group")

# Statistics.
st <- read_delivery(stats_file)
if (!identical(names(st)[1], "File.Name")) stop("the stats file does not start with File.Name", call. = FALSE)
st <- as_checked_numbers(st, names(st)[-1], "stats")

# Run mapping.
st_runs <- basename(gsub("\\\\", "/", st$File.Name))
pg_runs <- basename(gsub("\\\\", "/", runs))
if (!setequal(st_runs, pg_runs)) problems <- c(problems, "stats File.Name and matrix runs differ")
keys <- run_key(pg_runs)
if (anyNA(keys)) stop("runs without a sample/run-number key: ", paste(pg_runs[is.na(keys)], collapse = ", "), call. = FALSE)
if (anyDuplicated(keys)) stop("duplicated run keys in the matrix", call. = FALSE)
unknown <- setdiff(keys, meta_full$key)
if (length(unknown)) stop("runs not in the metadata: ", paste(pg_runs[keys %in% unknown], collapse = ", "), call. = FALSE)

delivered <- meta_full$key %in% keys
dropped <- meta_full$AnimalID %in% DROP_ANIMALS
excluded <- meta_full$exclude %in% TRUE
action <- ifelse(delivered & dropped, "dropped_not_in_design",
          ifelse(delivered & excluded, "delivered_excluded_by_metadata",
          ifelse(delivered, "analysed",
          ifelse(dropped, "not_delivered_not_in_design",
          ifelse(excluded, "not_delivered_excluded_by_metadata", "MISSING_from_delivery")))))
run_map <- data.frame(
  facility_run_name = pg_runs[match(meta_full$key, keys)],
  metadata_sample_id = meta_full$sample_id,
  sample_number = meta_full$sample_number,
  run_number = meta_full$run_order,
  AnimalID = meta_full$AnimalID,
  ExpGroup = meta_full$ExpGroup,
  ReplicateGroup = meta_full$ReplicateGroup,
  region = meta_full$region,
  metadata_layer = meta_full$layer,
  celltype_layer = meta_full$celltype_layer,
  exclude = meta_full$exclude,
  action = action,
  stringsAsFactors = FALSE
)
cat("Delivery:", basename(matrix_file), "(", attr(pg, "format"), ")", nrow(pg), "protein groups x", length(runs), "runs\n")
print(table(action))
missing <- run_map[run_map$action == "MISSING_from_delivery", c("metadata_sample_id", "AnimalID")]
if (nrow(missing)) cat("Runs in the design but not in the delivery:\n", paste0("  ", basename(gsub("\\\\", "/", missing$metadata_sample_id)), collapse = "\n"), "\n")

if (length(problems)) {
  cat("\nREFUSED:\n", paste0("  ", problems, collapse = "\n"), "\n")
  quit(status = 1)
}
cat("Format check passed.\n")
if (check_only) quit(status = 0)

# Outputs.
keep_meta <- meta_full[delivered & !dropped, ]
keep_meta$key <- NULL
keep_keys <- paste0(keep_meta$sample_number, "|", keep_meta$run_order)
run_cols <- runs[match(keep_keys, keys)]
out <- pg[, c(ANNOTATION, run_cols)]
names(out) <- c(ANNOTATION, keep_meta$sample_id)
for (cl in ANNOTATION) {
  if (any(grepl("[\t\r\n\"]", out[[cl]]))) stop("annotation column ", cl, " holds tabs, quotes or line breaks", call. = FALSE)
}

dir.create("data/raw/pg_matrix", recursive = TRUE, showWarnings = FALSE)
dir.create("data/raw/pg_matrix_rerun", recursive = TRUE, showWarnings = FALSE)
tsv_out <- "data/raw/pg_matrix/quicksearch.pg_matrix.tsv"
data.table::fwrite(out, tsv_out, sep = "\t", na = "", quote = FALSE)

# Round trip, read the way 01_impute.r and the joint QC read it.
back <- read_tsv(tsv_out, locale = locale(decimal_mark = "."), show_col_types = FALSE, progress = FALSE)
back_m <- as.matrix(back[, -(1:4)])
src_m <- as.matrix(out[, -(1:4)])
storage.mode(back_m) <- "double"
if (!identical(names(back), names(out)) || !identical(unname(back_m), unname(src_m)) ||
    !identical(back$Protein.Group, out$Protein.Group)) stop("round trip of the matrix failed (readr)", call. = FALSE)
back2 <- utils::read.delim(tsv_out, check.names = FALSE, stringsAsFactors = FALSE)
if (!identical(unname(as.matrix(back2[, -(1:4)])), unname(src_m))) stop("round trip of the matrix failed (read.delim)", call. = FALSE)

# Annotated statistics in metadata order.
st$File.Name <- NULL
st_keys <- run_key(st_runs)
st_out <- cbind(as.data.frame(keep_meta, stringsAsFactors = FALSE), File.Name = keep_meta$sample_id,
                st[match(keep_keys, st_keys), , drop = FALSE], stringsAsFactors = FALSE)
rownames(st_out) <- NULL
writexl::write_xlsx(list(quicksearch.stats = st_out), "data/raw/pg_matrix/quicksearch.stats.annotated.xlsx")

# Trimmed metadata, checked after reading back.
writexl::write_xlsx(list(Sheet1 = as.data.frame(keep_meta)), METADATA_REL)
meta_back <- read_excel(METADATA_REL)
if (!isTRUE(all.equal(as.data.frame(meta_back), as.data.frame(keep_meta), check.attributes = FALSE)) ||
    !identical(sapply(meta_back, class), sapply(keep_meta, class))) stop("metadata round trip failed", call. = FALSE)
utils::write.csv(run_map, "data/metadata/pg_matrix_rerun_run_map.csv", row.names = FALSE, na = "")

# Delivery copy and provenance.
copies <- file.path("data/raw/pg_matrix_rerun", basename(c(matrix_file, stats_file)))
if (!all(file.copy(c(matrix_file, stats_file), copies, overwrite = TRUE, copy.date = TRUE))) stop("copy failed", call. = FALSE)
written <- c(copies, tsv_out, "data/raw/pg_matrix/quicksearch.stats.annotated.xlsx", METADATA_REL,
             "data/metadata/pg_matrix_rerun_run_map.csv")
sums <- vapply(written, sha256, character(1))
if (!identical(unname(sums[1:2]), unname(vapply(c(matrix_file, stats_file), sha256, character(1))))) {
  stop("delivery copy differs from the source", call. = FALSE)
}
writeLines(paste0(sums, "  ", written), "data/raw/pg_matrix_rerun/SHA256SUMS")
tab <- table(action)
writeLines(c(
  "Raw inputs of the 3/3/2 design rerun",
  "",
  paste0("Converted ", format(Sys.time(), "%Y-%m-%d %H:%M"), " by tools/convert_pg_matrix_rerun.R from:"),
  paste0("  ", matrix_file, "  sha256 ", sha256(matrix_file)),
  paste0("  ", stats_file, "  sha256 ", sha256(stats_file)),
  "",
  "data/raw/pg_matrix/quicksearch.pg_matrix.tsv is this delivery, not the old quick search. The file name is kept",
  "because the pipeline and several scripts name it. Run columns carry the metadata sample_id; runs were matched by",
  "sample number plus acquisition run number (data/metadata/pg_matrix_rerun_run_map.csv).",
  "",
  paste0("Protein groups: ", nrow(out), "; runs in the delivery: ", length(runs), "; runs written: ", length(run_cols)),
  paste0("  ", names(tab), ": ", as.integer(tab)),
  if (nrow(missing)) c("Runs in the design but not in the delivery:",
                       paste0("  ", basename(gsub("\\\\", "/", missing$metadata_sample_id)))) else character()
), "data/raw/pg_matrix_rerun/README_PROVENANCE.txt")
cat("Written:", paste(written, collapse = ", "), "\n")

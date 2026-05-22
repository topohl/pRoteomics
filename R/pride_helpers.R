# Helper functions for PRIDE/ProteomeXchange package assembly.
#
# The functions are deliberately conservative: when information cannot be
# inferred, they keep NA values and let the validation report flag the problem.

standardize_names <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

read_delimited_auto <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    return(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE))
  }
  utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  invisible(path)
}

first_existing <- function(paths) {
  hits <- paths[file.exists(paths)]
  if (length(hits) == 0) return(NA_character_)
  hits[[1]]
}

find_metadata_file <- function(root = repo_root()) {
  candidates <- file.path(root, c(
    "results/analysis_ready/sample_metadata_clean.tsv",
    "results/analysis_ready/sample_metadata_clean.csv",
    "data/metadata/sample_metadata_clean.tsv",
    "data/metadata/sample_metadata_clean.csv",
    "01_preprocessing/output/sample_metadata_clean.tsv",
    "01_preprocessing/output/sample_metadata_clean.csv",
    "01_preprocessing/sample_metadata_clean.tsv",
    "01_preprocessing/sample_metadata_clean.csv"
  ))
  first_existing(candidates)
}

read_sample_metadata <- function(path = NULL) {
  if (is.null(path)) path <- find_metadata_file()
  if (is.na(path) || !file.exists(path)) {
    return(data.frame())
  }
  meta <- read_delimited_auto(path)
  names(meta) <- standardize_names(names(meta))
  meta
}

coalesce_column <- function(df, candidates, default = NA_character_) {
  candidates <- standardize_names(candidates)
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) {
    return(rep(default, nrow(df)))
  }
  out <- df[[hit[[1]]]]
  if (length(hit) > 1) {
    for (h in hit[-1]) {
      missing <- is.na(out) | out == ""
      out[missing] <- df[[h]][missing]
    }
  }
  out
}

classify_pride_file <- function(path) {
  p <- tolower(gsub("\\\\", "/", path))
  ext <- tolower(tools::file_ext(p))

  if (grepl("/01_raw/", p) || ext %in% c("raw", "wiff", "mzml", "mzxml", "d", "lcd", "baf", "tdf")) return("RAW")
  if (ext %in% c("mgf", "ms2", "mzdata")) return("PEAK")
  if (grepl("fasta|uniprot|database", p) || ext %in% c("fasta", "fa")) return("FASTA")
  if (grepl("parameter|settings|mqpar|search", p) && ext %in% c("txt", "xml", "json", "tsv", "csv")) return("PARAMETERS_FILE")
  if (grepl("sdrf|sample_metadata|experimental_design|manifest", p)) return("EXPERIMENTAL_DESIGN")
  if (grepl("search_results|proteinGroups|peptides|evidence|psm|msms|diann|spectronaut|maxquant", p, ignore.case = TRUE)) return("SEARCH")
  if (grepl("processed|quantification|analysis_ready|matrix|results|enrichment|ewce|wgcna|module|network", p)) return("RESULT")
  "OTHER"
}

list_pride_files <- function(package_dir = pride_package_dir()) {
  if (!dir.exists(package_dir)) return(character())
  files <- list.files(package_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE)
  files[file.info(files)$isdir == FALSE]
}

relative_to <- function(path, root) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  sub(paste0("^", gsub("([\\^$.|?*+(){}])", "\\\\\\1", root), "/?"), "", path)
}

make_md5_table <- function(files, root = pride_package_dir()) {
  if (length(files) == 0) {
    return(data.frame(md5 = character(), relative_path = character()))
  }
  sums <- tools::md5sum(files)
  data.frame(
    md5 = unname(sums),
    relative_path = relative_to(names(sums), root),
    stringsAsFactors = FALSE
  )
}

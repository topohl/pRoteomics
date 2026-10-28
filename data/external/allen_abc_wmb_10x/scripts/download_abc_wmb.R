# Download the ABC Atlas whole-mouse-brain (WMB-10X, Yao et al. 2023) files needed for
# per-subclass / per-supertype mean expression, verify them, and write the download manifest.
#
# Rerunnable: files already present with the expected MD5 are not downloaded again.
# Usage: Rscript download_abc_wmb.R            (from any working directory)
#
# Expression source = MapMyCells precomputed statistics for the WMB taxonomy (per-cluster
# n_cells, sum and sum of squares of log2(CPM+1), gt0/gt1/ge1 counts, all genes). It is the
# smallest official file with aggregate expression per cluster for all genes (1.38 GB);
# no raw cell-by-gene matrix is downloaded. The h5 stores the AIT21.0 cluster means
# (cl.means.freeze, n_cells = 1 per cluster), so cluster cell counts come from the taxonomy
# tables; cell_metadata.csv gives per-cluster dissection region and 10x chemistry, and
# example_genes_all_cells_expression.csv (log2(CPM+1) of a few genes in every cell) is used
# only to check how the cluster means were computed.

suppressPackageStartupMessages({
  library(jsonlite)
  library(digest)
})

args_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
script_dir <- if (length(args_file)) dirname(normalizePath(args_file)) else getwd()
abc_root <- normalizePath(file.path(script_dir, ".."), winslash = "/")
raw_dir <- file.path(abc_root, "raw")
dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)

bucket <- "https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com"
release <- "20261031"          # ABC Atlas release manifest used to locate the files
options(timeout = 7200)

# relative S3 key, source_name, version label
files <- data.frame(
  key = c(
    sprintf("releases/%s/manifest.json", release),
    "mapmycells/WMB-taxonomy/20240831/precomputed_stats_ABC_revision_230821.h5",
    "metadata/WMB-10X/20241115/gene.csv",
    "metadata/WMB-10X/20241115/cell_metadata.csv",
    "metadata/WMB-10X/20241115/views/example_genes_all_cells_expression.csv",
    "metadata/WMB-taxonomy/20231215/cluster.csv",
    "metadata/WMB-taxonomy/20231215/cluster_annotation_term.csv",
    "metadata/WMB-taxonomy/20231215/cluster_annotation_term_set.csv",
    "metadata/WMB-taxonomy/20231215/cluster_to_cluster_annotation_membership.csv",
    "metadata/WMB-taxonomy/20231215/views/cluster_to_cluster_annotation_membership_pivoted.csv",
    "metadata/WMB-taxonomy/20231215/views/cluster_annotation_term_with_counts.csv",
    "metadata/WMB-taxonomy/20231215/cl.df_CCN202307220.xlsx"
  ),
  source_name = c(
    "ABC_release_manifest",
    "ABC_WMB_mapmycells_precomputed_stats",
    "ABC_WMB10X_gene",
    "ABC_WMB10X_cell_metadata",
    "ABC_WMB10X_example_genes_all_cells",
    "ABC_WMB_taxonomy_cluster",
    "ABC_WMB_taxonomy_cluster_annotation_term",
    "ABC_WMB_taxonomy_cluster_annotation_term_set",
    "ABC_WMB_taxonomy_cluster_membership",
    "ABC_WMB_taxonomy_cluster_membership_pivoted",
    "ABC_WMB_taxonomy_term_with_counts",
    "ABC_WMB_taxonomy_cluster_table_CCN20230722"
  ),
  version = c(release, "20240831 (ABC_revision_230821; taxonomy CCN20230722)",
              "20241115", "20241115", "20241115", rep("20231215 (taxonomy CCN20230722)", 7)),
  stringsAsFactors = FALSE
)

download_one <- function(url, dest) {
  dir.create(dirname(dest), showWarnings = FALSE, recursive = TRUE)
  tmp <- paste0(dest, ".part")
  status <- utils::download.file(url, tmp, mode = "wb", quiet = TRUE, method = "libcurl")
  if (status != 0) stop("download failed: ", url)
  invisible(file.rename(tmp, dest))
}

# Release manifest first: it carries the official size and MD5 of every other file.
man_dest <- file.path(raw_dir, files$key[1])
if (!file.exists(man_dest)) download_one(file.path(bucket, files$key[1]), man_dest)
abc_man <- fromJSON(man_dest, simplifyVector = FALSE)

# collect url -> (size, md5) from the release manifest
official <- list()
walk <- function(x) {
  if (is.list(x)) {
    if (!is.null(x$url) && !is.null(x$file_hash)) {
      official[[x$url]] <<- list(size = x$size, md5 = x$file_hash)
    }
    for (el in x) walk(el)
  }
}
walk(abc_man$file_listing)

rows <- vector("list", nrow(files))
for (i in seq_len(nrow(files))) {
  url <- file.path(bucket, files$key[i])
  dest <- file.path(raw_dir, files$key[i])
  off <- official[[url]]
  # every file is fetched by this script; on reruns it is only re-verified, and
  # downloaded_at stays the file's write time (= time of the original download)
  status <- "downloaded"
  fresh <- !file.exists(dest) ||
    (!is.null(off) && unname(tools::md5sum(dest)) != off$md5)
  if (fresh) {
    message("downloading ", files$key[i])
    download_one(url, dest)
  }
  size <- file.info(dest)$size
  md5 <- unname(tools::md5sum(dest))
  md5_ok <- if (is.null(off)) NA else identical(md5, off$md5)
  size_ok <- if (is.null(off)) NA else isTRUE(size == off$size)
  if (isFALSE(md5_ok) || isFALSE(size_ok)) stop("verification failed for ", files$key[i])
  sha <- digest(file = dest, algo = "sha256")
  rows[[i]] <- data.frame(
    source_name = files$source_name[i],
    url = url,
    file = normalizePath(dest, winslash = "/"),
    status = status,
    downloaded_at = format(file.info(dest)$mtime, "%Y-%m-%dT%H:%M:%S%z"),
    hash = sha,
    version = files$version[i],
    size_bytes = format(size, scientific = FALSE),
    md5 = md5,
    md5_matches_abc_manifest = md5_ok,
    stringsAsFactors = FALSE
  )
  message(sprintf("%-45s %12s bytes  md5_ok=%s", files$source_name[i], size, md5_ok))
}
manifest <- do.call(rbind, rows)
write.csv(manifest, file.path(abc_root, "download_manifest.csv"), row.names = FALSE)
message("manifest written: ", file.path(abc_root, "download_manifest.csv"))

# Shared helpers for the gene annotation review (tools/gene_annotation_review). Read-only with respect to the
# pipeline: these scripts read DIA-NN outputs and public database records and write only to --out.
suppressPackageStartupMessages({ library(data.table); library(jsonlite) })

`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a

# --key value arguments; flags listed in `required` must be present
review_args <- function(required = character(), defaults = list()) {
  a <- commandArgs(trailingOnly = TRUE)
  out <- defaults
  i <- 1
  while (i <= length(a)) {
    if (startsWith(a[i], "--") && i < length(a) && !startsWith(a[i + 1], "--")) { out[[sub("^--", "", a[i])]] <- a[i + 1]; i <- i + 2 }
    else { out[[sub("^--", "", a[i])]] <- TRUE; i <- i + 1 }
  }
  miss <- setdiff(required, names(out))
  if (length(miss)) stop("missing argument(s): ", paste0("--", miss, collapse = ", "), call. = FALSE)
  out
}

# isoleucine and leucine have the same mass, so peptides are compared with I read as L
il <- function(x) gsub("I", "L", x, fixed = TRUE)
split_ids <- function(x) strsplit(as.character(x), ";", fixed = TRUE)

# facility run names end in _S<sample>_S<plate>-<well>_<n>_<run>.d; sample + run number identify a run
run_key <- function(x) sub(".*_(S[0-9]{3})_S[0-9]-[A-H][0-9]{1,2}_[0-9]+_([0-9]+)[.]d$", "\\1|\\2", basename(gsub("\\\\", "/", x)))

pr_meta_cols <- c("Protein.Group", "Protein.Ids", "Protein.Names", "Genes", "First.Protein.Description",
                  "Proteotypic", "Stripped.Sequence", "Modified.Sequence", "Precursor.Charge", "Precursor.Id")
pg_meta_cols <- c("Protein.Group", "Protein.Ids", "Protein.Names", "Genes", "First.Protein.Description")

# cached download with a retrieval log (URL, file, time, status); an existing non-empty file is reused
fetch_cached <- function(url, file, log_file, pause = 0.4) {
  if (file.exists(file) && file.size(file) > 0) return(invisible(file))
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  t0 <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  st <- tryCatch({ utils::download.file(url, file, quiet = TRUE, method = "libcurl",
    headers = c(`User-Agent` = "proteomics-annotation-check")); "ok" }, error = function(e) conditionMessage(e))
  fwrite(data.table(url = url, file = basename(file), retrieved = t0, status = st), log_file, append = file.exists(log_file))
  Sys.sleep(pause)
  invisible(file)
}

# MaxLFQ: median pairwise log2 ratios over shared precursors, then least squares for the run levels; levels are
# anchored to the mean log2 precursor sum so that values are on a comparable scale. X: precursors x runs, log2.
maxlfq <- function(X) {
  n <- ncol(X); idx <- which(colSums(!is.na(X)) > 0)
  if (!length(idx)) return(rep(NA_real_, n))
  pairs <- list(); rhs <- numeric(); k <- 0
  for (a in seq_along(idx)) for (b in seq_along(idx)) if (a < b) {
    d <- X[, idx[a]] - X[, idx[b]]
    if (all(is.na(d))) next
    k <- k + 1; pairs[[k]] <- c(a, b); rhs[k] <- stats::median(d, na.rm = TRUE)
  }
  out <- rep(NA_real_, n)
  if (!k) { out[idx] <- log2(colSums(2^X[, idx, drop = FALSE], na.rm = TRUE)); return(out) }
  A <- matrix(0, k + 1, length(idx))
  for (i in seq_len(k)) { A[i, pairs[[i]][1]] <- 1; A[i, pairs[[i]][2]] <- -1 }
  A[k + 1, ] <- 1
  sol <- qr.solve(A, c(rhs, 0))
  out[idx] <- sol - mean(sol) + mean(log2(colSums(2^X[, idx, drop = FALSE], na.rm = TRUE)))
  out
}

read_fasta <- function(path) {
  x <- readLines(path, warn = FALSE); h <- grep("^>", x); e <- c(h[-1] - 1, length(x))
  data.table(header = sub("^>", "", x[h]),
             accession = sub("^(sp|tr)\\|([^|]+)\\|.*", "\\2", sub("^>", "", x[h])),
             gene = ifelse(grepl(" GN=", x[h]), sub(".* GN=([^ ]+).*", "\\1", x[h]), NA_character_),
             sequence = vapply(seq_along(h), function(i) paste(x[(h[i] + 1):e[i]], collapse = ""), ""))
}

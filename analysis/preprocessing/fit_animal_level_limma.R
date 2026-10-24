# ================================================================
# Script: analysis/preprocessing/fit_animal_level_limma.R
# Stage: core
# Scope: dataset_specific
# Consumes: required data/processed/01_preprocessing/protigy_input_animal_level/<dataset>/<dataset>_animal_level.gct
# Produces: data/processed/01_preprocessing/protigy_output_animal_level/<dataset>/stat_results_for_ssGSEA_<label>_proteome.gct; results/preprocessing/fit_animal_level_limma/<dataset>/manifests/run_manifest.csv
# Notes: Scripted replacement for the former manual ProTigy step. One limma fit per dataset (~0 + region_layer_ExpGroup),
#        within-unit contrasts 2/1, 3/2 and 3/1, eBayes(robust = TRUE, trend = FALSE), BH per contrast. Reproduces the
#        stored ProTigy outputs to ~1e-13 on the same input (checked 2026-10-19). Writes ProTigy's stat-GCT layout so
#        extract_protigy_contrasts.R and the animal_level_protigy_da_v1 contract are unchanged.
# Dataset behavior: runs for neuron_neuropil,neuron_soma,microglia according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# ================================================================

suppressPackageStartupMessages(library(limma))
paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "protigy_stat_gct_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = "") { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[[i + 1L]] }
dataset_arg <- arg_value("--dataset", Sys.getenv("PROTEOMICS_DATASET", unset = ""))
if (!nzchar(dataset_arg)) dataset_arg <- current_dataset_from_cli()
datasets <- if (identical(dataset_arg, "all")) c("neuron_neuropil", "neuron_soma", "microglia") else validate_dataset(dataset_arg)
dry_run <- "--dry-run" %in% args
output_root <- arg_value("--output-root", repo_path("data", "processed", "01_preprocessing", "protigy_output_animal_level"))

GCT_LABEL <- c(neuron_neuropil = "neuropil", neuron_soma = "soma", microglia = "microglia")
PAIRS <- list(c("2", "1"), c("3", "2"), c("3", "1"))

read_input_gct <- function(path) {
  lines <- readLines(path, warn = FALSE)
  dims <- as.integer(strsplit(lines[[2]], "\t", fixed = TRUE)[[1]][1:4])
  nrd <- dims[[3]]; ncd <- dims[[4]]
  hdr <- strsplit(lines[[3]], "\t", fixed = TRUE)[[1]]
  desc <- do.call(rbind, strsplit(lines[4:(3 + ncd)], "\t", fixed = TRUE))
  body <- do.call(rbind, strsplit(lines[(4 + ncd):(3 + ncd + dims[[1]])], "\t", fixed = TRUE))
  sample_cols <- (nrd + 2L):length(hdr)
  values <- body[, sample_cols, drop = FALSE]
  m <- matrix(as.numeric(values), nrow = nrow(values), dimnames = list(body[, 1], hdr[sample_cols]))
  if (anyNA(m)) stop("Input GCT contains missing or non-numeric values: ", path, call. = FALSE)
  coldesc <- setNames(lapply(seq_len(nrow(desc)), function(i) desc[i, sample_cols]), desc[, 1])
  list(id = body[, 1], description = body[, 2], sample_names = hdr[sample_cols], sample_strings = values, m = m, coldesc = coldesc)
}

fmt <- function(x) as.character(x)

fit_dataset <- function(dataset) {
  in_path <- repo_path("data", "processed", "01_preprocessing", "protigy_input_animal_level", dataset, paste0(dataset, "_animal_level.gct"))
  if (!file.exists(in_path)) stop("Missing animal-level input: ", in_path, call. = FALSE)
  g <- read_input_gct(in_path)
  cell <- g$coldesc[["region_layer_ExpGroup"]]
  if (is.null(cell)) stop("Input GCT has no region_layer_ExpGroup descriptor: ", in_path, call. = FALSE)
  f <- factor(cell, levels = unique(cell))
  design <- model.matrix(~0 + f); colnames(design) <- levels(f)
  units <- unique(sub("_[0-9]+$", "", cell))
  contrasts <- unlist(lapply(PAIRS, function(pr) paste0(units, "_", pr[1], "_over_", units, "_", pr[2])))
  exprs <- unlist(lapply(PAIRS, function(pr) paste0(units, "_", pr[1], "-", units, "_", pr[2])))
  stopifnot(all(unlist(strsplit(gsub("-", " ", exprs), " ")) %in% colnames(design)))
  cm <- makeContrasts(contrasts = exprs, levels = design); colnames(cm) <- contrasts
  eb <- eBayes(contrasts.fit(lmFit(g$m, design), cm), robust = TRUE, trend = FALSE)

  stat_cols <- list(); sign_cols <- list()
  for (cc in contrasts) {
    tt <- topTable(eb, coef = cc, number = Inf, sort.by = "none", adjust.method = "BH")
    stopifnot(identical(rownames(tt), g$id))
    lp <- -log10(tt$P.Value)
    stat_cols[[cc]] <- cbind(fmt(tt$logFC), fmt(tt$AveExpr), fmt(tt$t), fmt(tt$P.Value), fmt(tt$adj.P.Val), fmt(tt$B),
                             ifelse(tt$adj.P.Val < 0.05, "TRUE", "FALSE"), fmt(lp))
    colnames(stat_cols[[cc]]) <- paste0(c("logFC.", "AveExpr.", "t.", "P.Value.", "adj.P.Val.", "B.", "significant.", "Log.P.Value."), cc)
    sign_cols[[cc]] <- fmt(round(sign(tt$logFC) * lp, 4))
  }
  stats <- do.call(cbind, stat_cols)
  signs <- do.call(cbind, sign_cols); colnames(signs) <- paste0("sign.logP.", contrasts)
  header <- c("id", "Description", colnames(stats), g$sample_names, colnames(signs))
  nrhd <- 1L + ncol(stats) + length(g$sample_names)
  body <- cbind(g$id, g$description, stats, g$sample_strings, signs)
  out_lines <- c("#1.3", paste(nrow(body), length(contrasts), nrhd, 0L, sep = "\t"), paste(header, collapse = "\t"),
                 apply(body, 1, paste, collapse = "\t"))

  out_dir <- file.path(output_root, dataset)
  out_path <- file.path(out_dir, paste0("stat_results_for_ssGSEA_", GCT_LABEL[[dataset]], "_proteome.gct"))
  if (dry_run) {
    cat(sprintf("[dry-run] %s: %d proteins, %d contrasts -> %s\n", dataset, nrow(body), length(contrasts), out_path))
    return(invisible(NULL))
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  existing <- list.files(out_dir, pattern = "\\.gct$")
  others <- existing[tolower(existing) != tolower(basename(out_path))]
  if (length(others)) stop("Other GCT files in ", out_dir, " would make the extractor ambiguous: ", paste(basename(others), collapse = ", "), call. = FALSE)
  con <- file(out_path, open = "wb"); writeLines(out_lines, con, sep = "\r\n"); close(con)

  chk <- validate_corrected_protigy_gct_contract(read_protigy_stat_gct(out_path, dataset, strict_primary = TRUE))
  if (nzchar(arg_value("--output-root"))) {   # validation runs into a scratch root write no repository manifest
    cat(sprintf("%s: %d proteins, %d contrasts, contract %s -> %s\n", dataset, nrow(body), length(contrasts), chk$status, out_path))
    return(invisible(NULL))
  }
  man_dir <- repo_path("results", "preprocessing", "fit_animal_level_limma", dataset, "manifests")
  dir.create(man_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(data.frame(dataset = dataset, input_gct = in_path, input_md5 = unname(tools::md5sum(in_path)),
                       output_gct = out_path, output_md5 = unname(tools::md5sum(out_path)), n_proteins = nrow(body),
                       n_contrasts = length(contrasts), model = "limma lmFit ~0 + region_layer_ExpGroup; eBayes(robust = TRUE, trend = FALSE); BH per contrast",
                       limma_version = as.character(packageVersion("limma")), R_version = R.version.string,
                       contract_status = chk$status, run_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                       sample_centring_stamp = as.character(normalisation_provenance_safe()$stamp_sample_centring %||% NA_character_),
                       sample_centring_setting = as.character(normalisation_provenance_safe()$sample_centring_setting %||% NA_character_)),
            file.path(man_dir, "run_manifest.csv"), row.names = FALSE)
  cat(sprintf("%s: %d proteins, %d contrasts, contract %s -> %s\n", dataset, nrow(body), length(contrasts), chk$status, out_path))
}

for (ds in datasets) fit_dataset(ds)

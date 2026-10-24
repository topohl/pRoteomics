# Evidence that the scripted limma reproduces the archived ProTigy statistics (normalisation review 2026-10-21).
# Fits the model of analysis/preprocessing/fit_animal_level_limma.R (limma ~0 + region_layer_ExpGroup; within-unit
# contrasts 2/1, 3/2, 3/1; eBayes(robust = TRUE, trend = FALSE); BH per contrast) on the ARCHIVED animal-level input and
# compares logFC, t, P and adjusted P with the ARCHIVED ProTigy output, contrast by contrast. Read-only; writes
# validation/limma_vs_protigy.csv.  USAGE  Rscript tools/normalisation_review/compare_limma_protigy.R
suppressPackageStartupMessages({ library(limma); library(dplyr); library(readr) })
ARCHIVE <- if (dir.exists("Q:/data")) "Q:" else "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/_archive_proteomics_pre_renormalisation_2026-10-19"
OUT <- file.path("results", "reviewer_audit", "normalisation_choice_2026-10-21", "validation")
LABEL <- c(neuron_neuropil = "neuropil", neuron_soma = "soma", microglia = "microglia")
read_gct <- function(path) {
  lines <- readLines(path, warn = FALSE)
  dims <- as.integer(strsplit(lines[[2]], "\t", fixed = TRUE)[[1]][1:4])
  hdr <- strsplit(lines[[3]], "\t", fixed = TRUE)[[1]]
  desc <- if (dims[4] > 0) do.call(rbind, strsplit(lines[4:(3 + dims[4])], "\t", fixed = TRUE)) else NULL
  body <- do.call(rbind, strsplit(lines[(4 + dims[4]):(3 + dims[4] + dims[1])], "\t", fixed = TRUE))
  colnames(body) <- hdr
  list(body = body, dims = dims, coldesc = if (!is.null(desc)) setNames(lapply(seq_len(nrow(desc)), function(i) desc[i, ]), desc[, 1]) else NULL, hdr = hdr)
}
rows <- list()
for (ds in names(LABEL)) {
  inp <- read_gct(file.path(ARCHIVE, "data/processed/01_preprocessing/protigy_input_animal_level", ds, paste0(ds, "_animal_level.gct")))
  sc <- (inp$dims[3] + 2L):length(inp$hdr)
  m <- matrix(as.numeric(inp$body[, sc]), nrow = nrow(inp$body), dimnames = list(inp$body[, 1], inp$hdr[sc]))
  cell <- inp$coldesc[["region_layer_ExpGroup"]][sc]
  f <- factor(cell, levels = unique(cell)); design <- model.matrix(~0 + f); colnames(design) <- levels(f)
  units <- unique(sub("_[0-9]+$", "", cell)); pairs <- list(c("2", "1"), c("3", "2"), c("3", "1"))
  nm <- unlist(lapply(pairs, function(p) paste0(units, "_", p[1], "_over_", units, "_", p[2])))
  ex <- unlist(lapply(pairs, function(p) paste0(units, "_", p[1], "-", units, "_", p[2])))
  cm <- makeContrasts(contrasts = ex, levels = design); colnames(cm) <- nm
  eb <- eBayes(contrasts.fit(lmFit(m, design), cm), robust = TRUE, trend = FALSE)
  pt <- read_gct(Sys.glob(file.path(ARCHIVE, "data/processed/01_preprocessing/protigy_output_animal_level", ds, "*.gct"))[1])
  for (cc in nm) {
    tt <- topTable(eb, coef = cc, number = Inf, sort.by = "none", adjust.method = "BH")
    k <- match(rownames(tt), pt$body[, "id"])
    get <- function(stat) as.numeric(pt$body[k, paste0(stat, ".", cc)])
    rows[[length(rows) + 1L]] <- tibble(dataset = ds, contrast = cc, n_proteins = sum(!is.na(k)),
      max_abs_diff_logFC = max(abs(tt$logFC - get("logFC")), na.rm = TRUE), max_abs_diff_t = max(abs(tt$t - get("t")), na.rm = TRUE),
      max_abs_diff_P = max(abs(tt$P.Value - get("P.Value")), na.rm = TRUE), max_abs_diff_adjP = max(abs(tt$adj.P.Val - get("adj.P.Val")), na.rm = TRUE),
      same_FDR05_calls = identical(tt$adj.P.Val < 0.05, get("adj.P.Val") < 0.05))
  }
}
res <- bind_rows(rows)
write_csv(res, file.path(OUT, "limma_vs_protigy.csv"))
print(res |> group_by(dataset) |> summarise(contrasts = n(), max_logFC = max(max_abs_diff_logFC), max_t = max(max_abs_diff_t),
                                            max_adjP = max(max_abs_diff_adjP), all_calls_equal = all(same_FDR05_calls)))

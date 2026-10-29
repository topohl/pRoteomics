# Shared by the stress-effect analyses added on 2026-10-29 (test_stress_exposure, test_unit_specificity).
# Both start from the inputs of the Figure 3 per-unit model (analysis/preprocessing/fit_animal_level_limma.R):
# the animal-level matrices (hemispheres averaged, imputed) and the canonical protein-group identities of the
# mapped per-unit tables.

STRESS_DATASETS <- c("neuron_neuropil", "neuron_soma", "microglia")
STRESS_GROUPS <- c("1" = "CON", "2" = "RES", "3" = "SUS")

animal_level_gct_path <- function(dataset) {
  repo_path("data", "processed", "01_preprocessing", "protigy_input_animal_level", dataset,
            paste0(dataset, "_animal_level.gct"))
}

# The reader of fit_animal_level_limma.R, returning the matrix with one column per animal x unit.
read_animal_level_gct <- function(dataset) {
  path <- animal_level_gct_path(dataset)
  lines <- readLines(path, warn = FALSE)
  dims <- as.integer(strsplit(lines[[2]], "\t", fixed = TRUE)[[1]][1:4])
  nrd <- dims[[3]]; ncd <- dims[[4]]
  hdr <- strsplit(lines[[3]], "\t", fixed = TRUE)[[1]]
  desc <- do.call(rbind, strsplit(lines[4:(3 + ncd)], "\t", fixed = TRUE))
  body <- do.call(rbind, strsplit(lines[(4 + ncd):(3 + ncd + dims[[1]])], "\t", fixed = TRUE))
  sample_cols <- (nrd + 2L):length(hdr)
  m <- matrix(as.numeric(body[, sample_cols, drop = FALSE]), nrow = nrow(body), dimnames = list(body[, 1], hdr[sample_cols]))
  if (anyNA(m)) stop("Animal-level GCT contains missing or non-numeric values: ", path, call. = FALSE)
  cd <- setNames(lapply(seq_len(nrow(desc)), function(i) desc[i, sample_cols]), desc[, 1])
  for (k in c("AnimalID", "ExpGroup", "region", "layer")) if (is.null(cd[[k]])) stop("GCT lacks ", k, ": ", path, call. = FALSE)
  group <- unname(STRESS_GROUPS[as.character(cd$ExpGroup)])
  if (anyNA(group)) stop("Unknown ExpGroup codes in ", path, call. = FALSE)
  list(m = m, animal = cd$AnimalID, group = group, unit = paste0(cd$region, cd$layer),
       cell = paste0(cd$region, cd$layer, "_", group), path = path)
}

# The canonical identity columns of a dataset's proteins, in GCT row order (identical in every per-unit table).
mapped_unit_table_path <- function(dataset, unit) {
  repo_path("results", "preprocessing", "map_protein_identifiers", dataset, "tables", "mapped", "forward", "per_file",
            paste0(unit, "sus_", unit, "res.csv"))
}
canonical_protein_identity <- function(dataset, unit, n_rows) {
  t <- utils::read.csv(mapped_unit_table_path(dataset, unit), stringsAsFactors = FALSE, check.names = FALSE)
  t <- t[order(t$source_row_id), , drop = FALSE]
  if (nrow(t) != n_rows || !identical(t$source_row_id, seq_len(n_rows)))
    stop("Mapped table does not follow the GCT rows: ", mapped_unit_table_path(dataset, unit), call. = FALSE)
  t
}

# Gene-level t of the Figure 3 model for every unit and contrast under any labelling of the nine animals
# (compare_external_stress_studies). One protein group per official symbol (claim-eligible, highest mean),
# fixed for every labelling.
STRESS_CONTRASTS <- c("SUS_vs_CON", "RES_vs_CON", "SUS_vs_RES", "stressed_vs_CON")
stress_gene_frame <- function(dataset) {
  g <- read_animal_level_gct(dataset)
  units <- unique(g$unit)
  ident <- canonical_protein_identity(dataset, units[[1]], nrow(g$m))
  sym <- ident$official_gene_symbol
  keep <- which(!is.na(sym) & nzchar(sym) & ident$gene_level_claim_allowed %in% c(TRUE, "TRUE"))
  keep <- keep[order(sym[keep], -rowMeans(g$m)[keep])]
  keep <- keep[!duplicated(sym[keep])]
  animals <- sort(unique(g$animal))
  list(dataset = dataset, m = g$m, animal = g$animal, unit = g$unit, units = units, rows = keep, genes = sym[keep],
       animals = animals, true_group = stats::setNames(g$group[match(animals, g$animal)], animals))
}
stress_contrast_t <- function(d, labels) {
  grp <- unname(labels[d$animal])
  cell <- factor(paste0(d$unit, "_", grp))
  design <- stats::model.matrix(~ 0 + cell); colnames(design) <- levels(cell)
  cm <- do.call(cbind, lapply(d$units, function(u) {
    w <- matrix(0, ncol(design), 4L, dimnames = list(colnames(design), paste0(u, "|", STRESS_CONTRASTS)))
    s <- paste0(u, "_SUS"); r <- paste0(u, "_RES"); c0 <- paste0(u, "_CON")
    w[c(s, c0), 1] <- c(1, -1); w[c(r, c0), 2] <- c(1, -1); w[c(s, r), 3] <- c(1, -1); w[c(s, r, c0), 4] <- c(0.5, 0.5, -1)
    w }))
  eb <- limma::eBayes(limma::contrasts.fit(limma::lmFit(d$m, design), cm), robust = TRUE, trend = FALSE)
  t <- eb$t[d$rows, , drop = FALSE]
  rownames(t) <- d$genes
  t
}
# All 1,680 assignments of nine animals to three groups of three (rows), in a fixed order.
stress_all_labellings <- function(animals) {
  stopifnot(length(animals) == 9L)
  out <- list()
  s1 <- utils::combn(9, 3)
  for (i in seq_len(ncol(s1))) {
    s2 <- utils::combn(setdiff(1:9, s1[, i]), 3)
    for (j in seq_len(ncol(s2))) {
      lab <- rep("CON", 9); lab[s1[, i]] <- "SUS"; lab[s2[, j]] <- "RES"
      out[[length(out) + 1L]] <- stats::setNames(lab, animals)
    }
  }
  do.call(rbind, out)
}

# Exact label permutation over animals for a difference of group means (two-sided on |difference|).
exact_label_test <- function(values, in_group) {
  ok <- is.finite(values)
  v <- values[ok]; g <- in_group[ok]
  n1 <- sum(g); n <- length(v)
  if (n1 < 2L || n - n1 < 2L) return(list(estimate = NA_real_, p = NA_real_, n_labellings = 0L, n1 = n1, n0 = n - n1))
  obs <- mean(v[g]) - mean(v[!g])
  idx <- utils::combn(n, n1)
  d <- apply(idx, 2, function(j) mean(v[j]) - mean(v[-j]))
  list(estimate = obs, p = mean(abs(d) >= abs(obs) - 1e-12), n_labellings = ncol(idx), n1 = n1, n0 = n - n1)
}

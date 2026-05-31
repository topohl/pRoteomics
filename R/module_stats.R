# Reusable module/overlap statistics.

jaccard_index <- function(a, b) {
  a <- unique(stats::na.omit(as.character(a)))
  b <- unique(stats::na.omit(as.character(b)))
  denom <- length(union(a, b))
  if (!denom) return(NA_real_)
  length(intersect(a, b)) / denom
}

fisher_overlap_p <- function(module_set, hit_set, universe = union(module_set, hit_set)) {
  module_set <- unique(stats::na.omit(as.character(module_set)))
  hit_set <- unique(stats::na.omit(as.character(hit_set)))
  universe <- unique(stats::na.omit(as.character(universe)))
  if (!length(module_set) || !length(hit_set) || !length(universe)) return(NA_real_)
  a <- length(intersect(module_set, hit_set))
  b <- length(setdiff(module_set, hit_set))
  c <- length(setdiff(hit_set, module_set))
  d <- max(length(universe) - a - b - c, 0)
  stats::fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value
}

bh_fdr <- function(p) stats::p.adjust(suppressWarnings(as.numeric(p)), method = "BH")

conservative_strength <- function(fdr = NA_real_, effect = NA_real_, n = NA_integer_, support = NA_real_) {
  if (exists("interpretation_strength", mode = "function")) {
    return(vapply(seq_along(fdr), function(i) {
      interpretation_strength(
        fdr = fdr[[i]],
        effect_size = if (length(effect) >= i) effect[[i]] else NA_real_,
        n = if (length(n) >= i) n[[i]] else NA_real_,
        bootstrap_support = if (length(support) >= i) support[[i]] else NA_real_
      )
    }, character(1)))
  }
  rep("exploratory", length(fdr))
}

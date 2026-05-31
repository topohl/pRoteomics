# Reusable spatial-network helpers.

make_edge_id <- function(source, target, sep = "--") {
  source <- as.character(source)
  target <- as.character(target)
  paste(pmin(source, target), pmax(source, target), sep = sep)
}

network_interpretation_strength <- function(fdr = NA_real_, stability_score = NA_real_, n_animals = NA_integer_) {
  fdr <- suppressWarnings(as.numeric(fdr))
  stability_score <- suppressWarnings(as.numeric(stability_score))
  n_animals <- suppressWarnings(as.numeric(n_animals))
  if (!is.na(fdr) && fdr <= 0.05 && !is.na(stability_score) && stability_score >= 0.7 && (is.na(n_animals) || n_animals >= 6)) {
    return("strong")
  }
  if ((!is.na(fdr) && fdr <= 0.10) || (!is.na(stability_score) && stability_score >= 0.6)) return("moderate")
  "exploratory"
}

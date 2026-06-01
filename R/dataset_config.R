# Shared biological dataset-family selection helpers.

valid_datasets <- function() {
  c("neuron_neuropil", "neuron_soma", "microglia")
}

dataset_contracts <- function() {
  list(
    neuron_neuropil = list(
      label = "Neuron neuropil",
      region = TRUE,
      layer = TRUE,
      celltype_roi = FALSE,
      purified_celltype = FALSE,
      interpretation = "Region/layer-resolved neuron neuropil proteomics."
    ),
    neuron_soma = list(
      label = "Neuron soma",
      region = TRUE,
      layer = TRUE,
      celltype_roi = FALSE,
      purified_celltype = FALSE,
      interpretation = "Region/layer-resolved neuronal soma-enriched proteomics."
    ),
    microglia = list(
      label = "Microglia-enriched ROI",
      region = TRUE,
      layer = FALSE,
      celltype_roi = TRUE,
      purified_celltype = FALSE,
      interpretation = "Region-resolved microglia-enriched ROI/local microenvironment proteomics; not purified microglia."
    )
  )
}

normalize_dataset <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[[:space:]-]+", "_", x)
  out <- x
  out[x %in% c("neuropil", "neuron_neuropil")] <- "neuron_neuropil"
  out[x %in% c("soma", "neuron_soma")] <- "neuron_soma"
  out[x %in% c("microglia", "microglial")] <- "microglia"
  out
}

validate_dataset <- function(dataset, source = "dataset") {
  dataset <- normalize_dataset(dataset)
  if (length(dataset) != 1L || is.na(dataset) || !nzchar(dataset)) {
    stop("Missing ", source, ". Expected one of: ", paste(valid_datasets(), collapse = ", "), call. = FALSE)
  }
  if (!dataset %in% valid_datasets()) {
    stop(
      "Unsupported ", source, ": '", dataset, "'. Expected one of: ",
      paste(valid_datasets(), collapse = ", "), ". Set PROTEOMICS_DATASET to a valid dataset family.",
      call. = FALSE
    )
  }
  dataset
}

current_dataset <- function(default = "neuron_neuropil") {
  candidates <- c(
    PROTEOMICS_DATASET = Sys.getenv("PROTEOMICS_DATASET", unset = ""),
    PROTEOMICS_COMPARISON = Sys.getenv("PROTEOMICS_COMPARISON", unset = ""),
    PROTEOMICS_GCT_COMPARISON = Sys.getenv("PROTEOMICS_GCT_COMPARISON", unset = "")
  )
  hit <- candidates[nzchar(candidates)][1]
  if (length(hit) && !is.na(hit)) {
    return(validate_dataset(unname(hit), source = names(hit)))
  }
  validate_dataset(default, source = "default dataset")
}

dataset_capabilities <- function(dataset = current_dataset()) {
  dataset <- validate_dataset(dataset)
  dataset_contracts()[[dataset]]
}

dataset_has_capability <- function(dataset = current_dataset(), capability) {
  caps <- dataset_capabilities(dataset)
  if (!capability %in% names(caps)) {
    stop("Unknown dataset capability: ", capability, call. = FALSE)
  }
  isTRUE(caps[[capability]])
}

assert_dataset_capability <- function(dataset = current_dataset(), capability, analysis = "analysis") {
  if (!dataset_has_capability(dataset, capability)) {
    stop(
      "Dataset '", dataset, "' does not support ", capability, "-level ", analysis,
      ". Dataset interpretation: ", dataset_capabilities(dataset)$interpretation,
      call. = FALSE
    )
  }
  invisible(TRUE)
}

dataset_interpretation <- function(dataset = current_dataset()) {
  dataset_capabilities(dataset)$interpretation
}

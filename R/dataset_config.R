# Shared biological dataset-family selection helpers.

valid_datasets <- function() {
  c("neuron_neuropil", "neuron_soma", "microglia")
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

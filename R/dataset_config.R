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
      region_layer = TRUE,
      spatial_unit = "region_layer",
      celltype_roi = FALSE,
      purified_celltype = FALSE,
      interpretation = "Region/layer-resolved neuron neuropil proteomics."
    ),
    neuron_soma = list(
      label = "Neuron soma",
      region = TRUE,
      layer = FALSE,
      region_layer = FALSE,
      spatial_unit = "region",
      celltype_roi = FALSE,
      purified_celltype = FALSE,
      interpretation = "Region-resolved neuronal soma-enriched proteomics."
    ),
    microglia = list(
      label = "Microglia-enriched ROI",
      region = TRUE,
      layer = FALSE,
      region_layer = FALSE,
      spatial_unit = "region",
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

cli_arg_value <- function(flag, default = "", args = commandArgs(trailingOnly = TRUE)) {
  hit <- which(args == flag)
  if (!length(hit) || hit[[1]] == length(args)) return(default)
  args[[hit[[1]] + 1L]]
}

current_dataset_from_cli <- function(default = "neuron_neuropil",
                                     flag = "--dataset",
                                     allow_all = FALSE,
                                     args = commandArgs(trailingOnly = TRUE)) {
  dataset_cli <- cli_arg_value(flag, default = "", args = args)
  if (nzchar(dataset_cli)) {
    dataset_cli <- normalize_dataset(dataset_cli)
    if (isTRUE(allow_all) && identical(dataset_cli, "all")) {
      Sys.setenv(PROTEOMICS_DATASET = dataset_cli)
      return(dataset_cli)
    }
    dataset_cli <- validate_dataset(dataset_cli, source = flag)
    Sys.setenv(PROTEOMICS_DATASET = dataset_cli)
    return(dataset_cli)
  }
  current_dataset(default = default)
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

dataset_spatial_unit <- function(dataset = current_dataset()) {
  dataset_capabilities(dataset)$spatial_unit
}

dataset_supports_region_layer <- function(dataset = current_dataset()) {
  dataset_has_capability(dataset, "region_layer")
}

# Narrow Figure 2 control-anatomy helpers.  This file intentionally has no
# general pipeline abstraction: it protects the fixed analysis requested here.

control_spatial_hemisphere <- function(sample_id) {
  x <- as.character(sample_id)
  out <- sub(".*_([LR])_[^_]+_[^_]+_Neuron_.*", "\\1", x)
  out[!out %in% c("L", "R")] <- NA_character_
  if (anyNA(out)) stop("Could not resolve hemisphere from canonical sample identifier.", call. = FALSE)
  factor(out, levels = c("L", "R"))
}

control_spatial_prepare_metadata <- function(metadata) {
  needed <- c("Sample", "StressGroup", "Region", "Layer", "SpatialUnit", "AnimalID")
  missing <- setdiff(needed, names(metadata))
  if (length(missing)) stop("Canonical metadata is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  x <- metadata[metadata$StressGroup == "CON" & (is.na(metadata$Exclude) | !metadata$Exclude), , drop = FALSE]
  if (!nrow(x)) stop("No CON samples in canonical metadata.", call. = FALSE)
  if (anyNA(x$AnimalID) || any(!nzchar(as.character(x$AnimalID)))) stop("CON samples have missing AnimalID.", call. = FALSE)
  x$hemisphere <- control_spatial_hemisphere(x$Sample)
  x$anatomical_unit <- ifelse(x$SpatialUnit == "region_layer", paste(toupper(x$Region), toupper(x$Layer), sep = "_"), toupper(x$Region))
  x
}

control_spatial_target_rest_weights <- function(levels, target) {
  levels <- as.character(levels); target <- as.character(target)
  if (!target %in% levels || length(levels) < 2L) stop("Target-versus-rest contrast is not estimable.", call. = FALSE)
  w <- stats::setNames(rep(-1 / (length(levels) - 1L), length(levels)), levels); w[[target]] <- 1
  if (!isTRUE(all.equal(sum(w), 0, tolerance = 1e-12))) stop("Contrast weights do not sum to zero.", call. = FALSE)
  w
}

control_spatial_region_mean_weights <- function(units, regions, target_region = "DG") {
  units <- as.character(units); regions <- as.character(regions)
  if (length(units) != length(regions) || !target_region %in% regions) stop("Requested neuropil contrast is not estimable.", call. = FALSE)
  region_units <- split(units, regions)
  other <- setdiff(names(region_units), target_region)
  if (!length(other)) stop("Requested neuropil contrast has no non-DG regions.", call. = FALSE)
  w <- stats::setNames(rep(0, length(units)), units)
  w[region_units[[target_region]]] <- 1 / length(region_units[[target_region]])
  for (r in other) w[region_units[[r]]] <- -1 / (length(other) * length(region_units[[r]]))
  if (!isTRUE(all.equal(sum(w), 0, tolerance = 1e-12))) stop("Contrast weights do not sum to zero.", call. = FALSE)
  w
}

control_spatial_design <- function(metadata) {
  unit <- factor(metadata$anatomical_unit)
  base <- stats::model.matrix(~ 0 + unit)
  use_hemi <- nlevels(droplevels(metadata$hemisphere)) == 2L
  design <- if (use_hemi) cbind(base, hemisphere_R = as.numeric(metadata$hemisphere == "R")) else base
  if (qr(design)$rank < ncol(design)) {
    if (!use_hemi) stop("Rank-deficient anatomical-unit design.", call. = FALSE)
    design <- base; use_hemi <- FALSE; reason <- "hemisphere perfectly aliased with anatomical unit or produced rank deficiency"
  } else reason <- if (use_hemi) NA_character_ else "only one hemisphere level present"
  colnames(design) <- sub("^unit", "anatomical_unit_", colnames(design))
  list(design = design, hemisphere_included = use_hemi, hemisphere_omission_reason = reason)
}

control_spatial_read_kaulich_sheet <- function(workbook, sheet) {
  x <- as.data.frame(readxl::read_excel(workbook, sheet = sheet, skip = 4, .name_repair = "minimal"), check.names = FALSE)
  names(x) <- make.unique(trimws(names(x)))
  x
}

control_spatial_gene_match_key <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x)] <- ""
  toupper(x)
}

control_spatial_parse_gene_cells <- function(x, source_row_id = seq_along(x)) {
  if (length(x) != length(source_row_id)) stop("Gene cells and source-row identifiers differ in length.", call. = FALSE)
  parsed <- lapply(seq_along(x), function(i) {
    raw <- as.character(x[[i]])
    genes <- if (is.na(raw)) character() else trimws(strsplit(raw, ";", fixed = TRUE)[[1]])
    genes <- genes[nzchar(genes)]
    keys <- control_spatial_gene_match_key(genes)
    keep <- nzchar(keys) & !duplicated(keys)
    genes <- genes[keep]; keys <- keys[keep]
    if (!length(genes)) {
      return(data.frame(
        source_row_id = source_row_id[[i]],
        source_gene_symbol_raw = raw,
        parsed_source_gene_symbol = NA_character_,
        gene_match_key = "",
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      source_row_id = rep(source_row_id[[i]], length(genes)),
      source_gene_symbol_raw = rep(raw, length(genes)),
      parsed_source_gene_symbol = genes,
      gene_match_key = keys,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, parsed)
}

control_spatial_key_audit <- function(symbols, keys = control_spatial_gene_match_key(symbols),
                                      count_name = "n_distinct_official_symbols",
                                      symbols_name = "official_symbols") {
  if (length(symbols) != length(keys)) stop("Symbols and matching keys differ in length.", call. = FALSE)
  x <- unique(data.frame(symbol = as.character(symbols), gene_match_key = as.character(keys), stringsAsFactors = FALSE))
  x <- x[nzchar(x$gene_match_key) & !is.na(x$symbol) & nzchar(trimws(x$symbol)), , drop = FALSE]
  if (!nrow(x)) {
    out <- data.frame(gene_match_key = character(), n = integer(), symbols = character(), stringsAsFactors = FALSE)
  } else {
    by_key <- split(x$symbol, x$gene_match_key)
    out <- data.frame(
      gene_match_key = names(by_key),
      n = vapply(by_key, function(z) length(unique(z)), integer(1)),
      symbols = vapply(by_key, function(z) paste(sort(unique(z)), collapse = ";"), character(1)),
      stringsAsFactors = FALSE
    )
  }
  names(out)[names(out) == "n"] <- count_name
  names(out)[names(out) == "symbols"] <- symbols_name
  out
}

control_spatial_kaulich_signatures <- function(workbook, sheet, workbook_sha256) {
  x <- control_spatial_read_kaulich_sheet(workbook, sheet)
  gene_col <- names(x)[tolower(names(x)) == "genes"][1]
  protein_col <- names(x)[tolower(names(x)) %in% c("protein.id", "protein_id")][1]
  if (is.na(gene_col) || is.na(protein_col)) stop("Kaulich sheet lacks Protein.ID or Genes: ", sheet, call. = FALSE)
  x$.source_row_id <- seq_len(nrow(x))
  sig_cols <- grep("^sig_", names(x), value = TRUE); out <- list()
  for (sig in sig_cols) {
    suffix <- sub("^sig_", "", sig); effect <- paste0("log2fc_", suffix); q <- paste0("qvalue_", suffix)
    keep <- as.logical(x[[sig]]) %in% TRUE & is.finite(suppressWarnings(as.numeric(x[[effect]]))) & suppressWarnings(as.numeric(x[[effect]])) > 0 & suppressWarnings(as.numeric(x[[q]])) < 0.05
    selected <- x[keep, , drop = FALSE]
    genes <- control_spatial_parse_gene_cells(selected[[gene_col]], selected$.source_row_id)
    source_rows <- match(genes$source_row_id, selected$.source_row_id)
    out[[length(out)+1L]] <- data.frame(
      source_sheet = sheet,
      external_signature = suffix,
      source_row_id = genes$source_row_id,
      source_protein_identifier = as.character(selected[[protein_col]][source_rows]),
      source_gene_symbol_raw = genes$source_gene_symbol_raw,
      parsed_source_gene_symbol = genes$parsed_source_gene_symbol,
      gene_match_key = genes$gene_match_key,
      source_effect = as.numeric(selected[[effect]][source_rows]),
      source_q_value = as.numeric(selected[[q]][source_rows]),
      workbook_sha256 = workbook_sha256,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

control_spatial_region_signature_membership <- function(
    ca1_ca3_effect, ca1_ca3_q, ca1_ca3_significant,
    ca1_dg_effect, ca1_dg_q, ca1_dg_significant,
    ca3_dg_effect, ca3_dg_q, ca3_dg_significant) {
  directional <- function(effect, q_value, significant, direction) {
    effect <- suppressWarnings(as.numeric(effect))
    q_value <- suppressWarnings(as.numeric(q_value))
    as.logical(significant) %in% TRUE &
      is.finite(effect) & is.finite(q_value) & q_value < 0.05 &
      direction * effect > 0
  }
  ca1_ca3_positive <- directional(ca1_ca3_effect, ca1_ca3_q, ca1_ca3_significant, 1)
  ca1_ca3_negative <- directional(ca1_ca3_effect, ca1_ca3_q, ca1_ca3_significant, -1)
  ca1_dg_positive <- directional(ca1_dg_effect, ca1_dg_q, ca1_dg_significant, 1)
  ca1_dg_negative <- directional(ca1_dg_effect, ca1_dg_q, ca1_dg_significant, -1)
  ca3_dg_positive <- directional(ca3_dg_effect, ca3_dg_q, ca3_dg_significant, 1)
  ca3_dg_negative <- directional(ca3_dg_effect, ca3_dg_q, ca3_dg_significant, -1)
  data.frame(
    CA1 = ca1_ca3_positive & ca1_dg_positive,
    `CA2/3` = ca1_ca3_negative & ca3_dg_positive,
    DG = ca1_dg_negative & ca3_dg_negative,
    check.names = FALSE
  )
}

control_spatial_kaulich_region_signatures <- function(workbook, sheet, workbook_sha256) {
  x <- control_spatial_read_kaulich_sheet(workbook, sheet)
  gene_col <- names(x)[tolower(names(x)) == "genes"][1]
  protein_col <- names(x)[tolower(names(x)) %in% c("protein.id", "protein_id")][1]
  if (is.na(gene_col) || is.na(protein_col)) stop("Kaulich sheet lacks Protein.ID or Genes: ", sheet, call. = FALSE)
  suffixes <- if (all(c("log2fc_CA1vCA3", "log2fc_CA1vDG", "log2fc_CA3vDG") %in% names(x))) {
    c(ca1_ca3 = "CA1vCA3", ca1_dg = "CA1vDG", ca3_dg = "CA3vDG")
  } else if (all(c("log2fc_CA1vsCA3", "log2fc_CA1vsDG", "log2fc_CA3vsDG") %in% names(x))) {
    c(ca1_ca3 = "CA1vsCA3", ca1_dg = "CA1vsDG", ca3_dg = "CA3vsDG")
  } else {
    stop("Kaulich subregion sheet has unsupported comparison orientations: ", sheet, call. = FALSE)
  }
  required <- unlist(lapply(suffixes, function(z) paste0(c("log2fc_", "qvalue_", "sig_"), z)))
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Kaulich subregion sheet is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  x$.source_row_id <- seq_len(nrow(x))
  membership <- control_spatial_region_signature_membership(
    x[[paste0("log2fc_", suffixes[["ca1_ca3"]])]],
    x[[paste0("qvalue_", suffixes[["ca1_ca3"]])]],
    x[[paste0("sig_", suffixes[["ca1_ca3"]])]],
    x[[paste0("log2fc_", suffixes[["ca1_dg"]])]],
    x[[paste0("qvalue_", suffixes[["ca1_dg"]])]],
    x[[paste0("sig_", suffixes[["ca1_dg"]])]],
    x[[paste0("log2fc_", suffixes[["ca3_dg"]])]],
    x[[paste0("qvalue_", suffixes[["ca3_dg"]])]],
    x[[paste0("sig_", suffixes[["ca3_dg"]])]]
  )
  out <- lapply(names(membership), function(signature_name) {
    selected <- x[membership[[signature_name]], , drop = FALSE]
    genes <- control_spatial_parse_gene_cells(selected[[gene_col]], selected$.source_row_id)
    source_rows <- match(genes$source_row_id, selected$.source_row_id)
    data.frame(
      source_sheet = sheet,
      external_signature = signature_name,
      source_row_id = genes$source_row_id,
      source_protein_identifier = as.character(selected[[protein_col]][source_rows]),
      source_gene_symbol_raw = genes$source_gene_symbol_raw,
      parsed_source_gene_symbol = genes$parsed_source_gene_symbol,
      gene_match_key = genes$gene_match_key,
      source_effect_CA1_vs_CA3 = suppressWarnings(as.numeric(selected[[paste0("log2fc_", suffixes[["ca1_ca3"]])]][source_rows])),
      source_q_value_CA1_vs_CA3 = suppressWarnings(as.numeric(selected[[paste0("qvalue_", suffixes[["ca1_ca3"]])]][source_rows])),
      source_significant_CA1_vs_CA3 = as.logical(selected[[paste0("sig_", suffixes[["ca1_ca3"]])]][source_rows]),
      source_effect_CA1_vs_DG = suppressWarnings(as.numeric(selected[[paste0("log2fc_", suffixes[["ca1_dg"]])]][source_rows])),
      source_q_value_CA1_vs_DG = suppressWarnings(as.numeric(selected[[paste0("qvalue_", suffixes[["ca1_dg"]])]][source_rows])),
      source_significant_CA1_vs_DG = as.logical(selected[[paste0("sig_", suffixes[["ca1_dg"]])]][source_rows]),
      source_effect_CA3_vs_DG = suppressWarnings(as.numeric(selected[[paste0("log2fc_", suffixes[["ca3_dg"]])]][source_rows])),
      source_q_value_CA3_vs_DG = suppressWarnings(as.numeric(selected[[paste0("qvalue_", suffixes[["ca3_dg"]])]][source_rows])),
      source_significant_CA3_vs_DG = as.logical(selected[[paste0("sig_", suffixes[["ca3_dg"]])]][source_rows]),
      workbook_sha256 = workbook_sha256,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

control_spatial_signature_interpretation <- function(dataset, internal_contrast, external_signature) {
  if (dataset == "neuron_soma") {
    domain <- "Soma region — tissue reference"
    compartment <- "tissue"
  } else if (grepl("_vs_mean_other_CA1_strata$", internal_contrast)) {
    domain <- "CA1 strata — synaptosome reference"
    compartment <- "synaptosome"
  } else {
    domain <- "Neuropil region — synaptosome reference"
    compartment <- "synaptosome"
  }
  expected <- FALSE
  match_type <- "specificity_comparison"
  note <- "Alternative anatomical signature used as a specificity comparison."
  if (dataset == "neuron_soma" &&
      ((grepl("^CA1_", internal_contrast) && external_signature == "CA1") ||
       (grepl("^DG_", internal_contrast) && external_signature == "DG"))) {
    expected <- TRUE; match_type <- "exact"
    note <- "Expected exact regional correspondence."
  } else if (dataset == "neuron_soma" &&
             grepl("^CA[23]_", internal_contrast) && external_signature == "CA2/3") {
    expected <- TRUE; match_type <- "approximate"
    note <- "External reference combines CA2/3; this is not independent CA3-specific validation."
  } else if (internal_contrast == "DG_neuropil_vs_mean_non_DG_regions" &&
             external_signature == "DG") {
    expected <- TRUE; match_type <- "exact_directional_match"
    note <- "Expected positive DG direction against non-DG neuropil."
  } else if (internal_contrast == "CA1_SO_vs_CA3_SO" && external_signature == "CA1") {
    expected <- TRUE; match_type <- "exact_CA1_direction"
    note <- "Expected positive CA1 direction."
  } else if (internal_contrast == "CA1_SO_vs_CA3_SO" && external_signature == "CA2/3") {
    expected <- TRUE; match_type <- "approximate_opposing_direction"
    note <- "Expected opposing direction; external reference combines CA2/3."
  } else if (grepl("_vs_mean_other_CA1_strata$", internal_contrast)) {
    target <- sub("^CA1_([^_]+)_.*$", "\\1", internal_contrast)
    if (identical(target, external_signature)) {
      expected <- TRUE; match_type <- "exact"
      note <- "Expected exact target-stratum correspondence."
    }
  }
  data.frame(
    validation_domain = domain,
    external_source_compartment = compartment,
    expected_match = expected,
    match_type = match_type,
    interpretation_note = note,
    stringsAsFactors = FALSE
  )
}

control_spatial_signature_family <- function(validation_domain) {
  unname(c(
    "Soma region — tissue reference" = "soma_tissue",
    "Neuropil region — synaptosome reference" = "neuropil_subregion",
    "CA1 strata — synaptosome reference" = "ca1_strata"
  )[as.character(validation_domain)])
}

control_spatial_signature_family_sizes <- function() {
  c(soma_tissue = 12L, neuropil_subregion = 6L, ca1_strata = 12L)
}

control_spatial_signature_max_gs_size <- function(largest_mapped_signature_size) {
  largest_mapped_signature_size <- as.integer(largest_mapped_signature_size)
  if (length(largest_mapped_signature_size) != 1L || is.na(largest_mapped_signature_size) ||
      largest_mapped_signature_size < 0L) {
    stop("Largest mapped signature size must be one non-negative integer.", call. = FALSE)
  }
  max(500L, largest_mapped_signature_size)
}

control_spatial_apply_signature_fdr <- function(x, family_sizes = control_spatial_signature_family_sizes()) {
  required <- c("signature_fdr_family", "status", "p_value")
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Signature GSEA results are missing: ", paste(missing, collapse = ", "), call. = FALSE)
  unknown <- setdiff(unique(x$signature_fdr_family), names(family_sizes))
  if (length(unknown)) stop("Unknown signature FDR families: ", paste(unknown, collapse = ", "), call. = FALSE)
  observed <- table(factor(x$signature_fdr_family, levels = names(family_sizes)))
  if (any(observed != family_sizes)) {
    stop(
      "Prespecified signature family sizes differ from observed tests: ",
      paste(names(family_sizes), observed, family_sizes, sep = "=", collapse = "; "),
      call. = FALSE
    )
  }
  x$signature_family_size <- unname(family_sizes[x$signature_fdr_family])
  x$signature_FDR <- NA_real_
  if ("p_adjust" %in% names(x) && !"single_set_p_adjust" %in% names(x)) {
    x$single_set_p_adjust <- x$p_adjust
  }
  for (family in names(family_sizes)) {
    idx <- which(
      x$signature_fdr_family == family &
        x$status == "completed" &
        is.finite(suppressWarnings(as.numeric(x$p_value)))
    )
    if (length(idx)) {
      x$signature_FDR[idx] <- stats::p.adjust(
        as.numeric(x$p_value[idx]), method = "BH", n = family_sizes[[family]]
      )
    }
  }
  x
}

control_spatial_contrast_label <- function(x) {
  labels <- c(
    CA1_vs_mean_other_soma_regions = "CA1 soma vs other regions",
    CA2_vs_mean_other_soma_regions = "CA2 soma vs other regions",
    CA3_vs_mean_other_soma_regions = "CA3 soma vs other regions",
    DG_vs_mean_other_soma_regions = "DG soma vs other regions",
    DG_neuropil_vs_mean_non_DG_regions = "DG neuropil vs non-DG neuropil",
    DG_MO_vs_mean_other_DG_layers = "DG molecular layer vs other DG layer",
    DG_PO_vs_mean_other_DG_layers = "DG polymorphic layer vs other DG layer",
    CA1_SO_vs_CA3_SO = "CA1 SO vs CA3 SO",
    CA1_SLM_vs_mean_other_CA1_strata = "CA1 SLM vs other CA1 strata",
    CA1_SO_vs_mean_other_CA1_strata = "CA1 SO vs other CA1 strata",
    CA1_SR_vs_mean_other_CA1_strata = "CA1 SR vs other CA1 strata"
  )
  out <- unname(labels[as.character(x)])
  out[is.na(out)] <- as.character(x)[is.na(out)]
  out
}

control_spatial_figure2e_external_contrasts <- function() {
  c(
    "CA1_vs_mean_other_soma_regions", "CA2_vs_mean_other_soma_regions",
    "CA3_vs_mean_other_soma_regions", "DG_vs_mean_other_soma_regions",
    "DG_neuropil_vs_mean_non_DG_regions", "CA1_SO_vs_CA3_SO",
    "CA1_SLM_vs_mean_other_CA1_strata", "CA1_SO_vs_mean_other_CA1_strata",
    "CA1_SR_vs_mean_other_CA1_strata"
  )
}

control_spatial_figure2f_display_contrasts <- function() {
  c(
    "CA1_vs_mean_other_soma_regions", "CA2_vs_mean_other_soma_regions",
    "CA3_vs_mean_other_soma_regions", "DG_vs_mean_other_soma_regions",
    "DG_MO_vs_mean_other_DG_layers", "DG_PO_vs_mean_other_DG_layers"
  )
}

control_spatial_figure2f_regions_ca1layers_display_contrasts <- function() {
  c(
    "CA1_vs_mean_other_soma_regions", "CA2_vs_mean_other_soma_regions",
    "CA3_vs_mean_other_soma_regions", "DG_vs_mean_other_soma_regions",
    "CA1_SLM_vs_mean_other_CA1_strata",
    "CA1_SO_vs_mean_other_CA1_strata",
    "CA1_SR_vs_mean_other_CA1_strata"
  )
}

control_spatial_figure2f_grouped_layout <- function() {
  contrasts <- control_spatial_figure2f_regions_ca1layers_display_contrasts()
  data.frame(
    contrast = contrasts,
    group = c(rep("Soma region", 4L), rep("CA1 neuropil layer", 3L)),
    short_label = c("CA1", "CA2", "CA3", "DG", "SLM", "SO", "SR"),
    stringsAsFactors = FALSE
  )
}

# This is deliberately separate from the GSEA result: it records the semantic
# match contract used by Figure 2e.  In particular, the Kaulich SP signature is
# useful reference/specificity context but cannot validate an internal CA1-SP
# contrast because that spatial unit is not present in this experiment.
control_spatial_figure2e_matching_audit <- function(x) {
  required <- c("dataset", "internal_contrast", "external_signature", "match_type", "interpretation_note")
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Figure 2e mapping audit is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  out <- unique(x[, required, drop = FALSE])
  out$matched_exactly <- out$match_type == "exact"
  out$match_reason <- ifelse(
    out$matched_exactly,
    "Internal target and external reference stratum are identical.",
    ifelse(out$external_signature == "SP" & grepl("_vs_mean_other_CA1_strata$", out$internal_contrast),
      "External-reference-only CA1 SP signature: no internal CA1-SP-versus-rest contrast exists.",
      out$interpretation_note
    )
  )
  out$external_reference_context <- out$external_signature == "SP" & !out$matched_exactly
  out[order(out$internal_contrast, out$external_signature), , drop = FALSE]
}

control_spatial_figure2e_plot_data <- function(x) {
  x <- control_spatial_successful_gsea_rows(x)
  if (!nrow(x)) return(x)
  if (!"signature_FDR" %in% names(x)) stop("Figure 2e results are missing signature_FDR.", call. = FALSE)
  x$internal_contrast_label <- control_spatial_contrast_label(x$internal_contrast)
  x$outline_status <- ifelse(
    is.finite(x$signature_FDR) & x$signature_FDR < 0.05,
    "FDR < 0.05", "Not significant"
  )
  x$outline_status <- factor(x$outline_status, levels = c("FDR < 0.05", "Not significant"))
  domain_levels <- c(
    "Soma region — tissue reference",
    "Neuropil region — synaptosome reference",
    "CA1 strata — synaptosome reference"
  )
  x$validation_domain <- factor(x$validation_domain, levels = domain_levels)
  x$external_signature <- factor(
    x$external_signature,
    levels = rev(c("CA1", "CA2/3", "DG", "SLM", "SO", "SP", "SR"))
  )
  x
}

control_spatial_select_go_display <- function(x, max_terms = 2L, contrasts = NULL) {
  required <- c("dataset", "contrast", "status", "ID", "Description", "NES", "p_adjust")
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("GO-BP display input is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  max_terms <- as.integer(max_terms)
  if (length(max_terms) != 1L || is.na(max_terms) || max_terms < 1L) stop("max_terms must be positive.", call. = FALSE)
  keep <- x$status == "completed" &
    is.finite(suppressWarnings(as.numeric(x$NES))) & x$NES > 0 &
    is.finite(suppressWarnings(as.numeric(x$p_adjust))) & x$p_adjust < 0.05
  if (!is.null(contrasts)) keep <- keep & x$contrast %in% contrasts
  eligible <- x[keep, , drop = FALSE]
  if (!nrow(eligible)) return(eligible)
  keys <- paste(eligible$dataset, eligible$contrast, sep = "\r")
  key_order <- unique(keys)
  selected <- lapply(key_order, function(key) {
    z <- eligible[keys == key, , drop = FALSE]
    z <- z[order(z$p_adjust, -z$NES, as.character(z$ID)), , drop = FALSE]
    z[seq_len(min(max_terms, nrow(z))), , drop = FALSE]
  })
  out <- do.call(rbind, selected)
  rownames(out) <- NULL
  out
}

control_spatial_complete_go_display_grid <- function(complete, selected, contrasts) {
  contrasts <- as.character(contrasts)
  universe <- unique(as.character(complete$contrast))
  missing_analytical <- setdiff(contrasts, universe)
  if (length(missing_analytical)) {
    stop(
      "Completed GO-BP table is missing display contrast(s): ",
      paste(missing_analytical, collapse = ", "), call. = FALSE
    )
  }
  selected$display_evidence_status <- "fdr_supported_positive_go_bp_term"
  selected$display_term_available <- TRUE
  missing_evidence <- setdiff(contrasts, unique(as.character(selected$contrast)))
  if (length(missing_evidence)) {
    template_index <- match(missing_evidence, as.character(complete$contrast))
    placeholders <- complete[template_index, , drop = FALSE]
    placeholders$contrast <- missing_evidence
    placeholders$ID <- "__no_fdr_supported_positive_go_bp_term__"
    placeholders$Description <- "No FDR-supported positive GO-BP term"
    placeholders$NES <- NA_real_
    placeholders$p_adjust <- NA_real_
    placeholders$display_evidence_status <- "no_fdr_supported_positive_go_bp_term"
    placeholders$display_term_available <- FALSE
    selected <- dplyr::bind_rows(selected, placeholders)
  }
  counts <- table(factor(selected$contrast, levels = contrasts))
  if (!setequal(unique(as.character(selected$contrast)), contrasts) ||
      any(counts < 1L | counts > 2L)) {
    stop(
      "GO-BP display must retain every requested contrast with zero to two supported terms.",
      call. = FALSE
    )
  }
  selected
}

control_spatial_map_signature <- function(signature, canonical_official_symbols, mapped_gene_threshold = 5L) {
  required <- c("source_row_id", "source_gene_symbol_raw", "parsed_source_gene_symbol", "gene_match_key")
  missing <- setdiff(required, names(signature))
  if (length(missing)) stop("Kaulich signature mapping input is missing: ", paste(missing, collapse = ", "), call. = FALSE)

  canonical <- unique(data.frame(
    canonical_official_gene_symbol = as.character(canonical_official_symbols),
    gene_match_key = control_spatial_gene_match_key(canonical_official_symbols),
    stringsAsFactors = FALSE
  ))
  canonical_audit <- control_spatial_key_audit(
    canonical$canonical_official_gene_symbol, canonical$gene_match_key,
    "n_distinct_official_symbols", "official_symbols"
  )
  canonical <- merge(canonical, canonical_audit, by = "gene_match_key", all.x = TRUE, sort = FALSE)
  canonical$canonical_key_ambiguous <- canonical$n_distinct_official_symbols > 1L
  eligible <- canonical[nzchar(canonical$gene_match_key) & !canonical$canonical_key_ambiguous, , drop = FALSE]

  external_audit <- control_spatial_key_audit(
    signature$parsed_source_gene_symbol, signature$gene_match_key,
    "n_distinct_external_symbols", "external_symbols"
  )
  candidates <- merge(signature, external_audit, by = "gene_match_key", all.x = TRUE, sort = FALSE)
  candidates$external_key_ambiguous <- !is.na(candidates$n_distinct_external_symbols) &
    candidates$n_distinct_external_symbols > 1L
  candidates <- merge(
    candidates,
    canonical_audit,
    by = "gene_match_key", all.x = TRUE, sort = FALSE
  )
  candidates$canonical_key_ambiguous <- !is.na(candidates$n_distinct_official_symbols) &
    candidates$n_distinct_official_symbols > 1L
  candidates <- merge(
    candidates,
    eligible[, c("gene_match_key", "canonical_official_gene_symbol"), drop = FALSE],
    by = "gene_match_key", all.x = TRUE, sort = FALSE
  )
  candidates$gene_mapping_status <- ifelse(
    !nzchar(candidates$gene_match_key), "missing_external_symbol",
    ifelse(candidates$external_key_ambiguous, "ambiguous_external_key",
      ifelse(!is.na(candidates$canonical_key_ambiguous) & candidates$canonical_key_ambiguous,
        "ambiguous_canonical_key",
        ifelse(is.na(candidates$canonical_official_gene_symbol), "unmapped", "mapped")
      )
    )
  )

  external_keys <- unique(signature$gene_match_key[nzchar(signature$gene_match_key)])
  external_ambiguous <- external_audit$gene_match_key[external_audit$n_distinct_external_symbols > 1L]
  canonical_ambiguous <- canonical_audit$gene_match_key[canonical_audit$n_distinct_official_symbols > 1L]
  ambiguous_keys <- intersect(external_keys, union(external_ambiguous, canonical_ambiguous))
  mapped_keys <- unique(candidates$gene_match_key[candidates$gene_mapping_status == "mapped"])
  mapped_symbols <- unique(candidates$canonical_official_gene_symbol[candidates$gene_mapping_status == "mapped"])
  source_size <- length(unique(signature$source_row_id))
  unique_external <- length(external_keys)
  mapped_n <- length(mapped_keys)
  summary <- data.frame(
    source_signature_size_before_parsing = source_size,
    unique_external_genes_after_parsing = unique_external,
    canonical_eligible_gene_universe_size = length(unique(eligible$gene_match_key)),
    mapped_unique_genes = mapped_n,
    mapping_coverage = if (unique_external) mapped_n / unique_external else 0,
    unmapped_unique_genes = length(setdiff(external_keys, mapped_keys)),
    ambiguous_key_count = length(ambiguous_keys),
    mapping_status = if (mapped_n < mapped_gene_threshold) "below_mapping_threshold" else "mapped",
    stringsAsFactors = FALSE
  )
  list(
    candidates = candidates,
    summary = summary,
    mapped_official_symbols = mapped_symbols,
    canonical_audit = canonical_audit,
    external_audit = external_audit
  )
}

control_spatial_bind_rows_fill <- function(x) {
  x <- Filter(function(z) is.data.frame(z) && nrow(z), x)
  if (!length(x)) return(data.frame())
  columns <- unique(unlist(lapply(x, names), use.names = FALSE))
  x <- lapply(x, function(z) {
    for (column in setdiff(columns, names(z))) z[[column]] <- NA
    z[, columns, drop = FALSE]
  })
  do.call(rbind, x)
}

control_spatial_successful_gsea_rows <- function(x) {
  if (!is.data.frame(x) || !nrow(x)) return(data.frame())
  if (!"status" %in% names(x)) stop("GSEA source data is missing status.", call. = FALSE)
  successful <- x$status %in% c("completed", "completed_zero_display_terms")
  if (!any(successful)) return(x[FALSE, , drop = FALSE])
  if (!"NES" %in% names(x)) stop("Successful GSEA rows are missing NES.", call. = FALSE)
  x[successful & is.finite(suppressWarnings(as.numeric(x$NES))), , drop = FALSE]
}

control_spatial_write_empty_state_svg <- function(path, message) {
  svglite::svglite(path, width = 8, height = 5)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::grid.text(message, gp = grid::gpar(col = "#4D4D4D", fontsize = 12))
  invisible(path)
}

control_spatial_empty_status <- function(dataset, contrast, status, message) data.frame(dataset=dataset, contrast=contrast, status=status, message=message, stringsAsFactors=FALSE)

control_spatial_direct_execution <- function(nframe = sys.nframe()) {
  identical(as.integer(nframe), 0L)
}

control_spatial_validate_output_bundle <- function(paths) {
  paths <- unique(as.character(paths))
  missing <- paths[!file.exists(paths)]
  empty <- paths[file.exists(paths) & (is.na(file.info(paths)$size) | file.info(paths)$size <= 0)]
  if (length(missing) || length(empty)) {
    problems <- c(
      if (length(missing)) paste0("missing: ", missing),
      if (length(empty)) paste0("empty: ", empty)
    )
    stop(
      "Control spatial identity output bundle is incomplete: ",
      paste(problems, collapse = "; "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

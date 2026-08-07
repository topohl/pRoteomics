# Focused helpers for the manuscript-facing SUS-RES spatial DAP atlas.

if (!exists("repo_path", mode = "function")) {
  paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
  source(paths_file)
}
if (!exists("dataset_capabilities", mode = "function")) source(repo_path("R", "dataset_config.R"))
if (!exists("build_enrichment_gene_inputs", mode = "function")) source(repo_path("R", "protein_group_enrichment_utils.R"))
if (!exists("validate_clusterprofiler_manifest_contract", mode = "function")) source(repo_path("R", "enrichment_io.R"))

SUS_RES_DAP_FDR_THRESHOLD <- 0.05
MIN_ORA_DAP_GENES <- 10L

sus_res_extract_phenotype <- function(side) {
  side <- tolower(trimws(as.character(side)))
  hit <- regmatches(side, regexpr("(sus|res|con)$", side, perl = TRUE))
  if (length(hit) != 1L || is.na(hit) || !nzchar(hit)) {
    stop("Comparison side does not end in a recognized phenotype token (SUS, RES, CON): ", side, call. = FALSE)
  }
  toupper(hit)
}

sus_res_resolve_orientation <- function(comparison) {
  comparison <- as.character(comparison)
  sides <- strsplit(comparison, "_", fixed = TRUE)[[1]]
  if (length(sides) != 2L || any(!nzchar(sides))) {
    stop("SUS-RES comparison must contain exactly two non-empty sides separated by one underscore: ", comparison, call. = FALSE)
  }
  left <- sus_res_extract_phenotype(sides[[1]])
  right <- sus_res_extract_phenotype(sides[[2]])
  if (!setequal(c(left, right), c("SUS", "RES"))) {
    stop("Comparison is not a phenotype-within-unit SUS versus RES contrast: ", comparison, call. = FALSE)
  }
  multiplier <- if (identical(c(left, right), c("SUS", "RES"))) 1 else -1
  data.frame(
    comparison = comparison,
    serialized_left_phenotype = left,
    serialized_right_phenotype = right,
    serialized_effect_definition = paste(left, "-", right),
    formal_effect_definition = "SUS - RES",
    formal_effect_multiplier = multiplier,
    sign_was_flipped = multiplier == -1,
    stringsAsFactors = FALSE
  )
}

sus_res_publication_dataset_label <- function(dataset) {
  dataset <- validate_dataset(dataset)
  labels <- c(
    neuron_neuropil = "Neuropil",
    neuron_soma = "Soma",
    microglia = dataset_capabilities("microglia")$label
  )
  unname(labels[[dataset]])
}

sus_res_spatial_unit <- function(dataset, route_unit) {
  dataset <- validate_dataset(dataset)
  route_unit <- as.character(route_unit)
  if (identical(dataset, "neuron_neuropil")) return(route_unit)
  sub("_.*$", "", route_unit)
}

sus_res_resolve_manifest_input <- function(manifest_input, dataset, repository_root = repo_root()) {
  manifest_input <- as.character(manifest_input)
  if (length(manifest_input) != 1L || is.na(manifest_input) || !nzchar(manifest_input)) {
    stop("Manifest input_gene_file is missing.", call. = FALSE)
  }
  if (file.exists(manifest_input)) {
    return(normalizePath(manifest_input, winslash = "/", mustWork = TRUE))
  }
  # Historical manifests may retain an absolute path from the producing checkout.
  # Relocation is deterministic: keep the manifest-declared basename and place it
  # only in the current canonical dataset-specific mapped-contrast directory.
  relocated <- file.path(
    repository_root, "data", "processed", "02_id_mapping", "mapped", dataset,
    "forward", "per_file", basename(manifest_input)
  )
  if (!file.exists(relocated)) {
    stop("Manifest-selected mapped contrast is missing at both its recorded and canonical relocated paths: ",
      manifest_input, " | ", relocated, call. = FALSE)
  }
  normalizePath(relocated, winslash = "/", mustWork = TRUE)
}

sus_res_manifest_contrasts <- function(dataset, manifest_path, repository_root = repo_root()) {
  dataset <- validate_dataset(dataset)
  if (!file.exists(manifest_path)) stop("Canonical clusterProfiler manifest not found: ", manifest_path, call. = FALSE)
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
  validate_clusterprofiler_manifest_contract(manifest, strict = TRUE, require_files = FALSE)
  manifest <- manifest[
    as.character(manifest$dataset) == dataset &
      as.character(manifest$route_category) == "phenotype_within_unit",
    , drop = FALSE
  ]
  comparisons <- sort(unique(as.character(manifest$comparison)), method = "radix")
  keep <- vapply(comparisons, function(comparison) {
    tryCatch({ sus_res_resolve_orientation(comparison); TRUE }, error = function(e) FALSE)
  }, logical(1))
  comparisons <- comparisons[keep]
  if (!length(comparisons)) stop("No manifest-selected phenotype-within-unit SUS versus RES comparisons for ", dataset, ".", call. = FALSE)

  rows <- lapply(comparisons, function(comparison) {
    x <- manifest[manifest$comparison == comparison, , drop = FALSE]
    invariant <- c("dataset", "comparison", "route_category", "route_unit", "input_gene_file", "input_hash")
    distinct <- unique(x[, invariant, drop = FALSE])
    if (nrow(distinct) != 1L) {
      stop("Manifest rows disagree on the canonical input identity for comparison ", comparison, ".", call. = FALSE)
    }
    expected_name <- paste0(comparison, ".csv")
    if (!identical(basename(distinct$input_gene_file[[1]]), expected_name)) {
      stop("Manifest input basename does not match comparison identity for ", comparison, ".", call. = FALSE)
    }
    resolved <- sus_res_resolve_manifest_input(distinct$input_gene_file[[1]], dataset, repository_root)
    current_hash <- file_hash(resolved)
    recorded_hash <- as.character(distinct$input_hash[[1]])
    if (!is.na(recorded_hash) && nzchar(recorded_hash) && !identical(recorded_hash, current_hash)) {
      stop("Mapped contrast hash differs from the clusterProfiler manifest for ", comparison, ".", call. = FALSE)
    }
    orientation <- sus_res_resolve_orientation(comparison)
    data.frame(
      dataset = dataset,
      comparison = comparison,
      route_category = "phenotype_within_unit",
      route_unit = as.character(distinct$route_unit[[1]]),
      spatial_unit = sus_res_spatial_unit(dataset, distinct$route_unit[[1]]),
      dataset_label = sus_res_publication_dataset_label(dataset),
      manifest_file = normalizePath(manifest_path, winslash = "/", mustWork = TRUE),
      input_file_manifest = as.character(distinct$input_gene_file[[1]]),
      input_file = resolved,
      input_hash = recorded_hash,
      current_input_hash = current_hash,
      input_hash_matches_manifest = is.na(recorded_hash) || !nzchar(recorded_hash) || identical(recorded_hash, current_hash),
      orientation[, setdiff(names(orientation), "comparison"), drop = FALSE],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

sus_res_normalize_contrast <- function(df, comparison, fdr_threshold = SUS_RES_DAP_FDR_THRESHOLD) {
  require_canonical_enrichment_input(df, strict = TRUE)
  if (anyNA(df$ProteinGroupID) || any(!nzchar(trimws(as.character(df$ProteinGroupID))))) {
    stop("Missing ProteinGroupID in mapped SUS-RES contrast.", call. = FALSE)
  }
  if (anyDuplicated(df$ProteinGroupID)) {
    stop("Duplicate ProteinGroupID in mapped SUS-RES contrast: ", comparison, call. = FALSE)
  }
  fdr_col <- find_enrichment_column(df, c("padj", "adj.p.val", "fdr", "qvalue"))
  effect_col <- find_enrichment_column(df, c("log2fc", "logfc", "log2foldchange", "avg_log2fc", "avg_logfc"))
  p_col <- find_enrichment_column(df, c("pval", "p.value", "pvalue"))
  if (is.na(fdr_col)) stop("Mapped contrast lacks a supported BH-adjusted protein-level FDR column: ", comparison, call. = FALSE)
  if (is.na(effect_col)) stop("Mapped contrast lacks a supported signed log2 fold-change column: ", comparison, call. = FALSE)
  orientation <- sus_res_resolve_orientation(comparison)
  fdr <- suppressWarnings(as.numeric(df[[fdr_col]]))
  serialized_effect <- suppressWarnings(as.numeric(df[[effect_col]]))
  formal_effect <- serialized_effect * orientation$formal_effect_multiplier[[1]]
  p_value <- if (!is.na(p_col)) suppressWarnings(as.numeric(df[[p_col]])) else rep(NA_real_, nrow(df))
  tested <- is.finite(fdr) & is.finite(formal_effect)
  dap <- tested & fdr <= fdr_threshold & formal_effect != 0
  direction <- rep(NA_character_, nrow(df))
  direction[dap & formal_effect > 0] <- "Higher in SUS"
  direction[dap & formal_effect < 0] <- "Higher in RES"

  df$comparison <- comparison
  df$protein_fdr_column <- fdr_col
  df$serialized_effect_column <- effect_col
  df$serialized_effect <- serialized_effect
  df$formal_SUS_minus_RES_effect <- formal_effect
  df$protein_FDR <- fdr
  df$protein_p_value <- p_value
  df$is_tested <- tested
  df$is_DAP_FDR05 <- dap
  df$DAP_direction <- direction
  df$is_nominal_p_lt_0_05 <- is.finite(p_value) & p_value < 0.05 & is.finite(formal_effect) & formal_effect != 0
  df$is_DAP_FDR10_audit_only <- tested & fdr < 0.10 & formal_effect != 0
  df
}

sus_res_attach_gene_mapping <- function(normalized_df) {
  gene_inputs <- build_enrichment_gene_inputs(normalized_df, strict = TRUE)
  mapping <- gene_inputs$transformation
  idx <- match(normalized_df$ProteinGroupID, mapping$ProteinGroupID)
  if (anyNA(idx)) stop("Canonical gene transformation omitted a ProteinGroupID.", call. = FALSE)
  normalized_df$eligible_official_gene_symbol <- mapping$GeneSymbol[idx]
  normalized_df$eligible_official_entrez_id <- mapping$EntrezID[idx]
  normalized_df$gene_mapping_eligibility_status <- mapping$eligibility_status[idx]
  normalized_df$gene_mapping_exclusion_reason <- mapping$exclusion_reason[idx]
  list(data = normalized_df, gene_inputs = gene_inputs)
}

sus_res_count_summary <- function(membership, min_ora_genes = MIN_ORA_DAP_GENES) {
  units <- unique(membership[, c("dataset", "dataset_label", "comparison", "route_category", "route_unit", "spatial_unit"), drop = FALSE])
  rows <- lapply(seq_len(nrow(units)), function(i) {
    key <- units[i, , drop = FALSE]
    x <- membership[membership$dataset == key$dataset & membership$comparison == key$comparison, , drop = FALSE]
    query_genes <- function(direction) sort(unique(x$eligible_official_gene_symbol[
      x$is_DAP_FDR05 & x$DAP_direction == direction &
        x$gene_mapping_eligibility_status == "eligible" & !is.na(x$eligible_official_gene_symbol)
    ]), method = "radix")
    sus_genes <- query_genes("Higher in SUS")
    res_genes <- query_genes("Higher in RES")
    data.frame(
      key,
      n_tested_ProteinGroupID = sum(x$is_tested),
      n_DAP_FDR05 = sum(x$is_DAP_FDR05),
      n_higher_in_SUS = sum(x$is_DAP_FDR05 & x$DAP_direction == "Higher in SUS"),
      n_higher_in_RES = sum(x$is_DAP_FDR05 & x$DAP_direction == "Higher in RES"),
      fraction_DAP_of_tested = if (sum(x$is_tested)) sum(x$is_DAP_FDR05) / sum(x$is_tested) else NA_real_,
      n_nominal_p_lt_0_05_audit_only = sum(x$is_nominal_p_lt_0_05),
      n_DAP_FDR10_audit_only = sum(x$is_DAP_FDR10_audit_only),
      n_gene_level_eligible_DAP_genes_higher_in_SUS = length(sus_genes),
      n_gene_level_eligible_DAP_genes_higher_in_RES = length(res_genes),
      higher_in_SUS_meets_MIN_ORA_DAP_GENES = length(sus_genes) >= min_ora_genes,
      higher_in_RES_meets_MIN_ORA_DAP_GENES = length(res_genes) >= min_ora_genes,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

sus_res_pairwise_similarity <- function(a, b, unit_a, unit_b) {
  # A protein is pairwise-testable only when it has finite FDR and formal effect
  # in both units. Direction-aware DAP elements are ProteinGroupID|direction.
  common <- intersect(as.character(a$ProteinGroupID[a$is_tested]), as.character(b$ProteinGroupID[b$is_tested]))
  a_common <- a[a$ProteinGroupID %in% common, , drop = FALSE]
  b_common <- b[b$ProteinGroupID %in% common, , drop = FALSE]
  a_dap <- a_common[a_common$is_DAP_FDR05, c("ProteinGroupID", "DAP_direction"), drop = FALSE]
  b_dap <- b_common[b_common$is_DAP_FDR05, c("ProteinGroupID", "DAP_direction"), drop = FALSE]
  a_elements <- paste(a_dap$ProteinGroupID, a_dap$DAP_direction, sep = "|")
  b_elements <- paste(b_dap$ProteinGroupID, b_dap$DAP_direction, sep = "|")
  a_ids <- unique(as.character(a_dap$ProteinGroupID))
  b_ids <- unique(as.character(b_dap$ProteinGroupID))
  shared_ids <- intersect(a_ids, b_ids)
  direction_a <- setNames(as.character(a_dap$DAP_direction), as.character(a_dap$ProteinGroupID))
  direction_b <- setNames(as.character(b_dap$DAP_direction), as.character(b_dap$ProteinGroupID))
  same <- sum(direction_a[shared_ids] == direction_b[shared_ids])
  opposite <- sum(direction_a[shared_ids] != direction_b[shared_ids])
  directional_union <- union(a_elements, b_elements)
  unsigned_union <- union(a_ids, b_ids)
  data.frame(
    spatial_unit_a = as.character(unit_a),
    spatial_unit_b = as.character(unit_b),
    n_tested_common = length(common),
    n_DAP_A_in_common_universe = length(a_ids),
    n_DAP_B_in_common_universe = length(b_ids),
    n_shared_same_direction = same,
    n_shared_opposite_direction = opposite,
    direction_aware_jaccard = if (length(directional_union)) length(intersect(a_elements, b_elements)) / length(directional_union) else NA_real_,
    unsigned_protein_jaccard = if (length(unsigned_union)) length(intersect(a_ids, b_ids)) / length(unsigned_union) else NA_real_,
    direction_concordance_among_shared = if (same + opposite) same / (same + opposite) else NA_real_,
    stringsAsFactors = FALSE
  )
}

sus_res_pairwise_by_dataset <- function(membership) {
  datasets <- unique(as.character(membership$dataset))
  rows <- list()
  for (dataset in datasets) {
    d <- membership[membership$dataset == dataset, , drop = FALSE]
    units <- unique(d[, c("comparison", "spatial_unit", "dataset_label"), drop = FALSE])
    units <- units[order(units$spatial_unit, method = "radix"), , drop = FALSE]
    if (nrow(units) < 2L) next
    pairs <- utils::combn(seq_len(nrow(units)), 2L, simplify = FALSE)
    for (pair in pairs) {
      a_key <- units[pair[[1]], , drop = FALSE]
      b_key <- units[pair[[2]], , drop = FALSE]
      x <- sus_res_pairwise_similarity(
        d[d$comparison == a_key$comparison, , drop = FALSE],
        d[d$comparison == b_key$comparison, , drop = FALSE],
        a_key$spatial_unit, b_key$spatial_unit
      )
      x$dataset <- dataset
      x$dataset_label <- a_key$dataset_label
      x$comparison_a <- a_key$comparison
      x$comparison_b <- b_key$comparison
      rows[[length(rows) + 1L]] <- x[, c("dataset", "dataset_label", "comparison_a", "comparison_b", setdiff(names(x), c("dataset", "dataset_label", "comparison_a", "comparison_b"))), drop = FALSE]
    }
  }
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

sus_res_prepare_ranked_program_panel <- function(source_df, expected_datasets = NULL) {
  required <- c(
    "dataset_compartment", "dataset_compartment_label", "comparison",
    "program", "spatial_unit", "spatial_unit_label", "mean_NES",
    "significant_term_count", "min_fdr", "publication_include"
  )
  missing <- setdiff(required, names(source_df))
  if (length(missing)) {
    stop("Ranked GO/GSEA Panel C source is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  out <- source_df[as.character(source_df$comparison) == "SUS_vs_RES", , drop = FALSE]
  if (!nrow(out)) {
    stop("Ranked GO/GSEA Panel C source contains no phenotype_contrast/comparison == SUS_vs_RES rows.", call. = FALSE)
  }
  out$dataset_compartment <- as.character(out$dataset_compartment)
  if (!is.null(expected_datasets)) {
    expected_datasets <- as.character(expected_datasets)
    missing_datasets <- setdiff(expected_datasets, unique(out$dataset_compartment))
    unexpected_datasets <- setdiff(unique(out$dataset_compartment), expected_datasets)
    if (length(missing_datasets) || length(unexpected_datasets)) {
      stop(
        "Ranked GO/GSEA Panel C dataset contract failed; missing: ",
        if (length(missing_datasets)) paste(missing_datasets, collapse = ", ") else "none",
        "; unexpected: ",
        if (length(unexpected_datasets)) paste(unexpected_datasets, collapse = ", ") else "none",
        call. = FALSE
      )
    }
  }

  key <- paste(out$dataset_compartment, out$spatial_unit, out$program, sep = "|")
  if (anyDuplicated(key)) {
    stop("Ranked GO/GSEA Panel C source has duplicate dataset + spatial_unit + program rows.", call. = FALSE)
  }
  if (any(!is.finite(suppressWarnings(as.numeric(out$mean_NES))))) {
    stop("Ranked GO/GSEA Panel C source contains non-finite mean_NES values.", call. = FALSE)
  }
  if (any(!is.finite(suppressWarnings(as.numeric(out$significant_term_count))))) {
    stop("Ranked GO/GSEA Panel C source contains non-finite significant_term_count values.", call. = FALSE)
  }

  label_map <- setNames(
    vapply(unique(out$dataset_compartment), sus_res_publication_dataset_label, character(1)),
    unique(out$dataset_compartment)
  )
  out$dataset_compartment_label <- unname(label_map[out$dataset_compartment])
  out$panel <- "c"
  out$statistical_view <- "ranked proteome-wide GO/GSEA program shifts"
  out$panel_data_basis <- "ranked_proteome_wide_GO_GSEA"
  out$panel_c_derived_from_DAP_sets <- FALSE
  out$NES_interpretation <- "positive = higher in SUS; negative = higher in RES"
  rownames(out) <- NULL
  out
}

sus_res_build_ora_plan <- function(membership, min_ora_genes = MIN_ORA_DAP_GENES) {
  contrasts <- unique(membership[, c("dataset", "dataset_label", "comparison", "route_unit", "spatial_unit"), drop = FALSE])
  directions <- c("Higher in SUS", "Higher in RES")
  rows <- list()
  queries <- list()
  universes <- list()
  for (i in seq_len(nrow(contrasts))) {
    meta <- contrasts[i, , drop = FALSE]
    x <- membership[membership$dataset == meta$dataset & membership$comparison == meta$comparison, , drop = FALSE]
    eligible <- x$gene_mapping_eligibility_status == "eligible" & !is.na(x$eligible_official_gene_symbol) & nzchar(x$eligible_official_gene_symbol)
    universe <- sort(unique(as.character(x$eligible_official_gene_symbol[x$is_tested & eligible])), method = "radix")
    for (direction in directions) {
      dap <- x$is_DAP_FDR05 & x$DAP_direction == direction
      query <- sort(unique(as.character(x$eligible_official_gene_symbol[dap & eligible])), method = "radix")
      eligible_pg <- unique(as.character(x$ProteinGroupID[dap & eligible]))
      key <- paste(meta$dataset, meta$spatial_unit, direction, sep = "|")
      queries[[key]] <- query
      universes[[key]] <- universe
      rows[[length(rows) + 1L]] <- data.frame(
        meta,
        direction = direction,
        n_DAP_ProteinGroupID = sum(dap),
        n_gene_level_eligible_ProteinGroupID = length(eligible_pg),
        n_unique_query_genes = length(query),
        n_universe_genes = length(universe),
        MIN_ORA_DAP_GENES = as.integer(min_ora_genes),
        ORA_status = if (length(query) < min_ora_genes) "insufficient_DAP_genes_for_ORA" else "eligible_for_ORA",
        stringsAsFactors = FALSE
      )
    }
  }
  summary <- do.call(rbind, rows)
  rownames(summary) <- NULL
  list(summary = summary, queries = queries, universes = universes)
}

sus_res_empty_go_display <- function() {
  data.frame(
    dataset = character(), dataset_label = character(), spatial_unit = character(),
    direction = character(), ID = character(), Description = character(),
    GeneRatio = character(), BgRatio = character(), pvalue = numeric(),
    p.adjust = numeric(), qvalue = numeric(), geneID = character(), Count = integer(),
    selected_for_set = logical(), FDR_supported = logical(),
    display_selection_rule = character(), stringsAsFactors = FALSE
  )
}

sus_res_select_go_display <- function(ora_results, max_terms = 2L, fdr_threshold = 0.05) {
  required <- c("dataset", "dataset_label", "spatial_unit", "direction", "ID", "Description", "p.adjust", "Count")
  if (!nrow(ora_results)) return(sus_res_empty_go_display())
  missing <- setdiff(required, names(ora_results))
  if (length(missing)) stop("GO ORA display input is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  supported <- ora_results[is.finite(ora_results$p.adjust) & ora_results$p.adjust <= fdr_threshold, , drop = FALSE]
  if (!nrow(supported)) return(sus_res_empty_go_display())
  set_key <- paste(supported$dataset, supported$spatial_unit, supported$direction, sep = "|")
  selected_rows <- lapply(split(seq_len(nrow(supported)), set_key), function(idx) {
    x <- supported[idx, , drop = FALSE]
    x <- x[order(x$p.adjust, -x$Count, x$ID, method = "radix"), , drop = FALSE]
    head(x, max_terms)
  })
  selected <- do.call(rbind, selected_rows)
  union_keys <- unique(selected[, c("dataset", "direction", "ID"), drop = FALSE])
  display <- merge(ora_results, union_keys, by = c("dataset", "direction", "ID"), all = FALSE, sort = FALSE)
  selected_key <- paste(selected$dataset, selected$spatial_unit, selected$direction, selected$ID, sep = "|")
  display_key <- paste(display$dataset, display$spatial_unit, display$direction, display$ID, sep = "|")
  display$selected_for_set <- display_key %in% selected_key
  display$FDR_supported <- is.finite(display$p.adjust) & display$p.adjust <= fdr_threshold
  display$display_selection_rule <- "up to two GO-BP terms with BH FDR <= 0.05 per eligible spatial-unit/direction set; union compared within dataset/direction"
  display <- display[order(display$dataset, display$direction, display$Description, display$spatial_unit, method = "radix"), , drop = FALSE]
  rownames(display) <- NULL
  display
}

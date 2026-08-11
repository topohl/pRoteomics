# Focused helpers for the downstream three-contrast stress-response biological audit.

if (!exists("repo_path", mode = "function")) {
  paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
  source(paths_file)
}
if (!exists("validate_dataset", mode = "function")) source(repo_path("R", "dataset_config.R"))
if (!exists("canonical_clusterprofiler_manifest_path", mode = "function")) source(repo_path("R", "enrichment_io.R"))
if (!exists("sus_res_resolve_manifest_input", mode = "function")) source(repo_path("R", "sus_res_spatial_dap_atlas_utils.R"))

STRESS_RESPONSE_CONTRASTS <- c("RES_vs_CON", "SUS_vs_CON", "SUS_vs_RES")
STRESS_RESPONSE_CONTROL_CONTRASTS <- c("RES_vs_CON", "SUS_vs_CON")
STRESS_RESPONSE_FDR_THRESHOLD <- 0.05
# Neuropil log2FC exports retain a final decimal that yields residuals of 1e-9
# when three serialized contrasts are combined; 2e-9 is an export-precision QC
# tolerance, not a biological effect threshold.
STRESS_RESPONSE_ALGEBRA_TOLERANCE <- 2e-9
STRESS_RESPONSE_DIRECT_SIGN_GEOMETRY <- c(
  "RES_up_SUS_down", "RES_down_SUS_up", "both_up_vs_CON", "both_down_vs_CON",
  "one_or_both_exact_zero"
)

stress_response_extract_phenotype <- function(side) {
  side <- tolower(trimws(as.character(side)))
  hit <- regmatches(side, regexpr("(con|res|sus)$", side, perl = TRUE))
  if (length(hit) != 1L || is.na(hit) || !nzchar(hit)) {
    stop("Comparison side does not end in CON, RES, or SUS: ", side, call. = FALSE)
  }
  toupper(hit)
}

stress_response_parse_comparison <- function(comparison) {
  comparison <- as.character(comparison)
  sides <- strsplit(comparison, "_", fixed = TRUE)[[1]]
  if (length(sides) != 2L || any(!nzchar(sides))) {
    stop("Comparison must contain exactly two non-empty sides: ", comparison, call. = FALSE)
  }
  left <- stress_response_extract_phenotype(sides[[1]])
  right <- stress_response_extract_phenotype(sides[[2]])
  pair <- sort(c(left, right), method = "radix")
  formal <- if (identical(pair, c("CON", "RES"))) {
    c(left = "RES", right = "CON", contrast = "RES_vs_CON")
  } else if (identical(pair, c("CON", "SUS"))) {
    c(left = "SUS", right = "CON", contrast = "SUS_vs_CON")
  } else if (identical(pair, c("RES", "SUS"))) {
    c(left = "SUS", right = "RES", contrast = "SUS_vs_RES")
  } else {
    stop("Comparison is outside the bounded three-contrast family: ", comparison, call. = FALSE)
  }
  multiplier <- if (identical(c(left, right), unname(formal[c("left", "right")]))) 1 else -1
  data.frame(
    comparison = comparison,
    contrast = unname(formal[["contrast"]]),
    serialized_left_phenotype = left,
    serialized_right_phenotype = right,
    serialized_effect_definition = paste(left, "-", right),
    formal_effect_definition = paste(formal[["left"]], "-", formal[["right"]]),
    formal_effect_multiplier = multiplier,
    sign_was_flipped = multiplier == -1,
    stringsAsFactors = FALSE
  )
}

stress_response_spatial_unit <- function(dataset, route_unit) {
  dataset <- validate_dataset(dataset)
  route_unit <- as.character(route_unit)
  if (dataset == "microglia") return(sub("_microglia$", "", route_unit))
  route_unit
}

stress_response_manifest_specs <- function(dataset, manifest_path = canonical_clusterprofiler_manifest_path(dataset),
                                           repository_root = repo_root()) {
  dataset <- validate_dataset(dataset)
  if (!file.exists(manifest_path)) stop("Missing canonical clusterProfiler manifest: ", manifest_path, call. = FALSE)
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
  validate_clusterprofiler_manifest_contract(manifest, strict = TRUE, require_files = FALSE)
  manifest <- manifest[
    as.character(manifest$dataset) == dataset &
      as.character(manifest$route_category) == "phenotype_within_unit",
    , drop = FALSE
  ]
  comparisons <- sort(unique(as.character(manifest$comparison)), method = "radix")
  parsed <- lapply(comparisons, function(x) tryCatch(stress_response_parse_comparison(x), error = function(e) NULL))
  keep <- !vapply(parsed, is.null, logical(1))
  comparisons <- comparisons[keep]
  parsed <- parsed[keep]
  if (!length(comparisons)) stop("No bounded stress-response contrasts found for ", dataset, ".", call. = FALSE)
  rows <- lapply(seq_along(comparisons), function(i) {
    comparison <- comparisons[[i]]
    x <- manifest[manifest$comparison == comparison, , drop = FALSE]
    invariant <- unique(x[c("dataset", "comparison", "route_category", "route_unit", "input_gene_file", "input_hash")])
    if (nrow(invariant) != 1L) stop("Manifest input identity disagrees for ", comparison, ".", call. = FALSE)
    input_file <- sus_res_resolve_manifest_input(invariant$input_gene_file[[1]], dataset, repository_root)
    current_hash <- file_hash(input_file)
    recorded_hash <- as.character(invariant$input_hash[[1]])
    if (!is.na(recorded_hash) && nzchar(recorded_hash) && !identical(recorded_hash, current_hash)) {
      stop("Mapped contrast hash differs from the canonical manifest for ", comparison, ".", call. = FALSE)
    }
    data.frame(
      dataset = dataset,
      route_unit = as.character(invariant$route_unit[[1]]),
      spatial_unit = stress_response_spatial_unit(dataset, invariant$route_unit[[1]]),
      manifest_file = normalizePath(manifest_path, winslash = "/", mustWork = TRUE),
      input_file_manifest = as.character(invariant$input_gene_file[[1]]),
      input_file = input_file,
      manifest_input_hash = recorded_hash,
      current_input_hash = current_hash,
      input_hash_matches_manifest = is.na(recorded_hash) || !nzchar(recorded_hash) || identical(recorded_hash, current_hash),
      parsed[[i]],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  out <- do.call(rbind, rows)
  key <- paste(out$dataset, out$spatial_unit, out$contrast, sep = "|")
  if (anyDuplicated(key)) stop("Duplicate dataset + spatial_unit + formal contrast in manifest selection.", call. = FALSE)
  unit_counts <- table(paste(out$dataset, out$spatial_unit, sep = "|"))
  if (any(unit_counts != length(STRESS_RESPONSE_CONTRASTS))) {
    stop("Every selected spatial context must contain exactly the three required contrasts.", call. = FALSE)
  }
  out[order(out$dataset, out$spatial_unit, match(out$contrast, STRESS_RESPONSE_CONTRASTS), method = "radix"), , drop = FALSE]
}

stress_response_find_column <- function(df, candidates) {
  nm <- tolower(names(df))
  hit <- match(tolower(candidates), nm)
  hit <- hit[!is.na(hit)]
  if (length(hit)) names(df)[hit[[1]]] else NA_character_
}

stress_response_normalize_protein_contrast <- function(df, spec) {
  if (!is.data.frame(df) || nrow(spec) != 1L) stop("Protein normalization requires one data frame and one manifest spec row.", call. = FALSE)
  required <- "ProteinGroupID"
  if (!required %in% names(df)) stop("Mapped protein contrast lacks ProteinGroupID.", call. = FALSE)
  if (anyNA(df$ProteinGroupID) || any(!nzchar(trimws(as.character(df$ProteinGroupID))))) stop("Mapped protein contrast has missing ProteinGroupID.", call. = FALSE)
  if (anyDuplicated(df$ProteinGroupID)) stop("Mapped protein contrast has duplicate ProteinGroupID: ", spec$comparison, call. = FALSE)
  p_col <- stress_response_find_column(df, c("pval", "p.value", "pvalue"))
  fdr_col <- stress_response_find_column(df, c("padj", "adj.p.val", "fdr", "qvalue"))
  effect_col <- stress_response_find_column(df, c("log2fc", "logfc", "log2foldchange", "avg_log2fc", "avg_logfc"))
  if (any(is.na(c(p_col, fdr_col, effect_col)))) {
    stop("Mapped protein contrast must contain raw p, canonical FDR, and signed log2FC: ", spec$comparison, call. = FALSE)
  }
  p <- suppressWarnings(as.numeric(df[[p_col]]))
  fdr <- suppressWarnings(as.numeric(df[[fdr_col]]))
  effect <- suppressWarnings(as.numeric(df[[effect_col]])) * as.numeric(spec$formal_effect_multiplier[[1]])
  data.frame(
    dataset = spec$dataset[[1]],
    route_unit = spec$route_unit[[1]],
    spatial_unit = spec$spatial_unit[[1]],
    contrast = spec$contrast[[1]],
    comparison = spec$comparison[[1]],
    ProteinGroupID = as.character(df$ProteinGroupID),
    log2FC = effect,
    p_value = p,
    canonical_contrast_FDR = fdr,
    canonical_FDR_support = is.finite(fdr) & fdr < STRESS_RESPONSE_FDR_THRESHOLD,
    raw_p_input_field = p_col,
    canonical_FDR_field = fdr_col,
    effect_input_field = effect_col,
    formal_effect_multiplier = as.numeric(spec$formal_effect_multiplier[[1]]),
    formal_effect_definition = spec$formal_effect_definition[[1]],
    source_file = spec$input_file[[1]],
    source_hash = spec$current_input_hash[[1]],
    stringsAsFactors = FALSE
  )
}

stress_response_bh_validation <- function(df, group_columns, p_column, fdr_column, tolerance = 1e-12) {
  keys <- interaction(df[group_columns], drop = TRUE, lex.order = TRUE)
  groups <- split(seq_len(nrow(df)), keys)
  rows <- lapply(groups, function(idx) {
    p <- suppressWarnings(as.numeric(df[[p_column]][idx]))
    canonical <- suppressWarnings(as.numeric(df[[fdr_column]][idx]))
    reconstructed <- stats::p.adjust(p, method = "BH")
    finite <- is.finite(canonical) & is.finite(reconstructed)
    max_diff <- if (any(finite)) max(abs(canonical[finite] - reconstructed[finite])) else NA_real_
    data.frame(
      df[idx[[1]], group_columns, drop = FALSE],
      n_rows = length(idx),
      n_valid_raw_p = sum(is.finite(p)),
      max_abs_FDR_difference = max_diff,
      tolerance = tolerance,
      validation_pass = is.finite(max_diff) && max_diff <= tolerance,
      p_adjust_method = "BH",
      raw_p_input_field = p_column,
      adjusted_p_used_as_input = FALSE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

stress_response_build_protein_outputs <- function(protein_long, algebra_tolerance = STRESS_RESPONSE_ALGEBRA_TOLERANCE) {
  required <- c("dataset", "spatial_unit", "contrast", "ProteinGroupID", "log2FC", "p_value", "canonical_contrast_FDR")
  missing <- setdiff(required, names(protein_long))
  if (length(missing)) stop("Protein input is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  context_key <- paste(protein_long$dataset, protein_long$spatial_unit, sep = "|")
  contexts <- split(protein_long, context_key)
  geometry_rows <- list()
  summary_rows <- list()
  common_audit_rows <- list()
  algebra_rows <- list()
  family_rows <- list()
  k <- 0L
  for (context in contexts) {
    by_contrast <- split(context, context$contrast)
    if (!all(STRESS_RESPONSE_CONTRASTS %in% names(by_contrast))) stop("Incomplete three-contrast protein context.", call. = FALSE)
    controls <- lapply(STRESS_RESPONSE_CONTROL_CONTRASTS, function(z) by_contrast[[z]])
    names(controls) <- STRESS_RESPONSE_CONTROL_CONTRASTS
    valid_ids <- lapply(controls, function(x) x$ProteinGroupID[is.finite(x$p_value)])
    common_ids <- sort(Reduce(intersect, valid_ids), method = "radix")
    if (!length(common_ids)) stop("No common valid-p ProteinGroupID universe for ", context$dataset[[1]], " / ", context$spatial_unit[[1]], ".", call. = FALSE)
    control_common <- lapply(controls, function(x) x[match(common_ids, x$ProteinGroupID), , drop = FALSE])
    joint_p <- c(control_common[["RES_vs_CON"]]$p_value, control_common[["SUS_vs_CON"]]$p_value)
    if (any(!is.finite(joint_p))) stop("Control-pair common universe contains non-finite raw p-values.", call. = FALSE)
    joint_fdr <- stats::p.adjust(joint_p, method = "BH")
    m <- length(common_ids)
    control_common[["RES_vs_CON"]]$control_pair_joint_FDR <- joint_fdr[seq_len(m)]
    control_common[["SUS_vs_CON"]]$control_pair_joint_FDR <- joint_fdr[m + seq_len(m)]
    family_id <- paste(context$dataset[[1]], context$spatial_unit[[1]], "protein_control_pair", sep = "::")
    for (z in STRESS_RESPONSE_CONTROL_CONTRASTS) {
      control_common[[z]]$control_pair_joint_FDR05 <- control_common[[z]]$control_pair_joint_FDR < STRESS_RESPONSE_FDR_THRESHOLD
      control_common[[z]]$control_pair_family_id <- family_id
      control_common[[z]]$control_pair_family_size <- 2L * m
      control_common[[z]]$common_control_contrast_universe_size <- m
    }
    direct <- by_contrast[["SUS_vs_RES"]]
    geometry_ids <- sort(intersect(common_ids, direct$ProteinGroupID[is.finite(direct$p_value)]), method = "radix")
    res <- control_common[["RES_vs_CON"]][match(geometry_ids, control_common[["RES_vs_CON"]]$ProteinGroupID), , drop = FALSE]
    sus <- control_common[["SUS_vs_CON"]][match(geometry_ids, control_common[["SUS_vs_CON"]]$ProteinGroupID), , drop = FALSE]
    sr <- direct[match(geometry_ids, direct$ProteinGroupID), , drop = FALSE]
    g <- data.frame(
      dataset = context$dataset[[1]], spatial_unit = context$spatial_unit[[1]], ProteinGroupID = geometry_ids,
      log2FC_RES_vs_CON = res$log2FC, p_RES_vs_CON = res$p_value,
      canonical_FDR_RES_vs_CON = res$canonical_contrast_FDR,
      canonical_FDR_support_RES_vs_CON = res$canonical_FDR_support,
      control_pair_joint_FDR_RES_vs_CON = res$control_pair_joint_FDR,
      control_pair_joint_FDR05_RES_vs_CON = res$control_pair_joint_FDR05,
      log2FC_SUS_vs_CON = sus$log2FC, p_SUS_vs_CON = sus$p_value,
      canonical_FDR_SUS_vs_CON = sus$canonical_contrast_FDR,
      canonical_FDR_support_SUS_vs_CON = sus$canonical_FDR_support,
      control_pair_joint_FDR_SUS_vs_CON = sus$control_pair_joint_FDR,
      control_pair_joint_FDR05_SUS_vs_CON = sus$control_pair_joint_FDR05,
      log2FC_SUS_vs_RES = sr$log2FC, p_SUS_vs_RES = sr$p_value,
      canonical_FDR_SUS_vs_RES = sr$canonical_contrast_FDR,
      canonical_FDR_support_SUS_vs_RES = sr$canonical_FDR_support,
      control_pair_family_id = family_id, control_pair_family_size = 2L * m,
      common_control_contrast_universe_size = m,
      stringsAsFactors = FALSE
    )
    g$abs_log2FC_RES_vs_CON <- abs(g$log2FC_RES_vs_CON)
    g$abs_log2FC_SUS_vs_CON <- abs(g$log2FC_SUS_vs_CON)
    g$delta_abs_log2FC <- g$abs_log2FC_RES_vs_CON - g$abs_log2FC_SUS_vs_CON
    g$signed_RES_minus_SUS_control_response <- g$log2FC_RES_vs_CON - g$log2FC_SUS_vs_CON
    g$same_control_response_direction <- sign(g$log2FC_RES_vs_CON) == sign(g$log2FC_SUS_vs_CON)
    g$RES_farther_from_CON <- g$abs_log2FC_RES_vs_CON > g$abs_log2FC_SUS_vs_CON
    g$SUS_farther_from_CON <- g$abs_log2FC_SUS_vs_CON > g$abs_log2FC_RES_vs_CON
    g$contrast_algebra_residual <- g$log2FC_SUS_vs_RES - (g$log2FC_SUS_vs_CON - g$log2FC_RES_vs_CON)
    g$contrast_algebra_absolute_residual <- abs(g$contrast_algebra_residual)
    g$derived_geometry_inference_role <- "descriptive_only_no_p_value"
    k <- k + 1L
    geometry_rows[[k]] <- g

    full_counts <- lapply(STRESS_RESPONSE_CONTRASTS, function(z) {
      x <- by_contrast[[z]]
      c(tested = sum(is.finite(x$p_value)), canonical = sum(is.finite(x$canonical_contrast_FDR) & x$canonical_contrast_FDR < STRESS_RESPONSE_FDR_THRESHOLD))
    })
    names(full_counts) <- STRESS_RESPONSE_CONTRASTS
    ratio <- sum(res$control_pair_joint_FDR05) / sum(sus$control_pair_joint_FDR05)
    if (sum(sus$control_pair_joint_FDR05) == 0L) ratio <- NA_real_
    summary_rows[[k]] <- data.frame(
      dataset = context$dataset[[1]], spatial_unit = context$spatial_unit[[1]],
      n_common_tested_proteins = m,
      n_three_contrast_geometry_proteins = nrow(g),
      n_full_tested_RES_vs_CON = full_counts$RES_vs_CON[["tested"]],
      n_full_tested_SUS_vs_CON = full_counts$SUS_vs_CON[["tested"]],
      n_full_tested_SUS_vs_RES = full_counts$SUS_vs_RES[["tested"]],
      n_full_canonical_FDR05_RES_vs_CON = full_counts$RES_vs_CON[["canonical"]],
      n_full_canonical_FDR05_SUS_vs_CON = full_counts$SUS_vs_CON[["canonical"]],
      n_full_canonical_FDR05_SUS_vs_RES = full_counts$SUS_vs_RES[["canonical"]],
      n_canonical_FDR05_RES_vs_CON = sum(res$canonical_FDR_support),
      n_canonical_FDR05_SUS_vs_CON = sum(sus$canonical_FDR_support),
      n_canonical_FDR05_SUS_vs_RES = sum(sr$canonical_FDR_support),
      n_joint_FDR05_RES_vs_CON = sum(res$control_pair_joint_FDR05),
      n_joint_FDR05_SUS_vs_CON = sum(sus$control_pair_joint_FDR05),
      delta_n_joint_DAP = sum(res$control_pair_joint_FDR05) - sum(sus$control_pair_joint_FDR05),
      ratio_n_joint_DAP_RES_to_SUS = ratio,
      median_abs_log2FC_RES_vs_CON = stats::median(g$abs_log2FC_RES_vs_CON),
      median_abs_log2FC_SUS_vs_CON = stats::median(g$abs_log2FC_SUS_vs_CON),
      median_delta_abs_log2FC = stats::median(g$delta_abs_log2FC),
      fraction_proteins_RES_farther_from_CON = mean(g$RES_farther_from_CON),
      fraction_proteins_SUS_farther_from_CON = mean(g$SUS_farther_from_CON),
      fraction_proteins_equal_distance_from_CON = mean(!g$RES_farther_from_CON & !g$SUS_farther_from_CON),
      median_log2FC_RES_vs_CON = stats::median(g$log2FC_RES_vs_CON),
      median_log2FC_SUS_vs_CON = stats::median(g$log2FC_SUS_vs_CON),
      spearman_control_response_log2FC = suppressWarnings(stats::cor(g$log2FC_RES_vs_CON, g$log2FC_SUS_vs_CON, method = "spearman")),
      burden_summary_inference_role = "descriptive_only_no_p_value",
      stringsAsFactors = FALSE
    )
    common_audit_rows[[k]] <- data.frame(
      dataset = context$dataset[[1]], spatial_unit = context$spatial_unit[[1]],
      n_rows_RES_vs_CON = nrow(controls$RES_vs_CON), n_valid_raw_p_RES_vs_CON = length(valid_ids$RES_vs_CON),
      n_rows_SUS_vs_CON = nrow(controls$SUS_vs_CON), n_valid_raw_p_SUS_vs_CON = length(valid_ids$SUS_vs_CON),
      n_common_valid_raw_p_ProteinGroupID = m,
      n_only_RES_vs_CON = length(setdiff(valid_ids$RES_vs_CON, valid_ids$SUS_vs_CON)),
      n_only_SUS_vs_CON = length(setdiff(valid_ids$SUS_vs_CON, valid_ids$RES_vs_CON)),
      matched_universe_rule = "ProteinGroupID intersection with finite raw pval in both control contrasts",
      validation_pass = m > 0L && length(joint_p) == 2L * m,
      stringsAsFactors = FALSE
    )
    algebra_rows[[k]] <- data.frame(
      dataset = context$dataset[[1]], spatial_unit = context$spatial_unit[[1]], n_proteins = nrow(g),
      expected_identity = "SUS_vs_RES = SUS_vs_CON - RES_vs_CON",
      max_absolute_residual = max(g$contrast_algebra_absolute_residual),
      median_absolute_residual = stats::median(g$contrast_algebra_absolute_residual),
      tolerance = algebra_tolerance,
      tolerance_pass = max(g$contrast_algebra_absolute_residual) <= algebra_tolerance,
      stringsAsFactors = FALSE
    )
    family_rows[[k]] <- data.frame(
      dataset = context$dataset[[1]], spatial_unit = context$spatial_unit[[1]], family_id = family_id,
      common_universe_size = m, family_size = length(joint_p),
      raw_p_input_fields = paste(unique(c(res$raw_p_input_field, sus$raw_p_input_field)), collapse = ";"),
      adjusted_p_used_as_input = FALSE, adjustment_method = "BH",
      family_definition = "RES_vs_CON and SUS_vs_CON raw protein p-values on the matched ProteinGroupID universe; SUS_vs_RES excluded",
      validation_pass = length(joint_p) == 2L * m && all(is.finite(joint_fdr)),
      stringsAsFactors = FALSE
    )
  }
  geometry <- do.call(rbind, geometry_rows)
  summary <- do.call(rbind, summary_rows)
  common_audit <- do.call(rbind, common_audit_rows)
  algebra <- do.call(rbind, algebra_rows)
  families <- do.call(rbind, family_rows)
  rownames(geometry) <- rownames(summary) <- rownames(common_audit) <- rownames(algebra) <- rownames(families) <- NULL
  list(geometry = geometry, summary = summary, common_universe_audit = common_audit, algebra_audit = algebra, joint_family_audit = families)
}

stress_response_direct_dap_control_geometry <- function(protein_geometry) {
  required <- c(
    "dataset", "spatial_unit", "ProteinGroupID", "log2FC_RES_vs_CON", "log2FC_SUS_vs_CON",
    "abs_log2FC_RES_vs_CON", "abs_log2FC_SUS_vs_CON", "delta_abs_log2FC",
    "canonical_FDR_support_SUS_vs_RES"
  )
  missing <- setdiff(required, names(protein_geometry))
  if (length(missing)) stop("Protein geometry is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  context_columns <- unique(protein_geometry[c("dataset", "spatial_unit")])
  context_columns <- context_columns[order(context_columns$dataset, context_columns$spatial_unit, method = "radix"), , drop = FALSE]
  selected <- protein_geometry[protein_geometry$canonical_FDR_support_SUS_vs_RES %in% TRUE, , drop = FALSE]
  if (nrow(selected) && any(!is.finite(selected$log2FC_RES_vs_CON) | !is.finite(selected$log2FC_SUS_vs_CON))) {
    stop("Direct SUS-RES DAP geometry requires finite RES-vs-CON and SUS-vs-CON log2FC values.", call. = FALSE)
  }
  selected$direct_DAP_selection_basis <- "canonical SUS_vs_RES FDR < 0.05"
  selected$control_geometry_inference_role <- "descriptive decomposition of primary direct DAP; not independent control-contrast validation"
  selected$sign_geometry <- ifelse(
    selected$log2FC_RES_vs_CON == 0 | selected$log2FC_SUS_vs_CON == 0,
    "one_or_both_exact_zero",
    ifelse(
      selected$log2FC_RES_vs_CON > 0 & selected$log2FC_SUS_vs_CON < 0,
      "RES_up_SUS_down",
      ifelse(
        selected$log2FC_RES_vs_CON < 0 & selected$log2FC_SUS_vs_CON > 0,
        "RES_down_SUS_up",
        ifelse(
          selected$log2FC_RES_vs_CON > 0 & selected$log2FC_SUS_vs_CON > 0,
          "both_up_vs_CON",
          "both_down_vs_CON"
        )
      )
    )
  )
  selected$opposite_side_of_CON <- selected$sign_geometry %in% c("RES_up_SUS_down", "RES_down_SUS_up")
  selected$same_side_of_CON <- selected$sign_geometry %in% c("both_up_vs_CON", "both_down_vs_CON")
  selected$control_side_relationship <- ifelse(
    selected$opposite_side_of_CON,
    "opposite_side_of_CON",
    ifelse(selected$same_side_of_CON, "same_side_of_CON", "one_or_both_exact_zero")
  )
  selected$relative_displacement <- ifelse(
    selected$abs_log2FC_RES_vs_CON > selected$abs_log2FC_SUS_vs_CON,
    "RES_farther",
    ifelse(selected$abs_log2FC_SUS_vs_CON > selected$abs_log2FC_RES_vs_CON, "SUS_farther", "equal_distance")
  )
  selected <- selected[order(selected$dataset, selected$spatial_unit, selected$ProteinGroupID, method = "radix"), , drop = FALSE]
  rownames(selected) <- NULL
  detail_columns <- c(
    "dataset", "spatial_unit", "ProteinGroupID",
    "log2FC_RES_vs_CON", "canonical_FDR_RES_vs_CON", "canonical_FDR_support_RES_vs_CON",
    "log2FC_SUS_vs_CON", "canonical_FDR_SUS_vs_CON", "canonical_FDR_support_SUS_vs_CON",
    "log2FC_SUS_vs_RES", "canonical_FDR_SUS_vs_RES", "canonical_FDR_support_SUS_vs_RES",
    "abs_log2FC_RES_vs_CON", "abs_log2FC_SUS_vs_CON", "delta_abs_log2FC",
    "direct_DAP_selection_basis", "control_geometry_inference_role", "sign_geometry",
    "opposite_side_of_CON", "same_side_of_CON", "control_side_relationship", "relative_displacement"
  )
  detail <- selected[intersect(detail_columns, names(selected))]

  safe_fraction <- function(n, denominator) if (denominator > 0L) n / denominator else NA_real_
  safe_median <- function(x) if (length(x)) stats::median(x) else NA_real_
  summarize <- function(x, global = FALSE, dataset = NA_character_, spatial_unit = NA_character_) {
    n <- nrow(x)
    counts <- setNames(vapply(STRESS_RESPONSE_DIRECT_SIGN_GEOMETRY, function(z) sum(x$sign_geometry == z), integer(1)), STRESS_RESPONSE_DIRECT_SIGN_GEOMETRY)
    if (!global) {
      out <- data.frame(dataset = dataset, spatial_unit = spatial_unit, n_direct_DAP = n, stringsAsFactors = FALSE)
    } else {
      out <- data.frame(n_direct_SUS_RES_DAP = n, stringsAsFactors = FALSE)
    }
    for (z in STRESS_RESPONSE_DIRECT_SIGN_GEOMETRY) {
      out[[paste0("n_", z)]] <- unname(counts[[z]])
      out[[paste0("fraction_", z)]] <- safe_fraction(unname(counts[[z]]), n)
    }
    out$n_opposite_side_of_CON <- sum(x$opposite_side_of_CON)
    out$fraction_opposite_side_of_CON <- safe_fraction(out$n_opposite_side_of_CON, n)
    out$n_same_side_of_CON <- sum(x$same_side_of_CON)
    out$fraction_same_side_of_CON <- safe_fraction(out$n_same_side_of_CON, n)
    out$n_RES_farther <- sum(x$relative_displacement == "RES_farther")
    out$fraction_RES_farther <- safe_fraction(out$n_RES_farther, n)
    out$n_SUS_farther <- sum(x$relative_displacement == "SUS_farther")
    out$fraction_SUS_farther <- safe_fraction(out$n_SUS_farther, n)
    out$n_equal_distance <- sum(x$relative_displacement == "equal_distance")
    out$fraction_equal_distance <- safe_fraction(out$n_equal_distance, n)
    out$median_delta_abs_log2FC <- safe_median(x$delta_abs_log2FC)
    out$median_log2FC_RES_vs_CON <- safe_median(x$log2FC_RES_vs_CON)
    out$median_log2FC_SUS_vs_CON <- safe_median(x$log2FC_SUS_vs_CON)
    out$direct_DAP_selection_basis <- "canonical SUS_vs_RES FDR < 0.05"
    out$control_geometry_inference_role <- "descriptive decomposition of primary direct DAP; not independent control-contrast validation"
    out
  }
  context_summary <- do.call(rbind, lapply(seq_len(nrow(context_columns)), function(i) {
    hit <- selected$dataset == context_columns$dataset[[i]] & selected$spatial_unit == context_columns$spatial_unit[[i]]
    summarize(
      selected[hit, , drop = FALSE],
      dataset = context_columns$dataset[[i]], spatial_unit = context_columns$spatial_unit[[i]]
    )
  }))
  rownames(context_summary) <- NULL
  global_summary <- summarize(selected, global = TRUE)
  list(detail = detail, context_summary = context_summary, global_summary = global_summary)
}

stress_response_add_joint_go_fdr <- function(go_long) {
  required <- c("dataset", "spatial_unit", "phenotype_contrast", "ID", "pvalue", "p.adjust")
  missing <- setdiff(required, names(go_long))
  if (length(missing)) stop("GO input is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  if (anyDuplicated(paste(go_long$dataset, go_long$spatial_unit, go_long$phenotype_contrast, go_long$ID, sep = "|"))) {
    stop("GO source has duplicate dataset + spatial_unit + contrast + GO ID rows.", call. = FALSE)
  }
  contexts <- split(go_long, paste(go_long$dataset, go_long$spatial_unit, sep = "|"))
  output <- list()
  audit <- list()
  k <- 0L
  for (context in contexts) {
    controls <- lapply(STRESS_RESPONSE_CONTROL_CONTRASTS, function(z) context[context$phenotype_contrast == z, , drop = FALSE])
    names(controls) <- STRESS_RESPONSE_CONTROL_CONTRASTS
    valid <- lapply(controls, function(x) as.character(x$ID[is.finite(x$pvalue)]))
    common <- sort(Reduce(intersect, valid), method = "radix")
    if (!length(common)) stop("No common valid-p GO universe for ", context$dataset[[1]], " / ", context$spatial_unit[[1]], ".", call. = FALSE)
    aligned <- lapply(controls, function(x) x[match(common, x$ID), , drop = FALSE])
    p <- c(aligned$RES_vs_CON$pvalue, aligned$SUS_vs_CON$pvalue)
    if (any(!is.finite(p))) stop("GO control-pair common universe contains non-finite nominal p-values.", call. = FALSE)
    q <- stats::p.adjust(p, method = "BH")
    m <- length(common)
    family_id <- paste(context$dataset[[1]], context$spatial_unit[[1]], "GO_BP_control_pair", sep = "::")
    for (z in STRESS_RESPONSE_CONTROL_CONTRASTS) {
      x <- aligned[[z]]
      idx <- if (z == "RES_vs_CON") seq_len(m) else m + seq_len(m)
      x$control_pair_joint_GO_FDR <- q[idx]
      x$control_pair_joint_GO_FDR05 <- q[idx] < STRESS_RESPONSE_FDR_THRESHOLD
      x$control_pair_GO_family_id <- family_id
      x$control_pair_GO_family_size <- 2L * m
      x$common_control_contrast_GO_universe_size <- m
      x$joint_GO_adjustment_input_field <- "pvalue"
      x$adjusted_p_used_as_joint_GO_input <- FALSE
      k <- k + 1L
      output[[k]] <- x
    }
    audit[[length(audit) + 1L]] <- data.frame(
      dataset = context$dataset[[1]], spatial_unit = context$spatial_unit[[1]], family_id = family_id,
      n_valid_nominal_p_RES_vs_CON = length(valid$RES_vs_CON), n_valid_nominal_p_SUS_vs_CON = length(valid$SUS_vs_CON),
      common_GO_ID_universe_size = m, family_size = length(p),
      n_only_RES_vs_CON = length(setdiff(valid$RES_vs_CON, valid$SUS_vs_CON)),
      n_only_SUS_vs_CON = length(setdiff(valid$SUS_vs_CON, valid$RES_vs_CON)),
      nominal_p_input_field = "pvalue", adjusted_p_used_as_input = FALSE, adjustment_method = "BH",
      family_definition = "RES_vs_CON and SUS_vs_CON nominal ranked-GSEA GO-BP p-values on the matched GO-ID universe; SUS_vs_RES excluded",
      validation_pass = length(p) == 2L * m && all(is.finite(q)), stringsAsFactors = FALSE
    )
  }
  joint <- do.call(rbind, output)
  audit <- do.call(rbind, audit)
  rownames(joint) <- rownames(audit) <- NULL
  list(joint = joint, audit = audit)
}

stress_response_build_supported_go_audits <- function(go_long, mapping) {
  required <- c("dataset", "spatial_unit", "phenotype_contrast", "ID", "Description", "pvalue", "p.adjust", "NES")
  missing <- setdiff(required, names(go_long))
  if (length(missing)) stop("Supported GO audit input is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  supported <- go_long[is.finite(go_long$p.adjust) & go_long$p.adjust < STRESS_RESPONSE_FDR_THRESHOLD, , drop = FALSE]
  supported$supported_occurrence_id <- paste(
    supported$dataset, supported$spatial_unit, supported$phenotype_contrast, supported$ID, sep = "|"
  )
  if (anyDuplicated(supported$supported_occurrence_id)) {
    stop("FDR-supported GO occurrence key is not unique before theme assignment expansion.", call. = FALSE)
  }
  status <- mapping$term_status
  status_columns <- c(
    "GO_ID", "assignment_status", "manuscript_themes", "theme_roles", "registry_version", "mapping_method"
  )
  occurrence <- merge(
    supported, status[status_columns], by.x = "ID", by.y = "GO_ID", all.x = TRUE, sort = FALSE
  )
  if (anyNA(occurrence$assignment_status)) stop("Every supported GO occurrence must have an assignment status.", call. = FALSE)
  names(occurrence)[names(occurrence) == "p.adjust"] <- "canonical_contrast_GO_FDR"
  names(occurrence)[names(occurrence) == "pvalue"] <- "nominal_GSEA_p_value"
  occurrence$canonical_GO_FDR_support <- TRUE
  occurrence$canonical_GO_FDR_family <- "original ranked-GSEA contrast-specific GO-BP BH family"
  occurrence$joint_control_GO_FDR_role <- ifelse(
    occurrence$phenotype_contrast %in% STRESS_RESPONSE_CONTROL_CONTRASTS,
    "secondary_control_burden_comparison_only", "not_applicable_to_SUS_vs_RES"
  )
  occurrence$row_grain <- "FDR_supported_GO_occurrence"
  occurrence <- occurrence[order(
    occurrence$phenotype_contrast, occurrence$dataset, occurrence$spatial_unit,
    occurrence$canonical_contrast_GO_FDR, occurrence$ID, method = "radix"
  ), , drop = FALSE]
  if (anyDuplicated(occurrence$supported_occurrence_id)) stop("Unique supported GO occurrence table contains duplicate occurrence IDs.", call. = FALSE)

  assignment <- collapse_go_theme_assignment_audit(mapping)
  exploded <- merge(
    occurrence,
    assignment[c("GO_ID", "theme_id", "manuscript_theme", "theme_role", "anchor_GO_IDs", "anchor_labels")],
    by.x = "ID", by.y = "GO_ID", all.x = TRUE, sort = FALSE
  )
  assignment_key <- ifelse(is.na(exploded$theme_id) | !nzchar(exploded$theme_id), "unclassified", exploded$theme_id)
  exploded$theme_assignment_row_id <- paste(exploded$supported_occurrence_id, assignment_key, sep = "|")
  exploded$row_grain <- "FDR_supported_GO_occurrence_x_theme_assignment"
  if (anyDuplicated(exploded$theme_assignment_row_id)) stop("Supported GO theme-assignment row ID is not unique.", call. = FALSE)
  exploded <- exploded[order(
    exploded$phenotype_contrast, exploded$dataset, exploded$spatial_unit,
    exploded$canonical_contrast_GO_FDR, exploded$ID, assignment_key, method = "radix"
  ), , drop = FALSE]
  rownames(occurrence) <- rownames(exploded) <- NULL
  list(occurrence = occurrence, exploded = exploded)
}

stress_response_protected_reference_artifacts <- function(paths, producer, note, repository_root = repo_root()) {
  paths <- as.character(paths)
  if (!length(paths) || anyNA(paths) || any(!nzchar(paths))) stop("Protected reference paths must be non-empty.", call. = FALSE)
  if (length(producer) != 1L || is.na(producer) || !nzchar(producer)) stop("Protected reference producer must be one non-empty value.", call. = FALSE)
  if (length(note) != 1L || is.na(note) || !nzchar(note)) stop("Protected reference note must be one non-empty value.", call. = FALSE)
  out <- data.frame(
    artifact_path = vapply(paths, relative_to, character(1), root = repository_root),
    exists = file.exists(paths),
    sha256_at_stage11_run = vapply(paths, file_hash_sha256, character(1)),
    role = "protected_reference_not_consumed",
    producer = producer,
    note = note,
    stringsAsFactors = FALSE
  )
  out[order(out$artifact_path, method = "radix"), , drop = FALSE]
}

stress_response_build_theme_detail <- function(go_long, mapping, joint_go) {
  assignment <- collapse_go_theme_assignment_audit(mapping)
  theme_meta <- unique(mapping$registry[c("theme_id", "display_label", "theme_role", "display_order")])
  names(theme_meta)[names(theme_meta) == "display_label"] <- "manuscript_theme"
  context_meta <- unique(go_long[c("dataset", "dataset_label", "compartment", "region", "layer", "spatial_unit", "phenotype_contrast")])
  grid <- merge(context_meta, theme_meta, by = NULL, sort = FALSE)
  term_theme <- merge(
    go_long,
    assignment[c("GO_ID", "theme_id", "manuscript_theme", "theme_role", "display_order", "anchor_GO_IDs", "anchor_labels", "mapping_method", "registry_version", "GO_db_package_version", "GO_source_name", "GO_source_url", "GO_source_date", "relationship_types_approved")],
    by.x = "ID", by.y = "GO_ID", all = FALSE, sort = FALSE
  )
  joint_key <- paste(joint_go$dataset, joint_go$spatial_unit, joint_go$phenotype_contrast, joint_go$ID, sep = "|")
  term_key <- paste(term_theme$dataset, term_theme$spatial_unit, term_theme$phenotype_contrast, term_theme$ID, sep = "|")
  joint_idx <- match(term_key, joint_key)
  term_theme$control_pair_joint_GO_FDR <- joint_go$control_pair_joint_GO_FDR[joint_idx]
  term_theme$control_pair_joint_GO_FDR05 <- joint_go$control_pair_joint_GO_FDR05[joint_idx]
  groups <- split(term_theme, paste(term_theme$dataset, term_theme$spatial_unit, term_theme$phenotype_contrast, term_theme$theme_id, sep = "|"))
  rows <- lapply(seq_len(nrow(grid)), function(i) {
    key <- paste(grid$dataset[[i]], grid$spatial_unit[[i]], grid$phenotype_contrast[[i]], grid$theme_id[[i]], sep = "|")
    x <- groups[[key]]
    if (is.null(x)) x <- term_theme[FALSE, , drop = FALSE]
    nes <- suppressWarnings(as.numeric(x$NES)); nes <- nes[is.finite(nes)]
    supported <- x[is.finite(x$p.adjust) & x$p.adjust < STRESS_RESPONSE_FDR_THRESHOLD, , drop = FALSE]
    representative <- if (nrow(supported)) supported[order(supported$p.adjust, -abs(supported$NES), supported$ID, supported$Description, method = "radix"), , drop = FALSE][1, , drop = FALSE] else NULL
    joint_supported <- x[is.finite(x$control_pair_joint_GO_FDR) & x$control_pair_joint_GO_FDR < STRESS_RESPONSE_FDR_THRESHOLD, , drop = FALSE]
    median_nes <- if (length(nes)) stats::median(nes) else NA_real_
    data.frame(
      dataset = grid$dataset[[i]], dataset_label = grid$dataset_label[[i]], compartment = grid$compartment[[i]],
      region = grid$region[[i]], layer = grid$layer[[i]], spatial_unit = grid$spatial_unit[[i]],
      contrast = grid$phenotype_contrast[[i]], theme_id = grid$theme_id[[i]], manuscript_theme = grid$manuscript_theme[[i]],
      theme_role = grid$theme_role[[i]], display_order = grid$display_order[[i]],
      n_theme_terms_tested = length(nes), median_NES_all_theme_terms = median_nes,
      q25_NES = if (length(nes)) stats::quantile(nes, 0.25, names = FALSE) else NA_real_,
      q75_NES = if (length(nes)) stats::quantile(nes, 0.75, names = FALSE) else NA_real_,
      n_positive_NES = sum(nes > 0), n_negative_NES = sum(nes < 0),
      fraction_positive_NES = if (length(nes)) mean(nes > 0) else NA_real_,
      fraction_negative_NES = if (length(nes)) mean(nes < 0) else NA_real_,
      descriptive_direction = if (!length(nes)) "not_evaluable" else if (median_nes > 0) "positive_formal_contrast_tendency" else if (median_nes < 0) "negative_formal_contrast_tendency" else "neutral",
      min_canonical_GO_FDR = if (nrow(x) && any(is.finite(x$p.adjust))) min(x$p.adjust, na.rm = TRUE) else NA_real_,
      n_canonical_FDR05_GO_terms = nrow(supported), canonical_GO_FDR_support = nrow(supported) > 0L,
      representative_supported_GO_ID = if (is.null(representative)) NA_character_ else representative$ID[[1]],
      representative_supported_GO_term = if (is.null(representative)) NA_character_ else representative$Description[[1]],
      representative_supported_NES = if (is.null(representative)) NA_real_ else representative$NES[[1]],
      representative_supported_FDR = if (is.null(representative)) NA_real_ else representative$p.adjust[[1]],
      representative_leading_edge_proteins = if (is.null(representative)) NA_character_ else representative$core_enrichment[[1]],
      representative_leading_edge_genes = if (is.null(representative)) NA_character_ else representative$core_enrichment_gene[[1]],
      min_control_pair_joint_GO_FDR = if (nrow(x) && any(is.finite(x$control_pair_joint_GO_FDR))) min(x$control_pair_joint_GO_FDR, na.rm = TRUE) else NA_real_,
      n_control_pair_joint_FDR05_GO_terms = nrow(joint_supported),
      control_pair_joint_GO_FDR_support = nrow(joint_supported) > 0L,
      descriptive_summary_inference_role = "descriptive_only_no_theme_level_test",
      canonical_GO_FDR_family = "original ranked-GSEA contrast-specific GO-BP BH family",
      joint_GO_FDR_family_role = if (grid$phenotype_contrast[[i]] %in% STRESS_RESPONSE_CONTROL_CONTRASTS) "secondary_control_burden_comparison_only" else "not_applicable_to_SUS_vs_RES",
      mapping_method = "go_id_ontology", registry_version = unique(mapping$registry$registry_version)[[1]],
      GO_db_package_version = mapping$provenance$go_db_package_version,
      GO_source_name = mapping$provenance$go_source_name, GO_source_url = mapping$provenance$go_source_url,
      GO_source_date = mapping$provenance$go_source_date, relationship_types_approved = mapping$provenance$approved_relationships,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- out[order(out$display_order, out$dataset, out$spatial_unit, match(out$contrast, STRESS_RESPONSE_CONTRASTS), method = "radix"), , drop = FALSE]
  rownames(out) <- NULL
  list(detail = out, term_theme = term_theme)
}

stress_response_build_theme_trajectory <- function(theme_detail) {
  primary <- theme_detail[theme_detail$theme_role == "primary", , drop = FALSE]
  groups <- split(primary, paste(primary$dataset, primary$spatial_unit, primary$theme_id, sep = "|"))
  rows <- lapply(groups, function(x) {
    get <- function(contrast, column) {
      hit <- x[x$contrast == contrast, column, drop = TRUE]
      if (length(hit)) hit[[1]] else NA
    }
    res <- as.numeric(get("RES_vs_CON", "median_NES_all_theme_terms"))
    sus <- as.numeric(get("SUS_vs_CON", "median_NES_all_theme_terms"))
    sr <- as.numeric(get("SUS_vs_RES", "median_NES_all_theme_terms"))
    same <- if (is.finite(res) && is.finite(sus) && res != 0 && sus != 0) sign(res) == sign(sus) else NA
    res_farther <- is.finite(res) && is.finite(sus) && abs(res) > abs(sus)
    sus_farther <- is.finite(res) && is.finite(sus) && abs(sus) > abs(res)
    label <- if (isFALSE(same)) {
      "opposite-direction stress response (descriptive)"
    } else if (isTRUE(same) && res_farther) {
      "RES-biased remodeling (descriptive)"
    } else if (isTRUE(same) && sus_farther) {
      "SUS-biased remodeling (descriptive)"
    } else if (isTRUE(same)) {
      "shared-direction stress response (descriptive)"
    } else {
      "mixed/weak descriptive pattern"
    }
    data.frame(
      dataset = x$dataset[[1]], dataset_label = x$dataset_label[[1]], compartment = x$compartment[[1]],
      region = x$region[[1]], layer = x$layer[[1]], spatial_unit = x$spatial_unit[[1]],
      theme_id = x$theme_id[[1]], manuscript_theme = x$manuscript_theme[[1]],
      median_NES_RES_vs_CON = res, median_NES_SUS_vs_CON = sus, median_NES_SUS_vs_RES = sr,
      abs_median_NES_RES_vs_CON = abs(res), abs_median_NES_SUS_vs_CON = abs(sus),
      delta_abs_theme_NES = abs(res) - abs(sus),
      same_direction_RES_and_SUS_vs_CON = same,
      RES_farther_from_CON_descriptively = res_farther,
      SUS_farther_from_CON_descriptively = sus_farther,
      descriptive_geometry_label = label,
      canonical_GO_FDR_support_RES_vs_CON = as.logical(get("RES_vs_CON", "canonical_GO_FDR_support")),
      canonical_GO_FDR_support_SUS_vs_CON = as.logical(get("SUS_vs_CON", "canonical_GO_FDR_support")),
      canonical_GO_FDR_support_SUS_vs_RES = as.logical(get("SUS_vs_RES", "canonical_GO_FDR_support")),
      min_canonical_GO_FDR_RES_vs_CON = as.numeric(get("RES_vs_CON", "min_canonical_GO_FDR")),
      min_canonical_GO_FDR_SUS_vs_CON = as.numeric(get("SUS_vs_CON", "min_canonical_GO_FDR")),
      min_canonical_GO_FDR_SUS_vs_RES = as.numeric(get("SUS_vs_RES", "min_canonical_GO_FDR")),
      n_canonical_FDR05_GO_terms_RES_vs_CON = as.integer(get("RES_vs_CON", "n_canonical_FDR05_GO_terms")),
      n_canonical_FDR05_GO_terms_SUS_vs_CON = as.integer(get("SUS_vs_CON", "n_canonical_FDR05_GO_terms")),
      n_canonical_FDR05_GO_terms_SUS_vs_RES = as.integer(get("SUS_vs_RES", "n_canonical_FDR05_GO_terms")),
      control_pair_joint_GO_FDR_support_RES_vs_CON = as.logical(get("RES_vs_CON", "control_pair_joint_GO_FDR_support")),
      control_pair_joint_GO_FDR_support_SUS_vs_CON = as.logical(get("SUS_vs_CON", "control_pair_joint_GO_FDR_support")),
      min_control_pair_joint_GO_FDR_RES_vs_CON = as.numeric(get("RES_vs_CON", "min_control_pair_joint_GO_FDR")),
      min_control_pair_joint_GO_FDR_SUS_vs_CON = as.numeric(get("SUS_vs_CON", "min_control_pair_joint_GO_FDR")),
      trajectory_inference_role = "descriptive_geometry; direct support only from original SUS_vs_RES GO BH FDR",
      NES_algebra_warning = "NES values are not algebraically additive across contrasts",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- out[order(out$theme_id, out$dataset, out$spatial_unit, method = "radix"), , drop = FALSE]
  rownames(out) <- NULL
  out
}

stress_response_theme_overview <- function(theme_detail, trajectory) {
  meta <- unique(theme_detail[c("theme_id", "manuscript_theme", "theme_role", "display_order")])
  rows <- lapply(seq_len(nrow(meta)), function(i) {
    d <- theme_detail[theme_detail$theme_id == meta$theme_id[[i]], , drop = FALSE]
    t <- trajectory[trajectory$theme_id == meta$theme_id[[i]], , drop = FALSE]
    if (!nrow(t)) {
      # QC themes are deliberately not included in the trajectory contract; calculate their control geometry locally.
      wide <- split(d, paste(d$dataset, d$spatial_unit, sep = "|"))
      t <- do.call(rbind, lapply(wide, function(x) {
        res <- x$median_NES_all_theme_terms[x$contrast == "RES_vs_CON"][[1]]
        sus <- x$median_NES_all_theme_terms[x$contrast == "SUS_vs_CON"][[1]]
        data.frame(abs_median_NES_RES_vs_CON = abs(res), abs_median_NES_SUS_vs_CON = abs(sus),
          delta_abs_theme_NES = abs(res) - abs(sus),
          same_direction_RES_and_SUS_vs_CON = if (res != 0 && sus != 0) sign(res) == sign(sus) else NA,
          canonical_GO_FDR_support_SUS_vs_RES = x$canonical_GO_FDR_support[x$contrast == "SUS_vs_RES"][[1]], stringsAsFactors = FALSE)
      }))
    }
    support_count <- function(z) sum(d$contrast == z & d$canonical_GO_FDR_support, na.rm = TRUE)
    joint_count <- function(z) sum(d$contrast == z & d$control_pair_joint_GO_FDR_support, na.rm = TRUE)
    res_more <- sum(t$abs_median_NES_RES_vs_CON > t$abs_median_NES_SUS_vs_CON, na.rm = TRUE)
    sus_more <- sum(t$abs_median_NES_SUS_vs_CON > t$abs_median_NES_RES_vs_CON, na.rm = TRUE)
    same <- sum(t$same_direction_RES_and_SUS_vs_CON %in% TRUE, na.rm = TRUE)
    opposite <- sum(t$same_direction_RES_and_SUS_vs_CON %in% FALSE, na.rm = TRUE)
    deltas <- t$delta_abs_theme_NES[is.finite(t$delta_abs_theme_NES)]
    median_delta <- if (length(deltas)) stats::median(deltas) else NA_real_
    interpretation <- paste0(
      res_more, "/", length(deltas), " units were RES-farther descriptively; median delta absolute NES = ",
      if (is.finite(median_delta)) formatC(median_delta, digits = 3L, format = "f") else "NA",
      ". Count imbalance and median magnitude are descriptive and spatially resolved."
    )
    data.frame(
      domain = if (exists("sus_res_theme_domain", mode = "function")) sus_res_theme_domain(meta$theme_id[[i]]) else "Manuscript GO theme",
      theme_id = meta$theme_id[[i]], theme = meta$manuscript_theme[[i]], role = if (meta$theme_role[[i]] == "qc_review") "QC" else meta$theme_role[[i]],
      n_evaluable_spatial_units = length(unique(paste(d$dataset, d$spatial_unit, sep = "|"))),
      n_canonical_supported_RES_vs_CON_units = support_count("RES_vs_CON"),
      n_canonical_supported_SUS_vs_CON_units = support_count("SUS_vs_CON"),
      n_canonical_supported_SUS_vs_RES_units = support_count("SUS_vs_RES"),
      n_joint_control_supported_RES_vs_CON_units = joint_count("RES_vs_CON"),
      n_joint_control_supported_SUS_vs_CON_units = joint_count("SUS_vs_CON"),
      n_units_abs_RES_theme_NES_gt_abs_SUS = res_more,
      n_units_abs_SUS_theme_NES_gt_abs_RES = sus_more,
      n_units_same_response_direction = same,
      n_units_opposite_response_direction = opposite,
      median_delta_abs_theme_NES = if (length(deltas)) stats::median(deltas) else NA_real_,
      q25_delta_abs_theme_NES = if (length(deltas)) stats::quantile(deltas, 0.25, names = FALSE) else NA_real_,
      q75_delta_abs_theme_NES = if (length(deltas)) stats::quantile(deltas, 0.75, names = FALSE) else NA_real_,
      median_abs_delta_theme_NES = if (length(deltas)) stats::median(abs(deltas)) else NA_real_,
      min_delta_abs_theme_NES = if (meta$theme_role[[i]] == "primary" && length(deltas)) min(deltas) else NA_real_,
      max_delta_abs_theme_NES = if (meta$theme_role[[i]] == "primary" && length(deltas)) max(deltas) else NA_real_,
      n_RES_farther_units_with_direct_SUS_RES_FDR_support = sum(t$abs_median_NES_RES_vs_CON > t$abs_median_NES_SUS_vs_CON & t$canonical_GO_FDR_support_SUS_vs_RES, na.rm = TRUE),
      n_SUS_farther_units_with_direct_SUS_RES_FDR_support = sum(t$abs_median_NES_SUS_vs_CON > t$abs_median_NES_RES_vs_CON & t$canonical_GO_FDR_support_SUS_vs_RES, na.rm = TRUE),
      neutral_interpretation = interpretation,
      display_order = meta$display_order[[i]], stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- out[order(out$display_order, method = "radix"), , drop = FALSE]
  out$display_order <- NULL
  rownames(out) <- NULL
  out
}

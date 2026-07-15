# Canonical ProteinGroupID-to-gene transformation for differential enrichment.

canonical_enrichment_contract_version <- function() "protein_group_enrichment_v1"

find_enrichment_column <- function(df, candidates) {
  normalized <- tolower(gsub("[^a-z0-9]", "", names(df)))
  wanted <- tolower(gsub("[^a-z0-9]", "", candidates))
  hit <- match(wanted, normalized)
  hit <- hit[!is.na(hit)]
  if (length(hit)) names(df)[hit[[1]]] else NA_character_
}

require_canonical_enrichment_input <- function(df, strict = TRUE) {
  required <- c("ProteinGroupID", "original_identifier", "member_accessions",
    "member_gene_symbols", "representative_accession", "representative_gene_symbol",
    "representative_selection_rule", "protein_group_ambiguity_class",
    "gene_level_claim_allowed", "protein_level_claim_allowed", "mapping_status")
  missing <- setdiff(required, names(df))
  if (length(missing) && isTRUE(strict)) {
    stop("Canonical ProteinGroupID enrichment input is required; missing: ", paste(missing, collapse = ", "),
      ". Legacy gene_symbol input is permitted only with strict = FALSE.", call. = FALSE)
  }
  invisible(missing)
}

select_rank_statistic <- function(df) {
  moderated <- find_enrichment_column(df, c("t", "statistic", "stat", "moderated_t", "t_statistic"))
  if (!is.na(moderated)) return(list(column = moderated, type = "moderated_or_signed_inferential", fallback_used = FALSE,
    reason = "preferred moderated or signed inferential statistic"))
  fallback <- find_enrichment_column(df, c("log2fc", "logfc", "log2foldchange", "avg_log2fc", "avg_logfc"))
  if (is.na(fallback)) stop("No supported signed rank statistic (t/statistic/log2fc/logFC) found.", call. = FALSE)
  list(column = fallback, type = "log_fold_change", fallback_used = TRUE, reason = "no moderated statistic available; using signed log fold-change")
}

split_gene_symbols <- function(x) {
  x <- as.character(x); x[is.na(x)] <- ""
  out <- trimws(unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE))
  unique(out[nzchar(out)])
}

protein_group_gene_transform <- function(df, statistic = NULL, strict = TRUE) {
  missing <- require_canonical_enrichment_input(df, strict = strict)
  if (length(missing)) {
    # Explicit legacy compatibility: a row remains a single observation but is not publication-facing.
    gene_col <- find_enrichment_column(df, c("gene_symbol", "gene", "genesymbol"))
    if (is.na(gene_col)) stop("Legacy compatibility input lacks gene_symbol.", call. = FALSE)
    df$ProteinGroupID <- paste0("LEGACY:", seq_len(nrow(df)))
    df$original_identifier <- as.character(df[[gene_col]])
    df$member_accessions <- NA_character_; df$member_gene_symbols <- as.character(df[[gene_col]])
    df$representative_accession <- NA_character_; df$representative_gene_symbol <- NA_character_
    df$representative_selection_rule <- "legacy_compatibility_only"
    df$protein_group_ambiguity_class <- "single_accession_single_gene"
    df$gene_level_claim_allowed <- TRUE; df$protein_level_claim_allowed <- FALSE; df$mapping_status <- "legacy_compatibility"
  }
  if (anyNA(df$ProteinGroupID) || any(!nzchar(as.character(df$ProteinGroupID))) || anyDuplicated(df$ProteinGroupID))
    stop("Missing or duplicate ProteinGroupID; enrichment identity cannot be repaired.", call. = FALSE)
  if (is.null(statistic)) statistic <- select_rank_statistic(df)
  values <- suppressWarnings(as.numeric(df[[statistic$column]]))
  p_col <- find_enrichment_column(df, c("pvalue", "p.value", "pval")); fdr_col <- find_enrichment_column(df, c("padj", "adj.p.val", "fdr", "qvalue"))
  fc_col <- find_enrichment_column(df, c("log2fc", "logfc", "log2foldchange"))
  allowed_classes <- c("single_accession_single_gene", "multi_accession_same_gene")
  rows <- lapply(seq_len(nrow(df)), function(i) {
    cls <- as.character(df$protein_group_ambiguity_class[[i]])
    allowed <- isTRUE(df$gene_level_claim_allowed[[i]]) && cls %in% allowed_classes
    genes <- split_gene_symbols(df$member_gene_symbols[[i]])
    eligible <- allowed && length(genes) == 1L
    reason <- if (eligible) "eligible_primary_gene_level" else if (!cls %in% allowed_classes) paste0("excluded_", cls) else if (!isTRUE(df$gene_level_claim_allowed[[i]])) "gene_level_claim_not_allowed" else "same_gene_group_without_exactly_one_mapped_gene"
    data.frame(ProteinGroupID = as.character(df$ProteinGroupID[[i]]), GeneSymbol = if (eligible) genes[[1]] else NA_character_,
      source_statistic = values[[i]], source_statistic_column = statistic$column,
      source_p_value = if (!is.na(p_col)) suppressWarnings(as.numeric(df[[p_col]][[i]])) else NA_real_,
      source_fdr = if (!is.na(fdr_col)) suppressWarnings(as.numeric(df[[fdr_col]][[i]])) else NA_real_,
      source_logfc = if (!is.na(fc_col)) suppressWarnings(as.numeric(df[[fc_col]][[i]])) else NA_real_,
      protein_group_ambiguity_class = cls, gene_level_claim_allowed = isTRUE(df$gene_level_claim_allowed[[i]]),
      eligibility_status = if (eligible) "eligible" else "excluded", exclusion_reason = reason,
      original_identifier = as.character(df$original_identifier[[i]]), member_accessions = as.character(df$member_accessions[[i]]),
      member_gene_symbols = as.character(df$member_gene_symbols[[i]]), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

collapse_protein_group_genes <- function(transform) {
  eligible <- transform[transform$eligibility_status == "eligible" & is.finite(transform$source_statistic), , drop = FALSE]
  if (!nrow(eligible)) return(data.frame())
  eligible <- eligible[order(eligible$GeneSymbol, eligible$ProteinGroupID, method = "radix"), , drop = FALSE]
  groups <- split(eligible, eligible$GeneSymbol)
  out <- lapply(groups, function(x) {
    stats <- x$source_statistic
    median_optional <- function(column) {
      if (!column %in% names(x)) return(NA_real_)
      values <- x[[column]]
      if (!any(is.finite(values))) return(NA_real_)
      stats::median(values[is.finite(values)])
    }
    signs <- sign(stats[is.finite(stats)])
    data.frame(GeneSymbol = x$GeneSymbol[[1]], gene_collapse_rule = "median_finite_statistics",
      n_protein_groups_for_gene = nrow(x), contributing_ProteinGroupIDs = paste(x$ProteinGroupID, collapse = ";"),
      individual_statistics = paste(format(stats, trim = TRUE, scientific = FALSE), collapse = ";"),
      collapsed_statistic = stats::median(stats[is.finite(stats)]),
      collapsed_p_value = median_optional("source_p_value"),
      collapsed_fdr = median_optional("source_fdr"),
      collapsed_logfc = median_optional("source_logfc"),
      direction_pattern = paste(sort(unique(ifelse(signs > 0, "up", ifelse(signs < 0, "down", "zero"))), method = "radix"), collapse = ";"),
      discordant_direction = any(signs > 0) && any(signs < 0),
      contributing_original_identifiers = paste(x$original_identifier, collapse = ";"),
      contributing_member_accessions = paste(x$member_accessions, collapse = ";"),
      contributing_member_gene_symbols = paste(x$member_gene_symbols, collapse = ";"), stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

build_enrichment_gene_inputs <- function(df, strict = TRUE) {
  statistic <- select_rank_statistic(df)
  transform <- protein_group_gene_transform(df, statistic, strict)
  collapsed <- collapse_protein_group_genes(transform)
  ranked <- if (nrow(collapsed)) stats::setNames(collapsed$collapsed_statistic, collapsed$GeneSymbol)[order(collapsed$collapsed_statistic, decreasing = TRUE)] else numeric()
  sensitivity <- if (nrow(collapsed)) {
    eligible <- transform[transform$eligibility_status == "eligible", , drop = FALSE]
    maxabs <- vapply(split(eligible$source_statistic, eligible$GeneSymbol), function(x) x[which.max(abs(x))], numeric(1))
    list(median = ranked, max_abs_signed = maxabs[order(maxabs, decreasing = TRUE)], unique_groups_only = ranked[collapsed$n_protein_groups_for_gene == 1L])
  } else list(median = numeric(), max_abs_signed = numeric(), unique_groups_only = numeric())
  list(statistic = statistic, transformation = transform, collapse = collapsed, ranked = ranked, sensitivity = sensitivity)
}

build_ora_inputs <- function(collapse, p_values = NULL, fdr_values = NULL, logfc = NULL, p_threshold = 0.05, fdr_threshold = 0.05, fc_threshold = 1) {
  if (!nrow(collapse)) return(list(universe = character(), up = character(), down = character(), annotation_coverage = FALSE))
  sig <- rep(TRUE, nrow(collapse))
  if (is.null(fdr_values) && "collapsed_fdr" %in% names(collapse)) fdr_values <- collapse$collapsed_fdr
  if (is.null(p_values) && "collapsed_p_value" %in% names(collapse)) p_values <- collapse$collapsed_p_value
  if (is.null(logfc) && "collapsed_logfc" %in% names(collapse)) logfc <- collapse$collapsed_logfc
  if (!is.null(fdr_values) && any(is.finite(fdr_values))) sig <- is.finite(fdr_values) & fdr_values <= fdr_threshold
  else if (!is.null(p_values) && any(is.finite(p_values))) sig <- is.finite(p_values) & p_values <= p_threshold
  stat <- collapse$collapsed_statistic
  if (!is.null(logfc)) sig <- sig & abs(logfc) >= fc_threshold else sig <- sig & abs(stat) >= fc_threshold
  list(universe = sort(unique(collapse$GeneSymbol)), up = sort(unique(collapse$GeneSymbol[sig & stat > 0])),
    down = sort(unique(collapse$GeneSymbol[sig & stat < 0])), annotation_coverage = FALSE)
}

make_ora_input_audit <- function(genes, direction) {
  genes <- as.character(genes)
  data.frame(
    GeneSymbol = genes,
    ORA_direction = rep(as.character(direction), length(genes)),
    stringsAsFactors = FALSE
  )
}

gsea_input_diagnostics <- function(gene_list, symbol_keys) {
  genes <- names(gene_list)
  non_empty <- !is.na(genes) & nzchar(genes)
  matched <- unique(genes[non_empty & genes %in% symbol_keys])
  values <- suppressWarnings(as.numeric(gene_list))
  data.frame(
    total_ranked_genes = length(gene_list),
    non_empty_unique_gene_names = length(unique(genes[non_empty])),
    orgdb_symbol_matches = length(matched),
    orgdb_symbol_match_fraction = if (sum(non_empty) == 0) 0 else length(matched) / length(unique(genes[non_empty])),
    finite_rank_statistics = sum(is.finite(values)),
    unique_rank_statistic_values = length(unique(values[is.finite(values)])),
    duplicated_gene_names = sum(duplicated(genes[non_empty])),
    min_rank_statistic = if (any(is.finite(values))) min(values[is.finite(values)]) else NA_real_,
    max_rank_statistic = if (any(is.finite(values))) max(values[is.finite(values)]) else NA_real_,
    stringsAsFactors = FALSE
  )
}

validate_gsea_input <- function(gene_list, diagnostics, minimum_symbol_matches = 10L) {
  genes <- names(gene_list)
  values <- suppressWarnings(as.numeric(gene_list))
  valid <- is.numeric(gene_list) && !is.null(genes) && all(!is.na(genes) & nzchar(genes)) &&
    !anyDuplicated(genes) && all(is.finite(values)) && diagnostics$orgdb_symbol_matches >= minimum_symbol_matches
  if (!valid) {
    stop("GO GSEA input precondition failed: ranked_genes=", diagnostics$total_ranked_genes,
      ", non_empty_unique_gene_names=", diagnostics$non_empty_unique_gene_names,
      ", orgdb_symbol_matches=", diagnostics$orgdb_symbol_matches,
      ", finite_rank_statistics=", diagnostics$finite_rank_statistics,
      ", duplicated_gene_names=", diagnostics$duplicated_gene_names,
      ", required_orgdb_symbol_matches=", minimum_symbol_matches, call. = FALSE)
  }
  invisible(TRUE)
}

gsea_result_table <- function(gse, diagnostics) {
  details <- paste0("ranked_genes=", diagnostics$total_ranked_genes,
    ", orgdb_symbol_matches=", diagnostics$orgdb_symbol_matches)
  if (is.null(gse)) stop("GO GSEA returned NULL: ", details, call. = FALSE)
  if (!inherits(gse, "gseaResult") && !inherits(gse, "enrichResult")) {
    stop("GO GSEA returned invalid result class: ", paste(class(gse), collapse = "/"), "; ", details, call. = FALSE)
  }
  result <- if (isS4(gse)) methods::slot(gse, "result") else gse$result
  if (!is.data.frame(result)) stop("GO GSEA returned result without a usable result data.frame: ", details, call. = FALSE)
  result
}

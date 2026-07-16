# Canonical ProteinGroupID-to-gene transformation for differential enrichment.

canonical_enrichment_contract_version <- function() "protein_group_enrichment_v3_term_gene_provenance"

canonical_enrichment_gene_namespace <- function() "SYMBOL"

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
    "gene_level_claim_allowed", "protein_level_claim_allowed", "mapping_status",
    "official_gene_symbol", "official_entrez_id", "protein_group_gene_annotation_status",
    "gene_annotation_contract_version", "uniprot_mapping_file_hash", "orgdb_package_version")
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
    df$official_gene_symbol <- as.character(df[[gene_col]])
    df$official_entrez_id <- NA_character_
    df$protein_group_gene_annotation_status <- "legacy_enrichment_compatibility"
    df$gene_annotation_contract_version <- "legacy_enrichment_compatibility"
    df$uniprot_mapping_file_hash <- NA_character_
    df$orgdb_package_version <- NA_character_
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
    gene <- as.character(df$official_gene_symbol[[i]])
    annotation_status <- as.character(df$protein_group_gene_annotation_status[[i]])
    annotation_eligible <- identical(annotation_status, "concordant_official_gene") || (!isTRUE(strict) && identical(annotation_status, "legacy_enrichment_compatibility"))
    eligible <- allowed && annotation_eligible && !is.na(gene) && nzchar(gene)
    reason <- if (eligible) "eligible_primary_gene_level" else if (!cls %in% allowed_classes) paste0("excluded_", cls) else if (!isTRUE(df$gene_level_claim_allowed[[i]])) "gene_level_claim_not_allowed" else if (!annotation_eligible) paste0("excluded_gene_annotation_", annotation_status) else "missing_official_gene_symbol"
    data.frame(ProteinGroupID = as.character(df$ProteinGroupID[[i]]), GeneSymbol = if (eligible) gene else NA_character_,
      EntrezID = if (eligible) as.character(df$official_entrez_id[[i]]) else NA_character_,
      source_statistic = values[[i]], source_statistic_column = statistic$column,
      source_p_value = if (!is.na(p_col)) suppressWarnings(as.numeric(df[[p_col]][[i]])) else NA_real_,
      source_fdr = if (!is.na(fdr_col)) suppressWarnings(as.numeric(df[[fdr_col]][[i]])) else NA_real_,
      source_logfc = if (!is.na(fc_col)) suppressWarnings(as.numeric(df[[fc_col]][[i]])) else NA_real_,
      protein_group_ambiguity_class = cls, gene_level_claim_allowed = isTRUE(df$gene_level_claim_allowed[[i]]),
      protein_group_gene_annotation_status = annotation_status,
      gene_annotation_contract_version = as.character(df$gene_annotation_contract_version[[i]]),
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
    entrez <- if ("EntrezID" %in% names(x)) sort(unique(x$EntrezID[!is.na(x$EntrezID) & nzchar(x$EntrezID)]), method = "radix") else character()
    median_optional <- function(column) {
      if (!column %in% names(x)) return(NA_real_)
      values <- x[[column]]
      if (!any(is.finite(values))) return(NA_real_)
      stats::median(values[is.finite(values)])
    }
    signs <- sign(stats[is.finite(stats)])
    data.frame(GeneSymbol = x$GeneSymbol[[1]], EntrezID = if (length(entrez) == 1L) entrez[[1]] else NA_character_,
      gene_collapse_rule = "median_finite_statistics",
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

empty_enrichment_term_gene_provenance <- function() {
  data.frame(
    dataset = character(), comparison = character(), result_type = character(),
    ontology = character(), term_id = character(), term_description = character(),
    official_gene_symbol = character(), official_entrez_id = character(),
    ProteinGroupID = character(), member_accessions = character(),
    protein_group_gene_annotation_status = character(), gene_level_claim_allowed = logical(),
    rank_statistic = numeric(), core_enrichment_member = logical(),
    enrichment_contract_version = character(), gene_annotation_contract_version = character(),
    stringsAsFactors = FALSE
  )
}

build_enrichment_term_gene_provenance <- function(gsea_table, collapsed_gene_input,
                                                  transformation, dataset, comparison,
                                                  result_type, ontology, strict = TRUE,
                                                  core_identifier = c("SYMBOL", "ENTREZID")) {
  core_identifier <- match.arg(core_identifier)
  required_result <- c("ID", "Description", "core_enrichment")
  missing_result <- setdiff(required_result, names(gsea_table))
  if (length(missing_result)) {
    stop("GSEA result table is missing term-provenance columns: ",
      paste(missing_result, collapse = ", "), call. = FALSE)
  }
  required_collapse <- c("GeneSymbol", "EntrezID", "collapsed_statistic", "contributing_ProteinGroupIDs")
  missing_collapse <- setdiff(required_collapse, names(collapsed_gene_input))
  if (length(missing_collapse)) {
    stop("Collapsed gene provenance is missing: ", paste(missing_collapse, collapse = ", "), call. = FALSE)
  }
  required_transform <- c("ProteinGroupID", "GeneSymbol", "EntrezID", "member_accessions",
    "protein_group_gene_annotation_status", "gene_level_claim_allowed",
    "gene_annotation_contract_version", "eligibility_status")
  missing_transform <- setdiff(required_transform, names(transformation))
  if (length(missing_transform)) {
    stop("Protein-group transformation provenance is missing: ",
      paste(missing_transform, collapse = ", "), call. = FALSE)
  }
  if (!nrow(gsea_table)) return(empty_enrichment_term_gene_provenance())

  collapse <- collapsed_gene_input
  collapse$GeneSymbol <- as.character(collapse$GeneSymbol)
  collapse$EntrezID <- as.character(collapse$EntrezID)
  transform <- transformation[
    transformation$eligibility_status == "eligible" &
      transformation$gene_level_claim_allowed %in% TRUE,
    , drop = FALSE
  ]
  identifier_column <- if (identical(core_identifier, "SYMBOL")) "GeneSymbol" else "EntrezID"
  collapse_ids <- as.character(collapse[[identifier_column]])

  rows <- list()
  for (term_index in seq_len(nrow(gsea_table))) {
    core_values <- trimws(unlist(strsplit(as.character(gsea_table$core_enrichment[[term_index]]), "/", fixed = TRUE)))
    core_values <- sort(unique(core_values[!is.na(core_values) & nzchar(core_values)]), method = "radix")
    missing_core <- setdiff(core_values, collapse_ids)
    if (length(missing_core) && isTRUE(strict)) {
      stop("Core-enrichment identifier(s) absent from collapsed official-gene provenance for ",
        dataset, "/", comparison, "/", result_type, ": ", paste(missing_core, collapse = ", "),
        call. = FALSE)
    }
    core_values <- intersect(core_values, collapse_ids)
    for (core_value in core_values) {
      collapse_rows <- collapse[collapse_ids == core_value, , drop = FALSE]
      for (collapse_index in seq_len(nrow(collapse_rows))) {
        collapsed_row <- collapse_rows[collapse_index, , drop = FALSE]
        contributor_ids <- trimws(unlist(strsplit(
          as.character(collapsed_row$contributing_ProteinGroupIDs[[1]]), ";", fixed = TRUE
        )))
        contributor_ids <- sort(unique(contributor_ids[nzchar(contributor_ids)]), method = "radix")
        contributor_rows <- transform[transform$ProteinGroupID %in% contributor_ids, , drop = FALSE]
        missing_contributors <- setdiff(contributor_ids, contributor_rows$ProteinGroupID)
        if (length(missing_contributors) && isTRUE(strict)) {
          stop("Collapsed official gene references ineligible or missing ProteinGroupID contributor(s): ",
            paste(missing_contributors, collapse = ", "), call. = FALSE)
        }
        contributor_rows <- contributor_rows[order(contributor_rows$ProteinGroupID, method = "radix"), , drop = FALSE]
        for (contributor_index in seq_len(nrow(contributor_rows))) {
          contributor <- contributor_rows[contributor_index, , drop = FALSE]
          rows[[length(rows) + 1L]] <- data.frame(
            dataset = as.character(dataset), comparison = as.character(comparison),
            result_type = as.character(result_type), ontology = as.character(ontology),
            term_id = as.character(gsea_table$ID[[term_index]]),
            term_description = as.character(gsea_table$Description[[term_index]]),
            official_gene_symbol = as.character(collapsed_row$GeneSymbol[[1]]),
            official_entrez_id = as.character(collapsed_row$EntrezID[[1]]),
            ProteinGroupID = as.character(contributor$ProteinGroupID[[1]]),
            member_accessions = as.character(contributor$member_accessions[[1]]),
            protein_group_gene_annotation_status = as.character(contributor$protein_group_gene_annotation_status[[1]]),
            gene_level_claim_allowed = isTRUE(contributor$gene_level_claim_allowed[[1]]),
            rank_statistic = as.numeric(collapsed_row$collapsed_statistic[[1]]),
            core_enrichment_member = TRUE,
            enrichment_contract_version = canonical_enrichment_contract_version(),
            gene_annotation_contract_version = as.character(contributor$gene_annotation_contract_version[[1]]),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  if (!length(rows)) return(empty_enrichment_term_gene_provenance())
  out <- do.call(rbind, rows)
  out <- unique(out)
  ordering <- order(out$dataset, out$comparison, out$result_type, out$ontology,
    out$term_id, out$official_gene_symbol, out$ProteinGroupID, method = "radix")
  rownames(out) <- NULL
  out[ordering, names(empty_enrichment_term_gene_provenance()), drop = FALSE]
}

# Backward-compatibility resolver only. Active strict enrichment consumes Stage 02 official symbols.
resolve_enrichment_symbols <- function(collapse, symbol_keys) {
  required <- c("GeneSymbol", "collapsed_statistic", "contributing_ProteinGroupIDs")
  missing <- setdiff(required, names(collapse))
  if (length(missing)) stop("Symbol resolution input is missing: ", paste(missing, collapse = ", "), call. = FALSE)

  keys <- sort(unique(as.character(symbol_keys[!is.na(symbol_keys) & nzchar(symbol_keys)])), method = "radix")
  submitted <- as.character(collapse$GeneSymbol)
  lower_keys <- tolower(keys)
  resolved <- rep(NA_character_, length(submitted))
  status <- rep("unmatched_symbol", length(submitted))
  exact <- !is.na(submitted) & nzchar(submitted) & submitted %in% keys
  resolved[exact] <- submitted[exact]
  status[exact] <- "exact_symbol_match"

  for (i in which(!exact & !is.na(submitted) & nzchar(submitted))) {
    hits <- unique(keys[lower_keys == tolower(submitted[[i]])])
    if (length(hits) == 1L) {
      resolved[[i]] <- hits[[1]]
      status[[i]] <- "case_normalized_symbol_match"
    } else if (length(hits) > 1L) {
      status[[i]] <- "ambiguous_case_insensitive_match"
    }
  }

  included <- !is.na(resolved) & nzchar(resolved)
  audit <- data.frame(
    submitted_GeneSymbol = submitted,
    resolved_GeneSymbol = resolved,
    symbol_resolution_status = status,
    exact_match = exact,
    case_normalized = status == "case_normalized_symbol_match",
    ambiguous_case_match = status == "ambiguous_case_insensitive_match",
    included_in_enrichment = included,
    contributing_ProteinGroupIDs = as.character(collapse$contributing_ProteinGroupIDs),
    stringsAsFactors = FALSE
  )
  provenance_columns <- intersect(c("individual_statistics", "collapsed_statistic", "gene_collapse_rule",
    "contributing_original_identifiers", "contributing_member_accessions", "contributing_member_gene_symbols",
    "discordant_direction"), names(collapse))
  if (length(provenance_columns)) audit <- cbind(audit, collapse[, provenance_columns, drop = FALSE])

  resolved_collapse <- collapse[included, , drop = FALSE]
  if (!nrow(resolved_collapse)) return(list(audit = audit, collapse = resolved_collapse, ranked = numeric()))
  resolved_collapse$submitted_GeneSymbol <- resolved_collapse$GeneSymbol
  resolved_collapse$GeneSymbol <- resolved[included]
  resolved_collapse <- resolved_collapse[order(resolved_collapse$GeneSymbol, resolved_collapse$submitted_GeneSymbol, method = "radix"), , drop = FALSE]

  if (anyDuplicated(resolved_collapse$GeneSymbol)) {
    combine_text <- function(x) paste(sort(unique(as.character(x[!is.na(x) & nzchar(as.character(x))])), method = "radix"), collapse = ";")
    combine_numeric <- function(x) if (any(is.finite(x))) stats::median(x[is.finite(x)]) else NA_real_
    resolved_collapse <- do.call(rbind, lapply(split(resolved_collapse, resolved_collapse$GeneSymbol), function(x) {
      out <- x[1, , drop = FALSE]
      out$submitted_GeneSymbol <- combine_text(x$submitted_GeneSymbol)
      out$collapsed_statistic <- combine_numeric(x$collapsed_statistic)
      for (nm in intersect(c("collapsed_p_value", "collapsed_fdr", "collapsed_logfc"), names(x))) out[[nm]] <- combine_numeric(x[[nm]])
      if ("n_protein_groups_for_gene" %in% names(x)) out$n_protein_groups_for_gene <- sum(x$n_protein_groups_for_gene)
      for (nm in intersect(c("contributing_ProteinGroupIDs", "individual_statistics", "contributing_original_identifiers",
        "contributing_member_accessions", "contributing_member_gene_symbols"), names(x))) out[[nm]] <- combine_text(x[[nm]])
      if ("direction_pattern" %in% names(x)) out$direction_pattern <- combine_text(x$direction_pattern)
      if ("discordant_direction" %in% names(x)) out$discordant_direction <- any(x$discordant_direction)
      out
    }))
    rownames(resolved_collapse) <- NULL
  }

  ranked <- stats::setNames(resolved_collapse$collapsed_statistic, resolved_collapse$GeneSymbol)
  ranked <- ranked[order(ranked, decreasing = TRUE)]
  list(audit = audit, collapse = resolved_collapse, ranked = ranked)
}

validate_precomputed_enrichment_symbols <- function(collapse, symbol_keys) {
  required <- c("GeneSymbol", "contributing_ProteinGroupIDs")
  missing <- setdiff(required, names(collapse))
  if (length(missing)) stop("Official symbol validation input is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  submitted <- as.character(collapse$GeneSymbol)
  exact <- !is.na(submitted) & nzchar(submitted) & submitted %in% symbol_keys
  data.frame(
    submitted_GeneSymbol = submitted,
    resolved_GeneSymbol = ifelse(exact, submitted, NA_character_),
    symbol_resolution_status = ifelse(exact, "precomputed_exact_symbol_validated", "invalid_precomputed_official_symbol"),
    exact_match = exact,
    case_normalized = FALSE,
    ambiguous_case_match = FALSE,
    included_in_enrichment = exact,
    contributing_ProteinGroupIDs = as.character(collapse$contributing_ProteinGroupIDs),
    stringsAsFactors = FALSE
  )
}

prepare_kegg_symbol_ranks <- function(collapse, symbol_to_entrez) {
  required <- c("GeneSymbol", "collapsed_statistic", "contributing_ProteinGroupIDs")
  missing <- setdiff(required, names(collapse))
  if (length(missing)) stop("KEGG rank input is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  if (!all(c("SYMBOL", "ENTREZID") %in% names(symbol_to_entrez))) {
    stop("SYMBOL-to-ENTREZ mapping must contain SYMBOL and ENTREZID columns.", call. = FALSE)
  }

  mapping <- unique(data.frame(SYMBOL = as.character(symbol_to_entrez$SYMBOL), ENTREZID = as.character(symbol_to_entrez$ENTREZID), stringsAsFactors = FALSE))
  mapping <- mapping[!is.na(mapping$SYMBOL) & nzchar(mapping$SYMBOL) & !is.na(mapping$ENTREZID) & nzchar(mapping$ENTREZID), , drop = FALSE]
  rows <- lapply(seq_len(nrow(collapse)), function(i) {
    ids <- sort(unique(mapping$ENTREZID[mapping$SYMBOL == collapse$GeneSymbol[[i]]]), method = "radix")
    data.frame(GeneSymbol = collapse$GeneSymbol[[i]], ENTREZID = if (length(ids) == 1L) ids[[1]] else NA_character_,
      entrez_mapping_status = if (!length(ids)) "unmapped_symbol" else if (length(ids) == 1L) "unique_entrez_match" else "ambiguous_entrez_match",
      n_entrez_matches = length(ids), included_in_kegg = length(ids) == 1L,
      collapsed_statistic = collapse$collapsed_statistic[[i]],
      contributing_ProteinGroupIDs = collapse$contributing_ProteinGroupIDs[[i]], stringsAsFactors = FALSE)
  })
  audit <- if (length(rows)) do.call(rbind, rows) else data.frame(GeneSymbol = character(), ENTREZID = character(), entrez_mapping_status = character(),
    n_entrez_matches = integer(), included_in_kegg = logical(), collapsed_statistic = numeric(), contributing_ProteinGroupIDs = character(), stringsAsFactors = FALSE)
  included <- audit[audit$included_in_kegg, , drop = FALSE]
  if (!nrow(included)) return(list(ranked = numeric(), audit = audit, collapse = data.frame()))

  collapsed_entrez <- do.call(rbind, lapply(split(included, included$ENTREZID), function(x) data.frame(
    ENTREZID = x$ENTREZID[[1]],
    collapsed_statistic = stats::median(x$collapsed_statistic[is.finite(x$collapsed_statistic)]),
    gene_collapse_rule = "median_finite_statistics_for_duplicate_entrez",
    n_symbols_for_entrez = nrow(x),
    contributing_GeneSymbols = paste(sort(unique(x$GeneSymbol), method = "radix"), collapse = ";"),
    contributing_ProteinGroupIDs = paste(sort(unique(x$contributing_ProteinGroupIDs), method = "radix"), collapse = ";"),
    stringsAsFactors = FALSE)))
  rownames(collapsed_entrez) <- NULL
  ranked <- stats::setNames(collapsed_entrez$collapsed_statistic, collapsed_entrez$ENTREZID)
  ranked <- ranked[is.finite(ranked)]
  ranked <- ranked[order(ranked, decreasing = TRUE)]
  list(ranked = ranked, audit = audit, collapse = collapsed_entrez)
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

enrichment_result_table <- function(result, analysis_label, diagnostics = NULL, required = TRUE) {
  details <- if (is.null(diagnostics)) "" else paste0(": ranked_genes=", diagnostics$total_ranked_genes,
    ", orgdb_symbol_matches=", diagnostics$orgdb_symbol_matches)
  fail <- function(message) {
    full <- paste0(analysis_label, " ", message, details)
    if (isTRUE(required)) stop(full, call. = FALSE)
    warning(full, call. = FALSE)
    out <- data.frame()
    attr(out, "skipped_reason") <- full
    out
  }
  if (is.null(result)) return(fail("returned NULL"))
  if (!inherits(result, "gseaResult") && !inherits(result, "enrichResult")) {
    return(fail(paste0("returned invalid result class: ", paste(class(result), collapse = "/"))))
  }
  table <- if (isS4(result)) methods::slot(result, "result") else result$result
  if (!is.data.frame(table)) return(fail("returned no usable result data.frame"))
  table
}

gsea_result_table <- function(gse, diagnostics) {
  enrichment_result_table(gse, "GO GSEA", diagnostics = diagnostics, required = TRUE)
}

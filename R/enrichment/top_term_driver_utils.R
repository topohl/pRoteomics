# Top GO terms per comparison and the genes that drive them.
#
# Used by analysis/differential_abundance/summarize_top_term_drivers.R, which
# reads the canonical compareGO bundle (GSEA_GO, BP, phenotype_within_unit).
# The selection and the driver table follow the former compareGO tail
# (top-term selection and its section 10), now on canonical inputs: the
# leading edge of a term is its claim-safe core, i.e. the term-gene provenance
# rows with core_enrichment_member and gene_level_claim_allowed both TRUE.
#
# The protein-group helpers at the end are shared with the two WGCNA readers of
# the driver table (compare_recurrent_module_proteins.R and
# build_curated_overlap_programs.R).

TOP_TERM_DRIVER_RESULT_TYPE <- "GSEA_GO"
TOP_TERM_DRIVER_ONTOLOGY <- "BP"
TOP_TERM_DRIVER_ROUTE <- "phenotype_within_unit"
TOP_TERM_DRIVER_FDR <- 0.05

# Selection parameters. The four selection fields come from
# config/compareGO_config.yml; the per-term gene cap is the former tail's
# top_driver_genes_per_term.
top_term_driver_parameters <- function(config = list()) {
  defaults <- list(
    significant_only = TRUE,
    target_n_terms = 5L,
    redundancy_threshold = 0.7,
    min_set_size = 10L,
    top_driver_genes_per_term = 100L
  )
  config <- config[intersect(names(config), names(defaults))]
  p <- utils::modifyList(defaults, config)
  out <- list(
    significant_only = isTRUE(as.logical(p$significant_only)),
    target_n_terms = as.integer(p$target_n_terms),
    redundancy_threshold = as.numeric(p$redundancy_threshold),
    min_set_size = as.integer(p$min_set_size),
    top_driver_genes_per_term = as.integer(p$top_driver_genes_per_term),
    fdr_threshold = TOP_TERM_DRIVER_FDR
  )
  bad <- names(out)[vapply(out, function(x) length(x) != 1L || is.na(x), logical(1))]
  if (length(bad)) {
    stop("Invalid top-term driver parameter(s): ", paste(bad, collapse = ", "), call. = FALSE)
  }
  out
}

# Jaccard index of two gene sets; 0 when either set is empty, as in the former
# tail.
top_term_jaccard <- function(a, b) {
  a <- unique(as.character(a))
  b <- unique(as.character(b))
  if (!length(a) || !length(b)) return(0)
  length(intersect(a, b)) / length(union(a, b))
}

# Claim-safe core rows of the term-gene provenance, one per term, gene and
# protein group.
top_term_core_rows <- function(provenance) {
  required <- c("comparison", "term_id", "official_gene_symbol", "ProteinGroupID",
    "member_accessions", "core_enrichment_member", "gene_level_claim_allowed")
  missing <- setdiff(required, names(provenance))
  if (length(missing)) {
    stop("Term-gene provenance is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  keep <- provenance$core_enrichment_member %in% TRUE &
    provenance$gene_level_claim_allowed %in% TRUE &
    !is.na(provenance$official_gene_symbol) & nzchar(trimws(provenance$official_gene_symbol))
  if ("result_type" %in% names(provenance)) {
    keep <- keep & provenance$result_type == TOP_TERM_DRIVER_RESULT_TYPE
  }
  if ("ontology" %in% names(provenance)) {
    keep <- keep & toupper(provenance$ontology) == TOP_TERM_DRIVER_ONTOLOGY
  }
  out <- provenance[keep, c("comparison", "term_id", "official_gene_symbol",
    "ProteinGroupID", "member_accessions"), drop = FALSE]
  out$comparison <- as.character(out$comparison)
  out$term_id <- as.character(out$term_id)
  out$official_gene_symbol <- trimws(as.character(out$official_gene_symbol))
  out$ProteinGroupID <- as.character(out$ProteinGroupID)
  out$member_accessions <- as.character(out$member_accessions)
  rownames(out) <- NULL
  out
}

# Leading-edge gene sets keyed by comparison and term.
top_term_leading_edges <- function(core_rows) {
  key <- paste(core_rows$comparison, core_rows$term_id, sep = "\r")
  lapply(split(core_rows$official_gene_symbol, key), function(x) sort(unique(x), method = "radix"))
}

# Term table in the shape the selection works on.
top_term_table <- function(terms) {
  required <- c("comparison", "term_id", "term_description", "NES", "p.adjust", "setSize")
  missing <- setdiff(required, names(terms))
  if (length(missing)) {
    stop("compareGO term comparison is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  data.frame(
    comparison = as.character(terms$comparison),
    term_id = as.character(terms$term_id),
    Description = as.character(terms$term_description),
    NES = suppressWarnings(as.numeric(terms$NES)),
    p.adjust = suppressWarnings(as.numeric(terms$p.adjust)),
    setSize = suppressWarnings(as.integer(terms$setSize)),
    stringsAsFactors = FALSE
  )
}

# The instance with the largest |NES| of each term description in each
# comparison. Ties go to the first term_id.
top_term_best_instances <- function(term_table) {
  if (!nrow(term_table)) return(term_table)
  ord <- order(term_table$comparison, term_table$Description, -abs(term_table$NES),
    term_table$term_id, method = "radix")
  x <- term_table[ord, , drop = FALSE]
  x <- x[!duplicated(x[c("comparison", "Description")]), , drop = FALSE]
  rownames(x) <- NULL
  x
}

# Walk one comparison and direction in NES order and keep a term unless its
# leading edge has Jaccard > redundancy_threshold with a term already kept; stop
# at n_target kept terms. Terms below min_set_size are not considered.
top_term_select_nonredundant <- function(df, leading_edges, direction,
                                         n_target, min_set_size, redundancy_threshold) {
  direction <- match.arg(direction, c("up", "down"))
  empty <- df[0, , drop = FALSE]
  empty$direction <- character()
  empty$selection_rank <- integer()
  empty$n_leading_edge_genes <- integer()
  empty$max_jaccard_with_kept <- numeric()
  df <- df[!is.na(df$setSize) & df$setSize >= min_set_size, , drop = FALSE]
  if (!nrow(df)) return(empty)
  ord <- if (identical(direction, "up")) {
    order(-df$NES, df$Description, df$term_id, method = "radix")
  } else {
    order(df$NES, df$Description, df$term_id, method = "radix")
  }
  df <- df[ord, , drop = FALSE]
  kept <- integer()
  kept_genes <- list()
  max_jaccard <- numeric()
  for (i in seq_len(nrow(df))) {
    genes <- leading_edges[[paste(df$comparison[i], df$term_id[i], sep = "\r")]] %||% character()
    sims <- vapply(kept_genes, function(k) top_term_jaccard(genes, k), numeric(1))
    if (length(sims) && any(sims > redundancy_threshold)) next
    kept <- c(kept, i)
    kept_genes[[length(kept_genes) + 1L]] <- genes
    max_jaccard <- c(max_jaccard, if (length(sims)) max(sims) else NA_real_)
    if (length(kept) >= n_target) break
  }
  if (!length(kept)) return(empty)
  out <- df[kept, , drop = FALSE]
  out$direction <- direction
  out$selection_rank <- seq_along(kept)
  out$n_leading_edge_genes <- vapply(kept_genes, length, integer(1))
  out$max_jaccard_with_kept <- max_jaccard
  rownames(out) <- NULL
  out
}

# Per comparison, the top n up (NES > 0) and top n down (NES < 0) terms from
# the candidate pool (best instances with p.adjust < 0.05 when significant_only).
top_term_selection <- function(best, leading_edges, params) {
  pool <- best
  if (isTRUE(params$significant_only)) {
    pool <- pool[!is.na(pool$p.adjust) & pool$p.adjust < params$fdr_threshold, , drop = FALSE]
  }
  pool <- pool[!is.na(pool$NES), , drop = FALSE]
  rows <- list()
  for (comp in sort(unique(pool$comparison), method = "radix")) {
    comp_df <- pool[pool$comparison == comp, , drop = FALSE]
    for (direction in c("up", "down")) {
      side <- if (identical(direction, "up")) comp_df[comp_df$NES > 0, , drop = FALSE] else comp_df[comp_df$NES < 0, , drop = FALSE]
      rows[[length(rows) + 1L]] <- top_term_select_nonredundant(
        side, leading_edges, direction, params$target_n_terms,
        params$min_set_size, params$redundancy_threshold)
    }
  }
  if (!length(rows)) {
    return(top_term_select_nonredundant(best[0, , drop = FALSE], leading_edges, "up",
      params$target_n_terms, params$min_set_size, params$redundancy_threshold))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

empty_top_term_drivers <- function() {
  data.frame(
    Description = character(), term_id = character(), Gene = character(),
    ProteinGroupIDs = character(), member_accessions = character(),
    Freq = integer(), n_significant_comparisons = integer(), Mean_NES = numeric(),
    Comparisons = character(), stringsAsFactors = FALSE
  )
}

# The former tail's section 10. For every selected term, across all
# comparisons in which the term has a finite NES (significant or not), count
# the comparisons in which each gene is in the term's claim-safe core:
# Freq = number of comparisons, Mean_NES = mean |NES| of the term over those
# comparisons, Comparisons = those comparisons. n_significant_comparisons is how
# many of them have the term at p.adjust < 0.05; it does not change Freq. Each
# term keeps its top_driver_genes_per_term genes by Freq, then Mean_NES, then
# gene symbol.
top_term_driver_genes <- function(term_table, core_rows, top_terms, params) {
  if (!length(top_terms) || !nrow(core_rows)) return(empty_top_term_drivers())
  terms <- term_table[term_table$Description %in% top_terms & is.finite(term_table$NES),
    c("comparison", "term_id", "Description", "NES", "p.adjust"), drop = FALSE]
  core <- merge(core_rows, terms, by = c("comparison", "term_id"), sort = FALSE)
  if (!nrow(core)) return(empty_top_term_drivers())
  core$significant <- !is.na(core$p.adjust) & core$p.adjust < params$fdr_threshold
  names(core)[names(core) == "official_gene_symbol"] <- "Gene"
  core <- core[order(core$Description, core$Gene, core$comparison, core$ProteinGroupID,
    method = "radix"), , drop = FALSE]

  ## one row per term, gene and comparison: a gene carried by two protein
  ## groups still counts once per comparison
  per_comparison <- core[!duplicated(core[c("Description", "Gene", "comparison")]), , drop = FALSE]
  stats <- per_comparison |>
    dplyr::group_by(.data$Description, .data$Gene) |>
    dplyr::summarise(
      term_id = paste(sort(unique(.data$term_id), method = "radix"), collapse = ";"),
      Freq = dplyr::n_distinct(.data$comparison),
      n_significant_comparisons = dplyr::n_distinct(.data$comparison[.data$significant]),
      Mean_NES = mean(abs(.data$NES)),
      Comparisons = paste(sort(unique(.data$comparison), method = "radix"), collapse = "; "),
      .groups = "drop"
    )
  groups <- core |>
    dplyr::group_by(.data$Description, .data$Gene) |>
    dplyr::summarise(
      encoded = list(top_term_collapse_protein_groups(.data$ProteinGroupID, .data$member_accessions)),
      .groups = "drop"
    )
  groups$ProteinGroupIDs <- vapply(groups$encoded, `[[`, character(1), "ProteinGroupIDs")
  groups$member_accessions <- vapply(groups$encoded, `[[`, character(1), "member_accessions")
  out <- merge(as.data.frame(stats), as.data.frame(groups[c("Description", "Gene", "ProteinGroupIDs", "member_accessions")]),
    by = c("Description", "Gene"), sort = FALSE)
  out <- out[c("Description", "term_id", "Gene", "ProteinGroupIDs", "member_accessions", "Freq",
    "n_significant_comparisons", "Mean_NES", "Comparisons")]
  out <- out[order(out$Description, -out$Freq, -out$Mean_NES, out$Gene, method = "radix"), , drop = FALSE]
  rank_in_term <- stats::ave(seq_len(nrow(out)), out$Description, FUN = seq_along)
  out <- out[rank_in_term <= params$top_driver_genes_per_term, , drop = FALSE]
  out$Freq <- as.integer(out$Freq)
  out$n_significant_comparisons <- as.integer(out$n_significant_comparisons)
  rownames(out) <- NULL
  out
}

# ---- protein groups ---------------------------------------------------------
#
# A gene can be carried by more than one protein group. The driver table keeps
# them all: ProteinGroupIDs lists the groups (";"-separated, sorted) and
# member_accessions lists each group's member accessions in the same order,
# groups separated by " | " and accessions within a group by ";". A gene with
# one protein group therefore reads exactly as in the provenance.

TOP_TERM_PROTEIN_GROUP_SEPARATOR <- " | "

top_term_split_accessions <- function(x) {
  x <- trimws(unlist(strsplit(as.character(x), ";", fixed = TRUE), use.names = FALSE))
  x[!is.na(x) & nzchar(x)]
}

# Decode one or more encoded (ProteinGroupIDs, member_accessions) pairs into
# one row per protein group.
top_term_protein_group_pairs <- function(protein_group_ids, member_accessions) {
  protein_group_ids <- as.character(protein_group_ids)
  member_accessions <- as.character(member_accessions)
  if (length(protein_group_ids) != length(member_accessions)) {
    stop("ProteinGroupIDs and member_accessions differ in length.", call. = FALSE)
  }
  rows <- lapply(seq_along(protein_group_ids), function(i) {
    pg <- top_term_split_accessions(protein_group_ids[i])
    if (!length(pg)) return(NULL)
    acc <- if (is.na(member_accessions[i])) character() else
      strsplit(member_accessions[i], TOP_TERM_PROTEIN_GROUP_SEPARATOR, fixed = TRUE)[[1]]
    if (length(acc) > length(pg)) {
      stop("member_accessions lists more protein groups than ProteinGroupIDs: ",
        protein_group_ids[i], call. = FALSE)
    }
    acc <- c(acc, rep("", length(pg) - length(acc)))
    data.frame(ProteinGroupID = pg, member_accessions = trimws(acc), stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) {
    return(data.frame(ProteinGroupID = character(), member_accessions = character(),
      stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

# Encode protein groups (one entry per row, duplicates allowed) as the aligned
# ProteinGroupIDs / member_accessions pair.
top_term_collapse_protein_groups <- function(protein_group_id, member_accessions) {
  pg <- as.character(protein_group_id)
  acc <- as.character(member_accessions)
  keep <- !is.na(pg) & nzchar(pg)
  pg <- pg[keep]
  acc <- acc[keep]
  if (!length(pg)) return(c(ProteinGroupIDs = NA_character_, member_accessions = NA_character_))
  first <- !duplicated(pg)
  pg <- pg[first]
  acc <- acc[first]
  ord <- order(pg, method = "radix")
  acc <- vapply(acc[ord], function(a) paste(top_term_split_accessions(a), collapse = ";"),
    character(1), USE.NAMES = FALSE)
  c(ProteinGroupIDs = paste(pg[ord], collapse = ";"),
    member_accessions = paste(acc, collapse = TOP_TERM_PROTEIN_GROUP_SEPARATOR))
}

# Merge several encoded pairs (for example all rows of one gene) into one.
top_term_merge_protein_groups <- function(protein_group_ids, member_accessions) {
  pairs <- top_term_protein_group_pairs(protein_group_ids, member_accessions)
  top_term_collapse_protein_groups(pairs$ProteinGroupID, pairs$member_accessions)
}

# Representative accession of a protein group: its first member accession in
# canonical (sorted) order, the fallback rule of select_representative_member()
# in R/data_contracts/protein_mapping_utils.R.
top_term_representative_accession <- function(member_accessions) {
  vapply(as.character(member_accessions), function(a) {
    acc <- top_term_split_accessions(a)
    if (!length(acc)) NA_character_ else sort(unique(acc), method = "radix")[[1]]
  }, character(1), USE.NAMES = FALSE)
}

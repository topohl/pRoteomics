# Helpers for representative WGCNA module and supermodule GO comparisons.
#
# Module heatmaps display FDR-significant ORA terms as -log10(BH FDR).
# Supermodule values are recurrence-weighted member-module evidence, not a
# pooled enrichment test on the union of module proteins.

wgcna_go_score <- function(p_adjust, fdr_cutoff = 0.05, score_cap = 6) {
  p_adjust <- suppressWarnings(as.numeric(p_adjust))
  significant <- is.finite(p_adjust) & p_adjust <= fdr_cutoff
  score <- rep(0, length(p_adjust))
  score[significant] <- pmin(-log10(pmax(p_adjust[significant], .Machine$double.xmin)), score_cap)
  score
}

validate_wgcna_module_go_enrichment <- function(x, artifact = "WGCNA module GO enrichment") {
  required <- c("ModuleID", "ModuleProteinSetType", "Ontology", "ID", "Description", "p.adjust")
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop(artifact, " is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

select_representative_go_terms <- function(
    x,
    entity_col = "EntityID",
    top_terms_per_entity = 2L,
    max_terms = 30L
) {
  required <- c(entity_col, "TermID", "Description", "EnrichmentScore", "Significant")
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Representative GO selection is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  if (!nrow(x)) return(data.frame())

  x <- x[!is.na(x$TermID) & nzchar(x$TermID) & !is.na(x$Description) & nzchar(x$Description), , drop = FALSE]
  if (!nrow(x)) return(data.frame())
  x$EnrichmentScore <- suppressWarnings(as.numeric(x$EnrichmentScore))
  x$EnrichmentScore[!is.finite(x$EnrichmentScore)] <- 0
  x$Significant <- as.logical(x$Significant)
  selection_basis <- if (any(x$Significant, na.rm = TRUE)) "FDR_significant" else "lowest_FDR_fallback"
  candidate <- if (identical(selection_basis, "FDR_significant")) x[x$Significant, , drop = FALSE] else x
  if (!nrow(candidate)) return(data.frame())

  candidate <- candidate[order(candidate[[entity_col]], -candidate$EnrichmentScore, candidate$TermID, candidate$Description), , drop = FALSE]
  candidate$entity_term_rank <- ave(candidate$EnrichmentScore, candidate[[entity_col]], FUN = function(z) seq_along(z))
  candidate <- candidate[candidate$entity_term_rank <= as.integer(top_terms_per_entity), , drop = FALSE]
  if (!nrow(candidate)) return(data.frame())

  key <- paste(candidate$TermID, candidate$Description, sep = "\r")
  summary_rows <- lapply(split(candidate, key), function(one) {
    data.frame(
      TermID = one$TermID[[1]], Description = one$Description[[1]],
      n_supporting_entities = length(unique(one[[entity_col]])),
      max_enrichment_score = max(one$EnrichmentScore, na.rm = TRUE),
      mean_enrichment_score = mean(one$EnrichmentScore, na.rm = TRUE),
      best_entity_rank = min(one$entity_term_rank, na.rm = TRUE),
      selection_basis = selection_basis, stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, summary_rows)
  out <- out[order(-out$n_supporting_entities, -out$max_enrichment_score, -out$mean_enrichment_score, out$best_entity_rank, out$Description, out$TermID), , drop = FALSE]
  out <- utils::head(out, as.integer(max_terms))
  out$selection_rank <- seq_len(nrow(out))
  rownames(out) <- NULL
  out
}

make_full_module_go_term_matrix <- function(
    enrichment,
    module_universe,
    ontology = "BP",
    fdr_cutoff = 0.05,
    score_cap = 6
) {
  validate_wgcna_module_go_enrichment(enrichment)
  if (!"ModuleID" %in% names(module_universe)) stop("Module universe must contain ModuleID.", call. = FALSE)
  modules <- unique(as.character(module_universe$ModuleID))
  modules <- modules[!is.na(modules) & nzchar(modules)]
  if (!length(modules)) stop("Module universe contains no ModuleID values.", call. = FALSE)

  go <- enrichment[
    as.character(enrichment$ModuleProteinSetType) == "all" &
      as.character(enrichment$Ontology) == ontology &
      as.character(enrichment$ModuleID) %in% modules,
    , drop = FALSE
  ]
  if (!nrow(go)) return(list(matrix = data.frame(), status = "no_terms_for_ontology"))
  go$ModuleID <- as.character(go$ModuleID)
  go$TermID <- as.character(go$ID)
  go$Description <- as.character(go$Description)
  go$p_adjust <- suppressWarnings(as.numeric(go$p.adjust))
  go$EnrichmentScore <- wgcna_go_score(go$p_adjust, fdr_cutoff = fdr_cutoff, score_cap = score_cap)
  go$Significant <- is.finite(go$p_adjust) & go$p_adjust <= fdr_cutoff
  go <- go[order(go$ModuleID, go$TermID, go$p_adjust), , drop = FALSE]
  go <- go[!duplicated(go[c("ModuleID", "TermID")]), , drop = FALSE]

  terms <- unique(go[c("TermID", "Description")])
  grid <- merge(
    expand.grid(TermID = terms$TermID, ModuleID = modules, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE),
    terms,
    by = "TermID", all.x = TRUE, sort = FALSE
  )
  retained <- c("ModuleID", "TermID", "p_adjust", "EnrichmentScore", "Significant", "GeneRatio", "Count", "MappedModuleSize", "ModuleSize")
  for (column in setdiff(retained, names(go))) go[[column]] <- NA
  grid <- merge(grid, go[retained], by = c("TermID", "ModuleID"), all.x = TRUE, sort = FALSE)
  grid$p_adjust[is.na(grid$p_adjust)] <- 1
  grid$EnrichmentScore[is.na(grid$EnrichmentScore)] <- 0
  grid$Significant[is.na(grid$Significant)] <- FALSE
  grid$Ontology <- ontology
  grid$FDRCutoff <- fdr_cutoff
  grid$ScoreCap <- score_cap
  grid$enrichment_measure <- "capped_-log10_BH_FDR; zero = not FDR-significant"
  list(matrix = grid, status = "ok")
}

make_module_go_heatmap_data <- function(
    enrichment,
    module_universe,
    ontology = "BP",
    fdr_cutoff = 0.05,
    score_cap = 6,
    top_terms_per_module = 2L,
    max_terms = 30L
) {
  full <- make_full_module_go_term_matrix(
    enrichment, module_universe, ontology = ontology,
    fdr_cutoff = fdr_cutoff, score_cap = score_cap
  )
  if (!nrow(full$matrix)) return(list(matrix = data.frame(), selection = data.frame(), status = full$status))
  selection <- select_representative_go_terms(
    data.frame(
      EntityID = full$matrix$ModuleID, TermID = full$matrix$TermID,
      Description = full$matrix$Description, EnrichmentScore = full$matrix$EnrichmentScore,
      Significant = full$matrix$Significant, stringsAsFactors = FALSE
    ),
    top_terms_per_entity = top_terms_per_module, max_terms = max_terms
  )
  if (!nrow(selection)) return(list(matrix = data.frame(), selection = selection, status = "no_selectable_terms"))
  selected_keys <- paste(selection$TermID, selection$Description, sep = "\r")
  matrix_keys <- paste(full$matrix$TermID, full$matrix$Description, sep = "\r")
  matrix <- full$matrix[matrix_keys %in% selected_keys, , drop = FALSE]
  matrix$selection_rank <- selection$selection_rank[match(paste(matrix$TermID, matrix$Description, sep = "\r"), selected_keys)]
  matrix <- matrix[order(matrix$selection_rank, matrix$ModuleID), , drop = FALSE]
  rownames(matrix) <- NULL
  list(matrix = matrix, selection = selection, status = "ok")
}

make_supermodule_go_heatmap_data <- function(
    module_heatmap,
    module_supermodule_map,
    ontology = "BP",
    fdr_cutoff = 0.05,
    score_cap = 6,
    top_terms_per_supermodule = 3L,
    max_terms = 30L
) {
  required_map <- c("ModuleID", "SupermoduleID")
  missing_map <- setdiff(required_map, names(module_supermodule_map))
  if (length(missing_map)) stop("Module-to-supermodule map is missing: ", paste(missing_map, collapse = ", "), call. = FALSE)
  if (is.null(module_heatmap) || !nrow(module_heatmap)) return(list(matrix = data.frame(), selection = data.frame(), status = "no_module_heatmap_terms"))
  map <- module_supermodule_map[, unique(c(required_map, intersect("SupermoduleDisplayLabel", names(module_supermodule_map)))), drop = FALSE]
  map$ModuleID <- as.character(map$ModuleID)
  map$SupermoduleID <- as.character(map$SupermoduleID)
  if (anyDuplicated(map$ModuleID)) stop("Module-to-supermodule map has duplicated ModuleID values.", call. = FALSE)
  x <- merge(module_heatmap, map, by = "ModuleID", all.x = FALSE, sort = FALSE)
  if (!nrow(x)) return(list(matrix = data.frame(), selection = data.frame(), status = "no_modules_with_supermodule_membership"))
  rows <- lapply(split(x, paste(x$SupermoduleID, x$TermID, sep = "\r")), function(one) {
    n_members <- length(unique(one$ModuleID)); n_sig <- sum(one$Significant, na.rm = TRUE)
    data.frame(
      SupermoduleID = one$SupermoduleID[[1]],
      SupermoduleDisplayLabel = if ("SupermoduleDisplayLabel" %in% names(one)) as.character(one$SupermoduleDisplayLabel[[1]]) else one$SupermoduleID[[1]],
      TermID = one$TermID[[1]], Description = one$Description[[1]],
      n_member_modules = n_members, n_modules_FDR_significant = n_sig,
      fraction_member_modules_FDR_significant = n_sig / n_members,
      mean_member_module_enrichment_score = mean(one$EnrichmentScore, na.rm = TRUE),
      max_member_module_enrichment_score = max(one$EnrichmentScore, na.rm = TRUE),
      min_member_module_FDR = if (any(one$Significant, na.rm = TRUE)) min(one$p_adjust[one$Significant], na.rm = TRUE) else 1,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$EnrichmentScore <- out$mean_member_module_enrichment_score
  out$Significant <- out$n_modules_FDR_significant > 0
  out$Ontology <- ontology
  out$FDRCutoff <- fdr_cutoff
  out$ScoreCap <- score_cap
  out$enrichment_measure <- "mean capped -log10(BH FDR) across member modules; zero = no member FDR-significant"
  selection <- select_representative_go_terms(
    data.frame(EntityID = out$SupermoduleID, TermID = out$TermID, Description = out$Description,
               EnrichmentScore = out$EnrichmentScore, Significant = out$Significant, stringsAsFactors = FALSE),
    top_terms_per_entity = top_terms_per_supermodule, max_terms = max_terms
  )
  if (!nrow(selection)) return(list(matrix = data.frame(), selection = selection, status = "no_selectable_terms"))
  selected_keys <- paste(selection$TermID, selection$Description, sep = "\r")
  keys <- paste(out$TermID, out$Description, sep = "\r")
  out <- out[keys %in% selected_keys, , drop = FALSE]
  out$selection_rank <- selection$selection_rank[match(paste(out$TermID, out$Description, sep = "\r"), selected_keys)]
  out <- out[order(out$selection_rank, out$SupermoduleID), , drop = FALSE]
  rownames(out) <- NULL
  list(matrix = out, selection = selection, status = "ok")
}

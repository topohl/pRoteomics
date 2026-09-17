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
    score_cap = 6,
    retain_gene_membership = FALSE
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
  if (isTRUE(retain_gene_membership) && "geneID" %in% names(go)) retained <- c(retained, "geneID")
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

wgcna_canonical_supermodule_order <- function(supermodule_ids) {
  ids <- unique(as.character(supermodule_ids))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  numeric_id <- suppressWarnings(as.integer(sub("^SM", "", ids)))
  ids[order(is.na(numeric_id), numeric_id, ids)]
}

summarize_supermodule_go_evidence <- function(
    module_go,
    module_supermodule_map,
    ontology = "BP",
    fdr_cutoff = 0.05,
    score_cap = 6
) {
  required_map <- c("ModuleID", "SupermoduleID")
  missing_map <- setdiff(required_map, names(module_supermodule_map))
  if (length(missing_map)) stop("Module-to-supermodule map is missing: ", paste(missing_map, collapse = ", "), call. = FALSE)
  if (is.null(module_go) || !nrow(module_go)) return(data.frame())
  required_go <- c("ModuleID", "TermID", "Description", "p_adjust", "EnrichmentScore", "Significant")
  missing_go <- setdiff(required_go, names(module_go))
  if (length(missing_go)) stop("Module GO matrix is missing: ", paste(missing_go, collapse = ", "), call. = FALSE)

  map <- module_supermodule_map[, unique(c(required_map, intersect("SupermoduleDisplayLabel", names(module_supermodule_map)))), drop = FALSE]
  map$ModuleID <- as.character(map$ModuleID)
  map$SupermoduleID <- as.character(map$SupermoduleID)
  if (anyDuplicated(map$ModuleID)) stop("Module-to-supermodule map has duplicated ModuleID values.", call. = FALSE)
  x <- merge(module_go, map, by = "ModuleID", all.x = FALSE, sort = FALSE)
  if (!nrow(x)) return(data.frame())

  term_descriptions <- split(as.character(x$Description), as.character(x$TermID))
  conflicting_terms <- names(term_descriptions)[vapply(term_descriptions, function(z) length(unique(z[!is.na(z)])) > 1L, logical(1))]
  if (length(conflicting_terms)) {
    stop("GO TermID values have conflicting Description values: ", paste(utils::head(conflicting_terms, 5L), collapse = ", "), call. = FALSE)
  }

  rows <- lapply(split(x, paste(x$SupermoduleID, x$TermID, sep = "\r")), function(one) {
    n_members <- length(unique(one$ModuleID))
    n_sig <- sum(one$Significant, na.rm = TRUE)
    gene_membership_available <- "geneID" %in% names(one)
    supporting_gene_ids <- character(0)
    if (gene_membership_available && n_sig > 0L) {
      supporting_rows <- as.character(one$geneID[one$Significant %in% TRUE])
      supporting_rows <- supporting_rows[!is.na(supporting_rows) & nzchar(supporting_rows)]
      if (length(supporting_rows)) {
        supporting_gene_ids <- unique(trimws(unlist(strsplit(supporting_rows, "/", fixed = TRUE))))
        supporting_gene_ids <- sort(supporting_gene_ids[!is.na(supporting_gene_ids) & nzchar(supporting_gene_ids)])
      }
    }
    data.frame(
      SupermoduleID = one$SupermoduleID[[1]],
      SupermoduleDisplayLabel = if ("SupermoduleDisplayLabel" %in% names(one)) as.character(one$SupermoduleDisplayLabel[[1]]) else one$SupermoduleID[[1]],
      TermID = one$TermID[[1]], Description = one$Description[[1]],
      n_member_modules = n_members, n_modules_FDR_significant = n_sig,
      fraction_member_modules_FDR_significant = n_sig / n_members,
      mean_member_module_enrichment_score = mean(one$EnrichmentScore, na.rm = TRUE),
      max_member_module_enrichment_score = max(one$EnrichmentScore, na.rm = TRUE),
      min_member_module_FDR = if (any(one$Significant, na.rm = TRUE)) min(one$p_adjust[one$Significant], na.rm = TRUE) else 1,
      is_singleton_supermodule = n_members == 1L,
      n_supporting_genes = if (gene_membership_available) length(supporting_gene_ids) else NA_integer_,
      supporting_gene_ids = if (gene_membership_available) paste(supporting_gene_ids, collapse = "/") else NA_character_,
      supporting_gene_id_delimiter = if (gene_membership_available) "/" else NA_character_,
      supporting_gene_membership_available = gene_membership_available && n_sig > 0L && length(supporting_gene_ids) > 0L,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$EnrichmentScore <- out$mean_member_module_enrichment_score
  out$Significant <- out$n_modules_FDR_significant > 0L
  out$Ontology <- ontology
  out$FDRCutoff <- fdr_cutoff
  out$ScoreCap <- score_cap
  out$enrichment_measure <- "mean capped -log10(BH FDR) across member modules; zero = no member FDR-significant"
  out$supermodule_summary_scope <- ifelse(
    out$is_singleton_supermodule,
    "member-module evidence; singleton reflects constituent module; not pooled supermodule ORA",
    "member-module evidence; not pooled supermodule ORA"
  )
  supermodule_order <- wgcna_canonical_supermodule_order(map$SupermoduleID)
  out <- out[order(match(out$SupermoduleID, supermodule_order), out$TermID), , drop = FALSE]
  rownames(out) <- NULL
  out
}

wgcna_go_hierarchy <- function(term_ids, ontology = "BP") {
  ontology <- toupper(as.character(ontology[[1]]))
  map_name <- switch(
    ontology,
    BP = "GOBPANCESTOR",
    MF = "GOMFANCESTOR",
    CC = "GOCCANCESTOR",
    NA_character_
  )
  unavailable <- function(reason) {
    list(
      available = FALSE, ancestors = list(), ontology = ontology,
      source = NA_character_, version = NA_character_, reason = reason
    )
  }
  if (is.na(map_name)) return(unavailable("unsupported_ontology"))
  if (!requireNamespace("GO.db", quietly = TRUE) || !requireNamespace("AnnotationDbi", quietly = TRUE)) {
    return(unavailable("GO.db_or_AnnotationDbi_not_installed"))
  }
  term_ids <- unique(as.character(term_ids))
  term_ids <- term_ids[!is.na(term_ids) & nzchar(term_ids)]
  if (!length(term_ids)) return(unavailable("no_GO_term_ids"))
  ancestors <- tryCatch(
    suppressWarnings(as.list(getExportedValue("GO.db", map_name)[term_ids])),
    error = function(e) NULL
  )
  if (is.null(ancestors)) return(unavailable("GO_ancestor_lookup_failed"))
  ancestors <- lapply(ancestors, function(one) {
    one <- as.character(one)
    sort(unique(one[grepl("^GO:[0-9]+$", one)]))
  })
  names(ancestors) <- term_ids
  list(
    available = TRUE, ancestors = ancestors, ontology = ontology,
    source = paste0("GO.db::", map_name),
    version = as.character(utils::packageVersion("GO.db")), reason = NA_character_
  )
}

wgcna_go_pair_relationship <- function(candidate_term_id, selected_term_id, go_hierarchy) {
  if (identical(candidate_term_id, selected_term_id)) return("identical_GO_ID")
  if (is.null(go_hierarchy) || !isTRUE(go_hierarchy$available)) return(NA_character_)
  ancestors <- go_hierarchy$ancestors
  candidate_ancestors <- ancestors[[candidate_term_id]]
  selected_ancestors <- ancestors[[selected_term_id]]
  if (!is.null(candidate_ancestors) && selected_term_id %in% candidate_ancestors) {
    return("candidate_descendant_of_selected")
  }
  if (!is.null(selected_ancestors) && candidate_term_id %in% selected_ancestors) {
    return("candidate_ancestor_of_selected")
  }
  "none"
}

wgcna_gene_jaccard <- function(gene_ids_a, gene_ids_b, delimiter = "/") {
  parse_one <- function(x) {
    if (length(x) != 1L || is.na(x) || !nzchar(x)) return(character(0))
    out <- trimws(strsplit(as.character(x), delimiter, fixed = TRUE)[[1]])
    sort(unique(out[!is.na(out) & nzchar(out)]))
  }
  a <- parse_one(gene_ids_a)
  b <- parse_one(gene_ids_b)
  union_size <- length(union(a, b))
  if (!length(a) || !length(b) || union_size == 0L) return(NA_real_)
  length(intersect(a, b)) / union_size
}

make_focused_supermodule_go_comparison <- function(
    module_go,
    module_supermodule_map,
    ontology = "BP",
    fdr_cutoff = 0.05,
    score_cap = 6,
    focused_terms_per_supermodule = 3L,
    go_hierarchy = NULL,
    hierarchy_gene_jaccard_threshold = 0.50,
    near_identical_gene_jaccard_threshold = 0.80
) {
  focused_terms_per_supermodule <- suppressWarnings(as.integer(focused_terms_per_supermodule))
  if (length(focused_terms_per_supermodule) != 1L || !is.finite(focused_terms_per_supermodule) || focused_terms_per_supermodule < 1L) {
    stop("focused_terms_per_supermodule must be a positive integer.", call. = FALSE)
  }
  summary <- summarize_supermodule_go_evidence(
    module_go, module_supermodule_map, ontology = ontology,
    fdr_cutoff = fdr_cutoff, score_cap = score_cap
  )
  empty_result <- function(status, redundancy_mode = NA_character_) {
    list(matrix = data.frame(), selection = data.frame(), audit = data.frame(), status = status,
         redundancy_mode = redundancy_mode)
  }
  if (!nrow(summary)) return(empty_result("no_supermodule_GO_terms"))

  supermodule_order <- wgcna_canonical_supermodule_order(module_supermodule_map$SupermoduleID)
  supported <- summary[summary$n_modules_FDR_significant > 0L, , drop = FALSE]
  hierarchy_available <- !is.null(go_hierarchy) && isTRUE(go_hierarchy$available)
  gene_membership_available <- nrow(supported) > 0L &&
    all(supported$supporting_gene_membership_available %in% TRUE) &&
    all(is.finite(supported$n_supporting_genes) & supported$n_supporting_genes > 0L)
  redundancy_mode <- if (hierarchy_available && gene_membership_available) {
    "GO_hierarchy_plus_gene_overlap"
  } else if (hierarchy_available) {
    "GO_hierarchy_only"
  } else {
    "evidence_rank_only"
  }
  if (!nrow(supported)) return(empty_result("no_FDR_supported_terms", redundancy_mode))
  supported <- supported[order(
    match(supported$SupermoduleID, supermodule_order),
    -supported$fraction_member_modules_FDR_significant,
    -supported$mean_member_module_enrichment_score,
    supported$min_member_module_FDR,
    supported$TermID,
    supported$Description
  ), , drop = FALSE]
  supported$initial_evidence_rank <- ave(
    seq_len(nrow(supported)), supported$SupermoduleID,
    FUN = function(z) seq_along(z)
  )
  supported$selected <- FALSE
  supported$focused_representative_rank <- NA_integer_
  supported$selection_status <- NA_character_
  supported$redundant_with_TermID <- NA_character_
  supported$redundant_with_Description <- NA_character_
  supported$go_relationship <- NA_character_
  supported$gene_jaccard <- NA_real_
  supported$redundancy_reason <- NA_character_
  supported$redundancy_mode <- redundancy_mode

  for (supermodule_id in supermodule_order) {
    candidate_indices <- which(supported$SupermoduleID == supermodule_id)
    selected_indices <- integer(0)
    for (candidate_index in candidate_indices) {
      if (length(selected_indices) >= focused_terms_per_supermodule) {
        supported$selection_status[[candidate_index]] <- "not_reached_after_three_selected"
        next
      }
      redundant <- FALSE
      if (length(selected_indices) && !identical(redundancy_mode, "evidence_rank_only")) {
        for (selected_index in selected_indices) {
          relationship <- wgcna_go_pair_relationship(
            supported$TermID[[candidate_index]], supported$TermID[[selected_index]], go_hierarchy
          )
          hierarchical_pair <- !is.na(relationship) && relationship %in% c(
            "identical_GO_ID", "candidate_descendant_of_selected", "candidate_ancestor_of_selected"
          )
          gene_jaccard <- if (identical(redundancy_mode, "GO_hierarchy_plus_gene_overlap")) {
            wgcna_gene_jaccard(
              supported$supporting_gene_ids[[candidate_index]],
              supported$supporting_gene_ids[[selected_index]]
            )
          } else {
            NA_real_
          }
          hierarchy_redundant <- if (identical(redundancy_mode, "GO_hierarchy_only")) {
            hierarchical_pair
          } else {
            hierarchical_pair && is.finite(gene_jaccard) && gene_jaccard >= hierarchy_gene_jaccard_threshold
          }
          near_identical_redundant <- identical(redundancy_mode, "GO_hierarchy_plus_gene_overlap") &&
            is.finite(gene_jaccard) && gene_jaccard >= near_identical_gene_jaccard_threshold
          if (hierarchy_redundant || near_identical_redundant) {
            redundant <- TRUE
            supported$selection_status[[candidate_index]] <- if (hierarchy_redundant) {
              "skipped_redundant_GO_ancestor_descendant"
            } else {
              "skipped_near_identical_gene_set"
            }
            supported$redundant_with_TermID[[candidate_index]] <- supported$TermID[[selected_index]]
            supported$redundant_with_Description[[candidate_index]] <- supported$Description[[selected_index]]
            supported$go_relationship[[candidate_index]] <- relationship
            supported$gene_jaccard[[candidate_index]] <- gene_jaccard
            supported$redundancy_reason[[candidate_index]] <- if (hierarchy_redundant) {
              if (identical(redundancy_mode, "GO_hierarchy_only")) {
                paste0("exact_GO_ancestor_descendant; gene_membership_unavailable; hierarchy_only_fallback")
              } else {
                paste0("exact_GO_ancestor_descendant_and_gene_jaccard>=", format(hierarchy_gene_jaccard_threshold, nsmall = 2L))
              }
            } else {
              paste0("gene_jaccard>=", format(near_identical_gene_jaccard_threshold, nsmall = 2L), "; hierarchy_not_required")
            }
            break
          }
        }
      }
      if (!redundant) {
        selected_indices <- c(selected_indices, candidate_index)
        supported$selected[[candidate_index]] <- TRUE
        supported$focused_representative_rank[[candidate_index]] <- length(selected_indices)
        supported$selection_status[[candidate_index]] <- "selected"
      }
    }
  }

  audit <- supported
  selection <- audit[audit$selected, , drop = FALSE]
  selection$selected_for_supermodule <- selection$SupermoduleID
  selection$representative_rank <- selection$focused_representative_rank
  selection$selection_basis <- paste0(
    "FDR_supported_member_module_evidence_ranked_by_fraction_desc_mean_score_desc_min_FDR_asc_TermID_asc; max_",
    focused_terms_per_supermodule, "_per_supermodule; sequential_conservative_redundancy_pruning; mode_",
    redundancy_mode, "; no_nonsignificant_fallback"
  )

  provenance_rows <- lapply(split(selection, selection$TermID), function(one) {
    one <- one[order(match(one$selected_for_supermodule, supermodule_order), one$focused_representative_rank), , drop = FALSE]
    data.frame(
      TermID = one$TermID[[1]],
      selecting_supermodules = paste(one$selected_for_supermodule, collapse = ";"),
      n_selecting_supermodules = nrow(one),
      selector_provenance = paste0(one$selected_for_supermodule, ":rank", one$focused_representative_rank, collapse = ";"),
      row_group_supermodule = one$selected_for_supermodule[[1]],
      row_group_representative_rank = one$focused_representative_rank[[1]],
      stringsAsFactors = FALSE
    )
  })
  provenance <- do.call(rbind, provenance_rows)
  provenance <- provenance[order(
    match(provenance$row_group_supermodule, supermodule_order),
    provenance$row_group_representative_rank,
    provenance$TermID
  ), , drop = FALSE]
  provenance$row_order <- seq_len(nrow(provenance))

  matrix <- summary[summary$TermID %in% provenance$TermID, , drop = FALSE]
  matrix <- merge(matrix, provenance, by = "TermID", all.x = TRUE, sort = FALSE)
  matrix$SupermodulePlotLabel <- ifelse(matrix$is_singleton_supermodule, paste0(matrix$SupermoduleID, "\u2020"), matrix$SupermoduleID)
  matrix$display_supported_dot <- matrix$n_modules_FDR_significant > 0L
  matrix$dot_colour_measure <- "mean_member_module_enrichment_score"
  matrix$dot_size_measure <- "fraction_member_modules_FDR_significant"
  matrix <- matrix[order(matrix$row_order, match(matrix$SupermoduleID, supermodule_order)), , drop = FALSE]
  selection$row_order <- provenance$row_order[match(selection$TermID, provenance$TermID)]
  selection <- selection[order(selection$row_order, match(selection$selected_for_supermodule, supermodule_order)), , drop = FALSE]
  rownames(matrix) <- NULL
  rownames(selection) <- NULL
  rownames(audit) <- NULL
  list(
    matrix = matrix, selection = selection, audit = audit, status = "ok",
    redundancy_mode = redundancy_mode,
    hierarchy_gene_jaccard_threshold = hierarchy_gene_jaccard_threshold,
    near_identical_gene_jaccard_threshold = near_identical_gene_jaccard_threshold
  )
}

wgcna_focused_dataset_order <- function() {
  c("neuron_neuropil", "neuron_soma", "microglia")
}

wgcna_focused_dataset_display_labels <- function() {
  c(
    neuron_neuropil = "Neuropil",
    neuron_soma = "Soma",
    microglia = "Microglia-enriched ROI"
  )
}

combine_focused_supermodule_go_sources <- function(
    source_tables,
    selection_tables,
    dataset_order = wgcna_focused_dataset_order()
) {
  dataset_order <- unique(as.character(dataset_order))
  if (!length(dataset_order) || any(is.na(dataset_order) | !nzchar(dataset_order))) {
    stop("dataset_order must contain non-empty dataset identifiers.", call. = FALSE)
  }
  if (!is.list(source_tables) || !is.list(selection_tables) ||
      is.null(names(source_tables)) || is.null(names(selection_tables))) {
    stop("Combined focused GO inputs must be named source and selection lists.", call. = FALSE)
  }
  missing_source <- setdiff(dataset_order, names(source_tables))
  missing_selection <- setdiff(dataset_order, names(selection_tables))
  if (length(missing_source) || length(missing_selection)) {
    stop(
      "Combined focused GO inputs are missing dataset(s): ",
      paste(unique(c(missing_source, missing_selection)), collapse = ", "),
      call. = FALSE
    )
  }

  source_required <- c(
    "dataset", "Ontology", "TermID", "Description", "SupermoduleID",
    "SupermodulePlotLabel", "row_order", "n_modules_FDR_significant",
    "fraction_member_modules_FDR_significant", "mean_member_module_enrichment_score",
    "is_singleton_supermodule", "supermodule_summary_scope", "FDRCutoff", "ScoreCap",
    "display_supported_dot", "dot_colour_measure", "dot_size_measure"
  )
  selection_required <- c(
    "dataset", "Ontology", "TermID", "Description", "SupermoduleID",
    "selected_for_supermodule", "representative_rank", "selection_basis"
  )
  validate_one <- function(x, dataset, required, artifact) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    missing <- setdiff(required, names(x))
    if (length(missing)) stop(artifact, " for ", dataset, " is missing: ", paste(missing, collapse = ", "), call. = FALSE)
    observed <- unique(as.character(x$dataset))
    observed <- observed[!is.na(observed) & nzchar(observed)]
    if (!identical(observed, dataset)) stop(artifact, " dataset column does not match ", dataset, ".", call. = FALSE)
    x
  }
  source_tables <- lapply(dataset_order, function(dataset) {
    validate_one(source_tables[[dataset]], dataset, source_required, "Focused source")
  })
  names(source_tables) <- dataset_order
  selection_tables <- lapply(dataset_order, function(dataset) {
    validate_one(selection_tables[[dataset]], dataset, selection_required, "Focused selection")
  })
  names(selection_tables) <- dataset_order
  source <- do.call(rbind, source_tables)
  selection <- do.call(rbind, selection_tables)
  rownames(source) <- NULL
  rownames(selection) <- NULL

  if (!nrow(source) || !nrow(selection)) stop("Combined focused GO inputs contain no selected terms.", call. = FALSE)
  if (!all(source$dot_colour_measure == "mean_member_module_enrichment_score")) {
    stop("Combined focused GO colour must use mean_member_module_enrichment_score.", call. = FALSE)
  }
  if (!all(source$dot_size_measure == "fraction_member_modules_FDR_significant")) {
    stop("Combined focused GO size must use fraction_member_modules_FDR_significant.", call. = FALSE)
  }
  support_fraction <- suppressWarnings(as.numeric(source$fraction_member_modules_FDR_significant))
  if (any(!is.finite(support_fraction) | support_fraction < 0 | support_fraction > 1)) {
    stop("Combined focused GO support fractions must remain on the unnormalized 0-1 scale.", call. = FALSE)
  }
  score_caps <- unique(suppressWarnings(as.numeric(source$ScoreCap)))
  if (length(score_caps) != 1L || !is.finite(score_caps) || score_caps <= 0) {
    stop("Combined focused GO sources must share one finite positive ScoreCap.", call. = FALSE)
  }
  fdr_cutoffs <- unique(suppressWarnings(as.numeric(source$FDRCutoff)))
  if (length(fdr_cutoffs) != 1L || !is.finite(fdr_cutoffs) || fdr_cutoffs <= 0 || fdr_cutoffs >= 1) {
    stop("Combined focused GO sources must share one valid FDRCutoff.", call. = FALSE)
  }
  ontologies <- unique(as.character(source$Ontology))
  if (length(ontologies) != 1L) stop("Combined focused GO sources must contain one ontology.", call. = FALSE)

  term_descriptions <- split(
    as.character(c(source$Description, selection$Description)),
    as.character(c(source$TermID, selection$TermID))
  )
  conflicting_terms <- names(term_descriptions)[vapply(
    term_descriptions,
    function(z) length(unique(z[!is.na(z) & nzchar(z)])) > 1L,
    logical(1)
  )]
  if (length(conflicting_terms)) {
    stop("Exact GO TermID values have conflicting canonical descriptions: ", paste(utils::head(conflicting_terms, 5L), collapse = ", "), call. = FALSE)
  }

  selected_presence <- unique(selection[c("dataset", "TermID")])
  recurrence_rows <- lapply(split(selected_presence, selected_presence$TermID), function(one) {
    datasets <- dataset_order[dataset_order %in% unique(as.character(one$dataset))]
    data.frame(
      TermID = one$TermID[[1]],
      n_datasets_selected = length(datasets),
      datasets_selected = paste(datasets, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  recurrence <- do.call(rbind, recurrence_rows)
  rownames(recurrence) <- NULL

  attach_recurrence <- function(x) {
    x$.combined_input_order <- seq_len(nrow(x))
    x <- merge(x, recurrence, by = "TermID", all.x = TRUE, sort = FALSE)
    x <- x[order(x$.combined_input_order), , drop = FALSE]
    x$.combined_input_order <- NULL
    if (any(is.na(x$n_datasets_selected))) stop("Combined focused GO provenance is incomplete.", call. = FALSE)
    display_labels <- wgcna_focused_dataset_display_labels()
    x$dataset_display_label <- unname(display_labels[match(x$dataset, names(display_labels))])
    x$panel_label <- letters[match(x$dataset, dataset_order)]
    x$cross_dataset_recurrence_exact_TermID <- x$n_datasets_selected >= 2L
    x$cross_dataset_recurrence_marker <- ifelse(x$cross_dataset_recurrence_exact_TermID, "\u25c6", "")
    x$cross_dataset_recurrence_scope <- "descriptive exact TermID selection only; no cross-dataset inference"
    x$combined_colour_measure <- "mean_member_module_enrichment_score"
    x$combined_colour_scale_min <- 0
    x$combined_colour_scale_max <- score_caps[[1]]
    x$combined_colour_normalization <- "none; shared unnormalized member-module summary scale"
    x$combined_size_measure <- "fraction_member_modules_FDR_significant"
    x$combined_size_scale_min <- 0
    x$combined_size_scale_max <- 1
    x$combined_size_normalization <- "none; shared member-module support fraction"
    x$combined_selection_FDR_cutoff <- fdr_cutoffs[[1]]
    rownames(x) <- NULL
    x
  }

  list(
    matrix = attach_recurrence(source),
    selection = attach_recurrence(selection),
    colour_limits = c(0, score_caps[[1]]),
    size_limits = c(0, 1),
    fdr_cutoff = fdr_cutoffs[[1]],
    ontology = ontologies[[1]],
    status = "ok"
  )
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
  if (is.null(module_heatmap) || !nrow(module_heatmap)) return(list(matrix = data.frame(), selection = data.frame(), status = "no_module_heatmap_terms"))
  out <- summarize_supermodule_go_evidence(
    module_heatmap, module_supermodule_map, ontology = ontology,
    fdr_cutoff = fdr_cutoff, score_cap = score_cap
  )
  if (!nrow(out)) return(list(matrix = data.frame(), selection = data.frame(), status = "no_modules_with_supermodule_membership"))
  # Keep the established broad supermodule heatmap source schema unchanged;
  # focused-only singleton/scope/gene annotations remain in companion outputs.
  out <- out[, setdiff(names(out), c(
    "is_singleton_supermodule", "supermodule_summary_scope", "n_supporting_genes",
    "supporting_gene_ids", "supporting_gene_id_delimiter",
    "supporting_gene_membership_available"
  )), drop = FALSE]
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

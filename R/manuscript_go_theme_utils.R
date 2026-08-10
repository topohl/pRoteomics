# Ontology-aware manuscript GO-theme mapping and descriptive redundancy QA.
# The broad text-regex program mapper in R/enrichment_io.R remains a legacy
# technical classifier. Valid GO IDs are never publication-classified by text.

manuscript_go_allowed_relationships <- function() c("is_a", "part_of")

normalize_go_relationship <- function(x) {
  key <- tolower(trimws(as.character(x)))
  key <- gsub("[[:space:]-]+", "_", key)
  key[key == "isa"] <- "is_a"
  key
}

require_go_ontology_dependencies <- function(include_semantic = FALSE) {
  required <- c("GO.db", "AnnotationDbi")
  if (isTRUE(include_semantic)) required <- c(required, "GOSemSim", "org.Mm.eg.db")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing required installed ontology package(s): ", paste(missing, collapse = ", "),
      ". Install packages explicitly; this manuscript mapper does not install dependencies.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

go_ontology_provenance <- function() {
  require_go_ontology_dependencies(FALSE)
  info <- GO.db::GO_dbInfo()
  values <- stats::setNames(as.character(info$value), as.character(info$name))
  list(
    go_db_package_version = as.character(utils::packageVersion("GO.db")),
    go_source_name = unname(values[["GOSOURCENAME"]]),
    go_source_url = unname(values[["GOSOURCEURL"]]),
    go_source_date = unname(values[["GOSOURCEDATE"]]),
    approved_relationships = paste(manuscript_go_allowed_relationships(), collapse = ";")
  )
}

read_manuscript_go_theme_registry <- function(path) {
  require_go_ontology_dependencies(FALSE)
  if (!file.exists(path)) stop("Missing manuscript GO-theme registry: ", path, call. = FALSE)
  registry <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE, quote = "")
  required <- c(
    "theme_id", "display_label", "anchor_go_id", "anchor_label", "theme_role",
    "display_order", "match_scope", "rationale", "registry_version"
  )
  missing <- setdiff(required, names(registry))
  if (length(missing)) stop("Manuscript GO-theme registry is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  registry <- registry[, required, drop = FALSE]
  character_columns <- setdiff(required, "display_order")
  registry[character_columns] <- lapply(registry[character_columns], function(x) trimws(as.character(x)))
  registry$display_order <- suppressWarnings(as.integer(registry$display_order))
  if (!nrow(registry) || any(!nzchar(registry$theme_id)) || any(!nzchar(registry$display_label))) {
    stop("Manuscript GO-theme registry must contain non-empty theme IDs and display labels.", call. = FALSE)
  }
  if (anyDuplicated(paste(registry$theme_id, registry$anchor_go_id, registry$match_scope, sep = "|"))) {
    stop("Manuscript GO-theme registry contains duplicate theme + anchor + scope rows.", call. = FALSE)
  }
  if (any(!grepl("^GO:[0-9]{7}$", registry$anchor_go_id))) {
    stop("Manuscript GO-theme registry contains invalid GO IDs.", call. = FALSE)
  }
  allowed_roles <- c("primary", "supporting", "qc_review")
  if (any(!registry$theme_role %in% allowed_roles)) {
    stop("Manuscript GO-theme registry theme_role must be primary, supporting, or qc_review.", call. = FALSE)
  }
  allowed_scopes <- c("anchor_and_descendants", "exact_go_id")
  if (any(!registry$match_scope %in% allowed_scopes)) {
    stop("Manuscript GO-theme registry match_scope must be anchor_and_descendants or exact_go_id.", call. = FALSE)
  }
  if (anyNA(registry$display_order) || any(registry$display_order < 1L)) {
    stop("Manuscript GO-theme registry display_order values must be positive integers.", call. = FALSE)
  }
  versions <- unique(registry$registry_version[nzchar(registry$registry_version)])
  if (length(versions) != 1L) stop("Manuscript GO-theme registry must declare exactly one registry_version.", call. = FALSE)
  theme_meta <- unique(registry[c("theme_id", "display_label", "theme_role", "display_order")])
  if (anyDuplicated(theme_meta$theme_id)) {
    stop("Each manuscript theme must have one display label, role, and display order.", call. = FALSE)
  }
  for (i in seq_len(nrow(registry))) {
    term <- GO.db::GOTERM[[registry$anchor_go_id[[i]]]]
    if (is.null(term)) stop("Registry anchor is absent from GO.db: ", registry$anchor_go_id[[i]], call. = FALSE)
    if (!identical(AnnotationDbi::Ontology(term), "BP")) {
      stop("Registry anchor is not a GO biological-process term: ", registry$anchor_go_id[[i]], call. = FALSE)
    }
    live_label <- AnnotationDbi::Term(term)
    if (!identical(tolower(trimws(live_label)), tolower(registry$anchor_label[[i]]))) {
      stop(
        "Registry anchor label differs from GO.db for ", registry$anchor_go_id[[i]],
        ": registry='", registry$anchor_label[[i]], "', GO.db='", live_label, "'.",
        call. = FALSE
      )
    }
  }
  registry
}

go_bp_allowed_ancestry <- function(go_id,
                                   allowed_relationships = manuscript_go_allowed_relationships()) {
  require_go_ontology_dependencies(FALSE)
  allowed_relationships <- unique(normalize_go_relationship(allowed_relationships))
  if (is.null(GO.db::GOTERM[[go_id]])) {
    return(data.frame(
      ancestor_GO_ID = character(), path_length = integer(),
      relationship_path = character(), GO_path = character(), stringsAsFactors = FALSE
    ))
  }
  queue <- list(list(node = go_id, relationships = character(), nodes = go_id))
  seen_depth <- stats::setNames(0L, go_id)
  result <- list(data.frame(
    ancestor_GO_ID = go_id, path_length = 0L, relationship_path = "anchor",
    GO_path = go_id, stringsAsFactors = FALSE
  ))
  while (length(queue)) {
    current <- queue[[1L]]
    queue <- queue[-1L]
    parents <- GO.db::GOBPPARENTS[[current$node]]
    if (is.null(parents) || !length(parents)) next
    relations <- normalize_go_relationship(names(parents))
    keep <- relations %in% allowed_relationships
    if (!any(keep)) next
    edge_df <- data.frame(
      parent = unname(parents[keep]),
      relationship = relations[keep],
      stringsAsFactors = FALSE
    )
    edge_df <- edge_df[order(edge_df$parent, edge_df$relationship, method = "radix"), , drop = FALSE]
    for (i in seq_len(nrow(edge_df))) {
      parent <- edge_df$parent[[i]]
      next_relationships <- c(current$relationships, edge_df$relationship[[i]])
      next_nodes <- c(current$nodes, parent)
      next_depth <- length(next_relationships)
      prior_depth <- unname(seen_depth[parent])
      if (!length(prior_depth) || is.na(prior_depth) || next_depth < prior_depth) {
        seen_depth[[parent]] <- next_depth
        result[[length(result) + 1L]] <- data.frame(
          ancestor_GO_ID = parent,
          path_length = as.integer(next_depth),
          relationship_path = paste(next_relationships, collapse = ">"),
          GO_path = paste(next_nodes, collapse = ">"),
          stringsAsFactors = FALSE
        )
        queue[[length(queue) + 1L]] <- list(
          node = parent, relationships = next_relationships, nodes = next_nodes
        )
      }
    }
  }
  out <- do.call(rbind, result)
  out <- out[order(out$path_length, out$ancestor_GO_ID, out$relationship_path, method = "radix"), , drop = FALSE]
  out <- out[!duplicated(out$ancestor_GO_ID), , drop = FALSE]
  rownames(out) <- NULL
  out
}

go_bp_shortest_allowed_path <- function(go_id, anchor_go_id,
                                        allowed_relationships = manuscript_go_allowed_relationships()) {
  ancestry <- go_bp_allowed_ancestry(go_id, allowed_relationships)
  hit <- ancestry[ancestry$ancestor_GO_ID == anchor_go_id, , drop = FALSE]
  if (!nrow(hit)) {
    return(list(matched = FALSE, path_length = NA_integer_, relationship_path = NA_character_, go_path = NA_character_))
  }
  list(
    matched = TRUE,
    path_length = hit$path_length[[1]],
    relationship_path = hit$relationship_path[[1]],
    go_path = hit$GO_path[[1]]
  )
}

map_go_terms_to_manuscript_themes <- function(go_terms, registry,
                                              id_col = "ID", description_col = "Description") {
  require_go_ontology_dependencies(FALSE)
  if (!is.data.frame(go_terms)) stop("go_terms must be a data frame.", call. = FALSE)
  missing <- setdiff(c(id_col, description_col), names(go_terms))
  if (length(missing)) stop("GO term input is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  if (is.character(registry) && length(registry) == 1L) registry <- read_manuscript_go_theme_registry(registry)
  required_registry <- c(
    "theme_id", "display_label", "anchor_go_id", "anchor_label", "theme_role",
    "display_order", "match_scope", "rationale", "registry_version"
  )
  if (!is.data.frame(registry) || length(setdiff(required_registry, names(registry)))) {
    stop("registry must be a validated manuscript GO-theme registry.", call. = FALSE)
  }
  terms <- unique(data.frame(
    GO_ID = trimws(as.character(go_terms[[id_col]])),
    GO_description = as.character(go_terms[[description_col]]),
    stringsAsFactors = FALSE
  ))
  terms <- terms[order(terms$GO_ID, terms$GO_description, method = "radix"), , drop = FALSE]
  provenance <- go_ontology_provenance()
  valid_bp <- vapply(terms$GO_ID, function(go_id) {
    term_obj <- if (grepl("^GO:[0-9]{7}$", go_id)) GO.db::GOTERM[[go_id]] else NULL
    !is.null(term_obj) && identical(AnnotationDbi::Ontology(term_obj), "BP")
  }, logical(1))
  ancestry_cache <- stats::setNames(
    lapply(terms$GO_ID[valid_bp], go_bp_allowed_ancestry),
    terms$GO_ID[valid_bp]
  )
  assignments <- vector("list", 0L)
  k <- 0L
  for (i in seq_len(nrow(terms))) {
    go_id <- terms$GO_ID[[i]]
    if (!valid_bp[[i]]) next
    ancestry <- ancestry_cache[[go_id]]
    for (j in seq_len(nrow(registry))) {
      if (identical(registry$match_scope[[j]], "exact_go_id")) {
        path <- list(
          matched = identical(go_id, registry$anchor_go_id[[j]]), path_length = 0L,
          relationship_path = "anchor", go_path = go_id
        )
      } else {
        hit <- ancestry[ancestry$ancestor_GO_ID == registry$anchor_go_id[[j]], , drop = FALSE]
        path <- if (nrow(hit)) {
          list(
            matched = TRUE, path_length = hit$path_length[[1]],
            relationship_path = hit$relationship_path[[1]], go_path = hit$GO_path[[1]]
          )
        } else {
          list(matched = FALSE, path_length = NA_integer_, relationship_path = NA_character_, go_path = NA_character_)
        }
      }
      if (!isTRUE(path$matched)) next
      k <- k + 1L
      assignments[[k]] <- data.frame(
        GO_ID = go_id,
        GO_description = terms$GO_description[[i]],
        theme_id = registry$theme_id[[j]],
        manuscript_theme = registry$display_label[[j]],
        theme_role = registry$theme_role[[j]],
        display_order = registry$display_order[[j]],
        anchor_GO_ID = registry$anchor_go_id[[j]],
        anchor_label = registry$anchor_label[[j]],
        match_scope = registry$match_scope[[j]],
        match_type = if (path$path_length == 0L) "exact_anchor" else "approved_descendant",
        path_length = path$path_length,
        relationship_path = path$relationship_path,
        GO_path = path$go_path,
        relationship_types_approved = provenance$approved_relationships,
        mapping_method = "go_id_ontology",
        registry_version = registry$registry_version[[j]],
        GO_db_package_version = provenance$go_db_package_version,
        GO_source_name = provenance$go_source_name,
        GO_source_url = provenance$go_source_url,
        GO_source_date = provenance$go_source_date,
        stringsAsFactors = FALSE
      )
    }
  }
  assignment_columns <- c(
    "GO_ID", "GO_description", "theme_id", "manuscript_theme", "theme_role", "display_order",
    "anchor_GO_ID", "anchor_label", "match_scope", "match_type", "path_length",
    "relationship_path", "GO_path", "relationship_types_approved", "mapping_method",
    "registry_version", "GO_db_package_version", "GO_source_name", "GO_source_url", "GO_source_date"
  )
  assignment_df <- if (length(assignments)) do.call(rbind, assignments) else {
    out <- as.data.frame(stats::setNames(replicate(length(assignment_columns), character(), simplify = FALSE), assignment_columns), stringsAsFactors = FALSE)
    out$display_order <- integer()
    out$path_length <- integer()
    out
  }
  if (nrow(assignment_df)) {
    assignment_df <- assignment_df[order(
      assignment_df$GO_ID, assignment_df$display_order, assignment_df$theme_id,
      assignment_df$anchor_GO_ID, assignment_df$path_length, method = "radix"
    ), , drop = FALSE]
    rownames(assignment_df) <- NULL
  }
  status_rows <- lapply(seq_len(nrow(terms)), function(i) {
    hit <- assignment_df[assignment_df$GO_ID == terms$GO_ID[[i]], , drop = FALSE]
    themes <- unique(hit$theme_id)
    roles <- unique(hit$theme_role)
    status <- if (!length(themes)) {
      "unclassified"
    } else if (all(roles == "qc_review")) {
      "qc_review"
    } else if (length(themes) > 1L || ("qc_review" %in% roles && any(roles != "qc_review"))) {
      "multi_theme"
    } else {
      "single_theme"
    }
    data.frame(
      GO_ID = terms$GO_ID[[i]],
      GO_description = terms$GO_description[[i]],
      assignment_status = status,
      n_manuscript_themes = length(themes),
      manuscript_theme_ids = if (length(themes)) paste(sort(themes), collapse = ";") else NA_character_,
      manuscript_themes = if (nrow(hit)) paste(unique(hit$manuscript_theme[order(hit$display_order)]), collapse = ";") else NA_character_,
      theme_roles = if (length(roles)) paste(sort(roles), collapse = ";") else NA_character_,
      anchor_GO_IDs = if (nrow(hit)) paste(unique(hit$anchor_GO_ID), collapse = ";") else NA_character_,
      mapping_method = if (nrow(hit)) "go_id_ontology" else "go_id_unclassified",
      registry_version = unique(registry$registry_version)[[1]],
      GO_db_package_version = provenance$go_db_package_version,
      GO_source_name = provenance$go_source_name,
      GO_source_url = provenance$go_source_url,
      GO_source_date = provenance$go_source_date,
      relationship_types_approved = provenance$approved_relationships,
      stringsAsFactors = FALSE
    )
  })
  status_df <- do.call(rbind, status_rows)
  rownames(status_df) <- NULL
  list(assignments = assignment_df, term_status = status_df, registry = registry, provenance = provenance)
}

collapse_go_theme_assignment_audit <- function(mapping) {
  assignments <- mapping$assignments
  if (!nrow(assignments)) return(assignments)
  groups <- split(assignments, paste(assignments$GO_ID, assignments$theme_id, sep = "|"))
  out <- do.call(rbind, lapply(groups, function(x) {
    data.frame(
      GO_ID = x$GO_ID[[1]],
      GO_description = x$GO_description[[1]],
      theme_id = x$theme_id[[1]],
      manuscript_theme = x$manuscript_theme[[1]],
      theme_role = x$theme_role[[1]],
      display_order = x$display_order[[1]],
      anchor_GO_IDs = collapse_distinct_values(x$anchor_GO_ID),
      anchor_labels = collapse_distinct_values(x$anchor_label),
      match_scopes = collapse_distinct_values(x$match_scope),
      match_types = collapse_distinct_values(x$match_type),
      shortest_path_length = min(x$path_length),
      relationship_paths = collapse_distinct_values(x$relationship_path),
      GO_paths = collapse_distinct_values(x$GO_path),
      relationship_types_approved = x$relationship_types_approved[[1]],
      mapping_method = x$mapping_method[[1]],
      registry_version = x$registry_version[[1]],
      GO_db_package_version = x$GO_db_package_version[[1]],
      GO_source_name = x$GO_source_name[[1]],
      GO_source_url = x$GO_source_url[[1]],
      GO_source_date = x$GO_source_date[[1]],
      stringsAsFactors = FALSE
    )
  }))
  out <- out[order(out$GO_ID, out$display_order, out$theme_id, method = "radix"), , drop = FALSE]
  rownames(out) <- NULL
  out
}

collapse_distinct_values <- function(x, separator = ";") {
  x <- unique(trimws(as.character(x)))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x)) paste(x, collapse = separator) else NA_character_
}

build_sus_res_supported_go_theme_audit <- function(enrichment_df, mapping, fdr_threshold = 0.05) {
  required <- c(
    "dataset", "compartment", "region", "layer", "spatial_unit", "phenotype_contrast",
    "Comparison", "ID", "Description", "NES", "p.adjust", "core_enrichment",
    "core_enrichment_gene", "program_class", "source_supplementary_file",
    "source_manifest_file", "evidence_source_family"
  )
  missing <- setdiff(required, names(enrichment_df))
  if (length(missing)) stop("Supported GO audit input is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  supported <- enrichment_df[
    as.character(enrichment_df$phenotype_contrast) == "SUS_vs_RES" &
      !is.na(enrichment_df$p.adjust) & enrichment_df$p.adjust < fdr_threshold,
    , drop = FALSE
  ]
  supported$supported_occurrence_id <- sprintf("SUS_RES_GO_%04d", seq_len(nrow(supported)))
  status <- mapping$term_status
  idx <- match(as.character(supported$ID), status$GO_ID)
  if (anyNA(idx)) stop("Every supported GO term must have an ontology assignment status.", call. = FALSE)
  out <- data.frame(
    supported_occurrence_id = supported$supported_occurrence_id,
    dataset = as.character(supported$dataset),
    compartment = as.character(supported$compartment),
    region = as.character(supported$region),
    layer = as.character(supported$layer),
    spatial_unit = as.character(supported$spatial_unit),
    original_comparison = as.character(supported$Comparison),
    phenotype_contrast = as.character(supported$phenotype_contrast),
    GO_ID = as.character(supported$ID),
    GO_description = as.character(supported$Description),
    NES = suppressWarnings(as.numeric(supported$NES)),
    p.adjust = suppressWarnings(as.numeric(supported$p.adjust)),
    leading_edge_proteins = as.character(supported$core_enrichment),
    leading_edge_genes = as.character(supported$core_enrichment_gene),
    legacy_program_class = as.character(supported$program_class),
    manuscript_theme_ids = status$manuscript_theme_ids[idx],
    manuscript_themes = status$manuscript_themes[idx],
    theme_roles = status$theme_roles[idx],
    anchor_GO_IDs = status$anchor_GO_IDs[idx],
    assignment_status = status$assignment_status[idx],
    mapping_method = status$mapping_method[idx],
    registry_version = status$registry_version[idx],
    GO_db_package_version = status$GO_db_package_version[idx],
    GO_source_name = status$GO_source_name[idx],
    GO_source_url = status$GO_source_url[idx],
    GO_source_date = status$GO_source_date[idx],
    relationship_types_approved = status$relationship_types_approved[idx],
    source_supplementary_file = as.character(supported$source_supplementary_file),
    source_manifest_file = as.character(supported$source_manifest_file),
    evidence_source_family = as.character(supported$evidence_source_family),
    stringsAsFactors = FALSE
  )
  out <- out[order(out$dataset, out$spatial_unit, out$p.adjust, -abs(out$NES), out$GO_ID, method = "radix"), , drop = FALSE]
  rownames(out) <- NULL
  out
}

go_semantic_redundancy_qa <- function(supported_audit, cutoff = 0.70) {
  require_go_ontology_dependencies(TRUE)
  if (!is.numeric(cutoff) || length(cutoff) != 1L || !is.finite(cutoff) || cutoff <= 0 || cutoff >= 1) {
    stop("Semantic redundancy cutoff must be a finite number between 0 and 1.", call. = FALSE)
  }
  ids <- sort(unique(as.character(supported_audit$GO_ID)), method = "radix")
  if (!length(ids)) stop("Semantic redundancy QA requires at least one supported GO ID.", call. = FALSE)
  sem_data <- GOSemSim::godata(annoDb = "org.Mm.eg.db", ont = "BP", computeIC = FALSE)
  similarity <- GOSemSim::mgoSim(ids, ids, semData = sem_data, measure = "Wang", combine = NULL)
  similarity <- as.matrix(similarity)
  similarity[!is.finite(similarity)] <- 0
  diag(similarity) <- 1
  if (length(ids) == 1L) {
    clusters <- stats::setNames(1L, ids)
  } else {
    clustering <- stats::hclust(stats::as.dist(1 - similarity), method = "average")
    clusters <- stats::cutree(clustering, h = 1 - cutoff)
  }
  min_fdr <- tapply(supported_audit$p.adjust, supported_audit$GO_ID, min, na.rm = TRUE)
  max_abs_nes <- tapply(abs(supported_audit$NES), supported_audit$GO_ID, max, na.rm = TRUE)
  description <- stats::setNames(
    vapply(ids, function(id) supported_audit$GO_description[match(id, supported_audit$GO_ID)], character(1)), ids
  )
  cluster_ids <- sort(unique(clusters))
  cluster_map <- stats::setNames(sprintf("WANG_BP_%03d", seq_along(cluster_ids)), cluster_ids)
  term_rows <- lapply(ids, function(id) {
    members <- names(clusters)[clusters == clusters[[id]]]
    ordered <- members[order(min_fdr[members], -max_abs_nes[members], members, description[members], method = "radix")]
    representative <- ordered[[1L]]
    data.frame(
      GO_ID = id,
      GO_description = unname(description[[id]]),
      semantic_cluster_id = unname(cluster_map[as.character(clusters[[id]])]),
      semantic_cluster_size = length(members),
      semantic_cluster_representative_GO_ID = representative,
      semantic_cluster_representative_term = unname(description[[representative]]),
      semantic_similarity_method = "Wang_BP",
      semantic_similarity_cutoff = cutoff,
      semantic_clustering_method = "average_linkage_hclust",
      semantic_inference_role = "descriptive_redundancy_QA_only",
      GOSemSim_package_version = as.character(utils::packageVersion("GOSemSim")),
      orgdb_package_version = as.character(utils::packageVersion("org.Mm.eg.db")),
      stringsAsFactors = FALSE
    )
  })
  term_audit <- do.call(rbind, term_rows)
  pair_index <- which(upper.tri(similarity, diag = FALSE), arr.ind = TRUE)
  pair_audit <- data.frame(
    GO_ID_1 = rownames(similarity)[pair_index[, 1]],
    GO_ID_2 = colnames(similarity)[pair_index[, 2]],
    Wang_similarity = similarity[pair_index],
    at_or_above_descriptive_cutoff = similarity[pair_index] >= cutoff,
    semantic_similarity_method = "Wang_BP",
    semantic_similarity_cutoff = cutoff,
    semantic_inference_role = "descriptive_redundancy_QA_only",
    stringsAsFactors = FALSE
  )
  pair_audit <- pair_audit[order(-pair_audit$Wang_similarity, pair_audit$GO_ID_1, pair_audit$GO_ID_2, method = "radix"), , drop = FALSE]
  rownames(term_audit) <- NULL
  rownames(pair_audit) <- NULL
  list(term_audit = term_audit, pair_audit = pair_audit)
}

build_sus_res_manuscript_theme_summary <- function(supported_audit, mapping, semantic_term_audit = NULL) {
  assignments <- mapping$assignments
  if (!nrow(assignments)) return(data.frame())
  term_theme <- unique(assignments[c(
    "GO_ID", "theme_id", "manuscript_theme", "theme_role", "display_order",
    "anchor_GO_ID", "anchor_label"
  )])
  collapse_by_term_theme <- split(term_theme, paste(term_theme$GO_ID, term_theme$theme_id, sep = "|"))
  term_theme <- do.call(rbind, lapply(collapse_by_term_theme, function(x) {
    data.frame(
      GO_ID = x$GO_ID[[1]], theme_id = x$theme_id[[1]], manuscript_theme = x$manuscript_theme[[1]],
      theme_role = x$theme_role[[1]], display_order = x$display_order[[1]],
      theme_anchor_GO_IDs = collapse_distinct_values(x$anchor_GO_ID),
      theme_anchor_labels = collapse_distinct_values(x$anchor_label),
      stringsAsFactors = FALSE
    )
  }))
  expanded <- merge(supported_audit, term_theme, by = "GO_ID", all = FALSE, sort = FALSE)
  if (!is.null(semantic_term_audit) && nrow(semantic_term_audit)) {
    expanded <- merge(expanded, semantic_term_audit[c("GO_ID", "semantic_cluster_id")], by = "GO_ID", all.x = TRUE, sort = FALSE)
  } else {
    expanded$semantic_cluster_id <- NA_character_
  }
  grouping <- c(
    "dataset", "compartment", "region", "layer", "spatial_unit", "phenotype_contrast",
    "theme_id", "manuscript_theme", "theme_role", "display_order"
  )
  grouping_values <- lapply(expanded[grouping], function(x) {
    x <- as.character(x)
    x[is.na(x)] <- "<NA>"
    x
  })
  keys <- do.call(paste, c(grouping_values, sep = "\u001f"))
  groups <- split(expanded, keys)
  summary_rows <- lapply(groups, function(x) {
    unique_terms <- x[!duplicated(x$GO_ID), , drop = FALSE]
    representative_order <- order(unique_terms$p.adjust, -abs(unique_terms$NES), unique_terms$GO_ID, unique_terms$GO_description, method = "radix")
    representative <- unique_terms[representative_order[[1]], , drop = FALSE]
    n_positive <- length(unique(unique_terms$GO_ID[unique_terms$NES > 0]))
    n_negative <- length(unique(unique_terms$GO_ID[unique_terms$NES < 0]))
    direction_consistency <- if (n_positive > 0L && n_negative > 0L) "mixed_direction" else if (n_positive > 0L) "positive_only" else if (n_negative > 0L) "negative_only" else "neutral_or_undirected"
    data.frame(
      dataset_compartment = x$dataset[[1]],
      compartment = x$compartment[[1]],
      region = x$region[[1]],
      layer = x$layer[[1]],
      spatial_unit = x$spatial_unit[[1]],
      comparison = x$phenotype_contrast[[1]],
      theme_id = x$theme_id[[1]],
      manuscript_theme = x$manuscript_theme[[1]],
      theme_role = x$theme_role[[1]],
      display_order = x$display_order[[1]],
      representative_term_id = representative$GO_ID[[1]],
      representative_term = representative$GO_description[[1]],
      representative_NES = representative$NES[[1]],
      representative_FDR = representative$p.adjust[[1]],
      direction = if (representative$NES[[1]] > 0) "higher in SUS" else if (representative$NES[[1]] < 0) "higher in RES" else "neutral",
      direction_consistency = direction_consistency,
      significant_term_count = nrow(unique_terms),
      n_positive_sig_terms = n_positive,
      n_negative_sig_terms = n_negative,
      n_semantic_clusters = length(unique(stats::na.omit(unique_terms$semantic_cluster_id))),
      representative_semantic_cluster_id = representative$semantic_cluster_id[[1]],
      representative_leading_edge_genes = representative$leading_edge_genes[[1]],
      representative_leading_edge_proteins = representative$leading_edge_proteins[[1]],
      representative_source_comparison = representative$original_comparison[[1]],
      representative_source_key = paste(representative$dataset[[1]], representative$original_comparison[[1]], representative$GO_ID[[1]], sep = "|"),
      representative_anchor_GO_IDs = representative$theme_anchor_GO_IDs[[1]],
      representative_anchor_labels = representative$theme_anchor_labels[[1]],
      source_supplementary_file = representative$source_supplementary_file[[1]],
      source_manifest_file = representative$source_manifest_file[[1]],
      evidence_source_family = "ranked_GSEA",
      registry_version = representative$registry_version[[1]],
      GO_db_package_version = representative$GO_db_package_version[[1]],
      GO_source_name = representative$GO_source_name[[1]],
      GO_source_url = representative$GO_source_url[[1]],
      GO_source_date = representative$GO_source_date[[1]],
      relationship_types_approved = representative$relationship_types_approved[[1]],
      mapping_method = "go_id_ontology",
      representative_selection_rule = "p.adjust < 0.05; lowest p.adjust; largest abs(NES); GO ID; Description",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, summary_rows)
  out$unit_key <- paste(out$dataset_compartment, out$spatial_unit, sep = "|")
  theme_groups <- split(out, out$theme_id)
  recurrence <- do.call(rbind, lapply(theme_groups, function(x) {
    positive_units <- unique(x$unit_key[x$direction_consistency == "positive_only"])
    negative_units <- unique(x$unit_key[x$direction_consistency == "negative_only"])
    mixed_units <- unique(x$unit_key[x$direction_consistency == "mixed_direction"])
    supported_units <- unique(x$unit_key)
    annotation <- if (length(supported_units) == 1L) {
      "spatially_specific"
    } else if (length(mixed_units) || (length(positive_units) && length(negative_units))) {
      "recurrent_mixed_direction"
    } else if (length(positive_units)) {
      "recurrent_consistent_higher_SUS"
    } else if (length(negative_units)) {
      "recurrent_consistent_higher_RES"
    } else {
      "recurrent_mixed_direction"
    }
    data.frame(
      theme_id = x$theme_id[[1]],
      n_supported_units = length(supported_units),
      n_supported_datasets = length(unique(x$dataset_compartment)),
      n_higher_SUS_units = length(positive_units),
      n_higher_RES_units = length(negative_units),
      n_mixed_direction_units = length(mixed_units),
      n_higher_SUS_datasets = length(unique(x$dataset_compartment[x$unit_key %in% positive_units])),
      n_higher_RES_datasets = length(unique(x$dataset_compartment[x$unit_key %in% negative_units])),
      directional_recurrence = annotation,
      stringsAsFactors = FALSE
    )
  }))
  out <- merge(out, recurrence, by = "theme_id", all.x = TRUE, sort = FALSE)
  out$sus_res_recurrent_units <- out$n_supported_units
  out$sus_res_recurrent_datasets <- out$n_supported_datasets
  out$sus_res_is_recurrent <- out$n_supported_units >= 2L
  out$recurrence_annotation <- out$directional_recurrence
  out$publication_include <- out$theme_role == "primary"
  out$panel_data_basis <- "ranked_proteome_wide_GO_GSEA"
  out$recurrence_inference_role <- "descriptive_only_not_a_significance_gate"
  out$NES_interpretation <- "positive = higher in SUS; negative = higher in RES"
  out$unit_key <- NULL
  out <- out[order(out$display_order, out$dataset_compartment, out$spatial_unit, out$representative_FDR, method = "radix"), , drop = FALSE]
  rownames(out) <- NULL
  key <- paste(out$dataset_compartment, out$spatial_unit, out$theme_id, sep = "|")
  if (anyDuplicated(key)) stop("Manuscript theme summary has duplicate dataset + spatial unit + theme rows.", call. = FALSE)
  out
}

summarize_go_theme_assignment_status <- function(supported_audit) {
  total <- nrow(supported_audit)
  statuses <- c("single_theme", "multi_theme", "unclassified", "qc_review")
  counts <- vapply(statuses, function(status) sum(supported_audit$assignment_status == status), integer(1))
  data.frame(
    metric = c("total_supported_term_occurrences", "unique_supported_GO_IDs", paste0(statuses, "_occurrences"), paste0(statuses, "_proportion")),
    value = c(total, length(unique(supported_audit$GO_ID)), counts, if (total) counts / total else rep(NA_real_, length(counts))),
    stringsAsFactors = FALSE
  )
}

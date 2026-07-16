# Shared enrichment IO and conservative biological-program mapping helpers.

if (!exists("repo_path", mode = "function")) {
  paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
  source(paths_file)
}
if (!exists("safe_name", mode = "function")) {
  source(repo_path("R", "validation_utils.R"))
}

canonical_clusterprofiler_manifest_contract_version <- function() {
  "clusterProfiler_manifest_v3_term_gene_provenance"
}

canonical_comparego_manifest_contract_version <- function() {
  "compareGO_manifest_v3_term_gene_provenance"
}

canonical_comparego_result_types <- function() "GSEA_GO"

clusterprofiler_manifest_columns <- function() c(
  "analysis_id", "dataset", "run_id", "ontology", "result_type", "contrast", "comparison",
  "route_category", "route_unit", "condition", "direction", "simplified", "plot_suffix",
  "used_for_plot", "input_gene_file", "gene_input_file", "input_hash",
  "collapsed_gene_input_file", "collapsed_gene_provenance_file", "term_gene_provenance_file",
  "enrichment_contract_version", "gene_annotation_contract_version", "protein_group_contract_version",
  "gene_mapping_policy", "primary_gene_level_eligibility_rule", "ambiguous_group_policy",
  "duplicate_gene_collapse_rule", "rank_statistic_column", "rank_statistic_type",
  "rank_statistic_fallback_used", "ORA_direction", "universe_definition", "config_file",
  "config_hash", "output_table", "output_plot", "n_genes", "n_terms", "analysis_status",
  "empty_result", "error_message", "checkpoint_status", "created_at"
)

comparego_manifest_columns <- function() c(
  clusterprofiler_manifest_columns(), "input_manifest", "comparego_contract_version",
  "comparego_analysis_status", "term_comparison_file", "term_gene_provenance_output_file",
  "analysis_status_summary_file"
)

clusterprofiler_output_path_audit <- function(paths, safe_limit = 240L) {
  paths <- unique(as.character(paths))
  paths <- paths[!is.na(paths) & nzchar(paths)]
  normalized <- normalizePath(paths, winslash = "/", mustWork = FALSE)
  out <- data.frame(
    path = normalized,
    path_length = nchar(normalized, type = "chars"),
    safe_limit = as.integer(safe_limit),
    within_safe_limit = nchar(normalized, type = "chars") <= as.integer(safe_limit),
    stringsAsFactors = FALSE
  )
  out[order(out$path_length, decreasing = TRUE, out$path, method = "radix"), , drop = FALSE]
}

validate_clusterprofiler_output_path_lengths <- function(paths, safe_limit = 240L) {
  audit <- clusterprofiler_output_path_audit(paths, safe_limit = safe_limit)
  excessive <- audit[!audit$within_safe_limit, , drop = FALSE]
  if (nrow(excessive)) {
    stop(
      "Expected clusterProfiler output path exceeds the safe Windows/R path limit of ",
      as.integer(safe_limit), " characters (maximum expected length: ", max(excessive$path_length), "). ",
      "Run the repository from a shorter project root such as P:\\ before launching workers. ",
      "Longest path: ", excessive$path[[1]],
      call. = FALSE
    )
  }
  audit
}

clusterprofiler_worker_error <- function(result, fallback) {
  if (!is.list(result) || is.null(result$error) || !length(result$error) ||
      is.na(result$error[[1]]) || !nzchar(as.character(result$error[[1]]))) {
    return(fallback)
  }
  as.character(result$error[[1]])
}

assess_clusterprofiler_worker_result <- function(result, expected_comparison) {
  failed <- function(message, worker_status = "MALFORMED", manifest_has_failure = FALSE) {
    data.frame(
      comparison = as.character(expected_comparison), worker_status = worker_status,
      analysis_status = "failed", error = as.character(message),
      manifest_has_failure = isTRUE(manifest_has_failure), stringsAsFactors = FALSE
    )
  }
  if (!is.list(result)) return(failed("Worker returned a missing or malformed result object."))
  required <- c("status", "comparison", "manifest")
  if (length(setdiff(required, names(result)))) {
    return(failed("Worker result object is missing required fields: status, comparison, manifest."))
  }
  worker_status <- as.character(result$status)[1]
  comparison <- as.character(result$comparison)[1]
  if (is.na(worker_status) || !nzchar(worker_status) || is.na(comparison) || !nzchar(comparison)) {
    return(failed("Worker result contains an empty status or comparison."))
  }
  if (!identical(comparison, as.character(expected_comparison))) {
    return(failed(
      paste0("Worker comparison identity mismatch: expected ", expected_comparison, ", received ", comparison, "."),
      worker_status = worker_status
    ))
  }
  manifest <- result$manifest
  manifest_is_table <- is.data.frame(manifest)
  manifest_has_failure <- manifest_is_table && all(c("result_type", "analysis_status") %in% names(manifest)) &&
    any(manifest$result_type == "GSEA_GO" & manifest$analysis_status == "failed", na.rm = TRUE)
  if (!worker_status %in% c("SUCCESS", "SKIPPED")) {
    return(failed(
      clusterprofiler_worker_error(result, paste0("Worker returned status ", worker_status, ".")),
      worker_status = worker_status, manifest_has_failure = manifest_has_failure
    ))
  }
  if (!manifest_is_table || !all(c("result_type", "analysis_status", "n_terms") %in% names(manifest))) {
    return(failed("Successful worker returned a missing or malformed manifest.", worker_status = worker_status))
  }
  primary <- manifest[manifest$result_type == "GSEA_GO", , drop = FALSE]
  if (nrow(primary) != 1L) {
    return(failed("Successful worker must return exactly one primary GSEA_GO manifest row.", worker_status = worker_status))
  }
  primary_status <- as.character(primary$analysis_status[[1]])
  n_terms <- suppressWarnings(as.integer(primary$n_terms[[1]]))
  if (identical(primary_status, "success_with_terms") && !is.na(n_terms) && n_terms > 0L) {
    category <- "success_with_terms"
  } else if (identical(primary_status, "success_zero_terms") && identical(n_terms, 0L)) {
    category <- "success_zero_terms"
  } else {
    return(failed(
      paste0("Primary GSEA_GO manifest has inconsistent analysis_status/n_terms: ",
        primary_status, "/", ifelse(is.na(n_terms), "NA", n_terms), "."),
      worker_status = worker_status, manifest_has_failure = manifest_has_failure
    ))
  }
  data.frame(
    comparison = comparison, worker_status = worker_status, analysis_status = category,
    error = NA_character_, manifest_has_failure = FALSE, stringsAsFactors = FALSE
  )
}

assess_clusterprofiler_worker_results <- function(results, expected_comparisons) {
  if (!is.list(results)) results <- list(results)
  expected_comparisons <- as.character(expected_comparisons)
  n_rows <- max(length(results), length(expected_comparisons))
  if (!n_rows) {
    return(data.frame(
      comparison = character(), worker_status = character(), analysis_status = character(),
      error = character(), manifest_has_failure = logical(), stringsAsFactors = FALSE
    ))
  }
  rows <- lapply(seq_len(n_rows), function(i) {
    if (i > length(expected_comparisons)) {
      return(data.frame(
        comparison = paste0("unexpected_worker_result_", i), worker_status = "MALFORMED",
        analysis_status = "failed", error = "Received an unexpected extra worker result object.",
        manifest_has_failure = FALSE, stringsAsFactors = FALSE
      ))
    }
    result <- if (i <= length(results)) results[[i]] else NULL
    assess_clusterprofiler_worker_result(result, expected_comparisons[[i]])
  })
  do.call(rbind, rows)
}

clusterprofiler_master_status_counts <- function(assessment) {
  statuses <- c("success_with_terms", "success_zero_terms", "failed")
  counts <- setNames(integer(length(statuses)), statuses)
  observed <- table(factor(assessment$analysis_status, levels = statuses))
  counts[] <- as.integer(observed)
  counts
}

clusterprofiler_master_exit_status <- function(assessment) {
  if (any(assessment$analysis_status == "failed")) 1L else 0L
}

clusterprofiler_master_summary_lines <- function(assessment) {
  counts <- clusterprofiler_master_status_counts(assessment)
  final <- if (counts[["failed"]] > 0L) {
    paste0("RUN FAILED: ", counts[["failed"]], " comparison(s) failed.")
  } else {
    "ALL COMPARISONS COMPLETED SUCCESSFULLY."
  }
  c(
    paste0("success_with_terms: ", counts[["success_with_terms"]]),
    paste0("success_zero_terms: ", counts[["success_zero_terms"]]),
    paste0("failed: ", counts[["failed"]]),
    final
  )
}

write_csv_strict <- function(x, path, label = "CSV output") {
  parent <- dirname(path)
  if (!dir.exists(parent)) {
    stop(label, " parent directory does not exist: ", parent, call. = FALSE)
  }
  tryCatch(
    utils::write.csv(x, path, row.names = FALSE),
    error = function(e) stop("Failed to write ", label, " to ", path, ": ", conditionMessage(e), call. = FALSE)
  )
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size <= 0) {
    stop("Failed to verify written ", label, ": ", path, call. = FALSE)
  }
  written <- tryCatch(
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) stop("Failed to read back written ", label, " at ", path, ": ", conditionMessage(e), call. = FALSE)
  )
  if (!identical(names(written), names(x)) || nrow(written) != nrow(x)) {
    stop("Written ", label, " failed row/column verification: ", path, call. = FALSE)
  }
  invisible(path)
}

term_gene_provenance_columns <- function() c(
  "dataset", "comparison", "result_type", "ontology", "term_id", "term_description",
  "official_gene_symbol", "official_entrez_id", "ProteinGroupID", "member_accessions",
  "protein_group_gene_annotation_status", "gene_level_claim_allowed", "rank_statistic",
  "core_enrichment_member", "enrichment_contract_version", "gene_annotation_contract_version"
)

validate_term_gene_provenance_contract <- function(x, strict = TRUE) {
  missing <- setdiff(term_gene_provenance_columns(), names(x))
  if (length(missing)) {
    stop("Term-gene provenance is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!is.character(x$official_entrez_id)) {
    stop("Term-gene provenance official_entrez_id must remain character.", call. = FALSE)
  }
  if (isTRUE(strict) && any(x$gene_level_claim_allowed %in% FALSE, na.rm = TRUE)) {
    stop("Ineligible protein groups are not permitted in strict term-gene provenance.", call. = FALSE)
  }
  invisible(TRUE)
}

validate_gsea_result_table_contract <- function(x, context = "GSEA result") {
  required <- c("ID", "Description", "NES", "p.adjust", "setSize", "core_enrichment")
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop(context, " is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

validate_clusterprofiler_manifest_contract <- function(manifest, strict = TRUE, require_files = TRUE) {
  required <- c(
    "dataset", "comparison", "result_type", "ontology", "analysis_status", "n_terms",
    "output_table", "collapsed_gene_input_file", "collapsed_gene_provenance_file",
    "term_gene_provenance_file", "enrichment_contract_version",
    "gene_annotation_contract_version"
  )
  missing <- setdiff(required, names(manifest))
  if (length(missing)) {
    stop("clusterProfiler manifest is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  expected_version <- canonical_clusterprofiler_manifest_contract_version()
  if (isTRUE(strict) && any(is.na(manifest$enrichment_contract_version) |
      manifest$enrichment_contract_version != expected_version)) {
    stop("Stale clusterProfiler manifest contract; expected ", expected_version, ".", call. = FALSE)
  }
  allowed_status <- c("success_with_terms", "success_zero_terms", "failed")
  invalid_status <- setdiff(unique(as.character(manifest$analysis_status)), allowed_status)
  if (length(invalid_status)) {
    stop("clusterProfiler manifest has unsupported analysis_status: ", paste(invalid_status, collapse = ", "), call. = FALSE)
  }
  supported_results <- c("GSEA_GO", "GSEA_KEGG")
  invalid_results <- setdiff(unique(as.character(manifest$result_type)), supported_results)
  if (length(invalid_results)) {
    stop("clusterProfiler manifest has unsupported result_type: ", paste(invalid_results, collapse = ", "), call. = FALSE)
  }
  success <- manifest$analysis_status %in% c("success_with_terms", "success_zero_terms")
  zero <- manifest$analysis_status == "success_zero_terms"
  with_terms <- manifest$analysis_status == "success_with_terms"
  if (any(zero & (is.na(manifest$n_terms) | manifest$n_terms != 0))) {
    stop("success_zero_terms rows must record n_terms = 0.", call. = FALSE)
  }
  if (any(with_terms & (is.na(manifest$n_terms) | manifest$n_terms <= 0))) {
    stop("success_with_terms rows must record n_terms > 0.", call. = FALSE)
  }
  if (isTRUE(require_files) && any(success)) {
    file_columns <- c("output_table", "collapsed_gene_input_file", "collapsed_gene_provenance_file", "term_gene_provenance_file")
    for (column in file_columns) {
      paths <- as.character(manifest[[column]][success])
      missing_paths <- is.na(paths) | !nzchar(paths) | !file.exists(paths)
      if (any(missing_paths)) {
        stop("Successful clusterProfiler manifest row references missing ", column, ": ",
          paste(paths[missing_paths], collapse = ", "), call. = FALSE)
      }
    }
  }
  invisible(TRUE)
}

validate_comparego_manifest_contract <- function(manifest, require_files = TRUE) {
  required <- c(
    "dataset", "comparison", "result_type", "ontology", "analysis_status",
    "comparego_analysis_status", "input_manifest", "term_comparison_file",
    "term_gene_provenance_output_file", "analysis_status_summary_file",
    "enrichment_contract_version", "comparego_contract_version"
  )
  missing <- setdiff(required, names(manifest))
  if (length(missing)) {
    stop("compareGO manifest is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (any(manifest$comparego_contract_version != canonical_comparego_manifest_contract_version())) {
    stop("Stale compareGO manifest contract.", call. = FALSE)
  }
  if (any(!manifest$result_type %in% canonical_comparego_result_types())) {
    stop("compareGO manifest contains unsupported result_type.", call. = FALSE)
  }
  allowed_actions <- c("included", "completed_zero_terms", "recorded_failed")
  if (any(!manifest$comparego_analysis_status %in% allowed_actions)) {
    stop("compareGO manifest contains unsupported comparego_analysis_status.", call. = FALSE)
  }
  if (isTRUE(require_files)) {
    file_columns <- c("input_manifest", "term_comparison_file", "term_gene_provenance_output_file", "analysis_status_summary_file")
    for (column in file_columns) {
      paths <- unique(as.character(manifest[[column]]))
      if (any(is.na(paths) | !nzchar(paths) | !file.exists(paths))) {
        stop("compareGO manifest references missing ", column, ".", call. = FALSE)
      }
    }
  }
  invisible(TRUE)
}

read_csv_contract <- function(path, character_columns = character()) {
  x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
    colClasses = if (length(character_columns)) c(setNames(rep("character", length(character_columns)), character_columns)) else NA)
  x
}

collect_canonical_comparego_outputs <- function(manifest, strict = TRUE, require_files = TRUE) {
  validate_clusterprofiler_manifest_contract(manifest, strict = strict, require_files = require_files)
  supported <- canonical_comparego_result_types()
  unsupported <- manifest[!manifest$result_type %in% supported, , drop = FALSE]
  selected <- manifest[manifest$result_type %in% supported, , drop = FALSE]
  if (!nrow(selected)) stop("No supported canonical compareGO result types were selected.", call. = FALSE)
  selected <- selected[order(selected$dataset, selected$comparison, selected$ontology, selected$result_type, method = "radix"), , drop = FALSE]

  term_rows <- list()
  provenance_rows <- list()
  status_rows <- list()
  for (i in seq_len(nrow(selected))) {
    row <- selected[i, , drop = FALSE]
    status_rows[[i]] <- data.frame(
      dataset = as.character(row$dataset), comparison = as.character(row$comparison),
      result_type = as.character(row$result_type), ontology = as.character(row$ontology),
      analysis_status = as.character(row$analysis_status), n_terms = as.integer(row$n_terms),
      comparego_action = if (row$analysis_status == "failed") "recorded_failed" else if (row$analysis_status == "success_zero_terms") "completed_zero_terms" else "included",
      stringsAsFactors = FALSE
    )
    if (row$analysis_status == "failed") next
    terms <- read_csv_contract(as.character(row$output_table))
    validate_gsea_result_table_contract(terms, paste0(row$dataset, "/", row$comparison))
    if (nrow(terms) != as.integer(row$n_terms)) {
      stop("Manifest n_terms does not match result table for ", row$dataset, "/", row$comparison, ".", call. = FALSE)
    }
    if (nrow(terms)) {
      terms$term_id <- as.character(terms$ID)
      terms$term_description <- as.character(terms$Description)
      terms$dataset <- as.character(row$dataset)
      terms$comparison <- as.character(row$comparison)
      terms$result_type <- as.character(row$result_type)
      terms$ontology <- as.character(row$ontology)
      term_rows[[length(term_rows) + 1L]] <- terms
    }
    provenance <- read_csv_contract(as.character(row$term_gene_provenance_file), character_columns = "official_entrez_id")
    validate_term_gene_provenance_contract(provenance, strict = strict)
    if (!nrow(terms) && nrow(provenance)) {
      stop("Zero-term GSEA result has non-empty term-gene provenance for ",
        row$dataset, "/", row$comparison, ".", call. = FALSE)
    }
    if (nrow(provenance) && any(!provenance$term_id %in% as.character(terms$ID))) {
      stop("Term-gene provenance contains term IDs absent from the declared GSEA result.", call. = FALSE)
    }
    if (nrow(provenance)) {
      expected <- c(dataset = as.character(row$dataset), comparison = as.character(row$comparison),
        result_type = as.character(row$result_type), ontology = as.character(row$ontology))
      for (column in names(expected)) {
        if (any(as.character(provenance[[column]]) != expected[[column]])) {
          stop("Term-gene provenance metadata does not match manifest field ", column, ".", call. = FALSE)
        }
      }
      provenance_rows[[length(provenance_rows) + 1L]] <- provenance
    }
  }
  empty_terms <- data.frame(
    ID = character(), Description = character(), NES = numeric(), p.adjust = numeric(),
    setSize = integer(), core_enrichment = character(), term_id = character(),
    term_description = character(), dataset = character(), comparison = character(),
    result_type = character(), ontology = character(), stringsAsFactors = FALSE
  )
  empty_provenance <- as.data.frame(setNames(lapply(term_gene_provenance_columns(), function(column) {
    if (column %in% c("gene_level_claim_allowed", "core_enrichment_member")) logical()
    else if (column == "rank_statistic") numeric() else character()
  }), term_gene_provenance_columns()), stringsAsFactors = FALSE)
  terms <- if (length(term_rows)) do.call(rbind, term_rows) else empty_terms
  provenance <- if (length(provenance_rows)) do.call(rbind, provenance_rows) else empty_provenance
  status <- do.call(rbind, status_rows)
  if (nrow(unsupported)) {
    unsupported_status <- data.frame(
      dataset = as.character(unsupported$dataset), comparison = as.character(unsupported$comparison),
      result_type = as.character(unsupported$result_type), ontology = as.character(unsupported$ontology),
      analysis_status = as.character(unsupported$analysis_status), n_terms = as.integer(unsupported$n_terms),
      comparego_action = "skipped_unsupported_result_type", stringsAsFactors = FALSE
    )
    status <- rbind(status, unsupported_status)
  }
  if (nrow(terms)) terms <- terms[order(terms$dataset, terms$comparison, terms$ontology, terms$ID, method = "radix"), , drop = FALSE]
  if (nrow(provenance)) provenance <- provenance[order(provenance$dataset, provenance$comparison, provenance$ontology,
    provenance$term_id, provenance$official_gene_symbol, provenance$ProteinGroupID, method = "radix"), , drop = FALSE]
  status <- status[order(status$dataset, status$comparison, status$result_type, status$ontology, method = "radix"), , drop = FALSE]
  rownames(terms) <- rownames(provenance) <- rownames(status) <- NULL
  list(input_manifest = selected, terms = terms, provenance = provenance, status = status)
}

biological_program_patterns <- function() {
  data.frame(
    biological_program = c(
      "RNA_RNP_processing",
      "Ribosome_Translation",
      "Mitochondria_OXPHOS_Metabolism",
      "Proteostasis_Ubiquitin_Folding",
      "Synapse_Vesicle_Organization",
      "Cytoskeleton_Motility",
      "Development_Patterning",
      "HPA_Glucocorticoid_Response",
      "Neuroimmune_Complement_Phagosome",
      "ECM_Vascular_Barrier",
      "Oxidative_Redox_Stress",
      "Lipid_Myelin_Membrane",
      "Autophagy_Lysosome"
    ),
    pattern = c(
      "rna|ribonucleoprotein|rnp|splice|splicing|mrna|ncrna|rrna|trna|nucleolus|ribonucle|rna processing",
      "translation|ribosom|peptide biosynthetic|cytoplasmic translation|translational initiation|elongation factor|initiation factor",
      "mitochond|oxidative phosphorylation|oxphos|electron transport|respiratory chain|atp synthesis|oxidoreduct|metabol|glycolys|tricarboxylic|acetyl.coa|energy",
      "proteas|ubiquitin|folding|chaperone|heat shock|proteostasis|protein quality",
      "synap|vesicle|neurotransmitter|axon|dendrit|postsynap|presynap|exocytosis|endocytosis|synaptic organization",
      "cytoskeleton|actin|tubulin|microtubule|motility|adhesion|migration|extracellular matrix",
      "develop|pattern|morphogen|differentiation|neurogenesis|gliogenesis|axon guidance|cell fate|regionalization|dorsal.ventral|anterior.posterior",
      "glucocorticoid|corticosterone|cortisol|steroid hormone|hpa axis|stress hormone|nr3c1|nuclear receptor subfamily 3 group c",
      "microglia|immune|inflamm|cytokine|chemokine|complement|phagocyt|phagosome|lysosomal engulfment|antigen presentation|mhc",
      "extracellular matrix|collagen|laminin|basement membrane|vascular|blood vessel|endothelial|pericyte|blood.brain barrier|barrier|integrin",
      "oxidative stress|redox|reactive oxygen|ros|peroxid|glutathione|superoxide|oxidant|antioxidant",
      "lipid|fatty acid|cholesterol|membrane|myelin|oligodendro|sphingolipid|phospholipid",
      "autophag|lysosom|endosom|phagolysosom|vacuolar|proteolysis"
    ),
    stringsAsFactors = FALSE
  )
}

map_terms_to_programs <- function(df, description_col = "Description") {
  if (!description_col %in% names(df)) {
    df$biological_program <- NA_character_
    return(df)
  }
  patterns <- biological_program_patterns()
  desc <- tolower(as.character(df[[description_col]]))
  hits <- lapply(seq_len(nrow(patterns)), function(i) grepl(patterns$pattern[[i]], desc, ignore.case = TRUE))
  hit_mat <- do.call(cbind, hits)
  program <- rep(NA_character_, length(desc))
  first_hit <- max.col(hit_mat, ties.method = "first")
  any_hit <- rowSums(hit_mat, na.rm = TRUE) > 0
  program[any_hit] <- patterns$biological_program[first_hit[any_hit]]
  df$biological_program <- program
  df
}

read_csv_if_exists <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  readr::read_csv(path, show_col_types = FALSE)
}

first_existing_path <- function(paths) {
  paths <- paths[!is.na(paths) & nzchar(paths)]
  if (!length(paths)) return(NA_character_)
  paths <- unique(normalizePath(paths, winslash = "/", mustWork = FALSE))
  hit <- paths[file.exists(paths)]
  if (!length(hit)) return(NA_character_)
  hit[[1]]
}

latest_file <- function(root, pattern) {
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  if (!dir.exists(root)) return(NA_character_)
  files <- list.files(root, pattern = pattern, full.names = TRUE, recursive = TRUE)
  files <- files[file.exists(files)]
  if (!length(files)) return(NA_character_)
  info <- file.info(files)
  normalizePath(rownames(info)[order(info$mtime, decreasing = TRUE)[1]], winslash = "/", mustWork = FALSE)
}

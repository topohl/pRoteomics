#!/usr/bin/env Rscript
# ================================================================
# Script: 04_differential_expression_enrichment/08_external_stress_disease_signature_overlap.r
# Stage: enrichment
# Scope: global
# Consumes: optional external cached stress/disease signatures and enrichment program summaries.
# Produces: external stress/disease signature overlap tables and source data.
# ================================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "integration_utils.R"))

SCRIPT_ID <- "04_differential_expression_enrichment/08_external_stress_disease_signature_overlap.r"
paths <- create_module_dirs("04_differential_expression_enrichment", "external_stress_disease_signature_overlap/global")
signature_dir <- Sys.getenv("PROTEOMICS_EXTERNAL_SIGNATURE_DIR", unset = path_external("stress_disease_signatures"))
signature_files <- Sys.glob(file.path(signature_dir, "*.csv"))

inputs <- c(signature_dir = signature_dir)
for (ds in valid_datasets()) {
  inputs[[paste0("program_summary_", ds)]] <- path_results("tables", "04_differential_expression_enrichment", "biological_program_summary", ds, "program_summary.csv")
}

if (is_dry_run()) {
  dry_run_inputs(SCRIPT_ID, inputs)
  dry_run_line("Signature CSV count", length(signature_files), if (length(signature_files)) "PASS" else "WARN")
  quit(status = 0, save = "no")
}

read_signatures <- function(files) {
  if (!length(files)) return(list(data = NULL, status = data.frame(
    dataset = "global", evidence_domain = "external_signature", input_type = "stress_disease_signature_dir",
    path = normalizePath(signature_dir, winslash = "/", mustWork = FALSE), required = FALSE,
    status = "missing_optional", message = "No cached external signatures found; no live downloads attempted.",
    n_rows = 0L, stringsAsFactors = FALSE
  )))
  rows <- list()
  status <- empty_status()
  for (f in files) {
    loaded <- read_csv_optional(f, "global", "external_signature", basename(f), required = FALSE)
    status <- rbind(status, loaded$status)
    df <- loaded$data
    if (is.null(df) || !nrow(df)) next
    gene_col <- first_col(df, c("GeneSymbol", "gene_symbol", "Gene", "gene", "symbol", "protein", "Protein"))
    sig_col <- first_col(df, c("signature", "Signature", "signature_name", "program", "disease", "term"))
    if (is.na(gene_col)) next
    if (is.na(sig_col)) df$signature <- tools::file_path_sans_ext(basename(f)) else df$signature <- as.character(df[[sig_col]])
    rows[[length(rows) + 1L]] <- data.frame(
      signature = chr_clean(df$signature, tools::file_path_sans_ext(basename(f))),
      gene_token = normalize_token(df[[gene_col]]),
      source_file = normalizePath(f, winslash = "/", mustWork = FALSE),
      stringsAsFactors = FALSE
    )
  }
  list(data = do.call(rbind, rows), status = status)
}

program_gene_rows <- function(ds) {
  f <- inputs[[paste0("program_summary_", ds)]]
  loaded <- read_csv_optional(f, ds, "enrichment_program", "biological_program_summary", required = FALSE)
  df <- loaded$data
  if (is.null(df) || !nrow(df)) {
    return(list(data = availability_evidence(ds, "external_signature_overlap", f, "program summary unavailable"), status = loaded$status))
  }
  gene_col <- first_col(df, c("key_genes", "core_genes", "key_proteins_genes", "leading_edge", "genes"))
  prog_col <- first_col(df, c("biological_program", "program", "top_term", "Description"))
  contrast_col <- first_col(df, c("comparison", "contrast", "phenotype_contrast"))
  if (is.na(gene_col) || is.na(prog_col)) {
    return(list(data = availability_evidence(ds, "external_signature_overlap", f, "program summary lacks gene/program columns"), status = loaded$status))
  }
  rows <- list()
  for (i in seq_len(nrow(df))) {
    toks <- split_gene_tokens(df[[gene_col]][[i]])
    if (!length(toks)) next
    rows[[length(rows) + 1L]] <- data.frame(
      dataset = ds,
      program_label = as.character(df[[prog_col]][[i]]),
      contrast = if (!is.na(contrast_col)) as.character(df[[contrast_col]][[i]]) else NA_character_,
      gene_token = toks,
      source_file = f,
      stringsAsFactors = FALSE
    )
  }
  list(data = do.call(rbind, rows), status = loaded$status)
}

sig <- read_signatures(signature_files)
status <- sig$status
programs <- lapply(valid_datasets(), program_gene_rows)
status <- rbind(status, do.call(rbind, lapply(programs, `[[`, "status")))
program_genes <- do.call(rbind, lapply(programs, `[[`, "data"))

if (is.null(sig$data) || !nrow(sig$data) || is.null(program_genes) || !nrow(program_genes) || !"gene_token" %in% names(program_genes)) {
  overlap <- availability_evidence("global", "external_signature_overlap", signature_dir, "Optional cached signature or program-gene evidence unavailable.")
} else {
  merged <- merge(program_genes, sig$data, by = "gene_token")
  if (!nrow(merged)) {
    overlap <- availability_evidence("global", "external_signature_overlap", signature_dir, "No overlapping gene/protein tokens detected.")
  } else {
    overlap <- aggregate(
      gene_token ~ dataset + program_label + contrast + signature + source_file.x + source_file.y,
      merged,
      function(x) paste(sort(unique(x)), collapse = ";")
    )
    names(overlap)[names(overlap) == "gene_token"] <- "overlap_genes"
    overlap$support_count <- lengths(strsplit(overlap$overlap_genes, ";", fixed = TRUE))
    overlap$evidence_domain <- "external_stress_disease_signature_overlap"
    overlap$evidence_id <- paste(overlap$dataset, overlap$signature, overlap$program_label, overlap$contrast, sep = "::")
    overlap$entity_type <- "signature"
    overlap$entity_id <- overlap$signature
    overlap$spatial_unit <- NA_character_
    overlap$effect_direction <- NA_character_
    overlap$effect_size <- NA_real_
    overlap$p_value <- NA_real_
    overlap$fdr <- NA_real_
    overlap$source_file <- paste(overlap$source_file.x, overlap$source_file.y, sep = ";")
    overlap$evidence_status <- ifelse(overlap$support_count >= 3, "supportive_overlap", "weak_overlap")
    overlap$interpretation_note <- paste0("Cached external signature overlap; no live downloads. Overlap genes: ", overlap$overlap_genes)
    overlap$qc_flag <- "PASS"
    overlap <- standardize_evidence(overlap)
  }
}

invisible(write_integration_table(overlap, paths, "external_stress_disease_signature_overlap.csv"))
write_csv_safe(status, file.path(paths$reports, "input_status.csv"))
write_csv_safe(status, file.path(paths$source_data, "external_stress_disease_signature_overlap_input_status.csv"))
write_integration_manifest(
  paths,
  inputs = as.list(inputs),
  outputs = list(tables = paths$tables, source_data = paths$source_data, reports = paths$reports),
  parameters = list(signature_dir = signature_dir, live_downloads = FALSE),
  notes = "Optional cached external signatures are used when present. Missing signatures produce status/evidence rows and do not fail the pipeline."
)

message("External stress/disease signature overlap complete: ", paths$tables)

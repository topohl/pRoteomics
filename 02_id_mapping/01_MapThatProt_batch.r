# ================================================================
# Script: 02_id_mapping/01_MapThatProt_batch.r
# Stage: core
# Scope: dataset_specific
# Consumes: required PROTEOMICS_GCT_EXTRACT_ROOT/<dataset>/<direction>/*.csv (historical default data/processed/01_preprocessing/gct_extractR); data/external/MOUSE_10090_idmapping.dat; optional manual protein and gene-annotation mapping tables.
# Produces: PROTEOMICS_MAPPING_OUTPUT_ROOT/{mapped,unmapped,member_bridge}/<dataset>/<direction>/per_file/ (historical default data/processed/02_id_mapping) plus namespaced QC outputs.
# Dataset behavior: runs for neuron_neuropil,neuron_soma,microglia according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Maps GCT-derived protein IDs to UniProt/gene symbols; enrichment consumes these mapped files.
# ================================================================

#' Batch UniProt ID Mapping Script for Proteomics Data
#'
#' This script maps protein identifiers to canonical UniProt accessions using the
#' local mouse UniProt dictionary and optional manual overrides. Official mouse
#' SYMBOL and ENTREZ annotations are then resolved from org.Mm.eg.db.
#'
#' Key features:
#' - Parallel processing of multiple CSV files using doParallel
#' - Deterministic local accession, entry-name, and gene-name/synonym mapping
#' - Accession-first official mouse SYMBOL/ENTREZ annotation from org.Mm.eg.db
#' - Manual mapping override support via Excel file
#' - Comprehensive mapping statistics and unmapped protein tracking
#' - Species-specific filtering (_MOUSE suffix enforcement)
#' 
#' Dependencies: dplyr, stringr, tidyr, purrr, readr, R.utils, foreach, doParallel, readxl, AnnotationDbi, org.Mm.eg.db
#' 
#' @title MapThatProt_batch
#' 
#' @author Tobias Pohl
#' 
#' @date 2025-12-15
#'
#' Consumes:
#'   - contrast CSVs from data/processed/01_preprocessing/gct_extractR/<comparison>/<forward|reverse>/
#'   - UniProt idmapping from data/external/MOUSE_10090_idmapping.dat
#'   - optional manual mapping workbook from data/metadata/manual_mapping.xlsx
#' Produces:
#'   - mapped contrast CSVs under data/processed/02_id_mapping/mapped/<comparison>/<forward|reverse>/per_file/
#'   - unmapped tracking under data/processed/02_id_mapping/unmapped/<comparison>/<forward|reverse>/per_file/
#'   - mapping summaries/reports/logs under results/{tables,reports,logs}/02_id_mapping/MapThatProt_batch/<comparison>/
#' File contract:
#'   - docs/file_contracts.tsv object mapped_contrast_csv

cat("====================================================\n")
cat("Starting MapThatProt_batch execution...\n")
cat("====================================================\n")

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "protein_mapping_utils.R"))
source(repo_path("R", "mapping_branch_utils.R"))
MODULE_ID <- "02_id_mapping"
SUBSTEP_ID <- "MapThatProt_batch"

load_required_packages <- function(pkgs) {
    missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing) > 0) {
        stop("Missing required R package(s): ", paste(missing, collapse = ", "),
             ". Install them explicitly before running this script.", call. = FALSE)
    }
    invisible(lapply(pkgs, library, character.only = TRUE))
}

# --- Configuration & Experimental Settings ---
mapped_comparisons <- current_dataset_from_cli()
map_direction <- Sys.getenv("PROTEOMICS_MAP_DIRECTION", unset = "forward")
if (!map_direction %in% c("forward", "reverse")) stop("PROTEOMICS_MAP_DIRECTION must be 'forward' or 'reverse'.", call. = FALSE)
map_reverse <- identical(map_direction, "reverse")

truthy_env <- function(name, default = FALSE) {
    value <- Sys.getenv(name, unset = if (isTRUE(default)) "true" else "false")
    tolower(value) %in% c("1", "true", "yes", "y")
}

has_cli_flag <- function(flags) {
    any(commandArgs(trailingOnly = TRUE) %in% flags) ||
        any(commandArgs(trailingOnly = FALSE) %in% flags)
}

force_rerun <- truthy_env("PROTEOMICS_FORCE_RERUN") ||
    truthy_env("PROTEOMICS_MAPTHATPROT_FORCE_RERUN") ||
    truthy_env("PROTEOMICS_RECOMPUTE") ||
    truthy_env("PROTEOMICS_MAPTHATPROT_RECOMPUTE") ||
    has_cli_flag(c("--force-rerun", "--recompute"))

cat("Target Comparison:", mapped_comparisons, "\n")
cat("Mapping Direction:", map_direction, "\n")
cat("Recompute existing tables:", force_rerun, "\n")

mapping_roots <- resolve_mapthatprot_roots()
mapping_paths <- resolve_mapthatprot_paths(mapped_comparisons, map_direction, mapping_roots)
CANONICAL_PATHS <- list(
    tables = mapping_paths$tables_root,
    logs = mapping_paths$logs_root,
    reports = mapping_paths$reports_root
)
raw_dir <- mapping_paths$raw_dir

# Define output directories for mapped datasets, unmapped trackers, and QC info
info_dir <- file.path(CANONICAL_PATHS$logs, mapped_comparisons, "mapped", map_direction, "info")
mapped_dir <- mapping_paths$mapped_dir
mapped_summary_dir <- file.path(CANONICAL_PATHS$tables, mapped_comparisons, "mapped", map_direction, "summaries")
unmapped_dir <- mapping_paths$unmapped_dir
unmapped_summary_dir <- file.path(CANONICAL_PATHS$tables, mapped_comparisons, "unmapped", map_direction, "summaries")
member_bridge_dir <- mapping_paths$member_bridge_dir
protein_group_audit_dir <- file.path(CANONICAL_PATHS$tables, mapped_comparisons, "protein_groups", map_direction, "audits")
accession_annotation_audit_dir <- file.path(CANONICAL_PATHS$tables, mapped_comparisons, "gene_annotation", map_direction, "accessions")
protein_group_annotation_audit_dir <- file.path(CANONICAL_PATHS$tables, mapped_comparisons, "gene_annotation", map_direction, "protein_groups")
report_dir <- file.path(CANONICAL_PATHS$reports, mapped_comparisons, "mapping_reports", map_direction)

# --- Reference Databases ---
# Define path for central UniProt species-specific knowledgebase flatfile
uniprot_mapping_file_path <- path_external("MOUSE_10090_idmapping.dat")

csv_files <- list_mapthatprot_input_files(raw_dir)
input_count_validation <- tryCatch(
    validate_mapthatprot_input_count(mapped_comparisons, map_direction, mapping_roots$gct_extract_root, csv_files),
    error = function(e) e
)
if (is_dry_run()) {
    expected_mapped <- file.path(mapped_dir, basename(csv_files))
    expected_unmapped <- file.path(unmapped_dir, basename(csv_files))
    expected_info <- file.path(
        info_dir,
        paste0(tools::file_path_sans_ext(basename(csv_files)), "_mapping_info.csv")
    )
    expected_mapped_summary <- file.path(
        mapped_summary_dir,
        paste0(tools::file_path_sans_ext(basename(csv_files)), "_summary.csv")
    )
    expected_unmapped_summary <- file.path(
        unmapped_summary_dir,
        paste0(tools::file_path_sans_ext(basename(csv_files)), "_summary.csv")
    )
    existing_complete <- file.exists(expected_mapped) &
        file.exists(expected_unmapped) &
        file.exists(expected_info) &
        file.exists(expected_mapped_summary) &
        file.exists(expected_unmapped_summary)
    resolution_report <- if (!inherits(input_count_validation, "error")) {
        mapthatprot_resolution_report(mapping_paths, input_count_validation)
    } else {
        c(
            "Resolved dataset" = mapped_comparisons,
            "Mapping direction" = map_direction,
            "GCT extract input root" = mapping_roots$gct_extract_root,
            "Mapping output root" = mapping_roots$mapping_output_root,
            "Results namespace" = mapping_paths$analysis_namespace
        )
    }
    dry_run_line("Script", "02_id_mapping/01_MapThatProt_batch.r")
    for (label in names(resolution_report)[seq_len(min(5L, length(resolution_report)))]) {
        dry_run_line(label, resolution_report[[label]])
    }
    dry_run_line("Raw contrast directory", raw_dir, if (dir.exists(raw_dir)) "PASS" else "FAIL")
    expected_count_label <- if (inherits(input_count_validation, "error") || is.na(input_count_validation$expected)) "at least 1" else input_count_validation$expected
    dry_run_line("Expected raw contrast CSV count", expected_count_label)
    dry_run_line("Raw contrast CSV count", length(csv_files), if (!inherits(input_count_validation, "error")) "PASS" else "FAIL")
    dry_run_line("Already mapped table sets", sum(existing_complete))
    dry_run_line("Remaining table sets", sum(!existing_complete))
    dry_run_line("Recompute existing tables", force_rerun)
    if (!dir.exists(raw_dir) || length(csv_files) == 0) {
        dry_run_line("Required upstream step", "Rscript 01_preprocessing/03_gct_extractR.r without --dry-run")
    }
    dry_run_line("UniProt mapping file", uniprot_mapping_file_path, if (file.exists(uniprot_mapping_file_path)) "PASS" else "FAIL")
    dry_run_line("Manual gene annotation overrides", Sys.getenv("PROTEOMICS_MANUAL_GENE_ANNOTATION_FILE", unset = path_metadata("manual_gene_annotation_overrides.csv")), "optional")
    dry_run_line("Mapped output directory", mapped_dir)
    dry_run_line("Unmapped output directory", unmapped_dir)
    dry_run_line("Member bridge output directory", member_bridge_dir)
    dry_run_line("Protein-group audit directory", protein_group_audit_dir)
    dry_run_line("Reports directory", report_dir)
    if (inherits(input_count_validation, "error")) dry_run_line("Input-count validation", conditionMessage(input_count_validation), "FAIL")
    quit(status = if (dir.exists(raw_dir) && !inherits(input_count_validation, "error") && file.exists(uniprot_mapping_file_path)) 0 else 1, save = "no")
}
if (!dir.exists(raw_dir)) stop("Raw contrast directory not found: ", raw_dir, call. = FALSE)
if (inherits(input_count_validation, "error")) stop(conditionMessage(input_count_validation), call. = FALSE)

# Initialize output folder structure only after all dry-run and input checks pass.
cat("Initializing output directories...\n")
for (directory in c(
    info_dir, mapped_dir, mapped_summary_dir, unmapped_dir, unmapped_summary_dir,
    member_bridge_dir, protein_group_audit_dir, accession_annotation_audit_dir,
    protein_group_annotation_audit_dir, report_dir
)) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

# Auto-download the UniProt mapping file if missing (useful for reproducibility)
if (!file.exists(uniprot_mapping_file_path)) {
    auto_download <- identical(tolower(Sys.getenv("AUTO_DOWNLOAD_UNIPROT_MAPPING", "false")), "true")
    if (!auto_download) {
        stop("UniProt mapping file not found: ", uniprot_mapping_file_path,
             ". Place MOUSE_10090_idmapping.dat under data/external or set AUTO_DOWNLOAD_UNIPROT_MAPPING=true.", call. = FALSE)
    }
    cat("UniProt mapping file not found at:", uniprot_mapping_file_path, "\nAttempting to download...\n")
    options(timeout = 3600)
    gz_url <- "https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/idmapping/by_organism/MOUSE_10090_idmapping.dat.gz"
    gz_file <- paste0(uniprot_mapping_file_path, ".gz")
    tryCatch({
        download.file(gz_url, gz_file, mode = "wb")
        R.utils::gunzip(gz_file, destname = uniprot_mapping_file_path, remove = TRUE)
        cat("Downloaded and unzipped the UniProt mapping file successfully.\n")
    }, error = function(e) stop("Failed to download/unzip UniProt mapping file: ", e$message))
}

cat("Checking and loading required libraries...\n")
load_required_packages(c("dplyr", "stringr", "tidyr", "purrr", "readr", "R.utils", "foreach", "doParallel", "readxl", "openxlsx", "AnnotationDbi", "org.Mm.eg.db"))
write_session_info(file.path(CANONICAL_PATHS$logs, "sessionInfo.txt"))
utils::write.csv(
    data.frame(
        comparison_family = mapped_comparisons,
        map_direction = map_direction,
        input_type = c("gct_extract_root", "mapping_output_root", "raw_contrast_directory", "uniprot_mapping", "manual_protein_mapping", "manual_gene_annotation"),
        path = c(mapping_roots$gct_extract_root, mapping_roots$mapping_output_root, raw_dir, uniprot_mapping_file_path,
            Sys.getenv("PROTEOMICS_MANUAL_MAPPING_FILE", unset = path_metadata("manual_mapping.xlsx")),
            Sys.getenv("PROTEOMICS_MANUAL_GENE_ANNOTATION_FILE", unset = path_metadata("manual_gene_annotation_overrides.csv"))),
        md5 = c(NA_character_, NA_character_, NA_character_, file_hash(uniprot_mapping_file_path),
            file_hash(Sys.getenv("PROTEOMICS_MANUAL_MAPPING_FILE", unset = path_metadata("manual_mapping.xlsx"))),
            file_hash(Sys.getenv("PROTEOMICS_MANUAL_GENE_ANNOTATION_FILE", unset = path_metadata("manual_gene_annotation_overrides.csv")))),
        stringsAsFactors = FALSE
    ),
    file.path(CANONICAL_PATHS$logs, mapped_comparisons, paste0("MapThatProt_batch_input_manifest_", map_direction, ".csv")),
    row.names = FALSE
)

# Parse the UniProt mapping dictionary natively into memory
cat("Parsing UniProt idmapping dictionary into memory... (This may take a moment)\n")
uniprot_mapping <- load_mouse_idmapping(uniprot_mapping_file_path)
mouse_maps <- build_mouse_maps(uniprot_mapping)
entry_map <- mouse_maps$entry_map
gene_map <- mouse_maps$gene_map
accession_gene_map <- mouse_maps$accession_gene_map
reviewed_map <- mouse_maps$reviewed_map
uniprot_mapping_file_hash <- file_hash(uniprot_mapping_file_path)
gene_annotation_maps <- build_mouse_gene_annotation_maps(
    org.Mm.eg.db,
    accessions = unique(accession_gene_map$UNIPROT)
)

cat("Found", length(csv_files), "CSV files to process in", raw_dir, "\n")

# --- Manual Override Configuration ---
# Optional manual curation file. This is crucial for resolving heavily ambiguous 
# protein groups or unannotated gene symbols common in exploratory proteomics.
cat("Checking for manual mapping override file...\n")
manual_mapping_path <- Sys.getenv("PROTEOMICS_MANUAL_MAPPING_FILE", unset = path_metadata("manual_mapping.xlsx"))
manual_override <- TRUE  # TRUE enforces curation over algorithmic mapping

manual_mapping <- read_manual_mapping_table(manual_mapping_path)
manual_gene_annotation_path <- Sys.getenv("PROTEOMICS_MANUAL_GENE_ANNOTATION_FILE", unset = path_metadata("manual_gene_annotation_overrides.csv"))
manual_gene_annotation_overrides <- read_manual_gene_annotation_overrides(manual_gene_annotation_path)
if (is.null(manual_mapping)) {
    cat("Notice: Manual mapping Excel not found at:", manual_mapping_path, "\n")
}

if (!is.null(manual_mapping)) {
    needed <- c("gene_symbol", "mapped_gene_symbol")
    missing_cols <- setdiff(needed, names(manual_mapping))
    if (length(missing_cols)) {
        cat("Warning: Manual mapping is missing expected columns:", paste(missing_cols, collapse = ", "), "\n")
    } else {
        cat("Successfully loaded", nrow(manual_mapping), "rows from manual mapping configuration.\n")
    }
}

load_existing_processed_result <- function(data_path, mapped_file, unmapped_file, info_table_file) {
    base <- tools::file_path_sans_ext(basename(data_path))
    cat("Skipping existing mapped tables for:", basename(data_path), "\n")

    mapping_info <- readr::read_csv(info_table_file, show_col_types = FALSE)
    if (!"comparison_family" %in% names(mapping_info)) mapping_info$comparison_family <- mapped_comparisons
    if (!"map_direction" %in% names(mapping_info)) mapping_info$map_direction <- map_direction

    mapping_table <- mapping_info %>%
        dplyr::mutate(
            comparison_family = mapped_comparisons,
            map_direction = map_direction,
            source_file = basename(data_path),
            mapped_status = ifelse(!is.na(final_accession) & nzchar(final_accession), "mapped", "unmapped")
        )

    unmapped_table <- readr::read_csv(unmapped_file, show_col_types = FALSE)
    if (!"gene_symbol" %in% names(unmapped_table)) {
        unmapped_table <- tibble::tibble(gene_symbol = character())
    }
    unmapped_table <- unmapped_table %>%
        dplyr::mutate(
            comparison_family = mapped_comparisons,
            map_direction = map_direction,
            source_file = basename(data_path)
        )

    invisible(list(
        mapping_table = mapping_table,
        unmapped_table = unmapped_table,
        multi_protein_log_table = tibble::tibble(
            source_file = character(),
            original_row_id = integer(),
            original_entry = character(),
            kept_protein = character(),
            dropped_proteins = character()
        ),
        manual_mapping_audit = tibble::tibble(
            comparison_family = character(),
            map_direction = character(),
            source_file = character(),
            original_token = character(),
            token_base = character(),
            resolved_uniprot = character(),
            mapping_strategy = character(),
            manual_mapping_used = logical()
        ),
        checkpoint_status = "skipped_existing",
        base = base
    ))
}

# --- Core Mapping Function ---
# This function is executed asynchronously for each input proteomics file.
process_file <- function(data_path) {
    base <- tools::file_path_sans_ext(basename(data_path))
    mapped_file <- file.path(mapped_dir, paste0(base, ".csv"))
    unmapped_file <- file.path(unmapped_dir, paste0(base, ".csv"))
    bridge_file <- file.path(member_bridge_dir, paste0(base, "_member_bridge.csv"))
    audit_file <- file.path(protein_group_audit_dir, paste0(base, "_protein_group_summary.csv"))
    collision_file <- file.path(protein_group_audit_dir, paste0(base, "_ProteinGroupID_collision_audit.csv"))
    accession_annotation_file <- file.path(accession_annotation_audit_dir, paste0(base, "_accession_gene_annotation_audit.csv"))
    protein_group_annotation_file <- file.path(protein_group_annotation_audit_dir, paste0(base, "_protein_group_gene_annotation_audit.csv"))
    info_table_file <- file.path(info_dir, paste0(base, "_mapping_info.csv"))
    info_summary_file <- file.path(info_dir, paste0(base, "_info.txt"))
    mapped_summary_file <- file.path(mapped_summary_dir, paste0(base, "_summary.csv"))
    unmapped_summary_file <- file.path(unmapped_summary_dir, paste0(base, "_summary.csv"))

    expected_tables <- c(mapped_file, unmapped_file, bridge_file, audit_file, collision_file,
        accession_annotation_file, protein_group_annotation_file,
        info_table_file, mapped_summary_file, unmapped_summary_file)
    if (!isTRUE(force_rerun) && all(file.exists(expected_tables))) {
        return(load_existing_processed_result(data_path, mapped_file, unmapped_file, info_table_file))
    }

    # Suppress internal messages for cleaner parallel console, file-level info handles writing
    # Ingest generic CSV arrays, supporting delimiter fallbacks (e.g. European CSVs)
    df_raw <- tryCatch(
        readr::read_csv(data_path, col_names = TRUE, show_col_types = FALSE, trim_ws = TRUE, quote = "\""),
        error = function(e) {
            readr::read_csv2(data_path, col_names = TRUE, show_col_types = FALSE, trim_ws = TRUE)
        }
    )

    # Coerce assumed key column to 'gene_symbol' representation
    if (!"gene_symbol" %in% names(df_raw)) {
        names(df_raw)[1] <- "gene_symbol"
    }

    canonical <- build_canonical_protein_group_tables(
        df_raw = df_raw,
        dataset = mapped_comparisons,
        source_file = data_path,
        entry_map = entry_map,
        gene_map = gene_map,
        accession_gene_map = accession_gene_map,
        reviewed_map = reviewed_map,
        manual_mapping = manual_mapping,
        manual_override = manual_override,
        gene_annotation_maps = gene_annotation_maps,
        manual_gene_annotation_overrides = manual_gene_annotation_overrides,
        uniprot_mapping_file_hash = uniprot_mapping_file_hash,
        strict = TRUE
    )

    df_mapped <- canonical$wide
    member_bridge <- canonical$bridge
    protein_group_summary <- canonical$summary
    collision_audit <- canonical$collision_audit
    accession_annotation_audit <- canonical$accession_annotation_audit
    protein_group_annotation_audit <- canonical$protein_group_annotation_audit

    mapping_info <- df_mapped %>%
        dplyr::transmute(
            comparison_family = mapped_comparisons,
            map_direction = map_direction,
            source_file = basename(data_path),
            original_row_id = .data$source_row_id,
            original_symbol = .data$original_identifier,
            base_name = .data$member_identifiers_canonical,
            ProteinGroupID = .data$ProteinGroupID,
            final_accession = .data$representative_accession,
            matched_by = .data$representative_selection_rule,
            manual_mapping_used = .data$ProteinGroupID %in% unique(member_bridge$ProteinGroupID[member_bridge$manual_mapping_used]),
            mapped_status = .data$mapping_status,
            protein_group_ambiguity_class = .data$protein_group_ambiguity_class,
            gene_level_claim_allowed = .data$gene_level_claim_allowed,
            protein_level_claim_allowed = .data$protein_level_claim_allowed,
            gene_symbol_compatibility_status = .data$gene_symbol_compatibility_status
        )

    manual_mapping_audit <- member_bridge %>%
        dplyr::filter(.data$manual_mapping_used) %>%
        dplyr::transmute(
            comparison_family = mapped_comparisons,
            map_direction = map_direction,
            source_file = basename(data_path),
            original_token = .data$member_identifier_original,
            token_base = .data$member_identifier_normalized,
            resolved_uniprot = .data$member_accession,
            mapping_strategy = .data$mapping_strategy,
            manual_mapping_used = TRUE
        )

    unmapped_proteins <- member_bridge %>%
        dplyr::filter(.data$member_mapping_status != "mapped") %>%
        dplyr::transmute(
            gene_symbol = .data$member_identifier_original,
            ProteinGroupID = .data$ProteinGroupID,
            original_row_id = .data$source_row_id,
            member_mapping_status = .data$member_mapping_status,
            mapping_strategy = .data$mapping_strategy
        ) %>%
        dplyr::distinct()

    multi_protein_log <- canonical$deprecated_dropped_log

    readr::write_csv(df_mapped, mapped_file)
    readr::write_csv(unmapped_proteins, unmapped_file)
    readr::write_csv(member_bridge, bridge_file)
    readr::write_csv(protein_group_summary, audit_file)
    readr::write_csv(collision_audit, collision_file)
    readr::write_csv(accession_annotation_audit, accession_annotation_file)
    readr::write_csv(protein_group_annotation_audit, protein_group_annotation_file)
    readr::write_csv(mapping_info, info_table_file)
    readr::write_tsv(manual_mapping_audit, file.path(info_dir, paste0(base, "_manual_mapping_audit.tsv")))

    mapped_summary <- df_mapped %>%
        dplyr::count(
            .data$protein_group_ambiguity_class,
            .data$mapping_status,
            .data$gene_level_claim_allowed,
            .data$protein_level_claim_allowed,
            name = "n"
        ) %>%
        dplyr::mutate(comparison_family = mapped_comparisons, map_direction = map_direction, .before = 1)
    readr::write_csv(mapped_summary, mapped_summary_file)

    unmapped_summary <- unmapped_proteins %>%
        dplyr::mutate(comparison_family = mapped_comparisons, map_direction = map_direction) %>%
        dplyr::group_by(.data$gene_symbol, .data$member_mapping_status) %>%
        dplyr::summarise(occurrences = dplyr::n(), .groups = "drop")
    readr::write_csv(unmapped_summary, unmapped_summary_file)

    summary_lines <- c(
        paste0("file: ", basename(data_path)),
        paste0("comparison_family: ", mapped_comparisons),
        paste0("map_direction: ", map_direction),
        paste0("total_quantitative_groups: ", protein_group_summary$total_quantitative_groups),
        paste0("mapped_groups: ", sum(df_mapped$mapping_status == "mapped")),
        paste0("partially_mapped_groups: ", protein_group_summary$partially_mapped_groups),
        paste0("unresolved_groups: ", protein_group_summary$unresolved_groups),
        paste0("ProteinGroupID_collisions: ", protein_group_summary$ProteinGroupID_collisions),
        paste0("groups_with_unstable_or_missing_source_feature_identifiers: ", protein_group_summary$groups_with_unstable_or_missing_source_feature_identifiers),
        "legacy gene_symbol is deprecated display-only; ProteinGroupID is the quantitative key."
    )
    writeLines(summary_lines, info_summary_file)

    invisible(list(
        mapping_table = mapping_info,
        unmapped_table = unmapped_proteins %>%
            dplyr::mutate(comparison_family = mapped_comparisons, map_direction = map_direction, source_file = basename(data_path)),
        multi_protein_log_table = multi_protein_log,
        manual_mapping_audit = manual_mapping_audit
    ))
}

# -------------------------
# Parallel Execution Engine
# -------------------------
# Distribute file-level mapping across available processor cores
n_files <- length(csv_files)
available_cores <- parallel::detectCores(logical = FALSE)
workers <- max(1, min(available_cores - 1, n_files))

cat("Setting up parallel processing backend using", workers, "out of", available_cores, "available cores...\n")
cl <- parallel::makeCluster(workers)
doParallel::registerDoParallel(cl)

cat("Initiating parallel mapping cascade for all files...\n")
results <- foreach(i = seq_along(csv_files),
                   .packages = c("dplyr", "stringr", "tidyr", "purrr", "readr", "R.utils", "tibble"),
                   .export = c(
                       "entry_map", "gene_map",
                       "accession_gene_map", "reviewed_map", "mapped_dir", "unmapped_dir",
                       "member_bridge_dir", "protein_group_audit_dir", "info_dir", "mapped_summary_dir",
                       "unmapped_summary_dir", "accession_annotation_audit_dir", "protein_group_annotation_audit_dir",
                       "gene_annotation_maps", "manual_gene_annotation_overrides", "uniprot_mapping_file_hash", "mapped_comparisons",
                       "map_direction", "force_rerun", "manual_mapping", "manual_override",
                       "load_existing_processed_result", "process_file",
                       "normalize_token", "to_base_no_iso_mouse", "is_uniprot_ac", "extract_ac", "extract_entry",
                       "split_protein_group_members", "normalize_member_identifier", "canonical_member_set",
                       "detect_source_feature_columns", "stable_pg_hash", "parse_member_identifier",
                       "lookup_accession_gene", "normalize_manual_gene_annotation_overrides",
                       "resolve_mouse_gene_annotation", "assess_protein_group_gene_annotation",
                       "mouse_gene_annotation_contract_version", "resolve_protein_group_member", "classify_protein_group",
                       "protein_group_claim_rules", "select_representative_member",
                       "validate_protein_group_id_collisions", "build_canonical_protein_group_tables"
                   )) %dopar% {
    process_file(csv_files[i])
}

parallel::stopCluster(cl)
cat("Batch parallel mapping successfully completed for", n_files, "files.\n")

# -------------------------
# Global Mapping Workbooks
# -------------------------

cat("Aggregating overall biology summaries and computing mapping strategy statistics...\n")

# Aggregate outputs back from processing clusters
all_mapping_tables <- purrr::map(results, "mapping_table") %>% dplyr::bind_rows()
all_unmapped_tables <- purrr::map(results, "unmapped_table") %>% dplyr::bind_rows()
all_dropped_proteins <- purrr::map(results, "multi_protein_log_table") %>% dplyr::bind_rows()
all_manual_mapping_audit <- purrr::map(results, "manual_mapping_audit") %>% dplyr::bind_rows()


# Consolidate standard mappings
global_mapping_summary <- all_mapping_tables %>%
    dplyr::distinct() %>%
    dplyr::arrange(source_file, original_symbol)

global_mapping_file <- file.path(info_dir, "GLOBAL_mapping_summary.csv")
readr::write_csv(global_mapping_summary, global_mapping_file)

# Backup in summaries dir
global_mapping_summary_file2 <- file.path(mapped_summary_dir, "GLOBAL_mapping_summary.csv")
readr::write_csv(global_mapping_summary, global_mapping_summary_file2)

cat("Saved global mapping summary to:", info_dir, "and", mapped_summary_dir, "\n")

readr::write_tsv(all_manual_mapping_audit, file.path(info_dir, "manual_mapping_audit.tsv"))

# Consolidate ID unmapped dropouts
global_unmapped_summary <- all_unmapped_tables %>%
    dplyr::group_by(comparison_family, map_direction, gene_symbol) %>%
    dplyr::summarise(
        occurrences = dplyr::n(),
        files_present = paste(unique(source_file), collapse = "; "),
        .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(occurrences))

global_unmapped_file <- file.path(unmapped_dir, "GLOBAL_unmapped_proteins.csv")
readr::write_csv(global_unmapped_summary, global_unmapped_file)

global_unmapped_summary_file2 <- file.path(unmapped_summary_dir, "GLOBAL_unmapped_proteins.csv")
readr::write_csv(global_unmapped_summary, global_unmapped_summary_file2)

cat("Saved global unmapped tracking to:", unmapped_dir, "and", unmapped_summary_dir, "\n")

# -------------------------
# Strategy Quality Control
# -------------------------

strategy_stats <- global_mapping_summary %>%
    dplyr::filter(mapped_status == "mapped") %>%
    dplyr::count(matched_by, name = "total_mapped") %>%
    dplyr::arrange(dplyr::desc(total_mapped))

strategy_file <- file.path(info_dir, "GLOBAL_strategy_statistics.csv")
readr::write_csv(strategy_stats, strategy_file)
cat("Saved global strategy statistics to:", strategy_file, "\n")

mapping_full <- all_mapping_tables %>%
    dplyr::distinct() %>%
    dplyr::arrange(source_file, original_symbol)

# Produce unique biological entity tracking
protein_summary <- mapping_full %>%
    dplyr::group_by(comparison_family, map_direction, original_symbol) %>%
    dplyr::summarise(
        mapped_to = paste(unique(na.omit(final_accession)), collapse = "; "),
        strategies = paste(unique(na.omit(matched_by)), collapse = "; "),
        files_present = paste(unique(source_file), collapse = "; "),
        mapped = any(!is.na(final_accession)),
        .groups = "drop"
    )

unmapped_proteins_global <- protein_summary %>%
    dplyr::filter(mapped == FALSE)


# Tally relative efficiency across mapping resolution methods
strategy_stats <- mapping_full %>%
    dplyr::filter(!is.na(final_accession)) %>%
    dplyr::count(matched_by, name = "count") %>%
    dplyr::arrange(desc(count))

strategy_stats_unique <- mapping_full %>%
    dplyr::filter(!is.na(final_accession)) %>%
    dplyr::distinct(matched_by, original_symbol) %>%
    dplyr::count(matched_by, name = "unique_proteins") %>%
    dplyr::arrange(desc(unique_proteins))

# Calculate dataset coverage % 
coverage_stats <- mapping_full %>%
    dplyr::group_by(comparison_family, map_direction, source_file) %>%
    dplyr::summarise(
        total = dplyr::n(),
        mapped = sum(!is.na(final_accession)),
        coverage = mapped / total,
        .groups = "drop"
    )

# --- Generate Comprehensive Excel QC Document ---
cat("Generating comprehensive Excel QC Report...\n")

report_file <- file.path(report_dir, "Mapping_QC_Report.xlsx")

wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "Full_Mapping_Table")
openxlsx::writeData(wb, "Full_Mapping_Table", mapping_full)

openxlsx::addWorksheet(wb, "Protein_Summary")
openxlsx::writeData(wb, "Protein_Summary", protein_summary)

openxlsx::addWorksheet(wb, "Unmapped_Proteins")
openxlsx::writeData(wb, "Unmapped_Proteins", unmapped_proteins_global)

openxlsx::addWorksheet(wb, "Dropped_Proteins")
openxlsx::writeData(wb, "Dropped_Proteins", all_dropped_proteins)

openxlsx::addWorksheet(wb, "Strategy_Stats")
openxlsx::writeData(wb, "Strategy_Stats", strategy_stats)

openxlsx::addWorksheet(wb, "Unique_Strategy_Stats")
openxlsx::writeData(wb, "Unique_Strategy_Stats", strategy_stats_unique)

openxlsx::addWorksheet(wb, "Coverage_Stats")
openxlsx::writeData(wb, "Coverage_Stats", coverage_stats)

openxlsx::saveWorkbook(wb, report_file, overwrite = TRUE)

cat("Saved master QC workbook to:", report_file, "\n")
cat("====================================================\n")
cat("MapThatProt_batch execution completed successfully!\n")
cat("====================================================\n")

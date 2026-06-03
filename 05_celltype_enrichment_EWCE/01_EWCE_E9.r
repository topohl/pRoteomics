# ==========================================
# Consumes:
#   - processed proteomics matrix from data/processed/01_preprocessing/ or local config override
#   - sample metadata from data/metadata/ or local config override
# Produces:
#   - EWCE tables, figures, source data, processed objects, cache and logs under canonical module folders
# File contract:
#   - docs/active_script_io_audit.tsv object 05_celltype_enrichment_EWCE/01_EWCE_E9.r
# ==========================================
# EWCE E9: Publication-oriented proteomics workflow
# ==========================================

# This script keeps the original EWCE analysis intent, but adds:
# - sample-level limma contrasts rather than mean-only differences
# - measured-proteome background for EWCE
# - up/down differential signatures
# - hit-list size and annotation-level sensitivity analyses
# - per-target and global FDR
# - panel-specific source data and reproducibility logs

options(repos = c(CRAN = "https://cloud.r-project.org"))
set.seed(42)

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "dataset_inputs.R"))
source(repo_path("R", "validation_utils.R"))
MODULE_ID <- "05_celltype_enrichment_EWCE"

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = "") {
  hit <- which(args == flag)
  if (!length(hit) || hit[1] == length(args)) return(default)
  args[[hit[1] + 1]]
}

dataset_cli <- arg_value("--dataset", default = "")
if (nzchar(dataset_cli)) {
  Sys.setenv(PROTEOMICS_DATASET = validate_dataset(dataset_cli, source = "--dataset"))
}

infer_dataset_from_path <- function(path) {
  if (is.null(path) || !length(path) || is.na(path) || !nzchar(path)) return(NA_character_)
  normalized_path <- normalize_dataset(path)
  hits <- valid_datasets()[vapply(
    valid_datasets(),
    function(dataset) grepl(dataset, normalized_path, fixed = TRUE),
    logical(1)
  )]
  if (length(hits) == 1L) hits else NA_character_
}

dataset_env_names <- c("PROTEOMICS_DATASET", "PROTEOMICS_COMPARISON", "PROTEOMICS_GCT_COMPARISON")
has_dataset_override <- any(nzchar(Sys.getenv(dataset_env_names, unset = "")))
if (!has_dataset_override) {
  explicit_matrix_paths <- c(
    Sys.getenv("PROTEOMICS_EWCE_MATRIX", unset = ""),
    Sys.getenv("PROTEOMICS_WGCNA_UPSTREAM_XLSX", unset = "")
  )
  inferred_dataset <- unique(stats::na.omit(vapply(explicit_matrix_paths, infer_dataset_from_path, character(1))))
  if (length(inferred_dataset) == 1L) {
    Sys.setenv(PROTEOMICS_DATASET = validate_dataset(inferred_dataset, source = "input path"))
  }
}

EWCE_DATASET <- current_dataset()
SUBSTEP_ID <- file.path("EWCE_E9", EWCE_DATASET)
CANONICAL_PATHS <- create_module_dirs(MODULE_ID, SUBSTEP_ID)

resolve_ewce_pg_matrix <- function() {
  env_path <- Sys.getenv("PROTEOMICS_EWCE_MATRIX", unset = "")
  if (nzchar(env_path)) {
    return(normalizePath(env_path, winslash = "/", mustWork = FALSE))
  }
  dataset_inputs <- resolve_dataset_inputs(EWCE_DATASET, purpose = "wgcna")
  if (!is.na(dataset_inputs$expression_file) && file.exists(dataset_inputs$expression_file)) {
    return(dataset_inputs$expression_file)
  }

  imputation_qc_path <- path_processed("01_preprocessing", "impute", "imputation_qc.csv")
  if (file.exists(imputation_qc_path)) {
    imputation_qc <- utils::read.csv(imputation_qc_path, stringsAsFactors = FALSE, check.names = FALSE)
    required_cols <- c("celltype_layer", "n_samples", "output_path")
    if (all(required_cols %in% colnames(imputation_qc))) {
      candidates <- imputation_qc[
        imputation_qc$celltype_layer == EWCE_DATASET &
          !is.na(imputation_qc$output_path) &
          nzchar(imputation_qc$output_path),
        ,
        drop = FALSE
      ]
      candidates <- candidates[file.exists(candidates$output_path), , drop = FALSE]
      if (nrow(candidates) > 0) {
        candidates$mtime <- file.info(candidates$output_path)$mtime
        candidates <- candidates[order(candidates$mtime, decreasing = TRUE), , drop = FALSE]
        return(normalizePath(candidates$output_path[1], winslash = "/", mustWork = FALSE))
      }
    }
  }

  legacy_path <- path_processed("01_preprocessing", paste0("20260218_pgmatrix_imputed_", EWCE_DATASET, "_180samples_missing70pct.xlsx"))
  if (file.exists(legacy_path)) {
    return(normalizePath(legacy_path, winslash = "/", mustWork = FALSE))
  }

  discovered <- list.files(
    path_processed("01_preprocessing", "impute"),
    pattern = paste0("pgmatrix_imputed_", EWCE_DATASET, "_[0-9]+samples_missing70pct\\.xlsx$"),
    full.names = TRUE,
    recursive = FALSE
  )
  discovered <- discovered[file.exists(discovered)]
  if (length(discovered) > 0) {
    discovered <- discovered[order(file.info(discovered)$mtime, decreasing = TRUE)]
    return(normalizePath(discovered[1], winslash = "/", mustWork = FALSE))
  }

  legacy_path
}

# Local machines can override these by defining the variables before source()
# or through PROTEOMICS_PROJECT_ROOT / PROTEOMICS_EWCE_MATRIX. Do not commit machine-specific paths.
if (!exists("data_path", inherits = FALSE)) {
  data_path <- resolve_ewce_pg_matrix()
}
if (!exists("sample_metadata_path", inherits = FALSE)) {
  sample_metadata_path <- resolve_dataset_inputs(EWCE_DATASET, purpose = "wgcna")$metadata_file
}
base_results <- CANONICAL_PATHS$reports
dirs <- list(
  plots   = CANONICAL_PATHS$figures,
  svgs    = file.path(CANONICAL_PATHS$figures, "SVG_Editable"),
  tables  = CANONICAL_PATHS$tables,
  qc      = CANONICAL_PATHS$logs,
  data    = CANONICAL_PATHS$processed,
  source  = CANONICAL_PATHS$source_data,
  cache   = file.path(CANONICAL_PATHS$processed, "cache")
)

if (is_dry_run()) {
  sample_match <- "not checked"
  sample_match_status <- "WARN"
  if (file.exists(data_path) && file.exists(sample_metadata_path) && requireNamespace("readxl", quietly = TRUE)) {
    sample_match <- tryCatch({
      expr_head <- readxl::read_excel(data_path, n_max = 3)
      md_head <- readxl::read_excel(sample_metadata_path, n_max = 1000)
      sample_col <- detect_column(md_head, c("sample_id", "SampleID", "SampleColumn", "row.names", "Sample"))
      expr_samples <- setdiff(names(expr_head), c("gene_symbol", "T: Protein.Names", "Genes", "Protein.Group"))
      if (is.na(sample_col)) {
        sample_match_status <<- "FAIL"
        "metadata sample column not detected"
      } else {
        s <- sample_overlap_summary(expr_samples, md_head[[sample_col]])
        sample_match_status <<- if (s$n_overlap > 0) "PASS" else "FAIL"
        paste0(s$n_overlap, " matching sample IDs")
      }
    }, error = function(e) {
      sample_match_status <<- "WARN"
      conditionMessage(e)
    })
  }
  dry_run_line("Script", "05_celltype_enrichment_EWCE/01_EWCE_E9.r")
  dry_run_line("Dataset", EWCE_DATASET)
  dry_run_line(
    "Dataset source",
    if (nzchar(dataset_cli)) "--dataset" else if (has_dataset_override) "environment" else "input path/default"
  )
  dry_run_line("Proteomics matrix", data_path, if (file.exists(data_path)) "PASS" else "FAIL")
  dry_run_line("Sample metadata", sample_metadata_path, if (file.exists(sample_metadata_path)) "PASS" else "FAIL")
  dry_run_line("Cheap sample matching", sample_match, sample_match_status)
  dry_run_line("Output folders", paste(unlist(dirs), collapse = "; "))
  quit(status = if (file.exists(data_path) && file.exists(sample_metadata_path)) 0 else 1, save = "no")
}

if (!file.exists(data_path)) stop("Proteomics matrix not found: ", data_path, call. = FALSE)
if (!file.exists(sample_metadata_path)) stop("Sample metadata not found: ", sample_metadata_path, call. = FALSE)

cran_packages <- c(
  "readxl", "dplyr", "tibble", "tidyr", "ggplot2", "pheatmap",
  "svglite", "ggrepel", "ggsci", "viridis",
  "openxlsx", "stringr", "patchwork", "future", "future.apply",
  "digest"
)

bioc_packages <- c(
  "limma", "EWCE", "ewceData", "org.Mm.eg.db", "AnnotationDbi"
)

load_required <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing required R package: ", pkg, ". Install it explicitly before running this script.", call. = FALSE)
  }
  library(pkg, character.only = TRUE)
}

invisible(lapply(c(cran_packages, bioc_packages), load_required))

# ==========================================
# 1. CONFIG
# ==========================================

analysis_params <- list(
  seed = 42,
  reps = 10000,
  annot_levels = c(1, 2),
  primary_annot_level = 2,
  top_n_values = c(100, 250, 500),
  primary_top_n = 250,
  fdr_alpha = 0.05,
  marker_top_n = 200,
  regions = c("CA1", "CA2", "CA3", "DG"),
  layers = c("so", "sp", "sr", "slm", "mo", "po", "sg"),
  conditions = c("con", "res", "sus"),
  expgroup_condition_map = c("1" = "con", "2" = "res", "3" = "sus")
)

lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
write_session_info(file.path(dirs$qc, "sessionInfo.txt"))

message("Step 1: Loading CTD and proteomics matrix...")
ctd <- ewceData::ctd()
ref_genes_by_level <- lapply(ctd, function(x) rownames(x$specificity))
ref_genes <- unique(unlist(ref_genes_by_level, use.names = FALSE))

raw_df <- readxl::read_excel(data_path)
input_hash <- if (file.exists(data_path)) digest::digest(file = data_path, algo = "sha256") else NA_character_
required_raw_cols <- c("Genes")
missing_raw_cols <- setdiff(required_raw_cols, colnames(raw_df))
if (length(missing_raw_cols) > 0) {
  stop("EWCE input matrix missing required columns: ", paste(missing_raw_cols, collapse = ", "), call. = FALSE)
}

# ==========================================
# 2. HELPERS
# ==========================================

theme_publication <- function() {
  ggplot2::theme_bw(base_size = 7) +
    ggplot2::theme(
      text = ggplot2::element_text(family = "sans", color = "black"),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.3),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(linewidth = 0.2, color = "black"),
      axis.text = ggplot2::element_text(size = 6, color = "black"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 7),
      legend.key.size = grid::unit(3, "mm"),
      legend.text = ggplot2::element_text(size = 5),
      legend.title = ggplot2::element_text(size = 6, face = "bold"),
      plot.title = ggplot2::element_text(size = 8, face = "bold", hjust = 0)
    )
}

safe_sheet <- function(x) {
  x <- stringr::str_replace_all(x, "[^A-Za-z0-9_]", "_")
  substr(x, 1, 31)
}

safe_file_stem <- function(x) {
  x <- stringr::str_replace_all(x, "[^A-Za-z0-9_.-]", "_")
  substr(x, 1, 180)
}

clean_gene_token <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_split_i(";", 1) %>%
    stringr::str_trim()
}

extract_sample_token <- function(sample_names, tokens, case = c("asis", "upper", "lower")) {
  case <- match.arg(case)
  pattern <- paste0(
    "(?i)(^|[^A-Za-z0-9])(",
    paste(tokens, collapse = "|"),
    ")(?=$|[^A-Za-z0-9])"
  )
  out <- stringr::str_match(sample_names, pattern)[, 3]
  if (case == "upper") out <- stringr::str_to_upper(out)
  if (case == "lower") out <- stringr::str_to_lower(out)
  out
}

normalize_sample_id <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_replace_all("\\\\", "/") %>%
    stringr::str_trim()
}

normalize_condition <- function(x) {
  x <- stringr::str_to_lower(as.character(x))
  mapped <- unname(analysis_params$expgroup_condition_map[x])
  dplyr::if_else(!is.na(mapped), mapped, x)
}

load_sample_condition_lookup <- function(metadata_path) {
  if (is.null(metadata_path) || !file.exists(metadata_path)) {
    return(tibble::tibble(Sample = character(), AnimalID = character(), Cond_Metadata = character()))
  }

  meta_df <- readxl::read_excel(metadata_path, sheet = readxl::excel_sheets(metadata_path)[1]) %>%
    tibble::as_tibble()

  sample_col <- intersect(c("sample_id", "Sample", "sample", "SampleID"), colnames(meta_df))[1]
  animal_col <- intersect(c("AnimalID", "animal_id", "Animal", "animal"), colnames(meta_df))[1]
  cond_col <- intersect(c("Cond", "Condition", "condition", "ExpGroup", "Group", "group"), colnames(meta_df))[1]

  if (is.na(cond_col) || (is.na(sample_col) && is.na(animal_col))) {
    warning("Sample metadata file found but no usable sample/animal condition columns were detected: ", metadata_path)
    return(tibble::tibble(Sample = character(), AnimalID = character(), Cond_Metadata = character()))
  }

  meta_df %>%
    dplyr::transmute(
      Sample = if (!is.na(sample_col)) normalize_sample_id(.data[[sample_col]]) else NA_character_,
      AnimalID = if (!is.na(animal_col)) as.character(.data[[animal_col]]) else stringr::str_match(Sample, "(?i)(?:^|_)(A[0-9]+)(?=_)")[, 2],
      Cond_Metadata = normalize_condition(.data[[cond_col]])
    ) %>%
    dplyr::filter(!is.na(Cond_Metadata), Cond_Metadata %in% analysis_params$conditions) %>%
    dplyr::distinct()
}

parse_sample_metadata <- function(sample_names, condition_lookup = NULL) {
  direct_meta <- tibble::tibble(
    Sample = sample_names,
    SampleKey = normalize_sample_id(sample_names),
    AnimalID = stringr::str_match(sample_names, "(?i)(?:^|_)(A[0-9]+)(?=_)")[, 2]
  )

  if (!is.null(condition_lookup) && nrow(condition_lookup) > 0) {
    by_sample <- condition_lookup %>%
      dplyr::filter(!is.na(Sample)) %>%
      dplyr::mutate(SampleKey = normalize_sample_id(Sample)) %>%
      dplyr::distinct(SampleKey, Cond_Metadata)

    by_animal <- condition_lookup %>%
      dplyr::filter(!is.na(AnimalID)) %>%
      dplyr::distinct(AnimalID, Cond_Animal = Cond_Metadata)

    direct_meta <- direct_meta %>%
      dplyr::left_join(by_sample, by = "SampleKey") %>%
      dplyr::left_join(by_animal, by = "AnimalID")
  } else {
    direct_meta <- direct_meta %>%
      dplyr::mutate(Cond_Metadata = NA_character_, Cond_Animal = NA_character_)
  }

  resolved_meta <- direct_meta %>%
    dplyr::select(Sample, AnimalID, Cond_Metadata, Cond_Animal) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      Cond_From_Name = extract_sample_token(Sample, analysis_params$conditions, case = "lower"),
      Cond_From_Name = dplyr::if_else(
        is.na(Cond_From_Name),
        normalize_condition(stringr::str_extract(stringr::str_to_lower(Sample), paste(analysis_params$conditions, collapse = "|"))),
        Cond_From_Name
      ),
      Cond_Resolved = dplyr::coalesce(Cond_From_Name, Cond_Metadata, Cond_Animal)
    ) %>%
    dplyr::select(Sample, AnimalID, Cond_Resolved) %>%
    dplyr::right_join(tibble::tibble(Sample = sample_names), by = "Sample")

  resolved_meta %>%
    dplyr::mutate(
      Region = extract_sample_token(Sample, analysis_params$regions, case = "upper"),
      Layer = extract_sample_token(Sample, analysis_params$layers, case = "lower"),
      Stratum = dplyr::if_else(is.na(Layer), Region, paste(Region, Layer, sep = "_")),
      Cond = factor(Cond_Resolved, levels = analysis_params$conditions),
      Batch = stringr::str_extract(Sample, "(?i)(batch|plate|run)[-_]?[A-Za-z0-9]+")
    )
}

make_expr_mat <- function(df, sample_cols) {
  df %>%
    dplyr::select(Gene, dplyr::all_of(sample_cols)) %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(sample_cols), ~ mean(as.numeric(.x), na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    tibble::column_to_rownames("Gene") %>%
    as.matrix()
}

run_limma_stratum <- function(expr_mat, sample_meta, stratum) {
  stratum_meta <- sample_meta %>%
    dplyr::filter(Stratum == stratum, !is.na(Cond))
  stratum_samples <- stratum_meta %>%
    dplyr::pull(Sample)

  if (length(stratum_samples) < 4) {
    return(tibble::tibble())
  }

  meta <- sample_meta %>%
    dplyr::filter(Sample %in% stratum_samples) %>%
    dplyr::arrange(match(Sample, stratum_samples))

  x <- expr_mat[, meta$Sample, drop = FALSE]
  keep <- rowSums(is.finite(x)) >= 2
  x <- x[keep, , drop = FALSE]
  x[!is.finite(x)] <- NA

  if (anyNA(x)) {
    row_medians <- apply(x, 1, stats::median, na.rm = TRUE)
    idx <- which(is.na(x), arr.ind = TRUE)
    x[idx] <- row_medians[idx[, 1]]
  }

  meta$Cond <- droplevels(factor(meta$Cond, levels = analysis_params$conditions))
  present_conditions <- levels(meta$Cond)[levels(meta$Cond) %in% meta$Cond]

  if (!"con" %in% present_conditions || length(present_conditions) < 2) {
    return(tibble::tibble())
  }

  design <- stats::model.matrix(~ 0 + Cond, data = meta)
  colnames(design) <- sub("^Cond", "", colnames(design))

  contrast_names <- c()
  if ("sus" %in% colnames(design)) contrast_names <- c(contrast_names, "sus - con")
  if ("res" %in% colnames(design)) contrast_names <- c(contrast_names, "res - con")
  if (length(contrast_names) == 0) return(tibble::tibble())

  contrast_matrix <- limma::makeContrasts(contrasts = contrast_names, levels = design)
  colnames(contrast_matrix) <- c(
    if ("sus - con" %in% contrast_names) "Sus_vs_Con",
    if ("res - con" %in% contrast_names) "Res_vs_Con"
  )

  fit <- limma::lmFit(x, design)
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, contrast_matrix), trend = TRUE, robust = TRUE)

  dplyr::bind_rows(lapply(colnames(contrast_matrix), function(coef_name) {
    limma::topTable(fit2, coef = coef_name, number = Inf, sort.by = "none") %>%
      tibble::rownames_to_column("Gene") %>%
      tibble::as_tibble() %>%
      dplyr::mutate(
        Stratum = stratum,
        Region = dplyr::first(meta$Region),
        Layer = dplyr::first(meta$Layer),
        Contrast = coef_name,
        Direction_for_EWCE = dplyr::case_when(
          logFC > 0 ~ "up",
          logFC < 0 ~ "down",
          TRUE ~ "flat"
        ),
        RankingStat = t
      )
  }))
}

make_baseline_targets <- function(baseline_mat, top_n_values) {
  dplyr::bind_rows(lapply(colnames(baseline_mat), function(target) {
    vals <- baseline_mat[, target]
    ord <- names(sort(vals, decreasing = TRUE, na.last = NA))
    target_parts <- stringr::str_split_fixed(target, "__", 2)
    stratum <- target_parts[, 1]
    metric <- target_parts[, 2]
    region <- stringr::str_extract(stratum, paste(analysis_params$regions, collapse = "|"))
    layer <- stringr::str_match(stratum, paste0("_(?:", paste(analysis_params$layers, collapse = "|"), ")$"))[, 1]
    layer <- stringr::str_remove(layer, "^_")
    dplyr::bind_rows(lapply(top_n_values, function(top_n) {
      tibble::tibble(
        Target = target,
        AnalysisType = "Baseline",
        Stratum = stratum,
        Region = region,
        Layer = layer,
        Contrast = NA_character_,
        Direction = "abundant",
        Metric = metric,
        TopN = top_n,
        Gene = head(ord, top_n),
        Rank = seq_along(head(ord, top_n)),
        RankValue = vals[head(ord, top_n)]
      )
    }))
  }))
}

make_differential_targets <- function(de_tbl, top_n_values) {
  dplyr::bind_rows(lapply(split(de_tbl, list(de_tbl$Stratum, de_tbl$Contrast), drop = TRUE), function(tbl) {
    if (nrow(tbl) == 0) return(tibble::tibble())
    stratum <- tbl$Stratum[1]
    region <- tbl$Region[1]
    layer <- tbl$Layer[1]
    contrast <- tbl$Contrast[1]

    dplyr::bind_rows(lapply(c("up", "down"), function(direction) {
      ranked <- if (direction == "up") {
        tbl %>% dplyr::arrange(dplyr::desc(RankingStat), dplyr::desc(logFC))
      } else {
        tbl %>% dplyr::arrange(RankingStat, logFC)
      }

      dplyr::bind_rows(lapply(top_n_values, function(top_n) {
        take <- ranked %>% dplyr::slice_head(n = top_n)
        tibble::tibble(
          Target = paste(stratum, contrast, direction, paste0("top", top_n), sep = "_"),
          AnalysisType = "Differential",
          Stratum = stratum,
          Region = region,
          Layer = layer,
          Contrast = contrast,
          Direction = direction,
          Metric = paste(contrast, direction, sep = "_"),
          TopN = top_n,
          Gene = take$Gene,
          Rank = seq_len(nrow(take)),
          RankValue = take$RankingStat,
          logFC = take$logFC,
          P.Value = take$P.Value,
          adj.P.Val = take$adj.P.Val
        )
      }))
    }))
  }))
}

run_ewce_once <- function(hits, bg, annot_level) {
  hits <- unique(stats::na.omit(hits))
  bg <- unique(stats::na.omit(bg))
  hits <- intersect(hits, bg)

  if (length(hits) < 10) {
    stop("EWCE hit list has fewer than 10 genes after background intersection.")
  }

  res <- try(
    EWCE::bootstrap_enrichment_test(
      sct_data = ctd,
      hits = hits,
      bg = bg,
      reps = analysis_params$reps,
      annotLevel = annot_level,
      genelistSpecies = "mouse",
      sctSpecies = "mouse"
    ),
    silent = TRUE
  )

  if (inherits(res, "try-error") && exists("bootstrap.enrichment.test", envir = asNamespace("EWCE"))) {
    bootstrap_enrichment_test_legacy <- get("bootstrap.enrichment.test", envir = asNamespace("EWCE"))
    res <- try(
      bootstrap_enrichment_test_legacy(
        sct_data = ctd[[annot_level]],
        mouse.hits = hits,
        mouse.bg = bg,
        reps = analysis_params$reps
      ),
      silent = TRUE
    )
  }

  if (inherits(res, "try-error")) {
    res <- try(
      EWCE::bootstrap_enrichment_test(
        sct_data = ctd,
        hits = hits,
        reps = analysis_params$reps,
        annotLevel = annot_level,
        genelistSpecies = "mouse",
        sctSpecies = "mouse"
      ),
      silent = TRUE
    )
  }

  if (inherits(res, "try-error")) stop(res)
  res$results
}

add_worksheet_safe <- function(wb, sheet, data) {
  sheet <- safe_sheet(sheet)
  if (sheet %in% names(wb)) {
    sheet <- safe_sheet(paste0(sheet, "_", length(names(wb)) + 1))
  }
  data <- as.data.frame(data, check.names = FALSE)
  col_names <- names(data)
  col_names <- as.character(col_names)
  col_names[is.na(col_names) | stringr::str_trim(col_names) == ""] <- "Column"
  col_names <- stringr::str_squish(col_names)
  col_names <- make.unique(col_names, sep = "_")
  names(data) <- col_names
  openxlsx::addWorksheet(wb, sheet)
  tryCatch(
    openxlsx::writeDataTable(wb, sheet, data),
    error = function(e) {
      if (!grepl("colNames must be a unique vector", conditionMessage(e), fixed = TRUE)) {
        stop(e)
      }
      openxlsx::writeData(wb, sheet, data, withFilter = TRUE)
    }
  )
}

celltype_marker_overlap <- function(results_tbl, target_gene_tbl, annot_level) {
  specificity <- ctd[[annot_level]]$specificity
  if (is.null(specificity) || nrow(specificity) == 0) return(tibble::tibble())

  marker_tbl <- dplyr::bind_rows(lapply(colnames(specificity), function(ct) {
    tibble::tibble(
      CellType = ct,
      MarkerGene = names(sort(specificity[, ct], decreasing = TRUE, na.last = NA))[
        seq_len(min(analysis_params$marker_top_n, nrow(specificity)))
      ]
    )
  }))

  results_tbl %>%
    dplyr::filter(
      AnnotLevel == annot_level,
      TopN == analysis_params$primary_top_n,
      q_global < analysis_params$fdr_alpha
    ) %>%
    dplyr::select(Target, CellType) %>%
    dplyr::distinct() %>%
    dplyr::left_join(
      target_gene_tbl %>%
        dplyr::filter(TopN == analysis_params$primary_top_n) %>%
        dplyr::distinct(Target, Gene),
      by = "Target",
      relationship = "many-to-many"
    ) %>%
    dplyr::inner_join(marker_tbl, by = c("CellType", "Gene" = "MarkerGene")) %>%
    dplyr::group_by(Target, CellType) %>%
    dplyr::summarise(
      N_Overlap_TopMarkers = dplyr::n_distinct(Gene),
      DriverGenes = paste(sort(unique(Gene)), collapse = ";"),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(N_Overlap_TopMarkers), Target, CellType)
}

# ==========================================
# 3. GENE MAPPING AND EXPRESSION MATRIX
# ==========================================

message("Step 2: Robust gene mapping and QC logs...")
sample_cols <- colnames(raw_df)[grep("^D:", colnames(raw_df))]
if (length(sample_cols) == 0) {
  stop("No sample columns detected. Expected sample columns starting with 'D:'.")
}

mapping_input <- raw_df %>%
  dplyr::mutate(
    RowID = dplyr::row_number(),
    Gene_Input = clean_gene_token(Genes)
  ) %>%
  dplyr::filter(!is.na(Gene_Input), Gene_Input != "0", Gene_Input != "")

mapping_selected <- AnnotationDbi::mapIds(
  org.Mm.eg.db,
  keys = mapping_input$Gene_Input,
  column = "SYMBOL",
  keytype = "ALIAS",
  multiVals = "first"
)

mapping_all <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = unique(mapping_input$Gene_Input),
  columns = "SYMBOL",
  keytype = "ALIAS"
) %>%
  tibble::as_tibble() %>%
  dplyr::group_by(ALIAS) %>%
  dplyr::summarise(
    All_SYMBOL_Matches = paste(sort(unique(stats::na.omit(SYMBOL))), collapse = ";"),
    N_SYMBOL_Matches = dplyr::n_distinct(stats::na.omit(SYMBOL)),
    .groups = "drop"
  )

clean_df <- mapping_input %>%
  dplyr::mutate(Gene = as.character(mapping_selected[Gene_Input])) %>%
  dplyr::left_join(mapping_all, by = c("Gene_Input" = "ALIAS")) %>%
  dplyr::mutate(
    MappingStatus = dplyr::case_when(
      is.na(Gene) ~ "unmapped",
      !Gene %in% ref_genes ~ "mapped_not_in_CTD",
      duplicated(Gene) | duplicated(Gene, fromLast = TRUE) ~ "mapped_duplicate_symbol",
      TRUE ~ "mapped_unique"
    )
  )

mapping_qc <- clean_df %>%
  dplyr::select(RowID, Genes, Gene_Input, Gene, MappingStatus, N_SYMBOL_Matches, All_SYMBOL_Matches)

openxlsx::write.xlsx(mapping_qc, file.path(dirs$qc, "gene_mapping_audit.xlsx"), overwrite = TRUE)
write.table(
  clean_df %>% dplyr::filter(MappingStatus == "mapped_not_in_CTD") %>% dplyr::pull(Gene) %>% unique(),
  file.path(dirs$qc, "dropped_genes_not_in_CTD.txt"),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

mapped_clean_df <- clean_df %>%
  dplyr::filter(!is.na(Gene), Gene %in% ref_genes)

expr_mat <- make_expr_mat(mapped_clean_df, sample_cols)
expr_mat <- expr_mat[rownames(expr_mat) %in% ref_genes, , drop = FALSE]
background_universe <- intersect(rownames(expr_mat), ref_genes)

condition_lookup <- load_sample_condition_lookup(sample_metadata_path)
sample_meta <- parse_sample_metadata(colnames(expr_mat), condition_lookup)
if (all(is.na(sample_meta$Cond))) {
  stop(
    "No sample conditions could be resolved. Provide condition tokens in sample names ",
    "or set sample_metadata_path to a metadata file with AnimalID/sample_id and ExpGroup/Condition columns."
  )
}
if (any(is.na(sample_meta$Cond))) {
  warning(sum(is.na(sample_meta$Cond)), " sample(s) have no resolved condition and will be excluded from contrasts.")
}
sample_meta_qc <- sample_meta %>%
  dplyr::count(Stratum, Region, Layer, Cond, name = "N_Samples") %>%
  dplyr::arrange(Region, Layer, Cond)

openxlsx::write.xlsx(
  list(SampleMetadata = sample_meta, SampleCounts = sample_meta_qc, ConditionLookup = condition_lookup),
  file.path(dirs$qc, "sample_metadata_qc.xlsx"),
  overwrite = TRUE
)

# ==========================================
# 4. TARGET SIGNATURES
# ==========================================

message("Step 3: Creating baseline and modeled differential target signatures...")

long_df <- expr_mat %>%
  as.data.frame(check.names = FALSE) %>%
  tibble::rownames_to_column("Gene") %>%
  tidyr::pivot_longer(-Gene, names_to = "Sample", values_to = "Exp") %>%
  dplyr::left_join(sample_meta, by = "Sample") %>%
  dplyr::filter(!is.na(Stratum), !is.na(Cond))

baseline_mat <- long_df %>%
  dplyr::group_by(Gene, Stratum, Cond) %>%
  dplyr::summarise(MeanExp = mean(Exp, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(ID = paste(Stratum, Cond, sep = "__")) %>%
  dplyr::select(Gene, ID, MeanExp) %>%
  tidyr::pivot_wider(names_from = ID, values_from = MeanExp) %>%
  tibble::column_to_rownames("Gene") %>%
  as.matrix()

analysis_strata <- sample_meta %>%
  dplyr::filter(!is.na(Stratum)) %>%
  dplyr::distinct(Stratum, Region, Layer) %>%
  dplyr::arrange(Region, Layer) %>%
  dplyr::pull(Stratum)

de_tbl <- dplyr::bind_rows(lapply(analysis_strata, function(stratum) {
  run_limma_stratum(expr_mat, sample_meta, stratum)
}))

if (nrow(de_tbl) == 0) {
  stop("No differential contrasts were created. Check sample names for region, optional layer, and condition labels.")
}

target_gene_tbl <- dplyr::bind_rows(
  make_baseline_targets(baseline_mat, analysis_params$top_n_values),
  make_differential_targets(de_tbl, analysis_params$top_n_values)
)

target_manifest <- target_gene_tbl %>%
  dplyr::group_by(Target, AnalysisType, Stratum, Region, Layer, Contrast, Direction, Metric, TopN) %>%
  dplyr::summarise(
    N_Hits = dplyr::n_distinct(Gene),
    .groups = "drop"
  )

openxlsx::write.xlsx(
  list(
    DifferentialModelResults = de_tbl,
    TargetManifest = target_manifest,
    TargetGenes = target_gene_tbl
  ),
  file.path(dirs$tables, "EWCE_input_signatures.xlsx"),
  overwrite = TRUE
)

# ==========================================
# 5. EWCE BOOTSTRAPPING
# ==========================================

message("Step 4: Running EWCE bootstraps across sensitivity grid...")
workers <- max(1, parallel::detectCores() - 1)
future::plan(future::multisession, workers = workers)

target_grid <- target_manifest %>%
  tidyr::crossing(AnnotLevel = analysis_params$annot_levels) %>%
  dplyr::mutate(TargetRun = paste(Target, paste0("top", TopN), paste0("annot", AnnotLevel), sep = "__"))

run_ewce_target <- function(target_run) {
  cache_file <- file.path(dirs$cache, paste0(safe_file_stem(target_run), ".rds"))
  if (file.exists(cache_file)) {
    return(readRDS(cache_file))
  }

  row <- target_grid %>% dplyr::filter(TargetRun == target_run) %>% dplyr::slice(1)
  legacy_cache_file <- file.path(dirs$cache, paste0(safe_file_stem(paste(row$Target, paste0("annot", row$AnnotLevel), sep = "__")), ".rds"))
  if (file.exists(legacy_cache_file)) {
    legacy_out <- readRDS(legacy_cache_file)
    if (all(legacy_out$TopN == row$TopN, na.rm = TRUE)) {
      saveRDS(legacy_out, cache_file)
      return(legacy_out)
    }
  }

  genes <- target_gene_tbl %>%
    dplyr::filter(Target == row$Target, TopN == row$TopN) %>%
    dplyr::pull(Gene)
  bg <- intersect(background_universe, ref_genes_by_level[[row$AnnotLevel]])

  res <- run_ewce_once(genes, bg, row$AnnotLevel)

  out <- res %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      Target = row$Target,
      AnalysisType = row$AnalysisType,
      Stratum = row$Stratum,
      Region = row$Region,
      Layer = row$Layer,
      Contrast = row$Contrast,
      Direction = row$Direction,
      Metric = row$Metric,
      TopN = row$TopN,
      AnnotLevel = row$AnnotLevel,
      N_Hits = row$N_Hits,
      N_Background = length(bg)
    )

  saveRDS(out, cache_file)
  out
}

results_all <- future.apply::future_lapply(
  target_grid$TargetRun,
  run_ewce_target,
  future.seed = analysis_params$seed
) %>%
  dplyr::bind_rows()

future::plan(future::sequential)

results_all <- results_all %>%
  dplyr::group_by(Target, AnnotLevel) %>%
  dplyr::mutate(q_target = p.adjust(p, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    q_global = p.adjust(p, method = "BH"),
    neglog10q_target = -log10(pmax(q_target, 1e-300)),
    neglog10q_global = -log10(pmax(q_global, 1e-300)),
    SignedSig_Target = sign(sd_from_mean) * neglog10q_target,
    SignedSig_Global = sign(sd_from_mean) * neglog10q_global,
    Significant_Global = q_global < analysis_params$fdr_alpha,
    Direction_EWCE = dplyr::case_when(
      sd_from_mean > 0 ~ "Positive",
      sd_from_mean < 0 ~ "Negative",
      TRUE ~ "Neutral"
    )
  )

primary_results <- results_all %>%
  dplyr::filter(
    AnnotLevel == analysis_params$primary_annot_level,
    TopN == analysis_params$primary_top_n
  )

driver_overlap_tbl <- celltype_marker_overlap(
  results_all,
  target_gene_tbl,
  analysis_params$primary_annot_level
)

sensitivity_tbl <- results_all %>%
  dplyr::group_by(AnalysisType, Stratum, Region, Layer, Metric, Direction, CellType, AnnotLevel) %>%
  dplyr::summarise(
    N_TopN_Tested = dplyr::n_distinct(TopN),
    N_TopN_GlobalSig = dplyr::n_distinct(TopN[Significant_Global]),
    Min_q_global = min(q_global, na.rm = TRUE),
    Max_abs_Z = max(abs(sd_from_mean), na.rm = TRUE),
    Direction_Consistent = dplyr::n_distinct(sign(sd_from_mean[is.finite(sd_from_mean)])) <= 1,
    RobustAcrossTopN = N_TopN_GlobalSig == N_TopN_Tested & Direction_Consistent,
    .groups = "drop"
  )

# ==========================================
# 5b. PUBLICATION-ORIENTED SYNTHESIS (CACHE-PRESERVING)
# ==========================================

message("Step 4b: Building Publication-oriented synthesis tables...")

stratum_order <- c(
  "CA1_so", "CA1_sp", "CA1_sr", "CA1_slm",
  "CA2_so", "CA2_sp", "CA2_sr", "CA2_slm",
  "CA3_so", "CA3_sp", "CA3_sr", "CA3_slm",
  "DG_mo", "DG_po", "DG_sg"
)
observed_strata <- unique(stats::na.omit(c(results_all$Stratum, sample_meta$Stratum)))
stratum_order <- c(stratum_order[stratum_order %in% observed_strata], sort(setdiff(observed_strata, stratum_order)))

format_contrast_label <- function(x) {
  dplyr::recode(
    as.character(x),
    "Sus_vs_Con" = "Sus vs Con",
    "Res_vs_Con" = "Res vs Con",
    .default = as.character(x)
  )
}

format_metric_label <- function(contrast, direction) {
  paste(format_contrast_label(contrast), stringr::str_to_title(as.character(direction)))
}

# Symmetric cap for signed -log10(FDR): keeps diverging heatmap readable when q -> 0.
SIGNED_LOG10Q_CAP <- 10L

cap_signed_value <- function(x, limit = SIGNED_LOG10Q_CAP) {
  pmax(pmin(x, limit), -limit)
}

# Cell types ranked for heatmaps: significant first, then best global FDR, then |Z|.
rank_celltypes_by_evidence <- function(
    df,
    celltype_col = "CellType",
    q_col = "q_global",
    z_col = "sd_from_mean",
    sig_col = "Significant_Global"
) {
  df %>%
    dplyr::group_by(.data[[celltype_col]]) %>%
    dplyr::summarise(
      best_q_global = min(.data[[q_col]], na.rm = TRUE),
      max_abs_z = max(abs(.data[[z_col]]), na.rm = TRUE),
      any_sig = any(.data[[sig_col]] %in% TRUE, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(
      dplyr::desc(any_sig),
      best_q_global,
      dplyr::desc(max_abs_z)
    ) %>%
    dplyr::pull(.data[[celltype_col]])
}

# Small dots on heatmap tiles: Significant_Global == TRUE; light/dark fill picks dot color.
geom_heatmap_sig_dots <- function(
    data,
    x,
    y,
    fill_col = "sd_from_mean",
    sig_col = "Significant_Global",
    fill_limit = NULL
) {
  sig_data <- data[data[[sig_col]] %in% TRUE, , drop = FALSE]
  if (nrow(sig_data) == 0) {
    return(ggplot2::geom_blank())
  }
  mid_threshold <- if (is.null(fill_limit)) {
    stats::median(abs(sig_data[[fill_col]]), na.rm = TRUE)
  } else {
    max(abs(fill_limit)) * 0.45
  }
  sig_data$.dot_fill <- ifelse(abs(sig_data[[fill_col]]) >= mid_threshold, "white", "black")
  sig_data$.x <- sig_data[[x]]
  sig_data$.y <- sig_data[[y]]
  ggplot2::geom_point(
    data = sig_data,
    ggplot2::aes(x = .x, y = .y),
    inherit.aes = FALSE,
    shape = 16,
    size = 0.55,
    color = sig_data$.dot_fill,
    show.legend = FALSE
  )
}

scale_fill_signed_z <- function(
    name = "Z-score",
    low = "#0072B2",
    mid = "white",
    high = "#D55E00",
    midpoint = 0,
    limits = NULL
) {
  ggplot2::scale_fill_gradient2(
    low = low,
    mid = mid,
    high = high,
    midpoint = midpoint,
    limits = limits,
    name = name
  )
}

primary_results <- primary_results %>%
  dplyr::mutate(
    Stratum = factor(Stratum, levels = stratum_order),
    ContrastLabel = format_contrast_label(Contrast),
    MetricLabel = format_metric_label(Contrast, Direction),
    MetricLabel = dplyr::if_else(AnalysisType == "Baseline", as.character(Metric), MetricLabel),
    MetricLabel = factor(
      MetricLabel,
      levels = c("con", "res", "sus", "Sus vs Con Up", "Sus vs Con Down", "Res vs Con Up", "Res vs Con Down")
    ),
    SignedSig_Global_Capped = cap_signed_value(SignedSig_Global)
  )

results_all <- results_all %>%
  dplyr::mutate(
    Stratum = factor(Stratum, levels = stratum_order),
    ContrastLabel = format_contrast_label(Contrast),
    MetricLabel = format_metric_label(Contrast, Direction),
    MetricLabel = dplyr::if_else(AnalysisType == "Baseline", as.character(Metric), MetricLabel),
    SignedSig_Global_Capped = cap_signed_value(SignedSig_Global)
  )

sensitivity_tbl <- sensitivity_tbl %>%
  dplyr::mutate(
    Stratum = factor(Stratum, levels = stratum_order),
    ContrastLabel = format_contrast_label(stringr::str_remove(as.character(Metric), "_(up|down)$")),
    RobustnessScore = N_TopN_GlobalSig / pmax(N_TopN_Tested, 1)
  )

annotation_consistency_tbl <- sensitivity_tbl %>%
  dplyr::group_by(AnalysisType, Stratum, Region, Layer, Metric, Direction, CellType) %>%
  dplyr::summarise(
    N_AnnotLevels_Tested = dplyr::n_distinct(AnnotLevel),
    N_AnnotLevels_With_GlobalSig = dplyr::n_distinct(AnnotLevel[N_TopN_GlobalSig > 0]),
    Mean_RobustnessScore = mean(RobustnessScore, na.rm = TRUE),
    Max_RobustnessScore = max(RobustnessScore, na.rm = TRUE),
    Min_q_global = min(Min_q_global, na.rm = TRUE),
    Max_abs_Z = max(Max_abs_Z, na.rm = TRUE),
    RobustInPrimaryAnnot = any(AnnotLevel == analysis_params$primary_annot_level & RobustAcrossTopN),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Min_q_global, dplyr::desc(Max_abs_Z))

high_confidence_hits <- primary_results %>%
  dplyr::filter(AnalysisType == "Differential", Significant_Global) %>%
  dplyr::left_join(
    sensitivity_tbl %>%
      dplyr::filter(AnnotLevel == analysis_params$primary_annot_level) %>%
      dplyr::select(
        AnalysisType, Stratum, Region, Layer, Metric, Direction, CellType,
        N_TopN_Tested, N_TopN_GlobalSig, RobustnessScore, RobustAcrossTopN,
        Sensitivity_Min_q_global = Min_q_global,
        Sensitivity_Max_abs_Z = Max_abs_Z
      ),
    by = c("AnalysisType", "Stratum", "Region", "Layer", "Metric", "Direction", "CellType")
  ) %>%
  dplyr::left_join(
    annotation_consistency_tbl %>%
      dplyr::select(
        AnalysisType, Stratum, Region, Layer, Metric, Direction, CellType,
        N_AnnotLevels_Tested, N_AnnotLevels_With_GlobalSig, Mean_RobustnessScore, Max_RobustnessScore
      ),
    by = c("AnalysisType", "Stratum", "Region", "Layer", "Metric", "Direction", "CellType")
  ) %>%
  dplyr::arrange(q_global, dplyr::desc(abs(sd_from_mean))) %>%
  dplyr::select(
    Stratum, Region, Layer, ContrastLabel, Direction, CellType,
    sd_from_mean, fold_change, p, q_target, q_global,
    RobustnessScore, RobustAcrossTopN, N_TopN_GlobalSig, N_TopN_Tested,
    N_AnnotLevels_With_GlobalSig, N_AnnotLevels_Tested,
    Target, Metric, TopN, AnnotLevel, N_Hits, N_Background
  )

publication_top_celltypes <- primary_results %>%
  dplyr::filter(AnalysisType == "Differential") %>%
  dplyr::group_by(CellType) %>%
  dplyr::summarise(
    Best_q_global = min(q_global, na.rm = TRUE),
    Max_abs_SignedSig = max(abs(SignedSig_Global), na.rm = TRUE),
    Max_abs_Z = max(abs(sd_from_mean), na.rm = TRUE),
    Any_GlobalSignificant = any(Significant_Global, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(!Any_GlobalSignificant, Best_q_global, dplyr::desc(Max_abs_Z)) %>%
  dplyr::slice_head(n = 35)

publication_diff_heatmap_tbl <- primary_results %>%
  dplyr::filter(AnalysisType == "Differential", CellType %in% publication_top_celltypes$CellType) %>%
  dplyr::mutate(
    # Rows: best global FDR, then max |Z| (most informative cell types at top of heatmap).
    CellType = factor(CellType, levels = rev(publication_top_celltypes$CellType)),
    Stratum = factor(Stratum, levels = stratum_order),
    # Plotted fill matches SignedSig_Global_Capped (±SIGNED_LOG10Q_CAP) in the figure legend.
    SignedSig_Global_Plot = SignedSig_Global_Capped
  )

publication_robust_tbl <- sensitivity_tbl %>%
  dplyr::filter(AnalysisType == "Differential", AnnotLevel == analysis_params$primary_annot_level) %>%
  dplyr::arrange(Min_q_global, dplyr::desc(Max_abs_Z)) %>%
  dplyr::mutate(
    FindingLabel = paste(Stratum, Metric, CellType, sep = " | "),
    FindingLabel = factor(FindingLabel, levels = rev(unique(FindingLabel)))
  ) %>%
  dplyr::slice_head(n = 35)

# ==========================================
# 6. PUBLICATION-STYLE VISUALIZATION
# ==========================================

message("Step 5: Building publication and extended-data figures...")

col_baseline <- "#4E79A7"
col_sus      <- "#D55E00"
col_res      <- "#0072B2"
col_down     <- "#009E73"

# Color mapping only: symmetric cap so a few extreme bootstrap Z values do not wash out mid-tones.
HEATMAP_Z_COLOR_CAP <- max(
  3,
  stats::quantile(abs(primary_results$sd_from_mean), 0.99, na.rm = TRUE, names = FALSE)
)

baseline_top_celltypes <- primary_results %>%
  dplyr::filter(AnalysisType == "Baseline") %>%
  dplyr::group_by(CellType) %>%
  dplyr::summarise(
    best_q_global = min(q_global, na.rm = TRUE),
    max_abs_z = max(abs(sd_from_mean), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(best_q_global, dplyr::desc(max_abs_z)) %>%
  dplyr::slice_head(n = 25)

baseline_celltype_order <- rank_celltypes_by_evidence(
  primary_results %>% dplyr::filter(AnalysisType == "Baseline", CellType %in% baseline_top_celltypes$CellType)
)

baseline_heatmap_tbl <- primary_results %>%
  dplyr::filter(AnalysisType == "Baseline", CellType %in% baseline_top_celltypes$CellType) %>%
  dplyr::mutate(
    Metric = factor(Metric, levels = analysis_params$conditions),
    CellType = factor(CellType, levels = rev(baseline_celltype_order)),
    sd_from_mean_plot = cap_signed_value(sd_from_mean, HEATMAP_Z_COLOR_CAP)
  )

# Fig 1A: only positive values, use sequential color scale
p1 <- ggplot2::ggplot(
  baseline_heatmap_tbl,
  ggplot2::aes(x = Metric, y = CellType, fill = sd_from_mean_plot)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.1) +
  geom_heatmap_sig_dots(
    baseline_heatmap_tbl,
    x = "Metric",
    y = "CellType",
    fill_col = "sd_from_mean_plot",
    fill_limit = c(0, HEATMAP_Z_COLOR_CAP)
  ) +
  ggplot2::facet_wrap(~Stratum, nrow = 1) +
  ggplot2::scale_fill_gradient(
    name = "Z-score",
    low = "white",
    high = col_sus,
    limits = c(0, HEATMAP_Z_COLOR_CAP)
  ) +
  theme_publication() +
  ggplot2::labs(x = NULL, y = "Cell type", title = "Baseline cell-type enrichment") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

diff_celltype_order <- rank_celltypes_by_evidence(
  primary_results %>% dplyr::filter(AnalysisType == "Differential")
)

diff_metric_order <- c(
  "Sus_vs_Con_up", "Sus_vs_Con_down", "Res_vs_Con_up", "Res_vs_Con_down"
)
diff_metric_order <- diff_metric_order[diff_metric_order %in% unique(primary_results$Metric[primary_results$AnalysisType == "Differential"])]
diff_metric_order <- c(diff_metric_order, setdiff(unique(primary_results$Metric[primary_results$AnalysisType == "Differential"]), diff_metric_order))

diff_heatmap_tbl <- primary_results %>%
  dplyr::filter(AnalysisType == "Differential") %>%
  dplyr::mutate(
    Metric = factor(Metric, levels = diff_metric_order),
    CellType = factor(CellType, levels = rev(diff_celltype_order)),
    sd_from_mean_plot = cap_signed_value(sd_from_mean, HEATMAP_Z_COLOR_CAP),
    SignedSig_Global_Plot = SignedSig_Global_Capped
  )

p2 <- ggplot2::ggplot(
  diff_heatmap_tbl,
  ggplot2::aes(x = Metric, y = CellType, fill = sd_from_mean_plot)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.1) +
  geom_heatmap_sig_dots(
    diff_heatmap_tbl,
    x = "Metric",
    y = "CellType",
    fill_col = "sd_from_mean_plot",
    fill_limit = c(-HEATMAP_Z_COLOR_CAP, HEATMAP_Z_COLOR_CAP)
  ) +
  ggplot2::facet_wrap(~Stratum, nrow = 1) +
  scale_fill_signed_z(
    limits = c(-HEATMAP_Z_COLOR_CAP, HEATMAP_Z_COLOR_CAP)
  ) +
  theme_publication() +
  ggplot2::labs(x = NULL, y = "Cell type", title = "Modeled stress contrasts") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

p3 <- ggplot2::ggplot(primary_results, ggplot2::aes(x = sd_from_mean, y = -log10(q_global))) +
  ggplot2::geom_point(ggplot2::aes(color = AnalysisType), alpha = 0.65, size = 1) +
  ggplot2::scale_color_manual(values = c("Baseline" = col_baseline, "Differential" = col_sus)) +
  ggplot2::geom_hline(yintercept = -log10(analysis_params$fdr_alpha), linetype = "dashed", linewidth = 0.2) +
  theme_publication() +
  ggplot2::labs(x = "Effect size (Z-score)", y = "-log10(global FDR)", title = "EWCE significance")

p4 <- ggplot2::ggplot(
  diff_heatmap_tbl,
  ggplot2::aes(x = Metric, y = CellType, fill = SignedSig_Global_Plot)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.1) +
  geom_heatmap_sig_dots(
    diff_heatmap_tbl,
    x = "Metric",
    y = "CellType",
    fill_col = "SignedSig_Global_Plot",
    fill_limit = c(-SIGNED_LOG10Q_CAP, SIGNED_LOG10Q_CAP)
  ) +
  ggplot2::facet_wrap(~Stratum, nrow = 1) +
  ggplot2::scale_fill_gradient2(
    low = col_res,
    mid = "white",
    high = col_sus,
    midpoint = 0,
    limits = c(-SIGNED_LOG10Q_CAP, SIGNED_LOG10Q_CAP),
    name = "Signed -log10(FDR)"
  ) +
  theme_publication() +
  ggplot2::labs(x = NULL, y = "Cell type", title = "Direction and global significance") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

top_celltypes <- primary_results %>%
  dplyr::group_by(CellType) %>%
  dplyr::summarise(
    max_abs_effect = max(abs(sd_from_mean), na.rm = TRUE),
    best_q_global = min(q_global, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(max_abs_effect)) %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::mutate(
    neglog10_q_global = -log10(pmax(best_q_global, 1e-300)),
    CellType = factor(CellType, levels = rev(CellType))
  )

p5 <- ggplot2::ggplot(
  top_celltypes,
  ggplot2::aes(x = max_abs_effect, y = CellType, color = neglog10_q_global)
) +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = max_abs_effect, yend = CellType),
    linewidth = 0.25,
    color = "grey55"
  ) +
  ggplot2::geom_point(size = 1.6, stroke = 0) +
  viridis::scale_color_viridis(option = "plasma", name = "-log10(global FDR)") +
  theme_publication() +
  ggplot2::labs(x = "Max |Z-score|", y = NULL, title = "Top cell-type effects")

p6_strip_tbl <- primary_results %>%
  dplyr::filter(AnalysisType == "Differential") %>%
  dplyr::mutate(
    Stratum = factor(Stratum, levels = stratum_order),
    Direction = factor(Direction, levels = c("up", "down"))
  )

# Jittered points (no KDE/ridges): avoids implying smooth densities in sparse strata.
p6 <- ggplot2::ggplot(
  p6_strip_tbl,
  ggplot2::aes(x = sd_from_mean, y = Stratum, color = Direction)
) +
  ggplot2::geom_jitter(
    position = ggplot2::position_jitter(width = 0, height = 0.22),
    size = 0.45,
    alpha = 0.55,
    stroke = 0
  ) +
  ggplot2::scale_color_manual(
    values = c("up" = col_sus, "down" = col_down),
    labels = c("up" = "Up", "down" = "Down"),
    name = NULL
  ) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.2, color = "grey35") +
  theme_publication() +
  ggplot2::labs(x = "Z-score", y = NULL, title = "Stratum-level EWCE effects (raw points)")

sensitivity_plot_tbl <- sensitivity_tbl %>%
  dplyr::filter(AnnotLevel == analysis_params$primary_annot_level) %>%
  dplyr::arrange(Min_q_global, dplyr::desc(Max_abs_Z)) %>%
  dplyr::slice_head(n = 30) %>%
  dplyr::mutate(
    FindingLabel = factor(
      paste(Stratum, Metric, CellType, sep = " | "),
      levels = rev(unique(paste(Stratum, Metric, CellType, sep = " | ")))
    )
  )

p7 <- ggplot2::ggplot(
  sensitivity_plot_tbl,
  ggplot2::aes(x = N_TopN_GlobalSig, y = FindingLabel, color = RobustAcrossTopN)
) +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = N_TopN_GlobalSig, yend = FindingLabel),
    linewidth = 0.25,
    color = "grey70"
  ) +
  ggplot2::geom_point(size = 1.5, stroke = 0) +
  ggplot2::scale_color_manual(values = c("FALSE" = "grey55", "TRUE" = col_sus), name = "Robust") +
  ggplot2::scale_x_continuous(breaks = sort(unique(sensitivity_plot_tbl$N_TopN_Tested))) +
  theme_publication() +
  ggplot2::labs(x = "Significant top-N settings", y = NULL, title = "Hit-list sensitivity")

publication_heatmap_fig <- ggplot2::ggplot(
  publication_diff_heatmap_tbl,
  ggplot2::aes(x = Stratum, y = CellType, fill = SignedSig_Global_Plot)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.12) +
  geom_heatmap_sig_dots(
    publication_diff_heatmap_tbl,
    x = "Stratum",
    y = "CellType",
    fill_col = "SignedSig_Global_Plot",
    fill_limit = c(-SIGNED_LOG10Q_CAP, SIGNED_LOG10Q_CAP)
  ) +
  ggplot2::facet_grid(. ~ MetricLabel, scales = "free_x", space = "free_x") +
  ggplot2::scale_fill_gradient2(
    low = col_res,
    mid = "white",
    high = col_sus,
    midpoint = 0,
    limits = c(-SIGNED_LOG10Q_CAP, SIGNED_LOG10Q_CAP),
    name = "Signed\n-log10(FDR)"
  ) +
  theme_publication() +
  ggplot2::labs(
    x = NULL,
    y = "Cell type",
    title = "Spatially resolved differential EWCE signatures"
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.spacing.x = grid::unit(1.5, "mm"),
    legend.position = "right"
  )

publication_robustness_fig <- ggplot2::ggplot(
  publication_robust_tbl,
  ggplot2::aes(x = RobustnessScore, y = FindingLabel)
) +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = RobustnessScore, yend = FindingLabel),
    linewidth = 0.3,
    color = "grey80"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      size = Max_abs_Z,
      fill = RobustAcrossTopN,
      color = -log10(pmax(Min_q_global, 1e-300))
    ),
    shape = 21,
    stroke = 0.2
  ) +
  ggplot2::scale_fill_manual(values = c("FALSE" = "grey85", "TRUE" = col_sus), name = "Robust\nall Top-N") +
  ggplot2::scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  viridis::scale_color_viridis(option = "plasma", name = "-log10\nmin FDR") +
  ggplot2::scale_size_area(max_size = 3.2, name = "Max |Z|") +
  theme_publication() +
  ggplot2::labs(
    x = "Fraction of Top-N settings significant",
    y = NULL,
    title = "High-confidence EWCE findings across hit-list sizes"
  )

main_fig <- (p1 / p2 | p3) +
  patchwork::plot_layout(widths = c(2, 1)) +
  patchwork::plot_annotation(tag_levels = "A") &
  ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 10))

supp_fig <- (p4 / p5 / p6 / p7) +
  patchwork::plot_layout(heights = c(1.1, 1, 0.8, 1)) +
  patchwork::plot_annotation(tag_levels = "A") &
  ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 10))

diff_heatmap_tbl <- primary_results %>%
  dplyr::filter(AnalysisType == "Differential") %>%
  dplyr::mutate(TargetLabel = paste(Stratum, Metric, sep = "_")) %>%
  dplyr::select(CellType, TargetLabel, SignedSig_Global) %>%
  dplyr::distinct() %>%
  tidyr::pivot_wider(names_from = TargetLabel, values_from = SignedSig_Global)

diff_heatmap_mat <- diff_heatmap_tbl %>%
  tibble::column_to_rownames("CellType") %>%
  as.matrix()

# ==========================================
# 7. EXPORTING
# ==========================================

message("Step 6: Exporting figures, source data, and reproducibility metadata...")

ggplot2::ggsave(file.path(dirs$plots, "Fig1_EWCE_Summary.pdf"), main_fig, width = 180, height = 200, units = "mm", device = grDevices::cairo_pdf)
ggplot2::ggsave(file.path(dirs$svgs, "Fig1_EWCE_Summary.svg"), main_fig, width = 180, height = 200, units = "mm")
ggplot2::ggsave(file.path(dirs$svgs, "Volcano_Panel.svg"), p3, width = 80, height = 80, units = "mm")
ggplot2::ggsave(file.path(dirs$plots, "FigS1_EWCE_Additional_Visuals.pdf"), supp_fig, width = 180, height = 280, units = "mm", device = grDevices::cairo_pdf)
ggplot2::ggsave(file.path(dirs$svgs, "FigS1_EWCE_Additional_Visuals.svg"), supp_fig, width = 180, height = 280, units = "mm")
ggplot2::ggsave(file.path(dirs$plots, "Fig2_Primary_Differential_EWCE_Heatmap.pdf"), publication_heatmap_fig, width = 183, height = 170, units = "mm", device = grDevices::cairo_pdf)
ggplot2::ggsave(file.path(dirs$svgs, "Fig2_Primary_Differential_EWCE_Heatmap.svg"), publication_heatmap_fig, width = 183, height = 170, units = "mm")
ggplot2::ggsave(file.path(dirs$plots, "Fig3_Robust_EWCE_Findings.pdf"), publication_robustness_fig, width = 183, height = 160, units = "mm", device = grDevices::cairo_pdf)
ggplot2::ggsave(file.path(dirs$svgs, "Fig3_Robust_EWCE_Findings.svg"), publication_robustness_fig, width = 183, height = 160, units = "mm")

if (nrow(diff_heatmap_mat) > 2 && ncol(diff_heatmap_mat) > 2) {
  grDevices::pdf(file.path(dirs$plots, "FigS2_EWCE_ClusteredHeatmap.pdf"), width = 8, height = 10, family = "sans")
  pheatmap::pheatmap(
    diff_heatmap_mat,
    color = grDevices::colorRampPalette(c(col_res, "white", col_sus))(101),
    border_color = NA,
    fontsize = 6,
    fontsize_row = 5,
    fontsize_col = 6,
    angle_col = 45,
    main = "Signed global significance by target"
  )
  grDevices::dev.off()
}

summary_tbl <- results_all %>%
  dplyr::group_by(AnalysisType, AnnotLevel, TopN, Stratum, Region, Layer, Metric) %>%
  dplyr::summarise(
    N_Tested = dplyr::n(),
    N_GlobalSignificant = sum(Significant_Global, na.rm = TRUE),
    Median_Abs_Effect = median(abs(sd_from_mean), na.rm = TRUE),
    Mean_Abs_Effect = mean(abs(sd_from_mean), na.rm = TRUE),
    Top_CellType = CellType[which.max(abs(sd_from_mean))],
    Top_Effect = max(abs(sd_from_mean), na.rm = TRUE),
    Best_q_target = min(q_target, na.rm = TRUE),
    Best_q_global = min(q_global, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(AnalysisType, AnnotLevel, TopN, Region, Layer, dplyr::desc(N_GlobalSignificant), dplyr::desc(Mean_Abs_Effect))

top_hits_tbl <- results_all %>%
  dplyr::group_by(Target, AnnotLevel) %>%
  dplyr::arrange(q_global, dplyr::desc(abs(sd_from_mean)), .by_group = TRUE) %>%
  dplyr::slice_head(n = 15) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    Target, AnalysisType, AnnotLevel, TopN, Stratum, Region, Layer, Metric, Direction,
    CellType, sd_from_mean, p, q_target, q_global, SignedSig_Global, N_Hits, N_Background
  )

sig_rank_tbl <- results_all %>%
  dplyr::mutate(Significant = q_global < analysis_params$fdr_alpha) %>%
  dplyr::group_by(CellType, AnalysisType, AnnotLevel) %>%
  dplyr::summarise(
    GlobalSignificant_Count = sum(Significant, na.rm = TRUE),
    Mean_SignedSig_Global = mean(SignedSig_Global, na.rm = TRUE),
    Max_Abs_Effect = max(abs(sd_from_mean), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(GlobalSignificant_Count), dplyr::desc(Max_Abs_Effect))

input_gene_stats <- tibble::tibble(
  Total_Input_Rows = nrow(raw_df),
  Rows_With_Gene_Input = nrow(mapping_input),
  Mapped_Rows = sum(!is.na(clean_df$Gene)),
  Unique_Mapped_Genes = dplyr::n_distinct(stats::na.omit(clean_df$Gene)),
  Unique_Mapped_In_CTD = length(background_universe),
  Dropped_Not_In_CTD = dplyr::n_distinct(clean_df$Gene[clean_df$MappingStatus == "mapped_not_in_CTD"]),
  Measured_Proteome_Background = length(background_universe),
  Input_SHA256 = input_hash
)

wb <- openxlsx::createWorkbook()
add_worksheet_safe(wb, "Full_Results", results_all)
add_worksheet_safe(wb, "Primary_Results", primary_results)
add_worksheet_safe(wb, "Significant_Global", results_all %>% dplyr::filter(q_global < analysis_params$fdr_alpha))
add_worksheet_safe(wb, "Summary_by_Contrast", summary_tbl)
add_worksheet_safe(wb, "Top_Hits_per_Target", top_hits_tbl)
add_worksheet_safe(wb, "Sensitivity", sensitivity_tbl)
add_worksheet_safe(wb, "High_Confidence", high_confidence_hits)
add_worksheet_safe(wb, "Annotation_Consistency", annotation_consistency_tbl)
add_worksheet_safe(wb, "CellType_Significance_Rank", sig_rank_tbl)
add_worksheet_safe(wb, "Driver_Marker_Overlap", driver_overlap_tbl)
add_worksheet_safe(wb, "Input_Gene_Stats", input_gene_stats)
add_worksheet_safe(wb, "Sample_Counts", sample_meta_qc)
openxlsx::saveWorkbook(wb, file.path(dirs$tables, "Supplementary_Table_EWCE.xlsx"), overwrite = TRUE)

publication_wb <- openxlsx::createWorkbook()
add_worksheet_safe(publication_wb, "High_Confidence_Findings", high_confidence_hits)
add_worksheet_safe(publication_wb, "Primary_Diff_Heatmap_Data", publication_diff_heatmap_tbl)
add_worksheet_safe(publication_wb, "Robustness_Data", publication_robust_tbl)
add_worksheet_safe(publication_wb, "Annotation_Consistency", annotation_consistency_tbl)
add_worksheet_safe(publication_wb, "Top_CellTypes", publication_top_celltypes)
openxlsx::saveWorkbook(publication_wb, file.path(dirs$tables, "High_Confidence_EWCE_Findings.xlsx"), overwrite = TRUE)

source_wb <- openxlsx::createWorkbook()
add_worksheet_safe(source_wb, "Fig1A_Baseline_Heatmap", baseline_heatmap_tbl)
add_worksheet_safe(source_wb, "Fig1B_Diff_Heatmap", diff_heatmap_tbl)
add_worksheet_safe(source_wb, "Fig1C_Volcano", primary_results)
add_worksheet_safe(
  source_wb,
  "FigS1A_Signed_Heatmap",
  diff_heatmap_tbl %>% dplyr::select(dplyr::any_of(c(
    "Stratum", "Region", "Layer", "Metric", "CellType", "Direction",
    "sd_from_mean", "SignedSig_Global", "SignedSig_Global_Plot", "Significant_Global", "q_global"
  )))
)
add_worksheet_safe(source_wb, "FigS1B_Top_CellTypes", top_celltypes)
add_worksheet_safe(source_wb, "FigS1C_Stratum_Strip", p6_strip_tbl)
add_worksheet_safe(source_wb, "FigS1D_Sensitivity", sensitivity_plot_tbl)
add_worksheet_safe(
  source_wb,
  "Fig2_Primary_Heatmap",
  publication_diff_heatmap_tbl %>% dplyr::select(dplyr::any_of(c(
    "Stratum", "MetricLabel", "CellType", "Direction", "sd_from_mean",
    "SignedSig_Global", "SignedSig_Global_Plot", "Significant_Global", "q_global"
  )))
)
add_worksheet_safe(source_wb, "Fig3_Robust_Findings", publication_robust_tbl)
add_worksheet_safe(
  source_wb,
  "Plot_Caps",
  tibble::tibble(
    Parameter = c("HEATMAP_Z_COLOR_CAP", "SIGNED_LOG10Q_CAP"),
    Value = c(HEATMAP_Z_COLOR_CAP, SIGNED_LOG10Q_CAP),
    Note = c(
      "Symmetric cap applied to sd_from_mean for heatmap fill only (Fig1A, Fig1B, FigS1 panels).",
      paste0("Symmetric cap on sign(q)*-log10(q) for signed FDR heatmaps (Fig2, FigS1A, Fig4); raw values in SignedSig_Global.")
    )
  )
)
openxlsx::saveWorkbook(source_wb, file.path(dirs$source, "Source_Data_EWCE_Figures.xlsx"), overwrite = TRUE)

saveRDS(
  list(
    results_all = results_all,
    primary_results = primary_results,
    differential_model_results = de_tbl,
    target_gene_tbl = target_gene_tbl,
    target_manifest = target_manifest,
    background_universe = background_universe,
    sample_metadata = sample_meta,
    condition_lookup = condition_lookup,
    mapping_qc = mapping_qc,
    high_confidence_hits = high_confidence_hits,
    annotation_consistency_tbl = annotation_consistency_tbl,
    publication_diff_heatmap_tbl = publication_diff_heatmap_tbl,
    publication_robust_tbl = publication_robust_tbl,
    analysis_params = analysis_params
  ),
  file.path(dirs$data, "EWCE_results_full.rds")
)

reproducibility_lines <- c(
  paste0("Run date: ", Sys.time()),
  paste0("Input file: ", data_path),
  paste0("Sample metadata file: ", sample_metadata_path),
  paste0("Input SHA256: ", input_hash),
  paste0("EWCE reps: ", analysis_params$reps),
  paste0("Top-N values: ", paste(analysis_params$top_n_values, collapse = ", ")),
  paste0("Annotation levels: ", paste(analysis_params$annot_levels, collapse = ", ")),
  paste0("Primary top-N: ", analysis_params$primary_top_n),
  paste0("Primary annotation level: ", analysis_params$primary_annot_level),
  paste0("Measured proteome background size: ", length(background_universe)),
  paste0("Heatmap Z color cap (99th percentile, min 3): ", HEATMAP_Z_COLOR_CAP),
  paste0("Signed -log10(FDR) color cap: +/-", SIGNED_LOG10Q_CAP),
  "",
  "Session info:",
  capture.output(utils::sessionInfo())
)

writeLines(reproducibility_lines, file.path(dirs$qc, "reproducibility_session_info.txt"))

cat("\nPipeline complete. Files organized in:", base_results, "\n")

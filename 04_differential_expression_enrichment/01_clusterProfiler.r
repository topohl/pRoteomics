#' ============================================================
#' Consumes:
#'   - mapped contrast CSVs from data/processed/02_id_mapping/... (configurable)
#'   - optional background/config inputs from config/clusterProfiler_config.yml
#' Produces:
#'   - canonical enrichment tables, figures, source data, reports and logs under
#'     data/processed/04_differential_expression_enrichment/clusterProfiler/<dataset> and
#'     results/*/04_differential_expression_enrichment/clusterProfiler/<dataset>
#'   - clusterProfiler_manifest.csv consumed by 02_compareGO.r
#' Contract:
#'   - docs/file_contracts.tsv object_id clusterProfiler_manifest
#' ============================================================
#' HIGH-THROUGHPUT Gene Set Enrichment Analysis (GSEA) Workflow
#' WITH PARALLEL PROCESSING AND GO TERM SIMPLIFICATION
#' ============================================================
#' 
#' This script performs comprehensive GSEA, KEGG, ORA, custom pathway analysis
#' for MULTIPLE cell type comparisons in parallel using future.
#' Includes GO term redundancy removal via simplify().
#' Celltype scoring runs separately at the end.
#'
#' @author Tobias Pohl
#' ============================================================

setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
MODULE_ID <- "04_differential_expression_enrichment"
SUBSTEP_ID <- "clusterProfiler"
CANONICAL_PATHS <- create_module_dirs(MODULE_ID, SUBSTEP_ID)

# ----------------------------------------------------
# 0. SIMPLIFICATION SETTINGS
# ----------------------------------------------------
# *** MAIN CONTROL: Set whether to perform simplification at all ***
PERFORM_SIMPLIFICATION <- FALSE  # Set to FALSE to skip simplification entirely

# Control which version to use for plots (only matters if PERFORM_SIMPLIFICATION = TRUE)
USE_SIMPLIFIED_FOR_PLOTS <- FALSE  # TRUE = simplified plots, FALSE = full results

# Semantic similarity cutoff for simplify() function (only used if simplification is enabled)
# Range: 0.4 (strict, removes more) to 0.9 (lenient, keeps more)
# Default: 0.7 (moderate redundancy removal)
SIMPLIFY_CUTOFF <- 0.7

# Package installation policy:
# FALSE (recommended): fail fast with a clear list of missing packages.
# TRUE: auto-install missing packages (can take a long time and look "stuck").
AUTO_INSTALL_MISSING_PACKAGES <- FALSE
DRY_RUN <- is_dry_run()

# Note: 
# - If PERFORM_SIMPLIFICATION = FALSE, only full results are computed and saved
# - If PERFORM_SIMPLIFICATION = TRUE, both versions are computed and saved

# ----------------------------------------------------
# 1. PACKAGE SETUP
# ----------------------------------------------------
checkBiocManager <- function() {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    stop("Missing required R package: BiocManager. Install it explicitly before running this script.", call. = FALSE)
  }
}

installBioC <- function(bioc_packages, auto_install = AUTO_INSTALL_MISSING_PACKAGES) {
  # Install only missing packages to avoid unloading issues, set update=FALSE
  missing <- bioc_packages[!sapply(bioc_packages, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing Bioconductor packages detected: ",
      paste(missing, collapse = ", "),
      "\nInstall them explicitly before running this script.",
      call. = FALSE
    )
  }
}

loadPkg <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) warning(paste("Could not load:", pkg))
}

setupPackages <- function() {
  setup_start <- Sys.time()
  if (isTRUE(DRY_RUN)) {
    message("[Package setup] Dry run: skipping analysis package checks.")
    return(invisible(TRUE))
  }
  message("[Package setup] Checking dependencies...")

  # Ensure a CRAN mirror is set for non-interactive sessions
  if (is.null(getOption("repos"))) {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
  }
  
  checkBiocManager()
  
  required_packages <- c("clusterProfiler", "pathview", "enrichplot", "DOSE", "ggplot2", "ggnewscale",
                         "cowplot", "ggridges", "europepmc", "ggpubr", "ggrepel", "ggsci", "ggthemes",
                         "ggExtra", "ggforce", "ggalluvial", "lattice", "latticeExtra", "BiocManager",
                         "org.Mm.eg.db", "ggplotify", "svglite", "tidyr", "dplyr", "pheatmap", "proxy",
                         "tibble", "openxlsx", "future", "future.apply", "GOSemSim", "yaml", "progressr", "withr")
  bioc_packages <- c("clusterProfiler", "pathview", "enrichplot", "DOSE", "org.Mm.eg.db", "GOSemSim")
  
  # 1. Install missing BioC packages first
  installBioC(bioc_packages)
  
  # 2. Identify and install missing CRAN packages
  cran_packages <- setdiff(required_packages, bioc_packages)
  # Only check strictly for missing namespaces to avoid touching loaded packages
  missing_cran <- cran_packages[!sapply(cran_packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_cran) > 0) {
    stop(
      "Missing CRAN packages detected: ",
      paste(missing_cran, collapse = ", "),
      "\nInstall them explicitly before running this script.",
      call. = FALSE
    )
  }
  
  # 3. Quietly load only packages needed in the master process.
  # Worker-specific packages are loaded inside analyze_comparison().
  master_packages <- c("future", "future.apply", "yaml", "progressr")
  suppressPackageStartupMessages({
    invisible(lapply(master_packages, function(pkg) {
      if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        warning(paste("Could not load:", pkg))
      }
    }))
  })

  message("[Package setup] Ready in ", round(as.numeric(difftime(Sys.time(), setup_start, units = "secs")), 2), " sec")
}

# ----------------------------------------------------
# 2. DIRECTORY ORGANIZATION FUNCTIONS
# ----------------------------------------------------
create_analysis_dirs <- function(base_dir, comparison_name, ontology) {
  route <- classify_comparison_route(comparison_name)
  results_root <- file.path(CANONICAL_PATHS$processed, route$category, route$unit_folder, comparison_name)
  plots_root <- file.path(CANONICAL_PATHS$figures, route$category, route$unit_folder, comparison_name)

  dirs <- list(
    results = results_root,
    go_ont = file.path(results_root, "GO", ontology),
    kegg = file.path(results_root, "KEGG"),
    ora = file.path(results_root, "ORA"),
    custom = file.path(results_root, "Custom", ontology),
    pathview = file.path(results_root, "pathview"),
    plots_go = file.path(plots_root, paste0("GO_", ontology)),
    plots_kegg = file.path(plots_root, "KEGG"),
    plots_ora = file.path(plots_root, "ORA"),
    plots_custom = file.path(plots_root, "Custom", ontology),
    core_enrich = file.path(CANONICAL_PATHS$source_data, ontology),
    core_enrich_routed = file.path(CANONICAL_PATHS$source_data, ontology, route$category, route$unit_folder),
    route_category = route$category,
    route_unit = route$unit_folder
  )
  # Create only actual directory paths, not metadata fields.
  dir_fields <- c("results", "go_ont", "kegg", "ora", "custom", "pathview", "plots_go", "plots_kegg", "plots_ora", "plots_custom", "core_enrich", "core_enrich_routed")
  lapply(dirs[dir_fields], function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE))
  return(dirs)
}

save_plot_organized <- function(plot, filename, directory) {
  # Wrap in tryCatch as specific plots sometimes fail to render
  tryCatch({
    ggsave(file.path(directory, filename), plot, units = "cm", dpi = 300)
  }, error = function(e) {
    warning("Plot save failed for ", filename, ": ", e$message)
  })
}

read_gct <- function(file_path) {
  gct_data <- read.delim(file_path, skip = 2, header = FALSE, check.names = FALSE)
  if ("Name" %in% colnames(gct_data)) {
    colnames(gct_data)[1] <- "Gene"
  }
  return(gct_data)
}

read_config <- function(config_path) {
  defaults <- list(
    simplification = list(
      perform = PERFORM_SIMPLIFICATION,
      use_for_plots = USE_SIMPLIFIED_FOR_PLOTS,
      cutoff = SIMPLIFY_CUTOFF
    ),
    analysis = list(
      dataset = Sys.getenv("PROTEOMICS_COMPARISON", unset = "neuron_neuropil"),
      organism = "org.Mm.eg.db",
      ontology = "BP",
      top_gene_abs_log2fc = 1,
      pvalue_cutoff = 1,
      qvalue_cutoff = 1,
      p_adjust_method = "BH",
      min_gs_size = 10,
      max_gs_size = 800
    ),
    runtime = list(
      workers = -1,
      future_globals_max_size_mb = 8000,
      resume_if_complete = TRUE,
      force_rerun = FALSE,
      prewarm_workers = TRUE,
      show_progress = TRUE,
      show_step_progress = TRUE,
      dry_run = FALSE,
      config_file = config_path
    ),
    paths = list(
      mapped_dir = "",
      working_base = repo_root(),
      mapped_data_base = "",
      background_universe_file = ""
    ),
    optional_inputs = list(
      nk3r_genes = character(0),
      selected_uniprot = character(0),
      path_ids = character(0)
    ),
    background_universe = list(
      enabled = FALSE,
      column = "UNIPROT"
    )
  )

  cfg <- defaults
  if (file.exists(config_path)) {
    yaml_cfg <- tryCatch(yaml::read_yaml(config_path), error = function(e) NULL)
    if (!is.null(yaml_cfg) && is.list(yaml_cfg)) {
      cfg <- utils::modifyList(defaults, yaml_cfg)
      message("Loaded config: ", config_path)
    } else {
      warning("Config file exists but could not be parsed. Using defaults: ", config_path)
    }
  } else {
    message("Config file not found; using script defaults: ", config_path)
  }
  cfg
}

manifest_columns <- c(
  "analysis_id", "dataset", "run_id", "ontology", "result_type", "contrast", "comparison",
  "route_category", "route_unit", "condition", "direction", "simplified",
  "plot_suffix", "used_for_plot", "input_gene_file", "input_hash",
  "config_file", "config_hash", "output_table", "output_plot",
  "n_genes", "n_terms", "empty_result", "checkpoint_status", "created_at"
)

make_manifest_row <- function(result_type, ontology, comparison_name, dirs, input_gene_file,
                              config_path, output_table, output_plot = NA_character_,
                              n_genes = NA_integer_, n_terms = NA_integer_,
                              simplified = FALSE, used_for_plot = FALSE,
                              plot_suffix = NA_character_, checkpoint_status = "computed",
                              condition = NA_character_, direction = NA_character_) {
  dataset <- get0("DATASET", ifnotfound = NA_character_)
  data.frame(
    analysis_id = paste("clusterProfiler", dataset, result_type, ontology, comparison_name, sep = "::"),
    dataset = dataset,
    run_id = run_id,
    ontology = ontology,
    result_type = result_type,
    contrast = comparison_name,
    comparison = comparison_name,
    route_category = dirs$route_category,
    route_unit = dirs$route_unit,
    condition = condition,
    direction = direction,
    simplified = isTRUE(simplified),
    plot_suffix = plot_suffix,
    used_for_plot = isTRUE(used_for_plot),
    input_gene_file = input_gene_file,
    input_hash = file_hash(input_gene_file),
    config_file = config_path,
    config_hash = file_hash(config_path),
    output_table = output_table,
    output_plot = output_plot,
    n_genes = n_genes,
    n_terms = n_terms,
    empty_result = isTRUE(is.na(n_terms) || n_terms == 0),
    checkpoint_status = checkpoint_status,
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    stringsAsFactors = FALSE
  )[manifest_columns]
}

count_table_rows <- function(path) {
  if (!file.exists(path)) return(NA_integer_)
  out <- tryCatch(nrow(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)), error = function(e) NA_integer_)
  as.integer(out)
}

reconstruct_checkpoint_manifest <- function(dirs, comparison_name, data_path, config_path, ont, plot_suffix) {
  rows <- list()
  gsea_table <- file.path(dirs$core_enrich_routed, paste0(comparison_name, plot_suffix, ".csv"))
  if (file.exists(gsea_table)) {
    gsea_plot <- file.path(dirs$plots_go, paste0("GSEA_", ont, "_dotplot", plot_suffix, ".svg"))
    rows[[length(rows) + 1]] <- make_manifest_row(
      result_type = "GSEA_GO",
      ontology = ont,
      comparison_name = comparison_name,
      dirs = dirs,
      input_gene_file = data_path,
      config_path = config_path,
      output_table = gsea_table,
      output_plot = if (file.exists(gsea_plot)) gsea_plot else NA_character_,
      n_genes = NA_integer_,
      n_terms = count_table_rows(gsea_table),
      simplified = identical(plot_suffix, "_simplified"),
      used_for_plot = TRUE,
      plot_suffix = plot_suffix,
      checkpoint_status = "reconstructed_from_checkpoint"
    )
  }

  kegg_table <- file.path(CANONICAL_PATHS$source_data, "KEGG", dirs$route_category, dirs$route_unit, paste0(comparison_name, "_KEGG.csv"))
  if (file.exists(kegg_table)) {
    kegg_plot <- file.path(dirs$plots_kegg, "KEGG_dotplot.svg")
    rows[[length(rows) + 1]] <- make_manifest_row(
      result_type = "GSEA_KEGG",
      ontology = "KEGG",
      comparison_name = comparison_name,
      dirs = dirs,
      input_gene_file = data_path,
      config_path = config_path,
      output_table = kegg_table,
      output_plot = if (file.exists(kegg_plot)) kegg_plot else NA_character_,
      n_genes = NA_integer_,
      n_terms = count_table_rows(kegg_table),
      simplified = FALSE,
      used_for_plot = TRUE,
      plot_suffix = "_KEGG",
      checkpoint_status = "reconstructed_from_checkpoint"
    )
  }

  dplyr::bind_rows(rows)
}

make_checkpoint_skip_result <- function(cell_types, working_base, mapped_data_base, ont, config_path) {
  comparison_name <- paste(cell_types, collapse = "_")
  dirs <- create_analysis_dirs(working_base, comparison_name, ont)
  data_path <- file.path(mapped_data_base, paste0(comparison_name, ".csv"))
  qc_path <- file.path(dirs$results, "QC_summary.csv")
  qc <- tryCatch(read.csv(qc_path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(qc) || nrow(qc) == 0) {
    qc <- data.frame(comparison = comparison_name, status = "SKIPPED", runtime_seconds = 0, stringsAsFactors = FALSE)
  } else {
    qc$status <- "SKIPPED"
    if ("runtime_seconds" %in% colnames(qc)) qc$runtime_seconds <- 0
  }
  manifest <- reconstruct_checkpoint_manifest(dirs, comparison_name, data_path, config_path, ont, checkpoint_plot_suffix())
  list(status = "SKIPPED", comparison = comparison_name, error = NA_character_, qc = qc, manifest = manifest)
}

write_log_line <- function(log_file, level = "INFO", comparison = "GLOBAL", step = "GENERAL", message_text = "") {
  if (is.null(log_file) || !nzchar(log_file)) return(invisible(NULL))
  dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
  line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, comparison, step, message_text, sep = " | ")
  cat(line, "\n", file = log_file, append = TRUE)
  invisible(NULL)
}

print_progress_step <- function(comparison, step, detail = "", enabled = TRUE) {
  if (!isTRUE(enabled)) return(invisible(NULL))
  stamp <- format(Sys.time(), "%H:%M:%S")
  msg <- paste0("[", stamp, "] [", comparison, "] ", step,
                if (nzchar(detail)) paste0(" - ", detail) else "")
  cat(msg, "\n")
  flush.console()
  invisible(NULL)
}

prewarm_workers <- function(n_workers, runtime_params) {
  if (!isTRUE(runtime_params$prewarm_workers) || n_workers < 1) return(invisible(NULL))

  cat("\n==============================================\n")
  cat("PREWARMING WORKERS (INITIAL PACKAGE LOAD)\n")
  cat("==============================================\n\n")

  worker_ids <- seq_len(n_workers)
  if (isTRUE(runtime_params$show_progress)) {
    progressr::handlers("txtprogressbar")
    progressr::with_progress({
      p <- progressr::progressor(steps = length(worker_ids))
      future_lapply(worker_ids, function(i) {
        suppressPackageStartupMessages({
          library(clusterProfiler)
          library(pathview)
          library(enrichplot)
          library(DOSE)
          library(ggplot2)
          library(tidyr)
          library(dplyr)
          library(openxlsx)
          library(GOSemSim)
        })
        options(clusterProfiler.worker_packages_loaded = TRUE)
        p(message = paste0("worker ", i, " ready"))
        TRUE
      }, future.seed = TRUE)
    })
  } else {
    future_lapply(worker_ids, function(i) {
      suppressPackageStartupMessages({
        library(clusterProfiler)
        library(pathview)
        library(enrichplot)
        library(DOSE)
        library(ggplot2)
        library(tidyr)
        library(dplyr)
        library(openxlsx)
        library(GOSemSim)
      })
      options(clusterProfiler.worker_packages_loaded = TRUE)
      TRUE
    }, future.seed = TRUE)
  }

  cat("Worker prewarm complete. Starting analysis...\n\n")
  invisible(NULL)
}

checkpoint_file_path <- function(dirs) {
  file.path(dirs$results, "ANALYSIS_COMPLETE.flag")
}

checkpoint_plot_suffix <- function() {
  if (PERFORM_SIMPLIFICATION && USE_SIMPLIFIED_FOR_PLOTS) "_simplified" else "_full"
}

expected_checkpoint_outputs <- function(dirs, comparison_name, ont, plot_suffix = checkpoint_plot_suffix()) {
  c(
    flag = checkpoint_file_path(dirs),
    qc = file.path(dirs$results, "QC_summary.csv"),
    gsea_core_table = file.path(dirs$core_enrich_routed, paste0(comparison_name, plot_suffix, ".csv")),
    gsea_full_table = file.path(dirs$go_ont, paste0("GSEA_", ont, "_results_full.csv")),
    enrichgo_table = file.path(dirs$go_ont, "enrichGO_ALL_results_full.csv")
  )
}

qc_reports_complete <- function(qc_path) {
  if (!file.exists(qc_path)) return(FALSE)
  qc <- tryCatch(read.csv(qc_path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(qc) || nrow(qc) == 0 || !"status" %in% colnames(qc)) return(FALSE)
  as.character(qc$status[[1]]) %in% c("SUCCESS", "SKIPPED")
}

is_completed_checkpoint <- function(dirs, comparison_name, ont, plot_suffix = checkpoint_plot_suffix()) {
  expected <- expected_checkpoint_outputs(dirs, comparison_name, ont, plot_suffix)
  all(file.exists(expected)) && qc_reports_complete(expected[["qc"]])
}

write_completed_checkpoint <- function(dirs) {
  flag <- checkpoint_file_path(dirs)
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), file = flag)
  invisible(flag)
}

prepare_gene_vectors <- function(df, top_gene_abs_log2fc = 1) {
  colnames(df)[1] <- "gene_symbol"
  original_gene_list <- df$log2fc
  names(original_gene_list) <- df$gene_symbol

  gene_list <- sort(na.omit(original_gene_list), decreasing = TRUE)
  gene_list <- gene_list[!duplicated(names(gene_list))]

  top_gene_list <- gene_list
  top_genes <- names(top_gene_list)[abs(top_gene_list) > top_gene_abs_log2fc]
  if (length(top_genes) > 0) {
    top_genes <- names(sort(top_gene_list[top_genes], decreasing = TRUE))
  } else {
    top_genes <- character(0)
  }

  list(
    original_gene_list = original_gene_list,
    gene_list = gene_list,
    fc_for_cnet = gene_list,
    top_genes = top_genes
  )
}

load_background_universe <- function(cfg) {
  if (!isTRUE(cfg$background_universe$enabled)) return(character(0))
  file_path <- cfg$paths$background_universe_file
  if (!nzchar(file_path) || !file.exists(file_path)) {
    warning("Background universe enabled but file missing: ", file_path)
    return(character(0))
  }

  universe_df <- tryCatch(read.csv(file_path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(universe_df) || nrow(universe_df) == 0) {
    warning("Background universe file is empty or unreadable: ", file_path)
    return(character(0))
  }

  column_name <- cfg$background_universe$column
  if (!column_name %in% colnames(universe_df)) {
    warning("Background universe column not found: ", column_name)
    return(character(0))
  }

  universe <- unique(stats::na.omit(universe_df[[column_name]]))
  as.character(universe)
}

pretty_unit_name <- function(unit) {
  unit <- gsub("[^A-Za-z0-9]", "", as.character(unit))
  if (!nzchar(unit)) return("")
  m <- regexec("^([A-Za-z]+[0-9]*)([A-Za-z0-9]+)$", unit)
  mm <- regmatches(unit, m)[[1]]
  if (length(mm) == 3 && nzchar(mm[3])) {
    return(paste(toupper(mm[2]), tolower(mm[3]), sep = "_"))
  }
  unit
}

parse_sample_token <- function(token, phenotypes = c("sus", "res", "con")) {
  token <- trimws(as.character(token))
  phenotype <- ""
  unit <- token

  for (ph in phenotypes) {
    pattern <- paste0("(^|[_\\.-])", ph, "$|", ph, "$")
    if (grepl(pattern, token, ignore.case = TRUE)) {
      phenotype <- tolower(ph)
      unit <- sub(paste0("([_\\.-]?)", ph, "$"), "", token, ignore.case = TRUE)
      break
    }
  }

  unit <- gsub("[^A-Za-z0-9]", "", unit)
  pretty_unit <- pretty_unit_name(unit)

  list(raw = token, unit = toupper(unit), pretty_unit = pretty_unit, phenotype = phenotype)
}

split_comparison_tokens <- function(comparison_name) {
  comparison_name <- sub("\\.csv$", "", basename(as.character(comparison_name)))
  chars <- gregexpr("_", comparison_name, fixed = TRUE)[[1]]
  split_positions <- if (identical(chars, -1L)) integer(0) else as.integer(chars)
  candidates <- lapply(split_positions, function(pos) {
    left <- substr(comparison_name, 1, pos - 1)
    right <- substr(comparison_name, pos + 1, nchar(comparison_name))
    a <- parse_sample_token(left)
    b <- parse_sample_token(right)
    score <- sum(nzchar(c(a$phenotype, b$phenotype, a$unit, b$unit)))
    list(left = left, right = right, a = a, b = b, score = score)
  })
  valid <- candidates[vapply(candidates, function(x) {
    nzchar(x$a$phenotype) && nzchar(x$b$phenotype) && nzchar(x$a$unit) && nzchar(x$b$unit)
  }, logical(1))]
  if (length(valid) > 0) {
    return(valid[[which.max(vapply(valid, `[[`, numeric(1), "score"))]])
  }
  if (length(candidates) > 0) {
    return(candidates[[which.max(vapply(candidates, `[[`, numeric(1), "score"))]])
  }
  empty <- parse_sample_token(comparison_name)
  list(left = comparison_name, right = "", a = empty, b = parse_sample_token(""), score = 0)
}

classify_comparison_route <- function(comparison_name) {
  tokens <- split_comparison_tokens(comparison_name)
  if (!nzchar(tokens$right)) {
    return(list(category = "unclassified", unit_folder = "unknown_unit"))
  }

  a <- tokens$a
  b <- tokens$b

  same_unit <- nzchar(a$unit) && nzchar(b$unit) && identical(a$unit, b$unit)
  same_pheno <- nzchar(a$phenotype) && nzchar(b$phenotype) && identical(a$phenotype, b$phenotype)

  category <- if (same_unit && !same_pheno) {
    "phenotype_within_unit"
  } else if (!same_unit && same_pheno) {
    "phenotype_between_unit"
  } else if (same_unit && same_pheno) {
    "within_unit_same_phenotype"
  } else if (!same_unit && !same_pheno) {
    "between_unit_and_phenotype"
  } else {
    "unclassified"
  }

  unit_folder <- if (same_unit && nzchar(a$pretty_unit)) {
    a$pretty_unit
  } else {
    unit_a <- if (nzchar(a$pretty_unit)) a$pretty_unit else "unitA"
    unit_b <- if (nzchar(b$pretty_unit)) b$pretty_unit else "unitB"
    paste(sort(c(unit_a, unit_b)), collapse = "_vs_")
  }

  list(
    category = category,
    unit_folder = unit_folder,
    left_token = tokens$left,
    right_token = tokens$right,
    left_unit = a$pretty_unit,
    right_unit = b$pretty_unit,
    left_phenotype = a$phenotype,
    right_phenotype = b$phenotype
  )
}

comparison_route_qc_row <- function(comparison_name) {
  route <- classify_comparison_route(comparison_name)
  expected_warning <- if (identical(route$category, "unclassified") ||
                          !nzchar(route$left_phenotype) ||
                          !nzchar(route$right_phenotype)) {
    "unclassified_or_missing_phenotype"
  } else {
    ""
  }
  data.frame(
    comparison = comparison_name,
    left_token = route$left_token %||% NA_character_,
    right_token = route$right_token %||% NA_character_,
    left_unit = route$left_unit %||% NA_character_,
    right_unit = route$right_unit %||% NA_character_,
    left_phenotype = route$left_phenotype %||% NA_character_,
    right_phenotype = route$right_phenotype %||% NA_character_,
    route_category = route$category,
    route_unit = route$unit_folder,
    expected_warning = expected_warning,
    stringsAsFactors = FALSE
  )
}

assert_comparison_route_examples <- function() {
  expected <- data.frame(
    comparison = c(
      "CA1slmres_CA1slmcon",
      "CA1slmsus_CA1slmres",
      "CA1slmsus_CA1slmcon",
      "CA1sores_CA1socon",
      "CA1slmcon_CA1socon",
      "CA1slmsus_CA1socon",
      "CA1_so_res_CA1_so_con",
      "CA1_slm_con_CA1_so_con"
    ),
    category = c(
      "phenotype_within_unit",
      "phenotype_within_unit",
      "phenotype_within_unit",
      "phenotype_within_unit",
      "phenotype_between_unit",
      "between_unit_and_phenotype",
      "phenotype_within_unit",
      "phenotype_between_unit"
    ),
    unit = c(
      "CA1_slm",
      "CA1_slm",
      "CA1_slm",
      "CA1_so",
      "CA1_slm_vs_CA1_so",
      "CA1_slm_vs_CA1_so",
      "CA1_so",
      "CA1_slm_vs_CA1_so"
    ),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(expected))) {
    route <- classify_comparison_route(expected$comparison[[i]])
    if (!identical(route$category, expected$category[[i]]) ||
        !identical(route$unit_folder, expected$unit[[i]])) {
      stop(
        "Comparison route assertion failed for ", expected$comparison[[i]],
        ": got ", route$category, " / ", route$unit_folder,
        "; expected ", expected$category[[i]], " / ", expected$unit[[i]],
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

# ----------------------------------------------------
# 3. DATA INPUT/COMPARISONS
# ----------------------------------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x
config_candidates <- c(
  file.path(getwd(), "clusterProfiler_config.yml"),
  file.path(getwd(), "clusterProfiler_config.local.yml"),
  file.path(getwd(), "config", "clusterProfiler_config.yml"),
  file.path(getwd(), "config", "clusterProfiler_config.local.yml"),
  repo_path("config", "clusterProfiler_config.local.yml"),
  repo_path("config", "clusterProfiler_config.yml"),
  file.path(getwd(), "Analysis", "clusterProfiler_config.yml")
)
config_path <- config_candidates[file.exists(config_candidates)][1] %||% config_candidates[1]
cfg <- read_config(config_path)
assert_comparison_route_examples()

as_repo_path <- function(path) {
  if (is.null(path) || !nzchar(path)) return(path)
  if (grepl("^([A-Za-z]:|/|~)", path)) return(path)
  repo_path(path)
}
DATASET <- as.character(cfg$analysis$dataset %||% Sys.getenv("PROTEOMICS_COMPARISON", unset = "neuron_neuropil"))
allowed_datasets <- trimws(strsplit(
  Sys.getenv("PROTEOMICS_ALLOWED_COMPARISONS", unset = "neuron_neuropil,neuron_soma,microglia"),
  ",",
  fixed = TRUE
)[[1]])
allowed_datasets <- allowed_datasets[nzchar(allowed_datasets)]
if (length(allowed_datasets) > 0 && !DATASET %in% allowed_datasets) {
  stop(
    "clusterProfiler dataset='", DATASET, "' is not in the allowed dataset families: ",
    paste(allowed_datasets, collapse = ", "),
    ". Extend PROTEOMICS_ALLOWED_COMPARISONS if this is an intentional new family.",
    call. = FALSE
  )
}
CANONICAL_PATHS <- lapply(CANONICAL_PATHS, function(path) file.path(path, DATASET))
invisible(lapply(CANONICAL_PATHS, dir.create, recursive = TRUE, showWarnings = FALSE))

dataset_mapped_default <- path_processed("02_id_mapping", "mapped", DATASET, "forward", "per_file")
if (is.null(cfg$paths$mapped_dir) || !nzchar(cfg$paths$mapped_dir)) {
  cfg$paths$mapped_dir <- dataset_mapped_default
}
if (is.null(cfg$paths$mapped_data_base) || !nzchar(cfg$paths$mapped_data_base)) {
  cfg$paths$mapped_data_base <- cfg$paths$mapped_dir
}
cfg$paths$mapped_dir <- as_repo_path(cfg$paths$mapped_dir)
cfg$paths$working_base <- as_repo_path(cfg$paths$working_base)
cfg$paths$mapped_data_base <- as_repo_path(cfg$paths$mapped_data_base)
cfg$paths$background_universe_file <- as_repo_path(cfg$paths$background_universe_file)
DRY_RUN <- is_dry_run(cfg)

legacy_mapped_dir <- normalizePath(path_processed("02_id_mapping", "mapped", "forward", "per_file"), winslash = "/", mustWork = FALSE)
configured_mapped_dir <- normalizePath(cfg$paths$mapped_dir, winslash = "/", mustWork = FALSE)
dataset_mapped_dir <- normalizePath(dataset_mapped_default, winslash = "/", mustWork = FALSE)
if (identical(configured_mapped_dir, legacy_mapped_dir)) {
  stop(
    "Refusing to scan legacy mixed mapped directory: ", cfg$paths$mapped_dir,
    "\nUse dataset-specific mapped_dir: ", dataset_mapped_default,
    "\nor set analysis.dataset/paths.mapped_dir explicitly for an intentional legacy migration.",
    call. = FALSE
  )
}

PERFORM_SIMPLIFICATION <- isTRUE(cfg$simplification$perform)
USE_SIMPLIFIED_FOR_PLOTS <- isTRUE(cfg$simplification$use_for_plots)
SIMPLIFY_CUTOFF <- as.numeric(cfg$simplification$cutoff)

mapped_dir <- cfg$paths$mapped_dir
mapped_dir_exists <- dir.exists(mapped_dir)
if (!mapped_dir_exists && !isTRUE(DRY_RUN)) {
  stop("Mapped data directory does not exist: ", mapped_dir)
}

comparison_files <- if (mapped_dir_exists) list.files(mapped_dir, pattern = "\\.csv$", full.names = FALSE) else character(0)
comparison_names <- sub("\\.csv$", "", comparison_files)
comparison_list <- lapply(comparison_names, function(comparison_name) {
  tokens <- split_comparison_tokens(comparison_name)
  if (nzchar(tokens$right)) c(tokens$left, tokens$right) else NULL
})
comparison_list <- Filter(Negate(is.null), comparison_list)

if (length(comparison_list) == 0 && !isTRUE(DRY_RUN)) {
  stop("No valid comparison files found in: ", mapped_dir)
}

working_base <- cfg$paths$working_base

working_dir <- working_base
mapped_data_base <- cfg$paths$mapped_data_base
#mapped_data_base <- file.path(working_base, "Datasets/mapped/baseline_cell_type_profiling/US")
organism <- cfg$analysis$organism
ont <- cfg$analysis$ontology

analysis_params <- list(
  pvalue_cutoff = as.numeric(cfg$analysis$pvalue_cutoff),
  qvalue_cutoff = as.numeric(cfg$analysis$qvalue_cutoff),
  p_adjust_method = as.character(cfg$analysis$p_adjust_method),
  min_gs_size = as.integer(cfg$analysis$min_gs_size),
  max_gs_size = as.integer(cfg$analysis$max_gs_size),
  top_gene_abs_log2fc = as.numeric(cfg$analysis$top_gene_abs_log2fc)
)

runtime_params <- list(
  workers = as.integer(cfg$runtime$workers),
  future_globals_max_size_mb = as.numeric(cfg$runtime$future_globals_max_size_mb),
  resume_if_complete = isTRUE(cfg$runtime$resume_if_complete),
  force_rerun = isTRUE(cfg$runtime$force_rerun),
  prewarm_workers = isTRUE(cfg$runtime$prewarm_workers),
  show_progress = isTRUE(cfg$runtime$show_progress),
  show_step_progress = isTRUE(cfg$runtime$show_step_progress),
  dry_run = isTRUE(cfg$runtime$dry_run)
)

if (isTRUE(DRY_RUN)) {
  expected_dirs <- unlist(CANONICAL_PATHS, use.names = TRUE)
  route_qc <- if (length(comparison_names) > 0) {
    dplyr::bind_rows(lapply(comparison_names, comparison_route_qc_row))
  } else {
    data.frame()
  }
  diagnostics <- data.frame(
    check = c("config_exists", "dataset", "mapped_dir_exists", "comparison_files_found", "canonical_dirs"),
    status = c(
      if (file.exists(config_path)) "PASS" else "WARN",
      "PASS",
      if (mapped_dir_exists) "PASS" else "FAIL",
      if (length(comparison_files) > 0) "PASS" else "FAIL",
      "PASS"
    ),
    detail = c(
      config_path,
      DATASET,
      mapped_dir,
      paste0(length(comparison_files), " csv file(s)"),
      paste(names(expected_dirs), expected_dirs, sep = "=", collapse = "; ")
    ),
    stringsAsFactors = FALSE
  )
  dry_run_line("Script", "01_clusterProfiler.r")
  dry_run_line("Config", config_path, diagnostics$status[1])
  dry_run_line("Dataset", DATASET)
  dry_run_line("Mapped directory", mapped_dir, diagnostics$status[3])
  dry_run_line("Comparison CSV count", length(comparison_files), diagnostics$status[4])
  if (length(comparison_files) > 0) {
    dry_run_line("First comparisons", paste(head(comparison_files, 10), collapse = ", "))
    route_counts <- as.data.frame(table(route_qc$route_category, route_qc$route_unit), stringsAsFactors = FALSE)
    names(route_counts) <- c("route_category", "route_unit", "n")
    route_counts <- route_counts[route_counts$n > 0, , drop = FALSE]
    dry_run_line("First routed comparisons", paste(utils::capture.output(print(head(route_qc, 10), row.names = FALSE)), collapse = " | "))
    dry_run_line("Route counts", paste(utils::capture.output(print(route_counts, row.names = FALSE)), collapse = " | "))
  }
  dry_run_line("Ontology", cfg$analysis$ontology)
  dry_run_line("Result manifest", file.path(CANONICAL_PATHS$processed, "clusterProfiler_manifest.csv"))
  dry_run_file <- file.path(CANONICAL_PATHS$reports, "clusterProfiler_dry_run_diagnostics.csv")
  write.csv(diagnostics, dry_run_file, row.names = FALSE)
  route_qc_file <- file.path(CANONICAL_PATHS$reports, "clusterProfiler_route_qc_dry_run.csv")
  write.csv(route_qc, route_qc_file, row.names = FALSE)
  dry_run_line("Route QC written", route_qc_file)
  dry_run_line("Diagnostics written", dry_run_file)
  quit(status = if (any(diagnostics$status == "FAIL")) 1 else 0, save = "no")
}

setupPackages()

run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_log_dir <- CANONICAL_PATHS$logs
dir.create(run_log_dir, recursive = TRUE, showWarnings = FALSE)
master_log <- file.path(run_log_dir, paste0("clusterProfiler_run_", run_id, ".log"))
write_log_line(master_log, "INFO", "GLOBAL", "CONFIG", paste0("Using config file: ", config_path))
write_log_line(master_log, "INFO", "GLOBAL", "DATASET", paste0("Dataset family: ", DATASET))
write_session_info(file.path(CANONICAL_PATHS$logs, "sessionInfo.txt"))
write_config_snapshot(cfg, file.path(CANONICAL_PATHS$logs, paste0("clusterProfiler_config_snapshot_", run_id, ".yml")))

route_qc <- dplyr::bind_rows(lapply(comparison_names, comparison_route_qc_row))
route_qc$dataset <- DATASET
route_qc_file <- file.path(CANONICAL_PATHS$reports, paste0("clusterProfiler_route_qc_", run_id, ".csv"))
write.csv(route_qc, route_qc_file, row.names = FALSE)
route_counts <- as.data.frame(table(route_qc$route_category, route_qc$route_unit), stringsAsFactors = FALSE)
names(route_counts) <- c("route_category", "route_unit", "n")
route_counts <- route_counts[route_counts$n > 0, , drop = FALSE]
cat("Route QC written to:", route_qc_file, "\n")
print(head(route_qc, 10), row.names = FALSE)
print(route_counts, row.names = FALSE)

#nk3r_genes <- c("P21279", "P21278", "P51432", "P11881", "P63318", "P68404", "P0DP26", "P0DP27", "P0DP28", "P11798", "P28652", "P47937", "P47713")
#selected_uniprot <- c("P21279", "P21278", "Q9Z1B3", "P51432", "P11881", "P68404", "P63318", "P0DP26", "P0DP27", "P11798", "P28652", "Q61411", "Q99N57", "P31938", "P63085", "Q63844", "Q8BWG8", "Q91YI4", "V9GXQ9")
#path_ids <- c("mmu04110", "mmu04115", "mmu04114", "mmu04113", "mmu04112", "mmu04111", "mmu04116", "mmu04117", "mmu04118", "mmu04119", "mmu04720", "mmu04721", "mmu04722", "mmu04725", "mmu04726", "mmu04727", "mmu04724", "mmu04080", "mmu00030", "mmu04151")

# Optional analysis inputs: default to empty vectors so downstream blocks can be skipped cleanly.
if (!exists("nk3r_genes", inherits = FALSE)) {
  nk3r_genes <- as.character(cfg$optional_inputs$nk3r_genes %||% character(0))
  if (length(nk3r_genes) > 0) {
    message("Loaded nk3r_genes from config (n=", length(nk3r_genes), ").")
  } else {
    message("No nk3r_genes provided; custom NK3R GSEA will be skipped.")
  }
}
if (!exists("selected_uniprot", inherits = FALSE)) {
  selected_uniprot <- as.character(cfg$optional_inputs$selected_uniprot %||% character(0))
  if (length(selected_uniprot) > 0) {
    message("Loaded selected_uniprot from config (n=", length(selected_uniprot), ").")
  } else {
    message("No selected_uniprot provided; predefined KEGG GSEA will be skipped.")
  }
}
if (!exists("path_ids", inherits = FALSE)) {
  path_ids <- as.character(cfg$optional_inputs$path_ids %||% character(0))
  if (length(path_ids) > 0) {
    message("Loaded path_ids from config (n=", length(path_ids), ").")
  } else {
    message("No path_ids provided; pathview rendering will be skipped.")
  }
}

background_universe <- load_background_universe(cfg)
if (length(background_universe) > 0) {
  write_log_line(master_log, "INFO", "GLOBAL", "BACKGROUND", paste0("Loaded background universe: n=", length(background_universe)))
} else {
  write_log_line(master_log, "INFO", "GLOBAL", "BACKGROUND", "No background universe loaded (disabled or unavailable).")
}

completed_results <- list()
if (isTRUE(runtime_params$resume_if_complete) && !isTRUE(runtime_params$force_rerun)) {
  complete_flags <- vapply(comparison_list, function(cell_types) {
    comparison_name <- paste(cell_types, collapse = "_")
    dirs <- create_analysis_dirs(working_base, comparison_name, ont)
    is_completed_checkpoint(dirs, comparison_name, ont, checkpoint_plot_suffix())
  }, logical(1))

  if (any(complete_flags)) {
    completed_results <- lapply(comparison_list[complete_flags], make_checkpoint_skip_result,
                                working_base = working_base,
                                mapped_data_base = mapped_data_base,
                                ont = ont,
                                config_path = config_path)
    skipped_names <- vapply(completed_results, function(x) x$comparison, character(1))
    cat("Resume check: skipping ", length(completed_results), " completed comparison(s).\n", sep = "")
    cat("Completed comparisons: ", paste(skipped_names, collapse = ", "), "\n", sep = "")
    write_log_line(master_log, "INFO", "GLOBAL", "RESUME",
                   paste0("Skipping completed comparisons: ", paste(skipped_names, collapse = ", ")))
  }

  comparison_list <- comparison_list[!complete_flags]
  cat("Resume check: ", length(comparison_list), " comparison(s) still need analysis.\n", sep = "")
} else if (isTRUE(runtime_params$force_rerun)) {
  cat("Resume check disabled: force_rerun=TRUE, all comparisons will be recomputed.\n")
}

if (length(comparison_list) > 0) {

# ----------------------------------------------------
# 4. SETUP PARALLEL PROCESSING WITH FUTURE
# ----------------------------------------------------
suppressPackageStartupMessages({
  library(future)
  library(future.apply)
})

# Ensure progress handlers are active in interactive runs
options(progressr.enable = TRUE)
progressr::handlers(global = TRUE)

# Set up multisession plan
# Set up multisession plan
detected_cores <- suppressWarnings(as.integer(availableCores()))
cat("[DEBUG PARALLEL] availableCores() returned:", detected_cores, "\n")
if (is.na(detected_cores) || detected_cores < 1L) detected_cores <- 1L
cat("[DEBUG PARALLEL] After validation, detected_cores =", detected_cores, "\n")
cat("[DEBUG PARALLEL] runtime_params$workers =", runtime_params$workers, "\n")

if (!is.na(runtime_params$workers) && runtime_params$workers > 0) {
  # Positive value: use specified workers, capped at detected cores
  n_cores <- min(runtime_params$workers, detected_cores)
  cat("[DEBUG PARALLEL] Using workers from config:", runtime_params$workers, "-> n_cores =", n_cores, "\n")
} else if (!is.na(runtime_params$workers) && runtime_params$workers < 0) {
  # Negative value (-1): use all cores except the specified number
  # e.g., -1 means use all but 1
  n_cores <- max(1L, detected_cores + runtime_params$workers)
  cat("[DEBUG PARALLEL] Using negative config value:", runtime_params$workers, "-> n_cores =", n_cores, "\n")
} else {
  # Default fallback: use detected cores - 1
  n_cores <- max(1L, detected_cores - 1L)
  cat("[DEBUG PARALLEL] Using fallback (detected - 1):", n_cores, "\n")
}

cat("[DEBUG PARALLEL] Calling plan(multisession, workers =", n_cores, ")\n")
plan(multisession, workers = n_cores)
cat("[DEBUG PARALLEL] Plan set. Current plan is:", toString(class(plan())), "\n")
on.exit(plan(sequential), add = TRUE)

# Increase global size limit for passing data to workers
options(future.globals.maxSize = runtime_params$future_globals_max_size_mb * 1024^2) 

cat("Setting up parallel processing with", n_cores, "cores using future\n")
write_log_line(master_log, "INFO", "GLOBAL", "PARALLEL", paste0("workers=", n_cores,
                                                                      ", future.globals.maxSize(MB)=", runtime_params$future_globals_max_size_mb))

prewarm_workers(n_cores, runtime_params)

# ----------------------------------------------------
# 5. ANALYSIS FUNCTION
# ----------------------------------------------------
analyze_comparison <- function(cell_types, working_base, mapped_data_base, organism, ont, 
                               nk3r_genes, selected_uniprot, path_ids,
                               analysis_params, runtime_params, run_log_dir,
                               background_universe = character(0),
                               config_path = NA_character_) {

  # Load required libraries in worker once per worker session.
  if (!isTRUE(getOption("clusterProfiler.worker_packages_loaded", FALSE))) {
    suppressPackageStartupMessages({
      library(clusterProfiler)
      library(pathview)
      library(enrichplot)
      library(DOSE)
      library(ggplot2)
      library(tidyr)
      library(dplyr)
      library(openxlsx)
      library(GOSemSim)
    })
    options(clusterProfiler.worker_packages_loaded = TRUE)
  }

  # === CRITICAL FIX FOR PARALLEL PROCESSING ===
  # To avoid "undefined columns selected" / connection errors:
  # We must ensure the OrgDb is loaded FRESH in this worker process.
  # If it was inherited from the master process, the SQLite connection is dead.
  comparison_name <- paste(cell_types, collapse = "_")
  comparison_log <- file.path(run_log_dir, paste0(comparison_name, ".log"))
  run_start <- Sys.time()
  qc <- data.frame(
    comparison = comparison_name,
    status = "STARTED",
    runtime_seconds = NA_real_,
    n_input_rows = NA_integer_,
    n_non_na_log2fc = NA_integer_,
    n_unique_gene_ids = NA_integer_,
    n_top_genes = NA_integer_,
    n_background_universe = length(background_universe),
    n_id_mapped_uniprot_to_entrez = NA_integer_,
    n_gsea_terms = NA_integer_,
    n_ora_terms = NA_integer_,
    n_kegg_terms = NA_integer_,
    stringsAsFactors = FALSE
  )
  manifest_rows <- list()

  cat("Analyzing comparison:", comparison_name, "\n")
  write_log_line(comparison_log, "INFO", comparison_name, "START", "Comparison analysis started")
  print_progress_step(comparison_name, "START", "Worker initialized", runtime_params$show_step_progress)

  tryCatch({
    # Load OrgDb namespace quietly (no attach chatter in console)
    if (!requireNamespace(organism, quietly = TRUE)) {
      stop("Could not load OrgDb namespace: ", organism)
    }

    # Get object reference from namespace
    org_db_obj <- get(organism, envir = asNamespace(organism))

    # Setup directories
    dirs <- create_analysis_dirs(working_base, comparison_name, ont)
    cat("Directories set up for comparison:", comparison_name, "\n")
    route_msg <- paste0("route=", dirs$route_category, " / ", dirs$route_unit)
    print_progress_step(comparison_name, "ROUTE", route_msg, runtime_params$show_step_progress)
    write_log_line(comparison_log, "INFO", comparison_name, "ROUTE", route_msg)
    print_progress_step(comparison_name, "SETUP", "Output directories ready", runtime_params$show_step_progress)
    qc_path <- file.path(dirs$results, "QC_summary.csv")

    # Define data file path before checkpoint handling so skipped comparisons
    # can still reconstruct manifest rows from existing canonical outputs.
    file_name <- paste0(comparison_name, ".csv")
    cat("Looking for data file:", file_name, "\n")
    data_path <- file.path(mapped_data_base, file_name)
    cat("Full data path:", data_path, "\n")

    if (isTRUE(runtime_params$resume_if_complete) &&
        !isTRUE(runtime_params$force_rerun) &&
        is_completed_checkpoint(dirs, comparison_name, ont, checkpoint_plot_suffix())) {
      write_log_line(comparison_log, "INFO", comparison_name, "CHECKPOINT", "Checkpoint found; skipping comparison")
      manifest_rows <- list(reconstruct_checkpoint_manifest(dirs, comparison_name, data_path, config_path, ont, plot_suffix = checkpoint_plot_suffix()))
      qc$status <- "SKIPPED"
      qc$runtime_seconds <- as.numeric(difftime(Sys.time(), run_start, units = "secs"))
      write.csv(qc, qc_path, row.names = FALSE)
      return(list(status = "SKIPPED", comparison = comparison_name, error = NA_character_, qc = qc, manifest = dplyr::bind_rows(manifest_rows)))
    }

    if (!file.exists(data_path)) {
      write_log_line(comparison_log, "ERROR", comparison_name, "INPUT", paste0("File not found: ", data_path))
      qc$status <- "FAILED"
      qc$runtime_seconds <- as.numeric(difftime(Sys.time(), run_start, units = "secs"))
      write.csv(qc, qc_path, row.names = FALSE)
      return(list(status = "FAILED", comparison = comparison_name, error = "File not found", qc = qc, manifest = dplyr::bind_rows(manifest_rows)))
    }

    # Load and prepare data
    df <- read.csv(data_path, header = TRUE, stringsAsFactors = FALSE)
    cat("Data loaded for comparison:", comparison_name, "\n")
    print_progress_step(comparison_name, "INPUT", paste0("Loaded rows=", nrow(df)), runtime_params$show_step_progress)

    if (nrow(df) == 0) {
      write_log_line(comparison_log, "ERROR", comparison_name, "INPUT", "Data file is empty")
      qc$status <- "FAILED"
      qc$runtime_seconds <- as.numeric(difftime(Sys.time(), run_start, units = "secs"))
      write.csv(qc, qc_path, row.names = FALSE)
      return(list(status = "FAILED", comparison = comparison_name, error = "Data file is empty", qc = qc, manifest = dplyr::bind_rows(manifest_rows)))
    }

    if (!"log2fc" %in% colnames(df)) {
      if ("logFC" %in% colnames(df)) {
        colnames(df)[colnames(df) == "logFC"] <- "log2fc"
      } else {
        write_log_line(comparison_log, "ERROR", comparison_name, "INPUT", "No log2fc/logFC column")
        qc$status <- "FAILED"
        qc$runtime_seconds <- as.numeric(difftime(Sys.time(), run_start, units = "secs"))
        write.csv(qc, qc_path, row.names = FALSE)
        return(list(status = "FAILED", comparison = comparison_name, error = "No log2fc column", qc = qc, manifest = dplyr::bind_rows(manifest_rows)))
      }
    }

    vectors <- prepare_gene_vectors(df, top_gene_abs_log2fc = analysis_params$top_gene_abs_log2fc)
    original_gene_list <- vectors$original_gene_list
    gene_list <- vectors$gene_list
    fc_for_cnet <- vectors$fc_for_cnet
    top_genes <- vectors$top_genes

    qc$n_input_rows <- nrow(df)
    qc$n_non_na_log2fc <- sum(!is.na(df$log2fc))
    qc$n_unique_gene_ids <- length(unique(df[[1]]))
    qc$n_top_genes <- length(top_genes)
    write_log_line(comparison_log, "INFO", comparison_name, "QC", paste0("rows=", qc$n_input_rows,
                                                                            ", nonNA_log2fc=", qc$n_non_na_log2fc,
                                                                            ", top_genes=", qc$n_top_genes))
    print_progress_step(comparison_name, "QC", paste0("top_genes=", qc$n_top_genes), runtime_params$show_step_progress)

    # ----------------------------------------------------
    # GSEA (GO) WITH OPTIONAL SIMPLIFICATION
    # ----------------------------------------------------
    # Ensure OrgDb is passed as object
    gse <- gseGO(geneList = gene_list, ont = ont, keyType = "UNIPROT", 
                 minGSSize = analysis_params$min_gs_size,
                 maxGSSize = analysis_params$max_gs_size,
                 pvalueCutoff = analysis_params$pvalue_cutoff,
                 verbose = FALSE,
                 OrgDb = org_db_obj,
                 pAdjustMethod = analysis_params$p_adjust_method)

    # Perform simplification ONLY if requested
    gse_simplified <- NULL  # Initialize as NULL

    if (PERFORM_SIMPLIFICATION && !is.null(gse) && nrow(gse@result) > 0) {
      gse_temp <- tryCatch({
        simplify(gse, cutoff = SIMPLIFY_CUTOFF, by = "p.adjust", select_fun = min)
      }, error = function(e) {
        warning("Simplify failed for ", comparison_name, ": ", e$message)
        NULL  # Return NULL on failure
      })

      # Validate S4 class
      if (!is.null(gse_temp) && (is(gse_temp, "gseaResult") || is(gse_temp, "enrichResult"))) {
        gse_simplified <- gse_temp
      }
    }

    # Determine which version to use for plotting
    if (PERFORM_SIMPLIFICATION && !is.null(gse_simplified)) {
      gse_for_plot <- if(USE_SIMPLIFIED_FOR_PLOTS) gse_simplified else gse
      plot_suffix <- if(USE_SIMPLIFIED_FOR_PLOTS) "_simplified" else "_full"
    } else {
      # No simplification performed - use full results only
      gse_for_plot <- gse
      plot_suffix <- "_full"
    }

    # Generate plots using selected version
    if (!is.null(gse_for_plot) && nrow(gse_for_plot@result) > 0) {
      plot_label <- if(PERFORM_SIMPLIFICATION && USE_SIMPLIFIED_FOR_PLOTS && !is.null(gse_simplified)) {
        "(Simplified)"
      } else {
        "(Full)"
      }

      tryCatch({
        gse_plot <- clusterProfiler::dotplot(gse_for_plot, showCategory = 10, split = ".sign") +
          facet_wrap(~ .sign, nrow = 1) +
          labs(title = paste("GSEA", ont, "of", comparison_name, plot_label), 
               x = "Gene Ratio", y = "Gene Set") +
          scale_fill_viridis_c(option = "cividis") + 
          theme_minimal(base_size = 12) +
          theme(
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = 10),
            strip.text = element_text(size = 12, face = "plain"),
            panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
            panel.grid.minor = element_blank()
          )
        save_plot_organized(gse_plot, paste0("GSEA_", ont, "_dotplot", plot_suffix, ".svg"), dirs$plots_go)
      }, error=function(e) warning("GSEA dotplot failed: ", e$message))

      # Other plots
      tryCatch({
        save_plot_organized(emapplot(pairwise_termsim(gse_for_plot), showCategory = 10), 
                            paste0("GSEA_", ont, "_emap", plot_suffix, ".svg"), dirs$plots_go)
      }, error=function(e) NULL)

      tryCatch({
        save_plot_organized(cnetplot(gse_for_plot, categorySize = "pvalue", foldChange = fc_for_cnet), 
                            paste0("GSEA_", ont, "_cnet", plot_suffix, ".svg"), dirs$plots_go)
      }, error = function(e) warning("cnetplot failed: ", e$message))

      tryCatch({
        ridgeplot_gse <- ridgeplot(gse_for_plot) + 
          labs(x = "Enrichment Distribution", 
               title = paste("GSEA", ont, "Ridgeplot", plot_label)) + 
          theme_minimal()
        save_plot_organized(ridgeplot_gse, paste0("GSEA_", ont, "_ridge", plot_suffix, ".svg"), dirs$plots_go)
      }, error=function(e) NULL)

      tryCatch({
        gseaplot_gse <- gseaplot(gse_for_plot, by = "all", 
                                 title = gse_for_plot@result$Description[1], geneSetID = 1)
        save_plot_organized(gseaplot_gse, paste0("GSEA_", ont, "_plot", plot_suffix, ".svg"), dirs$plots_go)
      }, error=function(e) NULL)

      tryCatch({
        top_terms <- head(gse_for_plot@result$Description, 3)
        pmcplot_gse <- pmcplot(top_terms, 2010:2025, proportion = FALSE) + 
          labs(title = paste("Publication Trends -", ont, plot_label))
        save_plot_organized(pmcplot_gse, paste0("GSEA_", ont, "_pubmed", plot_suffix, ".svg"), dirs$plots_go)
      }, error=function(e) NULL)
    }

    # Save results based on whether simplification was performed
    write.csv(gse@result, 
              file = file.path(dirs$go_ont, paste0("GSEA_", ont, "_results_full.csv")), 
              row.names = FALSE)

    if (PERFORM_SIMPLIFICATION && !is.null(gse_simplified)) {
      write.csv(gse_simplified@result, 
                file = file.path(dirs$go_ont, paste0("GSEA_", ont, "_results_simplified.csv")), 
                row.names = FALSE)
    }

    # Save the version used for plots to core_enrich
    gse_for_export <- if(PERFORM_SIMPLIFICATION && USE_SIMPLIFIED_FOR_PLOTS && !is.null(gse_simplified)) {
      gse_simplified
    } else {
      gse
    }

    write.csv(gse_for_export@result, 
              file = file.path(dirs$core_enrich, paste0(comparison_name, plot_suffix, ".csv")), 
              row.names = FALSE)
    write.csv(gse_for_export@result,
          file = file.path(dirs$core_enrich_routed, paste0(comparison_name, plot_suffix, ".csv")),
          row.names = FALSE)
    gsea_core_table <- file.path(dirs$core_enrich_routed, paste0(comparison_name, plot_suffix, ".csv"))
    gsea_plot <- file.path(dirs$plots_go, paste0("GSEA_", ont, "_dotplot", plot_suffix, ".svg"))
    manifest_rows[[length(manifest_rows) + 1]] <- make_manifest_row(
      result_type = "GSEA_GO",
      ontology = ont,
      comparison_name = comparison_name,
      dirs = dirs,
      input_gene_file = data_path,
      config_path = config_path,
      output_table = gsea_core_table,
      output_plot = if (file.exists(gsea_plot)) gsea_plot else NA_character_,
      n_genes = length(gene_list),
      n_terms = nrow(gse_for_export@result),
      simplified = identical(plot_suffix, "_simplified"),
      used_for_plot = TRUE,
      plot_suffix = plot_suffix
    )
    qc$n_gsea_terms <- nrow(gse@result)
    print_progress_step(comparison_name, "GSEA", paste0("terms=", qc$n_gsea_terms), runtime_params$show_step_progress)

    # ----------------------------------------------------
    # ORA WITH OPTIONAL SIMPLIFICATION
    # ----------------------------------------------------
    ora_universe <- if (length(background_universe) > 0) background_universe else names(gene_list)
    if(length(top_genes) > 0) {
      ora <- enrichGO(gene = top_genes, ont = ont, keyType = "UNIPROT", 
                      universe = ora_universe,
                      minGSSize = analysis_params$min_gs_size,
                      maxGSSize = analysis_params$max_gs_size,
                      pvalueCutoff = analysis_params$pvalue_cutoff,
                      OrgDb = org_db_obj,
                      pAdjustMethod = analysis_params$p_adjust_method)

      # Simplify ONLY if requested
      ora_simplified <- NULL

      if (PERFORM_SIMPLIFICATION && nrow(ora@result) > 0) {
        ora_simplified <- tryCatch({
          simplify(ora, cutoff = SIMPLIFY_CUTOFF, by = "p.adjust", select_fun = min)
        }, error = function(e) {
          warning("ORA simplify failed for ", comparison_name, ": ", e$message)
          NULL
        })

        # Validate S4 class
        if (!is.null(ora_simplified) && !is(ora_simplified, "enrichResult")) {
          ora_simplified <- NULL
        }
      }

      # Determine version for plotting
      if (PERFORM_SIMPLIFICATION && !is.null(ora_simplified)) {
        ora_for_plot <- if(USE_SIMPLIFIED_FOR_PLOTS) ora_simplified else ora
      } else {
        ora_for_plot <- ora
      }

      if (nrow(ora_for_plot@result) > 0) {
        plot_label <- if(PERFORM_SIMPLIFICATION && USE_SIMPLIFIED_FOR_PLOTS && !is.null(ora_simplified)) {
          "(Simplified)"
        } else {
          "(Full)"
        }

        ora_plot <- clusterProfiler::dotplot(ora_for_plot, showCategory = 10) + 
          labs(title = paste("ORA", ont, "- Top Regulated Genes", plot_label)) + 
          scale_x_continuous(limits = c(0, 1))
        save_plot_organized(ora_plot, paste0("ORA_", ont, "_dotplot", plot_suffix, ".svg"), dirs$plots_ora)
      }

      # Save results
      write.csv(ora@result, 
                file = file.path(dirs$ora, paste0("ORA_", ont, "_results_full.csv")), 
                row.names = FALSE)

      if (PERFORM_SIMPLIFICATION && !is.null(ora_simplified)) {
        write.csv(ora_simplified@result, 
                  file = file.path(dirs$ora, paste0("ORA_", ont, "_results_simplified.csv")), 
                  row.names = FALSE)
      }
      qc$n_ora_terms <- nrow(ora@result)
      print_progress_step(comparison_name, "ORA", paste0("terms=", qc$n_ora_terms), runtime_params$show_step_progress)
    }

    # ----------------------------------------------------
    # KEGG GSEA (No simplification - KEGG specific)
    # ----------------------------------------------------
    ids <- tryCatch({
      bitr(names(original_gene_list), fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org_db_obj)
    }, error = function(e) NULL)

    if(!is.null(ids)) {
      dedup_ids <- ids[!duplicated(ids$UNIPROT), ]
      qc$n_id_mapped_uniprot_to_entrez <- nrow(dedup_ids)
      df2 <- merge(df, dedup_ids, by.x = "gene_symbol", by.y = "UNIPROT")

      kegg_gene_list <- df2$log2fc
      names(kegg_gene_list) <- df2$ENTREZID
      kegg_gene_list <- kegg_gene_list[!duplicated(names(kegg_gene_list))]
      kegg_gene_list <- sort(na.omit(kegg_gene_list), decreasing = TRUE)

      kk2 <- gseKEGG(geneList = kegg_gene_list, organism = "mmu", 
                     minGSSize = analysis_params$min_gs_size,
                     maxGSSize = analysis_params$max_gs_size,
                     pvalueCutoff = analysis_params$pvalue_cutoff,
                     pAdjustMethod = analysis_params$p_adjust_method,
                     keyType = "ncbi-geneid", verbose = FALSE)

      if (!is.null(kk2) && nrow(kk2@result) > 0) {
        kegg_dot <- clusterProfiler::dotplot(kk2, showCategory = 10, split = ".sign") + 
          facet_wrap(~ .sign, nrow = 1) + 
          labs(title = "KEGG GSEA") + 
          theme_minimal()
        save_plot_organized(kegg_dot, "KEGG_dotplot.svg", dirs$plots_kegg)

        tryCatch({
          save_plot_organized(emapplot(pairwise_termsim(kk2), showCategory = 10), "KEGG_emap.svg", dirs$plots_kegg)
        }, error=function(e){})

        tryCatch({
          save_plot_organized(cnetplot(kk2, categorySize = "pvalue", foldChange = gene_list), "KEGG_cnet.svg", dirs$plots_kegg)
        }, error=function(e){})

        tryCatch({
          save_plot_organized(ridgeplot(kk2), "KEGG_ridge.svg", dirs$plots_kegg)
        }, error=function(e){})

        save_plot_organized(gseaplot(kk2, by = "all", title = kk2$Description[1], geneSetID = 1), "KEGG_plot.svg", dirs$plots_kegg)
        write.csv(kk2@result, file = file.path(dirs$kegg, "KEGG_GSEA_results.csv"), row.names = FALSE)
        # Also save the version used for plots to core_enrich
        write.csv(kk2@result, file = file.path(dirs$core_enrich, paste0(comparison_name, "_KEGG.csv")), row.names = FALSE)
        write.csv(kk2@result, file = file.path(dirs$core_enrich_routed, paste0(comparison_name, "_KEGG.csv")), row.names = FALSE)

        # Additionally, save to core_enrich with ontology set to 'KEGG'
        core_enrich_kegg_dir <- file.path(CANONICAL_PATHS$source_data, "KEGG")
        if (!dir.exists(core_enrich_kegg_dir)) dir.create(core_enrich_kegg_dir, recursive = TRUE)
        write.csv(kk2@result, file = file.path(core_enrich_kegg_dir, paste0(comparison_name, "_KEGG.csv")), row.names = FALSE)

        core_enrich_kegg_routed_dir <- file.path(CANONICAL_PATHS$source_data, "KEGG", dirs$route_category, dirs$route_unit)
        if (!dir.exists(core_enrich_kegg_routed_dir)) dir.create(core_enrich_kegg_routed_dir, recursive = TRUE)
        write.csv(kk2@result, file = file.path(core_enrich_kegg_routed_dir, paste0(comparison_name, "_KEGG.csv")), row.names = FALSE)
        kegg_core_table <- file.path(core_enrich_kegg_routed_dir, paste0(comparison_name, "_KEGG.csv"))
        kegg_plot <- file.path(dirs$plots_kegg, "KEGG_dotplot.svg")
        manifest_rows[[length(manifest_rows) + 1]] <- make_manifest_row(
          result_type = "GSEA_KEGG",
          ontology = "KEGG",
          comparison_name = comparison_name,
          dirs = dirs,
          input_gene_file = data_path,
          config_path = config_path,
          output_table = kegg_core_table,
          output_plot = if (file.exists(kegg_plot)) kegg_plot else NA_character_,
          n_genes = length(kegg_gene_list),
          n_terms = nrow(kk2@result),
          simplified = FALSE,
          used_for_plot = TRUE,
          plot_suffix = "_KEGG"
        )
        qc$n_kegg_terms <- nrow(kk2@result)
      } else {
        # Create empty file so we know it ran but found nothing
        write.csv(data.frame(), file = file.path(dirs$kegg, "KEGG_GSEA_results_EMPTY.csv"))
        qc$n_kegg_terms <- 0
      }
      print_progress_step(comparison_name, "KEGG", paste0("terms=", qc$n_kegg_terms), runtime_params$show_step_progress)

      # Pathview
      pathview_dir <- normalizePath(dirs$pathview, winslash = "/", mustWork = FALSE)
      if (length(path_ids) > 0) {
        tryCatch({
          if (!requireNamespace("withr", quietly = TRUE)) stop("Package 'withr' is required for pathview output isolation.")
          withr::with_dir(pathview_dir, {
            lapply(path_ids, function(pid) {
              try({
                pathview(gene.data = kegg_gene_list, pathway.id = pid, species = "mmu",
                         low = "#6698CC", mid = "white", high = "#F08C21", file.type = "svg")
              }, silent = TRUE)
            })
          })
        }, error = function(e) {
          warning("Pathview failed: ", e$message)
        })
      } else {
        message("No path_ids provided; skipping pathview for ", comparison_name)
      }
    }
    
    # ----------------------------------------------------
    # EnrichGO Analysis (ALL ontologies) WITH OPTIONAL SIMPLIFICATION
    # ----------------------------------------------------
    go_universe <- if (length(background_universe) > 0) background_universe else names(gene_list)
    go_enrich <- enrichGO(gene = names(gene_list), universe = go_universe, 
                          OrgDb = org_db_obj, keyType = 'UNIPROT', readable = TRUE, 
                ont = ont,
                pvalueCutoff = analysis_params$pvalue_cutoff,
                qvalueCutoff = analysis_params$qvalue_cutoff,
                pAdjustMethod = analysis_params$p_adjust_method,
                minGSSize = analysis_params$min_gs_size,
                maxGSSize = analysis_params$max_gs_size)
    
    
    # Initialize simplified version as NULL
    go_enrich_simplified <- NULL
    
    # Process simplification ONLY if requested
    if (PERFORM_SIMPLIFICATION) {
      simplified_results_list <- list()
      
      for (ont_type in c("BP", "CC", "MF")) {
        # Run enrichGO for this specific ontology
        go_single <- tryCatch({
          enrichGO(gene = names(gene_list), 
                   universe = go_universe, 
                   OrgDb = org_db_obj, 
                   keyType = 'UNIPROT', 
                   readable = TRUE, 
                   ont = ont_type,
                   pvalueCutoff = analysis_params$pvalue_cutoff,
                   qvalueCutoff = analysis_params$qvalue_cutoff,
                   pAdjustMethod = analysis_params$p_adjust_method,
                   minGSSize = analysis_params$min_gs_size,
                   maxGSSize = analysis_params$max_gs_size)
        }, error = function(e) {
          warning("enrichGO failed for ", ont_type, ": ", e$message)
          NULL
        })
        
        # Simplify if results exist
        if (!is.null(go_single) && nrow(go_single@result) > 0) {
          go_simplified <- tryCatch({
            simplify(go_single, cutoff = SIMPLIFY_CUTOFF, by = "p.adjust", select_fun = min)
          }, error = function(e) {
            warning("Simplify failed for ", ont_type, ": ", e$message)
            go_single
          })
          
          # Store results as data.frame
          if (is(go_simplified, "enrichResult")) {
            simplified_results_list[[ont_type]] <- go_simplified@result
          } else if (is.data.frame(go_simplified)) {
            simplified_results_list[[ont_type]] <- go_simplified
          } else {
            simplified_results_list[[ont_type]] <- go_single@result
          }
        }
      }
      
      # Combine simplified results
      if (length(simplified_results_list) > 0) {
        go_enrich_simplified <- go_enrich
        go_enrich_simplified@result <- do.call(rbind, simplified_results_list)
      }
    }
    
    
    # Choose version for plotting
    if (PERFORM_SIMPLIFICATION && !is.null(go_enrich_simplified)) {
      go_enrich_for_plot <- if(USE_SIMPLIFIED_FOR_PLOTS) go_enrich_simplified else go_enrich
    } else {
      go_enrich_for_plot <- go_enrich
    }
    
    
    # Plot with selected version
    if (nrow(go_enrich_for_plot@result) > 0) {
      plot_label <- if(PERFORM_SIMPLIFICATION && USE_SIMPLIFIED_FOR_PLOTS && !is.null(go_enrich_simplified)) {
        "(Simplified)"
      } else {
        "(Full)"
      }
      
      p4 <- clusterProfiler::dotplot(go_enrich_for_plot, showCategory = 20, split = "ONTOLOGY") +
        facet_grid(ONTOLOGY ~ ., scales = "free_y") +
        labs(title = paste("GO Enrichment Dotplot", plot_label), 
             x = "Gene Ratio", y = "GO Term", color = "p.adjust", size = "Count") +
        scale_color_viridis_c(option = "magma", direction = -1) +
        theme_minimal(base_size = 12) +
        theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
              axis.text.x = element_text(angle = 45, hjust = 1))
      
      save_plot_organized(p4, paste0("GOenrich_dotplot", plot_suffix, ".svg"), dirs$plots_go)
    }
    
    
    # Save results
    write.csv(go_enrich@result, 
              file = file.path(dirs$go_ont, paste0("enrichGO_ALL_results_full.csv")), 
              row.names = FALSE)
    print_progress_step(comparison_name, "enrichGO", paste0("terms=", nrow(go_enrich@result)), runtime_params$show_step_progress)
    
    if (PERFORM_SIMPLIFICATION && !is.null(go_enrich_simplified)) {
      write.csv(go_enrich_simplified@result, 
                file = file.path(dirs$go_ont, paste0("enrichGO_ALL_results_simplified.csv")), 
                row.names = FALSE)
    }
    
    # ----------------------------------------------------
    # KEGG GSEA with Predefined UniProt IDs
    # ----------------------------------------------------
     selected_entrez <- NULL
     if (length(selected_uniprot) > 0) {
      selected_entrez <- tryCatch({
        bitr(selected_uniprot, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org_db_obj)
      }, error=function(e) NULL)
     }
    
    if(!is.null(selected_entrez)) {
      selected_df <- merge(df, selected_entrez, by.x = "gene_symbol", by.y = "UNIPROT")
      if(nrow(selected_df) > 0) {
        selected_kegg_list <- selected_df$log2fc
        names(selected_kegg_list) <- selected_df$ENTREZID
        selected_kegg_list <- selected_kegg_list[!duplicated(names(selected_kegg_list))]
        selected_kegg_list <- sort(na.omit(selected_kegg_list), decreasing = TRUE)
        
        gsea_kegg_selected <- gseKEGG(geneList = selected_kegg_list, organism = "mmu", 
                                      minGSSize = analysis_params$min_gs_size,
                                      maxGSSize = analysis_params$max_gs_size,
                                      pvalueCutoff = analysis_params$pvalue_cutoff,
                                      pAdjustMethod = analysis_params$p_adjust_method,
                                      keyType = "ncbi-geneid", verbose = FALSE)
        
        if (!is.null(gsea_kegg_selected) && nrow(gsea_kegg_selected@result) > 0) {
          kegg_selected_dot <- clusterProfiler::dotplot(gsea_kegg_selected, showCategory = 10, split = ".sign") + 
            facet_wrap(~ .sign, nrow = 1) + 
            labs(title = "KEGG GSEA (Predefined)") + 
            theme_minimal()
          save_plot_organized(kegg_selected_dot, "KEGG_Predefined_dotplot.svg", dirs$plots_kegg)
          write.csv(gsea_kegg_selected@result, file = file.path(dirs$kegg, "KEGG_GSEA_Predefined_UniProt.csv"), row.names = FALSE)
          
          tryCatch({
            save_plot_organized(emapplot(pairwise_termsim(gsea_kegg_selected), showCategory = 10), "KEGG_Predefined_emap.svg", dirs$plots_kegg)
          }, error=function(e){})
          
          tryCatch({
            save_plot_organized(cnetplot(gsea_kegg_selected, categorySize = "pvalue", foldChange = selected_kegg_list), "KEGG_Predefined_cnet.svg", dirs$plots_kegg)
          }, error=function(e){})
          
          tryCatch({
            save_plot_organized(ridgeplot(gsea_kegg_selected), "KEGG_Predefined_ridge.svg", dirs$plots_kegg)
          }, error=function(e){})
          
          save_plot_organized(gseaplot(gsea_kegg_selected, by = "all", title = gsea_kegg_selected$Description[1], geneSetID = 1), "KEGG_Predefined_plot.svg", dirs$plots_kegg)
        }
      }
    }
    
    # ----------------------------------------------------
    # Custom GSEA: NK3R-signalling (No simplification needed)
    # ----------------------------------------------------
    if (length(nk3r_genes) > 0) {
      term2gene_nk3r <- data.frame(term = rep("NK3R-signalling", length(nk3r_genes)), gene = nk3r_genes)
      custom_gene_list <- df$log2fc
      names(custom_gene_list) <- df$gene_symbol
      custom_gene_list <- sort(na.omit(custom_gene_list), decreasing = TRUE)
      custom_gene_list <- custom_gene_list[!duplicated(names(custom_gene_list))]
      
      gsea_nk3r <- clusterProfiler::GSEA(geneList = custom_gene_list, TERM2GENE = term2gene_nk3r, 
                                         pvalueCutoff = analysis_params$pvalue_cutoff,
                                         minGSSize = 1,
                                         maxGSSize = 500,
                                         verbose = FALSE)
      
      if (!is.null(gsea_nk3r) && nrow(gsea_nk3r@result) > 0) {
        nk3r_dot <- clusterProfiler::dotplot(gsea_nk3r, showCategory = 10, split = ".sign") + 
          facet_wrap(~ .sign, nrow = 1) + 
          labs(title = "NK3R-signalling GSEA") + 
          theme_minimal()
        save_plot_organized(nk3r_dot, "NK3R_dotplot.svg", dirs$plots_custom)
        
        nk3r_plot <- gseaplot(gsea_nk3r, by = "all", title = "NK3R-signalling", geneSetID = 1)
        save_plot_organized(nk3r_plot, "NK3R_gsea_plot.svg", dirs$plots_custom)
        openxlsx::write.xlsx(gsea_nk3r@result, file = file.path(dirs$custom, "NK3R_GSEA_results.xlsx"))
      }
    } else {
      message("No nk3r_genes provided; skipping custom NK3R GSEA for ", comparison_name)
    }
    
    qc$status <- "SUCCESS"
    qc$runtime_seconds <- as.numeric(difftime(Sys.time(), run_start, units = "secs"))
    write.csv(qc, file.path(dirs$results, "QC_summary.csv"), row.names = FALSE)
    write_completed_checkpoint(dirs)
    write_log_line(comparison_log, "INFO", comparison_name, "DONE", paste0("Completed in ", round(qc$runtime_seconds, 2), " sec"))
    print_progress_step(comparison_name, "DONE", paste0("runtime_sec=", round(qc$runtime_seconds, 2)), runtime_params$show_step_progress)
    return(list(status = "SUCCESS", comparison = comparison_name, error = NA_character_, qc = qc, manifest = dplyr::bind_rows(manifest_rows)))
    
  }, error = function(e) {
    qc$status <- "ERROR"
    qc$runtime_seconds <- as.numeric(difftime(Sys.time(), run_start, units = "secs"))
    write_log_line(comparison_log, "ERROR", comparison_name, "UNHANDLED", conditionMessage(e))
    return(list(status = "ERROR", comparison = comparison_name, error = conditionMessage(e), qc = qc, manifest = dplyr::bind_rows(manifest_rows)))
  })
}

# ----------------------------------------------------
# 6. RUN PARALLEL ANALYSIS WITH FUTURE
# ----------------------------------------------------
cat("\n==============================================\n")
cat("STARTING PARALLEL GSEA ANALYSIS\n")
cat("==============================================\n\n")
cat("Launching", length(comparison_list), "comparisons across", n_cores, "workers...\n")
if (length(comparison_list) < n_cores) {
  cat("[WARNING] Number of comparisons (", length(comparison_list), ") < workers (", n_cores, ")\n")
  cat("[WARNING] Not all workers will be used; consider increasing comparison_list or reducing workers.\n")
}
cat("Progress updates now reflect both STARTED and FINISHED states.\n\n")

# Verify the plan before execution
cat("[DEBUG PARALLEL] Before future_lapply:\n")
cat("[DEBUG PARALLEL]   Current plan:", toString(class(plan())), "\n")
cat("[DEBUG PARALLEL]   Number of workers:", nbrOfWorkers(), "\n")
cat("[DEBUG PARALLEL]   Number of comparisons to process:", length(comparison_list), "\n")
cat("[DEBUG PARALLEL]   Future.seed setting: TRUE\n")

if (isTRUE(runtime_params$show_progress)) {
  progressr::handlers("txtprogressbar")
  results <- progressr::with_progress({
    p <- progressr::progressor(steps = max(1L, 2L * length(comparison_list)))
    future_lapply(comparison_list, function(cell_types) {
      comparison_name <- paste(cell_types, collapse = "_")
      p(message = paste0(comparison_name, " -> STARTED"))
      res <- analyze_comparison(cell_types, working_base, mapped_data_base, organism, ont,
                                nk3r_genes, selected_uniprot, path_ids,
                                analysis_params = analysis_params,
                                runtime_params = runtime_params,
                                run_log_dir = run_log_dir,
                                background_universe = background_universe,
                                config_path = config_path)
      p(message = paste0(res$comparison, " -> ", res$status))
      res
    }, future.seed = TRUE)
  })
} else {
  results <- future_lapply(comparison_list, function(cell_types) {
    analyze_comparison(cell_types, working_base, mapped_data_base, organism, ont,
                       nk3r_genes, selected_uniprot, path_ids,
                       analysis_params = analysis_params,
                       runtime_params = runtime_params,
                       run_log_dir = run_log_dir,
                       background_universe = background_universe,
                       config_path = config_path)
  }, future.seed = TRUE)
}

cat("[DEBUG PARALLEL] After future_lapply completed\n")
cat("[DEBUG PARALLEL]   Number of results:", length(results), "\n")

# Reset to sequential processing
plan(sequential)

results <- c(completed_results, results)

} else {
  cat("\n==============================================\n")
  cat("RESUME CHECK FOUND NO PENDING COMPARISONS\n")
  cat("==============================================\n\n")
  results <- completed_results
}

# Print summary
cat("\n==============================================\n")
cat("PARALLEL ANALYSIS SUMMARY\n")
cat("==============================================\n\n")

for (result in results) {
  if (result$status == "SUCCESS") {
    cat("✓", result$comparison, "- COMPLETED\n")
  } else if (result$status == "SKIPPED") {
    cat("○", result$comparison, "- SKIPPED (checkpoint exists)\n")
  } else {
    cat("✗", result$comparison, "- FAILED:", result$error, "\n")
  }
}

cat("\n==============================================\n")
cat("ALL COMPARISONS COMPLETED!\n")
cat("==============================================\n\n")

# ----------------------------------------------------
# 7. RUN-LEVEL SUMMARY OUTPUTS
# ----------------------------------------------------
summary_dir <- CANONICAL_PATHS$reports
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

status_vec <- vapply(results, function(x) x$status, character(1))
error_vec <- vapply(results, function(x) ifelse(is.null(x$error), NA_character_, x$error), character(1))
comparison_vec <- vapply(results, function(x) x$comparison, character(1))

run_summary <- data.frame(
  comparison = comparison_vec,
  status = status_vec,
  error = error_vec,
  stringsAsFactors = FALSE
)

run_summary_file <- file.path(summary_dir, paste0("clusterProfiler_run_summary_", run_id, ".csv"))
write.csv(run_summary, run_summary_file, row.names = FALSE)

qc_rows <- lapply(results, function(x) x$qc)
qc_rows <- qc_rows[!vapply(qc_rows, is.null, logical(1))]
if (length(qc_rows) > 0) {
  run_qc <- dplyr::bind_rows(qc_rows)
  run_qc_file <- file.path(summary_dir, paste0("clusterProfiler_qc_summary_", run_id, ".csv"))
  write.csv(run_qc, run_qc_file, row.names = FALSE)
}

manifest_rows <- lapply(results, function(x) x$manifest)
manifest_rows <- manifest_rows[!vapply(manifest_rows, is.null, logical(1))]
manifest <- if (length(manifest_rows) > 0) dplyr::bind_rows(manifest_rows) else {
  data.frame(matrix(ncol = length(manifest_columns), nrow = 0, dimnames = list(NULL, manifest_columns)))
}
manifest_file <- file.path(CANONICAL_PATHS$processed, "clusterProfiler_manifest.csv")
write.csv(manifest, manifest_file, row.names = FALSE)
write.csv(manifest, file.path(CANONICAL_PATHS$reports, paste0("clusterProfiler_manifest_", run_id, ".csv")), row.names = FALSE)

summary_txt <- file.path(summary_dir, paste0("clusterProfiler_run_summary_", run_id, ".txt"))
summary_lines <- c(
  "clusterProfiler parallel run summary",
  paste0("Run ID: ", run_id),
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Config: ", config_path),
  paste0("Total comparisons: ", length(results)),
  paste0("SUCCESS: ", sum(status_vec == "SUCCESS")),
  paste0("SKIPPED: ", sum(status_vec == "SKIPPED")),
  paste0("FAILED/ERROR: ", sum(status_vec %in% c("FAILED", "ERROR"))),
  paste0("Master log: ", master_log),
  paste0("Manifest: ", manifest_file),
  paste0("Run summary CSV: ", run_summary_file)
)
writeLines(summary_lines, con = summary_txt)
write_log_line(master_log, "INFO", "GLOBAL", "SUMMARY", paste0("Run summary written: ", run_summary_file))

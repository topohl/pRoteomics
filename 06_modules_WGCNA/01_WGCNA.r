# ================================ parallel-enabled
# WGCNA with profile-aware spatial traits, preservation, and condition panels
# Outputs organized into canonical module folders by default
# ================================ parallel-enabled

early_args <- commandArgs(trailingOnly = TRUE)
early_has_flag <- function(flag) flag %in% early_args
early_arg_value <- function(flag, default = "") {
  hit <- which(early_args == flag)
  if (!length(hit) || hit[1] == length(early_args)) return(default)
  early_args[[hit[1] + 1]]
}
if (early_has_flag("--dry-run") || tolower(Sys.getenv("PROTEOMICS_DRY_RUN", unset = "")) %in% c("1", "true", "yes")) {
  paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
  source(paths_file)
  source(repo_path("R", "dataset_config.R"))
  source(repo_path("R", "dataset_inputs.R"))
  dataset_cli_early <- early_arg_value("--dataset", default = "")
  if (nzchar(dataset_cli_early)) Sys.setenv(PROTEOMICS_DATASET = validate_dataset(dataset_cli_early, source = "--dataset"))
  dataset_profile_early <- {
    profile_override <- Sys.getenv("PROTEOMICS_WGCNA_DATASET_PROFILE", unset = "")
    if (nzchar(dataset_cli_early)) validate_dataset(dataset_cli_early, source = "--dataset")
    else if (nzchar(profile_override)) validate_dataset(profile_override, source = "PROTEOMICS_WGCNA_DATASET_PROFILE")
    else current_dataset()
  }
  canonical_paths_early <- module_paths("06_modules_WGCNA", file.path("01_WGCNA", dataset_profile_early))
  output_dir_env_early <- Sys.getenv("PROTEOMICS_WGCNA_OUTPUT_DIR", unset = "")
  output_root_early <- if (nzchar(output_dir_env_early)) file.path(output_dir_env_early, dataset_profile_early) else canonical_paths_early$reports
  subdirs_early <- if (nzchar(output_dir_env_early)) {
    list(
      figures_qc = file.path(output_root_early, "figures", "qc"),
      figures_network = file.path(output_root_early, "figures", "network"),
      figures_traits = file.path(output_root_early, "figures", "traits"),
      figures_main = file.path(output_root_early, "figures", "main"),
      tables_qc = file.path(output_root_early, "tables", "qc"),
      tables_mapping = file.path(output_root_early, "tables", "mapping"),
      tables_modules = file.path(output_root_early, "tables", "modules"),
      tables_pres = file.path(output_root_early, "tables", "preservation"),
      tables_traits = file.path(output_root_early, "tables", "traits"),
      tables_supermodules = file.path(output_root_early, "tables", "supermodules"),
      figures_supermodules = file.path(output_root_early, "figures", "supermodules"),
      source_data = file.path(output_root_early, "source_data"),
      state = file.path(output_root_early, "state"),
      logs = file.path(output_root_early, "logs")
    )
  } else {
    list(
      figures_qc = file.path(canonical_paths_early$figures, "qc"),
      figures_network = file.path(canonical_paths_early$figures, "network"),
      figures_traits = file.path(canonical_paths_early$figures, "traits"),
      figures_main = file.path(canonical_paths_early$figures, "main"),
      tables_qc = file.path(canonical_paths_early$tables, "qc"),
      tables_mapping = file.path(canonical_paths_early$tables, "mapping"),
      tables_modules = file.path(canonical_paths_early$tables, "modules"),
      tables_pres = file.path(canonical_paths_early$tables, "preservation"),
      tables_traits = file.path(canonical_paths_early$tables, "traits"),
      tables_supermodules = file.path(canonical_paths_early$tables, "supermodules"),
      figures_supermodules = file.path(canonical_paths_early$figures, "supermodules"),
      source_data = canonical_paths_early$source_data,
      state = canonical_paths_early$processed,
      logs = canonical_paths_early$logs
    )
  }
  invisible(lapply(c(output_root_early, unlist(subdirs_early)), dir_create))
  dataset_inputs_early <- resolve_dataset_inputs(dataset_profile_early, purpose = "wgcna")
  expr_xlsx_env_early <- Sys.getenv("PROTEOMICS_WGCNA_EXPR_XLSX", unset = "")
  meta_xlsx_env_early <- Sys.getenv("PROTEOMICS_WGCNA_META_XLSX", unset = "")
  expr_xlsx_early <- if (nzchar(expr_xlsx_env_early)) expr_xlsx_env_early else path_processed("06_modules_WGCNA", "01_WGCNA", dataset_profile_early, "inputs", "wgcna_expression.xlsx")
  meta_xlsx_early <- if (nzchar(meta_xlsx_env_early)) meta_xlsx_env_early else path_processed("06_modules_WGCNA", "01_WGCNA", dataset_profile_early, "inputs", "wgcna_sample_info.xlsx")
  idmap_dat_early <- Sys.getenv("PROTEOMICS_WGCNA_IDMAP_DAT", unset = dataset_inputs_early$idmap_file)
  wgcna_final_state_path_early <- file.path(subdirs_early$state, "wgcna_final_model_state.rds")
  reuse_completed_analysis_early <- tolower(Sys.getenv("PROTEOMICS_WGCNA_REUSE_STATE", unset = "true")) %in% c("1", "true", "yes", "y")
  force_full_analysis_early <- tolower(Sys.getenv("PROTEOMICS_WGCNA_FORCE_FULL", unset = "false")) %in% c("1", "true", "yes", "y")
  can_stage_inputs_early <- file.exists(dataset_inputs_early$expression_file) && file.exists(dataset_inputs_early$metadata_file)
  can_use_inputs_early <- file.exists(expr_xlsx_early) && file.exists(meta_xlsx_early)
  can_use_cache_early <- file.exists(wgcna_final_state_path_early) && !force_full_analysis_early
  can_write_outputs_early <- all(vapply(c(output_root_early, unlist(subdirs_early)), function(path) dir.exists(path) && file.access(path, 2) == 0, logical(1)))
  sample_check_early <- "not checked"
  sample_check_ok_early <- TRUE
  if (can_use_inputs_early) {
    sample_check_ok_early <- tryCatch({
      if (!requireNamespace("readxl", quietly = TRUE)) stop("readxl unavailable for cheap workbook sample check")
      first_existing_col_early <- function(df, candidates) {
        nms_clean <- tolower(gsub("[^a-z0-9]", "", names(df)))
        cand_clean <- tolower(gsub("[^a-z0-9]", "", candidates))
        hit <- match(cand_clean, nms_clean)
        hit <- hit[!is.na(hit)]
        if (!length(hit)) return(NA_character_)
        names(df)[hit[1]]
      }
      expr_head <- readxl::read_excel(expr_xlsx_early, n_max = 3)
      meta_head <- readxl::read_excel(meta_xlsx_early, n_max = 1000)
      sample_col <- first_existing_col_early(meta_head, dataset_inputs_early$sample_id_col_candidates)
      if (is.na(sample_col)) {
        sample_check_early <<- "metadata sample column not detected"
        FALSE
      } else {
        expr_samples <- setdiff(names(expr_head), dataset_inputs_early$protein_id_col_candidates)
        n_match <- length(intersect(expr_samples, as.character(meta_head[[sample_col]])))
        sample_check_early <<- paste0(n_match, " matching sample IDs")
        n_match > 0
      }
    }, error = function(e) {
      sample_check_early <<- conditionMessage(e)
      FALSE
    })
  } else if (can_stage_inputs_early) {
    sample_check_early <- "staged files missing; canonical upstream exists and can be staged"
  } else {
    sample_check_early <- "staged files missing and canonical upstream cannot be staged"
    sample_check_ok_early <- FALSE
  }
  downstream_contract_early <- file.path(subdirs_early$tables_modules, "WGCNA_module_definitions_for_downstream.csv")
  dry_run_line("Script", "06_modules_WGCNA/01_WGCNA.r")
  dry_run_line("Dataset", dataset_profile_early)
  dry_run_line("Resolved input diagnostics", paste(dataset_inputs_early$diagnostics, collapse = " | "))
  dry_run_line("Reuse cached final state", reuse_completed_analysis_early)
  dry_run_line("Force full analysis", force_full_analysis_early)
  dry_run_line("Cached state", wgcna_final_state_path_early, if (can_use_cache_early) "PASS" else if (reuse_completed_analysis_early) "WARN" else "INFO")
  dry_run_line("Explicit/staged expression workbook", expr_xlsx_early, if (file.exists(expr_xlsx_early)) "PASS" else "WARN")
  dry_run_line("Explicit/staged sample metadata workbook", meta_xlsx_early, if (file.exists(meta_xlsx_early)) "PASS" else "WARN")
  dry_run_line("Canonical upstream expression", dataset_inputs_early$expression_file, if (file.exists(dataset_inputs_early$expression_file)) "PASS" else "FAIL")
  dry_run_line("Canonical upstream metadata", dataset_inputs_early$metadata_file, if (file.exists(dataset_inputs_early$metadata_file)) "PASS" else "FAIL")
  dry_run_line("Staged input generation", if (can_use_inputs_early) "already staged" else if (can_stage_inputs_early) "can be generated" else "cannot be generated", if (can_use_inputs_early || can_stage_inputs_early) "PASS" else "FAIL")
  dry_run_line("Cheap sample matching", sample_check_early, if (sample_check_ok_early) "PASS" else "FAIL")
  dry_run_line("Mouse idmapping", idmap_dat_early, if (file.exists(idmap_dat_early)) "PASS" else "FAIL")
  dry_run_line("Output folders writable", paste(unlist(subdirs_early), collapse = "; "), if (can_write_outputs_early) "PASS" else "FAIL")
  dry_run_line("Downstream module contract", downstream_contract_early, if (file.exists(downstream_contract_early)) "PASS" else "WARN")
  dry_run_line("Optional DE/GSEA overlap bridge", repo_path("06_modules_WGCNA", "05_wgcna_de_gsea_overlap.r"), "WARN")
  quit(status = if (file.exists(idmap_dat_early) && can_write_outputs_early && sample_check_ok_early && (can_use_inputs_early || can_stage_inputs_early || can_use_cache_early)) 0 else 1, save = "no")
}

# Packages
required_pkgs <- c(
  "WGCNA", "flashClust", "curl", "readxl", "ggplot2", "svglite", "GO.db",
  "reshape2", "gtools", "patchwork", "cowplot", "pheatmap", "dplyr", "tidyr",
  "httr", "jsonlite", "purrr", "AnnotationDbi", "org.Mm.eg.db", "readr",
  "stringr", "tibble", "UniProt.ws", "RColorBrewer", "ggpubr", "broom", "grid",
  "clusterProfiler", "scales", "writexl"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop(
    "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
    ". Install them explicitly before running this manuscript pipeline."
  )
}
suppressWarnings(
  suppressPackageStartupMessages(
    invisible(lapply(required_pkgs, library, character.only = TRUE))
  )
)

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "dataset_inputs.R"))

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = "") {
  hit <- which(args == flag)
  if (!length(hit) || hit[1] == length(args)) return(default)
  args[[hit[1] + 1]]
}
has_flag <- function(flag) flag %in% args
dataset_cli <- arg_value("--dataset", default = "")
if (nzchar(dataset_cli)) Sys.setenv(PROTEOMICS_DATASET = validate_dataset(dataset_cli, source = "--dataset"))
wgcna_dry_run <- has_flag("--dry-run") || is_dry_run()

mm_to_in <- function(mm) mm / 25.4
figure_single_col <- mm_to_in(89)
figure_double_col <- mm_to_in(183)
graphics_font_works <- function(family) {
  tf_svg <- tempfile(fileext = ".svg")
  tf_pdf <- tempfile(fileext = ".pdf")
  dev_id <- grDevices::dev.cur()
  ok <- tryCatch({
    svglite::svglite(file = tf_svg, width = 1, height = 1)
    grid::grid.text("test", gp = grid::gpar(fontfamily = family, fontsize = 8))
    grDevices::dev.off()
    grDevices::pdf(file = tf_pdf, width = 1, height = 1, family = family, useDingbats = FALSE)
    grid::grid.newpage()
    TRUE
  }, error = function(e) FALSE)
  if (grDevices::dev.cur() != dev_id) grDevices::dev.off()
  unlink(c(tf_svg, tf_pdf), force = TRUE)
  ok
}

resolve_figure_font <- function(preferred = Sys.getenv("PROTEOMICS_FIGURE_FONT", unset = "sans")) {
  preferred <- trimws(preferred)
  if (!nzchar(preferred)) preferred <- "sans"

  font_available <- FALSE
  if (requireNamespace("systemfonts", quietly = TRUE)) {
    font_available <- tryCatch({
      fonts <- systemfonts::system_fonts()
      any(tolower(fonts$family) == tolower(preferred))
    }, error = function(e) FALSE)
  }

  if ((isTRUE(font_available) || identical(tolower(preferred), "sans")) && graphics_font_works(preferred)) {
    return(preferred)
  }

  warning(
    "Figure font '", preferred, "' is not available to the R graphics device; ",
    "falling back to 'sans'. Set PROTEOMICS_FIGURE_FONT to an installed family ",
    "to override.",
    call. = FALSE
  )
  "sans"
}
figure_font <- resolve_figure_font()
figure_base_size <- 7
figure_axis_size <- 6.2
figure_title_size <- 7
figure_line <- 0.25
figure_diverging <- c(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B")
figure_condition_cols <- c(con = "#4D4D4D", res = "#0072B2", sus = "#D55E00")
figure_condition_labels <- c(con = "CON", res = "RES", sus = "SUS")

theme_publication <- function(base_size = figure_base_size) {
  ggplot2::theme_classic(base_size = base_size, base_family = figure_font) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = figure_title_size, face = "plain", hjust = 0),
      plot.subtitle = ggplot2::element_text(size = figure_axis_size, color = "grey30", margin = ggplot2::margin(t = 1, b = 2)),
      axis.title = ggplot2::element_text(size = figure_axis_size),
      axis.text = ggplot2::element_text(size = figure_axis_size, color = "black"),
      axis.line = ggplot2::element_line(linewidth = figure_line, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = figure_line, color = "black"),
      legend.title = ggplot2::element_text(size = figure_axis_size),
      legend.text = ggplot2::element_text(size = figure_axis_size),
      legend.key.size = grid::unit(3, "mm"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = figure_axis_size, color = "black"),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(3, 3, 3, 3)
    )
}
ggplot2::theme_set(theme_publication())

# Parallel setup
nCores <- tryCatch({
  pc <- parallel::detectCores(logical = FALSE)
  if (is.na(pc) || pc < 2) 2 else pc
}, error = function(e) 2)
enableWGCNAThreads(nThreads = nCores)
WGCNAnThreads()

# --------------------------
# Paths and data load
# --------------------------
safe_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (file.access(path, 2) != 0) stop(sprintf("Not writable: %s", path))
  invisible(normalizePath(path))
}

wgcna_module <- "06_modules_WGCNA"
dataset_profile <- {
  profile_override <- Sys.getenv("PROTEOMICS_WGCNA_DATASET_PROFILE", unset = "")
  if (nzchar(dataset_cli)) validate_dataset(dataset_cli, source = "--dataset")
  else if (nzchar(profile_override)) validate_dataset(profile_override, source = "PROTEOMICS_WGCNA_DATASET_PROFILE")
  else current_dataset()
}
wgcna_substep <- file.path("01_WGCNA", dataset_profile)
output_dir_env <- Sys.getenv("PROTEOMICS_WGCNA_OUTPUT_DIR", unset = "")
if (nzchar(output_dir_env)) {
  output_dir <- file.path(output_dir_env, dataset_profile)
  subdirs <- list(
    figures_qc          = file.path(output_dir, "figures", "qc"),
    figures_network     = file.path(output_dir, "figures", "network"),
    figures_traits      = file.path(output_dir, "figures", "traits"),
    figures_main        = file.path(output_dir, "figures", "main"),
    tables_qc           = file.path(output_dir, "tables", "qc"),
    tables_mapping      = file.path(output_dir, "tables", "mapping"),
    tables_modules      = file.path(output_dir, "tables", "modules"),
    tables_pres         = file.path(output_dir, "tables", "preservation"),
    tables_traits       = file.path(output_dir, "tables", "traits"),
    tables_supermodules = file.path(output_dir, "tables", "supermodules"),
    figures_supermodules= file.path(output_dir, "figures", "supermodules"),
    source_data         = file.path(output_dir, "source_data"),
    state               = file.path(output_dir, "state"),
    logs                = file.path(output_dir, "logs")
  )
} else {
  canonical_paths <- module_paths(wgcna_module, wgcna_substep)
  output_dir <- canonical_paths$reports
  subdirs <- list(
    figures_qc          = file.path(canonical_paths$figures, "qc"),
    figures_network     = file.path(canonical_paths$figures, "network"),
    figures_traits      = file.path(canonical_paths$figures, "traits"),
    figures_main        = file.path(canonical_paths$figures, "main"),
    tables_qc           = file.path(canonical_paths$tables, "qc"),
    tables_mapping      = file.path(canonical_paths$tables, "mapping"),
    tables_modules      = file.path(canonical_paths$tables, "modules"),
    tables_pres         = file.path(canonical_paths$tables, "preservation"),
    tables_traits       = file.path(canonical_paths$tables, "traits"),
    tables_supermodules = file.path(canonical_paths$tables, "supermodules"),
    figures_supermodules= file.path(canonical_paths$figures, "supermodules"),
    source_data         = canonical_paths$source_data,
    state               = canonical_paths$processed,
    logs                = canonical_paths$logs
  )
}
invisible(lapply(c(output_dir, unlist(subdirs)), safe_dir))

# Path helpers
fp_qc       <- function(...) file.path(subdirs$figures_qc, ...)
fp_traits   <- function(...) file.path(subdirs$figures_traits, ...)
fp_mainfig  <- function(...) file.path(subdirs$figures_main, ...)
fp_net      <- function(...) file.path(subdirs$figures_network, ...)
fp_qctab    <- function(...) file.path(subdirs$tables_qc, ...)
fp_maptab   <- function(...) file.path(subdirs$tables_mapping, ...)
fp_modtab   <- function(...) file.path(subdirs$tables_modules, ...)
fp_traittab <- function(...) file.path(subdirs$tables_traits, ...)
fp_prestab  <- function(...) file.path(subdirs$tables_pres, ...)
fp_supertab <- function(...) file.path(subdirs$tables_supermodules, ...)
fp_superfig <- function(...) file.path(subdirs$figures_supermodules, ...)
fp_source   <- function(...) file.path(subdirs$source_data, ...)
fp_state    <- function(...) file.path(subdirs$state, ...)
fp_log      <- function(...) file.path(subdirs$logs, ...)

write_csv_safe <- function(x, path) {
  readr::write_csv(x, path, na = "")
  invisible(path)
}
write_tsv_safe <- function(x, path) {
  readr::write_tsv(x, path, na = "")
  invisible(path)
}
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x
}
ensure_module_label_schema <- function(df) {
  char_cols <- c(
    "ModuleLabel_Manual", "ModuleLabel_GO_BP", "ModuleLabel_GO_MF", "ModuleLabel_GO_CC",
    "ModuleLabel_Final", "ModuleLabel_Source",
    "best_GO_BP", "best_GO_MF", "best_GO_CC",
    "best_GO_qvalue_BP", "best_GO_qvalue_MF", "best_GO_qvalue_CC",
    "best_GO_gene_ratio_BP", "best_GO_gene_ratio_MF", "best_GO_gene_ratio_CC"
  )
  numeric_cols <- c("best_GO_padj_BP", "best_GO_padj_MF", "best_GO_padj_CC")
  for (col in char_cols) {
    if (!col %in% names(df)) df[[col]] <- NA_character_
  }
  for (col in numeric_cols) {
    if (!col %in% names(df)) df[[col]] <- NA_real_
  }
  df
}
safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  max(x)
}
log_session <- function() {
  writeLines(capture.output(utils::sessionInfo()), fp_log("session_info.txt"))
}
log_session()

# Analysis parameters reported with outputs
sample_tree_cut_height <- 80
sample_tree_plot_height <- 40
soft_threshold_rsquared <- 0.80
min_module_size <- 30
deep_split <- 2
merge_cut_height <- 0.25
module_preservation_permutations <- 100
reuse_completed_analysis <- tolower(Sys.getenv("PROTEOMICS_WGCNA_REUSE_STATE", unset = "true")) %in% c("1", "true", "yes", "y")
force_full_analysis <- tolower(Sys.getenv("PROTEOMICS_WGCNA_FORCE_FULL", unset = "false")) %in% c("1", "true", "yes", "y")
wgcna_final_state_path <- fp_state("wgcna_final_model_state.rds")

# Optional: safe svg helper
save_svg <- function(path, width, height, expr) {
  svglite::svglite(file = path, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  force(expr)
}

save_plot_publication <- function(plot, path_svg, width, height) {
  ggplot2::ggsave(path_svg, plot, device = svglite::svglite,
                  width = width, height = height, units = "in", limitsize = FALSE)
  ggplot2::ggsave(sub("\\.svg$", ".pdf", path_svg), plot, device = grDevices::pdf,
                  width = width, height = height, units = "in", limitsize = FALSE,
                  family = figure_font, useDingbats = FALSE)
  invisible(path_svg)
}

sig_dot <- function(fdr) {
  dplyr::case_when(
    is.na(fdr) ~ "",
    fdr < 0.01 ~ "\u2022",
    fdr < 0.05 ~ "\u00b7",
    TRUE ~ ""
  )
}

supermodule_levels <- c(
  "RNA processing / translation regulation",
  "Energy metabolism / mitochondrial support",
  "Synaptic and vesicle organization",
  "Proteostasis / post-translational regulation",
  "Developmental patterning / tissue-identity annotations",
  "Motility / structural remodeling",
  "Unassigned"
)

manual_supermodule_seed <- tibble::tribble(
  ~module_eigengene, ~ModuleColor, ~Supermodule, ~SupermoduleConfidence, ~SupermoduleRationale,
  "MEforestgreen", "forestgreen", "RNA processing / translation regulation", "high", "Top GO label is RNA processing; grouped with ncRNA and translation modules as an RNA expression-control program.",
  "MEkhaki", "khaki", "RNA processing / translation regulation", "moderate-high", "Top GO label is regulatory ncRNA-mediated gene silencing; mechanistically close to RNA regulation.",
  "MEcoral", "coral", "RNA processing / translation regulation", "high", "Top GO label is cytoplasmic translational initiation; complements RNA processing modules.",
  "MEmediumaquamarine", "mediumaquamarine", "Proteostasis / post-translational regulation", "high", "Top GO label is ubiquitin-dependent protein catabolic process, consistent with proteostasis.",
  "MElightpink", "lightpink", "Proteostasis / post-translational regulation", "high", "Top GO label is protein folding, consistent with chaperone/proteostasis biology.",
  "MEplum", "plum", "Proteostasis / post-translational regulation", "moderate-high", "Top GO label is protein dephosphorylation; interpreted as post-translational regulation.",
  "MElemonchiffon", "lemonchiffon", "Energy metabolism / mitochondrial support", "high", "Top GO label is generation of precursor metabolites and energy.",
  "MEaquamarine", "aquamarine", "Energy metabolism / mitochondrial support", "moderate-high", "Top GO label is acetyl-CoA biosynthetic process, linked to metabolic support.",
  "MEdeepskyblue", "deepskyblue", "Synaptic and vesicle organization", "high", "Top GO label is synapse organization.",
  "MElightblue", "lightblue", "Synaptic and vesicle organization", "high", "Top GO label is vesicle organization, linked to synaptic vesicle/trafficking biology.",
  "MEgold", "gold", "Developmental patterning / tissue-identity annotations", "moderate", "Top GO label is dorsal-ventral pattern formation; interpreted cautiously as patterning/tissue identity.",
  "MElavender", "lavender", "Developmental patterning / tissue-identity annotations", "low-moderate", "Top GO label is skin development; retained as a tissue-identity style annotation but flagged as lower confidence.",
  "MEsaddlebrown", "saddlebrown", "Motility / structural remodeling", "moderate", "Top GO label is negative regulation of cell motility, consistent with remodeling/motility annotation."
)

make_supermodule_annotation <- function(module_label_table = NULL, module_names = character()) {
  module_names <- unique(as.character(module_names))
  module_names <- module_names[nzchar(module_names)]
  present <- tibble::tibble(
    module_eigengene = module_names,
    ModuleColor = sub("^ME", "", module_names),
    present_in_dataset = TRUE
  )

  manual <- manual_supermodule_seed %>%
    dplyr::mutate(
      Supermodule = factor(.data$Supermodule, levels = supermodule_levels),
      manual_annotation = TRUE
    )

  out <- dplyr::full_join(present, manual, by = c("module_eigengene", "ModuleColor")) %>%
    dplyr::mutate(
      present_in_dataset = dplyr::coalesce(.data$present_in_dataset, FALSE),
      manual_annotation = dplyr::coalesce(.data$manual_annotation, FALSE),
      Supermodule = dplyr::coalesce(as.character(.data$Supermodule), "Unassigned"),
      SupermoduleConfidence = dplyr::coalesce(.data$SupermoduleConfidence, "unassigned"),
      SupermoduleRationale = dplyr::coalesce(.data$SupermoduleRationale, "No manual supermodule annotation for this dataset/module.")
    )

  if (!is.null(module_label_table) && nrow(module_label_table)) {
    label_cols <- intersect(
      c("ModuleColor", "ModuleID", "ModuleLabel_Final", "ModuleLabel_GO_BP", "best_GO_BP", "best_GO_padj_BP"),
      names(module_label_table)
    )
    out <- out %>%
      dplyr::left_join(module_label_table[, label_cols, drop = FALSE], by = "ModuleColor")
  }
  for (missing_col in c("ModuleLabel_Final", "ModuleLabel_GO_BP", "best_GO_BP", "best_GO_padj_BP")) {
    if (!missing_col %in% names(out)) out[[missing_col]] <- NA
  }

  out %>%
    dplyr::mutate(
      top_GO_label = dplyr::coalesce(
        .data$ModuleLabel_Final,
        .data$ModuleLabel_GO_BP,
        .data$best_GO_BP,
        paste("Module", .data$ModuleColor)
      ),
      Supermodule = factor(.data$Supermodule, levels = supermodule_levels)
    ) %>%
    dplyr::arrange(.data$Supermodule, .data$ModuleColor)
}

add_supermodule_cols <- function(df, annotation, module_col = "module", color_col = NULL) {
  if (!nrow(df)) return(df)
  ann <- annotation %>%
    dplyr::select(
      "module_eigengene", "ModuleColor", "Supermodule",
      "SupermoduleConfidence", "SupermoduleRationale", "top_GO_label",
      "present_in_dataset", "manual_annotation"
    )
  if (!is.null(color_col) && color_col %in% names(df)) {
    return(df %>% dplyr::left_join(ann, by = stats::setNames("ModuleColor", color_col)))
  }
  if (module_col %in% names(df)) {
    return(df %>% dplyr::left_join(ann, by = stats::setNames("module_eigengene", module_col)))
  }
  df
}

using_cached_final_state <- reuse_completed_analysis && !force_full_analysis && file.exists(wgcna_final_state_path)
if (!isTRUE(wgcna_dry_run) && using_cached_final_state) {
  message("Reusing completed WGCNA analysis from: ", wgcna_final_state_path)
  cached_state <- readRDS(wgcna_final_state_path)
  list2env(cached_state, envir = environment())
  if (exists("module_summary")) WGCNA_module_summary <- module_summary
  if (exists("module_preservation")) module_preservation_long <- module_preservation
  if (exists("module_label_table")) {
    module_label_table <- ensure_module_label_schema(module_label_table) %>%
      dplyr::mutate(
        ModuleLabel_Final = dplyr::coalesce(.data$ModuleLabel_Manual, .data$ModuleLabel_GO_BP, paste0("Module ", .data$ModuleColor)),
        ModuleLabel_Source = dplyr::case_when(
          !is.na(.data$ModuleLabel_Manual) & nzchar(.data$ModuleLabel_Manual) ~ "manual",
          !is.na(.data$ModuleLabel_GO_BP) & nzchar(.data$ModuleLabel_GO_BP) ~ "GO_BP_ORA_all_module",
          TRUE ~ "module_color"
        )
      )
  }
  if (exists("module_label_table")) {
    module_name_map <- stats::setNames(module_label_table$ModuleLabel_Final %||% module_label_table$ModuleLabel_GO_BP, module_label_table$ModuleColor)
  }
  dataset_profile_resolved <- cached_state$parameters$dataset_profile %||% dataset_profile
  if (!"condition" %in% names(sample_info)) sample_info$condition <- sample_info$ExpGroup
  if ("row.names" %in% names(sample_info)) {
    sample_info <- sample_info[match(rownames(expression.data), as.character(sample_info$row.names)), , drop = FALSE]
    rownames(sample_info) <- rownames(expression.data)
  }
  spatial_candidates <- intersect(c("celltype", "region", "layer"), names(sample_info))
  active_spatial_vars <- spatial_candidates[vapply(sample_info[spatial_candidates], function(x) {
    length(unique(stats::na.omit(as.character(x)))) > 1
  }, logical(1))]
  if (!length(active_spatial_vars) && "region" %in% names(sample_info)) active_spatial_vars <- "region"
  Samples <- rownames(expression.data)
  make_cached_trait_matrix <- function(sample_info, vars) {
    mats <- lapply(vars, function(v) {
      x <- droplevels(factor(sample_info[[v]]))
      if (length(unique(stats::na.omit(x))) <= 1) return(NULL)
      mm <- stats::model.matrix(~ 0 + x)
      colnames(mm) <- paste0(v, "_", sub("^x", "", colnames(mm)))
      mm
    })
    mats <- Filter(Negate(is.null), mats)
    if (!length(mats)) stop("No variable traits remain after filtering")
    as.data.frame(do.call(cbind, mats), stringsAsFactors = FALSE)
  }
  datTraits <- make_cached_trait_matrix(sample_info, c(active_spatial_vars, "condition"))
  keep_cols <- vapply(datTraits, function(x) stats::sd(as.numeric(x), na.rm = TRUE) > 0, logical(1))
  datTraits <- datTraits[, keep_cols, drop = FALSE]
  rownames(datTraits) <- Samples
  input_manifest <- tibble::tibble(
    role = "cached_wgcna_final_model_state",
    path = wgcna_final_state_path,
    md5 = file_hash(wgcna_final_state_path)
  )
} else {

expr_xlsx_env <- Sys.getenv("PROTEOMICS_WGCNA_EXPR_XLSX", unset = "")
meta_xlsx_env <- Sys.getenv("PROTEOMICS_WGCNA_META_XLSX", unset = "")
expr_xlsx_default <- path_processed(wgcna_module, "01_WGCNA", dataset_profile, "inputs", "wgcna_expression.xlsx")
meta_xlsx_default <- path_processed(wgcna_module, "01_WGCNA", dataset_profile, "inputs", "wgcna_sample_info.xlsx")
expr_xlsx <- if (nzchar(expr_xlsx_env)) expr_xlsx_env else expr_xlsx_default
meta_xlsx <- if (nzchar(meta_xlsx_env)) meta_xlsx_env else meta_xlsx_default
dataset_inputs <- resolve_dataset_inputs(dataset_profile, purpose = "wgcna")

# ================================
# Mouse-only mapping: robust idmapping parser + offline + SYMBOL/ALIAS + Entrez + UniProt gene_primary + QC
# ================================

# --------------------------
# Inputs
# --------------------------
expr_xlsx <- if (nzchar(expr_xlsx_env)) expr_xlsx_env else expr_xlsx
meta_xlsx <- if (nzchar(meta_xlsx_env)) meta_xlsx_env else meta_xlsx
idmap_dat <- Sys.getenv("PROTEOMICS_WGCNA_IDMAP_DAT", unset = path_external("MOUSE_10090_idmapping.dat"))

find_latest_wgcna_upstream <- function(profile = dataset_profile) {
  override <- Sys.getenv("PROTEOMICS_WGCNA_UPSTREAM_XLSX", unset = "")
  if (nzchar(override)) return(override)
  resolve_dataset_inputs(profile, purpose = "wgcna")$expression_file
}

first_existing_col <- function(df, candidates) {
  nms_clean <- tolower(gsub("[^a-z0-9]", "", names(df)))
  cand_clean <- tolower(gsub("[^a-z0-9]", "", candidates))
  hit <- match(cand_clean, nms_clean)
  hit <- hit[!is.na(hit)]
  if (!length(hit)) return(NA_character_)
  names(df)[hit[1]]
}

normalize_wgcna_group <- function(x) {
  x <- toupper(trimws(as.character(x)))
  dplyr::case_when(
    x %in% c("1", "CON", "CTRL", "CONTROL") ~ "con",
    x %in% c("2", "RES", "RESILIENT") ~ "res",
    x %in% c("3", "SUS", "SUSCEPTIBLE") ~ "sus",
    TRUE ~ tolower(x)
  )
}

profile_from_upstream_name <- function(path) {
  x <- tolower(basename(path))
  dplyr::case_when(
    grepl("microglia", x) ~ "microglia",
    grepl("neuron_soma", x) ~ "neuron_soma",
    grepl("neuron_neuropil", x) ~ "neuron_neuropil",
    TRUE ~ NA_character_
  )
}

auto_prepare_wgcna_inputs <- function(expr_path, meta_path) {
  if (file.exists(expr_path) && file.exists(meta_path)) return(invisible(FALSE))
  if (nzchar(expr_xlsx_env) || nzchar(meta_xlsx_env)) return(invisible(FALSE))

  upstream_xlsx <- find_latest_wgcna_upstream()
  upstream_meta <- dataset_inputs$metadata_file
  if (is.na(upstream_xlsx) || !file.exists(upstream_xlsx) || !file.exists(upstream_meta)) return(invisible(FALSE))

  message("Preparing missing WGCNA input workbooks from: ", upstream_xlsx)
  upstream <- as.data.frame(readxl::read_excel(upstream_xlsx))
  gene_col <- first_existing_col(upstream, c("T: Protein.Names", "gene_symbol", "Genes", "Protein.Group"))
  if (is.na(gene_col)) stop("Could not find a protein identifier column in upstream WGCNA input: ", upstream_xlsx, call. = FALSE)

  candidate_cols <- setdiff(names(upstream), c("Protein.Group", "T: Protein.Names", "Genes", "First.Protein.Description", gene_col))
  numeric_fraction <- vapply(candidate_cols, function(cc) {
    vals <- suppressWarnings(as.numeric(as.character(upstream[[cc]])))
    mean(!is.na(vals))
  }, numeric(1))
  sample_cols <- candidate_cols[numeric_fraction >= 0.70]
  if (length(sample_cols) < 4) stop("Could not detect expression sample columns in upstream WGCNA input: ", upstream_xlsx, call. = FALSE)

  expr_out <- upstream[, sample_cols, drop = FALSE]
  expr_out[] <- lapply(expr_out, function(x) suppressWarnings(as.numeric(as.character(x))))
  expr_out <- dplyr::bind_cols(tibble::tibble(gene_symbol = as.character(upstream[[gene_col]])), tibble::as_tibble(expr_out))

  metadata <- as.data.frame(readxl::read_excel(upstream_meta))
  sample_col <- first_existing_col(metadata, c("sample_id", "SampleID", "SampleColumn", "row.names"))
  if (is.na(sample_col)) stop("Could not find a sample identifier column in upstream metadata: ", upstream_meta, call. = FALSE)
  if ("exclude" %in% names(metadata)) metadata <- metadata[is.na(metadata$exclude) | metadata$exclude != TRUE, , drop = FALSE]

  metadata$.wgcna_sample <- as.character(metadata[[sample_col]])
  metadata <- metadata[metadata$.wgcna_sample %in% sample_cols, , drop = FALSE]
  metadata <- metadata[match(sample_cols, metadata$.wgcna_sample), , drop = FALSE]
  if (any(is.na(metadata$.wgcna_sample))) stop("Some WGCNA expression columns are missing from metadata.", call. = FALSE)

  source_profile <- profile_from_upstream_name(upstream_xlsx)
  metadata$celltype <- if (!is.na(source_profile)) {
    source_profile
  } else if ("celltype" %in% names(metadata)) {
    as.character(metadata$celltype)
  } else {
    NA_character_
  }
  region_col <- first_existing_col(metadata, c("region", "Region"))
  layer_col <- first_existing_col(metadata, c("layer", "Layer"))
  group_col <- first_existing_col(metadata, c("ExpGroup", "StressGroup", "group", "Group"))
  if (is.na(region_col)) stop("Could not find a region column in upstream metadata: ", upstream_meta, call. = FALSE)
  if (is.na(group_col)) stop("Could not find a condition/group column in upstream metadata: ", upstream_meta, call. = FALSE)
  meta_out <- tibble::tibble(
    row.names = metadata$.wgcna_sample,
    region = tolower(as.character(metadata[[region_col]])),
    layer = if (!is.na(layer_col)) tolower(as.character(metadata[[layer_col]])) else NA_character_,
    celltype = as.character(metadata$celltype),
    ExpGroup = normalize_wgcna_group(metadata[[group_col]])
  )

  dir.create(dirname(expr_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(meta_path), recursive = TRUE, showWarnings = FALSE)
  writexl::write_xlsx(expr_out, expr_path)
  writexl::write_xlsx(meta_out, meta_path)
  write_csv_safe(
    tibble::tibble(
      role = c("upstream_expression_matrix", "upstream_sample_metadata"),
      path = c(upstream_xlsx, upstream_meta),
      md5 = vapply(c(upstream_xlsx, upstream_meta), file_hash, character(1))
    ),
    fp_log("auto_prepared_input_manifest.csv")
  )
  invisible(TRUE)
}

stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Missing WGCNA input file: ", path,
      ". This script stages dataset-scoped WGCNA inputs under data/processed/06_modules_WGCNA/01_WGCNA/<dataset>/inputs. ",
      "If using custom inputs, set PROTEOMICS_WGCNA_EXPR_XLSX / PROTEOMICS_WGCNA_META_XLSX. Otherwise, ",
      "rerun 01_preprocessing/01_impute.r and ensure TPE9_sample_metadata_males.xlsx is available.",
      call. = FALSE
    )
  }
}

if (isTRUE(wgcna_dry_run)) {
  upstream_xlsx <- find_latest_wgcna_upstream()
  upstream_meta <- dataset_inputs$metadata_file
  can_stage_inputs <- file.exists(upstream_xlsx) && file.exists(upstream_meta)
  can_use_inputs <- file.exists(expr_xlsx) && file.exists(meta_xlsx)
  can_use_cache <- file.exists(wgcna_final_state_path) && !force_full_analysis
  can_write_outputs <- all(vapply(c(output_dir, unlist(subdirs)), function(path) dir.exists(path) && file.access(path, 2) == 0, logical(1)))
  downstream_contract <- fp_modtab("WGCNA_module_definitions_for_downstream.csv")
  sample_check <- "not checked"
  sample_check_ok <- TRUE
  if (can_use_inputs) {
    sample_check_ok <- tryCatch({
      expr_head <- readxl::read_excel(expr_xlsx, n_max = 3)
      meta_head <- readxl::read_excel(meta_xlsx, n_max = 1000)
      sample_col <- first_existing_col(meta_head, dataset_inputs$sample_id_col_candidates)
      if (is.na(sample_col)) {
        sample_check <<- "metadata sample column not detected"
        FALSE
      } else {
        expr_samples <- setdiff(names(expr_head), dataset_inputs$protein_id_col_candidates)
        meta_samples <- as.character(meta_head[[sample_col]])
        n_match <- length(intersect(expr_samples, meta_samples))
        sample_check <<- paste0(n_match, " matching sample IDs")
        n_match > 0
      }
    }, error = function(e) {
      sample_check <<- conditionMessage(e)
      FALSE
    })
  } else if (can_stage_inputs) {
    sample_check <- "staged files missing; canonical upstream exists and can be staged"
  } else {
    sample_check <- "staged files missing and canonical upstream cannot be staged"
    sample_check_ok <- FALSE
  }
  dry_run_line("Script", "06_modules_WGCNA/01_WGCNA.r")
  dry_run_line("Dataset", dataset_profile)
  dry_run_line("Resolved input diagnostics", paste(dataset_inputs$diagnostics, collapse = " | "))
  dry_run_line("Reuse cached final state", reuse_completed_analysis)
  dry_run_line("Force full analysis", force_full_analysis)
  dry_run_line("Cached state", wgcna_final_state_path, if (can_use_cache) "PASS" else if (reuse_completed_analysis) "WARN" else "INFO")
  dry_run_line("Explicit/staged expression workbook", expr_xlsx, if (file.exists(expr_xlsx)) "PASS" else "WARN")
  dry_run_line("Explicit/staged sample metadata workbook", meta_xlsx, if (file.exists(meta_xlsx)) "PASS" else "WARN")
  dry_run_line("Canonical upstream expression", upstream_xlsx, if (file.exists(upstream_xlsx)) "PASS" else "FAIL")
  dry_run_line("Canonical upstream metadata", upstream_meta, if (file.exists(upstream_meta)) "PASS" else "FAIL")
  dry_run_line("Staged input generation", if (can_use_inputs) "already staged" else if (can_stage_inputs) "can be generated" else "cannot be generated", if (can_use_inputs || can_stage_inputs) "PASS" else "FAIL")
  dry_run_line("Cheap sample matching", sample_check, if (sample_check_ok) "PASS" else "FAIL")
  dry_run_line("Mouse idmapping", idmap_dat, if (file.exists(idmap_dat)) "PASS" else "FAIL")
  dry_run_line("Output folders writable", paste(unlist(subdirs), collapse = "; "), if (can_write_outputs) "PASS" else "FAIL")
  dry_run_line("Downstream module contract", downstream_contract, if (file.exists(downstream_contract)) "PASS" else "WARN")
  dry_run_line("Optional DE/GSEA overlap bridge", repo_path("06_modules_WGCNA", "05_wgcna_de_gsea_overlap.r"), "WARN")
  quit(status = if (file.exists(idmap_dat) && can_write_outputs && sample_check_ok && (can_use_inputs || can_stage_inputs || can_use_cache)) 0 else 1, save = "no")
}

wgcna_inputs_auto_prepared <- auto_prepare_wgcna_inputs(expr_xlsx, meta_xlsx)

read_head <- function(path) {
  df <- readxl::read_excel(path)
  utils::write.table(utils::head(df, 10), fp_log(paste0(basename(path), "_head10.tsv")),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  df
}

norm_label <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[[:space:]-]+", "_", x)
  x[is.na(x) | !nzchar(x)] <- NA_character_
  x
}

infer_dataset_profile <- function(sample_info, expr_path, meta_path, requested = "auto") {
  requested <- tolower(requested)
  if (!identical(requested, "auto")) return(requested)
  path_hint <- tolower(paste(expr_path, meta_path, collapse = " "))
  if (grepl("microglia", path_hint)) return("microglia")
  if (grepl("neuropil", path_hint)) return("neuron_neuropil")
  if (grepl("soma", path_hint)) return("neuron_soma")
  cell_hint <- if ("celltype" %in% names(sample_info)) paste(unique(norm_label(sample_info$celltype)), collapse = " ") else ""
  if (grepl("microglia", cell_hint)) return("microglia")
  if (grepl("neuropil", cell_hint)) return("neuron_neuropil")
  if (grepl("soma", cell_hint)) return("neuron_soma")
  layer_vals <- if ("layer" %in% names(sample_info)) unique(stats::na.omit(norm_label(sample_info$layer))) else character()
  if (length(layer_vals) && any(layer_vals %in% c("so", "sr", "slm", "mo", "po"))) return("neuron_neuropil")
  if (!length(layer_vals) || all(layer_vals %in% c("sp", "sg", "none"))) return("neuron_soma")
  "region_only"
}

prepare_spatial_metadata <- function(sample_info, profile) {
  if (!"region" %in% names(sample_info)) stop("sample_info must contain a 'region' column")
  if (!"ExpGroup" %in% names(sample_info)) stop("sample_info must contain an 'ExpGroup' column")

  profile <- tolower(profile)
  allowed_profiles <- c("microglia", "neuron_soma", "neuron_neuropil", "region_only")
  if (!profile %in% allowed_profiles) {
    stop("Unsupported dataset_profile: ", profile, ". Use one of: ", paste(c("auto", allowed_profiles), collapse = ", "))
  }
  sample_info$region <- factor(norm_label(sample_info$region), levels = c("ca1", "ca2", "ca3", "dg"))
  sample_info$condition <- factor(norm_label(sample_info$ExpGroup), levels = c("con", "res", "sus"))
  sample_info$ExpGroup <- sample_info$condition

  raw_layer <- if ("layer" %in% names(sample_info)) norm_label(sample_info$layer) else rep(NA_character_, nrow(sample_info))
  raw_celltype <- if ("celltype" %in% names(sample_info)) norm_label(sample_info$celltype) else rep(NA_character_, nrow(sample_info))

  if (identical(profile, "microglia")) {
    sample_info$layer <- factor("none")
    sample_info$celltype <- factor("microglia")
    active_spatial_vars <- c("region")
  } else if (identical(profile, "neuron_soma")) {
    sample_info$soma_layer <- factor(raw_layer, levels = c("sp", "sg"))
    sample_info$layer <- factor("none")
    sample_info$celltype <- factor("neuron_soma")
    active_spatial_vars <- c("region")
  } else if (identical(profile, "neuron_neuropil")) {
    sample_info$layer <- factor(raw_layer, levels = c("so", "sr", "slm", "mo", "po"))
    sample_info$celltype <- factor("neuron_neuropil")
    active_spatial_vars <- c("region", "layer")
    invalid_layer <- (!is.na(sample_info$region) & !is.na(sample_info$layer)) & (
      (sample_info$region %in% c("ca1", "ca2", "ca3") & !sample_info$layer %in% c("so", "sr", "slm")) |
        (sample_info$region == "dg" & !sample_info$layer %in% c("mo", "po"))
    )
    if (any(invalid_layer)) {
      write_csv_safe(sample_info[invalid_layer, , drop = FALSE], fp_log("invalid_region_layer_combinations.csv"))
      warning("Invalid region/layer combinations found; see logs/invalid_region_layer_combinations.csv")
    }
  } else {
    sample_info$layer <- factor(ifelse(is.na(raw_layer), "none", raw_layer))
    sample_info$celltype <- factor(ifelse(is.na(raw_celltype), profile, raw_celltype))
    active_spatial_vars <- c("region")
    if (length(unique(stats::na.omit(sample_info$layer))) > 1) active_spatial_vars <- c(active_spatial_vars, "layer")
    if (length(unique(stats::na.omit(sample_info$celltype))) > 1) active_spatial_vars <- c(active_spatial_vars, "celltype")
  }

  required_vars <- c("condition", active_spatial_vars)
  missing_required <- Reduce(`|`, lapply(required_vars, function(v) is.na(sample_info[[v]])))
  if (any(missing_required)) {
    write_csv_safe(sample_info[missing_required, , drop = FALSE], fp_log("invalid_required_metadata.csv"))
    stop("Missing or invalid values in required metadata columns: ",
         paste(required_vars, collapse = ", "),
         ". See logs/invalid_required_metadata.csv")
  }

  write_csv_safe(
    tibble::tibble(
      dataset_profile = profile,
      active_spatial_vars = paste(active_spatial_vars, collapse = ","),
      n_samples = nrow(sample_info),
      regions = paste(sort(unique(stats::na.omit(as.character(sample_info$region)))), collapse = ","),
      layers = paste(sort(unique(stats::na.omit(as.character(sample_info$layer)))), collapse = ",")
    ),
    fp_log("dataset_profile.csv")
  )

  list(sample_info = sample_info, profile = profile, active_spatial_vars = active_spatial_vars)
}

male.data <- { stop_if_missing(expr_xlsx); read_head(expr_xlsx) } %>% dplyr::mutate(.row_id = dplyr::row_number())
meta.data <- { stop_if_missing(meta_xlsx); read_head(meta_xlsx) }

# --------------------------
# Robust UniProtKB-ID -> Accession map (offline)
# --------------------------
stop_if_missing(idmap_dat)
input_manifest <- tibble::tibble(
  role = c(
    "canonical_upstream_expression_matrix",
    "canonical_upstream_sample_metadata",
    "staged_wgcna_expression_matrix",
    "staged_wgcna_sample_metadata",
    "expression_matrix",
    "sample_metadata",
    "mouse_idmapping"
  ),
  path = c(
    dataset_inputs$expression_file,
    dataset_inputs$metadata_file,
    expr_xlsx_default,
    meta_xlsx_default,
    expr_xlsx,
    meta_xlsx,
    idmap_dat
  ),
  md5 = vapply(c(
    dataset_inputs$expression_file,
    dataset_inputs$metadata_file,
    expr_xlsx_default,
    meta_xlsx_default,
    expr_xlsx,
    meta_xlsx,
    idmap_dat
  ), file_hash, character(1))
)
if (isTRUE(wgcna_inputs_auto_prepared) && file.exists(fp_log("auto_prepared_input_manifest.csv"))) {
  input_manifest <- dplyr::bind_rows(
    input_manifest,
    readr::read_csv(fp_log("auto_prepared_input_manifest.csv"), show_col_types = FALSE)
  )
}
write_csv_safe(input_manifest, fp_log("input_manifest.csv"))
idmap_tbl <- readr::read_tsv(idmap_dat, col_names = c("ACC","DB","VAL"), col_types = "ccc", progress = FALSE, quote = "", comment = "")
idmap_uid <- idmap_tbl %>%
  dplyr::filter(DB == "UniProtKB-ID" & grepl("_MOUSE\\s*$", VAL) & nzchar(ACC)) %>%
  dplyr::transmute(
    UNIPROT    = toupper(trimws(ACC)),
    entry_full = toupper(trimws(VAL)),
    entry_base = toupper(gsub("_MOUSE$", "", trimws(VAL)))
  )
entry_map <- idmap_uid %>% dplyr::distinct(entry_base, .keep_all = TRUE)
if (!nrow(entry_map)) stop("entry_map is empty after robust parse")

# Sentinel sanity check
sentinels <- c("AIF1","AKAP2","ADCY1","AKAP1","AMPD3","ANXA3","1433S","ACK1","AIP","ADA10")
missing_sentinels <- setdiff(sentinels, entry_map$entry_base)
if (length(missing_sentinels)) {
  warning(sprintf("entry_map missing expected keys: %s", paste(missing_sentinels, collapse=", ")))
  write_tsv_safe(entry_map, fp_maptab("entry_map_debug.tsv"))
}

# --------------------------
# Mouse-only tokenization and classification
# --------------------------
normalize_token <- function(x) { x <- toupper(gsub("\\s+", "", x)); x <- gsub("\\u00A0", "", x); x <- gsub("\\.+", ".", x); x <- gsub("__+", "_", x); x }
to_base_no_iso_mouse <- function(x) { x <- gsub("-\\d+$", "", x); gsub("_MOUSE$", "", x) }

tokenize_mouse_only <- function(male_df) {
  tok <- male_df %>% tidyr::separate_rows(gene_symbol, sep = ";") %>% dplyr::mutate(token_raw = gene_symbol, token_up = normalize_token(gene_symbol))
  dropped_non_mouse <- tok %>% dplyr::filter(!grepl("_MOUSE$", token_up))
  if (nrow(dropped_non_mouse)) write_tsv_safe(dropped_non_mouse, fp_maptab("dropped_non_mouse_tokens.tsv"))
  tok %>%
    dplyr::filter(grepl("_MOUSE$", token_up)) %>%
    dplyr::mutate(
      token_base = to_base_no_iso_mouse(token_up),
      looks_ac   = grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$|^A0A[0-9A-Z]{7}$", token_base),
      looks_entry= grepl("^[A-Z0-9][A-Z0-9\\-\\.]+$", token_base),
      id_class = dplyr::case_when(
        looks_ac    ~ "UNIPROT_AC_MOUSE",
        looks_entry ~ "UNIPROT_ENTRY",
        TRUE        ~ "UNKNOWN"
      ),
      Resolved_UNIPROT = NA_character_,
      strategy = NA_character_
    )
}
resolved2 <- tokenize_mouse_only(male.data)

# --------------------------
# Mapping stack
# --------------------------
empty_id_map <- function() {
  tibble::tibble(input = character(), primaryAccession = character())
}

# 1) Accept accession-like bases
idx_ac <- which(resolved2$id_class == "UNIPROT_AC_MOUSE")
if (length(idx_ac)) {
  resolved2$Resolved_UNIPROT[idx_ac] <- resolved2$token_base[idx_ac]
  resolved2$strategy[idx_ac] <- "accept_accession_base"
}

# 2) Offline entry map (now robust)
idx_en <- which(resolved2$id_class == "UNIPROT_ENTRY" & (is.na(resolved2$Resolved_UNIPROT) | !nzchar(resolved2$Resolved_UNIPROT)))
if (length(idx_en)) {
  hit <- entry_map$UNIPROT[match(toupper(resolved2$token_base[idx_en]), entry_map$entry_base)]
  ok <- !is.na(hit) & nzchar(hit)
  if (any(ok)) { ii <- idx_en[ok]; resolved2$Resolved_UNIPROT[ii] <- hit[ok]; resolved2$strategy[ii] <- "entry_local_mouse" }
}

# 3) SYMBOL/ALIAS offline resolver (MGI-first)
need_idx <- which(is.na(resolved2$Resolved_UNIPROT) | !nzchar(resolved2$Resolved_UNIPROT))
need_ids <- toupper(unique(resolved2$token_base[need_idx]))
is_acc <- grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$|^A0A[0-9A-Z]{7}$", need_ids)
ids_ent <- unique(need_ids[!is_acc])

if (length(ids_ent)) {
  sel_sym <- try(AnnotationDbi::select(org.Mm.eg.db, keys = ids_ent, keytype = "SYMBOL", columns = c("MGIID","ENTREZID","UNIPROT","SYMBOL")), silent = TRUE)
  map_sym <- empty_id_map()
  if (!inherits(sel_sym, "try-error") && nrow(sel_sym)) {
    map_sym <- tibble::as_tibble(sel_sym) %>%
      dplyr::filter(!is.na(UNIPROT) & nzchar(UNIPROT)) %>%
      dplyr::group_by(SYMBOL) %>% dplyr::arrange(UNIPROT, .by_group = TRUE) %>% dplyr::slice_head(n=1) %>% dplyr::ungroup() %>%
      dplyr::transmute(input = toupper(SYMBOL), primaryAccession = toupper(UNIPROT))
  }
  kt <- try(AnnotationDbi::keytypes(org.Mm.eg.db), silent = TRUE)
  map_alias <- empty_id_map()
  if (!inherits(kt, "try-error") && "ALIAS" %in% kt) {
    sel_alias <- try(AnnotationDbi::select(org.Mm.eg.db, keys = ids_ent, keytype = "ALIAS", columns = c("UNIPROT","ALIAS")), silent = TRUE)
    if (!inherits(sel_alias, "try-error") && nrow(sel_alias)) {
      map_alias <- tibble::as_tibble(sel_alias) %>%
        dplyr::filter(!is.na(UNIPROT) & nzchar(UNIPROT)) %>%
        dplyr::group_by(ALIAS) %>% dplyr::arrange(UNIPROT, .by_group = TRUE) %>% dplyr::slice_head(n=1) %>% dplyr::ungroup() %>%
        dplyr::transmute(input = toupper(ALIAS), primaryAccession = toupper(UNIPROT))
    }
  }
  map_symall <- dplyr::bind_rows(map_sym, map_alias) %>% dplyr::distinct(input, .keep_all = TRUE)
  if (nrow(map_symall)) {
    need_idx_now <- which(is.na(resolved2$Resolved_UNIPROT) | !nzchar(resolved2$Resolved_UNIPROT))
    base_need <- toupper(resolved2$token_base[need_idx_now])
    hit <- map_symall$primaryAccession[match(base_need, map_symall$input)]
    ok <- !is.na(hit) & nzchar(hit)
    ii <- need_idx_now[ok]
    if (length(ii)) { resolved2$Resolved_UNIPROT[ii] <- hit[ok]; resolved2$strategy[ii] <- "orgdb_mgi_symbol_first" }
  }
}

# 4) SYMBOL -> Entrez -> UniProt (offline two-hop)
need_idx <- which(is.na(resolved2$Resolved_UNIPROT) | !nzchar(resolved2$Resolved_UNIPROT))
need_ids <- toupper(unique(resolved2$token_base[need_idx]))
sym_left <- unique(need_ids[grepl("^[A-Z0-9\\-]{2,}$", need_ids)])

if (length(sym_left)) {
  sym2eg <- try(AnnotationDbi::select(org.Mm.eg.db, keys = sym_left, keytype = "SYMBOL", columns = c("ENTREZID","SYMBOL")), silent = TRUE)
  eg2up  <- tibble::tibble(ENTREZID = character(), UNIPROT = character())
  if (!inherits(sym2eg, "try-error") && nrow(sym2eg)) {
    ekeys <- unique(na.omit(sym2eg$ENTREZID))
    if (length(ekeys)) {
      egsel <- try(AnnotationDbi::select(org.Mm.eg.db, keys = ekeys, keytype = "ENTREZID", columns = c("UNIPROT","ENTREZID")), silent = TRUE)
      if (!inherits(egsel, "try-error") && nrow(egsel)) {
        eg2up <- tibble::as_tibble(egsel) %>%
          dplyr::filter(!is.na(UNIPROT) & nzchar(UNIPROT)) %>%
          dplyr::group_by(ENTREZID) %>% dplyr::arrange(UNIPROT, .by_group = TRUE) %>% dplyr::slice_head(n=1) %>% dplyr::ungroup()
      }
    }
    if (nrow(eg2up)) {
      map_sym2up <- tibble::as_tibble(sym2eg) %>%
        dplyr::distinct(SYMBOL, ENTREZID) %>%
        dplyr::left_join(eg2up, by = "ENTREZID") %>%
        dplyr::filter(!is.na(UNIPROT) & nzchar(UNIPROT)) %>%
        dplyr::transmute(input = toupper(SYMBOL), primaryAccession = toupper(UNIPROT)) %>%
        dplyr::distinct(input, .keep_all = TRUE)
      need_idx2 <- which(is.na(resolved2$Resolved_UNIPROT) | !nzchar(resolved2$Resolved_UNIPROT))
      base_need <- toupper(resolved2$token_base[need_idx2])
      hit <- map_sym2up$primaryAccession[match(base_need, map_sym2up$input)]
      ok <- !is.na(hit) & nzchar(hit)
      ii <- need_idx2[ok]
      if (length(ii)) { resolved2$Resolved_UNIPROT[ii] <- hit[ok]; resolved2$strategy[ii] <- "orgdb_symbol_entrez_uniprot" }
    }
  }
}

# 5) UniProt gene_primary resolver (Mus musculus), batched with retry, prefer reviewed
need_idx <- which(is.na(resolved2$Resolved_UNIPROT) | !nzchar(resolved2$Resolved_UNIPROT))
need_ids <- unique(toupper(resolved2$token_base[need_idx]))
sym_left2 <- unique(need_ids[grepl("^[A-Z0-9\\-]{2,}$", need_ids)])

if (length(sym_left2)) {
  batch_vec <- split(sym_left2, ceiling(seq_along(sym_left2)/50))
  picks <- list()
  for (b in batch_vec) {
    q_list <- lapply(b, function(g) list(organism_id = 10090, gene_primary = g))
    query_once <- function(ql) try(UniProt.ws::queryUniProt(query = ql, fields = c("accession","id","gene_primary","reviewed"), collapse = "OR", n = 10, pageSize = 10), silent = TRUE)
    res_list <- lapply(q_list, function(ql) { out <- query_once(ql); if (inherits(out, "try-error") || !is.data.frame(out)) { Sys.sleep(0.8); out <- query_once(ql) }; out })
    ok <- res_list[!vapply(res_list, inherits, logical(1), "try-error")]
    if (length(ok)) {
      tbl <- dplyr::bind_rows(lapply(ok, tibble::as_tibble))
      if (nrow(tbl)) {
        tbl <- tbl %>% dplyr::mutate(gene_primary = toupper(.data$gene_primary), accession = toupper(.data$accession))
        pick <- tbl %>% dplyr::group_by(gene_primary) %>% dplyr::arrange(dplyr::desc(.data$reviewed), accession, .by_group = TRUE) %>% dplyr::slice_head(n=1) %>% dplyr::ungroup() %>% dplyr::transmute(input = gene_primary, primaryAccession = accession)
        picks[[length(picks)+1]] <- pick
      }
    }
  }
  if (length(picks)) {
    map_gene <- dplyr::bind_rows(c(list(empty_id_map()), picks)) %>% dplyr::distinct(input, .keep_all = TRUE)
    need_idx3 <- which(is.na(resolved2$Resolved_UNIPROT) | !nzchar(resolved2$Resolved_UNIPROT))
    base_need <- toupper(resolved2$token_base[need_idx3])
    hit <- map_gene$primaryAccession[match(base_need, map_gene$input)]
    ok <- !is.na(hit) & nzchar(hit)
    ii <- need_idx3[ok]
    if (length(ii)) { resolved2$Resolved_UNIPROT[ii] <- hit[ok]; resolved2$strategy[ii] <- "uniprot_gene_primary_retry" }
  }
}

# --------------------------
# Save unmapped lists and audits
# --------------------------
audit_mouse <- resolved2 %>%
  dplyr::mutate(mapped = !is.na(Resolved_UNIPROT) & nzchar(Resolved_UNIPROT)) %>%
  dplyr::count(id_class, strategy, mapped, name = "n") %>%
  dplyr::arrange(desc(n))
write_tsv_safe(audit_mouse, fp_maptab("mapping_audit_mouse_only_robust.tsv"))

unmapped_tokens <- resolved2 %>%
  dplyr::filter(is.na(Resolved_UNIPROT) | !nzchar(Resolved_UNIPROT)) %>%
  dplyr::transmute(
    .row_id, token_raw, token_up, token_base, id_class,
    reason = dplyr::case_when(
      grepl("[^A-Z0-9_\\-\\.]", token_base) ~ "illegal_chars",
      grepl("^[A-Z0-9\\-\\.]+$", token_base) ~ "entry_name_not_in_local_or_query",
      TRUE ~ "unexpected_format_or_na"
    )
  ) %>% dplyr::arrange(id_class, token_base)
write_tsv_safe(unmapped_tokens, fp_maptab("unmapped_mouse_tokens.tsv"))

unmapped_summary <- unmapped_tokens %>% dplyr::count(id_class, reason, token_base, name = "n") %>% dplyr::arrange(dplyr::desc(n))
write_tsv_safe(unmapped_summary, fp_maptab("unmapped_mouse_tokens_summary.tsv"))

# --------------------------
# Collapse to features and build expression matrix
# --------------------------
collapse_ids <- function(x) { x <- unique(x[!is.na(x) & nzchar(x)]); if (!length(x)) return(NA_character_); paste(x, collapse = ";") }
male.norm <- resolved2 %>%
  dplyr::group_by(.row_id) %>%
  dplyr::summarise(gene_symbol = collapse_ids(Resolved_UNIPROT), .groups = "drop") %>%
  dplyr::right_join(male.data %>% dplyr::select(-gene_symbol, .row_id), by = ".row_id") %>%
  dplyr::select(-.row_id) %>%
  dplyr::mutate(gene_symbol = dplyr::na_if(gene_symbol, ""))

to_numeric_matrix <- function(male_norm, qc_dir = subdirs$logs) {
  if (!"gene_symbol" %in% names(male_norm)) stop("male.norm must contain gene_symbol")
  expr <- as.data.frame(
    lapply(male_norm[, -1, drop = FALSE], function(x) suppressWarnings(as.numeric(x))),
    check.names = FALSE
  )
  if (!all(vapply(expr, is.numeric, logical(1)))) stop("Non-numeric columns remain after coercion")
  mat <- as.data.frame(t(expr))
  feat <- male_norm$gene_symbol
  empty <- which(!nzchar(ifelse(is.na(feat), "", feat)))
  if (length(empty)) feat[empty] <- paste0("UNMAPPED_", seq_along(empty))
  feat <- make.unique(feat, sep = "_")
  colnames(mat) <- feat
  if (any(!nzchar(colnames(mat)) | is.na(colnames(mat)))) stop("Empty/NA feature names after repair")
  utils::write.table(utils::head(mat[, 1:min(10, ncol(mat)), drop = FALSE]), file.path(qc_dir, "expression_head10.tsv"), sep = "\t", row.names = TRUE, quote = FALSE)
  mat
}
expression.data <- to_numeric_matrix(male.norm)

# Force single accession per feature name and ensure uniqueness
fix_feature_ids <- function(nms) {
  first <- sub(";.*$", "", nms)
  first <- toupper(trimws(first))
  first[is.na(first) | !nzchar(first)] <- "UNMAPPED"
  make.unique(first, sep = "_")
}

colnames(expression.data) <- fix_feature_ids(colnames(expression.data))

# Save core outputs
write_tsv_safe(resolved2, fp_maptab("resolved_tokens_mouse_only_robust.tsv"))
saveRDS(list(expression = expression.data, male.norm = male.norm, mapping = resolved2),
        file = fp_state("mouse_only_mapping_outputs_robust.rds"))

# --------------------------
# QC and sample clustering
# --------------------------
gsg <- goodSamplesGenes(expression.data)
sample_qc <- tibble::tibble(
  sample = rownames(expression.data),
  good_samples_genes = as.logical(gsg$goodSamples)
)
gene_qc <- tibble::tibble(
  gene = colnames(expression.data),
  good_samples_genes = as.logical(gsg$goodGenes)
)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) printFlush(paste("Removing genes:", paste(names(expression.data)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples) > 0) printFlush(paste("Removing samples:", paste(rownames(expression.data)[!gsg$goodSamples], collapse = ", ")))
  expression.data <- expression.data[gsg$goodSamples == TRUE, gsg$goodGenes == TRUE]
}

sampleTree <- hclust(dist(expression.data), method = "average")
svg(file = fp_qc("sample_clustering_outliers.svg"), width = figure_single_col * 1.4, height = 3.8, family = figure_font)
par(cex = 0.6, mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.main = 2)
abline(h = sample_tree_plot_height, col = "red")
dev.off()

cut.sampleTree <- cutreeStatic(sampleTree, cutHeight = sample_tree_cut_height, minSize = 10)
sample_tree_cluster <- stats::setNames(as.integer(cut.sampleTree), rownames(expression.data))
sample_qc <- sample_qc %>%
  dplyr::mutate(
    sample_tree_cluster = unname(sample_tree_cluster[.data$sample]),
    retained_after_sample_tree = .data$good_samples_genes & .data$sample_tree_cluster == 1
  )
write_csv_safe(sample_qc, fp_qctab("sample_filtering_qc.csv"))
write_csv_safe(gene_qc, fp_qctab("gene_filtering_qc.csv"))
expression.data <- expression.data[cut.sampleTree == 1, ]

# --------------------------
# Soft-threshold selection
# --------------------------
spt <- pickSoftThreshold(expression.data, networkType = "signed",
                         corFnc = "bicor",
                         corOptions = list(use = "p", maxPOutliers = 0.05))

svglite::svglite(file = fp_qc("soft_threshold_scale_independence.svg"), width = figure_single_col, height = 2.8)
par(mar = c(4,4,2,1))
plot(spt$fitIndices[,1], spt$fitIndices[,2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n", main = "Scale independence")
text(spt$fitIndices[,1], spt$fitIndices[,2], labels = spt$fitIndices[,1], col = "red")
abline(h = soft_threshold_rsquared, col = "red")
dev.off()

svglite::svglite(file = fp_qc("soft_threshold_mean_connectivity.svg"), width = figure_single_col, height = 2.8)
par(mar = c(4,4,2,1))
plot(spt$fitIndices[,1], spt$fitIndices[,5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = "Mean connectivity")
text(spt$fitIndices[,1], spt$fitIndices[,5], labels = spt$fitIndices[,1], col = "red")
dev.off()

fit_indices <- tibble::as_tibble(spt$fitIndices)
eligible_power <- fit_indices$Power[fit_indices$SFT.R.sq >= soft_threshold_rsquared]
softPower <- if (length(eligible_power)) {
  min(eligible_power)
} else {
  fit_indices$Power[which.max(fit_indices$SFT.R.sq)]
}
write_csv_safe(fit_indices, fp_qctab("soft_threshold_fit_indices.csv"))
write_csv_safe(
  tibble::tibble(
    parameter = c("soft_threshold_rsquared", "selected_soft_power"),
    value = c(soft_threshold_rsquared, softPower)
  ),
  fp_qctab("soft_threshold_selection.csv")
)

# --------------------------
# Network construction
# --------------------------
adjacency <- adjacency(expression.data, power = softPower, type = "signed",
                       corFnc = "bicor", corOptions = list(use="p", maxPOutliers=0.05))
TOM <- TOMsimilarity(adjacency)
TOM.dissimilarity <- 1 - TOM
geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average")

svg(file = fp_net("gene_dendrogram.svg"), width = figure_double_col, height = 5.2, family = figure_font)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

Modules <- cutreeDynamic(dendro = geneTree, distM = TOM.dissimilarity, deepSplit = deep_split,
                         pamRespectsDendro = FALSE, minClusterSize = min_module_size)

colorSeq <- c(
  "lemon" = "lemonchiffon", "sage" = "darkseagreen", "bluegray" = "steelblue",
  "mintblue" = "lightsteelblue", "azure" = "deepskyblue", "khaki" = "khaki",
  "skyblue" = "skyblue", "babyblue" = "lightblue", "amber" = "goldenrod",
  "tealgreen" = "darkcyan", "forestgreen" = "forestgreen", "gold" = "gold",
  "violet" = "violet", "seafoam" = "mediumaquamarine", "coral" = "coral",
  "salmonlight" = "lightsalmon", "peach" = "peachpuff", "mint" = "palegreen",
  "lime" = "limegreen", "mauve" = "plum", "freesia" = "lightpink", "cocoa" = "saddlebrown", 
  "lavender" = "lavender", "magenta" = "magenta", "salmon" = "salmon",
  "rose" = "mistyrose", "aquamarine" = "aquamarine", "tomato" = "tomato",
  "plum" = "plum", "hotpink" = "hotpink", "rust" = "sienna"
)

# Map numeric module labels -> colors from colorSeq (recycle palette if needed)
unique_mods <- sort(unique(Modules))
nmods <- length(unique_mods)
palette_vals <- unname(colorSeq)
if (length(palette_vals) < nmods) palette_vals <- rep(palette_vals, length.out = nmods)
mod_colors_map <- setNames(palette_vals[seq_len(nmods)], as.character(unique_mods))
ModuleColors <- as.character(mod_colors_map[as.character(Modules)])
stopifnot(length(ModuleColors) == length(Modules))

svg(file = fp_net("gene_dendrogram_module_colors.svg"), width = figure_double_col, height = 5.2, family = figure_font)
plotDendroAndColors(geneTree, ModuleColors, "Module", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

MElist <- moduleEigengenes(expression.data, colors = ModuleColors)
MEs <- MElist$eigengenes
ME.dissimilarity <- 1 - cor(MEs, use = "p", method = "pearson")
METree <- hclust(as.dist(ME.dissimilarity), method = "average")
merge <- mergeCloseModules(expression.data, ModuleColors, cutHeight = merge_cut_height)
mergedColors <- merge$colors
names(mergedColors) <- colnames(expression.data)
mergedMEs <- orderMEs(merge$newMEs)

write_csv_safe(
  tibble::tibble(
    parameter = c(
      "sample_tree_cut_height", "sample_tree_plot_height",
      "soft_threshold_rsquared", "selected_soft_power",
      "network_type", "correlation_function", "bicor_maxPOutliers",
      "deep_split", "min_module_size", "merge_cut_height",
      "module_preservation_permutations", "dataset_profile_requested",
      "output_layout"
    ),
    value = c(
      sample_tree_cut_height, sample_tree_plot_height,
      soft_threshold_rsquared, softPower,
      "signed", "bicor", 0.05,
      deep_split, min_module_size, merge_cut_height,
      module_preservation_permutations, dataset_profile,
      ifelse(nzchar(output_dir_env), "custom_bundle", "canonical_module_paths")
    )
  ),
  fp_log("analysis_parameters.csv")
)

svg(file = fp_net("gene_dendrogram_modules_merged.svg"), width = figure_double_col, height = 5.2, family = figure_font)
plotDendroAndColors(geneTree, cbind(ModuleColors, mergedColors),
                    c("Original Module","Merged Module"),
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors for original and merged modules")
dev.off()

# --------------------------
# Eigengene network plots
# --------------------------
MET <- orderMEs(mergedMEs)
svg(file = fp_net("eigengene_dendrogram.svg"), width = figure_single_col, height = 3.0, family = figure_font)
plotEigengeneNetworks(MET, "", plotHeatmaps = FALSE, marDendro = c(0, 4, 2, 0))
dev.off()
svg(file = fp_net("eigengene_adjacency_heatmap.svg"), width = figure_single_col, height = 3.0, family = figure_font)
par(mar = c(1,1,1,1))
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", plotDendrograms = FALSE,
                      marHeatmap = c(5, 5, 2, 2), xLabelsAngle = 90)
dev.off()

# --------------------------
# Spatial + condition trait matrix and heatmaps
# --------------------------
sample_info <- as.data.frame(readxl::read_excel(path = meta_xlsx))
stopifnot("row.names" %in% names(sample_info))

sample_key_for_match <- function(x) {
  x <- tolower(trimws(as.character(x)))
  gsub("[^a-z0-9]+", "", x)
}

Samples <- rownames(expression.data)
sample_match <- match(Samples, as.character(sample_info$row.names))
if (anyNA(sample_match)) {
  sample_key <- sample_key_for_match(Samples)
  metadata_key <- sample_key_for_match(sample_info$row.names)
  duplicated_metadata_key <- duplicated(metadata_key) | duplicated(metadata_key, fromLast = TRUE)
  if (any(duplicated_metadata_key)) {
    duplicated_metadata <- sample_info[duplicated(metadata_key) | duplicated(metadata_key, fromLast = TRUE), , drop = FALSE]
    write_csv_safe(duplicated_metadata, fp_log("metadata_duplicate_normalized_sample_keys.csv"))
  }
  metadata_key_unique <- metadata_key
  metadata_key_unique[duplicated_metadata_key] <- NA_character_
  fallback_match <- match(sample_key, metadata_key_unique)
  sample_match[is.na(sample_match)] <- fallback_match[is.na(sample_match)]
}
if (anyNA(sample_match)) {
  missing_samples <- Samples[is.na(sample_match)]
  write_csv_safe(
    tibble::tibble(sample = missing_samples),
    fp_log("samples_missing_from_metadata.csv")
  )
  stop(
    "Metadata is missing ", length(missing_samples), " expression samples. ",
    "See logs/samples_missing_from_metadata.csv",
    call. = FALSE
  )
}
sample_info <- sample_info[sample_match, , drop = FALSE]
rownames(sample_info) <- as.character(sample_info$row.names)
profile_info <- prepare_spatial_metadata(
  sample_info,
  infer_dataset_profile(sample_info, expr_xlsx, meta_xlsx, dataset_profile)
)
sample_info <- profile_info$sample_info
dataset_profile_resolved <- profile_info$profile
active_spatial_vars <- profile_info$active_spatial_vars

make_trait_matrix <- function(sample_info, vars) {
  mats <- lapply(vars, function(v) {
    x <- droplevels(factor(sample_info[[v]]))
    if (length(unique(stats::na.omit(x))) <= 1) return(NULL)
    mm <- model.matrix(~ 0 + x)
    colnames(mm) <- paste0(v, "_", sub("^x", "", colnames(mm)))
    mm
  })
  mats <- Filter(Negate(is.null), mats)
  if (!length(mats)) stop("No variable traits remain after filtering")
  as.data.frame(do.call(cbind, mats), stringsAsFactors = FALSE)
}

datTraits <- make_trait_matrix(sample_info, c(active_spatial_vars, "condition"))
keep_cols <- vapply(datTraits, function(x) sd(as.numeric(x), na.rm = TRUE) > 0, logical(1))
datTraits <- datTraits[, keep_cols, drop = FALSE]
rownames(datTraits) <- Samples

nSamples <- nrow(expression.data)
MEcorr <- cor(mergedMEs, datTraits, use = "p", method = "pearson")
MEp    <- corPvalueStudent(MEcorr, nSamples)
MEfdr  <- matrix(
  p.adjust(as.vector(MEp), method = "BH"),
  nrow = nrow(MEp), ncol = ncol(MEp), dimnames = dimnames(MEp)
)

plot_trait_heatmap <- function(matCorr, matFdr, cols, file) {
  if (length(cols) == 0) return(invisible(NULL))
  textMatrix <- paste(signif(matCorr[, cols, drop=FALSE], 2), "\nq=",
                      signif(matFdr[, cols, drop=FALSE], 1), sep = "")
  dim(textMatrix) <- dim(matCorr[, cols, drop=FALSE])
  svg(file = file, width = 6, height = max(4, nrow(matCorr) * 0.25 + 2))
  par(mar = c(6, 8.5, 3, 1))
  labeledHeatmap(Matrix = matCorr[, cols, drop=FALSE],
                 xLabels = colnames(matCorr)[cols],
                 yLabels = rownames(matCorr),
                 ySymbols = rownames(matCorr),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.6,
                 zlim = c(-1, 1),
                 main = "Module-trait relationships")
  dev.off()
}

trait_names <- colnames(datTraits)
groups <- list(
  celltype = grep("^celltype_", trait_names),
  layer    = grep("^layer_",   trait_names),
  region   = grep("^region_",  trait_names),
  condition = grep("^condition_", trait_names)
)
groups <- groups[c(intersect(c("celltype", "region", "layer"), active_spatial_vars), "condition")]
for (nm in names(groups)) {
  idx <- groups[[nm]]
  if (length(idx) > 0) {
    plot_trait_heatmap(MEcorr, MEfdr, idx, fp_traits(paste0("ME_trait_heatmap_", nm, ".svg")))
  }
}
write_csv_safe(
  reshape2::melt(MEcorr, varnames = c("module", "trait"), value.name = "r") |>
    dplyr::mutate(
      p = as.vector(MEp),
      fdr = as.vector(MEfdr)
    ),
  fp_source("ME_trait_correlations.csv")
)

# --------------------------
# Pairwise condition contrasts (optional)
# --------------------------
mk_contrast <- function(vec, a, b){v<-rep(NA_real_,length(vec));v[vec==a]<-0;v[vec==b]<-1;v}
grp <- as.character(sample_info$ExpGroup)
contrasts <- list(con_res = mk_contrast(grp, "con", "res"),
                  con_sus = mk_contrast(grp, "con", "sus"),
                  res_sus = mk_contrast(grp, "res", "sus"))
for (nm in names(contrasts)) {
  v <- contrasts[[nm]]; keep <- !is.na(v)
  cmat <- cor(mergedMEs[keep, , drop=FALSE], v[keep], use="p")
  pmat <- corPvalueStudent(cmat, sum(keep))
  qmat <- matrix(p.adjust(as.vector(pmat), method = "BH"),
                 nrow = nrow(pmat), ncol = ncol(pmat), dimnames = dimnames(pmat))
  txt <- paste(signif(cmat,2), "\nq=", signif(qmat,1), sep = "")
  dim(txt) <- dim(cmat)
  svg(file = fp_traits(paste0("ME_trait_heatmap_", nm, ".svg")),
      width = 3, height = max(4, ncol(mergedMEs) * 0.25 + 2))
  par(mar = c(6, 8.5, 3, 1))
  labeledHeatmap(Matrix = cmat, xLabels = nm,
                 yLabels = colnames(mergedMEs),
                 ySymbols = colnames(mergedMEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = txt,
                 setStdMargins = FALSE,
                  cex.text = 0.8, zlim = c(-1,1),
                  main = "Module-trait relationships")
  dev.off()
}

# --------------------------
# kME, GS, hubs
# --------------------------
modNames <- substring(colnames(mergedMEs), 3)
geneModuleMembership <- as.data.frame(cor(expression.data,
                                          mergedMEs, use = "p",
                                          method = "pearson"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste0("MM", modNames)
names(MMPvalue) <- paste0("p.MM", modNames)

condition_factor <- factor(sample_info$ExpGroup, levels = c("con", "res", "sus"))
gene_condition_p <- vapply(expression.data, function(x) {
  fit <- stats::lm(x ~ condition_factor)
  stats::anova(fit)[["Pr(>F)"]][1]
}, numeric(1))
gene_condition_fdr <- p.adjust(gene_condition_p, method = "BH")
geneTraitSignificance <- data.frame(
  GS.ExpGroup = -log10(pmax(gene_condition_fdr, .Machine$double.xmin)),
  row.names = colnames(expression.data)
)
GSPvalue <- data.frame(
  p.GS.ExpGroup = gene_condition_p,
  FDR.GS.ExpGroup = gene_condition_fdr,
  row.names = colnames(expression.data)
)
names(geneTraitSignificance) <- "GS.ExpGroup"

modules_of_interest <- unique(mergedColors)
for (module in modules_of_interest) {
  moduleGenes <- mergedColors == module
  if (!any(moduleGenes)) next
  mmcol <- paste0("MM", module); if (!mmcol %in% colnames(geneModuleMembership)) next
  gene_info <- data.frame(
    Gene = colnames(expression.data)[moduleGenes],
    Module = module,
    ModuleMembership = geneModuleMembership[moduleGenes, mmcol],
    GeneSignificance = geneTraitSignificance[moduleGenes, "GS.ExpGroup"],
    GeneSignificanceP = GSPvalue[moduleGenes, "p.GS.ExpGroup"],
    GeneSignificanceFDR = GSPvalue[moduleGenes, "FDR.GS.ExpGroup"]
  )
  write_csv_safe(gene_info, fp_modtab(paste0("genes_in_module_", module, ".csv")))
}

# ==========================
# WGCNA module contract and GO enrichment exports
# ==========================
all_feature_ids <- colnames(expression.data)
is_uniprot_accession <- function(x) {
  grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$|^A0A[0-9A-Z]{7}$", x)
}
feature_to_acc <- strsplit(all_feature_ids, ";", fixed = TRUE)
feature_to_acc <- lapply(feature_to_acc, function(v) {
  v <- toupper(trimws(v))
  v[is_uniprot_accession(v)]
})
names(feature_to_acc) <- all_feature_ids
feature_primary_acc <- vapply(feature_to_acc, function(x) if (length(x)) x[1] else NA_character_, character(1))

acc_universe <- unique(stats::na.omit(unlist(feature_to_acc, use.names = FALSE)))
if (length(acc_universe)) {
  map_df <- tryCatch(
    suppressWarnings(clusterProfiler::bitr(
      acc_universe,
      fromType = "UNIPROT",
      toType = c("ENTREZID", "SYMBOL"),
      OrgDb = org.Mm.eg.db
    )),
    error = function(e) NULL
  )
  map_df <- if (is.null(map_df) || !nrow(map_df)) {
    tibble::tibble(UNIPROT = character(), ENTREZID = character(), SYMBOL = character())
  } else {
    tibble::as_tibble(map_df) %>%
      dplyr::mutate(UNIPROT = toupper(.data$UNIPROT)) %>%
      dplyr::distinct(.data$UNIPROT, .keep_all = TRUE)
  }
} else {
  map_df <- tibble::tibble(UNIPROT = character(), ENTREZID = character(), SYMBOL = character())
}
acc_to_entrez <- stats::setNames(as.character(map_df$ENTREZID), map_df$UNIPROT)
acc_to_symbol <- stats::setNames(as.character(map_df$SYMBOL), map_df$UNIPROT)
universe_entrez <- unique(stats::na.omit(unname(acc_to_entrez[acc_universe])))
universe_entrez <- as.character(universe_entrez)

compact_term <- function(term) {
  term <- gsub("\\b(process|regulation|pathway|of|the|cellular)\\b", "", term, ignore.case = TRUE)
  term <- gsub("[[:punct:]]+", " ", term)
  term <- stringr::str_squish(term)
  stringr::str_to_title(term)
}
ratio_to_numeric <- function(x) {
  xy <- strsplit(as.character(x), "/", fixed = TRUE)[[1]]
  if (length(xy) != 2) return(NA_real_)
  suppressWarnings(as.numeric(xy[1]) / as.numeric(xy[2]))
}

ME_names_stable <- colnames(mergedMEs)
ME_colors <- sub("^ME", "", ME_names_stable)
color_to_MEcol <- stats::setNames(ME_names_stable, ME_colors)
module_colors <- sort(unique(as.character(mergedColors)))
module_ids <- stats::setNames(paste0("WGCNA_", module_colors), module_colors)

feature_module_tbl <- tibble::tibble(
  ProteinID = all_feature_ids,
  UniProt = unname(feature_primary_acc[all_feature_ids]),
  EntrezID = unname(acc_to_entrez[feature_primary_acc[all_feature_ids]]),
  GeneSymbol = unname(acc_to_symbol[feature_primary_acc[all_feature_ids]]),
  ModuleColor = as.character(mergedColors[all_feature_ids]),
  ModuleID = unname(module_ids[as.character(mergedColors[all_feature_ids])])
) %>%
  dplyr::mutate(
    EntrezID = as.character(.data$EntrezID),
    GeneSymbol = as.character(.data$GeneSymbol)
  )

kME_long <- purrr::map_dfr(module_colors, function(module_color) {
  mm_col <- paste0("MM", module_color)
  module_features <- feature_module_tbl$ProteinID[feature_module_tbl$ModuleColor == module_color]
  if (!mm_col %in% colnames(geneModuleMembership)) return(NULL)
  tibble::tibble(
    ProteinID = module_features,
    kME = as.numeric(geneModuleMembership[module_features, mm_col, drop = TRUE])
  )
})

make_module_sets <- function(module_color) {
  mod_tbl <- feature_module_tbl %>%
    dplyr::filter(.data$ModuleColor == module_color) %>%
    dplyr::left_join(kME_long, by = "ProteinID") %>%
    dplyr::mutate(abs_kME = abs(.data$kME))
  core_tbl <- mod_tbl %>% dplyr::filter(is.finite(.data$abs_kME), .data$abs_kME >= 0.6)
  hub_tbl <- mod_tbl %>% dplyr::arrange(dplyr::desc(.data$abs_kME)) %>% dplyr::slice_head(n = 25)
  list(all = mod_tbl, core_kME_0.6 = core_tbl, top_hub_25 = hub_tbl)
}

enrich_module_set <- function(module_color, set_name, set_tbl, ontology,
                              min_mapped_n = 5, min_universe_n = 100) {
  mapped_genes <- unique(stats::na.omit(as.character(set_tbl$EntrezID)))
  qc_base <- tibble::tibble(
    ModuleColor = module_color,
    ModuleID = unname(module_ids[module_color]),
    ModuleProteinSetType = set_name,
    Ontology = ontology,
    ModuleSize = nrow(set_tbl),
    MappedModuleSize = length(mapped_genes),
    UniverseSize = length(universe_entrez)
  )
  if (length(universe_entrez) < min_universe_n) {
    return(list(result = NULL, qc = dplyr::mutate(qc_base, status = "skipped", reason = "universe_too_small")))
  }
  if (length(mapped_genes) < min_mapped_n) {
    return(list(result = NULL, qc = dplyr::mutate(qc_base, status = "skipped", reason = "mapped_module_too_small")))
  }
  ego <- tryCatch(
    suppressWarnings(clusterProfiler::enrichGO(
      gene = mapped_genes,
      universe = universe_entrez,
      OrgDb = org.Mm.eg.db,
      keyType = "ENTREZID",
      ont = ontology,
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      readable = FALSE
    )),
    error = function(e) e
  )
  if (inherits(ego, "error")) {
    return(list(result = NULL, qc = dplyr::mutate(qc_base, status = "failed", reason = conditionMessage(ego))))
  }
  df <- as.data.frame(ego)
  if (!nrow(df)) {
    return(list(result = NULL, qc = dplyr::mutate(qc_base, status = "skipped", reason = "no_enriched_terms")))
  }
  df <- tibble::as_tibble(df) %>%
    dplyr::mutate(
      ModuleColor = module_color,
      ModuleID = unname(module_ids[module_color]),
      ModuleProteinSetType = set_name,
      Ontology = ontology,
      ModuleSize = nrow(set_tbl),
      MappedModuleSize = length(mapped_genes),
      UniverseSize = length(universe_entrez),
      .before = 1
    ) %>%
    dplyr::select(
      ModuleColor, ModuleID, ModuleProteinSetType, Ontology,
      ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue,
      geneID, Count, ModuleSize, MappedModuleSize, UniverseSize
    )
  list(result = df, qc = dplyr::mutate(qc_base, status = "ok", reason = NA_character_))
}

go_runs <- list()
go_qc <- list()
for (module_color in module_colors) {
  sets <- make_module_sets(module_color)
  for (set_name in names(sets)) {
    for (ontology in c("BP", "MF", "CC")) {
      run <- enrich_module_set(module_color, set_name, sets[[set_name]], ontology)
      go_runs[[length(go_runs) + 1]] <- run$result
      go_qc[[length(go_qc) + 1]] <- run$qc
    }
  }
}
GO_enrichment_long <- dplyr::bind_rows(go_runs)
GO_enrichment_QC <- dplyr::bind_rows(go_qc)
if (!nrow(GO_enrichment_long)) {
  GO_enrichment_long <- tibble::tibble(
    ModuleColor = character(), ModuleID = character(), ModuleProteinSetType = character(),
    Ontology = character(), ID = character(), Description = character(), GeneRatio = character(),
    BgRatio = character(), pvalue = numeric(), p.adjust = numeric(), qvalue = numeric(),
    geneID = character(), Count = integer(), ModuleSize = integer(), MappedModuleSize = integer(),
    UniverseSize = integer()
  )
}
write_csv_safe(GO_enrichment_long, fp_modtab("WGCNA_module_GO_enrichment_long.csv"))
write_csv_safe(GO_enrichment_QC, fp_modtab("WGCNA_module_GO_enrichment_QC.csv"))

best_go_labels <- GO_enrichment_long %>%
  dplyr::filter(.data$ModuleProteinSetType == "all") %>%
  dplyr::group_by(.data$ModuleColor, .data$Ontology) %>%
  dplyr::arrange(.data$p.adjust, .data$qvalue, .by_group = TRUE) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    ModuleColor,
    Ontology,
    Label = vapply(.data$Description, compact_term, character(1)),
    BestTerm = .data$Description,
    BestTermPAdjust = .data$p.adjust,
    BestTermQvalue = .data$qvalue,
    BestTermGeneRatio = .data$GeneRatio
  )
best_go_wide <- best_go_labels %>%
  dplyr::select(.data$ModuleColor, .data$Ontology, .data$Label) %>%
  tidyr::pivot_wider(names_from = .data$Ontology, values_from = .data$Label, names_prefix = "ModuleLabel_GO_")
best_term_wide <- best_go_labels %>%
  dplyr::select(.data$ModuleColor, .data$Ontology, .data$BestTerm) %>%
  tidyr::pivot_wider(names_from = .data$Ontology, values_from = .data$BestTerm, names_prefix = "best_GO_")
best_term_padj_wide <- best_go_labels %>%
  dplyr::select(.data$ModuleColor, .data$Ontology, .data$BestTermPAdjust) %>%
  tidyr::pivot_wider(names_from = .data$Ontology, values_from = .data$BestTermPAdjust, names_prefix = "best_GO_padj_")
best_term_qvalue_wide <- best_go_labels %>%
  dplyr::select(.data$ModuleColor, .data$Ontology, .data$BestTermQvalue) %>%
  tidyr::pivot_wider(names_from = .data$Ontology, values_from = .data$BestTermQvalue, names_prefix = "best_GO_qvalue_")
best_term_gene_ratio_wide <- best_go_labels %>%
  dplyr::select(.data$ModuleColor, .data$Ontology, .data$BestTermGeneRatio) %>%
  tidyr::pivot_wider(names_from = .data$Ontology, values_from = .data$BestTermGeneRatio, names_prefix = "best_GO_gene_ratio_")

module_label_table <- tibble::tibble(
  ModuleColor = module_colors,
  ModuleID = unname(module_ids[module_colors]),
  ModuleLabel_Manual = NA_character_
) %>%
  dplyr::left_join(best_go_wide, by = "ModuleColor") %>%
  dplyr::left_join(best_term_wide, by = "ModuleColor") %>%
  dplyr::left_join(best_term_padj_wide, by = "ModuleColor") %>%
  dplyr::left_join(best_term_qvalue_wide, by = "ModuleColor") %>%
  dplyr::left_join(best_term_gene_ratio_wide, by = "ModuleColor") %>%
  ensure_module_label_schema() %>%
  dplyr::mutate(
    ModuleLabel_GO_BP = dplyr::coalesce(.data$ModuleLabel_GO_BP, paste0("Module ", .data$ModuleColor)),
    ModuleLabel_GO_MF = dplyr::coalesce(.data$ModuleLabel_GO_MF, paste0("Module ", .data$ModuleColor)),
    ModuleLabel_GO_CC = dplyr::coalesce(.data$ModuleLabel_GO_CC, paste0("Module ", .data$ModuleColor)),
    ModuleLabel_Final = dplyr::coalesce(.data$ModuleLabel_Manual, .data$ModuleLabel_GO_BP, paste0("Module ", .data$ModuleColor)),
    ModuleLabel_Source = dplyr::case_when(
      !is.na(.data$ModuleLabel_Manual) & nzchar(.data$ModuleLabel_Manual) ~ "manual",
      !is.na(.data$best_GO_BP) & nzchar(.data$best_GO_BP) ~ "GO_BP_ORA_all_module",
      TRUE ~ "module_color"
    ),
    primary_label = .data$ModuleLabel_GO_BP,
    alternative_label_MF = .data$ModuleLabel_GO_MF,
    alternative_label_CC = .data$ModuleLabel_GO_CC,
    label_padj_BP = .data$best_GO_padj_BP,
    label_padj_MF = .data$best_GO_padj_MF,
    label_padj_CC = .data$best_GO_padj_CC,
    label_gene_ratio_BP = .data$best_GO_gene_ratio_BP,
    label_gene_ratio_MF = .data$best_GO_gene_ratio_MF,
    label_gene_ratio_CC = .data$best_GO_gene_ratio_CC,
    label_source = .data$ModuleLabel_Source,
    manual_label = .data$ModuleLabel_Manual,
    final_label = .data$ModuleLabel_Final
  )
module_name_map <- stats::setNames(module_label_table$ModuleLabel_Final, module_label_table$ModuleColor)

top_hub_flags <- feature_module_tbl %>%
  dplyr::left_join(kME_long, by = "ProteinID") %>%
  dplyr::mutate(abs_kME = abs(.data$kME)) %>%
  dplyr::group_by(.data$ModuleColor) %>%
  dplyr::arrange(dplyr::desc(.data$abs_kME), .by_group = TRUE) %>%
  dplyr::mutate(is_top_hub_25 = dplyr::row_number() <= pmin(25L, dplyr::n())) %>%
  dplyr::ungroup() %>%
  dplyr::select(.data$ProteinID, .data$is_top_hub_25)

WGCNA_modules_long <- feature_module_tbl %>%
  dplyr::left_join(module_label_table, by = c("ModuleColor", "ModuleID")) %>%
  dplyr::left_join(kME_long, by = "ProteinID") %>%
  dplyr::left_join(top_hub_flags, by = "ProteinID") %>%
  dplyr::mutate(
    ModuleSet = "WGCNA",
    abs_kME = abs(.data$kME),
    GeneSignificanceP = GSPvalue[.data$ProteinID, "p.GS.ExpGroup"],
    GeneSignificanceFDR = GSPvalue[.data$ProteinID, "FDR.GS.ExpGroup"],
    is_core_kME_0.6 = is.finite(.data$abs_kME) & .data$abs_kME >= 0.6,
    Source = "01_WGCNA.r"
  ) %>%
  dplyr::select(
    ModuleSet, ModuleID, ModuleColor,
    ModuleLabel_Final, ModuleLabel_Source, ModuleLabel_GO_BP, ModuleLabel_GO_MF, ModuleLabel_GO_CC, ModuleLabel_Manual,
    primary_label, alternative_label_MF, alternative_label_CC,
    label_padj_BP, label_padj_MF, label_padj_CC,
    label_gene_ratio_BP, label_gene_ratio_MF, label_gene_ratio_CC,
    label_source, manual_label, final_label,
    best_GO_BP, best_GO_MF, best_GO_CC, best_GO_padj_BP, best_GO_padj_MF, best_GO_padj_CC,
    ProteinID, UniProt, EntrezID, GeneSymbol, kME, abs_kME,
    GeneSignificanceP, GeneSignificanceFDR, is_core_kME_0.6, is_top_hub_25, Source
  )
write_csv_safe(WGCNA_modules_long, fp_modtab("WGCNA_modules_long.csv"))
writexl::write_xlsx(list(WGCNA_modules_long = WGCNA_modules_long), fp_modtab("WGCNA_modules_long.xlsx"))
write_csv_safe(module_label_table, fp_modtab("module_name_map.csv"))
write_tsv_safe(module_label_table, fp_modtab("module_name_map.tsv"))

make_wgcna_module_summary <- function(preservation_source = NULL) {
  pres_wide <- tibble::tibble(ModuleColor = module_colors)
  if (!is.null(preservation_source) && nrow(preservation_source)) {
    pres_wide <- preservation_source %>%
      dplyr::select(.data$module, .data$test_set, dplyr::contains("Zsummary")) %>%
      dplyr::rename(ModuleColor = .data$module) %>%
      tidyr::pivot_wider(
        names_from = .data$test_set,
        values_from = dplyr::contains("Zsummary"),
        names_glue = "preservation_{test_set}_{.value}"
      )
  }
  WGCNA_modules_long %>%
    dplyr::group_by(.data$ModuleID, .data$ModuleColor) %>%
    dplyr::summarise(
      n_features = dplyr::n(),
      n_mapped_entrez = dplyr::n_distinct(.data$EntrezID[!is.na(.data$EntrezID) & nzchar(.data$EntrezID)]),
      median_abs_kME = stats::median(.data$abs_kME, na.rm = TRUE),
      mean_abs_kME = mean(.data$abs_kME, na.rm = TRUE),
      top_hub_proteins = paste(utils::head(.data$ProteinID[order(.data$abs_kME, decreasing = TRUE)], 25), collapse = ";"),
      .groups = "drop"
    ) %>%
    dplyr::mutate(mapping_rate = .data$n_mapped_entrez / .data$n_features, .after = .data$n_mapped_entrez) %>%
    dplyr::left_join(best_term_wide, by = "ModuleColor") %>%
    dplyr::left_join(pres_wide, by = "ModuleColor")
}
WGCNA_module_summary <- make_wgcna_module_summary()
write_csv_safe(WGCNA_module_summary, fp_modtab("WGCNA_module_summary.csv"))

saveRDS(list(
  module_name_map = module_name_map,
  module_label_table = module_label_table,
  color_to_MEcol = color_to_MEcol,
  ME_names_stable = ME_names_stable
), file = fp_state("module_name_map.rds"))

# --------------------------
# Preservation across conditions
# --------------------------
idx_con <- which(sample_info$ExpGroup == "con")
idx_res <- which(sample_info$ExpGroup == "res")
idx_sus <- which(sample_info$ExpGroup == "sus")

multiExpr <- list(
  ALL = list(data = expression.data),
  CON = list(data = expression.data[idx_con, , drop = FALSE]),
  RES = list(data = expression.data[idx_res, , drop = FALSE]),
  SUS = list(data = expression.data[idx_sus, , drop = FALSE])
)

good_cols <- lapply(multiExpr, function(e) {
  gsg <- goodSamplesGenes(e$data, verbose = 3)
  which(gsg$goodGenes)
})
common_idx <- Reduce(intersect, good_cols)
common_genes <- colnames(multiExpr[[1]]$data)[common_idx]

has_na <- sapply(common_genes, function(g) any(vapply(multiExpr, function(e) any(is.na(e$data[, g])), logical(1))))
common_genes <- common_genes[!has_na]
stopifnot(length(common_genes) > 0)

multi_expr_clean <- lapply(multiExpr, function(e) {
  dat <- e$data[, common_genes, drop = FALSE]
  dat <- as.data.frame(dat)
  dat[] <- lapply(dat, function(col) as.numeric(col))
  keep_samp <- apply(dat, 1, function(r) sd(r, na.rm = TRUE) > 0)
  dat <- dat[keep_samp, , drop = FALSE]
  list(data = as.matrix(dat))
})
names(multi_expr_clean) <- names(multiExpr)

if (any(vapply(multi_expr_clean, function(e) nrow(e$data) < 2, logical(1)))) {
  stop("One or more sets have < 2 samples after cleaning; reduce filtering or combine sets.")
}

ref_colors <- mergedColors[match(common_genes, colnames(expression.data))]
stopifnot(length(ref_colors) == length(common_genes))
if (any(is.na(ref_colors) | !nzchar(ref_colors))) {
  stop("Missing module colors for one or more common genes before module preservation.")
}
ref_colors <- as.character(ref_colors)
names(ref_colors) <- common_genes

make_reserved_module_name <- function(base, existing) {
  out <- base
  i <- 1
  while (out %in% existing) {
    out <- paste0(base, "_", i)
    i <- i + 1
  }
  out
}

preservation_grey_name <- make_reserved_module_name("grey", ref_colors)
preservation_gold_name <- make_reserved_module_name("gold", c(ref_colors, preservation_grey_name))
preservation_min_module_size <- 3
ref_module_sizes <- table(ref_colors)
small_preservation_modules <- names(ref_module_sizes)[ref_module_sizes < preservation_min_module_size]
if (length(small_preservation_modules)) {
  ref_colors[ref_colors %in% small_preservation_modules] <- preservation_grey_name
}
ref_module_sizes <- table(ref_colors)
preserved_module_sizes <- ref_module_sizes[names(ref_module_sizes) != preservation_grey_name]
if (!length(preserved_module_sizes)) {
  stop("No modules have at least ", preservation_min_module_size,
       " shared genes for module preservation after filtering.")
}

preservation_input_summary <- dplyr::bind_rows(lapply(names(multi_expr_clean), function(set_name) {
  tibble::tibble(
    set = set_name,
    n_samples = nrow(multi_expr_clean[[set_name]]$data),
    n_genes = ncol(multi_expr_clean[[set_name]]$data),
    n_modules = length(setdiff(unique(ref_colors), preservation_grey_name)),
    smallest_module_size = min(preserved_module_sizes),
    largest_module_size = max(preserved_module_sizes),
    grey_name = preservation_grey_name,
    gold_name = preservation_gold_name,
    collapsed_small_modules = paste(small_preservation_modules, collapse = ";")
  )
}))
write_csv_safe(preservation_input_summary, fp_log("module_preservation_input_summary.csv"))

multi_color <- stats::setNames(
  rep(list(ref_colors), length(multi_expr_clean)),
  names(multi_expr_clean)
)

set.seed(12345)
mp <- modulePreservation(
  multi_expr_clean,
  multi_color,
  referenceNetworks = 1,
  nPermutations = module_preservation_permutations,
  networkType = "signed",
  corFnc = "bicor",
  corOptions = "use = 'p', maxPOutliers = 0.05",
  randomSeed = 12345,
  greyName = preservation_grey_name,
  goldName = preservation_gold_name,
  parallelCalculation = FALSE,
  verbose = 3
)

Z_CON <- mp$preservation$Z$ref.ALL$inColumnsAlsoPresentIn.CON
Z_RES <- mp$preservation$Z$ref.ALL$inColumnsAlsoPresentIn.RES
Z_SUS <- mp$preservation$Z$ref.ALL$inColumnsAlsoPresentIn.SUS

write_csv_safe(tibble::rownames_to_column(as.data.frame(Z_CON), "module"),
               fp_prestab("module_preservation_CON_vs_ALL_Zsummary.csv"))
write_csv_safe(tibble::rownames_to_column(as.data.frame(Z_RES), "module"),
               fp_prestab("module_preservation_RES_vs_ALL_Zsummary.csv"))
write_csv_safe(tibble::rownames_to_column(as.data.frame(Z_SUS), "module"),
               fp_prestab("module_preservation_SUS_vs_ALL_Zsummary.csv"))
module_preservation_long <- dplyr::bind_rows(
  tibble::rownames_to_column(as.data.frame(Z_CON), "module") |> dplyr::mutate(test_set = "CON"),
  tibble::rownames_to_column(as.data.frame(Z_RES), "module") |> dplyr::mutate(test_set = "RES"),
  tibble::rownames_to_column(as.data.frame(Z_SUS), "module") |> dplyr::mutate(test_set = "SUS")
)
write_csv_safe(module_preservation_long, fp_source("module_preservation.csv"))

WGCNA_module_summary <- make_wgcna_module_summary(module_preservation_long)
write_csv_safe(WGCNA_module_summary, fp_modtab("WGCNA_module_summary.csv"))
saveRDS(list(
  expression.data = expression.data,
  sample_info = sample_info,
  mergedColors = mergedColors,
  mergedMEs = mergedMEs,
  kME = WGCNA_modules_long %>%
    dplyr::select(.data$ModuleID, .data$ModuleColor, .data$ProteinID, .data$UniProt, .data$EntrezID, .data$GeneSymbol, .data$kME, .data$abs_kME),
  WGCNA_modules_long = WGCNA_modules_long,
  module_summary = WGCNA_module_summary,
  GO_enrichment = GO_enrichment_long,
  GO_enrichment_QC = GO_enrichment_QC,
  module_name_map = module_name_map,
  module_label_table = module_label_table,
  color_to_MEcol = color_to_MEcol,
  ME_names_stable = ME_names_stable,
  module_preservation = module_preservation_long,
  geneTree = geneTree,
  softPower = softPower,
  parameters = list(
    network_type = "signed",
    correlation_function = "bicor",
    bicor_maxPOutliers = 0.05,
    soft_threshold_rsquared = soft_threshold_rsquared,
    selected_soft_power = softPower,
    min_module_size = min_module_size,
    deep_split = deep_split,
    merge_cut_height = merge_cut_height,
    module_preservation_permutations = module_preservation_permutations,
    dataset_profile = dataset_profile_resolved
  )
), file = fp_state("wgcna_final_model_state.rds"))
}

# --------------------------
# Condition and active spatial trait panels
# --------------------------
combo_vars <- c("condition", active_spatial_vars)
combo_df <- sample_info[, combo_vars, drop = FALSE] |>
  tibble::rownames_to_column("Sample") |>
  dplyr::mutate(dplyr::across(dplyr::all_of(combo_vars), as.character))
combo <- factor(do.call(paste, c(combo_df[combo_vars], sep = "_")))

# One-hot and correlations
X_combo <- model.matrix(~ 0 + combo)
colnames(X_combo) <- levels(combo)  # no "comb" prefix at all
if (ncol(X_combo) == 0) stop("No combined strata present in X_combo (ncol == 0). Check combo label construction.")

MEcorr_combo <- cor(mergedMEs, X_combo, use = "p")
MEp_combo    <- corPvalueStudent(MEcorr_combo, nrow(expression.data))

df_combo <- reshape2::melt(MEcorr_combo, varnames = c("module","comb"), value.name = "r")
p_combo  <- reshape2::melt(MEp_combo,    varnames = c("module","comb"), value.name = "p")
df_combo$p <- p_combo$p
df_combo$fdr <- p.adjust(df_combo$p, method = "BH")
df_combo$sig <- sig_dot(df_combo$fdr)
df_combo$module <- factor(df_combo$module, levels = rownames(MEcorr_combo))

split_combo_labels <- function(labels, vars) {
  parts <- strsplit(labels, "_", fixed = TRUE)
  do.call(rbind, lapply(parts, function(p) {
    p <- as.character(p)
    if (length(p) < length(vars)) p <- c(p, rep("missing", length(vars) - length(p)))
    p[seq_along(vars)]
  }))
}
lab_parts <- split_combo_labels(as.character(df_combo$comb), combo_vars)
colnames(lab_parts) <- combo_vars
df_combo2 <- cbind(df_combo[, c("module","comb","r","p","fdr","sig")], as.data.frame(lab_parts, stringsAsFactors = FALSE))
df_combo2$spatial_trait <- do.call(paste, c(df_combo2[active_spatial_vars], sep = "_"))
df_combo2$trait <- do.call(paste, c(df_combo2[c("condition", active_spatial_vars)], sep = "_"))
write_csv_safe(df_combo2, fp_source("ME_condition_spatial_strata.csv"))

present_conds <- intersect(c("con","res","sus"), unique(df_combo2$condition))
by_cond <- split(df_combo2, df_combo2$condition)

panel_plot <- function(dfi, panel_title = "") {
  if (!nrow(dfi)) return(NULL)
  ord <- do.call(order, dfi[active_spatial_vars])
  dfi$spatial_trait <- factor(dfi$spatial_trait, levels = unique(dfi$spatial_trait[ord]))
  ggplot(dfi, aes(x = spatial_trait, y = module, fill = r)) +
    geom_tile(color = "white", linewidth = 0.15) +
    geom_text(aes(label = sig), size = 1.8, color = "black", na.rm = TRUE) +
    scale_fill_gradient2(limits = c(-1, 1), oob = scales::squish,
                         low = figure_diverging["low"], mid = figure_diverging["mid"], high = figure_diverging["high"]) +
    labs(title = panel_title, x = paste(active_spatial_vars, collapse = " / "), y = NULL, fill = "Pearson r") +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none")
}

plots <- lapply(present_conds, function(cn) panel_plot(by_cond[[cn]], toupper(cn)))
plots <- Filter(function(p) !is.null(p) && inherits(p, "ggplot"), plots)
if (length(plots) == 0) stop("No non-empty condition panels; check combo labels and present conditions.")

# Combine panels; if single, save directly
if (length(plots) == 1) {
  g_combined <- plots[[1]]
  save_plot_publication(g_combined, fp_traits("panel_ME_vs_condition_spatial_strata.svg"),
                   width = figure_double_col, height = 3.9)
} else {
  g_combined <- patchwork::wrap_plots(plots, nrow = 1)
  legend_df <- data.frame(x = 1:3, y = 1:3, r = c(-1, 0, 1))
  p_legend <- ggplot(legend_df, aes(x, y, fill = r)) +
    geom_tile() +
    scale_fill_gradient2(limits = c(-1, 1), oob = scales::squish,
                         low = figure_diverging["low"], mid = figure_diverging["mid"], high = figure_diverging["high"]) +
    theme_void() + theme(legend.position = "right") + labs(fill = "Pearson r")
  legend_only <- cowplot::get_legend(p_legend)
  g <- cowplot::plot_grid(g_combined, legend_only, rel_widths = c(1, 0.08))
  save_plot_publication(g, fp_traits("panel_ME_vs_condition_spatial_strata.svg"),
                   width = figure_double_col, height = 3.9)
}

# --------------------------
# Ensure required objects are defined for the pheatmap block
# --------------------------
if (!exists("df_all")) {
    if (exists("df_combo2")) {
        df_all <- df_combo2
        if (!"trait" %in% names(df_all)) {
            df_all$trait <- do.call(paste, c(df_all[c("condition", active_spatial_vars)], sep = "_"))
        }
    } else if (exists("MEcorr_combo") && exists("MEp_combo")) {
        df_all <- reshape2::melt(MEcorr_combo, varnames = c("module", "comb"), value.name = "r")
        p_tmp  <- reshape2::melt(MEp_combo,    varnames = c("module", "comb"), value.name = "p")
        df_all$p <- p_tmp$p
        parts <- split_combo_labels(as.character(df_all$comb), combo_vars)
        colnames(parts) <- combo_vars
        df_all <- cbind(df_all, as.data.frame(parts, stringsAsFactors = FALSE))
        df_all$trait <- do.call(paste, c(df_all[c("condition", active_spatial_vars)], sep = "_"))
    } else {
        stop("df_all is missing and cannot be reconstructed: provide df_combo2 or MEcorr_combo/MEp_combo.")
    }
}

# module_levels (ordered modules)
if (!exists("module_levels")) {
    if (exists("MEcorr_combo")) {
        module_levels <- rownames(MEcorr_combo)
    } else {
        module_levels <- unique(df_all$module)
    }
}

# traits_in_order (ordered trait columns)
if (!exists("traits_in_order")) {
    if (exists("traits_in_order") && length(traits_in_order) > 0) {
        # already provided
    } else {
        df_all$trait <- do.call(paste, c(df_all[c("condition", active_spatial_vars)], sep = "_"))
        ord <- do.call(order, df_all[c("condition", active_spatial_vars)])
        traits_in_order <- unique(df_all$trait[ord])
    }
}

# strip_df: module -> color mapping
if (!exists("strip_df")) {
    if (exists("mergedColors")) {
        # mergedColors is gene-level mapping; fall back to simple default
        strip_df <- data.frame(module = module_levels, mod_col = rep("grey80", length(module_levels)), stringsAsFactors = FALSE)
    } else {
        strip_df <- data.frame(module = module_levels, mod_col = rep("grey80", length(module_levels)), stringsAsFactors = FALSE)
    }
}

# 1) Wide matrices and order
r_mat  <- reshape2::acast(df_all, module ~ trait, value.var = "r")
p_mat  <- reshape2::acast(df_all, module ~ trait, value.var = "p")

# FDR matrix
fdr_vec <- p.adjust(as.vector(p_mat), method = "BH")
fdr_mat <- matrix(fdr_vec, nrow = nrow(p_mat), ncol = ncol(p_mat), dimnames = dimnames(p_mat))

stopifnot(all(module_levels %in% rownames(r_mat)), all(traits_in_order %in% colnames(r_mat)))
r_mat   <- r_mat[module_levels, traits_in_order, drop = FALSE]
p_mat   <- p_mat[module_levels, traits_in_order, drop = FALSE]
fdr_mat <- fdr_mat[module_levels, traits_in_order, drop = FALSE]

# 2) Column clustering and separators
if (ncol(r_mat) > 1) {
  trait_cor <- cor(r_mat, use = "pairwise.complete.obs")
  trait_cor[!is.finite(trait_cor)] <- 0
  diag(trait_cor) <- 1
  hc_cols <- hclust(as.dist(1 - trait_cor), method = "average")
  k_clusters <- min(6, ncol(r_mat))
  col_grp <- cutree(hc_cols, k = k_clusters)
  grp_ord <- col_grp[hc_cols$order]
  gap_pos <- which(grp_ord[-1] != head(grp_ord, -1))
  xlines <- gap_pos + 0.5
} else {
  hc_cols <- FALSE
  k_clusters <- 1
  xlines <- numeric(0)
}

# 3) Module color mapping for row strip
if (exists("mergedMEs")) {
  me_names <- colnames(mergedMEs)                      # e.g., "MEblue"
  row_me <- rownames(r_mat)
  row_mod_colors <- sub("^ME", "", row_me)
  row_mod_colors[is.na(row_mod_colors) | !nzchar(row_mod_colors)] <- sub("^ME", "", row_me[is.na(row_mod_colors) | !nzchar(row_mod_colors)])
  if (!exists("colorSeq")) {
    colorSeq <- c(
      "turquoise"="turquoise","blue"="blue","brown"="brown","yellow"="yellow",
      "green"="green","red"="red","black"="black","pink"="pink","magenta"="magenta",
      "purple"="purple","greenyellow"="greenyellow","tan"="tan","salmon"="salmon",
      "cyan"="cyan","midnightblue"="midnightblue","lightcyan"="lightcyan",
      "greenyellow2"="#ADFF2F","royalblue"="royalblue","darkred"="darkred",
      "skyblue"="skyblue","orange"="orange","grey"="grey"
    )
  }
  color_lookup <- colorSeq
  missing_cols <- setdiff(unique(row_mod_colors), names(color_lookup))
  if (length(missing_cols)) {
    add <- setNames(missing_cols, missing_cols)
    color_lookup <- c(color_lookup, add)
  }
  row_colors_vec <- unname(color_lookup[row_mod_colors])
  row_colors_vec[is.na(row_colors_vec)] <- "grey80"
  ann_row <- data.frame(ModuleColor = row_colors_vec, row.names = rownames(r_mat))
  ann_colors <- list(ModuleColor = setNames(unique(ann_row$ModuleColor), unique(ann_row$ModuleColor)))
} else {
  ann_row <- data.frame(ModuleColor = rep("grey80", nrow(r_mat)), row.names = rownames(r_mat))
  ann_colors <- list(ModuleColor = c("grey80" = "grey80"))
}

# Column condition annotation
col_condition <- sapply(strsplit(colnames(r_mat), "_", fixed = TRUE), `[`, 1)
ann_col <- data.frame(Condition = col_condition, row.names = colnames(r_mat))
ann_colors$Condition <- c(con = "#E6E6E6", res = "#D9ECF7", sus = "#F6D9CC")

# 4) Palette and breaks
r_lim <- 0.8
bk <- seq(-r_lim, r_lim, length.out = 201)
pal <- colorRampPalette(c(figure_diverging["low"], figure_diverging["mid"], figure_diverging["high"]))(200)

# 5) Significance dots
display_mat <- matrix("", nrow = nrow(r_mat), ncol = ncol(r_mat), dimnames = dimnames(r_mat))
display_mat[fdr_mat < 0.01] <- "\u2022"

# 6) Render pheatmap once and save
ph <- pheatmap::pheatmap(
  mat = r_mat,
  cluster_rows = FALSE,
  cluster_cols = hc_cols,
  treeheight_row = 0,
  treeheight_col = 50,
  gaps_col = NULL,
  annotation_row = ann_row,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  color = pal, breaks = bk, legend = TRUE,
  border_color = NA,
  cellwidth = 8.5, cellheight = 8.5,
  angle_col = 45,
  fontsize_col = 6.2, fontsize_row = 6.2,
  display_numbers = display_mat,
  number_color = "grey20",
  number_cex = 0.6,
  show_rownames = TRUE, show_colnames = TRUE,
  silent = TRUE
)

nr <- nrow(r_mat)

draw_trait_pheatmap <- function(ph, xlines = numeric(0), nr = nrow(r_mat), separator_lwd = 2) {
  grid::grid.newpage()
  grid::grid.draw(ph$gtable)
  if (!length(xlines)) return(invisible(NULL))
  try({
    panel_id <- grep("^matrix$", ph$gtable$layout$name)[1]
    if (is.na(panel_id)) stop("matrix panel not found", call. = FALSE)
    grid::seekViewport(ph$gtable$layout$name[panel_id])
    for (xl in xlines) {
      grid::grid.lines(
        x = grid::unit(c(xl, xl), "native"),
        y = grid::unit(c(0, nr), "native"),
        gp = grid::gpar(col = "white", lwd = separator_lwd)
      )
    }
    grid::upViewport(0)
  }, silent = TRUE)
  invisible(NULL)
}

png(fp_traits("ME_trait_pheatmap.png"), width = round(figure_double_col * 300), height = round(4.3 * 300), res = 300)
draw_trait_pheatmap(ph, xlines = xlines, nr = nr, separator_lwd = 3)
dev.off()

pdf(fp_traits("ME_trait_pheatmap.pdf"), width = figure_double_col, height = 4.3, family = figure_font, useDingbats = FALSE)
draw_trait_pheatmap(ph, xlines = xlines, nr = nr, separator_lwd = 2)
dev.off()

svglite::svglite(fp_traits("ME_trait_pheatmap.svg"), width = figure_double_col, height = 4.3)
draw_trait_pheatmap(ph, xlines = xlines, nr = nr, separator_lwd = 2)
dev.off()

block_list <- c(
  setNames(lapply(active_spatial_vars, function(v) grep(paste0("^", v, "_"), colnames(datTraits))), active_spatial_vars),
  list(condition = grep("^condition_", colnames(datTraits)))
)

block_cor_df <- function(block_idx, block_name) {
  if (length(block_idx) == 0) return(NULL)
  r <- cor(mergedMEs, datTraits[, block_idx, drop=FALSE], use="p")
  p <- corPvalueStudent(r, nrow(expression.data))
  df <- reshape2::melt(r, varnames = c("module","trait"), value.name = "r")
  p_long <- reshape2::melt(p, varnames = c("module","trait"), value.name = "p")
  df$p <- p_long$p
  df$fdr <- p.adjust(df$p, method = "BH")
  df$sig <- sig_dot(df$fdr)
  df$block <- block_name
  df$module <- factor(df$module, levels = rownames(r))
  df
}

dfs <- Filter(Negate(is.null),
              mapply(block_cor_df, block_list, names(block_list), SIMPLIFY = FALSE))
write_csv_safe(dplyr::bind_rows(dfs), fp_source("ME_trait_block_correlations.csv"))

panel_plot <- function(dfi, legend = "none") {
  if (nrow(dfi) == 0) return(NULL)
  dfi$trait <- factor(dfi$trait, levels = unique(dfi$trait))

  p_heat <- ggplot(dfi, aes(x = trait, y = module, fill = r)) +
    geom_tile(color = "white", linewidth = 0.15) +
    geom_text(aes(label = sig), size = 1.8, color = "black", na.rm = TRUE) +
    scale_fill_gradient2(limits = c(-1, 1), oob = scales::squish,
                         low = figure_diverging["low"], mid = figure_diverging["mid"], high = figure_diverging["high"]) +
    labs(x = NULL, y = NULL, fill = "Pearson r") +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks = element_blank(),
          legend.position = legend,
          plot.margin = margin(2, 2, 2, 2))

  if (!identical(legend, "none")) return(p_heat)

  modules <- levels(dfi$module)
  mod_cols <- sub("^ME", "", modules)
  strip_df <- data.frame(module = modules, mod_col = mod_cols, stringsAsFactors = FALSE)
  strip_df$module <- factor(strip_df$module, levels = modules)

  p_strip <- ggplot(strip_df, aes(x = 1, y = module, fill = mod_col)) +
    geom_tile() +
    scale_fill_identity() +
    theme_void() +
    theme(plot.margin = margin(2, 0, 2, 2))

  combined <- p_strip + p_heat + plot_layout(widths = c(0.05, 1))
  combined
}

plots <- lapply(dfs, panel_plot)
p_legend_heat <- panel_plot(dfs[[1]], legend = "right")
legend_only <- cowplot::get_legend(p_legend_heat)

combo <- wrap_plots(plots, nrow = 1, guides = "collect") +
  plot_annotation(title = paste("Module-trait relationships:", paste(c(active_spatial_vars, "condition"), collapse = ", ")))

main_trait_panel <- cowplot::plot_grid(combo, legend_only, rel_widths = c(1, 0.08))
save_plot_publication(main_trait_panel, fp_mainfig("panel_module_trait_relationships_spatial.svg"),
                 width = figure_double_col, height = 3.7)

# ==========================================================
# ME by condition plots with significance and exports
# ==========================================================

# Long-form eigengenes + metadata
stopifnot(nrow(mergedMEs) == nrow(sample_info))
ME_long <- mergedMEs %>%
  tibble::rownames_to_column("Sample") %>%
  mutate(
    condition = as.character(sample_info$ExpGroup),
    region = as.character(sample_info$region),
    layer = as.character(sample_info$layer),
    celltype = as.character(sample_info$celltype)
  ) %>%
  pivot_longer(cols = starts_with("ME"),
               names_to = "module", values_to = "ME")

# Ensure order con, res, sus
ME_long$condition <- factor(ME_long$condition, levels = c("con","res","sus"))
ME_long$region <- factor(ME_long$region)
ME_long$layer <- factor(ME_long$layer)
ME_long$celltype <- factor(ME_long$celltype)

supermodule_annotation <- make_supermodule_annotation(module_label_table, colnames(mergedMEs))
supermodule_matched <- supermodule_annotation %>% dplyr::filter(.data$present_in_dataset, .data$manual_annotation)
supermodule_missing <- supermodule_annotation %>% dplyr::filter(!.data$present_in_dataset, .data$manual_annotation)
supermodule_unassigned <- supermodule_annotation %>% dplyr::filter(.data$present_in_dataset, .data$Supermodule == "Unassigned")
message(
  "WGCNA supermodule annotation: matched=", nrow(supermodule_matched),
  "; manual-not-present=", nrow(supermodule_missing),
  "; unassigned-present=", nrow(supermodule_unassigned)
)
if (nrow(supermodule_missing)) {
  message("Manual supermodule modules not present in this dataset: ", paste(supermodule_missing$module_eigengene, collapse = ", "))
}
if (nrow(supermodule_unassigned)) {
  message("Present WGCNA modules without manual supermodule annotation: ", paste(supermodule_unassigned$module_eigengene, collapse = ", "))
}

write_csv_safe(supermodule_annotation, fp_supertab("wgcna_module_supermodule_annotation.csv"))
write_csv_safe(supermodule_annotation, fp_source("wgcna_module_supermodule_annotation.csv"))
writexl::write_xlsx(list(supermodule_annotation = supermodule_annotation), fp_supertab("wgcna_module_supermodule_annotation.xlsx"))

ME_long <- add_supermodule_cols(ME_long, supermodule_annotation, module_col = "module")
if (exists("df_combo2")) {
  df_combo2_supermodule <- add_supermodule_cols(df_combo2, supermodule_annotation, module_col = "module")
  write_csv_safe(df_combo2_supermodule, fp_supertab("wgcna_module_trait_correlations_with_supermodules.csv"))
  write_csv_safe(df_combo2_supermodule, fp_source("wgcna_module_trait_correlations_with_supermodules.csv"))
}

module_axis_labels <- function(modules, include_id = TRUE, width = 34) {
  module_colors <- sub("^ME", "", as.character(modules))
  labels <- unname(module_name_map[module_colors])
  labels[is.na(labels) | !nzchar(labels)] <- paste("Module", module_colors[is.na(labels) | !nzchar(labels)])
  labels <- stringr::str_to_sentence(labels)
  labels <- stringr::str_wrap(labels, width = width)
  if (isTRUE(include_id)) paste0(as.character(modules), "\n", labels) else labels
}

contrast_specs <- tibble::tibble(
  contrast = c("RES - CON", "SUS - CON", "SUS - RES"),
  group_a = c("con", "con", "res"),
  group_b = c("res", "sus", "sus")
)

lm_condition_contrast <- function(df, group_a, group_b, covars) {
  covars <- covars[vapply(df[covars], function(x) length(unique(stats::na.omit(x))) > 1, logical(1))]
  rhs <- paste(c("condition", covars), collapse = " + ")
  model_formula <- stats::as.formula(paste("ME ~", rhs))
  fit <- tryCatch(stats::lm(model_formula, data = df), error = function(e) NULL)
  if (is.null(fit)) {
    return(tibble::tibble(adjusted_delta = NA_real_, p = NA_real_, model_formula = deparse(model_formula)))
  }

  beta <- stats::coef(fit)
  vc <- stats::vcov(fit)
  contrast_vector <- rep(0, length(beta))
  names(contrast_vector) <- names(beta)
  condition_effect <- function(group) {
    out <- rep(0, length(beta))
    names(out) <- names(beta)
    coef_name <- paste0("condition", group)
    if (coef_name %in% names(out)) out[coef_name] <- 1
    out
  }
  contrast_vector <- condition_effect(group_b) - condition_effect(group_a)

  contrast_idx <- which(contrast_vector != 0)
  if (!length(contrast_idx) || any(is.na(beta[contrast_idx]))) {
    return(tibble::tibble(adjusted_delta = NA_real_, p = NA_real_, model_formula = deparse(model_formula)))
  }

  adjusted_delta <- unname(sum(contrast_vector[contrast_idx] * beta[contrast_idx]))
  se <- suppressWarnings(sqrt(as.numeric(t(contrast_vector[contrast_idx]) %*% vc[contrast_idx, contrast_idx, drop = FALSE] %*% contrast_vector[contrast_idx])))
  p <- if (is.finite(se) && se > 0) {
    2 * stats::pt(abs(adjusted_delta / se), df = stats::df.residual(fit), lower.tail = FALSE)
  } else {
    NA_real_
  }
  tibble::tibble(adjusted_delta = adjusted_delta, p = p, model_formula = deparse(model_formula))
}

ME_contrast_stats <- ME_long %>%
  dplyr::group_by(.data$module) %>%
  dplyr::group_map(~{
    dfm <- .x
    module_name <- .y$module
    purrr::map_dfr(seq_len(nrow(contrast_specs)), function(i) {
      a <- contrast_specs$group_a[i]
      b <- contrast_specs$group_b[i]
      adjusted <- lm_condition_contrast(dfm, a, b, active_spatial_vars)
      tibble::tibble(
        module = module_name,
        contrast = contrast_specs$contrast[i],
        group_a = a,
        group_b = b,
        n_a = sum(dfm$condition == a, na.rm = TRUE),
        n_b = sum(dfm$condition == b, na.rm = TRUE),
        mean_a = mean(dfm$ME[dfm$condition == a], na.rm = TRUE),
        mean_b = mean(dfm$ME[dfm$condition == b], na.rm = TRUE),
        delta_mean = mean_b - mean_a
      ) %>%
        dplyr::bind_cols(adjusted)
    })
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(
    fdr = p.adjust(.data$p, method = "BH"),
    sig = sig_dot(.data$fdr),
    contrast = factor(.data$contrast, levels = contrast_specs$contrast)
  )
ME_contrast_stats <- add_supermodule_cols(ME_contrast_stats, supermodule_annotation, module_col = "module")

write_csv_safe(ME_contrast_stats, fp_traittab("ME_by_condition_pairwise_contrasts.csv"))
write_csv_safe(ME_contrast_stats, fp_source("ME_by_condition_pairwise_contrasts.csv"))
write_csv_safe(ME_contrast_stats, fp_supertab("wgcna_group_contrasts_with_supermodules.csv"))

contrast_module_order <- ME_contrast_stats %>%
  dplyr::group_by(.data$module) %>%
  dplyr::summarise(max_abs_delta = max(abs(.data$adjusted_delta), na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(.data$max_abs_delta)) %>%
  dplyr::pull(.data$module)

ME_contrast_stats$module <- factor(ME_contrast_stats$module, levels = rev(contrast_module_order))
contrast_lim <- max(abs(ME_contrast_stats$adjusted_delta), na.rm = TRUE)
if (!is.finite(contrast_lim) || contrast_lim == 0) contrast_lim <- 1

contrast_panel <- ggplot(ME_contrast_stats, aes(x = contrast, y = module, fill = adjusted_delta)) +
  geom_tile(color = "white", linewidth = 0.15) +
  geom_text(aes(label = sig), size = 1.8, color = "black", na.rm = TRUE) +
  scale_fill_gradient2(
    limits = c(-contrast_lim, contrast_lim),
    oob = scales::squish,
    low = figure_diverging["low"],
    mid = figure_diverging["mid"],
    high = figure_diverging["high"]
  ) +
  labs(x = NULL, y = NULL, fill = "Adjusted delta") +
  theme_publication() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks = element_blank()
  )

save_plot_publication(
  contrast_panel,
  fp_traits("panel_ME_by_condition_pairwise_contrasts.svg"),
  width = figure_single_col,
  height = max(3.2, 0.18 * length(contrast_module_order) + 1.1)
)
save_plot_publication(
  contrast_panel,
  fp_mainfig("panel_ME_by_condition_pairwise_contrasts.svg"),
  width = figure_single_col,
  height = max(3.2, 0.18 * length(contrast_module_order) + 1.1)
)

contrast_panel_GO_labels <- contrast_panel +
  scale_y_discrete(labels = function(x) module_axis_labels(x, include_id = TRUE, width = 30))

save_plot_publication(
  contrast_panel_GO_labels,
  fp_traits("panel_ME_by_condition_pairwise_contrasts_GO_labels.svg"),
  width = figure_double_col,
  height = max(3.6, 0.27 * length(contrast_module_order) + 1.1)
)
save_plot_publication(
  contrast_panel_GO_labels,
  fp_mainfig("panel_ME_by_condition_pairwise_contrasts_GO_labels.svg"),
  width = figure_double_col,
  height = max(3.6, 0.27 * length(contrast_module_order) + 1.1)
)

ME_strata <- ME_long %>%
  dplyr::mutate(
    spatial_trait = do.call(paste, c(dplyr::pick(dplyr::all_of(active_spatial_vars)), sep = "_"))
  )

strata_levels <- ME_strata %>%
  dplyr::distinct(dplyr::across(dplyr::all_of(active_spatial_vars)), .data$spatial_trait) %>%
  dplyr::arrange(dplyr::across(dplyr::all_of(active_spatial_vars))) %>%
  dplyr::pull(.data$spatial_trait)

ME_strata_contrast_stats <- ME_strata %>%
  dplyr::group_by(.data$module, .data$spatial_trait) %>%
  dplyr::group_map(~{
    dfm <- .x
    module_name <- .y$module
    spatial_name <- .y$spatial_trait
    purrr::map_dfr(seq_len(nrow(contrast_specs)), function(i) {
      a <- contrast_specs$group_a[i]
      b <- contrast_specs$group_b[i]
      adjusted <- lm_condition_contrast(dfm, a, b, character(0))
      tibble::tibble(
        module = module_name,
        spatial_trait = spatial_name,
        contrast = contrast_specs$contrast[i],
        group_a = a,
        group_b = b,
        n_a = sum(dfm$condition == a, na.rm = TRUE),
        n_b = sum(dfm$condition == b, na.rm = TRUE),
        mean_a = mean(dfm$ME[dfm$condition == a], na.rm = TRUE),
        mean_b = mean(dfm$ME[dfm$condition == b], na.rm = TRUE),
        delta_mean = mean_b - mean_a
      ) %>%
        dplyr::bind_cols(adjusted)
    })
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(
    fdr = p.adjust(.data$p, method = "BH"),
    sig = sig_dot(.data$fdr),
    contrast = factor(.data$contrast, levels = contrast_specs$contrast),
    spatial_trait = factor(.data$spatial_trait, levels = strata_levels),
    module = factor(.data$module, levels = levels(ME_contrast_stats$module))
  )
ME_strata_contrast_stats <- add_supermodule_cols(ME_strata_contrast_stats, supermodule_annotation, module_col = "module")

write_csv_safe(ME_strata_contrast_stats, fp_traittab("ME_by_condition_spatial_strata_pairwise_contrasts.csv"))
write_csv_safe(ME_strata_contrast_stats, fp_source("ME_by_condition_spatial_strata_pairwise_contrasts.csv"))
write_csv_safe(ME_strata_contrast_stats, fp_supertab("wgcna_spatial_strata_group_contrasts_with_supermodules.csv"))

strata_contrast_lim <- max(abs(ME_strata_contrast_stats$adjusted_delta), na.rm = TRUE)
if (!is.finite(strata_contrast_lim) || strata_contrast_lim == 0) strata_contrast_lim <- 1

strata_contrast_plot <- ggplot(ME_strata_contrast_stats, aes(x = spatial_trait, y = module, fill = adjusted_delta)) +
  geom_tile(color = "white", linewidth = 0.15) +
  geom_text(aes(label = sig), size = 1.8, color = "black", na.rm = TRUE) +
  facet_wrap(~ contrast, nrow = 1) +
  scale_fill_gradient2(
    limits = c(-strata_contrast_lim, strata_contrast_lim),
    oob = scales::squish,
    low = figure_diverging["low"],
    mid = figure_diverging["mid"],
    high = figure_diverging["high"]
  ) +
  labs(x = paste(active_spatial_vars, collapse = " / "), y = NULL, fill = "Adjusted delta") +
  theme_publication() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks = element_blank(),
    panel.spacing.x = grid::unit(0.8, "lines")
  )

save_plot_publication(
  strata_contrast_plot,
  fp_traits("panel_ME_by_condition_spatial_strata_pairwise_contrasts.svg"),
  width = figure_double_col,
  height = max(3.6, 0.18 * length(contrast_module_order) + 1.2)
)
save_plot_publication(
  strata_contrast_plot,
  fp_mainfig("panel_ME_by_condition_spatial_strata_pairwise_contrasts.svg"),
  width = figure_double_col,
  height = max(3.6, 0.18 * length(contrast_module_order) + 1.2)
)

strata_contrast_plot_GO_labels <- strata_contrast_plot +
  scale_y_discrete(labels = function(x) module_axis_labels(x, include_id = TRUE, width = 30))

save_plot_publication(
  strata_contrast_plot_GO_labels,
  fp_traits("panel_ME_by_condition_spatial_strata_pairwise_contrasts_GO_labels.svg"),
  width = figure_double_col,
  height = max(4.2, 0.27 * length(contrast_module_order) + 1.2)
)
save_plot_publication(
  strata_contrast_plot_GO_labels,
  fp_mainfig("panel_ME_by_condition_spatial_strata_pairwise_contrasts_GO_labels.svg"),
  width = figure_double_col,
  height = max(4.2, 0.27 * length(contrast_module_order) + 1.2)
)

# Stats
do_stats <- function(df) {
  covars <- active_spatial_vars
  covars <- covars[vapply(df[covars], function(x) length(unique(stats::na.omit(x))) > 1, logical(1))]
  rhs_full <- paste(c("condition", covars), collapse = " + ")
  rhs_reduced <- if (length(covars)) paste(covars, collapse = " + ") else "1"
  full_formula <- stats::as.formula(paste("ME ~", rhs_full))
  reduced_formula <- stats::as.formula(paste("ME ~", rhs_reduced))
  full_fit <- stats::lm(full_formula, data = df)
  reduced_fit <- stats::lm(reduced_formula, data = df)
  condition_p <- stats::anova(reduced_fit, full_fit)[["Pr(>F)"]][2]
  list(condition_p = condition_p, model_formula = deparse(full_formula))
}

stat_list <- ME_long %>%
  group_by(module) %>%
  group_map(~{
    st <- do_stats(.x)
    tibble(module = .y$module,
           condition_model_p = st$condition_p,
           model_formula = st$model_formula)
  }) %>% bind_rows() %>%
  mutate(condition_model_fdr = p.adjust(condition_model_p, method = "BH")) %>%
  arrange(condition_model_fdr)
stat_list <- add_supermodule_cols(stat_list, supermodule_annotation, module_col = "module")

write_csv_safe(stat_list, fp_traittab("ME_by_condition_adjusted_lm_FDR.csv"))
write_csv_safe(stat_list, fp_supertab("wgcna_condition_omnibus_with_supermodules.csv"))
write_csv_safe(ME_long, fp_source("ME_by_condition.csv"))
write_csv_safe(ME_long, fp_supertab("wgcna_module_eigengenes_long_with_supermodules.csv"))

top_modules <- head(stat_list$module, 12)
comparisons <- list(c("con","res"), c("con","sus"), c("res","sus"))

# Color mapping (full circles)
cond_cols <- figure_condition_cols

# Helper: build one dotplot
plot_dot_mod <- function(dfm, mod, condition_fdr = NA_real_, show_subtitle = TRUE, errorbar = c("none","sem","sd")) {
  errorbar <- match.arg(errorbar)
  # Summary for error bars
  summ <- dfm %>%
    group_by(condition) %>%
    summarise(mean = mean(ME, na.rm = TRUE),
              sd = sd(ME, na.rm = TRUE),
              n = dplyr::n(),
              se = sd / sqrt(pmax(n, 1)), .groups = "drop")
  # Base: all sample dots (full circles)
  p <- ggplot(dfm, aes(x = condition, y = ME, color = condition)) +
    geom_point(position = position_jitter(width = 0.08, height = 0, seed = 1),
               size = 1.4, alpha = 0.75, shape = 16, stroke = 0) +
    scale_color_manual(values = cond_cols, labels = figure_condition_labels, guide = "none") +
    scale_x_discrete(labels = figure_condition_labels)

  # Add mean points (slightly larger) on top
  p <- p + geom_point(data = summ, aes(x = condition, y = mean),
                      inherit.aes = FALSE, size = 2.4, shape = 16, color = "black") +
           geom_point(data = summ, aes(x = condition, y = mean, color = condition),
                      inherit.aes = FALSE, size = 2.0, shape = 16)

  # Optional error bars
  if (errorbar != "none") {
    if (errorbar == "sem") {
      p <- p + geom_errorbar(data = summ,
                             aes(x = condition, ymin = mean - se, ymax = mean + se, color = condition),
                             inherit.aes = FALSE, width = 0.10, linewidth = 0.25, alpha = 0.9)
    } else if (errorbar == "sd") {
      p <- p + geom_errorbar(data = summ,
                             aes(x = condition, ymin = mean - sd, ymax = mean + sd, color = condition),
                             inherit.aes = FALSE, width = 0.10, linewidth = 0.25, alpha = 0.9)
    }
  }

  subtitle_txt <- if (isTRUE(show_subtitle))
    sprintf("adjusted LM FDR=%s", ifelse(is.na(condition_fdr), "NA", signif(condition_fdr, 3))) else NULL

  p <- p +
    labs(title = paste0(mod, " eigengene by condition"),
         subtitle = subtitle_txt,
         x = NULL, y = "Module eigengene") +
    theme_publication() +
    theme(panel.grid.major.x = element_blank())

  # Pairwise significance labels (Wilcoxon BH), shown as p.signif above groups
  p <- p + ggpubr::stat_compare_means(comparisons = comparisons,
                                      method = "wilcox.test",
                                      p.adjust.method = "BH",
                                      label = "p.signif",
                                      hide.ns = TRUE)
  p
}

# Multi-page PDF: all modules as dotplots
all_modules_pdf <- fp_traits("ME_by_condition_all_modules_dotplot.pdf")
all_modules_pdf_tmp <- tempfile(fileext = ".pdf")
grDevices::pdf(all_modules_pdf_tmp, width = 7, height = 4.5, family = figure_font, useDingbats = FALSE)
for (mod in unique(ME_long$module)) {
  dfm <- ME_long %>% filter(module == mod)
  condition_fdr <- stat_list$condition_model_fdr[stat_list$module == mod][1]
  p <- plot_dot_mod(dfm, mod, condition_fdr = condition_fdr, show_subtitle = TRUE, errorbar = "sem")
  print(p)
}
dev.off()
if (!file.copy(all_modules_pdf_tmp, all_modules_pdf, overwrite = TRUE)) {
  warning("Could not overwrite ", all_modules_pdf, "; it may be open in another application.", call. = FALSE)
}
unlink(all_modules_pdf_tmp, force = TRUE)

# SVG grid: top N most differential modules (compact)
plot_one_mod_dot <- function(mod) {
  dfm <- ME_long %>% filter(module == mod)
  condition_fdr <- stat_list$condition_model_fdr[stat_list$module == mod][1]
  summ <- dfm %>% group_by(condition) %>%
    summarise(mean = mean(ME, na.rm = TRUE), se = sd(ME, na.rm = TRUE)/sqrt(n()), .groups = "drop")

  ggplot(dfm, aes(x = condition, y = ME, color = condition)) +
    geom_point(position = position_jitter(width = 0.08, height = 0, seed = 1),
               size = 1.1, alpha = 0.75, shape = 16, stroke = 0) +
    geom_point(data = summ, aes(x = condition, y = mean),
               inherit.aes = FALSE, size = 2.2, shape = 16, color = "black") +
    geom_point(data = summ, aes(x = condition, y = mean, color = condition),
               inherit.aes = FALSE, size = 1.8, shape = 16) +
    # Optional tiny SEM bars for overview
    geom_errorbar(data = summ, aes(x = condition, ymin = mean - se, ymax = mean + se, color = condition),
                  inherit.aes = FALSE, width = 0.08, linewidth = 0.25, alpha = 0.8) +
    scale_color_manual(values = cond_cols, labels = figure_condition_labels, guide = "none") +
    scale_x_discrete(labels = figure_condition_labels) +
    labs(title = paste0(mod, "  (FDR=", signif(condition_fdr, 3), ")"),
         x = NULL, y = NULL) +
    theme_publication() +
    theme(legend.position = "none",
          plot.title = element_text(size = 6.2),
          panel.grid.major.x = element_blank())
}

if (length(top_modules) > 0) {
  plots_top <- lapply(top_modules, plot_one_mod_dot)
  ncol_grid <- min(4, ceiling(sqrt(length(plots_top))))
  nrow_grid <- ceiling(length(plots_top)/ncol_grid)
  g <- patchwork::wrap_plots(plots_top, ncol = ncol_grid)
  save_plot_publication(g, fp_mainfig("ME_by_condition_top_modules_dotplot.svg"),
                   width = min(figure_double_col, 1.75*ncol_grid),
                   height = 1.55*nrow_grid)
}

# Optional: export pairwise Wilcoxon BH summary across modules
pw_tables <- ME_long %>%
  group_by(module) %>%
  group_map(~{
    tt <- pairwise.wilcox.test(.x$ME, .x$condition, p.adjust.method = "BH", exact = FALSE)
    broom::tidy(tt) %>% mutate(module = .y$module)
  }) %>% bind_rows()
write_csv_safe(pw_tables, fp_traittab("ME_by_condition_pairwise_Wilcoxon_BH.csv"))

present_supermodule_annotation <- supermodule_annotation %>%
  dplyr::filter(.data$present_in_dataset) %>%
  dplyr::mutate(
    module_eigengene = factor(.data$module_eigengene, levels = colnames(mergedMEs)),
    Supermodule = factor(.data$Supermodule, levels = supermodule_levels)
  )
module_order_supermodule <- present_supermodule_annotation %>%
  dplyr::arrange(.data$Supermodule, .data$module_eigengene) %>%
  dplyr::pull(.data$module_eigengene) %>%
  as.character()
module_order_supermodule <- intersect(module_order_supermodule, colnames(mergedMEs))

ME_corr_mat <- stats::cor(mergedMEs[, module_order_supermodule, drop = FALSE], use = "pairwise.complete.obs")
ME_corr_export <- tibble::rownames_to_column(as.data.frame(ME_corr_mat), "module_eigengene")
write_csv_safe(ME_corr_export, fp_supertab("wgcna_module_eigengene_correlation_matrix.csv"))
write_csv_safe(ME_corr_export, fp_source("wgcna_module_eigengene_correlation_matrix.csv"))

ME_corr_pairs <- if (length(module_order_supermodule) >= 2) {
  pair_idx <- utils::combn(module_order_supermodule, 2, simplify = FALSE)
  purrr::map_dfr(pair_idx, function(pair) {
    r_val <- unname(ME_corr_mat[pair[[1]], pair[[2]]])
    tibble::tibble(
      module_a = pair[[1]],
      module_b = pair[[2]],
      r = r_val,
      abs_r = abs(r_val)
    )
  })
} else {
  tibble::tibble(module_a = character(), module_b = character(), r = numeric(), abs_r = numeric())
}

pair_ann <- present_supermodule_annotation %>%
  dplyr::select("module_eigengene", "ModuleColor", "Supermodule", "top_GO_label")
ME_corr_pairs <- ME_corr_pairs %>%
  dplyr::left_join(pair_ann, by = c("module_a" = "module_eigengene")) %>%
  dplyr::rename(ModuleColor_a = "ModuleColor", Supermodule_a = "Supermodule", top_GO_label_a = "top_GO_label") %>%
  dplyr::left_join(pair_ann, by = c("module_b" = "module_eigengene")) %>%
  dplyr::rename(ModuleColor_b = "ModuleColor", Supermodule_b = "Supermodule", top_GO_label_b = "top_GO_label") %>%
  dplyr::mutate(same_supermodule = .data$Supermodule_a == .data$Supermodule_b)
write_csv_safe(ME_corr_pairs, fp_supertab("wgcna_module_eigengene_correlation_pairs.csv"))

within_supermodule_corr_summary <- ME_corr_pairs %>%
  dplyr::filter(.data$same_supermodule, .data$Supermodule_a != "Unassigned") %>%
  dplyr::group_by(Supermodule = .data$Supermodule_a) %>%
  dplyr::summarise(
    n_module_pairs = dplyr::n(),
    median_abs_r = stats::median(.data$abs_r, na.rm = TRUE),
    max_abs_r = safe_max(.data$abs_r),
    mean_abs_r = mean(.data$abs_r, na.rm = TRUE),
    eigengene_supported = .data$median_abs_r >= 0.5 | .data$max_abs_r >= 0.6,
    .groups = "drop"
  )

between_supermodule_corr_summary <- ME_corr_pairs %>%
  dplyr::filter(!.data$same_supermodule) %>%
  dplyr::mutate(
    supermodule_pair = paste(
      pmin(as.character(.data$Supermodule_a), as.character(.data$Supermodule_b)),
      pmax(as.character(.data$Supermodule_a), as.character(.data$Supermodule_b)),
      sep = " vs "
    )
  ) %>%
  dplyr::group_by(.data$supermodule_pair) %>%
  dplyr::summarise(
    n_module_pairs = dplyr::n(),
    median_abs_r = stats::median(.data$abs_r, na.rm = TRUE),
    max_abs_r = safe_max(.data$abs_r),
    mean_abs_r = mean(.data$abs_r, na.rm = TRUE),
    .groups = "drop"
  )

write_csv_safe(within_supermodule_corr_summary, fp_supertab("wgcna_within_supermodule_eigengene_correlation_summary.csv"))
write_csv_safe(between_supermodule_corr_summary, fp_supertab("wgcna_between_supermodule_eigengene_correlation_summary.csv"))

hub_module_sets <- WGCNA_modules_long %>%
  add_supermodule_cols(supermodule_annotation, color_col = "ModuleColor") %>%
  dplyr::group_by(.data$module_eigengene, .data$ModuleColor, .data$Supermodule) %>%
  dplyr::arrange(dplyr::desc(.data$abs_kME), .by_group = TRUE) %>%
  dplyr::mutate(.hub_keep = if ("is_top_hub_25" %in% names(.)) .data$is_top_hub_25 | dplyr::row_number() <= 25 else dplyr::row_number() <= 25) %>%
  dplyr::filter(.data$.hub_keep) %>%
  dplyr::summarise(top_hub_proteins = list(unique(as.character(.data$ProteinID))), .groups = "drop")

hub_overlap_pairs <- if (nrow(hub_module_sets) >= 2) {
  pair_idx <- utils::combn(seq_len(nrow(hub_module_sets)), 2, simplify = FALSE)
  purrr::map_dfr(pair_idx, function(idx) {
    a <- hub_module_sets[idx[[1]], ]
    b <- hub_module_sets[idx[[2]], ]
    set_a <- a$top_hub_proteins[[1]]
    set_b <- b$top_hub_proteins[[1]]
    union_n <- length(union(set_a, set_b))
    tibble::tibble(
      module_a = a$module_eigengene,
      module_b = b$module_eigengene,
      Supermodule_a = as.character(a$Supermodule),
      Supermodule_b = as.character(b$Supermodule),
      same_supermodule = as.character(a$Supermodule) == as.character(b$Supermodule),
      hub_overlap_n = length(intersect(set_a, set_b)),
      hub_union_n = union_n,
      hub_jaccard = ifelse(union_n > 0, length(intersect(set_a, set_b)) / union_n, NA_real_)
    )
  })
} else {
  tibble::tibble(
    module_a = character(), module_b = character(), Supermodule_a = character(), Supermodule_b = character(),
    same_supermodule = logical(), hub_overlap_n = integer(), hub_union_n = integer(), hub_jaccard = numeric()
  )
}
write_csv_safe(hub_overlap_pairs, fp_supertab("wgcna_top_hub_overlap_pairs.csv"))

within_supermodule_hub_summary <- hub_overlap_pairs %>%
  dplyr::filter(.data$same_supermodule, .data$Supermodule_a != "Unassigned") %>%
  dplyr::group_by(Supermodule = .data$Supermodule_a) %>%
  dplyr::summarise(
    n_module_pairs = dplyr::n(),
    median_hub_jaccard = stats::median(.data$hub_jaccard, na.rm = TRUE),
    max_hub_jaccard = safe_max(.data$hub_jaccard),
    max_hub_overlap_n = safe_max(.data$hub_overlap_n),
    .groups = "drop"
  )
write_csv_safe(within_supermodule_hub_summary, fp_supertab("wgcna_within_supermodule_hub_overlap_summary.csv"))

supermodule_validation_summary <- present_supermodule_annotation %>%
  dplyr::filter(.data$Supermodule != "Unassigned") %>%
  dplyr::count(.data$Supermodule, name = "n_modules_present") %>%
  dplyr::left_join(within_supermodule_corr_summary, by = "Supermodule") %>%
  dplyr::left_join(within_supermodule_hub_summary, by = "Supermodule", suffix = c("_eigengene", "_hub")) %>%
  dplyr::mutate(
    eigengene_supported = dplyr::coalesce(.data$eigengene_supported, FALSE),
    hub_overlap_observed = dplyr::coalesce(.data$max_hub_overlap_n > 0, FALSE),
    support_call = dplyr::case_when(
      .data$n_modules_present < 2 ~ "single_module_interpretation_group",
      .data$eigengene_supported & .data$hub_overlap_observed ~ "supported_by_eigengenes_and_hub_overlap",
      .data$eigengene_supported ~ "supported_by_eigengene_correlation",
      .data$hub_overlap_observed ~ "supported_by_hub_overlap_only",
      TRUE ~ "weak_or_context_dependent_support"
    )
  ) %>%
  dplyr::arrange(factor(.data$Supermodule, levels = supermodule_levels))
write_csv_safe(supermodule_validation_summary, fp_supertab("wgcna_supermodule_validation_summary.csv"))
write_csv_safe(supermodule_validation_summary, fp_source("wgcna_supermodule_validation_summary.csv"))

ME_corr_long <- reshape2::melt(ME_corr_mat, varnames = c("module_x", "module_y"), value.name = "r") %>%
  dplyr::left_join(pair_ann, by = c("module_x" = "module_eigengene")) %>%
  dplyr::rename(Supermodule_x = "Supermodule") %>%
  dplyr::left_join(pair_ann, by = c("module_y" = "module_eigengene")) %>%
  dplyr::rename(Supermodule_y = "Supermodule") %>%
  dplyr::mutate(
    module_x = factor(.data$module_x, levels = module_order_supermodule),
    module_y = factor(.data$module_y, levels = rev(module_order_supermodule)),
    Supermodule_x = factor(.data$Supermodule_x, levels = supermodule_levels),
    Supermodule_y = factor(.data$Supermodule_y, levels = rev(supermodule_levels))
  )

corr_heatmap <- ggplot(ME_corr_long, aes(x = module_x, y = module_y, fill = r)) +
  geom_tile(color = "white", linewidth = 0.15) +
  facet_grid(Supermodule_y ~ Supermodule_x, scales = "free", space = "free") +
  scale_fill_gradient2(
    limits = c(-1, 1), oob = scales::squish,
    low = figure_diverging["low"], mid = figure_diverging["mid"], high = figure_diverging["high"]
  ) +
  labs(x = NULL, y = NULL, fill = "ME r") +
  theme_publication() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.spacing = grid::unit(0.25, "lines")
  )
save_plot_publication(
  corr_heatmap,
  fp_superfig("wgcna_module_eigengene_correlation_heatmap.svg"),
  width = figure_double_col,
  height = max(4.4, 0.22 * length(module_order_supermodule) + 2.0)
)

module_effects_supermodule_plot_df <- ME_strata_contrast_stats %>%
  dplyr::mutate(
    Supermodule = factor(.data$Supermodule, levels = supermodule_levels),
    module = factor(as.character(.data$module), levels = rev(module_order_supermodule))
  )

module_effects_by_supermodule <- ggplot(module_effects_supermodule_plot_df, aes(x = spatial_trait, y = module, fill = adjusted_delta)) +
  geom_tile(color = "white", linewidth = 0.15) +
  geom_text(aes(label = sig), size = 1.7, color = "black", na.rm = TRUE) +
  facet_grid(Supermodule ~ contrast, scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(x) module_axis_labels(x, include_id = TRUE, width = 26)) +
  scale_fill_gradient2(
    limits = c(-strata_contrast_lim, strata_contrast_lim), oob = scales::squish,
    low = figure_diverging["low"], mid = figure_diverging["mid"], high = figure_diverging["high"]
  ) +
  labs(x = paste(active_spatial_vars, collapse = " / "), y = NULL, fill = "Adjusted delta") +
  theme_publication() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.spacing = grid::unit(0.35, "lines")
  )
save_plot_publication(
  module_effects_by_supermodule,
  fp_superfig("wgcna_module_effects_by_supermodule.svg"),
  width = figure_double_col,
  height = max(5.0, 0.34 * length(module_order_supermodule) + 2.2)
)

if (exists("df_combo2_supermodule") && nrow(df_combo2_supermodule)) {
  trait_supermodule_plot_df <- df_combo2_supermodule %>%
    dplyr::mutate(
      Supermodule = factor(.data$Supermodule, levels = supermodule_levels),
      module = factor(as.character(.data$module), levels = rev(module_order_supermodule))
    )
  module_trait_supermodule_heatmap <- ggplot(trait_supermodule_plot_df, aes(x = spatial_trait, y = module, fill = r)) +
    geom_tile(color = "white", linewidth = 0.15) +
    geom_text(aes(label = sig), size = 1.7, color = "black", na.rm = TRUE) +
    facet_grid(Supermodule ~ condition, scales = "free_y", space = "free_y") +
    scale_y_discrete(labels = function(x) module_axis_labels(x, include_id = TRUE, width = 26)) +
    scale_fill_gradient2(
      limits = c(-1, 1), oob = scales::squish,
      low = figure_diverging["low"], mid = figure_diverging["mid"], high = figure_diverging["high"]
    ) +
    labs(x = paste(active_spatial_vars, collapse = " / "), y = NULL, fill = "Pearson r") +
    theme_publication() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      panel.spacing = grid::unit(0.35, "lines")
    )
  save_plot_publication(
    module_trait_supermodule_heatmap,
    fp_superfig("wgcna_module_trait_correlations_by_supermodule.svg"),
    width = figure_double_col,
    height = max(5.0, 0.34 * length(module_order_supermodule) + 2.2)
  )
}

ME_scaled <- as.data.frame(scale(mergedMEs[, module_order_supermodule, drop = FALSE]))
ME_scaled$Sample <- rownames(mergedMEs)
supermodule_summary_scores <- ME_scaled %>%
  tidyr::pivot_longer(cols = dplyr::all_of(module_order_supermodule), names_to = "module", values_to = "module_z") %>%
  dplyr::left_join(
    present_supermodule_annotation %>%
      dplyr::select("module_eigengene", "Supermodule", "top_GO_label"),
    by = c("module" = "module_eigengene")
  ) %>%
  dplyr::filter(.data$Supermodule != "Unassigned") %>%
  dplyr::group_by(.data$Sample, .data$Supermodule) %>%
  dplyr::summarise(
    aggregate_supermodule_score = mean(.data$module_z, na.rm = TRUE),
    n_member_modules = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    sample_info %>%
      tibble::rownames_to_column("Sample") %>%
      dplyr::select("Sample", "condition", dplyr::all_of(active_spatial_vars)),
    by = "Sample"
  ) %>%
  dplyr::mutate(
    Supermodule = factor(.data$Supermodule, levels = supermodule_levels),
    condition = factor(.data$condition, levels = c("con", "res", "sus"))
  )
write_csv_safe(supermodule_summary_scores, fp_supertab("wgcna_supermodule_summary_scores.csv"))
write_csv_safe(supermodule_summary_scores, fp_source("wgcna_supermodule_summary_scores.csv"))

supermodule_summary_plot_df <- supermodule_summary_scores %>%
  dplyr::group_by(.data$Supermodule, .data$condition) %>%
  dplyr::summarise(
    median_score = stats::median(.data$aggregate_supermodule_score, na.rm = TRUE),
    mean_score = mean(.data$aggregate_supermodule_score, na.rm = TRUE),
    se = stats::sd(.data$aggregate_supermodule_score, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  )
supermodule_summary_scores_plot <- ggplot(supermodule_summary_scores, aes(x = condition, y = aggregate_supermodule_score, color = condition)) +
  geom_point(position = position_jitter(width = 0.08, height = 0, seed = 1), size = 0.8, alpha = 0.45, shape = 16, stroke = 0) +
  geom_point(data = supermodule_summary_plot_df, aes(y = mean_score), size = 1.9, color = "black", inherit.aes = TRUE) +
  geom_errorbar(data = supermodule_summary_plot_df, aes(y = mean_score, ymin = mean_score - se, ymax = mean_score + se), width = 0.12, linewidth = 0.25, color = "black", inherit.aes = TRUE) +
  facet_wrap(~ Supermodule, scales = "free_y", ncol = 2) +
  scale_color_manual(values = figure_condition_cols, labels = figure_condition_labels, guide = "none") +
  scale_x_discrete(labels = figure_condition_labels) +
  labs(x = NULL, y = "Mean standardized module eigengene", title = "Aggregate supermodule score summary") +
  theme_publication() +
  theme(panel.grid.major.x = element_blank())
save_plot_publication(
  supermodule_summary_scores_plot,
  fp_superfig("wgcna_supermodule_summary_scores.svg"),
  width = figure_double_col,
  height = 5.5
)

module_label_export <- module_label_table %>%
  dplyr::mutate(
    ModuleLabel_Final = dplyr::coalesce(.data$ModuleLabel_Final, .data$ModuleLabel_GO_BP, paste0("Module ", .data$ModuleColor)),
    ModuleLabel_Source = dplyr::coalesce(.data$ModuleLabel_Source, "GO_BP_ORA_all_module")
  ) %>%
  dplyr::left_join(
    supermodule_annotation %>%
      dplyr::select(
        "ModuleColor", "module_eigengene", "Supermodule",
        "SupermoduleConfidence", "SupermoduleRationale",
        "present_in_dataset", "manual_annotation"
      ),
    by = "ModuleColor"
  )

module_condition_top <- ME_contrast_stats %>%
  dplyr::mutate(ModuleColor = sub("^ME", "", as.character(.data$module))) %>%
  dplyr::arrange(dplyr::desc(abs(.data$adjusted_delta))) %>%
  dplyr::group_by(.data$ModuleColor) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    ModuleColor,
    strongest_condition_contrast = as.character(.data$contrast),
    strongest_condition_adjusted_delta = .data$adjusted_delta,
    strongest_condition_fdr = .data$fdr
  )

module_strata_top <- ME_strata_contrast_stats %>%
  dplyr::mutate(ModuleColor = sub("^ME", "", as.character(.data$module))) %>%
  dplyr::arrange(dplyr::desc(abs(.data$adjusted_delta))) %>%
  dplyr::group_by(.data$ModuleColor) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    ModuleColor,
    strongest_spatial_stratum = as.character(.data$spatial_trait),
    strongest_spatial_strata_contrast = as.character(.data$contrast),
    strongest_spatial_strata_adjusted_delta = .data$adjusted_delta,
    strongest_spatial_strata_fdr = .data$fdr
  )

module_condition_omnibus <- stat_list %>%
  dplyr::mutate(ModuleColor = sub("^ME", "", as.character(.data$module))) %>%
  dplyr::transmute(
    ModuleColor,
    condition_model_p = .data$condition_model_p,
    condition_model_fdr = .data$condition_model_fdr,
    condition_model_formula = .data$model_formula
  )

WGCNA_module_priority_summary <- WGCNA_module_summary %>%
  dplyr::left_join(
    module_label_export %>%
      dplyr::select(
        "ModuleColor", "ModuleID", "module_eigengene", "ModuleLabel_Final", "ModuleLabel_Source",
        "Supermodule", "SupermoduleConfidence", "SupermoduleRationale"
      ),
    by = c("ModuleColor", "ModuleID")
  ) %>%
  dplyr::left_join(module_condition_omnibus, by = "ModuleColor") %>%
  dplyr::left_join(module_condition_top, by = "ModuleColor") %>%
  dplyr::left_join(module_strata_top, by = "ModuleColor") %>%
  dplyr::mutate(
    priority_score = -log10(pmax(dplyr::coalesce(.data$condition_model_fdr, 1), .Machine$double.xmin)) +
      abs(dplyr::coalesce(.data$strongest_condition_adjusted_delta, 0)) +
      abs(dplyr::coalesce(.data$strongest_spatial_strata_adjusted_delta, 0)),
    priority_rank = dplyr::min_rank(dplyr::desc(.data$priority_score))
  ) %>%
  dplyr::arrange(.data$priority_rank, .data$condition_model_fdr) %>%
  dplyr::select(
    "priority_rank", "priority_score",
    "ModuleID", "ModuleColor", "module_eigengene", "ModuleLabel_Final", "ModuleLabel_Source",
    "Supermodule", "SupermoduleConfidence", "SupermoduleRationale",
    dplyr::any_of(c(
      "n_features", "n_mapped_entrez", "mapping_rate", "median_abs_kME", "mean_abs_kME",
      "condition_model_p", "condition_model_fdr", "strongest_condition_contrast",
      "strongest_condition_adjusted_delta", "strongest_condition_fdr",
      "strongest_spatial_stratum", "strongest_spatial_strata_contrast",
      "strongest_spatial_strata_adjusted_delta", "strongest_spatial_strata_fdr",
      "best_GO_BP", "best_GO_MF", "best_GO_CC", "best_GO_padj_BP", "best_GO_padj_MF", "best_GO_padj_CC",
      "top_hub_proteins", "condition_model_formula"
    )),
    dplyr::contains("preservation_")
  )

WGCNA_module_definitions_for_downstream <- WGCNA_modules_long %>%
  dplyr::select(-dplyr::any_of(c("ModuleLabel_Final", "ModuleLabel_Source"))) %>%
  dplyr::left_join(
    module_label_export %>%
      dplyr::select(
        "ModuleColor", "ModuleID", "module_eigengene", "ModuleLabel_Final", "ModuleLabel_Source",
        "Supermodule", "SupermoduleConfidence", "SupermoduleRationale"
      ),
    by = c("ModuleColor", "ModuleID")
  ) %>%
  dplyr::left_join(
    WGCNA_module_priority_summary %>%
      dplyr::select(
        "ModuleColor", "ModuleID",
        dplyr::any_of(c(
          "condition_model_p", "condition_model_fdr",
          "strongest_condition_contrast", "strongest_condition_adjusted_delta", "strongest_condition_fdr",
          "strongest_spatial_stratum", "strongest_spatial_strata_contrast",
          "strongest_spatial_strata_adjusted_delta", "strongest_spatial_strata_fdr"
        )),
        dplyr::contains("preservation_")
      ) %>%
      dplyr::distinct(),
    by = c("ModuleColor", "ModuleID")
  ) %>%
  dplyr::select(
    "ModuleSet", "ModuleID", "ModuleColor", "module_eigengene", "ModuleLabel_Final", "ModuleLabel_Source",
    "Supermodule", "SupermoduleConfidence", "SupermoduleRationale",
    dplyr::everything()
  )

write_csv_safe(WGCNA_module_priority_summary, fp_modtab("WGCNA_module_priority_summary.csv"))
write_csv_safe(WGCNA_module_definitions_for_downstream, fp_modtab("WGCNA_module_definitions_for_downstream.csv"))
write_csv_safe(WGCNA_module_definitions_for_downstream, fp_supertab("wgcna_module_results_with_supermodules.csv"))
overlap_bridge_script <- repo_path("06_modules_WGCNA", "05_wgcna_de_gsea_overlap.r")
if (file.exists(overlap_bridge_script)) {
  tryCatch({
    source(overlap_bridge_script)
    overlap_summary <- run_wgcna_de_gsea_overlap(dataset_profile, dry_run = FALSE)
    if (file.exists(fp_modtab("WGCNA_module_priority_summary.csv"))) {
      WGCNA_module_priority_summary <- readr::read_csv(fp_modtab("WGCNA_module_priority_summary.csv"), show_col_types = FALSE)
    }
    WGCNA_module_priority_summary <- ensure_module_label_schema(WGCNA_module_priority_summary)
  }, error = function(e) {
    warning("Optional WGCNA-to-DE/GSEA overlap bridge skipped: ", conditionMessage(e), call. = FALSE)
  })
}

WGCNA_module_priority_summary <- ensure_module_label_schema(WGCNA_module_priority_summary)

WGCNA_module_preservation_summary <- WGCNA_module_priority_summary %>%
  dplyr::select("ModuleID", "ModuleColor", dplyr::contains("preservation_")) %>%
  dplyr::mutate(
    preservation_min_Zsummary = {
      z_cols <- grep("Zsummary", names(.), value = TRUE)
      if (length(z_cols)) do.call(pmin, c(dplyr::select(., dplyr::all_of(z_cols)), na.rm = TRUE)) else NA_real_
    },
    preservation_interpretation = dplyr::case_when(
      is.na(.data$preservation_min_Zsummary) ~ "not_estimable",
      .data$preservation_min_Zsummary >= 10 ~ "strong",
      .data$preservation_min_Zsummary >= 2 ~ "moderate",
      TRUE ~ "weak_or_unstable"
    ),
    preservation_warning = dplyr::case_when(
      is.na(.data$preservation_min_Zsummary) ~ "Preservation Zsummary not estimable for at least one comparison.",
      .data$preservation_min_Zsummary < 2 ~ "Low preservation support; interpret module transferability cautiously.",
      TRUE ~ NA_character_
    )
  )
write_csv_safe(WGCNA_module_preservation_summary, fp_modtab("WGCNA_module_preservation_summary.csv"))

gsea_overlap_file <- path_results("tables", "06_modules_WGCNA", "05_wgcna_de_gsea_overlap", dataset_profile, "WGCNA_vs_DE_GSEA_overlap.csv")
WGCNA_module_GSEA_coregene_overlap <- if (file.exists(gsea_overlap_file)) {
  readr::read_csv(gsea_overlap_file, show_col_types = FALSE) %>%
    dplyr::select(dplyr::any_of(c(
      "dataset", "ModuleID", "ModuleColor", "contrast", "route", "route_category",
      "n_DE_overlap", "n_leading_edge_overlap", "jaccard_DE", "fisher_p", "fisher_fdr",
      "top_overlap_proteins"
    )))
} else {
  tibble::tibble(
    dataset = dataset_profile,
    status = "skipped",
    reason = "Optional WGCNA-to-DE/GSEA overlap file was not available.",
    expected_file = gsea_overlap_file
  )
}
write_csv_safe(WGCNA_module_GSEA_coregene_overlap, fp_modtab("WGCNA_module_GSEA_coregene_overlap.csv"))

WGCNA_module_evidence_rank <- WGCNA_module_priority_summary %>%
  dplyr::left_join(
    WGCNA_module_preservation_summary %>%
      dplyr::select("ModuleID", "ModuleColor", "preservation_min_Zsummary", "preservation_interpretation"),
    by = c("ModuleID", "ModuleColor")
  ) %>%
  dplyr::mutate(
    condition_rank_component = -log10(pmax(dplyr::coalesce(.data$condition_model_fdr, 1), .Machine$double.xmin)),
    strongest_condition_rank_component = -log10(pmax(dplyr::coalesce(.data$strongest_condition_fdr, 1), .Machine$double.xmin)),
    effect_rank_component = abs(dplyr::coalesce(.data$strongest_condition_adjusted_delta, 0)),
    hub_rank_component = dplyr::coalesce(.data$median_abs_kME, .data$mean_abs_kME, 0),
    preservation_rank_component = dplyr::case_when(
      is.na(.data$preservation_min_Zsummary) ~ 0,
      .data$preservation_min_Zsummary >= 10 ~ 2,
      .data$preservation_min_Zsummary >= 2 ~ 1,
      TRUE ~ 0
    ),
    gsea_overlap_rank_component = dplyr::coalesce(.data$n_leading_edge_overlap, 0) + dplyr::coalesce(.data$n_DE_overlap, 0),
    label_confidence_component = dplyr::case_when(
      !is.na(.data$best_GO_padj_BP) & .data$best_GO_padj_BP <= 0.05 ~ 1,
      !is.na(.data$best_GO_padj_BP) & .data$best_GO_padj_BP <= 0.10 ~ 0.5,
      TRUE ~ 0
    ),
    evidence_score = .data$condition_rank_component +
      .data$strongest_condition_rank_component +
      .data$effect_rank_component +
      .data$hub_rank_component +
      .data$preservation_rank_component +
      log1p(.data$gsea_overlap_rank_component) +
      .data$label_confidence_component,
    evidence_rank = dplyr::min_rank(dplyr::desc(.data$evidence_score)),
    evidence_warning = dplyr::case_when(
      is.na(.data$condition_model_fdr) ~ "Condition model unavailable; rank is exploratory.",
      is.na(.data$preservation_min_Zsummary) ~ "Preservation unavailable; rank is exploratory.",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::arrange(.data$evidence_rank)
write_csv_safe(WGCNA_module_evidence_rank, fp_modtab("WGCNA_module_evidence_rank.csv"))

writexl::write_xlsx(
  list(
    priority_summary = WGCNA_module_priority_summary,
    downstream_definitions = WGCNA_module_definitions_for_downstream,
    module_name_map = module_label_export,
    evidence_rank = WGCNA_module_evidence_rank,
    preservation_summary = WGCNA_module_preservation_summary
  ),
  fp_modtab("WGCNA_module_contracts.xlsx")
)

manifest_category <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  roots <- vapply(subdirs, normalizePath, character(1), winslash = "/", mustWork = FALSE)
  hit <- names(roots)[vapply(roots, function(root) startsWith(path, root), logical(1))]
  if (length(hit)) hit[which.max(nchar(roots[hit]))] else "unclassified"
}
managed_files <- unique(unlist(lapply(unlist(subdirs), function(d) {
  if (!dir.exists(d)) character() else list.files(d, recursive = TRUE, full.names = TRUE)
}), use.names = FALSE))
manifest_paths <- unique(c(
  managed_files,
  fp_log("output_manifest.csv"),
  fp_log("run_manifest.yml")
))
output_manifest <- tibble::tibble(
  category = vapply(manifest_paths, manifest_category, character(1)),
  path = normalizePath(manifest_paths, winslash = "/", mustWork = FALSE),
  relative_path = relative_to(manifest_paths),
  exists = file.exists(manifest_paths),
  md5 = vapply(manifest_paths, file_hash, character(1))
) %>%
  dplyr::arrange(.data$category, .data$relative_path)
write_csv_safe(output_manifest, fp_log("output_manifest.csv"))
write_run_manifest(
  fp_log("run_manifest.yml"),
  inputs = as.list(stats::setNames(input_manifest$path, input_manifest$role)),
  outputs = list(
    output_manifest = fp_log("output_manifest.csv"),
    output_root = output_dir,
    managed_folders = subdirs
  ),
  parameters = list(
    soft_threshold_rsquared = soft_threshold_rsquared,
    selected_soft_power = if (exists("softPower")) softPower else NA,
    min_module_size = min_module_size,
    deep_split = deep_split,
    merge_cut_height = merge_cut_height,
    module_preservation_permutations = module_preservation_permutations,
    dataset_profile_requested = dataset_profile,
    dataset_profile_resolved = if (exists("dataset_profile_resolved")) dataset_profile_resolved else NA_character_,
    module_definition_contract = fp_modtab("WGCNA_module_definitions_for_downstream.csv"),
    module_priority_summary = fp_modtab("WGCNA_module_priority_summary.csv"),
    output_layout = ifelse(nzchar(output_dir_env), "custom_bundle", "canonical_module_paths")
  ),
  notes = "Dataset-scoped WGCNA module engine; see input_manifest.csv for exact inputs/hashes and WGCNA_module_contracts.xlsx for downstream module definitions."
)
output_manifest <- output_manifest %>%
  dplyr::mutate(
    exists = file.exists(.data$path),
    md5 = ifelse(
      basename(.data$path) == "output_manifest.csv",
      NA_character_,
      vapply(.data$path, file_hash, character(1))
    )
  )
write_csv_safe(output_manifest, fp_log("output_manifest.csv"))


# --------------------------
# End of script
# --------------------------

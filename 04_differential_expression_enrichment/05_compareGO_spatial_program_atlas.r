#' Spatial Program Atlas from compareGO outputs
#'
#' Consumes existing compareGO outputs across dataset families and synthesizes
#' manuscript-ready spatial program summaries and SVG figures. This script does
#' not rerun clusterProfiler or compareGO.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
dataset_config_file <- repo_path("R", "dataset_config.R")
if (!file.exists(dataset_config_file)) {
  stop("Missing required config: ", dataset_config_file, call. = FALSE)
}
source(dataset_config_file)

if (!exists("proteomics_dataset_order", inherits = TRUE)) {
  proteomics_dataset_order <- if (exists("valid_datasets", inherits = TRUE)) valid_datasets() else c("neuron_neuropil", "neuron_soma", "microglia")
}
if (!exists("dataset_label", inherits = TRUE)) {
  dataset_label <- function(dataset) {
    labels <- c(neuron_neuropil = "Neuron neuropil", neuron_soma = "Neuron soma", microglia = "Microglia")
    out <- unname(labels[dataset])
    ifelse(is.na(out), dataset, out)
  }
}
if (!exists("dataset_compartment", inherits = TRUE)) {
  dataset_compartment <- function(dataset) {
    compartments <- c(neuron_neuropil = "neuropil", neuron_soma = "soma", microglia = "microglia")
    out <- unname(compartments[dataset])
    ifelse(is.na(out), dataset, out)
  }
}

MODULE_ID <- "04_differential_expression_enrichment"
SUBSTEP_ID <- "compareGO_spatial_atlas"
CANONICAL_PATHS <- create_module_dirs(MODULE_ID, SUBSTEP_ID)

required_pkgs <- c("dplyr", "tidyr", "purrr", "readr", "readxl", "writexl", "ggplot2", "stringr", "tibble", "scales")

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

load_required_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing required R package(s): ", paste(missing, collapse = ", "),
         ". Install them explicitly before running this script.", call. = FALSE)
  }
  suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only = TRUE)))
}

# Clearly editable GO program dictionary. Patterns are case-insensitive.
program_dictionary <- list(
  RNA_RNP_processing = c("rna", "ribonucleoprotein", "splice", "splicing", "mrna", "rrna", "trna", "transcription", "ribonucl"),
  Translation_Ribosome = c("translation", "ribosom", "peptide biosynthetic", "cytoplasmic translation"),
  Mitochondria_OXPHOS = c("mitochond", "oxidative phosphorylation", "electron transport", "respiratory chain", "atp synthesis", "nadh", "proton transmembrane"),
  Proteostasis_Lysosome = c("proteas", "proteolysis", "ubiquitin", "lysosom", "autophag", "vacuolar", "folding", "chaperone"),
  Synapse_Plasticity = c("synap", "dendrit", "axon", "neurotransmitter", "vesicle", "plasticity", "glutamate", "postsynaptic", "presynaptic"),
  Immune_Microglia = c("immune", "inflamm", "complement", "microgl", "cytokine", "interferon", "antigen", "phagocyt", "leukocyte", "myeloid"),
  Cytoskeleton_Transport = c("cytoskeleton", "actin", "tubulin", "microtub", "transport", "traffick", "vesicle-mediated", "motor protein")
)

program_levels <- c(names(program_dictionary), "Other")
contrast_levels <- c("RES_vs_CON", "SUS_vs_CON", "SUS_vs_RES")
region_levels <- c("CA1", "CA2", "CA3", "DG")

is_truthy_env <- function(name) {
  tolower(Sys.getenv(name, unset = "")) %in% c("1", "true", "yes", "y")
}

is_script_dry_run <- function() {
  is_dry_run() || is_truthy_env("PROTEOMICS_SPATIAL_ATLAS_DRY_RUN")
}

atlas_theme <- function(base_size = 8) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      strip.background = ggplot2::element_rect(fill = "grey95", color = "grey70", linewidth = 0.3),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      panel.spacing = grid::unit(0.8, "lines")
    )
}

get_requested_datasets <- function() {
  raw <- Sys.getenv("PROTEOMICS_SPATIAL_ATLAS_DATASETS", unset = "neuron_neuropil,neuron_soma,microglia")
  datasets <- trimws(strsplit(raw, ",", fixed = TRUE)[[1]])
  datasets[nzchar(datasets)]
}

discover_dataset_files <- function(dataset) {
  table_dir <- path_results("tables", MODULE_ID, "compareGO", dataset, "BP", "phenotype_within_unit")
  processed_dir <- path_processed(MODULE_ID, "compareGO", dataset)
  list(
    dataset = dataset,
    supplementary = file.path(table_dir, "Supplementary_Data.xlsx"),
    drivers = file.path(table_dir, "07_Top_Genes_Driving_TopTerms.xlsx"),
    manifest = file.path(processed_dir, "compareGO_input_manifest.csv")
  )
}

diagnose_inputs <- function(files) {
  tibble::tibble(
    dataset = files$dataset,
    input_type = c("supplementary", "drivers", "manifest"),
    path = c(files$supplementary, files$drivers, files$manifest),
    exists = file.exists(c(files$supplementary, files$drivers, files$manifest))
  )
}

read_supplementary_enrichment <- function(path) {
  sheets <- readxl::excel_sheets(path)
  sheet <- dplyr::case_when(
    "GO_Enrichment_Results" %in% sheets ~ "GO_Enrichment_Results",
    "Sheet1" %in% sheets ~ "Sheet1",
    TRUE ~ sheets[[1]]
  )
  readxl::read_excel(path, sheet = sheet) %>% tibble::as_tibble()
}

read_driver_table <- function(path) {
  if (!file.exists(path)) {
    return(tibble::tibble(Description = character(), Gene = character(), Freq = integer(), Mean_NES = numeric(), Comparisons = character()))
  }
  readxl::read_excel(path, sheet = 1) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      Description = as.character(.data$Description),
      Gene = as.character(.data$Gene),
      Comparisons = as.character(.data$Comparisons)
    )
}

read_manifest <- function(path) {
  if (!file.exists(path)) {
    return(tibble::tibble(comparison = character(), route_unit = character()))
  }
  readr::read_csv(path, show_col_types = FALSE) %>%
    dplyr::distinct(comparison, route_unit)
}

parse_unit <- function(unit, dataset) {
  unit <- as.character(unit)
  unit <- gsub("_", "", unit)
  m <- stringr::str_match(unit, "^(CA1|CA2|CA3|DG)(.*)$")
  region <- ifelse(is.na(m[, 2]), NA_character_, m[, 2])
  layer <- ifelse(is.na(m[, 3]) | !nzchar(m[, 3]), NA_character_, m[, 3])
  if (identical(dataset, "microglia")) layer <- NA_character_
  tibble::tibble(region = region, layer = layer)
}

parse_comparison_name <- function(comparison, dataset, route_unit = NA_character_) {
  parts <- strsplit(as.character(comparison), "_", fixed = TRUE)[[1]]
  left <- parts[[1]]
  right <- if (length(parts) >= 2) parts[[2]] else NA_character_

  parse_side <- function(x) {
    m <- stringr::str_match(tolower(x), "^(.*?)(con|res|sus)$")
    tibble::tibble(unit = m[, 2], phenotype = toupper(m[, 3]))
  }

  left_parsed <- parse_side(left)
  right_parsed <- parse_side(right)
  unit <- if (!is.na(route_unit) && nzchar(as.character(route_unit))) as.character(route_unit) else left_parsed$unit
  unit_parsed <- parse_unit(unit, dataset)
  contrast <- paste0(left_parsed$phenotype, "_vs_", right_parsed$phenotype)

  tibble::tibble(
    Comparison = comparison,
    phenotype_contrast = contrast,
    region = unit_parsed$region,
    layer = unit_parsed$layer,
    compartment = dataset_compartment(dataset),
    spatial_unit = ifelse(is.na(unit_parsed$layer), unit_parsed$region, paste(unit_parsed$region, unit_parsed$layer, sep = "_"))
  )
}

classify_program <- function(description) {
  desc <- tolower(as.character(description))
  out <- rep("Other", length(desc))
  for (class_name in names(program_dictionary)) {
    pattern <- paste(program_dictionary[[class_name]], collapse = "|")
    hit <- grepl(pattern, desc, ignore.case = TRUE)
    out[hit & out == "Other"] <- class_name
  }
  factor(out, levels = program_levels)
}

split_genes <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  unique(trimws(unlist(strsplit(x, "/|;|,|\\s+"))))
}

collapse_top_terms <- function(df, n = 5) {
  df %>%
    dplyr::arrange(dplyr::desc(!is.na(.data$p.adjust) & .data$p.adjust < 0.05), .data$p.adjust, dplyr::desc(abs(.data$NES))) %>%
    dplyr::pull(.data$Description) %>%
    unique() %>%
    utils::head(n) %>%
    paste(collapse = "; ")
}

top_driver_genes_for_group <- function(driver_df, descriptions, comparisons, n = 10) {
  if (nrow(driver_df) == 0) return(NA_character_)
  comp_pattern <- paste(unique(comparisons), collapse = "|")
  out <- driver_df %>%
    dplyr::filter(.data$Description %in% descriptions) %>%
    dplyr::filter(is.na(.data$Comparisons) | grepl(comp_pattern, .data$Comparisons)) %>%
    dplyr::group_by(.data$Gene) %>%
    dplyr::summarise(
      total_freq = sum(suppressWarnings(as.numeric(.data$Freq)), na.rm = TRUE),
      mean_nes = mean(suppressWarnings(as.numeric(.data$Mean_NES)), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(.data$total_freq), dplyr::desc(.data$mean_nes)) %>%
    dplyr::slice_head(n = n) %>%
    dplyr::pull(.data$Gene)
  if (length(out) == 0) NA_character_ else paste(out, collapse = "; ")
}

load_dataset_atlas <- function(files) {
  enrichment <- read_supplementary_enrichment(files$supplementary)
  needed <- c("Description", "NES", "p.adjust", "core_enrichment", "Comparison")
  missing <- setdiff(needed, names(enrichment))
  if (length(missing) > 0) {
    stop("Supplementary enrichment table for ", files$dataset, " is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  manifest <- read_manifest(files$manifest)
  drivers <- read_driver_table(files$drivers)

  parsed <- enrichment %>%
    dplyr::mutate(
      dataset = files$dataset,
      dataset_label = dataset_label(files$dataset),
      NES = suppressWarnings(as.numeric(.data$NES)),
      p.adjust = suppressWarnings(as.numeric(.data$p.adjust)),
      Comparison = as.character(.data$Comparison),
      Description = as.character(.data$Description),
      core_enrichment = as.character(.data$core_enrichment)
    ) %>%
    dplyr::left_join(manifest, by = c("Comparison" = "comparison")) %>%
    dplyr::mutate(route_unit = dplyr::coalesce(as.character(.data$route_unit), NA_character_))

  comparison_meta <- purrr::pmap_dfr(
    list(parsed$Comparison, parsed$dataset, parsed$route_unit),
    parse_comparison_name
  ) %>%
    dplyr::distinct(.data$Comparison, .keep_all = TRUE)

  parsed <- parsed %>%
    dplyr::left_join(comparison_meta, by = "Comparison") %>%
    dplyr::filter(.data$phenotype_contrast %in% contrast_levels) %>%
    dplyr::mutate(
      program_class = classify_program(.data$Description),
      spatial_unit_order = factor(.data$spatial_unit, levels = order_spatial_units(unique(.data$spatial_unit))),
      phenotype_contrast = factor(.data$phenotype_contrast, levels = contrast_levels)
    )

  list(enrichment = parsed, drivers = drivers)
}

order_spatial_units <- function(units) {
  units <- unique(as.character(units))
  unit_df <- tibble::tibble(spatial_unit = units) %>%
    tidyr::separate(.data$spatial_unit, into = c("region", "layer"), sep = "_", fill = "right", remove = FALSE) %>%
    dplyr::mutate(
      region = factor(.data$region, levels = region_levels),
      layer = dplyr::coalesce(.data$layer, "")
    ) %>%
    dplyr::arrange(.data$region, .data$layer)
  unit_df$spatial_unit
}

calculate_program_summary <- function(enrichment_df, driver_by_dataset) {
  base_summary <- enrichment_df %>%
    dplyr::group_by(.data$dataset, .data$dataset_label, .data$region, .data$layer, .data$spatial_unit,
                    .data$compartment, .data$phenotype_contrast, .data$program_class) %>%
    dplyr::summarise(
      n_terms = dplyr::n_distinct(.data$Description),
      n_sig_terms = dplyr::n_distinct(.data$Description[!is.na(.data$p.adjust) & .data$p.adjust < 0.05]),
      mean_NES = mean(.data$NES, na.rm = TRUE),
      median_NES = median(.data$NES, na.rm = TRUE),
      mean_signed_log10FDR = mean(sign(.data$NES) * -log10(pmax(.data$p.adjust, .Machine$double.xmin)), na.rm = TRUE),
      leading_edge_union_size = length(unique(unlist(purrr::map(.data$core_enrichment, split_genes)))),
      top_GO_terms = collapse_top_terms(dplyr::pick(dplyr::everything())),
      comparisons = paste(unique(.data$Comparison), collapse = ";"),
      term_set = list(unique(.data$Description)),
      comparison_set = list(unique(.data$Comparison)),
      .groups = "drop"
    )

  base_summary$top_driver_genes <- purrr::pmap_chr(
    list(base_summary$dataset, base_summary$term_set, base_summary$comparison_set),
    function(dataset, term_set, comparison_set) {
      top_driver_genes_for_group(
        driver_by_dataset[[dataset]] %||% tibble::tibble(),
        unlist(term_set),
        unlist(comparison_set)
      ) %||% NA_character_
    }
  )

  base_summary %>%
    dplyr::select(-dplyr::all_of(c("term_set", "comparison_set")))
}

classify_program_behavior <- function(summary_df, min_abs_nes = 0.15) {
  wide <- summary_df %>%
    dplyr::group_by(.data$dataset, .data$dataset_label, .data$region, .data$layer, .data$spatial_unit,
                    .data$compartment, .data$program_class, .data$phenotype_contrast) %>%
    dplyr::summarise(mean_NES = mean(.data$mean_NES, na.rm = TRUE), n_sig_terms = sum(.data$n_sig_terms), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from = "phenotype_contrast",
      values_from = c("mean_NES", "n_sig_terms"),
      values_fill = list(mean_NES = NA_real_, n_sig_terms = 0)
    )

  get_col <- function(df, nm, default = NA_real_) {
    if (nm %in% names(df)) df[[nm]] else rep(default, nrow(df))
  }

  res <- get_col(wide, "mean_NES_RES_vs_CON")
  sus <- get_col(wide, "mean_NES_SUS_vs_CON")
  sus_res <- get_col(wide, "mean_NES_SUS_vs_RES")
  res_sig <- get_col(wide, "n_sig_terms_RES_vs_CON", 0)
  sus_sig <- get_col(wide, "n_sig_terms_SUS_vs_CON", 0)
  sr_sig <- get_col(wide, "n_sig_terms_SUS_vs_RES", 0)

  wide %>%
    dplyr::mutate(
      RES_vs_CON_mean_NES = res,
      SUS_vs_CON_mean_NES = sus,
      SUS_vs_RES_mean_NES = sus_res,
      behavior_class = dplyr::case_when(
        !is.na(res) & !is.na(sus) & abs(res) >= min_abs_nes & abs(sus) >= min_abs_nes & sign(res) != sign(sus) ~ "divergent",
        (sr_sig > 0 | (!is.na(sus_res) & abs(sus_res) >= min_abs_nes)) & !is.na(res) & !is.na(sus) ~ "phenotype_separating",
        res_sig > 0 & sus_sig > 0 & !is.na(res) & !is.na(sus) & sign(res) == sign(sus) ~ "shared_stress",
        res_sig > 0 & (sus_sig == 0 | is.na(sus) | abs(sus) < min_abs_nes) ~ "resilience_specific",
        sus_sig > 0 & (res_sig == 0 | is.na(res) | abs(res) < min_abs_nes) ~ "susceptibility_specific",
        TRUE ~ "phenotype_separating"
      )
    )
}

make_driver_recurrence <- function(enrichment_df, max_genes_per_class = 20) {
  long_genes <- enrichment_df %>%
    dplyr::filter(!is.na(.data$p.adjust), .data$p.adjust < 0.05) %>%
    dplyr::mutate(Gene = purrr::map(.data$core_enrichment, split_genes)) %>%
    tidyr::unnest("Gene") %>%
    dplyr::filter(!is.na(.data$Gene), nzchar(.data$Gene)) %>%
    dplyr::distinct(.data$dataset, .data$dataset_label, .data$spatial_unit, .data$program_class, .data$Gene, .data$Description)

  top_genes <- long_genes %>%
    dplyr::group_by(.data$program_class, .data$Gene) %>%
    dplyr::summarise(n_spatial_units = dplyr::n_distinct(paste(.data$dataset, .data$spatial_unit)), n_terms = dplyr::n_distinct(.data$Description), .groups = "drop") %>%
    dplyr::arrange(.data$program_class, dplyr::desc(.data$n_spatial_units), dplyr::desc(.data$n_terms)) %>%
    dplyr::group_by(.data$program_class) %>%
    dplyr::slice_head(n = max_genes_per_class) %>%
    dplyr::ungroup()

  long_genes %>%
    dplyr::semi_join(top_genes, by = c("program_class", "Gene")) %>%
    dplyr::group_by(.data$program_class, .data$Gene, .data$dataset_label, .data$spatial_unit) %>%
    dplyr::summarise(n_terms = dplyr::n_distinct(.data$Description), .groups = "drop") %>%
    dplyr::mutate(unit_label = paste(.data$dataset_label, .data$spatial_unit, sep = " | "))
}

plot_spatial_program_atlas <- function(summary_df, output_file) {
  plot_df <- summary_df %>%
    dplyr::mutate(
      spatial_unit = factor(.data$spatial_unit, levels = order_spatial_units(unique(.data$spatial_unit))),
      program_class = factor(.data$program_class, levels = rev(program_levels)),
      phenotype_contrast = factor(.data$phenotype_contrast, levels = contrast_levels)
    )
  max_abs <- max(abs(plot_df$mean_NES), na.rm = TRUE)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$spatial_unit, y = .data$program_class, color = .data$mean_NES, size = .data$n_sig_terms)) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::facet_grid(.data$dataset_label ~ .data$phenotype_contrast, scales = "free_x", space = "free_x") +
    ggplot2::scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-max_abs, max_abs), name = "Mean NES") +
    ggplot2::scale_size_continuous(range = c(0.6, 5), name = "Sig. terms") +
    ggplot2::labs(x = NULL, y = NULL, title = "Spatial GO Program Atlas") +
    atlas_theme(7)
  ggplot2::ggsave(output_file, p, width = 12, height = 7, device = "svg", limitsize = FALSE)
}

plot_res_sus_divergence <- function(behavior_df, output_file) {
  plot_df <- behavior_df %>%
    dplyr::filter(!is.na(.data$RES_vs_CON_mean_NES), !is.na(.data$SUS_vs_CON_mean_NES)) %>%
    dplyr::mutate(program_class = factor(.data$program_class, levels = program_levels))
  if (nrow(plot_df) == 0) return(invisible(FALSE))
  lim <- max(abs(c(plot_df$RES_vs_CON_mean_NES, plot_df$SUS_vs_CON_mean_NES)), na.rm = TRUE)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(.data$RES_vs_CON_mean_NES, .data$SUS_vs_CON_mean_NES, label = .data$spatial_unit, color = .data$behavior_class)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "grey75", linewidth = 0.3, linetype = "dashed") +
    ggplot2::geom_point(alpha = 0.85, size = 1.8) +
    ggplot2::geom_text(size = 2, check_overlap = TRUE, vjust = -0.5) +
    ggplot2::facet_grid(.data$dataset_label ~ .data$program_class) +
    ggplot2::coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    ggplot2::labs(x = "RES vs CON mean NES", y = "SUS vs CON mean NES", color = "Behavior", title = "Resilience and Susceptibility Program Divergence") +
    atlas_theme(7)
  ggplot2::ggsave(output_file, p, width = 13, height = 7, device = "svg", limitsize = FALSE)
}

plot_driver_recurrence <- function(driver_df, output_file) {
  if (nrow(driver_df) == 0) return(invisible(FALSE))
  driver_df <- driver_df %>%
    dplyr::mutate(
      gene_label = paste(.data$program_class, .data$Gene, sep = " | "),
      gene_label = factor(.data$gene_label, levels = rev(unique(.data$gene_label[order(as.character(.data$program_class), .data$Gene)]))),
      unit_label = factor(.data$unit_label, levels = unique(.data$unit_label))
    )
  p <- ggplot2::ggplot(driver_df, ggplot2::aes(.data$unit_label, .data$gene_label, fill = .data$n_terms)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradient(low = "white", high = "#3B7EA1", name = "Terms") +
    ggplot2::labs(x = NULL, y = "Program | leading-edge protein", title = "Recurrent Leading-Edge Driver Proteins") +
    atlas_theme(6) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 5))
  ggplot2::ggsave(output_file, p, width = 13, height = max(6, length(levels(driver_df$gene_label)) * 0.12), device = "svg", limitsize = FALSE)
}

plot_compartment_comparison <- function(summary_df, output_file) {
  plot_df <- summary_df %>%
    dplyr::group_by(.data$dataset_label, .data$compartment, .data$region, .data$phenotype_contrast, .data$program_class) %>%
    dplyr::summarise(mean_NES = mean(.data$mean_NES, na.rm = TRUE), n_sig_terms = sum(.data$n_sig_terms), .groups = "drop") %>%
    dplyr::mutate(
      region = factor(.data$region, levels = region_levels),
      program_class = factor(.data$program_class, levels = rev(program_levels))
    )
  max_abs <- max(abs(plot_df$mean_NES), na.rm = TRUE)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(.data$region, .data$program_class, color = .data$mean_NES, size = .data$n_sig_terms)) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::facet_grid(.data$dataset_label ~ .data$phenotype_contrast) +
    ggplot2::scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-max_abs, max_abs), name = "Mean NES") +
    ggplot2::scale_size_continuous(range = c(0.6, 5), name = "Sig. terms") +
    ggplot2::labs(x = NULL, y = NULL, title = "Compartment-Level Program Comparison by Region") +
    atlas_theme(7)
  ggplot2::ggsave(output_file, p, width = 10, height = 7, device = "svg", limitsize = FALSE)
}

write_outputs <- function(enrichment_df, summary_df, behavior_df, driver_df, diagnostics) {
  source_dir <- CANONICAL_PATHS$source_data
  table_dir <- CANONICAL_PATHS$tables
  report_dir <- CANONICAL_PATHS$reports

  readr::write_csv(diagnostics, file.path(report_dir, "input_diagnostics.csv"))
  readr::write_csv(enrichment_df, file.path(source_dir, "spatial_atlas_enrichment_long.csv"))
  readr::write_csv(summary_df, file.path(table_dir, "spatial_program_summary.csv"))
  readr::write_csv(behavior_df, file.path(table_dir, "spatial_program_behavior.csv"))
  readr::write_csv(driver_df, file.path(source_dir, "leading_edge_driver_recurrence.csv"))
  writexl::write_xlsx(
    list(
      input_diagnostics = diagnostics,
      spatial_program_summary = summary_df,
      spatial_program_behavior = behavior_df,
      leading_edge_driver_recurrence = driver_df
    ),
    file.path(table_dir, "compareGO_spatial_program_atlas_tables.xlsx")
  )
}

main <- function() {
  load_required_packages(required_pkgs)

  datasets <- get_requested_datasets()
  files <- purrr::map(datasets, discover_dataset_files)
  diagnostics <- purrr::map_dfr(files, diagnose_inputs)

  if (is_script_dry_run()) {
    message("[DRY-RUN] compareGO spatial atlas")
    print(diagnostics)
    missing <- diagnostics %>% dplyr::filter(!.data$exists)
    if (nrow(missing) > 0) {
      message("[DRY-RUN] Missing dataset outputs detected. Existing datasets will be used in non-dry-run mode only if all required files exist.")
    }
    return(invisible(diagnostics))
  }

  complete_datasets <- diagnostics %>%
    dplyr::group_by(.data$dataset) %>%
    dplyr::summarise(complete = all(.data$exists), .groups = "drop")
  missing_datasets <- complete_datasets$dataset[!complete_datasets$complete]
  if (length(missing_datasets) > 0) {
    warning("Skipping dataset(s) with missing compareGO atlas inputs: ", paste(missing_datasets, collapse = ", "))
  }

  files <- files[vapply(files, function(x) !(x$dataset %in% missing_datasets), logical(1))]
  if (length(files) == 0) {
    stop("No complete dataset compareGO outputs available for atlas synthesis.", call. = FALSE)
  }

  loaded <- purrr::map(files, load_dataset_atlas)
  enrichment_df <- purrr::map_dfr(loaded, "enrichment")
  driver_by_dataset <- stats::setNames(purrr::map(loaded, "drivers"), vapply(files, `[[`, character(1), "dataset"))

  program_summary <- calculate_program_summary(enrichment_df, driver_by_dataset)
  behavior_df <- classify_program_behavior(program_summary)
  driver_recurrence <- make_driver_recurrence(enrichment_df)

  write_outputs(enrichment_df, program_summary, behavior_df, driver_recurrence, diagnostics)

  plot_spatial_program_atlas(program_summary, file.path(CANONICAL_PATHS$figures, "Fig_SpatialProgramAtlas_dotheatmap.svg"))
  plot_res_sus_divergence(behavior_df, file.path(CANONICAL_PATHS$figures, "Fig_RES_SUS_divergence.svg"))
  plot_driver_recurrence(driver_recurrence, file.path(CANONICAL_PATHS$figures, "Fig_LeadingEdgeDriverRecurrence.svg"))
  plot_compartment_comparison(program_summary, file.path(CANONICAL_PATHS$figures, "Fig_Compartment_Comparison.svg"))

  writeLines(
    c(
      "compareGO spatial program atlas complete",
      paste0("Datasets requested: ", paste(datasets, collapse = ", ")),
      paste0("Datasets analyzed: ", paste(vapply(files, `[[`, character(1), "dataset"), collapse = ", ")),
      paste0("Program summary rows: ", nrow(program_summary)),
      paste0("Figures directory: ", CANONICAL_PATHS$figures)
    ),
    file.path(CANONICAL_PATHS$reports, "compareGO_spatial_program_atlas_summary.txt")
  )

  message("Spatial program atlas complete: ", CANONICAL_PATHS$figures)
  invisible(list(summary = program_summary, behavior = behavior_df, driver_recurrence = driver_recurrence))
}

if (sys.nframe() == 0) {
  main()
}

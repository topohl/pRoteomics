# ================================================================
# Consumes:
#   - spatial network RDS from canonical spatial-network output
#   - external behavior/movement tables from data/external/behavior or local config override
# Produces:
#   - network-behavior coupling tables, figures and logs in canonical folders
# File contract:
#   - docs/active_script_io_audit.tsv object 08_behavior_physio_coupling/02_network_behavior_coupling.r
# ================================================================
# Network-behavior coupling analysis for spatial proteomics
# ================================================================
# Purpose:
#   Link animal-level spatial proteomic network structure to behavioral and
#   physiological stress phenotypes.
#
# Core idea:
#   For each animal, compute the similarity between selected hippocampal
#   region/layer profiles, for example CA1_slm--CA2_slm. Then correlate this
#   animal-level network/edge score with movement, corticosterone reactivity,
#   sucrose preference, and composite stress burden (CombZ).
#
# Important interpretation:
#   These are spatial proteomic similarity/coupling metrics, not anatomical or
#   electrophysiological connectivity metrics.
#
# Required inputs:
#   1) network_spatial_relations_objects.rds from network_spatial_relations.r
#   2) GAMM movement AUC tables from MMMSociability
#   3) E9_Behavior_Data.xlsx with zScore sheet
# ================================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "validation_utils.R"))
MODULE_ID <- "08_behavior_physio_coupling"
SUBSTEP_ID <- "network_behavior_coupling"
CANONICAL_PATHS <- create_module_dirs(MODULE_ID, SUBSTEP_ID)
BEHAVIOR_DATASET <- current_dataset()
assert_dataset_capability(BEHAVIOR_DATASET, "layer", analysis = "network-behavior coupling")
behavior_spatial_unit <- if (BEHAVIOR_DATASET == "neuron_neuropil") "region_layer" else "region"

required_pkgs <- c(
  "dplyr", "tidyr", "stringr", "purrr", "tibble", "readr", "readxl",
  "ggplot2", "broom", "openxlsx", "svglite", "scales", "forcats"
)
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing required R package(s): ", paste(missing, collapse = ", "), ". Install them explicitly before running this script.", call. = FALSE)
invisible(lapply(required_pkgs, library, character.only = TRUE))

# -------------------------------
# 1) Parameters
# -------------------------------
resolve_spatial_rds <- function() {
  override <- Sys.getenv("PROTEOMICS_SPATIAL_NETWORK_OBJECT", unset = "")
  if (nzchar(override)) return(normalizePath(override, winslash = "/", mustWork = FALSE))
  scoped <- path_processed("07_spatial_networks", "network_spatial_relations", BEHAVIOR_DATASET, behavior_spatial_unit, "network_spatial_relations_objects.rds")
  if (file.exists(scoped)) return(scoped)
  path_processed("07_spatial_networks", "network_spatial_relations", "network_spatial_relations_objects.rds")
}

params <- list(
  spatial_rds = resolve_spatial_rds(),

  movement_auc_file = path_external("behavior", "auc_individual_animals_firstChangeActive.csv"),

  movement_auc_all_file = path_external("behavior", "auc_individual_animals_all.csv"),

  behavior_xlsx = path_external("behavior", "E9_Behavior_Data.xlsx"),
  behavior_sheet = "zScore",

  output_dir = CANONICAL_PATHS$reports,

  # Main edges to compute per animal.
  candidate_edges = tibble::tribble(
    ~Source, ~Target,
    "CA1_slm", "CA2_slm",
    "CA1_slm", "CA3_sr",
    "CA2_slm", "CA3_sr",
    "CA1_sr",  "DG_mo",
    "CA1_so",  "DG_mo",
    "CA1_sr",  "CA3_sr",
    "CA2_slm", "DG_po",
    "CA1_so",  "CA1_sr"
  ),

  # Higher values here mean more proteomic similarity/coupling between the two
  # region/layer profiles. Spearman is robust to scale and abundance effects.
  edge_correlation_method = "spearman",

  # Minimum number of proteins needed to estimate an animal-level edge.
  min_proteins_per_edge = 250,

  # Movement variables to correlate after z-scoring against controls.
  behavior_variables = c(
    "Movement_AUC_z_vs_CON",
    "CombZ",
    "delta_cort",
    "sucrose_pref"
  ),

  group_levels = c("CON", "RES", "SUS"),
  sex_levels = c("f", "m")
)

if (is_dry_run()) {
  dry_run_line("Script", "08_behavior_physio_coupling/02_network_behavior_coupling.r")
  dry_run_line("Spatial RDS", params$spatial_rds, if (file.exists(params$spatial_rds)) "PASS" else "FAIL")
  dry_run_line("Movement AUC", params$movement_auc_file, if (file.exists(params$movement_auc_file)) "PASS" else "FAIL")
  dry_run_line("Movement AUC all", params$movement_auc_all_file, if (file.exists(params$movement_auc_all_file)) "PASS" else "WARN")
  dry_run_line("Behavior workbook", params$behavior_xlsx, if (file.exists(params$behavior_xlsx)) "PASS" else "FAIL")
  dry_run_line("Output folders", paste(unlist(CANONICAL_PATHS), collapse = "; "))
  quit(status = if (file.exists(params$spatial_rds) && file.exists(params$movement_auc_file) && file.exists(params$behavior_xlsx)) 0 else 1, save = "no")
}
if (!file.exists(params$spatial_rds)) stop("spatial_rds not found: ", params$spatial_rds, call. = FALSE)
if (!file.exists(params$movement_auc_file)) stop("movement_auc_file not found: ", params$movement_auc_file, call. = FALSE)
if (!file.exists(params$behavior_xlsx)) stop("behavior_xlsx not found: ", params$behavior_xlsx, call. = FALSE)

# -------------------------------
# 2) Helpers
# -------------------------------
message2 <- function(...) message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ...)

safe_name <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_replace_all("[^A-Za-z0-9_\\-]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_|_$", "")
}

make_dirs <- function(base_dir) {
  dirs <- list(
    base = base_dir,
    tables = CANONICAL_PATHS$tables,
    figures = CANONICAL_PATHS$figures,
    logs = CANONICAL_PATHS$logs
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  dirs
}

theme_publication_coupling <- function(base_size = 7) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(family = "sans", colour = "black"),
      axis.text = ggplot2::element_text(colour = "black", size = base_size - 1),
      axis.title = ggplot2::element_text(colour = "black", face = "bold", size = base_size),
      axis.line = ggplot2::element_line(linewidth = 0.25, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.25, colour = "black"),
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = base_size),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = base_size + 1),
      plot.subtitle = ggplot2::element_text(hjust = 0, colour = "grey25", size = base_size - 1),
      plot.margin = ggplot2::margin(5, 6, 5, 6)
    )
}

group_colors <- c(
  "CON" = "#3E3C6F",
  "RES" = "#9E9A92",
  "SUS" = "#D7303F"
)

group_fills <- scales::alpha(group_colors, 0.35)

normalize_animal_id <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x <- stringr::str_extract(x, "A[0-9]{3,4}|[0-9]{3,4}")
  x <- ifelse(is.na(x), NA_character_, x)
  x <- ifelse(stringr::str_detect(x, "^A"), x, paste0("A", stringr::str_pad(x, 4, pad = "0")))
  x
}

first_existing_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (!length(hit)) return(rep(NA_character_, nrow(df)))
  as.character(df[[hit[[1]]]])
}

safe_cor <- function(data, x, y, method = "pearson") {
  d <- data %>%
    dplyr::filter(is.finite(.data[[x]]), is.finite(.data[[y]]))

  if (
    nrow(d) < 4 ||
    stats::sd(d[[x]], na.rm = TRUE) == 0 ||
    stats::sd(d[[y]], na.rm = TRUE) == 0
  ) {
    return(tibble::tibble(
      estimate = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      p.value = NA_real_,
      n = nrow(d)
    ))
  }

  ct <- suppressWarnings(stats::cor.test(d[[x]], d[[y]], method = method))
  tibble::tibble(
    estimate = unname(ct$estimate),
    conf.low = if (!is.null(ct$conf.int)) unname(ct$conf.int[1]) else NA_real_,
    conf.high = if (!is.null(ct$conf.int)) unname(ct$conf.int[2]) else NA_real_,
    p.value = ct$p.value,
    n = nrow(d)
  )
}

fit_lm_safe <- function(data, formula, model_name) {
  vars <- all.vars(formula)
  d <- data %>% tidyr::drop_na(dplyr::any_of(vars))

  if (nrow(d) < 6) {
    return(tibble::tibble(
      model = model_name,
      term = NA_character_,
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      n = nrow(d),
      r.squared = NA_real_,
      adj.r.squared = NA_real_
    ))
  }

  fit <- tryCatch(stats::lm(formula, data = d), error = function(e) NULL)
  if (is.null(fit)) {
    return(tibble::tibble(
      model = model_name,
      term = NA_character_,
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      n = nrow(d),
      r.squared = NA_real_,
      adj.r.squared = NA_real_
    ))
  }

  gl <- broom::glance(fit)
  broom::tidy(fit, conf.int = TRUE) %>%
    dplyr::mutate(
      model = model_name,
      n = nrow(d),
      r.squared = gl$r.squared,
      adj.r.squared = gl$adj.r.squared,
      .before = 1
    )
}

format_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}

# -------------------------------
# 3) Load spatial proteomics object
# -------------------------------
load_spatial_object <- function(params) {
  if (!file.exists(params$spatial_rds)) {
    stop("spatial_rds not found: ", params$spatial_rds)
  }
  obj <- readRDS(params$spatial_rds)

  if (!all(c("expression_matrix", "sample_metadata") %in% names(obj))) {
    stop("spatial_rds must contain expression_matrix and sample_metadata.")
  }

  expr <- obj$expression_matrix
  expr <- as.matrix(expr)
  storage.mode(expr) <- "numeric"
  expr[!is.finite(expr)] <- NA_real_

  sample_md_raw <- as.data.frame(obj$sample_metadata)
  animal_id_raw <- first_existing_col(sample_md_raw, c("AnimalID", "AnimalNum", "animal_id", "animal"))
  if (all(is.na(animal_id_raw)) && "SampleColumn" %in% names(sample_md_raw)) {
    animal_id_raw <- as.character(sample_md_raw$SampleColumn)
  }
  sample_column_raw <- first_existing_col(sample_md_raw, c("SampleColumn", "sample_id", "sample", "SampleID"))

  sample_md <- sample_md_raw %>%
    dplyr::mutate(
      AnimalID = .env$animal_id_raw,
      AnimalID = dplyr::if_else(
        is.na(.data$AnimalID) | .data$AnimalID == "",
        stringr::str_extract(.env$sample_column_raw, "A[0-9]{3,4}"),
        .data$AnimalID
      ),
      AnimalID = normalize_animal_id(.data$AnimalID),
      ExpGroup = toupper(as.character(.data$ExpGroup)),
      ExpGroup = factor(.data$ExpGroup, levels = params$group_levels),
      RegionLayer = as.character(.data$RegionLayer),
      SampleColumn = as.character(.env$sample_column_raw)
    ) %>%
    dplyr::filter(
      .data$SampleColumn %in% colnames(expr),
      !is.na(.data$AnimalID),
      !is.na(.data$RegionLayer),
      !is.na(.data$ExpGroup)
    )

  list(expr = expr, sample_md = sample_md)
}

# -------------------------------
# 4) Animal-level edge/network metrics
# -------------------------------
compute_animal_edge_scores <- function(expr, sample_md, candidate_edges, params) {
  animals <- sort(unique(sample_md$AnimalID))
  out <- list()

  for (animal in animals) {
    md_a <- sample_md %>% dplyr::filter(.data$AnimalID == .env$animal)
    group_a <- unique(as.character(md_a$ExpGroup))
    group_a <- group_a[!is.na(group_a)][1]

    for (i in seq_len(nrow(candidate_edges))) {
      source_rl <- candidate_edges$Source[i]
      target_rl <- candidate_edges$Target[i]

      source_cols <- md_a %>%
        dplyr::filter(.data$RegionLayer == .env$source_rl) %>%
        dplyr::pull(.data$SampleColumn)
      target_cols <- md_a %>%
        dplyr::filter(.data$RegionLayer == .env$target_rl) %>%
        dplyr::pull(.data$SampleColumn)

      source_cols <- source_cols[source_cols %in% colnames(expr)]
      target_cols <- target_cols[target_cols %in% colnames(expr)]

      if (length(source_cols) == 0 || length(target_cols) == 0) {
        out[[length(out) + 1]] <- tibble::tibble(
          AnimalID = animal,
          ExpGroup = group_a,
          Source = source_rl,
          Target = target_rl,
          Edge = paste(source_rl, target_rl, sep = " - "),
          EdgeR = NA_real_,
          EdgeAbsR = NA_real_,
          EdgeDistance = NA_real_,
          NProteins = 0L
        )
        next
      }

      source_profile <- rowMeans(expr[, source_cols, drop = FALSE], na.rm = TRUE)
      target_profile <- rowMeans(expr[, target_cols, drop = FALSE], na.rm = TRUE)
      keep <- is.finite(source_profile) & is.finite(target_profile)

      if (sum(keep) < params$min_proteins_per_edge ||
          stats::sd(source_profile[keep], na.rm = TRUE) == 0 ||
          stats::sd(target_profile[keep], na.rm = TRUE) == 0) {
        edge_r <- NA_real_
      } else {
        edge_r <- suppressWarnings(stats::cor(
          source_profile[keep],
          target_profile[keep],
          method = params$edge_correlation_method,
          use = "complete.obs"
        ))
      }

      out[[length(out) + 1]] <- tibble::tibble(
        AnimalID = animal,
        ExpGroup = group_a,
        Source = source_rl,
        Target = target_rl,
        Edge = paste(source_rl, target_rl, sep = " - "),
        EdgeR = as.numeric(edge_r),
        EdgeAbsR = abs(as.numeric(edge_r)),
        EdgeDistance = 1 - as.numeric(edge_r),
        NProteins = sum(keep)
      )
    }
  }

  dplyr::bind_rows(out) %>%
    dplyr::mutate(
      ExpGroup = factor(.data$ExpGroup, levels = params$group_levels),
      Edge = factor(.data$Edge, levels = unique(paste(candidate_edges$Source, candidate_edges$Target, sep = " - ")))
    )
}

compute_animal_global_network_metrics <- function(edge_scores) {
  edge_scores %>%
    dplyr::group_by(.data$AnimalID, .data$ExpGroup) %>%
    dplyr::summarise(
      MeanEdgeR = mean(.data$EdgeR, na.rm = TRUE),
      MeanAbsEdgeR = mean(.data$EdgeAbsR, na.rm = TRUE),
      SLM_CA1_CA2_EdgeR = .data$EdgeR[.data$Source == "CA1_slm" & .data$Target == "CA2_slm"][1],
      NEdgesEstimated = sum(is.finite(.data$EdgeR)),
      .groups = "drop"
    )
}

# -------------------------------
# 5) Load behavior/physiology
# -------------------------------
load_physiology <- function(params) {
  if (!file.exists(params$behavior_xlsx)) {
    warning("behavior_xlsx not found: ", params$behavior_xlsx)
    return(tibble::tibble())
  }

  readxl::read_excel(params$behavior_xlsx, sheet = params$behavior_sheet) %>%
    dplyr::rename(
      AnimalRaw = dplyr::any_of(c("ID", "AnimalNum", "AnimalID")),
      PhysioGroup = dplyr::any_of(c("Group", "ExpGroup"))
    ) %>%
    dplyr::mutate(
      AnimalID = normalize_animal_id(.data$AnimalRaw),
      Sex = as.character(.data$Sex),
      Sex = factor(.data$Sex, levels = params$sex_levels),
      PhysioGroup = toupper(as.character(.data$PhysioGroup)),
      CombZ = suppressWarnings(as.numeric(.data$CombZ)),
      delta_cort = suppressWarnings(as.numeric(.data$delta_cort)),
      sucrose_pref = suppressWarnings(as.numeric(.data$sucrose_pref))
    ) %>%
    dplyr::select(.data$AnimalID, .data$Sex, .data$PhysioGroup, .data$CombZ, .data$delta_cort, .data$sucrose_pref) %>%
    dplyr::distinct(.data$AnimalID, .keep_all = TRUE)
}

load_movement_auc <- function(auc_file, params, analysis_label) {
  if (is.null(auc_file) || !file.exists(auc_file)) {
    warning("Movement AUC file not found: ", auc_file)
    return(tibble::tibble())
  }

  auc_data <- readr::read_csv(auc_file, show_col_types = FALSE)

  needed <- c("AnimalNum", "Group", "Sex", "Phase", "metric", "AUC_norm")
  missing_needed <- setdiff(needed, names(auc_data))
  if (length(missing_needed) > 0) {
    warning("Movement AUC file is missing columns: ", paste(missing_needed, collapse = ", "))
    return(tibble::tibble())
  }

  auc_movement <- auc_data %>%
    dplyr::mutate(
      Analysis = analysis_label,
      AnimalID = normalize_animal_id(.data$AnimalNum),
      Batch = if ("Batch" %in% names(.)) as.character(.data$Batch) else NA_character_,
      Group = toupper(as.character(.data$Group)),
      Group = factor(.data$Group, levels = params$group_levels),
      Sex = as.character(.data$Sex),
      Sex = factor(.data$Sex, levels = params$sex_levels),
      Phase = as.character(.data$Phase),
      Phase = factor(.data$Phase, levels = c("Active", "Inactive")),
      Change = if ("Change" %in% names(.)) as.character(.data$Change) else NA_character_,
      window = if ("window" %in% names(.)) as.character(.data$window) else NA_character_,
      metric = as.character(.data$metric),
      MeanPredictedMovement = suppressWarnings(as.numeric(.data$AUC_norm))
    ) %>%
    dplyr::filter(
      .data$metric == "Movement",
      is.finite(.data$MeanPredictedMovement),
      !is.na(.data$AnimalID),
      !is.na(.data$Group),
      !is.na(.data$Sex),
      !is.na(.data$Phase)
    )

  # Batch-aware CON reference exactly mirrors the MMMSociability logic, with a
  # global CON fallback when the batch reference is too small.
  control_stats_batch <- auc_movement %>%
    dplyr::filter(.data$Group == "CON") %>%
    dplyr::group_by(.data$Batch, .data$Sex, .data$Phase, .data$Change, .data$window) %>%
    dplyr::summarise(
      con_mean_batch = mean(.data$MeanPredictedMovement, na.rm = TRUE),
      con_sd_batch = stats::sd(.data$MeanPredictedMovement, na.rm = TRUE),
      con_n_batch = sum(is.finite(.data$MeanPredictedMovement)),
      .groups = "drop"
    )

  control_stats_global <- auc_movement %>%
    dplyr::filter(.data$Group == "CON") %>%
    dplyr::group_by(.data$Sex, .data$Phase, .data$Change, .data$window) %>%
    dplyr::summarise(
      con_mean_global = mean(.data$MeanPredictedMovement, na.rm = TRUE),
      con_sd_global = stats::sd(.data$MeanPredictedMovement, na.rm = TRUE),
      con_n_global = sum(is.finite(.data$MeanPredictedMovement)),
      .groups = "drop"
    )

  auc_movement %>%
    dplyr::left_join(control_stats_batch, by = c("Batch", "Sex", "Phase", "Change", "window")) %>%
    dplyr::left_join(control_stats_global, by = c("Sex", "Phase", "Change", "window")) %>%
    dplyr::mutate(
      use_batch_reference = !is.na(.data$con_sd_batch) & .data$con_sd_batch > 0 & .data$con_n_batch >= 3,
      CON_mean_used = dplyr::if_else(.data$use_batch_reference, .data$con_mean_batch, .data$con_mean_global),
      CON_sd_used = dplyr::if_else(.data$use_batch_reference, .data$con_sd_batch, .data$con_sd_global),
      Movement_AUC_z_vs_CON = dplyr::if_else(
        is.na(.data$CON_sd_used) | .data$CON_sd_used == 0,
        NA_real_,
        (.data$MeanPredictedMovement - .data$CON_mean_used) / .data$CON_sd_used
      )
    ) %>%
    dplyr::filter(is.finite(.data$Movement_AUC_z_vs_CON)) %>%
    dplyr::select(
      .data$Analysis, .data$AnimalID, .data$Group, .data$Sex, .data$Phase,
      .data$Change, .data$window, .data$Batch,
      .data$MeanPredictedMovement, .data$Movement_AUC_z_vs_CON
    )
}

# -------------------------------
# 6) Coupling statistics
# -------------------------------
make_long_coupling_table <- function(edge_scores, global_metrics, movement_tbl, physio_tbl) {
  behavior_tbl <- movement_tbl %>%
    dplyr::left_join(physio_tbl, by = c("AnimalID", "Sex")) %>%
    dplyr::mutate(
      Group = dplyr::coalesce(as.character(.data$Group), as.character(.data$PhysioGroup)),
      Group = factor(.data$Group, levels = params$group_levels)
    )

  edge_joined <- edge_scores %>%
    dplyr::left_join(behavior_tbl, by = "AnimalID") %>%
    dplyr::mutate(
      Group = dplyr::coalesce(as.character(.data$Group), as.character(.data$ExpGroup)),
      Group = factor(.data$Group, levels = params$group_levels)
    )

  global_joined <- global_metrics %>%
    dplyr::left_join(behavior_tbl, by = "AnimalID") %>%
    dplyr::mutate(
      Group = dplyr::coalesce(as.character(.data$Group), as.character(.data$ExpGroup)),
      Group = factor(.data$Group, levels = params$group_levels)
    )

  list(edge_joined = edge_joined, global_joined = global_joined)
}

run_edge_behavior_correlations <- function(edge_joined, behavior_variables) {
  edge_joined %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(behavior_variables[behavior_variables %in% names(edge_joined)]),
      names_to = "Outcome",
      values_to = "OutcomeValue"
    ) %>%
    dplyr::group_by(.data$Analysis, .data$Phase, .data$Change, .data$window, .data$Edge, .data$Outcome) %>%
    dplyr::group_modify(~ safe_cor(.x, "EdgeR", "OutcomeValue", method = "pearson")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      p.adj_BH_within_outcome = ave(.data$p.value, .data$Outcome, FUN = function(p) p.adjust(p, method = "BH")),
      p.adj_BH_all_edge_phenotype_tests = p.adjust(.data$p.value, method = "BH"),
      interpretation_strength = vapply(seq_along(.data$p.value), function(i) {
        interpretation_strength(
          fdr = p.adj_BH_all_edge_phenotype_tests[[i]],
          effect_size = estimate[[i]],
          n = n[[i]]
        )
      }, character(1)),
      limitations = dplyr::if_else(.data$n < 6, "n < 6; exploratory only", NA_character_)
    ) %>%
    dplyr::arrange(.data$Outcome, .data$p.value)
}

run_sex_stratified_correlations <- function(edge_joined, behavior_variables) {
  if (!"Sex" %in% names(edge_joined)) return(tibble::tibble())
  edge_joined %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(behavior_variables[behavior_variables %in% names(edge_joined)]),
      names_to = "Outcome",
      values_to = "OutcomeValue"
    ) %>%
    dplyr::filter(!is.na(.data$Sex)) %>%
    dplyr::group_by(.data$Sex, .data$Analysis, .data$Phase, .data$Change, .data$window, .data$Edge, .data$Outcome) %>%
    dplyr::group_modify(~ safe_cor(.x, "EdgeR", "OutcomeValue", method = "pearson")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      p.adj_BH_all_sex_stratified_tests = p.adjust(.data$p.value, method = "BH"),
      interpretation_strength = vapply(seq_along(.data$p.value), function(i) {
        interpretation_strength(fdr = p.adj_BH_all_sex_stratified_tests[[i]], effect_size = estimate[[i]], n = n[[i]])
      }, character(1)),
      limitations = dplyr::if_else(.data$n < 6, "n < 6; exploratory only", NA_character_)
    )
}

run_group_adjusted_models <- function(edge_joined) {
  d <- edge_joined %>%
    dplyr::filter(is.finite(.data$EdgeR)) %>%
    dplyr::mutate(Group = droplevels(.data$Group))

  outcomes <- c("Movement_AUC_z_vs_CON", "CombZ", "delta_cort", "sucrose_pref")
  outcomes <- outcomes[outcomes %in% names(d)]

  purrr::map_dfr(outcomes, function(outcome) {
    d %>%
      dplyr::filter(is.finite(.data[[outcome]])) %>%
      dplyr::group_by(.data$Analysis, .data$Phase, .data$Change, .data$window, .data$Edge) %>%
      dplyr::group_split() %>%
      purrr::map_dfr(function(dd) {
        if (nrow(dd) < 6) return(tibble::tibble())
        edge_i <- as.character(unique(dd$Edge))[1]
        analysis_i <- as.character(unique(dd$Analysis))[1]
        phase_i <- as.character(unique(dd$Phase))[1]
        change_i <- as.character(unique(dd$Change))[1]
        window_i <- as.character(unique(dd$window))[1]

        covars <- c()
        if ("Sex" %in% names(dd) && length(unique(stats::na.omit(dd$Sex))) >= 2) covars <- c(covars, "Sex")
        if (length(unique(stats::na.omit(dd$Group))) >= 2) covars <- c(covars, "Group")
        formula_use <- if (length(covars)) {
          stats::as.formula(paste0(outcome, " ~ EdgeR + ", paste(covars, collapse = " + ")))
        } else {
          stats::as.formula(paste0(outcome, " ~ EdgeR"))
        }

        fit_lm_safe(
          dd,
          formula_use,
          paste(analysis_i, phase_i, change_i, window_i, edge_i, outcome, sep = "__")
        ) %>%
          dplyr::mutate(
            Analysis = analysis_i,
            Phase = phase_i,
            Change = change_i,
            window = window_i,
            Edge = edge_i,
            Outcome = outcome,
            fdr_all_model_terms = p.adjust(.data$p.value, method = "BH"),
            limitations = dplyr::if_else(.data$n < 6, "n < 6; exploratory only", NA_character_),
            .before = 1
          )
      })
  })
}

# -------------------------------
# 7) Figures
# -------------------------------
plot_edge_outcome_scatter <- function(data, edge_name, outcome, outfile, title = NULL) {
  d <- data %>%
    dplyr::filter(.data$Edge == edge_name, is.finite(.data$EdgeR), is.finite(.data[[outcome]]))

  if (nrow(d) < 4) return(invisible(NULL))

  stat <- safe_cor(d, "EdgeR", outcome, method = "pearson")
  subtitle <- if (is.finite(stat$estimate)) {
    sprintf("r = %.2f, p = %s, n = %d", stat$estimate, format_p(stat$p.value), stat$n)
  } else {
    sprintf("n = %d; correlation not estimated", stat$n)
  }

  if (is.null(title)) title <- paste(edge_name, "vs", outcome)

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$EdgeR, y = .data[[outcome]])) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, colour = "black", fill = "grey88", linewidth = 0.35, alpha = 0.45) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$Group), size = 2.1, alpha = 0.9) +
    ggplot2::scale_colour_manual(values = group_colors, drop = FALSE) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Animal-level edge similarity (Spearman rho)",
      y = outcome
    ) +
    theme_publication_coupling(base_size = 7)

  ggplot2::ggsave(outfile, p, width = 80, height = 70, units = "mm", device = svglite::svglite)
  invisible(p)
}

plot_correlation_forest <- function(cor_tbl, outfile, outcome_filter = NULL, top_n = 20) {
  dat <- cor_tbl %>%
    dplyr::filter(is.finite(.data$estimate))

  if (!is.null(outcome_filter)) {
    dat <- dat %>% dplyr::filter(.data$Outcome %in% outcome_filter)
  }

  dat <- dat %>%
    dplyr::arrange(.data$p.value) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::mutate(
      Label = paste(.data$Edge, .data$Outcome, sep = " | "),
      Label = forcats::fct_reorder(.data$Label, .data$estimate)
    )

  if (nrow(dat) == 0) return(invisible(NULL))

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$estimate, y = .data$Label)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, linetype = "dashed", colour = "grey45") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data$conf.low, xmax = .data$conf.high), height = 0, linewidth = 0.35) +
    ggplot2::geom_point(shape = 21, fill = "white", colour = "black", size = 2, stroke = 0.25) +
    ggplot2::facet_wrap(~ Analysis, scales = "free_y") +
    ggplot2::labs(
      title = "Edge-behavior coupling",
      x = "Pearson r",
      y = NULL
    ) +
    theme_publication_coupling(base_size = 7)

  ggplot2::ggsave(outfile, p, width = 110, height = 120, units = "mm", device = svglite::svglite)
  invisible(p)
}

plot_edge_group_distribution <- function(edge_scores, outfile) {
  dat <- edge_scores %>%
    dplyr::filter(is.finite(.data$EdgeR))

  if (nrow(dat) == 0) return(invisible(NULL))

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$ExpGroup, y = .data$EdgeR)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.25) +
    ggplot2::geom_violin(ggplot2::aes(fill = .data$ExpGroup), width = 0.72, alpha = 0.45, colour = NA, trim = FALSE) +
    ggplot2::geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", linewidth = 0.25) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$ExpGroup), position = ggplot2::position_jitter(width = 0.07, height = 0), size = 1.7, alpha = 0.85) +
    ggplot2::facet_wrap(~ Edge, scales = "free_y") +
    ggplot2::scale_fill_manual(values = group_fills, drop = FALSE) +
    ggplot2::scale_colour_manual(values = group_colors, drop = FALSE) +
    ggplot2::labs(
      title = "Animal-level spatial proteomic edge similarity",
      x = NULL,
      y = "Spearman rho"
    ) +
    theme_publication_coupling(base_size = 7) +
    ggplot2::theme(legend.position = "none")

  ggplot2::ggsave(outfile, p, width = 180, height = 120, units = "mm", device = svglite::svglite)
  invisible(p)
}

# -------------------------------
# 8) Main
# -------------------------------
dirs <- make_dirs(params$output_dir)
write_session_info(file.path(dirs$logs, "sessionInfo.txt"))
log_file <- file.path(dirs$logs, paste0("network_behavior_coupling_run_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

message2("Starting network-behavior coupling analysis")
print(params)

loaded <- load_spatial_object(params)
expr <- loaded$expr
sample_md <- loaded$sample_md

message2("Proteomics matrix: ", nrow(expr), " proteins x ", ncol(expr), " samples")
message2("Animals with proteomics metadata: ", length(unique(sample_md$AnimalID)))
message2("Region/layer units: ", paste(sort(unique(sample_md$RegionLayer)), collapse = ", "))

edge_scores <- compute_animal_edge_scores(expr, sample_md, params$candidate_edges, params)
global_metrics <- compute_animal_global_network_metrics(edge_scores)

readr::write_csv(edge_scores, file.path(dirs$tables, "animal_level_candidate_edge_scores.csv"))
readr::write_csv(global_metrics, file.path(dirs$tables, "animal_level_global_network_metrics.csv"))

physio_tbl <- load_physiology(params)
readr::write_csv(physio_tbl, file.path(dirs$tables, "physiology_traits_loaded.csv"))

movement_first <- load_movement_auc(params$movement_auc_file, params, "firstChangeActive")
movement_all <- load_movement_auc(params$movement_auc_all_file, params, "allPhases")
movement_tbl <- dplyr::bind_rows(movement_first, movement_all)

if (nrow(movement_tbl) == 0 && nrow(physio_tbl) == 0) {
  stop("No movement or physiology data could be loaded. Check params$movement_auc_file and params$behavior_xlsx.")
}

readr::write_csv(movement_tbl, file.path(dirs$tables, "movement_auc_z_loaded.csv"))

joined <- make_long_coupling_table(edge_scores, global_metrics, movement_tbl, physio_tbl)
edge_joined <- joined$edge_joined
global_joined <- joined$global_joined

proteomics_animals <- sort(unique(edge_scores$AnimalID))
behavior_animals <- sort(unique(c(movement_tbl$AnimalID, physio_tbl$AnimalID)))
joined_animals <- sort(unique(edge_joined$AnimalID[!is.na(edge_joined$Group)]))
join_diag <- dplyr::bind_rows(
  data.frame(metric = "proteomics_animals_before_join", value = length(proteomics_animals)),
  data.frame(metric = "behavior_animals_before_join", value = length(behavior_animals)),
  data.frame(metric = "animals_after_join", value = length(joined_animals)),
  data.frame(metric = "animals_lost_from_proteomics", value = length(setdiff(proteomics_animals, behavior_animals))),
  data.frame(metric = "animals_lost_from_behavior", value = length(setdiff(behavior_animals, proteomics_animals)))
)
readr::write_csv(join_diag, file.path(dirs$tables, "join_diagnostics_summary.csv"))
readr::write_csv(data.frame(source = "proteomics_only", AnimalID = setdiff(proteomics_animals, behavior_animals)), file.path(dirs$tables, "join_diagnostics_proteomics_only.csv"))
readr::write_csv(data.frame(source = "behavior_only", AnimalID = setdiff(behavior_animals, proteomics_animals)), file.path(dirs$tables, "join_diagnostics_behavior_only.csv"))
coverage_cols <- intersect(c("Sex", "Group", "Phase", "window"), names(edge_joined))
if (length(coverage_cols) > 0) {
  coverage <- edge_joined %>% dplyr::count(dplyr::across(dplyr::all_of(coverage_cols)), name = "n")
  readr::write_csv(coverage, file.path(dirs$tables, "join_diagnostics_coverage.csv"))
}

readr::write_csv(edge_joined, file.path(dirs$tables, "merged_edge_behavior_long.csv"))
readr::write_csv(global_joined, file.path(dirs$tables, "merged_global_network_behavior.csv"))

qc <- edge_joined %>%
  dplyr::count(.data$Analysis, .data$Phase, .data$Change, .data$window, .data$Edge, .data$Group, name = "n")
readr::write_csv(qc, file.path(dirs$tables, "qc_counts_edge_behavior.csv"))

cor_tbl <- run_edge_behavior_correlations(edge_joined, params$behavior_variables)
readr::write_csv(cor_tbl, file.path(dirs$tables, "edge_behavior_correlations.csv"))
readr::write_csv(cor_tbl, file.path(dirs$tables, "edge_behavior_correlations_fdr_all_tests.csv"))

sex_cor_tbl <- run_sex_stratified_correlations(edge_joined, params$behavior_variables)
readr::write_csv(sex_cor_tbl, file.path(dirs$tables, "edge_behavior_correlations_sex_stratified.csv"))

model_tbl <- run_group_adjusted_models(edge_joined)
if (nrow(model_tbl) && "p.value" %in% names(model_tbl)) {
  model_tbl <- model_tbl %>%
    dplyr::mutate(
      fdr_all_model_terms = p.adjust(.data$p.value, method = "BH"),
      interpretation_strength = vapply(seq_along(.data$p.value), function(i) {
        interpretation_strength(fdr = fdr_all_model_terms[[i]], effect_size = estimate[[i]], n = n[[i]])
      }, character(1)),
      limitations = dplyr::if_else(.data$n < 6, "n < 6; exploratory only", .data$limitations)
    )
}
readr::write_csv(model_tbl, file.path(dirs$tables, "edge_behavior_group_adjusted_models.csv"))
readr::write_csv(model_tbl, file.path(dirs$tables, "edge_behavior_model_summaries_fdr.csv"))

figure_ready_coupling <- cor_tbl %>%
  dplyr::filter(is.finite(.data$estimate)) %>%
  dplyr::arrange(.data$p.adj_BH_all_edge_phenotype_tests, dplyr::desc(abs(.data$estimate))) %>%
  dplyr::transmute(
    Analysis,
    Phase,
    Change,
    window,
    Edge,
    Outcome,
    estimate,
    p.value,
    fdr = .data$p.adj_BH_all_edge_phenotype_tests,
    n,
    interpretation_strength,
    limitations
  )
readr::write_csv(figure_ready_coupling, file.path(dirs$tables, "edge_behavior_figure_ready_table.csv"))

# Focused tables for the biologically central edge.
central_edge <- "CA1_slm - CA2_slm"
central_tbl <- edge_joined %>%
  dplyr::filter(.data$Edge == central_edge)
readr::write_csv(central_tbl, file.path(dirs$tables, "central_edge_CA1_slm_CA2_slm_behavior_table.csv"))

central_cor <- cor_tbl %>%
  dplyr::filter(.data$Edge == central_edge)
readr::write_csv(central_cor, file.path(dirs$tables, "central_edge_CA1_slm_CA2_slm_correlations.csv"))

# Figures
plot_edge_group_distribution(
  edge_scores,
  file.path(dirs$figures, "animal_level_edge_similarity_by_group.svg")
)

plot_correlation_forest(
  cor_tbl,
  file.path(dirs$figures, "edge_behavior_correlation_forest_top_hits.svg"),
  top_n = 24
)

# Candidate scatter plots for the central edge.
for (outcome in intersect(params$behavior_variables, names(edge_joined))) {
  plot_edge_outcome_scatter(
    central_tbl,
    central_edge,
    outcome,
    file.path(dirs$figures, paste0("central_edge_CA1_slm_CA2_slm_vs_", safe_name(outcome), ".svg")),
    title = paste0("CA1_slm-CA2_slm coupling vs ", outcome)
  )
}

saveRDS(
  list(
    params = params,
    edge_scores = edge_scores,
    global_metrics = global_metrics,
    physio_tbl = physio_tbl,
    movement_tbl = movement_tbl,
    edge_joined = edge_joined,
    global_joined = global_joined,
    cor_tbl = cor_tbl,
    model_tbl = model_tbl,
    sessionInfo = sessionInfo()
  ),
  file.path(dirs$logs, "network_behavior_coupling_objects.rds")
)
write_run_manifest(
  file.path(dirs$logs, "run_manifest.yml"),
  inputs = list(
    spatial_rds = params$spatial_rds,
    movement_auc_file = params$movement_auc_file,
    movement_auc_all_file = params$movement_auc_all_file,
    behavior_xlsx = params$behavior_xlsx
  ),
  outputs = list(tables = dirs$tables, figures = dirs$figures, logs = dirs$logs),
  parameters = params,
  notes = "Join diagnostics include animal counts, animals lost from each side, and sex/group/window coverage."
)

message2("Finished network-behavior coupling analysis")
message2("Output directory: ", params$output_dir)
message2("Interpretation caution: animal-level edge scores are molecular similarity metrics, not direct connectivity.")

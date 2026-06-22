#!/usr/bin/env Rscript
# ================================================================
# Script: 06_modules_WGCNA/09_microglia_neuropil_independence.R
# Stage: modules_downstream
# Scope: dataset_specific
# Consumes: microglia and neuron_neuropil WGCNA state/metadata plus optional marker traits.
# Produces: additive microglia neuropil-independence sensitivity and reviewer audit tables.
# Notes: Does not overwrite primary WGCNA group-effect outputs.
# ================================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "wgcna_downstream_utils.R"))
source(repo_path("R", "schema_validation.R"))

SCRIPT_ID <- "06_modules_WGCNA/09_microglia_neuropil_independence.R"
required_pkgs <- c("dplyr", "tidyr", "tibble", "readr", "yaml")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) && !is_dry_run()) stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
if (!length(missing_pkgs)) suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = "") {
  hit <- which(args == flag)
  if (!length(hit) || hit[[1]] == length(args)) return(default)
  args[[hit[[1]] + 1L]]
}
DATASET <- validate_dataset(value_after("--dataset", Sys.getenv("PROTEOMICS_DATASET", unset = "microglia")), source = "--dataset")
if (!identical(DATASET, "microglia") && !is_dry_run()) {
  stop("This sensitivity analysis is microglia-only. Use --dataset microglia.", call. = FALSE)
}

PATHS <- create_module_dirs("06_modules_WGCNA", file.path("microglia_neuropil_independence", "microglia"))
FILES_MICRO <- resolve_wgcna_files("microglia")
FILES_NEURO <- resolve_wgcna_files("neuron_neuropil")
CONFIG_FILE <- repo_path("config", "microglia_neuropil_independence.yml")
inputs <- list(
  independence_config = CONFIG_FILE,
  microglia_state = FILES_MICRO$state,
  neuron_neuropil_state = FILES_NEURO$state,
  microglia_definitions = FILES_MICRO$definitions,
  neuron_neuropil_definitions = FILES_NEURO$definitions,
  microglia_marker_traits = FILES_MICRO$marker_traits,
  neuron_neuropil_marker_traits = FILES_NEURO$marker_traits,
  microglia_interpretable = path_results("tables", "06_modules_WGCNA", "interpretable_summary", "microglia", "WGCNA_module_group_effects_interpretable.csv")
)

if (is_dry_run()) {
  dry_run_line("Script", SCRIPT_ID)
  dry_run_line("Dataset", DATASET)
  dry_run_line("Output tables", PATHS$tables)
  dry_run_line("Reviewer claim gate audit", path_results("reviewer_audit", "microglia_neuropil_independence_claim_gate.csv"))
  dry_run_line("Reviewer covariate selection audit", path_results("reviewer_audit", "microglia_neuropil_covariate_selection_audit.csv"))
  dry_run_line("Reviewer endpoint scope audit", path_results("reviewer_audit", "microglia_neuropil_independence_endpoint_scope_audit.csv"))
  for (nm in names(inputs)) dry_run_line(nm, inputs[[nm]], if (file.exists(inputs[[nm]])) "PASS" else "WARN")
  quit(status = 0, save = "no")
}

read_csv_optional2 <- function(path) safe_read_csv(path) %||% data.frame()

clean_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(trimws(x))] <- NA_character_
  x
}

endpoint_scope_for <- function(endpoint_type) {
  dplyr::case_when(
    endpoint_type == "module_eigengene" ~ "module",
    endpoint_type == "targeted_signature_score" ~ "targeted_signature",
    endpoint_type %in% c("supermodule_eigengene", "supermodule_score") ~ "supermodule",
    TRUE ~ "unavailable"
  )
}

as_int <- function(x, default) {
  y <- suppressWarnings(as.integer(x))
  ifelse(is.na(y), default, y)
}

has_repeats <- function(dat) {
  "AnimalID" %in% names(dat) &&
    any(!is.na(dat$AnimalID) & nzchar(as.character(dat$AnimalID))) &&
    any(table(dat$AnimalID[!is.na(dat$AnimalID) & nzchar(as.character(dat$AnimalID))]) > 1L)
}

valid_covariates <- function(dat, covars, protect = character()) {
  missing <- setdiff(covars, names(dat))
  present <- intersect(covars, names(dat))
  nonvarying <- present[vapply(present, function(v) dplyr::n_distinct(dat[[v]][!is.na(dat[[v]])]) < 2L, logical(1))]
  dropped <- unique(c(missing, setdiff(nonvarying, protect)))
  list(keep = setdiff(present, dropped), dropped = dropped)
}

fit_score_model <- function(dat, rhs_terms, use_random_animal = TRUE) {
  repeated <- has_repeats(dat)
  use_lmer <- isTRUE(use_random_animal) && repeated
  warning_text <- character()
  if (use_lmer && !requireNamespace("lmerTest", quietly = TRUE)) {
    warning_text <- c(warning_text, "AnimalID repeats detected but lmerTest is unavailable; used lm")
    use_lmer <- FALSE
  }
  rhs <- paste(rhs_terms, collapse = " + ")
  formula_requested <- paste0("score ~ ", rhs, if (use_lmer) " + (1 | AnimalID)" else "")
  fit <- tryCatch(
    if (use_lmer) lmerTest::lmer(stats::as.formula(formula_requested), data = dat, REML = FALSE) else stats::lm(stats::as.formula(formula_requested), data = dat),
    warning = function(w) {
      warning_text <<- c(warning_text, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(fit = NULL, model_type = if (use_lmer) "lmerTest_lmer" else "lm", warning = c(warning_text, conditionMessage(fit)), formula_requested = formula_requested))
  }
  list(
    fit = fit,
    model_type = if (use_lmer) "lmerTest_lmer" else "lm",
    warning = warning_text,
    formula_requested = formula_requested,
    formula_used = paste(deparse(stats::formula(fit)), collapse = "")
  )
}

extract_contrasts <- function(fit_info) {
  contrasts <- c("RES - CON", "SUS - CON", "SUS - RES")
  empty <- tibble::tibble(contrast = contrasts, estimate = NA_real_, SE = NA_real_, statistic = NA_real_, p_value = NA_real_)
  if (is.null(fit_info$fit) || !requireNamespace("emmeans", quietly = TRUE)) return(empty)
  emm <- tryCatch(emmeans::emmeans(fit_info$fit, specs = "StressGroup"), error = function(e) e)
  if (inherits(emm, "error")) return(empty)
  contr <- tryCatch(as.data.frame(emmeans::contrast(emm, method = list(
    "RES - CON" = c(-1, 1, 0),
    "SUS - CON" = c(-1, 0, 1),
    "SUS - RES" = c(0, -1, 1)
  ))), error = function(e) e)
  if (inherits(contr, "error")) return(empty)
  stat_col <- intersect(c("t.ratio", "z.ratio"), names(contr))[1]
  tibble::tibble(
    contrast = as.character(contr$contrast),
    estimate = as.numeric(contr$estimate),
    SE = as.numeric(contr$SE),
    statistic = if (!is.na(stat_col)) as.numeric(contr[[stat_col]]) else NA_real_,
    p_value = as.numeric(contr$p.value)
  ) |>
    dplyr::right_join(empty |> dplyr::select("contrast"), by = "contrast")
}

extract_covariate_term <- function(fit_info, term = "matched_neuropil_score") {
  if (is.null(fit_info$fit)) return(list(beta = NA_real_, p = NA_real_))
  cf <- tryCatch(as.data.frame(summary(fit_info$fit)$coefficients), error = function(e) NULL)
  if (is.null(cf) || !term %in% rownames(cf)) return(list(beta = NA_real_, p = NA_real_))
  p_col <- intersect(c("Pr(>|t|)", "Pr(>|z|)", "Pr(>|F|)"), names(cf))[1]
  list(beta = as.numeric(cf[term, "Estimate"]), p = if (!is.na(p_col)) as.numeric(cf[term, p_col]) else NA_real_)
}

make_module_map <- function(eig, definitions) {
  out <- tibble::tibble(endpoint_col = setdiff(names(eig), "Sample"), module_eigengene = endpoint_col, ModuleColor = module_col_to_id(endpoint_col))
  if (nrow(definitions)) {
    defs <- definitions[, intersect(c("ModuleColor", "ModuleID", "module_eigengene", "ModuleLabel_Final", "module_display_label"), names(definitions)), drop = FALSE] |>
      dplyr::distinct()
    if ("module_eigengene" %in% names(defs)) out <- out |> dplyr::left_join(defs, by = "module_eigengene")
  }
  for (nm in c("ModuleID", "ModuleColor.x", "ModuleColor.y", "ModuleLabel_Final", "module_display_label")) if (!nm %in% names(out)) out[[nm]] <- NA_character_
  out |>
    dplyr::mutate(
      endpoint_id = dplyr::coalesce(clean_chr(.data$ModuleID), clean_chr(.data$ModuleColor.x), clean_chr(.data$ModuleColor.y), clean_chr(.data$module_eigengene)),
      endpoint_label = dplyr::coalesce(clean_chr(.data$module_display_label), clean_chr(.data$ModuleLabel_Final), .data$endpoint_id),
      endpoint_type = "module_eigengene"
    )
}

metadata_with_scores <- function(state, dataset, score_df) {
  meta <- standardize_wgcna_metadata(state$sample_info, dataset)
  dplyr::inner_join(meta, score_df, by = "Sample")
}

marker_endpoint_candidates <- function(traits) {
  if (!nrow(traits)) return(tibble::tibble())
  score_cols <- grep("(^z_|^raw_|microglia_minus|microglia_to).*score$|ratio$", names(traits), value = TRUE)
  score_cols <- score_cols[vapply(traits[score_cols], function(x) is.numeric(x) || is.integer(x), logical(1))]
  score_cols <- score_cols[grepl("microglia|immune|dam|lysos|phago|complement|neuropil|neuronal", score_cols, ignore.case = TRUE)]
  tibble::tibble(endpoint_col = score_cols, endpoint_id = score_cols, endpoint_label = score_cols, endpoint_type = "targeted_signature_score")
}

state_micro <- load_wgcna_state(FILES_MICRO$state)
state_neuro <- load_wgcna_state(FILES_NEURO$state)
micro_eig <- extract_module_eigengenes(state_micro)
neuro_eig <- extract_module_eigengenes(state_neuro)
micro_defs <- read_csv_optional2(FILES_MICRO$definitions)
neuro_defs <- read_csv_optional2(FILES_NEURO$definitions)
micro_map <- make_module_map(micro_eig, micro_defs)
neuro_map <- make_module_map(neuro_eig, neuro_defs) |>
  dplyr::mutate(neuropil_covariate_source = "neuron_neuropil_module_eigengene")

config <- yaml::read_yaml(CONFIG_FILE)
predeclared_families <- config$predeclared_adjustments$families %||% list()
primary_family <- as.character(config$predeclared_adjustments$primary_family %||% "global_neuropil_score")
min_matched_animals_required <- as_int(config$predeclared_adjustments$minimum_matched_animals %||% 6L, 6L)
min_matched_samples_required <- as_int(config$predeclared_adjustments$minimum_matched_samples %||% 6L, 6L)
min_animals_per_group_required <- as_int(config$predeclared_adjustments$minimum_animals_per_group %||% 2L, 2L)
primary_effect_fdr_threshold <- suppressWarnings(as.numeric(config$predeclared_adjustments$primary_effect_fdr_threshold %||% 0.05))
if (!is.finite(primary_effect_fdr_threshold)) primary_effect_fdr_threshold <- 0.05
primary_effect_nominal_p_threshold <- suppressWarnings(as.numeric(config$predeclared_adjustments$primary_effect_nominal_p_threshold %||% 0.05))
if (!is.finite(primary_effect_nominal_p_threshold)) primary_effect_nominal_p_threshold <- 0.05
nominal_primary_effect_claim_relevant <- isTRUE(config$predeclared_adjustments$nominal_primary_effect_claim_relevant)
near_zero_effect_abs_threshold <- suppressWarnings(as.numeric(config$predeclared_adjustments$near_zero_effect_abs_threshold %||% 1e-6))
if (!is.finite(near_zero_effect_abs_threshold)) near_zero_effect_abs_threshold <- 1e-6
primary_effect_threshold_label <- paste0(
  "FDR<=", primary_effect_fdr_threshold,
  "; nominal_p<=", primary_effect_nominal_p_threshold,
  if (nominal_primary_effect_claim_relevant) " claim_relevant_by_config" else " diagnostic_only"
)

micro_dat <- metadata_with_scores(state_micro, "microglia", micro_eig)
neuro_dat <- metadata_with_scores(state_neuro, "neuron_neuropil", neuro_eig)
micro_traits <- read_csv_optional2(FILES_MICRO$marker_traits)
neuro_traits <- read_csv_optional2(FILES_NEURO$marker_traits)
if (nrow(micro_traits)) {
  trait_cols <- marker_endpoint_candidates(micro_traits)
  if (nrow(trait_cols)) {
    keep <- unique(c("Sample", trait_cols$endpoint_col))
    micro_dat <- micro_dat |> dplyr::left_join(micro_traits[, intersect(keep, names(micro_traits)), drop = FALSE], by = "Sample")
    micro_map <- dplyr::bind_rows(micro_map, trait_cols)
  }
}

neuro_covariates <- neuro_dat[, c("Sample", "AnimalID", "StressGroup", "Sex", "Batch", "Region", "SpatialLabel", setdiff(names(neuro_eig), "Sample")), drop = FALSE]
if (nrow(neuro_traits)) {
  neuro_score_cols <- grep("(^z_|^raw_).*score$|ratio$", names(neuro_traits), value = TRUE, ignore.case = TRUE)
  neuro_score_cols <- neuro_score_cols[vapply(neuro_traits[neuro_score_cols], function(x) is.numeric(x) || is.integer(x), logical(1))]
  if (length(neuro_score_cols)) {
    neuro_covariates <- neuro_covariates |>
      dplyr::left_join(neuro_traits[, c("Sample", neuro_score_cols), drop = FALSE], by = "Sample")
    neuro_map <- dplyr::bind_rows(
      neuro_map,
      tibble::tibble(
        endpoint_col = neuro_score_cols,
        module_eigengene = NA_character_,
        ModuleColor = NA_character_,
        ModuleID = NA_character_,
        endpoint_id = neuro_score_cols,
        endpoint_label = neuro_score_cols,
        endpoint_type = "neuron_neuropil_signature_score",
        neuropil_covariate_source = "neuron_neuropil_marker_trait"
      )
    )
  }
}

covariate_cols <- intersect(neuro_map$endpoint_col, names(neuro_covariates))
neuro_region <- neuro_covariates |>
  dplyr::mutate(
    AnimalID = as.character(.data$AnimalID),
    Region = as.character(.data$Region),
    StressGroup = as.character(.data$StressGroup),
    Sex = as.character(.data$Sex),
    Batch = as.character(.data$Batch)
  ) |>
  dplyr::filter(!is.na(.data$AnimalID), nzchar(.data$AnimalID), !is.na(.data$Region), nzchar(.data$Region)) |>
  dplyr::group_by(.data$AnimalID, .data$Region) |>
  dplyr::summarise(
    StressGroup_neuropil = dplyr::first(stats::na.omit(.data$StressGroup)) %||% NA_character_,
    Sex_neuropil = dplyr::first(stats::na.omit(.data$Sex)) %||% NA_character_,
    Batch_neuropil = paste(sort(unique(stats::na.omit(.data$Batch))), collapse = ";"),
    dplyr::across(dplyr::all_of(covariate_cols), ~ mean(suppressWarnings(as.numeric(.x)), na.rm = TRUE)),
    n_neuropil_regionlayer_samples = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::mutate(dplyr::across(dplyr::all_of(covariate_cols), ~ ifelse(is.nan(.x), NA_real_, .x)))

count_animal_region <- function(dat, source_label, count_col = "n_samples") {
  dat |>
    dplyr::mutate(AnimalID = as.character(.data$AnimalID), Region = as.character(.data$Region)) |>
    dplyr::filter(!is.na(.data$AnimalID), nzchar(.data$AnimalID), !is.na(.data$Region), nzchar(.data$Region)) |>
    dplyr::group_by(.data$AnimalID, .data$Region) |>
    dplyr::summarise("{count_col}" := dplyr::n(), .groups = "drop") |>
    dplyr::mutate(diagnostic_type = paste0(source_label, "_AnimalID_Region_counts"))
}

micro_counts <- count_animal_region(micro_dat, "microglia", "n_microglia_samples")
neuro_counts <- neuro_region |>
  dplyr::select("AnimalID", "Region", "n_neuropil_regionlayer_samples") |>
  dplyr::mutate(diagnostic_type = "neuropil_AnimalID_Region_counts")
matched_counts <- dplyr::inner_join(
  micro_counts |> dplyr::select("AnimalID", "Region", "n_microglia_samples"),
  neuro_counts |> dplyr::select("AnimalID", "Region", "n_neuropil_regionlayer_samples"),
  by = c("AnimalID", "Region")
) |>
  dplyr::mutate(diagnostic_type = "matched_AnimalID_Region_rows", detail = "microglia and neuron_neuropil both present for AnimalID + Region")
matched_preview <- micro_dat |>
  dplyr::mutate(AnimalID = as.character(.data$AnimalID), Region = as.character(.data$Region)) |>
  dplyr::filter(!is.na(.data$AnimalID), nzchar(.data$AnimalID), !is.na(.data$Region), nzchar(.data$Region)) |>
  dplyr::left_join(neuro_region, by = c("AnimalID", "Region"))
covariate_diagnostics <- dplyr::bind_rows(lapply(covariate_cols, function(cn) {
  vals <- if (cn %in% names(matched_preview)) suppressWarnings(as.numeric(matched_preview[[cn]])) else numeric()
  tibble::tibble(
    diagnostic_type = "matched_neuropil_covariate_variance",
    AnimalID = NA_character_,
    Region = NA_character_,
    neuropil_covariate = cn,
    n_matched_rows = nrow(matched_preview),
    n_finite_neuropil_values = sum(is.finite(vals)),
    variance = if (sum(is.finite(vals)) >= 2L) stats::var(vals, na.rm = TRUE) else NA_real_,
    detail = "finite variance is required before endpoint adjustment models can run"
  )
}))
matching_diagnostics <- dplyr::bind_rows(
  micro_counts |> dplyr::mutate(detail = "microglia rows by AnimalID + Region"),
  neuro_counts |> dplyr::mutate(detail = "neuron_neuropil rows aggregated by AnimalID + Region"),
  matched_counts,
  covariate_diagnostics
) |>
  dplyr::arrange(.data$diagnostic_type, .data$AnimalID, .data$Region, .data$neuropil_covariate)
write_table_and_source(matching_diagnostics, PATHS$tables, PATHS$source_data, "microglia_neuropil_matching_diagnostics.csv")

finite_covariates <- function(dat, candidates, min_n = 4L) {
  candidates <- candidates[candidates %in% names(dat)]
  candidates[vapply(candidates, function(cn) {
    vals <- suppressWarnings(as.numeric(dat[[cn]]))
    sum(is.finite(vals)) >= min_n && is.finite(stats::var(vals, na.rm = TRUE)) && stats::var(vals, na.rm = TRUE) > 0
  }, logical(1))]
}

family_covariates <- function(family, available_covariates) {
  preferred <- as.character(family$preferred_covariates %||% character())
  preferred <- preferred[preferred %in% available_covariates]
  regex <- as.character(family$regex %||% "")
  regex_hits <- if (nzchar(regex)) available_covariates[grepl(regex, available_covariates, ignore.case = TRUE, perl = TRUE)] else character()
  unique(c(preferred, regex_hits))
}

covariate_metadata <- function(cov_col, family, mode, rule, spearman = NA_real_) {
  cov_meta <- neuro_map[match(cov_col, neuro_map$endpoint_col), , drop = FALSE]
  list(
    covariate_col = cov_col,
    covariate_id = cov_meta$endpoint_id[[1]] %||% cov_col,
    covariate_label = cov_meta$endpoint_label[[1]] %||% cov_col,
    covariate_source = cov_meta$neuropil_covariate_source[[1]] %||% "neuron_neuropil_covariate",
    covariate_family = family,
    adjustment_mode = mode,
    covariate_selection_rule = rule,
    selection_spearman = spearman
  )
}

choose_exploratory_covariate <- function(dat) {
  candidates <- finite_covariates(dat, covariate_cols)
  if (!length(candidates)) return(NULL)
  cors <- vapply(candidates, function(cn) {
    ok <- is.finite(dat$score) & is.finite(suppressWarnings(as.numeric(dat[[cn]])))
    if (sum(ok) < 4L) return(NA_real_)
    suppressWarnings(stats::cor(dat$score[ok], as.numeric(dat[[cn]][ok]), method = "spearman"))
  }, numeric(1))
  if (!any(is.finite(cors))) return(NULL)
  cov_col <- names(which.max(abs(cors)))
  covariate_metadata(cov_col, "exploratory_best_spearman", "exploratory_best_spearman", "strongest_absolute_spearman_per_endpoint", unname(cors[[cov_col]]))
}

choose_predeclared_covariates <- function(dat) {
  available <- finite_covariates(dat, covariate_cols)
  out <- lapply(predeclared_families, function(fam) {
    fam_name <- as.character(fam$covariate_family %||% NA_character_)
    hits <- family_covariates(fam, available)
    if (!length(hits)) return(NULL)
    covariate_metadata(
      hits[[1]],
      fam_name,
      if (identical(fam_name, primary_family)) "predeclared_primary" else "predeclared_secondary",
      as.character(fam$selection_rule %||% "fixed_predeclared_family_priority")
    )
  })
  Filter(Negate(is.null), out)
}

selection_audit_row <- function(endpoint_row, dat0, predeclared_choices, exploratory_choice) {
  available <- finite_covariates(dat0, covariate_cols)
  predeclared_available <- unlist(lapply(predeclared_families, function(fam) {
    fam_name <- as.character(fam$covariate_family %||% NA_character_)
    hits <- family_covariates(fam, available)
    if (length(hits)) paste0(fam_name, "=", paste(hits, collapse = "|")) else character()
  }), use.names = FALSE)
  primary <- vapply(predeclared_choices, function(x) if (identical(x$adjustment_mode, "predeclared_primary")) x$covariate_id else NA_character_, character(1))
  secondary <- vapply(predeclared_choices, function(x) if (identical(x$adjustment_mode, "predeclared_secondary")) x$covariate_id else NA_character_, character(1))
  direct_tested <- length(predeclared_choices) > 0L || length(exploratory_choice) > 0L
  tibble::tibble(
    endpoint_id = endpoint_row$endpoint_id[[1]],
    module_or_supermodule_id = ifelse(endpoint_row$endpoint_type[[1]] %in% c("module_eigengene", "supermodule_eigengene", "supermodule_score"), endpoint_row$endpoint_id[[1]], NA_character_),
    endpoint_type = endpoint_row$endpoint_type[[1]],
    endpoint_scope = endpoint_scope_for(endpoint_row$endpoint_type[[1]]),
    source_level = endpoint_scope_for(endpoint_row$endpoint_type[[1]]),
    direct_independence_tested = direct_tested,
    direct_supermodule_test = direct_tested && endpoint_scope_for(endpoint_row$endpoint_type[[1]]) == "supermodule",
    candidate_covariates = paste(available, collapse = ";"),
    predeclared_covariates_available = paste(predeclared_available, collapse = ";"),
    selected_primary_covariate = paste(stats::na.omit(primary), collapse = ";"),
    selected_secondary_covariate = paste(stats::na.omit(secondary), collapse = ";"),
    selected_exploratory_covariate = exploratory_choice$covariate_id %||% NA_character_,
    selection_rule = "predeclared=fixed family priority from config; exploratory=strongest absolute Spearman per endpoint",
    primary_claim_gate_eligible = any(vapply(predeclared_choices, function(x) identical(x$adjustment_mode, "predeclared_primary"), logical(1))),
    secondary_claim_gate_eligible = any(vapply(predeclared_choices, function(x) identical(x$adjustment_mode, "predeclared_secondary"), logical(1))),
    exploratory_claim_gate_eligible = FALSE
  )
}

primary_effect_status_for <- function(p, fdr, endpoint_type = "module_eigengene") {
  p <- suppressWarnings(as.numeric(p))
  fdr <- suppressWarnings(as.numeric(fdr))
  dplyr::case_when(
    !endpoint_type %in% c("module_eigengene", "targeted_signature_score") ~ "not_applicable",
    is.na(p) & is.na(fdr) ~ "missing_effect_statistic",
    is.finite(fdr) & fdr <= primary_effect_fdr_threshold ~ "FDR_pass",
    is.finite(p) & p <= primary_effect_nominal_p_threshold ~ "nominal_only",
    TRUE ~ "no_primary_effect"
  )
}

primary_effect_claim_relevant_for <- function(status) {
  status == "FDR_pass" | (status == "nominal_only" & isTRUE(nominal_primary_effect_claim_relevant))
}

classify_adjustment <- function(effect_before, effect_after, fdr_after, n_animals, n_samples, min_group, mode,
                                primary_effect_claim_relevant, primary_effect_status,
                                percent_attenuation_reliable) {
  direction_preserved <- is.finite(effect_before) & is.finite(effect_after) & sign(effect_before) == sign(effect_after) & sign(effect_before) != 0
  attenuation <- ifelse(is.finite(effect_before) & abs(effect_before) > 1e-12, 100 * (abs(effect_before) - abs(effect_after)) / abs(effect_before), NA_real_)
  low_power <- is.na(n_animals) | is.na(n_samples) | is.na(min_group) |
    n_animals < min_matched_animals_required | n_samples < min_matched_samples_required | min_group < min_animals_per_group_required
  dplyr::case_when(
    identical(mode, "exploratory_best_spearman") ~ "exploratory_only",
    isTRUE(low_power) ~ "inconclusive_low_power",
    is.na(effect_before) | is.na(effect_after) ~ "inconclusive_missing_match",
    !isTRUE(primary_effect_claim_relevant) & primary_effect_status %in% c("missing_effect_statistic", "not_applicable") ~ "inconclusive_no_primary_effect",
    !isTRUE(primary_effect_claim_relevant) ~ "diagnostic_no_primary_effect",
    !direction_preserved ~ "neuropil_sensitive",
    isTRUE(percent_attenuation_reliable) & is.finite(attenuation) & attenuation >= 60 ~ "neuropil_sensitive",
    isTRUE(percent_attenuation_reliable) & is.finite(attenuation) & attenuation >= 25 ~ "partially_neuropil_adjusted",
    direction_preserved ~ "neuropil_independent",
    TRUE ~ "inconclusive_missing_match"
  )
}

fit_endpoint_adjustment <- function(endpoint_row, dat0, cov_choice) {
  dat <- dat0 |>
    dplyr::mutate(
      matched_neuropil_score = suppressWarnings(as.numeric(.data[[cov_choice$covariate_col]])),
      StressGroup = factor(.data$StressGroup, levels = c("CON", "RES", "SUS")),
      SpatialLabel = factor(.data$SpatialLabel),
      Sex = factor(.data$Sex),
      Batch = factor(.data$Batch)
    ) |>
    dplyr::filter(is.finite(.data$matched_neuropil_score))
  cov_base <- valid_covariates(dat, c("StressGroup", "SpatialLabel", "Sex", "Batch"), protect = c("StressGroup", "SpatialLabel"))
  cov_adj <- valid_covariates(dat, c("StressGroup", "matched_neuropil_score", "SpatialLabel", "Sex", "Batch"), protect = c("StressGroup", "matched_neuropil_score", "SpatialLabel"))
  before <- fit_score_model(dat, cov_base$keep)
  after <- fit_score_model(dat, cov_adj$keep)
  beta <- extract_covariate_term(after)
  b <- extract_contrasts(before) |> dplyr::rename(effect_before = "estimate", SE_before = "SE", statistic_before = "statistic", p_before = "p_value")
  a <- extract_contrasts(after) |> dplyr::rename(effect_after = "estimate", SE_after = "SE", statistic_after = "statistic", p_after = "p_value")
  group_counts <- dat |>
    dplyr::filter(!is.na(.data$AnimalID), nzchar(as.character(.data$AnimalID))) |>
    dplyr::distinct(.data$AnimalID, .data$StressGroup) |>
    dplyr::count(.data$StressGroup, name = "n_animals")
  min_animals_per_group <- if (nrow(group_counts)) min(group_counts$n_animals, na.rm = TRUE) else NA_integer_
  dplyr::left_join(b, a, by = "contrast") |>
    dplyr::mutate(
      dataset = "microglia",
      endpoint_type = endpoint_row$endpoint_type,
      endpoint_scope = endpoint_scope_for(endpoint_row$endpoint_type),
      source_level = endpoint_scope_for(endpoint_row$endpoint_type),
      direct_independence_tested = TRUE,
      direct_supermodule_test = endpoint_scope_for(endpoint_row$endpoint_type) == "supermodule",
      endpoint_id = endpoint_row$endpoint_id,
      endpoint_label = endpoint_row$endpoint_label,
      module_id = ifelse(endpoint_row$endpoint_type == "module_eigengene", endpoint_row$endpoint_id, NA_character_),
      adjustment_mode = cov_choice$adjustment_mode,
      covariate_family = cov_choice$covariate_family,
      covariate_selection_rule = cov_choice$covariate_selection_rule,
      covariate_source_dataset = "neuron_neuropil",
      covariate_source_module_or_score = cov_choice$covariate_id,
      matched_by = "AnimalID + Region",
      percent_attenuation = dplyr::if_else(is.finite(.data$effect_before) & abs(.data$effect_before) > 1e-12,
                                           100 * (abs(.data$effect_before) - abs(.data$effect_after)) / abs(.data$effect_before),
                                           NA_real_),
      direction_preserved = is.finite(.data$effect_before) & is.finite(.data$effect_after) & sign(.data$effect_before) == sign(.data$effect_after) & sign(.data$effect_before) != 0,
      neuropil_covariate_beta = beta$beta,
      neuropil_covariate_p = beta$p,
      matched_neuropil_covariate = cov_choice$covariate_id,
      matched_neuropil_label = cov_choice$covariate_label,
      neuropil_covariate_source = cov_choice$covariate_source,
      neuropil_selection_spearman = cov_choice$selection_spearman,
      n_samples = nrow(dat0),
      n_matched_samples = nrow(dat),
      n_matched_animals = dplyr::n_distinct(dat$AnimalID),
      n_animals = dplyr::n_distinct(dat$AnimalID),
      min_animals_per_group = min_animals_per_group,
      model_type_before = before$model_type,
      model_type_after = after$model_type,
      model_before_adjustment = before$formula_used %||% before$formula_requested,
      model_after_adjustment = after$formula_used %||% after$formula_requested,
      formula_before = .data$model_before_adjustment,
      formula_after = .data$model_after_adjustment,
      model_warning = paste(unique(c(before$warning, after$warning, cov_base$dropped, cov_adj$dropped)), collapse = ";")
    )
}

empty_adjustment_rows <- function(endpoint_row, reason, class = "inconclusive_missing_match") {
  tibble::tibble(
    dataset = "microglia",
    endpoint_type = endpoint_row$endpoint_type,
    endpoint_scope = endpoint_scope_for(endpoint_row$endpoint_type),
    source_level = endpoint_scope_for(endpoint_row$endpoint_type),
    direct_independence_tested = FALSE,
    direct_supermodule_test = FALSE,
    endpoint_id = endpoint_row$endpoint_id,
    endpoint_label = endpoint_row$endpoint_label,
    module_id = ifelse(endpoint_row$endpoint_type == "module_eigengene", endpoint_row$endpoint_id, NA_character_),
    contrast = c("RES - CON", "SUS - CON", "SUS - RES"),
    adjustment_mode = "predeclared_primary",
    covariate_family = primary_family,
    covariate_selection_rule = "fixed_predeclared_family_priority",
    covariate_source_dataset = "neuron_neuropil",
    covariate_source_module_or_score = NA_character_,
    matched_by = "AnimalID + Region",
    n_matched_animals = 0L,
    n_matched_samples = 0L,
    min_animals_per_group = NA_integer_,
    effect_before = NA_real_,
    effect_after = NA_real_,
    estimate_before = NA_real_,
    estimate_after = NA_real_,
    percent_attenuation = NA_real_,
    effect_before_abs = NA_real_,
    effect_before_near_zero = NA,
    percent_attenuation_reliable = NA,
    direction_preserved = NA,
    p_before = NA_real_,
    p_after = NA_real_,
    FDR_before = NA_real_,
    FDR_after = NA_real_,
    primary_effect_status = "missing_effect_statistic",
    primary_effect_claim_relevant = FALSE,
    primary_effect_threshold = primary_effect_threshold_label,
    primary_effect_p = NA_real_,
    primary_effect_FDR = NA_real_,
    neuropil_covariate_beta = NA_real_,
    neuropil_covariate_p = NA_real_,
    matched_neuropil_covariate = NA_character_,
    matched_neuropil_label = NA_character_,
    neuropil_covariate_source = NA_character_,
    neuropil_selection_spearman = NA_real_,
    n_samples = NA_integer_,
    n_animals = NA_integer_,
    model_type_before = NA_character_,
    model_type_after = NA_character_,
    formula_before = NA_character_,
    formula_after = NA_character_,
    model_before_adjustment = NA_character_,
    model_after_adjustment = NA_character_,
    independence_classification = class,
    claim_gate_eligible = FALSE,
    downgrade_reason = reason,
    classification = class,
    model_warning = reason
  )
}

fit_endpoint <- function(endpoint_row) {
  endpoint_col <- endpoint_row$endpoint_col[[1]]
  dat0 <- micro_dat |>
    dplyr::mutate(
      AnimalID = as.character(.data$AnimalID),
      Region = as.character(.data$Region),
      score = suppressWarnings(as.numeric(.data[[endpoint_col]]))
    ) |>
    dplyr::filter(is.finite(.data$score), !is.na(.data$StressGroup), .data$StressGroup %in% c("CON", "RES", "SUS")) |>
    dplyr::left_join(neuro_region, by = c("AnimalID", "Region"))
  if (nrow(dat0) < 6L || dplyr::n_distinct(dat0$StressGroup) < 2L) {
    return(empty_adjustment_rows(endpoint_row, "too few microglia samples/groups", "inconclusive_low_power"))
  }
  predeclared_choices <- choose_predeclared_covariates(dat0)
  exploratory_choice <- choose_exploratory_covariate(dat0)
  selection_audit_rows[[length(selection_audit_rows) + 1L]] <<- selection_audit_row(endpoint_row, dat0, predeclared_choices, exploratory_choice %||% list())
  choices <- c(predeclared_choices, if (!is.null(exploratory_choice)) list(exploratory_choice) else list())
  if (!length(choices)) return(empty_adjustment_rows(endpoint_row, "no matched neuron_neuropil covariate with finite variance", "inconclusive_missing_match"))
  dplyr::bind_rows(lapply(choices, function(choice) fit_endpoint_adjustment(endpoint_row, dat0, choice)))
}

endpoint_rows <- micro_map |>
  dplyr::filter(.data$endpoint_col %in% names(micro_dat)) |>
  dplyr::distinct(.data$endpoint_col, .keep_all = TRUE)
selection_audit_rows <- list()
if (!nrow(endpoint_rows)) {
  results <- empty_adjustment_rows(
    tibble::tibble(endpoint_type = NA_character_, endpoint_id = NA_character_, endpoint_label = NA_character_),
    "no microglia module or targeted-signature endpoints available",
    "inconclusive_missing_match"
  )
  selection_audit <- tibble::tibble(
    endpoint_id = character(), module_or_supermodule_id = character(), endpoint_type = character(),
    endpoint_scope = character(), source_level = character(), direct_independence_tested = logical(),
    direct_supermodule_test = logical(), candidate_covariates = character(),
    predeclared_covariates_available = character(), selected_primary_covariate = character(),
    selected_secondary_covariate = character(), selected_exploratory_covariate = character(),
    selection_rule = character(), primary_claim_gate_eligible = logical(),
    secondary_claim_gate_eligible = logical(), exploratory_claim_gate_eligible = logical()
  )
} else {
  results <- dplyr::bind_rows(lapply(seq_len(nrow(endpoint_rows)), function(i) fit_endpoint(endpoint_rows[i, , drop = FALSE])))
  selection_audit <- dplyr::bind_rows(selection_audit_rows)
}

if (!"estimate_before" %in% names(results)) results$estimate_before <- NA_real_
if (!"estimate_after" %in% names(results)) results$estimate_after <- NA_real_
if (!"effect_before" %in% names(results)) results$effect_before <- NA_real_
if (!"effect_after" %in% names(results)) results$effect_after <- NA_real_
results$estimate_before <- dplyr::coalesce(suppressWarnings(as.numeric(results$estimate_before)), suppressWarnings(as.numeric(results$effect_before)))
results$estimate_after <- dplyr::coalesce(suppressWarnings(as.numeric(results$estimate_after)), suppressWarnings(as.numeric(results$effect_after)))
results$effect_before <- dplyr::coalesce(suppressWarnings(as.numeric(results$effect_before)), suppressWarnings(as.numeric(results$estimate_before)))
results$effect_after <- dplyr::coalesce(suppressWarnings(as.numeric(results$effect_after)), suppressWarnings(as.numeric(results$estimate_after)))
results$FDR_before <- NA_real_
results$FDR_after <- NA_real_
ok_after <- is.finite(results$p_after)
results$FDR_after[ok_after] <- stats::p.adjust(results$p_after[ok_after], method = "BH")
primary_fdr_lookup <- results |>
  dplyr::distinct(.data$endpoint_type, .data$endpoint_id, .data$contrast, .data$p_before) |>
  dplyr::mutate(primary_FDR_before = NA_real_)
ok_primary <- is.finite(primary_fdr_lookup$p_before)
primary_fdr_lookup$primary_FDR_before[ok_primary] <- stats::p.adjust(primary_fdr_lookup$p_before[ok_primary], method = "BH")
results <- results |>
  dplyr::select(-"FDR_before") |>
  dplyr::left_join(
    primary_fdr_lookup |>
      dplyr::select("endpoint_type", "endpoint_id", "contrast", "p_before", "FDR_before" = "primary_FDR_before"),
    by = c("endpoint_type", "endpoint_id", "contrast", "p_before")
  )
results$primary_effect_p <- suppressWarnings(as.numeric(results$p_before))
results$primary_effect_FDR <- suppressWarnings(as.numeric(results$FDR_before))
results$primary_effect_status <- mapply(
  primary_effect_status_for,
  results$primary_effect_p,
  results$primary_effect_FDR,
  results$endpoint_type,
  SIMPLIFY = TRUE,
  USE.NAMES = FALSE
)
results$primary_effect_claim_relevant <- primary_effect_claim_relevant_for(results$primary_effect_status)
results$primary_effect_threshold <- primary_effect_threshold_label
results$effect_before_abs <- abs(suppressWarnings(as.numeric(results$effect_before)))
results$effect_before_near_zero <- is.finite(results$effect_before_abs) & results$effect_before_abs < near_zero_effect_abs_threshold
results$percent_attenuation_reliable <- is.finite(results$percent_attenuation) & !results$effect_before_near_zero

results <- results |>
  dplyr::mutate(
    independence_classification = mapply(
      classify_adjustment,
      .data$effect_before, .data$effect_after, .data$FDR_after,
      .data$n_matched_animals, .data$n_matched_samples, .data$min_animals_per_group,
      .data$adjustment_mode, .data$primary_effect_claim_relevant,
      .data$primary_effect_status, .data$percent_attenuation_reliable,
      SIMPLIFY = TRUE, USE.NAMES = FALSE
    ),
    claim_gate_eligible = .data$adjustment_mode %in% c("predeclared_primary", "predeclared_secondary") &
      .data$primary_effect_claim_relevant &
      .data$independence_classification %in% c("neuropil_independent", "partially_neuropil_adjusted") &
      .data$n_matched_animals >= min_matched_animals_required &
      .data$n_matched_samples >= min_matched_samples_required &
      .data$min_animals_per_group >= min_animals_per_group_required &
      (.data$percent_attenuation_reliable | .data$independence_classification == "neuropil_independent"),
    downgrade_reason = dplyr::case_when(
      .data$claim_gate_eligible ~ "none",
      .data$adjustment_mode == "exploratory_best_spearman" ~ "exploratory_best_spearman_diagnostic_only",
      !.data$primary_effect_claim_relevant ~ paste0(.data$independence_classification, "; primary_effect_status=", .data$primary_effect_status),
      .data$effect_before_near_zero ~ paste0(.data$independence_classification, "; unstable_attenuation_near_zero_effect"),
      !.data$percent_attenuation_reliable & .data$independence_classification != "neuropil_independent" ~ paste0(.data$independence_classification, "; attenuation_unreliable"),
      TRUE ~ .data$independence_classification
    ),
    classification = .data$independence_classification
  ) |>
  dplyr::relocate(
    "dataset", "endpoint_type", "endpoint_scope", "source_level", "direct_independence_tested",
    "direct_supermodule_test", "endpoint_id", "endpoint_label", "module_id", "contrast",
    "adjustment_mode", "covariate_family", "covariate_selection_rule", "covariate_source_dataset",
    "covariate_source_module_or_score", "matched_by", "n_matched_animals", "n_matched_samples",
    "min_animals_per_group", "model_before_adjustment", "model_after_adjustment",
    "primary_effect_status", "primary_effect_claim_relevant", "primary_effect_threshold",
    "primary_effect_p", "primary_effect_FDR",
    "effect_before", "effect_after", "estimate_before", "estimate_after",
    "effect_before_abs", "effect_before_near_zero", "percent_attenuation",
    "percent_attenuation_reliable", "direction_preserved", "p_before", "p_after",
    "FDR_before", "FDR_after", "neuropil_covariate_beta", "neuropil_covariate_p",
    "independence_classification", "claim_gate_eligible", "downgrade_reason", "classification"
  )

module_independence_classification <- function(primary_relevant, primary_status, classification) {
  primary_relevant <- suppressWarnings(as.logical(primary_relevant))
  primary_relevant[is.na(primary_relevant)] <- FALSE
  primary_status <- as.character(primary_status)
  classification <- as.character(classification)
  if (!any(primary_relevant)) {
    if (any(primary_status == "missing_effect_statistic", na.rm = TRUE)) return("inconclusive_no_primary_effect")
    return("diagnostic_no_primary_effect")
  }
  relevant_classes <- unique(stats::na.omit(classification[primary_relevant]))
  if ("neuropil_sensitive" %in% relevant_classes) return("neuropil_sensitive")
  if (any(relevant_classes %in% c("inconclusive_low_power", "inconclusive_missing_match", "inconclusive_no_primary_effect"))) return("inconclusive_missing_match")
  if (length(relevant_classes) > 1L) return("mixed_or_covariate_sensitive")
  if (length(relevant_classes) == 1L) return(relevant_classes[[1]])
  "inconclusive_no_primary_effect"
}

module_claim_gate_eligible <- function(primary_relevant, classification, eligible) {
  primary_relevant <- suppressWarnings(as.logical(primary_relevant))
  primary_relevant[is.na(primary_relevant)] <- FALSE
  if (!any(primary_relevant)) return(FALSE)
  relevant_classes <- unique(stats::na.omit(as.character(classification[primary_relevant])))
  length(relevant_classes) == 1L &&
    relevant_classes %in% c("neuropil_independent", "partially_neuropil_adjusted") &&
    all(suppressWarnings(as.logical(eligible[primary_relevant])), na.rm = TRUE)
}

module_classification <- results |>
  dplyr::filter(.data$endpoint_type == "module_eigengene", .data$adjustment_mode %in% c("predeclared_primary", "predeclared_secondary")) |>
  dplyr::group_by(.data$dataset, .data$module_id, .data$endpoint_id, .data$endpoint_label) |>
  dplyr::summarise(
    n_primary_effect_claim_relevant = sum(.data$primary_effect_claim_relevant, na.rm = TRUE),
    primary_effect_status_summary = paste(unique(stats::na.omit(.data$primary_effect_status)), collapse = ";"),
    neuropil_independence_classification = module_independence_classification(.data$primary_effect_claim_relevant, .data$primary_effect_status, .data$independence_classification),
    claim_gate_eligible = module_claim_gate_eligible(.data$primary_effect_claim_relevant, .data$independence_classification, .data$claim_gate_eligible),
    best_contrast = .data$contrast[which.min(dplyr::coalesce(.data$p_before, Inf))][[1]] %||% NA_character_,
    min_p_before = suppressWarnings(min(.data$p_before, na.rm = TRUE)),
    min_p_after = suppressWarnings(min(.data$p_after, na.rm = TRUE)),
    max_percent_attenuation = suppressWarnings(max(.data$percent_attenuation, na.rm = TRUE)),
    any_unstable_attenuation_near_zero_effect = any(.data$effect_before_near_zero, na.rm = TRUE),
    matched_neuropil_covariate = paste(unique(stats::na.omit(.data$matched_neuropil_covariate)), collapse = ";"),
    neuropil_covariate_source = paste(unique(stats::na.omit(.data$neuropil_covariate_source)), collapse = ";"),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    dplyr::across(c("min_p_before", "min_p_after", "max_percent_attenuation"), ~ ifelse(is.infinite(.x), NA_real_, .x)),
    endpoint_type = "module_eigengene",
    endpoint_scope = "module",
    source_level = "module",
    direct_independence_tested = TRUE,
    direct_supermodule_test = FALSE,
    neuropil_independence_note = "Classification from paired baseline versus AnimalID + Region matched neuron_neuropil-adjusted models. Only predeclared primary/secondary adjustments are claim-gate eligible."
  )

signature_results <- results |> dplyr::filter(.data$endpoint_type == "targeted_signature_score")
module_annotation <- read_csv_optional2(path_results("tables", "06_modules_WGCNA", "module_annotation", "microglia", "WGCNA_module_biological_annotation.csv"))
claim_gate_audit <- results |>
  dplyr::filter(.data$endpoint_type == "module_eigengene") |>
  dplyr::left_join(
    module_annotation |>
      dplyr::select(dplyr::any_of(c("ModuleID", "cleaned_biological_label", "safe_display_label", "module_biological_label", "microenvironment_class"))) |>
      dplyr::distinct(),
    by = c("module_id" = "ModuleID")
  ) |>
  dplyr::transmute(
    module_or_supermodule_id = .data$module_id,
    endpoint_type = .data$endpoint_type,
    endpoint_scope = .data$endpoint_scope,
    source_level = .data$source_level,
    direct_independence_tested = .data$direct_independence_tested,
    direct_supermodule_test = .data$direct_supermodule_test,
    contrast = .data$contrast,
    biological_program = dplyr::coalesce(.data$cleaned_biological_label, .data$safe_display_label, .data$module_biological_label, .data$endpoint_label),
    microenvironment_class = .data$microenvironment_class,
    adjustment_mode = .data$adjustment_mode,
    covariate_family = .data$covariate_family,
    primary_effect_status = .data$primary_effect_status,
    primary_effect_claim_relevant = .data$primary_effect_claim_relevant,
    primary_effect_threshold = .data$primary_effect_threshold,
    independence_classification = .data$independence_classification,
    claim_gate_eligible = .data$claim_gate_eligible,
    downgrade_reason = .data$downgrade_reason,
    n_matched_animals = .data$n_matched_animals,
    min_animals_per_group = .data$min_animals_per_group,
    effect_before = .data$effect_before,
    effect_after = .data$effect_after,
    effect_before_abs = .data$effect_before_abs,
    effect_before_near_zero = .data$effect_before_near_zero,
    percent_attenuation = .data$percent_attenuation,
    percent_attenuation_reliable = .data$percent_attenuation_reliable,
    direction_preserved = .data$direction_preserved
  )
claim_gate_counts <- claim_gate_audit |>
  dplyr::count(
    .data$adjustment_mode,
    .data$primary_effect_status,
    .data$independence_classification,
    .data$claim_gate_eligible,
    name = "audit_group_n"
  )
eligible_without_primary <- sum(
  claim_gate_audit$claim_gate_eligible & !claim_gate_audit$primary_effect_claim_relevant,
  na.rm = TRUE
)
claim_gate_audit <- claim_gate_audit |>
  dplyr::left_join(
    claim_gate_counts,
    by = c("adjustment_mode", "primary_effect_status", "independence_classification", "claim_gate_eligible")
  ) |>
  dplyr::mutate(eligible_without_primary_effect_count = eligible_without_primary)

endpoint_scope_audit <- results |>
  dplyr::filter(!is.na(.data$endpoint_type), .data$endpoint_scope != "unavailable") |>
  dplyr::distinct(.data$endpoint_type, .data$endpoint_scope, .data$source_level, .data$endpoint_id, .data$direct_independence_tested) |>
  dplyr::group_by(.data$endpoint_type, .data$endpoint_scope, .data$source_level, .data$direct_independence_tested) |>
  dplyr::summarise(n_endpoints = dplyr::n(), .groups = "drop") |>
  dplyr::mutate(
    consumed_by_claim_type = dplyr::case_when(
      .data$endpoint_scope == "module" ~ "WGCNA_module;WGCNA_module_group_effect",
      .data$endpoint_scope == "targeted_signature" ~ "targeted_signature_diagnostic",
      .data$endpoint_scope == "supermodule" ~ "WGCNA_supermodule_group_effect",
      TRUE ~ "none"
    ),
    notes = "Direct endpoint-level neuropil-independence model was run."
  ) |>
  dplyr::bind_rows(tibble::tibble(
    endpoint_type = "unavailable", endpoint_scope = "unavailable", source_level = "supermodule",
    n_endpoints = 0L, direct_independence_tested = FALSE,
    consumed_by_claim_type = "WGCNA_supermodule_group_effect",
    notes = "no_direct_supermodule_independence_test"
  )) |>
  dplyr::select("endpoint_type", "endpoint_scope", "source_level", "n_endpoints", "direct_independence_tested", "consumed_by_claim_type", "notes")

write_table_and_source(results, PATHS$tables, PATHS$source_data, "microglia_neuropil_independence_effects.csv")
write_table_and_source(module_classification, PATHS$tables, PATHS$source_data, "microglia_module_neuropil_independence_classification.csv")
write_table_and_source(signature_results, PATHS$tables, PATHS$source_data, "microglia_targeted_signature_neuropil_independence_effects.csv")

audit_dir <- path_results("reviewer_audit")
dir_create(audit_dir)
validate_table_schema(claim_gate_audit, "microglia_neuropil_independence_claim_gate", strict = TRUE)
validate_table_schema(selection_audit, "microglia_neuropil_covariate_selection_audit", strict = FALSE)
validate_table_schema(endpoint_scope_audit, "microglia_neuropil_independence_endpoint_scope_audit", strict = TRUE)
readr::write_csv(claim_gate_audit, file.path(audit_dir, "microglia_neuropil_independence_claim_gate.csv"), na = "")
readr::write_csv(selection_audit, file.path(audit_dir, "microglia_neuropil_covariate_selection_audit.csv"), na = "")
readr::write_csv(endpoint_scope_audit, file.path(audit_dir, "microglia_neuropil_independence_endpoint_scope_audit.csv"), na = "")

interp <- read_csv_optional2(inputs$microglia_interpretable)
if (nrow(interp) && nrow(module_classification)) {
  join_key <- intersect(c("module_id", "ModuleID"), names(interp))[1]
  class_key <- if (!is.na(join_key) && join_key == "ModuleID") "endpoint_id" else "module_id"
  joined <- interp |>
    dplyr::left_join(module_classification, by = stats::setNames(class_key, join_key))
  write_table_and_source(joined, PATHS$tables, PATHS$source_data, "WGCNA_module_group_effects_interpretable_with_neuropil_independence.csv")
}

if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(
    list(
      effects = results,
      module_classification = module_classification,
      targeted_signatures = signature_results,
      claim_gate_audit = claim_gate_audit,
      covariate_selection_audit = selection_audit,
      endpoint_scope_audit = endpoint_scope_audit
    ),
    file.path(PATHS$tables, "microglia_neuropil_independence.xlsx")
  )
}

write_run_manifest(
  file.path(PATHS$logs, "run_manifest.yml"),
  inputs = inputs,
  outputs = list(
    tables = PATHS$tables,
    source_data = PATHS$source_data,
    reviewer_audit = c(
      claim_gate = file.path(audit_dir, "microglia_neuropil_independence_claim_gate.csv"),
      covariate_selection = file.path(audit_dir, "microglia_neuropil_covariate_selection_audit.csv"),
      endpoint_scope = file.path(audit_dir, "microglia_neuropil_independence_endpoint_scope_audit.csv")
    )
  ),
  parameters = list(
    dataset = "microglia",
    baseline_model = "score ~ StressGroup + SpatialLabel + Sex + Batch + (1 | AnimalID)",
    adjusted_model = "score ~ StressGroup + matched_neuropil_score + SpatialLabel + Sex + Batch + (1 | AnimalID)",
    neuropil_matching = "neuron_neuropil layer-level values aggregated to AnimalID + Region; predeclared covariates selected by fixed config families; strongest absolute Spearman is retained as exploratory diagnostic only",
    independence_config = CONFIG_FILE
  ),
  notes = "Additive downstream sensitivity analysis. Primary group-effect outputs are not overwritten."
)

message("Microglia neuropil-independence sensitivity complete.")

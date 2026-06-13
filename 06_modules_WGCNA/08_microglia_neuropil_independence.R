#!/usr/bin/env Rscript
# ================================================================
# Script: 06_modules_WGCNA/08_microglia_neuropil_independence.R
# Stage: modules_downstream
# Scope: dataset_specific
# Consumes: microglia and neuron_neuropil WGCNA state/metadata plus optional marker traits.
# Produces: additive microglia neuropil-independence sensitivity tables.
# Notes: Does not overwrite primary WGCNA group-effect outputs.
# ================================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "wgcna_downstream_utils.R"))

SCRIPT_ID <- "06_modules_WGCNA/08_microglia_neuropil_independence.R"
required_pkgs <- c("dplyr", "tidyr", "tibble", "readr")
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
inputs <- list(
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
  for (nm in names(inputs)) dry_run_line(nm, inputs[[nm]], if (file.exists(inputs[[nm]])) "PASS" else "WARN")
  quit(status = 0, save = "no")
}

read_csv_optional2 <- function(path) safe_read_csv(path) %||% data.frame()

clean_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(trimws(x))] <- NA_character_
  x
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
  fit_formula <- stats::as.formula(formula_requested)
  fit <- tryCatch(
    if (use_lmer) lmerTest::lmer(fit_formula, data = dat, REML = FALSE) else stats::lm(fit_formula, data = dat),
    warning = function(w) {
      warning_text <<- c(warning_text, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(fit = NULL, model_type = if (use_lmer) "lmerTest_lmer" else "lm", warning = c(warning_text, conditionMessage(fit)), formula_requested = formula_requested))
  }
  rank_deficient <- tryCatch({
    X <- if (use_lmer && requireNamespace("lme4", quietly = TRUE)) lme4::getME(fit, "X") else stats::model.matrix(fit)
    qr(X)$rank < ncol(X)
  }, error = function(e) NA)
  list(
    fit = fit,
    model_type = if (use_lmer) "lmerTest_lmer" else "lm",
    warning = warning_text,
    formula_requested = formula_requested,
    formula_used = paste(deparse(stats::formula(fit)), collapse = ""),
    rank_deficient = rank_deficient
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
  list(
    beta = as.numeric(cf[term, "Estimate"]),
    p = if (!is.na(p_col)) as.numeric(cf[term, p_col]) else NA_real_
  )
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

state_micro <- load_wgcna_state(FILES_MICRO$state)
state_neuro <- load_wgcna_state(FILES_NEURO$state)
micro_eig <- extract_module_eigengenes(state_micro)
neuro_eig <- extract_module_eigengenes(state_neuro)
micro_defs <- read_csv_optional2(FILES_MICRO$definitions)
neuro_defs <- read_csv_optional2(FILES_NEURO$definitions)
micro_map <- make_module_map(micro_eig, micro_defs)
neuro_map <- make_module_map(neuro_eig, neuro_defs) |>
  dplyr::mutate(neuropil_covariate_source = "neuron_neuropil_module_eigengene")

micro_dat <- metadata_with_scores(state_micro, "microglia", micro_eig)
neuro_dat <- metadata_with_scores(state_neuro, "neuron_neuropil", neuro_eig)

marker_endpoint_candidates <- function(traits) {
  if (!nrow(traits)) return(tibble::tibble())
  score_cols <- grep("(^z_|^raw_|microglia_minus|microglia_to).*score$|ratio$", names(traits), value = TRUE)
  score_cols <- score_cols[vapply(traits[score_cols], function(x) is.numeric(x) || is.integer(x), logical(1))]
  score_cols <- score_cols[grepl("microglia|immune|dam|lysos|phago|complement|neuropil|neuronal", score_cols, ignore.case = TRUE)]
  tibble::tibble(
    endpoint_col = score_cols,
    endpoint_id = score_cols,
    endpoint_label = score_cols,
    endpoint_type = "targeted_signature_score"
  )
}

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
  neuro_score_cols <- grep("(^z_|^raw_).*neuropil.*score$|(^z_|^raw_).*neuronal.*score$", names(neuro_traits), value = TRUE, ignore.case = TRUE)
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

choose_covariate <- function(dat, endpoint_col) {
  candidates <- covariate_cols[covariate_cols %in% names(dat)]
  candidates <- candidates[vapply(dat[candidates], function(x) sum(is.finite(suppressWarnings(as.numeric(x)))) >= 4L && stats::var(suppressWarnings(as.numeric(x)), na.rm = TRUE) > 0, logical(1))]
  if (!length(candidates)) return(NULL)
  cors <- vapply(candidates, function(cn) {
    ok <- is.finite(dat$score) & is.finite(suppressWarnings(as.numeric(dat[[cn]])))
    if (sum(ok) < 4L) return(NA_real_)
    suppressWarnings(stats::cor(dat$score[ok], as.numeric(dat[[cn]][ok]), method = "spearman"))
  }, numeric(1))
  if (!any(is.finite(cors))) return(NULL)
  cov_col <- names(which.max(abs(cors)))
  cov_meta <- neuro_map[match(cov_col, neuro_map$endpoint_col), , drop = FALSE]
  list(
    covariate_col = cov_col,
    covariate_id = cov_meta$endpoint_id[[1]] %||% cov_col,
    covariate_label = cov_meta$endpoint_label[[1]] %||% cov_col,
    covariate_source = cov_meta$neuropil_covariate_source[[1]] %||% "neuron_neuropil_covariate",
    selection_spearman = unname(cors[[cov_col]])
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
    return(tibble::tibble(
      dataset = "microglia", endpoint_type = endpoint_row$endpoint_type, endpoint_id = endpoint_row$endpoint_id,
      endpoint_label = endpoint_row$endpoint_label, module_id = ifelse(endpoint_row$endpoint_type == "module_eigengene", endpoint_row$endpoint_id, NA_character_),
      contrast = c("RES - CON", "SUS - CON", "SUS - RES"), estimate_before = NA_real_, estimate_after = NA_real_,
      percent_attenuation = NA_real_, p_before = NA_real_, p_after = NA_real_, neuropil_covariate_beta = NA_real_,
      neuropil_covariate_p = NA_real_, classification = "unstable_or_unmatched", model_warning = "too few microglia samples/groups"
    ))
  }
  cov_choice <- choose_covariate(dat0, endpoint_col)
  if (is.null(cov_choice)) {
    rows <- tibble::tibble(contrast = c("RES - CON", "SUS - CON", "SUS - RES"))
    return(rows |>
      dplyr::mutate(
        dataset = "microglia", endpoint_type = endpoint_row$endpoint_type, endpoint_id = endpoint_row$endpoint_id,
        endpoint_label = endpoint_row$endpoint_label, module_id = ifelse(endpoint_row$endpoint_type == "module_eigengene", endpoint_row$endpoint_id, NA_character_),
        estimate_before = NA_real_, estimate_after = NA_real_, percent_attenuation = NA_real_, p_before = NA_real_,
        p_after = NA_real_, FDR_before = NA_real_, FDR_after = NA_real_, neuropil_covariate_beta = NA_real_,
        neuropil_covariate_p = NA_real_, matched_neuropil_covariate = NA_character_, matched_neuropil_label = NA_character_,
        neuropil_covariate_source = NA_character_, neuropil_selection_spearman = NA_real_, n_samples = nrow(dat0),
        n_matched_samples = 0L, n_animals = dplyr::n_distinct(dat0$AnimalID), model_type_before = NA_character_,
        model_type_after = NA_character_, formula_before = NA_character_, formula_after = NA_character_,
        classification = "unstable_or_unmatched", model_warning = "no matched neuron_neuropil covariate with finite variance"
      ))
  }
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
  b <- extract_contrasts(before) |> dplyr::rename(estimate_before = "estimate", SE_before = "SE", statistic_before = "statistic", p_before = "p_value")
  a <- extract_contrasts(after) |> dplyr::rename(estimate_after = "estimate", SE_after = "SE", statistic_after = "statistic", p_after = "p_value")
  dplyr::left_join(b, a, by = "contrast") |>
    dplyr::mutate(
      dataset = "microglia",
      endpoint_type = endpoint_row$endpoint_type,
      endpoint_id = endpoint_row$endpoint_id,
      endpoint_label = endpoint_row$endpoint_label,
      module_id = ifelse(endpoint_row$endpoint_type == "module_eigengene", endpoint_row$endpoint_id, NA_character_),
      percent_attenuation = dplyr::if_else(is.finite(.data$estimate_before) & abs(.data$estimate_before) > 1e-12,
                                           100 * (abs(.data$estimate_before) - abs(.data$estimate_after)) / abs(.data$estimate_before),
                                           NA_real_),
      neuropil_covariate_beta = beta$beta,
      neuropil_covariate_p = beta$p,
      matched_neuropil_covariate = cov_choice$covariate_id,
      matched_neuropil_label = cov_choice$covariate_label,
      neuropil_covariate_source = cov_choice$covariate_source,
      neuropil_selection_spearman = cov_choice$selection_spearman,
      n_samples = nrow(dat0),
      n_matched_samples = nrow(dat),
      n_animals = dplyr::n_distinct(dat$AnimalID),
      model_type_before = before$model_type,
      model_type_after = after$model_type,
      formula_before = before$formula_used %||% before$formula_requested,
      formula_after = after$formula_used %||% after$formula_requested,
      model_warning = paste(unique(c(before$warning, after$warning, cov_base$dropped, cov_adj$dropped)), collapse = ";")
    )
}

endpoint_rows <- micro_map |>
  dplyr::filter(.data$endpoint_col %in% names(micro_dat)) |>
  dplyr::distinct(.data$endpoint_col, .keep_all = TRUE)

if (!nrow(endpoint_rows)) {
  results <- tibble::tibble(
    dataset = "microglia",
    endpoint_type = NA_character_,
    endpoint_id = NA_character_,
    endpoint_label = NA_character_,
    module_id = NA_character_,
    contrast = c("RES - CON", "SUS - CON", "SUS - RES"),
    estimate_before = NA_real_,
    estimate_after = NA_real_,
    percent_attenuation = NA_real_,
    p_before = NA_real_,
    p_after = NA_real_,
    neuropil_covariate_beta = NA_real_,
    neuropil_covariate_p = NA_real_,
    classification = "unstable_or_unmatched",
    model_warning = "no microglia module or targeted-signature endpoints available"
  )
} else {
  results <- dplyr::bind_rows(lapply(seq_len(nrow(endpoint_rows)), function(i) fit_endpoint(endpoint_rows[i, , drop = FALSE])))
}
required_output_cols <- c(
  "dataset", "endpoint_type", "endpoint_id", "endpoint_label", "module_id", "contrast",
  "estimate_before", "estimate_after", "percent_attenuation", "p_before", "p_after",
  "FDR_before", "FDR_after", "neuropil_covariate_beta", "neuropil_covariate_p",
  "matched_neuropil_covariate", "matched_neuropil_label", "neuropil_covariate_source",
  "neuropil_selection_spearman", "n_samples", "n_matched_samples", "n_animals",
  "model_type_before", "model_type_after", "formula_before", "formula_after",
  "classification", "model_warning"
)
for (nm in required_output_cols) {
  if (!nm %in% names(results)) {
    results[[nm]] <- if (nm %in% c("n_samples", "n_matched_samples", "n_animals")) {
      NA_integer_
    } else if (nm %in% c(
      "estimate_before", "estimate_after", "percent_attenuation", "p_before", "p_after",
      "FDR_before", "FDR_after", "neuropil_covariate_beta", "neuropil_covariate_p",
      "neuropil_selection_spearman"
    )) {
      NA_real_
    } else {
      NA_character_
    }
  }
}
if (!"FDR_before" %in% names(results)) results$FDR_before <- NA_real_
if (!"FDR_after" %in% names(results)) results$FDR_after <- NA_real_
ok_before <- is.finite(results$p_before)
ok_after <- is.finite(results$p_after)
results$FDR_before[ok_before] <- stats::p.adjust(results$p_before[ok_before], method = "BH")
results$FDR_after[ok_after] <- stats::p.adjust(results$p_after[ok_after], method = "BH")

results <- results |>
  dplyr::mutate(
    classification = dplyr::case_when(
      .data$classification == "unstable_or_unmatched" ~ "unstable_or_unmatched",
      is.na(.data$estimate_before) | is.na(.data$estimate_after) | is.na(.data$matched_neuropil_covariate) ~ "unstable_or_unmatched",
      is.finite(.data$percent_attenuation) & .data$percent_attenuation >= 60 & (is.na(.data$p_after) | .data$p_after >= 0.05) ~ "neuropil_explained",
      is.finite(.data$percent_attenuation) & .data$percent_attenuation >= 25 ~ "partially_neuropil_adjusted",
      !is.na(.data$p_after) & .data$p_after < 0.05 & (is.na(.data$percent_attenuation) | .data$percent_attenuation < 25) ~ "neuropil_independent",
      is.finite(.data$percent_attenuation) & .data$percent_attenuation < 25 ~ "neuropil_independent",
      TRUE ~ "unstable_or_unmatched"
    )
  ) |>
  dplyr::relocate(
    "dataset", "endpoint_type", "endpoint_id", "endpoint_label", "module_id", "contrast",
    "estimate_before", "estimate_after", "percent_attenuation", "p_before", "p_after",
    "FDR_before", "FDR_after", "neuropil_covariate_beta", "neuropil_covariate_p",
    "classification"
  )

module_classification <- results |>
  dplyr::filter(.data$endpoint_type == "module_eigengene") |>
  dplyr::group_by(.data$dataset, .data$module_id, .data$endpoint_id, .data$endpoint_label) |>
  dplyr::summarise(
    neuropil_independence_classification = dplyr::case_when(
      any(.data$classification == "neuropil_independent", na.rm = TRUE) ~ "neuropil_independent",
      any(.data$classification == "partially_neuropil_adjusted", na.rm = TRUE) ~ "partially_neuropil_adjusted",
      any(.data$classification == "neuropil_explained", na.rm = TRUE) ~ "neuropil_explained",
      TRUE ~ "unstable_or_unmatched"
    ),
    best_contrast = .data$contrast[which.min(dplyr::coalesce(.data$p_before, Inf))][[1]] %||% NA_character_,
    min_p_before = suppressWarnings(min(.data$p_before, na.rm = TRUE)),
    min_p_after = suppressWarnings(min(.data$p_after, na.rm = TRUE)),
    max_percent_attenuation = suppressWarnings(max(.data$percent_attenuation, na.rm = TRUE)),
    matched_neuropil_covariate = paste(unique(stats::na.omit(.data$matched_neuropil_covariate)), collapse = ";"),
    neuropil_covariate_source = paste(unique(stats::na.omit(.data$neuropil_covariate_source)), collapse = ";"),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    dplyr::across(c("min_p_before", "min_p_after", "max_percent_attenuation"), ~ ifelse(is.infinite(.x), NA_real_, .x)),
    neuropil_independence_note = "Classification from paired baseline versus region-matched neuron_neuropil-adjusted mixed models."
  )

signature_results <- results |> dplyr::filter(.data$endpoint_type == "targeted_signature_score")

write_table_and_source(results, PATHS$tables, PATHS$source_data, "microglia_neuropil_independence_effects.csv")
write_table_and_source(module_classification, PATHS$tables, PATHS$source_data, "microglia_module_neuropil_independence_classification.csv")
write_table_and_source(signature_results, PATHS$tables, PATHS$source_data, "microglia_targeted_signature_neuropil_independence_effects.csv")

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
      targeted_signatures = signature_results
    ),
    file.path(PATHS$tables, "microglia_neuropil_independence.xlsx")
  )
}

write_run_manifest(
  file.path(PATHS$logs, "run_manifest.yml"),
  inputs = inputs,
  outputs = list(tables = PATHS$tables, source_data = PATHS$source_data),
  parameters = list(
    dataset = "microglia",
    baseline_model = "score ~ StressGroup + SpatialLabel + Sex + Batch + (1 | AnimalID)",
    adjusted_model = "score ~ StressGroup + matched_neuropil_score + SpatialLabel + Sex + Batch + (1 | AnimalID)",
    neuropil_matching = "neuron_neuropil layer-level values aggregated to AnimalID + Region; covariate selected by strongest absolute Spearman correlation per endpoint"
  ),
  notes = "Additive downstream sensitivity analysis. Primary group-effect outputs are not overwritten."
)

message("Microglia neuropil-independence sensitivity complete.")

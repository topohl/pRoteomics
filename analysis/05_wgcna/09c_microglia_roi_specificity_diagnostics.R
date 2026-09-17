#!/usr/bin/env Rscript
# ================================================================
# Script: 06_modules_WGCNA/09c_microglia_roi_specificity_diagnostics.R
# Stage: modules_downstream
# Scope: dataset_specific
# Consumes: Existing WGCNA states, module definitions, annotation, and marker traits.
# Produces: Contrast-free microglia ROI specificity diagnostics.
# Notes: No StressGroup contrast testing, no model refitting from 09, no WGCNA recomputation.
# ================================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "script_runtime.R"))
source(repo_path("R", "wgcna_downstream_utils.R"))
source(repo_path("R", "plotting_nature.R"))

SCRIPT_ID <- "06_modules_WGCNA/09c_microglia_roi_specificity_diagnostics.R"
runtime <- init_script_runtime(SCRIPT_ID, stage = "modules_downstream", default_dataset = "microglia")
if (!identical(runtime$dataset, "microglia") && !isTRUE(runtime$dry_run)) {
  stop("This diagnostic is microglia-only. Use --dataset microglia.", call. = FALSE)
}

bool_arg <- function(flag, default = FALSE, args = runtime$args) {
  value <- script_arg_value(flag, if (isTRUE(default)) "TRUE" else "FALSE", args = args)
  value <- tolower(as.character(value))
  if (value %in% c("1", "true", "t", "yes", "y")) return(TRUE)
  if (value %in% c("0", "false", "f", "no", "n")) return(FALSE)
  stop("Expected TRUE/FALSE for ", flag, ", got: ", value, call. = FALSE)
}

TARGET_MODULE <- script_arg_value("--target-module", "WGCNA_#4D4D4D", args = runtime$args)
METRIC <- script_arg_value("--metric", "spearman", args = runtime$args)
if (!METRIC %in% c("spearman", "delta_r2")) stop("--metric must be spearman or delta_r2.", call. = FALSE)
INCLUDE_EXPLORATORY <- bool_arg("--include-exploratory", FALSE, args = runtime$args)

thresholds <- list(
  high_abs_spearman = 0.5,
  moderate_abs_spearman = 0.30,
  high_jaccard = 0.15,
  high_overlap_fraction_microglia = 0.50,
  moderate_overlap_fraction_microglia = 0.25,
  low_abs_spearman = 0.25,
  low_jaccard = 0.05
)

required_pkgs <- c("dplyr", "tidyr", "tibble", "readr", "ggplot2", "yaml", "scales", "grid")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) && !isTRUE(runtime$dry_run)) {
  stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}
if (!length(missing_pkgs)) suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))

PATHS <- create_module_dirs("06_modules_WGCNA", file.path("microglia_roi_specificity", "microglia"))
REPORT_DIR <- path_results("reports", "06_modules_WGCNA", "microglia_roi_specificity", "microglia")
FILES_MICRO <- resolve_wgcna_files("microglia")
FILES_NEUROPIL <- resolve_wgcna_files("neuron_neuropil")
FILES_SOMA <- resolve_wgcna_files("neuron_soma")
CONFIG_FILE <- repo_path("config", "microglia_neuropil_independence.yml")
ANNOTATION_FILE <- path_results("tables", "06_modules_WGCNA", "module_annotation", "microglia", "WGCNA_module_biological_annotation.csv")
TARGETED_FILE <- path_results("tables", "06_modules_WGCNA", "module_annotation", "microglia", "WGCNA_module_targeted_signature_overlap_details.csv")
LABEL_LOOKUP_FILE <- path_results("tables", "06_modules_WGCNA", "interpretable_summary", "microglia", "WGCNA_final_label_lookup.csv")
soma_available <- file.exists(FILES_SOMA$state) && file.exists(FILES_SOMA$definitions)
INCLUDE_SOMA <- bool_arg("--include-soma", FALSE, args = runtime$args) && soma_available

inputs <- list(
  config = CONFIG_FILE,
  microglia_state = FILES_MICRO$state,
  microglia_definitions = FILES_MICRO$definitions,
  neuron_neuropil_state = FILES_NEUROPIL$state,
  neuron_neuropil_definitions = FILES_NEUROPIL$definitions,
  microglia_annotation = ANNOTATION_FILE,
  neuron_soma_state_optional = FILES_SOMA$state,
  neuron_soma_definitions_optional = FILES_SOMA$definitions,
  microglia_marker_traits_optional = FILES_MICRO$marker_traits,
  neuron_neuropil_marker_traits_optional = FILES_NEUROPIL$marker_traits,
  neuron_soma_marker_traits_optional = FILES_SOMA$marker_traits,
  targeted_signature_overlap_optional = TARGETED_FILE,
  final_label_lookup_optional = LABEL_LOOKUP_FILE
)
required_inputs <- inputs[c("config", "microglia_state", "microglia_definitions", "neuron_neuropil_state", "neuron_neuropil_definitions", "microglia_annotation")]

if (isTRUE(runtime$dry_run)) {
  dry_run_line("Script", SCRIPT_ID)
  dry_run_line("Dataset", runtime$dataset)
  dry_run_line("Target module", TARGET_MODULE)
  dry_run_line("Metric", METRIC)
  dry_run_line("Include neuron_soma", INCLUDE_SOMA)
  dry_run_line("Neuron_soma files available", soma_available)
  dry_run_line("Include exploratory in default figures", INCLUDE_EXPLORATORY)
  for (nm in names(required_inputs)) dry_run_line(nm, required_inputs[[nm]], if (file.exists(required_inputs[[nm]])) "PASS" else "FAIL")
  for (nm in setdiff(names(inputs), names(required_inputs))) dry_run_line(nm, inputs[[nm]], if (file.exists(inputs[[nm]])) "PASS" else "WARN")
  dry_run_line("Output tables", PATHS$tables)
  dry_run_line("Output source data", PATHS$source_data)
  dry_run_line("Output figures", PATHS$figures)
  dry_run_line("Output reports", REPORT_DIR)
  dry_run_line("Output logs", PATHS$logs)
  quit(status = if (all(file.exists(unlist(required_inputs, use.names = FALSE)))) 0 else 1, save = "no")
}

missing_required <- required_inputs[!file.exists(unlist(required_inputs, use.names = FALSE))]
if (length(missing_required)) {
  stop("Missing required input(s):\n", paste0(" - ", names(missing_required), ": ", unname(missing_required), collapse = "\n"), call. = FALSE)
}

read_csv_optional <- function(path) safe_read_csv(path) %||% data.frame()
clean_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(trimws(x))] <- NA_character_
  x
}
as_num <- function(x) suppressWarnings(as.numeric(x))
safe_cor_test <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4L || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) return(c(r = NA_real_, p = NA_real_))
  ct <- suppressWarnings(tryCatch(stats::cor.test(x[ok], y[ok], method = "spearman", exact = FALSE), error = function(e) NULL))
  if (is.null(ct)) c(r = NA_real_, p = NA_real_) else c(r = unname(ct$estimate), p = ct$p.value)
}

config <- yaml::read_yaml(CONFIG_FILE)
predeclared_families <- config$predeclared_adjustments$families %||% list()
family_order <- vapply(predeclared_families, function(fam) as.character(fam$covariate_family %||% NA_character_), character(1))
family_order <- family_order[family_order %in% c(
  "global_neuropil_score", "mitochondrial_neuropil_score", "synaptic_neuropil_score",
  "ECM_adhesion_neuropil_score", "RNA_translation_neuropil_score", "cytoskeletal_neuropil_score"
)]
family_display <- c(
  global_neuropil_score = "Global",
  mitochondrial_neuropil_score = "Mito",
  synaptic_neuropil_score = "Synaptic",
  ECM_adhesion_neuropil_score = "ECM",
  RNA_translation_neuropil_score = "RNA",
  cytoskeletal_neuropil_score = "Cytosk.",
  exploratory_best_spearman = "Best Spearman"
)

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
      module_id = dplyr::coalesce(clean_chr(.data$ModuleID), clean_chr(.data$ModuleColor.x), clean_chr(.data$ModuleColor.y), clean_chr(.data$module_eigengene)),
      endpoint_label = dplyr::coalesce(clean_chr(.data$module_display_label), clean_chr(.data$ModuleLabel_Final), .data$module_id)
    ) |>
    dplyr::distinct(.data$endpoint_col, .keep_all = TRUE)
}

metadata_with_scores <- function(state, dataset, eig) {
  standardize_wgcna_metadata(state$sample_info, dataset) |>
    dplyr::inner_join(eig, by = "Sample")
}

finite_covariates <- function(dat, candidates, min_n = 4L) {
  candidates <- candidates[candidates %in% names(dat)]
  candidates[vapply(candidates, function(cn) {
    vals <- as_num(dat[[cn]])
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

marker_endpoint_candidates <- function(traits) {
  if (!nrow(traits)) return(tibble::tibble(endpoint_col = character(), endpoint_label = character(), covariate_source = character()))
  score_cols <- grep("(^z_|^raw_|microglia_minus|microglia_to).*score$|ratio$", names(traits), value = TRUE)
  score_cols <- score_cols[vapply(traits[score_cols], function(x) is.numeric(x) || is.integer(x), logical(1))]
  tibble::tibble(endpoint_col = score_cols, endpoint_label = score_cols, covariate_source = "marker_trait")
}

prepare_compartment <- function(dataset, files, include_traits = TRUE) {
  state <- load_wgcna_state(files$state)
  eig <- extract_module_eigengenes(state)
  defs <- read_csv_optional(files$definitions)
  map <- make_module_map(eig, defs) |>
    dplyr::mutate(covariate_source = "module_eigengene")
  dat <- metadata_with_scores(state, dataset, eig)
  traits <- if (isTRUE(include_traits)) read_csv_optional(files$marker_traits) else data.frame()
  trait_map <- marker_endpoint_candidates(traits)
  if (nrow(traits) && nrow(trait_map) && "Sample" %in% names(traits)) {
    keep <- unique(c("Sample", trait_map$endpoint_col))
    dat <- dat |> dplyr::left_join(traits[, intersect(keep, names(traits)), drop = FALSE], by = "Sample")
    map <- dplyr::bind_rows(
      map,
      trait_map |>
        dplyr::transmute(endpoint_col = .data$endpoint_col, module_eigengene = NA_character_, ModuleColor = NA_character_, module_id = .data$endpoint_col, endpoint_label = .data$endpoint_label, covariate_source = .data$covariate_source)
    )
  }
  list(state = state, eig = eig, defs = defs, map = map, dat = dat, traits = traits)
}

aggregate_compartment <- function(comp) {
  cov_cols <- intersect(comp$map$endpoint_col, names(comp$dat))
  comp$dat |>
    dplyr::mutate(
      AnimalID = as.character(.data$AnimalID),
      Region = as.character(.data$Region)
    ) |>
    dplyr::filter(!is.na(.data$AnimalID), nzchar(.data$AnimalID), !is.na(.data$Region), nzchar(.data$Region)) |>
    dplyr::group_by(.data$AnimalID, .data$Region) |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(cov_cols), ~ mean(as_num(.x), na.rm = TRUE)),
      n_cross_compartment_region_samples = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(dplyr::across(dplyr::all_of(cov_cols), ~ ifelse(is.nan(.x), NA_real_, .x)))
}

valid_model_terms <- function(dat, terms) {
  terms <- terms[terms %in% names(dat)]
  terms[vapply(terms, function(v) dplyr::n_distinct(dat[[v]][!is.na(dat[[v]])]) >= 2L, logical(1))]
}

has_repeats <- function(dat) {
  "AnimalID" %in% names(dat) &&
    any(!is.na(dat$AnimalID) & nzchar(as.character(dat$AnimalID))) &&
    any(table(dat$AnimalID[!is.na(dat$AnimalID) & nzchar(as.character(dat$AnimalID))]) > 1L)
}

fit_coupling_model <- function(dat) {
  dat <- dat |>
    dplyr::mutate(module_score = as_num(.data$module_score), matched_cross_compartment_score = as_num(.data$matched_cross_compartment_score)) |>
    dplyr::filter(is.finite(.data$module_score), is.finite(.data$matched_cross_compartment_score))
  if (nrow(dat) < 4L) {
    return(list(beta = NA_real_, p = NA_real_, model_type = "not_fit", formula_used = NA_character_, delta_adjusted_r2 = NA_real_, warning = "fewer than four matched rows"))
  }
  covars <- valid_model_terms(dat, c("SpatialLabel", "Sex", "Batch"))
  rhs_full <- paste(c("matched_cross_compartment_score", covars), collapse = " + ")
  rhs_base <- if (length(covars)) paste(covars, collapse = " + ") else "1"
  use_lmer <- has_repeats(dat) && requireNamespace("lmerTest", quietly = TRUE)
  f_full <- stats::as.formula(paste0("module_score ~ ", rhs_full, if (use_lmer) " + (1 | AnimalID)" else ""))
  warn <- character()
  fit <- tryCatch(
    if (use_lmer) lmerTest::lmer(f_full, data = dat, REML = FALSE) else stats::lm(f_full, data = dat),
    warning = function(w) {
      warn <<- c(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(beta = NA_real_, p = NA_real_, model_type = if (use_lmer) "lmerTest_lmer_failed" else "lm_failed", formula_used = paste(deparse(f_full), collapse = ""), delta_adjusted_r2 = NA_real_, warning = conditionMessage(fit)))
  }
  cf <- tryCatch(as.data.frame(summary(fit)$coefficients), error = function(e) data.frame())
  p_col <- intersect(c("Pr(>|t|)", "Pr(>|z|)", "Pr(>|F|)"), names(cf))[1]
  beta <- if ("matched_cross_compartment_score" %in% rownames(cf)) as.numeric(cf["matched_cross_compartment_score", "Estimate"]) else NA_real_
  p <- if (!is.na(p_col) && "matched_cross_compartment_score" %in% rownames(cf)) as.numeric(cf["matched_cross_compartment_score", p_col]) else NA_real_
  delta <- NA_real_
  if (!use_lmer) {
    f_base <- stats::as.formula(paste0("module_score ~ ", rhs_base))
    fit_base <- tryCatch(stats::lm(f_base, data = dat), error = function(e) NULL)
    if (!is.null(fit_base)) {
      delta <- summary(fit)$adj.r.squared - summary(fit_base)$adj.r.squared
    }
  }
  list(beta = beta, p = p, model_type = if (use_lmer) "lmerTest_lmer" else "lm", formula_used = paste(deparse(stats::formula(fit)), collapse = ""), delta_adjusted_r2 = delta, warning = paste(warn, collapse = "; "))
}

choose_family_covariate <- function(comp_dat, comp_map, family) {
  available <- finite_covariates(comp_dat, comp_map$endpoint_col)
  hits <- family_covariates(family, available)
  if (!length(hits)) return(NULL)
  cov_col <- hits[[1]]
  meta <- comp_map[match(cov_col, comp_map$endpoint_col), , drop = FALSE]
  list(
    covariate_col = cov_col,
    covariate_family = as.character(family$covariate_family %||% NA_character_),
    covariate_label = meta$endpoint_label[[1]] %||% cov_col,
    covariate_source = meta$covariate_source[[1]] %||% "cross_compartment_covariate",
    adjustment_mode = if (identical(as.character(family$covariate_family %||% ""), "exploratory_best_spearman")) "exploratory_best_spearman" else "predeclared"
  )
}

module_protein_sets <- function(defs) {
  if (!nrow(defs) || !"ModuleID" %in% names(defs)) return(list())
  protein_col <- intersect(c("GeneSymbol", "UniProt", "ProteinID", "original_token"), names(defs))[1]
  if (is.na(protein_col)) return(list())
  defs <- defs |>
    dplyr::mutate(protein_token = normalize_gene_token(.data[[protein_col]])) |>
    dplyr::filter(!is.na(.data$ModuleID), nzchar(.data$ModuleID), !is.na(.data$protein_token), nzchar(.data$protein_token))
  split(defs$protein_token, defs$ModuleID) |>
    lapply(unique)
}

overlap_one_compartment <- function(micro_sets, other_sets, compartment) {
  dplyr::bind_rows(lapply(names(micro_sets), function(mid) {
    mset <- micro_sets[[mid]]
    dplyr::bind_rows(lapply(names(other_sets), function(oid) {
      oset <- other_sets[[oid]]
      ov <- intersect(mset, oset)
      union_n <- length(union(mset, oset))
      tibble::tibble(
        module_id = mid,
        other_compartment = compartment,
        other_module_id = oid,
        n_microglia_module_proteins = length(mset),
        n_other_module_proteins = length(oset),
        n_overlap = length(ov),
        jaccard = if (union_n > 0) length(ov) / union_n else NA_real_,
        overlap_fraction_microglia = if (length(mset)) length(ov) / length(mset) else NA_real_,
        overlap_fraction_other = if (length(oset)) length(ov) / length(oset) else NA_real_,
        overlap_proteins = paste(ov, collapse = ";")
      )
    }))
  }))
}

compact_module_label <- function(module_id, label) {
  lookup <- c(
    `WGCNA_#4D4D4D` = "#4D4D4D \u00b7 ECM/adhesion",
    `WGCNA_#3F007D` = "#3F007D \u00b7 mito/energy",
    `WGCNA_#7F2704` = "#7F2704 \u00b7 Acetyl-CoA",
    `WGCNA_#006D2C` = "#006D2C \u00b7 RNA/RNP",
    `WGCNA_#1F4E79` = "#1F4E79 \u00b7 mito gene expr.",
    `WGCNA_#08519C` = "#08519C \u00b7 syn/cytosk.",
    `WGCNA_#2B8CBE` = "#2B8CBE \u00b7 vesicle",
    `WGCNA_#737373` = "#737373 \u00b7 barrier",
    `WGCNA_#8C510A` = "#8C510A \u00b7 proteostasis",
    `WGCNA_#969696` = "#969696 \u00b7 actin/cytosk.",
    `WGCNA_#9E9AC8` = "#9E9AC8 \u00b7 translation",
    `WGCNA_#A6611A` = "#A6611A \u00b7 endomembrane",
    `WGCNA_#41AB5D` = "#41AB5D \u00b7 ion homeostasis"
  )
  out <- unname(lookup[as.character(module_id)])
  missing <- is.na(out) | !nzchar(out)
  fallback <- paste(sub("^WGCNA_", "", module_id), "\u00b7", label)
  out[missing] <- fallback[missing]
  out
}

state_micro <- load_wgcna_state(FILES_MICRO$state)
micro_eig <- extract_module_eigengenes(state_micro)
micro_defs <- read_csv_optional(FILES_MICRO$definitions)
micro_map <- make_module_map(micro_eig, micro_defs)
micro_dat <- metadata_with_scores(state_micro, "microglia", micro_eig)

comp_neuropil <- prepare_compartment("neuron_neuropil", FILES_NEUROPIL)
compartments <- list(neuron_neuropil = comp_neuropil)
if (isTRUE(INCLUDE_SOMA)) compartments$neuron_soma <- prepare_compartment("neuron_soma", FILES_SOMA)

coupling_rows <- list()
for (compartment_name in names(compartments)) {
  comp <- compartments[[compartment_name]]
  comp_region <- aggregate_compartment(comp)
  for (i in seq_len(nrow(micro_map))) {
    module_row <- micro_map[i, , drop = FALSE]
    endpoint_col <- module_row$endpoint_col[[1]]
    dat0 <- micro_dat |>
      dplyr::mutate(
        module_score = as_num(.data[[endpoint_col]]),
        AnimalID = as.character(.data$AnimalID),
        Region = as.character(.data$Region)
      ) |>
      dplyr::filter(is.finite(.data$module_score), !is.na(.data$AnimalID), nzchar(.data$AnimalID), !is.na(.data$Region), nzchar(.data$Region)) |>
      dplyr::left_join(comp_region, by = c("AnimalID", "Region"))

    family_choices <- lapply(predeclared_families, function(fam) choose_family_covariate(comp_region, comp$map, fam))
    family_choices <- Filter(Negate(is.null), family_choices)
    candidates <- finite_covariates(dat0, comp$map$endpoint_col)
    if (length(candidates)) {
      cors <- vapply(candidates, function(cn) safe_cor_test(dat0$module_score, as_num(dat0[[cn]]))[["r"]], numeric(1))
      if (any(is.finite(cors))) {
        best <- names(which.max(abs(cors)))
        meta <- comp$map[match(best, comp$map$endpoint_col), , drop = FALSE]
        family_choices[[length(family_choices) + 1L]] <- list(
          covariate_col = best,
          covariate_family = "exploratory_best_spearman",
          covariate_label = meta$endpoint_label[[1]] %||% best,
          covariate_source = meta$covariate_source[[1]] %||% "cross_compartment_covariate",
          adjustment_mode = "exploratory_best_spearman"
        )
      }
    }
    if (!length(family_choices)) next
    for (choice in family_choices) {
      cov_col <- choice$covariate_col
      dat_fit <- dat0 |>
        dplyr::mutate(matched_cross_compartment_score = as_num(.data[[cov_col]])) |>
        dplyr::filter(is.finite(.data$matched_cross_compartment_score))
      sp <- safe_cor_test(dat_fit$module_score, dat_fit$matched_cross_compartment_score)
      fit <- fit_coupling_model(dat_fit)
      group_counts <- if ("StressGroup" %in% names(dat_fit)) table(dat_fit$StressGroup[!is.na(dat_fit$StressGroup) & nzchar(dat_fit$StressGroup)]) else integer()
      coupling_rows[[length(coupling_rows) + 1L]] <- tibble::tibble(
        dataset = "microglia",
        module_id = module_row$module_id,
        endpoint_label = module_row$endpoint_label,
        covariate_source_compartment = compartment_name,
        covariate_family = choice$covariate_family,
        adjustment_mode = choice$adjustment_mode,
        covariate_col = cov_col,
        covariate_label = choice$covariate_label,
        covariate_source = choice$covariate_source,
        matched_by = "AnimalID + Region",
        n_matched_samples = nrow(dat_fit),
        n_matched_animals = dplyr::n_distinct(dat_fit$AnimalID[!is.na(dat_fit$AnimalID)]),
        min_animals_per_group = if (length(group_counts)) min(as.integer(group_counts)) else NA_integer_,
        matched_covariate_beta = fit$beta,
        matched_covariate_p = fit$p,
        matched_covariate_spearman_r = sp[["r"]],
        matched_covariate_spearman_p = sp[["p"]],
        abs_spearman = abs(sp[["r"]]),
        delta_adjusted_r2 = fit$delta_adjusted_r2,
        model_type = fit$model_type,
        formula_used = fit$formula_used,
        model_warning = fit$warning,
        uses_stressgroup_in_model = FALSE
      )
    }
  }
}
coupling_long <- dplyr::bind_rows(coupling_rows)

micro_sets <- module_protein_sets(micro_defs)
overlap_neuropil <- overlap_one_compartment(micro_sets, module_protein_sets(comp_neuropil$defs), "neuron_neuropil")
overlap_soma <- if (isTRUE(INCLUDE_SOMA)) overlap_one_compartment(micro_sets, module_protein_sets(compartments$neuron_soma$defs), "neuron_soma") else tibble::tibble()
overlap_long <- dplyr::bind_rows(overlap_neuropil, overlap_soma)

best_jaccard_overlap <- overlap_long |>
  dplyr::group_by(.data$module_id, .data$other_compartment) |>
  dplyr::slice_max(order_by = .data$jaccard, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::select(
    "module_id", "other_compartment", "other_module_id", "n_other_module_proteins",
    "n_overlap", "jaccard", "overlap_fraction_microglia", "overlap_fraction_other"
  ) |>
  tidyr::pivot_wider(
    names_from = "other_compartment",
    values_from = c(
      "other_module_id", "n_other_module_proteins", "n_overlap",
      "jaccard", "overlap_fraction_microglia", "overlap_fraction_other"
    ),
    names_glue = "{other_compartment}_best_jaccard_{.value}"
  )

best_microglia_fraction_overlap <- overlap_long |>
  dplyr::group_by(.data$module_id, .data$other_compartment) |>
  dplyr::slice_max(order_by = .data$overlap_fraction_microglia, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::select(
    "module_id", "other_compartment", "other_module_id", "n_other_module_proteins",
    "n_overlap", "jaccard", "overlap_fraction_microglia", "overlap_fraction_other"
  ) |>
  tidyr::pivot_wider(
    names_from = "other_compartment",
    values_from = c(
      "other_module_id", "n_other_module_proteins", "n_overlap",
      "jaccard", "overlap_fraction_microglia", "overlap_fraction_other"
    ),
    names_glue = "{other_compartment}_best_microglia_fraction_{.value}"
  )

annotation <- read_csv_optional(ANNOTATION_FILE)
targeted <- read_csv_optional(TARGETED_FILE)
label_lookup <- read_csv_optional(LABEL_LOOKUP_FILE)

annotation_keep <- c(
  "ModuleID", "cleaned_biological_label", "safe_display_label", "module_biological_label",
  "cleaned_biological_label_short", "module_biological_label_short", "microenvironment_class",
  "microenvironment_label", "annotation_confidence", "label_confidence", "marker_panels_supporting",
  "targeted_signature_primary_driver", "targeted_signature_driver_class",
  "targeted_signature_driver_signature", "targeted_signature_driver_padj",
  "targeted_signature_driver_overlap_proteins", "canonical_microglia_evidence",
  "empirical_microglia_roi_evidence", "canonical_neuropil_evidence",
  "empirical_neuropil_evidence", "empirical_shared_microenvironment_evidence",
  "microglia_state_or_activation_evidence", "interpretation_note"
)
annotation_slim <- annotation[, intersect(annotation_keep, names(annotation)), drop = FALSE] |>
  dplyr::distinct()

if (nrow(label_lookup) && all(c("entity_id", "level", "final_plot_label") %in% names(label_lookup))) {
  label_slim <- label_lookup |>
    dplyr::filter(.data$level == "module") |>
    dplyr::transmute(ModuleID = .data$entity_id, final_plot_label = .data$final_plot_label) |>
    dplyr::distinct()
  annotation_slim <- annotation_slim |> dplyr::left_join(label_slim, by = "ModuleID")
}

micro_base <- micro_map |>
  dplyr::transmute(module_id = .data$module_id, endpoint_label = .data$endpoint_label) |>
  dplyr::distinct()

best_coupling <- coupling_long |>
  dplyr::filter(.data$adjustment_mode != "exploratory_best_spearman") |>
  dplyr::group_by(.data$module_id, .data$covariate_source_compartment) |>
  dplyr::slice_max(order_by = .data$abs_spearman, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::select(
    "module_id", "covariate_source_compartment", "abs_spearman",
    "covariate_family", "covariate_label", "matched_covariate_spearman_p"
  ) |>
  tidyr::pivot_wider(
    names_from = "covariate_source_compartment",
    values_from = c("abs_spearman", "covariate_family", "covariate_label", "matched_covariate_spearman_p"),
    names_glue = "{covariate_source_compartment}_{.value}"
  )

marker_support_score <- function(row) {
  vals <- as_num(unlist(row[intersect(c("canonical_microglia_evidence", "empirical_microglia_roi_evidence", "microglia_state_or_activation_evidence"), names(row))]))
  any(is.finite(vals) & vals > 0, na.rm = TRUE)
}

shared_local_classes <- c(
  "shared_microenvironment", "vascular_basement_membrane_ecm", "astrocyte_endfoot",
  "astrocyte_or_endfoot_sensitive", "oligodendrocyte_myelin",
  "oligodendrocyte_or_myelin_sensitive", "neuropil_sensitive", "vascular_bbb_mural"
)

max_finite <- function(...) {
  vals <- as_num(unlist(list(...)))
  vals <- vals[is.finite(vals)]
  if (!length(vals)) NA_real_ else max(vals)
}

ensure_columns <- function(dat, cols, value = NA) {
  for (col in cols) {
    if (!col %in% names(dat)) dat[[col]] <- value
  }
  dat
}

fmt_fraction_label <- function(n_overlap, denom, frac) {
  n_overlap <- as_num(n_overlap)
  denom <- as_num(denom)
  frac <- as_num(frac)
  if (!is.finite(n_overlap) || !is.finite(denom) || !is.finite(frac)) return(NA_character_)
  paste0(as.integer(n_overlap), "/", as.integer(denom), " (", round(100 * frac), "%)")
}

source_of_max <- function(neuropil_value, soma_value, neuropil_source = "neuron_neuropil", soma_source = "neuron_soma") {
  neuropil_value <- as_num(neuropil_value)
  soma_value <- as_num(soma_value)
  if (!is.finite(neuropil_value) && !is.finite(soma_value)) return(NA_character_)
  if (!is.finite(soma_value) || (is.finite(neuropil_value) && neuropil_value >= soma_value)) return(neuropil_source)
  soma_source
}

classify_specificity <- function(row) {
  cls <- as.character(row$microenvironment_class %||% NA_character_)
  best_neuro <- as_num(row$best_neuropil_abs_spearman)
  best_soma <- as_num(row$best_soma_abs_spearman)
  jac_neuro <- as_num(row$best_jaccard_matching_neuropil_jaccard)
  jac_soma <- as_num(row$best_jaccard_matching_soma_jaccard)
  frac_neuro <- as_num(row$best_matching_neuropil_overlap_fraction_microglia)
  frac_soma <- as_num(row$best_matching_soma_overlap_fraction_microglia)
  coupling_vals <- c(best_neuro, best_soma)
  overlap_vals <- c(frac_neuro, frac_soma)
  jaccard_vals <- c(jac_neuro, jac_soma)
  strong_coupled <- any(coupling_vals >= thresholds$high_abs_spearman, na.rm = TRUE) ||
    any(jaccard_vals >= thresholds$high_jaccard, na.rm = TRUE) ||
    any(overlap_vals >= thresholds$high_overlap_fraction_microglia, na.rm = TRUE)
  low_coupled <- all(coupling_vals[is.finite(coupling_vals)] <= thresholds$low_abs_spearman) &&
    all(jaccard_vals[is.finite(jaccard_vals)] <= thresholds$low_jaccard) &&
    all(overlap_vals[is.finite(overlap_vals)] < thresholds$moderate_overlap_fraction_microglia)
  micro_support <- marker_support_score(row)
  if (!is.na(cls) && cls %in% shared_local_classes) return("local_microenvironment_or_shared")
  if (isTRUE(strong_coupled)) return("neuropil_coupled_or_shared")
  if (isTRUE(micro_support) && isTRUE(low_coupled)) return("microglia_roi_enriched")
  "ambiguous_diagnostic"
}

specificity_evidence_tier <- function(row) {
  cls <- as.character(row$microenvironment_class %||% NA_character_)
  shared_local <- !is.na(cls) && cls %in% shared_local_classes
  max_spearman <- max_finite(row$best_neuropil_abs_spearman, row$best_soma_abs_spearman)
  max_micro_frac <- max_finite(
    row$best_matching_neuropil_overlap_fraction_microglia,
    row$best_matching_soma_overlap_fraction_microglia
  )
  micro_support <- marker_support_score(row)
  if (isTRUE(shared_local) &&
      is.finite(max_spearman) && max_spearman >= thresholds$high_abs_spearman) {
    return("shared_local_program_high_score_coupling")
  }
  if (isTRUE(shared_local) &&
      is.finite(max_spearman) && max_spearman >= thresholds$moderate_abs_spearman &&
      max_spearman < thresholds$high_abs_spearman &&
      is.finite(max_micro_frac) && max_micro_frac >= thresholds$high_overlap_fraction_microglia) {
    return("shared_local_program_moderate_score_coupling_high_protein_sharing")
  }
  if (isTRUE(shared_local) &&
      is.finite(max_micro_frac) && max_micro_frac >= thresholds$high_overlap_fraction_microglia) {
    return("shared_local_program_high_protein_sharing")
  }
  if (isTRUE(micro_support) &&
      (!is.finite(max_spearman) || max_spearman < thresholds$low_abs_spearman) &&
      (!is.finite(max_micro_frac) || max_micro_frac < thresholds$moderate_overlap_fraction_microglia)) {
    return("microglia_roi_enriched_candidate")
  }
  "shared_or_ambiguous_diagnostic"
}

specificity_archetype <- function(row) {
  max_score <- as_num(row$max_score_coupling)
  max_share <- as_num(row$max_protein_sharing)
  cls <- as.character(row$microenvironment_class %||% NA_character_)
  micro_support <- marker_support_score(row)
  shared_local <- !is.na(cls) && cls %in% shared_local_classes
  if (is.finite(max_score) && is.finite(max_share) &&
      max_score >= thresholds$high_abs_spearman && max_share >= thresholds$high_overlap_fraction_microglia) {
    return("strongly_coupled_and_shared")
  }
  if (is.finite(max_score) && max_score >= thresholds$high_abs_spearman &&
      (!is.finite(max_share) || max_share < thresholds$high_overlap_fraction_microglia)) {
    return("score_coupled_shared_state")
  }
  if ((!is.finite(max_score) || max_score < thresholds$high_abs_spearman) &&
      is.finite(max_share) && max_share >= thresholds$high_overlap_fraction_microglia) {
    return("compositionally_shared_modest_coupling")
  }
  if (is.finite(max_score) && is.finite(max_share) &&
      max_score < thresholds$low_abs_spearman && max_share < thresholds$moderate_overlap_fraction_microglia &&
      isTRUE(micro_support) && !isTRUE(shared_local)) {
    return("microglia_roi_enriched_candidate")
  }
  "shared_or_ambiguous_diagnostic"
}

archetype_interpretation <- function(archetype) {
  dplyr::case_when(
    archetype == "strongly_coupled_and_shared" ~ "Module activity and protein composition are both shared across compartments; not microglia-restricted.",
    archetype == "score_coupled_shared_state" ~ "Module activity tracks matched cross-compartment scores; interpret as shared state/coupled ROI signal.",
    archetype == "compositionally_shared_modest_coupling" ~ "Module activity is not strongly coupled, but most proteins are shared across compartments; interpret as local/shared program.",
    archetype == "microglia_roi_enriched_candidate" ~ "Low cross-compartment coupling and low protein sharing with microglia evidence; candidate microglia ROI-enriched module.",
    TRUE ~ "Mixed or incomplete evidence; keep as diagnostic and avoid purified cell-type wording."
  )
}

interpret_note <- function(row) {
  class <- row$microglia_roi_specificity_class
  cls <- as.character(row$microenvironment_class %||% NA_character_)
  if (identical(row$module_id, TARGET_MODULE) && cls %in% c("shared_microenvironment", "vascular_basement_membrane_ecm", "vascular_bbb_mural")) {
    return("ECM/adhesion local microenvironment module; not purified microglia-restricted. Score coupling is modest/moderate, but protein sharing with large neuron_neuropil/soma modules is high.")
  }
  dplyr::case_when(
    class == "local_microenvironment_or_shared" ~ "Annotation and marker context indicate a local/shared ROI microenvironment diagnostic category.",
    class == "neuropil_coupled_or_shared" ~ "Matched cross-compartment coupling or protein sharing is high; interpret as shared/coupled diagnostic signal.",
    class == "microglia_roi_enriched" ~ "Marker context and low cross-compartment coupling support a microglia ROI-enriched diagnostic category.",
    TRUE ~ "Evidence is mixed or incomplete; keep as diagnostic without purified cell-type wording."
  )
}

summary_tbl <- micro_base |>
  dplyr::left_join(annotation_slim, by = c("module_id" = "ModuleID")) |>
  dplyr::left_join(best_coupling, by = "module_id") |>
  dplyr::left_join(best_jaccard_overlap, by = "module_id") |>
  dplyr::left_join(best_microglia_fraction_overlap, by = "module_id") |>
  ensure_columns(c(
    "cleaned_biological_label_short", "module_biological_label_short", "cleaned_biological_label",
    "module_biological_label", "microenvironment_class", "microenvironment_label",
    "neuron_neuropil_abs_spearman", "neuron_neuropil_covariate_family",
    "neuron_neuropil_covariate_label", "neuron_neuropil_matched_covariate_spearman_p",
    "neuron_soma_abs_spearman", "neuron_soma_covariate_family",
    "neuron_neuropil_best_microglia_fraction_other_module_id",
    "neuron_neuropil_best_microglia_fraction_jaccard",
    "neuron_neuropil_best_microglia_fraction_n_overlap",
    "neuron_neuropil_best_microglia_fraction_n_other_module_proteins",
    "neuron_neuropil_best_microglia_fraction_overlap_fraction_microglia",
    "neuron_neuropil_best_microglia_fraction_overlap_fraction_other",
    "neuron_neuropil_best_jaccard_other_module_id",
    "neuron_neuropil_best_jaccard_jaccard",
    "neuron_neuropil_best_jaccard_overlap_fraction_microglia",
    "neuron_soma_best_microglia_fraction_other_module_id",
    "neuron_soma_best_microglia_fraction_jaccard",
    "neuron_soma_best_microglia_fraction_n_overlap",
    "neuron_soma_best_microglia_fraction_n_other_module_proteins",
    "neuron_soma_best_microglia_fraction_overlap_fraction_microglia",
    "neuron_soma_best_microglia_fraction_overlap_fraction_other",
    "neuron_soma_best_jaccard_other_module_id",
    "neuron_soma_best_jaccard_jaccard",
    "neuron_soma_best_jaccard_overlap_fraction_microglia",
    "canonical_microglia_evidence", "empirical_microglia_roi_evidence",
    "microglia_state_or_activation_evidence"
  )) |>
  dplyr::mutate(
    compact_label = compact_module_label(.data$module_id, dplyr::coalesce(clean_chr(.data$cleaned_biological_label_short), clean_chr(.data$module_biological_label_short), clean_chr(.data$endpoint_label))),
    cleaned_biological_label = dplyr::coalesce(clean_chr(.data$cleaned_biological_label), clean_chr(.data$module_biological_label), clean_chr(.data$endpoint_label)),
    best_neuropil_abs_spearman = .data$neuron_neuropil_abs_spearman,
    best_neuropil_covariate_family = .data$neuron_neuropil_covariate_family,
    best_neuropil_covariate_label = .data$neuron_neuropil_covariate_label,
    best_neuropil_p = .data$neuron_neuropil_matched_covariate_spearman_p,
    best_soma_abs_spearman = .data$neuron_soma_abs_spearman %||% NA_real_,
    best_soma_covariate_family = .data$neuron_soma_covariate_family %||% NA_character_,
    best_matching_neuropil_module = .data$neuron_neuropil_best_microglia_fraction_other_module_id,
    best_matching_neuropil_jaccard = .data$neuron_neuropil_best_microglia_fraction_jaccard,
    best_matching_neuropil_n_overlap = .data$neuron_neuropil_best_microglia_fraction_n_overlap,
    best_matching_neuropil_n_other_module_proteins = .data$neuron_neuropil_best_microglia_fraction_n_other_module_proteins,
    best_matching_neuropil_overlap_fraction_microglia = .data$neuron_neuropil_best_microglia_fraction_overlap_fraction_microglia,
    best_matching_neuropil_overlap_fraction_other = .data$neuron_neuropil_best_microglia_fraction_overlap_fraction_other,
    best_jaccard_matching_neuropil_module = .data$neuron_neuropil_best_jaccard_other_module_id,
    best_jaccard_matching_neuropil_jaccard = .data$neuron_neuropil_best_jaccard_jaccard,
    best_jaccard_matching_neuropil_overlap_fraction_microglia = .data$neuron_neuropil_best_jaccard_overlap_fraction_microglia,
    best_matching_soma_module = .data$neuron_soma_best_microglia_fraction_other_module_id %||% NA_character_,
    best_matching_soma_jaccard = .data$neuron_soma_best_microglia_fraction_jaccard %||% NA_real_,
    best_matching_soma_n_overlap = .data$neuron_soma_best_microglia_fraction_n_overlap %||% NA_real_,
    best_matching_soma_n_other_module_proteins = .data$neuron_soma_best_microglia_fraction_n_other_module_proteins %||% NA_real_,
    best_matching_soma_overlap_fraction_microglia = .data$neuron_soma_best_microglia_fraction_overlap_fraction_microglia %||% NA_real_,
    best_matching_soma_overlap_fraction_other = .data$neuron_soma_best_microglia_fraction_overlap_fraction_other %||% NA_real_,
    best_jaccard_matching_soma_module = .data$neuron_soma_best_jaccard_other_module_id %||% NA_character_,
    best_jaccard_matching_soma_jaccard = .data$neuron_soma_best_jaccard_jaccard %||% NA_real_,
    best_jaccard_matching_soma_overlap_fraction_microglia = .data$neuron_soma_best_jaccard_overlap_fraction_microglia %||% NA_real_,
    max_abs_spearman = pmax(.data$best_neuropil_abs_spearman, .data$best_soma_abs_spearman, na.rm = TRUE),
    max_overlap_fraction_microglia = pmax(.data$best_matching_neuropil_overlap_fraction_microglia, .data$best_matching_soma_overlap_fraction_microglia, na.rm = TRUE),
    max_score_coupling = .data$max_abs_spearman,
    max_protein_sharing = .data$max_overlap_fraction_microglia
  ) |>
  dplyr::mutate(
    max_abs_spearman = ifelse(is.infinite(.data$max_abs_spearman), NA_real_, .data$max_abs_spearman),
    max_overlap_fraction_microglia = ifelse(is.infinite(.data$max_overlap_fraction_microglia), NA_real_, .data$max_overlap_fraction_microglia),
    max_score_coupling = ifelse(is.infinite(.data$max_score_coupling), NA_real_, .data$max_score_coupling),
    max_protein_sharing = ifelse(is.infinite(.data$max_protein_sharing), NA_real_, .data$max_protein_sharing)
  )

summary_tbl$microglia_roi_specificity_class <- vapply(seq_len(nrow(summary_tbl)), function(i) classify_specificity(summary_tbl[i, , drop = FALSE]), character(1))
summary_tbl$specificity_evidence_tier <- vapply(seq_len(nrow(summary_tbl)), function(i) specificity_evidence_tier(summary_tbl[i, , drop = FALSE]), character(1))
summary_tbl$max_score_coupling_source <- vapply(seq_len(nrow(summary_tbl)), function(i) {
  source_of_max(summary_tbl$best_neuropil_abs_spearman[[i]], summary_tbl$best_soma_abs_spearman[[i]])
}, character(1))
summary_tbl$max_protein_sharing_source <- vapply(seq_len(nrow(summary_tbl)), function(i) {
  source_of_max(summary_tbl$best_matching_neuropil_overlap_fraction_microglia[[i]], summary_tbl$best_matching_soma_overlap_fraction_microglia[[i]])
}, character(1))
summary_tbl$best_neuropil_sharing_label <- vapply(seq_len(nrow(summary_tbl)), function(i) {
  fmt_fraction_label(
    summary_tbl$best_matching_neuropil_n_overlap[[i]],
    overlap_long$n_microglia_module_proteins[match(summary_tbl$module_id[[i]], overlap_long$module_id)],
    summary_tbl$best_matching_neuropil_overlap_fraction_microglia[[i]]
  )
}, character(1))
summary_tbl$best_soma_sharing_label <- vapply(seq_len(nrow(summary_tbl)), function(i) {
  fmt_fraction_label(
    summary_tbl$best_matching_soma_n_overlap[[i]],
    overlap_long$n_microglia_module_proteins[match(summary_tbl$module_id[[i]], overlap_long$module_id)],
    summary_tbl$best_matching_soma_overlap_fraction_microglia[[i]]
  )
}, character(1))
summary_tbl$specificity_archetype <- vapply(seq_len(nrow(summary_tbl)), function(i) specificity_archetype(summary_tbl[i, , drop = FALSE]), character(1))
summary_tbl$interpretation <- archetype_interpretation(summary_tbl$specificity_archetype)
summary_tbl$boss_interpretation <- summary_tbl$interpretation
summary_tbl$boss_interpretation[summary_tbl$module_id == "WGCNA_#4D4D4D"] <- "ECM/adhesion local microenvironment module: modest score coupling but high protein sharing; not microglia-restricted."
summary_tbl$interpretation_note <- vapply(seq_len(nrow(summary_tbl)), function(i) interpret_note(summary_tbl[i, , drop = FALSE]), character(1))

summary_tbl <- summary_tbl |>
  dplyr::select(
    "module_id", "endpoint_label", "compact_label", dplyr::any_of(c(
      "cleaned_biological_label", "safe_display_label", "module_biological_label",
      "microenvironment_class", "microenvironment_label", "annotation_confidence", "label_confidence",
      "marker_panels_supporting", "targeted_signature_primary_driver",
      "targeted_signature_driver_class", "targeted_signature_driver_signature",
      "targeted_signature_driver_padj", "targeted_signature_driver_overlap_proteins"
    )),
    "best_neuropil_abs_spearman", "best_neuropil_covariate_family", "best_neuropil_covariate_label", "best_neuropil_p",
    "best_soma_abs_spearman", "best_soma_covariate_family",
    "best_matching_neuropil_module", "best_matching_neuropil_jaccard",
    "best_matching_neuropil_n_overlap", "best_matching_neuropil_n_other_module_proteins",
    "best_matching_neuropil_overlap_fraction_microglia", "best_matching_neuropil_overlap_fraction_other",
    "best_jaccard_matching_neuropil_module", "best_jaccard_matching_neuropil_jaccard",
    "best_jaccard_matching_neuropil_overlap_fraction_microglia",
    "best_matching_soma_module", "best_matching_soma_jaccard",
    "best_matching_soma_n_overlap", "best_matching_soma_n_other_module_proteins",
    "best_matching_soma_overlap_fraction_microglia", "best_matching_soma_overlap_fraction_other",
    "best_jaccard_matching_soma_module", "best_jaccard_matching_soma_jaccard",
    "best_jaccard_matching_soma_overlap_fraction_microglia",
    "max_abs_spearman", "max_overlap_fraction_microglia",
    "max_score_coupling", "max_score_coupling_source", "max_protein_sharing", "max_protein_sharing_source",
    "best_neuropil_sharing_label", "best_soma_sharing_label",
    "microglia_roi_specificity_class", "specificity_evidence_tier", "specificity_archetype",
    "interpretation", "boss_interpretation", "interpretation_note",
    dplyr::everything()
  )

plot_source <- summary_tbl |>
  dplyr::transmute(
    module_id = .data$module_id,
    compact_label = .data$compact_label,
    microglia_roi_specificity_class = .data$microglia_roi_specificity_class,
    specificity_evidence_tier = .data$specificity_evidence_tier,
    specificity_archetype = .data$specificity_archetype,
    max_score_coupling = .data$max_score_coupling,
    max_score_coupling_source = .data$max_score_coupling_source,
    max_protein_sharing = .data$max_protein_sharing,
    max_protein_sharing_source = .data$max_protein_sharing_source,
    best_neuropil_sharing_label = .data$best_neuropil_sharing_label,
    best_soma_sharing_label = .data$best_soma_sharing_label,
    microenvironment_class = .data$microenvironment_class,
    microenvironment_label = .data$microenvironment_label,
    matched_neuropil_coupling = .data$best_neuropil_abs_spearman,
    matched_soma_coupling = .data$best_soma_abs_spearman,
    best_neuropil_module_overlap = .data$best_matching_neuropil_jaccard,
    best_soma_module_overlap = .data$best_matching_soma_jaccard,
    best_neuropil_protein_sharing = .data$best_matching_neuropil_overlap_fraction_microglia,
    best_soma_protein_sharing = .data$best_matching_soma_overlap_fraction_microglia,
    microglia_support_score = pmax(as_num(.data$canonical_microglia_evidence), as_num(.data$empirical_microglia_roi_evidence), as_num(.data$microglia_state_or_activation_evidence), na.rm = TRUE),
    shared_ecm_evidence = .data$microenvironment_class %in% shared_local_classes,
    highlight_target = .data$module_id == TARGET_MODULE,
    interpretation = .data$interpretation,
    boss_interpretation = .data$boss_interpretation,
    interpretation_note = .data$interpretation_note
  ) |>
  dplyr::mutate(
    microglia_support_score = ifelse(is.infinite(.data$microglia_support_score), NA_real_, .data$microglia_support_score),
    microglia_roi_evidence = is.finite(.data$microglia_support_score) & .data$microglia_support_score > 0
  )

write_table_and_source(summary_tbl, PATHS$tables, PATHS$source_data, "microglia_roi_specificity_summary.csv")
write_table_and_source(coupling_long, PATHS$tables, PATHS$source_data, "microglia_matched_neuropil_coupling_long.csv")
write_table_and_source(overlap_long, PATHS$tables, PATHS$source_data, "microglia_cross_compartment_module_overlap.csv")
write_table_and_source(plot_source, PATHS$tables, PATHS$source_data, "microglia_roi_specificity_source_for_plot.csv")

save_svg_pdf <- function(plot, stem, width_mm = 100, height_mm = 70) {
  svg_path <- file.path(PATHS$figures, paste0(stem, ".svg"))
  pdf_path <- file.path(PATHS$figures, paste0(stem, ".pdf"))
  svg_device <- if (requireNamespace("svglite", quietly = TRUE)) svglite::svglite else "svg"
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggplot2::ggsave(svg_path, plot = plot, width = mm_to_in(width_mm), height = mm_to_in(height_mm), units = "in", device = svg_device, limitsize = FALSE)
  ggplot2::ggsave(pdf_path, plot = plot, width = mm_to_in(width_mm), height = mm_to_in(height_mm), units = "in", device = pdf_device, limitsize = FALSE)
  c(svg_path, pdf_path)
}

theme_compact <- ggplot2::theme_classic(base_size = 7, base_family = "sans") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 8.5),
    plot.caption = ggplot2::element_text(size = 5.5, hjust = 0, colour = "#5F6568"),
    axis.text = ggplot2::element_text(size = 6, colour = "#2D3436"),
    axis.title = ggplot2::element_text(size = 7, colour = "#2D3436"),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold", size = 6.5),
    legend.title = ggplot2::element_text(size = 6),
    legend.text = ggplot2::element_text(size = 5.5),
    plot.margin = ggplot2::margin(3, 4, 3, 3, "pt")
  )

decision_labels <- c(
  `WGCNA_#4D4D4D` = "#4D4D4D ECM",
  `WGCNA_#006D2C` = "#006D2C RNA",
  `WGCNA_#08519C` = "#08519C syn/cyto",
  `WGCNA_#969696` = "#969696 actin",
  `WGCNA_#1F4E79` = "#1F4E79 mito",
  `WGCNA_#7F2704` = "#7F2704 AcCoA",
  `WGCNA_#9E9AC8` = "#9E9AC8 transl."
)

decision_source <- summary_tbl |>
  dplyr::transmute(
    module_id = .data$module_id,
    compact_label = .data$compact_label,
    decision_label = unname(decision_labels[.data$module_id]),
    neuropil_score_coupling = .data$best_neuropil_abs_spearman,
    neuropil_protein_sharing = .data$best_matching_neuropil_overlap_fraction_microglia,
    max_score_coupling = .data$max_score_coupling,
    max_score_coupling_source = .data$max_score_coupling_source,
    max_protein_sharing = .data$max_protein_sharing,
    max_protein_sharing_source = .data$max_protein_sharing_source,
    specificity_archetype = .data$specificity_archetype,
    neuropil_specificity_archetype = dplyr::case_when(
      .data$best_neuropil_abs_spearman >= thresholds$high_abs_spearman &
        .data$best_matching_neuropil_overlap_fraction_microglia >= thresholds$high_overlap_fraction_microglia ~ "Strongly neuropil-coupled/shared",
      .data$best_neuropil_abs_spearman >= thresholds$high_abs_spearman &
        .data$best_matching_neuropil_overlap_fraction_microglia < thresholds$high_overlap_fraction_microglia ~ "Neuropil score-coupled",
      .data$best_neuropil_abs_spearman < thresholds$high_abs_spearman &
        .data$best_matching_neuropil_overlap_fraction_microglia >= thresholds$high_overlap_fraction_microglia ~ "Shared proteins, modest coupling",
      .data$best_neuropil_abs_spearman < thresholds$low_abs_spearman &
        .data$best_matching_neuropil_overlap_fraction_microglia < thresholds$moderate_overlap_fraction_microglia &
        pmax(as_num(.data$canonical_microglia_evidence), as_num(.data$empirical_microglia_roi_evidence), as_num(.data$microglia_state_or_activation_evidence), na.rm = TRUE) > 0 &
        !(.data$microenvironment_class %in% shared_local_classes) ~ "Microglia ROI-enriched candidate",
      TRUE ~ "Ambiguous low-neuropil diagnostic"
    ),
    microglia_roi_specificity_class = .data$microglia_roi_specificity_class,
    best_neuropil_sharing_label = .data$best_neuropil_sharing_label,
    best_soma_sharing_label = .data$best_soma_sharing_label,
    boss_interpretation = .data$boss_interpretation,
    label_key_module = .data$module_id %in% names(decision_labels)
  ) |>
  dplyr::filter(is.finite(.data$neuropil_score_coupling), is.finite(.data$neuropil_protein_sharing))
write_csv_safe2(decision_source, file.path(PATHS$source_data, "microglia_roi_specificity_decision_map_source.csv"))

decision_palette <- c(
  `Strongly neuropil-coupled/shared` = "#7A3B46",
  `Neuropil score-coupled` = "#6D8FB3",
  `Shared proteins, modest coupling` = "#176B57",
  `Microglia ROI-enriched candidate` = "#C4A24D",
  `Ambiguous low-neuropil diagnostic` = "#8A8F93"
)

decision_plot <- ggplot2::ggplot(decision_source, ggplot2::aes(x = .data$neuropil_score_coupling, y = .data$neuropil_protein_sharing)) +
  ggplot2::geom_hline(yintercept = thresholds$high_overlap_fraction_microglia, linetype = "dashed", linewidth = 0.18, colour = "#C6CCCD") +
  ggplot2::geom_vline(xintercept = thresholds$high_abs_spearman, linetype = "dashed", linewidth = 0.18, colour = "#C6CCCD") +
  ggplot2::annotate("label", x = 0.80, y = 0.93, label = "shared/coupled", size = 1.55, colour = "#6A7072", fill = "#FAFBFB", linewidth = 0, label.padding = grid::unit(0.8, "pt")) +
  ggplot2::annotate("label", x = 0.18, y = 0.61, label = "shared proteins\nmodest coupling", size = 1.55, colour = "#6A7072", fill = "#FAFBFB", linewidth = 0, label.padding = grid::unit(0.8, "pt"), lineheight = 0.86) +
  ggplot2::annotate("label", x = 0.80, y = 0.13, label = "coupled state", size = 1.55, colour = "#6A7072", fill = "#FAFBFB", linewidth = 0, label.padding = grid::unit(0.8, "pt")) +
  ggplot2::annotate("label", x = 0.22, y = 0.14, label = "low neuropil\ncoupling/sharing\ncandidate only with\nmicroglia evidence", size = 1.45, colour = "#6A7072", fill = "#FAFBFB", linewidth = 0, label.padding = grid::unit(0.8, "pt"), lineheight = 0.84) +
  ggplot2::geom_point(ggplot2::aes(fill = .data$neuropil_specificity_archetype), shape = 21, size = 1.75, stroke = 0.22, colour = "white", alpha = 0.95) +
  ggplot2::scale_fill_manual(values = decision_palette, name = NULL) +
  ggplot2::scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), expand = ggplot2::expansion(mult = c(0.02, 0.04))) +
  ggplot2::scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), expand = ggplot2::expansion(mult = c(0.03, 0.04))) +
  ggplot2::labs(
    title = "Microglia ROI specificity map",
    x = "Matched neuropil score coupling",
    y = "Protein sharing with neuropil modules",
    caption = "x: |Spearman rho| with matched neuron-neuropil scores.\ny: fraction of module proteins in the best matching neuropil module.\nLower-left needs microglia evidence plus no shared/local annotation.\nNo StressGroup contrasts used."
  ) +
  theme_compact +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 7.6),
    axis.title = ggplot2::element_text(size = 6.2, colour = "#2D3436"),
    axis.text = ggplot2::element_text(size = 5.2, colour = "#2D3436"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.size = grid::unit(2.8, "mm"),
    legend.text = ggplot2::element_text(size = 4.7),
    legend.margin = ggplot2::margin(t = -2, b = -2, unit = "pt"),
    legend.box.margin = ggplot2::margin(t = -2, b = -2, unit = "pt"),
    plot.caption = ggplot2::element_text(size = 4.65, hjust = 0, colour = "#5F6568", lineheight = 0.90, margin = ggplot2::margin(t = 2)),
    plot.margin = ggplot2::margin(5, 5, 5, 5, "pt")
  ) +
  ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 1.8)))

decision_label_source <- decision_source |> dplyr::filter(.data$label_key_module, !is.na(.data$decision_label))
if (requireNamespace("ggrepel", quietly = TRUE)) {
  decision_plot <- decision_plot +
    ggrepel::geom_text_repel(
      data = decision_label_source,
      ggplot2::aes(label = .data$decision_label),
      size = 1.65,
      segment.size = 0.14,
      segment.colour = "#8A8F93",
      box.padding = 0.12,
      point.padding = 0.06,
      force = 2.8,
      force_pull = 0.35,
      max.overlaps = Inf,
      seed = 9
    )
} else {
  decision_plot <- decision_plot +
    ggplot2::geom_text(
      data = decision_label_source,
      ggplot2::aes(label = .data$decision_label),
      size = 1.65,
      hjust = -0.04,
      vjust = -0.25,
      check_overlap = TRUE
    )
}
decision_files <- save_svg_pdf(decision_plot, "microglia_roi_specificity_decision_map", 55, 60)

row_levels <- rev(plot_source$compact_label[order(plot_source$matched_neuropil_coupling, na.last = TRUE)])
fingerprint_continuous <- plot_source |>
  dplyr::mutate(compact_label = factor(.data$compact_label, levels = row_levels)) |>
  tidyr::pivot_longer(
    cols = c("matched_neuropil_coupling", "best_neuropil_protein_sharing"),
    names_to = "diagnostic_axis",
    values_to = "value"
  ) |>
  dplyr::mutate(
    value_display = pmin(pmax(.data$value, 0), 1),
    value_type = "continuous",
    evidence_present = NA,
    diagnostic_axis = factor(
      .data$diagnostic_axis,
      levels = c("matched_neuropil_coupling", "best_neuropil_protein_sharing"),
      labels = c("Neuropil\nscore", "Neuropil\nprotein")
    )
  )

fingerprint_binary <- plot_source |>
  dplyr::mutate(compact_label = factor(.data$compact_label, levels = row_levels)) |>
  tidyr::pivot_longer(
    cols = c("microglia_roi_evidence"),
    names_to = "diagnostic_axis",
    values_to = "evidence_present"
  ) |>
  dplyr::mutate(
    value = NA_real_,
    value_display = NA_real_,
    value_type = "binary",
    diagnostic_axis = factor(
      .data$diagnostic_axis,
      levels = c("microglia_roi_evidence"),
      labels = c("Microglia\nROI ev.")
    )
  )

fingerprint_long <- dplyr::bind_rows(fingerprint_continuous, fingerprint_binary)
write_csv_safe2(fingerprint_long, file.path(PATHS$source_data, "microglia_roi_specificity_fingerprint_source.csv"))

fingerprint_axis_levels <- c("Neuropil\nscore", "Neuropil\nprotein", "Microglia\nROI ev.")
fingerprint_plot <- ggplot2::ggplot(fingerprint_long, ggplot2::aes(x = factor(.data$diagnostic_axis, levels = fingerprint_axis_levels), y = .data$compact_label)) +
  ggplot2::geom_point(
    data = fingerprint_continuous,
    ggplot2::aes(size = .data$value_display, fill = .data$value_display),
    shape = 21,
    colour = "white",
    stroke = 0.15
  ) +
  ggplot2::geom_point(
    data = fingerprint_binary |> dplyr::filter(.data$evidence_present %in% TRUE),
    shape = 3,
    colour = "#2D3436",
    size = 1.8,
    stroke = 0.35
  ) +
  ggplot2::geom_point(
    data = fingerprint_continuous |> dplyr::filter(.data$highlight_target),
    shape = 21,
    fill = NA,
    colour = "#176B57",
    size = 3.2,
    stroke = 0.35
  ) +
  ggplot2::scale_fill_gradient(low = "#E8ECEC", high = "#176B57", na.value = "grey94", limits = c(0, 1), oob = scales::squish, name = "Coupling / sharing") +
  ggplot2::scale_size(range = c(0.7, 3.0), limits = c(0, 1), guide = "none") +
  ggplot2::labs(title = "Microglia ROI module fingerprint", x = NULL, y = NULL, caption = "Score columns show |Spearman rho|; protein columns show fraction of microglia-module proteins found in best matching cross-compartment module. + indicates microglia ROI evidence.") +
  theme_compact +
  ggplot2::theme(
    legend.position = "right",
    legend.key.height = grid::unit(12, "mm"),
    legend.key.width = grid::unit(3, "mm"),
    legend.title = ggplot2::element_text(size = 5.5),
    legend.text = ggplot2::element_text(size = 5),
    axis.text.x = ggplot2::element_text(size = 5.5),
    plot.caption = ggplot2::element_text(size = 5.3, hjust = 0, colour = "#5F6568", lineheight = 0.92, margin = ggplot2::margin(t = 3)),
    plot.margin = ggplot2::margin(6, 8, 8, 6, "pt")
  ) +
  ggplot2::guides(fill = ggplot2::guide_colourbar(barheight = grid::unit(16, "mm"), barwidth = grid::unit(3, "mm")))
fingerprint_files <- save_svg_pdf(fingerprint_plot, "microglia_roi_specificity_fingerprint", 105, 75)

fingerprint_with_soma_files <- character()
if (isTRUE(INCLUDE_SOMA)) {
  fingerprint_with_soma_continuous <- plot_source |>
    dplyr::mutate(compact_label = factor(.data$compact_label, levels = row_levels)) |>
    tidyr::pivot_longer(
      cols = c("matched_neuropil_coupling", "matched_soma_coupling", "best_neuropil_protein_sharing", "best_soma_protein_sharing"),
      names_to = "diagnostic_axis",
      values_to = "value"
    ) |>
    dplyr::mutate(
      value_display = pmin(pmax(.data$value, 0), 1),
      value_type = "continuous",
      evidence_present = NA,
      diagnostic_axis = factor(
        .data$diagnostic_axis,
        levels = c("matched_neuropil_coupling", "matched_soma_coupling", "best_neuropil_protein_sharing", "best_soma_protein_sharing"),
        labels = c("Neuropil\nscore", "Soma\nscore", "Neuropil\nprotein", "Soma\nprotein")
      )
    )
  fingerprint_with_soma_binary <- fingerprint_binary
  fingerprint_with_soma_long <- dplyr::bind_rows(fingerprint_with_soma_continuous, fingerprint_with_soma_binary)
  write_csv_safe2(fingerprint_with_soma_long, file.path(PATHS$source_data, "microglia_roi_specificity_fingerprint_with_soma_source.csv"))
  fingerprint_with_soma_axis_levels <- c("Neuropil\nscore", "Soma\nscore", "Neuropil\nprotein", "Soma\nprotein", "Microglia\nROI ev.")
  fingerprint_with_soma_plot <- ggplot2::ggplot(fingerprint_with_soma_long, ggplot2::aes(x = factor(.data$diagnostic_axis, levels = fingerprint_with_soma_axis_levels), y = .data$compact_label)) +
    ggplot2::geom_point(
      data = fingerprint_with_soma_continuous,
      ggplot2::aes(size = .data$value_display, fill = .data$value_display),
      shape = 21,
      colour = "white",
      stroke = 0.15
    ) +
    ggplot2::geom_point(
      data = fingerprint_with_soma_binary |> dplyr::filter(.data$evidence_present %in% TRUE),
      shape = 3,
      colour = "#2D3436",
      size = 1.8,
      stroke = 0.35
    ) +
    ggplot2::scale_fill_gradient(low = "#E8ECEC", high = "#176B57", na.value = "grey94", limits = c(0, 1), oob = scales::squish, name = "Coupling / sharing") +
    ggplot2::scale_size(range = c(0.7, 3.0), limits = c(0, 1), guide = "none") +
    ggplot2::labs(title = "Microglia ROI module fingerprint with soma context", x = NULL, y = NULL, caption = "Score columns show |Spearman rho|; protein columns show fraction of microglia-module proteins found in best matching cross-compartment module. + indicates microglia ROI evidence.") +
    theme_compact +
    ggplot2::theme(
      legend.position = "right",
      legend.key.height = grid::unit(12, "mm"),
      legend.key.width = grid::unit(3, "mm"),
      legend.title = ggplot2::element_text(size = 5.5),
      legend.text = ggplot2::element_text(size = 5),
      axis.text.x = ggplot2::element_text(size = 5.5),
      plot.caption = ggplot2::element_text(size = 5.3, hjust = 0, colour = "#5F6568", lineheight = 0.92, margin = ggplot2::margin(t = 3)),
      plot.margin = ggplot2::margin(6, 8, 8, 6, "pt")
    ) +
    ggplot2::guides(fill = ggplot2::guide_colourbar(barheight = grid::unit(16, "mm"), barwidth = grid::unit(3, "mm")))
  fingerprint_with_soma_files <- save_svg_pdf(fingerprint_with_soma_plot, "microglia_roi_specificity_fingerprint_with_soma", 105, 75)
}

coupling_plot_source <- coupling_long |>
  dplyr::filter(.data$covariate_family %in% c(family_order, "exploratory_best_spearman")) |>
  dplyr::mutate(
    module_label = compact_module_label(.data$module_id, .data$endpoint_label),
    module_label = factor(.data$module_label, levels = row_levels),
    family_label = factor(dplyr::coalesce(unname(family_display[.data$covariate_family]), .data$covariate_family), levels = unname(family_display[c(family_order, "exploratory_best_spearman")])),
    plot_value = if (METRIC == "delta_r2") .data$delta_adjusted_r2 else .data$abs_spearman
  )
write_csv_safe2(coupling_plot_source, file.path(PATHS$source_data, "microglia_cross_compartment_coupling_heatmap_source.csv"))

coupling_plot_display <- coupling_plot_source |>
  dplyr::filter(isTRUE(INCLUDE_EXPLORATORY) | .data$covariate_family != "exploratory_best_spearman") |>
  dplyr::mutate(
    family_label = factor(
      as.character(.data$family_label),
      levels = unname(family_display[if (isTRUE(INCLUDE_EXPLORATORY)) c(family_order, "exploratory_best_spearman") else family_order])
    )
  )
coupling_facet_labels <- c(neuron_neuropil = "Neuropil", neuron_soma = "Soma")

coupling_plot <- ggplot2::ggplot(coupling_plot_display, ggplot2::aes(x = .data$family_label, y = .data$module_label, fill = .data$plot_value)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
  ggplot2::facet_wrap(~ covariate_source_compartment, nrow = 1, labeller = ggplot2::as_labeller(coupling_facet_labels)) +
  ggplot2::scale_fill_gradient(low = "#F7F8F8", high = "#176B57", na.value = "grey92", limits = c(0, if (METRIC == "delta_r2") max(coupling_plot_display$plot_value, na.rm = TRUE) else 1), oob = scales::squish, name = if (METRIC == "delta_r2") "Delta adj. R2" else "|Spearman r|") +
  ggplot2::labs(title = "Matched score coupling", x = NULL, y = NULL, caption = "Values show |Spearman rho| between microglia module scores and matched cross-compartment scores. No StressGroup contrasts used.") +
  theme_compact +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 8),
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
    legend.position = "right",
    plot.caption = ggplot2::element_text(size = 5.3, hjust = 0, colour = "#5F6568", lineheight = 0.92, margin = ggplot2::margin(t = 3)),
    plot.margin = ggplot2::margin(6, 8, 8, 6, "pt")
  )
coupling_files <- save_svg_pdf(coupling_plot, "microglia_cross_compartment_coupling_heatmap", 110, 75)

card_files <- character()
target_summary <- summary_tbl |> dplyr::filter(.data$module_id == TARGET_MODULE)
if (nrow(target_summary)) {
  target_coupling <- coupling_long |>
    dplyr::filter(.data$module_id == TARGET_MODULE, .data$covariate_family %in% family_order) |>
    dplyr::mutate(family_label = dplyr::coalesce(unname(family_display[.data$covariate_family]), .data$covariate_family))
  card_note <- "Shared local ECM signal in the microglia ROI; not microglia-restricted."
  card_title <- if (identical(TARGET_MODULE, "WGCNA_#4D4D4D")) "#4D4D4D \u00b7 ECM/adhesion" else compact_module_label(TARGET_MODULE, target_summary$endpoint_label[[1]])
  fmt_num <- function(x, digits = 2) ifelse(is.finite(as_num(x)), format(round(as_num(x), digits), nsmall = digits), "NA")
  fmt_pct <- function(x) ifelse(is.finite(as_num(x)), paste0(round(100 * as_num(x)), "%"), "NA")
  share_text <- function(prefix, n_overlap, denom, frac) {
    if (!is.finite(as_num(n_overlap)) || !is.finite(as_num(denom))) return(paste0(prefix, " NA"))
    paste0(prefix, " ", as.integer(n_overlap), "/", as.integer(denom), " (", fmt_pct(frac), ")")
  }
  n_microglia_proteins <- overlap_long |>
    dplyr::filter(.data$module_id == TARGET_MODULE) |>
    dplyr::pull(.data$n_microglia_module_proteins) |>
    unique()
  n_microglia_proteins <- n_microglia_proteins[is.finite(n_microglia_proteins)][1] %||% NA_real_
  card_metrics_all <- tibble::tibble(
    metric = c("Neuropil score coupling", "Soma score coupling", "Neuropil protein sharing", "Soma protein sharing"),
    value = c(
      target_summary$best_neuropil_abs_spearman[[1]],
      target_summary$best_soma_abs_spearman[[1]],
      target_summary$best_matching_neuropil_overlap_fraction_microglia[[1]],
      target_summary$best_matching_soma_overlap_fraction_microglia[[1]]
    ),
    display_value = c(
      fmt_num(target_summary$best_neuropil_abs_spearman[[1]]),
      fmt_num(target_summary$best_soma_abs_spearman[[1]]),
      target_summary$best_neuropil_sharing_label[[1]],
      target_summary$best_soma_sharing_label[[1]]
    ),
    matched_module = c(
      target_summary$best_neuropil_covariate_label[[1]],
      target_summary$best_soma_covariate_family[[1]] %||% NA_character_,
      target_summary$best_matching_neuropil_module[[1]],
      target_summary$best_matching_soma_module[[1]]
    ),
    n_overlap = c(NA_real_, NA_real_, target_summary$best_matching_neuropil_n_overlap[[1]], target_summary$best_matching_soma_n_overlap[[1]]),
    n_microglia_module_proteins = c(NA_real_, NA_real_, n_microglia_proteins, n_microglia_proteins),
    n_other_module_proteins = c(NA_real_, NA_real_, target_summary$best_matching_neuropil_n_other_module_proteins[[1]], target_summary$best_matching_soma_n_other_module_proteins[[1]]),
    jaccard = c(NA_real_, NA_real_, target_summary$best_matching_neuropil_jaccard[[1]], target_summary$best_matching_soma_jaccard[[1]]),
    interpretation_note = target_summary$boss_interpretation[[1]]
  ) |>
    dplyr::mutate(metric = factor(.data$metric, levels = rev(.data$metric)))
  card_metrics <- card_metrics_all |>
    dplyr::filter(as.character(.data$metric) %in% c("Neuropil score coupling", "Neuropil protein sharing")) |>
    dplyr::mutate(metric = factor(as.character(.data$metric), levels = rev(c("Neuropil score coupling", "Neuropil protein sharing"))))
  card_source <- dplyr::bind_rows(
    target_coupling |> dplyr::mutate(source_section = "predeclared_score_coupling"),
    card_metrics_all |> dplyr::mutate(source_section = "card_summary_metric")
  )
  write_csv_safe2(card_source, file.path(PATHS$source_data, "WGCNA_4D4D4D_roi_specificity_card_source.csv"))
  card_sub <- "Local ECM / shared microenvironment"
  card_plot <- ggplot2::ggplot(card_metrics, ggplot2::aes(x = .data$value, y = .data$metric)) +
    ggplot2::geom_col(width = 0.62, fill = "#176B57") +
    ggplot2::geom_text(ggplot2::aes(label = .data$display_value), hjust = -0.08, size = 2.0, colour = "#2D3436") +
    ggplot2::geom_vline(xintercept = thresholds$high_abs_spearman, linetype = "dashed", colour = "#C94A5A", linewidth = 0.25) +
    ggplot2::annotate("label", x = thresholds$high_abs_spearman + 0.055, y = 2.35, label = "high", size = 1.8, colour = "#C94A5A", fill = "white", linewidth = 0, label.padding = grid::unit(0.8, "pt"), hjust = 0) +
    ggplot2::scale_x_continuous(limits = c(0, 1.08), breaks = c(0, 0.5, 1), expand = ggplot2::expansion(mult = c(0, 0.01))) +
    ggplot2::labs(title = card_title, subtitle = card_sub, x = "Score coupling or protein sharing", y = NULL, caption = card_note) +
    theme_compact +
    ggplot2::theme(
      plot.subtitle = ggplot2::element_text(size = 6.2, colour = "#5F6568", margin = ggplot2::margin(b = 3)),
      plot.caption = ggplot2::element_text(size = 5.5, hjust = 0, colour = "#5F6568", margin = ggplot2::margin(t = 3)),
      plot.margin = ggplot2::margin(7, 12, 7, 7, "pt")
    )
  card_files <- save_svg_pdf(card_plot, "WGCNA_4D4D4D_roi_specificity_card", 105, 70)
  if (isTRUE(INCLUDE_SOMA)) {
    card_with_soma_plot <- ggplot2::ggplot(card_metrics_all, ggplot2::aes(x = .data$value, y = .data$metric)) +
      ggplot2::geom_col(width = 0.62, fill = "#176B57") +
      ggplot2::geom_text(ggplot2::aes(label = .data$display_value), hjust = -0.08, size = 2.0, colour = "#2D3436") +
      ggplot2::geom_vline(xintercept = thresholds$high_abs_spearman, linetype = "dashed", colour = "#C94A5A", linewidth = 0.25) +
      ggplot2::annotate("label", x = thresholds$high_abs_spearman + 0.055, y = 4.48, label = "high", size = 1.8, colour = "#C94A5A", fill = "white", linewidth = 0, label.padding = grid::unit(0.8, "pt"), hjust = 0) +
      ggplot2::scale_x_continuous(limits = c(0, 1.08), breaks = c(0, 0.5, 1), expand = ggplot2::expansion(mult = c(0, 0.01))) +
      ggplot2::labs(title = card_title, subtitle = paste0(card_sub, " + soma context"), x = "Score coupling or protein sharing", y = NULL, caption = card_note) +
      theme_compact +
      ggplot2::theme(
        plot.subtitle = ggplot2::element_text(size = 6.2, colour = "#5F6568", margin = ggplot2::margin(b = 3)),
        plot.caption = ggplot2::element_text(size = 5.5, hjust = 0, colour = "#5F6568", margin = ggplot2::margin(t = 3)),
        plot.margin = ggplot2::margin(7, 12, 7, 7, "pt")
      )
    card_files <- c(card_files, save_svg_pdf(card_with_soma_plot, "WGCNA_4D4D4D_roi_specificity_card_with_soma", 105, 70))
  }
}

dir_create(REPORT_DIR)
interpretation_guide <- file.path(REPORT_DIR, "microglia_roi_specificity_interpretation_guide.md")
writeLines(
  c(
    "# How to interpret microglia ROI specificity diagnostics",
    "",
    "## What the analysis asks",
    "",
    "This analysis does not test stress-group effects. It asks whether microglia WGCNA modules look specific to the microglia ROI dataset or are shared/coupled with matched neuron-neuropil and neuron-soma measurements.",
    "",
    "The primary specificity diagnostic compares microglia ROI modules to matched neuron-neuropil signals from the same AnimalID + Region.",
    "",
    "Soma comparisons are optional cross-compartment context and are not required for the main neuropil-specific interpretation.",
    "",
    "## Two independent axes",
    "",
    "- Score coupling = module activity moves with matched cross-compartment scores.",
    "- Protein sharing = the same proteins are present in cross-compartment modules.",
    "",
    "## How to read the specificity map",
    "",
    "- Low coupling + low sharing = best candidate enriched zone.",
    "- High coupling = shared/coupled state.",
    "- High protein sharing = shared/local composition.",
    "- High both = strongly shared/coupled.",
    "",
    "## Main result",
    "",
    "Most microglia WGCNA modules are best interpreted as shared/local ROI programs rather than purified microglia-restricted modules.",
    "",
    "Most microglia WGCNA modules are shared/local ROI programs rather than purified microglia-restricted modules. The ECM/adhesion module is not highly score-coupled to matched neuropil/soma, but its proteins are strongly shared across compartments, supporting a local ECM/shared microenvironment interpretation.",
    "",
    "## #4D4D4D interpretation",
    "",
    "The ECM/adhesion module shows modest score coupling but high protein sharing with neuron-neuropil/soma modules, supporting a local ECM/shared microenvironment interpretation rather than microglia specificity."
  ),
  interpretation_guide,
  useBytes = TRUE
)

write_run_manifest(
  file.path(PATHS$logs, "run_manifest.yml"),
  inputs = inputs,
  outputs = list(
    tables = PATHS$tables,
    source_data = PATHS$source_data,
    figures = c(decision_files, fingerprint_files, fingerprint_with_soma_files, coupling_files, card_files),
    reports = interpretation_guide
  ),
  parameters = list(
    dataset = "microglia",
    target_module = TARGET_MODULE,
    include_soma = INCLUDE_SOMA,
    include_exploratory = INCLUDE_EXPLORATORY,
    metric = METRIC,
    thresholds = thresholds,
    specificity_archetype_thresholds = list(
      high_score_coupling = thresholds$high_abs_spearman,
      high_protein_sharing = thresholds$high_overlap_fraction_microglia,
      low_score_coupling = thresholds$low_abs_spearman,
      low_protein_sharing = thresholds$moderate_overlap_fraction_microglia
    ),
    model = "module_score ~ matched_cross_compartment_score + SpatialLabel + Sex + Batch + optional (1 | AnimalID)",
    stressgroup_contrasts_used = FALSE,
    stressgroup_contrast_models_refit = FALSE,
    wgcna_recomputed = FALSE
  ),
  notes = "Contrast-free diagnostic of microglia ROI module specificity, matched neuropil/soma coupling, and cross-compartment module sharing. Does not test StressGroup effects and does not enable purified cell-type claims without supporting evidence."
)

message("Microglia ROI specificity diagnostics complete.")

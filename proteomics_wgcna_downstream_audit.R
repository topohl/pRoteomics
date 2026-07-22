#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, warn = 1)

report_path <- file.path(getwd(), "proteomics_wgcna_downstream_audit_output.txt")
con <- file(report_path, open = "wt", encoding = "UTF-8")
sink(con, split = TRUE)
on.exit({try(sink(), silent = TRUE); try(close(con), silent = TRUE)}, add = TRUE)

section <- function(x) {
  cat("\n", paste(rep("=", 72), collapse = ""), "\n", x, "\n",
      paste(rep("=", 72), collapse = ""), "\n", sep = "")
}
read_required <- function(path) {
  if (!file.exists(path)) stop("Required file missing: ", path, call. = FALSE)
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}
as_bool <- function(x) {
  y <- suppressWarnings(as.logical(x)); !is.na(y) & y
}
as_num <- function(x) suppressWarnings(as.numeric(x))
min_num <- function(x) {
  x <- as_num(x); x <- x[is.finite(x)]; if (length(x)) min(x) else NA_real_
}
count_le <- function(x, t) {
  x <- as_num(x); sum(is.finite(x) & x <= t)
}
print_df <- function(x) print(as.data.frame(x), row.names = FALSE)
ndistinct <- function(x) length(unique(as.character(x)[!is.na(x) & nzchar(as.character(x))]))
get_col <- function(df, candidates, default = NA) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit)) df[[hit[[1]]]] else rep(default, nrow(df))
}

summarise_stage05 <- function(dataset, level) {
  file <- if (level == "module") "module_group_effects.csv" else "supermodule_group_effects.csv"
  x <- read_required(file.path("results", "tables", "06_modules_WGCNA", "group_effects", dataset, file))
  id <- get_col(x, if (level == "module") c("module_id", "ModuleID", "endpoint_id") else c("supermodule_id", "SupermoduleID", "endpoint_id"))
  p <- as_num(get_col(x, c("p_value", "raw_p", "p")))
  fw <- as_num(get_col(x, c("FDR_within_dataset_level", "FDR", "fdr")))
  fg <- as_num(get_col(x, c("FDR_global")))
  est <- as_num(get_col(x, c("estimate", "effect_size")))
  claim <- as_bool(get_col(x, c("claim_allowed_model"), FALSE))
  stable <- as_bool(get_col(x, c("primary_model_stable"), FALSE))
  fallback <- as_bool(get_col(x, c("fallback_used"), FALSE))
  model_text <- tolower(paste(
    get_col(x, c("model_type", "model_family"), ""),
    get_col(x, c("fallback_type"), ""),
    get_col(x, c("evidence_status"), ""),
    get_col(x, c("model_warning"), "")
  ))
  fallback <- fallback | grepl("fallback|diagnostic|model_unstable", model_text)
  key_cols <- intersect(c(if (level == "module") "module_id" else "supermodule_id", "contrast", "effect_scope", "spatial_unit"), names(x))
  key <- if (length(key_cols)) do.call(paste, c(lapply(x[key_cols], as.character), sep = "||")) else rep(NA_character_, nrow(x))
  data.frame(
    dataset = dataset, level = level, n_rows = nrow(x), n_unique_ids = ndistinct(id),
    duplicate_endpoint_keys = sum(duplicated(key)), n_claim_allowed_model = sum(claim),
    n_primary_model_stable = sum(stable), n_fallback_or_diagnostic = sum(fallback),
    min_raw_p = min_num(p), min_FDR_within = min_num(fw), min_FDR_global = min_num(fg),
    n_FDR_within_10 = count_le(fw, 0.10), n_FDR_global_10 = count_le(fg, 0.10),
    invalid_p = sum(is.finite(p) & (p < 0 | p > 1)),
    invalid_FDR_within = sum(is.finite(fw) & (fw < 0 | fw > 1)),
    invalid_FDR_global = sum(is.finite(fg) & (fg < 0 | fg > 1)),
    nonfinite_estimates = sum(!is.na(est) & !is.finite(est))
  )
}

main <- function() {
  cat("Proteomics WGCNA downstream audit\n")
  cat("Working directory: ", normalizePath(getwd(), winslash = "/", mustWork = FALSE), "\n", sep = "")
  cat("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")

  section("A. ATLAS SEMANTICS")
  atlas <- read_required(file.path("results", "tables", "10_biological_integration", "cross_compartment_program_atlas", "global", "cross_compartment_program_atlas_long.csv"))
  ds <- sort(unique(atlas$dataset))
  atlas_counts <- do.call(rbind, lapply(ds, function(d) {
    z <- atlas[atlas$dataset == d, , drop = FALSE]
    data.frame(dataset = d, n_total = nrow(z), n_counting = sum(as_bool(z$counts_toward_convergence)), n_noncounting = sum(!as_bool(z$counts_toward_convergence)))
  }))
  print_df(atlas_counts)

  unavailable_mask <- grepl("Unavailable optional evidence|missing optional input|input unavailable", paste(atlas$program_label, atlas$interpretation_note, atlas$evidence_status), ignore.case = TRUE)
  unavailable <- atlas[unavailable_mask, c("dataset", "evidence_domain", "evidence_role", "evidence_status", "counts_toward_convergence", "program_label"), drop = FALSE]
  cat("\nUnavailable placeholders:\n"); print_df(unavailable)

  neuronal_wgcna <- atlas[atlas$dataset %in% c("neuron_soma", "neuron_neuropil") & atlas$evidence_domain == "wgcna_supermodule", , drop = FALSE]
  cat("\nNeuronal WGCNA effects:\n")
  print_df(do.call(rbind, lapply(sort(unique(neuronal_wgcna$dataset)), function(d) {
    z <- neuronal_wgcna[neuronal_wgcna$dataset == d, , drop = FALSE]
    data.frame(dataset = d, n_rows = nrow(z), n_counting = sum(as_bool(z$counts_toward_convergence)), evidence_roles = paste(sort(unique(z$evidence_role)), collapse = ";"), n_FDR10 = count_le(z$fdr, 0.10), min_fdr = min_num(z$fdr))
  })))

  section("B. MICROGLIA STAGE 13")
  readiness <- read_required(file.path("results", "tables", "06_modules_WGCNA", "claim_readiness", "microglia", "WGCNA_entity_claim_readiness.csv"))
  print_df(as.data.frame(table(readiness$claim_entity_role), stringsAsFactors = FALSE))
  eligible <- readiness[as_bool(readiness$separate_manuscript_claim_allowed), c("level", "entity_id", "claim_entity_role", "primary_architecture_status", "group_effect_status"), drop = FALSE]
  cat("\nEligible identities:\n"); print_df(eligible[order(eligible$level, eligible$entity_id), ])

  section("C. CLAIM GATING")
  claims <- read_required(file.path("results", "tables", "biological_claims_table.csv"))
  arch <- claims[claims$dataset %in% c("neuron_soma", "neuron_neuropil") & claims$claim_type == "wgcna_architecture", , drop = FALSE]
  print_df(do.call(rbind, lapply(sort(unique(arch$dataset)), function(d) {
    z <- arch[arch$dataset == d, , drop = FALSE]
    data.frame(dataset = d, n_rows = nrow(z), n_allowed = sum(as_bool(z$claim_allowed)), n_disallowed = sum(z$claim_gate_status == "disallowed", na.rm = TRUE), n_annotation_only = sum(z$claim_use_class == "annotation_only", na.rm = TRUE))
  })))
  effects <- claims[claims$dataset %in% c("neuron_soma", "neuron_neuropil") & claims$claim_type == "wgcna_group_effect", , drop = FALSE]
  cat("\nNeuronal WGCNA group effects:\n")
  print_df(do.call(rbind, lapply(sort(unique(effects$dataset)), function(d) {
    z <- effects[effects$dataset == d, , drop = FALSE]
    data.frame(dataset = d, n_rows = nrow(z), n_allowed = sum(as_bool(z$claim_allowed)), n_FDR10 = count_le(z$FDR, 0.10), min_FDR = min_num(z$FDR))
  })))

  section("D. STAGE 05 MATHEMATICS")
  stage05 <- do.call(rbind, list(
    summarise_stage05("neuron_soma", "module"),
    summarise_stage05("neuron_soma", "supermodule"),
    summarise_stage05("neuron_neuropil", "module"),
    summarise_stage05("neuron_neuropil", "supermodule")
  ))
  print_df(stage05)

  section("E. CIRCULAR ATLAS")
  segments <- read_required(file.path("results", "source_data", "10_biological_integration", "wgcna_circular_atlas", "global", "wgcna_circular_atlas_segments.csv"))
  print_df(do.call(rbind, lapply(sort(unique(segments$dataset)), function(d) {
    z <- segments[segments$dataset == d, , drop = FALSE]
    data.frame(dataset = d, n_segments = nrow(z), n_descriptive_selected = sum(as_bool(z$selected_for_descriptive_atlas)), n_manuscript_selected = sum(as_bool(z$selected_for_manuscript_claim)))
  })))
  cat("\nManuscript-selected identities:\n")
  print_df(segments[as_bool(segments$selected_for_manuscript_claim), c("dataset", "supermodule_id", "claim_display_status", "claim_entity_role"), drop = FALSE])

  duplicates <- read_required(file.path("results", "reports", "10_biological_integration", "wgcna_circular_atlas", "global", "wgcna_circular_atlas_duplicate_source_audit.csv"))
  cat("\nUnsafe duplicate source rows:\n")
  print_df(duplicates[!as_bool(duplicates$collapse_safe), c("dataset", "supermodule_id", "source_cleaned_label", "n_member_modules_source", "cleaned_label_differs", "n_member_modules_differs", "effect_estimate_differs", "FDR_global_differs"), drop = FALSE])

  callout <- read_required(file.path("results", "source_data", "10_biological_integration", "wgcna_circular_atlas", "global", "wgcna_supermodule_callout_source.csv"))
  cat("\nCallout metric sources:\n")
  print_df(as.data.frame(table(callout$dataset, callout$effect_source), stringsAsFactors = FALSE))

  section("F. FINAL BUNDLE")
  bundle <- file.path("results", "tables", "10_biological_integration", "final_evidence_bundle", "global")
  key_modules <- read_required(file.path(bundle, "wgcna_key_modules.csv"))
  key_supermodules <- read_required(file.path(bundle, "wgcna_key_supermodules.csv"))
  cat("\nKey modules by dataset/status:\n")
  print_df(as.data.frame(table(key_modules$dataset, key_modules$status, useNA = "ifany"), stringsAsFactors = FALSE))
  cat("\nKey supermodules by dataset/status:\n")
  print_df(as.data.frame(table(key_supermodules$dataset, key_supermodules$status, useNA = "ifany"), stringsAsFactors = FALSE))

  section("G. PASS/FAIL SUMMARY")
  count_for <- function(d) atlas_counts$n_counting[match(d, atlas_counts$dataset)]
  selected_keys <- paste(segments$dataset[as_bool(segments$selected_for_manuscript_claim)], segments$supermodule_id[as_bool(segments$selected_for_manuscript_claim)], sep = "::")
  checks <- data.frame(
    check = c(
      "atlas_microglia_counting_275", "atlas_neuropil_counting_660", "atlas_soma_counting_264",
      "all_unavailable_noncounting", "neuronal_wgcna_918_rows", "neuronal_wgcna_zero_counting",
      "stage13_22_identities", "stage13_13_modules", "stage13_3_blocks", "stage13_6_aliases",
      "stage13_10_eligible_modules", "stage13_SM09_only_eligible_block", "stage13_zero_FDR_supported_group_effects",
      "soma_8_architecture_rows", "neuropil_14_architecture_rows", "neuronal_architecture_all_disallowed",
      "neuronal_architecture_all_annotation_only", "neuronal_group_effects_zero_allowed",
      "circular_27_segments", "circular_SM09_only_manuscript_selected",
      "two_neuronal_module_status_rows", "two_neuronal_supermodule_status_rows",
      "microglia_10_unique_key_modules", "microglia_SM09_only_key_supermodule",
      "stage05_zero_claim_allowed_models", "stage05_zero_FDR10", "stage05_zero_duplicate_endpoint_keys", "stage05_probability_ranges_valid"
    ),
    pass = c(
      as.integer(count_for("microglia")) == 275L,
      as.integer(count_for("neuron_neuropil")) == 660L,
      as.integer(count_for("neuron_soma")) == 264L,
      nrow(unavailable) > 0 && all(!as_bool(unavailable$counts_toward_convergence)),
      nrow(neuronal_wgcna) == 918,
      !any(as_bool(neuronal_wgcna$counts_toward_convergence)),
      nrow(readiness) == 22,
      sum(readiness$claim_entity_role == "canonical_module") == 13,
      sum(readiness$claim_entity_role == "higher_order_block") == 3,
      sum(readiness$claim_entity_role == "compatibility_alias") == 6,
      sum(readiness$claim_entity_role == "canonical_module" & as_bool(readiness$separate_manuscript_claim_allowed)) == 10,
      setequal(readiness$entity_id[readiness$claim_entity_role == "higher_order_block" & as_bool(readiness$separate_manuscript_claim_allowed)], "SM09"),
      !any(readiness$group_effect_status == "FDR_supported"),
      sum(arch$dataset == "neuron_soma") == 8,
      sum(arch$dataset == "neuron_neuropil") == 14,
      !any(as_bool(arch$claim_allowed)),
      all(arch$claim_use_class == "annotation_only"),
      !any(as_bool(effects$claim_allowed)),
      nrow(segments) == 27,
      setequal(selected_keys, "microglia::SM09"),
      sum(key_modules$dataset %in% c("neuron_soma", "neuron_neuropil") & key_modules$status == "no_validated_neuronal_wgcna_key_rows") == 2,
      sum(key_supermodules$dataset %in% c("neuron_soma", "neuron_neuropil") & key_supermodules$status == "no_validated_neuronal_wgcna_key_rows") == 2,
      ndistinct(key_modules$module_id[key_modules$dataset == "microglia"]) == 10,
      setequal(unique(key_supermodules$supermodule_id[key_supermodules$dataset == "microglia"]), "SM09"),
      all(stage05$n_claim_allowed_model == 0),
      all(stage05$n_FDR_within_10 == 0 & stage05$n_FDR_global_10 == 0),
      all(stage05$duplicate_endpoint_keys == 0),
      all(stage05$invalid_p == 0 & stage05$invalid_FDR_within == 0 & stage05$invalid_FDR_global == 0 & stage05$nonfinite_estimates == 0)
    ),
    stringsAsFactors = FALSE
  )
  checks$status <- ifelse(checks$pass, "PASS", "FAIL")
  print_df(checks)
  cat("\nOverall: ", if (all(checks$pass)) "PASS" else "FAIL", "\n", sep = "")
}

status <- tryCatch({main(); 0L}, error = function(e) {
  section("FATAL AUDIT ERROR")
  cat(conditionMessage(e), "\n")
  1L
})
cat("\nReport written to: ", normalizePath(report_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
try(sink(), silent = TRUE)
try(close(con), silent = TRUE)
quit(status = status, save = "no")

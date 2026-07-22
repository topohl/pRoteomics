# Shared helpers for manuscript-level biological integration scripts.

if (!exists("repo_path", mode = "function")) {
  paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
  source(paths_file)
}
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "module_contracts.R"))

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}

integration_cli <- function(default_dataset = "all", allow_all = TRUE) {
  args <- commandArgs(trailingOnly = TRUE)
  value_after <- function(flag, default = "") {
    hit <- which(args == flag)
    if (!length(hit) || hit[[1]] == length(args)) return(default)
    args[[hit[[1]] + 1L]]
  }
  has_dataset_arg <- "--dataset" %in% args
  env_dataset <- Sys.getenv("PROTEOMICS_DATASET", unset = "")
  dataset_default <- if (isTRUE(allow_all) && identical(default_dataset, "all") && !has_dataset_arg) {
    "all"
  } else if (nzchar(env_dataset)) {
    env_dataset
  } else {
    default_dataset
  }
  dataset <- value_after("--dataset", dataset_default)
  dataset <- normalize_dataset(dataset)
  if (isTRUE(allow_all) && identical(dataset, "all")) {
    Sys.setenv(PROTEOMICS_DATASET = dataset)
  } else {
    dataset <- validate_dataset(dataset, source = "--dataset")
    Sys.setenv(PROTEOMICS_DATASET = dataset)
  }
  list(args = args, dataset = dataset, dry_run = is_dry_run())
}

integration_datasets <- function(dataset_arg) {
  if (identical(dataset_arg, "all")) valid_datasets() else validate_dataset(dataset_arg)
}

integration_paths <- function(substep, dataset = "global") {
  create_module_dirs("10_biological_integration", file.path(substep, dataset))
}

read_csv_optional <- function(path, dataset = "global", evidence_domain = "input",
                              input_type = basename(path), required = FALSE) {
  exists <- file.exists(path)
  record_input_resolution(
    script = Sys.getenv("PROTEOMICS_SCRIPT_ID", unset = NA_character_),
    dataset = dataset,
    stage = "integration",
    input_name = input_type,
    expected_path = path,
    resolved_path = path,
    resolution_mode = if (exists) "canonical" else if (required) "missing_required" else "missing_optional",
    strict_mode = strict_inputs_enabled(),
    allowed_in_strict_mode = TRUE,
    producer_script_or_artifact_id = evidence_domain,
    warning = if (!exists && isTRUE(required)) "Required integration input missing." else NA_character_
  )
  status <- data.frame(
    dataset = dataset,
    evidence_domain = evidence_domain,
    input_type = input_type,
    path = normalizePath(path, winslash = "/", mustWork = FALSE),
    required = isTRUE(required),
    status = if (exists) "present" else if (required) "missing_required" else "missing_optional",
    message = if (exists) "loaded" else "input not available; downstream rows marked unavailable",
    n_rows = 0L,
    stringsAsFactors = FALSE
  )
  if (!exists) return(list(data = NULL, status = status))
  data <- tryCatch(
    {
      if (requireNamespace("readr", quietly = TRUE)) {
        readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
      } else {
        utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
      }
    },
    error = function(e) {
      status$status <- "read_error"
      status$message <- conditionMessage(e)
      NULL
    }
  )
  if (!is.null(data)) status$n_rows <- nrow(data)
  list(data = data, status = status)
}

write_csv_safe <- function(x, path) {
  dir_create(dirname(path))
  if (requireNamespace("readr", quietly = TRUE)) {
    readr::write_csv(x, path, na = "")
  } else {
    utils::write.csv(x, path, row.names = FALSE, na = "")
  }
  invisible(path)
}

write_integration_table <- function(x, paths, filename) {
  list(
    table = write_csv_safe(x, file.path(paths$tables, filename)),
    source = write_csv_safe(x, file.path(paths$source_data, filename))
  )
}

first_col <- function(df, candidates) {
  if (is.null(df)) return(NA_character_)
  hit <- candidates[candidates %in% names(df)]
  if (length(hit)) hit[[1]] else NA_character_
}

col_or_na <- function(df, col, default = NA) {
  if (!is.null(df) && col %in% names(df)) df[[col]] else rep(default, if (is.null(df)) 0L else nrow(df))
}

num_or_na <- function(x) suppressWarnings(as.numeric(x))

chr_clean <- function(x, fallback = NA_character_) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(trimws(x))] <- fallback
  x
}

normalize_token <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("_MOUSE$", "", x, ignore.case = TRUE)
  gsub("[^A-Z0-9]", "", x)
}

split_gene_tokens <- function(x) {
  x <- as.character(x %||% "")
  toks <- unlist(strsplit(x, "[/;,|[:space:]]+"), use.names = FALSE)
  unique(normalize_token(toks[nzchar(toks)]))
}

empty_status <- function() {
  data.frame(
    dataset = character(), evidence_domain = character(), input_type = character(),
    path = character(), required = logical(), status = character(),
    message = character(), n_rows = integer(), stringsAsFactors = FALSE
  )
}

empty_evidence <- function() {
  data.frame(
    dataset = character(), evidence_domain = character(), evidence_id = character(),
    program_label = character(), entity_type = character(), entity_id = character(),
    wgcna_level = character(), canonical_claim_entity_id = character(),
    claim_entity_role = character(), separate_manuscript_claim_allowed = logical(),
    wgcna_architecture_status = character(), wgcna_group_effect_status = character(),
    wgcna_allowed_claim_scope = character(), wgcna_prohibited_claim_scope = character(),
    readiness_contract_version = character(), counts_toward_convergence = logical(),
    evidence_semantic_class = character(), evidence_role = character(), evidence_role_reason = character(),
    contrast = character(), spatial_unit = character(), effect_direction = character(),
    effect_size = numeric(), p_value = numeric(), fdr = numeric(),
    support_count = numeric(), source_file = character(), evidence_status = character(),
    interpretation_note = character(), qc_flag = character(), stringsAsFactors = FALSE
  )
}

standardize_evidence <- function(df) {
  cols <- names(empty_evidence())
  if (is.null(df) || !nrow(df)) return(empty_evidence())
  for (col in cols) if (!col %in% names(df)) df[[col]] <- NA
  df <- df[, cols, drop = FALSE]
  missing_semantic <- is.na(df$evidence_semantic_class) | !nzchar(trimws(as.character(df$evidence_semantic_class)))
  df$evidence_semantic_class[missing_semantic] <- ifelse(
    grepl("wgcna|module", as.character(df$evidence_domain[missing_semantic]), ignore.case = TRUE),
    "wgcna_stress_group_effect",
    "non_wgcna_evidence"
  )
  missing_count <- is.na(df$counts_toward_convergence)
  df$counts_toward_convergence[missing_count] <- !(
    as.character(df$evidence_status[missing_count]) %in% c("unavailable_optional_input", "missing_required") |
      as.character(df$entity_type[missing_count]) == "status"
  )
  for (col in c("effect_size", "p_value", "fdr", "support_count")) df[[col]] <- num_or_na(df[[col]])
  for (col in c("separate_manuscript_claim_allowed", "counts_toward_convergence")) df[[col]] <- as.logical(df[[col]])
  for (col in setdiff(cols, c("effect_size", "p_value", "fdr", "support_count", "separate_manuscript_claim_allowed", "counts_toward_convergence"))) df[[col]] <- as.character(df[[col]])
  df
}

# Apply after all atlas evidence has been combined.  This keeps source rows intact
# while separating independent evidence from annotation, diagnostic, and unavailable
# context.  Microglia Stage 13 countability is intentionally preserved as supplied.
assign_downstream_evidence_roles <- function(df) {
  df <- standardize_evidence(df)
  domain <- tolower(trimws(as.character(df$evidence_domain)))
  status <- tolower(trimws(as.character(df$evidence_status)))
  source <- tolower(as.character(df$source_file))
  program <- tolower(trimws(as.character(df$program_label)))
  note <- tolower(as.character(df$interpretation_note))
  status[is.na(status)] <- ""
  source[is.na(source)] <- ""
  program[is.na(program)] <- ""
  note[is.na(note)] <- ""
  neuronal <- as.character(df$dataset) %in% c("neuron_soma", "neuron_neuropil")
  unavailable <- as.character(df$entity_type) == "status" |
    grepl("unavailable|missing", status) |
    grepl("unavailable optional evidence|missing optional input", program) |
    grepl("input unavailable", note)
  annotation <- domain %in% c("microenvironment_marker", "complex_organelle", "wgcna_compatibility_provenance")
  diagnostic <- domain %in% c("module_robustness_sensitivity", "robustness_sensitivity", "qc_confounding") |
    domain == "module_behavior_coupling" |
    (domain == "behavior_coupling" & grepl("08_behavior_physio_coupling/(03_)?module_behavior_coupling", source))
  neuronal_wgcna <- neuronal & domain %in% c("wgcna_supermodule", "wgcna_stress_group_effect")
  microglia_architecture <- as.character(df$dataset) == "microglia" & domain == "wgcna_architecture"

  df$evidence_role <- "direct_biological_evidence"
  df$evidence_role_reason <- NA_character_
  df$evidence_role[annotation] <- "descriptive_annotation"
  df$evidence_role_reason[annotation] <- "annotation_or_compatibility_provenance_not_independent_convergence"
  df$evidence_role[diagnostic] <- "diagnostic_gate"
  df$evidence_role_reason[diagnostic] <- "diagnostic_or_context_input_not_independent_convergence"
  df$evidence_role[microglia_architecture] <- "validated_wgcna_architecture"
  df$evidence_role_reason[microglia_architecture] <- "microglia_stage13_claim_readiness_contract"
  df$evidence_role[neuronal_wgcna] <- "unvalidated_neuronal_wgcna"
  df$evidence_role_reason[neuronal_wgcna] <- "neuronal_wgcna_readiness_contract_not_available"
  df$evidence_role[unavailable] <- "unavailable_evidence"
  df$evidence_role_reason[unavailable] <- "missing_or_unavailable_optional_input"

  force_noncounting <- annotation | diagnostic | neuronal_wgcna | unavailable
  force_noncounting[is.na(force_noncounting)] <- FALSE
  df$counts_toward_convergence[force_noncounting] <- FALSE
  df
}

availability_evidence <- function(dataset, evidence_domain, source_file, message) {
  data.frame(
    dataset = dataset,
    evidence_domain = evidence_domain,
    evidence_id = paste(dataset, evidence_domain, "unavailable", sep = "::"),
    program_label = "Unavailable optional evidence",
    entity_type = "status",
    entity_id = NA_character_,
    contrast = NA_character_,
    spatial_unit = NA_character_,
    effect_direction = NA_character_,
    effect_size = NA_real_,
    p_value = NA_real_,
    fdr = NA_real_,
    support_count = 0,
    source_file = source_file,
    evidence_status = "unavailable_optional_input",
    interpretation_note = message,
    qc_flag = "WARN",
    stringsAsFactors = FALSE
  )
}

evidence_strength <- function(fdr, support_count = NA_real_) {
  fdr <- num_or_na(fdr)
  support_count <- num_or_na(support_count)
  if (!is.na(fdr) && fdr <= 0.05) return("strong")
  if (!is.na(fdr) && fdr <= 0.10) return("moderate")
  if (!is.na(support_count) && support_count >= 2) return("supportive")
  "exploratory"
}

program_key <- function(x) {
  z <- tolower(as.character(x))
  dplyr_case <- function(...) {
    if (requireNamespace("dplyr", quietly = TRUE)) return(dplyr::case_when(...))
    NULL
  }
  out <- dplyr_case(
    grepl("mitochond|respirat|oxidative|\\batp\\b|tca", z) ~ "Mitochondrial metabolism",
    grepl("\\brna\\b|ribosom|translation|splice|rnp", z) ~ "RNA/translation",
    grepl("synap|vesicle|postsynap|cytoskeleton|actin|microtub", z) ~ "Synaptic/cytoskeletal",
    grepl("ecm|adhesion|collagen|laminin|integrin|basement", z) ~ "Perivascular ECM/adhesion",
    grepl("microglia|immune|phago|lysosom|complement|inflamm", z) ~ "Microglia/immune state",
    grepl("vascular|bbb|endothelial|pericyte|blood vessel", z) ~ "Vascular/BBB",
    grepl("myelin|oligodendro", z) ~ "Myelin/oligodendrocyte",
    grepl("astrocyte|endfoot", z) ~ "Astrocyte/endfoot",
    TRUE ~ "Mixed/unresolved"
  )
  if (!is.null(out)) return(out)
  out <- rep("Mixed/unresolved", length(z))
  out[grepl("mitochond|respirat|oxidative|\\batp\\b|tca", z)] <- "Mitochondrial metabolism"
  out[grepl("\\brna\\b|ribosom|translation|splice|rnp", z)] <- "RNA/translation"
  out[grepl("synap|vesicle|postsynap|cytoskeleton|actin|microtub", z)] <- "Synaptic/cytoskeletal"
  out
}

write_integration_manifest <- function(paths, inputs, outputs, parameters, notes) {
  write_run_manifest(
    file.path(paths$logs, "run_manifest.yml"),
    inputs = inputs,
    outputs = outputs,
    parameters = parameters,
    notes = notes
  )
}

dry_run_inputs <- function(label, inputs) {
  dry_run_line("Script", label)
  old_script <- Sys.getenv("PROTEOMICS_SCRIPT_ID", unset = NA_character_)
  Sys.setenv(PROTEOMICS_SCRIPT_ID = label)
  on.exit({
    if (is.na(old_script)) Sys.unsetenv("PROTEOMICS_SCRIPT_ID") else Sys.setenv(PROTEOMICS_SCRIPT_ID = old_script)
  }, add = TRUE)
  nms <- names(inputs)
  for (i in seq_along(inputs)) {
    nm <- nms[[i]] %||% paste0("input_", i)
    path <- inputs[[i]]
    record_input_resolution(
      script = label,
      dataset = Sys.getenv("PROTEOMICS_DATASET", unset = "global"),
      stage = "integration",
      input_name = nm,
      expected_path = path,
      resolved_path = path,
      resolution_mode = if (file.exists(path) || dir.exists(path)) "canonical" else "missing_optional",
      strict_mode = strict_inputs_enabled(),
      allowed_in_strict_mode = TRUE,
      producer_script_or_artifact_id = "dry_run_inputs"
    )
    dry_run_line(nm, path, if (file.exists(path) || dir.exists(path)) "PASS" else "WARN")
  }
}

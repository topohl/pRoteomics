# Shared helpers for WGCNA downstream interpretation scripts.

if (!exists("repo_path", mode = "function")) {
  paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
  source(paths_file)
}
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "dataset_inputs.R"))
source(repo_path("R", "module_contracts.R"))

WGCNA_ROI_NOTE <- "microglia-enriched ROI/local microenvironment; annotation only, not purity correction."

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}

shorten_supermodule_label <- function(x, max_chars = 45) {
  vapply(as.character(x), function(z) {
    z <- trimws(z)
    if (is.na(z) || !nzchar(z)) return("Unresolved / mixed")
    explicit_multi_theme <- grepl("^mixed:\\s*", z, ignore.case = TRUE)
    z <- gsub("\\s+([0-9]+)\\s+[Mm]odules?$", "", z)
    z <- gsub("\\b[Mm]odules?\\b", "", z)
    z <- gsub("\\s*\\([^)]*modules?[^)]*\\)", "", z, ignore.case = TRUE)
    z <- gsub("\\b(process|regulation|pathway|of|the|cellular|biological|positive|negative)\\b", "", z, ignore.case = TRUE)
    z <- gsub("\\s+", " ", z)
    parts <- trimws(unlist(strsplit(z, "\\s*[;|]\\s*", perl = TRUE), use.names = FALSE))
    parts <- parts[nzchar(parts)]
    if (length(parts) && !explicit_multi_theme) z <- paste(utils::head(parts, 2), collapse = " / ")
    if (length(parts) && explicit_multi_theme) z <- paste(parts, collapse = "; ")
    z <- stringr::str_squish(z)
    if (!nzchar(z)) z <- "Unresolved / mixed"
    if (nchar(z) > max_chars && !explicit_multi_theme) {
      words <- unlist(strsplit(z, "\\s+"), use.names = FALSE)
      keep <- character()
      for (word in words) {
        cand <- paste(c(keep, word), collapse = " ")
        if (nchar(cand) > max_chars) break
        keep <- c(keep, word)
      }
      z <- if (length(keep)) paste(keep, collapse = " ") else substr(z, 1, max_chars)
    }
    z
  }, character(1))
}

wgcna_member_theme_priority <- function() {
  c(
    "mitochondrial / energy metabolism",
    "translation / proteostasis",
    "RNA/RNP regulatory module",
    "synaptic/cytoskeletal trafficking",
    "ECM/adhesion",
    "barrier / cell-junction structural",
    "neuronal signalling / regulatory",
    "myelin / oligodendrocyte-associated",
    "mixed / low-specificity"
  )
}

wgcna_order_member_themes <- function(themes, fractions = NULL) {
  themes <- as.character(themes)
  if (!length(themes)) return(themes)
  priority <- match(themes, wgcna_member_theme_priority())
  priority[is.na(priority)] <- length(wgcna_member_theme_priority()) + seq_len(sum(is.na(priority)))
  if (is.null(fractions)) fractions <- rep(0, length(themes))
  fractions <- suppressWarnings(as.numeric(fractions))
  fractions[is.na(fractions)] <- 0
  themes[order(-fractions, priority, themes)]
}

wgcna_member_theme_summary <- function(member_theme, display_threshold = 0.20, dominant_threshold = 0.60) {
  themes <- trimws(as.character(member_theme))
  themes <- themes[!is.na(themes) & nzchar(themes)]
  if (!length(themes)) {
    return(list(
      counts = setNames(integer(), character()),
      fractions = setNames(numeric(), character()),
      counts_label = "",
      fractions_label = "",
      n_distinct = 0L,
      is_multi = FALSE,
      themes_above_display_threshold = "",
      display_label = "mixed / low-specificity",
      themes_omitted_from_display_label = "",
      displayed_themes = character()
    ))
  }
  raw_counts <- table(themes)
  raw_fracs <- as.numeric(raw_counts) / length(themes)
  ord <- wgcna_order_member_themes(names(raw_counts), raw_fracs)
  counts <- raw_counts[ord]
  fractions <- setNames(as.numeric(counts) / length(themes), names(counts))
  counts_label <- paste(paste0(names(counts), "=", as.integer(counts)), collapse = "; ")
  fractions_label <- paste(paste0(names(fractions), "=", sprintf("%.2f", fractions)), collapse = "; ")
  informative <- names(fractions)[names(fractions) != "mixed / low-specificity"]
  above <- informative[fractions[informative] >= display_threshold]
  above <- wgcna_order_member_themes(above, fractions[above])
  dominant <- informative[fractions[informative] >= dominant_threshold]
  dominant <- wgcna_order_member_themes(dominant, fractions[dominant])

  if (length(dominant)) {
    displayed <- dominant[[1]]
    display_label <- paste0("mostly ", displayed)
  } else if (length(above) >= 2L && length(above) <= 3L) {
    displayed <- above
    display_label <- paste0("mixed: ", paste(displayed, collapse = "; "))
  } else if (length(above) > 3L) {
    displayed <- character()
    display_label <- "mixed multi-program"
  } else {
    displayed <- character()
    display_label <- "mixed / low-specificity"
  }

  omitted <- setdiff(above, displayed)
  list(
    counts = counts,
    fractions = fractions,
    counts_label = counts_label,
    fractions_label = fractions_label,
    n_distinct = length(counts),
    is_multi = length(counts) > 1L,
    themes_above_display_threshold = paste(above, collapse = "; "),
    display_label = display_label,
    themes_omitted_from_display_label = paste(omitted, collapse = "; "),
    displayed_themes = displayed
  )
}

macroprogram_display <- function(x) {
  vapply(as.character(x), function(z) {
    z0 <- tolower(trimws(z %||% ""))
    if (!nzchar(z0)) return("Unresolved / mixed")
    if (grepl("ecm|extracellular matrix|basement membrane|collagen|laminin|nidogen|\\bagrn\\b|hspg2|perlecan|\\bcol4a[12]?\\b|\\blama[0-9]?\\b|\\blamb[0-9]?\\b|\\blamc[0-9]?\\b|\\bnid[12]\\b|\\bbcam\\b|serpinh1", z0)) return("Perivascular ECM / adhesion")
    if (grepl("mitochondr|respiratory|oxidative|\\batp\\b|\\btca\\b|acetyl-coa", z0)) return("Mitochondrial metabolism")
    if (grepl("\\brna\\b|ribosome|translation|splice|\\brnp\\b|ncrna", z0)) return("RNA / translation")
    if (grepl("synapse|vesicle|postsynaptic|actin|cytoskeleton", z0)) return("Synaptic / cytoskeletal")
    if (grepl("microglia|phagolysosomal|immune", z0)) return("Microglia state")
    if (grepl("neuropil|neuronal", z0)) return("Neuropil / neuronal")
    if (grepl("\\bbbb\\b|endothelial|pericyte|vascular", z0)) return("Vascular / BBB")
    "Unresolved / mixed"
  }, character(1))
}

supermodule_microenvironment_label <- function(cls, dataset = current_dataset()) {
  vapply(as.character(cls), function(z) {
    z <- tolower(trimws(z %||% ""))
    if (!nzchar(z)) return(NA_character_)
    if (z %in% c("vascular_basement_membrane_ecm", "vascular/ecm")) return("Perivascular ECM")
    if (z %in% c("vascular_bbb_mural", "vascular")) return("Vascular / BBB")
    if (z %in% c("neuropil_sensitive", "neuropil")) return("Neuropil reference overlap")
    if (z %in% c("astrocyte_or_endfoot_sensitive", "astrocyte")) return("Astrocyte / endfoot")
    if (z %in% c("oligodendrocyte_or_myelin_sensitive", "oligodendrocyte/myelin")) return("Oligodendrocyte / myelin")
    if (z %in% c("ambiguous_or_mixed", "shared_microenvironment", "mixed")) return("Mixed / unresolved")
    if (identical(as.character(dataset), "microglia") && z %in% c("microglia_supported", "microglia_state_or_activation_supported")) return("Microglia-associated ROI")
    NA_character_
  }, character(1))
}

compose_supermodule_display_label <- function(supermodule_id, short_label) {
  id <- as.character(supermodule_id)
  id[is.na(id) | !nzchar(id)] <- "SM??"
  label <- shorten_supermodule_label(short_label, max_chars = 30)
  label[is.na(label) | !nzchar(label)] <- "Mixed / unresolved"
  paste0(id, " \u00b7 ", label)
}

classify_supermodule_label_confidence <- function(n_modules, go_class = "unresolved",
                                                  has_coherent_hubs = FALSE,
                                                  microenvironment_class = NA_character_,
                                                  high_unmapped_fraction = FALSE) {
  n_modules <- suppressWarnings(as.integer(n_modules))
  singleton <- is.na(n_modules) | n_modules <= 1L
  go_class <- as.character(go_class %||% "unresolved")
  micro_label <- supermodule_microenvironment_label(microenvironment_class)
  micro_supported <- !is.na(micro_label) & nzchar(micro_label) & !micro_label %in% "Mixed / unresolved"
  mixed_micro <- !is.na(micro_label) & micro_label == "Mixed / unresolved"
  go_supported <- go_class %in% c("GO_supported", "data_driven_GO")
  suggestive_go <- go_class %in% c("suggestive_GO", "manual_only")
  if (singleton || isTRUE(high_unmapped_fraction) || mixed_micro) return("low")
  if ((go_supported || micro_supported) && isTRUE(has_coherent_hubs)) return("high")
  if ((go_supported && micro_supported) || (suggestive_go && isTRUE(has_coherent_hubs)) || (micro_supported && isTRUE(has_coherent_hubs))) return("medium")
  if (isTRUE(has_coherent_hubs) || suggestive_go || go_supported || micro_supported) return("low")
  "unresolved"
}

wgcna_cli <- function(default_dataset = "neuron_neuropil", allow_all = FALSE) {
  args <- commandArgs(trailingOnly = TRUE)
  value_after <- function(flag, default = "") {
    hit <- which(args == flag)
    if (!length(hit) || hit[[1]] == length(args)) return(default)
    args[[hit[[1]] + 1L]]
  }
  dataset <- value_after("--dataset", Sys.getenv("PROTEOMICS_DATASET", unset = default_dataset))
  if (allow_all && identical(tolower(dataset), "all")) {
    dataset <- "all"
  } else {
    dataset <- validate_dataset(dataset, source = "--dataset")
    Sys.setenv(PROTEOMICS_DATASET = dataset)
  }
  list(
    args = args,
    dataset = dataset,
    level = value_after("--level", "both"),
    dry_run = is_dry_run()
  )
}

wgcna_downstream_paths <- function(substep, dataset) {
  create_module_dirs("06_modules_WGCNA", file.path(substep, dataset))
}

safe_read_csv <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  if (requireNamespace("readr", quietly = TRUE)) {
    return(tryCatch(readr::read_csv(path, show_col_types = FALSE, progress = FALSE), error = function(e) NULL))
  }
  tryCatch(utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
}

write_csv_safe2 <- function(x, path) {
  dir_create(dirname(path))
  if (requireNamespace("readr", quietly = TRUE)) readr::write_csv(x, path, na = "") else utils::write.csv(x, path, row.names = FALSE, na = "")
  invisible(path)
}

write_table_and_source <- function(x, table_dir, source_dir, filename) {
  out <- list(
    table = write_csv_safe2(x, file.path(table_dir, filename)),
    source = write_csv_safe2(x, file.path(source_dir, filename))
  )
  invisible(out)
}

first_present_col <- function(df, candidates) {
  if (is.null(df)) return(NA_character_)
  norm <- function(x) tolower(gsub("[^a-z0-9]", "", x))
  hit <- match(norm(candidates), norm(names(df)))
  hit <- hit[!is.na(hit)]
  if (!length(hit)) return(NA_character_)
  names(df)[hit[[1]]]
}

normalize_wgcna_animal_id <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "NA", "N/A", "NULL", "NONE", "NAN")] <- NA_character_
  out <- rep(NA_character_, length(x))
  has_a_id <- !is.na(x) & grepl("(^|[^A-Z0-9])A0*[0-9]+([^A-Z0-9]|$)", x, perl = TRUE)
  if (any(has_a_id)) {
    token <- regmatches(x[has_a_id], regexpr("A0*[0-9]+", x[has_a_id], perl = TRUE))
    digits <- sub("^A0*", "", token)
    digits[digits == ""] <- "0"
    out[has_a_id] <- paste0("A", digits)
  }
  remaining <- is.na(out) & !is.na(x)
  if (any(remaining)) {
    numeric_like <- grepl("^[0-9]+$", x[remaining])
    vals <- x[remaining]
    vals[numeric_like] <- paste0("A", sub("^0+", "", vals[numeric_like]))
    vals[vals == "A"] <- "A0"
    vals[!numeric_like] <- gsub("[^A-Z0-9]+", "", vals[!numeric_like])
    vals[!nzchar(vals)] <- NA_character_
    out[remaining] <- vals
  }
  out
}

infer_wgcna_animal_id <- function(meta, sample_values) {
  aliases <- c(
    "AnimalID", "AnimalNum", "AnimalNumber", "AnimalNo", "Animal",
    "MouseID", "Mouse", "MouseNum", "MouseNumber", "mouse_id",
    "animal_id", "animal_num", "animal_number", "subject", "SubjectID",
    "donor", "DonorID", "ID", "animalid", "animalnum"
  )
  animal_col <- first_present_col(meta, aliases)
  from_col <- if (!is.na(animal_col)) normalize_wgcna_animal_id(meta[[animal_col]]) else rep(NA_character_, nrow(meta))
  missing <- is.na(from_col) | !nzchar(from_col)
  if (any(missing)) {
    parsed <- normalize_wgcna_animal_id(sample_values)
    from_col[missing] <- parsed[missing]
  }
  from_col
}

normalize_gene_token <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("_MOUSE$", "", x, ignore.case = TRUE)
  gsub("[^A-Z0-9]", "", x)
}

split_tokens <- function(x) {
  x <- as.character(x)
  out <- unlist(strsplit(x, "[/;,|[:space:]]+"), use.names = FALSE)
  out <- toupper(trimws(out))
  unique(out[nzchar(out) & !is.na(out)])
}

wgcna_split_label_terms <- function(x) {
  vals <- unlist(strsplit(as.character(x), "\\s*[;|]\\s*"), use.names = FALSE)
  vals <- gsub("_", " ", vals)
  vals <- gsub("\\bGO\\b", "", vals, ignore.case = TRUE)
  vals <- gsub("\\s+", " ", vals)
  vals <- trimws(vals)
  unique(vals[nzchar(vals) & !is.na(vals)])
}

wgcna_collapse_terms <- function(x, n = 5L) {
  x <- unique(trimws(as.character(x)))
  x <- x[nzchar(x) & !is.na(x)]
  paste(utils::head(x, n), collapse = ";")
}

wgcna_marker_panel_summary <- function(marker_panel_fractions = NULL) {
  if (is.null(marker_panel_fractions)) return(list(label = "", values = numeric()))
  if (is.data.frame(marker_panel_fractions)) {
    nms <- names(marker_panel_fractions)
    vals <- suppressWarnings(as.numeric(marker_panel_fractions[1, , drop = TRUE]))
    names(vals) <- nms
  } else {
    vals <- suppressWarnings(as.numeric(marker_panel_fractions))
    names(vals) <- names(marker_panel_fractions)
  }
  vals <- vals[is.finite(vals) & !is.na(names(vals)) & nzchar(names(vals))]
  vals <- vals[order(vals, decreasing = TRUE)]
  vals <- vals[vals > 0]
  label <- if (length(vals)) paste0(names(utils::head(vals, 8L)), "=", sprintf("%.2f", utils::head(vals, 8L)), collapse = ";") else ""
  list(label = label, values = vals)
}

wgcna_assign_semantic_program <- function(raw_label = NA_character_, go_terms = character(),
                                          hubs = character(), marker_label = NA_character_,
                                          bp_terms = NULL, mf_terms = NULL, cc_terms = NULL,
                                          marker_panel_fractions = NULL) {
  bp <- if (is.null(bp_terms)) character() else wgcna_split_label_terms(bp_terms)
  mf <- if (is.null(mf_terms)) character() else wgcna_split_label_terms(mf_terms)
  cc <- if (is.null(cc_terms)) character() else wgcna_split_label_terms(cc_terms)
  all_go <- wgcna_split_label_terms(c(go_terms, bp, mf, cc))
  raw_terms <- wgcna_split_label_terms(raw_label)
  hub_tokens <- split_tokens(hubs)
  marker_summary <- wgcna_marker_panel_summary(marker_panel_fractions)
  marker_values <- marker_summary$values

  bp_text <- tolower(paste(bp, collapse = " "))
  mf_text <- tolower(paste(mf, collapse = " "))
  cc_text <- tolower(paste(cc, collapse = " "))
  evidence_text <- tolower(paste(raw_terms, all_go, hub_tokens, marker_label, marker_summary$label, collapse = " "))
  hub_text <- tolower(paste(hub_tokens, collapse = " "))
  marker_text <- tolower(paste(marker_label, marker_summary$label, collapse = " "))

  has <- function(pattern, text = evidence_text) grepl(pattern, text, perl = TRUE)
  panel_max <- function(pattern) {
    hit <- grepl(pattern, names(marker_values), ignore.case = TRUE)
    if (any(hit)) max(marker_values[hit], na.rm = TRUE) else NA_real_
  }
  ecm_panel <- panel_max("ecm|basement|matrix|collagen|laminin|perivascular")
  vascular_panel <- panel_max("vascular|endothelial|pericyte|bbb|mural")
  high_ecm_panel <- is.finite(ecm_panel) && ecm_panel >= 0.20
  vascular_context <- is.finite(vascular_panel) && vascular_panel >= 0.20

  explicit_ecm <- has("ecm|extracellular matrix|basement membrane|collagen|laminin|nidogen|\\bagrn\\b|hspg2|perlecan|\\bcol4a[12]?\\b|\\blama[0-9]?\\b|\\blamb[0-9]?\\b|\\blamc[0-9]?\\b|\\bnid[12]\\b|\\bbcam\\b|serpinh1")
  adhesion_evidence <- has("adhesion|cell-cell adhesion mediator activity|cell adhesion molecule")
  integrin_evidence <- has("integrin|\\bitga[0-9a-z]*\\b|\\bitgb[0-9a-z]*\\b")
  synaptic_cc <- has("presynaptic membrane|postsynaptic density|synaptic membrane|synapse|postsynap|presynap", cc_text)
  synaptic_evidence <- synaptic_cc || has("synap|vesicle|postsynap|presynap|neurotrans|axon|dendrit|scaffold|pdz|ankyrin|membrane raft|traffick|actin|cytoskeleton|microtubule")
  ion_regulatory <- has("monoatomic ion transport|ion transport|cation transport|channel activity|regulation of.*ion transport", bp_text)
  pdz_binding <- has("pdz domain binding", mf_text)

  warnings <- character()
  adhesion_interpretation <- NA_character_
  if (adhesion_evidence && synaptic_cc) {
    adhesion_interpretation <- "synaptic adhesion / membrane scaffold"
  } else if ((adhesion_evidence || integrin_evidence) && (explicit_ecm || high_ecm_panel || vascular_context)) {
    adhesion_interpretation <- "ECM-supported adhesion"
  } else if (adhesion_evidence || integrin_evidence) {
    adhesion_interpretation <- "ambiguous adhesion/integrin evidence without ECM context"
    warnings <- c(warnings, "adhesion_or_integrin_without_explicit_ecm_context")
  }

  scores <- c(
    `mitochondrial / energy metabolism` = as.integer(has("generation precursor metabolites|precursor metabolites|mitochond|respirat|oxidative|electron transport|\\batp\\b|\\btca\\b|acetyl[- ]coa|nad[hp]?|cytochrome")),
    `translation / proteostasis` = as.integer(has("ribosom|translation|protein folding|proteasom|chaperon|proteostasis|ubiquitin|er stress")),
    `RNA/RNP regulatory module` = as.integer(has("\\brna\\b|mrna|splice|\\brnp\\b|hnrnp|nucleolar|ribonucl|ncrna|ribonucleoprotein")),
    `synaptic/cytoskeletal trafficking` = as.integer(synaptic_evidence),
    `ECM/adhesion` = as.integer(explicit_ecm || high_ecm_panel || ((adhesion_evidence || integrin_evidence) && vascular_context)),
    `barrier / cell-junction structural` = as.integer(has("skin development|keratinocyte|epiderm|cornified envelope|desmosome|adherens junction|cell-cell junction|tight junction|barrier")),
    `neuronal signalling / regulatory` = as.integer(has("signal|kinase|phosph|redox|neuronal|neuron projection|stimulus|ion homeostasis")),
    `myelin / oligodendrocyte-associated` = as.integer(has("ensheathment neurons|myelin|oligodendro|\\bmbp\\b|\\bmag\\b|\\bmog\\b|\\bplp1\\b|\\bcnp\\b")),
    `microglia/immune-associated ROI` = as.integer(has("lysosom|phago|complement|immune|inflamm|antigen|interferon|microglia")),
    `vascular / BBB` = as.integer(has("vascular|blood vessel|endothelial|pericyte|\\bbbb\\b|blood-brain barrier|mural"))
  )
  if (explicit_ecm || high_ecm_panel) scores[["ECM/adhesion"]] <- scores[["ECM/adhesion"]] + 1L
  if (synaptic_cc) scores[["synaptic/cytoskeletal trafficking"]] <- scores[["synaptic/cytoskeletal trafficking"]] + 1L
  if (ion_regulatory && synaptic_cc) scores[["synaptic/cytoskeletal trafficking"]] <- scores[["synaptic/cytoskeletal trafficking"]] + 1L

  ontology_mismatch <- FALSE
  label <- NA_character_
  primary <- NA_character_
  secondary <- NA_character_
  rule <- "evidence_scoring_with_specificity_gates"
  confidence <- "low"

  if (synaptic_cc && (pdz_binding || adhesion_evidence) && ion_regulatory) {
    label <- "synaptic membrane / ion-transport regulatory scaffold"
    primary <- "synaptic/cytoskeletal trafficking"
    secondary <- "neuronal signalling / regulatory"
    confidence <- "high"
    rule <- "synaptic_CC_plus_PDZ_or_adhesion_MF_plus_ion_transport_BP"
  } else if (has("skin development|keratinocyte|epiderm|cornified envelope|desmosome|adherens junction|cell-cell junction")) {
    label <- "barrier / cell-junction structural module"
    primary <- "barrier / cell-junction structural"
    confidence <- "low"
    rule <- "brain_context_cleanup_for_epidermal_cell_junction_GO"
    ontology_mismatch <- TRUE
  } else if (has("binding sperm|fertilization|egg coat|zona pellucida", tolower(paste(raw_terms, collapse = " ")))) {
    ontology_mismatch <- TRUE
    if (synaptic_evidence) {
      label <- "synaptic / vesicle-neurosecretory module"
      primary <- "synaptic/cytoskeletal trafficking"
      confidence <- "low"
      rule <- "fertilization_GO_context_mismatch_cleaned_by_synaptic_evidence"
    } else if (explicit_ecm || high_ecm_panel) {
      label <- "ECM/adhesion"
      primary <- "ECM/adhesion"
      confidence <- "low"
      rule <- "fertilization_GO_context_mismatch_cleaned_by_ECM_evidence"
    } else {
      label <- "low-confidence ontology artefact"
      primary <- "mixed / low-specificity"
      confidence <- "low"
      rule <- "fertilization_GO_context_mismatch_without_specific_support"
    }
  } else if (has("dorsal ventral pattern", tolower(paste(raw_terms, collapse = " ")))) {
    ontology_mismatch <- TRUE
    if (scores[["RNA/RNP regulatory module"]] > 0L) {
      label <- "RNA/RNP regulatory module"
      primary <- "RNA/RNP regulatory module"
    } else if (scores[["neuronal signalling / regulatory"]] > 0L || synaptic_evidence) {
      label <- "neuronal signalling / regulatory module"
      primary <- "neuronal signalling / regulatory"
    } else {
      label <- "mixed regulatory module"
      primary <- "mixed / low-specificity"
    }
    confidence <- "low"
    rule <- "developmental_patterning_GO_context_mismatch_cleaned_by_supporting_evidence"
  }

  if (is.na(primary) || !nzchar(primary)) {
    ranked <- sort(scores, decreasing = TRUE)
    positive <- ranked[ranked > 0]
    if (length(positive)) {
      primary <- names(positive)[[1]]
      secondary <- if (length(positive) >= 2L) names(positive)[[2]] else NA_character_
      label <- primary
      confidence <- if (positive[[1]] >= 2L || (primary != "ECM/adhesion" && length(raw_terms) + length(all_go) > 0L)) "medium" else "low"
    } else if (length(raw_terms) || length(all_go)) {
      label <- shorten_supermodule_label(wgcna_collapse_terms(c(raw_terms, all_go), n = 4L), max_chars = 42)
      primary <- "mixed / low-specificity"
      confidence <- "low"
      rule <- "no_specific_semantic_gate_matched_raw_label_retained_for_audit"
      warnings <- c(warnings, "no_specific_semantic_program_classified")
    } else {
      label <- "mixed / low-specificity"
      primary <- "mixed / low-specificity"
      rule <- "no_semantic_evidence"
      warnings <- c(warnings, "no_semantic_evidence")
    }
  }
  if (is.na(secondary) || !nzchar(secondary)) {
    ranked <- sort(scores, decreasing = TRUE)
    positive <- names(ranked[ranked > 0])
    secondary_candidates <- positive[positive != primary]
    secondary <- if (length(secondary_candidates)) secondary_candidates[[1]] else NA_character_
  }

  flag <- dplyr::case_when(
    ontology_mismatch ~ "ontology_context_mismatch",
    rule %in% c("brain_context_cleanup_for_epidermal_cell_junction_GO") ~ "cleaned_go_display",
    TRUE ~ "direct_or_suggestive"
  )
  rationale <- paste0(
    "semantic classifier used BP/MF/CC/hub/marker-panel evidence with adhesion/ECM specificity gates; rule=", rule,
    if (nzchar(marker_summary$label)) paste0("; marker_panels=", marker_summary$label) else ""
  )

  list(
    module_specific_label = label,
    module_program_primary = primary,
    module_program_secondary = secondary,
    label_confidence = confidence,
    label_decision_rule = rule,
    label_warning = paste(unique(warnings), collapse = ";"),
    adhesion_interpretation = adhesion_interpretation %||% "",
    evidence_BP = wgcna_collapse_terms(bp, n = 20L),
    evidence_MF = wgcna_collapse_terms(mf, n = 20L),
    evidence_CC = wgcna_collapse_terms(cc, n = 20L),
    evidence_hubs = paste(utils::head(unique(hub_tokens), 30L), collapse = ";"),
    evidence_marker_panels = marker_summary$label,
    GO_label_relevance_flag = flag,
    GO_label_relevance_rationale = rationale,
    ManualReviewFromSemanticLabel = ontology_mismatch || length(warnings) > 0L
  )
}

normalize_wgcna_module_key <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^WGCNA[_-]?", "", x, ignore.case = TRUE)
  x <- sub("^ME", "", x, ignore.case = FALSE)
  x <- sub("^#", "", x)
  tolower(gsub("[^A-Za-z0-9]", "", x))
}

wgcna_go_terms_for_module <- function(go, module_color, n = 5L) {
  empty <- list(
    raw_GO_BP_terms = NA_character_,
    raw_GO_MF_terms = NA_character_,
    raw_GO_CC_terms = NA_character_,
    raw_top_GO_label = NA_character_,
    top_GO_BP_labels = NA_character_,
    top_GO_MF_labels = NA_character_,
    top_GO_CC_labels = NA_character_
  )
  if (is.null(go) || !nrow(go)) return(empty)
  color_col <- first_present_col(go, c("ModuleColor", "module", "ModuleID", "module_color", "module_eigengene"))
  desc_col <- first_present_col(go, c("Description", "description", "term_description", "ID"))
  ont_col <- first_present_col(go, c("Ontology", "ONTOLOGY", "ontology"))
  padj_col <- first_present_col(go, c("p.adjust", "p_adj", "padj", "qvalue", "FDR"))
  gene_col <- first_present_col(go, c("geneID", "gene_id", "genes", "core_enrichment", "overlap_proteins"))
  if (is.na(color_col) || is.na(desc_col)) return(empty)

  target <- normalize_wgcna_module_key(module_color)
  tab <- go[normalize_wgcna_module_key(go[[color_col]]) %in% target, , drop = FALSE]
  if (!nrow(tab)) return(empty)
  tab[[desc_col]] <- trimws(as.character(tab[[desc_col]]))
  tab <- tab[nzchar(tab[[desc_col]]) & !is.na(tab[[desc_col]]), , drop = FALSE]
  if (!nrow(tab)) return(empty)
  if (is.na(ont_col)) {
    tab$.__ontology <- "BP"
    ont_col <- ".__ontology"
  }
  tab$.__ontology_norm <- toupper(trimws(as.character(tab[[ont_col]])))
  tab$.__ontology_norm[!tab$.__ontology_norm %in% c("BP", "MF", "CC")] <- "BP"
  tab$.__gene_key <- if (!is.na(gene_col)) as.character(tab[[gene_col]]) else ""
  tab <- tab[!duplicated(paste(tab$.__ontology_norm, tab[[desc_col]], tab$.__gene_key, sep = "\r")), , drop = FALSE]
  if (!is.na(padj_col)) {
    tab$.__padj <- suppressWarnings(as.numeric(tab[[padj_col]]))
    tab <- tab[order(is.na(tab$.__padj), tab$.__padj), , drop = FALSE]
  }
  one <- function(ont) wgcna_collapse_terms(tab[[desc_col]][tab$.__ontology_norm == ont], n = n)
  bp <- one("BP")
  mf <- one("MF")
  cc <- one("CC")
  top <- wgcna_collapse_terms(c(wgcna_split_label_terms(bp), wgcna_split_label_terms(mf), wgcna_split_label_terms(cc)), n = 1L)
  list(
    raw_GO_BP_terms = if (nzchar(bp)) bp else NA_character_,
    raw_GO_MF_terms = if (nzchar(mf)) mf else NA_character_,
    raw_GO_CC_terms = if (nzchar(cc)) cc else NA_character_,
    raw_top_GO_label = if (nzchar(top)) top else NA_character_,
    top_GO_BP_labels = if (nzchar(bp)) bp else NA_character_,
    top_GO_MF_labels = if (nzchar(mf)) mf else NA_character_,
    top_GO_CC_labels = if (nzchar(cc)) cc else NA_character_
  )
}

wgcna_clean_semantic_label <- function(raw_label = NA_character_, go_terms = character(),
                                       hubs = character(), marker_label = NA_character_,
                                       bp_terms = NULL, mf_terms = NULL, cc_terms = NULL,
                                       marker_panel_fractions = NULL) {
  semantic <- wgcna_assign_semantic_program(
    raw_label = raw_label,
    go_terms = go_terms,
    hubs = hubs,
    marker_label = marker_label,
    bp_terms = bp_terms,
    mf_terms = mf_terms,
    cc_terms = cc_terms,
    marker_panel_fractions = marker_panel_fractions
  )
  raw <- wgcna_collapse_terms(c(raw_label, go_terms, bp_terms, mf_terms, cc_terms), n = 8L)
  label <- semantic$module_specific_label
  flag <- semantic$GO_label_relevance_flag
  rationale <- semantic$GO_label_relevance_rationale

  c(
    list(
      cleaned_biological_label = label,
      cleaned_biological_label_short = shorten_supermodule_label(label, max_chars = 34),
      cleaned_biological_label_source = dplyr::case_when(
        grepl("ontology_context_mismatch|cleaned_go_display", flag) ~ "cleaned_GO_or_hub_evidence",
        nzchar(raw) ~ "GO_or_module_label",
        TRUE ~ "unresolved"
      ),
      cleaned_biological_label_confidence = semantic$label_confidence,
      cleaned_biological_label_rationale = rationale,
      GO_label_relevance_flag = flag,
      GO_label_relevance_rationale = rationale,
      ManualReviewFromSemanticLabel = semantic$ManualReviewFromSemanticLabel
    ),
    semantic[c(
      "module_specific_label", "module_program_primary", "module_program_secondary",
      "label_confidence", "label_decision_rule", "label_warning", "adhesion_interpretation",
      "evidence_BP", "evidence_MF", "evidence_CC", "evidence_hubs", "evidence_marker_panels"
    )]
  )
}

wgcna_microenvironment_caution <- function(cls, label = NA_character_, rationale = NA_character_, dataset = current_dataset()) {
  cls <- as.character(cls %||% NA_character_)
  lab <- as.character(label %||% NA_character_)
  ds <- as.character(dataset %||% current_dataset())
  caution <- if (identical(ds, "microglia")) {
    dplyr::case_when(
      cls == "shared_microenvironment" ~ "shared microglia-neuropil ROI",
      cls == "neuropil_sensitive" ~ "neuropil-sensitive microglia ROI",
      cls == "vascular_basement_membrane_ecm" ~ "perivascular/ECM microglia ROI",
      cls == "vascular_bbb_mural" ~ "vascular/BBB microglia ROI",
      cls %in% c("microglia_supported", "microglia_state_or_activation_supported") ~ "microglia-supported ROI",
      cls == "ambiguous_or_mixed" ~ "curated/shared ROI caution; not microglia-specific by itself",
      TRUE ~ "curated/shared ROI caution; not microglia-specific by itself"
    )
  } else if (identical(ds, "neuron_soma")) {
    dplyr::case_when(
      cls %in% c("vascular_basement_membrane_ecm", "vascular_bbb_mural", "astrocyte_or_endfoot_sensitive", "oligodendrocyte_or_myelin_sensitive", "ambiguous_or_mixed") ~ "shared local tissue / marker-overlap caution",
      TRUE ~ "neuronal soma-enriched tissue program"
    )
  } else if (identical(ds, "neuron_neuropil")) {
    dplyr::case_when(
      cls %in% c("vascular_basement_membrane_ecm", "vascular_bbb_mural", "astrocyte_or_endfoot_sensitive", "oligodendrocyte_or_myelin_sensitive", "ambiguous_or_mixed") ~ "shared local tissue / marker-overlap caution",
      TRUE ~ "neuron-neuropil tissue program"
    )
  } else {
    "not_applicable"
  }
  list(
    microenvironment_caution_label = caution,
    microenvironment_caution_class = cls,
    microenvironment_caution_rationale = paste(c(na.omit(c(rationale, WGCNA_ROI_NOTE))), collapse = "; ")
  )
}

compose_supermodule_composition_display <- function(composition_label, caution_label = NA_character_, dataset = current_dataset()) {
  comp <- as.character(composition_label)
  caution <- as.character(caution_label)
  prefix <- dplyr::case_when(
    identical(as.character(dataset), "microglia") & grepl("^shared", caution, ignore.case = TRUE) ~ "shared ROI: ",
    identical(as.character(dataset), "microglia") & grepl("^neuropil", caution, ignore.case = TRUE) ~ "neuropil-sensitive ROI: ",
    identical(as.character(dataset), "microglia") & grepl("perivascular|vascular|BBB", caution, ignore.case = TRUE) ~ "vascular/ECM ROI: ",
    identical(as.character(dataset), "microglia") & grepl("microglia-supported", caution, ignore.case = TRUE) ~ "microglia-supported ROI: ",
    TRUE ~ ""
  )
  out <- paste0(prefix, comp)
  out[is.na(comp) | !nzchar(comp)] <- "mixed / low-specificity"
  out
}

wgcna_marker_sets <- function() {
  list(
    microglia = c("Aif1", "Tmem119", "P2ry12", "Cx3cr1", "Csf1r", "Hexb",
  "Fcrls", "Olfml3", "Sall1", "Siglech", "Gpr34", "P2ry13",
  "Tgfbr1", "Mertk", "Spi1", "C1qa", "C1qb", "C1qc", "Tyrobp", "Trem2", "Apoe",
  "Lpl", "Ctsb", "Ctsd", "Ctsz", "Lgals3", "Itgam",
  "Cd68", "Fcgr3", "Clec7a"),
    neuronal_synaptic_neuropil = c("Stxbp1", "Gpm6a", "Nptn", "Sh3gl2", "Atp6v1g2",
  "Snap25", "Snap91", "Syn1", "Syn2", "Syp", "Syt1",
  "Vamp2", "Vamp1", "Dlg4", "Dlg3", "Shank1", "Shank2",
  "Homer1", "Camk2a", "Camk2b", "Gria1", "Gria2",
  "Grin1", "Grin2a", "Map2", "Tubb3", "Nefl", "Nefm",
  "Nefh", "Rbfox3"),
    neuropil_synaptic_neuronal = c("Stxbp1", "Gpm6a", "Nptn", "Sh3gl2", "Atp6v1g2",
  "Snap25", "Snap91", "Syn1", "Syn2", "Syp", "Syt1",
  "Vamp2", "Vamp1", "Dlg4", "Dlg3", "Shank1", "Shank2",
  "Homer1", "Camk2a", "Camk2b", "Gria1", "Gria2",
  "Grin1", "Grin2a", "Map2", "Tubb3", "Nefl", "Nefm",
  "Nefh", "Rbfox3"),
    nuclear_soma = c("H2ac1", "H4c1", "H3-3a", "H1-4", "H1-3",
  "Matr3", "Srsf3", "Ddx39b", "Lmna", "Lmnb1", "Lmnb2",
  "Ncl", "Npm1", "Nono", "Sfpq", "Hnrnpm", "Hnrnpu",
  "Rbm39", "Top2b", "Hist1h1c"),
    astrocyte = c("Gfap", "Aqp4", "Aldh1l1", "Slc1a2", "Slc1a3",
  "Aldoc", "Glul", "Gja1", "S100b", "Fabp7",
  "Clu", "Sparcl1", "Fgfr3", "Pla2g7", "Mlc1"),
    oligodendrocyte_myelin = c("Mbp", "Mog", "Plp1", "Cnp", "Mag", "Mobp", "Cldn11",
  "Myrf", "Olig1", "Olig2", "Sox10", "Mal", "Fa2h",
  "Ugt8a", "Opalin", "Ermn", "Tspan2"),
    endothelial_pericyte_vascular = c("Pecam1", "Cldn5", "Kdr", "Flt1", "Tek", "Klf2",
  "Klf4", "Slco1a4", "Abcb1a", "Plvap", "Vwf", "Emcn"),
    mitochondrial_oxphos = c("Ndufs1", "Ndufs2", "Ndufv1", "Ndufa9", "Ndufb8",
  "Sdha", "Sdhb", "Uqcrc1", "Uqcrc2", "Cyc1",
  "Cox4i1", "Cox5a", "Cox6c", "Atp5f1a", "Atp5f1b",
  "Atp5f1c", "Atp5mc1", "Atp5mc2", "Atp5pb"),
    ribosomal_translation = c("Rpl3", "Rpl4", "Rpl5", "Rpl7", "Rpl10", "Rpl13",
  "Rps3", "Rps6", "Rps8", "Rps14", "Rps18",
  "Eef1a1", "Eef1a2", "Eef2", "Eif3a", "Eif4a1",
  "Eif4g1", "Eif5a"),
    rnp_rna_processing = c("Hnrnpa2b1", "Hnrnpc", "Hnrnpm", "Hnrnpu",
  "Sfpq", "Nono", "Snrnp70", "Snrpa", "Snrpb",
  "Ddx5", "Ddx17", "Ddx39b", "Pabpc1", "Pabpc4",
  "Rbm39", "Srsf1", "Srsf3", "Srsf7")
  )
}

wgcna_registry_required_columns <- c(
  "marker_set", "cell_class", "cell_state", "gene_symbol", "source_type",
  "source_name", "source_reference", "selection_rule", "confidence", "use_for", "notes"
)

read_wgcna_marker_registry <- function(path = Sys.getenv("PROTEOMICS_WGCNA_MARKER_REGISTRY_FILE", unset = "")) {
  if (!nzchar(path)) path <- repo_path("config", "marker_panels", "wgcna_reference_marker_sets.csv")
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  registry <- safe_read_csv(path)
  if (is.null(registry) || !nrow(registry)) return(NULL)
  missing <- setdiff(wgcna_registry_required_columns, names(registry))
  if (length(missing)) stop("Marker registry is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  registry$gene_symbol <- as.character(registry$gene_symbol)
  registry$gene_token <- normalize_gene_token(registry$gene_symbol)
  registry <- registry[nzchar(registry$gene_token) & !is.na(registry$gene_token), , drop = FALSE]
  attr(registry, "marker_registry_file") <- path
  attr(registry, "marker_registry_version") <- paste(unique(registry$source_name), collapse = ";")
  registry
}

marker_registry_to_sets <- function(registry) {
  if (is.null(registry) || !nrow(registry)) return(list())
  split(as.character(registry$gene_symbol), as.character(registry$marker_set)) |>
    lapply(function(x) unique(x[nzchar(normalize_gene_token(x))]))
}

read_empirical_roi_marker_sets <- function(path = Sys.getenv("PROTEOMICS_WGCNA_EMPIRICAL_MARKER_FILE", unset = "")) {
  if (!nzchar(path)) path <- path_results("tables", "03_qc_exploration", "05_empirical_roi_marker_discovery", "empirical_roi_marker_sets.csv")
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  empirical <- safe_read_csv(path)
  if (is.null(empirical) || !nrow(empirical)) return(NULL)
  required <- c("marker_set", "GeneSymbol")
  missing <- setdiff(required, names(empirical))
  if (length(missing)) stop("Empirical marker file is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  empirical$gene_symbol <- as.character(empirical$GeneSymbol)
  empirical$gene_token <- normalize_gene_token(empirical$gene_symbol)
  empirical <- empirical[nzchar(empirical$gene_token) & !is.na(empirical$gene_token), , drop = FALSE]
  attr(empirical, "empirical_marker_file") <- path
  attr(empirical, "empirical_marker_set_version") <- paste(unique(empirical$marker_source %||% "empirical_roi_marker_sets"), collapse = ";")
  empirical
}

load_wgcna_marker_sets <- function(include_empirical = TRUE, include_legacy_aliases = TRUE, quiet = FALSE) {
  registry <- read_wgcna_marker_registry()
  empirical <- if (isTRUE(include_empirical)) read_empirical_roi_marker_sets() else NULL
  sets <- list()
  metadata <- data.frame(marker_set = character(), marker_source = character(), source_file = character(), stringsAsFactors = FALSE)

  if (!is.null(registry)) {
    ref_sets <- marker_registry_to_sets(registry)
    sets <- c(sets, ref_sets)
    metadata <- rbind(metadata, data.frame(
      marker_set = names(ref_sets),
      marker_source = "reference_registry",
      source_file = attr(registry, "marker_registry_file") %||% NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  if (!is.null(empirical)) {
    emp_sets <- split(as.character(empirical$GeneSymbol), as.character(empirical$marker_set)) |>
      lapply(function(x) unique(x[nzchar(normalize_gene_token(x))]))
    sets <- c(sets, emp_sets)
    metadata <- rbind(metadata, data.frame(
      marker_set = names(emp_sets),
      marker_source = "empirical_roi_marker_sets",
      source_file = attr(empirical, "empirical_marker_file") %||% NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  if (!length(sets)) {
    if (!isTRUE(quiet)) warning("Falling back to legacy hard-coded WGCNA marker panels; run 03_qc_exploration/04b_import_reference_marker_sources.r to create the registry.", call. = FALSE)
    sets <- wgcna_marker_sets()
    metadata <- data.frame(marker_set = names(sets), marker_source = "legacy_hardcoded_fallback", source_file = NA_character_, stringsAsFactors = FALSE)
  } else if (isTRUE(include_legacy_aliases)) {
    legacy <- wgcna_marker_sets()
    alias_map <- c(
      microglia = "canonical_microglia_homeostatic",
      neuronal_synaptic_neuropil = "canonical_neuronal_synaptic_neuropil",
      neuropil_synaptic_neuronal = "canonical_neuronal_synaptic_neuropil",
      nuclear_soma = "canonical_neuronal_soma_nuclear",
      astrocyte = "canonical_astrocyte",
      oligodendrocyte_myelin = "canonical_oligodendrocyte_myelin",
      endothelial_pericyte_vascular = "canonical_endothelial_vascular",
      mitochondrial_oxphos = "canonical_mitochondrial_oxphos",
      ribosomal_translation = "canonical_ribosomal_translation",
      rnp_rna_processing = "canonical_rnp_rna_processing"
    )
    for (alias in names(alias_map)) {
      target <- unname(alias_map[[alias]])
      if (!alias %in% names(sets)) sets[[alias]] <- sets[[target]] %||% legacy[[alias]]
    }
  }

  sets <- sets[!duplicated(names(sets))]
  attr(sets, "marker_source_metadata") <- metadata
  attr(sets, "marker_registry_version") <- if (!is.null(registry)) attr(registry, "marker_registry_version") else NA_character_
  attr(sets, "empirical_marker_set_version") <- if (!is.null(empirical)) attr(empirical, "empirical_marker_set_version") else NA_character_
  sets
}

standardize_wgcna_metadata <- function(meta, dataset) {
  meta <- as.data.frame(meta, check.names = FALSE, stringsAsFactors = FALSE)
  sample_col <- first_present_col(meta, c("Sample", "sample", "SampleID", "sample_id", "row.names"))
  if (is.na(sample_col)) meta$Sample <- rownames(meta) else meta$Sample <- as.character(meta[[sample_col]])
  meta$AnimalID <- infer_wgcna_animal_id(meta, meta$Sample)
  for (target in c("Region", "Layer", "Sex", "Batch")) {
    col <- first_present_col(meta, c(target, tolower(target), if (target == "Batch") c("plate", "run", "batch_id") else character()))
    meta[[target]] <- if (!is.na(col)) as.character(meta[[col]]) else NA_character_
  }
  group_col <- first_present_col(meta, c("StressGroup", "ExpGroup", "condition", "Group", "group", "group2"))
  meta$StressGroup <- if (!is.na(group_col)) toupper(as.character(meta[[group_col]])) else NA_character_
  meta$StressGroup <- dplyr::case_when(
    grepl("^CON|^CTRL|CONTROL", meta$StressGroup, ignore.case = TRUE) ~ "CON",
    grepl("^RES", meta$StressGroup, ignore.case = TRUE) ~ "RES",
    grepl("^SUS", meta$StressGroup, ignore.case = TRUE) ~ "SUS",
    TRUE ~ meta$StressGroup
  )
  if (!"ExpGroup" %in% names(meta)) meta$ExpGroup <- meta$StressGroup
  meta$RegionLayer <- ifelse(!is.na(meta$Region) & nzchar(meta$Region) & !is.na(meta$Layer) & nzchar(meta$Layer), paste(meta$Region, meta$Layer, sep = "_"), NA_character_)
  spatial_col <- if (dataset == "neuron_neuropil" && any(!is.na(meta$RegionLayer))) "RegionLayer" else "Region"
  meta$SpatialUnit <- spatial_col
  meta$SpatialLabel <- as.character(meta[[spatial_col]])
  meta
}

resolve_wgcna_files <- function(dataset) {
  list(
    state = path_processed("06_modules_WGCNA", "01_WGCNA", dataset, "wgcna_final_model_state.rds"),
    definitions = path_results("tables", "06_modules_WGCNA", "01_WGCNA", dataset, "modules", "WGCNA_module_definitions_for_downstream.csv"),
    module_summary = path_results("tables", "06_modules_WGCNA", "01_WGCNA", dataset, "modules", "WGCNA_module_summary.csv"),
    go = path_results("tables", "06_modules_WGCNA", "01_WGCNA", dataset, "modules", "WGCNA_module_GO_enrichment_long.csv"),
    supermodule_annotation = path_results("tables", "06_modules_WGCNA", "01_WGCNA", dataset, "supermodules", "wgcna_module_supermodule_annotation.csv"),
    supermodule_summary = path_results("tables", "06_modules_WGCNA", "01_WGCNA", dataset, "supermodules", "wgcna_supermodule_summary.csv"),
    marker_traits = path_results("tables", "03_qc_exploration", "06_wgcna_marker_trait_export", dataset, "wgcna_marker_traits_by_sample.csv"),
    neuropil_annotation = path_results("tables", "04_differential_expression_enrichment", "neuropil_reference_annotation", "microglia", "microglia_neuropil_annotation_latest.csv")
  )
}

load_wgcna_state <- function(path) {
  if (!file.exists(path)) stop("Missing WGCNA final state: ", path, call. = FALSE)
  readRDS(path)
}

extract_module_eigengenes <- function(state) {
  MEs <- state$mergedMEs %||% state$MEs %||% state$moduleEigengenes
  if (is.null(MEs)) stop("WGCNA state does not contain mergedMEs/module eigengenes.", call. = FALSE)
  MEs <- as.data.frame(MEs, check.names = FALSE)
  if (!"Sample" %in% names(MEs)) MEs <- tibble::rownames_to_column(MEs, "Sample")
  MEs
}

module_col_to_id <- function(x) {
  out <- sub("^ME", "", as.character(x))
  out
}

make_supermodule_eigengenes <- function(module_eigengenes, super_map) {
  me_cols <- setdiff(names(module_eigengenes), "Sample")
  super_map <- super_map[super_map$module_eigengene %in% me_cols & !is.na(super_map$SupermoduleID), , drop = FALSE]
  rows <- list(Sample = module_eigengenes$Sample)
  comp <- list()
  for (sid in unique(super_map$SupermoduleID)) {
    members <- unique(super_map$module_eigengene[super_map$SupermoduleID == sid])
    vals <- module_eigengenes[, members, drop = FALSE]
    vals <- vals[, vapply(vals, function(z) stats::var(as.numeric(z), na.rm = TRUE) > 0, logical(1)), drop = FALSE]
    if (!ncol(vals)) next
    if (ncol(vals) == 1L) {
      score <- as.numeric(vals[[1]])
    } else {
      pc <- stats::prcomp(vals, center = TRUE, scale. = TRUE)$x[, 1L]
      mean_vec <- rowMeans(vals, na.rm = TRUE)
      score <- if (stats::cor(pc, mean_vec, use = "pairwise.complete.obs") < 0) -pc else pc
    }
    col <- paste0("SM__", safe_filename(sid))
    rows[[col]] <- as.numeric(score)
    comp[[length(comp) + 1L]] <- data.frame(
      supermodule_id = sid,
      supermodule_eigengene = col,
      n_member_modules = length(members),
      member_modules = paste(members, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
  list(eigengenes = as.data.frame(rows, check.names = FALSE), composition = dplyr::bind_rows(comp))
}

empty_group_effects <- function(dataset, level, reason) {
  data.frame(
    dataset = dataset, level = level, endpoint_id = NA_character_, endpoint_label = NA_character_,
    module_id = NA_character_, supermodule_id = NA_character_,
    module_label = NA_character_, supermodule_label = NA_character_, spatial_unit = NA_character_,
    effect_scope = NA_character_, SpatialUnitType = NA_character_, model_type = NA_character_,
    has_repeated_animals = NA, n_animals = NA_integer_,
    n_animals_total = NA_integer_, min_animals_per_group = NA_integer_,
    animal_level_status = "mixed_or_unclear",
    contrast = NA_character_, estimate = NA_real_, SE = NA_real_, statistic = NA_real_,
    p_value = NA_real_, FDR_within_dataset_level = NA_real_, FDR_global = NA_real_,
    evidence_status = "not_supported",
    direction = NA_character_, n_samples = 0L, formula_requested = NA_character_, formula_used = NA_character_,
    dropped_covariates = NA_character_, rank_deficient_model = NA, model_warning = reason,
    stringsAsFactors = FALSE
  )
}

required_group_effect_columns <- c(
  "dataset", "level", "endpoint_id", "endpoint_label", "module_id", "supermodule_id", "module_label", "supermodule_label",
  "spatial_unit", "effect_scope", "SpatialUnitType", "model_type", "has_repeated_animals",
  "n_animals", "n_animals_total", "min_animals_per_group", "animal_level_status",
  "contrast", "estimate", "SE", "statistic", "p_value",
  "FDR_within_dataset_level", "FDR_global", "evidence_status", "direction", "n_samples",
  "formula_requested", "formula_used", "dropped_covariates", "rank_deficient_model", "model_warning"
)

required_module_annotation_columns <- c(
  "dataset", "ModuleID", "ModuleColor", "n_proteins", "microenvironment_class", "interpretation_note"
)

required_interpretable_columns <- c(
  "dataset", "level", "contrast", "estimate", "p_value", "FDR_global", "interpretation_sentence"
)

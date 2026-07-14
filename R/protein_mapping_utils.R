# Shared protein identifier mapping helpers.

normalize_token <- function(x) {
  x <- as.character(x)
  x <- sub("\\s.*$", "", x)
  x <- gsub("\\|\\|+", "|", x)
  x <- toupper(gsub("\\s+", "", x))
  x <- gsub("\\u00A0", "", x)
  x <- gsub("\\.+", ".", x)
  x <- gsub("__+", "_", x)
  x
}

to_base_no_iso_mouse <- function(x) {
  x <- gsub("-\\d+$", "", as.character(x))
  gsub("_MOUSE$", "", x)
}

is_uniprot_ac <- function(x) {
  grepl("^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9]|A0A[0-9A-Z]{7})$", as.character(x))
}

extract_ac <- function(s) {
  s <- as.character(s)
  pattern <- "(?i)(?:^|\\|)([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9]|A0A[0-9A-Z]{7})(?:\\-|\\||$|[^A-Z0-9])"
  if (requireNamespace("stringr", quietly = TRUE)) {
    m <- stringr::str_match(s, pattern)
    out <- ifelse(is.na(m[, 2]), NA_character_, toupper(m[, 2]))
  } else {
    out <- regmatches(s, regexpr(pattern, s, perl = TRUE, ignore.case = TRUE))
    out <- sub(".*?([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9]|A0A[0-9A-Z]{7}).*", "\\1", out, ignore.case = TRUE)
    out[!nzchar(out)] <- NA_character_
    out <- toupper(out)
  }
  gsub("-\\d+$", "", out)
}

extract_entry <- function(s) {
  s <- as.character(s)
  if (requireNamespace("stringr", quietly = TRUE)) {
    m <- stringr::str_match(s, "(?i)(?:^|\\|)([A-Z0-9]+_MOUSE)(?:\\||$|\\s)")
    out <- ifelse(is.na(m[, 2]), NA_character_, m[, 2])
    if (all(is.na(out))) {
      m2 <- stringr::str_match(s, "(?i)\\b([A-Z0-9]+_MOUSE)\\b")
      out <- m2[, 2]
    }
    return(toupper(out))
  }
  out <- regmatches(s, regexpr("\\b[A-Z0-9]+_MOUSE\\b", toupper(s), perl = TRUE))
  out[!nzchar(out)] <- NA_character_
  toupper(out)
}

load_mouse_idmapping <- function(path = NULL, auto_download = FALSE) {
  if (is.null(path) || !nzchar(path)) {
    if (exists("path_external", mode = "function")) {
      path <- path_external("MOUSE_10090_idmapping.dat")
    } else {
      path <- "MOUSE_10090_idmapping.dat"
    }
  }
  if (!file.exists(path) && isTRUE(auto_download)) {
    gz_url <- "https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/idmapping/by_organism/MOUSE_10090_idmapping.dat.gz"
    gz_file <- paste0(path, ".gz")
    options(timeout = 3600)
    utils::download.file(gz_url, gz_file, mode = "wb")
    if (!requireNamespace("R.utils", quietly = TRUE)) stop("R.utils is required to unzip UniProt idmapping.", call. = FALSE)
    R.utils::gunzip(gz_file, destname = path, remove = TRUE)
  }
  if (!file.exists(path)) stop("UniProt mapping file not found: ", path, call. = FALSE)
  if (requireNamespace("readr", quietly = TRUE)) {
    readr::read_tsv(path, col_names = c("UniProt_Accession", "Type", "Value"), col_types = "ccc", quote = "", progress = FALSE)
  } else {
    utils::read.delim(path, header = FALSE, col.names = c("UniProt_Accession", "Type", "Value"), stringsAsFactors = FALSE, quote = "")
  }
}

build_mouse_maps <- function(uniprot_mapping) {
  entry_map <- uniprot_mapping |>
    dplyr::filter(.data$Type == "UniProtKB-ID") |>
    dplyr::transmute(
      UNIPROT = toupper(trimws(.data$UniProt_Accession)),
      entry_full = toupper(trimws(.data$Value)),
      entry_base = toupper(gsub("_MOUSE$", "", trimws(.data$Value)))
    ) |>
    dplyr::filter(grepl("_MOUSE\\s*$", .data$entry_full), nzchar(.data$UNIPROT)) |>
    dplyr::distinct(.data$entry_base, .keep_all = TRUE)

  gene_map <- uniprot_mapping |>
    dplyr::filter(.data$Type %in% c("Gene_Name", "Gene_Name(synonym)", "Gene_Synonym")) |>
    dplyr::transmute(
      primaryAccession = toupper(trimws(.data$UniProt_Accession)),
      input = toupper(trimws(.data$Value))
    ) |>
    dplyr::filter(nzchar(.data$input), nzchar(.data$primaryAccession)) |>
    dplyr::mutate(pref = !startsWith(.data$primaryAccession, "A0A")) |>
    dplyr::arrange(dplyr::desc(.data$pref), .data$primaryAccession, .data$input) |>
    dplyr::group_by(.data$input) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup() |>
    dplyr::select(dplyr::all_of(c("input", "primaryAccession")))

  gene_map <- dplyr::bind_rows(
    gene_map,
    entry_map |> dplyr::transmute(input = .data$entry_base, primaryAccession = .data$UNIPROT)
  ) |>
    dplyr::distinct(.data$input, .keep_all = TRUE)

  accession_gene_map <- uniprot_mapping |>
    dplyr::filter(.data$Type == "Gene_Name") |>
    dplyr::transmute(
      UNIPROT = toupper(trimws(.data$UniProt_Accession)),
      gene_symbol = toupper(trimws(.data$Value))
    ) |>
    dplyr::filter(nzchar(.data$UNIPROT), nzchar(.data$gene_symbol)) |>
    dplyr::distinct(.data$UNIPROT, .keep_all = TRUE)

  reviewed_map <- uniprot_mapping |>
    dplyr::filter(.data$Type %in% c("UniProtKB-ID", "Reviewed", "Status")) |>
    dplyr::transmute(
      UNIPROT = toupper(trimws(.data$UniProt_Accession)),
      reviewed_status = dplyr::case_when(
        .data$Type == "Reviewed" ~ as.character(.data$Value),
        startsWith(toupper(trimws(.data$UniProt_Accession)), "A0A") ~ "unreviewed_inferred",
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::filter(nzchar(.data$UNIPROT)) |>
    dplyr::group_by(.data$UNIPROT) |>
    dplyr::summarise(reviewed_status = dplyr::first(stats::na.omit(.data$reviewed_status)), .groups = "drop")

  list(
    entry_map = entry_map,
    gene_map = gene_map,
    accession_gene_map = accession_gene_map,
    reviewed_map = reviewed_map
  )
}

map_token_to_mouse_accession <- function(token, entry_map, gene_map) {
  if (is.na(token) || !nzchar(trimws(token))) return(NA_character_)
  token <- unlist(strsplit(as.character(token), ";", fixed = TRUE), use.names = FALSE)[1]
  token_up <- normalize_token(token)
  token_base <- to_base_no_iso_mouse(token_up)
  ac_guess <- extract_ac(token)
  if (!is.na(ac_guess) && nzchar(ac_guess)) return(ac_guess)
  if (is_uniprot_ac(token_base)) return(token_base)
  hit_entry <- entry_map$UNIPROT[match(toupper(token_base), entry_map$entry_base)]
  if (!is.na(hit_entry) && nzchar(hit_entry)) return(hit_entry)
  hit_gene <- gene_map$primaryAccession[match(toupper(token_base), gene_map$input)]
  if (!is.na(hit_gene) && nzchar(hit_gene)) return(hit_gene)
  NA_character_
}

read_manual_mapping_table <- function(path = Sys.getenv("PROTEOMICS_MANUAL_MAPPING_FILE", unset = "")) {
  if (!nzchar(path)) {
    path <- if (exists("path_metadata", mode = "function")) path_metadata("manual_mapping.xlsx") else file.path("data", "metadata", "manual_mapping.xlsx")
  }
  if (!file.exists(path)) {
    return(structure(NULL, path = path, status = "missing"))
  }
  ext <- tolower(tools::file_ext(path))
  mm <- tryCatch({
    if (ext %in% c("xlsx", "xls")) {
      if (!requireNamespace("readxl", quietly = TRUE)) stop("readxl is required for manual mapping workbooks.")
      readxl::read_excel(path, sheet = 1)
    } else if (requireNamespace("readr", quietly = TRUE)) {
      readr::read_csv(path, show_col_types = FALSE)
    } else {
      utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    }
  }, error = function(e) NULL)
  if (is.null(mm) || !is.data.frame(mm) || !nrow(mm)) {
    return(structure(NULL, path = path, status = "empty_or_unreadable"))
  }
  names(mm) <- tolower(gsub("[^a-z0-9]+", "_", trimws(names(mm))))
  input_col <- intersect(c("gene_symbol", "input", "source_id", "original", "original_symbol", "symbol", "token_raw"), names(mm))
  mapped_col <- intersect(c("mapped_gene_symbol", "mapped", "mapped_id", "final_accession", "accession", "uniprot", "uniprot_accession"), names(mm))
  if (!length(input_col) || !length(mapped_col)) {
    attr(mm, "path") <- path
    attr(mm, "status") <- "missing_required_columns"
    return(mm[0, , drop = FALSE])
  }
  out <- mm |>
    dplyr::transmute(
      gene_symbol = toupper(trimws(as.character(.data[[input_col[1]]]))),
      mapped_gene_symbol = toupper(trimws(as.character(.data[[mapped_col[1]]]))),
      manual_mapping_source_column = input_col[1],
      manual_mapping_target_column = mapped_col[1]
    ) |>
    dplyr::filter(!is.na(.data$gene_symbol), nzchar(.data$gene_symbol), !is.na(.data$mapped_gene_symbol), nzchar(.data$mapped_gene_symbol)) |>
    dplyr::distinct(.data$gene_symbol, .keep_all = TRUE)
  attr(out, "path") <- path
  attr(out, "status") <- "loaded"
  out
}

apply_manual_mapping_override <- function(resolved, manual_mapping, entry_map, gene_map,
                                          token_col = "token_raw", base_col = "token_base",
                                          resolved_col = "Resolved_UNIPROT", strategy_col = "strategy",
                                          override = TRUE) {
  if (!"manual_mapping_used" %in% names(resolved)) resolved$manual_mapping_used <- FALSE
  empty_audit <- resolved[0, intersect(c(token_col, base_col, resolved_col, strategy_col, "manual_mapping_used"), names(resolved)), drop = FALSE]
  if (is.null(manual_mapping) || !nrow(manual_mapping)) {
    return(list(data = resolved, audit = empty_audit))
  }
  map_to_acc <- function(vals) {
    vals <- toupper(trimws(as.character(vals)))
    out <- ifelse(is_uniprot_ac(vals), vals, NA_character_)
    need <- is.na(out) | !nzchar(out)
    if (any(need)) {
      base <- to_base_no_iso_mouse(vals[need])
      hit <- entry_map$UNIPROT[match(base, entry_map$entry_base)]
      ok <- !is.na(hit) & nzchar(hit)
      if (any(ok)) out[which(need)[ok]] <- hit[ok]
    }
    need <- is.na(out) | !nzchar(out)
    if (any(need) && nrow(gene_map)) {
      key <- toupper(to_base_no_iso_mouse(vals[need]))
      hit <- gene_map$primaryAccession[match(key, gene_map$input)]
      ok <- !is.na(hit) & nzchar(hit)
      if (any(ok)) out[which(need)[ok]] <- hit[ok]
    }
    out
  }

  mm <- manual_mapping |>
    dplyr::mutate(
      manual_input_norm = normalize_token(.data$gene_symbol),
      manual_input_base = to_base_no_iso_mouse(.data$manual_input_norm),
      manual_mapped_accession = map_to_acc(.data$mapped_gene_symbol)
    ) |>
    dplyr::filter(!is.na(.data$manual_mapped_accession), nzchar(.data$manual_mapped_accession))
  if (!nrow(mm)) return(list(data = resolved, audit = empty_audit))

  token_norm <- normalize_token(resolved[[token_col]])
  token_base <- if (base_col %in% names(resolved)) toupper(as.character(resolved[[base_col]])) else to_base_no_iso_mouse(token_norm)
  hit <- match(token_norm, mm$manual_input_norm)
  missing_hit <- is.na(hit)
  hit[missing_hit] <- match(token_base[missing_hit], mm$manual_input_base)
  idx <- which(!is.na(hit))
  if (!length(idx)) return(list(data = resolved, audit = empty_audit))
  if (!isTRUE(override)) {
    idx <- idx[is.na(resolved[[resolved_col]][idx]) | !nzchar(resolved[[resolved_col]][idx])]
  }
  if (!length(idx)) return(list(data = resolved, audit = empty_audit))

  previous_accession <- resolved[[resolved_col]][idx]
  previous_strategy <- resolved[[strategy_col]][idx]
  mapped <- mm$manual_mapped_accession[hit[idx]]
  resolved[[resolved_col]][idx] <- mapped
  resolved[[strategy_col]][idx] <- ifelse(
    is.na(previous_strategy) | !nzchar(previous_strategy),
    "manual_mapping",
    paste0(previous_strategy, "|manual_mapping")
  )
  resolved$manual_mapping_used[idx] <- TRUE
  audit <- data.frame(
    original_token = as.character(resolved[[token_col]][idx]),
    token_base = as.character(token_base[idx]),
    manual_input = mm$gene_symbol[hit[idx]],
    manual_mapped_gene_symbol = mm$mapped_gene_symbol[hit[idx]],
    previous_accession = previous_accession,
    resolved_uniprot = mapped,
    previous_strategy = previous_strategy,
    mapping_strategy = resolved[[strategy_col]][idx],
    manual_mapping_used = TRUE,
    stringsAsFactors = FALSE
  )
  list(data = resolved, audit = audit)
}

protein_group_ambiguity_levels <- function() {
  c(
    "single_accession_single_gene",
    "multi_accession_same_gene",
    "multi_gene_indistinguishable",
    "explicit_master_with_subordinate_members",
    "mixed_species_or_contaminant",
    "partially_mapped_group",
    "unresolved_group"
  )
}

split_protein_group_members <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  # Avoid splitting UniProt pipe notation such as sp|Q9CQH5|ENTRY_MOUSE.
  pieces <- unlist(strsplit(x, "\\s*(?:;|,(?=\\s)|\\s/\\s)\\s*", perl = TRUE), use.names = FALSE)
  pieces <- trimws(pieces)
  pieces[nzchar(pieces)]
}

normalize_member_identifier <- function(x) normalize_token(trimws(as.character(x)))

canonical_member_set <- function(members) {
  norm <- normalize_member_identifier(members)
  norm <- norm[nzchar(norm)]
  sort(unique(norm), method = "radix")
}

detect_source_feature_columns <- function(df) {
  nms <- names(df)
  first_match <- function(candidates) {
    hit <- intersect(candidates, nms)
    if (length(hit)) hit[1] else NA_character_
  }
  list(
    feature_id = first_match(c(
      "ProteinGroupID", "protein_group_id", "feature_id", "FeatureID",
      "Protein.Group", "Protein Group", "T: Protein.Group", "row_id", "id"
    )),
    original_identifier = first_match(c(
      "original_identifier", "gene_symbol", "Protein.Group", "Protein Group",
      "T: Protein.Group", "id"
    )),
    explicit_master = first_match(c(
      "master_protein", "master_protein_id", "Master.Protein", "Leading.Protein",
      "leading_protein", "Primary.Protein", "is_master_protein"
    )),
    unique_peptide_count = first_match(c(
      "unique_peptide_count", "Unique.Peptides", "Unique Peptides", "unique_peptides"
    )),
    razor_peptide_count = first_match(c(
      "razor_peptide_count", "Razor.Peptides", "Razor Peptides", "razor_peptides"
    )),
    sequence_coverage = first_match(c(
      "sequence_coverage", "Sequence.Coverage", "Sequence Coverage", "Coverage"
    )),
    identification_confidence = first_match(c(
      "identification_confidence", "Identification.Confidence", "Confidence", "Score"
    ))
  )
}

stable_pg_hash <- function(x) {
  f <- tempfile("protein_group_id_")
  on.exit(unlink(f), add = TRUE)
  writeLines(enc2utf8(as.character(x)), f, useBytes = TRUE)
  substr(toupper(unname(tools::md5sum(f))), 1, 16)
}

parse_member_identifier <- function(member_identifier) {
  token_norm <- normalize_member_identifier(member_identifier)
  entry <- extract_entry(member_identifier)
  accession <- extract_ac(member_identifier)
  species <- dplyr::case_when(
    grepl("_MOUSE\\b", token_norm) ~ "mouse",
    grepl("_HUMAN\\b", token_norm) ~ "human",
    grepl("_RAT\\b", token_norm) ~ "rat",
    TRUE ~ NA_character_
  )
  contaminant <- grepl("(^|[_|])CON(__|_)|(^|[_|])REV(__|_)|CONTAM|KERATIN|TRYPSIN", token_norm)
  list(
    member_identifier_normalized = token_norm,
    parsed_accession = accession,
    parsed_entry = entry,
    token_base = to_base_no_iso_mouse(ifelse(!is.na(entry), entry, token_norm)),
    member_species = species,
    contaminant_status = ifelse(contaminant, "contaminant", "not_flagged"),
    isoform = ifelse(grepl("-\\d+$", token_norm), sub(".*(-\\d+)$", "\\1", token_norm), NA_character_)
  )
}

lookup_accession_gene <- function(accession, accession_gene_map = NULL, gene_map = NULL) {
  accession <- toupper(as.character(accession))
  if (!is.null(accession_gene_map) && nrow(accession_gene_map)) {
    hit <- accession_gene_map$gene_symbol[match(accession, accession_gene_map$UNIPROT)]
    if (length(hit) && !is.na(hit) && nzchar(hit)) return(hit)
  }
  if (!is.null(gene_map) && nrow(gene_map)) {
    hit <- gene_map$input[match(accession, gene_map$primaryAccession)]
    if (length(hit) && !is.na(hit) && nzchar(hit)) return(hit)
  }
  NA_character_
}

resolve_protein_group_member <- function(member_identifier, entry_map, gene_map,
                                         accession_gene_map = NULL, reviewed_map = NULL,
                                         manual_mapping = NULL, manual_override = TRUE) {
  parsed <- parse_member_identifier(member_identifier)
  acc <- parsed$parsed_accession
  strategy <- NA_character_
  manual_used <- FALSE

  if (!is.na(acc) && nzchar(acc)) {
    strategy <- "accept_accession_in_member"
  }

  if ((is.na(acc) || !nzchar(acc)) && is_uniprot_ac(parsed$token_base)) {
    acc <- parsed$token_base
    strategy <- "accept_accession_base"
  }

  if ((is.na(acc) || !nzchar(acc)) && !is.na(parsed$parsed_entry) && nrow(entry_map)) {
    hit <- entry_map$UNIPROT[match(to_base_no_iso_mouse(parsed$parsed_entry), entry_map$entry_base)]
    if (!is.na(hit) && nzchar(hit)) {
      acc <- hit
      strategy <- "entry_local_mouse"
    }
  }

  if ((is.na(acc) || !nzchar(acc)) && nrow(gene_map)) {
    hit <- gene_map$primaryAccession[match(toupper(parsed$token_base), gene_map$input)]
    if (!is.na(hit) && nzchar(hit)) {
      acc <- hit
      strategy <- "gene_local_mouse"
    }
  }

  if (!is.null(manual_mapping) && nrow(manual_mapping)) {
    manual <- data.frame(
      token_raw = member_identifier,
      token_base = parsed$token_base,
      Resolved_UNIPROT = acc,
      strategy = strategy,
      manual_mapping_used = FALSE,
      stringsAsFactors = FALSE
    )
    manual_out <- apply_manual_mapping_override(manual, manual_mapping, entry_map, gene_map, override = manual_override)$data
    if (isTRUE(manual_out$manual_mapping_used[1])) {
      acc <- manual_out$Resolved_UNIPROT[1]
      strategy <- manual_out$strategy[1]
      manual_used <- TRUE
    }
  }

  gene <- if (!is.na(acc) && nzchar(acc)) lookup_accession_gene(acc, accession_gene_map, gene_map) else NA_character_
  reviewed_status <- NA_character_
  if (!is.null(reviewed_map) && nrow(reviewed_map) && !is.na(acc) && nzchar(acc)) {
    reviewed_status <- reviewed_map$reviewed_status[match(acc, reviewed_map$UNIPROT)]
  }

  data.frame(
    member_identifier_original = as.character(member_identifier),
    member_identifier_normalized = parsed$member_identifier_normalized,
    member_accession = ifelse(!is.na(acc) && nzchar(acc), toupper(acc), NA_character_),
    member_gene_symbol = gene,
    member_species = parsed$member_species,
    member_mapping_status = ifelse(!is.na(acc) && nzchar(acc), "mapped", "unmapped"),
    mapping_strategy = ifelse(!is.na(strategy) && nzchar(strategy), strategy, "unmapped"),
    manual_mapping_used = manual_used,
    reviewed_status = reviewed_status,
    contaminant_status = parsed$contaminant_status,
    stringsAsFactors = FALSE
  )
}

classify_protein_group <- function(member_bridge, explicit_master_present = FALSE) {
  mapped <- member_bridge[member_bridge$member_mapping_status == "mapped", , drop = FALSE]
  n_mapped <- length(unique(stats::na.omit(mapped$member_accession)))
  genes <- unique(stats::na.omit(mapped$member_gene_symbol))
  contaminant <- any(member_bridge$contaminant_status == "contaminant", na.rm = TRUE)
  non_mouse <- any(!is.na(member_bridge$member_species) & member_bridge$member_species != "mouse")
  partially_mapped <- n_mapped > 0 && any(member_bridge$member_mapping_status != "mapped")

  if (contaminant || non_mouse) return("mixed_species_or_contaminant")
  if (n_mapped == 0) return("unresolved_group")
  if (partially_mapped) return("partially_mapped_group")
  if (isTRUE(explicit_master_present) && n_mapped > 1) return("explicit_master_with_subordinate_members")
  if (n_mapped == 1 && length(genes) <= 1) return("single_accession_single_gene")
  if (length(genes) == 1) return("multi_accession_same_gene")
  "multi_gene_indistinguishable"
}

protein_group_claim_rules <- function(ambiguity_class) {
  gene_allowed <- ambiguity_class %in% c("single_accession_single_gene", "multi_accession_same_gene")
  protein_allowed <- ambiguity_class %in% c("single_accession_single_gene")
  data.frame(
    protein_group_ambiguity_class = ambiguity_class,
    gene_level_claim_allowed = gene_allowed,
    protein_level_claim_allowed = protein_allowed,
    stringsAsFactors = FALSE
  )
}

select_representative_member <- function(member_bridge, row = NULL, source_cols = list()) {
  mapped <- member_bridge[member_bridge$member_mapping_status == "mapped", , drop = FALSE]
  if (!nrow(mapped)) {
    return(list(accession = NA_character_, gene_symbol = NA_character_, rule = "none_unmapped"))
  }

  explicit_master_col <- source_cols$explicit_master
  if (!is.null(row) && length(explicit_master_col) && !is.na(explicit_master_col) && explicit_master_col %in% names(row)) {
    master_value <- as.character(row[[explicit_master_col]][1])
    if (!is.na(master_value) && nzchar(master_value)) {
      master_members <- split_protein_group_members(master_value)
      master_norm <- normalize_member_identifier(master_members)
      hit <- mapped[mapped$member_identifier_normalized %in% master_norm | mapped$member_accession %in% toupper(master_norm), , drop = FALSE]
      if (nrow(hit)) return(list(accession = hit$member_accession[1], gene_symbol = hit$member_gene_symbol[1], rule = "explicit_upstream_master"))
    }
  }

  reviewed <- mapped[!is.na(mapped$reviewed_status) & grepl("reviewed", mapped$reviewed_status, ignore.case = TRUE), , drop = FALSE]
  if (nrow(reviewed)) {
    return(list(accession = reviewed$member_accession[order(reviewed$member_accession)][1], gene_symbol = reviewed$member_gene_symbol[order(reviewed$member_accession)][1], rule = "reviewed_uniprot_status"))
  }

  ord <- order(mapped$member_accession)
  list(
    accession = mapped$member_accession[ord][1],
    gene_symbol = mapped$member_gene_symbol[ord][1],
    rule = "canonical_accession_order_display_only"
  )
}

validate_protein_group_id_collisions <- function(wide, strict = TRUE) {
  collisions <- wide |>
    dplyr::count(.data$ProteinGroupID, .data$source_file, name = "n_rows") |>
    dplyr::filter(.data$n_rows > 1)
  if (nrow(collisions) && isTRUE(strict)) {
    stop("ProteinGroupID collision or unstable identifier detected; refusing to repair with make.unique().", call. = FALSE)
  }
  collisions
}

build_canonical_protein_group_tables <- function(df_raw, dataset, source_file,
                                                 entry_map, gene_map,
                                                 accession_gene_map = NULL,
                                                 reviewed_map = NULL,
                                                 manual_mapping = NULL,
                                                 manual_override = TRUE,
                                                 strict = TRUE,
                                                 identifier_col = NULL,
                                                 feature_col = NULL) {
  source_cols <- detect_source_feature_columns(df_raw)
  if (is.null(identifier_col)) identifier_col <- source_cols$original_identifier
  if (is.na(identifier_col)) identifier_col <- names(df_raw)[1]
  if (is.null(feature_col)) feature_col <- source_cols$feature_id
  if (length(feature_col) == 0L) feature_col <- NA_character_
  if (!identifier_col %in% names(df_raw)) {
    stop("Canonical protein-group identifier column is missing: ", identifier_col, call. = FALSE)
  }

  rows <- vector("list", nrow(df_raw))
  bridges <- vector("list", nrow(df_raw))

  for (i in seq_len(nrow(df_raw))) {
    row <- df_raw[i, , drop = FALSE]
    original_identifier <- as.character(row[[identifier_col]][1])
    members_original <- split_protein_group_members(original_identifier)
    members_canonical <- canonical_member_set(members_original)
    source_feature_id <- if (!is.na(feature_col) && feature_col %in% names(row)) as.character(row[[feature_col]][1]) else NA_character_

    member_tbls <- lapply(members_original, resolve_protein_group_member,
                          entry_map = entry_map, gene_map = gene_map,
                          accession_gene_map = accession_gene_map,
                          reviewed_map = reviewed_map,
                          manual_mapping = manual_mapping,
                          manual_override = manual_override)
    bridge <- if (length(member_tbls)) dplyr::bind_rows(member_tbls) else data.frame()
    if (!nrow(bridge)) {
      bridge <- data.frame(
        member_identifier_original = NA_character_,
        member_identifier_normalized = NA_character_,
        member_accession = NA_character_,
        member_gene_symbol = NA_character_,
        member_species = NA_character_,
        member_mapping_status = "unmapped",
        mapping_strategy = "unmapped",
        manual_mapping_used = FALSE,
        reviewed_status = NA_character_,
        contaminant_status = "not_flagged",
        stringsAsFactors = FALSE
      )
    }
    bridge$member_rank_original <- seq_len(nrow(bridge))
    canonical_rank <- match(bridge$member_identifier_normalized, members_canonical)
    bridge$member_rank_canonical <- canonical_rank

    explicit_master_present <- !is.na(source_cols$explicit_master) &&
      source_cols$explicit_master %in% names(row) &&
      !is.na(row[[source_cols$explicit_master]][1]) &&
      nzchar(as.character(row[[source_cols$explicit_master]][1]))
    ambiguity_class <- classify_protein_group(bridge, explicit_master_present)
    claims <- protein_group_claim_rules(ambiguity_class)
    representative <- select_representative_member(bridge, row, source_cols)

    feature_is_membership <- !is.na(feature_col) && feature_col %in%
      c("Protein.Group", "Protein Group", "T: Protein.Group")
    identity_members <- ifelse(
      !is.na(bridge$member_accession) & nzchar(bridge$member_accession),
      bridge$member_accession,
      bridge$member_identifier_normalized
    )
    identity_members <- sort(unique(identity_members[!is.na(identity_members) & nzchar(identity_members)]), method = "radix")
    membership_basis <- paste(identity_members, collapse = ";")
    source_membership_basis <- if (feature_is_membership && !is.na(source_feature_id) && nzchar(source_feature_id)) {
      paste(canonical_member_set(split_protein_group_members(source_feature_id)), collapse = ";")
    } else {
      ""
    }
    id_basis <- if (!is.na(source_feature_id) && nzchar(source_feature_id) && !feature_is_membership) {
      paste(dataset, "source_feature_id", source_feature_id, sep = "|")
    } else if (feature_is_membership && nzchar(source_membership_basis) && !identical(source_membership_basis, membership_basis)) {
      paste(dataset, "canonical_members", membership_basis, "source_membership_discriminator", source_membership_basis, sep = "|")
    } else {
      paste(dataset, "canonical_members", membership_basis, sep = "|")
    }
    protein_group_id <- if (!is.na(feature_col) && identical(feature_col, "ProteinGroupID")) {
      source_feature_id
    } else {
      paste0("PG:", dataset, ":", stable_pg_hash(id_basis))
    }
    if (is.na(protein_group_id) || !nzchar(protein_group_id)) {
      stop("Missing ProteinGroupID; canonical feature identity cannot be repaired from row order.", call. = FALSE)
    }

    evidence_na <- list(
      unique_peptide_count = NA_real_,
      razor_peptide_count = NA_real_,
      sequence_coverage = NA_real_,
      identification_confidence = NA_character_
    )
    for (nm in names(evidence_na)) bridge[[nm]] <- evidence_na[[nm]]

    bridge$ProteinGroupID <- protein_group_id
    bridge$source_file <- basename(source_file)
    bridge$source_feature_id <- source_feature_id
    bridge$source_row_id <- i
    bridge$original_identifier <- original_identifier

    mapped_accessions <- unique(stats::na.omit(bridge$member_accession))
    mapped_genes <- unique(stats::na.omit(bridge$member_gene_symbol))
    stat_cols <- setdiff(names(row), unique(stats::na.omit(c(
      identifier_col, feature_col, source_cols$explicit_master,
      source_cols$unique_peptide_count, source_cols$razor_peptide_count,
      source_cols$sequence_coverage, source_cols$identification_confidence
    ))))
    rows[[i]] <- cbind(
      data.frame(
        ProteinGroupID = protein_group_id,
        source_file = basename(source_file),
        source_feature_id = source_feature_id,
        source_row_id = i,
        original_identifier = original_identifier,
        member_identifiers_original = paste(members_original, collapse = ";"),
        member_identifiers_canonical = paste(members_canonical, collapse = ";"),
        member_accessions = paste(mapped_accessions, collapse = ";"),
        member_gene_symbols = paste(mapped_genes, collapse = ";"),
        representative_accession = representative$accession,
        representative_gene_symbol = representative$gene_symbol,
        representative_selection_rule = representative$rule,
        n_members_original = length(members_original),
        n_members_canonical = length(members_canonical),
        n_mapped_accessions = length(mapped_accessions),
        n_unmapped_members = sum(bridge$member_mapping_status != "mapped"),
        n_gene_symbols = length(mapped_genes),
        same_gene_group = length(mapped_genes) == 1 && length(mapped_accessions) > 1,
        protein_group_ambiguity_class = ambiguity_class,
        protein_level_claim_allowed = claims$protein_level_claim_allowed,
        gene_level_claim_allowed = claims$gene_level_claim_allowed,
        mapping_status = dplyr::case_when(
          length(mapped_accessions) == 0 ~ "unmapped",
          sum(bridge$member_mapping_status != "mapped") > 0 ~ "partially_mapped",
          TRUE ~ "mapped"
        ),
        gene_symbol = ifelse(length(mapped_genes) == 1, mapped_genes, representative$gene_symbol),
        gene_symbol_compatibility_status = "deprecated_display_only",
        stringsAsFactors = FALSE
      ),
      row[, stat_cols, drop = FALSE]
    )
    bridges[[i]] <- bridge
  }

  wide <- dplyr::bind_rows(rows)
  bridge <- dplyr::bind_rows(bridges) |>
    dplyr::select(dplyr::all_of(c(
      "ProteinGroupID", "source_file", "source_feature_id", "source_row_id",
      "original_identifier",
      "member_identifier_original", "member_identifier_normalized",
      "member_rank_original", "member_rank_canonical",
      "member_accession", "member_gene_symbol", "member_species",
      "member_mapping_status", "mapping_strategy", "manual_mapping_used",
      "reviewed_status", "contaminant_status",
      "unique_peptide_count", "razor_peptide_count",
      "sequence_coverage", "identification_confidence"
    )))

  missing_feature_id <- is.na(wide$source_feature_id) | !nzchar(wide$source_feature_id)
  duplicate_member_sets <- wide |>
    dplyr::filter(missing_feature_id) |>
    dplyr::count(.data$source_file, .data$member_identifiers_canonical, name = "n_rows") |>
    dplyr::filter(.data$n_rows > 1)
  if (nrow(duplicate_member_sets) && isTRUE(strict)) {
    stop("Unstable ProteinGroupID: duplicate quantitative rows share the same member set and no stable source feature identifier.", call. = FALSE)
  }
  collisions <- validate_protein_group_id_collisions(wide, strict = strict)

  summary <- data.frame(
    source_file = basename(source_file),
    dataset = dataset,
    total_quantitative_groups = nrow(wide),
    single_accession_groups = sum(wide$protein_group_ambiguity_class == "single_accession_single_gene"),
    same_gene_multi_accession_groups = sum(wide$protein_group_ambiguity_class == "multi_accession_same_gene"),
    multi_gene_groups = sum(wide$protein_group_ambiguity_class == "multi_gene_indistinguishable"),
    explicit_master_groups = sum(wide$protein_group_ambiguity_class == "explicit_master_with_subordinate_members"),
    partially_mapped_groups = sum(wide$protein_group_ambiguity_class == "partially_mapped_group"),
    unresolved_groups = sum(wide$protein_group_ambiguity_class == "unresolved_group"),
    mixed_species_or_contaminant_groups = sum(wide$protein_group_ambiguity_class == "mixed_species_or_contaminant"),
    groups_allowed_gene_level = sum(wide$gene_level_claim_allowed),
    groups_allowed_protein_level = sum(wide$protein_level_claim_allowed),
    groups_requiring_downstream_exclusion_or_sensitivity = sum(!wide$gene_level_claim_allowed | !wide$protein_level_claim_allowed),
    ProteinGroupID_collisions = nrow(collisions),
    groups_with_unstable_or_missing_source_feature_identifiers = sum(missing_feature_id),
    stringsAsFactors = FALSE
  )

  deprecated_dropped_log <- wide |>
    dplyr::filter(.data$n_members_original > 1) |>
    dplyr::transmute(
      source_file = .data$source_file,
      original_row_id = .data$source_row_id,
      original_entry = .data$original_identifier,
      kept_protein = .data$representative_accession,
      dropped_proteins = NA_character_,
      compatibility_note = "deprecated audit: all members are retained in the member bridge"
    )

  list(wide = wide, bridge = bridge, summary = summary, collision_audit = collisions, deprecated_dropped_log = deprecated_dropped_log)
}

wgcna_feature_display_label <- function(feature_table) {
  vapply(seq_len(nrow(feature_table)), function(i) {
    row <- feature_table[i, , drop = FALSE]
    pgid <- as.character(row$ProteinGroupID)
    gene <- as.character(row$representative_gene_symbol)
    genes <- as.character(row$member_gene_symbols)
    n_acc <- suppressWarnings(as.integer(row$n_mapped_accessions))
    cls <- as.character(row$protein_group_ambiguity_class)
    if (identical(cls, "single_accession_single_gene") && !is.na(gene) && nzchar(gene)) return(gene)
    if (identical(cls, "multi_accession_same_gene") && !is.na(gene) && nzchar(gene)) {
      return(paste0(gene, " (+", max(0L, n_acc - 1L), " accessions, same gene)"))
    }
    if (identical(cls, "unresolved_group")) return(paste("Unresolved protein group", pgid))
    member_label <- if (!is.na(genes) && nzchar(genes)) genes else as.character(row$member_accessions)
    paste0("Protein group ", pgid, if (!is.na(member_label) && nzchar(member_label)) paste0(": ", member_label) else "")
  }, character(1))
}

build_wgcna_canonical_features <- function(df_raw, dataset, source_file,
                                           entry_map, gene_map,
                                           accession_gene_map = NULL,
                                           reviewed_map = NULL,
                                           manual_mapping = NULL,
                                           manual_override = TRUE,
                                           sample_columns = NULL,
                                           strict = TRUE) {
  identifier_col <- if ("original_identifier" %in% names(df_raw)) {
    "original_identifier"
  } else if ("T: Protein.Names" %in% names(df_raw)) {
    "T: Protein.Names"
  } else if ("gene_symbol" %in% names(df_raw)) {
    "gene_symbol"
  } else if ("Protein.Group" %in% names(df_raw)) {
    "Protein.Group"
  } else {
    names(df_raw)[1]
  }
  feature_col <- if ("ProteinGroupID" %in% names(df_raw)) {
    "ProteinGroupID"
  } else if ("source_feature_id" %in% names(df_raw)) {
    "source_feature_id"
  } else if ("Protein.Group" %in% names(df_raw)) {
    "Protein.Group"
  } else {
    NA_character_
  }
  canonical <- build_canonical_protein_group_tables(
    df_raw = df_raw,
    dataset = dataset,
    source_file = source_file,
    entry_map = entry_map,
    gene_map = gene_map,
    accession_gene_map = accession_gene_map,
    reviewed_map = reviewed_map,
    manual_mapping = manual_mapping,
    manual_override = manual_override,
    strict = strict,
    identifier_col = identifier_col,
    feature_col = feature_col
  )
  feature_table <- canonical$wide
  feature_table$FeatureDisplayLabel <- wgcna_feature_display_label(feature_table)
  feature_table$quantitative_eligibility <- !feature_table$protein_group_ambiguity_class %in%
    "mixed_species_or_contaminant"
  feature_table$wgcna_exclusion_reason <- ifelse(
    feature_table$quantitative_eligibility,
    NA_character_,
    "mixed_species_or_contaminant"
  )
  feature_table$included_in_wgcna <- feature_table$quantitative_eligibility

  if (is.null(sample_columns)) {
    annotation_columns <- unique(c(
      "ProteinGroupID", "Protein.Group", "T: Protein.Names", "Genes",
      "First.Protein.Description", "gene_symbol", "original_identifier",
      "source_feature_id", "source_row_id", ".row_id",
      "member_identifiers_original", "member_identifiers_canonical",
      "member_accessions", "member_gene_symbols", "representative_accession",
      "representative_gene_symbol", "representative_selection_rule",
      "protein_group_ambiguity_class", "gene_level_claim_allowed",
      "protein_level_claim_allowed", "mapping_status"
    ))
    candidates <- setdiff(names(df_raw), annotation_columns)
    sample_columns <- candidates[vapply(df_raw[candidates], function(x) {
      vals <- suppressWarnings(as.numeric(as.character(x)))
      mean(!is.na(vals)) >= 0.7
    }, logical(1))]
  }
  if (!length(sample_columns)) stop("No numeric WGCNA sample columns were detected.", call. = FALSE)
  if (nrow(feature_table) != nrow(df_raw)) stop("Canonical WGCNA feature table lost expression rows.", call. = FALSE)
  keep <- feature_table$included_in_wgcna
  expression_rows <- as.data.frame(lapply(df_raw[keep, sample_columns, drop = FALSE], function(x) {
    suppressWarnings(as.numeric(as.character(x)))
  }), check.names = FALSE)
  expression_data <- as.data.frame(t(expression_rows), check.names = FALSE)
  colnames(expression_data) <- feature_table$ProteinGroupID[keep]
  if (anyNA(colnames(expression_data)) || any(!nzchar(colnames(expression_data)))) {
    stop("Missing ProteinGroupID in WGCNA expression columns; row-order repair is forbidden.", call. = FALSE)
  }
  if (anyDuplicated(colnames(expression_data))) {
    stop("Duplicate ProteinGroupID in WGCNA expression columns; make.unique() repair is forbidden.", call. = FALSE)
  }
  if (!identical(colnames(expression_data), feature_table$ProteinGroupID[keep])) {
    stop("WGCNA expression columns are not aligned to the canonical feature table.", call. = FALSE)
  }
  list(
    expression_data = expression_data,
    feature_table = feature_table,
    member_bridge = canonical$bridge,
    collision_audit = canonical$collision_audit,
    sample_columns = sample_columns
  )
}

tokenize_wgcna_mouse_only <- function(male_df, dropped_non_mouse_path = NULL) {
  tok <- male_df |>
    dplyr::mutate(gene_symbol = as.character(.data$gene_symbol)) |>
    tidyr::separate_rows("gene_symbol", sep = ";") |>
    dplyr::mutate(
      token_raw = .data$gene_symbol,
      token_up = normalize_token(.data$gene_symbol)
    )
  dropped_non_mouse <- tok |>
    dplyr::filter(is.na(.data$token_up) | !grepl("_MOUSE$", .data$token_up))
  if (!is.null(dropped_non_mouse_path) && nrow(dropped_non_mouse)) {
    if (exists("write_tsv_safe", mode = "function")) {
      write_tsv_safe(dropped_non_mouse, dropped_non_mouse_path)
    } else if (requireNamespace("readr", quietly = TRUE)) {
      readr::write_tsv(dropped_non_mouse, dropped_non_mouse_path)
    } else {
      utils::write.table(dropped_non_mouse, dropped_non_mouse_path, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  }
  tok |>
    dplyr::filter(!is.na(.data$token_up), grepl("_MOUSE$", .data$token_up)) |>
    dplyr::mutate(
      token_base = to_base_no_iso_mouse(.data$token_up),
      looks_ac = grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$|^A0A[0-9A-Z]{7}$", .data$token_base),
      looks_entry = grepl("^[A-Z0-9][A-Z0-9\\-\\.]+$", .data$token_base),
      id_class = dplyr::case_when(
        .data$looks_ac ~ "UNIPROT_AC_MOUSE",
        .data$looks_entry ~ "UNIPROT_ENTRY",
        TRUE ~ "UNKNOWN"
      ),
      Resolved_UNIPROT = NA_character_,
      strategy = NA_character_
    )
}

collapse_wgcna_ids <- function(x) {
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (!length(x)) return(NA_character_)
  paste(x, collapse = ";")
}

collapse_wgcna_bool <- function(x) any(as.logical(x), na.rm = TRUE)

build_wgcna_feature_exclusion_audit <- function(male_data, resolved) {
  row_tokens <- male_data |>
    dplyr::transmute(.row_id = .data$.row_id, original_input_token = as.character(.data$gene_symbol)) |>
    tidyr::separate_rows("original_input_token", sep = ";") |>
    dplyr::mutate(token_up = normalize_token(.data$original_input_token))

  has_mouse <- row_tokens |>
    dplyr::group_by(.data$.row_id) |>
    dplyr::summarise(has_mouse_identifier = any(!is.na(.data$token_up) & grepl("_MOUSE$", .data$token_up)), .groups = "drop")

  mapped_rows <- resolved |>
    dplyr::group_by(.data$.row_id) |>
    dplyr::summarise(has_resolved_uniprot = any(!is.na(.data$Resolved_UNIPROT) & nzchar(.data$Resolved_UNIPROT)), .groups = "drop")

  male_data |>
    dplyr::transmute(.row_id = .data$.row_id, original_input_token = as.character(.data$gene_symbol)) |>
    dplyr::left_join(has_mouse, by = ".row_id") |>
    dplyr::left_join(mapped_rows, by = ".row_id") |>
    dplyr::mutate(
      has_mouse_identifier = dplyr::coalesce(.data$has_mouse_identifier, FALSE),
      has_resolved_uniprot = dplyr::coalesce(.data$has_resolved_uniprot, FALSE),
      is_blank_identifier = is.na(.data$original_input_token) | !nzchar(trimws(.data$original_input_token)),
      exclusion_category = dplyr::case_when(
        .data$is_blank_identifier ~ "blank_identifier",
        !.data$has_mouse_identifier ~ "non_mouse_identifier",
        !.data$has_resolved_uniprot ~ "unresolved_mouse_identifier",
        TRUE ~ NA_character_
      ),
      exclusion_reason = dplyr::case_when(
        .data$exclusion_category == "blank_identifier" ~ "Original identifier is blank or missing.",
        .data$exclusion_category == "non_mouse_identifier" ~ "No semicolon-separated token ends with _MOUSE.",
        .data$exclusion_category == "unresolved_mouse_identifier" ~ "At least one _MOUSE token was present but none resolved to a UniProt accession.",
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::filter(!is.na(.data$exclusion_category)) |>
    dplyr::select(".row_id", "original_input_token", "exclusion_category", "exclusion_reason") |>
    dplyr::arrange(.data$.row_id)
}

build_wgcna_input_tables <- function(male_data, resolved) {
  feature_mapping_pre <- resolved |>
    dplyr::group_by(.data$.row_id) |>
    dplyr::summarise(
      original_token = collapse_wgcna_ids(.data$token_raw),
      resolved_uniprot = collapse_wgcna_ids(.data$Resolved_UNIPROT),
      mapping_strategy = collapse_wgcna_ids(.data$strategy),
      manual_mapping_used = collapse_wgcna_bool(.data$manual_mapping_used),
      .groups = "drop"
    ) |>
    dplyr::right_join(male_data |> dplyr::transmute(.row_id = .data$.row_id, original_input_token = .data$gene_symbol), by = ".row_id") |>
    dplyr::mutate(
      original_token = dplyr::coalesce(.data$original_token, as.character(.data$original_input_token), NA_character_),
      resolved_uniprot = dplyr::coalesce(.data$resolved_uniprot, NA_character_),
      mapping_strategy = dplyr::coalesce(.data$mapping_strategy, NA_character_),
      manual_mapping_used = dplyr::coalesce(.data$manual_mapping_used, FALSE),
      mapping_status = dplyr::if_else(!is.na(.data$resolved_uniprot) & nzchar(.data$resolved_uniprot), "mapped", "unmapped")
    ) |>
    dplyr::arrange(.data$.row_id)

  included_row_ids <- feature_mapping_pre |>
    dplyr::filter(.data$mapping_status == "mapped") |>
    dplyr::select(".row_id")

  exclusion_audit <- build_wgcna_feature_exclusion_audit(male_data, resolved)

  male_norm <- resolved |>
    dplyr::group_by(.data$.row_id) |>
    dplyr::summarise(gene_symbol = collapse_wgcna_ids(.data$Resolved_UNIPROT), .groups = "drop") |>
    dplyr::semi_join(included_row_ids, by = ".row_id") |>
    dplyr::left_join(male_data |> dplyr::select(-dplyr::all_of("gene_symbol")), by = ".row_id") |>
    dplyr::arrange(.data$.row_id) |>
    dplyr::select(-dplyr::all_of(".row_id")) |>
    dplyr::mutate(gene_symbol = dplyr::na_if(.data$gene_symbol, ""))

  feature_mapping_final <- feature_mapping_pre |>
    dplyr::semi_join(included_row_ids, by = ".row_id") |>
    dplyr::arrange(.data$.row_id)

  list(
    male_norm = male_norm,
    feature_mapping_pre = feature_mapping_pre,
    feature_mapping_final = feature_mapping_final,
    exclusion_audit = exclusion_audit
  )
}

validate_wgcna_expression_inputs <- function(expression_data, feature_mapping_tbl, exclusion_audit) {
  if (any(startsWith(colnames(expression_data), "UNMAPPED"))) {
    stop("Final WGCNA expression matrix contains UNMAPPED feature IDs.", call. = FALSE)
  }
  missing_accession <- is.na(feature_mapping_tbl$resolved_uniprot) | !nzchar(feature_mapping_tbl$resolved_uniprot)
  if (any(missing_accession)) {
    stop("Final WGCNA feature mapping contains rows without resolved UniProt accessions.", call. = FALSE)
  }
  if (ncol(expression_data) != nrow(feature_mapping_tbl)) {
    stop("Final WGCNA expression matrix and feature mapping audit have different feature counts.", call. = FALSE)
  }
  if (nrow(exclusion_audit)) {
    excluded_final <- intersect(exclusion_audit$.row_id, feature_mapping_tbl$.row_id)
    if (length(excluded_final)) {
      stop("Excluded blank, non-mouse, or unresolved rows appear in expression.data.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

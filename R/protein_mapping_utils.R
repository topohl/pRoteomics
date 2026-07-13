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
    dplyr::select(.data$input, .data$primaryAccession)

  gene_map <- dplyr::bind_rows(
    gene_map,
    entry_map |> dplyr::transmute(input = .data$entry_base, primaryAccession = .data$UNIPROT)
  ) |>
    dplyr::distinct(.data$input, .keep_all = TRUE)

  list(entry_map = entry_map, gene_map = gene_map)
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

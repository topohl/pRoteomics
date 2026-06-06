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

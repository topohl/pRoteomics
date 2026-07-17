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

mouse_gene_annotation_contract_version <- function() "mouse_gene_annotation_v1"

build_mouse_gene_annotation_maps <- function(org_db_obj, accessions = NULL) {
  symbol_keys <- as.character(AnnotationDbi::keys(org_db_obj, keytype = "SYMBOL"))
  symbol_map <- suppressMessages(AnnotationDbi::select(
    org_db_obj, keys = symbol_keys, keytype = "SYMBOL", columns = "ENTREZID"
  ))
  symbol_map <- unique(symbol_map[!is.na(symbol_map$SYMBOL) & nzchar(symbol_map$SYMBOL), c("SYMBOL", "ENTREZID"), drop = FALSE])

  alias_keys <- as.character(AnnotationDbi::keys(org_db_obj, keytype = "ALIAS"))
  alias_map <- suppressMessages(AnnotationDbi::select(
    org_db_obj, keys = alias_keys, keytype = "ALIAS", columns = c("SYMBOL", "ENTREZID")
  ))
  alias_map <- unique(alias_map[!is.na(alias_map$ALIAS) & nzchar(alias_map$ALIAS) &
    !is.na(alias_map$SYMBOL) & nzchar(alias_map$SYMBOL), c("ALIAS", "SYMBOL", "ENTREZID"), drop = FALSE])

  valid_accessions <- as.character(AnnotationDbi::keys(org_db_obj, keytype = "UNIPROT"))
  if (!is.null(accessions)) valid_accessions <- intersect(toupper(as.character(accessions)), valid_accessions)
  uniprot_map <- if (length(valid_accessions)) suppressMessages(AnnotationDbi::select(
    org_db_obj, keys = valid_accessions, keytype = "UNIPROT", columns = c("SYMBOL", "ENTREZID")
  )) else data.frame(UNIPROT = character(), SYMBOL = character(), ENTREZID = character(), stringsAsFactors = FALSE)
  uniprot_map <- unique(uniprot_map[!is.na(uniprot_map$UNIPROT) & nzchar(uniprot_map$UNIPROT) &
    !is.na(uniprot_map$SYMBOL) & nzchar(uniprot_map$SYMBOL), c("UNIPROT", "SYMBOL", "ENTREZID"), drop = FALSE])

  make_index <- function(df, key) {
    split(df, as.character(df[[key]]), drop = TRUE)
  }

  list(
    uniprot_map = uniprot_map,
    symbol_map = symbol_map,
    alias_map = alias_map,
    uniprot_index = make_index(uniprot_map, "UNIPROT"),
    symbol_index = make_index(symbol_map, "SYMBOL"),
    symbol_ci_index = split(symbol_map, tolower(as.character(symbol_map$SYMBOL)), drop = TRUE),
    alias_index = make_index(alias_map, "ALIAS"),
    orgdb_package_version = as.character(utils::packageVersion("org.Mm.eg.db")),
    annotation_contract_version = mouse_gene_annotation_contract_version()
  )
}

normalize_manual_gene_annotation_overrides <- function(overrides) {
  empty <- data.frame(lookup_value = character(), official_gene_symbol = character(), official_entrez_id = character(), stringsAsFactors = FALSE)
  if (is.null(overrides) || !is.data.frame(overrides) || !nrow(overrides)) return(empty)
  nms <- tolower(gsub("[^a-z0-9]+", "_", names(overrides)))
  names(overrides) <- nms
  lookup_col <- intersect(c("lookup_value", "member_accession", "uniprot", "submitted_gene_symbol", "gene_symbol"), nms)
  symbol_col <- intersect(c("official_gene_symbol", "resolved_gene_symbol", "symbol"), nms)
  entrez_col <- intersect(c("official_entrez_id", "entrezid", "entrez_id"), nms)
  if (!length(lookup_col) || !length(symbol_col)) return(empty)
  data.frame(
    lookup_value = toupper(trimws(as.character(overrides[[lookup_col[[1]]]]))),
    official_gene_symbol = trimws(as.character(overrides[[symbol_col[[1]]]])),
    official_entrez_id = if (length(entrez_col)) trimws(as.character(overrides[[entrez_col[[1]]]])) else NA_character_,
    stringsAsFactors = FALSE
  )
}

read_manual_gene_annotation_overrides <- function(path = Sys.getenv("PROTEOMICS_MANUAL_GENE_ANNOTATION_FILE", unset = "")) {
  if (!nzchar(path)) {
    path <- if (exists("path_metadata", mode = "function")) path_metadata("manual_gene_annotation_overrides.csv") else file.path("data", "metadata", "manual_gene_annotation_overrides.csv")
  }
  if (!file.exists(path)) {
    out <- normalize_manual_gene_annotation_overrides(NULL)
    attr(out, "path") <- path
    attr(out, "status") <- "missing"
    return(out)
  }
  ext <- tolower(tools::file_ext(path))
  raw <- if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("readxl is required for manual gene annotation workbooks.", call. = FALSE)
    readxl::read_excel(path, sheet = 1)
  } else if (requireNamespace("readr", quietly = TRUE)) {
    readr::read_csv(path, show_col_types = FALSE)
  } else {
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
  out <- normalize_manual_gene_annotation_overrides(as.data.frame(raw))
  attr(out, "path") <- path
  attr(out, "status") <- if (nrow(out)) "loaded" else "empty_or_invalid"
  out
}

resolve_mouse_gene_annotation <- function(member_accession, submitted_gene_symbol,
                                           annotation_maps, manual_overrides = NULL) {
  accession <- toupper(trimws(as.character(member_accession)))
  submitted <- trimws(as.character(submitted_gene_symbol))
  if (is.na(accession)) accession <- ""
  if (is.na(submitted)) submitted <- ""
  manual <- normalize_manual_gene_annotation_overrides(manual_overrides)
  manual_hits <- manual[toupper(manual$lookup_value) %in% unique(c(accession, toupper(submitted))), , drop = FALSE]
  manual_symbols <- sort(unique(manual_hits$official_gene_symbol[nzchar(manual_hits$official_gene_symbol)]), method = "radix")

  make_result <- function(symbol = NA_character_, entrez = NA_character_, status, strategy,
                          candidates = character(), secondary_conflict = FALSE, manual_used = FALSE) {
    data.frame(
      submitted_gene_symbol = if (nzchar(submitted)) submitted else NA_character_,
      official_gene_symbol = symbol,
      official_entrez_id = entrez,
      gene_annotation_status = status,
      gene_annotation_strategy = strategy,
      gene_annotation_candidates = paste(sort(unique(candidates), method = "radix"), collapse = ";"),
      gene_annotation_secondary_conflict = secondary_conflict,
      gene_annotation_manual_override_used = manual_used,
      gene_annotation_contract_version = if (is.null(annotation_maps$annotation_contract_version)) mouse_gene_annotation_contract_version() else annotation_maps$annotation_contract_version,
      orgdb_package_version = if (is.null(annotation_maps$orgdb_package_version)) NA_character_ else annotation_maps$orgdb_package_version,
      stringsAsFactors = FALSE
    )
  }

  use_manual <- function() {
    if (length(manual_symbols) != 1L) return(NULL)
    rows <- manual_hits[manual_hits$official_gene_symbol == manual_symbols[[1]], , drop = FALSE]
    entrez <- sort(unique(rows$official_entrez_id[!is.na(rows$official_entrez_id) & nzchar(rows$official_entrez_id)]), method = "radix")
    make_result(manual_symbols[[1]], if (length(entrez) == 1L) entrez[[1]] else NA_character_,
      "resolved", "manual_annotation_override", manual_symbols, manual_used = TRUE)
  }

  accession_rows <- if (!is.null(annotation_maps$uniprot_index)) {
    indexed <- annotation_maps$uniprot_index[[accession]]
    if (is.null(indexed)) annotation_maps$uniprot_map[0, , drop = FALSE] else indexed
  } else annotation_maps$uniprot_map[annotation_maps$uniprot_map$UNIPROT == accession, , drop = FALSE]
  accession_symbols <- sort(unique(accession_rows$SYMBOL[!is.na(accession_rows$SYMBOL) & nzchar(accession_rows$SYMBOL)]), method = "radix")
  if (length(accession_symbols) == 1L) {
    symbol_rows <- accession_rows[accession_rows$SYMBOL == accession_symbols[[1]], , drop = FALSE]
    entrez <- sort(unique(as.character(symbol_rows$ENTREZID[!is.na(symbol_rows$ENTREZID) & nzchar(as.character(symbol_rows$ENTREZID))])), method = "radix")
    alias_symbols <- sort(unique(annotation_maps$alias_map$SYMBOL[annotation_maps$alias_map$ALIAS == submitted]), method = "radix")
    secondary_conflict <- length(alias_symbols) > 0L && any(alias_symbols != accession_symbols[[1]])
    return(make_result(accession_symbols[[1]], if (length(entrez) == 1L) entrez[[1]] else NA_character_,
      "resolved", "unique_uniprot_to_symbol", accession_symbols, secondary_conflict))
  }
  if (length(accession_symbols) > 1L) {
    override <- use_manual()
    if (!is.null(override)) return(override)
    return(make_result(status = "ambiguous", strategy = "ambiguous_uniprot_mapping", candidates = accession_symbols))
  }

  exact_rows <- if (!is.null(annotation_maps$symbol_index)) {
    indexed <- annotation_maps$symbol_index[[submitted]]
    if (is.null(indexed)) annotation_maps$symbol_map[0, , drop = FALSE] else indexed
  } else annotation_maps$symbol_map[annotation_maps$symbol_map$SYMBOL == submitted, , drop = FALSE]
  if (nrow(exact_rows)) {
    entrez <- sort(unique(as.character(exact_rows$ENTREZID[!is.na(exact_rows$ENTREZID)])), method = "radix")
    return(make_result(submitted, if (length(entrez) == 1L) entrez[[1]] else NA_character_, "resolved", "exact_symbol_match", submitted))
  }

  ci_rows <- if (!is.null(annotation_maps$symbol_ci_index)) {
    indexed <- annotation_maps$symbol_ci_index[[tolower(submitted)]]
    if (is.null(indexed)) annotation_maps$symbol_map[0, , drop = FALSE] else indexed
  } else annotation_maps$symbol_map[tolower(annotation_maps$symbol_map$SYMBOL) == tolower(submitted), , drop = FALSE]
  ci_symbols <- sort(unique(ci_rows$SYMBOL), method = "radix")
  if (length(ci_symbols) == 1L) {
    rows <- annotation_maps$symbol_map[annotation_maps$symbol_map$SYMBOL == ci_symbols[[1]], , drop = FALSE]
    entrez <- sort(unique(as.character(rows$ENTREZID[!is.na(rows$ENTREZID)])), method = "radix")
    return(make_result(ci_symbols[[1]], if (length(entrez) == 1L) entrez[[1]] else NA_character_, "resolved", "unique_case_insensitive_symbol", ci_symbols))
  }
  if (length(ci_symbols) > 1L) return(make_result(status = "ambiguous", strategy = "ambiguous_case_insensitive_symbol", candidates = ci_symbols))

  alias_rows <- if (!is.null(annotation_maps$alias_index)) {
    indexed <- annotation_maps$alias_index[[submitted]]
    if (is.null(indexed)) annotation_maps$alias_map[0, , drop = FALSE] else indexed
  } else annotation_maps$alias_map[annotation_maps$alias_map$ALIAS == submitted, , drop = FALSE]
  alias_symbols <- sort(unique(alias_rows$SYMBOL[!is.na(alias_rows$SYMBOL) & nzchar(alias_rows$SYMBOL)]), method = "radix")
  if (length(alias_symbols) == 1L) {
    entrez <- sort(unique(as.character(alias_rows$ENTREZID[alias_rows$SYMBOL == alias_symbols[[1]] & !is.na(alias_rows$ENTREZID)])), method = "radix")
    return(make_result(alias_symbols[[1]], if (length(entrez) == 1L) entrez[[1]] else NA_character_, "resolved", "unique_alias_to_symbol", alias_symbols))
  }
  if (length(alias_symbols) > 1L) return(make_result(status = "ambiguous", strategy = "ambiguous_alias_mapping", candidates = alias_symbols))

  override <- use_manual()
  if (!is.null(override)) return(override)
  make_result(status = "unresolved", strategy = "unresolved_gene_annotation")
}

assess_protein_group_gene_annotation <- function(member_bridge) {
  relevant <- member_bridge[!is.na(member_bridge$member_accession) & nzchar(member_bridge$member_accession), , drop = FALSE]
  if (!nrow(relevant)) {
    return(data.frame(official_gene_symbol = NA_character_, official_entrez_id = NA_character_,
      protein_group_gene_annotation_status = "no_mapped_accessions", all_member_accessions_gene_annotated = FALSE, stringsAsFactors = FALSE))
  }
  resolved <- !is.na(relevant$member_gene_symbol) & nzchar(relevant$member_gene_symbol) & relevant$gene_annotation_status == "resolved"
  symbols <- sort(unique(relevant$member_gene_symbol[resolved]), method = "radix")
  entrez <- sort(unique(relevant$member_entrez_id[resolved & !is.na(relevant$member_entrez_id) & nzchar(relevant$member_entrez_id)]), method = "radix")
  all_resolved <- all(resolved)
  status <- if (!all_resolved) "incomplete_or_ambiguous_member_annotation" else if (length(symbols) != 1L || length(entrez) > 1L) "conflicting_member_annotations" else "concordant_official_gene"
  data.frame(
    official_gene_symbol = if (identical(status, "concordant_official_gene")) symbols[[1]] else NA_character_,
    official_entrez_id = if (identical(status, "concordant_official_gene") && length(entrez) == 1L) entrez[[1]] else NA_character_,
    protein_group_gene_annotation_status = status,
    all_member_accessions_gene_annotated = all_resolved,
    stringsAsFactors = FALSE
  )
}

apply_enrichment_gene_annotation_fallback <- function(df, annotation_maps,
                                                       uniprot_mapping_file_hash = NA_character_) {
  required <- c("ProteinGroupID", "member_accessions", "member_gene_symbols")
  missing <- setdiff(required, names(df))
  if (length(missing)) stop("Compatibility annotation fallback is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  bridge_rows <- list()
  group_rows <- vector("list", nrow(df))
  for (i in seq_len(nrow(df))) {
    accessions <- split_protein_group_members(df$member_accessions[[i]])
    submitted <- split_protein_group_members(df$member_gene_symbols[[i]])
    member_rows <- lapply(seq_along(accessions), function(j) {
      submitted_j <- if (length(submitted) == 1L) submitted[[1]] else if (length(submitted) == length(accessions)) submitted[[j]] else NA_character_
      annotation <- resolve_mouse_gene_annotation(accessions[[j]], submitted_j, annotation_maps)
      data.frame(ProteinGroupID = as.character(df$ProteinGroupID[[i]]), member_accession = toupper(accessions[[j]]),
        member_gene_symbol_submitted = submitted_j, member_gene_symbol = annotation$official_gene_symbol,
        member_entrez_id = annotation$official_entrez_id, gene_annotation_status = annotation$gene_annotation_status,
        gene_annotation_strategy = paste0("enrichment_compatibility_fallback|", annotation$gene_annotation_strategy),
        gene_annotation_candidates = annotation$gene_annotation_candidates,
        gene_annotation_secondary_conflict = annotation$gene_annotation_secondary_conflict,
        gene_annotation_manual_override_used = annotation$gene_annotation_manual_override_used,
        gene_annotation_contract_version = "enrichment_compatibility_fallback_v1",
        uniprot_mapping_file_hash = uniprot_mapping_file_hash,
        orgdb_package_version = annotation$orgdb_package_version, stringsAsFactors = FALSE)
    })
    bridge <- if (length(member_rows)) do.call(rbind, member_rows) else data.frame(
      member_accession = character(), member_gene_symbol = character(), member_entrez_id = character(),
      gene_annotation_status = character(), stringsAsFactors = FALSE)
    group <- assess_protein_group_gene_annotation(bridge)
    group_rows[[i]] <- cbind(data.frame(ProteinGroupID = as.character(df$ProteinGroupID[[i]]), stringsAsFactors = FALSE), group)
    bridge_rows[[i]] <- bridge
  }
  groups <- do.call(rbind, group_rows)
  df$official_gene_symbol <- groups$official_gene_symbol
  df$official_entrez_id <- groups$official_entrez_id
  df$protein_group_gene_annotation_status <- groups$protein_group_gene_annotation_status
  df$all_member_accessions_gene_annotated <- groups$all_member_accessions_gene_annotated
  df$gene_level_claim_allowed <- as.logical(df$gene_level_claim_allowed) &
    df$protein_group_gene_annotation_status == "concordant_official_gene"
  df$gene_annotation_contract_version <- "enrichment_compatibility_fallback_v1"
  df$uniprot_mapping_file_hash <- uniprot_mapping_file_hash
  df$orgdb_package_version <- annotation_maps$orgdb_package_version
  list(data = df, accession_audit = do.call(rbind, bridge_rows), protein_group_audit = groups)
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
    out <- data.frame(gene_symbol = character(), mapped_gene_symbol = character(), stringsAsFactors = FALSE)
    attr(out, "path") <- path
    attr(out, "status") <- "missing"
    return(out)
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
    out <- data.frame(gene_symbol = character(), mapped_gene_symbol = character(), stringsAsFactors = FALSE)
    attr(out, "path") <- path
    attr(out, "status") <- "empty_or_unreadable"
    return(out)
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

detect_source_provenance_columns <- function(df) {
  candidates <- c(
    "ProteinGroupID", "protein_group_id", "source_feature_id", "feature_id", "FeatureID",
    "Protein.Group", "Protein Group", "T: Protein.Group",
    "original_identifier", "gene_symbol", "T: Protein.Names", "Protein.Names",
    "Genes", "ProteinID", "UniProt", "row_id", "id"
  )
  intersect(candidates, names(df))
}

classify_source_provenance <- function(row, provenance_columns = detect_source_provenance_columns(row)) {
  provenance_columns <- intersect(provenance_columns, names(row))
  values <- if (length(provenance_columns)) {
    vapply(provenance_columns, function(nm) as.character(row[[nm]][1]), character(1))
  } else character()
  present <- !is.na(values) & nzchar(trimws(values))
  values <- values[present]
  columns <- provenance_columns[present]
  normalized <- toupper(trimws(values))

  # Stable source identifiers use several common contaminant/decoy prefixes.
  # This deliberately does not inspect free-text protein descriptions: a valid
  # mouse keratin must not become a contaminant merely because of its name.
  contaminant_hit <- grepl(
    "(^|[;|,/[:space:]])(?:CON__|CONT?_|CONTAMINANT(?:_|:)?|REV__|REVERSE_|DECOY_)",
    normalized, perl = TRUE
  )
  non_mouse_hit <- grepl("_(?:HUMAN|RAT)(?:$|[;|,/[:space:]])", normalized, perl = TRUE)

  evidence <- function(hit) {
    if (!any(hit)) return("")
    paste(paste0(columns[hit], "=", values[hit]), collapse = ";")
  }
  reasons <- c(
    if (any(contaminant_hit)) "contaminant_source_provenance",
    if (any(non_mouse_hit)) "non_mouse_source_provenance"
  )
  list(
    source_provenance_columns = paste(columns, collapse = ";"),
    source_provenance_values = paste(values, collapse = ";"),
    contaminant_source_evidence = evidence(contaminant_hit),
    non_mouse_source_evidence = evidence(non_mouse_hit),
    source_provenance_exclusion_reason = paste(reasons, collapse = ";"),
    source_provenance_contaminant = any(contaminant_hit),
    source_provenance_non_mouse = any(non_mouse_hit)
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
  contaminant <- grepl(
    "(^|[|;,:/])(?:CON__|CONT?_|CONTAMINANT(?:_|:)?|REV__|REVERSE_|DECOY_)",
    token_norm, perl = TRUE
  )
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
                                         manual_mapping = NULL, manual_override = TRUE,
                                         gene_annotation_maps = NULL,
                                         manual_gene_annotation_overrides = NULL,
                                         uniprot_mapping_file_hash = NA_character_) {
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

  submitted_gene <- if (!is.na(acc) && nzchar(acc)) lookup_accession_gene(acc, accession_gene_map, gene_map) else NA_character_
  if (!is.null(gene_annotation_maps)) {
    annotation <- resolve_mouse_gene_annotation(acc, submitted_gene, gene_annotation_maps, manual_gene_annotation_overrides)
    gene <- annotation$official_gene_symbol[[1]]
    entrez <- annotation$official_entrez_id[[1]]
  } else {
    gene <- submitted_gene
    entrez <- NA_character_
    annotation <- data.frame(
      gene_annotation_status = ifelse(!is.na(gene) && nzchar(gene), "resolved", "unresolved"),
      gene_annotation_strategy = "legacy_uniprot_gene_name",
      gene_annotation_candidates = ifelse(!is.na(gene), gene, ""),
      gene_annotation_secondary_conflict = FALSE,
      gene_annotation_manual_override_used = FALSE,
      gene_annotation_contract_version = "legacy_uniprot_gene_name",
      orgdb_package_version = NA_character_, stringsAsFactors = FALSE
    )
  }
  reviewed_status <- NA_character_
  if (!is.null(reviewed_map) && nrow(reviewed_map) && !is.na(acc) && nzchar(acc)) {
    reviewed_status <- reviewed_map$reviewed_status[match(acc, reviewed_map$UNIPROT)]
  }

  data.frame(
    member_identifier_original = as.character(member_identifier),
    member_identifier_normalized = parsed$member_identifier_normalized,
    member_accession = ifelse(!is.na(acc) && nzchar(acc), toupper(acc), NA_character_),
    member_gene_symbol_submitted = submitted_gene,
    member_gene_symbol = gene,
    member_entrez_id = entrez,
    gene_annotation_status = annotation$gene_annotation_status[[1]],
    gene_annotation_strategy = annotation$gene_annotation_strategy[[1]],
    gene_annotation_candidates = annotation$gene_annotation_candidates[[1]],
    gene_annotation_secondary_conflict = annotation$gene_annotation_secondary_conflict[[1]],
    gene_annotation_manual_override_used = annotation$gene_annotation_manual_override_used[[1]],
    gene_annotation_contract_version = annotation$gene_annotation_contract_version[[1]],
    orgdb_package_version = annotation$orgdb_package_version[[1]],
    uniprot_mapping_file_hash = uniprot_mapping_file_hash,
    member_species = parsed$member_species,
    member_mapping_status = ifelse(!is.na(acc) && nzchar(acc), "mapped", "unmapped"),
    mapping_strategy = ifelse(!is.na(strategy) && nzchar(strategy), strategy, "unmapped"),
    manual_mapping_used = manual_used,
    reviewed_status = reviewed_status,
    contaminant_status = parsed$contaminant_status,
    stringsAsFactors = FALSE
  )
}

classify_protein_group <- function(member_bridge, explicit_master_present = FALSE,
                                   source_provenance = NULL) {
  mapped <- member_bridge[member_bridge$member_mapping_status == "mapped", , drop = FALSE]
  n_mapped <- length(unique(stats::na.omit(mapped$member_accession)))
  genes <- unique(stats::na.omit(mapped$member_gene_symbol))
  contaminant <- any(member_bridge$contaminant_status == "contaminant", na.rm = TRUE) ||
    isTRUE(source_provenance$source_provenance_contaminant)
  non_mouse <- any(!is.na(member_bridge$member_species) & member_bridge$member_species != "mouse") ||
    isTRUE(source_provenance$source_provenance_non_mouse)
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
                                                 gene_annotation_maps = NULL,
                                                 manual_gene_annotation_overrides = NULL,
                                                 uniprot_mapping_file_hash = NA_character_,
                                                 strict = TRUE,
                                                 identifier_col = NULL,
                                                 feature_col = NULL) {
  source_cols <- detect_source_feature_columns(df_raw)
  provenance_cols <- detect_source_provenance_columns(df_raw)
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
    source_feature_id <- if (!is.na(feature_col) && feature_col %in% names(row)) as.character(row[[feature_col]][1]) else NA_character_
    member_identifier_source <- identifier_col
    member_identifier_value <- original_identifier
    if ((!length(split_protein_group_members(member_identifier_value))) &&
        !is.na(feature_col) && feature_col %in% c("Protein.Group", "Protein Group", "T: Protein.Group") &&
        !is.na(source_feature_id) && nzchar(source_feature_id)) {
      member_identifier_source <- feature_col
      member_identifier_value <- source_feature_id
    }
    members_original <- split_protein_group_members(member_identifier_value)
    members_canonical <- canonical_member_set(members_original)
    source_provenance <- classify_source_provenance(row, provenance_cols)

    member_tbls <- lapply(members_original, resolve_protein_group_member,
                          entry_map = entry_map, gene_map = gene_map,
                          accession_gene_map = accession_gene_map,
                          reviewed_map = reviewed_map,
                          manual_mapping = manual_mapping,
                          manual_override = manual_override,
                          gene_annotation_maps = gene_annotation_maps,
                          manual_gene_annotation_overrides = manual_gene_annotation_overrides,
                          uniprot_mapping_file_hash = uniprot_mapping_file_hash)
    bridge <- if (length(member_tbls)) dplyr::bind_rows(member_tbls) else data.frame()
    if (!nrow(bridge)) {
      bridge <- data.frame(
        member_identifier_original = NA_character_,
        member_identifier_normalized = NA_character_,
        member_accession = NA_character_,
        member_gene_symbol_submitted = NA_character_,
        member_gene_symbol = NA_character_,
        member_entrez_id = NA_character_,
        gene_annotation_status = "unresolved",
        gene_annotation_strategy = "unresolved_gene_annotation",
        gene_annotation_candidates = "",
        gene_annotation_secondary_conflict = FALSE,
        gene_annotation_manual_override_used = FALSE,
        gene_annotation_contract_version = if (is.null(gene_annotation_maps)) "legacy_uniprot_gene_name" else gene_annotation_maps$annotation_contract_version,
        orgdb_package_version = if (is.null(gene_annotation_maps)) NA_character_ else gene_annotation_maps$orgdb_package_version,
        uniprot_mapping_file_hash = uniprot_mapping_file_hash,
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
    group_annotation <- if (is.null(gene_annotation_maps)) {
      genes <- unique(stats::na.omit(bridge$member_gene_symbol))
      data.frame(official_gene_symbol = if (length(genes) == 1L) genes[[1]] else NA_character_, official_entrez_id = NA_character_,
        protein_group_gene_annotation_status = if (length(genes) == 1L) "legacy_uniprot_gene_name" else "legacy_ambiguous_gene_name",
        all_member_accessions_gene_annotated = length(genes) == 1L, stringsAsFactors = FALSE)
    } else assess_protein_group_gene_annotation(bridge)
    ambiguity_class <- classify_protein_group(bridge, explicit_master_present, source_provenance)
    claims <- protein_group_claim_rules(ambiguity_class)
    if (!is.null(gene_annotation_maps)) {
      claims$gene_level_claim_allowed <- claims$gene_level_claim_allowed &&
        identical(group_annotation$protein_group_gene_annotation_status[[1]], "concordant_official_gene")
    }
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
    bridge$member_identifier_source <- member_identifier_source
    bridge$source_provenance_columns <- source_provenance$source_provenance_columns
    bridge$source_provenance_values <- source_provenance$source_provenance_values
    bridge$contaminant_source_evidence <- source_provenance$contaminant_source_evidence
    bridge$non_mouse_source_evidence <- source_provenance$non_mouse_source_evidence
    bridge$source_provenance_exclusion_reason <- source_provenance$source_provenance_exclusion_reason

    mapped_accessions <- unique(stats::na.omit(bridge$member_accession))
    submitted_genes <- unique(stats::na.omit(bridge$member_gene_symbol_submitted))
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
        member_identifier_source = member_identifier_source,
        member_identifiers_original = paste(members_original, collapse = ";"),
        member_identifiers_canonical = paste(members_canonical, collapse = ";"),
        member_accessions = paste(mapped_accessions, collapse = ";"),
        member_gene_symbols_submitted = paste(submitted_genes, collapse = ";"),
        member_gene_symbols = paste(mapped_genes, collapse = ";"),
        official_gene_symbol = group_annotation$official_gene_symbol[[1]],
        official_entrez_id = group_annotation$official_entrez_id[[1]],
        protein_group_gene_annotation_status = group_annotation$protein_group_gene_annotation_status[[1]],
        all_member_accessions_gene_annotated = group_annotation$all_member_accessions_gene_annotated[[1]],
        gene_annotation_contract_version = if (is.null(gene_annotation_maps)) "legacy_uniprot_gene_name" else gene_annotation_maps$annotation_contract_version,
        uniprot_mapping_file_hash = uniprot_mapping_file_hash,
        orgdb_package_version = if (is.null(gene_annotation_maps)) NA_character_ else gene_annotation_maps$orgdb_package_version,
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
        source_provenance_columns = source_provenance$source_provenance_columns,
        source_provenance_values = source_provenance$source_provenance_values,
        contaminant_source_evidence = source_provenance$contaminant_source_evidence,
        non_mouse_source_evidence = source_provenance$non_mouse_source_evidence,
        source_provenance_exclusion_reason = source_provenance$source_provenance_exclusion_reason,
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
      "original_identifier", "member_identifier_source",
      "source_provenance_columns", "source_provenance_values",
      "contaminant_source_evidence", "non_mouse_source_evidence",
      "source_provenance_exclusion_reason",
      "member_identifier_original", "member_identifier_normalized",
      "member_rank_original", "member_rank_canonical",
      "member_accession", "member_gene_symbol_submitted", "member_gene_symbol", "member_entrez_id", "member_species",
      "member_mapping_status", "mapping_strategy", "manual_mapping_used",
      "gene_annotation_status", "gene_annotation_strategy", "gene_annotation_candidates",
      "gene_annotation_secondary_conflict", "gene_annotation_manual_override_used",
      "gene_annotation_contract_version", "uniprot_mapping_file_hash", "orgdb_package_version",
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

  protein_group_annotation_audit <- wide |>
    dplyr::select(dplyr::all_of(c("ProteinGroupID", "source_file", "source_feature_id", "source_row_id",
      "member_accessions", "member_gene_symbols", "official_gene_symbol", "official_entrez_id",
      "protein_group_gene_annotation_status", "all_member_accessions_gene_annotated",
      "gene_level_claim_allowed", "gene_annotation_contract_version",
      "uniprot_mapping_file_hash", "orgdb_package_version")))
  accession_annotation_audit <- bridge |>
    dplyr::select(dplyr::all_of(c("ProteinGroupID", "source_file", "source_row_id",
      "member_identifier_original", "member_accession", "member_gene_symbol_submitted",
      "member_gene_symbol", "member_entrez_id", "gene_annotation_status",
      "gene_annotation_strategy", "gene_annotation_candidates",
      "gene_annotation_secondary_conflict", "gene_annotation_manual_override_used",
      "gene_annotation_contract_version", "uniprot_mapping_file_hash", "orgdb_package_version")))

  list(wide = wide, bridge = bridge, summary = summary, collision_audit = collisions,
    accession_annotation_audit = accession_annotation_audit,
    protein_group_annotation_audit = protein_group_annotation_audit,
    deprecated_dropped_log = deprecated_dropped_log)
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
                                           gene_annotation_maps = NULL,
                                           manual_gene_annotation_overrides = NULL,
                                           uniprot_mapping_file_hash = NA_character_,
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
    gene_annotation_maps = gene_annotation_maps,
    manual_gene_annotation_overrides = manual_gene_annotation_overrides,
    uniprot_mapping_file_hash = uniprot_mapping_file_hash,
    strict = strict,
    identifier_col = identifier_col,
    feature_col = feature_col
  )
  feature_table <- canonical$wide
  feature_table$FeatureDisplayLabel <- wgcna_feature_display_label(feature_table)
  feature_table$quantitative_eligibility <- !feature_table$protein_group_ambiguity_class %in%
    "mixed_species_or_contaminant"
  feature_table$wgcna_exclusion_reason <- dplyr::case_when(
    feature_table$quantitative_eligibility ~ NA_character_,
    nzchar(feature_table$source_provenance_exclusion_reason) ~ feature_table$source_provenance_exclusion_reason,
    TRUE ~ "mixed_species_or_contaminant_member_evidence"
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

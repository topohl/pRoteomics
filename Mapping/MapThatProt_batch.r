#' Batch UniProt ID Mapping for ClusterProfiler Analysis (parallelized with doParallel)
#' Processes all .csv files in Datasets/raw and maps UniProtKB IDs to Accessions.

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, stringr, tidyr, purrr, readr, R.utils, foreach, doParallel)

mapped_comparisons <- "unknown-comparison"  # specify the comparison folder to process
map_reverse <- FALSE

#working_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/"
working_dir <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler"
comparison_dir <- if (isTRUE(map_reverse)) file.path(mapped_comparisons, "reverse") else mapped_comparisons
raw_dir <- file.path(working_dir, "Datasets", "raw", comparison_dir)
mapped_dir <- file.path(working_dir, "Datasets", "mapped", comparison_dir)
unmapped_dir <- file.path(working_dir, "Datasets", "unmapped", comparison_dir)

dir.create(mapped_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(unmapped_dir, recursive = TRUE, showWarnings = FALSE)

setwd(working_dir)

uniprot_mapping_file_path <- file.path(working_dir, "Datasets", "MOUSE_10090_idmapping.dat")
#protein2ipr_mapping_file_path <- file.path(working_dir, "Datasets", "protein2ipr.dat")

# ensure mapping file is present (download if necessary)
if (!file.exists(uniprot_mapping_file_path)) {
    cat("UniProt mapping file not found at:", uniprot_mapping_file_path, "\nAttempting to download...\n")
    options(timeout = 3600)
    gz_url <- "https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/idmapping/by_organism/MOUSE_10090_idmapping.dat.gz"
    gz_file <- paste0(uniprot_mapping_file_path, ".gz")
    tryCatch({
        download.file(gz_url, gz_file, mode = "wb")
        R.utils::gunzip(gz_file, destname = uniprot_mapping_file_path, remove = TRUE)
        cat("Downloaded and unzipped the UniProt mapping file successfully.\n")
    }, error = function(e) stop("Failed to download/unzip UniProt mapping file: ", e$message))
}

# load mapping once
cat("Loading UniProt mapping file from:", uniprot_mapping_file_path, "\n")
uniprot_mapping <- readr::read_tsv(
    uniprot_mapping_file_path,
    col_names = c("UniProt_Accession", "Type", "Value"),
    col_types = "ccc",
    quote = ""
)
entry_name_to_accession <- uniprot_mapping %>%
    filter(Type == "UniProtKB-ID") %>%
    select(UniProt_Accession, UniProtKB_ID = Value) %>%
    distinct(UniProtKB_ID, .keep_all = TRUE)

if (nrow(entry_name_to_accession) == 0) stop("No UniProtKB-ID mappings found in mapping file.")

# load protein to IPR/domain mapping (6 columns as in sample file)
#protein2ipr_mapping <- readr::read_tsv(
#    protein2ipr_mapping_file_path,
#    col_names = c("UniProt_Accession","IPR_ID","IPR_Name","Source_DB_ID","Start","End"),
#    col_types = readr::cols(
#        UniProt_Accession = readr::col_character(),
#        IPR_ID            = readr::col_character(),
#        IPR_Name          = readr::col_character(),
#        Source_DB_ID      = readr::col_character(),
#        Start             = readr::col_integer(),
#        End               = readr::col_integer()
#    ),
#    quote = ""
#) %>%
#    dplyr::mutate(
#        UniProt_Accession = toupper(trimws(UniProt_Accession)),
#        IPR_ID            = toupper(trimws(IPR_ID)),
#        Source_DB_ID      = toupper(trimws(Source_DB_ID))
#    ) %>%
#    dplyr::filter(
#        nzchar(UniProt_Accession),
#        nzchar(IPR_ID)
#    )

# find CSV files
csv_files <- list.files(raw_dir, pattern = ".*_.*\\.csv$", full.names = TRUE)
if (length(csv_files) == 0) stop("No .csv files found in: ", raw_dir)

# Optional manual mapping file (user-provided), columns can include:
# original_symbol OR symbol OR input OR gene_symbol
# final_accession OR accession OR uniprot OR uniprot_accession
# Manual mapping is provided as an Excel file with columns:
# "gene_symbol" (unmapped) and "mapped_gene_symbol" (mapped)
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")

manual_mapping_path <- file.path(working_dir, "Datasets", "manual_mapping.xlsx")
manual_override <- TRUE  # if TRUE, manual map can overwrite previously resolved entries

read_manual_xlsx <- function(path) {
    mm <- try(readxl::read_excel(path, sheet = 1), silent = TRUE)
    if (inherits(mm, "try-error") || !is.data.frame(mm) || !nrow(mm)) {
        message("Manual mapping Excel unreadable or empty: ", path)
        return(NULL)
    }
    mm %>%
        dplyr::mutate(dplyr::across(dplyr::everything(), ~ toupper(trimws(as.character(.))))) %>%
        dplyr::rename_with(tolower)
}

manual_mapping <- if (file.exists(manual_mapping_path)) {
    read_manual_xlsx(manual_mapping_path)
} else {
    message("Manual mapping Excel not found at: ", manual_mapping_path)
    NULL
}

if (!is.null(manual_mapping)) {
    # Ensure expected columns exist for symbol-to-symbol mapping
    needed <- c("gene_symbol", "mapped_gene_symbol")
    missing_cols <- setdiff(needed, names(manual_mapping))
    if (length(missing_cols)) {
        message("Manual mapping is missing columns: ", paste(missing_cols, collapse = ", "))
    } else {
        message("Loaded manual mapping rows: ", nrow(manual_mapping))
    }
}

process_file <- function(data_path) {
    message("Processing file: ", data_path)

    df_raw <- tryCatch(
        readr::read_csv(data_path, col_names = TRUE, show_col_types = FALSE, trim_ws = TRUE, quote = "\""),
        error = function(e) {
            message("read_csv failed, trying read_csv2(): ", e$message)
            readr::read_csv2(data_path, col_names = TRUE, show_col_types = FALSE, trim_ws = TRUE)
        }
    )

    if (!"gene_symbol" %in% names(df_raw)) {
        names(df_raw)[1] <- "gene_symbol"
        message("Renamed first column to 'gene_symbol' for file: ", data_path)
    }

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
        x <- gsub("-\\d+$", "", x)
        gsub("_MOUSE$", "", x)
    }
    is_uniprot_ac <- function(x) {
        grepl("^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9]|A0A[0-9A-Z]{7})$", x)
    }
    extract_ac <- function(s) {
        s <- as.character(s)
        m <- stringr::str_match(s, "(?i)(?:^|\\|)([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9]|A0A[0-9A-Z]{7})(?:\\-|\\||$|[^A-Z0-9])")
        out <- ifelse(is.na(m[,2]), NA_character_, toupper(m[,2]))
        gsub("-\\d+$", "", out)
    }
    extract_entry <- function(s) {
        s <- as.character(s)
        m <- stringr::str_match(s, "(?i)(?:^|\\|)([A-Z0-9]+_MOUSE)(?:\\||$|\\s)")
        out <- ifelse(is.na(m[,2]), NA_character_, m[,2])
        if (all(is.na(out))) {
            m2 <- stringr::str_match(s, "(?i)\\b([A-Z0-9]+_MOUSE)\\b")
            out <- m2[,2]
        }
        toupper(out)
    }
    nz <- function(x) !is.na(x) & nzchar(x)

    df_tok <- df_raw %>%
        tidyr::separate_rows(gene_symbol, sep = ";") %>%
        dplyr::mutate(
            token_raw = trimws(as.character(gene_symbol)),
            token_up  = normalize_token(token_raw),
            acc_guess = extract_ac(token_raw),
            entry_guess_up = extract_entry(token_raw)
        ) %>%
        dplyr::filter(nzchar(token_up)) %>%
        dplyr::mutate(
            token_base = dplyr::case_when(
                nz(entry_guess_up) ~ to_base_no_iso_mouse(entry_guess_up),
                TRUE               ~ to_base_no_iso_mouse(token_up)
            ),
            token_kind = dplyr::case_when(
                nz(acc_guess)                          ~ "AC_GUESS",
                nz(entry_guess_up)                     ~ "ENTRY_GUESS",
                grepl("_MOUSE$", token_up)             ~ "ENTRY",
                is_uniprot_ac(token_base)              ~ "AC",
                TRUE                                   ~ "SYMBOL_OR_ALIAS"
            )
        ) %>%
        dplyr::distinct(token_up, .keep_all = TRUE) %>%
        dplyr::filter(grepl("_MOUSE$", token_up))  # enforce _MOUSE filter

    if (nrow(df_tok) == 0) {
        message("No _MOUSE entries remaining after exclusion in file: ", data_path)
        return(invisible(list(unmapped_table = tibble::tibble(gene_symbol = character()))))
    }

    entry_map <- uniprot_mapping %>%
        dplyr::filter(Type == "UniProtKB-ID") %>%
        dplyr::transmute(
            UNIPROT    = toupper(trimws(UniProt_Accession)),
            entry_full = toupper(trimws(Value)),
            entry_base = toupper(gsub("_MOUSE$", "", trimws(Value)))
        ) %>%
        dplyr::filter(grepl("_MOUSE\\s*$", entry_full), nzchar(UNIPROT)) %>%
        dplyr::distinct(entry_base, .keep_all = TRUE)

    gene_map <- uniprot_mapping %>%
        dplyr::filter(Type %in% c("Gene_Name", "Gene_Name(synonym)", "Gene_Synonym")) %>%
        dplyr::transmute(
            primaryAccession = toupper(trimws(UniProt_Accession)),
            input            = toupper(trimws(Value))
        ) %>%
        dplyr::filter(nzchar(input), nzchar(primaryAccession)) %>%
        dplyr::mutate(pref = !startsWith(primaryAccession, "A0A")) %>%
        dplyr::arrange(dplyr::desc(pref), primaryAccession, input) %>%
        dplyr::group_by(input) %>% dplyr::slice_head(n = 1) %>% dplyr::ungroup() %>%
        dplyr::select(input, primaryAccession)

    gene_map <- dplyr::bind_rows(
        gene_map,
        entry_map %>% dplyr::transmute(input = entry_base, primaryAccession = UNIPROT)
    ) %>% dplyr::distinct(input, .keep_all = TRUE)

    if (!nrow(entry_map)) stop("entry_map is empty after robust parse")

    resolved <- df_tok %>%
        dplyr::transmute(
            token_raw, token_up, token_base, token_kind,
            acc_guess = toupper(acc_guess),
            entry_guess_up = toupper(entry_guess_up),
            id_class = token_kind,
            Resolved_UNIPROT = NA_character_,
            strategy = NA_character_
        )

    idx0 <- which(nz(resolved$acc_guess))
    if (length(idx0)) {
        resolved$Resolved_UNIPROT[idx0] <- resolved$acc_guess[idx0]
        resolved$strategy[idx0] <- "accept_accession_in_token"
    }

    need <- which(is.na(resolved$Resolved_UNIPROT) | !nzchar(resolved$Resolved_UNIPROT))
    if (length(need)) {
        idx_ac <- need[is_uniprot_ac(resolved$token_base[need])]
        if (length(idx_ac)) {
            resolved$Resolved_UNIPROT[idx_ac] <- resolved$token_base[idx_ac]
            resolved$strategy[idx_ac] <- "accept_accession_base"
        }
    }

    need <- which(is.na(resolved$Resolved_UNIPROT) | !nzchar(resolved$Resolved_UNIPROT))
    idx_en_guess <- need[nz(resolved$entry_guess_up[need])]
    if (length(idx_en_guess)) {
        key <- to_base_no_iso_mouse(resolved$entry_guess_up[idx_en_guess])
        hit <- entry_map$UNIPROT[match(key, entry_map$entry_base)]
        ok <- nz(hit)
        if (any(ok)) {
            ii <- idx_en_guess[ok]
            resolved$Resolved_UNIPROT[ii] <- hit[ok]
            resolved$strategy[ii] <- "entry_from_token"
        }
    }

    need <- which(is.na(resolved$Resolved_UNIPROT) | !nzchar(resolved$Resolved_UNIPROT))
    idx_en <- need[resolved$id_class[need] %in% c("ENTRY")]
    if (length(idx_en)) {
        hit <- entry_map$UNIPROT[match(toupper(resolved$token_base[idx_en]), entry_map$entry_base)]
        ok <- nz(hit)
        if (any(ok)) {
            ii <- idx_en[ok]
            resolved$Resolved_UNIPROT[ii] <- hit[ok]
            resolved$strategy[ii] <- "entry_local_mouse"
        }
    }

    need <- which(is.na(resolved$Resolved_UNIPROT) | !nzchar(resolved$Resolved_UNIPROT))
    idx_sym <- need[resolved$id_class[need] == "SYMBOL_OR_ALIAS"]
    if (length(idx_sym) && nrow(gene_map)) {
        base_need <- toupper(resolved$token_base[idx_sym])
        hit <- gene_map$primaryAccession[match(base_need, gene_map$input)]
        ok <- nz(hit)
        if (any(ok)) {
            ii <- idx_sym[ok]
            resolved$Resolved_UNIPROT[ii] <- hit[ok]
            resolved$strategy[ii] <- "gene_local_mouse"
        }
    }

    need <- which(is.na(resolved$Resolved_UNIPROT) | !nzchar(resolved$Resolved_UNIPROT))
    idx_entry_to_gene <- need[resolved$id_class[need] %in% c("ENTRY","ENTRY_GUESS")]
    if (length(idx_entry_to_gene) && nrow(gene_map)) {
        entry_basis <- ifelse(nz(resolved$entry_guess_up[idx_entry_to_gene]),
                                                    to_base_no_iso_mouse(resolved$entry_guess_up[idx_entry_to_gene]),
                                                    resolved$token_base[idx_entry_to_gene])
        gene_keys <- toupper(entry_basis)
        hit <- gene_map$primaryAccession[match(gene_keys, gene_map$input)]
        ok <- nz(hit)
        if (any(ok)) {
            ii <- idx_entry_to_gene[ok]
            resolved$Resolved_UNIPROT[ii] <- hit[ok]
            resolved$strategy[ii] <- "entry_prefix_as_gene_local"
        }
    }

    need <- which(is.na(resolved$Resolved_UNIPROT) | !nzchar(resolved$Resolved_UNIPROT))
    idx_family <- need[resolved$id_class[need] %in% c("SYMBOL_OR_ALIAS","ENTRY","ENTRY_GUESS")]
    if (length(idx_family) && nrow(entry_map)) {
        bases <- unique(toupper(resolved$token_base[idx_family]))
        pick_idx <- lapply(bases, function(b) {
            ix <- which(grepl(paste0("^", b, "[A-Z0-9]+$"), entry_map$entry_base))
            if (!length(ix)) return(NA_integer_)
            em <- entry_map[ix, , drop = FALSE]
            em$suffix <- sub(paste0("^", b), "", em$entry_base)
            em$score_review <- ifelse(startsWith(em$UNIPROT, "A0A"), 1L, 0L)
            em$score_suffixA <- ifelse(grepl("^A", em$suffix), 0L, 1L)
            ord <- order(em$score_review, em$score_suffixA, em$suffix, em$UNIPROT)
            ix[ord[1]]
        })
        sel <- !is.na(pick_idx)
        if (any(sel)) {
            b_ok <- bases[sel]
            acc_ok <- entry_map$UNIPROT[unlist(pick_idx[sel])]
            map_vec <- stats::setNames(acc_ok, b_ok)
            key_now <- toupper(resolved$token_base[idx_family])
            hit <- unname(map_vec[key_now])
            ok <- nz(hit)
            ii <- idx_family[ok]
            if (length(ii)) {
                resolved$Resolved_UNIPROT[ii] <- hit[ok]
                resolved$strategy[ii] <- "entry_family_prefix_pick"
            }
        }
    }

    # OrgDb strategies
    if (requireNamespace("AnnotationDbi", quietly = TRUE) && requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
        need_idx <- which((is.na(resolved$Resolved_UNIPROT) | !nzchar(resolved$Resolved_UNIPROT)) & resolved$id_class == "SYMBOL_OR_ALIAS")
        need_ids <- toupper(unique(resolved$token_base[need_idx]))
        ids_ent <- unique(need_ids[!is_uniprot_ac(need_ids)])
        if (length(ids_ent)) {
            sel_sym <- try(AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys = ids_ent, keytype = "SYMBOL", columns = c("MGIID","ENTREZID","UNIPROT","SYMBOL")), silent = TRUE)
            map_sym <- tibble::tibble(input = character(), primaryAccession = character())
            if (!inherits(sel_sym, "try-error") && nrow(sel_sym)) {
                map_sym <- tibble::as_tibble(sel_sym) %>%
                    dplyr::filter(!is.na(UNIPROT) & nzchar(UNIPROT)) %>%
                    dplyr::group_by(SYMBOL) %>% dplyr::arrange(UNIPROT, .by_group = TRUE) %>% dplyr::slice_head(n = 1) %>% dplyr::ungroup() %>%
                    dplyr::transmute(input = toupper(SYMBOL), primaryAccession = toupper(UNIPROT))
            }
            kt <- try(AnnotationDbi::keytypes(org.Mm.eg.db::org.Mm.eg.db), silent = TRUE)
            map_alias <- tibble::tibble(input = character(), primaryAccession = character())
            if (!inherits(kt, "try-error") && "ALIAS" %in% kt) {
                sel_alias <- try(AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys = ids_ent, keytype = "ALIAS", columns = c("UNIPROT","ALIAS")), silent = TRUE)
                if (!inherits(sel_alias, "try-error") && nrow(sel_alias)) {
                    map_alias <- tibble::as_tibble(sel_alias) %>%
                        dplyr::filter(!is.na(UNIPROT) & nzchar(UNIPROT)) %>%
                        dplyr::group_by(ALIAS) %>% dplyr::arrange(UNIPROT, .by_group = TRUE) %>% dplyr::slice_head(n = 1) %>% dplyr::ungroup() %>%
                        dplyr::transmute(input = toupper(ALIAS), primaryAccession = toupper(UNIPROT))
                }
            }
            map_symall <- dplyr::bind_rows(map_sym, map_alias) %>% dplyr::distinct(input, .keep_all = TRUE)
            if (nrow(map_symall)) {
                need_idx_now <- which((is.na(resolved$Resolved_UNIPROT) | !nzchar(resolved$Resolved_UNIPROT)) & resolved$id_class == "SYMBOL_OR_ALIAS")
                base_need <- toupper(resolved$token_base[need_idx_now])
                hit <- map_symall$primaryAccession[match(base_need, map_symall$input)]
                ok <- !is.na(hit) & nzchar(hit)
                ii <- need_idx_now[ok]
                if (length(ii)) {
                    resolved$Resolved_UNIPROT[ii] <- hit[ok]
                    resolved$strategy[ii] <- "orgdb_mgi_symbol_first"
                }
            }
        }

        need_idx <- which((is.na(resolved$Resolved_UNIPROT) | !nzchar(resolved$Resolved_UNIPROT)) & resolved$id_class == "SYMBOL_OR_ALIAS")
        need_ids <- toupper(unique(resolved$token_base[need_idx]))
        sym_left <- unique(need_ids[grepl("^[A-Z0-9\\-]{2,}$", need_ids)])
        if (length(sym_left)) {
            sym2eg <- try(AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys = sym_left, keytype = "SYMBOL", columns = c("ENTREZID","SYMBOL")), silent = TRUE)
            eg2up  <- tibble::tibble()
            if (!inherits(sym2eg, "try-error") && nrow(sym2eg)) {
                ekeys <- unique(na.omit(sym2eg$ENTREZID))
                if (length(ekeys)) {
                    egsel <- try(AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys = ekeys, keytype = "ENTREZID", columns = c("UNIPROT","ENTREZID")), silent = TRUE)
                    if (!inherits(egsel, "try-error") && nrow(egsel)) {
                        eg2up <- tibble::as_tibble(egsel) %>%
                            dplyr::filter(!is.na(UNIPROT) & nzchar(UNIPROT)) %>%
                            dplyr::group_by(ENTREZID) %>% dplyr::arrange(UNIPROT, .by_group = TRUE) %>% dplyr::slice_head(n = 1) %>% dplyr::ungroup()
                    }
                }
                if (nrow(eg2up)) {
                    map_sym2up <- tibble::as_tibble(sym2eg) %>%
                        dplyr::distinct(SYMBOL, ENTREZID) %>%
                        dplyr::left_join(eg2up, by = "ENTREZID") %>%
                        dplyr::filter(!is.na(UNIPROT) & nzchar(UNIPROT)) %>%
                        dplyr::transmute(input = toupper(SYMBOL), primaryAccession = toupper(UNIPROT)) %>%
                        dplyr::distinct(input, .keep_all = TRUE)
                    need_idx2 <- which((is.na(resolved$Resolved_UNIPROT) | !nzchar(resolved$Resolved_UNIPROT)) & resolved$id_class == "SYMBOL_OR_ALIAS")
                    base_need <- toupper(resolved$token_base[need_idx2])
                    hit <- map_sym2up$primaryAccession[match(base_need, map_sym2up$input)]
                    ok <- !is.na(hit) & nzchar(hit)
                    ii <- need_idx2[ok]
                    if (length(ii)) {
                        resolved$Resolved_UNIPROT[ii] <- hit[ok]
                        resolved$strategy[ii] <- "orgdb_symbol_entrez_uniprot"
                    }
                }
            }
        }
    } else {
        message("org.Mm.eg.db/AnnotationDbi not available; skipping OrgDb strategies")
    }

    need_idx <- which((is.na(resolved$Resolved_UNIPROT) | !nzchar(resolved$Resolved_UNIPROT)) & resolved$id_class == "SYMBOL_OR_ALIAS")
    need_ids <- unique(toupper(resolved$token_base[need_idx]))
    sym_left2 <- unique(need_ids[grepl("^[A-Z0-9\\-]{2,}$", need_ids)])
    if (length(sym_left2) && requireNamespace("UniProt.ws", quietly = TRUE)) {
        batch_vec <- split(sym_left2, ceiling(seq_along(sym_left2) / 50))
        picks <- list()
        for (b in batch_vec) {
            q_list <- lapply(b, function(g) list(organism_id = 10090, gene_primary = g))
            query_once <- function(ql) try(UniProt.ws::queryUniProt(query = ql, fields = c("accession","id","gene_primary","reviewed"), collapse = "OR", n = 10, pageSize = 10), silent = TRUE)
            res_list <- lapply(q_list, function(ql) { out <- query_once(ql); if (inherits(out, "try-error") || !is.data.frame(out)) { Sys.sleep(0.8); out <- query_once(ql) }; out })
            ok <- res_list[!vapply(res_list, inherits, logical(1), "try-error")]
            if (length(ok)) {
                tbl <- dplyr::bind_rows(lapply(ok, tibble::as_tibble))
                if (nrow(tbl)) {
                    tbl <- tbl %>% dplyr::mutate(gene_primary = toupper(.data$gene_primary), accession = toupper(.data$accession))
                    pick <- tbl %>% dplyr::group_by(gene_primary) %>% dplyr::arrange(dplyr::desc(.data$reviewed), accession, .by_group = TRUE) %>% dplyr::slice_head(n = 1) %>% dplyr::ungroup() %>% dplyr::transmute(input = gene_primary, primaryAccession = accession)
                    picks[[length(picks) + 1]] <- pick
                }
            }
        }
        if (length(picks)) {
            map_gene <- dplyr::bind_rows(picks) %>% dplyr::distinct(input, .keep_all = TRUE)
            need_idx3 <- which((is.na(resolved$Resolved_UNIPROT) | !nzchar(resolved$Resolved_UNIPROT)) & resolved$id_class == "SYMBOL_OR_ALIAS")
            base_need <- toupper(resolved$token_base[need_idx3])
            hit <- map_gene$primaryAccession[match(base_need, map_gene$input)]
            ok <- !is.na(hit) & nzchar(hit)
            ii <- need_idx3[ok]
            if (length(ii)) {
                resolved$Resolved_UNIPROT[ii] <- hit[ok]
                resolved$strategy[ii] <- "uniprot_gene_primary_retry"
            }
        }
    } else if (length(sym_left2)) {
        message("UniProt.ws not available; skipping online gene strategy")
    }

    need <- which(is.na(resolved$Resolved_UNIPROT) | !nzchar(resolved$Resolved_UNIPROT))
    idx_entry_ws <- need[resolved$id_class[need] %in% c("ENTRY","ENTRY_GUESS")]
    if (length(idx_entry_ws) && requireNamespace("UniProt.ws", quietly = TRUE)) {
        left_ids <- unique(toupper(ifelse(nz(resolved$entry_guess_up[idx_entry_ws]),
                                                                            resolved$entry_guess_up[idx_entry_ws],
                                                                            paste0(resolved$token_base[idx_entry_ws], "_MOUSE"))))
        left_ids <- gsub("-\\d+$", "", left_ids)
        batch_vec <- split(left_ids, ceiling(seq_along(left_ids) / 50))
        picks <- list()
        for (b in batch_vec) {
            q_list <- lapply(b, function(ua) list(organism_id = 10090, id = ua))
            query_once <- function(ql) try(UniProt.ws::queryUniProt(query = ql, fields = c("accession","id","reviewed"), collapse = "OR", n = 5, pageSize = 5), silent = TRUE)
            res_list <- lapply(q_list, function(ql) { out <- query_once(ql); if (inherits(out, "try-error") || !is.data.frame(out)) { Sys.sleep(0.8); out <- query_once(ql) }; out })
            ok <- res_list[!vapply(res_list, inherits, logical(1), "try-error")]
            if (length(ok)) {
                tbl <- dplyr::bind_rows(lapply(ok, tibble::as_tibble))
                if (nrow(tbl)) {
                    tbl <- tbl %>% dplyr::mutate(id = toupper(.data$id), accession = toupper(.data$accession))
                    pick <- tbl %>% dplyr::group_by(id) %>% dplyr::arrange(dplyr::desc(.data$reviewed), accession, .by_group = TRUE) %>% dplyr::slice_head(n = 1) %>% dplyr::ungroup() %>% dplyr::transmute(input = id, primaryAccession = accession)
                    picks[[length(picks) + 1]] <- pick
                }
            }
        }
        if (length(picks)) {
            map_id <- dplyr::bind_rows(picks) %>% dplyr::distinct(input, .keep_all = TRUE)
            keys_now <- toupper(ifelse(nz(resolved$entry_guess_up[idx_entry_ws]), resolved$entry_guess_up[idx_entry_ws], paste0(resolved$token_base[idx_entry_ws], "_MOUSE")))
            hit <- map_id$primaryAccession[match(keys_now, map_id$input)]
            ok <- nz(hit)
            ii <- idx_entry_ws[ok]
            if (length(ii)) {
                resolved$Resolved_UNIPROT[ii] <- hit[ok]
                resolved$strategy[ii] <- "uniprot_id_retry"
            }
        }
    }

    # Manual mapping application supporting "gene_symbol" -> "mapped_gene_symbol" (symbol-to-symbol) or accession
    if (!is.null(manual_mapping)) {
        symbol_cols <- intersect(names(manual_mapping), c("original_symbol","symbol","input","gene_symbol","token_raw"))
        mapped_cols <- intersect(names(manual_mapping), c("final_accession","accession","uniprot","uniprot_accession","mapped_gene_symbol"))
        base_cols   <- intersect(names(manual_mapping), c("base_name","token_base","base","symbol_base"))

        if (length(mapped_cols) && (length(symbol_cols) || length(base_cols))) {
            mm_symbol <- if (length(symbol_cols)) symbol_cols[1] else NA
            mm_mapped <- mapped_cols[1]
            mm_base   <- if (length(base_cols)) base_cols[1] else NA

            mm_clean <- manual_mapping
            if (!is.na(mm_symbol)) mm_clean <- mm_clean %>% dplyr::filter(nzchar(.data[[mm_symbol]]))
            if (!is.na(mm_base))   mm_clean <- mm_clean %>% dplyr::filter(nzchar(.data[[mm_base]]))
            mm_clean <- mm_clean %>% dplyr::filter(nzchar(.data[[mm_mapped]]))

            if (nrow(mm_clean)) {
                # Convert mapped values to UniProt accessions if needed
                map_to_acc <- function(vals) {
                    vals <- toupper(trimws(as.character(vals)))
                    out <- ifelse(is_uniprot_ac(vals), vals, NA_character_)
                    need <- is.na(out) | !nzchar(out)
                    if (any(need)) {
                        base <- to_base_no_iso_mouse(vals[need])
                        hit <- entry_map$UNIPROT[match(base, entry_map$entry_base)]
                        ok <- nz(hit)
                        if (any(ok)) out[need][ok] <- hit[ok]
                    }
                    need <- is.na(out) | !nzchar(out)
                    if (any(need) && nrow(gene_map)) {
                        key <- toupper(to_base_no_iso_mouse(vals[need]))
                        hit <- gene_map$primaryAccession[match(key, gene_map$input)]
                        ok <- nz(hit)
                        if (any(ok)) out[need][ok] <- hit[ok]
                    }
                    out
                }

                mm_clean$.__acc__ <- map_to_acc(mm_clean[[mm_mapped]])
                mm_clean <- mm_clean %>% dplyr::filter(!is.na(.__acc__) & nzchar(.__acc__))

                if (nrow(mm_clean)) {
                    if (!is.na(mm_symbol)) {
                        sym_map <- stats::setNames(mm_clean$.__acc__, mm_clean[[mm_symbol]])
                        tgt <- toupper(resolved$token_raw)
                        hit <- toupper(unname(sym_map[tgt]))
                        ok <- nz(hit)
                        idx <- which(ok)
                        if (length(idx)) {
                            if (isTRUE(manual_override)) {
                                resolved$Resolved_UNIPROT[idx] <- hit[idx]
                                resolved$strategy[idx] <- ifelse(is.na(resolved$strategy[idx]), "manual_symbol", paste0(resolved$strategy[idx], "|manual_symbol"))
                            } else {
                                need_idx <- idx[is.na(resolved$Resolved_UNIPROT[idx]) | !nzchar(resolved$Resolved_UNIPROT[idx])]
                                if (length(need_idx)) {
                                    resolved$Resolved_UNIPROT[need_idx] <- hit[need_idx]
                                    resolved$strategy[need_idx] <- "manual_symbol"
                                }
                            }
                        }
                    }
                    if (!is.na(mm_base)) {
                        base_map <- stats::setNames(mm_clean$.__acc__, mm_clean[[mm_base]])
                        tgt <- toupper(resolved$token_base)
                        hit <- toupper(unname(base_map[tgt]))
                        ok <- nz(hit)
                        idx <- which(ok)
                        if (length(idx)) {
                            if (isTRUE(manual_override)) {
                                resolved$Resolved_UNIPROT[idx] <- hit[idx]
                                resolved$strategy[idx] <- ifelse(is.na(resolved$strategy[idx]), "manual_base", paste0(resolved$strategy[idx], "|manual_base"))
                            } else {
                                need_idx <- idx[is.na(resolved$Resolved_UNIPROT[idx]) | !nzchar(resolved$Resolved_UNIPROT[idx])]
                                if (length(need_idx)) {
                                    resolved$Resolved_UNIPROT[need_idx] <- hit[need_idx]
                                    resolved$strategy[need_idx] <- "manual_base"
                                }
                            }
                        }
                    }
                } else {
                    message("Manual mapping present but mapped values could not be converted to UniProt accessions; skipping.")
                }
            }
        } else {
            message("Manual mapping file lacks required columns; expected one of symbol columns and 'mapped_gene_symbol' or accession column.")
        }
    }

    mapping_info <- resolved %>%
        dplyr::transmute(
            original_symbol = token_raw,
            base_name = token_base,
            final_accession = Resolved_UNIPROT,
            matched_by = strategy
        )

    df_joined <- df_tok %>%
        dplyr::left_join(mapping_info, by = c("token_raw" = "original_symbol", "token_base" = "base_name"))

    df_mapped <- df_joined %>%
        dplyr::filter(!is.na(final_accession) & nzchar(final_accession)) %>%
        dplyr::mutate(gene_symbol = final_accession)

    keep_cols <- intersect(c("gene_symbol", "pval", "padj", "log2fc", "P.Value", "adj.P.Val", "logFC"), names(df_mapped))
    if (!length(keep_cols)) keep_cols <- "gene_symbol"
    df_mapped <- df_mapped %>% dplyr::select(dplyr::all_of(keep_cols)) %>% dplyr::distinct()

    unmapped_proteins <- df_joined %>%
        dplyr::filter(is.na(final_accession) | !nzchar(final_accession)) %>%
        dplyr::transmute(gene_symbol = token_raw) %>%
        dplyr::distinct()

    base <- tools::file_path_sans_ext(basename(data_path))
    mapped_file <- file.path(mapped_dir, paste0(base, ".csv"))
    unmapped_file <- file.path(unmapped_dir, paste0(base, ".csv"))
    info_table_file <- file.path(mapped_dir, paste0(base, "_mapping_info.csv"))
    info_summary_file <- file.path(mapped_dir, paste0(base, "_info.txt"))

    readr::write_csv(df_mapped, mapped_file)
    readr::write_csv(unmapped_proteins, unmapped_file)
    readr::write_csv(mapping_info, info_table_file)

    total_in <- nrow(df_tok)
    total_mapped <- sum(!is.na(mapping_info$final_accession) & nzchar(mapping_info$final_accession))
    total_unmapped <- total_in - total_mapped
    strategy_counts <- mapping_info %>%
        dplyr::filter(!is.na(matched_by) & nzchar(matched_by)) %>%
        dplyr::count(matched_by, name = "n") %>%
        dplyr::arrange(dplyr::desc(n))

    summary_lines <- c(
        paste0("file: ", basename(data_path)),
        paste0("total_valid_input_after__MOUSE_filter: ", total_in),
        paste0("mapped: ", total_mapped),
        paste0("unmapped: ", total_unmapped),
        "strategy_counts:"
    )
    if (nrow(strategy_counts) > 0) {
        strategy_str <- paste0("  - ", strategy_counts$matched_by, ": ", strategy_counts$n)
        summary_lines <- c(summary_lines, strategy_str)
    } else {
        summary_lines <- c(summary_lines, "  - none")
    }
    writeLines(summary_lines, info_summary_file)

    message("Saved mapped -> ", mapped_file)
    message("Saved unmapped -> ", unmapped_file, " (", nrow(unmapped_proteins), " entries )")
    message("Saved mapping info -> ", info_table_file)
    message("Saved summary info -> ", info_summary_file)

    invisible(list(
        mapped = mapped_file,
        unmapped = unmapped_file,
        info = info_table_file,
        summary = info_summary_file,
        unmapped_table = unmapped_proteins  # collected per-file
    ))
}
# NOTE: Add "manual_mapping","manual_override" to .export in foreach if running in parallel.

# -------------------------
# Parallel execution block
# -------------------------
n_files <- length(csv_files)
available_cores <- parallel::detectCores(logical = FALSE)
# leave one core free when possible
workers <- max(1, min(available_cores - 1, n_files))

cat("Starting parallel processing with", workers, "workers for", n_files, "files...\n")
cl <- parallel::makeCluster(workers)
doParallel::registerDoParallel(cl)

# Export the necessary objects and ensure packages are loaded on workers
results <- foreach(i = seq_along(csv_files),
                   .packages = c("dplyr", "stringr", "tidyr", "purrr", "readr", "R.utils"),
                   .export = c("uniprot_mapping", "entry_name_to_accession", "mapped_dir", "unmapped_dir", "process_file")) %dopar% {
    process_file(csv_files[i])
}

parallel::stopCluster(cl)
cat("Batch mapping completed for", length(csv_files), "files.\n")
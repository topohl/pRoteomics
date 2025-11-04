#' Batch UniProt ID Mapping for ClusterProfiler Analysis (parallelized with doParallel)
#' Processes all .csv files in Datasets/raw and maps UniProtKB IDs to Accessions.

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, stringr, tidyr, purrr, readr, R.utils, foreach, doParallel)

mapped_comparisons <- "unknown-comparison"  # specify the comparison folder to process

#working_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/"
working_dir <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler"
raw_dir <- file.path(working_dir, "Datasets", "raw", mapped_comparisons, "reverse")
mapped_dir <- file.path(working_dir, "Datasets", "mapped", mapped_comparisons, "reverse")
unmapped_dir <- file.path(working_dir, "Datasets", "unmapped", mapped_comparisons, "reverse")

dir.create(mapped_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(unmapped_dir, recursive = TRUE, showWarnings = FALSE)

setwd(working_dir)

uniprot_mapping_file_path <- file.path(working_dir, "Datasets", "MOUSE_10090_idmapping.dat")

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

# find CSV files
csv_files <- list.files(raw_dir, pattern = ".*_.*\\.csv$", full.names = TRUE)
if (length(csv_files) == 0) stop("No .csv files found in: ", raw_dir)

process_file <- function(data_path) {
    message("Processing file: ", data_path)

    # read CSV (try comma first, fall back to semicolon)
    df_raw <- tryCatch(
        readr::read_csv(data_path, col_names = TRUE, show_col_types = FALSE, trim_ws = TRUE, quote = "\""),
        error = function(e) {
            message("read_csv failed, trying read_csv2(): ", e$message)
            readr::read_csv2(data_path, col_names = TRUE, show_col_types = FALSE, trim_ws = TRUE)
        }
    )

    # ensure a gene_symbol column exists (assume first column otherwise)
    if (!"gene_symbol" %in% names(df_raw)) {
        names(df_raw)[1] <- "gene_symbol"
        message("Renamed first column to 'gene_symbol' for file: ", data_path)
    }

    df <- df_raw %>%
        mutate(gene_symbol = as.character(gene_symbol)) %>%
        mutate(gene_symbol = purrr::map_chr(stringr::str_split(gene_symbol, ";"), ~ stringr::str_trim(.x[1]))) %>%
        filter(gene_symbol != "" & !is.na(gene_symbol)) %>%
        distinct(gene_symbol, .keep_all = TRUE)

    if (nrow(df) == 0) {
        message("No valid gene symbols in file: ", data_path)
        return(invisible(NULL))
    }

    # build synonym map (Gene_Name and Gene_Synonym), normalized to lowercase for joins
    synonym_map <- uniprot_mapping %>%
        filter(Type %in% c("Gene_Name", "Gene_Synonym")) %>%
        tidyr::separate_rows(Value, sep = ";") %>%
        mutate(Value = stringr::str_trim(Value)) %>%
        filter(Value != "" & !is.na(Value)) %>%
        transmute(
            Syn_Accession = UniProt_Accession,
            GeneName = Value,
            GeneName_lower = stringr::str_to_lower(Value)
        ) %>%
        distinct(Syn_Accession, GeneName_lower, .keep_all = TRUE)

    # prepare entry name map (do not overwrite global variable)
    entry_map <- entry_name_to_accession %>%
        mutate(UniProtKB_ID_lower = stringr::str_to_lower(UniProtKB_ID))

    # prepare dataframe for matching: derive base_name if _MOUSE suffix present
    df_prepped <- df %>%
        mutate(
            base_name = if_else(stringr::str_detect(gene_symbol, "_MOUSE$"),
                                                    stringr::str_remove(gene_symbol, "_MOUSE$"),
                                                    NA_character_),
            gene_symbol_lower = stringr::str_to_lower(gene_symbol),
            base_name_lower = if_else(!is.na(base_name), stringr::str_to_lower(base_name), NA_character_)
        )

    # join attempts:
    # 1) match full entry name (case-insensitive) using entry_map
    # 2) match base_name to synonyms / gene names
    joined <- df_prepped %>%
        left_join(entry_map, by = c("gene_symbol_lower" = "UniProtKB_ID_lower")) %>%
        left_join(synonym_map %>% select(Syn_Accession, GeneName_lower),
                            by = c("base_name_lower" = "GeneName_lower"))

    # determine final accession and whether row is considered an accession
    resolved <- joined %>%
        mutate(
            final_accession = dplyr::coalesce(UniProt_Accession, Syn_Accession),
            is_accession = if_else(
                !is.na(final_accession) & final_accession != "",
                TRUE,
                stringr::str_detect(gene_symbol, "^[A-Za-z][A-Za-z0-9]{5,9}$")
            )
        )

    # mapped proteins: keep only rows resolved as accessions, prefer resolved accession as gene_symbol
    df_mapped <- resolved %>%
        filter(is_accession) %>%
        mutate(gene_symbol = if_else(!is.na(final_accession) & final_accession != "", final_accession, gene_symbol))

    # select only the common expected output columns if they exist
    keep_cols <- intersect(c("gene_symbol", "pval", "padj", "log2fc", "P.Value", "adj.P.Val", "logFC"), names(df_mapped))
    df_mapped <- df_mapped %>%
        select(all_of(keep_cols)) %>%
        distinct()

    # unmapped proteins: original gene_symbol values that could not be resolved
    unmapped_proteins <- resolved %>%
        filter(!is_accession) %>%
        select(gene_symbol) %>%
        distinct()

    # write outputs
    base <- tools::file_path_sans_ext(basename(data_path))
    mapped_file <- file.path(mapped_dir, paste0(base, ".csv"))
    unmapped_file <- file.path(unmapped_dir, paste0(base, ".csv"))

    readr::write_csv(df_mapped, mapped_file)
    readr::write_csv(unmapped_proteins, unmapped_file)

    message("Saved mapped -> ", mapped_file)
    message("Saved unmapped -> ", unmapped_file, " (", nrow(unmapped_proteins), " entries )")

    invisible(list(mapped = mapped_file, unmapped = unmapped_file))
}

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
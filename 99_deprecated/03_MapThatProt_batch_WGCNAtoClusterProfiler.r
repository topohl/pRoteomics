# -------------------------------
# Map WGCNA modules to UniProtKB Accessions
# Reuses local UniProt idmapping.dat (MOUSE_10090_idmapping.dat)
# Produces accession-mapped module CSVs compatible with clusterProfiler keyType="UNIPROT"
# -------------------------------

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, stringr, tidyr, purrr, readr, R.utils, tools, tibble)

working_dir   <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/"
modules_dir   <- file.path(working_dir, "wgcna", "output", "tables_modules")
out_root      <- file.path(working_dir, "wgcna", "output", "wgcna_modules_mapped")
unmapped_root <- file.path(working_dir, "wgcna", "output", "wgcna_modules_unmapped")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(unmapped_root, recursive = TRUE, showWarnings = FALSE)

uniprot_mapping_file_path <- file.path(working_dir, "Datasets", "MOUSE_10090_idmapping.dat")
if (!file.exists(uniprot_mapping_file_path)) {
  options(timeout = 3600)
  gz_url  <- "https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/idmapping/by_organism/MOUSE_10090_idmapping.dat.gz"
  gz_file <- paste0(uniprot_mapping_file_path, ".gz")
  download.file(gz_url, gz_file, mode = "wb")
  R.utils::gunzip(gz_file, destname = uniprot_mapping_file_path, remove = TRUE)
}

uniprot_mapping <- readr::read_tsv(
  uniprot_mapping_file_path,
  col_names = c("UniProt_Accession", "Type", "Value"),
  col_types = "ccc",
  quote = ""
) |> tibble::as_tibble()

entry_name_to_accession <- uniprot_mapping |>
  dplyr::filter(Type == "UniProtKB-ID") |>
  dplyr::select(UniProt_Accession, UniProtKB_ID = Value) |>
  dplyr::distinct(UniProtKB_ID, .keep_all = TRUE) |>
  dplyr::mutate(UniProtKB_ID_lower = stringr::str_to_lower(UniProtKB_ID))

stopifnot(nrow(entry_name_to_accession) > 0)

synonym_map <- uniprot_mapping |>
  dplyr::filter(Type %in% c("Gene_Name", "Gene_Synonym")) |>
  tidyr::separate_rows(Value, sep = ";") |>
  dplyr::mutate(Value = stringr::str_trim(Value)) |>
  dplyr::filter(Value != "" & !is.na(Value)) |>
  dplyr::transmute(Syn_Accession = UniProt_Accession,
                   GeneName_lower = stringr::str_to_lower(Value)) |>
  dplyr::distinct(GeneName_lower, .keep_all = TRUE)

looks_like_accession_core <- function(x) stringr::str_detect(x, "^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$|^[OPQ][0-9][A-Z0-9]{3}[0-9]$")
looks_like_accession_loose <- function(x) stringr::str_detect(x, "^[A-Z0-9]{6,10}$")
looks_like_accession_with_mouse <- function(x) stringr::str_detect(x, "^[A-Z0-9]{6,10}_MOUSE$")

strip_mouse_suffix_if_accession <- function(x) {
  y <- ifelse(looks_like_accession_with_mouse(x), stringr::str_replace(x, "_MOUSE$", ""), NA_character_)
  y
}

map_module_file <- function(fp, out_dir = out_root, unmapped_dir = unmapped_root) {
  message("Mapping module file: ", fp)
  df_raw <- utils::read.csv(fp, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  if (!("Gene" %in% names(df_raw))) stop("File missing 'Gene' column: ", fp)

  df <- df_raw |>
    dplyr::mutate(Gene = as.character(Gene)) |>
    dplyr::filter(!is.na(Gene), Gene != "")

  if (nrow(df) == 0) {
    message("No valid Gene entries in: ", fp)
    return(invisible(NULL))
  }

  df_prep <- df |>
    dplyr::mutate(
      Gene_trim   = stringr::str_trim(Gene),
      Gene_lower  = stringr::str_to_lower(Gene_trim),
      base_name   = dplyr::if_else(stringr::str_ends(Gene_trim, "_MOUSE"),
                                   stringr::str_remove(Gene_trim, "_MOUSE$"),
                                   NA_character_),
      base_lower  = dplyr::if_else(!is.na(base_name), stringr::str_to_lower(base_name), NA_character_),
      accession_from_suffix = strip_mouse_suffix_if_accession(Gene_trim),
      is_accession_core  = looks_like_accession_core(Gene_trim),
      is_accession_loose = looks_like_accession_loose(Gene_trim)
    )

  entry_map <- entry_name_to_accession |>
    dplyr::select(UniProt_Accession, UniProtKB_ID_lower)

  joined <- df_prep |>
    dplyr::left_join(entry_map, by = c("Gene_lower" = "UniProtKB_ID_lower"))

  joined2 <- joined |>
    dplyr::left_join(synonym_map, by = c("base_lower" = "GeneName_lower"))

  joined3 <- joined2 |>
    dplyr::mutate(
      final_accession = dplyr::coalesce(
        UniProt_Accession,
        Syn_Accession,
        accession_from_suffix,
        dplyr::if_else(is_accession_core | is_accession_loose, Gene_trim, NA_character_)
      )
    )

  mapped <- joined3 |>
    dplyr::filter(!is.na(final_accession) & final_accession != "") |>
    dplyr::mutate(Gene = final_accession)

  # Reorder to original df_raw column order if present
  common_cols <- intersect(names(df_raw), names(mapped))
  mapped <- dplyr::select(mapped, dplyr::all_of(common_cols)) |>
    dplyr::distinct()

  unmapped <- joined3 |>
    dplyr::filter(is.na(final_accession) | final_accession == "") |>
    dplyr::select(Gene) |>
    dplyr::distinct()

  base <- tools::file_path_sans_ext(basename(fp))
  out_file <- file.path(out_dir, paste0(base, "_mapped.csv"))
  unmapped_file <- file.path(unmapped_dir, paste0(base, "_unmapped.csv"))
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(unmapped_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(mapped, out_file)
  readr::write_csv(unmapped, unmapped_file)

  tibble::tibble(
    file = basename(fp),
    n_input = nrow(df),
    n_mapped = nrow(mapped),
    n_unmapped = nrow(unmapped),
    pct_mapped = round(100 * nrow(mapped) / max(1, nrow(df)), 1)
  )
}

module_files <- list.files(modules_dir, pattern = "\\.csv$", full.names = TRUE)
stopifnot(length(module_files) > 0)

stats_list <- lapply(module_files, map_module_file)
stats <- dplyr::bind_rows(stats_list)
readr::write_csv(stats, file.path(out_root, "module_mapping_stats.csv"))
print(stats)

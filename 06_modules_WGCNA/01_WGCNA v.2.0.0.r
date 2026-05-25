# ================================ parallel-enabled
# WGCNA with profile-aware spatial traits, preservation, and condition panels
# Outputs organized into subfolders under output_dir
# ================================ parallel-enabled

# Packages
required_pkgs <- c(
  "WGCNA", "flashClust", "curl", "readxl", "ggplot2", "svglite", "GO.db",
  "reshape2", "gtools", "patchwork", "cowplot", "pheatmap", "dplyr", "tidyr",
  "httr", "jsonlite", "purrr", "AnnotationDbi", "org.Mm.eg.db", "readr",
  "stringr", "tibble", "UniProt.ws", "RColorBrewer", "ggpubr", "broom", "grid",
  "clusterProfiler", "scales"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop(
    "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
    ". Install them explicitly before running this manuscript pipeline."
  )
}
suppressPackageStartupMessages(
  invisible(lapply(required_pkgs, library, character.only = TRUE))
)

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)

mm_to_in <- function(mm) mm / 25.4
nature_single_col <- mm_to_in(89)
nature_double_col <- mm_to_in(183)
nature_font <- "Arial"
nature_base_size <- 7
nature_axis_size <- 6.2
nature_title_size <- 7
nature_line <- 0.25
nature_diverging <- c(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B")
nature_condition_cols <- c(con = "#4D4D4D", res = "#0072B2", sus = "#D55E00")
nature_condition_labels <- c(con = "CON", res = "RES", sus = "SUS")

theme_nature <- function(base_size = nature_base_size) {
  ggplot2::theme_classic(base_size = base_size, base_family = nature_font) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = nature_title_size, face = "plain", hjust = 0),
      plot.subtitle = ggplot2::element_text(size = nature_axis_size, color = "grey30", margin = ggplot2::margin(t = 1, b = 2)),
      axis.title = ggplot2::element_text(size = nature_axis_size),
      axis.text = ggplot2::element_text(size = nature_axis_size, color = "black"),
      axis.line = ggplot2::element_line(linewidth = nature_line, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = nature_line, color = "black"),
      legend.title = ggplot2::element_text(size = nature_axis_size),
      legend.text = ggplot2::element_text(size = nature_axis_size),
      legend.key.size = grid::unit(3, "mm"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = nature_axis_size, color = "black"),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(3, 3, 3, 3)
    )
}
ggplot2::theme_set(theme_nature())

# Parallel setup
nCores <- tryCatch({
  pc <- parallel::detectCores(logical = FALSE)
  if (is.na(pc) || pc < 2) 2 else pc
}, error = function(e) 2)
enableWGCNAThreads(nThreads = nCores)
WGCNAnThreads()

# --------------------------
# Paths and data load
# --------------------------
output_dir <- Sys.getenv("PROTEOMICS_WGCNA_OUTPUT_DIR", unset = path_results("06_modules_WGCNA", "wgcna_output"))

subdirs <- list(
  figures_qc          = file.path(output_dir, "figures", "qc"),
  figures_network     = file.path(output_dir, "figures", "network"),
  figures_traits      = file.path(output_dir, "figures", "traits"),
  figures_main        = file.path(output_dir, "figures", "main"),
  tables_qc           = file.path(output_dir, "tables", "qc"),
  tables_mapping      = file.path(output_dir, "tables", "mapping"),
  tables_modules      = file.path(output_dir, "tables", "modules"),
  tables_pres         = file.path(output_dir, "tables", "preservation"),
  tables_traits       = file.path(output_dir, "tables", "traits"),
  source_data         = file.path(output_dir, "source_data"),
  state               = file.path(output_dir, "state"),
  logs                = file.path(output_dir, "logs")
)

safe_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (file.access(path, 2) != 0) stop(sprintf("Not writable: %s", path))
  invisible(normalizePath(path))
}
invisible(lapply(c(output_dir, unlist(subdirs)), safe_dir))

# Path helpers
fp_qc       <- function(...) file.path(subdirs$figures_qc, ...)
fp_traits   <- function(...) file.path(subdirs$figures_traits, ...)
fp_mainfig  <- function(...) file.path(subdirs$figures_main, ...)
fp_net      <- function(...) file.path(subdirs$figures_network, ...)
fp_qctab    <- function(...) file.path(subdirs$tables_qc, ...)
fp_maptab   <- function(...) file.path(subdirs$tables_mapping, ...)
fp_modtab   <- function(...) file.path(subdirs$tables_modules, ...)
fp_traittab <- function(...) file.path(subdirs$tables_traits, ...)
fp_prestab  <- function(...) file.path(subdirs$tables_pres, ...)
fp_source   <- function(...) file.path(subdirs$source_data, ...)
fp_state    <- function(...) file.path(subdirs$state, ...)
fp_log      <- function(...) file.path(subdirs$logs, ...)

write_csv_safe <- function(x, path) {
  readr::write_csv(x, path, na = "")
  invisible(path)
}
write_tsv_safe <- function(x, path) {
  readr::write_tsv(x, path, na = "")
  invisible(path)
}
log_session <- function() {
  writeLines(capture.output(utils::sessionInfo()), fp_log("session_info.txt"))
}
log_session()

# Analysis parameters reported with outputs
sample_tree_cut_height <- 80
sample_tree_plot_height <- 40
soft_threshold_rsquared <- 0.80
min_module_size <- 30
deep_split <- 2
merge_cut_height <- 0.25
module_preservation_permutations <- 1000
dataset_profile <- "auto"  # one of: auto, microglia, neuron_soma, neuron_neuropil

# Optional: safe svg helper
save_svg <- function(path, width, height, expr) {
  svglite::svglite(file = path, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  force(expr)
}

save_plot_nature <- function(plot, path_svg, width, height) {
  ggplot2::ggsave(path_svg, plot, device = svglite::svglite,
                  width = width, height = height, units = "in", limitsize = FALSE)
  ggplot2::ggsave(sub("\\.svg$", ".pdf", path_svg), plot, device = grDevices::pdf,
                  width = width, height = height, units = "in", limitsize = FALSE,
                  family = nature_font, useDingbats = FALSE)
  invisible(path_svg)
}

sig_dot <- function(fdr) {
  dplyr::case_when(
    is.na(fdr) ~ "",
    fdr < 0.01 ~ "\u2022",
    fdr < 0.05 ~ "\u00b7",
    TRUE ~ ""
  )
}

expr_xlsx <- Sys.getenv("PROTEOMICS_WGCNA_EXPR_XLSX", unset = path_processed("variancePartition", "data", "male.data.xlsx"))
meta_xlsx <- Sys.getenv("PROTEOMICS_WGCNA_META_XLSX", unset = path_processed("variancePartition", "data", "sample_info.xlsx"))

# ================================
# Mouse-only mapping: robust idmapping parser + offline + SYMBOL/ALIAS + Entrez + UniProt gene_primary + QC
# ================================

# --------------------------
# Inputs
# --------------------------
expr_xlsx <- Sys.getenv("PROTEOMICS_WGCNA_EXPR_XLSX", unset = expr_xlsx)
meta_xlsx <- Sys.getenv("PROTEOMICS_WGCNA_META_XLSX", unset = meta_xlsx)
idmap_dat <- Sys.getenv("PROTEOMICS_WGCNA_IDMAP_DAT", unset = path_external("MOUSE_10090_idmapping.dat"))

stop_if_missing <- function(path) if (!file.exists(path)) stop(sprintf("Missing file: %s", path))
read_head <- function(path) {
  df <- readxl::read_excel(path)
  utils::write.table(utils::head(df, 10), fp_log(paste0(basename(path), "_head10.tsv")),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  df
}

norm_label <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[[:space:]-]+", "_", x)
  x[is.na(x) | !nzchar(x)] <- NA_character_
  x
}

infer_dataset_profile <- function(sample_info, expr_path, meta_path, requested = "auto") {
  requested <- tolower(requested)
  if (!identical(requested, "auto")) return(requested)
  path_hint <- tolower(paste(expr_path, meta_path, collapse = " "))
  if (grepl("microglia", path_hint)) return("microglia")
  if (grepl("neuropil", path_hint)) return("neuron_neuropil")
  if (grepl("soma", path_hint)) return("neuron_soma")
  cell_hint <- if ("celltype" %in% names(sample_info)) paste(unique(norm_label(sample_info$celltype)), collapse = " ") else ""
  if (grepl("microglia", cell_hint)) return("microglia")
  if (grepl("neuropil", cell_hint)) return("neuron_neuropil")
  if (grepl("soma", cell_hint)) return("neuron_soma")
  layer_vals <- if ("layer" %in% names(sample_info)) unique(stats::na.omit(norm_label(sample_info$layer))) else character()
  if (length(layer_vals) && any(layer_vals %in% c("so", "sr", "slm", "mo", "po"))) return("neuron_neuropil")
  if (!length(layer_vals) || all(layer_vals %in% c("sp", "sg", "none"))) return("neuron_soma")
  "region_only"
}

prepare_spatial_metadata <- function(sample_info, profile) {
  if (!"region" %in% names(sample_info)) stop("sample_info must contain a 'region' column")
  if (!"ExpGroup" %in% names(sample_info)) stop("sample_info must contain an 'ExpGroup' column")

  profile <- tolower(profile)
  allowed_profiles <- c("microglia", "neuron_soma", "neuron_neuropil", "region_only")
  if (!profile %in% allowed_profiles) {
    stop("Unsupported dataset_profile: ", profile, ". Use one of: ", paste(c("auto", allowed_profiles), collapse = ", "))
  }
  sample_info$region <- factor(norm_label(sample_info$region), levels = c("ca1", "ca2", "ca3", "dg"))
  sample_info$condition <- factor(norm_label(sample_info$ExpGroup), levels = c("con", "res", "sus"))
  sample_info$ExpGroup <- sample_info$condition

  raw_layer <- if ("layer" %in% names(sample_info)) norm_label(sample_info$layer) else rep(NA_character_, nrow(sample_info))
  raw_celltype <- if ("celltype" %in% names(sample_info)) norm_label(sample_info$celltype) else rep(NA_character_, nrow(sample_info))

  if (identical(profile, "microglia")) {
    sample_info$layer <- factor("none")
    sample_info$celltype <- factor("microglia")
    active_spatial_vars <- c("region")
  } else if (identical(profile, "neuron_soma")) {
    sample_info$soma_layer <- factor(raw_layer, levels = c("sp", "sg"))
    sample_info$layer <- factor("none")
    sample_info$celltype <- factor("neuron_soma")
    active_spatial_vars <- c("region")
  } else if (identical(profile, "neuron_neuropil")) {
    sample_info$layer <- factor(raw_layer, levels = c("so", "sr", "slm", "mo", "po"))
    sample_info$celltype <- factor("neuron_neuropil")
    active_spatial_vars <- c("region", "layer")
    invalid_layer <- (!is.na(sample_info$region) & !is.na(sample_info$layer)) & (
      (sample_info$region %in% c("ca1", "ca2", "ca3") & !sample_info$layer %in% c("so", "sr", "slm")) |
        (sample_info$region == "dg" & !sample_info$layer %in% c("mo", "po"))
    )
    if (any(invalid_layer)) {
      write_csv_safe(sample_info[invalid_layer, , drop = FALSE], fp_log("invalid_region_layer_combinations.csv"))
      warning("Invalid region/layer combinations found; see logs/invalid_region_layer_combinations.csv")
    }
  } else {
    sample_info$layer <- factor(ifelse(is.na(raw_layer), "none", raw_layer))
    sample_info$celltype <- factor(ifelse(is.na(raw_celltype), profile, raw_celltype))
    active_spatial_vars <- c("region")
    if (length(unique(stats::na.omit(sample_info$layer))) > 1) active_spatial_vars <- c(active_spatial_vars, "layer")
    if (length(unique(stats::na.omit(sample_info$celltype))) > 1) active_spatial_vars <- c(active_spatial_vars, "celltype")
  }

  required_vars <- c("condition", active_spatial_vars)
  missing_required <- Reduce(`|`, lapply(required_vars, function(v) is.na(sample_info[[v]])))
  if (any(missing_required)) {
    write_csv_safe(sample_info[missing_required, , drop = FALSE], fp_log("invalid_required_metadata.csv"))
    stop("Missing or invalid values in required metadata columns: ",
         paste(required_vars, collapse = ", "),
         ". See logs/invalid_required_metadata.csv")
  }

  write_csv_safe(
    tibble::tibble(
      dataset_profile = profile,
      active_spatial_vars = paste(active_spatial_vars, collapse = ","),
      n_samples = nrow(sample_info),
      regions = paste(sort(unique(stats::na.omit(as.character(sample_info$region)))), collapse = ","),
      layers = paste(sort(unique(stats::na.omit(as.character(sample_info$layer)))), collapse = ",")
    ),
    fp_log("dataset_profile.csv")
  )

  list(sample_info = sample_info, profile = profile, active_spatial_vars = active_spatial_vars)
}

male.data <- { stop_if_missing(expr_xlsx); read_head(expr_xlsx) } %>% dplyr::mutate(.row_id = dplyr::row_number())
meta.data <- { stop_if_missing(meta_xlsx); read_head(meta_xlsx) }

# --------------------------
# Robust UniProtKB-ID -> Accession map (offline)
# --------------------------
stop_if_missing(idmap_dat)
idmap_tbl <- readr::read_tsv(idmap_dat, col_names = c("ACC","DB","VAL"), col_types = "ccc", progress = FALSE, quote = "", comment = "")
idmap_uid <- idmap_tbl %>%
  dplyr::filter(DB == "UniProtKB-ID" & grepl("_MOUSE\\s*$", VAL) & nzchar(ACC)) %>%
  dplyr::transmute(
    UNIPROT    = toupper(trimws(ACC)),
    entry_full = toupper(trimws(VAL)),
    entry_base = toupper(gsub("_MOUSE$", "", trimws(VAL)))
  )
entry_map <- idmap_uid %>% dplyr::distinct(entry_base, .keep_all = TRUE)
if (!nrow(entry_map)) stop("entry_map is empty after robust parse")

# Sentinel sanity check
sentinels <- c("AIF1","AKAP2","ADCY1","AKAP1","AMPD3","ANXA3","1433S","ACK1","AIP","ADA10")
missing_sentinels <- setdiff(sentinels, entry_map$entry_base)
if (length(missing_sentinels)) {
  warning(sprintf("entry_map missing expected keys: %s", paste(missing_sentinels, collapse=", ")))
  write_tsv_safe(entry_map, fp_maptab("entry_map_debug.tsv"))
}

# --------------------------
# Mouse-only tokenization and classification
# --------------------------
normalize_token <- function(x) { x <- toupper(gsub("\\s+", "", x)); x <- gsub("\\u00A0", "", x); x <- gsub("\\.+", ".", x); x <- gsub("__+", "_", x); x }
to_base_no_iso_mouse <- function(x) { x <- gsub("-\\d+$", "", x); gsub("_MOUSE$", "", x) }

tokenize_mouse_only <- function(male_df) {
  tok <- male_df %>% tidyr::separate_rows(gene_symbol, sep = ";") %>% dplyr::mutate(token_raw = gene_symbol, token_up = normalize_token(gene_symbol))
  dropped_non_mouse <- tok %>% dplyr::filter(!grepl("_MOUSE$", token_up))
  if (nrow(dropped_non_mouse)) write_tsv_safe(dropped_non_mouse, fp_maptab("dropped_non_mouse_tokens.tsv"))
  tok %>%
    dplyr::filter(grepl("_MOUSE$", token_up)) %>%
    dplyr::mutate(
      token_base = to_base_no_iso_mouse(token_up),
      looks_ac   = grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$|^A0A[0-9A-Z]{7}$", token_base),
      looks_entry= grepl("^[A-Z0-9][A-Z0-9\\-\\.]+$", token_base),
      id_class = dplyr::case_when(
        looks_ac    ~ "UNIPROT_AC_MOUSE",
        looks_entry ~ "UNIPROT_ENTRY",
        TRUE        ~ "UNKNOWN"
      ),
      Resolved_UNIPROT = NA_character_,
      strategy = NA_character_
    )
}
resolved2 <- tokenize_mouse_only(male.data)

# --------------------------
# Mapping stack
# --------------------------

# 1) Accept accession-like bases
idx_ac <- which(resolved2$id_class == "UNIPROT_AC_MOUSE")
if (length(idx_ac)) {
  resolved2$Resolved_UNIPROT[idx_ac] <- resolved2$token_base[idx_ac]
  resolved2$strategy[idx_ac] <- "accept_accession_base"
}

# 2) Offline entry map (now robust)
idx_en <- which(resolved2$id_class == "UNIPROT_ENTRY" & (is.na(resolved2$Resolved_UNIPROT) | !nzchar(resolved2$Resolved_UNIPROT)))
if (length(idx_en)) {
  hit <- entry_map$UNIPROT[match(toupper(resolved2$token_base[idx_en]), entry_map$entry_base)]
  ok <- !is.na(hit) & nzchar(hit)
  if (any(ok)) { ii <- idx_en[ok]; resolved2$Resolved_UNIPROT[ii] <- hit[ok]; resolved2$strategy[ii] <- "entry_local_mouse" }
}

# 3) SYMBOL/ALIAS offline resolver (MGI-first)
need_idx <- which(is.na(resolved2$Resolved_UNIPROT) | !nzchar(resolved2$Resolved_UNIPROT))
need_ids <- toupper(unique(resolved2$token_base[need_idx]))
is_acc <- grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$|^A0A[0-9A-Z]{7}$", need_ids)
ids_ent <- unique(need_ids[!is_acc])

if (length(ids_ent)) {
  sel_sym <- try(AnnotationDbi::select(org.Mm.eg.db, keys = ids_ent, keytype = "SYMBOL", columns = c("MGIID","ENTREZID","UNIPROT","SYMBOL")), silent = TRUE)
  map_sym <- tibble::tibble()
  if (!inherits(sel_sym, "try-error") && nrow(sel_sym)) {
    map_sym <- tibble::as_tibble(sel_sym) %>%
      dplyr::filter(!is.na(UNIPROT) & nzchar(UNIPROT)) %>%
      dplyr::group_by(SYMBOL) %>% dplyr::arrange(UNIPROT, .by_group = TRUE) %>% dplyr::slice_head(n=1) %>% dplyr::ungroup() %>%
      dplyr::transmute(input = toupper(SYMBOL), primaryAccession = toupper(UNIPROT))
  }
  kt <- try(AnnotationDbi::keytypes(org.Mm.eg.db), silent = TRUE)
  map_alias <- tibble::tibble()
  if (!inherits(kt, "try-error") && "ALIAS" %in% kt) {
    sel_alias <- try(AnnotationDbi::select(org.Mm.eg.db, keys = ids_ent, keytype = "ALIAS", columns = c("UNIPROT","ALIAS")), silent = TRUE)
    if (!inherits(sel_alias, "try-error") && nrow(sel_alias)) {
      map_alias <- tibble::as_tibble(sel_alias) %>%
        dplyr::filter(!is.na(UNIPROT) & nzchar(UNIPROT)) %>%
        dplyr::group_by(ALIAS) %>% dplyr::arrange(UNIPROT, .by_group = TRUE) %>% dplyr::slice_head(n=1) %>% dplyr::ungroup() %>%
        dplyr::transmute(input = toupper(ALIAS), primaryAccession = toupper(UNIPROT))
    }
  }
  map_symall <- dplyr::bind_rows(map_sym, map_alias) %>% dplyr::distinct(input, .keep_all = TRUE)
  if (nrow(map_symall)) {
    need_idx_now <- which(is.na(resolved2$Resolved_UNIPROT) | !nzchar(resolved2$Resolved_UNIPROT))
    base_need <- toupper(resolved2$token_base[need_idx_now])
    hit <- map_symall$primaryAccession[match(base_need, map_symall$input)]
    ok <- !is.na(hit) & nzchar(hit)
    ii <- need_idx_now[ok]
    if (length(ii)) { resolved2$Resolved_UNIPROT[ii] <- hit[ok]; resolved2$strategy[ii] <- "orgdb_mgi_symbol_first" }
  }
}

# 4) SYMBOL -> Entrez -> UniProt (offline two-hop)
need_idx <- which(is.na(resolved2$Resolved_UNIPROT) | !nzchar(resolved2$Resolved_UNIPROT))
need_ids <- toupper(unique(resolved2$token_base[need_idx]))
sym_left <- unique(need_ids[grepl("^[A-Z0-9\\-]{2,}$", need_ids)])

if (length(sym_left)) {
  sym2eg <- try(AnnotationDbi::select(org.Mm.eg.db, keys = sym_left, keytype = "SYMBOL", columns = c("ENTREZID","SYMBOL")), silent = TRUE)
  eg2up  <- tibble::tibble()
  if (!inherits(sym2eg, "try-error") && nrow(sym2eg)) {
    ekeys <- unique(na.omit(sym2eg$ENTREZID))
    if (length(ekeys)) {
      egsel <- try(AnnotationDbi::select(org.Mm.eg.db, keys = ekeys, keytype = "ENTREZID", columns = c("UNIPROT","ENTREZID")), silent = TRUE)
      if (!inherits(egsel, "try-error") && nrow(egsel)) {
        eg2up <- tibble::as_tibble(egsel) %>%
          dplyr::filter(!is.na(UNIPROT) & nzchar(UNIPROT)) %>%
          dplyr::group_by(ENTREZID) %>% dplyr::arrange(UNIPROT, .by_group = TRUE) %>% dplyr::slice_head(n=1) %>% dplyr::ungroup()
      }
    }
    if (nrow(eg2up)) {
      map_sym2up <- tibble::as_tibble(sym2eg) %>%
        dplyr::distinct(SYMBOL, ENTREZID) %>%
        dplyr::left_join(eg2up, by = "ENTREZID") %>%
        dplyr::filter(!is.na(UNIPROT) & nzchar(UNIPROT)) %>%
        dplyr::transmute(input = toupper(SYMBOL), primaryAccession = toupper(UNIPROT)) %>%
        dplyr::distinct(input, .keep_all = TRUE)
      need_idx2 <- which(is.na(resolved2$Resolved_UNIPROT) | !nzchar(resolved2$Resolved_UNIPROT))
      base_need <- toupper(resolved2$token_base[need_idx2])
      hit <- map_sym2up$primaryAccession[match(base_need, map_sym2up$input)]
      ok <- !is.na(hit) & nzchar(hit)
      ii <- need_idx2[ok]
      if (length(ii)) { resolved2$Resolved_UNIPROT[ii] <- hit[ok]; resolved2$strategy[ii] <- "orgdb_symbol_entrez_uniprot" }
    }
  }
}

# 5) UniProt gene_primary resolver (Mus musculus), batched with retry, prefer reviewed
need_idx <- which(is.na(resolved2$Resolved_UNIPROT) | !nzchar(resolved2$Resolved_UNIPROT))
need_ids <- unique(toupper(resolved2$token_base[need_idx]))
sym_left2 <- unique(need_ids[grepl("^[A-Z0-9\\-]{2,}$", need_ids)])

if (length(sym_left2)) {
  batch_vec <- split(sym_left2, ceiling(seq_along(sym_left2)/50))
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
        pick <- tbl %>% dplyr::group_by(gene_primary) %>% dplyr::arrange(dplyr::desc(.data$reviewed), accession, .by_group = TRUE) %>% dplyr::slice_head(n=1) %>% dplyr::ungroup() %>% dplyr::transmute(input = gene_primary, primaryAccession = accession)
        picks[[length(picks)+1]] <- pick
      }
    }
  }
  if (length(picks)) {
    map_gene <- dplyr::bind_rows(picks) %>% dplyr::distinct(input, .keep_all = TRUE)
    need_idx3 <- which(is.na(resolved2$Resolved_UNIPROT) | !nzchar(resolved2$Resolved_UNIPROT))
    base_need <- toupper(resolved2$token_base[need_idx3])
    hit <- map_gene$primaryAccession[match(base_need, map_gene$input)]
    ok <- !is.na(hit) & nzchar(hit)
    ii <- need_idx3[ok]
    if (length(ii)) { resolved2$Resolved_UNIPROT[ii] <- hit[ok]; resolved2$strategy[ii] <- "uniprot_gene_primary_retry" }
  }
}

# --------------------------
# Save unmapped lists and audits
# --------------------------
audit_mouse <- resolved2 %>%
  dplyr::mutate(mapped = !is.na(Resolved_UNIPROT) & nzchar(Resolved_UNIPROT)) %>%
  dplyr::count(id_class, strategy, mapped, name = "n") %>%
  dplyr::arrange(desc(n))
write_tsv_safe(audit_mouse, fp_maptab("mapping_audit_mouse_only_robust.tsv"))

unmapped_tokens <- resolved2 %>%
  dplyr::filter(is.na(Resolved_UNIPROT) | !nzchar(Resolved_UNIPROT)) %>%
  dplyr::transmute(
    .row_id, token_raw, token_up, token_base, id_class,
    reason = dplyr::case_when(
      grepl("[^A-Z0-9_\\-\\.]", token_base) ~ "illegal_chars",
      grepl("^[A-Z0-9\\-\\.]+$", token_base) ~ "entry_name_not_in_local_or_query",
      TRUE ~ "unexpected_format_or_na"
    )
  ) %>% dplyr::arrange(id_class, token_base)
write_tsv_safe(unmapped_tokens, fp_maptab("unmapped_mouse_tokens.tsv"))

unmapped_summary <- unmapped_tokens %>% dplyr::count(id_class, reason, token_base, name = "n") %>% dplyr::arrange(dplyr::desc(n))
write_tsv_safe(unmapped_summary, fp_maptab("unmapped_mouse_tokens_summary.tsv"))

# --------------------------
# Collapse to features and build expression matrix
# --------------------------
collapse_ids <- function(x) { x <- unique(x[!is.na(x) & nzchar(x)]); if (!length(x)) return(NA_character_); paste(x, collapse = ";") }
male.norm <- resolved2 %>%
  dplyr::group_by(.row_id) %>%
  dplyr::summarise(gene_symbol = collapse_ids(Resolved_UNIPROT), .groups = "drop") %>%
  dplyr::right_join(male.data %>% dplyr::select(-gene_symbol, .row_id), by = ".row_id") %>%
  dplyr::select(-.row_id) %>%
  dplyr::mutate(gene_symbol = dplyr::na_if(gene_symbol, ""))

to_numeric_matrix <- function(male_norm, qc_dir = subdirs$logs) {
  if (!"gene_symbol" %in% names(male_norm)) stop("male.norm must contain gene_symbol")
  expr <- as.data.frame(lapply(male_norm[, -1, drop = FALSE], function(x) suppressWarnings(as.numeric(x))))
  if (!all(vapply(expr, is.numeric, logical(1)))) stop("Non-numeric columns remain after coercion")
  mat <- as.data.frame(t(expr))
  feat <- male_norm$gene_symbol
  empty <- which(!nzchar(ifelse(is.na(feat), "", feat)))
  if (length(empty)) feat[empty] <- paste0("UNMAPPED_", seq_along(empty))
  feat <- make.unique(feat, sep = "_")
  colnames(mat) <- feat
  if (any(!nzchar(colnames(mat)) | is.na(colnames(mat)))) stop("Empty/NA feature names after repair")
  utils::write.table(utils::head(mat[, 1:min(10, ncol(mat)), drop = FALSE]), file.path(qc_dir, "expression_head10.tsv"), sep = "\t", row.names = TRUE, quote = FALSE)
  mat
}
expression.data <- to_numeric_matrix(male.norm)

# Force single accession per feature name and ensure uniqueness
fix_feature_ids <- function(nms) {
  first <- sub(";.*$", "", nms)
  first <- toupper(trimws(first))
  first[is.na(first) | !nzchar(first)] <- "UNMAPPED"
  make.unique(first, sep = "_")
}

colnames(expression.data) <- fix_feature_ids(colnames(expression.data))

# Save core outputs
write_tsv_safe(resolved2, fp_maptab("resolved_tokens_mouse_only_robust.tsv"))
saveRDS(list(expression = expression.data, male.norm = male.norm, mapping = resolved2),
        file = fp_state("mouse_only_mapping_outputs_robust.rds"))

# --------------------------
# QC and sample clustering
# --------------------------
gsg <- goodSamplesGenes(expression.data)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) printFlush(paste("Removing genes:", paste(names(expression.data)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples) > 0) printFlush(paste("Removing samples:", paste(rownames(expression.data)[!gsg$goodSamples], collapse = ", ")))
  expression.data <- expression.data[gsg$goodSamples == TRUE, gsg$goodGenes == TRUE]
}

sampleTree <- hclust(dist(expression.data), method = "average")
svg(file = fp_qc("sample_clustering_outliers.svg"), width = nature_single_col * 1.4, height = 3.8, family = nature_font)
par(cex = 0.6, mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.main = 2)
abline(h = sample_tree_plot_height, col = "red")
dev.off()

cut.sampleTree <- cutreeStatic(sampleTree, cutHeight = sample_tree_cut_height, minSize = 10)
expression.data <- expression.data[cut.sampleTree == 1, ]

# --------------------------
# Soft-threshold selection
# --------------------------
spt <- pickSoftThreshold(expression.data, networkType = "signed",
                         corFnc = "bicor",
                         corOptions = list(use = "p", maxPOutliers = 0.05))

svglite::svglite(file = fp_qc("soft_threshold_scale_independence.svg"), width = nature_single_col, height = 2.8)
par(mar = c(4,4,2,1))
plot(spt$fitIndices[,1], spt$fitIndices[,2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n", main = "Scale independence")
text(spt$fitIndices[,1], spt$fitIndices[,2], labels = spt$fitIndices[,1], col = "red")
abline(h = soft_threshold_rsquared, col = "red")
dev.off()

svglite::svglite(file = fp_qc("soft_threshold_mean_connectivity.svg"), width = nature_single_col, height = 2.8)
par(mar = c(4,4,2,1))
plot(spt$fitIndices[,1], spt$fitIndices[,5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = "Mean connectivity")
text(spt$fitIndices[,1], spt$fitIndices[,5], labels = spt$fitIndices[,1], col = "red")
dev.off()

fit_indices <- tibble::as_tibble(spt$fitIndices)
eligible_power <- fit_indices$Power[fit_indices$SFT.R.sq >= soft_threshold_rsquared]
softPower <- if (length(eligible_power)) {
  min(eligible_power)
} else {
  fit_indices$Power[which.max(fit_indices$SFT.R.sq)]
}
write_csv_safe(fit_indices, fp_qctab("soft_threshold_fit_indices.csv"))
write_csv_safe(
  tibble::tibble(
    parameter = c("soft_threshold_rsquared", "selected_soft_power"),
    value = c(soft_threshold_rsquared, softPower)
  ),
  fp_qctab("soft_threshold_selection.csv")
)

# --------------------------
# Network construction
# --------------------------
adjacency <- adjacency(expression.data, power = softPower, type = "signed",
                       corFnc = "bicor", corOptions = list(use="p", maxPOutliers=0.05))
TOM <- TOMsimilarity(adjacency)
TOM.dissimilarity <- 1 - TOM
geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average")

svg(file = fp_net("gene_dendrogram.svg"), width = nature_double_col, height = 5.2, family = nature_font)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

Modules <- cutreeDynamic(dendro = geneTree, distM = TOM.dissimilarity, deepSplit = deep_split,
                         pamRespectsDendro = FALSE, minClusterSize = min_module_size)

colorSeq <- c(
  "lemon" = "lemonchiffon", "sage" = "darkseagreen", "bluegray" = "steelblue",
  "mintblue" = "lightsteelblue", "azure" = "deepskyblue", "khaki" = "khaki",
  "skyblue" = "skyblue", "babyblue" = "lightblue", "amber" = "goldenrod",
  "tealgreen" = "darkcyan", "forestgreen" = "forestgreen", "gold" = "gold",
  "violet" = "violet", "seafoam" = "mediumaquamarine", "coral" = "coral",
  "salmonlight" = "lightsalmon", "peach" = "peachpuff", "mint" = "palegreen",
  "lime" = "limegreen", "mauve" = "plum", "freesia" = "lightpink", "cocoa" = "saddlebrown", 
  "lavender" = "lavender", "magenta" = "magenta", "salmon" = "salmon",
  "rose" = "mistyrose", "aquamarine" = "aquamarine", "tomato" = "tomato",
  "plum" = "plum", "hotpink" = "hotpink", "rust" = "sienna"
)

# Map numeric module labels -> colors from colorSeq (recycle palette if needed)
unique_mods <- sort(unique(Modules))
nmods <- length(unique_mods)
palette_vals <- unname(colorSeq)
if (length(palette_vals) < nmods) palette_vals <- rep(palette_vals, length.out = nmods)
mod_colors_map <- setNames(palette_vals[seq_len(nmods)], as.character(unique_mods))
ModuleColors <- as.character(mod_colors_map[as.character(Modules)])
stopifnot(length(ModuleColors) == length(Modules))

svg(file = fp_net("gene_dendrogram_module_colors.svg"), width = nature_double_col, height = 5.2, family = nature_font)
plotDendroAndColors(geneTree, ModuleColors, "Module", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

MElist <- moduleEigengenes(expression.data, colors = ModuleColors)
MEs <- MElist$eigengenes
ME.dissimilarity <- 1 - cor(MEs, use = "p", method = "pearson")
METree <- hclust(as.dist(ME.dissimilarity), method = "average")
merge <- mergeCloseModules(expression.data, ModuleColors, cutHeight = merge_cut_height)
mergedColors <- merge$colors
names(mergedColors) <- colnames(expression.data)
mergedMEs <- orderMEs(merge$newMEs)

write_csv_safe(
  tibble::tibble(
    parameter = c(
      "sample_tree_cut_height", "sample_tree_plot_height",
      "soft_threshold_rsquared", "selected_soft_power",
      "network_type", "correlation_function", "bicor_maxPOutliers",
      "deep_split", "min_module_size", "merge_cut_height",
      "module_preservation_permutations"
    ),
    value = c(
      sample_tree_cut_height, sample_tree_plot_height,
      soft_threshold_rsquared, softPower,
      "signed", "bicor", 0.05,
      deep_split, min_module_size, merge_cut_height,
      module_preservation_permutations
    )
  ),
  fp_log("analysis_parameters.csv")
)

svg(file = fp_net("gene_dendrogram_modules_merged.svg"), width = nature_double_col, height = 5.2, family = nature_font)
plotDendroAndColors(geneTree, cbind(ModuleColors, mergedColors),
                    c("Original Module","Merged Module"),
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors for original and merged modules")
dev.off()

# --------------------------
# Eigengene network plots
# --------------------------
MET <- orderMEs(mergedMEs)
svg(file = fp_net("eigengene_dendrogram.svg"), width = nature_single_col, height = 3.0, family = nature_font)
plotEigengeneNetworks(MET, "", plotHeatmaps = FALSE, marDendro = c(0, 4, 2, 0))
dev.off()
svg(file = fp_net("eigengene_adjacency_heatmap.svg"), width = nature_single_col, height = 3.0, family = nature_font)
par(mar = c(1,1,1,1))
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", plotDendrograms = FALSE,
                      marHeatmap = c(5, 5, 2, 2), xLabelsAngle = 90)
dev.off()

# --------------------------
# Spatial + condition trait matrix and heatmaps
# --------------------------
sample_info <- readxl::read_excel(path = meta_xlsx)
stopifnot("row.names" %in% names(sample_info))
rownames(sample_info) <- as.character(sample_info$row.names)

Samples <- rownames(expression.data)
sample_info <- sample_info[Samples, , drop = FALSE]
profile_info <- prepare_spatial_metadata(
  sample_info,
  infer_dataset_profile(sample_info, expr_xlsx, meta_xlsx, dataset_profile)
)
sample_info <- profile_info$sample_info
dataset_profile_resolved <- profile_info$profile
active_spatial_vars <- profile_info$active_spatial_vars

make_trait_matrix <- function(sample_info, vars) {
  mats <- lapply(vars, function(v) {
    x <- droplevels(factor(sample_info[[v]]))
    if (length(unique(stats::na.omit(x))) <= 1) return(NULL)
    mm <- model.matrix(~ 0 + x)
    colnames(mm) <- paste0(v, "_", sub("^x", "", colnames(mm)))
    mm
  })
  mats <- Filter(Negate(is.null), mats)
  if (!length(mats)) stop("No variable traits remain after filtering")
  as.data.frame(do.call(cbind, mats), stringsAsFactors = FALSE)
}

datTraits <- make_trait_matrix(sample_info, c(active_spatial_vars, "condition"))
keep_cols <- vapply(datTraits, function(x) sd(as.numeric(x), na.rm = TRUE) > 0, logical(1))
datTraits <- datTraits[, keep_cols, drop = FALSE]
rownames(datTraits) <- Samples

nSamples <- nrow(expression.data)
MEcorr <- cor(mergedMEs, datTraits, use = "p", method = "pearson")
MEp    <- corPvalueStudent(MEcorr, nSamples)
MEfdr  <- matrix(
  p.adjust(as.vector(MEp), method = "BH"),
  nrow = nrow(MEp), ncol = ncol(MEp), dimnames = dimnames(MEp)
)

plot_trait_heatmap <- function(matCorr, matFdr, cols, file) {
  if (length(cols) == 0) return(invisible(NULL))
  textMatrix <- paste(signif(matCorr[, cols, drop=FALSE], 2), "\nq=",
                      signif(matFdr[, cols, drop=FALSE], 1), sep = "")
  dim(textMatrix) <- dim(matCorr[, cols, drop=FALSE])
  svg(file = file, width = 6, height = max(4, nrow(matCorr) * 0.25 + 2))
  par(mar = c(6, 8.5, 3, 1))
  labeledHeatmap(Matrix = matCorr[, cols, drop=FALSE],
                 xLabels = colnames(matCorr)[cols],
                 yLabels = rownames(matCorr),
                 ySymbols = rownames(matCorr),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.6,
                 zlim = c(-1, 1),
                 main = "Module–trait relationships")
  dev.off()
}

trait_names <- colnames(datTraits)
groups <- list(
  celltype = grep("^celltype_", trait_names),
  layer    = grep("^layer_",   trait_names),
  region   = grep("^region_",  trait_names),
  condition = grep("^condition_", trait_names)
)
groups <- groups[c(intersect(c("celltype", "region", "layer"), active_spatial_vars), "condition")]
for (nm in names(groups)) {
  idx <- groups[[nm]]
  if (length(idx) > 0) {
    plot_trait_heatmap(MEcorr, MEfdr, idx, fp_traits(paste0("ME_trait_heatmap_", nm, ".svg")))
  }
}
write_csv_safe(
  reshape2::melt(MEcorr, varnames = c("module", "trait"), value.name = "r") |>
    dplyr::mutate(
      p = as.vector(MEp),
      fdr = as.vector(MEfdr)
    ),
  fp_source("ME_trait_correlations.csv")
)

# --------------------------
# Pairwise condition contrasts (optional)
# --------------------------
mk_contrast <- function(vec, a, b){v<-rep(NA_real_,length(vec));v[vec==a]<-0;v[vec==b]<-1;v}
grp <- as.character(sample_info$ExpGroup)
contrasts <- list(con_res = mk_contrast(grp, "con", "res"),
                  con_sus = mk_contrast(grp, "con", "sus"),
                  res_sus = mk_contrast(grp, "res", "sus"))
for (nm in names(contrasts)) {
  v <- contrasts[[nm]]; keep <- !is.na(v)
  cmat <- cor(mergedMEs[keep, , drop=FALSE], v[keep], use="p")
  pmat <- corPvalueStudent(cmat, sum(keep))
  qmat <- matrix(p.adjust(as.vector(pmat), method = "BH"),
                 nrow = nrow(pmat), ncol = ncol(pmat), dimnames = dimnames(pmat))
  txt <- paste(signif(cmat,2), "\nq=", signif(qmat,1), sep = "")
  dim(txt) <- dim(cmat)
  svg(file = fp_traits(paste0("ME_trait_heatmap_", nm, ".svg")),
      width = 3, height = max(4, ncol(mergedMEs) * 0.25 + 2))
  par(mar = c(6, 8.5, 3, 1))
  labeledHeatmap(Matrix = cmat, xLabels = nm,
                 yLabels = colnames(mergedMEs),
                 ySymbols = colnames(mergedMEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = txt,
                 setStdMargins = FALSE,
                 cex.text = 0.8, zlim = c(-1,1),
                 main = "Module–trait relationships")
  dev.off()
}

# --------------------------
# kME, GS, hubs
# --------------------------
modNames <- substring(colnames(mergedMEs), 3)
geneModuleMembership <- as.data.frame(cor(expression.data,
                                          mergedMEs, use = "p",
                                          method = "pearson"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste0("MM", modNames)
names(MMPvalue) <- paste0("p.MM", modNames)

condition_factor <- factor(sample_info$ExpGroup, levels = c("con", "res", "sus"))
gene_condition_p <- vapply(expression.data, function(x) {
  fit <- stats::lm(x ~ condition_factor)
  stats::anova(fit)[["Pr(>F)"]][1]
}, numeric(1))
gene_condition_fdr <- p.adjust(gene_condition_p, method = "BH")
geneTraitSignificance <- data.frame(
  GS.ExpGroup = -log10(pmax(gene_condition_fdr, .Machine$double.xmin)),
  row.names = colnames(expression.data)
)
GSPvalue <- data.frame(
  p.GS.ExpGroup = gene_condition_p,
  FDR.GS.ExpGroup = gene_condition_fdr,
  row.names = colnames(expression.data)
)
names(geneTraitSignificance) <- "GS.ExpGroup"

modules_of_interest <- unique(mergedColors)
for (module in modules_of_interest) {
  moduleGenes <- mergedColors == module
  if (!any(moduleGenes)) next
  mmcol <- paste0("MM", module); if (!mmcol %in% colnames(geneModuleMembership)) next
  gene_info <- data.frame(
    Gene = colnames(expression.data)[moduleGenes],
    Module = module,
    ModuleMembership = geneModuleMembership[moduleGenes, mmcol],
    GeneSignificance = geneTraitSignificance[moduleGenes, "GS.ExpGroup"],
    GeneSignificanceP = GSPvalue[moduleGenes, "p.GS.ExpGroup"],
    GeneSignificanceFDR = GSPvalue[moduleGenes, "FDR.GS.ExpGroup"]
  )
  write_csv_safe(gene_info, fp_modtab(paste0("genes_in_module_", module, ".csv")))
}

# ==========================
# Module naming via GO terms (robust UNIPROT→ENTREZ + consistent renaming)
# ==========================
# Feature IDs may contain multiple UniProt accessions separated by ';'
all_feature_ids <- colnames(expression.data)
feature_to_acc <- strsplit(all_feature_ids, ";", fixed = TRUE)
feature_to_acc <- lapply(feature_to_acc, function(v){
  v <- toupper(trimws(v))
  # Keep only accession-like IDs (incl. A0A…)
  v[grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$|^A0A[0-9A-Z]{7}$", v)]
})
names(feature_to_acc) <- all_feature_ids

# Accession universe from all features used in WGCNA
acc_universe <- unique(unlist(feature_to_acc, use.names = FALSE))

# Map UNIPROT accessions → ENTREZ
if (length(acc_universe)) {
  map_df <- suppressWarnings(clusterProfiler::bitr(
    acc_universe,
    fromType = "UNIPROT",
    toType   = "ENTREZID",
    OrgDb    = org.Mm.eg.db
  ))
  map_df <- dplyr::distinct(map_df, UNIPROT, .keep_all = TRUE)
} else {
  map_df <- tibble::tibble(UNIPROT = character(), ENTREZID = character())
}
sym_to_entrez <- setNames(map_df$ENTREZID, map_df$UNIPROT)

# Universe of ENTREZ IDs present in the network (character)
universe_entrez <- unique(na.omit(unname(sym_to_entrez[acc_universe])))
universe_entrez <- as.character(universe_entrez)

# Helper: compact a GO term description
compact_term <- function(term) {
  term <- gsub("\\b(process|regulation|pathway|of|the|cellular)\\b", "", term, ignore.case = TRUE)
  term <- gsub("[[:punct:]]+", " ", term)
  term <- str_squish(term)
  str_to_title(term)
}

# Simple redundancy pruning by Jaccard on tokens
prune_terms <- function(df, top_n = 5) {
  if (nrow(df) <= 1) return(utils::head(df, top_n))
  keep <- rep(TRUE, nrow(df))
  names_list <- stringr::str_split(tolower(df$Description), " ")
  for (i in seq_len(nrow(df))) {
    if (!keep[i]) next
    for (j in seq.int(i + 1, nrow(df))) {
      if (!keep[j]) next
      s1 <- names_list[[i]]; s2 <- names_list[[j]]
      inter <- length(intersect(s1, s2))
      uni   <- length(union(s1, s2))
      jac   <- if (uni == 0) 0 else inter / uni
      if (jac >= 0.60) keep[j] <- FALSE
    }
  }
  utils::head(df[keep, , drop = FALSE], top_n)
}

# Per-module enrichment
enrich_one_module <- function(mod_color, universe_entrez,
                              min_mod_n = 10, min_univ_n = 100,
                              pcut = 0.1, qcut = 0.1) {

  # Features in module (by color)
  features_in_mod <- names(mergedColors)[mergedColors == mod_color]
  if (!length(features_in_mod)) return(NULL)

  # Accessions for those features
  acc_in_mod <- unique(unlist(feature_to_acc[features_in_mod], use.names = FALSE))

  # ENTREZ IDs for the module
  entrez_in_mod <- unique(na.omit(unname(sym_to_entrez[acc_in_mod])))
  entrez_in_mod <- as.character(entrez_in_mod)

  # Guard rails
  if (length(universe_entrez) < min_univ_n || length(entrez_in_mod) < min_mod_n) return(NULL)

  # Enrichment
  ego <- suppressWarnings(clusterProfiler::enrichGO(
    gene          = entrez_in_mod,
    universe      = universe_entrez,
    OrgDb         = org.Mm.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = pcut,
    qvalueCutoff  = qcut,
    readable      = TRUE
  ))
  if (is.null(ego)) return(NULL)
  df <- as.data.frame(ego)
  if (!nrow(df)) return(NULL)

  # Score = -log10(qvalue) * GeneRatio
  gr <- vapply(df$GeneRatio, function(x){
    xy <- strsplit(x, "/", fixed = TRUE)[[1]]
    as.numeric(xy[1]) / as.numeric(xy[2])
  }, numeric(1))
  df$score <- (-log10(pmax(df$qvalue, 1e-300))) * gr

  df <- df |>
    dplyr::arrange(dplyr::desc(score)) |>
    dplyr::filter(qvalue <= qcut)

  if (!nrow(df)) return(NULL)

  df_top <- prune_terms(df, top_n = 5)
  df_top$Compact <- vapply(df_top$Description, compact_term, character(1))
  df_top$Module  <- mod_color
  df_top
}

# Run enrichment over module colors
modules_vec <- unique(mergedColors)
annot_list <- lapply(modules_vec, enrich_one_module, universe_entrez = universe_entrez)
annot_df <- dplyr::bind_rows(annot_list)

# Fallback if no terms
mods_no_hits <- setdiff(modules_vec, unique(annot_df$Module))
if (length(mods_no_hits)) {
  fallback <- data.frame(
    Module      = mods_no_hits,
    ID          = NA_character_,
    Description = NA_character_,
    qvalue      = NA_real_,
    GeneRatio   = NA_character_,
    score       = NA_real_,
    Compact     = NA_character_,
    stringsAsFactors = FALSE
  )
  annot_df <- dplyr::bind_rows(annot_df, fallback)
}

# Best name per module color
best_by_mod <- annot_df |>
  dplyr::group_by(Module) |>
  dplyr::slice_max(order_by = score, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::transmute(
    Module,
    BestTerm    = ifelse(is.na(Description), paste0("module_", Module), Description),
    BestCompact = ifelse(is.na(Compact),    paste0("Module ", Module), Compact),
    qvalue,
    score
  )

# Persist annotations
write_csv_safe(annot_df |> dplyr::arrange(Module, dplyr::desc(score)),
               fp_modtab("module_annotations_GO.csv"))
write_csv_safe(best_by_mod, fp_modtab("module_name_map.csv"))

# Color -> Name map
module_name_map <- setNames(best_by_mod$BestCompact, best_by_mod$Module)

# Rename ME columns from "ME<color>" to "ME<Compact>"
ME_names_old <- colnames(mergedMEs)
ME_colors    <- sub("^ME", "", ME_names_old)
ME_new_suffix <- ifelse(ME_colors %in% names(module_name_map),
                        make.names(module_name_map[ME_colors]),
                        ME_colors)
ME_names_new <- paste0("ME", ME_new_suffix)
colnames(mergedMEs) <- ME_names_new

# Provide a robust color->MEcol lookup for downstream MM/GS code
color_to_MEcol <- setNames(ME_names_new, ME_colors)

# Save gene→module with color and compact name
gene_module_named <- data.frame(
  Gene          = colnames(expression.data),
  ModuleColor   = mergedColors[match(colnames(expression.data), names(mergedColors))],
  ModuleName    = unname(module_name_map[mergedColors[match(colnames(expression.data), names(mergedColors))]]),
  stringsAsFactors = FALSE
)
write_csv_safe(gene_module_named, fp_modtab("genes_module_with_names.csv"))

# Save maps
module_name_map_df <- data.frame(
  ModuleColor = names(module_name_map),
  ModuleName  = unname(module_name_map),
  stringsAsFactors = FALSE
)
write_tsv_safe(module_name_map_df, fp_modtab("module_name_map.tsv"))
me_rename_df <- data.frame(
  ME_old = ME_names_old,
  ME_new = ME_names_new,
  stringsAsFactors = FALSE
)
write_tsv_safe(me_rename_df, fp_modtab("ME_column_rename_map.tsv"))

saveRDS(list(
  module_name_map = module_name_map,
  color_to_MEcol  = color_to_MEcol,
  ME_names_old    = ME_names_old,
  ME_names_new    = ME_names_new
), file = fp_state("module_name_map.rds"))

# --------------------------
# Preservation across conditions
# --------------------------
idx_con <- which(sample_info$ExpGroup == "con")
idx_res <- which(sample_info$ExpGroup == "res")
idx_sus <- which(sample_info$ExpGroup == "sus")

multiExpr <- list(
  ALL = list(data = expression.data),
  CON = list(data = expression.data[idx_con, , drop = FALSE]),
  RES = list(data = expression.data[idx_res, , drop = FALSE]),
  SUS = list(data = expression.data[idx_sus, , drop = FALSE])
)

good_cols <- lapply(multiExpr, function(e) {
  gsg <- goodSamplesGenes(e$data, verbose = 3)
  which(gsg$goodGenes)
})
common_idx <- Reduce(intersect, good_cols)
common_genes <- colnames(multiExpr[[1]]$data)[common_idx]

has_na <- sapply(common_genes, function(g) any(vapply(multiExpr, function(e) any(is.na(e$data[, g])), logical(1))))
common_genes <- common_genes[!has_na]
stopifnot(length(common_genes) > 0)

multi_expr_clean <- lapply(multiExpr, function(e) {
  dat <- e$data[, common_genes, drop = FALSE]
  dat <- as.data.frame(dat)
  dat[] <- lapply(dat, function(col) as.numeric(col))
  keep_samp <- apply(dat, 1, function(r) sd(r, na.rm = TRUE) > 0)
  dat <- dat[keep_samp, , drop = FALSE]
  list(data = as.matrix(dat))
})
names(multi_expr_clean) <- names(multiExpr)

if (any(vapply(multi_expr_clean, function(e) nrow(e$data) < 2, logical(1)))) {
  stop("One or more sets have < 2 samples after cleaning; reduce filtering or combine sets.")
}

ref_colors <- mergedColors[match(common_genes, colnames(expression.data))]
stopifnot(length(ref_colors) == length(common_genes))

multi_color <- list(
  ALL = ref_colors,
  CON = rep("grey", length(common_genes)),
  RES = rep("grey", length(common_genes)),
  SUS = rep("grey", length(common_genes))
)

set.seed(12345)
mp <- modulePreservation(
  multi_expr_clean,
  multi_color,
  referenceNetworks = 1,
  nPermutations = module_preservation_permutations,
  networkType = "signed",
  corFnc = "bicor",
  corOptions = "use = 'p', maxPOutliers = 0.05",
  parallelCalculation = TRUE,  # when supported by your WGCNA build
  verbose = 3
)

Z_CON <- mp$preservation$Z$ref.ALL$inColumnsAlsoPresentIn.CON
Z_RES <- mp$preservation$Z$ref.ALL$inColumnsAlsoPresentIn.RES
Z_SUS <- mp$preservation$Z$ref.ALL$inColumnsAlsoPresentIn.SUS

write_csv_safe(tibble::rownames_to_column(as.data.frame(Z_CON), "module"),
               fp_prestab("module_preservation_CON_vs_ALL_Zsummary.csv"))
write_csv_safe(tibble::rownames_to_column(as.data.frame(Z_RES), "module"),
               fp_prestab("module_preservation_RES_vs_ALL_Zsummary.csv"))
write_csv_safe(tibble::rownames_to_column(as.data.frame(Z_SUS), "module"),
               fp_prestab("module_preservation_SUS_vs_ALL_Zsummary.csv"))
write_csv_safe(
  dplyr::bind_rows(
    tibble::rownames_to_column(as.data.frame(Z_CON), "module") |> dplyr::mutate(test_set = "CON"),
    tibble::rownames_to_column(as.data.frame(Z_RES), "module") |> dplyr::mutate(test_set = "RES"),
    tibble::rownames_to_column(as.data.frame(Z_SUS), "module") |> dplyr::mutate(test_set = "SUS")
  ),
  fp_source("module_preservation.csv")
)

# --------------------------
# Condition and active spatial trait panels
# --------------------------
combo_vars <- c("condition", active_spatial_vars)
combo_df <- sample_info[, combo_vars, drop = FALSE] |>
  tibble::rownames_to_column("Sample") |>
  dplyr::mutate(dplyr::across(dplyr::all_of(combo_vars), as.character))
combo <- factor(do.call(paste, c(combo_df[combo_vars], sep = "_")))

# One-hot and correlations
X_combo <- model.matrix(~ 0 + combo)
colnames(X_combo) <- levels(combo)  # no "comb" prefix at all
if (ncol(X_combo) == 0) stop("No combined strata present in X_combo (ncol == 0). Check combo label construction.")

MEcorr_combo <- cor(mergedMEs, X_combo, use = "p")
MEp_combo    <- corPvalueStudent(MEcorr_combo, nrow(expression.data))

df_combo <- reshape2::melt(MEcorr_combo, varnames = c("module","comb"), value.name = "r")
p_combo  <- reshape2::melt(MEp_combo,    varnames = c("module","comb"), value.name = "p")
df_combo$p <- p_combo$p
df_combo$fdr <- p.adjust(df_combo$p, method = "BH")
df_combo$sig <- sig_dot(df_combo$fdr)
df_combo$module <- factor(df_combo$module, levels = rownames(MEcorr_combo))

split_combo_labels <- function(labels, vars) {
  parts <- strsplit(labels, "_", fixed = TRUE)
  do.call(rbind, lapply(parts, function(p) {
    p <- as.character(p)
    if (length(p) < length(vars)) p <- c(p, rep("missing", length(vars) - length(p)))
    p[seq_along(vars)]
  }))
}
lab_parts <- split_combo_labels(as.character(df_combo$comb), combo_vars)
colnames(lab_parts) <- combo_vars
df_combo2 <- cbind(df_combo[, c("module","comb","r","p","fdr","sig")], as.data.frame(lab_parts, stringsAsFactors = FALSE))
df_combo2$spatial_trait <- do.call(paste, c(df_combo2[active_spatial_vars], sep = "_"))
df_combo2$trait <- do.call(paste, c(df_combo2[c("condition", active_spatial_vars)], sep = "_"))
write_csv_safe(df_combo2, fp_source("ME_condition_spatial_strata.csv"))

present_conds <- intersect(c("con","res","sus"), unique(df_combo2$condition))
by_cond <- split(df_combo2, df_combo2$condition)

panel_plot <- function(dfi, panel_title = "") {
  if (!nrow(dfi)) return(NULL)
  ord <- do.call(order, dfi[active_spatial_vars])
  dfi$spatial_trait <- factor(dfi$spatial_trait, levels = unique(dfi$spatial_trait[ord]))
  ggplot(dfi, aes(x = spatial_trait, y = module, fill = r)) +
    geom_tile(color = "white", linewidth = 0.15) +
    geom_text(aes(label = sig), size = 1.8, color = "black", na.rm = TRUE) +
    scale_fill_gradient2(limits = c(-1, 1), oob = scales::squish,
                         low = nature_diverging["low"], mid = nature_diverging["mid"], high = nature_diverging["high"]) +
    labs(title = panel_title, x = paste(active_spatial_vars, collapse = " / "), y = NULL, fill = "Pearson r") +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none")
}

plots <- lapply(present_conds, function(cn) panel_plot(by_cond[[cn]], toupper(cn)))
plots <- Filter(function(p) !is.null(p) && inherits(p, "ggplot"), plots)
if (length(plots) == 0) stop("No non-empty condition panels; check combo labels and present conditions.")

# Combine panels; if single, save directly
if (length(plots) == 1) {
  g_combined <- plots[[1]]
  save_plot_nature(g_combined, fp_traits("panel_ME_vs_condition_spatial_strata.svg"),
                   width = nature_double_col, height = 3.9)
} else {
  g_combined <- patchwork::wrap_plots(plots, nrow = 1)
  legend_df <- data.frame(x = 1:3, y = 1:3, r = c(-1, 0, 1))
  p_legend <- ggplot(legend_df, aes(x, y, fill = r)) +
    geom_tile() +
    scale_fill_gradient2(limits = c(-1, 1), oob = scales::squish,
                         low = nature_diverging["low"], mid = nature_diverging["mid"], high = nature_diverging["high"]) +
    theme_void() + theme(legend.position = "right") + labs(fill = "Pearson r")
  legend_only <- cowplot::get_legend(p_legend)
  g <- cowplot::plot_grid(g_combined, legend_only, rel_widths = c(1, 0.08))
  save_plot_nature(g, fp_traits("panel_ME_vs_condition_spatial_strata.svg"),
                   width = nature_double_col, height = 3.9)
}

# --------------------------
# Ensure required objects are defined for the pheatmap block
# --------------------------
if (!exists("df_all")) {
    if (exists("df_combo2")) {
        df_all <- df_combo2
        if (!"trait" %in% names(df_all)) {
            df_all$trait <- do.call(paste, c(df_all[c("condition", active_spatial_vars)], sep = "_"))
        }
    } else if (exists("MEcorr_combo") && exists("MEp_combo")) {
        df_all <- reshape2::melt(MEcorr_combo, varnames = c("module", "comb"), value.name = "r")
        p_tmp  <- reshape2::melt(MEp_combo,    varnames = c("module", "comb"), value.name = "p")
        df_all$p <- p_tmp$p
        parts <- split_combo_labels(as.character(df_all$comb), combo_vars)
        colnames(parts) <- combo_vars
        df_all <- cbind(df_all, as.data.frame(parts, stringsAsFactors = FALSE))
        df_all$trait <- do.call(paste, c(df_all[c("condition", active_spatial_vars)], sep = "_"))
    } else {
        stop("df_all is missing and cannot be reconstructed: provide df_combo2 or MEcorr_combo/MEp_combo.")
    }
}

# module_levels (ordered modules)
if (!exists("module_levels")) {
    if (exists("MEcorr_combo")) {
        module_levels <- rownames(MEcorr_combo)
    } else {
        module_levels <- unique(df_all$module)
    }
}

# traits_in_order (ordered trait columns)
if (!exists("traits_in_order")) {
    if (exists("traits_in_order") && length(traits_in_order) > 0) {
        # already provided
    } else {
        df_all$trait <- do.call(paste, c(df_all[c("condition", active_spatial_vars)], sep = "_"))
        ord <- do.call(order, df_all[c("condition", active_spatial_vars)])
        traits_in_order <- unique(df_all$trait[ord])
    }
}

# strip_df: module -> color mapping
if (!exists("strip_df")) {
    if (exists("mergedColors")) {
        # mergedColors is gene-level mapping; fall back to simple default
        strip_df <- data.frame(module = module_levels, mod_col = rep("grey80", length(module_levels)), stringsAsFactors = FALSE)
    } else {
        strip_df <- data.frame(module = module_levels, mod_col = rep("grey80", length(module_levels)), stringsAsFactors = FALSE)
    }
}

# 1) Wide matrices and order
r_mat  <- reshape2::acast(df_all, module ~ trait, value.var = "r")
p_mat  <- reshape2::acast(df_all, module ~ trait, value.var = "p")

# FDR matrix
fdr_vec <- p.adjust(as.vector(p_mat), method = "BH")
fdr_mat <- matrix(fdr_vec, nrow = nrow(p_mat), ncol = ncol(p_mat), dimnames = dimnames(p_mat))

stopifnot(all(module_levels %in% rownames(r_mat)), all(traits_in_order %in% colnames(r_mat)))
r_mat   <- r_mat[module_levels, traits_in_order, drop = FALSE]
p_mat   <- p_mat[module_levels, traits_in_order, drop = FALSE]
fdr_mat <- fdr_mat[module_levels, traits_in_order, drop = FALSE]

# 2) Column clustering and separators
hc_cols <- hclust(as.dist(1 - cor(r_mat, use = "pairwise.complete.obs")), method = "average")
k_clusters <- 6
col_grp <- cutree(hc_cols, k = k_clusters)
grp_ord <- col_grp[hc_cols$order]
gap_pos <- which(grp_ord[-1] != head(grp_ord, -1))
xlines <- gap_pos + 0.5

# 3) Module color mapping for row strip
if (exists("mergedMEs")) {
  me_names <- colnames(mergedMEs)                      # e.g., "MEblue"
  row_me <- rownames(r_mat)
  row_mod_colors <- if (exists("ME_names_new") && exists("ME_colors")) {
    unname(setNames(ME_colors, ME_names_new)[row_me])
  } else {
    sub("^ME", "", row_me)
  }
  row_mod_colors[is.na(row_mod_colors) | !nzchar(row_mod_colors)] <- sub("^ME", "", row_me[is.na(row_mod_colors) | !nzchar(row_mod_colors)])
  if (!exists("colorSeq")) {
    colorSeq <- c(
      "turquoise"="turquoise","blue"="blue","brown"="brown","yellow"="yellow",
      "green"="green","red"="red","black"="black","pink"="pink","magenta"="magenta",
      "purple"="purple","greenyellow"="greenyellow","tan"="tan","salmon"="salmon",
      "cyan"="cyan","midnightblue"="midnightblue","lightcyan"="lightcyan",
      "greenyellow2"="#ADFF2F","royalblue"="royalblue","darkred"="darkred",
      "skyblue"="skyblue","orange"="orange","grey"="grey"
    )
  }
  color_lookup <- colorSeq
  missing_cols <- setdiff(unique(row_mod_colors), names(color_lookup))
  if (length(missing_cols)) {
    add <- setNames(missing_cols, missing_cols)
    color_lookup <- c(color_lookup, add)
  }
  row_colors_vec <- unname(color_lookup[row_mod_colors])
  row_colors_vec[is.na(row_colors_vec)] <- "grey80"
  ann_row <- data.frame(ModuleColor = row_colors_vec, row.names = rownames(r_mat))
  ann_colors <- list(ModuleColor = setNames(unique(ann_row$ModuleColor), unique(ann_row$ModuleColor)))
} else {
  ann_row <- data.frame(ModuleColor = rep("grey80", nrow(r_mat)), row.names = rownames(r_mat))
  ann_colors <- list(ModuleColor = c("grey80" = "grey80"))
}

# Column condition annotation
col_condition <- sapply(strsplit(colnames(r_mat), "_", fixed = TRUE), `[`, 1)
ann_col <- data.frame(Condition = col_condition, row.names = colnames(r_mat))
ann_colors$Condition <- c(con = "#E6E6E6", res = "#D9ECF7", sus = "#F6D9CC")

# 4) Palette and breaks
r_lim <- 0.8
bk <- seq(-r_lim, r_lim, length.out = 201)
pal <- colorRampPalette(c(nature_diverging["low"], nature_diverging["mid"], nature_diverging["high"]))(200)

# 5) Significance dots
display_mat <- matrix("", nrow = nrow(r_mat), ncol = ncol(r_mat), dimnames = dimnames(r_mat))
display_mat[fdr_mat < 0.01] <- "•"

# 6) Render pheatmap once and save
ph <- pheatmap::pheatmap(
  mat = r_mat,
  cluster_rows = FALSE,
  cluster_cols = hc_cols,
  treeheight_row = 0,
  treeheight_col = 50,
  gaps_col = NULL,
  annotation_row = ann_row,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  color = pal, breaks = bk, legend = TRUE,
  border_color = NA,
  cellwidth = 8.5, cellheight = 8.5,
  angle_col = 45,
  fontsize_col = 6.2, fontsize_row = 6.2,
  display_numbers = display_mat,
  number_color = "grey20",
  number_cex = 0.6,
  show_rownames = TRUE, show_colnames = TRUE,
  silent = TRUE
)

nr <- nrow(r_mat)

png(fp_traits("ME_trait_pheatmap.png"), width = round(nature_double_col * 300), height = round(4.3 * 300), res = 300)
grid::grid.newpage(); grid::grid.draw(ph$gtable)
panel_id <- grep("matrix", ph$gtable$layout$name)[1]
seekViewport(ph$gtable$layout$name[panel_id])
for (xl in xlines) {
  grid::grid.lines(x = unit(c(xl, xl), "native"),
                   y = unit(c(0, nr), "native"),
                   gp = gpar(col = "white", lwd = 3))
}
upViewport(0); dev.off()

pdf(fp_traits("ME_trait_pheatmap.pdf"), width = nature_double_col, height = 4.3, family = nature_font, useDingbats = FALSE)
grid::grid.newpage(); grid::grid.draw(ph$gtable)
panel_id <- grep("matrix", ph$gtable$layout$name)[1]
seekViewport(ph$gtable$layout$name[panel_id])
for (xl in xlines) {
  grid::grid.lines(x = unit(c(xl, xl), "native"),
                   y = unit(c(0, nr), "native"),
                   gp = gpar(col = "white", lwd = 2))
}
upViewport(0); dev.off()

svglite::svglite(fp_traits("ME_trait_pheatmap.svg"), width = nature_double_col, height = 4.3)
grid::grid.newpage(); grid::grid.draw(ph$gtable)
panel_id <- grep("matrix", ph$gtable$layout$name)[1]
seekViewport(ph$gtable$layout$name[panel_id])
for (xl in xlines) {
  grid::grid.lines(x = unit(c(xl, xl), "native"),
                   y = unit(c(0, nr), "native"),
                   gp = gpar(col = "white", lwd = 2))
}
upViewport(0); dev.off()

block_list <- c(
  setNames(lapply(active_spatial_vars, function(v) grep(paste0("^", v, "_"), colnames(datTraits))), active_spatial_vars),
  list(condition = grep("^condition_", colnames(datTraits)))
)

block_cor_df <- function(block_idx, block_name) {
  if (length(block_idx) == 0) return(NULL)
  r <- cor(mergedMEs, datTraits[, block_idx, drop=FALSE], use="p")
  p <- corPvalueStudent(r, nrow(expression.data))
  df <- reshape2::melt(r, varnames = c("module","trait"), value.name = "r")
  p_long <- reshape2::melt(p, varnames = c("module","trait"), value.name = "p")
  df$p <- p_long$p
  df$fdr <- p.adjust(df$p, method = "BH")
  df$sig <- sig_dot(df$fdr)
  df$block <- block_name
  df$module <- factor(df$module, levels = rownames(r))
  df
}

dfs <- Filter(Negate(is.null),
              mapply(block_cor_df, block_list, names(block_list), SIMPLIFY = FALSE))
write_csv_safe(dplyr::bind_rows(dfs), fp_source("ME_trait_block_correlations.csv"))

panel_plot <- function(dfi, legend = "none") {
  if (nrow(dfi) == 0) return(NULL)
  dfi$trait <- factor(dfi$trait, levels = unique(dfi$trait))

  p_heat <- ggplot(dfi, aes(x = trait, y = module, fill = r)) +
    geom_tile(color = "white", linewidth = 0.15) +
    geom_text(aes(label = sig), size = 1.8, color = "black", na.rm = TRUE) +
    scale_fill_gradient2(limits = c(-1, 1), oob = scales::squish,
                         low = nature_diverging["low"], mid = nature_diverging["mid"], high = nature_diverging["high"]) +
    labs(x = NULL, y = NULL, fill = "Pearson r") +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks = element_blank(),
          legend.position = legend,
          plot.margin = margin(2, 2, 2, 2))

  if (!identical(legend, "none")) return(p_heat)

  modules <- levels(dfi$module)
  mod_cols <- sub("^ME", "", modules)
  strip_df <- data.frame(module = modules, mod_col = mod_cols, stringsAsFactors = FALSE)
  strip_df$module <- factor(strip_df$module, levels = modules)

  p_strip <- ggplot(strip_df, aes(x = 1, y = module, fill = mod_col)) +
    geom_tile() +
    scale_fill_identity() +
    theme_void() +
    theme(plot.margin = margin(2, 0, 2, 2))

  combined <- p_strip + p_heat + plot_layout(widths = c(0.05, 1))
  combined
}

plots <- lapply(dfs, panel_plot)
p_legend_heat <- panel_plot(dfs[[1]], legend = "right")
legend_only <- cowplot::get_legend(p_legend_heat)

combo <- wrap_plots(plots, nrow = 1, guides = "collect") +
  plot_annotation(title = paste("Module-trait relationships:", paste(c(active_spatial_vars, "condition"), collapse = ", ")))

main_trait_panel <- cowplot::plot_grid(combo, legend_only, rel_widths = c(1, 0.08))
save_plot_nature(main_trait_panel, fp_mainfig("panel_module_trait_relationships_spatial.svg"),
                 width = nature_double_col, height = 3.7)

# ==========================================================
# ME by condition plots with significance and exports
# ==========================================================

# Long-form eigengenes + metadata
stopifnot(nrow(mergedMEs) == nrow(sample_info))
ME_long <- mergedMEs %>%
  tibble::rownames_to_column("Sample") %>%
  mutate(
    condition = as.character(sample_info$ExpGroup),
    region = as.character(sample_info$region),
    layer = as.character(sample_info$layer),
    celltype = as.character(sample_info$celltype)
  ) %>%
  pivot_longer(cols = starts_with("ME"),
               names_to = "module", values_to = "ME")

# Ensure order con, res, sus
ME_long$condition <- factor(ME_long$condition, levels = c("con","res","sus"))
ME_long$region <- factor(ME_long$region)
ME_long$layer <- factor(ME_long$layer)
ME_long$celltype <- factor(ME_long$celltype)

# Stats
do_stats <- function(df) {
  covars <- active_spatial_vars
  covars <- covars[vapply(df[covars], function(x) length(unique(stats::na.omit(x))) > 1, logical(1))]
  rhs_full <- paste(c("condition", covars), collapse = " + ")
  rhs_reduced <- if (length(covars)) paste(covars, collapse = " + ") else "1"
  full_formula <- stats::as.formula(paste("ME ~", rhs_full))
  reduced_formula <- stats::as.formula(paste("ME ~", rhs_reduced))
  full_fit <- stats::lm(full_formula, data = df)
  reduced_fit <- stats::lm(reduced_formula, data = df)
  condition_p <- stats::anova(reduced_fit, full_fit)[["Pr(>F)"]][2]
  list(condition_p = condition_p, model_formula = deparse(full_formula))
}

stat_list <- ME_long %>%
  group_by(module) %>%
  group_map(~{
    st <- do_stats(.x)
    tibble(module = unique(.x$module),
           condition_model_p = st$condition_p,
           model_formula = st$model_formula)
  }) %>% bind_rows() %>%
  mutate(condition_model_fdr = p.adjust(condition_model_p, method = "BH")) %>%
  arrange(condition_model_fdr)

write_csv_safe(stat_list, fp_traittab("ME_by_condition_adjusted_lm_FDR.csv"))
write_csv_safe(ME_long, fp_source("ME_by_condition.csv"))

top_modules <- head(stat_list$module, 12)
comparisons <- list(c("con","res"), c("con","sus"), c("res","sus"))

# Color mapping (full circles)
cond_cols <- nature_condition_cols

# Helper: build one dotplot
plot_dot_mod <- function(dfm, mod, condition_fdr = NA_real_, show_subtitle = TRUE, errorbar = c("none","sem","sd")) {
  errorbar <- match.arg(errorbar)
  # Summary for error bars
  summ <- dfm %>%
    group_by(condition) %>%
    summarise(mean = mean(ME, na.rm = TRUE),
              sd = sd(ME, na.rm = TRUE),
              n = dplyr::n(),
              se = sd / sqrt(pmax(n, 1)), .groups = "drop")
  # Base: all sample dots (full circles)
  p <- ggplot(dfm, aes(x = condition, y = ME, color = condition)) +
    geom_point(position = position_jitter(width = 0.08, height = 0, seed = 1),
               size = 1.4, alpha = 0.75, shape = 16, stroke = 0) +
    scale_color_manual(values = cond_cols, labels = nature_condition_labels, guide = "none") +
    scale_x_discrete(labels = nature_condition_labels)

  # Add mean points (slightly larger) on top
  p <- p + geom_point(data = summ, aes(x = condition, y = mean),
                      inherit.aes = FALSE, size = 2.4, shape = 16, color = "black") +
           geom_point(data = summ, aes(x = condition, y = mean, color = condition),
                      inherit.aes = FALSE, size = 2.0, shape = 16)

  # Optional error bars
  if (errorbar != "none") {
    if (errorbar == "sem") {
      p <- p + geom_errorbar(data = summ,
                             aes(x = condition, ymin = mean - se, ymax = mean + se, color = condition),
                             inherit.aes = FALSE, width = 0.10, linewidth = 0.25, alpha = 0.9)
    } else if (errorbar == "sd") {
      p <- p + geom_errorbar(data = summ,
                             aes(x = condition, ymin = mean - sd, ymax = mean + sd, color = condition),
                             inherit.aes = FALSE, width = 0.10, linewidth = 0.25, alpha = 0.9)
    }
  }

  subtitle_txt <- if (isTRUE(show_subtitle))
    sprintf("adjusted LM FDR=%s", ifelse(is.na(condition_fdr), "NA", signif(condition_fdr, 3))) else NULL

  p <- p +
    labs(title = paste0(mod, " eigengene by condition"),
         subtitle = subtitle_txt,
         x = NULL, y = "Module eigengene") +
    theme_nature() +
    theme(panel.grid.major.x = element_blank())

  # Pairwise significance labels (Wilcoxon BH), shown as p.signif above groups
  p <- p + ggpubr::stat_compare_means(comparisons = comparisons,
                                      method = "wilcox.test",
                                      p.adjust.method = "BH",
                                      label = "p.signif",
                                      hide.ns = TRUE)
  p
}

# Multi-page PDF: all modules as dotplots
pdf(fp_traits("ME_by_condition_all_modules_dotplot.pdf"), width = 7, height = 4.5)
for (mod in unique(ME_long$module)) {
  dfm <- ME_long %>% filter(module == mod)
  condition_fdr <- stat_list$condition_model_fdr[stat_list$module == mod][1]
  p <- plot_dot_mod(dfm, mod, condition_fdr = condition_fdr, show_subtitle = TRUE, errorbar = "sem")
  print(p)
}
dev.off()

# SVG grid: top N most differential modules (compact)
plot_one_mod_dot <- function(mod) {
  dfm <- ME_long %>% filter(module == mod)
  condition_fdr <- stat_list$condition_model_fdr[stat_list$module == mod][1]
  summ <- dfm %>% group_by(condition) %>%
    summarise(mean = mean(ME, na.rm = TRUE), se = sd(ME, na.rm = TRUE)/sqrt(n()), .groups = "drop")

  ggplot(dfm, aes(x = condition, y = ME, color = condition)) +
    geom_point(position = position_jitter(width = 0.08, height = 0, seed = 1),
               size = 1.1, alpha = 0.75, shape = 16, stroke = 0) +
    geom_point(data = summ, aes(x = condition, y = mean),
               inherit.aes = FALSE, size = 2.2, shape = 16, color = "black") +
    geom_point(data = summ, aes(x = condition, y = mean, color = condition),
               inherit.aes = FALSE, size = 1.8, shape = 16) +
    # Optional tiny SEM bars for overview
    geom_errorbar(data = summ, aes(x = condition, ymin = mean - se, ymax = mean + se, color = condition),
                  inherit.aes = FALSE, width = 0.08, linewidth = 0.25, alpha = 0.8) +
    scale_color_manual(values = cond_cols, labels = nature_condition_labels, guide = "none") +
    scale_x_discrete(labels = nature_condition_labels) +
    labs(title = paste0(mod, "  (FDR=", signif(condition_fdr, 3), ")"),
         x = NULL, y = NULL) +
    theme_nature() +
    theme(legend.position = "none",
          plot.title = element_text(size = 6.2),
          panel.grid.major.x = element_blank())
}

if (length(top_modules) > 0) {
  plots_top <- lapply(top_modules, plot_one_mod_dot)
  ncol_grid <- min(4, ceiling(sqrt(length(plots_top))))
  nrow_grid <- ceiling(length(plots_top)/ncol_grid)
  g <- patchwork::wrap_plots(plots_top, ncol = ncol_grid)
  save_plot_nature(g, fp_mainfig("ME_by_condition_top_modules_dotplot.svg"),
                   width = min(nature_double_col, 1.75*ncol_grid),
                   height = 1.55*nrow_grid)
}

# Optional: export pairwise Wilcoxon BH summary across modules
pw_tables <- ME_long %>%
  group_by(module) %>%
  group_map(~{
    tt <- pairwise.wilcox.test(.x$ME, .x$condition, p.adjust.method = "BH", exact = FALSE)
    broom::tidy(tt) %>% mutate(module = unique(.x$module))
  }) %>% bind_rows()
write_csv_safe(pw_tables, fp_traittab("ME_by_condition_pairwise_Wilcoxon_BH.csv"))

output_manifest <- tibble::tibble(
  category = c(
    "logs", "logs", "qc_table", "mapping_table", "state",
    "network_figure", "trait_figure", "main_figure", "main_figure",
    "source_data", "source_data", "source_data", "source_data",
    "module_table", "preservation_table", "trait_table"
  ),
  path = c(
    fp_log("session_info.txt"),
    fp_log("analysis_parameters.csv"),
    fp_qctab("soft_threshold_selection.csv"),
    fp_maptab("mapping_audit_mouse_only_robust.tsv"),
    fp_state("mouse_only_mapping_outputs_robust.rds"),
    fp_net("gene_dendrogram_modules_merged.svg"),
    fp_traits("ME_trait_pheatmap.svg"),
    fp_mainfig("panel_module_trait_relationships_spatial.svg"),
    fp_mainfig("ME_by_condition_top_modules_dotplot.svg"),
    fp_source("ME_trait_correlations.csv"),
    fp_source("ME_condition_spatial_strata.csv"),
    fp_source("ME_trait_block_correlations.csv"),
    fp_source("ME_by_condition.csv"),
    fp_modtab("module_name_map.csv"),
    fp_prestab("module_preservation_CON_vs_ALL_Zsummary.csv"),
    fp_traittab("ME_by_condition_adjusted_lm_FDR.csv")
  )
)
write_csv_safe(output_manifest, fp_log("output_manifest.csv"))


# --------------------------
# End of script
# --------------------------

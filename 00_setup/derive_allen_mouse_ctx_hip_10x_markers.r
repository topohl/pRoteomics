#!/usr/bin/env Rscript

# derive_allen_mouse_ctx_hip_10x_markers.R
#
# Purpose:
#   Derive broad cell-type marker panels from Allen Mouse Whole Cortex and
#   Hippocampus 10x "Gene Expression by Cluster, trimmed means", using
#   dend.RData taxonomy where possible.
#
# Inputs:
#   data/external/reference_markers/allen_mouse_ctx_hip_10x/
#     allen_mouse_ctx_hip_10x_trimmed_means.csv
#     dend.RData
#
# Outputs:
#   data/external/reference_markers/derived/
#     allen_mouse_ctx_hip_10x_cluster_taxonomy_from_dend.csv
#     allen_mouse_ctx_hip_10x_cluster_broad_label_map.csv
#     allen_mouse_ctx_hip_10x_broad_expression.csv
#     allen_mouse_ctx_hip_10x_marker_candidates.csv
#     allen_mouse_ctx_hip_10x_marker_panels_top.csv
#     allen_mouse_ctx_hip_10x_marker_panels_top_proteomics_intersect.csv
#     allen_mouse_ctx_hip_10x_source_manifest.csv
#
# Notes:
#   - This script does not use the full cell-by-gene matrix.
#   - It uses aggregated trimmed means per Allen cluster.
#   - It uses dend.RData if available.
#   - If dend.RData cannot be parsed sufficiently, it falls back to
#     cluster-name-based broad labels.
#   - Interpretation should remain broad: cell-type support / contamination /
#     compartment signal, not fine subtype claims.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(purrr)
  library(tibble)
})

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

message2 <- function(...) {
  message(sprintf(...))
}

normalize_gene_symbol <- function(x) {
  x %>%
    as.character() %>%
    str_trim()
}

clean_text <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("[\r\n\t]+", " ") %>%
    str_squish()
}

make_cluster_id <- function(x) {
  str_extract(as.character(x), "^[0-9]+")
}

make_cluster_label <- function(x) {
  str_remove(as.character(x), "^[0-9]+_")
}

assign_broad_label_from_text <- function(text) {
  text <- clean_text(text)

  case_when(
    str_detect(text, regex("Micro|PVM|Macrophage|Myeloid", ignore_case = TRUE)) ~
      "microglia_pvm",

    str_detect(text, regex("Astro", ignore_case = TRUE)) ~
      "astrocyte",

    str_detect(text, regex("Oligo", ignore_case = TRUE)) ~
      "oligodendrocyte",

    str_detect(text, regex("\\bOPC\\b|Oligodendrocyte precursor", ignore_case = TRUE)) ~
      "opc",

    str_detect(text, regex("Endo|SMC|Peri|VLMC|Vascular|Pericyte|Mural", ignore_case = TRUE)) ~
      "vascular",

    str_detect(text, regex("CA1|CA2|CA3|DG|SUB|ProS|Mossy|HATA|HPF|Hippoc", ignore_case = TRUE)) ~
      "hippocampal_excitatory_neuron",

    str_detect(text, regex("IT|PT|CT|NP|L2|L3|L4|L5|L6|Car3|ENT|RSP|ACA|PFC|CTX|Cortex", ignore_case = TRUE)) ~
      "cortical_excitatory_neuron",

    str_detect(text, regex("Pvalb|Sst|Vip|Lamp5|Sncg|Pax6|Meis2|GABA|Gad1|Gad2|Interneuron|Inh", ignore_case = TRUE)) ~
      "inhibitory_interneuron",

    TRUE ~
      "other"
  )
}

# -------------------------------------------------------------------------
# dend.RData parsing helpers
# -------------------------------------------------------------------------

get_rdata_objects <- function(path) {
  e <- new.env(parent = emptyenv())
  obj_names <- load(path, envir = e)
  objs <- mget(obj_names, envir = e)
  objs
}

is_dendrogram_like <- function(x) {
  inherits(x, "dendrogram")
}

flatten_dendrogram <- function(dend) {
  rows <- list()

  walk_node <- function(node, path_labels = character(), depth = 0L) {
    node_label <- attr(node, "label")
    node_label <- if (!is.null(node_label)) clean_text(node_label) else NA_character_

    next_path <- path_labels
    if (!is.na(node_label) && nzchar(node_label)) {
      next_path <- c(next_path, node_label)
    }

    if (is.leaf(node)) {
      leaf_label <- node_label

      rows[[length(rows) + 1L]] <<- tibble(
        leaf_label = leaf_label,
        path = paste(next_path, collapse = " | "),
        depth = depth
      )
    } else {
      for (child in node) {
        walk_node(child, next_path, depth + 1L)
      }
    }
  }

  walk_node(dend)

  bind_rows(rows) %>%
    mutate(
      allen_cluster = leaf_label,
      cluster_id = make_cluster_id(allen_cluster),
      cluster_label = make_cluster_label(allen_cluster),
      taxonomy_text = paste(allen_cluster, path, sep = " | "),
      broad_label_from_dend = assign_broad_label_from_text(taxonomy_text)
    ) %>%
    distinct()
}

# Recursive list/data.frame flattening for unknown Allen RData object shapes.
# This makes the script robust if dend.RData stores a list/table rather than
# a base R dendrogram.
flatten_unknown_taxonomy_object <- function(x, object_name = "dend_object") {
  out <- list()

  if (is.data.frame(x)) {
    df <- as_tibble(x)

    possible_cluster_cols <- names(df)[
      str_detect(
        names(df),
        regex("cluster|label|name|alias|cell|type|node", ignore_case = TRUE)
      )
    ]

    if (length(possible_cluster_cols) == 0) {
      possible_cluster_cols <- names(df)
    }

    df2 <- df %>%
      mutate(.row_id = row_number()) %>%
      mutate(
        taxonomy_text = pmap_chr(
          select(., any_of(possible_cluster_cols)),
          ~ paste(clean_text(c(...)), collapse = " | ")
        ),
        allen_cluster = str_extract(taxonomy_text, "\\b[0-9]+_[^|]+"),
        cluster_id = make_cluster_id(allen_cluster),
        cluster_label = make_cluster_label(allen_cluster),
        broad_label_from_dend = assign_broad_label_from_text(taxonomy_text),
        source_object = object_name
      ) %>%
      filter(!is.na(allen_cluster), allen_cluster != "") %>%
      select(
        source_object,
        allen_cluster,
        cluster_id,
        cluster_label,
        taxonomy_text,
        broad_label_from_dend
      ) %>%
      distinct()

    return(df2)
  }

  if (is.list(x)) {
    # Try to find useful labels in list attributes/names recursively.
    walk_list <- function(obj, path = character(), depth = 0L) {
      nm <- names(obj)

      attr_label <- attr(obj, "label")
      attr_name <- attr(obj, "name")
      attr_members <- attr(obj, "members")

      this_bits <- c(
        path,
        clean_text(attr_label %||% NA_character_),
        clean_text(attr_name %||% NA_character_)
      )
      this_bits <- this_bits[!is.na(this_bits) & nzchar(this_bits)]

      if (!is.null(nm) && length(nm) > 0) {
        this_bits <- c(this_bits, clean_text(nm))
      }

      text <- paste(unique(this_bits), collapse = " | ")

      possible_cluster <- str_extract(text, "\\b[0-9]+_[^|]+")

      if (!is.na(possible_cluster)) {
        out[[length(out) + 1L]] <<- tibble(
          source_object = object_name,
          allen_cluster = clean_text(possible_cluster),
          cluster_id = make_cluster_id(possible_cluster),
          cluster_label = make_cluster_label(possible_cluster),
          taxonomy_text = text,
          broad_label_from_dend = assign_broad_label_from_text(text)
        )
      }

      if (is.list(obj) && length(obj) > 0) {
        for (i in seq_along(obj)) {
          child_name <- names(obj)[i]
          child_name <- if (!is.null(child_name) && nzchar(child_name)) child_name else paste0("child_", i)
          child <- obj[[i]]
          if (is.list(child) || is.data.frame(child)) {
            walk_list(child, c(path, child_name), depth + 1L)
          }
        }
      }
    }

    walk_list(x)

    if (length(out) > 0) {
      return(bind_rows(out) %>% distinct())
    }
  }

  tibble()
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

load_dend_taxonomy <- function(dend_file) {
  if (!file.exists(dend_file)) {
    warning("dend.RData not found: ", dend_file)
    return(tibble())
  }

  message2("Loading dend.RData:")
  message2("  %s", dend_file)

  objs <- get_rdata_objects(dend_file)

  message2("Objects in dend.RData:")
  for (nm in names(objs)) {
    message2("  %s : %s", nm, paste(class(objs[[nm]]), collapse = ", "))
  }

  parsed <- list()

  for (nm in names(objs)) {
    obj <- objs[[nm]]

    if (is_dendrogram_like(obj)) {
      message2("Parsing dendrogram object: %s", nm)
      parsed[[length(parsed) + 1L]] <- flatten_dendrogram(obj) %>%
        mutate(source_object = nm) %>%
        select(
          source_object,
          allen_cluster,
          cluster_id,
          cluster_label,
          taxonomy_text,
          broad_label_from_dend
        )
    } else {
      message2("Trying generic parser for object: %s", nm)
      parsed[[length(parsed) + 1L]] <- flatten_unknown_taxonomy_object(obj, nm)
    }
  }

  taxonomy <- bind_rows(parsed) %>%
    filter(!is.na(allen_cluster), allen_cluster != "") %>%
    mutate(
      cluster_id = make_cluster_id(allen_cluster),
      cluster_label = make_cluster_label(allen_cluster),
      taxonomy_text = clean_text(taxonomy_text),
      broad_label_from_dend = assign_broad_label_from_text(taxonomy_text)
    ) %>%
    distinct()

  taxonomy
}

# -------------------------------------------------------------------------
# Optional proteomics detected genes
# -------------------------------------------------------------------------

read_optional_detected_genes <- function(paths) {
  existing <- paths[file.exists(paths)]

  if (length(existing) == 0) {
    return(tibble())
  }

  map_dfr(existing, function(path) {
    dataset <- basename(path) %>%
      str_remove("\\.csv$") %>%
      str_remove("\\.tsv$") %>%
      str_remove("\\.txt$") %>%
      str_remove("_detected_genes$") %>%
      str_remove("_detected_proteins$")

    delim <- ifelse(str_detect(path, "\\.tsv$|\\.txt$"), "\t", ",")

    x <- suppressMessages(read_delim(path, delim = delim, show_col_types = FALSE))

    gene_col <- names(x)[
      str_to_lower(names(x)) %in% c(
        "gene", "gene_symbol", "symbol", "genesymbol",
        "external_gene_name", "protein_gene_symbol"
      )
    ][1]

    if (is.na(gene_col)) {
      gene_col <- names(x)[1]
      warning("No obvious gene-symbol column found in ", path,
              ". Using first column: ", gene_col)
    }

    tibble(
      dataset = dataset,
      gene_symbol = normalize_gene_symbol(x[[gene_col]])
    ) %>%
      filter(!is.na(gene_symbol), gene_symbol != "") %>%
      distinct()
  })
}

# -------------------------------------------------------------------------
# Paths
# -------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

repo_root <- if (length(args) >= 1) {
  args[[1]]
} else {
  "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics"
}

repo_root <- normalizePath(repo_root, winslash = "/", mustWork = FALSE)

reference_dir <- file.path(
  repo_root,
  "data/external/reference_markers/allen_mouse_ctx_hip_10x"
)

derived_dir <- file.path(
  repo_root,
  "data/external/reference_markers/derived"
)

dir.create(reference_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)

trimmed_file <- file.path(
  reference_dir,
  "allen_mouse_ctx_hip_10x_trimmed_means.csv"
)

dend_file <- file.path(
  reference_dir,
  "dend.RData"
)

if (!file.exists(trimmed_file)) {
  stop(
    "Missing Allen trimmed-means file:\n  ", trimmed_file, "\n\n",
    "Expected file structure:\n",
    "  feature,108_Pvalb,229_L6 IT CTX,26_Ntng1 HPF,...\n\n",
    "Rename the Allen trimmed-means CSV to:\n",
    "  allen_mouse_ctx_hip_10x_trimmed_means.csv\n"
  )
}

detected_gene_files <- c(
  file.path(repo_root, "data/processed/01_preprocessing/detected_genes_microglia.csv"),
  file.path(repo_root, "data/processed/01_preprocessing/detected_genes_neuron_neuropil.csv"),
  file.path(repo_root, "data/processed/01_preprocessing/detected_genes_neuron_soma.csv"),
  file.path(repo_root, "data/external/reference_markers/detected_genes_microglia.csv"),
  file.path(repo_root, "data/external/reference_markers/detected_genes_neuron_neuropil.csv"),
  file.path(repo_root, "data/external/reference_markers/detected_genes_neuron_soma.csv")
)

# -------------------------------------------------------------------------
# Parameters
# -------------------------------------------------------------------------

min_mean_expr <- 1.0
min_specificity_ratio <- 1.5
top_n_per_group <- 100

valid_broad_labels <- c(
  "microglia_pvm",
  "astrocyte",
  "oligodendrocyte",
  "opc",
  "vascular",
  "hippocampal_excitatory_neuron",
  "cortical_excitatory_neuron",
  "inhibitory_interneuron"
)

# -------------------------------------------------------------------------
# Load taxonomy from dend.RData
# -------------------------------------------------------------------------

taxonomy <- load_dend_taxonomy(dend_file)

taxonomy_out <- file.path(
  derived_dir,
  "allen_mouse_ctx_hip_10x_cluster_taxonomy_from_dend.csv"
)

if (nrow(taxonomy) > 0) {
  write_csv(taxonomy, taxonomy_out)
  message2("Wrote parsed dend taxonomy:")
  message2("  %s", taxonomy_out)
  message2("Parsed taxonomy rows: %s", format(nrow(taxonomy), big.mark = ","))
} else {
  warning(
    "Could not parse useful cluster taxonomy from dend.RData. ",
    "The script will fall back to broad labels inferred from trimmed-means column names."
  )
}

# -------------------------------------------------------------------------
# Read trimmed means
# -------------------------------------------------------------------------

message2("Reading Allen trimmed-means file:")
message2("  %s", trimmed_file)

expr_wide <- fread(trimmed_file, check.names = FALSE)

if (!"feature" %in% names(expr_wide)) {
  stop(
    "Expected first gene column named 'feature'. Found columns:\n",
    paste(head(names(expr_wide), 20), collapse = ", ")
  )
}

expr_wide <- as_tibble(expr_wide) %>%
  rename(gene_symbol = feature) %>%
  mutate(gene_symbol = normalize_gene_symbol(gene_symbol)) %>%
  filter(!is.na(gene_symbol), gene_symbol != "")

cluster_cols <- setdiff(names(expr_wide), "gene_symbol")

if (length(cluster_cols) < 10) {
  stop("Too few cluster columns detected. Check whether the file was read correctly.")
}

message2("Detected genes in trimmed-means table: %s", format(nrow(expr_wide), big.mark = ","))
message2("Detected Allen cluster columns: %s", format(length(cluster_cols), big.mark = ","))

# -------------------------------------------------------------------------
# Build cluster map from trimmed-means columns and join dend taxonomy
# -------------------------------------------------------------------------

cluster_map <- tibble(
  allen_cluster = cluster_cols,
  cluster_id = make_cluster_id(allen_cluster),
  cluster_label = make_cluster_label(allen_cluster),
  broad_label_from_name = assign_broad_label_from_text(cluster_label)
)

if (nrow(taxonomy) > 0) {
  taxonomy_slim <- taxonomy %>%
    select(
      allen_cluster,
      cluster_id,
      cluster_label,
      taxonomy_text,
      broad_label_from_dend
    ) %>%
    distinct()

  cluster_map <- cluster_map %>%
    left_join(
      taxonomy_slim %>%
        select(cluster_id, taxonomy_text, broad_label_from_dend) %>%
        distinct(),
      by = "cluster_id"
    ) %>%
    mutate(
      broad_label = case_when(
        !is.na(broad_label_from_dend) &
          broad_label_from_dend != "other" ~ broad_label_from_dend,
        TRUE ~ broad_label_from_name
      ),
      broad_label_source = case_when(
        !is.na(broad_label_from_dend) &
          broad_label_from_dend != "other" ~ "dend.RData",
        TRUE ~ "cluster_column_name"
      )
    )
} else {
  cluster_map <- cluster_map %>%
    mutate(
      taxonomy_text = NA_character_,
      broad_label_from_dend = NA_character_,
      broad_label = broad_label_from_name,
      broad_label_source = "cluster_column_name"
    )
}

cluster_map_out <- file.path(
  derived_dir,
  "allen_mouse_ctx_hip_10x_cluster_broad_label_map.csv"
)

cluster_map %>%
  arrange(broad_label, as.numeric(cluster_id), cluster_label) %>%
  write_csv(cluster_map_out)

message2("Wrote cluster broad-label map:")
message2("  %s", cluster_map_out)

message2("Cluster counts by broad_label:")
cluster_map %>%
  count(broad_label, broad_label_source, name = "n_clusters") %>%
  arrange(broad_label, broad_label_source) %>%
  print(n = Inf)

# -------------------------------------------------------------------------
# Reshape wide gene x cluster matrix to long format
# -------------------------------------------------------------------------

message2("Reshaping wide gene x cluster matrix to long format...")

expr_long <- expr_wide %>%
  pivot_longer(
    cols = all_of(cluster_cols),
    names_to = "allen_cluster",
    values_to = "trimmed_mean_expr"
  ) %>%
  mutate(
    trimmed_mean_expr = suppressWarnings(as.numeric(trimmed_mean_expr))
  ) %>%
  left_join(
    cluster_map %>%
      select(
        allen_cluster,
        cluster_id,
        cluster_label,
        taxonomy_text,
        broad_label,
        broad_label_source
      ),
    by = "allen_cluster"
  )

long_out <- file.path(
  derived_dir,
  "allen_mouse_ctx_hip_10x_expression_long.csv.gz"
)

message2("Writing long expression table:")
message2("  %s", long_out)

write_csv(expr_long, long_out)

# -------------------------------------------------------------------------
# Aggregate expression by broad biological group
# -------------------------------------------------------------------------

message2("Aggregating expression by broad_label...")

broad_expr <- expr_long %>%
  filter(
    !is.na(trimmed_mean_expr),
    !is.na(broad_label),
    broad_label != "other"
  ) %>%
  group_by(gene_symbol, broad_label) %>%
  summarise(
    mean_expr = mean(trimmed_mean_expr, na.rm = TRUE),
    median_expr = median(trimmed_mean_expr, na.rm = TRUE),
    max_expr = max(trimmed_mean_expr, na.rm = TRUE),
    n_clusters = n(),
    n_clusters_expr_gt_0 = sum(trimmed_mean_expr > 0, na.rm = TRUE),
    frac_clusters_expr_gt_0 = n_clusters_expr_gt_0 / n_clusters,
    .groups = "drop"
  )

broad_expr_out <- file.path(
  derived_dir,
  "allen_mouse_ctx_hip_10x_broad_expression.csv"
)

write_csv(broad_expr, broad_expr_out)

message2("Wrote broad expression table:")
message2("  %s", broad_expr_out)

# -------------------------------------------------------------------------
# Derive specificity metrics
# -------------------------------------------------------------------------

message2("Calculating specificity metrics...")

marker_candidates <- broad_expr %>%
  group_by(gene_symbol) %>%
  mutate(
    max_other_mean_expr = map_dbl(
      broad_label,
      ~ {
        other_vals <- mean_expr[broad_label != .x]
        if (length(other_vals) == 0 || all(is.na(other_vals))) {
          NA_real_
        } else {
          max(other_vals, na.rm = TRUE)
        }
      }
    ),
    mean_other_mean_expr = map_dbl(
      broad_label,
      ~ {
        other_vals <- mean_expr[broad_label != .x]
        if (length(other_vals) == 0 || all(is.na(other_vals))) {
          NA_real_
        } else {
          mean(other_vals, na.rm = TRUE)
        }
      }
    ),
    specificity_ratio = (mean_expr + 0.1) / (max_other_mean_expr + 0.1),
    specificity_delta = mean_expr - max_other_mean_expr
  ) %>%
  ungroup() %>%
  mutate(
    marker_rank_score =
      log2(specificity_ratio) *
      log2(mean_expr + 1) *
      pmax(frac_clusters_expr_gt_0, 0.05)
  ) %>%
  arrange(
    broad_label,
    desc(marker_rank_score),
    desc(specificity_ratio),
    desc(mean_expr)
  )

marker_candidates_out <- file.path(
  derived_dir,
  "allen_mouse_ctx_hip_10x_marker_candidates.csv"
)

write_csv(marker_candidates, marker_candidates_out)

message2("Wrote full marker candidate table:")
message2("  %s", marker_candidates_out)

# -------------------------------------------------------------------------
# Conservative top marker panels
# -------------------------------------------------------------------------

message2("Selecting conservative top marker panels...")

top_markers <- marker_candidates %>%
  filter(
    broad_label %in% valid_broad_labels,
    mean_expr >= min_mean_expr,
    specificity_ratio >= min_specificity_ratio
  ) %>%
  group_by(broad_label) %>%
  arrange(
    desc(marker_rank_score),
    desc(specificity_ratio),
    desc(mean_expr),
    .by_group = TRUE
  ) %>%
  slice_head(n = top_n_per_group) %>%
  ungroup() %>%
  mutate(
    reference_source = "Allen Mouse Whole Cortex and Hippocampus 10x trimmed means + dend.RData taxonomy",
    panel_group = "celltype",
    panel_name = broad_label,
    evidence_type = "allen_trimmed_mean_specificity",
    use_for = "qc_module_annotation"
  ) %>%
  select(
    reference_source,
    panel_group,
    panel_name,
    gene_symbol,
    evidence_type,
    use_for,
    mean_expr,
    median_expr,
    max_expr,
    max_other_mean_expr,
    mean_other_mean_expr,
    specificity_ratio,
    specificity_delta,
    marker_rank_score,
    n_clusters,
    n_clusters_expr_gt_0,
    frac_clusters_expr_gt_0
  )

top_markers_out <- file.path(
  derived_dir,
  "allen_mouse_ctx_hip_10x_marker_panels_top.csv"
)

write_csv(top_markers, top_markers_out)

message2("Wrote top marker panels:")
message2("  %s", top_markers_out)

message2("Top marker count by panel:")
top_markers %>%
  count(panel_name, name = "n_markers") %>%
  arrange(panel_name) %>%
  print(n = Inf)

# -------------------------------------------------------------------------
# Optional: intersect with detected proteomics genes
# -------------------------------------------------------------------------

detected_genes <- read_optional_detected_genes(detected_gene_files)

if (nrow(detected_genes) > 0) {
  message2("Detected proteomics gene files found. Intersecting markers...")

  detected_wide <- detected_genes %>%
    mutate(value = TRUE) %>%
    distinct(dataset, gene_symbol, value) %>%
    pivot_wider(
      names_from = dataset,
      values_from = value,
      values_fill = FALSE,
      names_prefix = "detected_"
    )

  intersected_wide <- top_markers %>%
    left_join(detected_wide, by = "gene_symbol") %>%
    mutate(across(starts_with("detected_"), ~ if_else(is.na(.x), FALSE, .x)))

  intersected_out <- file.path(
    derived_dir,
    "allen_mouse_ctx_hip_10x_marker_panels_top_proteomics_intersect.csv"
  )

  write_csv(intersected_wide, intersected_out)

  message2("Wrote proteomics-intersected marker panels:")
  message2("  %s", intersected_out)

} else {
  message2("No detected proteomics gene files found. Skipping proteomics intersection.")
  message2("Expected optional files such as:")
  for (p in detected_gene_files) {
    message2("  %s", p)
  }
}

# -------------------------------------------------------------------------
# Source manifest
# -------------------------------------------------------------------------

manifest <- tibble(
  source_id = "allen_mouse_ctx_hip_10x_trimmed_means_dend",
  source_name = "Allen Mouse Whole Cortex and Hippocampus 10x",
  source_file_trimmed_means = trimmed_file,
  source_file_taxonomy = dend_file,
  derived_from = "Gene Expression by Cluster, trimmed means; cluster taxonomy dend.RData",
  use = "Broad reference marker derivation for proteomics QC and WGCNA/module annotation",
  limitations = paste(
    "Markers are derived from transcriptomic trimmed-mean expression, not protein abundance.",
    "Broad labels are assigned using dend.RData taxonomy where parseable, with cluster-name fallback.",
    "Use for broad cell-type support, contamination checks, and module annotation only.",
    "Do not treat these as definitive fine cell-subtype markers in proteomics."
  )
)

manifest_out <- file.path(
  derived_dir,
  "allen_mouse_ctx_hip_10x_source_manifest.csv"
)

write_csv(manifest, manifest_out)

message2("Wrote source manifest:")
message2("  %s", manifest_out)

message2("Done.")
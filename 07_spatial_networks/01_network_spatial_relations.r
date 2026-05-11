# ================================================================
# Spatial layer-region network analysis for proteomics matrices
# ================================================================
# Purpose:
#   Quantify relationships between hippocampal region/layer units from
#   protein expression data and export network-ready tables and figures.
#
# Main outputs:
#   1) Region/layer expression similarity network
#   2) Top-protein Jaccard-overlap network
#   3) Group-specific spatial similarity networks, if ExpGroup is available
#   4) Network centrality tables
#   5) Heatmaps, network plots, and run logs
#
# Notes:
#   - This script is deliberately separate from WGCNA.r. WGCNA infers
#     protein co-expression modules; this script infers relationships among
#     anatomical sampling units such as CA1_sr, CA3_sr, DG_mo, etc.
#   - Edge-level interpretation is unstable at low n. Treat single edges as
#     exploratory unless bootstrap/sensitivity support is added.
# ================================================================

# -------------------------------
# 0) Packages
# -------------------------------
required_pkgs <- c(
  "readxl", "dplyr", "tidyr", "stringr", "purrr", "tibble", "ggplot2",
  "pheatmap", "igraph", "ggraph", "openxlsx", "scales", "svglite"
)

install_if_missing <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}
install_if_missing(required_pkgs)

# -------------------------------
# 1) User parameters
# -------------------------------
params <- list(
  # Typical Exp9 neuron-neuropil matrix. Change as needed.
  #protein_file = "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Datasets/morpheus/20260218_pgmatrix_imputed_neuron_neuropil_180samples_missing70pct_with_metadata.xlsx",
  protein_file = "/Users/tobiaspohl/Documents/pRoteomics/20260218_pgmatrix_imputed_neuron_neuropil_180samples_missing70pct.xlsx",

  # EWCE/sample metadata file. The script uses sample_id/region/layer/celltype_layer/ExpGroup when available.
  #metadata_file = "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Results/module_scores/sample_metadata_merged_clean_for_module_scores.xlsx",
  metadata_file = "/Users/tobiaspohl/Documents/Data/proteomics/TPE9_sample_metadata_males.xlsx",

  #output_dir = "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Results/network_spatial_relations",
  output_dir = "/Users/tobiaspohl/Documents/pRoteomics/Results/network_spatial_relations",

  # Minimum fraction of non-missing values required for a protein after sample-column selection.
  min_nonmissing_fraction = 0.5,

  # Top proteins used for Jaccard overlap per region/layer.
  top_n_proteins = 100,

  # Edge filters for exported network tables and visual networks.
  min_abs_spearman_r = 0.60,
  min_jaccard = 0.10,

  # If TRUE, z-score each protein across all samples before aggregation.
  # Recommended for spatial relation networks, because it avoids high-abundance proteins dominating.
  zscore_proteins = TRUE,

  # If TRUE, keep only neuron neuropil samples from metadata column celltype_layer.
  # This prevents soma/microglia samples from entering the spatial neuropil network.
  filter_neuron_neuropil_only = TRUE,

  # Numeric ExpGroup recoding used by the EWCE metadata. Adjust here if EWCE_E9 used another coding.
  # Common Exp9 convention assumed here: 1 = CON, 2 = RES, 3 = SUS.
  expgroup_numeric_map = c("1" = "CON", "2" = "RES", "3" = "SUS"),

  # If TRUE, also creates group-specific spatial similarity networks where ExpGroup is available.
  run_group_specific_networks = TRUE,

  # Group level order for plots/tables.
  group_levels = c("CON", "RES", "SUS")
)

# -------------------------------
# 2) Helpers
# -------------------------------
message2 <- function(...) message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ...)

safe_name <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_replace_all("[^A-Za-z0-9_\\-]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_|_$", "")
}

make_dirs <- function(base_dir) {
  dirs <- list(
    base = base_dir,
    tables = file.path(base_dir, "01_Tables"),
    figures = file.path(base_dir, "02_Figures"),
    networks = file.path(base_dir, "03_Network_Files"),
    logs = file.path(base_dir, "04_Logs")
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  dirs
}

theme_nature_network <- function(base_size = 8) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(family = "sans", colour = "black"),
      axis.text = ggplot2::element_text(colour = "black"),
      axis.line = ggplot2::element_line(linewidth = 0.3),
      axis.ticks = ggplot2::element_line(linewidth = 0.3),
      legend.key = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

find_first_col <- function(df, candidates) {
  nms <- names(df)
  nms_clean <- tolower(gsub("[^a-z0-9]", "", nms))
  cand_clean <- tolower(gsub("[^a-z0-9]", "", candidates))
  idx <- match(cand_clean, nms_clean)
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0) return(NA_character_)
  nms[idx[1]]
}

coerce_numeric_safely <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

row_zscore <- function(mat) {
  m <- rowMeans(mat, na.rm = TRUE)
  s <- apply(mat, 1, sd, na.rm = TRUE)
  s[is.na(s) | s == 0] <- 1
  sweep(sweep(mat, 1, m, "-"), 1, s, "/")
}

normalize_sample_key <- function(x) {
  x %>%
    as.character() %>%
    basename() %>%
    stringr::str_replace("\\.d$", "") %>%
    stringr::str_to_lower()
}

normalize_expgroup <- function(x, numeric_map = params$expgroup_numeric_map) {
  x_chr <- as.character(x)
  x_chr <- stringr::str_trim(x_chr)
  x_up <- toupper(x_chr)

  mapped_numeric <- unname(numeric_map[x_chr])

  dplyr::case_when(
    !is.na(mapped_numeric) ~ mapped_numeric,
    x_up %in% c("CONTROL", "CTRL", "CON") ~ "CON",
    x_up %in% c("RESILIENT", "RES") ~ "RES",
    x_up %in% c("SUSCEPTIBLE", "SUS") ~ "SUS",
    x_up %in% c("NA", "NAN", "") ~ NA_character_,
    TRUE ~ x_up
  )
}

parse_sample_metadata_from_names <- function(sample_names) {
  sample_key <- basename(sample_names)

  tibble::tibble(
    SampleColumn = sample_names,
    SampleKey = sample_key,
    Region = stringr::str_extract(SampleKey, regex("CA1|CA2|CA3|DG", ignore_case = TRUE)) %>% toupper(),
    Layer = stringr::str_extract(SampleKey, regex("slm|sr|so|mo|po|sp|sg", ignore_case = TRUE)) %>% tolower(),
    ExpGroup = stringr::str_extract(SampleKey, regex("CON|RES|SUS|control|resilient|susceptible", ignore_case = TRUE)) %>% normalize_expgroup()
  ) %>%
    dplyr::mutate(
      RegionLayer = ifelse(!is.na(Region) & !is.na(Layer), paste(Region, Layer, sep = "_"), NA_character_)
    )
}

standardize_metadata <- function(metadata_df, sample_names, numeric_map = params$expgroup_numeric_map) {
  md <- as.data.frame(metadata_df)

  sample_col <- find_first_col(md, c(
    "sample_id", "SampleID", "SampleColumn", "Raw_Sample", "Sample", "SampleName",
    "sample_name", "row.names", "filename", "File", "Run", "Label"
  ))
  region_col <- find_first_col(md, c("region", "Region", "BrainRegion"))
  layer_col  <- find_first_col(md, c("layer", "Layer", "Stratum", "stratum"))
  group_col  <- find_first_col(md, c("ExpGroup", "EWCE_Group", "group", "Group", "Phenotype", "Condition"))
  celltype_layer_col <- find_first_col(md, c("celltype_layer", "CelltypeLayer", "cell_type_layer"))
  celltype_col <- find_first_col(md, c("celltype", "Celltype", "CellType", "cell_type"))
  exclude_col <- find_first_col(md, c("exclude", "Exclude", "excluded", "Excluded"))

  out <- tibble::tibble(
    SampleColumn = sample_names,
    .JoinKey = normalize_sample_key(sample_names)
  )

  if (!is.na(sample_col)) {
    md2 <- md %>%
      dplyr::mutate(
        .SampleColumn = as.character(.data[[sample_col]]),
        .JoinKey = normalize_sample_key(.SampleColumn)
      ) %>%
      dplyr::select(.JoinKey, .SampleColumn, dplyr::everything()) %>%
      dplyr::distinct(.JoinKey, .keep_all = TRUE)

    out <- out %>%
      dplyr::left_join(md2, by = ".JoinKey")
  } else {
    warning("No sample identifier column found in metadata. Falling back to parsing sample names.")
  }

  parsed <- parse_sample_metadata_from_names(sample_names)

  out <- out %>%
    dplyr::mutate(
      Region = if (!is.na(region_col) && region_col %in% names(out)) as.character(.data[[region_col]]) else NA_character_,
      Layer = if (!is.na(layer_col) && layer_col %in% names(out)) as.character(.data[[layer_col]]) else NA_character_,
      ExpGroup = if (!is.na(group_col) && group_col %in% names(out)) as.character(.data[[group_col]]) else NA_character_,
      CelltypeLayer = if (!is.na(celltype_layer_col) && celltype_layer_col %in% names(out)) as.character(.data[[celltype_layer_col]]) else NA_character_,
      Celltype = if (!is.na(celltype_col) && celltype_col %in% names(out)) as.character(.data[[celltype_col]]) else NA_character_,
      Exclude = if (!is.na(exclude_col) && exclude_col %in% names(out)) as.character(.data[[exclude_col]]) else NA_character_
    ) %>%
    dplyr::mutate(
      Region = ifelse(is.na(Region) | Region == "", parsed$Region, Region),
      Layer = ifelse(is.na(Layer) | Layer == "" | toupper(Layer) == "NA", parsed$Layer, Layer),
      ExpGroup = ifelse(is.na(ExpGroup) | ExpGroup == "" | toupper(ExpGroup) == "NA", parsed$ExpGroup, ExpGroup),
      Region = toupper(Region),
      Layer = tolower(Layer),
      ExpGroup = normalize_expgroup(ExpGroup, numeric_map = numeric_map),
      CelltypeLayer = tolower(CelltypeLayer),
      Celltype = tolower(Celltype),
      Exclude = dplyr::case_when(
        tolower(Exclude) %in% c("true", "t", "1", "yes", "y") ~ TRUE,
        tolower(Exclude) %in% c("false", "f", "0", "no", "n", "") ~ FALSE,
        is.na(Exclude) ~ FALSE,
        TRUE ~ NA
      ),
      RegionLayer = ifelse(!is.na(Region) & !is.na(Layer), paste(Region, Layer, sep = "_"), NA_character_)
    ) %>%
    dplyr::select(
      SampleColumn, Region, Layer, RegionLayer, ExpGroup,
      CelltypeLayer, Celltype, Exclude, dplyr::everything()
    )

  out
}

select_expression_columns <- function(df, protein_id_col) {
  candidate_cols <- setdiff(names(df), protein_id_col)

  numeric_fraction <- vapply(candidate_cols, function(cc) {
    vals <- coerce_numeric_safely(df[[cc]])
    mean(!is.na(vals))
  }, numeric(1))

  # Keep columns with mostly numeric expression values.
  expr_cols <- candidate_cols[numeric_fraction >= 0.70]

  # Prefer path-like/sample-like expression columns if present, but do not require this.
  path_like <- grepl("(^[A-Za-z]:)|(/)|(_CON|_RES|_SUS|CON|RES|SUS|CA1|CA2|CA3|DG)", expr_cols, ignore.case = TRUE)
  if (sum(path_like) >= 4) expr_cols <- expr_cols[path_like]

  expr_cols
}

load_expression_matrix <- function(protein_file, metadata_file = NULL, min_nonmissing_fraction = 0.5, zscore = TRUE) {
  message2("Reading protein matrix: ", protein_file)
  raw <- readxl::read_excel(protein_file)
  raw <- as.data.frame(raw)

  protein_col <- find_first_col(raw, c(
    "Gene", "Genes", "gene_symbol", "GeneSymbol", "Majority protein IDs",
    "Protein IDs", "Protein.IDs", "Accession", "UniprotID", "UniProt", "ID"
  ))
  if (is.na(protein_col)) {
    protein_col <- names(raw)[1]
    warning("No known protein/gene ID column detected; using first column: ", protein_col)
  }

  raw[[protein_col]] <- as.character(raw[[protein_col]])
  raw[[protein_col]][is.na(raw[[protein_col]]) | raw[[protein_col]] == ""] <- paste0("missing_id_", seq_len(sum(is.na(raw[[protein_col]]) | raw[[protein_col]] == "")))
  raw[[protein_col]] <- make.unique(raw[[protein_col]])

  expr_cols <- select_expression_columns(raw, protein_col)
  if (length(expr_cols) < 4) {
    stop("Fewer than 4 expression columns detected. Check protein_file and column naming.")
  }

  expr <- raw[, expr_cols, drop = FALSE]
  expr[] <- lapply(expr, coerce_numeric_safely)
  expr <- as.matrix(expr)
  rownames(expr) <- raw[[protein_col]]

  keep <- rowMeans(!is.na(expr)) >= min_nonmissing_fraction
  expr <- expr[keep, , drop = FALSE]

  if (zscore) expr <- row_zscore(expr)

  if (!is.null(metadata_file) && file.exists(metadata_file)) {
    message2("Reading metadata: ", metadata_file)
    metadata_raw <- readxl::read_excel(metadata_file)
    sample_md <- standardize_metadata(metadata_raw, colnames(expr), numeric_map = params$expgroup_numeric_map)
  } else {
    message2("No metadata file found; parsing region/layer/group from sample names.")
    sample_md <- parse_sample_metadata_from_names(colnames(expr))
  }

  sample_md <- sample_md %>%
    dplyr::filter(SampleColumn %in% colnames(expr)) %>%
    dplyr::distinct(SampleColumn, .keep_all = TRUE)

  if (isTRUE(params$filter_neuron_neuropil_only) && "CelltypeLayer" %in% names(sample_md)) {
    n_before <- nrow(sample_md)
    sample_md <- sample_md %>%
      dplyr::filter(is.na(Exclude) | Exclude == FALSE) %>%
      dplyr::filter(CelltypeLayer == "neuron_neuropil")
    message2("Filtered to neuron_neuropil non-excluded samples: ", nrow(sample_md), " / ", n_before)
  }

  message2("Sample metadata group counts:")
  print(sample_md %>% dplyr::count(ExpGroup, sort = TRUE))
  message2("Sample metadata region/layer counts:")
  print(sample_md %>% dplyr::count(RegionLayer, sort = TRUE))

  usable_samples <- sample_md %>% dplyr::filter(!is.na(RegionLayer)) %>% dplyr::pull(SampleColumn)
  if (length(usable_samples) < 4) {
    stop("Could not identify enough samples with RegionLayer metadata. Check metadata_file or sample names.")
  }

  expr <- expr[, usable_samples, drop = FALSE]
  sample_md <- sample_md %>% dplyr::filter(SampleColumn %in% usable_samples)

  list(expr = expr, sample_metadata = sample_md, protein_id_col = protein_col)
}

aggregate_region_layer_expression <- function(expr, sample_md) {
  sample_order <- match(colnames(expr), sample_md$SampleColumn)
  md <- sample_md[sample_order, ]

  long <- as.data.frame(expr) %>%
    tibble::rownames_to_column("Protein") %>%
    tidyr::pivot_longer(-Protein, names_to = "SampleColumn", values_to = "Expression") %>%
    dplyr::left_join(md %>% dplyr::select(SampleColumn, Region, Layer, RegionLayer, ExpGroup), by = "SampleColumn")

  agg <- long %>%
    dplyr::group_by(Protein, RegionLayer) %>%
    dplyr::summarise(MeanExpression = mean(Expression, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = RegionLayer, values_from = MeanExpression) %>%
    tibble::column_to_rownames("Protein") %>%
    as.matrix()

  # Remove proteins with all missing after aggregation.
  agg <- agg[rowSums(!is.na(agg)) > 1, , drop = FALSE]
  agg
}

compute_similarity_edges <- function(unit_matrix, min_abs_r = 0.6) {
  cor_mat <- suppressWarnings(cor(unit_matrix, method = "spearman", use = "pairwise.complete.obs"))
  units <- colnames(cor_mat)

  edges <- expand.grid(Source = units, Target = units, stringsAsFactors = FALSE) %>%
    tibble::as_tibble() %>%
    dplyr::filter(Source < Target) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(SpearmanR = cor_mat[Source, Target]) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      AbsSpearmanR = abs(SpearmanR),
      Direction = ifelse(SpearmanR >= 0, "positive", "negative")
    ) %>%
    dplyr::arrange(dplyr::desc(AbsSpearmanR))

  list(
    cor_mat = cor_mat,
    all_edges = edges,
    filtered_edges = edges %>% dplyr::filter(AbsSpearmanR >= min_abs_r)
  )
}

compute_top_protein_jaccard_edges <- function(unit_matrix, top_n = 100, min_jaccard = 0.1) {
  top_sets <- lapply(colnames(unit_matrix), function(unit) {
    vals <- unit_matrix[, unit]
    names(sort(vals, decreasing = TRUE, na.last = NA))[seq_len(min(top_n, sum(!is.na(vals))))]
  })
  names(top_sets) <- colnames(unit_matrix)

  units <- names(top_sets)
  edges <- expand.grid(Source = units, Target = units, stringsAsFactors = FALSE) %>%
    tibble::as_tibble() %>%
    dplyr::filter(Source < Target) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      IntersectionN = length(intersect(top_sets[[Source]], top_sets[[Target]])),
      UnionN = length(union(top_sets[[Source]], top_sets[[Target]])),
      Jaccard = ifelse(UnionN > 0, IntersectionN / UnionN, NA_real_),
      OverlapCoefficient = IntersectionN / min(length(top_sets[[Source]]), length(top_sets[[Target]]))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(Jaccard))

  list(
    top_sets = top_sets,
    all_edges = edges,
    filtered_edges = edges %>% dplyr::filter(!is.na(Jaccard), Jaccard >= min_jaccard)
  )
}

make_node_table <- function(sample_md, edge_tbl = NULL) {
  nodes <- sample_md %>%
  dplyr::filter(!is.na(RegionLayer)) %>%
  dplyr::count(Region, Layer, RegionLayer, name = "NSamples") %>%
  dplyr::rename(name = RegionLayer) %>%
  dplyr::select(name, Region, Layer, NSamples)

  if (!is.null(edge_tbl) && nrow(edge_tbl) > 0) {
    g <- igraph::graph_from_data_frame(edge_tbl %>% dplyr::select(Source, Target), directed = FALSE, vertices = nodes)
    cent <- tibble::tibble(
      name = igraph::V(g)$name,
      Degree = igraph::degree(g),
      Strength = igraph::strength(g, weights = rep(1, igraph::ecount(g))),
      Betweenness = igraph::betweenness(g, normalized = TRUE)
    )
    nodes <- nodes %>% dplyr::left_join(cent, by = "name")
  } else {
    nodes <- nodes %>% dplyr::mutate(Degree = 0, Strength = 0, Betweenness = 0)
  }

  nodes
}

plot_similarity_heatmap <- function(cor_mat, outfile, title = "Spatial expression similarity") {
  pheatmap::pheatmap(
    cor_mat,
    color = colorRampPalette(c("#457B9D", "white", "#E63946"))(101),
    breaks = seq(-1, 1, length.out = 102),
    border_color = NA,
    fontsize = 7,
    main = title,
    filename = outfile,
    width = 6.5,
    height = 5.8
  )
}

plot_network <- function(nodes, edges, outfile, edge_weight_col = "AbsSpearmanR", title = "Spatial network") {
  if (nrow(edges) == 0) {
    warning("No edges available for network plot: ", outfile)
    return(invisible(NULL))
  }

  edges2 <- edges %>%
    dplyr::rename(from = Source, to = Target) %>%
    dplyr::mutate(weight_for_plot = .data[[edge_weight_col]])

  g <- igraph::graph_from_data_frame(edges2, directed = FALSE, vertices = nodes)

  p <- ggraph::ggraph(g, layout = "fr") +
    ggraph::geom_edge_link(aes(width = weight_for_plot, alpha = weight_for_plot), colour = "grey45") +
    ggraph::geom_node_point(aes(size = NSamples), shape = 21, fill = "white", colour = "black", stroke = 0.4) +
    ggraph::geom_node_text(aes(label = name), repel = TRUE, size = 2.7) +
    ggraph::scale_edge_width(range = c(0.2, 1.5), guide = "none") +
    ggraph::scale_edge_alpha(range = c(0.25, 0.9), guide = "none") +
    ggplot2::scale_size_continuous(range = c(2, 6), name = "n samples") +
    ggplot2::labs(title = title) +
    ggplot2::theme_void(base_size = 8) +
    ggplot2::theme(
      text = element_text(family = "sans", colour = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  ggplot2::ggsave(outfile, p, width = 150, height = 120, units = "mm", device = svglite::svglite)
  invisible(p)
}

write_network_files <- function(edges, nodes, prefix, dirs) {
  readr_write_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE)
  readr_write_csv(nodes, file.path(dirs$networks, paste0(prefix, "_nodes.csv")))
  readr_write_csv(edges, file.path(dirs$networks, paste0(prefix, "_edges.csv")))
}

run_group_specific <- function(expr, sample_md, dirs, params) {
  groups_available <- sample_md %>%
    dplyr::filter(!is.na(ExpGroup), ExpGroup %in% params$group_levels) %>%
    dplyr::count(ExpGroup)

  if (nrow(groups_available) == 0) {
    message2("No usable ExpGroup metadata detected; skipping group-specific networks.")
    return(NULL)
  }

  outputs <- list()

  for (grp in groups_available$ExpGroup) {
    grp_samples <- sample_md %>% dplyr::filter(ExpGroup == grp, SampleColumn %in% colnames(expr)) %>% dplyr::pull(SampleColumn)
    if (length(grp_samples) < 4) next

    message2("Running group-specific spatial network: ", grp)
    expr_grp <- expr[, grp_samples, drop = FALSE]
    md_grp <- sample_md %>% dplyr::filter(SampleColumn %in% grp_samples)

    unit_mat <- aggregate_region_layer_expression(expr_grp, md_grp)
    if (ncol(unit_mat) < 3) next

    sim <- compute_similarity_edges(unit_mat, params$min_abs_spearman_r)
    jac <- compute_top_protein_jaccard_edges(unit_mat, params$top_n_proteins, params$min_jaccard)
    nodes <- make_node_table(md_grp, sim$filtered_edges)

    prefix <- paste0("group_", safe_name(grp))
    utils::write.csv(sim$all_edges, file.path(dirs$tables, paste0(prefix, "_spearman_all_edges.csv")), row.names = FALSE)
    utils::write.csv(sim$filtered_edges, file.path(dirs$tables, paste0(prefix, "_spearman_filtered_edges.csv")), row.names = FALSE)
    utils::write.csv(jac$all_edges, file.path(dirs$tables, paste0(prefix, "_jaccard_all_edges.csv")), row.names = FALSE)
    utils::write.csv(nodes, file.path(dirs$tables, paste0(prefix, "_nodes_centrality.csv")), row.names = FALSE)

    plot_similarity_heatmap(
      sim$cor_mat,
      file.path(dirs$figures, paste0(prefix, "_spearman_heatmap.pdf")),
      paste0("Spatial expression similarity — ", grp)
    )
    plot_network(
      nodes,
      sim$filtered_edges,
      file.path(dirs$figures, paste0(prefix, "_spearman_network.svg")),
      "AbsSpearmanR",
      paste0("Spatial similarity network — ", grp)
    )

    write_network_files(sim$filtered_edges, nodes, paste0(prefix, "_spearman"), dirs)

    outputs[[grp]] <- list(unit_matrix = unit_mat, spearman = sim, jaccard = jac, nodes = nodes)
  }

  outputs
}

# -------------------------------
# 3) Main analysis
# -------------------------------
dirs <- make_dirs(params$output_dir)
log_file <- file.path(dirs$logs, paste0("network_spatial_relations_run_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
sink(log_file, split = TRUE)
on.exit({ sink(); }, add = TRUE)

message2("Starting spatial layer-region network analysis")
print(params)

loaded <- load_expression_matrix(
  protein_file = params$protein_file,
  metadata_file = params$metadata_file,
  min_nonmissing_fraction = params$min_nonmissing_fraction,
  zscore = params$zscore_proteins
)

expr <- loaded$expr
sample_md <- loaded$sample_metadata

message2("Expression matrix retained: ", nrow(expr), " proteins x ", ncol(expr), " samples")
message2("Region/layer units: ", paste(sort(unique(sample_md$RegionLayer)), collapse = ", "))

# Save standardized metadata used by the script.
utils::write.csv(sample_md, file.path(dirs$tables, "standardized_sample_metadata_used.csv"), row.names = FALSE)

unit_matrix <- aggregate_region_layer_expression(expr, sample_md)
utils::write.csv(
  as.data.frame(unit_matrix) %>% tibble::rownames_to_column("Protein"),
  file.path(dirs$tables, "region_layer_mean_expression_matrix.csv"),
  row.names = FALSE
)

# Overall spatial similarity network.
message2("Computing overall region/layer Spearman similarity network")
sim <- compute_similarity_edges(unit_matrix, params$min_abs_spearman_r)
nodes <- make_node_table(sample_md, sim$filtered_edges)

utils::write.csv(sim$all_edges, file.path(dirs$tables, "overall_spearman_all_edges.csv"), row.names = FALSE)
utils::write.csv(sim$filtered_edges, file.path(dirs$tables, "overall_spearman_filtered_edges.csv"), row.names = FALSE)
utils::write.csv(nodes, file.path(dirs$tables, "overall_nodes_centrality.csv"), row.names = FALSE)
utils::write.csv(as.data.frame(sim$cor_mat) %>% tibble::rownames_to_column("RegionLayer"), file.path(dirs$tables, "overall_spearman_correlation_matrix.csv"), row.names = FALSE)

plot_similarity_heatmap(
  sim$cor_mat,
  file.path(dirs$figures, "overall_spearman_heatmap.pdf"),
  "Spatial expression similarity across region/layer units"
)
plot_network(
  nodes,
  sim$filtered_edges,
  file.path(dirs$figures, "overall_spearman_network.svg"),
  "AbsSpearmanR",
  "Spatial expression similarity network"
)
write_network_files(sim$filtered_edges, nodes, "overall_spearman", dirs)

# Top-protein overlap network.
message2("Computing top-protein Jaccard overlap network")
jac <- compute_top_protein_jaccard_edges(unit_matrix, params$top_n_proteins, params$min_jaccard)
jac_nodes <- make_node_table(sample_md, jac$filtered_edges %>% dplyr::rename(SpearmanR = Jaccard))
utils::write.csv(jac$all_edges, file.path(dirs$tables, "overall_top_protein_jaccard_all_edges.csv"), row.names = FALSE)
utils::write.csv(jac$filtered_edges, file.path(dirs$tables, "overall_top_protein_jaccard_filtered_edges.csv"), row.names = FALSE)
utils::write.csv(jac_nodes, file.path(dirs$tables, "overall_jaccard_nodes_centrality.csv"), row.names = FALSE)

# Export top-protein set membership in long format.
top_set_long <- purrr::imap_dfr(jac$top_sets, ~ tibble::tibble(RegionLayer = .y, Protein = .x))
utils::write.csv(top_set_long, file.path(dirs$tables, "top_protein_sets_by_region_layer.csv"), row.names = FALSE)

plot_network(
  jac_nodes,
  jac$filtered_edges,
  file.path(dirs$figures, "overall_top_protein_jaccard_network.svg"),
  "Jaccard",
  paste0("Top-", params$top_n_proteins, " protein overlap network")
)
write_network_files(jac$filtered_edges, jac_nodes, "overall_top_protein_jaccard", dirs)

# Optional group-specific networks.
group_outputs <- NULL
if (isTRUE(params$run_group_specific_networks)) {
  group_outputs <- run_group_specific(expr, sample_md, dirs, params)
}

# Workbook summary.
message2("Writing Excel summary workbook")
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Sample_Metadata")
openxlsx::writeData(wb, "Sample_Metadata", sample_md)
openxlsx::addWorksheet(wb, "Nodes")
openxlsx::writeData(wb, "Nodes", nodes)
openxlsx::addWorksheet(wb, "Spearman_All")
openxlsx::writeData(wb, "Spearman_All", sim$all_edges)
openxlsx::addWorksheet(wb, "Spearman_Filtered")
openxlsx::writeData(wb, "Spearman_Filtered", sim$filtered_edges)
openxlsx::addWorksheet(wb, "Jaccard_All")
openxlsx::writeData(wb, "Jaccard_All", jac$all_edges)
openxlsx::addWorksheet(wb, "Jaccard_Filtered")
openxlsx::writeData(wb, "Jaccard_Filtered", jac$filtered_edges)
openxlsx::addWorksheet(wb, "Top_Protein_Sets")
openxlsx::writeData(wb, "Top_Protein_Sets", top_set_long)
openxlsx::saveWorkbook(wb, file.path(dirs$tables, "network_spatial_relations_summary.xlsx"), overwrite = TRUE)

# Save R object for later figure assembly or sensitivity analyses.
saveRDS(
  list(
    params = params,
    expression_matrix = expr,
    sample_metadata = sample_md,
    region_layer_matrix = unit_matrix,
    overall_spearman = sim,
    overall_jaccard = jac,
    nodes = nodes,
    group_specific = group_outputs,
    sessionInfo = sessionInfo()
  ),
  file.path(dirs$logs, "network_spatial_relations_objects.rds")
)

message2("Finished spatial layer-region network analysis")
message2("Output directory: ", params$output_dir)
message2("Important caution: edges represent similarity/overlap, not causal protein-protein interactions.")

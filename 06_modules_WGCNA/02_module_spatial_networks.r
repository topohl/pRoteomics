# ================================================================
# Module-level spatial network analysis
# ================================================================
# Purpose:
#   Build region/layer networks from biologically interpretable protein modules
#   instead of all proteins. This links spatial proteomics to pathway/module
#   programs such as RNP/RNA processing, translation, mitochondrial
#   bioenergetics, and endolysosomal/proteostasis signals.
#
# Required first step:
#   Run Analysis/network_spatial_relations.r to create:
#   network_spatial_relations_objects.rds
#
# Inputs:
#   1) spatial_rds from network_spatial_relations.r
#   2) optional module_file with columns Module and Protein
#      If missing, conservative regex modules are created from protein names.
#
# Outputs:
#   - module scores per sample and region/layer
#   - module spatial similarity networks
#   - module-by-region heatmaps
#   - group-specific module spatial profiles
#   - Cytoscape/GraphML network files
# ================================================================

required_pkgs <- c(
  "dplyr", "tidyr", "stringr", "purrr", "tibble", "ggplot2", "pheatmap",
  "igraph", "ggraph", "openxlsx", "svglite", "readxl"
)
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) install.packages(missing, repos = "https://cloud.r-project.org")
invisible(lapply(required_pkgs, library, character.only = TRUE))

params <- list(
  spatial_rds = "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Results/network_spatial_relations/04_Logs/network_spatial_relations_objects.rds",
  module_file = "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Results/module_scores/curated_neuropil_modules.xlsx",
  output_dir = "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Results/module_spatial_networks",
  min_module_size = 5,
  min_abs_module_similarity = 0.50,
  group_levels = c("CON", "RES", "SUS"),
  run_group_specific_profiles = TRUE
)

message2 <- function(...) message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ...)
safe_name <- function(x) stringr::str_replace_all(as.character(x), "[^A-Za-z0-9_\\-]+", "_")

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

normalise_id <- function(x) {
  x %>%
    as.character() %>%
    str_replace(";.*$", "") %>%
    str_replace("\\|.*$", "") %>%
    str_replace("_MOUSE$", "") %>%
    str_to_lower()
}

row_zscore <- function(mat) {
  m <- rowMeans(mat, na.rm = TRUE)
  s <- apply(mat, 1, sd, na.rm = TRUE)
  s[is.na(s) | s == 0] <- 1
  sweep(sweep(mat, 1, m, "-"), 1, s, "/")
}

load_module_file <- function(module_file) {
  if (is.null(module_file) || !file.exists(module_file)) return(NULL)
  message2("Reading curated module file: ", module_file)
  ext <- tools::file_ext(module_file)
  if (tolower(ext) %in% c("xlsx", "xls")) {
    df <- readxl::read_excel(module_file)
  } else {
    df <- read.csv(module_file, stringsAsFactors = FALSE)
  }
  df <- as.data.frame(df)
  nms <- names(df)
  module_col <- nms[tolower(gsub("[^a-z0-9]", "", nms)) %in% c("module", "modulename", "set", "geneset")][1]
  protein_col <- nms[tolower(gsub("[^a-z0-9]", "", nms)) %in% c("protein", "gene", "genes", "accession", "uniprot", "uniprotid")][1]
  if (is.na(module_col) || is.na(protein_col)) {
    warning("Module file found but required columns Module and Protein/Gene could not be detected. Using regex fallback.")
    return(NULL)
  }
  df %>%
    transmute(Module = as.character(.data[[module_col]]), Protein = as.character(.data[[protein_col]])) %>%
    filter(!is.na(Module), !is.na(Protein), Module != "", Protein != "") %>%
    mutate(Module = safe_name(Module), ProteinNorm = normalise_id(Protein)) %>%
    distinct(Module, ProteinNorm, .keep_all = TRUE)
}

regex_modules_from_proteins <- function(proteins) {
  prot_norm <- normalise_id(proteins)
  tibble(Protein = proteins, ProteinNorm = prot_norm) %>%
    mutate(
      RNP_RNA_processing = str_detect(ProteinNorm, paste(c(
        "^hnrnp", "^snrnp", "^srsf", "^sf3", "^u2af", "^prpf", "^ddx", "^dhx",
        "^rbm", "^nono", "^pcbp", "^fus", "^matrin", "^pabp", "^pabpc", "^pabpn",
        "^ncl", "^nol", "^nop", "^fbl"
      ), collapse = "|")),
      Ribosome_translation = str_detect(ProteinNorm, paste(c(
        "^rps", "^rpl", "^eif", "^eef", "^mrpl", "^mrps"
      ), collapse = "|")),
      Mito_bioenergetics = str_detect(ProteinNorm, paste(c(
        "^nduf", "^cox", "^uqcr", "^atp5", "^sdh", "^cyc1", "^mt-", "^mdh", "^idha", "^id hb", "^ogdh"
      ), collapse = "|")),
      Endolysosomal_proteostasis = str_detect(ProteinNorm, paste(c(
        "^lamp", "^ct s", "^ctsb", "^ctsd", "^ctsl", "^psm", "^ub", "^hspa", "^hspb", "^vps", "^rab", "^atg"
      ), collapse = "|"))
    ) %>%
    pivot_longer(cols = c(RNP_RNA_processing, Ribosome_translation, Mito_bioenergetics, Endolysosomal_proteostasis), names_to = "Module", values_to = "Member") %>%
    filter(Member) %>%
    select(Module, Protein, ProteinNorm) %>%
    distinct()
}

match_modules_to_expression <- function(modules, expr) {
  expr_norm <- tibble(ProteinExpr = rownames(expr), ProteinNorm = normalise_id(rownames(expr)))
  modules %>%
    left_join(expr_norm, by = "ProteinNorm") %>%
    filter(!is.na(ProteinExpr)) %>%
    distinct(Module, ProteinExpr, .keep_all = TRUE)
}

compute_module_scores <- function(expr, sample_md, module_map, min_module_size = 5) {
  module_sizes <- module_map %>% count(Module, name = "NMatched")
  keep_modules <- module_sizes %>% filter(NMatched >= min_module_size) %>% pull(Module)
  module_map <- module_map %>% filter(Module %in% keep_modules)
  if (nrow(module_map) == 0) stop("No modules retained after min_module_size filter.")

  expr_z <- row_zscore(expr)

  score_mat <- sapply(split(module_map$ProteinExpr, module_map$Module), function(prots) {
    colMeans(expr_z[intersect(prots, rownames(expr_z)), , drop = FALSE], na.rm = TRUE)
  })

  score_df <- as.data.frame(score_mat) %>%
    rownames_to_column("SampleColumn") %>%
    left_join(sample_md %>% select(SampleColumn, Region, Layer, RegionLayer, ExpGroup), by = "SampleColumn") %>%
    pivot_longer(cols = all_of(colnames(score_mat)), names_to = "Module", values_to = "ModuleScore")

  list(score_df = score_df, module_sizes = module_sizes, retained_module_map = module_map)
}

aggregate_module_region_layer <- function(score_df) {
  score_df %>%
    filter(!is.na(RegionLayer), !is.na(ModuleScore)) %>%
    group_by(Module, RegionLayer) %>%
    summarise(MeanModuleScore = mean(ModuleScore, na.rm = TRUE), SEModuleScore = sd(ModuleScore, na.rm = TRUE) / sqrt(n()), N = n(), .groups = "drop")
}

make_module_matrix <- function(module_region_df) {
  module_region_df %>%
    select(Module, RegionLayer, MeanModuleScore) %>%
    pivot_wider(names_from = RegionLayer, values_from = MeanModuleScore) %>%
    column_to_rownames("Module") %>%
    as.matrix()
}

compute_region_similarity_from_modules <- function(module_matrix, min_abs_r = 0.5) {
  cor_mat <- suppressWarnings(cor(module_matrix, method = "spearman", use = "pairwise.complete.obs"))
  units <- colnames(cor_mat)
  edges <- expand.grid(Source = units, Target = units, stringsAsFactors = FALSE) %>%
    as_tibble() %>%
    filter(Source < Target) %>%
    rowwise() %>%
    mutate(ModuleProfileR = cor_mat[Source, Target]) %>%
    ungroup() %>%
    mutate(AbsModuleProfileR = abs(ModuleProfileR), Direction = ifelse(ModuleProfileR >= 0, "positive", "negative")) %>%
    arrange(desc(AbsModuleProfileR))
  list(cor_mat = cor_mat, all_edges = edges, filtered_edges = edges %>% filter(AbsModuleProfileR >= min_abs_r))
}

plot_module_heatmap <- function(module_matrix, outfile) {
  pheatmap::pheatmap(
    module_matrix,
    scale = "row",
    color = colorRampPalette(c("#457B9D", "white", "#E63946"))(101),
    border_color = NA,
    fontsize = 7,
    main = "Module spatial profiles",
    filename = outfile,
    width = 7.2,
    height = 4.8
  )
}

plot_module_network <- function(edges, outfile, title) {
  if (nrow(edges) == 0) {
    warning("No module-profile edges passed threshold for ", title)
    return(invisible(NULL))
  }
  nodes <- unique(c(edges$Source, edges$Target)) %>% tibble(name = .)
  g <- igraph::graph_from_data_frame(edges %>% transmute(from = Source, to = Target, AbsModuleProfileR, ModuleProfileR, Direction), directed = FALSE, vertices = nodes)
  p <- ggraph::ggraph(g, layout = "fr") +
    ggraph::geom_edge_link(aes(width = AbsModuleProfileR, alpha = AbsModuleProfileR), colour = "grey35") +
    ggraph::geom_node_point(shape = 21, fill = "white", colour = "black", size = 4, stroke = 0.4) +
    ggraph::geom_node_text(aes(label = name), repel = TRUE, size = 2.8) +
    scale_edge_width(range = c(0.2, 1.6), guide = "none") +
    scale_edge_alpha(range = c(0.25, 0.9), guide = "none") +
    labs(title = title) +
    theme_void(base_size = 8) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(outfile, p, width = 150, height = 120, units = "mm", device = svglite::svglite)
  invisible(p)
}

plot_group_module_profiles <- function(group_df, outfile) {
  p <- ggplot(group_df, aes(x = RegionLayer, y = MeanModuleScore, group = ExpGroup, shape = ExpGroup)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey70") +
    geom_point(size = 1.8, position = position_dodge(width = 0.45)) +
    geom_errorbar(aes(ymin = MeanModuleScore - SEModuleScore, ymax = MeanModuleScore + SEModuleScore), width = 0.15, linewidth = 0.25, position = position_dodge(width = 0.45)) +
    facet_wrap(~ Module, scales = "free_y", ncol = 2) +
    labs(x = NULL, y = "Module score", title = "Group-specific spatial module profiles") +
    theme_classic(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"), axis.text.y = element_text(colour = "black"), strip.background = element_blank(), strip.text = element_text(face = "bold"))
  ggsave(outfile, p, width = 180, height = 140, units = "mm", device = svglite::svglite)
  invisible(p)
}

write_graphml <- function(edges, outfile) {
  if (nrow(edges) == 0) return(invisible(NULL))
  nodes <- unique(c(edges$Source, edges$Target)) %>% tibble(name = .)
  g <- igraph::graph_from_data_frame(edges %>% rename(from = Source, to = Target), directed = FALSE, vertices = nodes)
  igraph::write_graph(g, outfile, format = "graphml")
}

# -------------------------------
# Main
# -------------------------------
dirs <- make_dirs(params$output_dir)
log_file <- file.path(dirs$logs, paste0("module_spatial_networks_run_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

message2("Starting module-level spatial network analysis")
if (!file.exists(params$spatial_rds)) stop("spatial_rds not found. Run network_spatial_relations.r first: ", params$spatial_rds)
obj <- readRDS(params$spatial_rds)
expr <- obj$expression_matrix
sample_md <- obj$sample_metadata

modules <- load_module_file(params$module_file)
if (is.null(modules)) {
  message2("Using regex-based fallback module definitions.")
  modules <- regex_modules_from_proteins(rownames(expr))
}

matched_modules <- match_modules_to_expression(modules, expr)
utils::write.csv(matched_modules, file.path(dirs$tables, "matched_module_membership.csv"), row.names = FALSE)

scores <- compute_module_scores(expr, sample_md, matched_modules, params$min_module_size)
score_df <- scores$score_df
utils::write.csv(score_df, file.path(dirs$tables, "module_scores_per_sample.csv"), row.names = FALSE)
utils::write.csv(scores$module_sizes, file.path(dirs$tables, "module_sizes_before_filter.csv"), row.names = FALSE)
utils::write.csv(scores$retained_module_map, file.path(dirs$tables, "retained_module_membership.csv"), row.names = FALSE)

module_region <- aggregate_module_region_layer(score_df)
utils::write.csv(module_region, file.path(dirs$tables, "module_scores_by_region_layer.csv"), row.names = FALSE)
module_matrix <- make_module_matrix(module_region)
utils::write.csv(as.data.frame(module_matrix) %>% rownames_to_column("Module"), file.path(dirs$tables, "module_region_layer_matrix.csv"), row.names = FALSE)

plot_module_heatmap(module_matrix, file.path(dirs$figures, "module_spatial_profile_heatmap.pdf"))

sim <- compute_region_similarity_from_modules(module_matrix, params$min_abs_module_similarity)
utils::write.csv(sim$all_edges, file.path(dirs$tables, "module_profile_region_similarity_all_edges.csv"), row.names = FALSE)
utils::write.csv(sim$filtered_edges, file.path(dirs$tables, "module_profile_region_similarity_filtered_edges.csv"), row.names = FALSE)
utils::write.csv(as.data.frame(sim$cor_mat) %>% rownames_to_column("RegionLayer"), file.path(dirs$tables, "module_profile_region_similarity_matrix.csv"), row.names = FALSE)
plot_module_network(sim$filtered_edges, file.path(dirs$figures, "module_profile_region_similarity_network.svg"), "Region/layer similarity from module profiles")

utils::write.csv(sim$filtered_edges, file.path(dirs$networks, "module_profile_similarity_edges.csv"), row.names = FALSE)
utils::write.csv(unique(c(sim$filtered_edges$Source, sim$filtered_edges$Target)) %>% tibble(name = .), file.path(dirs$networks, "module_profile_similarity_nodes.csv"), row.names = FALSE)
write_graphml(sim$filtered_edges, file.path(dirs$networks, "module_profile_similarity.graphml"))

# Group-specific module spatial profiles.
group_module_region <- NULL
if (isTRUE(params$run_group_specific_profiles) && "ExpGroup" %in% names(score_df)) {
  group_module_region <- score_df %>%
    filter(!is.na(ExpGroup), ExpGroup %in% params$group_levels, !is.na(RegionLayer)) %>%
    group_by(ExpGroup, Module, RegionLayer) %>%
    summarise(MeanModuleScore = mean(ModuleScore, na.rm = TRUE), SEModuleScore = sd(ModuleScore, na.rm = TRUE) / sqrt(n()), N = n(), .groups = "drop") %>%
    mutate(ExpGroup = factor(ExpGroup, levels = params$group_levels))
  utils::write.csv(group_module_region, file.path(dirs$tables, "group_specific_module_scores_by_region_layer.csv"), row.names = FALSE)
  plot_group_module_profiles(group_module_region, file.path(dirs$figures, "group_specific_module_spatial_profiles.svg"))
}

wb <- openxlsx::createWorkbook()
for (nm in c("Matched_Modules", "Scores_Per_Sample", "RegionLayer", "SimilarityEdges", "SimilarityFiltered")) openxlsx::addWorksheet(wb, nm)
openxlsx::writeData(wb, "Matched_Modules", matched_modules)
openxlsx::writeData(wb, "Scores_Per_Sample", score_df)
openxlsx::writeData(wb, "RegionLayer", module_region)
openxlsx::writeData(wb, "SimilarityEdges", sim$all_edges)
openxlsx::writeData(wb, "SimilarityFiltered", sim$filtered_edges)
if (!is.null(group_module_region)) {
  openxlsx::addWorksheet(wb, "Group_Profiles")
  openxlsx::writeData(wb, "Group_Profiles", group_module_region)
}
openxlsx::saveWorkbook(wb, file.path(dirs$tables, "module_spatial_networks_summary.xlsx"), overwrite = TRUE)

saveRDS(
  list(params = params, matched_modules = matched_modules, scores = score_df, module_region = module_region, module_matrix = module_matrix, similarity = sim, group_profiles = group_module_region, sessionInfo = sessionInfo()),
  file.path(dirs$logs, "module_spatial_networks_objects.rds")
)

message2("Finished module-level spatial network analysis")
message2("Caution: regex fallback modules are exploratory. Prefer curated module_file for publication analyses.")

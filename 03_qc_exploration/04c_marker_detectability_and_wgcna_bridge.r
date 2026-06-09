# ================================================================
# Script: 03_qc_exploration/04c_marker_detectability_and_wgcna_bridge.r
# Stage: qc
# Scope: dataset_specific
# Consumes: required data/processed/02_id_mapping/mapped/<dataset>/forward/per_file/*.csv; optional config/marker_panels/wgcna_reference_marker_sets.csv; results/tables/03_qc_exploration/05_empirical_roi_marker_discovery/empirical_roi_marker_sets.csv; +1 more.
# Produces: results/tables/03_qc_exploration/04c_marker_detectability_and_wgcna_bridge/<dataset>/.
# Dataset behavior: runs for neuron_neuropil,neuron_soma,microglia according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: QC bridge into WGCNA marker interpretation; optional WGCNA state is used when present.
# ================================================================

# Marker-panel detectability and optional WGCNA marker-enrichment bridge.
# This script is QC/interpreter support, not a formal purity or deconvolution analysis.
# It reports individual marker-panel scores, broad marker-compartment scores,
# and a narrow soma/neuropil/microglia marker-fidelity layer.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "dataset_inputs.R"))
source(repo_path("R", "qc_exploration_utils.R"))

run <- qc_args()
DATASET <- run$dataset
SUBSTEP_ID <- "04c_marker_detectability_and_wgcna_bridge"
PATHS <- qc_paths(SUBSTEP_ID, DATASET)
matrix_file <- qc_resolve_matrix(DATASET, env = "PROTEOMICS_MARKER_DETECTABILITY_MATRIX_FILE")
metadata_file <- qc_resolve_metadata(DATASET, env = "PROTEOMICS_MARKER_DETECTABILITY_METADATA_FILE")
marker_file <- Sys.getenv(
  "PROTEOMICS_MARKER_PANEL_FILE",
  unset = repo_path("config", "marker_panels", "wgcna_reference_marker_sets.csv")
)

if (run$dry_run) {
  status <- qc_dry_run_contract(
    "03_qc_exploration/04c_marker_detectability_and_wgcna_bridge.r",
    DATASET,
    matrix_file = matrix_file,
    metadata_file = metadata_file,
    paths = PATHS,
    extra = c(
      paste0("Marker panel file: ", marker_file, if (file.exists(marker_file)) " [PASS]" else " [fallback internal panels]"),
      "Writes marker-panel, marker-compartment, and soma/neuropil/microglia fidelity detectability tables and SVG figures.",
      "Optionally writes WGCNA module x marker-panel, marker-compartment, and fidelity-class enrichment if module assignments exist."
    )
  )
  quit(status = status, save = "no")
}

required_pkgs <- c("dplyr", "tidyr", "tibble", "ggplot2", "svglite", "scales")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))

if (!file.exists(matrix_file)) {
  stop("Marker-detectability matrix not found for dataset '", DATASET, "': ", matrix_file,
       ". Set PROTEOMICS_MARKER_DETECTABILITY_MATRIX_FILE or PROTEOMICS_QC_MATRIX_FILE.", call. = FALSE)
}

expr <- qc_read_expression(matrix_file, metadata_file, DATASET)
mat <- expr$mat
meta <- expr$meta

# -----------------------------------------------------------------------------
# Marker registry handling
# -----------------------------------------------------------------------------

gene_norm <- function(x) toupper(gsub("[^A-Za-z0-9]", "", sub("_MOUSE$", "", as.character(x), ignore.case = TRUE)))
collapse_unique <- function(x) paste(sort(unique(x[!is.na(x) & nzchar(x)])), collapse = ";")
first_nonempty <- function(x) {
  x <- unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))
  if (length(x)) x[[1]] else NA_character_
}

classify_marker_compartment <- function(marker_panel, cell_class = NA_character_, cell_state = NA_character_) {
  panel_l <- tolower(as.character(marker_panel))
  class_l <- tolower(as.character(cell_class))
  state_l <- tolower(as.character(cell_state))
  all_l <- paste(panel_l, class_l, state_l)

  dplyr::case_when(
    grepl("peripheral_myeloid|myeloid|monocyte|macrophage", all_l) ~ "myeloid_caution",

    grepl("microglia|microglia_pvm", all_l) ~ "microglia_core_or_pvm",

    grepl("neuron|neuronal|excitatory|inhibitory|interneuron", all_l) &
      grepl("soma|nuclear|histone", all_l) ~ "neuronal_soma_nuclear",

    grepl("neuron|neuronal|excitatory|inhibitory|interneuron", all_l) ~ "neuronal_neuropil_synaptic",

    grepl("astrocyte", all_l) ~ "astrocyte",

    grepl("opc", all_l) ~ "opc",

    grepl("oligodendrocyte|myelin", all_l) ~ "oligodendrocyte_myelin",

    grepl("endothelial|pericyte|vascular", all_l) ~ "vascular",

    grepl("extracellular|ecm|basement|collagen|laminin", all_l) ~ "ecm_basement_membrane",

    grepl("mitochondrial|oxphos", all_l) ~ "mitochondrial_oxphos",

    grepl("ribosomal|translation", all_l) ~ "ribosomal_translation",

    grepl("rnp|rna", all_l) ~ "rnp_rna_processing",

    TRUE ~ "other_unclassified"
  )
}

# Narrow marker set used for the cross-compartment fidelity QC question:
# Do neuron_soma, neuron_neuropil, and microglia samples show the expected
# relative abundance of soma, neuropil, and microglia marker proteins?
#
# This intentionally excludes broad Allen neuronal cell-class reference panels
# such as reference_cortical_excitatory_neuron, because those indicate neuronal
# identity/background rather than soma-vs-neuropil spatial compartment identity.
classify_fidelity_marker_class <- function(marker_panel) {
  panel_l <- tolower(as.character(marker_panel))

  dplyr::case_when(
    panel_l %in% c("canonical_neuronal_soma_nuclear", "nuclear_soma") ~ "Soma markers",
    panel_l %in% c("canonical_neuronal_synaptic_neuropil", "neuronal_synaptic_neuropil") ~ "Neuropil markers",
    panel_l %in% c(
      "canonical_microglia_homeostatic",
      "canonical_microglia_phagolysosomal_state",
      "reference_microglia_pvm",
      "microglia_homeostatic",
      "microglia_phagolysosomal"
    ) ~ "Microglia/PVM markers",
    TRUE ~ NA_character_
  )
}

fidelity_marker_class_order <- c("Soma markers", "Neuropil markers", "Microglia/PVM markers")

internal_marker_registry <- function() {
  panels <- list(
    neuronal_synaptic_neuropil = list(
      cell_class = "neuron", cell_state = "synaptic_neuropil",
      genes = c("Stxbp1", "Gpm6a", "Nptn", "Sh3gl2", "Atp6v1g2", "Snap25", "Syn1", "Syp", "Dlg4", "Camk2a")
    ),
    nuclear_soma = list(
      cell_class = "neuron", cell_state = "soma_nuclear",
      genes = c("H2ac1", "H4c1", "H3-3a", "H1-4", "H1-3", "Matr3", "Srsf3", "Ddx39b")
    ),
    microglia_homeostatic = list(
      cell_class = "microglia", cell_state = "homeostatic_identity",
      genes = c("Aif1", "Tmem119", "P2ry12", "Cx3cr1", "Csf1r", "Hexb", "Fcrls", "Gpr34", "Mertk")
    ),
    microglia_phagolysosomal = list(
      cell_class = "microglia", cell_state = "phagolysosomal_complement_state",
      genes = c("Tyrobp", "Trem2", "Apoe", "Lpl", "Cst7", "Ctsb", "Ctsd", "Lgals3", "Itgax", "Axl", "C1qa", "C1qb", "C1qc")
    ),
    astrocyte = list(
      cell_class = "astrocyte", cell_state = "canonical",
      genes = c("Gfap", "Aqp4", "Aldh1l1", "Slc1a2", "Slc1a3", "Aldoc")
    ),
    oligodendrocyte_myelin = list(
      cell_class = "oligodendrocyte", cell_state = "myelin",
      genes = c("Mbp", "Mog", "Plp1", "Cnp", "Mag", "Mobp")
    ),
    endothelial_pericyte_vascular = list(
      cell_class = "vascular", cell_state = "endothelial_pericyte",
      genes = c("Pecam1", "Cldn5", "Kdr", "Rgs5", "Pdgfrb", "Vtn")
    ),
    extracellular_matrix_basement_membrane = list(
      cell_class = "ecm", cell_state = "basement_membrane",
      genes = c("Agrn", "Lamb2", "Lamc1", "Col4a1", "Col4a2", "Nid1", "Nid2", "Lama2", "Lama4", "Lama5", "Itga1", "Itgb1")
    ),
    mitochondrial_oxphos = list(
      cell_class = "mitochondrial", cell_state = "oxphos",
      genes = c("Ndufs1", "Ndufa9", "Sdha", "Uqcrc2", "Cox4i1", "Atp5f1a", "Atp5f1b")
    ),
    ribosomal_translation = list(
      cell_class = "ribosomal", cell_state = "translation",
      genes = c("Rpl3", "Rpl4", "Rpl5", "Rps3", "Rps6", "Eef1a1", "Eef2")
    ),
    rnp_rna_processing = list(
      cell_class = "rnp", cell_state = "rna_processing",
      genes = c("Hnrnpa2b1", "Hnrnpc", "Sfpq", "Snrnp70", "Ddx5", "Ddx17", "Pabpc1")
    )
  )

  dplyr::bind_rows(lapply(names(panels), function(panel) {
    x <- panels[[panel]]
    data.frame(
      marker_panel = panel,
      cell_class = x$cell_class,
      cell_state = x$cell_state,
      requested_marker = unique(x$genes),
      source_type = "internal_fallback",
      source_name = "04c_internal_marker_sets",
      confidence = "curated_conservative",
      use_for = "qc_interpretation",
      stringsAsFactors = FALSE
    )
  }))
}

read_marker_registry <- function(path) {
  if (!file.exists(path)) return(internal_marker_registry())
  df <- tryCatch(qc_read_table(path), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(internal_marker_registry())

  panel_col <- qc_first_col(
    df,
    c("marker_set", "marker_panel", "panel", "celltype", "cell_type", "signature", "set", "source_panel")
  )
  gene_col <- qc_first_col(df, c("gene", "gene_symbol", "marker", "symbol", "protein", "protein_id", "Gene"))
  if (is.na(panel_col) || is.na(gene_col)) return(internal_marker_registry())

  optional_col <- function(candidates) {
    col <- qc_first_col(df, candidates)
    if (is.na(col)) rep(NA_character_, nrow(df)) else as.character(df[[col]])
  }

  registry <- data.frame(
    marker_panel = as.character(df[[panel_col]]),
    cell_class = optional_col(c("cell_class", "celltype_class", "class", "cell_type_class")),
    cell_state = optional_col(c("cell_state", "state", "subclass", "cell_subclass", "signature_state")),
    registry_fidelity_marker_class = optional_col(c("fidelity_marker_class", "marker_class")),
    requested_marker = as.character(df[[gene_col]]),
    source_type = optional_col(c("source_type", "source_category")),
    source_name = optional_col(c("source_name", "source", "reference_source")),
    source_reference = optional_col(c("source_reference", "reference", "url", "citation")),
    source_term_or_category = optional_col(c("source_term_or_category", "source_term", "category", "term")),
    evidence_level = optional_col(c("evidence_level", "evidence_code")),
    confidence = optional_col(c("confidence", "evidence")),
    use_for = optional_col(c("use_for", "intended_use")),
    fidelity_subpanel = optional_col(c("fidelity_subpanel", "subpanel")),
    include_in_fidelity_score = optional_col(c("include_in_fidelity_score")),
    marker_notes = optional_col(c("notes", "marker_notes")),
    stringsAsFactors = FALSE
  )

  registry |>
    dplyr::filter(!is.na(marker_panel), nzchar(marker_panel), !is.na(requested_marker), nzchar(requested_marker)) |>
    dplyr::mutate(
      marker_key = gene_norm(requested_marker),
      marker_compartment = classify_marker_compartment(marker_panel, cell_class, cell_state),
      registry_fidelity_marker_class = dplyr::na_if(registry_fidelity_marker_class, ""),
      fidelity_marker_class = dplyr::coalesce(registry_fidelity_marker_class, classify_fidelity_marker_class(marker_panel))
    ) |>
    dplyr::select(-"registry_fidelity_marker_class") |>
    dplyr::filter(!is.na(marker_key), nzchar(marker_key)) |>
    dplyr::distinct(marker_panel, requested_marker, marker_key, .keep_all = TRUE)
}

marker_lookup <- read_marker_registry(marker_file)
marker_sets <- split(marker_lookup$requested_marker, marker_lookup$marker_panel)
protein_ids <- rownames(mat)
protein_key <- gene_norm(protein_ids)
universe <- unique(protein_key[!is.na(protein_key) & nzchar(protein_key)])

panel_metadata <- marker_lookup |>
  dplyr::group_by(marker_panel) |>
  dplyr::summarise(
    marker_compartment = first_nonempty(marker_compartment),
    fidelity_marker_class = first_nonempty(fidelity_marker_class),
    cell_class = first_nonempty(cell_class),
    cell_state = first_nonempty(cell_state),
    source_type = first_nonempty(source_type),
    source_name = first_nonempty(source_name),
    source_reference = first_nonempty(source_reference),
    source_term_or_category = first_nonempty(source_term_or_category),
    evidence_level = first_nonempty(evidence_level),
    confidence = first_nonempty(confidence),
    use_for = first_nonempty(use_for),
    fidelity_subpanel = first_nonempty(fidelity_subpanel),
    include_in_fidelity_score = first_nonempty(include_in_fidelity_score),
    marker_notes = first_nonempty(marker_notes),
    .groups = "drop"
  )

# -----------------------------------------------------------------------------
# Dataset abundance/rank summaries
# -----------------------------------------------------------------------------

long <- as.data.frame(mat, check.names = FALSE) |>
  tibble::rownames_to_column("Protein") |>
  tidyr::pivot_longer(-Protein, names_to = "Sample", values_to = "Log2Intensity") |>
  dplyr::left_join(meta, by = "Sample")

rank_group_cols <- intersect(c("Region", "region", "Layer", "layer", "Group", "group", "ExpGroup", "Sex", "sex", "plate", "batch"), names(meta))
rank_group_cols <- rank_group_cols[!duplicated(tolower(rank_group_cols))]
if (!length(rank_group_cols)) rank_group_cols <- "Sample"
rank_data <- long |>
  dplyr::mutate(RankGroup = do.call(paste, c(dplyr::across(dplyr::all_of(rank_group_cols)), sep = " | "))) |>
  dplyr::filter(!is.na(Log2Intensity), !is.na(RankGroup), nzchar(RankGroup)) |>
  dplyr::group_by(RankGroup, Protein) |>
  dplyr::summarise(MeanLog2 = mean(Log2Intensity, na.rm = TRUE), .groups = "drop") |>
  dplyr::group_by(RankGroup) |>
  dplyr::arrange(dplyr::desc(MeanLog2), .by_group = TRUE) |>
  dplyr::mutate(Rank = dplyr::row_number()) |>
  dplyr::ungroup() |>
  dplyr::mutate(protein_key = gene_norm(Protein))

protein_stats <- data.frame(Protein = protein_ids, protein_key = protein_key, stringsAsFactors = FALSE) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    n_samples = ncol(mat),
    n_nonmissing = sum(is.finite(mat[Protein, ])),
    fraction_nonmissing = n_nonmissing / n_samples,
    mean_log2_abundance = ifelse(n_nonmissing > 0, mean(mat[Protein, ], na.rm = TRUE), NA_real_),
    median_log2_abundance = ifelse(n_nonmissing > 0, median(mat[Protein, ], na.rm = TRUE), NA_real_),
    sd_log2_abundance = ifelse(n_nonmissing > 1, stats::sd(mat[Protein, ], na.rm = TRUE), NA_real_),
    cv_like = ifelse(is.finite(sd_log2_abundance) && is.finite(mean_log2_abundance) && abs(mean_log2_abundance) > .Machine$double.eps,
                     sd_log2_abundance / abs(mean_log2_abundance), NA_real_)
  ) |>
  dplyr::ungroup()

rank_stats <- rank_data |>
  dplyr::group_by(protein_key) |>
  dplyr::summarise(mean_rank = mean(Rank, na.rm = TRUE), median_rank = median(Rank, na.rm = TRUE), .groups = "drop")

# -----------------------------------------------------------------------------
# Marker-panel detectability
# -----------------------------------------------------------------------------

marker_detectability_by_protein <- marker_lookup |>
  dplyr::left_join(protein_stats, by = c("marker_key" = "protein_key")) |>
  dplyr::group_by(
    dataset = DATASET,
    marker_panel,
    marker_compartment,
    fidelity_marker_class,
    fidelity_subpanel,
    cell_class,
    cell_state,
    source_type,
    source_name,
    source_reference,
    source_term_or_category,
    evidence_level,
    confidence,
    use_for,
    include_in_fidelity_score,
    marker_notes,
    requested_marker,
    marker_key
  ) |>
  dplyr::summarise(
    matched_protein_id = collapse_unique(Protein),
    detected = any(!is.na(Protein)),
    n_matched_protein_ids = sum(!is.na(Protein)),
    n_samples = dplyr::first(n_samples[!is.na(n_samples)]) %||% ncol(mat),
    n_nonmissing = ifelse(all(is.na(n_nonmissing)), 0L, sum(n_nonmissing, na.rm = TRUE)),
    fraction_nonmissing = ifelse(n_matched_protein_ids > 0, max(fraction_nonmissing, na.rm = TRUE), 0),
    mean_log2_abundance = ifelse(n_matched_protein_ids > 0, mean(mean_log2_abundance, na.rm = TRUE), NA_real_),
    median_log2_abundance = ifelse(n_matched_protein_ids > 0, median(median_log2_abundance, na.rm = TRUE), NA_real_),
    sd_log2_abundance = ifelse(n_matched_protein_ids > 0, mean(sd_log2_abundance, na.rm = TRUE), NA_real_),
    cv_like = ifelse(n_matched_protein_ids > 0, mean(cv_like, na.rm = TRUE), NA_real_),
    .groups = "drop"
  ) |>
  dplyr::left_join(rank_stats, by = c("marker_key" = "protein_key")) |>
  dplyr::mutate(
    matched_protein_id = ifelse(nzchar(matched_protein_id), matched_protein_id, NA_character_),
    fraction_nonmissing = ifelse(is.finite(fraction_nonmissing), fraction_nonmissing, 0),
    mean_log2_abundance = ifelse(is.nan(mean_log2_abundance), NA_real_, mean_log2_abundance),
    median_log2_abundance = ifelse(is.nan(median_log2_abundance), NA_real_, median_log2_abundance),
    sd_log2_abundance = ifelse(is.nan(sd_log2_abundance), NA_real_, sd_log2_abundance),
    cv_like = ifelse(is.nan(cv_like), NA_real_, cv_like)
  )

sample_scores <- lapply(names(marker_sets), function(panel) {
  panel_meta <- panel_metadata[panel_metadata$marker_panel == panel, , drop = FALSE]
  genes <- gene_norm(marker_sets[[panel]])
  idx <- which(protein_key %in% genes)
  if (!length(idx)) {
    out <- data.frame(
      dataset = DATASET,
      Sample = colnames(mat),
      marker_panel = panel,
      marker_score = NA_real_,
      n_detected_markers = 0L,
      fraction_detected_markers = 0,
      stringsAsFactors = FALSE
    )
  } else {
    out <- data.frame(
      dataset = DATASET,
      Sample = colnames(mat),
      marker_panel = panel,
      marker_score = colMeans(mat[idx, , drop = FALSE], na.rm = TRUE),
      n_detected_markers = length(unique(protein_key[idx])),
      fraction_detected_markers = length(unique(protein_key[idx])) / length(unique(genes)),
      stringsAsFactors = FALSE
    )
  }
  dplyr::left_join(out, panel_meta, by = "marker_panel")
}) |>
  dplyr::bind_rows() |>
  dplyr::left_join(meta, by = "Sample")

panel_score_summary <- sample_scores |>
  dplyr::group_by(dataset, marker_panel) |>
  dplyr::summarise(
    mean_marker_score = mean(marker_score, na.rm = TRUE),
    median_marker_score = median(marker_score, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    mean_marker_score = ifelse(is.nan(mean_marker_score), NA_real_, mean_marker_score),
    median_marker_score = ifelse(is.nan(median_marker_score), NA_real_, median_marker_score)
  )

marker_detectability_by_panel <- marker_detectability_by_protein |>
  dplyr::group_by(dataset, marker_panel, marker_compartment, fidelity_marker_class, fidelity_subpanel, cell_class, cell_state, source_type, source_name, source_reference, evidence_level, confidence, use_for, include_in_fidelity_score) |>
  dplyr::summarise(
    n_markers_requested = dplyr::n_distinct(requested_marker),
    n_marker_keys_requested = dplyr::n_distinct(marker_key),
    n_markers_detected = sum(detected),
    fraction_markers_detected = n_markers_detected / n_markers_requested,
    median_fraction_nonmissing = median(fraction_nonmissing, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::left_join(panel_score_summary, by = c("dataset", "marker_panel"))

marker_panels_not_detected <- marker_detectability_by_panel |>
  dplyr::filter(n_markers_detected == 0) |>
  dplyr::arrange(marker_compartment, marker_panel)

# -----------------------------------------------------------------------------
# Marker-compartment detectability using deduplicated marker proteins
# -----------------------------------------------------------------------------

compartment_keys <- marker_lookup |>
  dplyr::filter(!is.na(marker_compartment), nzchar(marker_compartment)) |>
  dplyr::group_by(marker_compartment) |>
  dplyr::summarise(
    marker_panels = collapse_unique(marker_panel),
    requested_markers = collapse_unique(requested_marker),
    requested_keys = list(unique(marker_key)),
    .groups = "drop"
  )

sample_compartment_scores <- lapply(seq_len(nrow(compartment_keys)), function(i) {
  comp <- compartment_keys$marker_compartment[[i]]
  genes <- unique(compartment_keys$requested_keys[[i]])
  idx <- which(protein_key %in% genes)
  if (!length(idx)) {
    data.frame(
      dataset = DATASET,
      Sample = colnames(mat),
      marker_compartment = comp,
      marker_compartment_score = NA_real_,
      n_detected_markers = 0L,
      fraction_detected_markers = 0,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      dataset = DATASET,
      Sample = colnames(mat),
      marker_compartment = comp,
      marker_compartment_score = colMeans(mat[idx, , drop = FALSE], na.rm = TRUE),
      n_detected_markers = length(unique(protein_key[idx])),
      fraction_detected_markers = length(unique(protein_key[idx])) / length(unique(genes)),
      stringsAsFactors = FALSE
    )
  }
}) |>
  dplyr::bind_rows() |>
  dplyr::left_join(meta, by = "Sample")

compartment_score_summary <- sample_compartment_scores |>
  dplyr::group_by(dataset, marker_compartment) |>
  dplyr::summarise(
    mean_marker_compartment_score = mean(marker_compartment_score, na.rm = TRUE),
    median_marker_compartment_score = median(marker_compartment_score, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    mean_marker_compartment_score = ifelse(is.nan(mean_marker_compartment_score), NA_real_, mean_marker_compartment_score),
    median_marker_compartment_score = ifelse(is.nan(median_marker_compartment_score), NA_real_, median_marker_compartment_score)
  )

marker_detectability_by_compartment <- marker_detectability_by_protein |>
  dplyr::group_by(dataset, marker_compartment) |>
  dplyr::summarise(
    n_marker_panels = dplyr::n_distinct(marker_panel),
    marker_panels = collapse_unique(marker_panel),
    n_markers_requested = dplyr::n_distinct(marker_key),
    n_markers_detected = dplyr::n_distinct(marker_key[detected]),
    fraction_markers_detected = n_markers_detected / n_markers_requested,
    median_fraction_nonmissing = median(fraction_nonmissing, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::left_join(compartment_score_summary, by = c("dataset", "marker_compartment")) |>
  dplyr::arrange(marker_compartment)

# -----------------------------------------------------------------------------
# Narrow soma/neuropil/microglia marker-fidelity layer
# -----------------------------------------------------------------------------

fidelity_lookup <- marker_lookup |>
  dplyr::filter(!is.na(fidelity_marker_class), nzchar(fidelity_marker_class)) |>
  dplyr::mutate(
    fidelity_marker_class = factor(fidelity_marker_class, levels = fidelity_marker_class_order, ordered = TRUE)
  )

if (nrow(fidelity_lookup)) {
  fidelity_keys <- fidelity_lookup |>
    dplyr::group_by(fidelity_marker_class) |>
    dplyr::summarise(
      marker_panels = collapse_unique(marker_panel),
      requested_markers = collapse_unique(requested_marker),
      requested_keys = list(unique(marker_key)),
      .groups = "drop"
    ) |>
    dplyr::arrange(fidelity_marker_class)

  fidelity_sample_scores <- lapply(seq_len(nrow(fidelity_keys)), function(i) {
    cls <- as.character(fidelity_keys$fidelity_marker_class[[i]])
    genes <- unique(fidelity_keys$requested_keys[[i]])
    idx <- which(protein_key %in% genes)

    if (!length(idx)) {
      data.frame(
        dataset = DATASET,
        Sample = colnames(mat),
        fidelity_marker_class = cls,
        fidelity_marker_score = NA_real_,
        n_detected_markers = 0L,
        fraction_detected_markers = 0,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        dataset = DATASET,
        Sample = colnames(mat),
        fidelity_marker_class = cls,
        fidelity_marker_score = colMeans(mat[idx, , drop = FALSE], na.rm = TRUE),
        n_detected_markers = length(unique(protein_key[idx])),
        fraction_detected_markers = length(unique(protein_key[idx])) / length(unique(genes)),
        stringsAsFactors = FALSE
      )
    }
  }) |>
    dplyr::bind_rows() |>
    dplyr::mutate(
      fidelity_marker_class = factor(fidelity_marker_class, levels = fidelity_marker_class_order, ordered = TRUE)
    ) |>
    dplyr::left_join(meta, by = "Sample")

  fidelity_score_summary <- fidelity_sample_scores |>
    dplyr::group_by(dataset, fidelity_marker_class) |>
    dplyr::summarise(
      mean_fidelity_marker_score = mean(fidelity_marker_score, na.rm = TRUE),
      median_fidelity_marker_score = median(fidelity_marker_score, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      mean_fidelity_marker_score = ifelse(is.nan(mean_fidelity_marker_score), NA_real_, mean_fidelity_marker_score),
      median_fidelity_marker_score = ifelse(is.nan(median_fidelity_marker_score), NA_real_, median_fidelity_marker_score)
    )

  marker_fidelity_by_class <- marker_detectability_by_protein |>
    dplyr::filter(!is.na(fidelity_marker_class), nzchar(fidelity_marker_class)) |>
    dplyr::mutate(
      fidelity_marker_class = factor(fidelity_marker_class, levels = fidelity_marker_class_order, ordered = TRUE)
    ) |>
    dplyr::group_by(dataset, fidelity_marker_class) |>
    dplyr::summarise(
      n_marker_panels = dplyr::n_distinct(marker_panel),
      marker_panels = collapse_unique(marker_panel),
      n_markers_requested = dplyr::n_distinct(marker_key),
      n_markers_detected = dplyr::n_distinct(marker_key[detected]),
      fraction_markers_detected = n_markers_detected / n_markers_requested,
      median_fraction_nonmissing = median(fraction_nonmissing, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::left_join(fidelity_score_summary, by = c("dataset", "fidelity_marker_class")) |>
    dplyr::arrange(fidelity_marker_class)
} else {
  fidelity_sample_scores <- data.frame(
    dataset = character(),
    Sample = character(),
    fidelity_marker_class = character(),
    fidelity_marker_score = numeric(),
    n_detected_markers = integer(),
    fraction_detected_markers = numeric(),
    stringsAsFactors = FALSE
  )

  marker_fidelity_by_class <- data.frame(
    dataset = character(),
    fidelity_marker_class = character(),
    n_marker_panels = integer(),
    marker_panels = character(),
    n_markers_requested = integer(),
    n_markers_detected = integer(),
    fraction_markers_detected = numeric(),
    median_fraction_nonmissing = numeric(),
    mean_fidelity_marker_score = numeric(),
    median_fidelity_marker_score = numeric(),
    stringsAsFactors = FALSE
  )
}


# Explicit protein-level provenance for the narrow fidelity layer. This is the
# audit trail for 04d: it shows exactly which requested marker genes matched
# proteins in the current dataset and therefore contributed to Soma/Neuropil/
# Microglia-PVM fidelity scores.
fidelity_marker_proteins_used <- marker_detectability_by_protein |>
  dplyr::filter(!is.na(fidelity_marker_class), nzchar(as.character(fidelity_marker_class))) |>
  dplyr::mutate(
    fidelity_marker_class = factor(fidelity_marker_class, levels = fidelity_marker_class_order, ordered = TRUE),
    used_in_fidelity_score = detected & !is.na(matched_protein_id) & nzchar(as.character(matched_protein_id))
  ) |>
  dplyr::transmute(
    dataset,
    fidelity_marker_class = as.character(fidelity_marker_class),
    marker_panel,
    marker_compartment,
    cell_class,
    cell_state,
    source_type,
    source_name,
    source_reference,
    source_term_or_category,
    evidence_level,
    confidence,
    use_for,
    fidelity_subpanel,
    include_in_fidelity_score,
    marker_notes,
    requested_marker,
    marker_key,
    matched_protein_id,
    used_in_fidelity_score,
    n_matched_protein_ids,
    n_samples,
    n_nonmissing,
    fraction_nonmissing,
    mean_log2_abundance,
    median_log2_abundance,
    sd_log2_abundance,
    cv_like,
    mean_rank,
    median_rank
  ) |>
  dplyr::arrange(factor(fidelity_marker_class, levels = fidelity_marker_class_order, ordered = TRUE), marker_panel, requested_marker)


# Deduplicated class-level audit matching the actual fidelity-score logic.
# The fidelity score is calculated from unique marker_key values within each
# fidelity_marker_class, so repeated markers across microglia/PVM panels do not
# overweight the score. This table exposes that exact scoring-level provenance:
# one row per dataset x fidelity class x marker_key.
fidelity_marker_proteins_used_deduplicated <- fidelity_marker_proteins_used |>
  dplyr::group_by(dataset, fidelity_marker_class, marker_key) |>
  dplyr::summarise(
    requested_markers = collapse_unique(requested_marker),
    matched_protein_id = collapse_unique(matched_protein_id),
    used_in_fidelity_score = any(used_in_fidelity_score, na.rm = TRUE),
    contributing_marker_panels = collapse_unique(marker_panel),
    contributing_fidelity_subpanels = collapse_unique(fidelity_subpanel),
    source_types = collapse_unique(source_type),
    source_names = collapse_unique(source_name),
    source_references = collapse_unique(source_reference),
    source_terms_or_categories = collapse_unique(source_term_or_category),
    evidence_levels = collapse_unique(evidence_level),
    confidences = collapse_unique(confidence),
    include_in_fidelity_score_values = collapse_unique(include_in_fidelity_score),
    n_marker_panels_containing_marker = dplyr::n_distinct(marker_panel),
    n_matched_protein_ids = max(n_matched_protein_ids, na.rm = TRUE),
    n_samples = dplyr::first(n_samples),
    n_nonmissing = max(n_nonmissing, na.rm = TRUE),
    fraction_nonmissing = max(fraction_nonmissing, na.rm = TRUE),
    mean_log2_abundance = mean(mean_log2_abundance, na.rm = TRUE),
    median_log2_abundance = median(median_log2_abundance, na.rm = TRUE),
    sd_log2_abundance = mean(sd_log2_abundance, na.rm = TRUE),
    cv_like = mean(cv_like, na.rm = TRUE),
    mean_rank = mean(mean_rank, na.rm = TRUE),
    median_rank = median(median_rank, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    n_matched_protein_ids = ifelse(is.infinite(n_matched_protein_ids), 0L, n_matched_protein_ids),
    n_nonmissing = ifelse(is.infinite(n_nonmissing), 0L, n_nonmissing),
    fraction_nonmissing = ifelse(is.infinite(fraction_nonmissing), 0, fraction_nonmissing),
    mean_log2_abundance = ifelse(is.nan(mean_log2_abundance), NA_real_, mean_log2_abundance),
    median_log2_abundance = ifelse(is.nan(median_log2_abundance), NA_real_, median_log2_abundance),
    sd_log2_abundance = ifelse(is.nan(sd_log2_abundance), NA_real_, sd_log2_abundance),
    cv_like = ifelse(is.nan(cv_like), NA_real_, cv_like),
    mean_rank = ifelse(is.nan(mean_rank), NA_real_, mean_rank),
    median_rank = ifelse(is.nan(median_rank), NA_real_, median_rank)
  ) |>
  dplyr::arrange(
    factor(fidelity_marker_class, levels = fidelity_marker_class_order, ordered = TRUE),
    marker_key
  )

fidelity_duplicate_markers <- fidelity_marker_proteins_used_deduplicated |>
  dplyr::filter(n_marker_panels_containing_marker > 1) |>
  dplyr::arrange(
    factor(fidelity_marker_class, levels = fidelity_marker_class_order, ordered = TRUE),
    marker_key
  )

marker_fidelity_protein_summary <- fidelity_marker_proteins_used |>
  dplyr::group_by(dataset, fidelity_marker_class, marker_panel) |>
  dplyr::summarise(
    source_types = collapse_unique(source_type),
    source_names = collapse_unique(source_name),
    source_references = collapse_unique(source_reference),
    evidence_levels = collapse_unique(evidence_level),
    fidelity_subpanels = collapse_unique(fidelity_subpanel),
    include_in_fidelity_score_values = collapse_unique(include_in_fidelity_score),
    n_requested_markers = dplyr::n_distinct(marker_key),
    n_used_markers = dplyr::n_distinct(marker_key[used_in_fidelity_score]),
    fraction_used_markers = ifelse(n_requested_markers > 0, n_used_markers / n_requested_markers, NA_real_),
    used_requested_markers = collapse_unique(requested_marker[used_in_fidelity_score]),
    used_matched_protein_ids = collapse_unique(matched_protein_id[used_in_fidelity_score]),
    missing_requested_markers = collapse_unique(requested_marker[!used_in_fidelity_score]),
    median_fraction_nonmissing_used = ifelse(any(used_in_fidelity_score), median(fraction_nonmissing[used_in_fidelity_score], na.rm = TRUE), NA_real_),
    mean_log2_abundance_used = ifelse(any(used_in_fidelity_score), mean(mean_log2_abundance[used_in_fidelity_score], na.rm = TRUE), NA_real_),
    median_log2_abundance_used = ifelse(any(used_in_fidelity_score), median(median_log2_abundance[used_in_fidelity_score], na.rm = TRUE), NA_real_),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    median_fraction_nonmissing_used = ifelse(is.nan(median_fraction_nonmissing_used), NA_real_, median_fraction_nonmissing_used),
    mean_log2_abundance_used = ifelse(is.nan(mean_log2_abundance_used), NA_real_, mean_log2_abundance_used),
    median_log2_abundance_used = ifelse(is.nan(median_log2_abundance_used), NA_real_, median_log2_abundance_used)
  ) |>
  dplyr::arrange(factor(fidelity_marker_class, levels = fidelity_marker_class_order, ordered = TRUE), marker_panel)

marker_source_detectability_summary <- marker_detectability_by_protein |>
  dplyr::group_by(dataset, source_name, source_type, confidence, evidence_level, use_for, include_in_fidelity_score) |>
  dplyr::summarise(
    n_marker_panels = dplyr::n_distinct(marker_panel),
    marker_panels = collapse_unique(marker_panel),
    marker_compartments = collapse_unique(marker_compartment),
    fidelity_marker_classes = collapse_unique(fidelity_marker_class),
    n_requested_markers = dplyr::n_distinct(marker_key),
    n_detected_markers = dplyr::n_distinct(marker_key[detected]),
    fraction_detected_markers = ifelse(n_requested_markers > 0, n_detected_markers / n_requested_markers, NA_real_),
    n_used_fidelity_markers = dplyr::n_distinct(marker_key[!is.na(fidelity_marker_class) & nzchar(as.character(fidelity_marker_class)) & detected]),
    source_terms_or_categories = collapse_unique(source_term_or_category),
    source_references = collapse_unique(source_reference),
    .groups = "drop"
  ) |>
  dplyr::arrange(source_name, source_type, confidence, evidence_level)

# -----------------------------------------------------------------------------
# Write tables
# -----------------------------------------------------------------------------

qc_write_csv(marker_detectability_by_protein, file.path(PATHS$tables, "marker_detectability_by_protein.csv"))
qc_write_csv(marker_detectability_by_panel, file.path(PATHS$tables, "marker_detectability_by_panel.csv"))
qc_write_csv(sample_scores, file.path(PATHS$tables, "marker_detectability_by_sample.csv"))
qc_write_csv(marker_detectability_by_compartment, file.path(PATHS$tables, "marker_detectability_by_compartment.csv"))
qc_write_csv(sample_compartment_scores, file.path(PATHS$tables, "marker_detectability_by_sample_compartment.csv"))
qc_write_csv(marker_panels_not_detected, file.path(PATHS$tables, "marker_panels_not_detected.csv"))
qc_write_csv(fidelity_sample_scores, file.path(PATHS$tables, "marker_fidelity_by_sample.csv"))
qc_write_csv(marker_fidelity_by_class, file.path(PATHS$tables, "marker_fidelity_by_class.csv"))
qc_write_csv(fidelity_marker_proteins_used, file.path(PATHS$tables, "marker_fidelity_proteins_used.csv"))
qc_write_csv(fidelity_marker_proteins_used_deduplicated, file.path(PATHS$tables, "marker_fidelity_proteins_used_deduplicated.csv"))
qc_write_csv(fidelity_duplicate_markers, file.path(PATHS$tables, "marker_fidelity_duplicate_marker_membership.csv"))
qc_write_csv(marker_fidelity_protein_summary, file.path(PATHS$tables, "marker_fidelity_protein_summary.csv"))
qc_write_csv(marker_source_detectability_summary, file.path(PATHS$tables, "marker_source_detectability_summary.csv"))
qc_write_xlsx(
  list(
    panel_sample = sample_scores,
    comp_sample = sample_compartment_scores,
    fidelity_sample = fidelity_sample_scores,
    protein_detect = marker_detectability_by_protein,
    panel_detect = marker_detectability_by_panel,
    comp_detect = marker_detectability_by_compartment,
    absent_panels = marker_panels_not_detected,
    fidelity_class = marker_fidelity_by_class,
    fidelity_proteins = fidelity_marker_proteins_used,
    fidelity_dedup = fidelity_marker_proteins_used_deduplicated,
    fidelity_dup = fidelity_duplicate_markers,
    fidelity_prot_sum = marker_fidelity_protein_summary,
    source_summary = marker_source_detectability_summary
  ),
  file.path(PATHS$tables, "marker_detectability_qc.xlsx")
)

# -----------------------------------------------------------------------------
# Figures
# -----------------------------------------------------------------------------

p_detection <- ggplot(marker_detectability_by_panel, aes(x = reorder(marker_panel, fraction_markers_detected), y = fraction_markers_detected)) +
  geom_col(fill = "grey55") +
  geom_text(aes(label = paste0(n_markers_detected, "/", n_markers_requested)), hjust = -0.05, size = 2) +
  coord_flip(ylim = c(0, 1.05)) +
  labs(x = NULL, y = "Detected marker fraction") +
  theme_classic(base_size = 8)
ggsave(file.path(PATHS$figures, "marker_panel_detection_barplot.svg"), p_detection, width = 140, height = 100, units = "mm", device = svglite::svglite)

p_comp_detection <- ggplot(marker_detectability_by_compartment, aes(x = reorder(marker_compartment, fraction_markers_detected), y = fraction_markers_detected)) +
  geom_col(fill = "grey55") +
  geom_text(aes(label = paste0(n_markers_detected, "/", n_markers_requested)), hjust = -0.05, size = 2) +
  coord_flip(ylim = c(0, 1.05)) +
  labs(x = NULL, y = "Detected marker fraction") +
  theme_classic(base_size = 8)
ggsave(file.path(PATHS$figures, "marker_compartment_detection_barplot.svg"), p_comp_detection, width = 140, height = 90, units = "mm", device = svglite::svglite)

if (nrow(marker_fidelity_by_class)) {
  p_fidelity_detection <- ggplot(marker_fidelity_by_class, aes(x = fidelity_marker_class, y = fraction_markers_detected)) +
    geom_col(fill = "grey55") +
    geom_text(aes(label = paste0(n_markers_detected, "/", n_markers_requested)), vjust = -0.25, size = 2) +
    scale_x_discrete(drop = FALSE) +
    coord_cartesian(ylim = c(0, 1.05)) +
    labs(x = NULL, y = "Detected marker fraction") +
    theme_classic(base_size = 8) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  ggsave(file.path(PATHS$figures, "marker_fidelity_detection_barplot.svg"), p_fidelity_detection, width = 110, height = 75, units = "mm", device = svglite::svglite)
}

plot_missing <- marker_detectability_by_protein |>
  dplyr::mutate(marker_label = ifelse(detected, requested_marker, paste0(requested_marker, " (not detected)")))
p_missing <- ggplot(plot_missing, aes(x = marker_panel, y = marker_label, size = fraction_nonmissing, fill = mean_log2_abundance)) +
  geom_point(shape = 21, color = "grey30", alpha = 0.85) +
  scale_size_continuous(range = c(0.4, 3.2), limits = c(0, 1)) +
  labs(x = NULL, y = NULL, size = "Nonmissing fraction", fill = "Mean log2 abundance") +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
ggsave(file.path(PATHS$figures, "marker_missingness_dotplot.svg"), p_missing, width = 180, height = max(90, min(260, 3.2 * nrow(plot_missing))), units = "mm", device = svglite::svglite)

group_col <- intersect(c("Dataset", "dataset", "RegionLayer", "Region", "region", "Layer", "layer", "Group", "group", "ExpGroup", "Sex", "sex", "Batch", "batch", "plate"), names(sample_scores))
group_col <- group_col[!duplicated(tolower(group_col))]
if (length(group_col)) {
  sample_scores_plot <- sample_scores |>
    dplyr::group_by(marker_panel) |>
    dplyr::filter(any(is.finite(marker_score))) |>
    dplyr::ungroup()

  if (nrow(sample_scores_plot)) {
    p_box <- ggplot(sample_scores_plot, aes(x = marker_panel, y = marker_score, color = .data[[group_col[[1]]]])) +
      geom_boxplot(outlier.shape = NA, color = "grey35") +
      geom_point(position = position_jitter(width = 0.15), alpha = 0.75, size = 0.9) +
      coord_flip() +
      labs(x = NULL, y = "Mean marker abundance", color = group_col[[1]]) +
      theme_classic(base_size = 8) +
      theme(legend.position = "bottom")
    ggsave(file.path(PATHS$figures, "marker_score_sample_boxplot.svg"), p_box, width = 150, height = 105, units = "mm", device = svglite::svglite)
  }

  heat <- sample_scores_plot |>
    dplyr::mutate(marker_group = as.character(.data[[group_col[[1]]]])) |>
    dplyr::filter(!is.na(marker_group), nzchar(marker_group), is.finite(marker_score)) |>
    dplyr::group_by(marker_panel, marker_group) |>
    dplyr::summarise(median_marker_score = median(marker_score, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(median_marker_score = ifelse(is.nan(median_marker_score), NA_real_, median_marker_score))
  qc_write_csv(heat, file.path(PATHS$tables, "marker_detectability_heatmap_source.csv"))
  if (nrow(heat)) {
    p_heat <- ggplot(heat, aes(x = marker_group, y = marker_panel, fill = median_marker_score)) +
      geom_tile(color = "white", linewidth = 0.2) +
      labs(x = group_col[[1]], y = NULL, fill = "Median marker score") +
      theme_classic(base_size = 8) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
    ggsave(file.path(PATHS$figures, "marker_detectability_heatmap.svg"), p_heat, width = 150, height = 95, units = "mm", device = svglite::svglite)
  }

  sample_compartment_scores_plot <- sample_compartment_scores |>
    dplyr::group_by(marker_compartment) |>
    dplyr::filter(any(is.finite(marker_compartment_score))) |>
    dplyr::ungroup()

  if (nrow(sample_compartment_scores_plot)) {
    p_comp_box <- ggplot(sample_compartment_scores_plot, aes(x = marker_compartment, y = marker_compartment_score, color = .data[[group_col[[1]]]])) +
      geom_boxplot(outlier.shape = NA, color = "grey35") +
      geom_point(position = position_jitter(width = 0.15), alpha = 0.75, size = 0.9) +
      coord_flip() +
      labs(x = NULL, y = "Mean compartment marker abundance", color = group_col[[1]]) +
      theme_classic(base_size = 8) +
      theme(legend.position = "bottom")
    ggsave(file.path(PATHS$figures, "marker_compartment_score_sample_boxplot.svg"), p_comp_box, width = 150, height = 95, units = "mm", device = svglite::svglite)
  }

  comp_heat <- sample_compartment_scores_plot |>
    dplyr::mutate(marker_group = as.character(.data[[group_col[[1]]]])) |>
    dplyr::filter(!is.na(marker_group), nzchar(marker_group), is.finite(marker_compartment_score)) |>
    dplyr::group_by(marker_compartment, marker_group) |>
    dplyr::summarise(median_marker_compartment_score = median(marker_compartment_score, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(median_marker_compartment_score = ifelse(is.nan(median_marker_compartment_score), NA_real_, median_marker_compartment_score))
  qc_write_csv(comp_heat, file.path(PATHS$tables, "marker_compartment_detectability_heatmap_source.csv"))
  if (nrow(comp_heat)) {
    p_comp_heat <- ggplot(comp_heat, aes(x = marker_group, y = marker_compartment, fill = median_marker_compartment_score)) +
      geom_tile(color = "white", linewidth = 0.2) +
      labs(x = group_col[[1]], y = NULL, fill = "Median compartment score") +
      theme_classic(base_size = 8) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
    ggsave(file.path(PATHS$figures, "marker_compartment_detectability_heatmap.svg"), p_comp_heat, width = 150, height = 90, units = "mm", device = svglite::svglite)
  }

  fidelity_sample_scores_plot <- fidelity_sample_scores |>
    dplyr::group_by(fidelity_marker_class) |>
    dplyr::filter(any(is.finite(fidelity_marker_score))) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      fidelity_marker_class = factor(fidelity_marker_class, levels = fidelity_marker_class_order, ordered = TRUE)
    )

  if (nrow(fidelity_sample_scores_plot)) {
    p_fidelity_box <- ggplot(fidelity_sample_scores_plot, aes(x = fidelity_marker_class, y = fidelity_marker_score, color = .data[[group_col[[1]]]])) +
      geom_boxplot(outlier.shape = NA, color = "grey35") +
      geom_point(position = position_jitter(width = 0.15), alpha = 0.75, size = 0.9) +
      scale_x_discrete(drop = FALSE) +
      labs(x = NULL, y = "Mean fidelity marker abundance", color = group_col[[1]]) +
      theme_classic(base_size = 8) +
      theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom")
    ggsave(file.path(PATHS$figures, "marker_fidelity_score_sample_boxplot.svg"), p_fidelity_box, width = 120, height = 85, units = "mm", device = svglite::svglite)
  }

  fidelity_heat <- fidelity_sample_scores_plot |>
    dplyr::mutate(marker_group = as.character(.data[[group_col[[1]]]])) |>
    dplyr::filter(!is.na(marker_group), nzchar(marker_group), is.finite(fidelity_marker_score)) |>
    dplyr::group_by(fidelity_marker_class, marker_group) |>
    dplyr::summarise(median_fidelity_marker_score = median(fidelity_marker_score, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(median_fidelity_marker_score = ifelse(is.nan(median_fidelity_marker_score), NA_real_, median_fidelity_marker_score))
  qc_write_csv(fidelity_heat, file.path(PATHS$tables, "marker_fidelity_heatmap_source.csv"))
  if (nrow(fidelity_heat)) {
    p_fidelity_heat <- ggplot(fidelity_heat, aes(x = marker_group, y = fidelity_marker_class, fill = median_fidelity_marker_score)) +
      geom_tile(color = "white", linewidth = 0.2) +
      scale_y_discrete(drop = FALSE) +
      labs(x = group_col[[1]], y = NULL, fill = "Median fidelity score") +
      theme_classic(base_size = 8) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
    ggsave(file.path(PATHS$figures, "marker_fidelity_heatmap.svg"), p_fidelity_heat, width = 120, height = 75, units = "mm", device = svglite::svglite)
  }
}

# -----------------------------------------------------------------------------
# Optional WGCNA marker enrichment bridge
# -----------------------------------------------------------------------------

find_wgcna_modules <- function(dataset) {
  override <- Sys.getenv("PROTEOMICS_WGCNA_MODULE_ASSIGNMENT_FILE", unset = "")
  if (nzchar(override)) return(override)
  roots <- c(
    path_results("tables", "06_modules_WGCNA", "01_WGCNA", dataset, "modules"),
    path_results("tables", "06_modules_WGCNA", "01_WGCNA", dataset),
    path_results("tables", "06_modules_WGCNA", dataset),
    path_processed("06_modules_WGCNA", "01_WGCNA", dataset, "modules")
  )
  patterns <- c("module.*assignment.*\\.(csv|tsv)$", "module.*membership.*\\.(csv|tsv)$", "protein.*module.*\\.(csv|tsv)$", ".*modules.*\\.(csv|tsv)$")
  files <- unique(unlist(lapply(roots[dir.exists(roots)], function(root) {
    unique(unlist(lapply(patterns, function(pat) list.files(root, pattern = pat, recursive = TRUE, full.names = TRUE, ignore.case = TRUE))))
  })))
  if (length(files)) files[[1]] else NA_character_
}

read_wgcna_modules <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  df <- tryCatch(qc_read_table(path), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(NULL)
  protein_col <- qc_first_col(df, c("Protein", "protein", "ProteinID", "protein_id", "gene", "gene_symbol", "Gene", "feature", "Feature", "id"))
  module_col <- qc_first_col(df, c("module", "Module", "module_label", "module_color", "color", "WGCNA_module"))
  if (is.na(protein_col) || is.na(module_col)) return(NULL)
  df |>
    dplyr::transmute(module = as.character(.data[[module_col]]), Protein = as.character(.data[[protein_col]]), protein_key = gene_norm(.data[[protein_col]])) |>
    dplyr::filter(!is.na(module), nzchar(module), !is.na(protein_key), nzchar(protein_key)) |>
    dplyr::distinct(module, Protein, protein_key)
}

run_enrichment <- function(feature_keys, feature_col, modules, universe, dataset) {
  dplyr::bind_rows(lapply(names(modules), function(mod) {
    mod_keys <- unique(modules[[mod]][modules[[mod]] %in% universe])
    dplyr::bind_rows(lapply(seq_len(nrow(feature_keys)), function(i) {
      feature <- feature_keys[[feature_col]][[i]]
      keys <- unique(feature_keys$keys[[i]])
      overlap <- intersect(mod_keys, keys)
      a <- length(overlap)
      b <- length(setdiff(mod_keys, keys))
      c <- length(setdiff(keys, mod_keys))
      d <- length(setdiff(universe, union(mod_keys, keys)))
      ft <- tryCatch(stats::fisher.test(matrix(c(a, b, c, d), nrow = 2)), error = function(e) NULL)
      data.frame(
        dataset = dataset,
        module = mod,
        feature = feature,
        module_size = length(mod_keys),
        feature_size_detected = length(keys),
        overlap_n = a,
        overlap_genes = paste(sort(overlap), collapse = ";"),
        odds_ratio = if (is.null(ft)) NA_real_ else unname(ft$estimate),
        p_value = if (is.null(ft)) NA_real_ else ft$p.value,
        stringsAsFactors = FALSE
      )
    }))
  })) |>
    dplyr::group_by(feature) |>
    dplyr::mutate(padj_BH = stats::p.adjust(p_value, method = "BH")) |>
    dplyr::ungroup() |>
    dplyr::mutate(enrichment_status = dplyr::case_when(
      !is.na(padj_BH) & padj_BH < 0.05 ~ "FDR_significant",
      !is.na(p_value) & p_value < 0.05 ~ "nominal",
      TRUE ~ "not_supported"
    ))
}

module_file <- find_wgcna_modules(DATASET)
module_df <- read_wgcna_modules(module_file)
if (!is.null(module_df) && nrow(module_df)) {
  modules <- split(module_df$protein_key, module_df$module)

  marker_panel_keys <- marker_lookup |>
    dplyr::filter(marker_key %in% universe) |>
    dplyr::group_by(marker_panel) |>
    dplyr::summarise(keys = list(unique(marker_key)), .groups = "drop")

  enrich <- run_enrichment(marker_panel_keys, "marker_panel", modules, universe, DATASET) |>
    dplyr::rename(marker_panel = feature, panel_size_detected = feature_size_detected) |>
    dplyr::left_join(panel_metadata, by = "marker_panel")
  qc_write_csv(enrich, file.path(PATHS$tables, "wgcna_module_marker_panel_enrichment.csv"))

  heat_wgcna <- enrich |>
    dplyr::mutate(score = ifelse(is.finite(padj_BH) & padj_BH > 0, -log10(padj_BH), NA_real_),
                  score = ifelse(is.finite(odds_ratio) & odds_ratio < 1, -score, score))
  p_wgcna <- ggplot(heat_wgcna, aes(x = marker_panel, y = module, fill = score)) +
    geom_tile(color = "white", linewidth = 0.2) +
    labs(x = NULL, y = "WGCNA module", fill = "signed -log10 FDR") +
    theme_classic(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
  ggsave(file.path(PATHS$figures, "wgcna_module_marker_panel_enrichment_heatmap.svg"), p_wgcna,
         width = 170, height = max(90, min(240, 4 * length(unique(heat_wgcna$module)))), units = "mm", device = svglite::svglite)

  marker_compartment_keys <- marker_lookup |>
    dplyr::filter(marker_key %in% universe, !is.na(marker_compartment), nzchar(marker_compartment)) |>
    dplyr::group_by(marker_compartment) |>
    dplyr::summarise(keys = list(unique(marker_key)), .groups = "drop")

  enrich_comp <- run_enrichment(marker_compartment_keys, "marker_compartment", modules, universe, DATASET) |>
    dplyr::rename(marker_compartment = feature, compartment_size_detected = feature_size_detected)
  qc_write_csv(enrich_comp, file.path(PATHS$tables, "wgcna_module_marker_compartment_enrichment.csv"))

  heat_wgcna_comp <- enrich_comp |>
    dplyr::mutate(score = ifelse(is.finite(padj_BH) & padj_BH > 0, -log10(padj_BH), NA_real_),
                  score = ifelse(is.finite(odds_ratio) & odds_ratio < 1, -score, score))
  p_wgcna_comp <- ggplot(heat_wgcna_comp, aes(x = marker_compartment, y = module, fill = score)) +
    geom_tile(color = "white", linewidth = 0.2) +
    labs(x = NULL, y = "WGCNA module", fill = "signed -log10 FDR") +
    theme_classic(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
  ggsave(file.path(PATHS$figures, "wgcna_module_marker_compartment_enrichment_heatmap.svg"), p_wgcna_comp,
         width = 150, height = max(90, min(240, 4 * length(unique(heat_wgcna_comp$module)))), units = "mm", device = svglite::svglite)

  fidelity_marker_keys <- marker_lookup |>
    dplyr::filter(marker_key %in% universe, !is.na(fidelity_marker_class), nzchar(fidelity_marker_class)) |>
    dplyr::mutate(fidelity_marker_class = factor(fidelity_marker_class, levels = fidelity_marker_class_order, ordered = TRUE)) |>
    dplyr::group_by(fidelity_marker_class) |>
    dplyr::summarise(keys = list(unique(marker_key)), .groups = "drop") |>
    dplyr::arrange(fidelity_marker_class)

  if (nrow(fidelity_marker_keys)) {
    enrich_fidelity <- run_enrichment(fidelity_marker_keys, "fidelity_marker_class", modules, universe, DATASET) |>
      dplyr::rename(fidelity_marker_class = feature, fidelity_class_size_detected = feature_size_detected)
    qc_write_csv(enrich_fidelity, file.path(PATHS$tables, "wgcna_module_marker_fidelity_enrichment.csv"))

    heat_wgcna_fidelity <- enrich_fidelity |>
      dplyr::mutate(score = ifelse(is.finite(padj_BH) & padj_BH > 0, -log10(padj_BH), NA_real_),
                    score = ifelse(is.finite(odds_ratio) & odds_ratio < 1, -score, score))
    p_wgcna_fidelity <- ggplot(heat_wgcna_fidelity, aes(x = fidelity_marker_class, y = module, fill = score)) +
      geom_tile(color = "white", linewidth = 0.2) +
      scale_x_discrete(drop = FALSE) +
      labs(x = NULL, y = "WGCNA module", fill = "signed -log10 FDR") +
      theme_classic(base_size = 8) +
      theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom")
    ggsave(file.path(PATHS$figures, "wgcna_module_marker_fidelity_enrichment_heatmap.svg"), p_wgcna_fidelity,
           width = 120, height = max(90, min(240, 4 * length(unique(heat_wgcna_fidelity$module)))), units = "mm", device = svglite::svglite)
  }
} else {
  message("No WGCNA module assignment table found for dataset ", DATASET, "; skipping WGCNA marker-panel, marker-compartment, and marker-fidelity enrichment bridge.")
}

# -----------------------------------------------------------------------------
# Report and manifest
# -----------------------------------------------------------------------------


# Additional human-readable audit report for protein provenance in the narrow
# fidelity layer. The CSV files are the authoritative machine-readable outputs;
# this markdown is for quick inspection.
fidelity_audit_lines <- c(
  "# Soma/neuropil/microglia fidelity protein audit",
  "",
  paste0("Dataset: ", DATASET),
  paste0("Expression matrix: ", matrix_file),
  paste0("Marker panel source: ", if (file.exists(marker_file)) marker_file else "internal fallback panels"),
  "",
  "This audit lists the marker panels and proteins that actually contributed to the narrow fidelity scores used downstream by 04d.",
  "A protein contributes when the requested marker matches at least one protein ID in the active matrix.",
  "The actual fidelity score uses deduplicated marker keys within each fidelity class; repeated markers across panels are not overweighted.",
  "Broad Allen neuronal reference panels are intentionally excluded from the soma/neuropil fidelity layer.",
  ""
)

if (nrow(marker_fidelity_protein_summary)) {
  for (cls in fidelity_marker_class_order) {
    cls_tbl <- marker_fidelity_protein_summary |>
      dplyr::filter(.data$fidelity_marker_class == cls)
    if (!nrow(cls_tbl)) next

    class_used <- sum(cls_tbl$n_used_markers, na.rm = TRUE)
    class_requested <- sum(cls_tbl$n_requested_markers, na.rm = TRUE)
    fidelity_audit_lines <- c(
      fidelity_audit_lines,
      paste0("## ", cls),
      paste0("Detected/used marker keys: ", class_used, "/", class_requested),
      ""
    )

    for (i in seq_len(nrow(cls_tbl))) {
      row <- cls_tbl[i, , drop = FALSE]
      fidelity_audit_lines <- c(
        fidelity_audit_lines,
        paste0("### ", row$marker_panel),
        paste0("Detected/used: ", row$n_used_markers, "/", row$n_requested_markers),
        paste0("Used markers: ", ifelse(nzchar(row$used_requested_markers), row$used_requested_markers, "none")),
        paste0("Matched protein IDs: ", ifelse(nzchar(row$used_matched_protein_ids), row$used_matched_protein_ids, "none")),
        paste0("Missing requested markers: ", ifelse(nzchar(row$missing_requested_markers), row$missing_requested_markers, "none")),
        ""
      )
    }
  }
} else {
  fidelity_audit_lines <- c(fidelity_audit_lines, "No fidelity marker proteins were available for this dataset.", "")
}



if (nrow(fidelity_marker_proteins_used_deduplicated)) {
  fidelity_audit_lines <- c(
    fidelity_audit_lines,
    "# Deduplicated scoring-level marker list",
    "",
    "This section lists one row per unique marker key per fidelity class, matching the actual scoring logic.",
    "If a marker occurs in more than one panel, it is still counted once in the score.",
    ""
  )

  for (cls in fidelity_marker_class_order) {
    cls_tbl <- fidelity_marker_proteins_used_deduplicated |>
      dplyr::filter(.data$fidelity_marker_class == cls)
    if (!nrow(cls_tbl)) next

    used_tbl <- cls_tbl |>
      dplyr::filter(.data$used_in_fidelity_score)
    missing_tbl <- cls_tbl |>
      dplyr::filter(!.data$used_in_fidelity_score)

    fidelity_audit_lines <- c(
      fidelity_audit_lines,
      paste0("## Deduplicated ", cls),
      paste0("Used unique marker keys: ", nrow(used_tbl), "/", nrow(cls_tbl)),
      paste0("Used marker keys: ", ifelse(nrow(used_tbl), paste(sort(unique(used_tbl$marker_key)), collapse = ";"), "none")),
      paste0("Matched protein IDs: ", ifelse(nrow(used_tbl), collapse_unique(used_tbl$matched_protein_id), "none")),
      paste0("Missing marker keys: ", ifelse(nrow(missing_tbl), paste(sort(unique(missing_tbl$marker_key)), collapse = ";"), "none")),
      ""
    )
  }
}

if (nrow(fidelity_duplicate_markers)) {
  fidelity_audit_lines <- c(
    fidelity_audit_lines,
    "# Repeated marker membership across panels",
    "",
    "These markers occur in more than one marker panel within the same fidelity class. They are reported here for transparency but are counted only once in the fidelity score.",
    ""
  )

  for (i in seq_len(nrow(fidelity_duplicate_markers))) {
    row <- fidelity_duplicate_markers[i, , drop = FALSE]
    fidelity_audit_lines <- c(
      fidelity_audit_lines,
      paste0("- ", row$fidelity_marker_class, ": ", row$marker_key,
             " in ", row$n_marker_panels_containing_marker, " panels [",
             row$contributing_marker_panels, "]" )
    )
  }
  fidelity_audit_lines <- c(fidelity_audit_lines, "")
}

writeLines(fidelity_audit_lines, file.path(PATHS$reports, "marker_fidelity_protein_audit.md"))

notes <- c(
  "# Marker detectability and WGCNA marker-panel bridge", "",
  "This QC report evaluates whether canonical and reference marker proteins are detectable and consistently quantified in the active dataset.", "",
  "Interpretation constraints:",
  "- Marker-panel scoring is QC and interpretive support, not formal cell-type purity estimation or deconvolution.",
  "- marker_panel refers to an individual marker_set/signature from the registry, e.g. canonical_microglia_homeostatic or reference_microglia_pvm.",
  "- marker_compartment is a broad collapsed biological class used only for QC interpretation and robust downstream plotting.",
  "- fidelity_marker_class is a narrow soma/neuropil/microglia marker class for comparing neuron_soma, neuron_neuropil, and microglia samples.",
  "- The fidelity layer intentionally uses only curated soma/neuropil panels and conservative microglia/PVM panels; broad Allen neuronal cell-class panels are not treated as soma or neuropil markers.",
  "- Compartment-level and fidelity-level scores use deduplicated detected marker proteins, not simple averages of panel scores.",
  "- marker_fidelity_proteins_used.csv is the panel-level protein audit trail for the 04d soma/neuropil/microglia fidelity plot.",
  "- marker_fidelity_proteins_used_deduplicated.csv is the score-level audit trail: one marker key per fidelity class, matching the deduplicated scoring logic.",
  "- marker_fidelity_duplicate_marker_membership.csv lists markers that occur in multiple panels but are counted only once in the fidelity score.",
  "- marker_fidelity_protein_summary.csv summarizes which requested markers matched proteins and which were missing for each fidelity marker panel.",
  "- marker_source_detectability_summary.csv summarizes requested and detected markers by upstream source_name/source_type/evidence_level, so SynGO, GO/MGI, curated seed, and generic reference sources can be checked directly.",
  "- Low detectability or high missingness of canonical microglia markers can explain why WGCNA modules are driven by better-detected ECM, synaptic, vascular, mitochondrial, ribosomal, or microenvironment-associated proteins.",
  "- Marker abundance and WGCNA module assignment answer different questions: sample identity/abundance versus covariance/module topology.",
  "- Agreement in directionality between WGCNA eigengenes and previous GSEA/leading-edge results is robustness evidence even when module labels differ from marker-panel labels.", "",
  paste0("Dataset: ", DATASET),
  paste0("Expression matrix: ", matrix_file),
  paste0("Marker panel source: ", if (file.exists(marker_file)) marker_file else "internal fallback panels"),
  paste0("WGCNA module assignment source: ", if (!is.na(module_file) && file.exists(module_file)) module_file else "not found"),
  paste0("Marker panels scored: ", length(unique(marker_lookup$marker_panel))),
  paste0("Marker compartments scored: ", length(unique(marker_lookup$marker_compartment))),
  paste0("Fidelity marker classes scored: ", length(unique(marker_lookup$fidelity_marker_class[!is.na(marker_lookup$fidelity_marker_class)])))
)
writeLines(notes, file.path(PATHS$reports, "marker_detectability_interpretation_notes.md"))

write_run_manifest(
  file.path(PATHS$logs, "run_manifest.yml"),
  inputs = list(matrix = matrix_file, metadata = metadata_file, marker_panels = marker_file, wgcna_modules = module_file),
  outputs = list(figures = PATHS$figures, tables = PATHS$tables, reports = PATHS$reports),
  parameters = list(
    dataset = DATASET,
    marker_sets = sort(unique(marker_lookup$marker_panel)),
    marker_compartments = sort(unique(marker_lookup$marker_compartment)),
    fidelity_marker_classes = sort(unique(marker_lookup$fidelity_marker_class[!is.na(marker_lookup$fidelity_marker_class)]))
  ),
  notes = "Marker-panel, marker-compartment, and soma/neuropil/microglia fidelity detectability/missingness QC plus optional WGCNA enrichment bridge; not formal purity/deconvolution."
)

message("Marker detectability QC complete for dataset: ", DATASET)

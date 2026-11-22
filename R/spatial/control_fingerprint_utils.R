# Control spatial fingerprint: the CON baseline profile of every protein, the
# external signature scores and the named proteins of Extended Data Fig. 6c/e.
#
# Used by analysis/differential_abundance/build_control_fingerprint.R. Ported
# from two retired producers: figures/spatial_v6_fingerprint_selection.R
# (Exp9_manuscript, removed in 8e22dd9) and the CON baseline builder kept as
# audits/spatial/spatial_v6_baseline_profile.R. Three CON animals only;
# descriptive: no model, no test, no p-value.
#
# Rows are chosen without any stress information: external hippocampal
# signatures (Kaulich 2025), external cell-class panels (config/marker_panels/)
# or the top proteins of the prespecified CON-only anatomical contrasts. Values
# are within-protein z-scores of the CON mean.
#
# Gene symbols: every join uses the official gene symbol of the canonical
# feature table, the symbol the anatomical contrasts and the Kaulich mapping
# carry. The retired producer took the baseline symbol from the WGCNA shortlist,
# which then spelled histone H1.3 "H1-3" while the contrasts spelled it "H1f3":
# H1F3 was selected but dropped from the proteins table.
# cfp_check_selection_in_proteins() now stops on any such gap.

CFP_MIN_SIGNATURE_MEMBERS <- 5L
CFP_TOP_PROTEINS_PER_CONTRAST <- 2L

# External cell-class panels shown as signature rows, with their row labels.
CFP_REFERENCE_SETS <- c(
  reference_hippocampal_excitatory_neuron = "Excitatory neuron",
  reference_inhibitory_interneuron = "Inhibitory interneuron",
  reference_astrocyte = "Astrocyte",
  reference_oligodendrocyte = "Oligodendrocyte",
  reference_microglia_pvm = "Microglia / PVM",
  reference_vascular = "Vascular",
  canonical_neuronal_synaptic_neuropil = "Synaptic neuropil"
)

CFP_EXCLUDED_SOURCE_NOTE <- paste0(
  "protein_baseline_spatial_profile.csv (build_protein_spatial_atlas) was NOT ",
  "used for selection: its scope is is_sus_res_fdr_supported | ",
  "is_wgcna_candidate, so its rows were chosen by SUS - RES support or WGCNA ",
  "candidacy. Its VALUES are CON-only and are reproduced exactly here on the ",
  "full protein universe, but its ROW MEMBERSHIP is not phenotype-blind.")

# Upper-case symbol used as the join key; blank and missing become NA.
cfp_gene_key <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[!is.na(x) & !nzchar(x)] <- NA_character_
  x
}

cfp_require_columns <- function(df, cols, what) {
  missing <- setdiff(cols, names(df))
  if (length(missing)) {
    stop(what, " is missing column(s): ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  invisible(df)
}

# ---- spatial units (config/manuscript_spatial_order.yml) -------------------

# The 18 canonical units and their spellings. A unit is resolved within its
# dataset, because soma and microglia both spell their units CA1/CA2/CA3/DG.
cfp_spatial_order <- function(path = repo_path("config", "manuscript_spatial_order.yml")) {
  if (!file.exists(path)) stop("Spatial order contract not found: ", path, call. = FALSE)
  y <- yaml::read_yaml(path)
  chr <- function(x) if (is.null(x)) NA_character_ else as.character(x)
  units <- do.call(rbind, lapply(y$units, function(z) data.frame(
    dataset = chr(z$dataset),
    order = as.integer(z$order),
    region = chr(z$region),
    layer = chr(z$layer),
    layer_is_resolution = !identical(z$layer_is_resolution, FALSE) && !is.null(z$layer),
    analysis_key = chr(z$analysis_key),
    display = chr(z$display),
    stringsAsFactors = FALSE)))
  units <- units[order(units$order), , drop = FALSE]
  rownames(units) <- NULL
  aliases <- do.call(rbind, lapply(y$units, function(z) {
    keys <- unique(as.character(c(z$analysis_key, z$atlas_key, unlist(z$aliases))))
    keys <- keys[!is.na(keys) & nzchar(keys)]
    data.frame(dataset = chr(z$dataset), alias = tolower(keys),
               analysis_key = chr(z$analysis_key), stringsAsFactors = FALSE)
  }))
  aliases <- unique(aliases)
  clash <- duplicated(paste(aliases$dataset, aliases$alias))
  if (any(clash)) {
    stop("Spatial order contract uses an alias twice within a dataset: ",
         paste(unique(paste(aliases$dataset[clash], aliases$alias[clash])), collapse = ", "),
         call. = FALSE)
  }
  compartments <- do.call(rbind, lapply(y$compartments, function(z) data.frame(
    id = chr(z$id), order = as.integer(z$order), short = chr(z$short),
    stringsAsFactors = FALSE)))
  layers <- do.call(rbind, lapply(y$layers, function(z) data.frame(
    id = chr(z$id), display = chr(z$display), stringsAsFactors = FALSE)))
  list(units = units, aliases = aliases,
       compartments = compartments[order(compartments$order), , drop = FALSE],
       layers = layers)
}

# (dataset, unit) to the canonical unit key. An unknown pair is an error,
# never a silently dropped row.
cfp_resolve_unit <- function(unit, dataset, spatial_order) {
  u <- tolower(trimws(as.character(unit)))
  d <- as.character(dataset)
  if (length(d) == 1L) d <- rep(d, length(u))
  if (length(d) != length(u)) {
    stop("cfp_resolve_unit: dataset must have length 1 or length(unit).", call. = FALSE)
  }
  a <- spatial_order$aliases
  out <- a$analysis_key[match(paste(d, u), paste(a$dataset, a$alias))]
  bad <- is.na(out) & !is.na(u) & nzchar(u)
  if (any(bad)) {
    stop("Unrecognised (dataset, spatial unit) pair(s): ",
         paste(unique(paste0(d[bad], "/", unit[bad])), collapse = ", "), call. = FALSE)
  }
  out
}

# The sg_* columns the figure readers use, with the names and meanings of the
# manuscript's spatial grammar.
cfp_annotate_units <- function(df, spatial_order, unit_col = "spatial_unit",
                               dataset_col = "dataset") {
  cfp_require_columns(df, c(unit_col, dataset_col), "cfp_annotate_units input")
  key <- cfp_resolve_unit(df[[unit_col]], df[[dataset_col]], spatial_order)
  u <- spatial_order$units
  i <- match(paste(df[[dataset_col]], key), paste(u$dataset, u$analysis_key))
  cs <- spatial_order$compartments
  df$sg_unit <- key
  df$sg_compartment <- as.character(df[[dataset_col]])
  df$sg_compartment_label <- cs$short[match(df[[dataset_col]], cs$id)]
  df$sg_region <- u$region[i]
  df$sg_layer <- u$layer[i]
  df$sg_layer_display <- spatial_order$layers$display[match(u$layer[i], spatial_order$layers$id)]
  df$sg_layer_is_resolution <- u$layer_is_resolution[i]
  df$sg_unit_label <- u$display[i]
  df$sg_unit_f <- key
  df
}

# ---- CON baseline ----------------------------------------------------------

# CON baseline of one dataset: per protein and spatial unit, the mean over CON
# animals of the level-2 log2 abundance (animal x unit, hemispheres combined),
# and its within-protein z-score across the dataset's units. `mat` and `meta`
# are lv$level2$mat and lv$level2$meta from sps_build_spatial_levels().
cfp_con_profile <- function(mat, meta, dataset, spatial_order) {
  cfp_require_columns(meta, c("AnimalID", "StressGroup", "SpatialUnit"), "Level-2 metadata")
  if (ncol(mat) != nrow(meta)) stop("Level-2 matrix and metadata differ in length.", call. = FALSE)
  con <- meta$StressGroup %in% "CON"
  if (!any(con)) stop("No CON samples for ", dataset, ".", call. = FALSE)
  mat <- mat[, con, drop = FALSE]
  meta <- meta[con, , drop = FALSE]
  units <- sort(unique(as.character(meta$SpatialUnit)))
  prof <- vapply(units, function(u)
    rowMeans(mat[, meta$SpatialUnit == u, drop = FALSE], na.rm = TRUE),
    numeric(nrow(mat)))
  if (is.null(dim(prof))) prof <- matrix(prof, nrow = nrow(mat))
  dimnames(prof) <- list(rownames(mat), units)

  # within-protein standardisation, so the panel shows spatial pattern rather
  # than which proteins happen to be abundant
  m <- rowMeans(prof, na.rm = TRUE)
  s <- apply(prof, 1L, stats::sd, na.rm = TRUE)
  z <- (prof - m) / s
  z[!(is.finite(s) & s > 0), ] <- NA_real_

  raw_unit <- rep(colnames(prof), each = nrow(prof))
  data.frame(
    dataset = dataset,
    ProteinGroupID = rep(rownames(prof), times = ncol(prof)),
    raw_unit = raw_unit,
    con_mean_log2 = as.vector(prof),
    spatial_unit = cfp_resolve_unit(raw_unit, dataset, spatial_order),
    con_z = as.vector(z),
    protein_mean_log2 = rep(unname(m), times = ncol(prof)),
    protein_sd_log2 = rep(unname(s), times = ncol(prof)),
    n_con_animals = length(unique(meta$AnimalID)),
    stringsAsFactors = FALSE)
}

# The baseline must equal the atlas CON baseline on every cell they share: same
# code path, the atlas only keeps fewer proteins. `atlas` is
# protein_baseline_spatial_profile.csv (one column per unit, atlas spelling).
cfp_reconcile_atlas <- function(base, atlas, spatial_order, tolerance = 1e-8) {
  cfp_require_columns(atlas, c("dataset", "ProteinGroupID"), "Atlas baseline")
  unit_cols <- setdiff(names(atlas), c("dataset", "ProteinGroupID", "n_con_animals"))
  long <- do.call(rbind, lapply(unit_cols, function(u) data.frame(
    dataset = atlas$dataset, ProteinGroupID = atlas$ProteinGroupID,
    raw_unit = u, atlas_value = atlas[[u]], stringsAsFactors = FALSE)))
  long <- long[!is.na(long$atlas_value), , drop = FALSE]
  long$spatial_unit <- cfp_resolve_unit(long$raw_unit, long$dataset, spatial_order)
  j <- merge(long, base[, c("dataset", "ProteinGroupID", "spatial_unit", "con_mean_log2")],
             by = c("dataset", "ProteinGroupID", "spatial_unit"))
  d <- abs(j$atlas_value - j$con_mean_log2)
  max_d <- if (nrow(j)) max(d, na.rm = TRUE) else NA_real_
  if (nrow(j) && is.finite(max_d) && max_d > tolerance) {
    stop("The CON baseline does not reproduce the atlas CON baseline on the cells ",
         "they share (max abs difference ", format(max_d, scientific = TRUE), ").",
         call. = FALSE)
  }
  list(shared_cells = nrow(j), max_abs_difference = max_d)
}

# ---- form B: external signature scores -------------------------------------

# Signature membership: the Kaulich 2025 signatures (mapped genes only) and the
# external cell-class panels of config/marker_panels/.
cfp_signature_sets <- function(kaulich_mapping, reference_panels,
                               reference_sets = CFP_REFERENCE_SETS) {
  cfp_require_columns(kaulich_mapping, c("external_signature", "canonical_official_gene_symbol"),
                      "Kaulich signature mapping")
  cfp_require_columns(reference_panels, c("marker_set", "gene_symbol", "source_name", "use_for"),
                      "Reference marker panels")
  kau <- kaulich_mapping
  kau$gene <- cfp_gene_key(kau$canonical_official_gene_symbol)
  kau <- kau[!is.na(kau$gene), , drop = FALSE]
  k <- unique(kau[, c("external_signature", "gene")])
  k$signature <- paste0("Kaulich ", k$external_signature)
  k$signature_family <- ifelse(
    k$external_signature %in% c("CA1", "CA2/3", "DG"),
    "external hippocampal subregion signature",
    "external hippocampal strata signature")
  k$source <- "Kaulich 2025 external published hippocampal proteome"
  k$selection_rule <- paste0(
    "all genes of the external published signature that map to a canonical ",
    "official mouse symbol (kaulich_signature_mapping.csv, mapping_status = mapped)")
  k$declared_allowed_use <- "external spatial validation reference"

  ref <- reference_panels
  ref$gene <- cfp_gene_key(ref$gene_symbol)
  ref <- ref[ref$marker_set %in% names(reference_sets) & !is.na(ref$gene), , drop = FALSE]
  r <- unique(ref[, c("marker_set", "gene", "source_name", "use_for")])
  r$signature <- unname(reference_sets[r$marker_set])
  r$signature_family <- "external curated cell-class reference panel"
  r$source <- paste0("config/marker_panels/wgcna_reference_marker_sets.csv (",
                     r$source_name, ")")
  r$selection_rule <- paste0(
    "every gene of the prespecified marker_set, fixed in a config file before ",
    "any stress contrast was computed")
  r$declared_allowed_use <- r$use_for

  cols <- c("signature", "signature_family", "gene", "source", "selection_rule",
            "declared_allowed_use")
  out <- unique(rbind(k[, cols, drop = FALSE], r[, cols, drop = FALSE]))
  rownames(out) <- NULL
  out
}

# Score = mean within-protein z of a signature's measured member proteins per
# (dataset, spatial unit); NA when fewer than `min_members` are measured.
cfp_score_signatures <- function(base, signatures, min_members = CFP_MIN_SIGNATURE_MEMBERS) {
  b <- base[!is.na(base$gene), c("dataset", "spatial_unit", "gene", "con_z", "ProteinGroupID"),
            drop = FALSE]
  members <- unique(signatures[, c("signature", "signature_family", "gene"), drop = FALSE])
  j <- merge(b, members, by = "gene")
  empty <- data.frame(signature = character(), signature_family = character(),
                      dataset = character(), spatial_unit = character(),
                      n_measured_members = integer(), score = numeric(),
                      stringsAsFactors = FALSE)
  if (!nrow(j)) return(empty)
  key <- paste(j$signature, j$dataset, j$spatial_unit, sep = "\r")
  out <- do.call(rbind, lapply(split(seq_len(nrow(j)), key), function(ix) {
    z <- j[ix, , drop = FALSE]
    data.frame(signature = z$signature[1], signature_family = z$signature_family[1],
               dataset = z$dataset[1], spatial_unit = z$spatial_unit[1],
               n_measured_members = length(unique(z$ProteinGroupID)),
               score = mean(z$con_z, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out$score[out$n_measured_members < min_members] <- NA_real_
  out[order(out$signature, out$dataset, out$spatial_unit, method = "radix"), , drop = FALSE]
}

# ---- form A: named proteins from the CON-only anatomical contrasts ---------

# The top `top_n` proteins of the enriched side (logFC > 0) of each contrast,
# by BH-adjusted p, then |logFC|, then symbol; one row per symbol and contrast.
cfp_select_contrast_proteins <- function(contrasts, top_n = CFP_TOP_PROTEINS_PER_CONTRAST) {
  cfp_require_columns(contrasts, c("dataset", "contrast", "ProteinGroupID",
                                   "official_gene_symbol", "logFC", "adj.P.Val"),
                      "Anatomical protein contrasts")
  a <- contrasts
  a$gene <- cfp_gene_key(a$official_gene_symbol)
  a <- a[!is.na(a$gene) & is.finite(a$adj.P.Val) & is.finite(a$logFC) & a$logFC > 0, ,
         drop = FALSE]
  if (!nrow(a)) stop("No anatomical contrast has an enriched protein with a gene symbol.", call. = FALSE)
  pick <- do.call(rbind, lapply(split(seq_len(nrow(a)), paste(a$dataset, a$contrast)), function(ix) {
    z <- a[ix, , drop = FALSE]
    z <- z[order(z$adj.P.Val, -abs(z$logFC), z$gene, method = "radix"), , drop = FALSE]
    z <- z[!duplicated(z$gene), , drop = FALSE]
    z <- utils::head(z, top_n)
    data.frame(gene = z$gene, contrast = z$contrast, dataset = z$dataset,
               ProteinGroupID = z$ProteinGroupID, logFC = z$logFC,
               adj_P_Val = z$adj.P.Val, stringsAsFactors = FALSE)
  }))
  rownames(pick) <- NULL
  pick
}

cfp_contrast_selection_rule <- function(top_n = CFP_TOP_PROTEINS_PER_CONTRAST) {
  paste0(
    "the top ", top_n, " genes by BH-adjusted p among the ENRICHED side ",
    "(logFC > 0) of each prespecified CON-only anatomical contrast; the ",
    "contrast set was fixed in advance and is fitted on CON animals only")
}

# Every selected symbol across all spatial units of all compartments where it
# is measured. selected_dataset and selected_ProteinGroupID name the protein
# group the contrast picked.
cfp_fingerprint_proteins <- function(base, pick, top_n = CFP_TOP_PROTEINS_PER_CONTRAST) {
  b <- base[!is.na(base$gene), c("gene", "dataset", "spatial_unit", "con_z",
                                 "con_mean_log2", "ProteinGroupID"), drop = FALSE]
  sel <- unique(pick[, c("gene", "contrast", "dataset", "ProteinGroupID"), drop = FALSE])
  names(sel) <- c("gene", "contrast", "selected_dataset", "selected_ProteinGroupID")
  prot <- merge(b, sel, by = "gene")
  prot$signature_family <- rep("prespecified CON-only anatomical contrast", nrow(prot))
  prot$source <- rep("validate_control_spatial_identity/anatomical_protein_contrasts.csv", nrow(prot))
  prot$selection_rule <- rep(cfp_contrast_selection_rule(top_n), nrow(prot))
  prot$declared_allowed_use <- rep("control spatial identity validation", nrow(prot))
  cols <- c("gene", "dataset", "spatial_unit", "con_z", "con_mean_log2", "ProteinGroupID",
            "contrast", "signature_family", "source", "selection_rule",
            "declared_allowed_use", "selected_dataset", "selected_ProteinGroupID")
  prot <- prot[order(prot$gene, prot$contrast, prot$dataset, prot$spatial_unit,
                     method = "radix"), cols, drop = FALSE]
  rownames(prot) <- NULL
  prot
}

# Stops unless every selected protein is in the proteins table: its symbol and
# contrast, and the selected protein group itself in its own compartment. Also
# stops on two rows for one symbol, contrast and spatial unit, which would draw
# two tiles on one cell.
cfp_check_selection_in_proteins <- function(pick, proteins) {
  cfp_require_columns(pick, c("gene", "contrast", "dataset", "ProteinGroupID"), "Selection")
  cfp_require_columns(proteins, c("gene", "contrast", "dataset", "ProteinGroupID", "spatial_unit"),
                      "Proteins table")
  problems <- character()
  sel <- unique(pick[, c("gene", "contrast"), drop = FALSE])
  have <- paste(proteins$gene, proteins$contrast, sep = "\r")
  gone <- sel[!paste(sel$gene, sel$contrast, sep = "\r") %in% have, , drop = FALSE]
  if (nrow(gone)) {
    problems <- c(problems, paste0(
      "selected but absent from the proteins table: ",
      paste0(gone$gene, " (", gone$contrast, ")", collapse = ", ")))
  }
  own <- unique(pick[, c("gene", "dataset", "ProteinGroupID"), drop = FALSE])
  have_own <- paste(proteins$gene, proteins$dataset, proteins$ProteinGroupID, sep = "\r")
  lost <- own[!paste(own$gene, own$dataset, own$ProteinGroupID, sep = "\r") %in% have_own, ,
              drop = FALSE]
  if (nrow(lost)) {
    problems <- c(problems, paste0(
      "selected protein group not shown under its symbol: ",
      paste0(lost$gene, " = ", lost$ProteinGroupID, " in ", lost$dataset, collapse = ", ")))
  }
  cell <- paste(proteins$gene, proteins$contrast, proteins$dataset, proteins$spatial_unit, sep = "\r")
  dup <- unique(proteins[duplicated(cell), c("gene", "dataset", "spatial_unit"), drop = FALSE])
  if (nrow(dup)) {
    problems <- c(problems, paste0(
      "more than one protein group per cell: ",
      paste0(dup$gene, " in ", dup$dataset, "/", dup$spatial_unit, collapse = ", ")))
  }
  if (length(problems)) {
    stop("Fingerprint selection and proteins table disagree; ",
         paste(problems, collapse = "; "), call. = FALSE)
  }
  invisible(TRUE)
}

# ---- selection provenance --------------------------------------------------

cfp_selection_table <- function(pick, signatures, top_n = CFP_TOP_PROTEINS_PER_CONTRAST) {
  sel_a <- unique(data.frame(
    display_form = "A_individual_proteins",
    row_id = pick$gene,
    source = "validate_control_spatial_identity/anatomical_protein_contrasts.csv",
    selection_rule = cfp_contrast_selection_rule(top_n),
    phenotype_blind = "YES",
    phenotype_blind_evidence = paste0(
      "contrast fitted on StressGroup == CON only ",
      "(R/spatial/control_spatial_identity_utils.R); no stress group, contrast, ",
      "effect or FDR enters selection"),
    dataset = pick$dataset,
    biological_interpretation = pick$contrast,
    selected_from_external_or_prespecified_source = "prespecified CON-only contrast",
    ProteinGroupID = pick$ProteinGroupID,
    stringsAsFactors = FALSE))
  sel_b <- unique(data.frame(
    display_form = "B_signature_scores",
    row_id = signatures$signature,
    source = signatures$source,
    selection_rule = signatures$selection_rule,
    phenotype_blind = "YES",
    phenotype_blind_evidence = ifelse(
      grepl("^Kaulich", signatures$signature),
      "membership comes from an external published dataset with no contact with these animals",
      "membership fixed in a repository config file, external curated panel"),
    dataset = NA_character_,
    biological_interpretation = signatures$signature_family,
    selected_from_external_or_prespecified_source = "external",
    ProteinGroupID = NA_character_,
    stringsAsFactors = FALSE))
  sel_b <- sel_b[!duplicated(sel_b$row_id), , drop = FALSE]
  sel <- rbind(sel_a, sel_b)
  sel$excluded_source_note <- CFP_EXCLUDED_SOURCE_NOTE
  rownames(sel) <- NULL
  sel[, c("display_form", "row_id", "source", "selection_rule", "phenotype_blind",
          "phenotype_blind_evidence", "dataset", "biological_interpretation",
          "selected_from_external_or_prespecified_source", "excluded_source_note",
          "ProteinGroupID")]
}

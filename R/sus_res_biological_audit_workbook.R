# Human-facing XLSX presentation of the canonical SUS-RES ontology-theme audits.
# CSV files remain the machine-readable computational contracts.

sus_res_audit_required_columns <- function() {
  list(
    overview = c(
      "biological_domain", "manuscript_theme_id", "manuscript_theme", "role",
      "supported_GO_occurrences", "unique_supported_GO_IDs", "n_units_evaluable",
      "n_units_with_FDR_support", "descriptive_higher_SUS_units",
      "descriptive_higher_RES_units", "FDR_supported_SUS_units",
      "FDR_supported_RES_units", "spatial_concordance_classification",
      "strongest_supported_representative_GO_term", "strongest_original_FDR", "short_audit_note"
    ),
    detail = c(
      "dataset", "compartment", "region", "layer", "spatial_unit",
      "manuscript_theme_id", "manuscript_theme", "theme_role",
      "n_theme_terms_tested", "median_NES_all_theme_terms", "q25_NES_all_theme_terms",
      "q75_NES_all_theme_terms", "fraction_positive_NES", "fraction_negative_NES",
      "descriptive_direction", "low_theme_term_coverage", "min_original_GO_FDR",
      "FDR_support_present", "n_FDR_supported_GO_terms", "n_positive_supported_terms",
      "n_negative_supported_terms", "representative_GO_ID", "representative_GO_term",
      "representative_NES", "representative_FDR", "direction_consistency",
      "n_semantic_clusters", "leading_edge_genes", "leading_edge_proteins",
      "n_FDR_supported_SUS_RES_DAPs"
    ),
    go_audit = c(
      "dataset", "compartment", "region", "layer", "spatial_unit", "GO_ID",
      "GO_description", "NES", "p.adjust", "assignment_status", "manuscript_themes",
      "theme_roles", "semantic_cluster_id", "semantic_cluster_representative_term",
      "leading_edge_genes", "leading_edge_proteins"
    ),
    dap = c(
      "dataset", "dataset_label", "spatial_unit", "n_tested_ProteinGroupID",
      "n_DAP_FDR05", "n_higher_in_SUS", "n_higher_in_RES", "fraction_DAP_of_tested"
    )
  )
}

assert_sus_res_audit_columns <- function(df, required, label) {
  missing <- setdiff(required, names(df))
  if (length(missing)) stop(label, " is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

sus_res_audit_spatial_order <- function(detail) {
  datasets <- c("neuron_neuropil", "neuron_soma", "microglia")
  unlist(lapply(datasets, function(dataset) {
    units <- unique(as.character(detail$spatial_unit[detail$dataset == dataset]))
    units <- anatomical_spatial_unit_levels(units)
    paste(dataset, units, sep = "|")
  }), use.names = FALSE)
}

sus_res_audit_spatial_label <- function(key) {
  dataset <- sub("\\|.*$", "", key)
  unit <- sub("^[^|]+\\|", "", key)
  prefix <- c(neuron_neuropil = "Neuropil", neuron_soma = "Soma", microglia = "Microglia ROI")[[dataset]]
  paste(prefix, clean_spatial_unit_label(unit), sep = " | ")
}

sus_res_audit_diverging_color <- function(value, limit) {
  if (!is.finite(value)) return("#E7E9ED")
  t <- max(-1, min(1, value / limit))
  low <- grDevices::col2rgb("#3B4CC0")[, 1]
  mid <- grDevices::col2rgb("#FFFFFF")[, 1]
  high <- grDevices::col2rgb("#B40426")[, 1]
  if (t < 0) rgb <- mid + (-t) * (low - mid) else rgb <- mid + t * (high - mid)
  grDevices::rgb(rgb[[1]], rgb[[2]], rgb[[3]], maxColorValue = 255)
}

write_sus_res_biological_audit_workbook <- function(
    overview, detail, supported_go_audit, dap_counts, registry, output_file,
    protected_hashes = NULL,
    source_scripts = c(
      "04_differential_expression_enrichment/07_compareGO_spatial_program_atlas.r",
      "04_differential_expression_enrichment/10_sus_res_spatial_dap_atlas.r"
    )) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Writing sus_res_biological_audit.xlsx requires the installed openxlsx package; packages are not installed automatically.", call. = FALSE)
  }
  req <- sus_res_audit_required_columns()
  assert_sus_res_audit_columns(overview, req$overview, "Theme overview")
  assert_sus_res_audit_columns(detail, req$detail, "Theme detail")
  assert_sus_res_audit_columns(supported_go_audit, req$go_audit, "Supported GO audit")
  assert_sus_res_audit_columns(dap_counts, req$dap, "DAP summary")

  overview_human <- overview[, req$overview, drop = FALSE]
  names(overview_human) <- c(
    "Biological domain", "Theme ID", "Manuscript theme", "Role",
    "Supported GO occurrences", "Unique supported GO IDs", "Evaluable spatial units",
    "FDR-supported spatial units", "Descriptive higher-SUS units",
    "Descriptive higher-RES units", "FDR-supported SUS units", "FDR-supported RES units",
    "Descriptive spatial concordance", "Strongest supported representative GO term",
    "Strongest original BH FDR", "Audit note"
  )

  detail_human <- detail[, req$detail, drop = FALSE]
  names(detail_human) <- c(
    "Dataset", "Compartment", "Region", "Layer", "Spatial unit", "Theme ID", "Manuscript theme", "Theme role",
    "N tested theme terms", "Median NES (all theme terms)", "NES Q25", "NES Q75",
    "Positive NES fraction", "Negative NES fraction", "Descriptive direction", "Low term coverage",
    "Min original GO BH FDR", "FDR support present", "N FDR-supported GO terms",
    "N positive supported terms", "N negative supported terms", "Representative supported GO ID",
    "Representative supported GO description", "Representative supported NES",
    "Representative supported BH FDR", "Supported-term direction consistency",
    "Semantic cluster count", "Leading-edge genes", "Leading-edge proteins", "FDR-supported DAP count"
  )

  go_human <- supported_go_audit[, req$go_audit, drop = FALSE]
  names(go_human) <- c(
    "Dataset", "Compartment", "Region", "Layer", "Spatial unit", "GO ID", "GO description",
    "Original NES", "Original BH FDR", "Assignment status", "Manuscript theme(s)", "Theme role(s)",
    "Semantic cluster ID", "Semantic cluster representative", "Leading-edge genes", "Leading-edge proteins"
  )

  dap_human <- dap_counts[, req$dap, drop = FALSE]
  dap_human$proportion_DAP_higher_SUS <- ifelse(dap_human$n_DAP_FDR05 > 0, dap_human$n_higher_in_SUS / dap_human$n_DAP_FDR05, NA_real_)
  dap_human$proportion_DAP_higher_RES <- ifelse(dap_human$n_DAP_FDR05 > 0, dap_human$n_higher_in_RES / dap_human$n_DAP_FDR05, NA_real_)
  names(dap_human) <- c(
    "Dataset", "Dataset label", "Spatial unit", "Tested ProteinGroupIDs", "FDR-supported DAPs",
    "Higher in SUS", "Higher in RES", "Fraction DAP", "Proportion DAP higher SUS", "Proportion DAP higher RES"
  )
  dap_human <- dap_human[order(-dap_human$`FDR-supported DAPs`, dap_human$Dataset, dap_human$`Spatial unit`, method = "radix"), , drop = FALSE]

  wb <- openxlsx::createWorkbook(creator = "Codex downstream presentation layer")
  sheets <- c("Overview", "Spatial Map", "Theme Detail", "GO Term Audit", "DAP Summary", "Provenance")
  tab_colours <- rep("#5B6770", length(sheets))
  invisible(Map(function(sheet, tab_colour) openxlsx::addWorksheet(wb, sheet, gridLines = FALSE, tabColour = tab_colour), sheets, tab_colours))

  title_style <- openxlsx::createStyle(fontName = "Arial", fontSize = 14, fontColour = "#1F2933", textDecoration = "bold", halign = "left", valign = "center")
  note_style <- openxlsx::createStyle(fontName = "Arial", fontSize = 9, fontColour = "#5B6770", textDecoration = "italic", wrapText = TRUE, valign = "top")
  header_style <- openxlsx::createStyle(fontName = "Arial", fontSize = 9, fontColour = "#1F2933", fgFill = "#E9EDF0", textDecoration = "bold", wrapText = TRUE, halign = "center", valign = "center", border = "Bottom", borderColour = "#5B6770", borderStyle = "thin")
  primary_style <- openxlsx::createStyle(fontName = "Arial", fontSize = 9)
  qc_style <- openxlsx::createStyle(fontName = "Arial", fontSize = 9, fgFill = "#F2F3F4")
  supported_border <- openxlsx::createStyle(border = "TopBottomLeftRight", borderColour = "#5B6770", borderStyle = "thin", textDecoration = "bold")
  body_wrap <- openxlsx::createStyle(fontName = "Arial", fontSize = 9, wrapText = TRUE, valign = "top")
  descriptive_block <- openxlsx::createStyle(fontName = "Arial", fontSize = 9, fgFill = "#F7FAFC")
  inferential_block <- openxlsx::createStyle(fontName = "Arial", fontSize = 9)
  supported_fill <- openxlsx::createStyle(fontName = "Arial", fontSize = 9, fgFill = "#E6F0E8", textDecoration = "bold", halign = "center")
  descriptive_only_fill <- openxlsx::createStyle(fontName = "Arial", fontSize = 9, fgFill = "#F7F7F7", fontColour = "#6B737A", halign = "center")
  body_base <- openxlsx::createStyle(fontName = "Arial", fontSize = 9, valign = "top")
  fdr_style <- openxlsx::createStyle(numFmt = "0.00E+00")
  nes_style <- openxlsx::createStyle(numFmt = "0.00")
  fraction_style <- openxlsx::createStyle(numFmt = "0.0%")

  write_table_sheet <- function(sheet, title, note, data, table_name) {
    openxlsx::mergeCells(wb, sheet, cols = 1:max(1, ncol(data)), rows = 1)
    openxlsx::writeData(wb, sheet, title, startRow = 1, startCol = 1)
    openxlsx::addStyle(wb, sheet, title_style, rows = 1, cols = 1:max(1, ncol(data)), gridExpand = TRUE, stack = TRUE)
    openxlsx::setRowHeights(wb, sheet, rows = 1, heights = 23)
    openxlsx::mergeCells(wb, sheet, cols = 1:max(1, ncol(data)), rows = 2)
    openxlsx::writeData(wb, sheet, note, startRow = 2, startCol = 1)
    openxlsx::addStyle(wb, sheet, note_style, rows = 2, cols = 1:max(1, ncol(data)), gridExpand = TRUE, stack = TRUE)
    openxlsx::setRowHeights(wb, sheet, rows = 2, heights = 34)
    openxlsx::writeDataTable(wb, sheet, data, startRow = 4, startCol = 1, tableName = table_name, tableStyle = "TableStyleLight9", withFilter = TRUE)
    if (nrow(data)) openxlsx::addStyle(wb, sheet, body_base, rows = 5:(nrow(data) + 4L), cols = seq_len(ncol(data)), gridExpand = TRUE, stack = TRUE)
    openxlsx::addStyle(wb, sheet, header_style, rows = 4, cols = seq_len(ncol(data)), gridExpand = TRUE, stack = TRUE)
    openxlsx::setRowHeights(wb, sheet, rows = 4, heights = 36)
    openxlsx::freezePane(wb, sheet, firstActiveRow = 5, firstActiveCol = 2)
  }

  write_table_sheet(
    "Overview", "SUS-RES biological audit: ontology themes",
    "Original GO-term BH FDR is the only inferential pathway support. Spatial direction counts summarize median NES across all mapped tested GO-BP terms and are descriptive, not independent replications or theme-level tests.",
    overview_human, "OverviewTable"
  )
  role_col <- match("Role", names(overview_human))
  primary_rows <- which(overview_human$Role == "primary") + 4L
  qc_rows <- which(overview_human$Role == "QC") + 4L
  if (length(primary_rows)) openxlsx::addStyle(wb, "Overview", primary_style, rows = primary_rows, cols = seq_len(ncol(overview_human)), gridExpand = TRUE, stack = TRUE)
  if (length(qc_rows)) openxlsx::addStyle(wb, "Overview", qc_style, rows = qc_rows, cols = seq_len(ncol(overview_human)), gridExpand = TRUE, stack = TRUE)
  openxlsx::setRowHeights(wb, "Overview", rows = 5:(nrow(overview_human) + 4L), heights = 32)
  openxlsx::addStyle(wb, "Overview", fdr_style, rows = 5:(nrow(overview_human) + 4L), cols = match("Strongest original BH FDR", names(overview_human)), gridExpand = TRUE, stack = TRUE)
  openxlsx::setColWidths(wb, "Overview", cols = 1:ncol(overview_human), widths = c(24, 25, 39, 10, 13, 13, 12, 13, 13, 13, 12, 12, 20, 44, 15, 48))
  openxlsx::addStyle(wb, "Overview", body_wrap, rows = 5:(nrow(overview_human) + 4L), cols = c(1:4, 13:16), gridExpand = TRUE, stack = TRUE)

  spatial_keys <- sus_res_audit_spatial_order(detail)
  theme_meta <- unique(detail[c("manuscript_theme_id", "manuscript_theme", "theme_role")])
  registry_order <- unique(registry[c("theme_id", "display_order")])
  theme_meta$order <- registry_order$display_order[match(theme_meta$manuscript_theme_id, registry_order$theme_id)]
  theme_meta <- theme_meta[order(theme_meta$order, method = "radix"), , drop = FALSE]
  spatial_values <- matrix(NA_real_, nrow = nrow(theme_meta), ncol = length(spatial_keys), dimnames = list(theme_meta$manuscript_theme_id, spatial_keys))
  spatial_support <- matrix(FALSE, nrow = nrow(theme_meta), ncol = length(spatial_keys), dimnames = list(theme_meta$manuscript_theme_id, spatial_keys))
  for (i in seq_len(nrow(detail))) {
    r <- match(detail$manuscript_theme_id[[i]], rownames(spatial_values))
    c <- match(paste(detail$dataset[[i]], detail$spatial_unit[[i]], sep = "|"), colnames(spatial_values))
    if (!is.na(r) && !is.na(c)) {
      spatial_values[r, c] <- detail$median_NES_all_theme_terms[[i]]
      spatial_support[r, c] <- isTRUE(detail$FDR_support_present[[i]])
    }
  }
  spatial_table <- data.frame(
    `Manuscript theme` = theme_meta$manuscript_theme,
    Role = ifelse(theme_meta$theme_role == "qc_review", "QC", theme_meta$theme_role),
    spatial_values,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(spatial_table)[-(1:2)] <- vapply(spatial_keys, sus_res_audit_spatial_label, character(1))
  write_table_sheet(
    "Spatial Map", "SUS-RES theme × spatial-unit map",
    "Cell value/background = descriptive median NES across all mapped tested GO-BP terms. Medium black border = at least one constituent GO term with original BH FDR < 0.05. Unbordered colored cells are descriptive only. Grey = not evaluable. QC rows are separated and excluded from primary Panel C.",
    spatial_table, "SpatialMapTable"
  )
  finite_values <- abs(spatial_values[is.finite(spatial_values)])
  limit <- if (length(finite_values)) max(finite_values) else 1
  for (r in seq_len(nrow(spatial_values))) {
    for (c in seq_len(ncol(spatial_values))) {
      excel_row <- r + 4L
      excel_col <- c + 2L
      fill_style <- openxlsx::createStyle(fgFill = sus_res_audit_diverging_color(spatial_values[r, c], limit), numFmt = "0.00", halign = "center")
      openxlsx::addStyle(wb, "Spatial Map", fill_style, rows = excel_row, cols = excel_col, stack = TRUE)
      if (spatial_support[r, c]) openxlsx::addStyle(wb, "Spatial Map", supported_border, rows = excel_row, cols = excel_col, stack = TRUE)
    }
  }
  qc_map_rows <- which(theme_meta$theme_role == "qc_review") + 4L
  if (length(qc_map_rows)) openxlsx::addStyle(wb, "Spatial Map", qc_style, rows = qc_map_rows, cols = 1:2, gridExpand = TRUE, stack = TRUE)
  openxlsx::setRowHeights(wb, "Spatial Map", rows = 5:(nrow(spatial_table) + 4L), heights = 27)
  openxlsx::setColWidths(wb, "Spatial Map", cols = 1, widths = 40)
  openxlsx::setColWidths(wb, "Spatial Map", cols = 2, widths = 9)
  openxlsx::setColWidths(wb, "Spatial Map", cols = 3:ncol(spatial_table), widths = 14)
  openxlsx::addStyle(wb, "Spatial Map", body_wrap, rows = 5:(nrow(spatial_table) + 4L), cols = 1, gridExpand = TRUE, stack = TRUE)

  write_table_sheet(
    "Theme Detail", "Theme detail: descriptive direction and original GO-term support",
    "Descriptive columns use all tested GO-BP terms mapped to the theme. Representative GO fields are populated only when at least one original constituent term has BH FDR < 0.05; no nonsignificant representative is selected.",
    detail_human, "ThemeDetailTable"
  )
  detail_rows <- 5:(nrow(detail_human) + 4L)
  openxlsx::addStyle(wb, "Theme Detail", nes_style, rows = detail_rows, cols = match(c("Median NES (all theme terms)", "NES Q25", "NES Q75", "Representative supported NES"), names(detail_human)), gridExpand = TRUE, stack = TRUE)
  openxlsx::addStyle(wb, "Theme Detail", fraction_style, rows = detail_rows, cols = match(c("Positive NES fraction", "Negative NES fraction"), names(detail_human)), gridExpand = TRUE, stack = TRUE)
  openxlsx::addStyle(wb, "Theme Detail", fdr_style, rows = detail_rows, cols = match(c("Min original GO BH FDR", "Representative supported BH FDR"), names(detail_human)), gridExpand = TRUE, stack = TRUE)
  descriptive_cols <- match(c("N tested theme terms", "Median NES (all theme terms)", "NES Q25", "NES Q75", "Positive NES fraction", "Negative NES fraction", "Descriptive direction", "Low term coverage"), names(detail_human))
  inferential_cols <- match(c("Min original GO BH FDR", "FDR support present", "N FDR-supported GO terms", "N positive supported terms", "N negative supported terms", "Representative supported GO ID", "Representative supported GO description", "Representative supported NES", "Representative supported BH FDR", "Supported-term direction consistency", "Semantic cluster count", "FDR-supported DAP count"), names(detail_human))
  openxlsx::addStyle(wb, "Theme Detail", descriptive_block, rows = detail_rows, cols = descriptive_cols, gridExpand = TRUE, stack = TRUE)
  openxlsx::addStyle(wb, "Theme Detail", inferential_block, rows = detail_rows, cols = inferential_cols, gridExpand = TRUE, stack = TRUE)
  support_col <- match("FDR support present", names(detail_human))
  supported_detail_rows <- detail_rows[isTRUE(detail_human$`FDR support present`) | detail_human$`FDR support present` %in% c(TRUE, "TRUE")]
  unsupported_detail_rows <- setdiff(detail_rows, supported_detail_rows)
  if (length(supported_detail_rows)) openxlsx::addStyle(wb, "Theme Detail", supported_fill, rows = supported_detail_rows, cols = support_col, stack = TRUE)
  if (length(unsupported_detail_rows)) openxlsx::addStyle(wb, "Theme Detail", descriptive_only_fill, rows = unsupported_detail_rows, cols = support_col, stack = TRUE)
  qc_detail_rows <- detail_rows[detail_human$`Theme role` == "qc_review"]
  if (length(qc_detail_rows)) openxlsx::addStyle(wb, "Theme Detail", qc_style, rows = qc_detail_rows, cols = c(1:8), gridExpand = TRUE, stack = TRUE)
  openxlsx::setColWidths(wb, "Theme Detail", cols = 1:ncol(detail_human), widths = 14)
  openxlsx::setColWidths(wb, "Theme Detail", cols = c(7, 23, 28, 29), widths = c(38, 48, 58, 58))
  openxlsx::addStyle(wb, "Theme Detail", body_wrap, rows = detail_rows, cols = c(7, 23, 28, 29), gridExpand = TRUE, stack = TRUE)
  openxlsx::setRowHeights(wb, "Theme Detail", rows = detail_rows, heights = 34)

  write_table_sheet(
    "GO Term Audit", "Complete FDR-supported SUS-RES GO-term audit",
    "Every original FDR-supported SUS-RES GO occurrence is retained. Unclassified terms remain visible; QC-associated terms are highlighted separately from primary manuscript themes.",
    go_human, "GOTermAuditTable"
  )
  go_rows <- 5:(nrow(go_human) + 4L)
  openxlsx::addStyle(wb, "GO Term Audit", nes_style, rows = go_rows, cols = match("Original NES", names(go_human)), gridExpand = TRUE, stack = TRUE)
  openxlsx::addStyle(wb, "GO Term Audit", fdr_style, rows = go_rows, cols = match("Original BH FDR", names(go_human)), gridExpand = TRUE, stack = TRUE)
  unclassified_rows <- which(go_human$`Assignment status` == "unclassified") + 4L
  qc_go_rows <- which(go_human$`Assignment status` == "qc_review") + 4L
  primary_go_rows <- setdiff(seq_len(nrow(go_human)) + 4L, c(unclassified_rows, qc_go_rows))
  if (length(unclassified_rows)) openxlsx::addStyle(wb, "GO Term Audit", openxlsx::createStyle(fgFill = "#F2F2F2"), rows = unclassified_rows, cols = 1:ncol(go_human), gridExpand = TRUE, stack = TRUE)
  if (length(qc_go_rows)) openxlsx::addStyle(wb, "GO Term Audit", qc_style, rows = qc_go_rows, cols = 1:ncol(go_human), gridExpand = TRUE, stack = TRUE)
  if (length(primary_go_rows)) openxlsx::addStyle(wb, "GO Term Audit", primary_style, rows = primary_go_rows, cols = 1:ncol(go_human), gridExpand = TRUE, stack = TRUE)
  openxlsx::setColWidths(wb, "GO Term Audit", cols = 1:ncol(go_human), widths = 14)
  openxlsx::setColWidths(wb, "GO Term Audit", cols = c(7, 11, 14, 15, 16), widths = c(48, 38, 46, 58, 58))
  openxlsx::addStyle(wb, "GO Term Audit", body_wrap, rows = go_rows, cols = c(7, 11, 14, 15, 16), gridExpand = TRUE, stack = TRUE)
  openxlsx::setRowHeights(wb, "GO Term Audit", rows = go_rows, heights = 38)

  write_table_sheet(
    "DAP Summary", "SUS-RES protein/DAP summary",
    "Protein evidence is based on individual ProteinGroupIDs at BH FDR <= 0.05 with no added fold-change threshold. DAP-set ORA and overlap summaries remain supporting/descriptive and are not the source of ranked-GSEA Panel C.",
    dap_human, "DAPSummaryTable"
  )
  dap_rows <- 5:(nrow(dap_human) + 4L)
  openxlsx::addStyle(wb, "DAP Summary", fraction_style, rows = dap_rows, cols = match(c("Fraction DAP", "Proportion DAP higher SUS", "Proportion DAP higher RES"), names(dap_human)), gridExpand = TRUE, stack = TRUE)
  high_dap_rows <- dap_rows[as.numeric(dap_human$`FDR-supported DAPs`) >= 10]
  if (length(high_dap_rows)) openxlsx::addStyle(wb, "DAP Summary", openxlsx::createStyle(fgFill = "#E8F2FA", textDecoration = "bold"), rows = high_dap_rows, cols = 1:ncol(dap_human), gridExpand = TRUE, stack = TRUE)
  openxlsx::setColWidths(wb, "DAP Summary", cols = 1:ncol(dap_human), widths = c(18, 20, 16, 18, 18, 14, 14, 14, 18, 18))
  openxlsx::setRowHeights(wb, "DAP Summary", rows = dap_rows, heights = 22)

  ontology <- go_ontology_provenance()
  git_head <- tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)[[1]], error = function(e) NA_character_)
  metadata <- data.frame(
    Field = c(
      "Workbook purpose", "Inferential support definition", "Descriptive score definition",
      "GO.db version", "GO source name", "GO source URL", "GO source date",
      "Approved ontology relationships", "Registry version", "Source scripts",
      "Repository HEAD at generation", "Workbook generated at", "Workbook generator"
    ),
    Value = c(
      "Human-facing scientific audit; canonical CSVs remain computational contracts",
      "At least one original ranked-GSEA GO term with BH p.adjust < 0.05",
      "Median original NES across all mapped tested GO-BP terms; no theme-level p-value or FDR",
      ontology$go_db_package_version, ontology$go_source_name, ontology$go_source_url,
      ontology$go_source_date, ontology$approved_relationships,
      unique(registry$registry_version)[[1]], paste(source_scripts, collapse = "; "),
      git_head, format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
      "R/sus_res_biological_audit_workbook.R using openxlsx fallback"
    ),
    stringsAsFactors = FALSE
  )
  registry_human <- registry[c("theme_id", "display_label", "theme_role", "anchor_go_id", "anchor_label", "match_scope", "rationale", "registry_version")]
  names(registry_human) <- c("Theme ID", "Display label", "Theme role", "Anchor GO ID", "Anchor label", "Match scope", "Rationale", "Registry version")
  openxlsx::mergeCells(wb, "Provenance", cols = 1:8, rows = 1)
  openxlsx::writeData(wb, "Provenance", "Ontology, registry, source, and hash provenance", startRow = 1, startCol = 1)
  openxlsx::addStyle(wb, "Provenance", title_style, rows = 1, cols = 1:8, gridExpand = TRUE, stack = TRUE)
  openxlsx::writeDataTable(wb, "Provenance", metadata, startRow = 3, startCol = 1, tableName = "ProvenanceMetadataTable", tableStyle = "TableStyleLight9")
  registry_start <- nrow(metadata) + 6L
  openxlsx::writeData(wb, "Provenance", "Manuscript GO-theme registry", startRow = registry_start, startCol = 1)
  openxlsx::addStyle(wb, "Provenance", openxlsx::createStyle(textDecoration = "bold", fontSize = 11, fontColour = "#23384D"), rows = registry_start, cols = 1)
  openxlsx::writeDataTable(wb, "Provenance", registry_human, startRow = registry_start + 1L, startCol = 1, tableName = "ThemeRegistryTable", tableStyle = "TableStyleLight9")
  hash_start <- registry_start + nrow(registry_human) + 4L
  if (!is.null(protected_hashes) && nrow(protected_hashes)) {
    openxlsx::writeData(wb, "Provenance", "Protected artifact hash comparison", startRow = hash_start, startCol = 1)
    openxlsx::writeDataTable(wb, "Provenance", protected_hashes, startRow = hash_start + 1L, startCol = 1, tableName = "ProtectedHashTable", tableStyle = "TableStyleLight9")
  }
  openxlsx::freezePane(wb, "Provenance", firstActiveRow = 3)
  openxlsx::setColWidths(wb, "Provenance", cols = 1:8, widths = c(28, 58, 16, 16, 38, 24, 70, 24))
  openxlsx::addStyle(wb, "Provenance", body_wrap, rows = 3:(hash_start + if (is.null(protected_hashes)) 1L else nrow(protected_hashes) + 2L), cols = 1:8, gridExpand = TRUE, stack = TRUE)

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  openxlsx::saveWorkbook(wb, output_file, overwrite = TRUE)
  invisible(output_file)
}

validate_sus_res_biological_audit_workbook <- function(output_file, detail_source) {
  if (!requireNamespace("readxl", quietly = TRUE)) stop("Workbook validation requires installed readxl.", call. = FALSE)
  expected_sheets <- c("Overview", "Spatial Map", "Theme Detail", "GO Term Audit", "DAP Summary", "Provenance")
  actual_sheets <- readxl::excel_sheets(output_file)
  if (!identical(actual_sheets, expected_sheets)) stop("Workbook sheet contract failed.", call. = FALSE)
  detail <- readxl::read_excel(output_file, sheet = "Theme Detail", skip = 3)
  required <- c("Dataset", "Spatial unit", "Theme ID", "Median NES (all theme terms)", "FDR support present")
  missing <- setdiff(required, names(detail))
  if (length(missing)) stop("Workbook Theme Detail is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  source_key <- paste(detail_source$dataset, detail_source$spatial_unit, detail_source$manuscript_theme_id, sep = "|")
  workbook_key <- paste(detail$Dataset, detail$`Spatial unit`, detail$`Theme ID`, sep = "|")
  idx <- match(workbook_key, source_key)
  if (anyNA(idx) || !isTRUE(all.equal(as.numeric(detail$`Median NES (all theme terms)`), detail_source$median_NES_all_theme_terms[idx], tolerance = 1e-12, check.attributes = FALSE))) {
    stop("Workbook descriptive spatial values do not match the canonical summary CSV.", call. = FALSE)
  }
  error_tokens <- c("#REF!", "#VALUE!", "#DIV/0!", "#NAME?", "#N/A")
  token_hits <- unlist(lapply(expected_sheets, function(sheet) {
    x <- readxl::read_excel(output_file, sheet = sheet, col_names = FALSE, .name_repair = "minimal")
    any(vapply(x, function(column) any(as.character(column) %in% error_tokens, na.rm = TRUE), logical(1)))
  }))
  if (any(token_hits)) stop("Workbook contains spreadsheet error tokens.", call. = FALSE)
  data.frame(
    sheet = expected_sheets,
    data_rows = c(
      nrow(readxl::read_excel(output_file, sheet = "Overview", skip = 3)),
      nrow(readxl::read_excel(output_file, sheet = "Spatial Map", skip = 3)),
      nrow(detail),
      nrow(readxl::read_excel(output_file, sheet = "GO Term Audit", skip = 3)),
      nrow(readxl::read_excel(output_file, sheet = "DAP Summary", skip = 3)),
      nrow(readxl::read_excel(output_file, sheet = "Provenance", skip = 2, .name_repair = "minimal"))
    ),
    validation_status = "passed",
    stringsAsFactors = FALSE
  )
}

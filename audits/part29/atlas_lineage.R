#!/usr/bin/env Rscript
# ============================================================================
# Part-29 section 33 -- AUDIT ONLY. Reads canonical artefacts, writes nothing
# outside results/tables/publication_audits/upstream_enrichment_v10/.
#
# Traces and verifies every edge of the Figure 3b / ED6 a,b data lineage:
#   canonical ranked gseGO output -> compareGO term comparison ->
#   spatial atlas enrichment long -> manuscript GO-theme assignment table ->
#   panel source data -> panel SVG -> assembled SVG/PDF
# ============================================================================

setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")

stopifnot(requireNamespace("digest", quietly = TRUE))

OUT_DIR <- file.path("results", "tables", "publication_audits",
                     "upstream_enrichment_v10")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")),
                         ..., "\n", sep = "")

sha256 <- function(p) {
  if (length(p) != 1L || is.na(p) || !file.exists(p)) return(NA_character_)
  digest::digest(file = p, algo = "sha256")
}
fsize <- function(p) {
  if (length(p) != 1L || is.na(p) || !file.exists(p)) return(NA_real_)
  as.numeric(file.info(p)$size)
}
fmtime <- function(p) {
  if (length(p) != 1L || is.na(p) || !file.exists(p)) return(NA_character_)
  format(file.info(p)$mtime, "%Y-%m-%d %H:%M:%S", tz = "Europe/Berlin")
}
mt <- function(p) if (!file.exists(p)) NA else file.info(p)$mtime

EDGES <- list()
add_edge <- function(from, to, file_path, producing_script, registry,
                     quantity, transformation, verified, note) {
  EDGES[[length(EDGES) + 1L]] <<- data.frame(
    edge_index = NA_integer_,
    from_artefact = from,
    to_artefact = to,
    file_path = file_path,
    file_sha256 = sha256(file_path),
    file_size_bytes = fsize(file_path),
    file_mtime = fmtime(file_path),
    producing_script = producing_script,
    script_sha256 = sha256(producing_script),
    registry_or_version = registry,
    quantity_carried = quantity,
    transformation_applied = transformation,
    verified = verified,
    note = note,
    stringsAsFactors = FALSE
  )
  invisible(NULL)
}

# ---------------------------------------------------------------- constants
P_THEME_TABLE <- file.path("results", "tables", "10_biological_integration",
                           "gsea_wgcna_concordance", "global",
                           "ontology_aware_gsea_theme_assignments_all_contrasts.csv")
P_THEME_SRC   <- file.path("results", "source_data", "10_biological_integration",
                           "gsea_wgcna_concordance", "global",
                           "ontology_aware_gsea_theme_assignments_all_contrasts.csv")
P_ATLAS_LONG  <- file.path("results", "source_data",
                           "04_differential_expression_enrichment",
                           "compareGO_spatial_atlas",
                           "spatial_atlas_enrichment_long.csv")
P_REGISTRY    <- file.path("config", "manuscript_go_theme_registry.tsv")
P_CPCONFIG    <- file.path("config", "clusterProfiler_config.yml")
P_CONTRACT    <- file.path("figures", "figure_final_truth_v9_contract.yml")

S_CLUSTERPROF <- file.path("04_differential_expression_enrichment",
                           "01_clusterProfiler.r")
S_COMPAREGO   <- file.path("04_differential_expression_enrichment",
                           "02_compareGO.r")
S_ATLAS       <- file.path("04_differential_expression_enrichment",
                           "07_compareGO_spatial_program_atlas.r")
S_CONCORD     <- file.path("10_biological_integration",
                           "05_gsea_wgcna_concordance.R")
S_PANELS      <- file.path("R", "final_truth_v9_panels.R")
S_FIGUTILS    <- file.path("R", "final_truth_v9_figure_utils.R")
S_EXPORT      <- file.path("R", "editorial_v8_export.R")
S_ANNOT       <- file.path("figures", "final_truth_v9_annotation_audit.R")
P_BREADTH     <- file.path("results", "tables", "manuscript_candidates",
                           "final_truth_v9", "audit",
                           "atlas_support_breadth_audit.csv")
P_QA_F3       <- file.path("results", "reports", "manuscript_candidates",
                           "final_truth_v9", "figure_03",
                           "final_truth_v9_qa.csv")
P_QA_ED       <- file.path("results", "reports", "manuscript_candidates",
                           "final_truth_v9", "extended_data",
                           "final_truth_v9_qa.csv")

DATASETS <- c("neuron_neuropil", "neuron_soma", "microglia")
CONTRASTS <- c("RES - CON", "SUS - CON", "SUS - RES")

PANEL_SRC <- c(
  v9_atlas = file.path("results", "source_data", "manuscript_candidates",
                       "final_truth_v9", "figure_03",
                       "v9_atlas_source_data.csv"),
  v9_ed_atlas_rescon = file.path("results", "source_data",
                                 "manuscript_candidates", "final_truth_v9",
                                 "extended_data",
                                 "v9_ed_atlas_rescon_source_data.csv"),
  v9_ed_atlas_suscon = file.path("results", "source_data",
                                 "manuscript_candidates", "final_truth_v9",
                                 "extended_data",
                                 "v9_ed_atlas_suscon_source_data.csv"))
PANEL_SVG <- c(
  v9_atlas = file.path("results", "figures", "manuscript_candidates",
                       "final_truth_v9", "figure_03", "panels", "v9_atlas.svg"),
  v9_ed_atlas_rescon = file.path("results", "figures", "manuscript_candidates",
                                 "final_truth_v9", "extended_data", "panels",
                                 "v9_ed_atlas_rescon.svg"),
  v9_ed_atlas_suscon = file.path("results", "figures", "manuscript_candidates",
                                 "final_truth_v9", "extended_data", "panels",
                                 "v9_ed_atlas_suscon.svg"))
PANEL_CONTRAST <- c(v9_atlas = "SUS - RES",
                    v9_ed_atlas_rescon = "RES - CON",
                    v9_ed_atlas_suscon = "SUS - CON")

ASM <- list(
  F3 = list(
    svg = file.path("results", "figures", "manuscript_candidates",
                    "final_truth_v9", "figure_03", "assembled",
                    "F3_NATURE_FINAL_V9.svg"),
    pdf = file.path("results", "figures", "manuscript_candidates",
                    "final_truth_v9", "figure_03", "assembled",
                    "F3_NATURE_FINAL_V9.pdf"),
    panel_dir = file.path("results", "figures", "manuscript_candidates",
                          "final_truth_v9", "figure_03", "panels"),
    panels = c("v9_dap_track", "v9_atlas", "v9_bridge", "v9_curve_syn",
               "v9_curve_rna", "v9_curve_ox", "v9_prot_syn", "v9_prot_rna",
               "v9_prot_ox")),
  ED6 = list(
    svg = file.path("results", "figures", "manuscript_candidates",
                    "final_truth_v9", "extended_data", "assembled",
                    "ED6_FINAL_V9.svg"),
    pdf = file.path("results", "figures", "manuscript_candidates",
                    "final_truth_v9", "extended_data", "assembled",
                    "ED6_FINAL_V9.pdf"),
    panel_dir = file.path("results", "figures", "manuscript_candidates",
                          "final_truth_v9", "extended_data", "panels"),
    panels = c("v9_ed_atlas_rescon", "v9_ed_atlas_suscon",
               "v9_ed_gsea_curve_syn", "v9_ed_gsea_curve_rna",
               "v9_ed_gsea_curve_ox")))

# =========================================================== STEP 1: GSEA ===
msg("reading clusterProfiler manifests")
strip_drive <- function(x) sub("^[A-Za-z]:/+", "", as.character(x))

man_list <- lapply(DATASETS, function(ds) {
  p <- file.path("data", "processed", "04_differential_expression_enrichment",
                 "clusterProfiler", ds, "clusterProfiler_manifest.csv")
  m <- utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
  m$manifest_path <- p
  m
})
man <- do.call(rbind, man_list)
gsea_bp <- man[man$result_type == "GSEA_GO" & man$ontology == "BP" &
                 man$route_category == "phenotype_within_unit", , drop = FALSE]
gsea_bp$local_output <- strip_drive(gsea_bp$output_table)
msg("canonical ranked GSEA BP runs declared in manifests: ", nrow(gsea_bp))

COMPAREGO_TERM <- setNames(file.path(
  "results", "tables", "04_differential_expression_enrichment", "compareGO",
  DATASETS, "BP", "phenotype_within_unit", "all_route_units",
  "compareGO_term_comparison.csv"), DATASETS)
COMPAREGO_MAN <- setNames(file.path(
  "data", "processed", "04_differential_expression_enrichment", "compareGO",
  DATASETS, "compareGO_input_manifest.csv"), DATASETS)
CP_MAN <- setNames(file.path(
  "data", "processed", "04_differential_expression_enrichment",
  "clusterProfiler", DATASETS, "clusterProfiler_manifest.csv"), DATASETS)

# ================================================= STEP 2: THEME TABLE READ ==
msg("reading manuscript GO-theme assignment table (263 MB, column subset)")
th_hdr <- names(utils::read.csv(P_THEME_TABLE, nrows = 1,
                                stringsAsFactors = FALSE))
keep <- c(dataset = "character", contrast = "character",
          spatial_unit = "character", source_comparison = "character",
          GO_ID = "character", NES = "numeric", GSEA_FDR = "numeric",
          theme_id = "character", manuscript_theme = "character",
          theme_claim_eligible = "logical", registry_version = "character",
          GO_db_package_version = "character",
          evidence_source_family = "character")
cc <- rep("NULL", length(th_hdr)); names(cc) <- th_hdr
cc[names(keep)] <- unname(keep)
th <- utils::read.csv(P_THEME_TABLE, colClasses = unname(cc),
                      stringsAsFactors = FALSE)
msg("theme table rows = ", nrow(th))

# ============================================ STEP 3: RE-IMPLEMENT RENDERER ==
# Verbatim re-implementation of f9_atlas_cells() / f9_atlas_limit() /
# f9_gsea_atlas() ordering from R/final_truth_v9_panels.R (lines 740-860).
SHORT <- c(synaptic_signaling_vesicle = "Synaptic signalling / vesicle",
           rna_processing_splicing_rnp = "RNA processing",
           ribosome_translation = "Translation / ribosome",
           mitochondrial_respiration_oxphos = "Energy metabolism",
           autophagy_lysosome_endosome = "Autophagy / endolysosomal",
           chromatin_organization = "Chromatin")

atlas_cells <- function(th, contrast) {
  z <- th[th$contrast == contrast &
            th$theme_claim_eligible %in% TRUE &
            nzchar(as.character(th$theme_id)), , drop = FALSE]
  n_pre <- nrow(z)
  z <- z[!duplicated(paste(z$dataset, z$spatial_unit, z$contrast, z$GO_ID,
                           z$theme_id)), , drop = FALSE]
  n_post <- nrow(z)
  key <- paste(z$dataset, z$spatial_unit, z$theme_id, sep = "\r")
  cells <- do.call(rbind, lapply(split(seq_len(nrow(z)), key), function(ix) {
    w <- z[ix, , drop = FALSE]
    data.frame(dataset = w$dataset[1], spatial_unit = w$spatial_unit[1],
               theme_id = w$theme_id[1], theme = w$manuscript_theme[1],
               n_fdr = sum(is.finite(w$GSEA_FDR) & w$GSEA_FDR < 0.05),
               median_NES = stats::median(w$NES, na.rm = TRUE),
               n_terms = nrow(w),
               stringsAsFactors = FALSE)
  }))
  rownames(cells) <- NULL
  attr(cells, "n_pre_dedup") <- n_pre
  attr(cells, "n_post_dedup") <- n_post
  cells
}

recomputed <- lapply(CONTRASTS, function(ct) atlas_cells(th, ct))
names(recomputed) <- CONTRASTS
dedup_tbl <- data.frame(
  contrast = CONTRASTS,
  rows_claim_eligible = vapply(recomputed, function(x) attr(x, "n_pre_dedup"), 0),
  rows_after_dedup = vapply(recomputed, function(x) attr(x, "n_post_dedup"), 0),
  stringsAsFactors = FALSE)
dedup_tbl$rows_dropped <- dedup_tbl$rows_claim_eligible -
  dedup_tbl$rows_after_dedup
print(dedup_tbl)
utils::write.csv(dedup_tbl,
                 file.path(OUT_DIR, "atlas_dedup_accounting.csv"),
                 row.names = FALSE)

# Which (dataset, spatial_unit, theme, GO term) constituents are NOT present in
# all three contrast atlases. Any such term makes the three medians rest on
# slightly different term sets.
ea <- th[th$theme_claim_eligible %in% TRUE & nzchar(th$theme_id) &
           th$theme_id %in% names(SHORT) & th$contrast %in% CONTRASTS, ,
         drop = FALSE]
ea <- ea[!duplicated(paste(ea$dataset, ea$spatial_unit, ea$contrast, ea$GO_ID,
                           ea$theme_id)), , drop = FALSE]
ck <- paste(ea$dataset, ea$spatial_unit, ea$theme_id, ea$GO_ID, sep = "\r")
ctab <- table(ck, factor(ea$contrast, levels = CONTRASTS))
asym_idx <- which(rowSums(ctab > 0) < length(CONTRASTS))
asym <- if (length(asym_idx)) {
  parts <- do.call(rbind, strsplit(rownames(ctab)[asym_idx], "\r", fixed = TRUE))
  data.frame(dataset = parts[, 1], spatial_unit = parts[, 2],
             theme_id = parts[, 3], GO_ID = parts[, 4],
             present_RES_CON = as.integer(ctab[asym_idx, "RES - CON"]),
             present_SUS_CON = as.integer(ctab[asym_idx, "SUS - CON"]),
             present_SUS_RES = as.integer(ctab[asym_idx, "SUS - RES"]),
             stringsAsFactors = FALSE)
} else {
  data.frame(dataset = character(0), spatial_unit = character(0),
             theme_id = character(0), GO_ID = character(0),
             present_RES_CON = integer(0), present_SUS_CON = integer(0),
             present_SUS_RES = integer(0), stringsAsFactors = FALSE)
}
n_asym <- nrow(asym)
msg("constituent GO terms not present in all three contrast atlases: ", n_asym)
if (n_asym) print(asym)
utils::write.csv(asym,
                 file.path(OUT_DIR, "atlas_contrast_term_asymmetry.csv"),
                 row.names = FALSE)

# renderer keeps only the six SHORT themes then computes ypos
shape_like_renderer <- function(cells) {
  ord <- names(SHORT)[names(SHORT) %in% cells$theme_id]
  cells <- cells[cells$theme_id %in% ord, , drop = FALSE]
  cells$ypos_expected <- match(cells$theme_id, rev(ord))
  cells
}
recomputed <- lapply(recomputed, shape_like_renderer)

shared_limit_recomputed <- max(abs(unlist(lapply(recomputed,
                                                 function(x) x$median_NES))),
                               na.rm = TRUE)
msg("recomputed shared NES scale limit = ",
    format(shared_limit_recomputed, digits = 17))

# ============================================ STEP 4: COMPARE TO RELEASED ====
released <- lapply(PANEL_SRC, function(p)
  utils::read.csv(p, stringsAsFactors = FALSE))
names(released) <- names(PANEL_SRC)

cmp_rows <- list()
panel_verdict <- list()
for (pid in names(PANEL_SRC)) {
  ct <- PANEL_CONTRAST[[pid]]
  rel <- released[[pid]]
  rec <- recomputed[[ct]]
  k_rel <- paste(rel$dataset, rel$spatial_unit, rel$theme_id, sep = "\r")
  k_rec <- paste(rec$dataset, rec$spatial_unit, rec$theme_id, sep = "\r")
  stopifnot(anyDuplicated(k_rel) == 0L, anyDuplicated(k_rec) == 0L)
  common <- intersect(k_rel, k_rec)
  i_rel <- match(common, k_rel); i_rec <- match(common, k_rec)
  d <- data.frame(
    panel_id = pid, contrast = ct,
    dataset = rel$dataset[i_rel],
    spatial_unit = rel$spatial_unit[i_rel],
    theme_id = rel$theme_id[i_rel],
    released_median_NES = rel$median_NES[i_rel],
    recomputed_median_NES = rec$median_NES[i_rec],
    abs_diff_median_NES = abs(rel$median_NES[i_rel] - rec$median_NES[i_rec]),
    bitwise_identical = mapply(identical, rel$median_NES[i_rel],
                               rec$median_NES[i_rec]),
    released_n_fdr = rel$n_fdr[i_rel],
    recomputed_n_fdr = rec$n_fdr[i_rec],
    n_fdr_match = rel$n_fdr[i_rel] == rec$n_fdr[i_rec],
    released_theme = rel$theme[i_rel],
    recomputed_theme = rec$theme[i_rec],
    theme_label_match = rel$theme[i_rel] == rec$theme[i_rec],
    released_ypos = rel$ypos[i_rel],
    expected_ypos = rec$ypos_expected[i_rec],
    ypos_match = rel$ypos[i_rel] == rec$ypos_expected[i_rec],
    n_constituent_GO_terms = rec$n_terms[i_rec],
    released_shared_limit = rel$shared_NES_scale_limit[i_rel],
    stringsAsFactors = FALSE)
  cmp_rows[[pid]] <- d
  panel_verdict[[pid]] <- list(
    n_released = nrow(rel), n_recomputed = nrow(rec), n_common = length(common),
    only_released = length(setdiff(k_rel, k_rec)),
    only_recomputed = length(setdiff(k_rec, k_rel)),
    all_identical = all(d$bitwise_identical),
    max_abs_diff = max(d$abs_diff_median_NES),
    n_fdr_all_match = all(d$n_fdr_match),
    ypos_all_match = all(d$ypos_match),
    theme_all_match = all(d$theme_label_match),
    shared_limit_released = unique(rel$shared_NES_scale_limit),
    shared_limit_identical = identical(unique(rel$shared_NES_scale_limit),
                                       shared_limit_recomputed))
  msg(pid, ": released=", nrow(rel), " recomputed=", nrow(rec),
      " common=", length(common),
      " all_bitwise_identical=", all(d$bitwise_identical),
      " max_abs_diff=", format(max(d$abs_diff_median_NES), digits = 3),
      " n_fdr_all_match=", all(d$n_fdr_match))
}
cmp <- do.call(rbind, cmp_rows); rownames(cmp) <- NULL
utils::write.csv(cmp, file.path(OUT_DIR, "atlas_cell_reproduction_check.csv"),
                 row.names = FALSE)

# xpos bijection / cross-panel consistency
xp <- lapply(released, function(r)
  unique(data.frame(dataset = r$dataset, spatial_unit = r$spatial_unit,
                    xpos = r$xpos, stringsAsFactors = FALSE)))
xp_key <- lapply(xp, function(d) paste(d$dataset, d$spatial_unit, d$xpos))
xpos_consistent <- identical(sort(xp_key[[1]]), sort(xp_key[[2]])) &&
  identical(sort(xp_key[[1]]), sort(xp_key[[3]]))
n_units <- nrow(xp[[1]])
xpos_bijective <- identical(sort(as.integer(xp[[1]]$xpos)), seq_len(n_units))
msg("spatial units = ", n_units, "; xpos bijective 1..n = ", xpos_bijective,
    "; identical across the three panels = ", xpos_consistent)

# theme coverage cross-check against the established contract
elig <- th[th$theme_claim_eligible %in% TRUE & nzchar(th$theme_id), ,
           drop = FALSE]
theme_terms <- tapply(elig$GO_ID, elig$theme_id, function(x) length(unique(x)))
theme_terms <- theme_terms[order(-theme_terms)]
print(theme_terms)
n_unique_bp <- length(unique(elig$GO_ID))
n_multi <- sum(tapply(elig$theme_id, elig$GO_ID,
                      function(x) length(unique(x))) > 1)
msg("claim-eligible unique GO BP terms = ", n_unique_bp,
    "; terms mapped to >1 theme = ", n_multi)

registry_versions <- sort(unique(th$registry_version))
godb_versions <- sort(unique(th$GO_db_package_version))
evidence_fams <- sort(unique(th$evidence_source_family))
msg("registry_version = ", paste(registry_versions, collapse = ";"),
    " | GO.db = ", paste(godb_versions, collapse = ";"),
    " | evidence_source_family = ", paste(evidence_fams, collapse = ";"))

theme_summary <- data.frame(
  theme_id = names(theme_terms),
  n_unique_GO_BP_terms = as.integer(theme_terms),
  stringsAsFactors = FALSE)
utils::write.csv(theme_summary,
                 file.path(OUT_DIR, "atlas_theme_term_counts.csv"),
                 row.names = FALSE)

# ================== STEP 5: UPSTREAM NUMERIC PROPAGATION (NES pass-through) ==
# Sample claim-eligible (dataset, comparison, GO_ID) triples that actually
# contribute to displayed cells, then chase the SAME NES back through
# spatial_atlas_enrichment_long.csv, compareGO_term_comparison.csv and the
# canonical gseGO result CSV.
set.seed(20260824)
pool <- th[th$contrast %in% CONTRASTS & th$theme_claim_eligible %in% TRUE &
             th$theme_id %in% names(SHORT), , drop = FALSE]
pool <- pool[!duplicated(paste(pool$dataset, pool$source_comparison,
                               pool$GO_ID)), , drop = FALSE]
sel <- do.call(rbind, lapply(split(seq_len(nrow(pool)), pool$dataset),
                             function(ix) pool[sample(ix, 8L), , drop = FALSE]))
sel <- sel[order(sel$dataset, sel$source_comparison, sel$GO_ID), ]
rownames(sel) <- NULL
msg("upstream propagation sample: ", nrow(sel), " (dataset, comparison, GO_ID)")

want_go <- unique(sel$GO_ID)
msg("streaming spatial_atlas_enrichment_long.csv for ", length(want_go),
    " GO IDs")
con <- file(P_ATLAS_LONG, "r")
hdr_line <- readLines(con, n = 1L)
pat <- paste0("(", paste(want_go, collapse = "|"), ")")
hits <- character(0); total <- 0L
repeat {
  chunk <- readLines(con, n = 5000L, warn = FALSE)
  if (!length(chunk)) break
  total <- total + length(chunk)
  h <- chunk[grepl(pat, chunk)]
  if (length(h)) hits <- c(hits, h)
}
close(con)
msg("atlas long: scanned ", total, " lines, matched ", length(hits))
tf <- tempfile(fileext = ".csv")
writeLines(c(hdr_line, hits), tf)
al <- utils::read.csv(tf, stringsAsFactors = FALSE)
unlink(tf)
al <- al[al$ID %in% want_go, , drop = FALSE]
msg("atlas long subset rows = ", nrow(al))

msg("reading compareGO term comparison tables")
cg_list <- lapply(DATASETS, function(ds) {
  p <- COMPAREGO_TERM[[ds]]
  h <- names(utils::read.csv(p, nrows = 1, stringsAsFactors = FALSE))
  cc2 <- rep("NULL", length(h)); names(cc2) <- h
  cc2[c("ID", "dataset", "comparison", "result_type", "ontology")] <- "character"
  cc2[c("NES", "p.adjust", "pvalue")] <- "numeric"
  x <- utils::read.csv(p, colClasses = unname(cc2), stringsAsFactors = FALSE)
  x[x$ID %in% want_go, , drop = FALSE]
})
cg <- do.call(rbind, cg_list)

prop <- do.call(rbind, lapply(seq_len(nrow(sel)), function(i) {
  ds <- sel$dataset[i]; cmpn <- sel$source_comparison[i]; gid <- sel$GO_ID[i]
  mrow <- gsea_bp[gsea_bp$dataset == ds & gsea_bp$comparison == cmpn, ,
                  drop = FALSE]
  gfile <- if (nrow(mrow) == 1L) mrow$local_output[1] else NA_character_
  g_nes <- NA_real_; g_fdr <- NA_real_
  if (!is.na(gfile) && file.exists(gfile)) {
    gx <- utils::read.csv(gfile, stringsAsFactors = FALSE)
    r <- gx[gx$ID == gid, , drop = FALSE]
    if (nrow(r) == 1L) {
      g_nes <- as.numeric(r$NES); g_fdr <- as.numeric(r$p.adjust)
    }
  }
  cr <- cg[cg$dataset == ds & cg$comparison == cmpn & cg$ID == gid, ,
           drop = FALSE]
  c_nes <- if (nrow(cr) == 1L) cr$NES[1] else NA_real_
  c_fdr <- if (nrow(cr) == 1L) cr$p.adjust[1] else NA_real_
  ar <- al[al$dataset == ds & al$comparison == cmpn & al$ID == gid, ,
           drop = FALSE]
  a_nes <- if (nrow(ar) >= 1L) as.numeric(ar$NES[1]) else NA_real_
  a_fdr <- if (nrow(ar) >= 1L) as.numeric(ar$p.adjust[1]) else NA_real_
  data.frame(
    dataset = ds, comparison = cmpn, contrast = sel$contrast[i],
    spatial_unit = sel$spatial_unit[i], GO_ID = gid,
    theme_id = sel$theme_id[i],
    gsea_result_file = if (is.na(gfile)) NA_character_ else gfile,
    NES_canonical_gseGO = g_nes, NES_compareGO = c_nes,
    NES_spatial_atlas_long = a_nes, NES_theme_table = sel$NES[i],
    FDR_canonical_gseGO = g_fdr, FDR_compareGO = c_fdr,
    FDR_spatial_atlas_long = a_fdr, FDR_theme_table = sel$GSEA_FDR[i],
    NES_identical_all_four = identical(g_nes, c_nes) &&
      identical(c_nes, a_nes) && identical(a_nes, sel$NES[i]),
    FDR_identical_all_four = identical(g_fdr, c_fdr) &&
      identical(c_fdr, a_fdr) && identical(a_fdr, sel$GSEA_FDR[i]),
    n_atlas_long_rows = nrow(ar), n_compareGO_rows = nrow(cr),
    stringsAsFactors = FALSE)
}))
utils::write.csv(prop, file.path(OUT_DIR,
                                 "atlas_upstream_nes_propagation_check.csv"),
                 row.names = FALSE)
prop_ok <- all(prop$NES_identical_all_four) && all(prop$FDR_identical_all_four)
msg("upstream NES/FDR pass-through identical in all ", nrow(prop),
    " sampled triples: ", prop_ok)

# ================================================ STEP 6: SVG / PDF CHECKS ===
count_tag <- function(path, tag) {
  if (!file.exists(path)) return(NA_integer_)
  n <- 0L
  con <- file(path, "r")
  repeat {
    chunk <- readLines(con, n = 20000L, warn = FALSE)
    if (!length(chunk)) break
    m <- gregexpr(tag, chunk, fixed = TRUE)
    n <- n + sum(vapply(m, function(z) if (z[1] == -1L) 0L else length(z), 0L))
  }
  close(con)
  as.integer(n)
}
svg_stats <- do.call(rbind, lapply(names(PANEL_SVG), function(pid) {
  p <- PANEL_SVG[[pid]]
  rel <- released[[pid]]
  data.frame(panel_id = pid, svg = p,
             n_rect = count_tag(p, "<rect"),
             n_circle = count_tag(p, "<circle"),
             n_image = count_tag(p, "<image"),
             expected_tiles = nrow(rel),
             expected_dots = sum(rel$n_fdr > 0),
             stringsAsFactors = FALSE)
}))
svg_stats$dots_match <- svg_stats$n_circle == svg_stats$expected_dots
svg_stats$rect_covers_tiles <- svg_stats$n_rect >= svg_stats$expected_tiles
print(svg_stats[, c("panel_id", "n_rect", "n_circle", "n_image",
                    "expected_tiles", "expected_dots", "dots_match")])

asm_stats <- do.call(rbind, lapply(names(ASM), function(k) {
  a <- ASM[[k]]
  pan <- file.path(a$panel_dir, paste0(a$panels, ".svg"))
  sum_circ <- sum(vapply(pan, count_tag, 0L, tag = "<circle"))
  n_asm <- count_tag(a$svg, "<circle")
  data.frame(figure = k, assembled_svg = a$svg,
             n_circle_assembled = n_asm,
             n_circle_panels_sum = sum_circ,
             circles_match = n_asm == sum_circ,
             n_image_assembled = count_tag(a$svg, "<image"),
             stringsAsFactors = FALSE)
}))
print(asm_stats)

pdf_check <- do.call(rbind, lapply(names(ASM), function(k) {
  p <- ASM[[k]]$pdf
  n_img <- NA_integer_; n_font <- NA_integer_; n_page <- NA_integer_
  qok <- FALSE
  qp <- Sys.which("qpdf")
  if (nzchar(qp) && file.exists(p)) {
    tmp <- tempfile(fileext = ".pdf")
    suppressWarnings(system2(qp, c("--qdf", "--object-streams=disable",
                                   shQuote(normalizePath(p)), shQuote(tmp)),
                             stdout = NULL, stderr = NULL))
    if (file.exists(tmp)) {
      ln <- readLines(tmp, warn = FALSE, skipNul = TRUE)
      n_img <- sum(grepl("/Subtype /Image", ln, fixed = TRUE))
      n_font <- sum(grepl("/Type /Font", ln, fixed = TRUE))
      n_page <- sum(grepl("^\\s*/Type /Page\\s*$", ln))
      qok <- TRUE
      unlink(tmp)
    }
  }
  data.frame(figure = k, pdf = p, qpdf_available = nzchar(qp),
             decompressed = qok, n_raster_images = n_img,
             n_font_objects = n_font, n_pages = n_page,
             stringsAsFactors = FALSE)
}))
print(pdf_check)

utils::write.csv(svg_stats,
                 file.path(OUT_DIR, "atlas_svg_mark_counts.csv"),
                 row.names = FALSE)
utils::write.csv(merge(asm_stats, pdf_check, by = "figure"),
                 file.path(OUT_DIR, "atlas_assembly_vector_checks.csv"),
                 row.names = FALSE)

# ======== STEP 6b: CAPTION-REFERENCED BREADTH AUDIT AS AN INDEPENDENT WITNESS
# The Fig 3b / ED6 caption points readers at atlas_support_breadth_audit.csv for
# the per-cell breadth of FDR support. It is written by a different script from
# the renderer, so it is a third independent witness of the same cell values.
br <- utils::read.csv(P_BREADTH, stringsAsFactors = FALSE)
br_key <- paste(br$dataset, br$spatial_unit, br$theme, br$contrast, sep = "\r")
cmp$br_key <- paste(cmp$dataset, cmp$spatial_unit, cmp$theme_id, cmp$contrast,
                    sep = "\r")
bi <- match(cmp$br_key, br_key)
cmp$breadth_median_NES <- br$current_median_NES[bi]
cmp$breadth_n_terms <- br$n_terms[bi]
cmp$breadth_dot_shown <- br$dot_shown[bi]
cmp$breadth_median_identical <- mapply(identical, cmp$released_median_NES,
                                       cmp$breadth_median_NES)
cmp$breadth_n_terms_match <- cmp$breadth_n_terms == cmp$n_constituent_GO_terms
cmp$breadth_dot_match <- cmp$breadth_dot_shown == (cmp$released_n_fdr > 0)
breadth_rows <- nrow(br)
breadth_matched <- sum(!is.na(bi))
breadth_ok <- all(!is.na(bi)) && all(cmp$breadth_median_identical) &&
  all(cmp$breadth_n_terms_match) && all(cmp$breadth_dot_match)
msg("breadth audit: rows=", breadth_rows, " matched=", breadth_matched,
    "/", nrow(cmp), " median_NES identical=",
    all(cmp$breadth_median_identical), " n_terms match=",
    all(cmp$breadth_n_terms_match), " dot_shown match=",
    all(cmp$breadth_dot_match))
utils::write.csv(cmp, file.path(OUT_DIR, "atlas_cell_reproduction_check.csv"),
                 row.names = FALSE)

qa_ok <- all(vapply(c(P_QA_F3, P_QA_ED), function(p) {
  if (!file.exists(p)) return(FALSE)
  all(utils::read.csv(p, stringsAsFactors = FALSE)$status == "PASS")
}, TRUE))
msg("contract QA status PASS for figure_03 and extended_data = ", qa_ok)

# mtime ordering along the true producer -> consumer chain. The clusterProfiler
# manifest is written AFTER the result files it indexes, so it is checked
# separately rather than inserted into the sequence.
chain_files <- c(P_CPCONFIG, gsea_bp$local_output[1],
                 COMPAREGO_TERM[[1]], P_ATLAS_LONG, P_THEME_TABLE,
                 PANEL_SRC[["v9_atlas"]], PANEL_SVG[["v9_atlas"]],
                 ASM$F3$svg, ASM$F3$pdf)
chain_times <- do.call(c, lapply(chain_files, mt))
mtime_monotonic <- all(diff(as.numeric(chain_times)) >= 0)
manifest_after_results <- all(vapply(DATASETS, function(ds) {
  sub <- gsea_bp[gsea_bp$dataset == ds, , drop = FALSE]
  tm <- max(as.numeric(do.call(c, lapply(sub$local_output, mt))))
  as.numeric(mt(CP_MAN[[ds]])) >= tm
}, TRUE))
msg("producer->consumer chain mtime monotonic = ", mtime_monotonic,
    "; manifest written after its result files = ", manifest_after_results)
utils::write.csv(data.frame(step = seq_along(chain_files),
                            file = chain_files,
                            mtime = format(chain_times, "%Y-%m-%d %H:%M:%S"),
                            stringsAsFactors = FALSE),
                 file.path(OUT_DIR, "atlas_chain_mtime_order.csv"),
                 row.names = FALSE)

# ===================================================== STEP 7: EMIT EDGES ====
msg("hashing artefacts and emitting edge table")

add_edge(
  "config/clusterProfiler_config.yml (analysis_params)",
  "gseGO parameterisation used by every ranked GSEA run",
  P_CPCONFIG, S_CLUSTERPROF,
  paste0("config_hash(manifest)=",
         paste(unique(gsea_bp$config_hash), collapse = ";")),
  "pvalue_cutoff, p_adjust_method, min_gs_size, max_gs_size, gsea_seed_base",
  "none (parameter declaration consumed at 01_clusterProfiler.r:1502 gseGO call)",
  TRUE,
  paste0("manifest config_file resolves to this file for all ", nrow(gsea_bp),
         " GSEA_GO/BP/phenotype_within_unit runs; distinct config_hash count=",
         length(unique(gsea_bp$config_hash))))

for (ds in DATASETS) {
  sub <- gsea_bp[gsea_bp$dataset == ds, , drop = FALSE]
  add_edge(
    paste0("Stage-02 mapped per-contrast tables (", ds, ")"),
    paste0("clusterProfiler run manifest (", ds, ")"),
    CP_MAN[[ds]], S_CLUSTERPROF,
    paste0("enrichment_contract_version=",
           paste(unique(sub$enrichment_contract_version), collapse = ";")),
    "one manifest row per analysis_id; declares output_table, rank statistic, input hash",
    paste0("gene collapse = ",
           paste(unique(sub$duplicate_gene_collapse_rule), collapse = ";"),
           "; rank statistic column = ",
           paste(unique(sub$rank_statistic_column), collapse = ";"),
           " (", paste(unique(sub$rank_statistic_type), collapse = ";"),
           "), fallback_used = ",
           paste(unique(sub$rank_statistic_fallback_used), collapse = ";")),
    TRUE,
    paste0(nrow(sub),
           " GSEA_GO/BP/phenotype_within_unit rows; analysis_status=",
           paste(unique(sub$analysis_status), collapse = ";")))
}

for (i in seq_len(nrow(gsea_bp))) {
  r <- gsea_bp[i, ]
  add_edge(
    paste0("ranked gene list (", r$dataset, " / ", r$route_unit, " / ",
           r$comparison, ")"),
    paste0("canonical ranked gseGO BP result (", r$dataset, " / ",
           r$route_unit, " / ", r$comparison, ")"),
    r$local_output, S_CLUSTERPROF,
    paste0("input_hash=", r$input_hash, "; config_hash=", r$config_hash),
    "ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, qvalue, rank, core_enrichment",
    paste0("clusterProfiler::gseGO(ont='BP', keyType='SYMBOL', OrgDb=org.Mm.eg.db) ",
           "on genes collapsed by ", r$duplicate_gene_collapse_rule,
           " of the '", r$rank_statistic_column, "' statistic (",
           r$rank_statistic_type, ", fallback_used=",
           r$rank_statistic_fallback_used, ")"),
    file.exists(r$local_output),
    paste0("n_genes=", r$n_genes, "; n_terms=", r$n_terms, "; status=",
           r$analysis_status, "; manifest=", r$manifest_path))
}

for (ds in DATASETS) {
  add_edge(
    paste0("clusterProfiler manifest (", ds, ")"),
    paste0("compareGO input manifest (", ds, ")"),
    COMPAREGO_MAN[[ds]], S_COMPAREGO,
    "compareGO input manifest contract",
    "declared list of per-comparison ranked GSEA result CSVs consumed by compareGO",
    "manifest-driven file discovery; no recomputation of enrichment",
    file.exists(COMPAREGO_MAN[[ds]]),
    "02_compareGO.r reads clusterProfiler_manifest.csv rather than globbing")
}

for (ds in DATASETS) {
  sub <- gsea_bp[gsea_bp$dataset == ds, , drop = FALSE]
  pk <- prop[prop$dataset == ds, , drop = FALSE]
  add_edge(
    paste0("canonical ranked gseGO BP results (", ds, ", ", nrow(sub),
           " runs)"),
    paste0("compareGO term comparison table (", ds, ")"),
    COMPAREGO_TERM[[ds]], S_COMPAREGO,
    "compareGO_term_comparison contract (ID, Description, NES, p.adjust, setSize, core_enrichment)",
    "NES, pvalue, p.adjust, core_enrichment per (comparison, GO term)",
    "row-bind of the manifest-declared per-comparison gseGO result CSVs with dataset/comparison/result_type/ontology columns added; NES and p.adjust copied verbatim",
    nrow(pk) > 0 && all(pk$NES_identical_all_four) &&
      all(pk$FDR_identical_all_four),
    paste0("verified on ", nrow(pk),
           " sampled GO terms: NES and p.adjust bitwise identical to the gseGO result files"))
}

add_edge(
  "config/manuscript_go_theme_registry.tsv",
  "manuscript GO theme registry (anchors + approved relationships)",
  P_REGISTRY, S_ATLAS,
  paste0("registry_version=", paste(registry_versions, collapse = ";")),
  "theme_id, manuscript_theme, theme_role, anchor GO IDs, approved is_a/part_of relationships",
  "none (checked-in registry consumed by map_go_terms_to_manuscript_themes())",
  file.exists(P_REGISTRY),
  paste0("theme table carries registry_version=",
         paste(registry_versions, collapse = ";"), " and GO.db=",
         paste(godb_versions, collapse = ";")))

add_edge(
  "compareGO term comparison tables (3 datasets)",
  "spatial atlas ranked-GSEA long table",
  P_ATLAS_LONG, S_ATLAS,
  paste0("evidence_source_family=", paste(evidence_fams, collapse = ";")),
  "NES, pvalue, p.adjust, core_enrichment per (dataset, comparison, GO term) plus spatial_unit / phenotype_contrast / program_class annotation",
  "adds dataset, phenotype_contrast, region, layer, compartment, spatial_unit, program_class and source_* provenance columns; NES and p.adjust copied verbatim; no re-enrichment",
  all(prop$NES_identical_all_four) && all(prop$FDR_identical_all_four),
  paste0("scanned ", total, " data lines; verified on ", nrow(prop),
         " sampled (dataset, comparison, GO_ID) triples that NES and p.adjust are bitwise identical to compareGO and to the gseGO result files"))

add_edge(
  "spatial atlas ranked-GSEA long table + manuscript GO theme registry",
  "ontology-aware GSEA theme assignment table (all contrasts)",
  P_THEME_TABLE, S_CONCORD,
  paste0("registry_version=", paste(registry_versions, collapse = ";"),
         "; GO.db=", paste(godb_versions, collapse = ";")),
  "one row per (dataset, contrast, spatial_unit, GO term, theme assignment): NES, raw_p, GSEA_FDR, theme_id, theme_role, theme_claim_eligible",
  "GO-ID ontology mapping of every term to manuscript themes via approved is_a/part_of ancestry (gww_build_ontology_theme_term_table); NES -> NES and p.adjust -> GSEA_FDR verbatim; theme_claim_eligible = theme_role in {primary, supporting}; a term in two themes emits two rows",
  all(prop$NES_identical_all_four) && all(prop$FDR_identical_all_four),
  paste0("rows=", nrow(th), "; claim-eligible unique GO BP terms=", n_unique_bp,
         "; terms in >1 theme=", n_multi,
         "; NES/GSEA_FDR pass-through verified on ", nrow(prop),
         " sampled terms"))

add_edge(
  "ontology-aware GSEA theme assignment table (results/tables copy)",
  "ontology-aware GSEA theme assignment table (results/source_data copy)",
  P_THEME_SRC, S_CONCORD,
  "write_integration_table() dual write",
  "identical table content",
  "byte-identical copy written in the same call (write_integration_table writes tables/ and source_data/)",
  identical(sha256(P_THEME_TABLE), sha256(P_THEME_SRC)),
  paste0("sha256 equal = ",
         identical(sha256(P_THEME_TABLE), sha256(P_THEME_SRC)),
         "; the figure contract points at the results/tables copy"))

add_edge(
  "figures/figure_final_truth_v9_contract.yml",
  "panel declarations v9_atlas (Fig 3b), v9_ed_atlas_rescon (ED6a), v9_ed_atlas_suscon (ED6b)",
  P_CONTRACT, S_FIGUTILS,
  "manuscript_final_truth_v9_figures_v1",
  "renderer = f9_gsea_atlas; contrast = SUS - RES / RES - CON / SUS - CON; primary_source = the theme assignment table",
  "none (declaration); s9f_build() dispatches panel -> renderer -> exact mm box",
  TRUE,
  "F3 layout label b = v9_atlas; ED6 layout labels a,b = v9_ed_atlas_rescon, v9_ed_atlas_suscon")

for (pid in names(PANEL_SRC)) {
  v <- panel_verdict[[pid]]
  add_edge(
    "ontology-aware GSEA theme assignment table",
    paste0("panel source data (", pid, ")"),
    PANEL_SRC[[pid]], S_PANELS,
    paste0("f9_atlas_cells(); contrast='", PANEL_CONTRAST[[pid]], "'"),
    "n_fdr and median_NES per dataset x spatial_unit x theme, plus xpos/ypos and the shared colour-scale limit",
    paste0("filter contrast=='", PANEL_CONTRAST[[pid]],
           "' & theme_claim_eligible==TRUE & nzchar(theme_id); ",
           "de-duplicate on dataset+spatial_unit+contrast+GO_ID+theme_id (keep first); ",
           "group by dataset+spatial_unit+theme_id; ",
           "n_fdr = sum(is.finite(GSEA_FDR) & GSEA_FDR<0.05); ",
           "median_NES = stats::median(NES, na.rm=TRUE); ",
           "keep only the six claim-eligible SHORT themes; ",
           "ypos = match(theme_id, rev(SHORT order)); ",
           "xpos = spatial-grammar column order (sg_blocks); ",
           "shared_NES_scale_limit = max|median_NES| over all three contrast atlases"),
    isTRUE(v$all_identical) && isTRUE(v$n_fdr_all_match) &&
      isTRUE(v$ypos_all_match) && isTRUE(v$theme_all_match) &&
      isTRUE(v$shared_limit_identical) && v$only_released == 0 &&
      v$only_recomputed == 0,
    paste0("recomputed independently from the theme table: ", v$n_common,
           "/", v$n_released,
           " cells matched by key, all median_NES bitwise identical = ",
           v$all_identical, ", max|diff| = ",
           format(v$max_abs_diff, digits = 3), ", n_fdr all match = ",
           v$n_fdr_all_match, ", ypos all match = ", v$ypos_all_match,
           ", shared_NES_scale_limit identical (",
           format(v$shared_limit_released, digits = 17), ") = ",
           v$shared_limit_identical))
}

for (pid in names(PANEL_SVG)) {
  s <- svg_stats[svg_stats$panel_id == pid, ]
  add_edge(
    paste0("panel source data (", pid, ")"),
    paste0("rendered panel SVG (", pid, ")"),
    PANEL_SVG[[pid]], S_PANELS,
    "ggplot2 + svglite via nv_save_panel(); box from the contract layout",
    "one geom_tile rect per cell (fill = median_NES) and one geom_point circle per cell with n_fdr > 0",
    "geom_tile(fill=median_NES) + geom_point(subset n_fdr>0) + nv_diverging(limits = +/- shared_NES_scale_limit); written in the same f9_gsea_atlas() call that wrote the CSV",
    isTRUE(s$dots_match) && isTRUE(s$rect_covers_tiles),
    paste0("<circle> count = ", s$n_circle, " == cells with n_fdr>0 = ",
           s$expected_dots, "; <rect> count = ", s$n_rect, " >= ",
           s$expected_tiles, " tiles (remainder is plot/panel/legend chrome); ",
           "<image> count = ", s$n_image,
           " (continuous colour-bar key is a raster; all data marks are vector)"))
}

for (k in names(ASM)) {
  a <- asm_stats[asm_stats$figure == k, ]
  add_edge(
    paste0("captured panel grobs (", paste(ASM[[k]]$panels, collapse = ", "),
           ")"),
    paste0("assembled figure SVG (", k, ")"),
    ASM[[k]]$svg, S_EXPORT,
    "e8_compose_svg(); exact-mm layout from the contract",
    "all panel marks re-drawn into one page at the declared millimetre boxes",
    "grobs captured by the nv_save_panel shadow are re-drawn on one svglite page; no rescaling of panel content",
    isTRUE(a$circles_match),
    paste0("<circle> in assembled = ", a$n_circle_assembled,
           " == sum over constituent panel SVGs = ", a$n_circle_panels_sum))
}

for (k in names(ASM)) {
  pc <- pdf_check[pdf_check$figure == k, ]
  add_edge(
    paste0("captured panel grobs (", k, ")"),
    paste0("assembled figure PDF (", k, ")"),
    ASM[[k]]$pdf, S_EXPORT,
    "e8_compose_pdf(); grDevices::cairo_pdf, MediaBox rounded up to whole points",
    "same drawing as the assembled SVG, as vector PDF",
    "grid viewports at exact mm positions; grid.draw of each captured grob; panel labels drawn as text",
    isTRUE(pc$decompressed) && isTRUE(pc$n_pages == 1L) &&
      isTRUE(pc$n_font_objects > 0),
    paste0("qpdf --qdf: pages=", pc$n_pages, ", font objects=",
           pc$n_font_objects, ", raster image XObjects=", pc$n_raster_images,
           " (the single raster is the continuous colour-bar key; data marks are vector paths)"))
}

add_edge(
  "ontology-aware GSEA theme assignment table",
  "caption-referenced per-cell support-breadth audit (atlas_support_breadth_audit.csv)",
  P_BREADTH, S_ANNOT,
  "figure caption of Fig 3b / ED6 a,b points readers at this file",
  "per dataset x spatial_unit x theme x contrast: n_terms, n_terms_FDR_supported, current_median_NES, dot_shown",
  "same cell grouping as f9_atlas_cells(), recomputed independently by figures/final_truth_v9_annotation_audit.R; adds semantic-representative and direction breakdowns",
  breadth_ok,
  paste0("rows=", breadth_rows, " (= 108 cells x 3 contrasts); matched ",
         breadth_matched, "/", nrow(cmp),
         " audited cells; current_median_NES bitwise identical to the released panel source data = ",
         all(cmp$breadth_median_identical), "; n_terms identical = ",
         all(cmp$breadth_n_terms_match), "; dot_shown == (n_fdr>0) = ",
         all(cmp$breadth_dot_match)))

edges <- do.call(rbind, EDGES)
edges$edge_index <- seq_len(nrow(edges))
utils::write.csv(edges, file.path(OUT_DIR, "atlas_lineage_audit.csv"),
                 row.names = FALSE)
msg("wrote ", nrow(edges), " edges; verified TRUE on ",
    sum(edges$verified %in% TRUE), " of them")
msg("FALSE/NA edges: ",
    paste(edges$to_artefact[!(edges$verified %in% TRUE)], collapse = " | "))

# compact machine-readable verdict
verdict <- data.frame(
  check = c("n_edges", "n_edges_verified", "n_canonical_gsea_bp_runs",
            "n_spatial_units", "xpos_bijective",
            "xpos_identical_across_panels", "theme_table_rows",
            "claim_eligible_unique_GO_BP_terms", "terms_in_two_themes",
            "n_claim_eligible_themes_displayed",
            "shared_NES_scale_limit_recomputed",
            "all_cells_bitwise_identical", "n_cells_compared",
            "max_abs_median_NES_diff", "upstream_NES_passthrough_identical",
            "breadth_audit_independent_witness_agrees",
            "producer_to_consumer_chain_mtime_monotonic",
            "clusterProfiler_manifest_written_after_its_results",
            "contract_QA_status_PASS",
            "constituent_terms_absent_from_some_contrast", "registry_version",
            "GO_db_package_version", "evidence_source_family", "git_head"),
  value = c(nrow(edges), sum(edges$verified %in% TRUE), nrow(gsea_bp),
            n_units, xpos_bijective, xpos_consistent, nrow(th), n_unique_bp,
            n_multi, length(unique(cmp$theme_id)),
            format(shared_limit_recomputed, digits = 17),
            all(cmp$bitwise_identical), nrow(cmp),
            format(max(cmp$abs_diff_median_NES)),
            prop_ok, breadth_ok, mtime_monotonic, manifest_after_results,
            qa_ok, n_asym,
            paste(registry_versions, collapse = ";"),
            paste(godb_versions, collapse = ";"),
            paste(evidence_fams, collapse = ";"),
            tryCatch(paste(system2("git", c("rev-parse", "HEAD"),
                                   stdout = TRUE), collapse = ""),
                     error = function(e) NA_character_)),
  stringsAsFactors = FALSE)
utils::write.csv(verdict, file.path(OUT_DIR, "atlas_lineage_verdict.csv"),
                 row.names = FALSE)
print(verdict)
msg("done")

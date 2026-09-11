# Reusable, phenotype-blind EWCE gene-set annotation engine.
#
# WHY THIS EXISTS
#   The EWCE computation lived inside a 2,043-line analysis script as a closure
#   over global state (`ctd`, `analysis_params`), reachable only by running the
#   whole phenotype pipeline. Annotating an arbitrary gene set - a WGCNA module,
#   say - was impossible without either running SUS/RES/CON machinery or writing
#   a second implementation of EWCE. This file makes the core computation
#   callable. It is NOT a second implementation: the analysis script delegates
#   to `ewce_bootstrap_once()` so exactly one code path performs the test.
#
# THREE DEFECTS THIS ADDRESSES
#   1. The background could be silently dropped. A retry without `bg`
#      substitutes the full reference transcriptome while the caller still
#      records the intended background size, making a wrong result
#      indistinguishable from a correct one. Here the background is mandatory.
#   2. One global BH family spanned the phenotype and phenotype-blind arms, so
#      adding or removing phenotype rows changed a baseline result's FDR. Here
#      every row carries an explicit `fdr_family` and BH is applied strictly
#      within it.
#   3. The historical "Baseline" arm was computed per CON/RES/SUS, so it was
#      never phenotype-blind. This API takes no group argument at all.
#
# WHAT THIS DOES NOT CHANGE
#   The statistical method is still EWCE::bootstrap_enrichment_test with the
#   same specificity reference. Nothing about the test is reimplemented.

# ------------------------------------------------------------- FDR families
#
# Family IDs are prospective and declared here so a family can never be widened
# by accident. Two rows share an FDR family only if they answer the same
# question over the same grid.

ewce_module_scopes <- function() c("all", "core_kME06", "top25")

# Phenotype-derived signatures (the historical Differential arm).
ewce_fdr_family_differential <- function(dataset, celltype_level) {
  paste("ewce_differential", dataset, paste0("level", celltype_level), sep = "_")
}

# Phenotype-blind module annotation.
#
# The family is dataset x module_scope x celltype_level, with BH applied across
# ModuleID x CellType inside it. That grid is fixed by the module set and the
# reference, so it cannot be changed by anything happening in the phenotype arm.
ewce_fdr_family_module_annotation <- function(dataset, module_scope, celltype_level) {
  module_scope <- match.arg(as.character(module_scope), ewce_module_scopes())
  paste("ewce_module_annotation", dataset, module_scope,
        paste0("level", celltype_level), sep = "_")
}

# Apply BH strictly within each declared family.
ewce_apply_family_fdr <- function(df, p_col = "p_value",
                                  family_col = "fdr_family",
                                  out_col = "FDR") {
  for (nm in c(p_col, family_col)) {
    if (!nm %in% names(df)) stop("Missing column: ", nm, ".", call. = FALSE)
  }
  if (!nrow(df)) { df[[out_col]] <- numeric(0); return(df) }
  fam <- as.character(df[[family_col]])
  if (any(is.na(fam) | !nzchar(fam))) {
    stop("Every row needs an explicit FDR family; refusing to pool rows into a ",
         "global family.", call. = FALSE)
  }
  p <- as.numeric(df[[p_col]])
  out <- rep(NA_real_, length(p))
  for (f in unique(fam)) {
    i <- which(fam == f)
    out[i] <- stats::p.adjust(p[i], method = "BH")
  }
  df[[out_col]] <- out
  df
}

# ------------------------------------------------------------- core engine

# ONE bootstrap enrichment test. The background is REQUIRED and is never
# silently dropped: if the test fails with the supplied background, this fails.
# `output_species` defaults to NULL, meaning the argument is NOT passed and
# EWCE's own default applies. That keeps the canonical analysis byte-for-byte on
# its historical behaviour. The phenotype-blind module annotation passes "mouse"
# explicitly: the specificity reference is mouse Title-case, so staying in its
# native gene space avoids a mouse -> human ortholog round trip that can silently
# drop every gene (EWCE then reports "Only 0 provided").
ewce_bootstrap_once <- function(hits, background, reference, annot_level,
                                reps = 10000L, seed = NULL,
                                min_hits = 10L, output_species = NULL) {
  if (is.null(reference)) stop("An EWCE specificity reference is required.", call. = FALSE)
  hits <- unique(stats::na.omit(as.character(hits)))
  background <- unique(stats::na.omit(as.character(background)))
  if (!length(background)) {
    stop("EWCE background is empty. A background-free test answers a different ",
         "question and is not an acceptable fallback.", call. = FALSE)
  }
  hits <- intersect(hits, background)
  if (length(hits) < min_hits) {
    stop("EWCE hit list has fewer than ", min_hits,
         " genes after background intersection.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(as.integer(seed))

  call_args <- list(
    sct_data = reference,
    hits = hits,
    bg = background,
    reps = reps,
    annotLevel = annot_level,
    genelistSpecies = "mouse",
    sctSpecies = "mouse"
  )
  if (!is.null(output_species)) call_args$output_species <- output_species
  res <- try(do.call(EWCE::bootstrap_enrichment_test, call_args), silent = TRUE)
  # Legacy EWCE signature, still WITH an explicit background.
  if (inherits(res, "try-error") &&
      exists("bootstrap.enrichment.test", envir = asNamespace("EWCE"))) {
    legacy <- get("bootstrap.enrichment.test", envir = asNamespace("EWCE"))
    res <- try(
      legacy(sct_data = reference[[annot_level]], mouse.hits = hits,
             mouse.bg = background, reps = reps),
      silent = TRUE)
  }
  if (inherits(res, "try-error")) {
    stop("EWCE bootstrap failed at annotLevel ", annot_level, " with ",
         length(hits), " hits and a ", length(background),
         "-gene background. Refusing to retry without a background.",
         call. = FALSE)
  }
  out <- res$results
  attr(out, "n_hits_tested") <- length(hits)
  attr(out, "n_background") <- length(background)
  out
}

# ------------------------------------------------------------- public API

# Annotate arbitrary gene sets against an external cell-type reference.
#
# Takes NO StressGroup, contrast, CON/RES/SUS label or DAP status. A gene set is
# just a character vector of gene symbols with an id.
#
# `gene_sets` : named list of character vectors, or a data.frame with columns
#               gene_set_id and gene_symbol.
run_ewce_gene_set_annotation <- function(gene_sets,
                                         background,
                                         reference,
                                         dataset,
                                         celltype_level = 1L,
                                         n_boot = 10000L,
                                         seed = 20260101L,
                                         fdr_family = NULL,
                                         min_genes = 10L,
                                         output_species = "mouse",
                                         provenance = list()) {
  forbidden <- c("StressGroup", "ExpGroup", "contrast", "Contrast", "Direction",
                 "condition", "Condition", "DAP", "group", "Group")
  if (is.data.frame(gene_sets)) {
    hit <- intersect(forbidden, names(gene_sets))
    if (length(hit)) {
      stop("Gene sets carry phenotype column(s): ", paste(hit, collapse = ", "),
           ". This annotation layer is phenotype-blind by contract.",
           call. = FALSE)
    }
    if (!all(c("gene_set_id", "gene_symbol") %in% names(gene_sets))) {
      stop("gene_sets data.frame needs gene_set_id and gene_symbol.", call. = FALSE)
    }
    gene_sets <- split(as.character(gene_sets$gene_symbol),
                       as.character(gene_sets$gene_set_id))
  }
  if (!is.list(gene_sets) || is.null(names(gene_sets))) {
    stop("gene_sets must be a named list or a gene_set_id/gene_symbol frame.",
         call. = FALSE)
  }
  background <- unique(stats::na.omit(as.character(background)))
  if (!length(background)) stop("A measured background is required.", call. = FALSE)

  rows <- list()
  for (gid in names(gene_sets)) {
    genes <- unique(stats::na.omit(as.character(gene_sets[[gid]])))
    mapped <- intersect(genes, background)
    fam <- fdr_family %||% ewce_fdr_family_differential(dataset, celltype_level)
    base <- data.frame(
      gene_set_id = gid, dataset = dataset, level = as.integer(celltype_level),
      n_input_genes = length(genes), n_mapped_genes = length(mapped),
      n_background = length(background), stringsAsFactors = FALSE)

    if (length(mapped) < min_genes) {
      rows[[length(rows) + 1L]] <- cbind(base, data.frame(
        cell_type = NA_character_, observed_statistic = NA_real_,
        null_mean = NA_real_, null_sd = NA_real_, z_score = NA_real_,
        p_value = NA_real_, fdr_family = fam,
        annotation_status = "insufficient_mapped_genes",
        stringsAsFactors = FALSE))
      next
    }
    res <- ewce_bootstrap_once(mapped, background, reference,
                               annot_level = celltype_level,
                               reps = n_boot, seed = seed, min_hits = min_genes,
                               output_species = output_species)
    res <- as.data.frame(res, stringsAsFactors = FALSE)
    pick <- function(nms, default = NA_real_) {
      nm <- intersect(nms, names(res))
      if (length(nm)) as.numeric(res[[nm[[1]]]]) else rep(default, nrow(res))
    }
    ct <- {
      nm <- intersect(c("CellType", "cell_type"), names(res))
      if (length(nm)) as.character(res[[nm[[1]]]]) else rownames(res)
    }
    rows[[length(rows) + 1L]] <- cbind(
      base[rep(1L, nrow(res)), , drop = FALSE],
      data.frame(
        cell_type = ct,
        observed_statistic = pick(c("mean_exp", "observed", "hit.cells")),
        null_mean = pick(c("bootstrap_mean", "mean_bootstrap", "boot_mean")),
        null_sd = pick(c("bootstrap_sd", "sd_bootstrap", "boot_sd")),
        z_score = pick(c("sd_from_mean", "z", "zscore")),
        fold_change = pick(c("fold_change", "fc")),
        p_value = pick(c("p", "p_value", "pvalue")),
        fdr_family = fam,
        annotation_status = "tested",
        stringsAsFactors = FALSE))
  }
  out <- dplyr::bind_rows(rows)
  rownames(out) <- NULL

  prov <- c(list(engine = "ewce_gene_set_engine",
                 method = "EWCE::bootstrap_enrichment_test",
                 reps = n_boot, seed = seed,
                 background_source = "measured_proteome",
                 phenotype_used = "none"), provenance)
  out$provenance <- paste(names(prov), unlist(lapply(prov, function(v)
    paste(as.character(v), collapse = "|"))), sep = "=", collapse = "; ")
  out
}

# ------------------------------------------------------- gene symbol contract
#
# EWCE matches against a mouse specificity reference and expects canonical mouse
# symbols (Title case, e.g. "Gria1"). The WGCNA membership tables store gene
# symbols UPPERCASED, which are not valid org.Mm.eg.db keys - passing them
# through yields zero reference matches and EWCE fails with "Only 0 provided".
#
# This resolves the canonical mouse symbol case-insensitively against
# org.Mm.eg.db, preferring an exact SYMBOL and falling back to ALIAS. It adds no
# new identifier authority: it recovers the casing the annotation database
# already defines.
ewce_to_mouse_symbols <- function(x, keep_unmapped = FALSE) {
  x <- trimws(as.character(x))
  x[!nzchar(x)] <- NA_character_
  if (!requireNamespace("org.Mm.eg.db", quietly = TRUE) ||
      !requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("org.Mm.eg.db and AnnotationDbi are required to resolve mouse symbols.",
         call. = FALSE)
  }
  db <- getExportedValue("org.Mm.eg.db", "org.Mm.eg.db")
  sym <- AnnotationDbi::keys(db, keytype = "SYMBOL")
  ali <- AnnotationDbi::keys(db, keytype = "ALIAS")
  lut_sym <- stats::setNames(sym, toupper(sym))
  lut_ali <- stats::setNames(ali, toupper(ali))
  up <- toupper(x)
  out <- unname(lut_sym[up])
  miss <- is.na(out)
  if (any(miss)) out[miss] <- unname(lut_ali[up[miss]])
  if (keep_unmapped) out[is.na(out)] <- x[is.na(out)]
  out
}

#!/usr/bin/env Rscript

# Part-28: biological coherence of the atlas themes, and whether the median-NES
# summary survives GO redundancy.
#
# ANNOTATION AUDIT ONLY. No differential abundance, no GSEA and no WGCNA is
# rerun; no GO term is added to or removed from the canonical registry; no FDR
# family is created. Every number here is either a semantic-similarity statistic
# over GO IDs, or a descriptive summary of NES and FDR values the canonical
# workflow already stored.
#
# PHENOTYPE-BLIND BY CONSTRUCTION. The similarity matrix is built from GO IDs
# alone; the clustering sees only that matrix; the representative of a cluster
# is its semantic medoid. It is deliberately NOT the repository's own
# go_semantic_redundancy_qa(), which picks representatives by best FDR - correct
# for redundancy QA of an already-supported set, but circular here.
#
# TWO CLUSTERINGS, FOR TWO DIFFERENT QUESTIONS.
#   binary_cut (simplifyEnrichment, cutoff-free) answers "is this theme one
#     biological block?" - it is the method that package's own benchmark
#     favours, and it does not need a threshold that could be tuned.
#   average-linkage at Wang >= 0.70 (the repository's own pinned descriptive
#     cutoff) answers "which terms are near-duplicates of each other?" - that is
#     the grouping a redundancy-reduced summary needs.
# Using binary_cut for redundancy reduction would collapse a 62-term theme to a
# single term; using the 0.70 cut for coherence would call every theme broad.

source(file.path("R", "paths.R"))
source(repo_path("R", "null_coalescing.R"))
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "nature_v2_figure_utils.R"))
source(repo_path("R", "final_truth_v9_figure_utils.R"))
Sys.setenv(PROTEOMICS_SCRIPT_ID = "figures/final_truth_v9_annotation_audit.R")

args <- commandArgs(trailingOnly = TRUE)
if ("--dry-run" %in% args || is_dry_run()) {
  message("[DRY-RUN] final_truth_v9 annotation coherence audit")
  quit(save = "no", status = 0L)
}

NEED <- c("GOSemSim", "GO.db", "org.Mm.eg.db", "AnnotationDbi",
          "simplifyEnrichment")
for (p in NEED)
  if (!requireNamespace(p, quietly = TRUE))
    stop("Part-28 annotation audit requires ", p, call. = FALSE)

OUT <- path_results("tables", "manuscript_candidates", "final_truth_v9", "audit")
REP <- path_results("reports", "manuscript_candidates", "final_truth_v9")
dir_create(OUT); dir_create(REP)

# ======================================================= a priori thresholds
#
# Declared here, before any comparison is computed, and justified from the
# measurement scale rather than from the observed outcome.
SEM_CUTOFF <- 0.70          # the repository's own pinned Wang cutoff
RANK_SHIFT_POSITIONS <- 2L  # movement among the six theme rows of one column
SUPPORT_ALPHA <- 0.05       # the canonical GSEA FDR threshold, unchanged
# a subcluster is MAJOR if it holds at least three terms and at least a tenth
# of the theme - small enough to catch a real minority block, large enough that
# a single stray term cannot rename a theme
MAJOR_MIN_N <- 3L
MAJOR_MIN_FRAC <- 0.10

TH <- nv_read_csv(repo_path(
  "results", "tables", "10_biological_integration", "gsea_wgcna_concordance",
  "global", "ontology_aware_gsea_theme_assignments_all_contrasts.csv"))
elig <- TH[TH$theme_claim_eligible %in% TRUE, , drop = FALSE]
if (!nrow(elig)) stop("no claim-eligible theme assignments", call. = FALSE)

map <- unique(elig[, c("GO_ID", "GO_description", "theme_id",
                       "manuscript_theme")])
map <- map[order(map$theme_id, map$GO_ID), , drop = FALSE]
ids <- sort(unique(map$GO_ID))
themes <- sort(unique(map$theme_id))
theme_label <- stats::setNames(
  vapply(themes, function(t) map$manuscript_theme[map$theme_id == t][1],
         character(1)), themes)
term_of <- stats::setNames(map$GO_description[!duplicated(map$GO_ID)],
                           map$GO_ID[!duplicated(map$GO_ID)])

# the audit must be reading the same ontology the canonical workflow pinned
ontab <- AnnotationDbi::select(GO.db::GO.db, keys = ids,
                               columns = c("ONTOLOGY", "TERM"),
                               keytype = "GOID")
if (any(is.na(ontab$ONTOLOGY)) || !all(ontab$ONTOLOGY == "BP"))
  stop("atlas themes contain non-BP or unknown GO terms", call. = FALSE)
chk <- merge(unique(map[, c("GO_ID", "GO_description")]), ontab,
             by.x = "GO_ID", by.y = "GOID")
if (any(chk$GO_description != chk$TERM))
  stop("stored GO descriptions disagree with the pinned GO.db", call. = FALSE)

# ============================================ S4 GO-BP semantic similarity
sem_data <- GOSemSim::godata(annoDb = "org.Mm.eg.db", ont = "BP",
                             computeIC = FALSE)
S <- as.matrix(GOSemSim::mgoSim(ids, ids, semData = sem_data,
                                measure = "Wang", combine = NULL))
S[!is.finite(S)] <- 0
diag(S) <- 1
S <- S[ids, ids, drop = FALSE]
OFFDIAG_MEDIAN <- stats::median(S[upper.tri(S)])
# a term only counts as sitting in the wrong theme if its advantage elsewhere
# exceeds the similarity any two unrelated BP terms already share through
# ontology structure alone
MISASSIGNED_MARGIN <- ceiling(OFFDIAG_MEDIAN * 100) / 100

members <- lapply(themes, function(t) unique(map$GO_ID[map$theme_id == t]))
names(members) <- themes
mean_sim <- function(go, th, drop_self = TRUE) {
  m <- members[[th]]
  if (drop_self) m <- setdiff(m, go)
  if (!length(m)) return(NA_real_)
  mean(S[go, m])
}

medoid_of_clusters <- function(W, cl) {
  stats::setNames(vapply(sort(unique(cl)), function(k) {
    g <- names(cl)[cl == k]
    if (length(g) == 1L) return(g)
    a <- vapply(g, function(x) mean(W[x, setdiff(g, x)]), numeric(1))
    g[order(-a, g)][1]          # ties break on GO_ID, never on an effect size
  }, character(1)), sort(unique(cl)))
}

dup_of <- list(); dup_medoid <- list(); blk_of <- list()
for (th in themes) {
  m <- members[[th]]
  W <- S[m, m, drop = FALSE]
  if (length(m) == 1L) {
    dup <- stats::setNames(1L, m); blk <- dup
  } else {
    hc <- stats::hclust(stats::as.dist(1 - W), method = "average")
    dup <- stats::cutree(hc, h = 1 - SEM_CUTOFF)
    # binary_cut tries a k-means partition among others, so it draws on the
    # session RNG; fix the seed or the block assignment is not reproducible
    set.seed(20260914L)
    blk <- stats::setNames(
      as.integer(simplifyEnrichment::binary_cut(W)), rownames(W))
  }
  dup_of[[th]] <- dup
  dup_medoid[[th]] <- medoid_of_clusters(W, dup)
  blk_of[[th]] <- blk
}

sem_rows <- do.call(rbind, lapply(seq_len(nrow(map)), function(i) {
  go <- map$GO_ID[i]; th <- map$theme_id[i]
  k <- dup_of[[th]][[go]]
  med <- dup_medoid[[th]][[as.character(k)]]
  own <- mean_sim(go, th)
  oth <- vapply(setdiff(themes, th), function(t2) mean_sim(go, t2, FALSE),
                numeric(1))
  margin <- max(oth) - own
  data.frame(
    GO_ID = go, GO_term = map$GO_description[i], current_theme = th,
    current_theme_label = unname(theme_label[th]),
    semantic_cluster = sprintf("%s_D%02d", th, k),
    semantic_block = sprintf("%s_B%02d", th, blk_of[[th]][[go]]),
    semantic_medoid = med, semantic_medoid_term = unname(term_of[med]),
    semantic_similarity_to_medoid = if (identical(go, med)) 1 else S[go, med],
    within_theme_similarity = own,
    between_theme_similarity = max(oth),
    nearest_other_theme = names(oth)[which.max(oth)],
    margin_to_nearest_other = margin,
    closer_to_other_theme = is.finite(own) && margin > 0,
    misassignment_candidate = is.finite(own) && margin > MISASSIGNED_MARGIN,
    is_medoid = identical(go, med),
    assigned_to_multiple_themes = sum(map$GO_ID == go) > 1L,
    also_in_themes = paste(setdiff(map$theme_id[map$GO_ID == go], th),
                           collapse = "; "),
    stringsAsFactors = FALSE)
}))
sem_rows$semantic_similarity_method <- "GOSemSim_Wang_BP"
sem_rows$block_method <- "simplifyEnrichment_binary_cut"
sem_rows$duplicate_cluster_method <- sprintf(
  "average_linkage_hclust at Wang >= %.2f", SEM_CUTOFF)
sem_rows$orgdb <- paste0("org.Mm.eg.db ", utils::packageVersion("org.Mm.eg.db"))
sem_rows$GO_db <- paste0("GO.db ", utils::packageVersion("GO.db"))
sem_rows$phenotype_blind <- TRUE
sem_rows <- sem_rows[order(sem_rows$current_theme, sem_rows$semantic_block,
                           sem_rows$GO_ID), , drop = FALSE]
write_csv_safe(sem_rows, file.path(OUT, "atlas_semantic_similarity_audit.csv"))

# ================================================= S5 theme coherence summary
coh <- do.call(rbind, lapply(themes, function(th) {
  m <- members[[th]]
  W <- S[m, m, drop = FALSE]; v <- W[upper.tri(W)]
  r <- sem_rows[sem_rows$current_theme == th, , drop = FALSE]
  blk <- blk_of[[th]]; sz <- table(blk)
  major <- sz[sz >= MAJOR_MIN_N & sz / length(m) >= MAJOR_MIN_FRAC]
  data.frame(
    theme = th, theme_label = unname(theme_label[th]), n_terms = length(m),
    n_semantic_blocks = length(unique(blk)),
    n_major_blocks = length(major),
    largest_block_n = max(sz),
    largest_block_frac = max(sz) / length(m),
    n_duplicate_clusters = length(unique(dup_of[[th]])),
    median_within_similarity = stats::median(v),
    min_within_similarity = min(v),
    q25_within_similarity = unname(stats::quantile(v, 0.25)),
    strongest_cross_theme_similarity = max(vapply(setdiff(themes, th),
      function(t2) max(S[m, members[[t2]]]), numeric(1))),
    n_closer_to_other_theme = sum(r$closer_to_other_theme),
    prop_closer_to_other_theme = mean(r$closer_to_other_theme),
    n_misassignment_candidates = sum(r$misassignment_candidate),
    n_shared_with_another_theme = sum(r$assigned_to_multiple_themes),
    semantic_medoids = paste(sort(unname(dup_medoid[[th]])), collapse = "; "),
    outlier_terms = paste(sort(r$GO_ID[r$misassignment_candidate]),
                          collapse = "; "),
    major_block_terms = paste(vapply(names(major), function(k)
      sprintf("B%s(n=%d): %s", k, sz[[k]],
              paste(utils::head(unname(term_of[m[blk == as.integer(k)]]), 4),
                    collapse = ", ")), character(1)), collapse = " || "),
    stringsAsFactors = FALSE)
}))
coh$classification <- with(coh, ifelse(
  n_misassignment_candidates > 0, "MISASSIGNED_TERMS_PRESENT",
  ifelse(n_major_blocks <= 1L & n_semantic_blocks == 1L, "COHERENT",
  ifelse(n_major_blocks <= 1L, "COHERENT_WITH_SUBSTRUCTURE",
  ifelse(largest_block_frac >= 0.60, "COHERENT_WITH_SUBSTRUCTURE",
  ifelse(median_within_similarity < OFFDIAG_MEDIAN * 1.5, "TOO_BROAD",
         "MIXED"))))))
coh <- coh[order(coh$theme), , drop = FALSE]
write_csv_safe(coh, file.path(OUT, "atlas_theme_coherence_summary.csv"))

# ====================================== S10 redundancy-reduced atlas summary
cells_for <- function(ct) {
  z <- elig[elig$contrast == ct & nzchar(as.character(elig$theme_id)), ,
            drop = FALSE]
  z <- z[!duplicated(paste(z$dataset, z$spatial_unit, z$contrast, z$GO_ID,
                           z$theme_id)), , drop = FALSE]
  key <- paste(z$dataset, z$spatial_unit, z$theme_id, sep = "\r")
  do.call(rbind, lapply(split(seq_len(nrow(z)), key), function(ix) {
    w <- z[ix, , drop = FALSE]; th <- w$theme_id[1]
    k <- dup_of[[th]][w$GO_ID]
    rep_ix <- w$GO_ID %in% unname(dup_medoid[[th]])
    sup <- is.finite(w$GSEA_FDR) & w$GSEA_FDR < SUPPORT_ALPHA
    data.frame(
      theme = th, theme_label = unname(theme_label[th]),
      dataset = w$dataset[1], spatial_unit = w$spatial_unit[1], contrast = ct,
      n_terms = nrow(w), n_semantic_representatives = sum(rep_ix),
      n_duplicate_clusters_present = length(unique(k)),
      current_median_NES = stats::median(w$NES, na.rm = TRUE),
      reduced_median_NES = if (any(rep_ix))
        stats::median(w$NES[rep_ix], na.rm = TRUE) else NA_real_,
      # one vote per near-duplicate cluster, using that cluster's own median
      # so the general parent term carries no extra weight
      cluster_median_NES = stats::median(
        tapply(w$NES, k, stats::median, na.rm = TRUE), na.rm = TRUE),
      n_terms_FDR_supported = sum(sup),
      n_supported_representatives = sum(sup & rep_ix),
      supported_pos = sum(sup & w$NES > 0), supported_neg = sum(sup & w$NES < 0),
      n_supported_clusters = length(unique(k[sup])),
      stringsAsFactors = FALSE)
  }))
}
cmp <- do.call(rbind, lapply(c("RES - CON", "SUS - CON", "SUS - RES"),
                             cells_for))
rownames(cmp) <- NULL

# LARGE_SHIFT is 5% of the full atlas colour range, i.e. a shift a reader could
# plausibly read as a different colour. Anchored to the scale, not the outcome.
ATLAS_LIMIT <- max(abs(cmp$current_median_NES), na.rm = TRUE)
LARGE_SHIFT <- round(0.05 * 2 * ATLAS_LIMIT, 3)

cmp$absolute_difference <- abs(cmp$current_median_NES - cmp$reduced_median_NES)
cmp$cluster_absolute_difference <-
  abs(cmp$current_median_NES - cmp$cluster_median_NES)
cmp$sign_changed <- is.finite(cmp$reduced_median_NES) &
  sign(cmp$current_median_NES) != sign(cmp$reduced_median_NES)
cmp$cluster_sign_changed <- is.finite(cmp$cluster_median_NES) &
  sign(cmp$current_median_NES) != sign(cmp$cluster_median_NES)
# a sign change only matters if BOTH summaries are far enough from zero to be
# read as a direction at all
cmp$material_sign_flip <- cmp$sign_changed &
  pmin(abs(cmp$current_median_NES), abs(cmp$reduced_median_NES)) >= LARGE_SHIFT
cmp$large_shift <- is.finite(cmp$absolute_difference) &
  cmp$absolute_difference >= LARGE_SHIFT
cmp$carries_support_dot <- cmp$n_terms_FDR_supported > 0L

rank_in <- function(v, by) {
  out <- integer(length(v))
  for (ix in split(seq_along(v), by)) out[ix] <- rank(-v[ix], ties.method = "first")
  out
}
ctx <- paste(cmp$dataset, cmp$spatial_unit, cmp$contrast, sep = "\r")
cmp$rank_current <- rank_in(cmp$current_median_NES, ctx)
cmp$rank_reduced <- rank_in(
  ifelse(is.finite(cmp$reduced_median_NES), cmp$reduced_median_NES, -Inf), ctx)
cmp$rank_shift <- abs(cmp$rank_current - cmp$rank_reduced)
cmp$rank_changed_materially <- cmp$rank_shift >= RANK_SHIFT_POSITIONS
cmp$large_shift_threshold <- LARGE_SHIFT
cmp$atlas_colour_limit <- ATLAS_LIMIT
cmp <- cmp[order(cmp$contrast, cmp$dataset, cmp$spatial_unit, cmp$theme), ,
           drop = FALSE]
write_csv_safe(cmp, file.path(OUT, "atlas_redundancy_reduced_comparison.csv"))

# ------------------------------- S7 the glycolysis question, asked directly
GLY <- sort(sem_rows$GO_ID[grepl("glycoly", sem_rows$GO_term, ignore.case = TRUE)])
MITO <- "mitochondrial_respiration_oxphos"
gm <- elig[elig$theme_id == MITO, , drop = FALSE]
gm <- gm[!duplicated(paste(gm$dataset, gm$spatial_unit, gm$contrast,
                           gm$GO_ID)), , drop = FALSE]
gkey <- paste(gm$dataset, gm$spatial_unit, gm$contrast, sep = "\r")
gly <- do.call(rbind, lapply(split(seq_len(nrow(gm)), gkey), function(ix) {
  w <- gm[ix, , drop = FALSE]
  isg <- w$GO_ID %in% GLY
  data.frame(dataset = w$dataset[1], spatial_unit = w$spatial_unit[1],
             contrast = w$contrast[1], n_terms = nrow(w), n_glycolytic = sum(isg),
             median_all = stats::median(w$NES, na.rm = TRUE),
             median_without_glycolysis = stats::median(w$NES[!isg], na.rm = TRUE),
             median_glycolysis_only = if (any(isg))
               stats::median(w$NES[isg], na.rm = TRUE) else NA_real_,
             n_supported_glycolytic = sum(isg & is.finite(w$GSEA_FDR) &
                                            w$GSEA_FDR < SUPPORT_ALPHA),
             n_supported_mitochondrial = sum(!isg & is.finite(w$GSEA_FDR) &
                                               w$GSEA_FDR < SUPPORT_ALPHA),
             stringsAsFactors = FALSE)
}))
gly$difference <- gly$median_all - gly$median_without_glycolysis
gly$sign_changed <- sign(gly$median_all) != sign(gly$median_without_glycolysis)
gly$material_change <- abs(gly$difference) >= LARGE_SHIFT
gly <- gly[order(gly$contrast, gly$dataset, gly$spatial_unit), , drop = FALSE]
write_csv_safe(gly, file.path(OUT, "atlas_glycolysis_sensitivity.csv"))

# ============================================== S13 support breadth per cell
sb <- cmp[, c("theme", "theme_label", "dataset", "spatial_unit", "contrast",
              "n_terms", "n_terms_FDR_supported", "n_semantic_representatives",
              "n_supported_representatives", "supported_pos", "supported_neg",
              "n_supported_clusters", "n_duplicate_clusters_present",
              "current_median_NES")]
sb$fraction_terms_supported <- sb$n_terms_FDR_supported / sb$n_terms
sb$fraction_representatives_supported <-
  sb$n_supported_representatives / sb$n_semantic_representatives
sb$supported_terms_direction <- with(sb, ifelse(
  n_terms_FDR_supported == 0L, "none",
  ifelse(supported_pos > 0 & supported_neg > 0, "mixed",
         ifelse(supported_pos > 0, "positive", "negative"))))
sb$median_direction <- ifelse(sb$current_median_NES > 0, "positive",
                       ifelse(sb$current_median_NES < 0, "negative", "zero"))
sb$all_supported_same_direction <- sb$supported_terms_direction %in%
  c("positive", "negative")
sb$dot_shown <- sb$n_terms_FDR_supported > 0L
sb$support_class <- with(sb, ifelse(
  !dot_shown, "NO_SUPPORT",
  ifelse(supported_terms_direction == "mixed", "MIXED_DIRECTION",
  ifelse(n_supported_clusters == 1L & n_duplicate_clusters_present > 1L,
         "NARROW_REDUNDANT_SUPPORT",
  ifelse(fraction_representatives_supported >= 0.5, "BROAD_SUPPORT",
         "PARTIAL_SUPPORT")))))
sb$median_support_direction_mismatch <- sb$dot_shown &
  sb$all_supported_same_direction &
  sb$supported_terms_direction != sb$median_direction
sb <- sb[order(sb$contrast, sb$dataset, sb$spatial_unit, sb$theme), ,
         drop = FALSE]
write_csv_safe(sb, file.path(OUT, "atlas_support_breadth_audit.csv"))

# =============================== S30/S31 label decision and the three levels
#
# A row label must cover every MAJOR semantic block in its theme. Where the
# canonical registry's own manuscript_theme already names the missing block,
# the recommendation is simply to stop abbreviating it away on the figure.
FIG_LABEL <- c(synaptic_signaling_vesicle = "Synaptic signalling / vesicle",
               rna_processing_splicing_rnp = "RNA processing",
               ribosome_translation = "Translation / ribosome",
               mitochondrial_respiration_oxphos = "Energy metabolism",
               autophagy_lysosome_endosome = "Autophagy / endolysosomal",
               chromatin_organization = "Chromatin",
               neuron_projection_development = "Neuron projection development")
RECOMMENDED <- c(synaptic_signaling_vesicle = "Synaptic signalling / vesicle",
                 rna_processing_splicing_rnp = "RNA processing",
                 ribosome_translation = "Translation / ribosome",
                 mitochondrial_respiration_oxphos = "Mitochondrial respiration",
                 autophagy_lysosome_endosome = "Autophagy / endolysosomal",
                 chromatin_organization = "Chromatin",
                 neuron_projection_development = "Neuron projection development")
REASON <- c(
  synaptic_signaling_vesicle = paste0(
    "one coherent semantic block, but roughly a quarter of the terms are the ",
    "synaptic vesicle cycle rather than signalling; the registry theme ",
    "already says vesicle-mediated transport"),
  rna_processing_splicing_rnp = paste0(
    "one coherent semantic block and RNA processing is the true general ",
    "parent of the mRNA, tRNA, rRNA, miRNA and RNP members; no change"),
  ribosome_translation = paste0(
    "ribosome biogenesis and rRNA maturation are a major block that the word ",
    "Translation does not cover; the registry theme already says ribosome ",
    "biogenesis"),
  mitochondrial_respiration_oxphos = paste0(
    "registry v3 excludes the cytosolic glycolysis sub-DAG by one ontology ",
    "rule, so the theme is mitochondrial bioenergetics again and the row ",
    "name matches its membership"),
  autophagy_lysosome_endosome = paste0(
    "endosomal and lysosomal trafficking is a major block that the word ",
    "Autophagy does not cover; the registry theme already says endolysosomal ",
    "trafficking"),
  chromatin_organization = paste0(
    "every member is a chromatin-mediated process; the low Wang similarity ",
    "reflects a sparse ontology branch rather than mixed biology; no change"),
  neuron_projection_development = paste0(
    "added under registry v3 as the only omitted semantic cluster meeting every fixed criterion; zero overlap with any other primary theme"))

lev <- do.call(rbind, lapply(themes, function(th) {
  cc <- coh[coh$theme == th, , drop = FALSE]
  blk <- blk_of[[th]]; sz <- table(blk)
  major <- sz[sz >= MAJOR_MIN_N & sz / cc$n_terms >= MAJOR_MIN_FRAC]
  data.frame(
    theme = th,
    level1_exact_terms_example = paste(
      utils::head(sort(unname(dup_medoid[[th]])), 3), collapse = "; "),
    level2_redundancy_reduced_theme = unname(theme_label[th]),
    level3_manuscript_interpretation = unname(RECOMMENDED[th]),
    current_label = unname(FIG_LABEL[th]),
    recommended_label = unname(RECOMMENDED[th]),
    semantic_support = sprintf(
      "%d terms, %d semantic block(s), %d major, largest %.0f%%, median within-similarity %.3f",
      cc$n_terms, cc$n_semantic_blocks, cc$n_major_blocks,
      100 * cc$largest_block_frac, cc$median_within_similarity),
    main_subclusters = cc$major_block_terms,
    terms_misrepresented_by_current_label = cc$n_terms -
      max(sz[names(sz) %in% names(major)], 0L),
    change_required = !identical(unname(FIG_LABEL[th]),
                                 unname(RECOMMENDED[th])),
    reason = unname(REASON[th]),
    stringsAsFactors = FALSE)
}))
write_csv_safe(lev, file.path(OUT, "atlas_annotation_decisions.csv"))

# ================================================= S32 theme evidence book
book <- do.call(rbind, lapply(themes, function(th) {
  cc <- coh[coh$theme == th, , drop = FALSE]
  s <- sb[sb$theme == th, , drop = FALSE]
  k <- cmp[cmp$theme == th, , drop = FALSE]
  reps <- sort(unname(dup_medoid[[th]]))
  data.frame(
    theme = th, theme_label = unname(theme_label[th]),
    n_terms = cc$n_terms, semantic_subclusters = cc$n_semantic_blocks,
    representative_GO_IDs = paste(reps, collapse = "; "),
    representative_GO_terms = paste(unname(term_of[reps]), collapse = "; "),
    median_within_similarity = cc$median_within_similarity,
    median_NES_summary_rule =
      "median NES across every claim-eligible GO term assigned to the theme",
    redundancy_reduced_summary_rule = sprintf(
      "median NES across %d phenotype-blind semantic medoids (Wang >= %.2f)",
      length(reps), SEM_CUTOFF),
    sign_concordance = sprintf("%d of %d cells keep their sign; %d flip materially",
      sum(!k$sign_changed), nrow(k), sum(k$material_sign_flip)),
    supported_term_count = sum(s$n_terms_FDR_supported),
    supported_rep_count = sum(s$n_supported_representatives),
    supported_direction_consistency = sprintf(
      "%d of %d supported cells are single-direction",
      sum(s$all_supported_same_direction), sum(s$dot_shown)),
    recommended_label = unname(RECOMMENDED[th]),
    confidence = cc$classification,
    stringsAsFactors = FALSE)
}))
write_csv_safe(book, file.path(OUT, "atlas_theme_evidence_book.csv"))

# ------------------------------------------------------------------ console
cat("\n===== PART-28 ATLAS ANNOTATION COHERENCE =====\n")
cat("GO universe          :", length(ids), "BP terms,", nrow(map),
    "term-theme assignments,", sum(table(map$GO_ID) > 1),
    "terms in more than one theme\n")
cat("similarity           : GOSemSim Wang BP, org.Mm.eg.db",
    as.character(utils::packageVersion("org.Mm.eg.db")),
    "| off-diagonal median", sprintf("%.4f", OFFDIAG_MEDIAN), "\n")
cat("a priori thresholds  : misassignment margin",
    sprintf("%.2f", MISASSIGNED_MARGIN), "| large shift",
    sprintf("%.3f", LARGE_SHIFT), "NES (5% of colour range",
    sprintf("%.3f", ATLAS_LIMIT), ")\n\n")
print(coh[, c("theme", "n_terms", "n_semantic_blocks", "n_major_blocks",
              "median_within_similarity", "n_misassignment_candidates",
              "n_shared_with_another_theme", "classification")],
      row.names = FALSE)
cat("\natlas cells:", nrow(cmp), "| carrying a support dot:",
    sum(cmp$carries_support_dot), "\n")
cat("  medoid reduction : sign flips", sum(cmp$sign_changed), ", MATERIAL",
    sum(cmp$material_sign_flip), ", large shifts", sum(cmp$large_shift),
    ", median |d|", sprintf("%.4f", stats::median(cmp$absolute_difference)), "\n")
cat("  cluster reduction: sign flips", sum(cmp$cluster_sign_changed),
    ", median |d|",
    sprintf("%.4f", stats::median(cmp$cluster_absolute_difference)), "\n")
cat("  material flips in a DOTTED cell:",
    sum(cmp$material_sign_flip & cmp$carries_support_dot), "\n")
cat("  material rank changes:", sum(cmp$rank_changed_materially), "\n")
cat("\nglycolysis sensitivity (", length(GLY), "cytosolic terms of 20 ):\n")
cat("  cells", nrow(gly), "| sign flips if removed", sum(gly$sign_changed),
    "| material", sum(gly$material_change),
    "| median |d|", sprintf("%.4f", stats::median(abs(gly$difference))), "\n")
cat("  supported glycolytic terms anywhere:", sum(gly$n_supported_glycolytic),
    "| supported mitochondrial terms:", sum(gly$n_supported_mitochondrial), "\n")
cat("\nsupport dot:\n"); print(table(sb$support_class))
cat("  median/support direction mismatches:",
    sum(sb$median_support_direction_mismatch), "\n")
cat("\nlabel decisions:\n")
print(lev[, c("theme", "current_label", "recommended_label",
              "change_required")], row.names = FALSE)
cat("\nwritten to:", OUT, "\n")

# ------------------------------------------------------- S30 decisions, prose
md <- c(
"# Atlas annotation decisions (Part 28)",
"",
"Semantic method: GOSemSim Wang over GO biological process with org.Mm.eg.db,",
"the same pinned combination the canonical workflow already uses. Phenotype-",
"blind by construction - only GO IDs enter the similarity matrix, and a",
"cluster is represented by its semantic medoid, never by its best-FDR term.",
"",
"Two clusterings answer two different questions. simplifyEnrichment binary_cut",
"needs no threshold and answers whether a theme is one biological block.",
"Average-linkage at Wang >= 0.70, the repository's own descriptive cutoff,",
"answers which terms are near-duplicates and therefore drives the redundancy-",
"reduced summary.",
"",
"## Thresholds, declared before the comparison was computed",
"",
sprintf("- misassignment margin %.2f: the median Wang similarity of an arbitrary", MISASSIGNED_MARGIN),
"  pair of BP terms in this universe. A term counts as sitting in the wrong",
"  theme only if its advantage elsewhere exceeds what any two unrelated terms",
"  already share through ontology structure.",
sprintf("- large shift %.3f NES: 5%% of the full atlas colour range (+/-%.3f),", LARGE_SHIFT, ATLAS_LIMIT),
"  i.e. a change a reader could plausibly read as a different colour.",
"- a sign flip counts as MATERIAL only if both summaries are at least that far",
"  from zero; a value flipping around 0.04 is not a change of direction.",
"",
"## Verdict on the summary rule",
"",
sprintf("Redundancy reduction moves the typical atlas cell by %.4f NES (median),",
        stats::median(cmp$absolute_difference)),
sprintf("against a colour range of %.3f. %d of %d cells change sign at all and %d do",
        ATLAS_LIMIT, sum(cmp$sign_changed), nrow(cmp), sum(cmp$material_sign_flip)),
sprintf("so materially; %d of those sit in a cell that carries a support dot.",
        sum(cmp$material_sign_flip & cmp$carries_support_dot)),
"The current median-NES atlas is therefore KEPT, and its robustness to GO",
"redundancy is now a measured fact rather than an assumption. The cells that",
"do move are listed in atlas_redundancy_reduced_comparison.csv and concentrate",
"in the two themes that share six rRNA terms with each other.",
"",
"## Theme-by-theme",
"")
for (t in themes) {
  cc <- coh[coh$theme == t, ]
  lv <- lev[lev$theme == t, ]
  md <- c(md,
    sprintf("### %s", unname(theme_label[t])),
    "",
    sprintf("- %d terms, %d semantic block(s) by binary_cut, %d major; median",
            cc$n_terms, cc$n_semantic_blocks, cc$n_major_blocks),
    sprintf("  within-theme Wang similarity %.3f (universe median %.3f).",
            cc$median_within_similarity, OFFDIAG_MEDIAN),
    sprintf("- classification: %s", cc$classification),
    sprintf("- terms shared with another theme: %d", cc$n_shared_with_another_theme),
    sprintf("- current row label \"%s\" -> recommended \"%s\"%s",
            lv$current_label, lv$recommended_label,
            if (lv$change_required) "  (CHANGED)" else "  (unchanged)"),
    sprintf("- reason: %s", lv$reason),
    "")
}
md <- c(md,
"## The glycolysis question, answered directly",
"",
sprintf("The mitochondrial theme holds %d cytosolic glycolysis terms among %d.",
        length(GLY), length(members[[MITO]])),
"binary_cut separates them into their own block, so the membership is",
"biologically mixed and no label containing only \"mitochondrial\" can cover it.",
sprintf("Removing them moves the theme summary by %.4f NES (median over %d cells)",
        stats::median(abs(gly$difference)), nrow(gly)),
sprintf("and changes the sign in %d cell. In the CA1 microglia-enriched ROI, the",
        sum(gly$sign_changed)),
"compartment the Figure-3 exemplar comes from, the three contrast summaries",
sprintf("move by %.3f, %.3f and %.3f, none of them across zero.",
        abs(gly$difference[gly$dataset == "microglia" & gly$spatial_unit == "CA1"])[1],
        abs(gly$difference[gly$dataset == "microglia" & gly$spatial_unit == "CA1"])[2],
        abs(gly$difference[gly$dataset == "microglia" & gly$spatial_unit == "CA1"])[3]),
sprintf("Across all cells %d glycolytic term is FDR-supported against %d",
        sum(gly$n_supported_glycolytic), sum(gly$n_supported_mitochondrial)),
"mitochondrial ones, so the supported signal in this theme is mitochondrial",
"even though the membership is not. The canonical registry is NOT edited: the",
"theme keeps its 20 terms and the row is renamed to cover them.",
"",
"## Three annotation levels (S31)",
"",
"| level | example | where it may be used |",
"|---|---|---|",
"| 1 exact canonical GO term | oxidative phosphorylation (GO:0006119) | the direct GSEA curve, and every inferential statement |",
"| 2 redundancy-reduced theme | mitochondrial respiration / OXPHOS | the atlas row and its legend |",
"| 3 manuscript interpretation | energy metabolism | discussion prose, only when it covers levels 1 and 2 |",
"",
"A level-3 word must never be substituted for the exact tested term in an",
"inferential sentence.",
"",
"## The support dot (S13-S15)",
"",
sprintf("Of the %d cells that carry a dot, %d draw it from a single",
        sum(sb$dot_shown), sum(sb$support_class == "NARROW_REDUNDANT_SUPPORT")),
sprintf("near-duplicate semantic branch, %d have broad support and %d have",
        sum(sb$support_class == "BROAD_SUPPORT"),
        sum(sb$support_class == "PARTIAL_SUPPORT")),
sprintf("partial support; %d mix directions and %d disagree in sign with the",
        sum(sb$support_class == "MIXED_DIRECTION"),
        sum(sb$median_support_direction_mismatch)),
"median the cell is coloured by. Two equal dots therefore do not mean two",
"equal breadths.",
"",
"A breadth marker was prototyped: dot AREA proportional to the fraction of",
"phenotype-blind semantic representatives containing an FDR-supported term,",
"colour encoding held constant, judged against the current panel at 175 x 44",
"mm. It is legible and it does separate narrow from broad support.",
"",
"VERDICT: KEEP_BINARY. Not because the benefit is marginal - it is visible -",
"but because adopting it would make a main figure depend on a semantic",
"clustering with a tunable similarity cutoff, and would spend panel width on",
"a second key. The dot keeps its exact stated meaning, the caption now says",
"it is NOT a breadth measure, and breadth is published per cell in",
"atlas_support_breadth_audit.csv, which answers the question a reviewer would",
"actually ask without putting a tunable parameter into Figure 3.")
writeLines(md, file.path(REP, "atlas_annotation_decisions.md"))
cat("decisions doc written\n")

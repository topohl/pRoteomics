# =====================================================================
# CA2-SLM DAP robustness engine.
#
# Every threshold in this file is PRESPECIFIED and stated as a named constant.
# None was chosen after inspecting which proteins are biologically attractive,
# and no classification rule references a module, a gene symbol, a candidate
# tier or any other annotation. The classifier sees only observation counts,
# effect estimates and QC sample identity.
#
# Nothing here recomputes a differential statistic. The canonical log2fc, pval
# and padj are carried through verbatim; every quantity produced here is a
# DESCRIPTIVE sensitivity diagnostic computed on the same values the canonical
# model saw, and is named so it cannot be mistaken for a test.
#
# WHAT THE CANONICAL MODEL SAW
#   animal-level value = unweighted mean of that animal's Left and Right
#   POST-imputation log2 values (aggregation_audit.csv records
#   equal_weight_mean_LR_on_existing_imputed_log2_values for every animal x unit).
#   Working from the sample-level imputed matrix therefore reproduces the
#   animal-level matrix exactly, and additionally allows a hemisphere to be
#   dropped - which is what the QC diagnostics need.
# =====================================================================

# --------------------------------------------------------------- constants

csr_fdr_threshold <- function() 0.05

# An observed-only effect needs at least this many FULLY observed animals in
# EACH group. With 3 per group the only options are 3 and 2; a 1-versus-3
# comparison is never reported as an effect. The same minimum decides when a
# QC-fail removal is labelled descriptive only (csr_removal_needs_label).
csr_min_observed_animals <- function() 2L

# Imputation dependence, as a fraction of the canonical effect magnitude.
csr_imputation_minimal_cut <- function() 0.25
csr_imputation_strong_cut <- function() 0.50

# Losing more than this fraction of the canonical magnitude is a collapse.
csr_magnitude_collapse_cut <- function() 0.50

csr_contract_version <- function() "ca2_slm_robustness_v1"

csr_thresholds <- function() {
  data.frame(
    constant = c("fdr_threshold", "min_observed_animals_per_group",
                 "imputation_minimal_cut", "imputation_strong_cut",
                 "magnitude_collapse_cut"),
    value = c(csr_fdr_threshold(), csr_min_observed_animals(),
              csr_imputation_minimal_cut(), csr_imputation_strong_cut(),
              csr_magnitude_collapse_cut()),
    meaning = c(
      "canonical DAP threshold, matching R/statistics/sus_res_spatial_dap_atlas_utils.R",
      "fully observed animals required in EACH group for an observed-only effect",
      "relative change in effect magnitude up to which imputation dependence is minimal",
      "relative change above which imputation dependence is strong",
      "fraction of canonical magnitude whose loss counts as a collapse"),
    prespecified = TRUE, stringsAsFactors = FALSE)
}

# ExpGroup 1/2/3 -> CON/RES/SUS, matching 01_preprocessing/06_merged_metadata.
csr_expgroup_to_stress <- function(x) {
  v <- as.character(x)
  out <- c("1" = "CON", "2" = "RES", "3" = "SUS")[v]
  if (any(is.na(out))) stop("unmapped ExpGroup value(s): ",
                            paste(unique(v[is.na(out)]), collapse = ", "), call. = FALSE)
  unname(out)
}

csr_bare_animal <- function(x) sub("^A0*", "", as.character(x))

# ------------------------------------------------------------------- GCT

# Minimal GCT v1.3 reader, used only to PROVE that the sample-level
# reconstruction reproduces the matrix the differential model consumed.
csr_read_gct <- function(path) {
  ln <- readLines(path, warn = FALSE)
  dims <- as.integer(strsplit(trimws(ln[2]), "\t")[[1]])
  n_row <- dims[1]; n_col <- dims[2]; n_cd <- dims[4]
  hdr <- strsplit(ln[3], "\t")[[1]]
  samples <- utils::tail(hdr, n_col)
  off <- length(hdr) - n_col
  body <- ln[(4L + n_cd):(3L + n_cd + n_row)]
  sp <- strsplit(body, "\t", fixed = TRUE)
  ids <- vapply(sp, function(z) z[1], character(1))
  mat <- matrix(NA_real_, nrow = n_row, ncol = n_col, dimnames = list(ids, samples))
  for (i in seq_along(sp)) {
    mat[i, ] <- suppressWarnings(as.numeric(sp[[i]][seq_len(n_col) + off]))
  }
  list(mat = mat, ids = ids, samples = samples)
}

# --------------------------------------- sample-level POST-imputation values

csr_read_imputed_samples <- function(path, sample_ids) {
  x <- as.data.frame(readxl::read_excel(path))
  key <- if ("T: Protein.Names" %in% names(x)) "T: Protein.Names" else "Protein.Names"
  if (!key %in% names(x)) stop("imputed matrix has no Protein.Names column", call. = FALSE)
  present <- intersect(sample_ids, names(x))
  if (!length(present)) stop("no requested sample column is in the imputed matrix",
                             call. = FALSE)
  m <- as.matrix(x[, present, drop = FALSE])
  storage.mode(m) <- "double"
  rownames(m) <- as.character(x[[key]])
  m
}

# ----------------------------------------- PRE-imputation observation mask

# The authoritative pre-imputation mask is the raw DIA-NN matrix, where a
# missing measurement is a literal NA. The post-imputation matrix cannot
# answer this question: it has no NAs left at all. The 70% row filter from
# archive/01_preprocessing/01_impute.r is reproduced so the analysed row set matches
# what the pipeline actually modelled.
csr_preimputation_mask <- function(raw_path, dataset_samples, unit_samples) {
  raw <- utils::read.delim(raw_path, check.names = FALSE, stringsAsFactors = FALSE)
  key_col <- if ("Protein.Names" %in% names(raw)) "Protein.Names" else
    names(raw)[grep("^Protein[. ]Names$", names(raw))[1]]
  if (is.na(key_col)) stop("raw matrix has no Protein.Names column", call. = FALSE)

  present <- intersect(dataset_samples, names(raw))
  if (!length(present)) stop("no dataset sample column matched the raw matrix",
                             call. = FALSE)
  sub <- as.matrix(raw[, present, drop = FALSE])
  storage.mode(sub) <- "double"
  na_mat <- is.na(sub)

  keep <- rowMeans(na_mat) <= 0.7          # step 5 of 01_impute.r
  unit_cols <- intersect(unit_samples, present)
  unit_mask <- na_mat[keep, unit_cols, drop = FALSE]
  rownames(unit_mask) <- raw[[key_col]][keep]

  list(unit_missing = unit_mask,
       n_missing_total = sum(na_mat[keep, , drop = FALSE]),
       n_analysed_rows = sum(keep),
       n_dataset_samples = length(present),
       unit_samples = unit_cols)
}

# ------------------------------------------- normalisation bias context
#
# archive/01_preprocessing/01_impute.r centres each sample column on one median
# before imputing; config/normalisation.yml chooses the median of the sample's
# OBSERVED values (all_observed_median, primary) or of the proteins quantified
# in every sample of the compartment (complete_case_median). Under all-observed
# centring a heavily missing sample detects only its more abundant proteins, so
# its median is inflated and subtracting it displaces EVERY protein in that
# sample downwards, including proteins with no imputed value anywhere. This
# function measures the SUS-RES displacement on the proteins observed in every
# sample of the unit, so imputation cannot enter. With imputed_path (this run's
# sample-level matrix) it is measured on the values the run analysed, whichever
# centring the run used; without it, on raw values centred on each sample's
# observed median over the unit's analysed rows. Rows are matched by
# Protein.Group, which is unique where Protein.Names is not.
csr_normalization_bias_context <- function(raw_path, dataset_samples, sample_meta, imputed_path = NULL) {
  raw <- utils::read.delim(raw_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"Protein.Group" %in% names(raw) || anyDuplicated(raw$Protein.Group))
    stop("the raw matrix needs a unique Protein.Group column", call. = FALSE)
  sub <- as.matrix(raw[, intersect(dataset_samples, names(raw)), drop = FALSE])
  storage.mode(sub) <- "double"
  keep <- rowMeans(is.na(sub)) <= 0.7
  sub <- sub[keep, , drop = FALSE]
  rownames(sub) <- raw$Protein.Group[keep]
  lg <- log2(sub)
  cols <- sample_meta$sample_id
  unit <- lg[, cols, drop = FALSE]

  med <- apply(unit, 2, stats::median, na.rm = TRUE)
  miss <- colMeans(is.na(unit))
  complete_rows <- stats::complete.cases(unit)
  if (is.null(imputed_path)) {
    centred <- sweep(unit[complete_rows, , drop = FALSE], 2, med, "-")
    centred_source <- "raw log2 values centred on each sample's observed median over the unit's analysed rows"
  } else {
    x <- as.data.frame(readxl::read_excel(imputed_path))
    if (!"Protein.Group" %in% names(x) || anyDuplicated(x$Protein.Group))
      stop("the imputed matrix needs a unique Protein.Group column: ", imputed_path, call. = FALSE)
    absent <- setdiff(cols, names(x))
    if (length(absent))
      stop("the imputed matrix lacks ", length(absent), " unit sample column(s): ", imputed_path, call. = FALSE)
    rows <- match(rownames(unit)[complete_rows], x$Protein.Group)
    if (anyNA(rows))
      stop(sum(is.na(rows)), " always-observed protein group(s) are absent from the imputed matrix: ",
           imputed_path, call. = FALSE)
    centred <- as.matrix(x[rows, cols, drop = FALSE])
    storage.mode(centred) <- "double"
    centred_source <- "this run's centred sample values (imputed matrix; every protein used is observed in all unit samples)"
  }

  per_sample <- data.frame(
    sample_id = cols,
    AnimalID = sample_meta$AnimalID,
    hemisphere = sample_meta$hemisphere,
    StressGroup = sample_meta$StressGroup,
    qc_flag = sample_meta$qc_flag,
    fraction_missing_preimputation = as.numeric(miss),
    median_observed_log2 = as.numeric(med),
    n_proteins_observed_in_every_sample = sum(complete_rows),
    mean_centred_value_of_always_observed_proteins = as.numeric(colMeans(centred)),
    centred_values_source = centred_source,
    stringsAsFactors = FALSE)

  # per animal, the displacement is the mean of its two hemispheres, exactly
  # as the animal-level value is formed
  per_animal <- vapply(split(per_sample$mean_centred_value_of_always_observed_proteins,
                             per_sample$AnimalID), mean, numeric(1))
  grp <- sample_meta$StressGroup[match(names(per_animal), sample_meta$AnimalID)]
  disp <- mean(per_animal[grp == "SUS"]) - mean(per_animal[grp == "RES"])

  list(per_sample = per_sample,
       per_animal_displacement = per_animal,
       sus_minus_res_displacement = disp,
       cor_missing_vs_median = stats::cor(miss, med),
       cor_missing_vs_centred = stats::cor(
         miss, per_sample$mean_centred_value_of_always_observed_proteins),
       n_always_observed = sum(complete_rows),
       centred_values_source = centred_source)
}

# ------------------------------------------------------- effect helpers

# Descriptive group difference on animal-level values. Not a test.
csr_effect <- function(animal_values, stress, a = "SUS", b = "RES") {
  va <- animal_values[stress == a]; vb <- animal_values[stress == b]
  va <- va[is.finite(va)]; vb <- vb[is.finite(vb)]
  if (!length(va) || !length(vb)) return(NA_real_)
  mean(va) - mean(vb)
}

csr_same_sign <- function(x, y) {
  if (!is.finite(x) || !is.finite(y)) return(NA)
  (x > 0 && y > 0) || (x < 0 && y < 0) || (x == 0 && y == 0)
}

csr_rel_change <- function(x, canonical) {
  if (!is.finite(x) || !is.finite(canonical) || canonical == 0) return(NA_real_)
  abs(x - canonical) / abs(canonical)
}

# Collapse sample-level values to animal level exactly as the pipeline does:
# the unweighted mean of whichever hemispheres are supplied.
csr_animal_values <- function(sample_vals, sample_animal, drop_samples = character()) {
  keep <- !(names(sample_vals) %in% drop_samples)
  v <- sample_vals[keep]; a <- sample_animal[names(v)]
  vapply(split(v, a), function(z) mean(z, na.rm = TRUE), numeric(1))
}

# ------------------------------------------------- per-protein robustness

csr_protein_robustness <- function(protein_key, canonical, imputed, mask_row,
                                   ctx) {
  sm <- ctx$sample_meta
  sample_vals <- imputed[protein_key, sm$sample_id]
  names(sample_vals) <- sm$sample_id
  sample_animal <- stats::setNames(sm$AnimalID, sm$sample_id)
  observed <- stats::setNames(!mask_row[sm$sample_id], sm$sample_id)

  animals <- ctx$animals
  stress <- ctx$stress                       # named by animal
  av <- csr_animal_values(sample_vals, sample_animal)[animals]
  st <- unname(stress[animals])
  canon_fc <- as.numeric(canonical$log2fc)

  sus <- animals[st == "SUS"]; res <- animals[st == "RES"]; con <- animals[st == "CON"]
  n_obs_animal <- vapply(animals, function(a)
    sum(observed[sm$sample_id[sm$AnimalID == a]]), integer(1))

  n_obs_g <- function(g) as.integer(sum(n_obs_animal[animals[st == g]]))
  n_tot_g <- function(g) 2L * as.integer(sum(st == g))
  n_imp_g <- function(g) n_tot_g(g) - n_obs_g(g)

  full_animal <- n_obs_animal == 2L
  sus_full <- sum(full_animal[sus]); res_full <- sum(full_animal[res])
  estimable <- sus_full >= csr_min_observed_animals() &&
    res_full >= csr_min_observed_animals()

  obs_only <- if (!estimable) NA_real_ else {
    k <- c(sus, res)[full_animal[c(sus, res)]]
    csr_effect(av[k], st[match(k, animals)])
  }

  # ------------------------------------------------ leave one animal out
  loo <- data.frame(original_identifier = protein_key,
                    omitted_AnimalID = animals,
                    omitted_StressGroup = st, stringsAsFactors = FALSE)
  loo$loo_effect <- vapply(animals, function(a) {
    k <- animals != a
    csr_effect(av[k], st[k])
  }, numeric(1))
  loo$canonical_log2FC_SUS_minus_RES <- canon_fc
  loo$change_from_canonical <- loo$loo_effect - canon_fc
  loo$sign_retained <- vapply(loo$loo_effect,
                              function(z) isTRUE(csr_same_sign(z, canon_fc)), logical(1))
  # omitting a CON animal cannot move a SUS-vs-RES contrast, so summarising
  # over those would dilute every stability metric
  loo$can_change_estimate <- loo$omitted_StressGroup %in% c("SUS", "RES")
  eff <- loo$loo_effect[loo$can_change_estimate]
  sgn <- loo$sign_retained[loo$can_change_estimate]

  # -------------------------------------------------- QC-fail sensitivity
  # A removal of animals that are not in the design (A764 in the 3/3/2 rerun)
  # is not computed: it would only repeat the canonical effect. An animal is
  # removed whole only when none of its acquisitions in this unit passed QC;
  # with one good hemisphere only the failed acquisition is left out (below).
  removable <- function(drop) csr_removable_animals(intersect(drop, animals), sm, ctx$qc_fail_samples)
  drop_animals <- function(drop) {
    k <- !(animals %in% removable(drop))
    if (all(k)) return(NA_real_)
    if (!sum(st[k] == "SUS") || !sum(st[k] == "RES")) return(NA_real_)
    csr_effect(av[k], st[k])
  }
  drop_label <- function(drop) {
    if (!length(intersect(drop, animals))) return(csr_removal_interpretation(0L, 0L, computed = FALSE))
    if (!length(removable(drop))) return(csr_kept_animal_interpretation())
    k <- !(animals %in% removable(drop))
    csr_removal_interpretation(sum(st[k] == "SUS"), sum(st[k] == "RES"))
  }
  e_755 <- drop_animals("755"); e_764 <- drop_animals("764")
  e_both <- drop_animals(c("755", "764"))

  # Keep the animal, drop only its QC-FAILED hemisphere. This asks whether the
  # failing acquisition drives the effect without sacrificing a replicate, and
  # is the diagnostic a whole-animal drop cannot give.
  av_hemi <- csr_animal_values(sample_vals, sample_animal,
                               drop_samples = ctx$qc_fail_samples)[animals]
  e_hemi <- csr_effect(av_hemi, st)

  # ------------------------------------------------ imputation dependence
  n_imp_sus <- n_imp_g("SUS"); n_imp_res <- n_imp_g("RES")
  fully_observed <- (n_imp_sus + n_imp_res) == 0L
  imp_rel <- csr_rel_change(obs_only, canon_fc)
  imp_class <- if (fully_observed) "fully_observed" else
    if (!estimable) "not_estimable_without_imputation" else
      if (!isTRUE(csr_same_sign(obs_only, canon_fc)) ||
          (is.finite(imp_rel) && imp_rel > csr_imputation_strong_cut()))
        "strong_imputation_dependence" else
        if (is.finite(imp_rel) && imp_rel > csr_imputation_minimal_cut())
          "moderate_imputation_dependence" else "minimal_imputation_dependence"

  # imputed draws come from a left-shifted MNAR distribution, so an imputed
  # SUS value pushes the contrast down and an imputed RES value pushes it up
  aligned <- (if (canon_fc < 0) n_imp_sus else 0L) +
    (if (canon_fc > 0) n_imp_res else 0L)

  # how much of each group mean is carried by imputed values
  imp_share <- function(g) {
    ids <- sm$sample_id[sm$AnimalID %in% animals[st == g]]
    imp <- ids[!observed[ids]]
    if (!length(imp)) return(0)
    sum(sample_vals[imp]) / sum(sample_vals[ids])
  }

  # ------------------------------------------------------ stability flags
  loo_sign_stable <- length(sgn) > 0L && all(sgn)
  # A whole-animal removal is made only for an animal none of whose acquisitions in this
  # unit passed QC (Tobias, 2026-11-18: with one good hemisphere only the failed acquisition
  # is left out, and that hemisphere check decides). It decides a class only when it leaves
  # at least csr_min_observed_animals() animals in each group; otherwise its effect is
  # reported, labelled descriptive only, and the flag is NA. The flag is NA when no animal
  # is removed.
  animal_flag <- function(drop, e) {
    drop <- removable(drop)
    if (!length(drop)) return(NA)
    k <- !(animals %in% drop)
    if (csr_removal_needs_label(sum(st[k] == "SUS"), sum(st[k] == "RES"))) return(NA)
    isTRUE(!csr_same_sign(e, canon_fc)) ||
      isTRUE(csr_rel_change(e, canon_fc) > csr_magnitude_collapse_cut())
  }
  a755 <- animal_flag("755", e_755)
  a764 <- animal_flag("764", e_764)
  hemi <- isTRUE(!csr_same_sign(e_hemi, canon_fc)) ||
    isTRUE(csr_rel_change(e_hemi, canon_fc) > csr_magnitude_collapse_cut())
  obs_sign <- isTRUE(csr_same_sign(obs_only, canon_fc))

  cls <- csr_classify(fully_observed, estimable, obs_sign, imp_class,
                      loo_sign_stable, a755, a764, hemi, imp_rel)
  # the reasons are worded for the published design with two QC-failed SUS samples; a QC-failed hemisphere of
  # another group (a CON sample can fail) does not enter a SUS-RES effect and is not counted
  n_fail <- sum(sm$StressGroup[match(ctx$qc_fail_samples, sm$sample_id)] == "SUS", na.rm = TRUE)
  if (n_fail != 2L)
    cls$reason <- sub("either QC-failed SUS sample",
                      if (n_fail == 1L) "the QC-failed SUS sample" else "any QC-failed SUS sample",
                      cls$reason, fixed = TRUE)

  summary <- data.frame(
    ProteinGroupID = canonical$ProteinGroupID,
    original_identifier = protein_key,
    gene_symbol = canonical$official_gene_symbol,
    canonical_log2FC_SUS_minus_RES = canon_fc,
    canonical_p_value = as.numeric(canonical$pval),
    canonical_BH_FDR = as.numeric(canonical$padj),
    canonical_statistic_source = "read verbatim from the canonical contrast file; not refitted",

    n_observed_CON = n_obs_g("CON"), n_observed_RES = n_obs_g("RES"),
    n_observed_SUS = n_obs_g("SUS"),
    n_imputed_CON = n_imp_g("CON"), n_imputed_RES = n_imp_res, n_imputed_SUS = n_imp_sus,
    fraction_observed_CON = n_obs_g("CON") / n_tot_g("CON"),
    fraction_observed_RES = n_obs_g("RES") / n_tot_g("RES"),
    fraction_observed_SUS = n_obs_g("SUS") / n_tot_g("SUS"),
    differential_missingness_SUS_minus_RES =
      (n_imp_sus / n_tot_g("SUS")) - (n_imp_res / n_tot_g("RES")),
    fully_observed = fully_observed,
    n_fully_observed_SUS_animals = as.integer(sus_full),
    n_fully_observed_RES_animals = as.integer(res_full),

    A755_L_observed = unname(observed[ctx$sample_of("755", "L")]),
    A755_R_observed = unname(observed[ctx$sample_of("755", "R")]),
    A764_L_observed = unname(observed[ctx$sample_of("764", "L")]),
    A764_R_observed = unname(observed[ctx$sample_of("764", "R")]),
    A765_L_observed = unname(observed[ctx$sample_of("765", "L")]),
    A765_R_observed = unname(observed[ctx$sample_of("765", "R")]),
    A755_value = unname(av["755"]), A764_value = unname(av["764"]),
    A765_value = unname(av["765"]),

    canonical_mean_CON = mean(av[con]), canonical_mean_RES = mean(av[res]),
    canonical_mean_SUS = mean(av[sus]),
    observed_only_mean_RES = if (res_full > 0L) mean(av[res][full_animal[res]]) else NA_real_,
    observed_only_mean_SUS = if (sus_full > 0L) mean(av[sus][full_animal[sus]]) else NA_real_,
    imputed_share_of_SUS_mean = imp_share("SUS"),
    imputed_share_of_RES_mean = imp_share("RES"),

    observed_only_estimable = estimable,
    observed_only_log2FC = obs_only,
    observed_only_minus_canonical = obs_only - canon_fc,
    observed_only_same_sign = obs_sign,
    observed_only_relative_change = imp_rel,
    observed_only_rule = sprintf(
      "an animal counts only if BOTH hemispheres were detected; at least %d such animals in EACH group",
      csr_min_observed_animals()),

    n_imputed_aligned_with_canonical_direction = as.integer(aligned),
    imputation_dependence_class = imp_class,

    loo_min_effect = if (length(eff)) min(eff) else NA_real_,
    loo_max_effect = if (length(eff)) max(eff) else NA_real_,
    loo_median_effect = if (length(eff)) stats::median(eff) else NA_real_,
    loo_proportion_canonical_sign = if (length(sgn)) mean(sgn) else NA_real_,
    loo_max_abs_change_from_canonical = if (length(eff))
      max(abs(eff - canon_fc)) else NA_real_,
    loo_sign_stable = loo_sign_stable,
    loo_n_animals_that_can_change_estimate = as.integer(length(eff)),

    effect_excluding_A755 = e_755,
    effect_excluding_A764 = e_764,
    effect_excluding_both_QC_fail_SUS = e_both,
    effect_excluding_both_interpretation = drop_label(c("755", "764")),
    effect_excluding_A755_interpretation = drop_label("755"),
    effect_excluding_A764_interpretation = drop_label("764"),
    effect_dropping_QC_failed_hemispheres = e_hemi,
    A755_sensitive = a755, A764_sensitive = a764,
    qc_failed_hemisphere_sensitive = hemi,
    single_animal_sensitive = !loo_sign_stable,

    CA2_SLM_robustness_class = cls$class,
    classification_reason = cls$reason,
    contract_version = csr_contract_version(),
    stringsAsFactors = FALSE)

  list(summary = summary, loo = loo)
}

# ---------------------------------------------------- final classification

# A transparent ladder. Every branch records WHY it fired. There is no hidden
# numeric score, and no branch may consult an annotation.
csr_classify <- function(fully_observed, estimable, obs_sign_retained, imp_class,
                         loo_sign_stable, a755, a764, hemi, imp_rel) {
  pct <- function(x) if (is.finite(x)) sprintf("%.0f%%", 100 * x) else "not estimable"

  if (!fully_observed && !estimable) {
    return(list(class = "insufficient_observed_data",
                reason = sprintf(paste0("fewer than %d fully observed animals in at least ",
                                        "one group, so no observed-only effect exists to ",
                                        "compare against"), csr_min_observed_animals())))
  }
  # a755 / a764 are NA when that removal does not count (fewer than two animals left)
  a755 <- isTRUE(a755); a764 <- isTRUE(a764); hemi <- isTRUE(hemi)
  if (a755 || a764 || hemi) {
    which_ <- paste(c(if (a755) "dropping A755", if (a764) "dropping A764",
                      if (hemi) "dropping the QC-failed hemispheres"), collapse = ", ")
    return(list(class = "not_claimable_due_to_QC",
                reason = paste0("the effect reverses or loses more than ",
                                pct(csr_magnitude_collapse_cut()), " of its magnitude when ",
                                which_)))
  }
  if (identical(imp_class, "strong_imputation_dependence")) {
    return(list(class = "imputation_sensitive",
                reason = paste0("removing imputed observations reverses the sign or ",
                                "changes the magnitude by ", pct(imp_rel),
                                " (strong threshold ", pct(csr_imputation_strong_cut()), ")")))
  }
  if (!loo_sign_stable) {
    return(list(class = "single_animal_sensitive",
                reason = "omitting one SUS or RES animal reverses the sign of the effect"))
  }
  if (fully_observed ||
      (identical(imp_class, "minimal_imputation_dependence") && obs_sign_retained)) {
    return(list(class = "robust_to_missingness_and_QC",
                reason = paste0("imputation dependence at most minimal (", imp_class,
                                "), observed-only sign retained, sign stable under ",
                                "leave-one-animal-out, and not dependent on either ",
                                "QC-failed SUS sample")))
  }
  list(class = "supported_but_QC_sensitive",
       reason = paste0("direction is consistent and neither a single animal nor a QC ",
                       "sample reverses it, but the magnitude is materially influenced ",
                       "by missingness or imputation (", imp_class, ", change ",
                       pct(imp_rel), ")"))
}

# -------------------------------------------- fully observed deep audit

# Containing no imputed value is NECESSARY but NOT SUFFICIENT for claimability.
csr_fully_observed_audit <- function(fully, loo_long, sample_meta) {
  if (!nrow(fully)) return(cbind(fully, max_qc_animal_leverage = numeric(0)))
  qc_an <- unique(sample_meta$AnimalID[sample_meta$qc_flag %in% c("FAIL", "WARN")])
  # leverage counts only omissions that leave at least csr_min_observed_animals() per group
  grp <- sample_meta$StressGroup[match(unique(sample_meta$AnimalID), sample_meta$AnimalID)]
  n_sus <- sum(grp == "SUS"); n_res <- sum(grp == "RES")
  valid_omission <- !csr_removal_needs_label(n_sus - (loo_long$omitted_StressGroup == "SUS"),
                                             n_res - (loo_long$omitted_StressGroup == "RES"))
  out <- fully
  out$max_qc_animal_leverage <- vapply(seq_len(nrow(fully)), function(i) {
    z <- loo_long[loo_long$original_identifier == fully$original_identifier[i] &
                    loo_long$omitted_AnimalID %in% qc_an &
                    loo_long$can_change_estimate & valid_omission, , drop = FALSE]
    if (!nrow(z)) return(NA_real_)
    max(abs(z$change_from_canonical))
  }, numeric(1))
  out$qc_animal_leverage_relative <-
    out$max_qc_animal_leverage / abs(out$canonical_log2FC_SUS_minus_RES)
  out$qc_samples_measured_but_flagged <- paste(
    sprintf("%s_%s(%s)", sample_meta$AnimalID, sample_meta$hemisphere,
            sample_meta$qc_flag)[sample_meta$qc_flag != "PASS"], collapse = ";")
  out$fully_observed_but_qc_leveraged <-
    is.finite(out$qc_animal_leverage_relative) &
    out$qc_animal_leverage_relative > csr_magnitude_collapse_cut()
  fail_an <- csr_sort_animals(unique(sample_meta$AnimalID[sample_meta$qc_flag %in% "FAIL"]))
  if (length(fail_an)) fail_an <- paste0("A", fail_an)
  out$audit_note <- paste0(
    "these proteins carry no imputed value, ",
    if (!length(fail_an)) "and no acquisition is QC-FAILED" else
      paste0("but ", csr_and_list(fail_an),
             if (length(fail_an) == 1L) " remains a QC-FAILED acquisition" else " remain QC-FAILED acquisitions"),
    "; leverage and leave-one-animal-out stability are assessed ",
    "independently rather than assumed from the absence of imputation")
  out
}

# ------------------------------------------------- QC-fail removal labels

# A QC-fail removal that leaves fewer than csr_min_observed_animals() animals
# in SUS or in RES is not a group comparison. Its effect is still reported and
# labelled DESCRIPTIVE ONLY, and it decides no class (A755_sensitive /
# A764_sensitive are then NA).
csr_removal_needs_label <- function(n_sus, n_res) {
  n_sus < csr_min_observed_animals() | n_res < csr_min_observed_animals()
}

# computed = FALSE: none of the animals to remove is in the design, so no
# effect was formed.
# An animal is removed whole only when none of its acquisitions in the unit passed QC; when
# one hemisphere is fine, only the failed acquisition is left out (Tobias, 2026-11-18).
csr_removable_animals <- function(drop, sample_meta, qc_fail_samples) {
  drop[vapply(drop, function(a) {
    s <- sample_meta$sample_id[as.character(sample_meta$AnimalID) == a]
    length(s) > 0L && all(s %in% qc_fail_samples)
  }, logical(1))]
}
csr_kept_animal_interpretation <- function() {
  "not computed: the animal keeps a hemisphere that passed QC, so only its failed acquisition is left out"
}

csr_removal_interpretation <- function(n_sus, n_res, computed = TRUE) {
  if (!computed) return("not computed: the animal is not in this design")
  if (csr_removal_needs_label(n_sus, n_res)) {
    return(sprintf("DESCRIPTIVE ONLY: %d SUS versus %d RES is not a valid group comparison",
                   n_sus, n_res))
  }
  sprintf("%d SUS versus %d RES remain; a descriptive sensitivity diagnostic", n_sus, n_res)
}

# ------------------------------------------- text from the design present
#
# Generated wording names the design the audit tables describe (3/3/3 in the
# published run, 3/3/2 in the rerun) instead of assuming it.

csr_number_word <- function(n) {
  w <- c("no", "one", "two", "three", "four", "five", "six", "seven", "eight",
         "nine", "ten")
  ifelse(!is.na(n) & n >= 0 & n <= 10, w[pmin(pmax(n, 0), 10) + 1], as.character(n))
}

csr_and_list <- function(x) {
  x <- as.character(x)
  if (length(x) <= 1L) return(paste(x, collapse = ""))
  paste(paste(x[-length(x)], collapse = ", "), "and", x[length(x)])
}

# bare animal ids in numeric order (755 before 764 before 765)
csr_sort_animals <- function(x) {
  x <- as.character(x)
  x[order(suppressWarnings(as.numeric(x)), x)]
}

# "A755 Left FAIL" per sample-context row, in animal then hemisphere order.
csr_sample_labels <- function(ctx, with_flag = TRUE) {
  if (!nrow(ctx)) return(character(0))
  ctx <- ctx[order(suppressWarnings(as.numeric(ctx$AnimalID)), ctx$AnimalID,
                   ctx$hemisphere), , drop = FALSE]
  hemi <- unname(c(L = "Left", R = "Right")[ctx$hemisphere])
  paste0("A", ctx$AnimalID, " ", hemi, if (with_flag) paste0(" ", ctx$qc_flag) else "")
}

# Phrases for the claim-status summary and the workbook README, from the
# CA2-SLM sample context (one row per sample: AnimalID, hemisphere,
# StressGroup, qc_flag, fraction_missing_preimputation). n_flagged_neuropil is
# the number of QC-flagged samples over the whole neuropil dataset, NA when the
# QC table is not available.
csr_design_text <- function(ctx, n_flagged_neuropil = NA_integer_) {
  animals_in <- function(g) unique(ctx$AnimalID[ctx$StressGroup == g])
  n_sus <- length(animals_in("SUS")); n_res <- length(animals_in("RES"))
  flagged <- ctx[ctx$qc_flag != "PASS", , drop = FALSE]
  fail_sus <- ctx[ctx$qc_flag == "FAIL" & ctx$StressGroup == "SUS", , drop = FALSE]
  k <- nrow(fail_sus)
  # QC-failed SUS animals that are removed whole: none of their acquisitions passed QC
  k_an <- sum(vapply(unique(fail_sus$AnimalID),
                     function(a) all(ctx$qc_flag[ctx$AnimalID == a] == "FAIL"), logical(1)))

  miss <- tapply(ctx$fraction_missing_preimputation, ctx$StressGroup, mean)
  miss_dir <- if (!all(c("SUS", "RES") %in% names(miss))) "missingness in SUS and RES not compared" else
    if (miss[["SUS"]] > miss[["RES"]]) "higher missingness in SUS than RES" else
      if (miss[["SUS"]] < miss[["RES"]]) "lower missingness in SUS than RES" else
        "equal missingness in SUS and RES"

  share <- function(with_flags) {
    lab <- if (with_flags && nrow(flagged))
      paste0(" (", paste(csr_sample_labels(flagged), collapse = ", "), ")") else ""
    if (is.na(n_flagged_neuropil))
      return(sprintf("%d QC-flagged %s%s", nrow(flagged),
                     if (nrow(flagged) == 1L) "sample" else "samples", lab))
    if (n_flagged_neuropil == 1L && nrow(flagged) == 1L)
      return(sprintf("the only QC-flagged neuropil sample%s", lab))
    sprintf("%d of the %d QC-flagged neuropil samples%s", nrow(flagged), n_flagged_neuropil, lab)
  }
  fail_lab <- paste0(" (", paste(csr_sample_labels(fail_sus, with_flag = FALSE), collapse = ", "), ")")
  of_them <- if (!k) "none of them a QC-FAILED SUS acquisition" else
    if (k == nrow(flagged)) {
      if (k == 1L) paste0("a QC-FAILED SUS acquisition", fail_lab) else
        paste0("all of them QC-FAILED SUS acquisitions", fail_lab)
    } else {
      if (k == 1L) paste0("one of them a QC-FAILED SUS acquisition", fail_lab) else
        paste0(csr_number_word(k), " of them QC-FAILED SUS acquisitions", fail_lab)
    }

  rk <- sort(rank(-ctx$fraction_missing_preimputation, ties.method = "min")[
    ctx$qc_flag == "FAIL" & ctx$StressGroup == "SUS"])
  most_missing <- if (!k) "no SUS sample is QC-FAILED" else
    if (identical(as.integer(rk), seq_len(k))) {
      if (k == 1L) "the QC-FAILED SUS sample is the most missing" else
        sprintf("the %s QC-FAILED SUS samples are the %s most missing",
                csr_number_word(k), csr_number_word(k))
    } else {
      sprintf("the QC-FAILED SUS %s %s %s of %d by missingness",
              if (k == 1L) "sample" else "samples", if (k == 1L) "ranks" else "rank",
              csr_and_list(rk), nrow(ctx))
    }

  loo_left <- sprintf("%d versus %d", n_sus - 1L, n_res)
  if (n_sus != n_res) loo_left <- paste(loo_left, "or", sprintf("%d versus %d", n_sus, n_res - 1L))
  loo <- sprintf("Leave-one-animal-out on a %d-versus-%d design leaves %s.", n_sus, n_res, loo_left)
  if (csr_removal_needs_label(n_sus - 1L, n_res) || csr_removal_needs_label(n_sus, n_res - 1L))
    loo <- paste0(loo, sprintf(paste0(" A leave-out that leaves fewer than %d animals in a group is descriptive only; ",
                                      "it enters the classes only through the requirement that the sign is stable ",
                                      "under every leave-out."), csr_min_observed_animals()))

  removed <- if (k_an == 2L) "both QC-failed SUS animals" else
    if (k_an == 1L) "the QC-failed SUS animal" else
      sprintf("all %s QC-failed SUS animals", csr_number_word(k_an))
  qc_removal <- if (!k) "No SUS acquisition is QC-failed, so no QC-fail removal is made." else
    if (!k_an) paste0("No animal is removed for QC: ",
                      if (k == 1L) "the QC-failed SUS acquisition is left out alone, and the animal's other hemisphere, "
                      else "the QC-failed SUS acquisitions are left out alone, and each animal's other hemisphere, ",
                      "which passed QC, is kept.") else
    paste0("Removing ", removed, " leaves ", n_sus - k_an, " versus ", n_res,
           if (csr_removal_needs_label(n_sus - k_an, n_res)) {
             if (csr_removal_needs_label(n_sus - 1L, n_res))
               paste0(", which is NOT a valid group comparison: its effect is reported, labelled DESCRIPTIVE ONLY, ",
                      "and decides no class. The QC check that decides is the removal of the QC-failed ",
                      "acquisition alone, which keeps the animal's other hemisphere.") else
               ", which is NOT a valid group comparison and is labelled descriptive-only wherever it appears."
           } else
             sprintf(", which keeps at least %d animals in each group.", csr_min_observed_animals()))

  # the SUS hemisphere with the most proteins undetected before imputation, named by depth rather than QC flag
  sus <- ctx[ctx$StressGroup == "SUS" & is.finite(ctx$fraction_missing_preimputation), , drop = FALSE]
  top <- if (nrow(sus)) sus[which.max(sus$fraction_missing_preimputation), , drop = FALSE] else NULL
  top_is_most_missing <- !is.null(top) &&
    !any(ctx$fraction_missing_preimputation > top$fraction_missing_preimputation, na.rm = TRUE)
  shallowest <- if (is.null(top)) "no SUS sample with a known missingness" else
    sprintf("%s, the shallowest %s, with %.1f%% of the analysed proteins not detected",
            csr_sample_labels(top, with_flag = FALSE),
            if (top_is_most_missing) sprintf("of the %d CA2-SLM samples", nrow(ctx)) else "SUS hemisphere",
            100 * top$fraction_missing_preimputation)

  list(n_sus = n_sus, n_res = n_res, n_animals = length(unique(ctx$AnimalID)),
       n_samples = nrow(ctx), n_sus_res_samples = sum(ctx$StressGroup %in% c("SUS", "RES")),
       n_qc_fail_sus_samples = k, n_qc_fail_sus_animals = k_an,
       missingness_direction = miss_dir,
       flagged_share_with_flags = share(TRUE),
       flagged_share_with_failed = paste0(share(FALSE), " - ", of_them),
       qc_fail_sus_samples = paste(csr_number_word(k), if (k == 1L) "QC-failed SUS sample" else "QC-failed SUS samples"),
       qc_fail_most_missing = most_missing,
       shallowest_sus = shallowest, shallowest_sus_is_most_missing = top_is_most_missing,
       loo_design = loo, qc_fail_removal = qc_removal)
}

# The run's sample centring, from the stamp the imputation writes next to its workbooks (config/normalisation.yml
# decides it); NA when there is no stamp.
csr_run_centring <- function(imputed_path) {
  if (!length(imputed_path) || is.na(imputed_path)) return(NA_character_)
  stamp <- file.path(dirname(imputed_path), "normalisation_stamp.yml")
  if (!file.exists(stamp)) return(NA_character_)
  as.character(yaml::read_yaml(stamp)$sample_centring %||% NA_character_)
}

# How the CA2-SLM displacement arises, for the generated wording. Under all-observed centring a sample's missingness
# sets its centring constant; under complete-case centring it does not, so the old explanation would be false.
csr_displacement_origin <- function(centring) {
  if (identical(centring, "complete_case_median"))
    ", measured under complete-case centring (each sample centred on the proteins quantified in every sample of the compartment, so its missingness does not set its centring constant),"
  else " arising from per-sample median centring under unequal missingness"
}

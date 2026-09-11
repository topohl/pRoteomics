# Shared statistics for bilateral (left vs right) reproducibility.
#
# THE QUESTION THESE ANSWER
#   Not "is the effect significant on both sides?" - requiring independent
#   significance in each hemisphere just discards power and rewards large n.
#   The question is whether an effect REPRODUCES in magnitude and direction on
#   the opposite side of the same brain.
#
# WHAT THESE ARE NOT
#   Two hemispheres of one animal are one biological replicate. Agreement
#   between sides is evidence about measurement reliability and within-animal
#   symmetry. It is NOT independent replication and does not raise n.

# Pair left and right values on an exact key.
#
# A left value may only ever meet the right value carrying the SAME id columns.
# Returns one row per id combination with `left` and `right` columns; a side
# that is absent becomes NA rather than being silently dropped or recycled.
sps_pair_sides <- function(df, id_cols, side_col = "Hemisphere",
                           value_col = "value") {
  need <- c(id_cols, side_col, value_col)
  missing <- setdiff(need, names(df))
  if (length(missing)) {
    stop("Pairing table is missing column(s): ", paste(missing, collapse = ", "),
         ".", call. = FALSE)
  }
  side <- as.character(df[[side_col]])
  if (!all(side %in% c("L", "R"))) {
    stop("Side column must be normalized to L/R before pairing.", call. = FALSE)
  }
  key <- do.call(paste, c(df[id_cols], list(sep = "\037")))
  if (anyDuplicated(paste(key, side))) {
    stop("More than one value per id x side; collapse technical rows within a ",
         "side before pairing.", call. = FALSE)
  }
  ids <- df[!duplicated(key), id_cols, drop = FALSE]
  ukey <- key[!duplicated(key)]
  ord <- order(ukey, method = "radix")
  ids <- ids[ord, , drop = FALSE]; ukey <- ukey[ord]
  val <- function(s) {
    v <- rep(NA_real_, length(ukey))
    hit <- side == s
    v[match(key[hit], ukey)] <- as.numeric(df[[value_col]][hit])
    v
  }
  out <- cbind(ids, data.frame(left = val("L"), right = val("R"),
                               stringsAsFactors = FALSE))
  rownames(out) <- NULL
  out
}

# Correlation that returns NA instead of erroring on degenerate input.
sps_safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) return(NA_real_)
  if (stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(stats::cor(x[ok], y[ok], method = method))
}

# Agreement summary for one paired vector.
#
# Reports BOTH a signed mean difference (does one side sit systematically
# higher?) and unsigned error (how far apart are they?), because a constant
# offset and random disagreement are different findings.
sps_paired_agreement <- function(left, right) {
  ok <- is.finite(left) & is.finite(right)
  l <- left[ok]; r <- right[ok]
  n <- length(l)
  d <- l - r
  slope <- if (n >= 3L && stats::sd(r) > 0) {
    unname(stats::coef(stats::lm(l ~ r))[["r"]])
  } else NA_real_
  r2 <- if (n >= 3L) {
    pr <- sps_safe_cor(l, r, "pearson"); if (is.na(pr)) NA_real_ else pr^2
  } else NA_real_
  data.frame(
    n_pairs = n,
    pearson_r = sps_safe_cor(l, r, "pearson"),
    spearman_rho = sps_safe_cor(l, r, "spearman"),
    regression_slope = slope,
    r_squared = r2,
    sign_agreement_fraction = if (n) mean(sign(l) == sign(r)) else NA_real_,
    MAE = if (n) mean(abs(d)) else NA_real_,
    RMSE = if (n) sqrt(mean(d^2)) else NA_real_,
    median_abs_difference = if (n) stats::median(abs(d)) else NA_real_,
    p90_abs_difference = if (n) unname(stats::quantile(abs(d), 0.9, na.rm = TRUE)) else NA_real_,
    mean_signed_L_minus_R = if (n) mean(d) else NA_real_,
    stringsAsFactors = FALSE
  )
}

# Descriptive class naming WHICH reproducibility holds.
#
# Deliberately not a quality score and deliberately not a pass/fail: poor
# bilateral agreement may be genuine hemispheric asymmetry rather than a
# measurement problem, so the classes describe the pattern and stop there.
sps_reproducibility_class <- function(absolute_r, profile_r,
                                      mean_signed_difference,
                                      median_abs_difference,
                                      strong = 0.7, weak = 0.3) {
  n <- max(length(absolute_r), length(profile_r))
  rep_len_ <- function(x) if (length(x) == n) x else rep(x, length.out = n)
  a <- rep_len_(absolute_r); p <- rep_len_(profile_r)
  off <- abs(rep_len_(mean_signed_difference))
  mad <- rep_len_(median_abs_difference)
  offset_dominates <- is.finite(off) & is.finite(mad) & mad > 0 & off >= 0.5 * mad

  out <- rep("insufficient_data", n)
  known <- is.finite(a) | is.finite(p)
  out[known] <- "poor_bilateral_reproducibility"
  out[known & is.finite(p) & p >= strong & is.finite(a) & a >= strong] <-
    "reproducible_level_and_pattern"
  out[known & is.finite(p) & p >= strong & is.finite(a) & a < strong] <-
    "reproducible_spatial_pattern_only"
  out[known & is.finite(p) & p >= strong & is.finite(a) & a < strong &
        offset_dominates[known & is.finite(p) & p >= strong & is.finite(a) & a < strong]] <-
    "systematic_side_offset_pattern_preserved"
  out[known & is.finite(a) & a >= strong & (!is.finite(p) | p < weak)] <-
    "reproducible_level_only"
  out
}

# Sign agreement between two effect vectors, ignoring exact zeros.
sps_sign_agreement <- function(a, b) {
  ok <- is.finite(a) & is.finite(b) & a != 0 & b != 0
  if (!any(ok)) return(NA_real_)
  mean(sign(a[ok]) == sign(b[ok]))
}

# Rank-position summary for a marker set discovered on one side and evaluated on
# the other. `discovered` must come from the discovery side ONLY.
sps_rank_transfer <- function(discovered, evaluation_stat, all_ids,
                              higher_is_better = TRUE) {
  discovered <- intersect(as.character(discovered), as.character(all_ids))
  if (!length(discovered)) {
    return(data.frame(n_transferred = 0L, median_eval_rank_fraction = NA_real_,
                      auc_like = NA_real_, sign_retained_fraction = NA_real_,
                      stringsAsFactors = FALSE))
  }
  stat <- as.numeric(evaluation_stat)
  names(stat) <- as.character(all_ids)
  r <- rank(if (higher_is_better) -stat else stat, na.last = "keep")
  rf <- r / sum(!is.na(r))
  hit <- rf[discovered]
  # rank-fraction based AUC: probability a discovered marker outranks a random
  # non-discovered one.
  rest <- rf[setdiff(names(rf), discovered)]
  auc <- if (length(rest) && length(hit)) {
    mean(outer(hit, rest, "<"), na.rm = TRUE)
  } else NA_real_
  data.frame(
    n_transferred = length(discovered),
    median_eval_rank_fraction = stats::median(hit, na.rm = TRUE),
    auc_like = auc,
    sign_retained_fraction = mean(sign(stat[discovered]) ==
                                    sign(stats::median(stat[discovered], na.rm = TRUE)),
                                  na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

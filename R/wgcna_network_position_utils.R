# Pure helpers for the SUS - RES network-position audit.
#
# QUESTION
#   Are SUS - RES FDR-supported proteins distributed differently WITHIN their
#   frozen WGCNA modules than an equal number of DA-tested proteins drawn from
#   the same modules?
#
# WHY A PERMUTATION IS NEEDED
#   The A1/A2/D tier split cannot answer this. Those tiers are defined BY kME
#   (top-10 hub, |kME| >= 0.60), so their composition is largely fixed by module
#   size: only 10 of a ~335-member module can be top-10 hubs, so a small A1 count
#   is expected under the null and is not evidence of peripheral targeting.
#
# WHAT THIS FILE DOES NOT DO
#   No WGCNA and no differential abundance is recomputed. kME, ranks and BH FDR
#   are copied from frozen artifacts. Nothing here is causal or independent
#   evidence; it is a descriptive positional comparison against a
#   module-stratified null.

# `%||%` comes from the canonical R/null_coalescing.R, loaded via R/paths.R.

# ---------------------------------------------------------------- vocabulary

wnp_contract_version <- function() "wgcna_sus_res_network_position_v1"

wnp_default_permutations <- function() 50000L

# For each statistic, which direction of the observed value means the hits sit
# MORE CENTRALLY in their modules. Used to orient the one-sided p-values.
wnp_statistic_directions <- function() {
  c(
    median_rank_fraction = "low",
    mean_rank_fraction = "low",
    median_abs_kME = "high",
    fraction_core_kME_0.6 = "high",
    fraction_top25_hub = "high",
    fraction_top10_hub = "high"
  )
}

wnp_statistic_labels <- function() {
  c(
    median_rank_fraction = "Median within-module |kME| rank fraction",
    mean_rank_fraction = "Mean within-module |kME| rank fraction",
    median_abs_kME = "Median |kME|",
    `fraction_core_kME_0.6` = "Fraction with |kME| >= 0.60",
    fraction_top25_hub = "Fraction among module top-25 hubs",
    fraction_top10_hub = "Fraction among module top-10 hubs"
  )
}

.wnp_stop <- function(...) stop(..., call. = FALSE)

.wnp_require <- function(data, columns, label) {
  if (!is.data.frame(data)) .wnp_stop(label, " must be a data frame.")
  missing <- setdiff(columns, names(data))
  if (length(missing)) {
    .wnp_stop(label, " is missing required column(s): ",
              paste(missing, collapse = ", "), ".")
  }
  invisible(TRUE)
}

.wnp_is_true <- function(x) {
  if (is.logical(x)) return(x %in% TRUE)
  toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1")
}

# --------------------------------------------------------- rank fraction

# Midpoint rank fraction: (rank - 0.5) / n_members.
#
# The midpoint form is used rather than rank/n so that a uniformly drawn set has
# an expected rank fraction of exactly 0.5 regardless of module size, which
# keeps modules of very different sizes comparable when they are pooled.
# 0 = most central (highest |kME|), 1 = most peripheral.
wnp_rank_fraction <- function(rank, n_members) {
  rank <- suppressWarnings(as.numeric(rank))
  n_members <- suppressWarnings(as.numeric(n_members))
  out <- (rank - 0.5) / n_members
  out[!is.finite(out) | !is.finite(n_members) | n_members <= 0] <- NA_real_
  out
}

# ------------------------------------------------------------- statistics

# The six positional statistics, computed over a set of row indices.
# `rows` selects the observed hits, or one permuted draw.
wnp_position_statistics <- function(data, rows) {
  .wnp_require(
    data,
    c("rank_fraction", "abs_kME", "is_core_kME_0.6", "is_top_hub_25",
      "is_top10_module_hub"),
    "Network-position table"
  )
  if (!length(rows)) {
    return(stats::setNames(
      rep(NA_real_, length(wnp_statistic_directions())),
      names(wnp_statistic_directions())
    ))
  }
  rf <- suppressWarnings(as.numeric(data$rank_fraction[rows]))
  kme <- suppressWarnings(as.numeric(data$abs_kME[rows]))
  c(
    median_rank_fraction = stats::median(rf, na.rm = TRUE),
    mean_rank_fraction = mean(rf, na.rm = TRUE),
    median_abs_kME = stats::median(kme, na.rm = TRUE),
    `fraction_core_kME_0.6` = mean(.wnp_is_true(data[["is_core_kME_0.6"]][rows])),
    fraction_top25_hub = mean(.wnp_is_true(data$is_top_hub_25[rows])),
    fraction_top10_hub = mean(.wnp_is_true(data$is_top10_module_hub[rows]))
  )
}

# ------------------------------------------------ module-stratified null

# Draw one null set: within every module, sample exactly as many DA-tested
# proteins as that module actually contributed hits.
#
# Stratifying by module is what makes this a fair comparison. A global shuffle
# would confound module size, module kME distribution and how many hits each
# module produced; stratifying holds all three fixed and asks only whether the
# hits sit unusually within their own modules.
wnp_draw_null_rows <- function(module_rows, hits_per_module) {
  out <- vector("list", length(module_rows))
  nms <- names(module_rows)
  for (i in seq_along(module_rows)) {
    k <- hits_per_module[[nms[[i]]]]
    if (is.na(k) || k <= 0L) next
    pool <- module_rows[[i]]
    # sample() would treat a length-1 pool as 1:n; guard it explicitly.
    out[[i]] <- if (length(pool) == 1L) pool else sample(pool, k)
  }
  unlist(out, use.names = FALSE)
}

# Full module-stratified permutation. Returns the observed statistics and the
# B x n_statistics null matrix.
wnp_permute <- function(data, n_permutations = wnp_default_permutations(),
                        seed = 20260911L) {
  .wnp_require(data, c("ModuleID", "is_hit"), "Network-position table")
  n_permutations <- as.integer(n_permutations)
  if (is.na(n_permutations) || n_permutations < 1L) {
    .wnp_stop("n_permutations must be a positive integer.")
  }
  hit <- .wnp_is_true(data$is_hit)
  observed <- wnp_position_statistics(data, which(hit))

  module_rows <- split(seq_len(nrow(data)), as.character(data$ModuleID))
  hits_per_module <- vapply(module_rows, function(rows) sum(hit[rows]),
                            integer(1))
  # Modules that contributed no hit contribute nothing to either side.
  keep <- hits_per_module > 0L
  module_rows <- module_rows[keep]
  hits_per_module <- hits_per_module[keep]

  # A module can only supply a null draw if it has at least as many tested
  # members as it had hits. By construction it always does (hits are tested),
  # but assert it rather than silently under-drawing.
  short <- vapply(seq_along(module_rows), function(i) {
    length(module_rows[[i]]) < hits_per_module[[i]]
  }, logical(1))
  if (any(short)) {
    .wnp_stop("Module(s) ", paste(names(module_rows)[short], collapse = ", "),
              " have fewer DA-tested members than observed hits.")
  }

  null <- matrix(
    NA_real_, nrow = n_permutations, ncol = length(observed),
    dimnames = list(NULL, names(observed))
  )
  if (length(module_rows)) {
    old_seed <- if (exists(".Random.seed", envir = globalenv())) {
      get(".Random.seed", envir = globalenv())
    } else NULL
    set.seed(seed)
    on.exit({
      if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = globalenv())
    }, add = TRUE)
    for (b in seq_len(n_permutations)) {
      null[b, ] <- wnp_position_statistics(
        data, wnp_draw_null_rows(module_rows, hits_per_module)
      )
    }
  }
  list(observed = observed, null = null, n_permutations = n_permutations,
       n_hits = sum(hit), n_eligible = nrow(data),
       n_modules_with_hits = length(module_rows), seed = seed)
}

# ------------------------------------------------------- permutation summary

# Observed vs null, with (1 + count) / (B + 1) corrected one-sided p-values
# oriented by wnp_statistic_directions(), plus a two-sided p-value.
wnp_permutation_summary <- function(permutation, dataset = NA_character_,
                                    scope = NA_character_) {
  observed <- permutation$observed
  null <- permutation$null
  B <- permutation$n_permutations
  directions <- wnp_statistic_directions()
  labels <- wnp_statistic_labels()

  rows <- lapply(names(observed), function(nm) {
    obs <- observed[[nm]]
    draws <- null[, nm]
    draws <- draws[is.finite(draws)]
    central_high <- identical(unname(directions[[nm]]), "high")

    # "at least as central as observed" and "at least as peripheral as observed"
    if (central_high) {
      n_as_central <- sum(draws >= obs)
      n_as_peripheral <- sum(draws <= obs)
    } else {
      n_as_central <- sum(draws <= obs)
      n_as_peripheral <- sum(draws >= obs)
    }
    p_central <- (1 + n_as_central) / (B + 1)
    p_peripheral <- (1 + n_as_peripheral) / (B + 1)

    data.frame(
      dataset = dataset,
      scope = scope,
      statistic = nm,
      statistic_label = unname(labels[[nm]]),
      central_direction = unname(directions[[nm]]),
      observed = obs,
      null_median = if (length(draws)) stats::median(draws) else NA_real_,
      null_mean = if (length(draws)) mean(draws) else NA_real_,
      null_q025 = if (length(draws)) unname(stats::quantile(draws, 0.025)) else NA_real_,
      null_q975 = if (length(draws)) unname(stats::quantile(draws, 0.975)) else NA_real_,
      p_more_central = p_central,
      p_more_peripheral = p_peripheral,
      p_two_sided = min(1, 2 * min(p_central, p_peripheral)),
      n_permutations = B,
      n_hits = permutation$n_hits,
      n_eligible = permutation$n_eligible,
      n_modules_with_hits = permutation$n_modules_with_hits,
      permutation_seed = permutation$seed,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# ------------------------------------------------------- module-level summary

# Descriptive per-module hit topology plus an EXPLORATORY hypergeometric test of
# whether a module carries more hits than its share of the DA-tested universe.
#
# The test is exploratory only: it uses the DA-tested universe and is BH
# corrected across modules within a dataset, but it does not account for the
# differing number of spatial contexts, differing power, or the correlation
# between spatial contexts.
wnp_module_summary <- function(data) {
  .wnp_require(
    data,
    c("dataset", "ModuleID", "is_hit", "rank_fraction", "abs_kME"),
    "Network-position table"
  )
  if (!nrow(data)) {
    return(data.frame(
      dataset = character(), ModuleID = character(),
      n_da_eligible = integer(), n_hits = integer(), hit_fraction = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  hit <- .wnp_is_true(data$is_hit)
  key <- paste(as.character(data$dataset), as.character(data$ModuleID),
               sep = "\r")
  idx <- split(seq_len(nrow(data)), key)
  parts <- do.call(rbind, strsplit(names(idx), "\r", fixed = TRUE))

  total_tested <- vapply(split(hit, as.character(data$dataset)), length, integer(1))
  total_hits <- vapply(split(hit, as.character(data$dataset)), sum, integer(1))

  rows <- lapply(seq_along(idx), function(i) {
    rows_i <- idx[[i]]
    ds <- parts[i, 1]
    h <- hit[rows_i]
    rf <- suppressWarnings(as.numeric(data$rank_fraction[rows_i]))
    kme <- suppressWarnings(as.numeric(data$abs_kME[rows_i]))
    n <- length(rows_i)
    k <- sum(h)
    # exploratory hypergeometric: P(X >= k) for this module vs the dataset
    p_over <- if (k > 0L) {
      stats::phyper(k - 1L, m = n, n = total_tested[[ds]] - n,
                    k = total_hits[[ds]], lower.tail = FALSE)
    } else NA_real_
    data.frame(
      dataset = ds,
      ModuleID = parts[i, 2],
      n_da_eligible = n,
      n_hits = k,
      hit_fraction = k / n,
      n_tier_A1 = sum(.wnp_is_true(data[["is_tier_A1"]][rows_i]) %||% FALSE),
      n_tier_A2 = sum(.wnp_is_true(data[["is_tier_A2"]][rows_i]) %||% FALSE),
      n_tier_D = sum(.wnp_is_true(data[["is_tier_D"]][rows_i]) %||% FALSE),
      median_hit_rank_fraction = if (k) stats::median(rf[h], na.rm = TRUE) else NA_real_,
      median_nonhit_rank_fraction = if (n - k) stats::median(rf[!h], na.rm = TRUE) else NA_real_,
      median_hit_abs_kME = if (k) stats::median(kme[h], na.rm = TRUE) else NA_real_,
      median_nonhit_abs_kME = if (n - k) stats::median(kme[!h], na.rm = TRUE) else NA_real_,
      exploratory_hypergeometric_p = p_over,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  # BH across modules WITHIN a dataset
  out$exploratory_hypergeometric_BH <- NA_real_
  for (ds in unique(out$dataset)) {
    sel <- out$dataset == ds & is.finite(out$exploratory_hypergeometric_p)
    if (any(sel)) {
      out$exploratory_hypergeometric_BH[sel] <-
        stats::p.adjust(out$exploratory_hypergeometric_p[sel], method = "BH")
    }
  }
  out[order(out$dataset, out$ModuleID, method = "radix"), , drop = FALSE]
}

# ------------------------------------------------------------ interpretation

# One-line, deliberately conservative reading of a permutation summary row.
wnp_interpret <- function(summary_row, alpha = 0.05) {
  p_c <- suppressWarnings(as.numeric(summary_row$p_more_central))
  p_p <- suppressWarnings(as.numeric(summary_row$p_more_peripheral))
  ifelse(
    !is.finite(p_c) | !is.finite(p_p), "not_evaluable",
    ifelse(p_c <= alpha, "more_central_than_null",
           ifelse(p_p <= alpha, "more_peripheral_than_null",
                  "consistent_with_null"))
  )
}

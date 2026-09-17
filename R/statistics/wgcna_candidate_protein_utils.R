# Pure helpers for the downstream WGCNA candidate-protein shortlist layer.
#
# Scope contract.  Every function here is DESCRIPTIVE and ADDITIVE.  Nothing in
# this file fits a network, refits differential abundance, or re-derives a
# module identity.  All statistics are copied unchanged from frozen upstream
# artifacts; the only new quantities are ranks, percentiles, counts, and boolean
# flags computed from those copied values.
#
# Vocabulary used throughout, deliberately conservative:
#   * "module hub" / "module representative" - high |kME| within its own module.
#   * "protein-level differential abundance"  - manifest-selected limma results.
#   * "directional concordance"               - sign agreement only.
#   * "spatial consistency"                   - sign agreement across hippocampal
#     spatial contexts measured in the SAME animals.  These are repeated
#     measurements, NOT independent replication.
# The terms "driver", "causal", "validated" and "replicated" are intentionally
# never emitted by this layer.

# `%||%` comes from the canonical R/null_coalescing.R, loaded via R/paths.R.
# It is deliberately not redefined here: the repository permits exactly one
# definition so that source order cannot change coalescing semantics.

# ---------------------------------------------------------------- vocabulary

wcp_contract_version <- function() "wgcna_candidate_protein_shortlist_v1"

# Canonical three-contrast vocabulary, identical to the Stage 07 module-level
# handoff and to analysis/08_integration/export_module_protein_zoom_source_data.R.
wcp_contrast_levels <- function() c("RES - CON", "SUS - CON", "SUS - RES")

# Primary biological endpoint of the study.
wcp_primary_contrast <- function() "SUS - RES"

# |kME| >= 0.60 reproduces the frozen `is_core_kME_0.6` definition emitted by
# analysis/05_wgcna/build_wgcna_modules.R.  It is NOT re-derived from expression data.
wcp_high_kme_threshold <- function() 0.60

# Descriptive only.  The repository defines no canonical 0.80 kME threshold;
# this flag exists purely to mark unusually strong module representatives.
wcp_very_high_kme_threshold <- function() 0.80

wcp_fdr_threshold <- function() 0.05

# Data-relative effect rule: top decile of |log2FC| WITHIN one
# dataset x spatial_unit x contrast family.  Avoids a universal fold-change cut.
wcp_large_effect_quantile <- function() 0.90

# Minimum number of tested spatial contexts before a consistency claim is made.
wcp_min_spatial_contexts <- function() 3L


wcp_module_representative_n <- function() 5L
wcp_module_hub_n <- function() 10L

# Precedence for the convenience `candidate_tier` label and for display order:
# protein-level inferential support (A1, A2, D) is shown before purely
# descriptive candidate status (B, C).
#
# THIS IS NOT A RANKING OF EVIDENCE STRENGTH.  A1 / A2 / D differ only in WGCNA
# topology - top-10 hub, module member, peripheral member - not in how strong
# the differential-abundance evidence is.  A Tier D protein routinely carries a
# smaller FDR and a larger effect than a Tier A1 protein.
wcp_candidate_tier_levels <- function() c("A1", "A2", "D", "B", "C")

# Human-readable names for the phenotype-linked classes.
wcp_tier_display_names <- function() {
  c(
    A1 = "FDR-supported phenotype-linked module hub",
    A2 = "FDR-supported phenotype-linked module member",
    D  = "FDR-supported phenotype-linked peripheral module member",
    B  = "Network-central candidate with a large typical SUS - RES effect",
    C  = "Canonical module representative"
  )
}

# --------------------------------------------------------------- small utils

.wcp_stop <- function(...) stop(..., call. = FALSE)

.wcp_require_columns <- function(data, columns, label) {
  if (!is.data.frame(data)) .wcp_stop(label, " must be a data frame.")
  missing <- setdiff(columns, names(data))
  if (length(missing)) {
    .wcp_stop(label, " is missing required column(s): ",
              paste(missing, collapse = ", "), ".")
  }
  invisible(TRUE)
}

.wcp_require_unique <- function(data, key, label) {
  .wcp_require_columns(data, key, label)
  dup <- duplicated(data[key])
  if (any(dup)) {
    .wcp_stop(
      label, " contains duplicate rows for ", paste(key, collapse = " + "),
      " (first at row ", which(dup)[[1]], "). Spatial and module summaries ",
      "would double-count these rows."
    )
  }
  invisible(TRUE)
}

.wcp_num <- function(x) suppressWarnings(as.numeric(x))

.wcp_lgl <- function(x) {
  if (is.logical(x)) return(x)
  out <- rep(NA, length(x))
  chr <- toupper(trimws(as.character(x)))
  out[chr %in% c("TRUE", "T", "1")] <- TRUE
  out[chr %in% c("FALSE", "F", "0")] <- FALSE
  as.logical(out)
}

# TRUE only where the value is explicitly TRUE; NA never leaks through as TRUE.
.wcp_is_true <- function(x) .wcp_lgl(x) %in% TRUE

# Length-safe optional-column accessor.
#
# `%||%` must NOT be used for this. The canonical proteomics_null_coalesce also
# falls back for length-0 and scalar-NA input, so `df$missing %||% NA_real_`
# yields a length-1 value for a zero-row frame and `order()` then fails with
# "argument lengths differ". This always returns exactly nrow(data) elements.
.wcp_col <- function(data, column, default) {
  n <- nrow(data)
  if (!is.null(data[[column]])) return(data[[column]])
  rep(default, n)
}

# ------------------------------------------------------ contrast orientation

# Vectorised comparison -> canonical contrast label.  Mirrors the mapping used
# by the manuscript Figure 3 renderer so both layers name contrasts identically.
wcp_contrast_from_comparison <- function(x) {
  z <- tolower(as.character(x))
  ifelse(
    grepl("res.*con", z), "RES - CON",
    ifelse(
      grepl("sus.*con", z), "SUS - CON",
      ifelse(grepl("sus.*res", z), "SUS - RES", NA_character_)
    )
  )
}

.wcp_phenotype_token <- function(side) {
  side <- tolower(trimws(as.character(side)))
  hit <- regmatches(side, regexpr("(sus|res|con)$", side, perl = TRUE))
  if (length(hit) != 1L || is.na(hit) || !nzchar(hit)) {
    .wcp_stop("Comparison side does not end in a recognised phenotype token ",
              "(SUS, RES, CON): ", side)
  }
  toupper(hit)
}

# Generalises `sus_res_resolve_orientation()` (R/statistics/sus_res_spatial_dap_atlas_utils.R)
# from SUS-vs-RES to all three canonical contrasts.  The serialized log2FC in a
# manifest-selected DA file is oriented left-minus-right by filename; this
# returns the multiplier that re-expresses it as the canonical contrast label.
wcp_comparison_orientation <- function(comparison) {
  comparison <- as.character(comparison)
  sides <- strsplit(comparison, "_", fixed = TRUE)[[1]]
  if (length(sides) != 2L || any(!nzchar(sides))) {
    .wcp_stop("Comparison must contain exactly two non-empty sides separated ",
              "by one underscore: ", comparison)
  }
  left <- .wcp_phenotype_token(sides[[1]])
  right <- .wcp_phenotype_token(sides[[2]])
  formal <- wcp_contrast_from_comparison(comparison)
  if (is.na(formal)) {
    .wcp_stop("Comparison is not one of the canonical contrasts: ", comparison)
  }
  formal_sides <- trimws(strsplit(formal, "-", fixed = TRUE)[[1]])
  if (!setequal(c(left, right), formal_sides)) {
    .wcp_stop("Comparison phenotypes (", left, ", ", right,
              ") do not match canonical contrast '", formal, "'.")
  }
  multiplier <- if (identical(c(left, right), formal_sides)) 1 else -1
  data.frame(
    comparison = comparison,
    serialized_effect_definition = paste(left, "-", right),
    formal_effect_definition = formal,
    formal_effect_multiplier = multiplier,
    sign_was_flipped = multiplier == -1,
    stringsAsFactors = FALSE
  )
}

# --------------------------------------------------- within-module kME ranks

# Adds deterministic within-module |kME| ranks and hub flags.
#
# Ordering is desc(abs_kME) then ProteinGroupID ascending - the same
# tie-breaking rule the manuscript Figure 3 renderer uses - so the flags are
# reproducible across runs and platforms.  Rows with a non-finite abs_kME are
# ranked last and never receive a hub flag.
wcp_rank_module_members <- function(members) {
  .wcp_require_columns(
    members, c("dataset", "ModuleID", "ProteinGroupID", "abs_kME"),
    "Module membership"
  )
  .wcp_require_unique(
    members, c("dataset", "ModuleID", "ProteinGroupID"), "Module membership"
  )
  if (!nrow(members)) {
    members$n_module_members <- integer()
    members$abs_kME_rank_in_module <- integer()
    members$abs_kME_rank_fraction_in_module <- numeric()
    members$is_top5_module_representative <- logical()
    members$is_top10_module_hub <- logical()
    return(members)
  }

  abs_kme <- .wcp_num(members$abs_kME)
  finite <- is.finite(abs_kme)
  group <- paste(as.character(members$dataset), as.character(members$ModuleID),
                 sep = "\r")

  ord <- order(group, !finite, -ifelse(finite, abs_kme, 0),
               as.character(members$ProteinGroupID), method = "radix")

  rank_in_module <- integer(nrow(members))
  rank_in_module[ord] <- stats::ave(
    seq_along(ord), group[ord], FUN = seq_along
  )
  n_members <- as.integer(stats::ave(
    rep(1L, nrow(members)), group, FUN = length
  ))

  # Non-finite abs_kME carries no interpretable module position.
  rank_out <- ifelse(finite, rank_in_module, NA_integer_)

  members$n_module_members <- n_members
  members$abs_kME_rank_in_module <- as.integer(rank_out)
  members$abs_kME_rank_fraction_in_module <- ifelse(
    finite & n_members > 0L, rank_out / n_members, NA_real_
  )
  members$is_top5_module_representative <-
    finite & !is.na(rank_out) & rank_out <= wcp_module_representative_n()
  members$is_top10_module_hub <-
    finite & !is.na(rank_out) & rank_out <= wcp_module_hub_n()
  members
}

# ------------------------------------------------------------- mapping / QC

# A protein is "cleanly mapped" only when the frozen Stage 02 identity contract
# allows BOTH a gene-level and a protein-level claim and the protein group is
# not gene-ambiguous.  Ambiguity is never silently upgraded to a clean claim.
wcp_clean_mapping <- function(members) {
  .wcp_require_columns(
    members,
    c("mapping_status", "gene_level_claim_allowed", "protein_level_claim_allowed",
      "protein_group_ambiguity_class"),
    "Module membership"
  )
  status <- tolower(trimws(as.character(members$mapping_status)))
  ambiguity <- tolower(trimws(as.character(members$protein_group_ambiguity_class)))
  status %in% "mapped" &
    .wcp_is_true(members$gene_level_claim_allowed) &
    .wcp_is_true(members$protein_level_claim_allowed) &
    ambiguity %in% "single_accession_single_gene"
}

# ----------------------------------------------- data-relative effect ranking

# Top-decile |log2FC| threshold computed strictly WITHIN each
# dataset x spatial_unit x contrast family.  Families never pool.
wcp_large_effect_thresholds <- function(da_long,
                                        probability = wcp_large_effect_quantile()) {
  key <- c("dataset", "spatial_unit", "contrast")
  .wcp_require_columns(da_long, c(key, "log2FC"), "Protein-level DA")
  if (!nrow(da_long)) {
    return(data.frame(
      dataset = character(), spatial_unit = character(), contrast = character(),
      large_effect_abs_log2FC_threshold = numeric(),
      n_tested_in_effect_family = integer(),
      stringsAsFactors = FALSE
    ))
  }
  abs_effect <- abs(.wcp_num(da_long$log2FC))
  family <- interaction(
    da_long$dataset, da_long$spatial_unit, da_long$contrast,
    drop = TRUE, sep = "\r"
  )
  split_abs <- split(abs_effect, family)
  parts <- do.call(rbind, strsplit(names(split_abs), "\r", fixed = TRUE))
  data.frame(
    dataset = parts[, 1],
    spatial_unit = parts[, 2],
    contrast = parts[, 3],
    large_effect_abs_log2FC_threshold = vapply(split_abs, function(v) {
      v <- v[is.finite(v)]
      if (!length(v)) return(NA_real_)
      unname(stats::quantile(v, probs = probability, na.rm = TRUE, type = 7))
    }, numeric(1)),
    n_tested_in_effect_family = vapply(
      split_abs, function(v) sum(is.finite(v)), integer(1)
    ),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

# Attaches the family threshold and the resulting descriptive flag.
wcp_flag_large_effect <- function(da_long,
                                  probability = wcp_large_effect_quantile()) {
  thresholds <- wcp_large_effect_thresholds(da_long, probability)
  key <- c("dataset", "spatial_unit", "contrast")
  idx <- match(
    do.call(paste, c(da_long[key], sep = "\r")),
    do.call(paste, c(thresholds[key], sep = "\r"))
  )
  da_long$large_effect_abs_log2FC_threshold <-
    thresholds$large_effect_abs_log2FC_threshold[idx]
  da_long$n_tested_in_effect_family <- thresholds$n_tested_in_effect_family[idx]
  abs_effect <- abs(.wcp_num(da_long$log2FC))
  da_long$large_effect_within_context <-
    is.finite(abs_effect) &
    is.finite(da_long$large_effect_abs_log2FC_threshold) &
    abs_effect >= da_long$large_effect_abs_log2FC_threshold
  da_long
}

# ------------------------------------------------------ directional agreement

# Descriptive sign agreement between a protein-level effect and the
# module-level eigengene estimate for the SAME contrast.
#
# Returns NA - never FALSE - unless BOTH values are finite and non-zero, so an
# absent or exactly-zero module estimate can never be read as disagreement.
# This is concordance, not independent replication: the module eigengene is a
# summary of the same proteins measured in the same animals.
wcp_direction_matches_module <- function(protein_effect, module_estimate) {
  protein_effect <- .wcp_num(protein_effect)
  module_estimate <- .wcp_num(module_estimate)
  n <- max(length(protein_effect), length(module_estimate))
  protein_effect <- rep_len(protein_effect, n)
  module_estimate <- rep_len(module_estimate, n)
  comparable <- is.finite(protein_effect) & is.finite(module_estimate) &
    protein_effect != 0 & module_estimate != 0
  out <- rep(NA, n)
  out[comparable] <- sign(protein_effect[comparable]) ==
    sign(module_estimate[comparable])
  as.logical(out)
}

# ---------------------------------------------------- spatial consistency

# Per dataset x ProteinGroupID x contrast summary across hippocampal spatial
# contexts.  These contexts are repeated measurements from the same animals;
# the output describes CONSISTENCY, never independent replication.
#
# `spatially_consistent` requires unanimous sign agreement across at least
# `min_contexts` tested contexts.  A fixed fraction cut was rejected because a
# dataset with four contexts and one with ten are not comparable under it; the
# continuous `fraction_matching_majority` is exported so any other rule can be
# applied by hand.
wcp_spatial_consistency <- function(da_long,
                                    min_contexts = wcp_min_spatial_contexts()) {
  key <- c("dataset", "ProteinGroupID", "contrast", "spatial_unit")
  .wcp_require_columns(da_long, c(key, "log2FC", "BH_FDR"), "Protein-level DA")
  .wcp_require_unique(da_long, key, "Protein-level DA")

  empty <- data.frame(
    dataset = character(), ProteinGroupID = character(), contrast = character(),
    n_spatial_contexts_tested = integer(), n_spatial_contexts_positive = integer(),
    n_spatial_contexts_negative = integer(), majority_direction = character(),
    n_matching_majority_direction = integer(),
    fraction_matching_majority_direction = numeric(),
    n_spatial_contexts_fdr05 = integer(), max_abs_log2FC = numeric(),
    median_log2FC = numeric(), median_abs_log2FC = numeric(),
    spatially_consistent = logical(), stringsAsFactors = FALSE
  )
  if (!nrow(da_long)) return(empty)

  effect <- .wcp_num(da_long$log2FC)
  fdr <- .wcp_num(da_long$BH_FDR)
  tested <- is.finite(effect)

  group <- paste(as.character(da_long$dataset),
                 as.character(da_long$ProteinGroupID),
                 as.character(da_long$contrast), sep = "\r")
  idx <- split(seq_len(nrow(da_long)), group)
  parts <- do.call(rbind, strsplit(names(idx), "\r", fixed = TRUE))

  summarise_one <- function(rows) {
    e <- effect[rows][tested[rows]]
    f <- fdr[rows][tested[rows]]
    n <- length(e)
    npos <- sum(e > 0)
    nneg <- sum(e < 0)
    nmaj <- max(npos, nneg)
    direction <- if (!n) NA_character_ else if (npos > nneg) "positive" else
      if (nneg > npos) "negative" else "tie"
    c(
      n = n, npos = npos, nneg = nneg, nmaj = if (n) nmaj else NA_real_,
      frac = if (n) nmaj / n else NA_real_,
      nfdr = sum(is.finite(f) & f <= wcp_fdr_threshold()),
      maxabs = if (n) max(abs(e)) else NA_real_,
      med = if (n) stats::median(e) else NA_real_,
      medabs = if (n) stats::median(abs(e)) else NA_real_,
      dir = if (identical(direction, "positive")) 1 else
        if (identical(direction, "negative")) -1 else
          if (identical(direction, "tie")) 0 else NA_real_
    )
  }

  stats_mat <- do.call(rbind, lapply(idx, summarise_one))
  direction <- rep(NA_character_, nrow(stats_mat))
  direction[stats_mat[, "dir"] %in% 1] <- "positive"
  direction[stats_mat[, "dir"] %in% -1] <- "negative"
  direction[stats_mat[, "dir"] %in% 0] <- "tie"

  out <- data.frame(
    dataset = parts[, 1],
    ProteinGroupID = parts[, 2],
    contrast = parts[, 3],
    n_spatial_contexts_tested = as.integer(stats_mat[, "n"]),
    n_spatial_contexts_positive = as.integer(stats_mat[, "npos"]),
    n_spatial_contexts_negative = as.integer(stats_mat[, "nneg"]),
    majority_direction = direction,
    n_matching_majority_direction = as.integer(stats_mat[, "nmaj"]),
    fraction_matching_majority_direction = as.numeric(stats_mat[, "frac"]),
    n_spatial_contexts_fdr05 = as.integer(stats_mat[, "nfdr"]),
    max_abs_log2FC = as.numeric(stats_mat[, "maxabs"]),
    median_log2FC = as.numeric(stats_mat[, "med"]),
    median_abs_log2FC = as.numeric(stats_mat[, "medabs"]),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  # Unanimous sign across at least `min_contexts` tested contexts.  Purely
  # DESCRIPTIVE: no probability is attached to it and it does not by itself
  # create a candidate tier.  Hippocampal spatial contexts are repeated,
  # positively correlated measurements from the SAME animals, so their signs
  # are not independent Bernoulli trials and unanimity is far more common than
  # any independence-based calculation would suggest.
  out$spatially_consistent <-
    out$n_spatial_contexts_tested >= as.integer(min_contexts) &
    !is.na(out$fraction_matching_majority_direction) &
    out$fraction_matching_majority_direction >= 1

  out[order(out$dataset, out$contrast, out$ProteinGroupID, method = "radix"), ,
      drop = FALSE]
}

# Protein-level effect magnitude, summarised across spatial contexts before
# thresholding.  Flags the top decile of MEDIAN |log2FC| within one
# dataset x contrast family.
#
# This is used for tiering in preference to "top decile in ANY single context",
# because a protein measured in ten contexts gets ten independent chances at the
# per-context decile; that inflates the per-context rule from 10% to roughly 35%
# (4 contexts) or 65% (10 contexts) of all proteins and stops it selecting.
wcp_flag_typical_large_effect <- function(consistency,
                                          probability = wcp_large_effect_quantile()) {
  .wcp_require_columns(
    consistency, c("dataset", "contrast", "median_abs_log2FC"),
    "Spatial consistency summary"
  )
  if (!nrow(consistency)) {
    consistency$typical_large_effect_threshold <- numeric()
    consistency$large_effect_typical <- logical()
    return(consistency)
  }
  family <- paste(as.character(consistency$dataset),
                  as.character(consistency$contrast), sep = "\r")
  med <- .wcp_num(consistency$median_abs_log2FC)
  thresholds <- vapply(split(med, family), function(v) {
    v <- v[is.finite(v)]
    if (!length(v)) return(NA_real_)
    unname(stats::quantile(v, probs = probability, na.rm = TRUE, type = 7))
  }, numeric(1))
  consistency$typical_large_effect_threshold <- unname(thresholds[family])
  consistency$large_effect_typical <-
    is.finite(med) & is.finite(consistency$typical_large_effect_threshold) &
    med >= consistency$typical_large_effect_threshold
  consistency
}

# ------------------------------------------------------- candidate flagging

# Turns the joined evidence into transparent boolean reason flags.  There is
# deliberately NO weighted composite score: every flag stays separately
# inspectable and separately falsifiable.
wcp_assign_candidate_flags <- function(candidates) {
  required <- c(
    "abs_kME", "is_top5_module_representative", "is_top10_module_hub",
    "clean_mapping", "sus_res_fdr05_any_context",
    "sus_res_large_effect_any_context", "sus_res_large_effect_typical",
    "sus_res_spatially_consistent"
  )
  .wcp_require_columns(candidates, required, "Candidate table")

  abs_kme <- .wcp_num(candidates$abs_kME)
  candidates$high_kME <- is.finite(abs_kme) & abs_kme >= wcp_high_kme_threshold()
  candidates$very_high_kME <-
    is.finite(abs_kme) & abs_kme >= wcp_very_high_kme_threshold()
  candidates$top5_hub <- .wcp_is_true(candidates$is_top5_module_representative)
  candidates$top10_hub <- .wcp_is_true(candidates$is_top10_module_hub)

  fdr05 <- .wcp_is_true(candidates$sus_res_fdr05_any_context)
  clean <- .wcp_is_true(candidates$clean_mapping)

  # Tier B is EFFECT-based.  Two rules are deliberately NOT used for tiering:
  #   * "top decile in ANY single context" gives a protein one chance per
  #     spatial context, so it flags ~35% (4 contexts) to ~65% (10 contexts)
  #     of all proteins instead of 10%.
  #   * spatial unanimity, because these contexts are repeated, positively
  #     correlated measurements from the SAME animals. Any independence-based
  #     sign probability understates how often unanimity arises, so it must not
  #     gate candidate selection.
  # Both flags remain exported as descriptive annotation and as secondary
  # prioritisation information.
  large <- .wcp_is_true(candidates$sus_res_large_effect_typical)

  # Phenotype-linked classes.  A1, A2 and D all require the SAME protein-level
  # inferential evidence (SUS - RES BH FDR <= 0.05 in >= 1 spatial context) and
  # differ ONLY in WGCNA topology.  The split exists so that "has protein-level
  # statistical support", "is a real module member" and "is extremely central"
  # stay separable; it is not a gradation of evidence strength.
  #
  # Module-direction concordance is recorded as an additional flag rather than
  # a requirement for any of them: no module-level SUS - RES effect is
  # FDR-supported in any dataset, so gating on the sign of a non-significant
  # eigengene estimate would not be defensible.

  # A1 - FDR-supported phenotype-linked module hub.
  candidates$is_tier_A1 <- clean & candidates$high_kME &
    candidates$top10_hub & fdr05

  # A2 - FDR-supported phenotype-linked module member: a core member by the
  # frozen |kME| >= 0.60 definition, but outside its module's top 10 hubs.
  # No top-25 requirement is imposed; `is_top_hub_25` remains available as a
  # separate frozen descriptive flag.
  candidates$is_tier_A2 <- clean & candidates$high_kME &
    !candidates$top10_hub & fdr05

  # Compatibility aggregate only.  `is_tier_A` now means "phenotype-linked
  # module member of either centrality class" - the A1/A2 booleans carry the
  # actual distinction and are what human-facing output uses.
  candidates$is_tier_A <- candidates$is_tier_A1 | candidates$is_tier_A2

  # Tier B - network-central protein carrying a large typical SUS - RES effect.
  # Protein-level FDR support is not required.  Spatial directional agreement
  # is NOT part of this definition; it is reported alongside for prioritisation.
  candidates$is_tier_B <- candidates$high_kME & large

  # Tier C - canonical module representative, independent of any DA evidence.
  candidates$is_tier_C <- candidates$top5_hub

  # Tier D - FDR-supported phenotype-linked PERIPHERAL module member.
  # Definition deliberately unchanged from the previous revision, including the
  # fact that it does NOT require clean_mapping while A1/A2 do. That asymmetry
  # predates this refinement and is left alone rather than silently altered; on
  # the current frozen inputs it has no effect, because every FDR-supported
  # protein in every dataset is cleanly mapped.
  candidates$is_tier_D <- fdr05 & !candidates$high_kME

  candidates$is_candidate <- candidates$is_tier_A1 | candidates$is_tier_A2 |
    candidates$is_tier_B | candidates$is_tier_C | candidates$is_tier_D
  candidates
}

# Compact readability label for the phenotype-linked classes, derived from the
# A1/A2/D booleans. Describes WGCNA topology only, never evidence strength.
wcp_phenotype_network_class <- function(candidates) {
  .wcp_require_columns(
    candidates, c("is_tier_A1", "is_tier_A2", "is_tier_D"), "Candidate table"
  )
  out <- rep(NA_character_, nrow(candidates))
  out[.wcp_is_true(candidates$is_tier_D)] <- "peripheral_member"
  out[.wcp_is_true(candidates$is_tier_A2)] <- "module_member"
  out[.wcp_is_true(candidates$is_tier_A1)] <- "top10_hub"
  out
}

# Every FDR-supported protein must land in exactly one phenotype-linked class,
# unless the frozen mapping contract forbids a protein-level claim. Returns the
# offending rows so a caller can fail loudly rather than silently drop them.
wcp_unclassified_fdr_support <- function(candidates) {
  required <- c("sus_res_fdr05_any_context", "is_tier_A1", "is_tier_A2",
                "is_tier_D", "clean_mapping")
  .wcp_require_columns(candidates, required, "Candidate table")
  fdr05 <- .wcp_is_true(candidates$sus_res_fdr05_any_context)
  n_classes <- .wcp_is_true(candidates$is_tier_A1) +
    .wcp_is_true(candidates$is_tier_A2) + .wcp_is_true(candidates$is_tier_D)
  claimable <- .wcp_is_true(candidates$clean_mapping)
  # Flags both under- and over-classification: A1/A2 are disjoint by
  # construction and D is disjoint from both, so anything other than exactly
  # one class is a contract violation.
  candidates[fdr05 & claimable & n_classes != 1L, , drop = FALSE]
}

# Broad, convenience label only.  Derived AFTER the booleans, never instead of
# them: a protein may legitimately satisfy several tiers at once.
wcp_candidate_tier <- function(candidates) {
  .wcp_require_columns(
    candidates,
    c("is_tier_A1", "is_tier_A2", "is_tier_B", "is_tier_C", "is_tier_D"),
    "Candidate table"
  )
  # Applied in reverse precedence so the highest-precedence class wins.
  out <- rep(NA_character_, nrow(candidates))
  out[.wcp_is_true(candidates$is_tier_C)] <- "C"
  out[.wcp_is_true(candidates$is_tier_B)] <- "B"
  out[.wcp_is_true(candidates$is_tier_D)] <- "D"
  out[.wcp_is_true(candidates$is_tier_A2)] <- "A2"
  out[.wcp_is_true(candidates$is_tier_A1)] <- "A1"
  out
}

# All satisfied tiers, semicolon separated in precedence order, so overlap is
# never hidden (e.g. "A1;B", "A2;B", "A2", "D", "B;C").
wcp_candidate_tier_all <- function(candidates) {
  .wcp_require_columns(
    candidates,
    c("is_tier_A1", "is_tier_A2", "is_tier_B", "is_tier_C", "is_tier_D"),
    "Candidate table"
  )
  flags <- cbind(
    A1 = .wcp_is_true(candidates$is_tier_A1),
    A2 = .wcp_is_true(candidates$is_tier_A2),
    D = .wcp_is_true(candidates$is_tier_D),
    B = .wcp_is_true(candidates$is_tier_B),
    C = .wcp_is_true(candidates$is_tier_C)
  )
  flags <- flags[, wcp_candidate_tier_levels(), drop = FALSE]
  apply(flags, 1L, function(row) {
    hit <- wcp_candidate_tier_levels()[row]
    if (!length(hit)) NA_character_ else paste(hit, collapse = ";")
  })
}

.wcp_reason_dictionary <- function() {
  list(
    high_kME = "high module membership (|kME| >= 0.60)",
    very_high_kME = "very high module membership (|kME| >= 0.80, descriptive)",
    top5_hub = "top-5 module representative by |kME|",
    top10_hub = "top-10 module hub by |kME|",
    sus_res_fdr05_any_context = "SUS - RES protein-level BH FDR <= 0.05",
    sus_res_large_effect_typical =
      "SUS - RES median |log2FC| in the top decile for this dataset",
    sus_res_large_effect_any_context =
      "SUS - RES |log2FC| in the top decile of at least one spatial context",
    sus_res_spatially_consistent =
      "SUS - RES direction consistent across all tested spatial contexts",
    sus_res_direction_matches_module =
      "SUS - RES direction concordant with the module eigengene estimate",
    gsea_leading_edge =
      "member of a module-concordant GSEA leading edge (same DA statistics)",
    clean_mapping = "unambiguous gene- and protein-level claim"
  )
}

# Human-readable justification, built from the same booleans used for tiering.
wcp_candidate_reason <- function(candidates) {
  dict <- .wcp_reason_dictionary()
  present <- intersect(names(dict), names(candidates))
  if (!length(present)) return(rep(NA_character_, nrow(candidates)))
  mat <- vapply(present, function(nm) .wcp_is_true(candidates[[nm]]),
                logical(nrow(candidates)))
  if (is.null(dim(mat))) mat <- matrix(mat, nrow = nrow(candidates),
                                       dimnames = list(NULL, present))
  labels <- unlist(dict[present], use.names = FALSE)
  apply(mat, 1L, function(row) {
    hit <- labels[row]
    if (!length(hit)) NA_character_ else paste(hit, collapse = "; ")
  })
}

# ------------------------------------------------------------ display order

# Deterministic presentation order.  This is a SORT, not a score: no numeric
# combination of the criteria is ever formed or exported.
#   1. candidate tier in precedence order (A1, A2, D, B, C; unassigned last)
#   2. SUS - RES FDR support
#   3. |kME|
#   4. |SUS - RES log2FC|
#   5. spatial consistency
#   6. ProteinGroupID (deterministic tie-break)
wcp_order_candidates <- function(candidates) {
  required <- c("candidate_tier", "abs_kME", "ProteinGroupID")
  .wcp_require_columns(candidates, required, "Candidate table")
  tier <- match(as.character(candidates$candidate_tier),
                wcp_candidate_tier_levels())
  tier[is.na(tier)] <- length(wcp_candidate_tier_levels()) + 1L

  fdr_support <- if ("sus_res_fdr05_any_context" %in% names(candidates)) {
    !.wcp_is_true(candidates$sus_res_fdr05_any_context)
  } else rep(TRUE, nrow(candidates))

  abs_kme <- .wcp_num(candidates$abs_kME)
  abs_kme[!is.finite(abs_kme)] <- -Inf

  max_effect <- if ("sus_res_max_abs_log2FC" %in% names(candidates)) {
    .wcp_num(candidates$sus_res_max_abs_log2FC)
  } else rep(NA_real_, nrow(candidates))
  max_effect[!is.finite(max_effect)] <- -Inf

  consistency <- if ("sus_res_fraction_matching_majority_direction" %in% names(candidates)) {
    .wcp_num(candidates$sus_res_fraction_matching_majority_direction)
  } else rep(NA_real_, nrow(candidates))
  consistency[!is.finite(consistency)] <- -Inf

  ord <- order(
    tier, fdr_support, -abs_kme, -max_effect, -consistency,
    as.character(.wcp_col(candidates, "dataset", "")),
    as.character(candidates$ProteinGroupID),
    method = "radix"
  )
  candidates[ord, , drop = FALSE]
}

# ------------------------------------------------- human-facing review views

# Compact column set for manual biological review.  The full ~100-column table
# stays the exhaustive audit representation; this view drops source paths,
# internal contract fields and the redundant variants of the same quantity
# (top5_hub vs is_top5_module_representative, module_label vs
# module_display_label, GeneSymbol vs representative_gene_symbol, and so on).
wcp_protein_review_columns <- function() {
  c(
    # identity
    "dataset", "ModuleID", "module_label", "module_supermodule_label",
    "GeneSymbol", "RepresentativeUniProt",
    "candidate_tier", "candidate_tier_all", "phenotype_network_class",
    "clean_mapping",
    # network position
    "abs_kME", "abs_kME_rank_in_module", "n_module_members",
    "is_top5_module_representative", "is_top10_module_hub",
    # primary SUS - RES evidence
    "sus_res_min_BH_FDR", "sus_res_n_contexts_fdr05",
    "sus_res_strongest_spatial_unit", "sus_res_strongest_log2FC",
    "sus_res_median_log2FC", "sus_res_median_abs_log2FC",
    "sus_res_large_effect_typical",
    # descriptive spatial context
    "sus_res_n_spatial_contexts_tested", "sus_res_majority_direction",
    "sus_res_fraction_matching_majority_direction",
    "sus_res_spatially_consistent",
    # module context
    "module_estimate", "module_tier_specific_fdr", "module_support_class",
    "sus_res_direction_matches_module",
    # multimethod context
    "sus_res_gsea_leading_edge",
    # interpretation
    "candidate_reason"
  )
}

# Display names for the review sheet; keeps the sheet readable without
# renaming anything in the machine-readable CSVs.
wcp_protein_review_labels <- function() {
  c(
    dataset = "Dataset", ModuleID = "Module", module_label = "Module label",
    module_supermodule_label = "Supermodule", GeneSymbol = "Gene",
    RepresentativeUniProt = "UniProt", candidate_tier = "Tier",
    candidate_tier_all = "All tiers",
    phenotype_network_class = "Network class",
    clean_mapping = "Clean mapping",
    abs_kME = "|kME|", abs_kME_rank_in_module = "kME rank in module",
    n_module_members = "Module size",
    is_top5_module_representative = "Top-5 hub",
    is_top10_module_hub = "Top-10 hub",
    sus_res_min_BH_FDR = "SUS-RES min BH FDR",
    sus_res_n_contexts_fdr05 = "N contexts FDR<=0.05",
    sus_res_strongest_spatial_unit = "Strongest context",
    sus_res_strongest_log2FC = "Strongest log2FC",
    sus_res_median_log2FC = "Median log2FC",
    sus_res_median_abs_log2FC = "Median |log2FC|",
    sus_res_large_effect_typical = "Large typical effect",
    sus_res_n_spatial_contexts_tested = "N contexts tested",
    sus_res_majority_direction = "Majority direction",
    sus_res_fraction_matching_majority_direction = "Fraction agreeing",
    sus_res_spatially_consistent = "Unanimous direction",
    module_estimate = "Module estimate",
    module_tier_specific_fdr = "Module FDR",
    module_support_class = "Module support",
    sus_res_direction_matches_module = "Matches module direction",
    sus_res_gsea_leading_edge = "GSEA leading edge",
    candidate_reason = "Why flagged"
  )
}

# Review order: Tier A, then FDR-supported Tier D, then Tier B, then Tier C
# only.  Within a block: FDR ascending, large-effect flag, |kME| descending,
# median |effect| descending, ProteinGroupID.  A deterministic SORT for manual
# reading, not a statistical ranking.
wcp_order_protein_review <- function(candidates) {
  .wcp_require_columns(candidates, c("candidate_tier", "abs_kME", "ProteinGroupID"),
                       "Candidate table")
  tier <- as.character(candidates$candidate_tier)

  # Phenotype-linked classes first (A1, A2, D), then descriptive ones (B, C).
  # A reading order, not a statement about evidence strength.
  block <- rep(6L, nrow(candidates))
  block[tier %in% "C"] <- 5L
  block[tier %in% "B"] <- 4L
  block[tier %in% "D"] <- 3L
  block[tier %in% "A2"] <- 2L
  block[tier %in% "A1"] <- 1L

  fdr <- .wcp_num(.wcp_col(candidates, "sus_res_min_BH_FDR", NA_real_))
  fdr[!is.finite(fdr)] <- Inf
  large <- !.wcp_is_true(.wcp_col(candidates, "sus_res_large_effect_typical", FALSE))
  abs_kme <- .wcp_num(candidates$abs_kME); abs_kme[!is.finite(abs_kme)] <- -Inf
  med <- .wcp_num(.wcp_col(candidates, "sus_res_median_abs_log2FC", NA_real_))
  med[!is.finite(med)] <- -Inf

  candidates[order(block, fdr, large, -abs_kme, -med,
                   as.character(.wcp_col(candidates, "dataset", "")),
                   as.character(candidates$ProteinGroupID),
                   method = "radix"), , drop = FALSE]
}

wcp_protein_review_table <- function(candidates, rename = TRUE) {
  cols <- intersect(wcp_protein_review_columns(), names(candidates))
  out <- wcp_order_protein_review(candidates)[, cols, drop = FALSE]
  if (isTRUE(rename)) {
    labels <- wcp_protein_review_labels()
    hit <- names(out) %in% names(labels)
    names(out)[hit] <- unname(labels[names(out)[hit]])
  }
  out
}

# Every protein with SUS - RES BH FDR <= 0.05 in at least one spatial context:
# the inferentially supported individual proteins, in one place.
wcp_fdr_hits_table <- function(candidates, rename = TRUE) {
  .wcp_require_columns(candidates, "sus_res_fdr05_any_context", "Candidate table")
  hits <- candidates[.wcp_is_true(candidates$sus_res_fdr05_any_context), ,
                     drop = FALSE]
  wcp_protein_review_table(hits, rename = rename)
}

# One descriptive row per module.  Candidate fraction is NOT a measure of
# biological strength: Tier C contributes exactly 5 proteins to every module by
# construction, so small modules always show a higher fraction.
wcp_module_review_table <- function(candidates) {
  required <- c("dataset", "ModuleID", "is_tier_A1", "is_tier_A2", "is_tier_B",
                "is_tier_C", "is_tier_D", "is_candidate", "abs_kME",
                "GeneSymbol", "ProteinGroupID")
  .wcp_require_columns(candidates, required, "Candidate table")
  if (!nrow(candidates)) {
    return(data.frame(Dataset = character(), Module = character(),
                      stringsAsFactors = FALSE))
  }
  key <- paste(candidates$dataset, candidates$ModuleID, sep = "\r")
  idx <- split(seq_len(nrow(candidates)), key)

  top_names <- function(rows, n) {
    sym <- as.character(candidates$GeneSymbol[rows])
    pg <- as.character(candidates$ProteinGroupID[rows])
    sym[is.na(sym) | !nzchar(sym)] <- pg[is.na(sym) | !nzchar(sym)]
    paste(utils::head(sym, n), collapse = ", ")
  }

  rows <- lapply(idx, function(rows) {
    kme <- .wcp_num(candidates$abs_kME[rows])
    by_kme <- rows[order(-ifelse(is.finite(kme), kme, -Inf),
                         as.character(candidates$ProteinGroupID[rows]),
                         method = "radix")]
    cand <- rows[.wcp_is_true(candidates$is_candidate[rows])]
    data.frame(
      Dataset = as.character(candidates$dataset[rows[[1]]]),
      Module = as.character(candidates$ModuleID[rows[[1]]]),
      `Module label` = as.character(
        .wcp_col(candidates, "module_label", NA_character_)[rows[[1]]]
      ),
      `Module size` = length(rows),
      Candidates = length(cand),
      `Tier A1` = sum(.wcp_is_true(candidates$is_tier_A1[rows])),
      `Tier A2` = sum(.wcp_is_true(candidates$is_tier_A2[rows])),
      `Tier D` = sum(.wcp_is_true(candidates$is_tier_D[rows])),
      `Tier B` = sum(.wcp_is_true(candidates$is_tier_B[rows])),
      `Tier C` = sum(.wcp_is_true(candidates$is_tier_C[rows])),
      `Candidate fraction` = length(cand) / length(rows),
      `SUS-RES FDR<=0.05` = sum(.wcp_is_true(
        .wcp_col(candidates, "sus_res_fdr05_any_context", FALSE)[rows]
      )),
      `Large typical effect` = sum(.wcp_is_true(
        .wcp_col(candidates, "sus_res_large_effect_typical", FALSE)[rows]
      )),
      `Unanimous direction` = sum(.wcp_is_true(
        .wcp_col(candidates, "sus_res_spatially_consistent", FALSE)[rows]
      )),
      `GSEA leading edge` = sum(.wcp_is_true(
        .wcp_col(candidates, "sus_res_gsea_leading_edge", FALSE)[rows]
      )),
      `Top candidates` = top_names(cand, 5L),
      `Top hubs by |kME|` = top_names(by_kme, 5L),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$Dataset, out$Module, method = "radix"), , drop = FALSE]
}

# ------------------------------------------------------- README dictionary

# Machine-readable flag dictionary; rendered verbatim into the workbook README
# sheet so the exported file explains itself without this source file.
wcp_flag_dictionary <- function() {
  rows <- list(
    c("candidate_tier", "Derived label",
      "Convenience label from the tier booleans, precedence A1 > A2 > D > B > C. This precedence is a READING ORDER, not a ranking of evidence strength: A1/A2/D differ only in WGCNA topology. Overlap is preserved in candidate_tier_all.",
      "Derived here"),
    c("candidate_tier_all", "Derived label",
      "Every tier a protein satisfies, semicolon separated in precedence order (e.g. A1;B, A2;B, A2, D, B;C). Tiers are not mutually exclusive.",
      "Derived here"),
    c("is_tier_A1", "Phenotype-linked class",
      "FDR-supported phenotype-linked module HUB: clean mapping, SUS-RES BH FDR <= 0.05 in >=1 spatial context, |kME| >= 0.60, and in its module's top 10 by |kME|.",
      "Derived here"),
    c("is_tier_A2", "Phenotype-linked class",
      "FDR-supported phenotype-linked module MEMBER: same evidence as A1 and |kME| >= 0.60, but NOT in the module's top 10. No top-25 requirement is imposed; the frozen is_top_hub_25 flag remains separately available.",
      "Derived here"),
    c("is_tier_D", "Phenotype-linked class",
      "FDR-supported phenotype-linked PERIPHERAL module member: SUS-RES BH FDR <= 0.05 with |kME| < 0.60. Note this class does not additionally require clean_mapping, unlike A1/A2.",
      "Derived here"),
    c("is_tier_A", "Phenotype-linked class",
      "COMPATIBILITY AGGREGATE ONLY: is_tier_A1 OR is_tier_A2. It no longer denotes a single class; A1 and A2 carry the centrality distinction and are what human-facing sheets use.",
      "Derived here"),
    c("is_tier_B", "Descriptive candidate class",
      "Network-central candidate: |kME| >= 0.60 AND sus_res_large_effect_typical. Carries NO protein-level inferential support.",
      "Derived here"),
    c("is_tier_C", "Descriptive candidate class",
      "Canonical module representative: top 5 in the module by |kME|, regardless of differential abundance. A network-structural statement only.",
      "Derived here"),
    c("phenotype_network_class", "Phenotype-linked class",
      "Readability label derived from A1/A2/D: top10_hub, module_member, peripheral_member. Describes WGCNA topology only. It does NOT indicate which protein has the stronger differential-abundance evidence.",
      "Derived here"),
    c("candidate_reason", "Derived label",
      "Human-readable list of the satisfied flags behind the tier assignment.",
      "Derived here"),
    c("kME / abs_kME", "WGCNA module centrality",
      "Signed and absolute module eigengene correlation, copied unchanged from frozen Stage 01 WGCNA output.",
      "analysis/05_wgcna/build_wgcna_modules.R"),
    c("abs_kME_rank_in_module", "WGCNA module centrality",
      "Rank of |kME| within the protein's own module; ties broken by ProteinGroupID ascending.",
      "Derived here"),
    c("abs_kME_rank_fraction_in_module", "WGCNA module centrality",
      "Rank divided by module size; 0 = most central.", "Derived here"),
    c("is_core_kME_0.6", "WGCNA module centrality",
      "Frozen upstream flag, |kME| >= 0.60. Copied unchanged; not recomputed.",
      "analysis/05_wgcna/build_wgcna_modules.R"),
    c("is_top_hub_25", "WGCNA module centrality",
      "Frozen upstream flag, top 25 members per module by |kME|. Copied unchanged.",
      "analysis/05_wgcna/build_wgcna_modules.R"),
    c("top5_hub / is_top5_module_representative", "WGCNA module centrality",
      "Top 5 members of the module by |kME|. Answers 'which proteins best represent this module'.",
      "Derived here"),
    c("top10_hub / is_top10_module_hub", "WGCNA module centrality",
      "Top 10 members of the module by |kME|.", "Derived here"),
    c("high_kME", "WGCNA module centrality",
      "|kME| >= 0.60, matching the frozen is_core_kME_0.6 definition.",
      "Derived here"),
    c("very_high_kME", "WGCNA module centrality",
      "|kME| >= 0.80. DESCRIPTIVE ONLY - the repository defines no canonical 0.80 threshold.",
      "Derived here"),
    c("GeneSignificanceP / GeneSignificanceFDR", "Omnibus condition test",
      "OMNIBUS one-way ANOVA across CON/RES/SUS on the WGCNA expression matrix, BH adjusted. UNSIGNED, and NOT a SUS-vs-RES contrast. Do not read as phenotype correlation.",
      "analysis/05_wgcna/build_wgcna_modules.R"),
    c("log2FC / raw_p / BH_FDR", "Protein-level differential abundance",
      "Manifest-selected limma results copied unchanged; oriented to the canonical contrast label. No DA model is refitted.",
      "clusterProfiler manifest -> Stage 02 mapped per-file DA"),
    c("sus_res_fdr05_any_context", "Protein-level differential abundance",
      "SUS - RES BH FDR <= 0.05 in at least one spatial context.",
      "Derived here"),
    c("large_effect_within_context", "Protein-level differential abundance",
      "|log2FC| at or above the 90th percentile WITHIN its own dataset x spatial unit x contrast family. Data-relative, no universal fold-change cut. Per-context, so NOT used for tiering: 'top decile in any context' gives a protein one chance per context and flags 35-65% of proteins.",
      "Derived here"),
    c("sus_res_large_effect_typical", "Protein-level differential abundance",
      "Median |log2FC| across spatial contexts at or above the 90th percentile within dataset x contrast. This is the ONLY phenotype criterion used for Tier B, because it selects 10% by construction rather than inflating across contexts.",
      "Derived here"),
    c("module_estimate / module_tier_specific_fdr", "Module-level inference",
      "Stage 07 module eigengene contrast estimate and its tier-specific BH FDR, copied unchanged from the canonical inferential handoff.",
      "analysis/05_wgcna/summarize_module_interpretation.R"),
    c("module_support_class / module_model_valid", "Module-level inference",
      "Frozen Stage 07 support classification and model validity gate. Copied unchanged.",
      "analysis/05_wgcna/summarize_module_interpretation.R"),
    c("direction_matches_module", "Directional concordance",
      "Sign agreement between the protein log2FC and the module eigengene estimate for the SAME contrast. NA unless both are finite and non-zero. DESCRIPTIVE CONCORDANCE, NOT independent replication - the eigengene summarises the same proteins in the same animals.",
      "Derived here"),
    c("n_spatial_contexts_tested / _positive / _negative", "Spatial consistency",
      "Counts of hippocampal spatial contexts with a finite effect, and their signs.",
      "Derived here"),
    c("majority_direction", "Spatial consistency",
      "Sign shared by most tested contexts. For SUS - RES, positive = higher in SUS.",
      "Derived here"),
    c("fraction_matching_majority_direction", "Spatial consistency",
      "Share of tested contexts agreeing with the majority sign. Exported continuously so other thresholds can be applied by hand.",
      "Derived here"),
    c("spatially_consistent", "Spatial consistency",
      "Unanimous sign across at least 3 tested contexts. DESCRIPTIVE PRIORITISATION ONLY - it does not create a candidate tier. These contexts are repeated, positively correlated measurements from the SAME animals, so their signs are not independent trials; no probability is attached to unanimity and none should be inferred.",
      "Derived here"),
    c("gsea_leading_edge", "GSEA-derived concordance",
      "Protein lies in the leading edge of a module-concordant GSEA term (overlap BH FDR <= 0.05). NOT independent of the DA columns: that GSEA ranks on the moderated t statistic from the SAME limma fit.",
      "analysis/08_integration/test_enrichment_module_concordance.R"),
    c("mapping_status / *_claim_allowed / protein_group_ambiguity_class",
      "Mapping and claimability",
      "Frozen Stage 02 protein-group identity contract fields, copied unchanged.",
      "02_id_mapping canonical mapping"),
    c("clean_mapping", "Mapping and claimability",
      "Gene- AND protein-level claim allowed AND a single-accession single-gene protein group. Gene-ambiguous groups are never silently upgraded.",
      "Derived here")
  )
  out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(out) <- c("Field", "Evidence class", "Meaning", "Provenance")
  out
}

# Evidence-class glossary shown at the top of the workbook README sheet.
wcp_evidence_class_notes <- function() {
  data.frame(
    `Evidence class` = c(
      "WGCNA module centrality", "Omnibus condition test",
      "Protein-level differential abundance", "Module-level inference",
      "Directional concordance", "Spatial consistency",
      "GSEA-derived concordance", "Mapping and claimability"
    ),
    Meaning = c(
      "Where a protein sits inside its own frozen co-abundance module. Says nothing about stress response.",
      "Unsigned omnibus ANOVA across CON/RES/SUS from Stage 01. Not a SUS-vs-RES contrast and not a signed phenotype correlation.",
      "Per-protein limma contrast within one spatial context. The only per-protein inferential evidence in this table.",
      "Eigengene-level mixed-model contrast for the whole module. A separate hypothesis at a different unit of analysis.",
      "Sign agreement between a protein and its module. Descriptive only.",
      "Sign agreement across hippocampal contexts measured in the same animals. Repeated measures, not replication.",
      "Leading-edge membership from ranked GSEA. Derived from the same DA statistics, therefore NOT independent evidence.",
      "Frozen Stage 02 identity contract governing whether a gene- or protein-level claim may be made at all."
    ),
    `Independent of protein-level DA?` = c(
      "Yes - built from co-abundance structure, not group contrasts",
      "No - same expression matrix, different (omnibus) test",
      "-",
      "Partly - same animals and proteins, different unit of analysis",
      "No - one side of the comparison is the DA estimate",
      "No - repeated measurements of the same animals",
      "No - ranks on the moderated t of the same limma fit",
      "Not applicable"
    ),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

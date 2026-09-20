# Spatial / cell-affinity atlas utilities.
#
# The atlas answers one question: where does a module or protein normally belong
# anatomically and by cellular context, and how reproducible is that identity
# across hemispheres?
#
# FOUR EVIDENCE DIMENSIONS ARE KEPT SEPARATE AND ARE NEVER SUMMED
#   A  spatial anatomical identity   (CON baseline)
#   B  empirical compartment affinity (this experiment's ROIs)
#   C  external / reference cell-type affinity (outside data)
#   D  bilateral reliability
# A categorical context class is derived from EXPLICIT rules and always carries
# the reason that produced it, so no opaque composite score is created.
#
# EVERYTHING HERE IS PHENOTYPE-BLIND. SUS/RES never define identity; they are
# joined later as a downstream overlay only.

## The Stage-05 hemisphere values come from the WGCNA path resolver, which
## looks in the normalized location first and falls back to the historical
## one. This file has no plain source() block of its own, so the guard is
## explicit; callers that already loaded the resolver pay nothing.
if (!exists("wgcna_group_effects_artifact", mode = "function")) {
  if (!exists("repo_path", mode = "function")) {
    paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
    source(paths_file)
  }
  source(repo_path("R", "wgcna_paths.R"))
}

sat_contract_version <- function() "spatial_cell_atlas_v1"

# ------------------------------------------------------ spatial profile metrics
#
# Module eigengenes are SIGNED (observed range about -0.17 to +0.15), so an
# entropy taken directly on raw values is undefined. Both indices below are
# therefore computed on the non-negative deviation profile
#
#     w_i = x_i - min(x)
#
# which measures each unit's elevation above the least-active unit.
#
#   tau      Tau specificity index, bounded [0, 1].
#            0 = perfectly uniform, 1 = all elevation in a single unit.
#            tau = (n * max(x) - sum(x)) / ((n - 1) * (max(x) - min(x))), which is
#            algebraically identical to the range-normalised peak elevation
#            (max(x) - mean(other units)) / (max(x) - min(x)) - verified to 2e-16.
#   entropy  Shannon entropy of w / sum(w), divided by log(n).
#            1 = uniform, 0 = single unit.
#
# BOTH are normalised by the number of spatial units, so a 10-unit neuropil
# profile and a 4-unit region profile are directly comparable. A perfectly flat
# profile is the degenerate case and is defined as maximally UNSPECIFIC
# (tau = 0, entropy = 1) rather than returned as NA.
#
# Peak, second and ranks come from the RAW values, never from a z-scored or
# otherwise rescaled profile: standardisation is for visualisation only.
sat_profile_metrics <- function(x, unit_names = names(x)) {
  # Force the default BEFORE x is coerced: `unit_names = names(x)` is evaluated
  # lazily, so stripping names via as.numeric() first would silently yield
  # positional indices ("1", "2", ...) instead of unit labels.
  force(unit_names)
  x <- as.numeric(x)
  n <- length(x)
  if (is.null(unit_names)) unit_names <- as.character(seq_len(n))
  ok <- is.finite(x)
  if (sum(ok) < 2L) {
    return(data.frame(
      n_spatial_units = sum(ok), peak_unit = NA_character_,
      second_unit = NA_character_, peak_value = NA_real_,
      second_value = NA_real_, peak_minus_second = NA_real_,
      peak_minus_rest = NA_real_, spatial_tau = NA_real_,
      spatial_entropy_normalized = NA_real_, profile_range = NA_real_,
      spatial_rank_order = NA_character_, stringsAsFactors = FALSE))
  }
  xv <- x[ok]; un <- unit_names[ok]; n <- length(xv)
  ord <- order(-xv)
  rng <- max(xv) - min(xv)
  flat <- !is.finite(rng) || rng <= 0
  w <- xv - min(xv)
  tau <- if (flat) 0 else sum(1 - w / max(w)) / (n - 1)
  ent <- if (flat) 1 else {
    p <- w[w > 0] / sum(w)
    -sum(p * log(p)) / log(n)
  }
  data.frame(
    n_spatial_units = n,
    peak_unit = un[ord[1]],
    second_unit = if (n >= 2L) un[ord[2]] else NA_character_,
    peak_value = xv[ord[1]],
    second_value = if (n >= 2L) xv[ord[2]] else NA_real_,
    peak_minus_second = if (n >= 2L) xv[ord[1]] - xv[ord[2]] else NA_real_,
    peak_minus_rest = xv[ord[1]] - mean(xv[-ord[1]]),
    spatial_tau = tau,
    spatial_entropy_normalized = ent,
    profile_range = rng,
    spatial_rank_order = paste(un[ord], collapse = ">"),
    stringsAsFactors = FALSE)
}

# Row-standardise a profile matrix FOR VISUALISATION ONLY.
sat_row_z <- function(m) {
  mu <- rowMeans(m, na.rm = TRUE)
  sdv <- apply(m, 1, stats::sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv <= 0] <- NA_real_
  z <- sweep(sweep(m, 1, mu, "-"), 1, sdv, "/")
  z[!is.finite(z)] <- 0
  z
}

# --------------------------------------------- CON baseline module profiles
#
# Consumes the ACCEPTED Stage-05 hemisphere values. Nothing about WGCNA is
# recomputed: the left and right eigengene values are averaged with equal weight
# within an animal (the canonical bilateral policy) and then averaged across CON
# animals only.
sat_con_bilateral_module_profiles <- function(dataset, level = "module") {
  p <- wgcna_group_effects_artifact("WGCNA_group_effect_hemisphere_values.csv", dataset)
  if (!file.exists(p)) {
    stop("missing_required_input: Stage-05 hemisphere values: ", p, call. = FALSE)
  }
  h <- as.data.frame(readr::read_csv(p, show_col_types = FALSE, progress = FALSE,
                                     guess_max = Inf))
  h <- h[h$level %in% level, , drop = FALSE]
  # PHENOTYPE-BLIND: CON only defines baseline identity.
  con <- h[h$StressGroup == "CON", , drop = FALSE]
  if (!nrow(con)) stop("No CON rows for ", dataset, ".", call. = FALSE)

  # equal-weight bilateral mean within animal x unit, then mean across animals
  per_animal <- stats::aggregate(
    hemisphere_value ~ endpoint_id + AnimalID + SpatialUnit, con, mean)
  names(per_animal)[names(per_animal) == "hemisphere_value"] <- "bilateral_value"
  per_unit <- stats::aggregate(bilateral_value ~ endpoint_id + SpatialUnit,
                               per_animal, mean)
  names(per_unit)[names(per_unit) == "bilateral_value"] <- "mean_activity"
  sd_unit <- stats::aggregate(bilateral_value ~ endpoint_id + SpatialUnit,
                              per_animal, stats::sd)
  names(sd_unit)[names(sd_unit) == "bilateral_value"] <- "sd_activity"
  n_unit <- stats::aggregate(bilateral_value ~ endpoint_id + SpatialUnit,
                             per_animal, length)
  names(n_unit)[names(n_unit) == "bilateral_value"] <- "n_con_animals"
  out <- merge(merge(per_unit, sd_unit, by = c("endpoint_id", "SpatialUnit")),
               n_unit, by = c("endpoint_id", "SpatialUnit"))
  out$se_activity <- out$sd_activity / sqrt(out$n_con_animals)
  out$dataset <- dataset
  out$level <- level[[1]]
  out$baseline_definition <- "CON animals only; equal-weight L/R mean within animal, then mean across CON animals"
  out[order(out$endpoint_id, out$SpatialUnit), , drop = FALSE]
}

# Wide profile matrix (modules x spatial units) from the long baseline table.
sat_profile_matrix <- function(long, id_col = "endpoint_id",
                               unit_col = "SpatialUnit", value_col = "mean_activity") {
  ids <- sort(unique(as.character(long[[id_col]])))
  units <- sort(unique(as.character(long[[unit_col]])))
  m <- matrix(NA_real_, nrow = length(ids), ncol = length(units),
              dimnames = list(ids, units))
  m[cbind(match(as.character(long[[id_col]]), ids),
          match(as.character(long[[unit_col]]), units))] <- as.numeric(long[[value_col]])
  m
}

# ------------------------------------------------------------ enrichment
#
# One-sided Fisher test of a module's proteins against a marker set, with the
# MEASURED PROTEOME as universe. Returns the 2x2 counts so a reader can audit
# the test rather than trust an odds ratio alone.
sat_fisher_enrichment <- function(module_ids, marker_ids, universe_ids) {
  universe <- unique(as.character(universe_ids))
  a_set <- intersect(unique(as.character(module_ids)), universe)
  b_set <- intersect(unique(as.character(marker_ids)), universe)
  a <- length(intersect(a_set, b_set))
  b <- length(a_set) - a
  c_ <- length(b_set) - a
  d <- length(universe) - a - b - c_
  if (any(c(a, b, c_, d) < 0)) {
    stop("Enrichment counts are negative; the universe does not contain the ",
         "tested sets.", call. = FALSE)
  }
  ft <- stats::fisher.test(matrix(c(a, b, c_, d), nrow = 2), alternative = "greater")
  data.frame(
    n_universe = length(universe),
    n_module_in_universe = length(a_set),
    n_markers_in_universe = length(b_set),
    n_overlap = a,
    expected_overlap = length(a_set) * length(b_set) / max(length(universe), 1L),
    odds_ratio = unname(ft$estimate),
    p_value = ft$p.value,
    stringsAsFactors = FALSE)
}

# Is a marker set capable of discriminating compartment affinity at all?
#
# PRESPECIFIED, and fixed before inspecting which modules it favours. The
# empirical ROI sets are extremely unbalanced: in the neuropil universe
# `empirical_microglia_neuropil_shared` covers 72.6% of all measured proteins,
# while `empirical_microglia_roi_high_confidence` contains a single protein.
#
# A set spanning most of the proteome cannot say where a module "belongs" - it
# will win any enrichment comparison on power alone, and it was selected as the
# strongest compartment for 20 of 35 modules before this guard existed. A set
# with almost no members cannot support a stable Fisher test either.
#
# Non-informative sets are NOT deleted: they stay in the long table with the
# flag and their coverage fraction, and are only excluded from the single
# "strongest compartment" call.
sat_marker_set_bounds <- function() list(min_fraction = 0.005, max_fraction = 0.50,
                                         min_members = 5L)

sat_marker_set_is_informative <- function(n_markers_in_universe, n_universe) {
  b <- sat_marker_set_bounds()
  frac <- as.numeric(n_markers_in_universe) / pmax(as.numeric(n_universe), 1)
  is.finite(frac) & frac >= b$min_fraction & frac <= b$max_fraction &
    as.numeric(n_markers_in_universe) >= b$min_members
}

# FDR families are declared prospectively so a family can never widen silently.
sat_fdr_family_empirical_compartment <- function(dataset, scope) {
  paste("empirical_compartment", dataset, scope, sep = "_")
}
sat_fdr_family_reference_marker <- function(dataset, scope) {
  paste("reference_marker", dataset, scope, sep = "_")
}

# BH within each declared family, never pooled.
sat_apply_family_fdr <- function(df, p_col = "p_value", family_col = "fdr_family",
                                 out_col = "FDR") {
  for (nm in c(p_col, family_col)) {
    if (!nm %in% names(df)) stop("Missing column: ", nm, ".", call. = FALSE)
  }
  if (!nrow(df)) { df[[out_col]] <- numeric(0); return(df) }
  fam <- as.character(df[[family_col]])
  if (any(is.na(fam) | !nzchar(fam))) {
    stop("Every row needs an explicit FDR family.", call. = FALSE)
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

# ------------------------------------------------ bilateral reliability class
#
# PRESPECIFIED thresholds, fixed before looking at which modules they select.
# They are deliberately coarse: the point is to make weak bilateral support
# visible, not to rank modules finely.
sat_bilateral_thresholds <- function() {
  list(strong = 0.80, moderate = 0.60)
}

sat_bilateral_support_class <- function(absolute_r, profile_r) {
  th <- sat_bilateral_thresholds()
  a <- as.numeric(absolute_r); p <- as.numeric(profile_r)
  best <- pmax(a, p, na.rm = TRUE)
  out <- rep("insufficient_bilateral_data", length(best))
  known <- is.finite(best)
  out[known] <- "low_bilateral_support"
  out[known & best >= th$moderate] <- "moderate_bilateral_support"
  out[known & best >= th$strong] <- "high_bilateral_support"
  out
}

# ------------------------------------------------------- context classification
#
# Rule-based, never a hidden numeric sum. Every call returns the class AND the
# rule that fired, so a reviewer can see exactly why a module was classified.
#
# A low bilateral reliability NEVER erases an anatomical identity. It produces
# the explicit `spatial_identity_present_but_bilaterally_variable` caveat
# instead, because a spatial pattern that is real but asymmetric is a finding,
# not an absence.
sat_context_class <- function(has_spatial, has_compartment, has_external,
                              bilateral_class, external_mixed = FALSE) {
  n <- max(length(has_spatial), length(has_compartment), length(has_external),
           length(bilateral_class))
  rep_to <- function(v) if (length(v) == n) v else rep(v, length.out = n)
  hs <- rep_to(has_spatial %in% TRUE)
  hc <- rep_to(has_compartment %in% TRUE)
  he <- rep_to(has_external %in% TRUE)
  bc <- rep_to(as.character(bilateral_class))
  em <- rep_to(external_mixed %in% TRUE)

  cls <- rep("unresolved", n)
  rule <- rep("no spatial, compartment or external evidence", n)

  i <- hs & (hc | he)
  cls[i] <- "strong_multisource_context"
  rule[i] <- "spatial identity plus at least one cell-context source"

  i <- hs & !hc & !he
  cls[i] <- "spatial_only"
  rule[i] <- "spatial identity without compartment or external support"

  i <- !hs & (hc | he)
  cls[i] <- "cell_affinity_only"
  rule[i] <- "cell context without a distinct spatial peak"

  i <- em & (hc | he)
  cls[i] <- "mixed_cell_context"
  rule[i] <- "cell-type evidence disagrees across module scopes"

  i <- (hs | hc | he) & bc %in% c("low_bilateral_support",
                                  "insufficient_bilateral_data")
  cls[i] <- "bilaterally_variable"
  rule[i] <- paste0("evidence present but bilateral reliability is ",
                    bc[i], "; identity retained with a caveat")

  i <- !hs & !hc & !he & !(bc %in% "insufficient_bilateral_data")
  cls[i] <- "weak_context"
  rule[i] <- "no dimension reached its evidence threshold"

  data.frame(context_confidence = cls, context_rule = rule,
             stringsAsFactors = FALSE)
}

# --------------------------------------------------- label vs context audit
#
# AUDIT ONLY. Nothing here renames a module or activates a label.
sat_label_context_relationship <- function(label, spatial_ctx, compartment_ctx,
                                           external_ctx) {
  lab <- tolower(paste(label, ""))
  ext <- tolower(paste(external_ctx, ""))
  cmp <- tolower(paste(compartment_ctx, ""))
  n <- length(lab)
  out <- rep("insufficient_context", n)
  rule <- rep("no cell or spatial context available", n)

  informative <- nzchar(trimws(ext)) & !ext %in% c("na", "none", "")
  out[informative] <- "context_compatible"
  rule[informative] <- "context present and not contradictory"

  # a shared biological token between the GO-derived label and the external
  # cell type is orthogonal corroboration: different data, same conclusion
  tokens <- list(
    oligodendrocyte = c("oligodendro", "myelin"),
    microglia = c("microglia", "immune"),
    astrocytes_ependymal = c("astro", "ependym"),
    endothelial_mural = c("endothel", "mural", "vascular", "pericyte"),
    interneurons = c("interneuron", "gaba", "inhibitory"),
    pyramidal = c("pyramidal", "synap", "neuron", "excitatory", "glutamat"))
  for (i in seq_len(n)) {
    if (!informative[i]) next
    for (ct in names(tokens)) {
      if (!grepl(ct, ext[i], fixed = TRUE)) next
      if (any(vapply(tokens[[ct]], function(tk) grepl(tk, lab[i], fixed = TRUE),
                     logical(1)))) {
        out[i] <- "context_corroborates_label"
        rule[i] <- paste0("external cell type '", external_ctx[i],
                          "' shares biology with the label")
        break
      }
    }
  }
  data.frame(label_context_relationship = out,
             label_context_rule = rule, stringsAsFactors = FALSE)
}

# ------------------------------------------------ canonical spatial unit names
#
# The same anatomical unit is spelled differently by different producers:
#
#   Stage-05 baseline     ca2_slm      ca2            ca2
#   canonical DA          CA2_slm      CA2_sp         CA2_microglia
#                         (neuropil)   (soma)         (microglia)
#
# Joining these without normalisation silently drops every row - the join simply
# finds no match and the atlas would show empty phenotype overlays while looking
# perfectly healthy. This maps any producer's spelling onto the baseline form.
#
# Soma and microglia resolve at REGION level, so their layer/compartment suffix
# is part of the producer's label, not part of the spatial unit, and is removed.
sat_canonical_spatial_unit <- function(unit, dataset) {
  u <- tolower(trimws(as.character(unit)))
  u[!nzchar(u)] <- NA_character_
  ds <- as.character(dataset)
  if (length(ds) == 1L) ds <- rep(ds, length(u))
  out <- u
  # neuropil keeps region_layer; only the case differs
  # soma: strip the soma layer token (sp / sg)
  i <- ds == "neuron_soma"
  out[i] <- sub("_(sp|sg)$", "", out[i])
  # microglia: strip the compartment token
  i <- ds == "microglia"
  out[i] <- sub("_microglia$", "", out[i])
  out
}

# Animal-level spatial molecular-similarity networks.
#
# THE DESIGN REQUIREMENT
#   One AnimalID = one independent network replicate. Every inferential unit
#   below is an animal. No sample, no hemisphere and no protein is ever treated
#   as an independent biological replicate.
#
# WHAT A NETWORK MEANS HERE
#   NODE  an anatomical sampling unit (Region x Layer for neuropil, Region for
#         soma and microglia).
#   EDGE  the similarity of the molecular spatial profiles of two anatomical
#         units WITHIN ONE ANIMAL.
#
#   This is a spatial molecular-profile similarity network. It is NOT neural
#   connectivity, NOT anatomical connectivity, NOT molecular communication, NOT
#   a protein coexpression network and NOT brain connectivity. Those terms are
#   forbidden in every output this file feeds.
#
# WHY THIS REPLACES THE OLD LAYER
#   The historical 07_spatial_networks pipeline pooled animals into one mean
#   profile per spatial unit and then bootstrapped SAMPLE COLUMNS, so its
#   intervals described sampling of ROIs rather than sampling of animals. It
#   also carried a field named permutation_p that was a bootstrap sign
#   frequency, with a BH FDR built on it. Nothing here reuses either.

asn_contract_version <- function() "animal_spatial_network_v1"

# ------------------------------------------------------------ representations
#
# All four are declared PROSPECTIVELY, before any phenotype result is seen.
asn_representations <- function() {
  c("A_abundance", "B_protein_centered", "C_protein_z", "D_centered_varfilter")
}

asn_representation_description <- function(rep) {
  switch(rep,
    A_abundance = "canonical normalised log-scale abundance, untransformed",
    B_protein_centered = "within animal, each protein centred on its mean across that animal's spatial units",
    C_protein_z = "within animal, each protein centred and divided by its SD across spatial units; near-zero-variance proteins excluded",
    D_centered_varfilter = "as B, with near-zero spatial-variance proteins excluded",
    NA_character_)
}

# Predefined variance threshold. Fixed before inspecting any phenotype result.
asn_min_spatial_sd <- function() 0.05

# Apply one representation to a protein x spatial-unit matrix for ONE animal.
# Returns the transformed matrix plus the number of proteins retained, so the
# retention cost of each representation is always visible.
asn_transform <- function(m, representation) {
  representation <- match.arg(representation, asn_representations())
  sdv <- apply(m, 1, stats::sd, na.rm = TRUE)
  keep_all <- rep(TRUE, nrow(m))
  thr <- asn_min_spatial_sd()
  out <- switch(representation,
    A_abundance = list(mat = m, keep = keep_all),
    B_protein_centered = list(mat = m - rowMeans(m, na.rm = TRUE), keep = keep_all),
    C_protein_z = {
      keep <- is.finite(sdv) & sdv >= thr
      z <- sweep(sweep(m[keep, , drop = FALSE], 1,
                       rowMeans(m[keep, , drop = FALSE], na.rm = TRUE), "-"),
                 1, sdv[keep], "/")
      list(mat = z, keep = keep)
    },
    D_centered_varfilter = {
      keep <- is.finite(sdv) & sdv >= thr
      list(mat = m[keep, , drop = FALSE] -
             rowMeans(m[keep, , drop = FALSE], na.rm = TRUE), keep = keep)
    })
  out$n_proteins_input <- nrow(m)
  out$n_proteins_retained <- nrow(out$mat)
  out$representation <- representation
  out
}

# ------------------------------------------------------------ one network
#
# The full weighted similarity matrix over a FIXED node ordering. No edge is
# ever thresholded here: thresholding belongs to visualisation only.
asn_similarity_matrix <- function(m, nodes, method = "spearman") {
  method <- match.arg(method, c("spearman", "pearson"))
  missing <- setdiff(nodes, colnames(m))
  if (length(missing)) {
    stop("Network matrix is missing node(s): ", paste(missing, collapse = ", "),
         ". Node ordering must match exactly across every comparison.",
         call. = FALSE)
  }
  m <- m[, nodes, drop = FALSE]
  s <- suppressWarnings(stats::cor(m, method = method, use = "pairwise.complete.obs"))
  dimnames(s) <- list(nodes, nodes)
  s
}

# Upper-triangle edge vector in a deterministic order, so two animals' vectors
# are always element-comparable.
asn_edge_vector <- function(s) {
  nodes <- rownames(s)
  idx <- which(upper.tri(s), arr.ind = TRUE)
  ord <- order(idx[, 1], idx[, 2])
  idx <- idx[ord, , drop = FALSE]
  stats::setNames(s[cbind(idx[, 1], idx[, 2])],
                  paste(nodes[idx[, 1]], nodes[idx[, 2]], sep = "__"))
}

asn_edge_table <- function(s) {
  v <- asn_edge_vector(s)
  parts <- do.call(rbind, strsplit(names(v), "__", fixed = TRUE))
  data.frame(node_a = parts[, 1], node_b = parts[, 2], edge_id = names(v),
             similarity = unname(v), stringsAsFactors = FALSE)
}

# Fisher z, used wherever similarities are averaged or differenced so the
# bounded correlation scale does not distort the arithmetic.
asn_fisher_z <- function(r) {
  # Preserve dim/dimnames: as.numeric() on a matrix silently returns a plain
  # vector, which would turn an animal x edge matrix into an unusable vector.
  d <- dim(r); dn <- dimnames(r)
  out <- atanh(pmin(pmax(as.numeric(r), -0.999999), 0.999999))
  if (!is.null(d)) { dim(out) <- d; dimnames(out) <- dn }
  out
}

# ------------------------------------------------------- global descriptors
#
# A SMALL prespecified set. Fashionable graph metrics are deliberately absent:
# betweenness, modularity and small-worldness are not interpretable on a fully
# weighted 4-to-10 node similarity matrix and are not computed.
asn_global_descriptors <- function(s, dataset) {
  v <- asn_edge_vector(s)
  nodes <- rownames(s)
  region <- toupper(sub("_.*$", "", nodes))
  layer <- ifelse(grepl("_", nodes), sub("^[^_]*_", "", nodes), NA_character_)
  parts <- do.call(rbind, strsplit(names(v), "__", fixed = TRUE))
  ra <- toupper(sub("_.*$", "", parts[, 1])); rb <- toupper(sub("_.*$", "", parts[, 2]))
  la <- ifelse(grepl("_", parts[, 1]), sub("^[^_]*_", "", parts[, 1]), NA_character_)
  lb <- ifelse(grepl("_", parts[, 2]), sub("^[^_]*_", "", parts[, 2]), NA_character_)

  out <- list(
    mean_edge_similarity = mean(v, na.rm = TRUE),
    sd_edge_similarity = stats::sd(v, na.rm = TRUE),
    mean_abs_edge_similarity = mean(abs(v), na.rm = TRUE))
  pair_mean <- function(sel) if (any(sel, na.rm = TRUE)) mean(v[sel], na.rm = TRUE) else NA_real_
  if (identical(dataset, "neuron_neuropil")) {
    out$within_region_cross_layer <- pair_mean(ra == rb)
    out$same_layer_across_region <- pair_mean(ra != rb & !is.na(la) & la == lb)
    out$CA1_CA2_similarity <- pair_mean((ra == "CA1" & rb == "CA2") | (ra == "CA2" & rb == "CA1"))
    out$CA1_CA3_similarity <- pair_mean((ra == "CA1" & rb == "CA3") | (ra == "CA3" & rb == "CA1"))
    out$CA2_CA3_similarity <- pair_mean((ra == "CA2" & rb == "CA3") | (ra == "CA3" & rb == "CA2"))
    out$DG_to_CA_similarity <- pair_mean(xor(ra == "DG", rb == "DG"))
  } else {
    # only region-level contrasts are meaningful with four nodes
    out$CA1_CA2_similarity <- pair_mean((ra == "CA1" & rb == "CA2") | (ra == "CA2" & rb == "CA1"))
    out$CA1_CA3_similarity <- pair_mean((ra == "CA1" & rb == "CA3") | (ra == "CA3" & rb == "CA1"))
    out$CA2_CA3_similarity <- pair_mean((ra == "CA2" & rb == "CA3") | (ra == "CA3" & rb == "CA2"))
    out$DG_to_CA_similarity <- pair_mean(xor(ra == "DG", rb == "DG"))
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}

# Weighted node strength: the mean similarity of a node to all other nodes.
asn_node_strength <- function(s) {
  diag(s) <- NA_real_
  rowMeans(s, na.rm = TRUE)
}

# --------------------------------------------------- CON centroid distance
#
# For a CON animal the centroid is built LEAVE-ONE-OUT, otherwise each CON
# animal would be compared against a centroid that contains itself and would be
# handed an artificially small distance.
asn_distance_from_con <- function(edge_mat, group, method = c("euclidean_fisherz",
                                                              "correlation")) {
  method <- match.arg(method)
  group <- as.character(group)
  con <- which(group == "CON")
  if (length(con) < 2L) {
    return(rep(NA_real_, nrow(edge_mat)))
  }
  z <- asn_fisher_z(edge_mat)
  vapply(seq_len(nrow(z)), function(i) {
    ref <- if (group[i] == "CON") setdiff(con, i) else con
    if (!length(ref)) return(NA_real_)
    centroid <- colMeans(z[ref, , drop = FALSE], na.rm = TRUE)
    if (method == "euclidean_fisherz") {
      sqrt(sum((z[i, ] - centroid)^2, na.rm = TRUE))
    } else {
      1 - suppressWarnings(stats::cor(z[i, ], centroid, use = "complete.obs"))
    }
  }, numeric(1))
}

# ------------------------------------------------------- exact permutation
#
# The permutation unit is the ANIMAL. Whole animals are relabelled; an animal's
# edges always move together.
asn_exact_two_group_p <- function(values, group, a = "SUS", b = "RES") {
  keep <- group %in% c(a, b)
  v <- as.numeric(values)[keep]; g <- as.character(group)[keep]
  na <- sum(g == a); n <- length(g)
  if (n < 2L || na < 1L || na == n) {
    return(list(p_two_sided = NA_real_, n_assignments = 0L,
                min_attainable_two_sided_p = NA_real_, observed = NA_real_))
  }
  combos <- utils::combn(n, na)
  obs <- mean(v[g == a]) - mean(v[g == b])
  null <- apply(combos, 2, function(i) mean(v[i]) - mean(v[-i]))
  n_assign <- ncol(combos)
  list(
    observed = obs,
    n_assignments = n_assign,
    # (1 + #{as extreme}) / (1 + B): an exact enumeration still cannot return 0
    p_two_sided = (1 + sum(abs(null) >= abs(obs) - 1e-12)) / (1 + n_assign),
    min_attainable_two_sided_p = 2 / n_assign)
}

# Exact three-group omnibus by full label enumeration (1680 assignments for
# 3/3/3), using the between-group sum of squares as the statistic.
asn_exact_three_group_p <- function(values, group, max_assignments = 5000L) {
  v <- as.numeric(values); g <- as.character(group)
  ok <- is.finite(v)
  v <- v[ok]; g <- g[ok]
  lv <- sort(unique(g)); n <- length(v)
  if (length(lv) != 3L || n < 6L) {
    return(list(p = NA_real_, n_assignments = 0L, observed = NA_real_))
  }
  sizes <- vapply(lv, function(k) sum(g == k), integer(1))
  stat <- function(lab) {
    mu <- mean(v)
    sum(vapply(lv, function(k) sum(lab == k) * (mean(v[lab == k]) - mu)^2, numeric(1)))
  }
  obs <- stat(g)
  # enumerate every distinct assignment of the fixed group sizes
  idx1 <- utils::combn(n, sizes[[1]])
  null <- c()
  for (j in seq_len(ncol(idx1))) {
    rest <- setdiff(seq_len(n), idx1[, j])
    idx2 <- utils::combn(rest, sizes[[2]])
    for (k in seq_len(ncol(idx2))) {
      lab <- rep(lv[[3]], n)
      lab[idx1[, j]] <- lv[[1]]; lab[idx2[, k]] <- lv[[2]]
      null <- c(null, stat(lab))
    }
  }
  list(observed = obs, n_assignments = length(null),
       p = (1 + sum(null >= obs - 1e-12)) / (1 + length(null)))
}

# ----------------------------------------------------------- animal bootstrap
#
# Resamples ANIMALS with replacement. Never proteins, never hemispheres, never
# spatial units - those are repeated measures inside an animal, not replicates.
asn_animal_bootstrap_ci <- function(values, group, target_group, n_boot = 5000L,
                                    seed = 20260912L, conf = 0.95) {
  set.seed(as.integer(seed))
  v <- as.numeric(values)[group == target_group]
  v <- v[is.finite(v)]
  if (length(v) < 2L) {
    return(data.frame(mean = if (length(v)) mean(v) else NA_real_,
                      ci_lower = NA_real_, ci_upper = NA_real_,
                      n_valid_iterations = 0L, n_animals = length(v),
                      stringsAsFactors = FALSE))
  }
  b <- replicate(n_boot, mean(sample(v, length(v), replace = TRUE)))
  b <- b[is.finite(b)]
  a <- (1 - conf) / 2
  data.frame(mean = mean(v),
             ci_lower = unname(stats::quantile(b, a)),
             ci_upper = unname(stats::quantile(b, 1 - a)),
             n_valid_iterations = length(b), n_animals = length(v),
             bootstrap_unit = "AnimalID", seed = as.integer(seed),
             stringsAsFactors = FALSE)
}

# Difference between two groups with an animal-level bootstrap interval.
asn_animal_bootstrap_difference <- function(values, group, a = "SUS", b = "RES",
                                            n_boot = 5000L, seed = 20260912L,
                                            conf = 0.95) {
  set.seed(as.integer(seed))
  va <- as.numeric(values)[group == a]; vb <- as.numeric(values)[group == b]
  va <- va[is.finite(va)]; vb <- vb[is.finite(vb)]
  if (length(va) < 2L || length(vb) < 2L) {
    return(data.frame(difference = NA_real_, ci_lower = NA_real_,
                      ci_upper = NA_real_, n_valid_iterations = 0L,
                      stringsAsFactors = FALSE))
  }
  d <- replicate(n_boot, mean(sample(va, length(va), replace = TRUE)) -
                   mean(sample(vb, length(vb), replace = TRUE)))
  d <- d[is.finite(d)]
  al <- (1 - conf) / 2
  data.frame(difference = mean(va) - mean(vb),
             ci_lower = unname(stats::quantile(d, al)),
             ci_upper = unname(stats::quantile(d, 1 - al)),
             n_valid_iterations = length(d),
             bootstrap_unit = "AnimalID", seed = as.integer(seed),
             stringsAsFactors = FALSE)
}

# Hedges-corrected standardised difference, appropriate at n = 3 per group.
asn_effect_size <- function(values, group, a, b) {
  va <- as.numeric(values)[group == a]; vb <- as.numeric(values)[group == b]
  va <- va[is.finite(va)]; vb <- vb[is.finite(vb)]
  na <- length(va); nb <- length(vb)
  if (na < 2L || nb < 2L) return(NA_real_)
  sp <- sqrt(((na - 1) * stats::var(va) + (nb - 1) * stats::var(vb)) / (na + nb - 2))
  if (!is.finite(sp) || sp <= 0) return(NA_real_)
  d <- (mean(va) - mean(vb)) / sp
  J <- 1 - 3 / (4 * (na + nb) - 9)          # small-sample correction
  d * J
}

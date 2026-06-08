# Dataset-aware variance partitioning for spatial proteomics.
#
# Key design choices:
# - variancePartition::fitExtractVarPartModel() requires either all or no
#   categorical variables as random effects.
# - Therefore, categorical metadata terms are modeled as random effects.
# - Numeric metadata terms are modeled as fixed effects.
# - Samples with missing formula metadata are dropped before fitting.
# - Near-redundant categorical terms are screened with Cramer's V.
# - Outputs include ranked summary, distribution plot, and top-protein heatmap.

paths_file <- if (file.exists(file.path("R", "paths.R"))) {
  file.path("R", "paths.R")
} else {
  file.path("..", "R", "paths.R")
}

source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "dataset_inputs.R"))
source(repo_path("R", "qc_exploration_utils.R"))

run <- qc_args()
DATASET <- run$dataset

PATHS <- qc_paths("06_variance_partitioning", DATASET)

matrix_file <- path_or_env(
  "PROTEOMICS_VARPART_MATRIX_FILE",
  qc_resolve_matrix(DATASET),
  must_exist = FALSE
)

metadata_file <- qc_resolve_metadata(
  DATASET,
  env = "PROTEOMICS_VARPART_METADATA_FILE"
)

if (run$dry_run) {
  status <- qc_dry_run_contract(
    "03_qc_exploration/06_variance_partitioning.r",
    DATASET,
    matrix_file = matrix_file,
    metadata_file = metadata_file,
    paths = PATHS,
    extra = c(
      "Formula is adapted to available metadata terms.",
      "Categorical metadata terms are modeled as random effects for variancePartition.",
      "Numeric metadata terms are modeled as fixed effects.",
      "Uses (1|AnimalID) only when repeated samples per animal exist.",
      "Drops samples with missing formula metadata before fitting.",
      "Writes ranked variance, distribution, and top-protein heatmap figures."
    )
  )
  quit(status = status, save = "no")
}

required_pkgs <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "svglite",
  "variancePartition",
  "reformulas",
  "grid"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs)) {
  stop(
    "Missing required R package(s): ",
    paste(missing_pkgs, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages(
  invisible(lapply(required_pkgs, library, character.only = TRUE))
)

if (!file.exists(matrix_file)) {
  stop("VariancePartition matrix not found: ", matrix_file, call. = FALSE)
}

expr <- qc_read_expression(matrix_file, metadata_file, DATASET)
mat <- qc_impute_for_pca(expr$mat)

meta <- expr$meta[colnames(mat), , drop = FALSE]
meta <- as.data.frame(meta, stringsAsFactors = FALSE)

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

non_missing_nonempty <- function(x) {
  !is.na(x) & nzchar(as.character(x))
}

drop_if_present <- function(x, preferred, redundant) {
  if (preferred %in% x && redundant %in% x) {
    setdiff(x, redundant)
  } else {
    x
  }
}

cramers_v <- function(x, y) {
  ok <- non_missing_nonempty(x) & non_missing_nonempty(y)
  x <- as.factor(x[ok])
  y <- as.factor(y[ok])

  if (length(x) < 3L) return(NA_real_)
  if (nlevels(x) < 2L || nlevels(y) < 2L) return(NA_real_)

  tab <- table(x, y)

  if (any(dim(tab) < 2L)) return(NA_real_)

  ch <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))

  n <- sum(tab)
  if (!is.finite(ch$statistic) || n <= 0L) return(NA_real_)

  phi2 <- as.numeric(ch$statistic) / n
  r <- nrow(tab)
  k <- ncol(tab)

  denom <- min(k - 1L, r - 1L)
  if (denom <= 0L) return(NA_real_)

  sqrt(phi2 / denom)
}

numeric_eta2_by_group <- function(x, group) {
  ok <- is.finite(x) & non_missing_nonempty(group)
  x <- x[ok]
  group <- as.factor(group[ok])

  if (length(x) < 3L) return(NA_real_)
  if (nlevels(group) < 2L) return(NA_real_)

  grand_mean <- mean(x, na.rm = TRUE)
  ss_total <- sum((x - grand_mean)^2, na.rm = TRUE)

  if (!is.finite(ss_total) || ss_total <= 0) return(NA_real_)

  group_means <- tapply(x, group, mean, na.rm = TRUE)
  group_ns <- table(group)

  ss_between <- sum(group_ns[names(group_means)] * (group_means - grand_mean)^2)

  as.numeric(ss_between / ss_total)
}

make_term_summary <- function(meta, terms) {
  dplyr::bind_rows(lapply(terms, function(term) {
    x <- meta[[term]]
    ok <- non_missing_nonempty(x)

    data.frame(
      term = term,
      storage_mode = typeof(x),
      class = paste(class(x), collapse = ";"),
      n_total = length(x),
      n_nonmissing_nonempty = sum(ok),
      n_missing_or_empty = sum(!ok),
      n_unique_nonmissing = length(unique(x[ok])),
      example_values = paste(utils::head(unique(as.character(x[ok])), 8), collapse = "; "),
      stringsAsFactors = FALSE
    )
  }))
}

safe_filename <- function(x) {
  gsub("[^A-Za-z0-9_\\-]+", "_", x)
}

# -------------------------------------------------------------------------
# Candidate metadata terms
# -------------------------------------------------------------------------

candidate_terms <- c(
  "Group",
  "group",
  "ExpGroup",
  "Sex",
  "sex",
  "Region",
  "region",
  "Layer",
  "layer",
  "ReplicateGroup",
  "plate",
  "batch",
  "sample_prep",
  "run",
  "order"
)

candidate_terms <- intersect(candidate_terms, names(meta))

usable <- candidate_terms[vapply(candidate_terms, function(term) {
  x <- meta[[term]]
  ok <- non_missing_nonempty(x)

  sum(ok) >= 3L &&
    length(unique(x[ok])) >= 2L
}, logical(1))]

# Remove case-only duplicates.
usable <- usable[!duplicated(tolower(usable))]

# Prefer biologically more informative / standardized terms.
# If ExpGroup exists, it usually carries CON/RES/SUS and is preferable to simple group.
if ("ExpGroup" %in% usable) {
  usable <- setdiff(usable, c("Group", "group"))
}

usable <- drop_if_present(usable, "Sex", "sex")
usable <- drop_if_present(usable, "Region", "region")
usable <- drop_if_present(usable, "Layer", "layer")

# Add AnimalID only if repeated samples per animal exist.
if ("AnimalID" %in% names(meta)) {
  animal_ok <- non_missing_nonempty(meta$AnimalID)

  if (sum(animal_ok) >= 3L && any(table(meta$AnimalID[animal_ok]) > 1L)) {
    usable <- unique(c(usable, "AnimalID"))
  }
}

if (!length(usable)) {
  stop("No usable metadata terms for variance partitioning.", call. = FALSE)
}

term_summary_pre <- make_term_summary(meta, usable)

qc_write_csv(
  term_summary_pre,
  file.path(PATHS$tables, "metadata_terms_considered_pre_filter.csv")
)

# -------------------------------------------------------------------------
# Remove rows with missing formula metadata
# -------------------------------------------------------------------------

complete_for_formula <- Reduce(
  `&`,
  lapply(usable, function(term) non_missing_nonempty(meta[[term]]))
)

n_before <- ncol(mat)
n_after <- sum(complete_for_formula)

if (n_after < 3L) {
  stop(
    "Too few samples remain after filtering missing metadata for candidate terms: ",
    n_after,
    call. = FALSE
  )
}

if (n_after < n_before) {
  message(
    "Dropping ",
    n_before - n_after,
    " sample(s) with missing/empty metadata among candidate variance terms."
  )
}

mat <- mat[, complete_for_formula, drop = FALSE]
meta <- meta[complete_for_formula, , drop = FALSE]

# Re-check terms after sample filtering.
usable <- usable[vapply(usable, function(term) {
  x <- meta[[term]]
  ok <- non_missing_nonempty(x)

  sum(ok) >= 3L &&
    length(unique(x[ok])) >= 2L
}, logical(1))]

if (!length(usable)) {
  stop("No usable metadata terms remain after missing-data filtering.", call. = FALSE)
}

# -------------------------------------------------------------------------
# Split into numeric fixed effects and categorical random effects
# -------------------------------------------------------------------------

is_numeric_term <- vapply(usable, function(term) {
  is.numeric(meta[[term]]) || is.integer(meta[[term]])
}, logical(1))

numeric_terms <- usable[is_numeric_term]
categorical_terms <- usable[!is_numeric_term]

# Keep numeric terms only if they have enough distinct numeric values.
numeric_terms <- numeric_terms[vapply(numeric_terms, function(term) {
  x <- suppressWarnings(as.numeric(meta[[term]]))
  sum(is.finite(x)) >= 3L && length(unique(x[is.finite(x)])) >= 3L
}, logical(1))]

# Keep categorical terms only if they have enough levels.
categorical_terms <- categorical_terms[vapply(categorical_terms, function(term) {
  x <- as.factor(meta[[term]])
  nlevels(x) >= 2L
}, logical(1))]

# Remove categorical terms with one sample per every level, except AnimalID.
# Such terms often behave poorly as random effects and explain residual identity
# rather than structured variance.
categorical_terms <- categorical_terms[vapply(categorical_terms, function(term) {
  tab <- table(meta[[term]])
  if (term == "AnimalID") {
    any(tab > 1L)
  } else {
    length(tab) < nrow(meta) && any(tab > 1L)
  }
}, logical(1))]

if (!length(numeric_terms) && !length(categorical_terms)) {
  stop(
    "No usable numeric or categorical terms remain for variance partitioning.",
    call. = FALSE
  )
}

# Convert categorical terms to factors explicitly.
for (term in categorical_terms) {
  meta[[term]] <- droplevels(as.factor(meta[[term]]))
}

# Numeric terms should be numeric explicitly.
for (term in numeric_terms) {
  meta[[term]] <- as.numeric(meta[[term]])
}

# -------------------------------------------------------------------------
# Confounding / redundancy diagnostics
# -------------------------------------------------------------------------

if (length(categorical_terms) >= 2L) {
  cat_pairs <- utils::combn(categorical_terms, 2, simplify = FALSE)

  categorical_confounding <- dplyr::bind_rows(lapply(cat_pairs, function(pair) {
    data.frame(
      term_1 = pair[[1]],
      term_2 = pair[[2]],
      cramers_v = cramers_v(meta[[pair[[1]]]], meta[[pair[[2]]]]),
      stringsAsFactors = FALSE
    )
  }))

  categorical_confounding <- categorical_confounding |>
    dplyr::arrange(dplyr::desc(.data$cramers_v))

  qc_write_csv(
    categorical_confounding,
    file.path(PATHS$tables, "categorical_term_cramers_v.csv")
  )
} else {
  categorical_confounding <- data.frame()
}

if (length(numeric_terms) && length(categorical_terms)) {
  numeric_categorical_confounding <- dplyr::bind_rows(lapply(numeric_terms, function(nm) {
    dplyr::bind_rows(lapply(categorical_terms, function(cat) {
      data.frame(
        numeric_term = nm,
        categorical_term = cat,
        eta2_numeric_by_categorical = numeric_eta2_by_group(meta[[nm]], meta[[cat]]),
        stringsAsFactors = FALSE
      )
    }))
  }))

  numeric_categorical_confounding <- numeric_categorical_confounding |>
    dplyr::arrange(dplyr::desc(.data$eta2_numeric_by_categorical))

  qc_write_csv(
    numeric_categorical_confounding,
    file.path(PATHS$tables, "numeric_categorical_eta2.csv")
  )
} else {
  numeric_categorical_confounding <- data.frame()
}

# Optional conservative pruning of near-redundant categorical terms.
# Threshold 0.98 removes only near-deterministic redundancy.
pruned_terms <- categorical_terms

if (nrow(categorical_confounding)) {
  high_pairs <- categorical_confounding |>
    dplyr::filter(!is.na(.data$cramers_v), .data$cramers_v >= 0.98)

  if (nrow(high_pairs)) {
    for (i in seq_len(nrow(high_pairs))) {
      a <- high_pairs$term_1[[i]]
      b <- high_pairs$term_2[[i]]

      if (!(a %in% pruned_terms && b %in% pruned_terms)) next

      drop_term <- NULL

      if ("ExpGroup" %in% c(a, b) && any(c("Group", "group") %in% c(a, b))) {
        drop_term <- intersect(c("Group", "group"), c(a, b))[1]
      } else if ("AnimalID" %in% c(a, b) && any(c("plate", "batch") %in% c(a, b))) {
        drop_term <- intersect(c("plate", "batch"), c(a, b))[1]
      } else if (any(c("Region", "region", "Layer", "layer") %in% c(a, b)) &&
                 any(c("plate", "batch") %in% c(a, b))) {
        drop_term <- intersect(c("plate", "batch"), c(a, b))[1]
      } else {
        drop_term <- b
      }

      pruned_terms <- setdiff(pruned_terms, drop_term)
    }
  }
}

categorical_terms <- pruned_terms

term_summary_final <- make_term_summary(meta, c(numeric_terms, categorical_terms))

qc_write_csv(
  term_summary_final,
  file.path(PATHS$tables, "metadata_terms_used_final.csv")
)

# -------------------------------------------------------------------------
# Build variancePartition formula
# -------------------------------------------------------------------------

random_terms <- paste0("(1|", categorical_terms, ")")
formula_terms <- c(numeric_terms, random_terms)

if (!length(formula_terms)) {
  stop("No usable metadata terms for variance partitioning.", call. = FALSE)
}

form <- as.formula(paste("~", paste(formula_terms, collapse = " + ")))

message("VariancePartition formula: ", deparse(form))

# -------------------------------------------------------------------------
# Canonical correlation diagnostics
# -------------------------------------------------------------------------

cca <- try(
  variancePartition::canCorPairs(form, meta),
  silent = TRUE
)

if (!inherits(cca, "try-error")) {
  cca_df <- as.data.frame(as.table(as.matrix(cca)))
  names(cca_df) <- c("term_1", "term_2", "canonical_correlation")

  cca_df <- cca_df |>
    dplyr::filter(.data$term_1 != .data$term_2) |>
    dplyr::arrange(dplyr::desc(.data$canonical_correlation))

  qc_write_csv(
    cca_df,
    file.path(PATHS$tables, "metadata_term_canonical_correlations.csv")
  )
} else {
  cca_df <- data.frame()

  writeLines(
    c(
      "canCorPairs() failed.",
      "",
      as.character(cca)
    ),
    file.path(PATHS$logs, "canCorPairs_failed.txt")
  )
}

# -------------------------------------------------------------------------
# Fit variance partitioning model
# -------------------------------------------------------------------------

vp <- try(
  variancePartition::fitExtractVarPartModel(mat, form, meta),
  silent = TRUE
)

if (inherits(vp, "try-error")) {
  writeLines(
    c(
      "# Variance partitioning failed",
      "",
      paste0("Dataset: ", DATASET),
      paste0("Formula: ", deparse(form)),
      "",
      "Error:",
      as.character(vp)
    ),
    file.path(PATHS$reports, "variance_partitioning_failed.md")
  )

  stop(
    "Variance partitioning failed. See report: ",
    file.path(PATHS$reports, "variance_partitioning_failed.md"),
    call. = FALSE
  )
}

vp <- variancePartition::sortCols(vp)

vp_df <- as.data.frame(vp)
vp_df$Protein <- rownames(vp_df)

vp_long <- vp_df |>
  tidyr::pivot_longer(
    -Protein,
    names_to = "term",
    values_to = "variance_fraction"
  )

median_vp <- vp_long |>
  dplyr::group_by(.data$term) |>
  dplyr::summarise(
    median_variance_fraction = median(.data$variance_fraction, na.rm = TRUE),
    mean_variance_fraction = mean(.data$variance_fraction, na.rm = TRUE),
    max_variance_fraction = max(.data$variance_fraction, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(.data$median_variance_fraction))

top_driven <- vp_long |>
  dplyr::group_by(.data$term) |>
  dplyr::slice_max(
    .data$variance_fraction,
    n = 50,
    with_ties = FALSE
  ) |>
  dplyr::ungroup()

qc_write_csv(
  vp_long,
  file.path(PATHS$tables, "variance_fraction_by_protein.csv")
)

qc_write_csv(
  median_vp,
  file.path(PATHS$tables, "median_variance_by_term.csv")
)

qc_write_csv(
  top_driven,
  file.path(PATHS$tables, "top_term_driven_proteins.csv")
)

qc_write_xlsx(
  list(
    variance_fraction = vp_long,
    median_by_term = median_vp,
    top_term_driven = top_driven,
    terms_considered_pre_filter = term_summary_pre,
    terms_used_final = term_summary_final,
    categorical_cramers_v = categorical_confounding,
    numeric_categorical_eta2 = numeric_categorical_confounding,
    canonical_correlations = cca_df
  ),
  file.path(PATHS$tables, "variance_partitioning.xlsx")
)

# -------------------------------------------------------------------------
# Plot helpers
# -------------------------------------------------------------------------

plot_palette <- c(
  "Residuals" = "#B8B8B8",
  "ExpGroup" = "#3B5B92",
  "Group" = "#3B5B92",
  "group" = "#3B5B92",
  "Sex" = "#8C4E70",
  "sex" = "#8C4E70",
  "Region" = "#4F7F6F",
  "region" = "#4F7F6F",
  "Layer" = "#9A6A3A",
  "layer" = "#9A6A3A",
  "AnimalID" = "#6B6B8A",
  "plate" = "#7A7A7A",
  "batch" = "#7A7A7A",
  "ReplicateGroup" = "#7A7A7A",
  "sample_prep" = "#7A7A7A",
  "run" = "#7A7A7A",
  "order" = "#7A7A7A"
)

fallback_palette <- c(
  "#3B5B92",
  "#4F7F6F",
  "#8C4E70",
  "#9A6A3A",
  "#6B6B8A",
  "#7A7A7A",
  "#A05A4F",
  "#5F758E"
)

all_terms <- sort(unique(vp_long$term))
missing_cols <- setdiff(all_terms, names(plot_palette))

if (length(missing_cols)) {
  extra_cols <- rep(fallback_palette, length.out = length(missing_cols))
  names(extra_cols) <- missing_cols
  plot_palette <- c(plot_palette, extra_cols)
}

term_order <- median_vp |>
  dplyr::arrange(.data$median_variance_fraction) |>
  dplyr::pull(.data$term)

vp_long <- vp_long |>
  dplyr::mutate(
    term = factor(.data$term, levels = term_order),
    variance_percent = 100 * .data$variance_fraction
  )

median_vp_plot <- median_vp |>
  dplyr::mutate(
    term = factor(.data$term, levels = term_order),
    median_percent = 100 * .data$median_variance_fraction,
    mean_percent = 100 * .data$mean_variance_fraction,
    max_percent = 100 * .data$max_variance_fraction
  )

theme_varpart <- function(base_size = 7.5, base_family = "Arial") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.25, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.25, colour = "black"),
      axis.ticks.length = grid::unit(1.5, "mm"),
      axis.text = ggplot2::element_text(colour = "black", size = base_size),
      axis.title = ggplot2::element_text(colour = "black", size = base_size + 0.5),
      plot.title = ggplot2::element_text(
        colour = "black",
        size = base_size + 1.5,
        face = "bold",
        hjust = 0
      ),
      plot.subtitle = ggplot2::element_text(
        colour = "grey25",
        size = base_size,
        hjust = 0,
        margin = ggplot2::margin(b = 4)
      ),
      plot.caption = ggplot2::element_text(
        colour = "grey35",
        size = base_size - 1,
        hjust = 0,
        margin = ggplot2::margin(t = 4)
      ),
      panel.grid.major.x = ggplot2::element_line(
        linewidth = 0.2,
        colour = "grey88"
      ),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "none",
      plot.margin = ggplot2::margin(4, 6, 4, 4)
    )
}

# -------------------------------------------------------------------------
# Plot 1: ranked median variance explained
# -------------------------------------------------------------------------

p_median <- ggplot(
  median_vp_plot,
  aes(
    x = .data$term,
    y = .data$median_percent,
    fill = .data$term
  )
) +
  geom_col(
    width = 0.68,
    linewidth = 0,
    alpha = 0.95
  ) +
  geom_text(
    aes(label = sprintf("%.1f", .data$median_percent)),
    hjust = -0.15,
    size = 2.1,
    family = "Arial",
    colour = "black"
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = plot_palette, drop = FALSE) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.12)),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    title = "Median variance explained",
    subtitle = paste0(DATASET, " proteome"),
    x = NULL,
    y = "Median variance explained"
  ) +
  theme_varpart(base_size = 7.5) +
  theme(
    panel.grid.major.x = element_line(linewidth = 0.2, colour = "grey88")
  )

ranked_height <- max(55, 8 + 7 * length(term_order))

ggsave(
  file.path(PATHS$figures, "variance_partition_median_ranked.svg"),
  p_median,
  width = 95,
  height = ranked_height,
  units = "mm",
  device = svglite::svglite
)

ggsave(
  file.path(PATHS$figures, "variance_partition_median_ranked.pdf"),
  p_median,
  width = 95,
  height = ranked_height,
  units = "mm",
  device = grDevices::cairo_pdf
)

# -------------------------------------------------------------------------
# Plot 2: distribution across proteins
# -------------------------------------------------------------------------

p_distribution <- ggplot(
  vp_long,
  aes(
    x = .data$term,
    y = .data$variance_percent,
    fill = .data$term
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.25,
    alpha = 0.55,
    colour = "grey30"
  ) +
  geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    linewidth = 0.25,
    alpha = 0.85,
    colour = "black"
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    size = 0.8,
    colour = "black"
  ) +
  coord_flip() +
  scale_fill_manual(values = plot_palette, drop = FALSE) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0.01, 0.04))
  ) +
  labs(
    title = "Variance explained across proteins",
    subtitle = paste0(DATASET, " proteome"),
    x = NULL,
    y = "Variance explained"
  ) +
  theme_varpart(base_size = 7.5)

distribution_height <- max(60, 8 + 7 * length(term_order))

ggsave(
  file.path(PATHS$figures, "variance_partition_distribution.svg"),
  p_distribution,
  width = 110,
  height = distribution_height,
  units = "mm",
  device = svglite::svglite
)

ggsave(
  file.path(PATHS$figures, "variance_partition_distribution.pdf"),
  p_distribution,
  width = 110,
  height = distribution_height,
  units = "mm",
  device = grDevices::cairo_pdf
)

# Backward-compatible filename for existing pipeline expectations.
ggsave(
  file.path(PATHS$figures, "variance_partition_summary.svg"),
  p_distribution,
  width = 110,
  height = distribution_height,
  units = "mm",
  device = svglite::svglite
)

# -------------------------------------------------------------------------
# Plot 3: strongest term-driven proteins
# -------------------------------------------------------------------------

n_heatmap_terms <- min(8L, nrow(median_vp))

top_heatmap_terms <- median_vp |>
  dplyr::arrange(dplyr::desc(.data$median_variance_fraction)) |>
  dplyr::slice_head(n = n_heatmap_terms) |>
  dplyr::pull(.data$term)

protein_rank_df <- vp_long |>
  dplyr::filter(as.character(.data$term) %in% top_heatmap_terms) |>
  dplyr::group_by(.data$Protein) |>
  dplyr::summarise(
    max_variance_percent = max(.data$variance_percent, na.rm = TRUE),
    dominant_term = as.character(.data$term[which.max(.data$variance_percent)]),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(.data$max_variance_percent))

n_heatmap_proteins <- min(40L, nrow(protein_rank_df))

top_heatmap_proteins <- protein_rank_df |>
  dplyr::slice_head(n = n_heatmap_proteins) |>
  dplyr::pull(.data$Protein)

if (length(top_heatmap_terms) && length(top_heatmap_proteins)) {
  heatmap_df <- vp_long |>
    dplyr::filter(
      .data$Protein %in% top_heatmap_proteins,
      as.character(.data$term) %in% top_heatmap_terms
    ) |>
    dplyr::mutate(
      Protein = factor(
        .data$Protein,
        levels = rev(top_heatmap_proteins)
      ),
      term = factor(
        as.character(.data$term),
        levels = rev(top_heatmap_terms)
      )
    )

  p_heatmap <- ggplot(
    heatmap_df,
    aes(
      x = .data$term,
      y = .data$Protein,
      fill = .data$variance_percent
    )
  ) +
    geom_tile(
      linewidth = 0.15,
      colour = "white"
    ) +
    scale_fill_gradient(
      low = "#F3F3F3",
      high = "#2F4F6F",
      name = "Variance\nexplained",
      labels = function(x) paste0(x, "%")
    ) +
    labs(
      title = "Top term-driven proteins",
      subtitle = paste0(DATASET, " proteome"),
      x = NULL,
      y = NULL
    ) +
    theme_classic(base_size = 6.5, base_family = "Arial") +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(
        colour = "black",
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      axis.text.y = element_text(colour = "black", size = 5.5),
      plot.title = element_text(
        colour = "black",
        size = 8.5,
        face = "bold",
        hjust = 0
      ),
      plot.subtitle = element_text(
        colour = "grey25",
        size = 6.5,
        hjust = 0,
        margin = margin(b = 4)
      ),
      legend.position = "right",
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 5.5),
      legend.key.height = grid::unit(12, "mm"),
      plot.margin = margin(4, 6, 4, 4)
    )

  ggsave(
    file.path(PATHS$figures, "variance_partition_top_proteins_heatmap.svg"),
    p_heatmap,
    width = 115,
    height = 120,
    units = "mm",
    device = svglite::svglite
  )

  ggsave(
    file.path(PATHS$figures, "variance_partition_top_proteins_heatmap.pdf"),
    p_heatmap,
    width = 115,
    height = 120,
    units = "mm",
    device = grDevices::cairo_pdf
  )
}

# -------------------------------------------------------------------------
# Group / technical confounding screen
# -------------------------------------------------------------------------

group_terms <- intersect(c("ExpGroup", "Group", "group", "group2"), names(meta))
technical_terms <- intersect(c("plate", "batch", "sample_prep", "run", "order"), names(meta))

confound_rows <- dplyr::bind_rows(lapply(group_terms, function(g) {
  dplyr::bind_rows(lapply(technical_terms, function(t) {
    data.frame(
      group_term = g,
      technical_term = t,
      cramers_v_group_by_technical = cramers_v(meta[[g]], meta[[t]]),
      stringsAsFactors = FALSE
    )
  }))
}))

if (nrow(confound_rows)) {
  confound_rows <- confound_rows |>
    dplyr::arrange(dplyr::desc(.data$cramers_v_group_by_technical))

  qc_write_csv(
    confound_rows,
    file.path(PATHS$tables, "group_technical_confounding_screen.csv")
  )
}

flag <- if (
  nrow(confound_rows) &&
    any(confound_rows$cramers_v_group_by_technical >= 0.8, na.rm = TRUE)
) {
  "WARN"
} else {
  "PASS"
}

# -------------------------------------------------------------------------
# Report and manifest
# -------------------------------------------------------------------------

report_lines <- c(
  "# Variance Partitioning",
  "",
  paste0("Dataset: ", DATASET),
  paste0("Input matrix: ", matrix_file),
  paste0("Input metadata: ", metadata_file),
  "",
  "## Model",
  "",
  paste0("Formula: ", deparse(form)),
  "",
  "## Terms",
  "",
  paste0(
    "Numeric fixed-effect terms: ",
    if (length(numeric_terms)) paste(numeric_terms, collapse = ", ") else "none"
  ),
  paste0(
    "Categorical random-effect terms: ",
    if (length(categorical_terms)) paste(categorical_terms, collapse = ", ") else "none"
  ),
  "",
  "## Samples",
  "",
  paste0("Samples before metadata filtering: ", n_before),
  paste0("Samples after metadata filtering: ", n_after),
  "",
  "## Confounding screen",
  "",
  paste0(
    flag,
    ": Group/technical confounding screen written where metadata allowed."
  ),
  "",
  "## Figures",
  "",
  "- variance_partition_median_ranked.svg",
  "- variance_partition_median_ranked.pdf",
  "- variance_partition_distribution.svg",
  "- variance_partition_distribution.pdf",
  "- variance_partition_summary.svg",
  "- variance_partition_top_proteins_heatmap.svg",
  "- variance_partition_top_proteins_heatmap.pdf"
)

if (nrow(categorical_confounding)) {
  top_cat <- categorical_confounding |>
    dplyr::filter(!is.na(.data$cramers_v)) |>
    dplyr::slice_head(n = 10)

  if (nrow(top_cat)) {
    report_lines <- c(
      report_lines,
      "",
      "Top categorical metadata associations by Cramer's V:",
      "",
      utils::capture.output(print(top_cat, row.names = FALSE))
    )
  }
}

if (exists("cca_df") && nrow(cca_df)) {
  top_cca <- cca_df |>
    dplyr::filter(!is.na(.data$canonical_correlation)) |>
    dplyr::slice_head(n = 10)

  if (nrow(top_cca)) {
    report_lines <- c(
      report_lines,
      "",
      "Top canonical correlations:",
      "",
      utils::capture.output(print(top_cca, row.names = FALSE))
    )
  }
}

writeLines(
  report_lines,
  file.path(PATHS$reports, "variance_partitioning_summary.md")
)

write_run_manifest(
  file.path(PATHS$logs, "run_manifest.yml"),
  inputs = list(
    matrix = matrix_file,
    metadata = metadata_file
  ),
  outputs = list(
    tables = PATHS$tables,
    figures = PATHS$figures,
    reports = PATHS$reports
  ),
  parameters = list(
    dataset = DATASET,
    formula = deparse(form),
    numeric_fixed_terms = numeric_terms,
    categorical_random_terms = categorical_terms,
    samples_before_metadata_filter = n_before,
    samples_after_metadata_filter = n_after
  ),
  notes = paste(
    "Canonical variance partitioning with adaptive formula.",
    "Categorical metadata terms are modeled as random effects.",
    "Numeric metadata terms are modeled as fixed effects.",
    "Missing formula metadata samples are dropped before model fitting.",
    "Near-redundant categorical terms are screened using Cramer's V.",
    "Figures use restrained colors, compact labels, and vector output."
  )
)

message("Variance partitioning complete for dataset: ", DATASET)
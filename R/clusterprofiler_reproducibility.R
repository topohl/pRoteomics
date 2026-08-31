# Deterministic RNG, fgsea controls, and targeted-run helpers for clusterProfiler.
#
# Reproducibility contract (measured, not assumed)
# ------------------------------------------------
# What is pinned: the ranked input, the per-comparison derived seed, RNGkind
# ("L'Ecuyer-CMRG"/"Inversion"/"Rejection"), nPermSimple, by = "fgsea", and
# clusterProfiler's own logical `seed` flag (kept FALSE so the local RNG scope
# below governs). Software versions are pinned by renv.
#
# What that guarantees: within one execution context, the same seed and input
# reproduce results field-identically. Verified for
# neuron_soma/CA1_vs_mean_other_soma_regions (gseGO BP, seed 1751020239,
# nPermSimple 100000) across separate clean R --vanilla sessions.
#
# What it does NOT guarantee: bit-identical p-values across *different*
# execution contexts. clusterProfiler routes GSEA through
# DOSE:::GSEA_fgsea(), which hardcodes nproc = 0, so fgsea:::setUpBPPARAM()
# falls back to the ambient bpparam(). Backend choice was measured NOT to be
# the source of cross-context drift: runs under the default SnowParam (30
# workers, unseeded), an explicit SnowParam with a pinned RNGseed at 8 and at
# 4 workers, and a run preceded by another GSEA call in the same session were
# all field-identical to one another. Pinning BPPARAM therefore buys no
# additional determinism here and is deliberately not done.
#
# Observed cross-context tolerance: enrichmentScore agrees to ~1.6e-15 (last
# bit), setSize/rank/leading_edge/core_enrichment agree exactly, and that
# last-bit difference propagates through fgseaMultilevel's adaptive estimator
# to <= ~2.4e-05 in NES and <= ~2.3e-05 in p-value/FDR. In the audited
# comparison this produced zero FDR-0.05 crossings and an identical Figure-2f
# display selection, i.e. no inferential or presentational consequence.
# Reproducibility is therefore asserted numerically within that tolerance,
# not as floating-point bit identity.

clusterprofiler_integer_scalar <- function(value, name, min_value = 1L,
                                           max_value = .Machine$integer.max) {
  if (length(value) != 1L || !is.numeric(value) || is.na(value) ||
      !is.finite(value) || value != floor(value) ||
      value < min_value || value > max_value) {
    stop(
      name, " must be one finite whole number between ", min_value,
      " and ", max_value, ".",
      call. = FALSE
    )
  }
  as.integer(value)
}

validate_clusterprofiler_gsea_config <- function(cfg) {
  if (!is.list(cfg) || !is.list(cfg$analysis)) {
    stop("clusterProfiler config must contain an analysis mapping.", call. = FALSE)
  }
  cfg$analysis$gsea_seed_base <- clusterprofiler_integer_scalar(
    cfg$analysis$gsea_seed_base,
    "analysis.gsea_seed_base"
  )
  cfg$analysis$n_perm_simple <- clusterprofiler_integer_scalar(
    cfg$analysis$n_perm_simple,
    "analysis.n_perm_simple"
  )
  cfg
}

clusterprofiler_gsea_analysis_types <- function(ontology = "BP") {
  ontology <- toupper(trimws(as.character(ontology)))
  if (length(ontology) != 1L || is.na(ontology) || !nzchar(ontology)) {
    stop("ontology must be one non-empty value.", call. = FALSE)
  }
  c(
    go = paste0("gseGO_", ontology),
    kegg = "gseKEGG",
    kegg_predefined = "gseKEGG_predefined",
    nk3r = "GSEA_NK3R"
  )
}

derive_clusterprofiler_gsea_seed <- function(base_seed, comparison, analysis_type) {
  base_seed <- clusterprofiler_integer_scalar(base_seed, "base_seed")
  comparison <- trimws(as.character(comparison))
  analysis_type <- trimws(as.character(analysis_type))
  if (length(comparison) != 1L || is.na(comparison) || !nzchar(comparison)) {
    stop("comparison must be one non-empty value.", call. = FALSE)
  }
  if (length(analysis_type) != 1L || is.na(analysis_type) || !nzchar(analysis_type)) {
    stop("analysis_type must be one non-empty value.", call. = FALSE)
  }

  # A dependency-free, platform-stable rolling hash. Length prefixes prevent
  # separator ambiguity, and all arithmetic stays below exact-double limits.
  modulus <- 2147483646
  key <- paste0(
    nchar(comparison, type = "bytes"), ":", comparison, "|",
    nchar(analysis_type, type = "bytes"), ":", analysis_type
  )
  hash <- as.double(base_seed) %% modulus
  for (code in utf8ToInt(enc2utf8(key))) {
    hash <- (hash * 131 + as.double(code)) %% modulus
  }
  seed <- as.integer(hash)
  if (seed < 1L) seed <- 1L
  seed
}

clusterprofiler_fgsea_control_args <- function(n_perm_simple) {
  list(
    by = "fgsea",
    # clusterProfiler exposes a logical seed flag. The numeric per-comparison
    # seed is established by the local RNG scope below, so its internal flag
    # must remain FALSE rather than replacing that state.
    seed = FALSE,
    nPermSimple = clusterprofiler_integer_scalar(n_perm_simple, "n_perm_simple")
  )
}

run_with_stable_gsea_rng <- function(gsea_fun, gsea_seed, ...) {
  if (!is.function(gsea_fun)) stop("gsea_fun must be a function.", call. = FALSE)
  gsea_seed <- clusterprofiler_integer_scalar(gsea_seed, "gsea_seed")
  if (!requireNamespace("withr", quietly = TRUE)) {
    stop("Package 'withr' is required for deterministic GSEA RNG scoping.", call. = FALSE)
  }
  dots <- list(...)
  withr::with_preserve_seed({
    previous_kind <- RNGkind()
    tryCatch({
      RNGkind(
        kind = "L'Ecuyer-CMRG", normal.kind = "Inversion",
        sample.kind = "Rejection"
      )
      set.seed(gsea_seed)
      do.call(gsea_fun, dots)
    }, finally = {
      do.call(RNGkind, as.list(previous_kind))
    })
  })
}

run_seeded_clusterprofiler_gsea <- function(gsea_fun, gsea_seed, n_perm_simple, ...) {
  if (!is.function(gsea_fun)) stop("gsea_fun must be a function.", call. = FALSE)
  gsea_seed <- clusterprofiler_integer_scalar(gsea_seed, "gsea_seed")
  dots <- list(...)
  reserved <- intersect(names(dots), c("by", "seed", "nPermSimple"))
  if (length(reserved)) {
    stop(
      "GSEA RNG controls must be supplied by run_seeded_clusterprofiler_gsea, not via ...: ",
      paste(reserved, collapse = ", "),
      call. = FALSE
    )
  }
  call_args <- c(dots, clusterprofiler_fgsea_control_args(n_perm_simple))
  run_with_stable_gsea_rng(
    function() do.call(gsea_fun, call_args), gsea_seed = gsea_seed
  )
}

clusterprofiler_gsea_reproducibility_table <- function(comparison, ontology,
                                                       gsea_seed_base,
                                                       n_perm_simple,
                                                       requested = NULL) {
  analysis_types <- clusterprofiler_gsea_analysis_types(ontology)
  if (is.null(requested)) requested <- rep(TRUE, length(analysis_types))
  if (is.null(names(requested))) names(requested) <- names(analysis_types)
  if (!setequal(names(requested), names(analysis_types))) {
    stop("requested must name every clusterProfiler GSEA analysis type.", call. = FALSE)
  }
  requested <- as.logical(requested[names(analysis_types)])
  n_perm_simple <- clusterprofiler_integer_scalar(n_perm_simple, "n_perm_simple")
  gsea_seed_base <- clusterprofiler_integer_scalar(gsea_seed_base, "gsea_seed_base")
  data.frame(
    comparison = rep(as.character(comparison), length(analysis_types)),
    analysis_type = unname(analysis_types),
    requested = unname(requested),
    gsea_seed_base = rep(gsea_seed_base, length(analysis_types)),
    gsea_seed = vapply(
      unname(analysis_types),
      function(type) derive_clusterprofiler_gsea_seed(gsea_seed_base, comparison, type),
      integer(1)
    ),
    n_perm_simple = rep(n_perm_simple, length(analysis_types)),
    rng_kind = rep("L'Ecuyer-CMRG/Inversion/Rejection", length(analysis_types)),
    clusterprofiler_by = rep("fgsea", length(analysis_types)),
    clusterprofiler_seed_argument = rep(FALSE, length(analysis_types)),
    stringsAsFactors = FALSE
  )
}

clusterprofiler_reproducibility_audit_matches <- function(path, comparison, ontology,
                                                          gsea_seed_base,
                                                          n_perm_simple) {
  if (!file.exists(path)) return(FALSE)
  observed <- tryCatch(
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(observed)) return(FALSE)
  expected <- clusterprofiler_gsea_reproducibility_table(
    comparison = comparison,
    ontology = ontology,
    gsea_seed_base = gsea_seed_base,
    n_perm_simple = n_perm_simple
  )
  required <- setdiff(names(expected), "requested")
  if (!all(required %in% names(observed)) || nrow(observed) != nrow(expected)) return(FALSE)
  observed <- observed[order(observed$analysis_type, method = "radix"), required, drop = FALSE]
  expected <- expected[order(expected$analysis_type, method = "radix"), required, drop = FALSE]
  rownames(observed) <- rownames(expected) <- NULL
  identical(observed, expected)
}

resolve_clusterprofiler_comparison_filter <- function(
    args = commandArgs(trailingOnly = TRUE),
    env_value = Sys.getenv("PROTEOMICS_CLUSTERPROFILER_COMPARISON", unset = "")) {
  hits <- which(args == "--comparison")
  if (length(hits) > 1L) stop("--comparison may be supplied only once.", call. = FALSE)
  if (length(hits) == 1L && hits[[1]] == length(args)) {
    stop("--comparison requires an exact comparison name.", call. = FALSE)
  }
  cli_value <- if (length(hits)) args[[hits[[1]] + 1L]] else ""
  value <- trimws(if (nzchar(cli_value)) cli_value else env_value)
  if (!nzchar(value)) return("")
  if (grepl("[/\\\\]", value)) {
    stop("The comparison filter must be a comparison name, not a path.", call. = FALSE)
  }
  value <- sub("\\.csv$", "", value, ignore.case = TRUE)
  if (!grepl("^[A-Za-z0-9._-]+$", value)) {
    stop("The comparison filter contains unsupported characters.", call. = FALSE)
  }
  value
}

filter_clusterprofiler_comparisons <- function(comparison_names, comparison_filter = "") {
  comparison_names <- as.character(comparison_names)
  if (!nzchar(comparison_filter)) return(seq_along(comparison_names))
  hits <- which(comparison_names == comparison_filter)
  if (length(hits) != 1L) {
    stop(
      "Requested clusterProfiler comparison was not found exactly once: ", comparison_filter,
      call. = FALSE
    )
  }
  hits
}

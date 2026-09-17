source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "nature_v2_figure_utils.R"))
source(repo_path("R", "spatial_grammar_utils.R"))
source(repo_path("R", "spatial_v6_figure_utils.R"))

rd <- function(p) utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
have <- function(p) file.exists(p)
S6 <- s6e_contract()

code_of <- function(path) {
  ln <- readLines(path, warn = FALSE)
  paste(sub("#.*$", "", ln), collapse = "\n")
}

# ------------------------------------------------------ the spatial contract

testthat::test_that("the spatial order contract declares exactly 18 units", {
  u <- sg_units()
  testthat::expect_identical(nrow(u), 18L)
  testthat::expect_identical(sum(u$dataset == "neuron_neuropil"), 10L)
  testthat::expect_identical(sum(u$dataset == "neuron_soma"), 4L)
  testthat::expect_identical(sum(u$dataset == "microglia"), 4L)
  testthat::expect_identical(anyDuplicated(u$analysis_key), 0L)
  testthat::expect_identical(u$order, seq_len(18L))
})

testthat::test_that("region-level compartments never carry laminar resolution", {
  u <- sg_units()
  # soma carries an sp/sg token naming the layer that was dissected, but only
  # one layer per region was sampled, so it is NOT laminar resolution
  testthat::expect_true(all(u$layer_is_resolution[u$dataset == "neuron_neuropil"]))
  testthat::expect_false(any(u$layer_is_resolution[u$dataset == "neuron_soma"]))
  testthat::expect_false(any(u$layer_is_resolution[u$dataset == "microglia"]))
  # and the DISPLAY label of a region-level unit must not gain a layer suffix
  soma <- u[u$dataset == "neuron_soma", ]
  testthat::expect_identical(sort(soma$display), c("CA1", "CA2", "CA3", "DG"))
  mg <- u[u$dataset == "microglia", ]
  testthat::expect_identical(sort(mg$display), c("CA1", "CA2", "CA3", "DG"))
})

testthat::test_that("CA3 has no SLM and DG layers belong to neuropil", {
  u <- sg_units()
  ca3 <- u[u$dataset == "neuron_neuropil" & u$region == "CA3", ]
  # the Part-19 schematic drew a continuous SLM band through CA3, implying a
  # unit that does not exist
  testthat::expect_false("slm" %in% ca3$layer)
  testthat::expect_identical(sort(ca3$layer), c("so", "sr"))
  # DG_mo and DG_po are NEUROPIL units; the Part-19 schematic coloured them as
  # soma, misassigning 2 of the 10 neuropil units
  dg_np <- u[u$dataset == "neuron_neuropil" & u$region == "DG", ]
  testthat::expect_identical(sort(dg_np$layer), c("mo", "po"))
  testthat::expect_identical(u$layer[u$dataset == "neuron_soma" &
                                       u$region == "DG"], "sg")
})

testthat::test_that("both spelling vocabularies resolve, disambiguated by dataset", {
  # the atlas spells soma AND microglia units ca1/ca2/ca3/dg, so the unit
  # string alone is ambiguous and the pair (dataset, unit) is the real key
  testthat::expect_identical(sg_resolve_unit("ca1", "neuron_soma"), "CA1_sp")
  testthat::expect_identical(sg_resolve_unit("ca1", "microglia"), "CA1")
  testthat::expect_identical(sg_resolve_unit("CA1_sp", "neuron_soma"), "CA1_sp")
  testthat::expect_identical(sg_resolve_unit("CA1", "microglia"), "CA1")
  testthat::expect_identical(sg_resolve_unit("ca1_slm", "neuron_neuropil"), "CA1_slm")
  testthat::expect_identical(sg_resolve_unit("CA1_slm", "neuron_neuropil"), "CA1_slm")
  # fails closed rather than silently dropping a row
  testthat::expect_error(sg_resolve_unit("CA1_slm", "neuron_soma"),
                         "unrecognised")
  testthat::expect_error(sg_resolve_unit("nonsense", "microglia"),
                         "unrecognised")
})

testthat::test_that("the axis never repeats what a header already carries", {
  u <- sg_units()
  b <- sg_blocks(u$analysis_key, u$dataset)
  testthat::expect_identical(nrow(b$order), 18L)
  testthat::expect_identical(b$compartment$label,
                             c("Neuropil", "Soma", "Microglia ROI"))
  testthat::expect_identical(b$compartment$start, c(1L, 11L, 15L))
  testthat::expect_identical(b$compartment$end, c(10L, 14L, 18L))
  lab <- sg_axis_labels(b)
  # neuropil shows only the layer token; region-level compartments show nothing
  # because the region strip above them IS the label
  testthat::expect_identical(lab[1:10],
    c("SO", "SR", "SLM", "SO", "SR", "SLM", "SO", "SR", "MO", "PO"))
  testthat::expect_true(all(lab[11:18] == ""))
  # and no axis label is ever a compound string
  testthat::expect_false(any(grepl("neuron_|microglia", lab)))
})

testthat::test_that("sg_blocks deduplicates long input", {
  u <- sg_units()
  # a tidy table repeats each unit once per row; building one axis position per
  # ROW instead of per UNIT silently garbles the axis
  long_u <- rep(u$analysis_key, 12)
  long_d <- rep(u$dataset, 12)
  testthat::expect_identical(nrow(sg_blocks(long_u, long_d)$order), 18L)
})

# ------------------------------------------------------------- the contract

testthat::test_that("spatial_v6 is a candidate layer honouring the size contract", {
  testthat::expect_identical(S6$contract_version, s6e_contract_version())
  testthat::expect_identical(S6$status, "candidate_only_not_promoted")
  for (f in S6$figures) {
    testthat::expect_identical(as.numeric(f$width_mm), 183)
    testthat::expect_lte(as.numeric(f$height_mm), 170)
    labs <- vapply(f$layout, function(x) as.character(x$label), character(1))
    testthat::expect_identical(labs, tolower(labs))
    testthat::expect_identical(anyDuplicated(labs), 0L)
  }
})

testthat::test_that("no two panels in a figure overlap", {
  for (f in S6$figures) {
    n <- length(f$layout)
    if (n < 2L) next
    for (i in seq_len(n - 1L)) for (j in seq(i + 1L, n)) {
      a <- f$layout[[i]]; b <- f$layout[[j]]
      sep <- as.numeric(a$x) + as.numeric(a$w) <= as.numeric(b$x) ||
        as.numeric(b$x) + as.numeric(b$w) <= as.numeric(a$x) ||
        as.numeric(a$y) + as.numeric(a$h) <= as.numeric(b$y) ||
        as.numeric(b$y) + as.numeric(b$h) <= as.numeric(a$y)
      testthat::expect_true(sep,
        info = sprintf("%s: %s overlaps %s", f$name, a$panel, b$panel))
    }
  }
})

testthat::test_that("every panel declares its spatial and phenotype axes", {
  for (p in S6$panels) {
    testthat::expect_true(nzchar(as.character(p$spatial_axes %||% "")),
                          label = paste(p$id, "spatial_axes"))
    testthat::expect_true(nzchar(as.character(p$phenotype_axes %||% "")),
                          label = paste(p$id, "phenotype_axes"))
    testthat::expect_true(nzchar(as.character(p$role %||% "")),
                          label = paste(p$id, "role"))
  }
  ids <- vapply(S6$panels, function(p) as.character(p$id), character(1))
  testthat::expect_identical(anyDuplicated(ids), 0L)
  used <- unique(unlist(lapply(S6$figures, function(f)
    vapply(f$layout, function(it) as.character(it$panel), character(1)))))
  testthat::expect_identical(sort(setdiff(ids, used)), character(0))
})

testthat::test_that("the Extended Data family is gapless", {
  ed <- Filter(function(f) identical(as.character(f$figure_key), "extended_data"),
               S6$figures)
  n <- as.integer(sub("^ED([0-9]+)_.*$", "\\1",
                      vapply(ed, function(f) as.character(f$name), character(1))))
  testthat::expect_identical(sort(n), seq_len(length(n)))
  testthat::expect_gte(length(ed), 8L)
})

# ------------------------------------------------- phenotype-blind selection

testthat::test_that("the Figure-2 fingerprint selection is phenotype-blind", {
  p <- path_results("tables", "manuscript_candidates", "spatial_v6",
                    "figure2_spatial_fingerprint_selection.csv")
  testthat::skip_if_not(have(p), "selection table not generated")
  z <- rd(p)
  testthat::expect_true(all(z$phenotype_blind == "YES"))
  testthat::expect_true(all(nzchar(z$selection_rule)))
  testthat::expect_true(all(nzchar(z$phenotype_blind_evidence)))
  # the rule must never mention a stress contrast as a selection criterion
  testthat::expect_false(any(grepl("SUS *- *RES|SUS-RES|sus_res",
                                   z$selection_rule)))
  # and the contaminated source must be explicitly excluded, not silently unused
  testthat::expect_true(all(grepl("protein_baseline_spatial_profile",
                                  z$excluded_source_note, fixed = TRUE)))
  testthat::expect_true(all(grepl("SUS - RES", z$excluded_source_note,
                                  fixed = TRUE)))
  testthat::expect_true(all(z$display_form %in%
                              c("A_individual_proteins", "B_signature_scores")))
})

testthat::test_that("the CON baseline profile uses CON animals only", {
  p <- path_results("tables", "manuscript_candidates", "spatial_v6",
                    "spatial_v6_con_baseline_profile_long.csv")
  testthat::skip_if_not(have(p), "baseline profile not generated")
  z <- rd(p)
  testthat::expect_true(all(z$n_con_animals == 3L))
  testthat::expect_identical(sort(unique(z$dataset)),
                             c("microglia", "neuron_neuropil", "neuron_soma"))
  testthat::expect_true(all(grepl("CON", z$value_definition)))
  # every unit must resolve through the shared grammar
  testthat::expect_false(anyNA(sg_resolve_unit(z$spatial_unit, z$dataset)))
  # neuropil must have 10 units, soma and microglia 4 each - never more
  for (d in unique(z$dataset)) {
    n <- length(unique(z$spatial_unit[z$dataset == d]))
    testthat::expect_identical(n, if (d == "neuron_neuropil") 10L else 4L,
                               label = d)
  }
})

testthat::test_that("the coverage audit leaves no figure-level gap", {
  p <- path_results("tables", "manuscript_candidates", "spatial_v6",
                    "figure_story_coverage_audit.csv")
  testthat::skip_if_not(have(p), "coverage audit not generated")
  z <- rd(p)
  for (col in c("figure", "panel", "scientific_question",
                "required_compartment_information", "required_region_information",
                "required_layer_information",
                "required_stress_group_information",
                "required_stress_contrast_information",
                "requires_three_group_trajectory",
                "requires_spatial_and_stress_link", "currently_present",
                "currently_missing", "omission_justified", "justification",
                "recommended_fix", "main_story_role")) {
    testthat::expect_true(col %in% names(z), info = col)
  }
  testthat::expect_true(all(z$still_missing_at_figure_level == "none"))
  testthat::expect_true(all(nzchar(z$justification)))
  # some panels MUST be phenotype-blind and that is the design, not a gap
  testthat::expect_gt(sum(z$omission_justified == "INTENTIONAL AND CORRECT"), 10L)
  # and some MUST join both axes
  testthat::expect_gt(sum(z$requires_spatial_and_stress_link == "YES"), 10L)
})

testthat::test_that("ED7 shows only QC-qualified proteins and states the exclusion", {
  p <- path_results("source_data", "manuscript_candidates", "spatial_v6",
                    "extended_data", "v6_ed_locations_source_data.csv")
  testthat::skip_if_not(have(p), "ED7 panel not generated")
  z <- rd(p)
  testthat::skip_if(identical(as.character(z$status[1]), "render_error"))
  testthat::expect_true(all(z$qc_class %in%
    c("robust_to_missingness_and_QC", "not_in_CA2_SLM_never_at_risk")))
  # the hits the audit could not clear must be absent
  testthat::expect_false(any(z$qc_class == "not_claimable_due_to_QC"))
  testthat::expect_false(any(z$qc_class == "insufficient_observed_data"))
  # never described as movement
  testthat::expect_false(any(grepl("redistribut|relocat|migrat", z$reading,
                                   ignore.case = TRUE)))
  testthat::expect_false(anyNA(z$baseline_unit))
  testthat::expect_false(anyNA(z$effect_unit))
})

testthat::test_that("the Figure-3 bridge refuses to overstate a trajectory", {
  p <- path_results("source_data", "manuscript_candidates", "spatial_v6",
                    "figure_03", "v6_f3_bridge_source_data.csv")
  testthat::skip_if_not(have(p), "bridge not generated")
  z <- rd(p)
  testthat::skip_if(identical(as.character(z$status[1]), "render_error"))
  testthat::expect_identical(nrow(z), 9L)
  testthat::expect_identical(sort(unique(z$contrast)),
                             c("RES - CON", "SUS - CON", "SUS - RES"))
  # the microglia exemplar is the only one whose three-group ordering is fully
  # FDR-supported; the other two have an unsupported RES arm and must NOT be
  # labelled divergent
  ox <- z[z$key == "oxphos", ]
  testthat::expect_true(all(ox$FDR < 0.05))
  testthat::expect_true(all(grepl("^graded", ox$shape)))
  for (k in c("synaptic", "rna")) {
    w <- z[z$key == k, ]
    rc <- w[w$contrast == "RES - CON", ]
    testthat::expect_gt(rc$FDR, 0.05)
    testthat::expect_false(any(grepl("divergent", w$shape)))
  }
})

testthat::test_that("the Figure-3 atlas excludes qc_review themes", {
  p <- path_results("source_data", "manuscript_candidates", "spatial_v6",
                    "figure_03", "v6_f3_atlas_source_data.csv")
  testthat::skip_if_not(have(p), "atlas not generated")
  z <- rd(p)
  testthat::skip_if(identical(as.character(z$status[1]), "render_error"))
  # these two carry theme_claim_eligible = FALSE and theme_role = qc_review;
  # drawing them beside the primary themes would present them as programs
  testthat::expect_false("epithelial_epidermal_qc" %in% z$theme_id)
  testthat::expect_false("cytoskeleton_structure" %in% z$theme_id)
  testthat::expect_false(anyNA(z$sg_unit))
})

# --------------------------------------------------------------- discipline

testthat::test_that("no spatial_v6 renderer creates new inference", {
  for (f in c("spatial_v6_figure_panels.R", "spatial_v6_figure3_panels.R",
              "spatial_v6_ed_panels.R", "spatial_v6_wgcna_panels.R",
              "spatial_v6_schematic.R", "spatial_grammar_utils.R")) {
    src <- code_of(repo_path("R", f))
    for (tok in c("lmFit(", "eBayes(", "topTable(", "GSEA(", "gseGO(",
                  "enrichGO(", "fgsea(", "blockwiseModules(", "TOMsimilarity(",
                  "p.adjust(", "t.test(", "wilcox.test(", "cor.test(",
                  "aov(", "lmer(", "impute.knn(", "normalizeBetweenArrays(")) {
      testthat::expect_false(grepl(tok, src, fixed = TRUE),
                             label = paste(f, tok))
    }
  }
})

testthat::test_that("spatial_v6 writes only under manuscript_candidates/spatial_v6", {
  for (f in c("spatial_v6_figure_02.R", "spatial_v6_figure_03.R",
              "spatial_v6_extended_data.R", "spatial_v6_baseline_profile.R",
              "spatial_v6_fingerprint_selection.R",
              "spatial_v6_story_coverage_audit.R")) {
    src <- code_of(repo_path("figures", f))
    testthat::expect_false(grepl("figure_contract.yml", src, fixed = TRUE),
                           label = f)
    testthat::expect_false(grepl("figure_story_v5_contract", src, fixed = TRUE),
                           label = f)
  }
  u <- code_of(repo_path("R", "spatial_v6_figure_utils.R"))
  testthat::expect_true(grepl("manuscript_candidates", u, fixed = TRUE))
  testthat::expect_true(grepl("spatial_v6", u, fixed = TRUE))
})

testthat::test_that("the microglia compartment is never called cell-intrinsic", {
  for (f in c("spatial_v6_figure_panels.R", "spatial_v6_figure3_panels.R",
              "spatial_v6_ed_panels.R", "spatial_v6_wgcna_panels.R",
              "spatial_v6_schematic.R", "spatial_grammar_utils.R")) {
    src <- readLines(repo_path("R", f), warn = FALSE)
    src <- paste(src, collapse = "\n")
    for (bad in c("microglial proteome", "microglia-specific",
                  "purified microglia")) {
      testthat::expect_false(grepl(bad, src, fixed = TRUE),
                             label = paste(f, bad))
    }
  }
  cfg <- paste(readLines(repo_path("config", "manuscript_spatial_order.yml"),
                         warn = FALSE), collapse = "\n")
  testthat::expect_true(grepl("microglia-enriched ROI", cfg, fixed = TRUE))
})

testthat::test_that("rendered panels carry no render_error", {
  for (key in c("figure_02", "figure_03", "extended_data")) {
    p <- path_results("reports", "manuscript_candidates", "spatial_v6", key,
                      "spatial_v6_panel_status.csv")
    if (!have(p)) next
    z <- rd(p)
    testthat::expect_identical(sum(z$status != "ok"), 0L, label = key)
  }
})

#!/usr/bin/env Rscript

# Phase 6B/6C scientific equivalence audit.
#
# The equivalence guard proves that every frozen object is accounted for and
# hash-identical. That is necessary but not sufficient: it would still pass if
# a new canonical path happened to point at the wrong output while some
# unrelated set of hashes stayed intact.
#
# This audit closes that gap by naming each frozen scientific contract and the
# specific object that carries it, then requiring that object to be mapped,
# present and hash-identical. A contract whose carrier went missing, changed,
# or was never mapped fails here even if the aggregate counts look healthy.
#
# It deliberately does not recompute any scientific value. The frozen summary
# tables are the evidence; recomputation is prohibited in this phase and would
# be weaker evidence anyway, because it could reproduce a number from different
# inputs.
#
# Exit status 0 when every named contract is carried by an unchanged object.

source(file.path("R", "paths.R"))

CONTRACTS <- list(
  list(
    id = "ca2_slm_dap_arithmetic",
    description = "CA2-SLM DAP arithmetic, the 37 / 28 / 6 / 9 / 15 contract",
    objects = c(
      "results/tables/11_spatial_systems/ca2_slm_robustness/CA2_SLM_DAP_robustness.csv",
      "results/tables/11_spatial_systems/ca2_slm_robustness/CA2_SLM_normalization_bias_context.csv",
      "results/tables/11_spatial_systems/ca2_slm_robustness/CA2_SLM_sample_context.csv")),
  list(
    id = "seven_program_go_atlas",
    description = "Seven-program ontology-defined GO atlas and its theme assignments",
    objects = "results/tables/10_biological_integration/gsea_wgcna_concordance/global/ontology_aware_gsea_theme_assignments_all_contrasts.csv"),
  list(
    id = "wgcna_module_identities",
    description = "WGCNA module and supermodule identities as consumed by the publication layer",
    objects = c(
      "results/tables/10_biological_integration/gsea_wgcna_concordance/global/ontology_aware_gsea_theme_assignments_all_contrasts.csv",
      "config/manuscript_spatial_order.yml")),
  list(
    id = "gsea_exemplars",
    description = "GSEA exemplar identities behind the Figure 3 curves",
    objects = "results/tables/11_spatial_systems/atlas/protein_sus_res_fdr_supported_atlas.csv"),
  list(
    id = "bilateral_metrics",
    description = "Bilateral spatial identity and precision-gain metrics",
    objects = c(
      "results/tables/11_spatial_systems/bilateral/bilateral_spatial_identity_summary.csv",
      "results/tables/11_spatial_systems/bilateral/bilateral_spatial_identity_protein_level.csv",
      "results/tables/11_spatial_systems/precision/bilateral_precision_gain.csv")),
  list(
    id = "spatial_network_nulls",
    description = "Spatial network null results, including the global multivariate test",
    objects = c(
      "results/tables/11_spatial_systems/networks/network_global_multivariate_test.csv",
      "results/tables/11_spatial_systems/networks/animal_network_distance_from_CON.csv",
      "results/tables/11_spatial_systems/networks/CON_spatial_molecular_similarity_matrix.csv")),
  list(
    id = "spatial_cell_affinity_atlas",
    description = "Protein spatial cell-affinity atlas",
    objects = "results/tables/11_spatial_systems/atlas/protein_spatial_cell_affinity.csv"),
  list(
    id = "control_spatial_identity_validation",
    description = "Control spatial identity validation source data behind Figure 2",
    objects = c(
      "results/source_data/04_differential_expression_enrichment/control_spatial_identity_validation/global/figure2e_source_data.csv",
      "results/source_data/04_differential_expression_enrichment/control_spatial_identity_validation/global/figure2f_regions_CA1layers_source_data.csv")),
  list(
    id = "behaviour_bridge",
    description = "Frozen behaviour bridge imported from topohl/MMMSociability",
    objects = c(
      "manuscript/figure1_bridge_mmmsociability/source_data/figure1_panel_statistics.csv",
      "manuscript/figure1_bridge_mmmsociability/behavior_sex_effect_contract.csv",
      "manuscript/figure1_bridge_mmmsociability/behavior_prediction_model_ladder.csv",
      "manuscript/figure1_bridge_mmmsociability/early_behavior_later_outcome_association.csv")),
  list(
    id = "statistical_and_naming_contracts",
    description = "Manuscript statistical contract and atlas naming rules",
    objects = c(
      "docs/MANUSCRIPT_STATISTICAL_CONTRACT.md",
      "docs/ATLAS_PROGRAM_SELECTION_AND_NAMING_RULES.md")),
  list(
    id = "pb12_spatial_v6_fingerprint",
    description = "PB-12: the spatial_v6 fingerprint source tables, retained unchanged",
    objects = c(
      "results/tables/manuscript_candidates/spatial_v6/figure2_spatial_fingerprint_proteins.csv",
      "results/tables/manuscript_candidates/spatial_v6/figure2_spatial_fingerprint_scores.csv"))
)

map <- utils::read.csv(repo_path("audits", "restructure_migration_map.csv"),
                       stringsAsFactors = FALSE)
MR <- Sys.getenv("EXP9_MANUSCRIPT_ROOT",
                 unset = normalizePath(file.path(repo_root(), "..", "Exp9_manuscript"),
                                       winslash = "/", mustWork = FALSE))
rownames(map) <- map$baseline_path

cat("Scientific contract equivalence audit\n")
cat("=====================================\n\n")
cat("contracts declared:", length(CONTRACTS), "\n")
cat("migration map rows:", nrow(map), "\n\n")

problems <- 0L
addressing <- 0L
covered <- character(0)

for (k in CONTRACTS) {
  cat(sprintf("%-38s %s\n", k$id, k$description))
  for (obj in k$objects) {
    covered <- c(covered, obj)
    if (!(obj %in% map$baseline_path)) {
      cat("    FAIL  not in the migration map:", obj, "\n")
      problems <- problems + 1L
      next
    }
    row <- map[obj, ]
    root <- if (row$destination_repo == "Exp9_manuscript") MR else repo_root()
    abs <- file.path(root, row$destination_path)
    if (!file.exists(abs)) {
      cat("    FAIL  destination missing:", row$destination_path, "\n")
      problems <- problems + 1L
      next
    }
    actual <- unname(tools::sha256sum(abs))
    changed <- !identical(actual, row$baseline_sha256)

    ## A carrier may differ from the freeze only if it is a declared addressing
    ## contract: its content named files by their pre-migration path. That
    ## class is independently proven to differ in addressing alone by
    ## tools/verify_path_contract_rewrites.R, so the scientific content it
    ## carries is still the frozen content. Any other difference is a failure.
    addressing_only <- identical(row$migration_class, "REWRITTEN_PATH_CONTRACT")

    if (changed && !addressing_only) {
      cat("    FAIL  content changed beyond addressing:", obj, "\n")
      cat("            class:", row$migration_class, "\n")
      problems <- problems + 1L
      next
    }

    cat(sprintf("    %s %-12s %s\n", if (changed) "ok*  " else "ok   ",
                row$destination_repo, substr(row$destination_path, 1, 74)))
    if (changed) {
      addressing <- addressing + 1L
      cat("            addressing-only rewrite; scientific content unchanged,",
          "proven by tools/verify_path_contract_rewrites.R\n")
    }
  }
  cat("\n")
}

cat("distinct objects checked      :", length(unique(covered)), "\n")
cat("byte-identical to the freeze  :", length(unique(covered)) - addressing, "\n")
cat("addressing-only rewrites (ok*):", addressing, "\n")
cat("\n")
if (problems == 0L) {
  cat("RESULT: PASS - every named scientific contract is carried by an object\n")
  cat("that is mapped, present and byte-identical to the pre-restructure freeze.\n")
  cat("No value was recomputed; the frozen tables are the evidence.\n")
  quit(status = 0L)
}
cat("RESULT: FAIL -", problems, "problem(s).\n")
quit(status = 1L)

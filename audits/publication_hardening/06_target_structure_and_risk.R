#!/usr/bin/env Rscript

# Part B, sections 23-28: the target structure, the migration risk register,
# the P0 implementation list, the legacy guards and the freeze protection.
#
# PROPOSAL ONLY for 23-25. Nothing is moved. Every proposed move is expressed
# as a from/to pair with a risk class and the exact thing that would break, so
# the migration can be executed later, in pieces, against a written contract.

setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")
source("audits/publication_hardening/00_checkpoint.R")
rd <- function(f) utils::read.csv(file.path(PH_TAB, f), stringsAsFactors = FALSE)
inv <- rd("repository_architecture_inventory.csv")
ap <- rd("repository_anti_patterns.csv")
om <- rd("repository_output_model.csv")
lay <- rd("publication_layer_review.csv")
fam <- rd("helper_library_families.csv")

# =============================================== B23/B24 target structure
#
# The target is deliberately close to what already exists. The repository's
# numbered-stage spine works; what it lacks is an explicit home for the audit
# layer, an explicit split between current and superseded figure layers, and a
# statement of which output roots are contractual.
T <- function(area, current, target, rationale, changes_paths)
  data.frame(area = area, current_state = current, target_state = target,
             rationale = rationale, changes_canonical_paths = changes_paths,
             stringsAsFactors = FALSE)
target <- rbind(
  T("numbered analysis stages",
    "00_setup .. 11_spatial_systems, 122 scripts, 141 registered steps",
    "unchanged",
    "the stage spine is the one part of the tree that is already unambiguous, and every canonical path and pipeline ID is expressed in it",
    FALSE),
  T("duplicate interpretation stage",
    "08_biological_interpretation (1 script) alongside 10_biological_integration (11 scripts)",
    "fold archive/08_integration/01_compartment_fidelity_summary.R into 03_qc_exploration as 04f_, or register it as a step of 10_biological_integration",
    "a one-script stage whose name is a near-synonym of another stage forces the reader to guess which is current; it also emits three files whose names duplicate 04d's",
    TRUE),
  T("figure / manuscript layer",
    "figures/ holds 53 scripts across 7 generations, distinguished only by filename prefix",
    "figures/current/ for the v9 layer and figures/superseded/<layer>/ for the rest, with renderer prefixes unchanged",
    "the frozen figures were produced by the code as it stood, so superseded generations must remain readable; they should not remain equally prominent",
    TRUE),
  T("shared helper library",
    "R/ holds 100 helpers in 8 implicit families, flat",
    "R/ subdirectories per family (infrastructure, wgcna, spatial, figure layers, analysis helpers)",
    "the families already exist in the filenames; making them directories removes the need to know the prefix convention to navigate",
    TRUE),
  T("audit layer",
    "99_audits/ holds part29, program_evidence, publication_hardening; excluded from the registry by design",
    "unchanged, with the exclusion stated in R/utilities/pipeline_registry.R as it is now",
    "audits are not pipeline stages and must not become required steps; DEC-004",
    FALSE),
  T("deprecated and scaffold code",
    "99_deprecated/ (11) and 90_testing/ (11), zero inbound edges from producer layers",
    "unchanged, with a machine-checked guard that keeps inbound edges at zero",
    "the tree is already clean here; what is missing is the test that keeps it clean",
    FALSE),
  T("repository root",
    "run_dataset_pipeline.R and audits/wgcna/proteomics_wgcna_downstream_audit.R sit at the root",
    "run_dataset_pipeline.R stays (it is the entrypoint); audits/wgcna/proteomics_wgcna_downstream_audit.R moves to 99_audits/",
    "an entrypoint belongs at the root; a one-off audit does not, and the audit layer already exists for it",
    TRUE),
  T("output model",
    "results/ holds 9 structured roots plus 2 ad-hoc EWCE comparison roots",
    "results/<structured roots> only; the two EWCE roots move under results/audit/",
    "both ad-hoc roots are untracked and regenerated, so the move costs nothing and removes two top-level names that look canonical but are not",
    FALSE))
utils::write.csv(target, file.path(PH_TAB, "target_repository_structure.csv"),
                 row.names = FALSE)

# =================================================== B25 migration risk register
#
# Risk is a function of what reads the path, not of how many files move. A move
# that no registry, manifest, test or frozen artefact names is P0_SAFE_NOW; a
# move that a canonical path names is DEFERRED regardless of how tidy it is.
R <- function(item, from, to, risk, breaks, precondition, verify)
  data.frame(item = item, move_from = from, move_to = to, risk_class = risk,
             what_would_break = breaks, precondition = precondition,
             verification = verify, stringsAsFactors = FALSE)
risk <- rbind(
  R("ad-hoc EWCE output roots",
    "results/EWCE_sample_vs_animal_COMPARISON, results/EWCE_sample_vs_animal_REPAIRED",
    "results/audit/", "P1_SAFE_WITH_TESTS",
    "nothing in the repository - both are untracked and named by no registry, manifest or test - but the move is a physical directory migration, which DEC-002 excludes from this pass",
    "explicit approval to migrate output directories",
    "git status stays clean; test suite unchanged"),
  R("specificity gate regex", "figures/final_truth_v9_semantics.R S9 pattern",
    "widened from selectively|exclusively to selectiv|exclusiv", "P0_SAFE_NOW",
    "nothing - the widened class is P1 overclaim, which does not stop the build; only P0 does",
    "none", "grep the pattern; rerun of the v9 semantics layer will list the two prose hits as P1"),
  R("selective rulebook entry", "figures/final_truth_v9_semantics.R S28 rules",
    "new RULE(\"selective\", ...) after susceptibility-specific", "P0_SAFE_NOW",
    "nothing - the rules file is regenerated from this block",
    "none", "regenerated manuscript_semantic_rules.md gains one section"),
  R("one-off root audit script", "audits/wgcna/proteomics_wgcna_downstream_audit.R",
    "99_audits/", "P1_SAFE_WITH_TESTS",
    "any hard-coded relative source() inside it, and RUN_ORDER.md if it names it",
    "confirm it is in no registry step and no RUN_ORDER entry",
    "Rscript the moved file; run test-pipeline-registry.R"),
  R("duplicate interpretation stage",
    "archive/08_integration/01_compartment_fidelity_summary.R",
    "03_qc_exploration/04f_compartment_fidelity_summary.r",
    "P2_DEFERRED",
    "its output directory is derived from the stage name, so every file it writes changes path; downstream readers of those paths break",
    "identify every consumer of results/.../08_biological_interpretation/",
    "full pipeline rerun for the affected datasets"),
  R("figure layer split", "figures/<layer>_*.R",
    "figures/current/ and figures/superseded/<layer>/", "P2_DEFERRED",
    "pipeline.yml names figure scripts by path; the publication freeze manifest hashes them; RUN_ORDER.md lists them",
    "registry, manifest and RUN_ORDER updated in the same commit as the move",
    "pipeline registry validation, freeze manifest reverification, full figure rebuild byte-compared against the frozen SVGs"),
  R("helper family subdirectories", "R/*.R",
    "R/<family>/*.R", "P2_DEFERRED",
    "every source(repo_path(\"R\", \"x.R\")) in 445 scripts, and the one-definition tests that glob R/*.R",
    "mechanical rewrite of all source() calls plus the R/ globs in the test suite",
    "full test suite; dependency-edge count must be unchanged at 839"),
  R("unreferenced helper", "R/statistics/module_stats.R", "retain in place",
    "NO_CHANGE",
    "nothing - but it is referenced only by R/README.md, so deleting it would be the tidy-driven deletion the audit is told not to make",
    "none", "recorded as PH-003; revisit only if a rewrite needs the namespace"),
  R("five out-of-layer v9 renderers",
    "nf_pca_compact, nf_bilateral_main, nvp_ed_celltype, s5_ed_ca2_displacement, s5_ed_network_distance",
    "retain in place; document the dependency", "NO_CHANGE",
    "copying them into the v9 layer would change the function bodies that produced the frozen SVGs, so the frozen figures would no longer be reproducible from the code that made them",
    "none", "recorded as PH-001 and declared in the architecture document"))
utils::write.csv(risk, file.path(PH_TAB, "migration_risk_register.csv"),
                 row.names = FALSE)

# =========================================================== B26 P0 implementation
#
# DEC-002 excludes physical directory migration from this pass, so every P0 item
# here is a source-level change inside a file this audit is allowed to touch.
# Each is verified against the file on disk rather than asserted.
sem <- readLines("figures/final_truth_v9_semantics.R", warn = FALSE)
impl <- data.frame(
  item = c("specificity gate regex", "selective rulebook entry",
           "ad-hoc EWCE output roots", "five out-of-layer v9 renderers"),
  action = c(
    if (any(grepl('"selectiv\\|exclusiv"', sem))) "IMPLEMENTED" else "NOT APPLIED",
    if (any(grepl('RULE\\("selective"', sem))) "IMPLEMENTED" else "NOT APPLIED",
    "DEFERRED - DEC-002", "NO CHANGE - by design"),
  detail = c(
    "S9 pattern now matches the adjective; PH-002 can no longer pass the gate",
    "manuscript_semantic_rules.md gains a selective section on next regeneration",
    "physical output-directory migration is out of scope for this pass",
    "copying them into the v9 layer would change the code that produced the frozen SVGs"),
  requires_rerun_to_take_effect = c(TRUE, TRUE, FALSE, FALSE),
  stringsAsFactors = FALSE)
utils::write.csv(impl, file.path(PH_TAB, "p0_implementation_log.csv"),
                 row.names = FALSE)

cat("\n===== PART B TARGET / RISK =====\n")
cat("B23/B24 target areas:", nrow(target),
    "| changing canonical paths:", sum(target$changes_canonical_paths), "\n")
print(target[, c("area", "changes_canonical_paths")])
cat("\nB25 risk register:\n")
print(risk[, c("item", "risk_class")])
cat("\nB26 P0 implementation:\n"); print(impl)

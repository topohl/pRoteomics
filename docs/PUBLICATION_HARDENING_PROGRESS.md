# Publication Hardening Progress

## Run identity

- **Starting commit:** b26cd6f (Codify the atlas selection and naming rules, and rename the Chromatin row)
- **Date:** 2026-09-15
- **Audit version:** publication_hardening_v1
- **Local repo authoritative:** yes
- **Status:** all 28 sections COMPLETE in run 1

## Scope

A manuscript-wide statistical/terminology consistency audit (Part A) and a
repository architecture review (Part B).

**May change:** documentation, audit artefacts, registries of entrypoints,
deprecation banners, validation tests, low-risk publication-hardening fixes.

**May NOT change:** frozen atlas biology, canonical statistics, source-data
values, figure bytes (except to correct a factual defect), canonical paths,
run order, pipeline IDs. No physical repository migration in this pass.

## Completed sections

Machine state: `results/reports/publication_hardening/checkpoint_state.csv`
(28 COMPLETE, 0 DEFERRED).

| # | Section | Output |
|---|---|---|
| A01 | Manuscript-reachable output inventory | `manuscript_reachability_inventory.csv` (42 panels) |
| A02 | Effect / contrast contract | `manuscript_effect_contract.csv` |
| A03 | Statistic identity | `manuscript_statistic_identity_audit.csv` (10 quantities, 0 FAIL) |
| A04 | Biological n / replication | `manuscript_n_replication_audit.csv` |
| A05 | Multiple-testing contract | `manuscript_multiple_testing_audit.csv` (44 mentions, 0 FAIL) |
| A06 | Null-language audit | `manuscript_claim_language_audit.csv` |
| A07 | Interaction / specificity claims | same |
| A08 | Claim-language audit | same (55 hits, 2 genuine) |
| A09 | Dataset-specific interpretation rules | `manuscript_dataset_interpretation_rules.csv`, `manuscript_specificity_gate_defect.csv` |
| A10 | Zero / NA / not-evaluable | `manuscript_zero_na_audit.csv` (39 files, 0 REVIEW) |
| A11 | Numerical-floor disclosure | verified present in corpus |
| A12 | Methods contract table | `manuscript_methods_contract.csv` (11 analyses) |
| B13 | Architecture inventory | `repository_architecture_inventory.csv` (445 scripts) |
| B14 | Active dependency graph | `repository_dependency_edges.csv`, `figures/repository_dependency_graph.svg` (839 edges, 0 unresolved) |
| B15 | Structural anti-patterns | `repository_anti_patterns.csv` |
| B16 | Output-model review | `repository_output_model.csv` |
| B17 | Publication layer review | `publication_layer_review.csv` |
| B18 | Code organization review | `helper_library_families.csv` |
| B19 | Script naming audit | `repository_script_naming_audit.csv` |
| B20 | Identity vs display label | `identity_display_label_contract.csv` (7/7 PASS) |
| B21 | One authoritative script rule | `repository_contested_outputs.csv` (0 contested paths) |
| B22 | Source-data-first contract | `source_data_first_contract.csv` (3/4 PASS) |
| B23 | Target repository structure | `target_repository_structure.csv` |
| B24 | Output target structure | same + `repository_output_model.csv` |
| B25 | Migration risk classification | `migration_risk_register.csv` |
| B26 | P0 implementation | `p0_implementation_log.csv` |
| B27 | Legacy guards | `tests/testthat/test-publication-hardening.R` |
| B28 | Publication freeze protection | same |

Documents produced: `docs/MANUSCRIPT_STATISTICAL_CONTRACT.md`,
`docs/REPOSITORY_ARCHITECTURE.md`, `docs/CANONICAL_ANALYSIS_ENTRYPOINTS.md`.

## Findings ledger

Findings are numbered PH-001 onward and are **never renumbered**.

### PH-001 — five v9 publication panels are rendered from superseded layers
**Severity:** P1 reachability. **Disposition:** ACCEPTED, guarded.

Of the 42 panels in the frozen v9 contract, 37 are rendered by `f9_`/`s9f_`
functions. Five are not:

| Figure | Panel | Renderer | Home layer |
|---|---|---|---|
| F2_NATURE_FINAL_V9 | c | `nf_pca_compact` | nature_final_v7 |
| F2_NATURE_FINAL_V9 | f | `nf_bilateral_main` | nature_final_v7 |
| ED3_FINAL_V9 | c | `s5_ed_ca2_displacement` | story_v5 |
| ED_WGCNA_FINAL_V9 | c | `nvp_ed_celltype` | nature_v2 |
| ED8_FINAL_V9 | b | `s5_ed_network_distance` | story_v5 |

Not repaired: the frozen SVGs were produced by these exact function bodies, so
copying them into the v9 layer would mean the released figures are no longer
reproducible from the code that made them. Two tests now fix the set at exactly
these five and assert every contract renderer stays defined.

### PH-002 — the specificity gate matched the adverb but not the adjective
**Severity:** P1 overclaim. **Disposition:** gate FIXED; prose flagged, not rewritten.

`figures/final_truth_v9_semantics.R` §S9 exists to catch specificity wording that
no test in the study licenses. Its pattern was `selectively|exclusively` — the
adverb only. The adjective passed, and it passed in the two highest-visibility
sentences in the package:

- `final_figure_story_v9.md:25` — the Figure 3 title claim: "stronger spatially
  **selective** coordinated molecular programs"
- `final_figure_story_v9.md:51` — "The defensible conclusion": "**selective**
  local molecular-program differences"

"Spatially selective" is the spatial form of the `susceptibility-specific` claim
the same rulebook already bans with *"ONLY USE WHEN: a specificity analysis was
actually performed."* The only heterogeneity test executed anywhere in the
package is the WGCNA stress × spatial-unit omnibus: **0 of 35 cells
FDR-supported, smallest FDR 0.27**. No equivalent test exists at GO-program
level.

Implemented (both P0_SAFE_NOW, source-level):
1. S9 pattern widened to `selectiv|exclusiv`. The class is P1, and only P0 stops
   the build, so this cannot break a rebuild.
2. A `selective` entry added to the S28 rulebook generator, naming the test that
   would license the word.

**Not implemented:** the two sentences themselves. Rewriting the Figure 3 title
claim is a change to the scientific story, not a hardening fix, and is the
user's call. Suggested replacement wording — "spatially restricted" or "program
differences detectable in some spatial units and not others" — states the same
observation as a count rather than as selectivity. Both fixes take effect only
when the v9 semantics layer is rerun, which this pass did not do.

### PH-003 — `R/module_stats.R` is unreferenced
**Severity:** P2. **Disposition:** RETAINED deliberately.

The only helper of 100 with no `source()` edge and no indirect reference; it is
named only in `R/README.md`. Deleting it would be exactly the tidiness-driven
deletion this audit was told not to make. (`R/renv_lock_audit.R` also has no
`source()` edge but is reached through a variable and through
`testthat::test_path` — static edge counting alone would misreport it.)

### PH-004 — `08_biological_interpretation` is a one-script near-duplicate stage
**Severity:** P2 clarity. **Disposition:** DEFERRED (P2_DEFERRED).

One script, a name that is a near-synonym of `10_biological_integration`, not a
registry step, and it emits three files whose basenames duplicate those of
`03_qc_exploration/04d_compartment_marker_fidelity.r`. The paths differ, so
nothing is overwritten — but a reader handed
`compartment_marker_fidelity_scores.csv` cannot tell which stage produced it.
Folding it in would change its output directory, which derives from the stage
name, so every path it writes would move.

### PH-005 — 20 scripts inside numbered stages are not registry steps
**Severity:** P1 reachability. **Disposition:** documented, no change.

Eight are `01_preprocessing` predecessors of the animal-level contract, five are
in explicit `legacy/` subdirectories, and the rest are audits or wrappers
(`09_export_pride_journal/RUN_EXPORT.R`). None is reachable from a publication
artefact, and no producer layer sources any of them (0 edges, now tested). The
residual risk is that each can still be run by hand while never being validated
by the registry.

### PH-006 — two ad-hoc output roots sit beside the structured ones
**Severity:** P2 clarity. **Disposition:** DEFERRED (DEC-002).

`results/EWCE_sample_vs_animal_COMPARISON` and
`results/EWCE_sample_vs_animal_REPAIRED`. Both untracked and named by no
registry, manifest or test; they belong under `results/audit/`. Moving them is a
physical directory migration, which this pass excludes.

### PH-007 — naming conventions inside numbered stages are mixed
**Severity:** P2 cosmetic. **Disposition:** no change.

19 scripts use a sub-step suffix (`02a_`, `04d_`) rather than a bare `NN_`; this
is a deliberate insertion convention. Numbered stages mix `.R` and `.r`
extensions while every other layer is internally consistent. Renaming would
churn every registry path for no functional gain.

## Verified clean

- **0** producer-layer edges into `90_testing/` or `99_deprecated/` (now tested).
- **0** unresolved `source()` targets outside the test suite.
- **0** registry steps naming a missing file.
- **0** output files written by two scripts at the same relative path.
- **0** statistic-identity violations across 10 displayed quantities.
- **0** multiple-testing statements asserting a theme-level p-value or FDR.
- **7 of 7** primary atlas themes carry exactly one display label, in registry order.
- **42 of 42** panels have released source data and a vector SVG.
- Full suite after all changes: **11,364 passing, 0 failures, 0 errors, 1 skip.**

## Decisions

- **DEC-001:** Frozen atlas biology is out of scope unless a factual correctness
  defect is found. *(No such defect was found.)*
- **DEC-002:** No physical directory migration in this pass. Structural changes
  are proposed with a migration map and risk class only.
- **DEC-003:** Only `P0_SAFE_NOW` items may be implemented automatically. Anything
  touching canonical paths, hashes, manifests or publication bytes is PROPOSED.
- **DEC-004:** `99_audits/` is excluded from the pipeline registry by design
  (established at commit bab450a) — it is an isolated audit layer, not a stage.
- **DEC-005:** PH-002's detector is fixed but its two prose sentences are not
  rewritten. Correcting an overclaim in the Figure 3 title claim changes the
  scientific story, which is the author's call, not the audit's.
- **DEC-006:** The five out-of-layer renderers (PH-001) are frozen in place, not
  migrated. Reproducibility of the released figures outranks layer tidiness.

## Deferred architecture changes

See `results/tables/publication_hardening/migration_risk_register.csv` and §5 of
`docs/REPOSITORY_ARCHITECTURE.md`. Nothing on that list is blocking; each entry
carries its precondition and its verification step.

## Open questions

1. **PH-002 prose.** Rewrite the two "selective" sentences, or keep them and
   accept the P1 flag? Requires an author decision, then a rerun of
   `figures/final_truth_v9_semantics.R`.
2. **Rerun timing.** Both PH-002 fixes are source-level and take effect only when
   the v9 semantics layer is regenerated. That rerun rewrites frozen report
   bytes and was therefore not performed here.

## Resume point

**LAST COMPLETED SECTION:** B28 — all 28 sections COMPLETE.

**NEXT EXACT ACTION:** none within this audit. The two follow-ups are the open
questions above, both of which need an author decision first.

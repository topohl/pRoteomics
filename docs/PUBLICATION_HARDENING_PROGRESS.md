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
**Severity:** P1 overclaim. **Disposition:** **RESOLVED** (gate fixed at b392977;
prose rewritten and the layer rebuilt in this pass).

`figures/final_truth_v9_semantics.R` §S9 exists to catch specificity wording that
no test in the study licenses. Its pattern was `selectively|exclusively` — the
adverb only. The adjective passed, and it passed in the two highest-visibility
sentences in the package.

**The licensing test does not exist.** Independently recomputed from the result
tables during this pass: every `interaction_omnibus` FDR family in the repository
is WGCNA — `{neuron_neuropil, neuron_soma, microglia}` × `{module,
higher_order_multimodule}`, six families and no others. Module level 15+13+7 =
**35 tests, 0 below 0.05, smallest FDR 0.2741**; supermodule level 4+3+2 = **9
tests, 0 below 0.05, smallest FDR 0.1924**. There is no GO-program-level
equivalent, which is the level both sentences describe.

**Sentence A — the Figure 3 title claim.** Source
`figures/final_truth_v9_vector_audit.R:123-125`.

> **OLD:** Later resilient and susceptible outcomes are associated with sparse
> protein-level effects but **stronger spatially selective** coordinated molecular
> programs across distinct hippocampal compartments.

> **NEW:** Later resilient and susceptible outcomes are associated with sparse
> protein-level effects **and with spatially resolved** coordinated
> **molecular-program differences** across distinct hippocampal compartments.

Two unlicensed words, not one. `stronger` is a cross-family comparison: protein-
and program-level results sit in separate BH families that were never placed on a
common scale. The repository had **already** excised the identical comparison —
`claim_chain_audit.md` rejects the synonym `richer` as "a loose comparison between
two different statistical objects", `core_story_corrected.md` states "neither is
claimed to be the stronger", and lists "that program-level evidence is stronger
than protein-level evidence" under what it deliberately does not say. The Figure 3
title was the one artefact still asserting it. Removing it **resolves** an existing
internal contradiction rather than creating one. `but` became `and with` for the
same reason: the adversative carries the ranking that the conjunction does not.

**Sentence B — the defensible conclusion.** Source
`figures/final_truth_v9_vector_audit.R:150-153`.

> **OLD:** Stress outcome is associated with **selective local**
> molecular-program differences superimposed on a hippocampal spatial molecular
> architecture that remained evident across groups: bilaterally reproducible, and
> with no whole-network group difference detectable at three animals per group.

> **NEW:** Stress outcome is associated with **coordinated** molecular-program
> differences **across hippocampal spatial contexts**, superimposed on a spatial
> molecular architecture that remained evident across groups: bilaterally
> reproducible, and with no whole-network group difference detectable at three
> animals per group.

`local` was dropped with `selective`: it carries the same restriction claim.
`coordinated` is retained — the claim chain licenses it, because coordination is
what a ranked gene-set statistic measures. The second `hippocampal` was dropped to
avoid repeating it within one sentence; the sentence still names the hippocampus.

`spatially resolved` was chosen over the rulebook's other sanctioned substitutes
because it describes the analytical resolution — 18 prespecified spatial units —
rather than asserting a detection boundary.

**Verification.** Semantic scan after rebuild: **0 unresolved (0 P1, 0 P2, 0 P0)**.
Specificity gate rows for the selective class: **0**. No `selectiv` string remains
anywhere in the reader-facing corpus outside the two files that exist to quote
banned wording. Delta audit: **3 files changed, all TEXT_EXPECTED**; 29,131 numeric
values compared across 455 columns, **0 changed**; 29 canonical statistic files
hashed, only the edited generator differs. No SVG or PDF changed, and vector
integrity held at 51 of 51.

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

### PH-008 — "spatially restricted" describes the central claim and was reviewed, not changed
**Severity:** P2 / open question. **Disposition:** REVIEWED, deliberately unchanged.

Raised during the PH-002 adjudication. `spatially restricted` describes the
CENTRAL claim in two further reader-facing artefacts —
`core_story_corrected.md:18` and `nature_reviewer_vulnerabilities.md:112` — and is
the spatial-restriction cousin of the wording PH-002 removed.

It was left alone because it does **not** violate the existing statistical
contract; it is the contract's own prescribed alternative:

- `final_truth_v9_semantics.R` `RULE("selective")` lists "spatially restricted" in
  the **USE** column, as a sanctioned substitute for the banned word.
- The claim chain certifies it explicitly, with its evidence named:
  `CS(2, "spatially restricted", TRUE, "the theme atlas and ED6: FDR-supported
  terms occur in a subset of units", "keep")`.

**The argument against it, recorded for an author decision.** At three animals per
group, what is restricted is the *detection*, not the effect; absence of FDR
support in a unit is not absence of an effect there. Two independent reviewers
argued it should be reclassified PROHIBITED_UNSUPPORTED on that basis. Two others
argued the opposite, on the contract grounds above. Changing it would mean
amending `docs/MANUSCRIPT_STATISTICAL_CONTRACT.md` and the rulebook, not just
rewording a sentence — a contract change, not a hardening fix.

**Note on vocabulary.** The Figure 3 title now says "spatially resolved" while the
corrected core story says "spatially restricted". Both are licensed and they are
compatible — one describes measurement resolution, the other describes where
support was found — but the package now uses two vocabularies for one claim.
Harmonising them is the natural follow-up if PH-008 is taken up.

### PH-009 — "specificity comparisons" in the external-validation legend
**Severity:** P2 / unverified. **Disposition:** RECORDED, not changed, not verified.

`final_figure_legends_v9.md:61` and `figure_final_truth_v9_contract.yml:118,286`
use "specificity"/"specificity comparisons" for the external-validation inventory.
Three independent reviewers argued this must change, on the grounds that the
phrase names a matched-versus-mismatched discrimination that no statistic in the
repository computes, and one asserted that the plotted `p_adjust` is an
uncorrected single-set p rather than a BH value within the signature families.

Not changed, and deliberately **not verified**, for three reasons: it concerns
external-anatomical-validation specificity rather than the spatial selectivity
PH-002 addresses; confirming it would be a new statistical audit, which this pass
is scoped out of; and two of the three locations are in the frozen figure
contract. If the p-value claim is correct it is a genuine defect and would need
its own pass.

### PH-010 — a duplicate `%||%` definition was masked by the pre-commit test run
**Severity:** P1 correctness (test-visible). **Disposition:** FIXED in this pass.

`99_audits/publication_hardening/01_manuscript_contract.R:66` defined
``` `%||%` <- function(a, b) ... ```, which violates the repository's rule that
`R/null_coalescing.R` holds the only definition
(`tests/testthat/test-null-coalescing.R`).

**Why it was not caught at b392977.** The test enumerates candidate files with
`git ls-files`. The full suite for b392977 was run *before* the commit, while
`01_manuscript_contract.R` was still untracked — so the test could not see the
file and the run reported 11,364 passing. The commit then made the file tracked,
and the defect became visible on the next run. **The b392977 benchmark was
therefore measured on a tree that did not match what was committed, and the
committed state was already failing this test.**

Fixed by deleting the redundant definition rather than by sourcing
`R/null_coalescing.R`: `%||%` has been in base R since 4.4, the two other audit
scripts in this directory already rely on that, and the definition sat *after*
its first use in the file, so it was never doing anything.

**Process lesson.** A benchmark taken before `git add` is not a benchmark of the
commit. Any future suite run used as a baseline must be taken against the staged
or committed tree, not the working tree.

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
- **DEC-005:** ~~PH-002's detector is fixed but its two prose sentences are not
  rewritten.~~ **Superseded by DEC-007.** The author authorised the rewrite; both
  sentences were corrected and the layer rebuilt.
- **DEC-006:** The five out-of-layer renderers (PH-001) are frozen in place, not
  migrated. Reproducibility of the released figures outranks layer tidiness.
- **DEC-007:** PH-002 is resolved by rewriting the two sentences at source and
  rerunning only `final_truth_v9_vector_audit.R` and `final_truth_v9_semantics.R`.
  No analysis was rerun and no figure re-rendered, because the wording lives in a
  generated markdown report and appears in no panel SVG.
- **DEC-008:** Wording that the repository's own contract sanctions is not changed
  by this pass, even where a reviewer argues it overclaims. `spatially restricted`
  (PH-008) is listed in the rulebook's USE column and carries a claim-chain
  certification with named evidence, so changing it would be a contract amendment,
  not a hardening fix.
- **DEC-009:** Findings outside PH-002's subject matter are recorded, not acted
  on. PH-009 concerns external-validation specificity and would require a new
  statistical audit to confirm; this pass is scoped out of that.

## Deferred architecture changes

See `results/tables/publication_hardening/migration_risk_register.csv` and §5 of
`docs/REPOSITORY_ARCHITECTURE.md`. Nothing on that list is blocking; each entry
carries its precondition and its verification step.

## Open questions

1. ~~**PH-002 prose.**~~ **Closed.** Both sentences rewritten; semantic scan
   returns 0 unresolved.
2. ~~**Rerun timing.**~~ **Closed.** The two generators were rerun; the delta
   audit confirms 3 files changed, all TEXT_EXPECTED, and 0 statistical values
   moved.
3. **PH-008 — vocabulary harmonisation.** The Figure 3 title now says "spatially
   resolved" while `core_story_corrected.md` says "spatially restricted". Both are
   licensed and compatible, but the package now uses two vocabularies for one
   claim. Harmonising needs an author decision, and dropping "restricted" would
   mean amending the statistical contract and the rulebook.
4. **PH-009 — external-validation specificity.** Three reviewers argued
   "specificity comparisons" overclaims, and one asserted the plotted `p_adjust`
   is an uncorrected single-set p. Unverified. If correct it is a real defect and
   needs its own pass.

## Resume point

**LAST COMPLETED SECTION:** B28 plus the PH-002 resolution pass.

**NEXT EXACT ACTION:** none. PH-002 is RESOLVED and verified. PH-001 is unchanged
and remains an accepted architecture/provenance condition. PH-008 and PH-009 are
recorded for an author decision and are not blocking. The next task is manuscript
drafting; no repository migration should begin from this state without reopening
the migration risk register.

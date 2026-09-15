# Manuscript drafting progress

## Run identity

- **Starting HEAD:** b339d58 (Harmonize spatial wording on "resolved" and add a
  verification-state contract)
- **Worktree at start:** clean
- **Phase:** 1 — skeleton + Figure 2 and Figure 3 Results
- **Date:** 2026-09-15

## Scope of this phase

**Drafting:** the manuscript skeleton, Results §2 (Figure 2) and Results §3
(Figure 3), with claim-by-claim provenance and a red-team review.

**Not drafting:** Introduction, Discussion, Methods, Figure 1 / behavioural
Results. Those are phase 2 and later, and are left as explicit placeholders.

**May write:** manuscript prose, provenance tables, Methods TODOs, this file.

**May NOT write:** source data, statistical outputs, figures, the atlas registry,
analysis scripts, canonical configs. If drafting reveals a possible scientific
defect, the claim is stopped and recorded — the analysis is not repaired here.

## Authoritative sources consulted

| Source | Used for |
|---|---|
| `figures/figure_final_truth_v9_contract.yml` | panel inventory, `primary_source` per panel |
| `results/reports/.../final_figure_legends_v9.md` | panel logic and the standing caveats |
| `results/reports/.../final_figure_story_v9.md` | frozen story wording |
| `results/reports/.../core_story_corrected.md` | audited preferred core story |
| `results/reports/.../claim_chain_audit.md` | clause-by-clause licensing |
| `results/reports/.../manuscript_semantic_rules.md` | 20 USE / AVOID / ONLY USE WHEN rules |
| `results/reports/.../known_issues_v9.md` | CA2-SLM QC, leading-edge FDR |
| `docs/MANUSCRIPT_STATISTICAL_CONTRACT.md` | n, families, statistic identity, §6 spatial wording |
| `docs/ATLAS_PROGRAM_SELECTION_AND_NAMING_RULES.md` | atlas selection/naming contracts |
| `docs/PUBLICATION_HARDENING_PROGRESS.md` | PH-001…PH-011 dispositions |
| `results/source_data/manuscript_candidates/final_truth_v9/**` | every number in the draft |

## Figure structure (from the frozen contract)

**F2_NATURE_FINAL_V9** — *Does the spatial proteomics experiment faithfully
resolve biologically meaningful hippocampal molecular architecture?*
a schematic · b acquisition depth · c global proteome structure · d CON-only
spatial fingerprint · e compartment marker identity · f bilateral reproducibility
· g **external** anatomical validation · h functional characterisation (**not**
validation)

**F3_NATURE_FINAL_V9** — *What molecular programs distinguish later resilient and
susceptible outcomes, where do they occur, and which proteins carry them?*
a differential-abundance burden · b theme-level atlas · c anatomical bridge ·
d/e/f ranked enrichment for the three exemplars · g/h/i selected leading-edge
proteins

## Sections completed

| Section | Status |
|---|---|
| Skeleton | **done** — `manuscript/manuscript_draft.md` |
| Results §2 (Figure 2) | **drafted and revised** — 786 words |
| Results §3 (Figure 3) | **drafted and revised** — 999 words |
| Provenance tables | **done** — 44 statement rows, 11 claim rows, 13 Methods TODOs |
| Red-team review | **done** — 18 logged, 15 real defects fixed, 3 clean |

## Red-team outcome

Seven independent reviewers checked every material sentence against the claim
contract, one per checklist question. **Fifteen real defects were found and
fixed; three dimensions came back clean.** Full log:
`manuscript/results_red_team_review.csv`.

**The most serious was mine and it was factually false.** I wrote *"The lower
tail is specific to CA1: … the three CA1 strata contrasts were the least
reproducible in the dataset."* Four reviewers flagged it independently, and
re-derivation from `bilateral_spatial_identity_summary.csv` refutes it:
CA1-SLM (r = 0.738) ranks **8th of 15, at the inventory median**, and four
non-CA1 contrasts fall below it. It also used a word reserved by contract §6 —
the same contract this project spent the previous pass enforcing.

The corrected finding is stronger and better sourced: **bilateral concordance
varies more by measurement compartment than by anatomical scale** — neuronal
soma median 0.80, neuropil 0.74, microglia-enriched 0.56. Anatomical scale is
not the discriminating factor; dentate-gyrus layer contrasts reproduce as well as
the better regional ones (r = 0.81).

Other fixes of substance: both PCA *P* values removed (F2c is declared
descriptive and their degrees of freedom come from 323 non-independent
acquisitions); "establishes"/"establish" downgraded; the adversative
"nevertheless" removed — I had reintroduced the exact construction PH-002 struck
from the Figure 3 title; "Program families differ across spatial contexts"
rewritten as a factual count (a heterogeneity assertion dodging the banned word);
"compartment-restricted" replaced; the CA2-SLM over-negation corrected (6 of 28
*did* qualify); a GSEA FDR caveat added; CAMERA scoped to the SUS-RES arm;
edge–behaviour scoped to the eight neuropil pairs; external validation scoped to
neuropil and soma with the microglia gap stated; and the unverifiable union count
5,694 removed.

**Clean on first pass:** CAMERA never described as validation; leading-edge
proteins never implied individually significant; all microglia and soma wording
correct throughout.

## Known limitation of this pass

The adjudication phase of the red-team workflow lost 51 of 59 agents to an
organisation spend limit. The seven reviews themselves completed in full, and the
load-bearing finding (RT-01) was independently re-derived by the author before
acting on it. Remaining dispositions in
`results_red_team_review.csv` are the author's adjudication of reviewer evidence,
not a second automated pass. **Open scope question:** whether the WGCNA
"0 of 45" module × contrast count is compartment-scoped (one reviewer asserts the
full set is 0 of 105, smallest FDR 0.16). The draft attributes the count to the
displayed Extended Data panel, which is defensible either way, but this should be
settled in Methods.

Evidence was extracted across nine domains and then **independently re-derived
from the cited sources: 108 checks, 0 wrong.** No number entered the draft
without a second agent reopening its source file and reproducing it.

## Findings raised during drafting

These are documentation defects found while writing, recorded and **not
repaired** — repairing them means editing generator scripts, which this phase
must not do (§31). None changes a scientific result.

### PH-012 — the OXPHOS rulebook entry describes the pre-v3 atlas
`RULE("OXPHOS")` in `manuscript_semantic_rules.md` says the theme is *"the
20-term atlas theme"* and *"also contains glycolytic terms"*. Verified against
`ontology_aware_gsea_theme_assignments_all_contrasts.csv` at
`registry_version = manuscript_go_themes_v3`: the theme has **16 GO terms and
zero glycolytic members**. 324 glycolytic GO rows exist but carry `theme_id = ""`
and `assignment_status = "unclassified"`, correctly excluded by
`GO:0006096 -> exclude_anchor_and_descendants`. PDH (GO:0006086) and TCA
(GO:0006099) are retained. `final_terminology_registry.csv` carries the same
stale claim in its `notes` field while its own `canonical_source` field
correctly says 16 — the row is internally inconsistent.
**Effect on drafting:** following the rulebook would have put glycolysis into the
manuscript. The draft uses the verified membership instead.

### PH-013 — the eps-floor disclosure count is stale
Legends state the floor applies to **"90 of the 851"** displayed FDR-supported
occurrences. Current source data gives **94 of 953** (RES-CON 13, SUS-CON 54,
SUS-RES 27). The draft does not quote either figure; `MT-09` requires the
current count in Methods.

### PH-014 — the synaptic exemplar rank "51st of 60" matches no computed scope
Hardcoded in `figures/final_truth_v9_semantics.R`. Computed values: 52nd of 65
(audit table), 51st of 64 unique FDR-supported terms after deduplication, 94th of
195, or 175th of 237 depending on scope. The draft avoids the number and says
only that the exemplars "are not the strongest result in their own units", which
is true under every scope.

### PH-015 — `p_adjust` in the external-validation tables is not adjusted
In `v9_ed_external_full_source_data.csv` the column named `p_adjust` is
numerically **identical to the raw GSEA p-value** across all 30 pairings. The
real BH correction is stored separately in three signature families (12, 6 and 12
tests). The external-validation conclusion is unaffected — all 10 expected
pairings survive under the stored signature-family FDR *and* under BH recomputed
across all 30 — but Methods must not call that column adjusted. Recorded as
`MT-04`, high priority.

### PH-016 — leading-edge `fdr_note` provenance string is wrong in two panels
The `fdr_note` column in the panel g and panel i source data prints *"smallest BH
FDR = 1.00"* while asserting cross-program scope; only panel h carries the
correct cross-panel value of 0.53. The renderer computes the minimum on the
panel-local table. The draft uses 0.53, which is the correct cross-panel value.

## PH-009 disposition

**PH-009_REWORD_RECOMMENDED**, with the technical meaning confirmed.

"Specificity comparisons" denotes 20 **off-target** contrast–signature pairings —
internal contrasts tested against signatures they are not expected to match. That
is a technical reference-signature check, so the phrase is literally accurate.

It is nonetheless misleading in effect: **18 of the 20 off-target pairings also
reach `p_adjust < 0.05`**, so the comparisons do not discriminate, and the label
invites a reader to infer that spatial specificity was demonstrated. A formal
expected-versus-off-target test does not exist (`NOT_FOUND`).

The drafted Results avoid the word "specificity" entirely in the external
validation paragraph and claim only what the ten expected pairings support.

## Wording decisions

- **`spatially resolved`** is the only spatial descriptor (DEC-010, statistical
  contract §6). `restricted` / `selective` / `specific` are reserved for a formal
  test or an explicitly factual support count.
- **No cross-family comparison.** Protein-level and program-level results sit in
  separate BH families that were never placed on a common scale; neither is
  claimed to be the stronger.
- **CAMERA is a sensitivity analysis**, never validation, replication or
  confirmation. The only external validation is the CON-only anatomical
  signature test (F2g, ED2b).
- **Leading-edge proteins are descriptive** and are not implied to be
  individually FDR-supported.
- **Microglia** wording is always `microglia-enriched ROI` or
  `microglia-enriched local microenvironment`; never a cell-intrinsic property.
- **Soma** is region-level; `CA2-SP`-style labels are source-localisation
  context, not a tested laminar hierarchy.
- **Nulls** are stated as "did not survive correction" / "no detectable … at
  three animals per group", never as absence of an effect.

## Unresolved factual questions

_(populated during drafting)_

## Resume point

**LAST COMPLETED:** Phase 1 complete. Skeleton, Results §2 and §3 drafted,
red-teamed and revised; all four provenance/review tables generated.

**NEXT PHASE (2):** Methods, plus Figure 1 / behavioural Results. Work the 13
rows of `manuscript/methods_todo.csv` in priority order; MT-04 (the `p_adjust`
column is not adjusted) and MT-13 (bilateral contrast registry and the DG mirror)
are the two that could otherwise produce a wrong Methods sentence.

**PHASE 3:** Introduction, Discussion, Abstract, title selection.

**DO NOT REPEAT:** the publication-hardening audit; atlas selection; any
statistical recomputation; the evidence extraction (108 checks, 0 wrong, all
recorded in the provenance tables).

**CARRY FORWARD INTO METHODS:** the six documentation findings PH-012…PH-016 and
PH-009 are recorded but unrepaired. Methods must be written from the verified
source data, not from the stale artefacts, in every one of those places.

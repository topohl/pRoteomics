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

## Phase 2 (HEAD cbe8463 → this commit)

### Figure 1 is BLOCKED — the analysis is not in this repository

> **SUPERSEDED BY PHASE 2C.** This section is the phase-2 record and is retained
> as written. Figure 1 is no longer blocked: the behavioural analysis was
> reconstructed in the upstream repository and imported here as a frozen,
> hash-verified evidence bundle. Read the Phase 2C section below for the current
> state. One conclusion recorded here was not merely unresolved but **wrong** —
> the predictor is not the GAMM-derived AUC (see `figure1_bridge_conflicts.csv`
> row FC-01).

This is the headline result of phase 2 and it is a finding, not a delay.
Item-by-item evidence: `manuscript/figure1_authoritative_source_inventory.csv`
(12 rows: 1 PRESENT, 2 PRESENT-as-external-input, 4 NOT IN REPOSITORY, 3 ABSENT,
1 EMPTY, 1 NOT RECOVERABLE).

**Present as finished external inputs:** `E9_Behavior_Data.xlsx` (34 sheets; the
`zScore` sheet has 117 animals and a **precomputed `CombZ` column**) and
`auc_individual_animals_*.csv` (322 animals, already carrying the RES/SUS label,
`prediction_type = subject_specific_gamm_observed_grid`).

**Absent:** the CombZ construction (read and renamed at
`analysis/preprocessing/build_module_score_metadata.R:399-427`; no z-scoring,
orientation, aggregation, sex or batch handling anywhere); the RES/SUS
cut-point; the GAMM specification; **any prediction or cross-validation analysis
at all** — no `glmnet`, `caret`, `pROC`, `randomForest`, `cv.glmnet`,
`trainControl`, `createFolds`, with every "leave-one-animal-out" match being a
*proteomic* stability analysis; any HMM or state model; any sex-stratified or
interaction model; and any record of ages or windows.

**Consequence:** "predicts" cannot be licensed, no AUC or permutation null can be
quoted, no sex claim of any strength is available, and temporal ordering cannot
be asserted. The §1 placeholder states this rather than substituting invented
content. The full checklist was still applied —
`manuscript/figure1_red_team_review.csv`, 8 questions, all unanswerable, which is
itself the auditable record.

### Methods drafted for Figures 2–3

Sixteen sections, every parameter traced to
`manuscript/methods_statement_provenance.csv` (18 rows: 15 VERIFIED, 1 PARTIAL,
2 NOT_VERIFIABLE). Three `[METHOD DETAIL UNRESOLVED]` markers rather than
conventional filler. Methods red-team:
`manuscript/methods_red_team_review.csv` — 12 checks, 7 clean, 2 defects found
and fixed, 1 avoided, 1 marked unresolved, 1 outstanding.

### WGCNA "0 of 45" resolved — it was neuropil-only

Recomputed from `results/tables/06_modules_WGCNA/group_effects/*/module_group_effects.csv`
(`manuscript/wgcna_45_contract.csv`). Each compartment carries its own BH
families.

| Family | Scope | n | supported | smallest FDR |
|---|---|---|---|---|
| primary + secondary | **neuropil only**, 15 modules × 3 group contrasts | **45** | 0 | **0.245** |
| primary + secondary | all three compartments, 35 modules × 3 contrasts | **105** | 0 | **0.156** |
| interaction omnibus | all three compartments, 1 test per module | **35** | 0 | **0.274** |

The reviewer who argued the full set is "0 of 105, smallest 0.16" was right. The
Phase-1 Results sentence has been given the minimum correction to state both
scopes; Methods states all three families.

### PH-009 and MT-04 resolved

See `docs/PUBLICATION_HARDENING_PROGRESS.md`. In short: "specificity comparisons"
is literally accurate but misleading (18 of 20 off-target pairings also clear the
threshold; no discrimination test exists), and `p_adjust` equals the raw *P*
**because each pairing is a single-signature GSEA and BH over a family of one is
a no-op** — a scope artefact, not a coding error. My earlier note was right about
the effect and wrong about the cause.

### PH-012 bypassed

Methods states the verified frozen membership — 16 mitochondrial GO terms,
glycolysis excluded, PDH and TCA retained — not the stale 20-term/glycolytic text
still carried by two reader-facing artefacts.

## Resume point

**LAST COMPLETED:** Phase 1 complete. Skeleton, Results §2 and §3 drafted,
red-teamed and revised; all four provenance/review tables generated.

**PHASE 2 COMPLETE** apart from Figure 1, which is blocked on an upstream
repository. Methods for Figures 2–3 drafted; WGCNA 0/45, PH-009, MT-04 and
PH-012 all resolved.

**BLOCKING QUESTION FOR THE AUTHOR:** where is the behavioural analysis
repository? Figure 1 Results and the behavioural half of Methods cannot be
written without it. Everything needed is listed in
`manuscript/figure1_authoritative_source_inventory.csv`,
`outcome_score_contract.csv` and `behavior_prediction_contract.csv` — the
critical items are the CombZ construction, the RES/SUS cut-point, the GAMM
specification, the experimental timeline, and whether any out-of-sample
prediction analysis exists at all. If no prediction analysis exists anywhere, the
manuscript must not use "predicts".

**PHASE 3:** Introduction, Discussion, Abstract, title selection. Also
outstanding: MT-01 (software versions, seeds, ontology release) is the last
mechanical gap in the Figures 2–3 Methods.

**DO NOT REPEAT:** the publication-hardening audit; atlas selection; any
statistical recomputation; the evidence extraction (108 checks, 0 wrong, all
recorded in the provenance tables).

**CARRY FORWARD INTO METHODS:** the six documentation findings PH-012…PH-016 and
PH-009 are recorded but unrepaired. Methods must be written from the verified
source data, not from the stale artefacts, in every one of those places.

## Phase 2C (HEAD e2fa49b → this commit)

### Figure 1 is unblocked — by import, not by analysis here

The behavioural analysis was reconstructed in `topohl/MMMSociability` and frozen
for manuscript use. Phase 2C consumes that freeze. **No behavioural statistic is
computed in this repository and no behavioural analysis code was copied into
it.** The evidence interface is five CSV contracts, mirrored byte-identically at
`manuscript/figure1_bridge_mmmsociability/`.

- **analysis commit:** `4b0f90f` · **bundle commit:** `53bc7e9` ·
  **verified source HEAD:** `a53d73f` · source tests 30/30, tested state HEAD
- **Integrity:** all 5 contract files byte-identical to the frozen export and
  provably unchanged from `53bc7e9` through `a53d73f`; all 8 upstream source
  tables re-hashed at import, 8/8 MD5 and 8/8 byte sizes matching.
- `.gitattributes` pins the bridge to `-text`. Without it the repo-wide
  `* text=auto` would LF-normalise these CRLF files and the recorded SHA-256
  would not reproduce on checkout.

### Drafted

**Results §1** — 819 words. Heading: *Early spontaneous home-cage activity
predicts later composite stress outcome*. Six paragraphs: experimental logic and
temporal ordering; outcome definition with the by-construction circularity stated;
the movement association led by ρ = −0.39; prospective prediction; sex; limitations
and the bridge into §2.

**Behavioural Methods** — 988 words, ten subsections, led by a provenance
subsection. `M-19` … `M-28`.

**Results §4** — moved from pending placeholder to **CLOSED**. No
behaviour–proteomics claim is supported in either repository (upstream BH-006;
the edge–behaviour null already reported in §3).

### Provenance

- `results_claim_provenance.csv` — 13 new rows `C1-1` … `C1-13`
- `results_statement_provenance.csv` — 13 new rows `S1-01` … `S1-13`
- `methods_statement_provenance.csv` — 10 new rows; `M-16`/`M-17` marked
  `SUPERSEDED_BY_FROZEN_BRIDGE`, `M-18` `SUPERSEDED_AND_CORRECTED`
- `figure1_bridge_provenance.csv`, `figure1_bridge_import_manifest.csv`,
  `figure1_bridge_conflicts.csv` — new
- `figure1_red_team_review.csv` — 8 UNANSWERABLE verdicts → 8 answered, phase-2B
  verdicts retained in a `phase2b_verdict` column; 2 questions added

Claim IDs deliberately do **not** reuse the bundle's `F1-xx` namespace, which
already means something else (`FC-03`). Every `C1-x` row names the `F1-xx` or
`FM-xx` row it resolves to.

### Twelve disagreements found and recorded

`manuscript/figure1_bridge_conflicts.csv`. The four that changed the prose:

| id | finding | resolved to |
|---|---|---|
| FC-01 | phase-2B called the predictor a GAMM-derived AUC | raw `Movement_mean` per bundle `FM-02` |
| FC-02 | bundle says repeated-CV `seed 123`; upstream code passes `seed = 521` | **521** — the bundle field is a transcription defect |
| FC-06 | "movement-only" is an upstream model name meaning `Sex + Group + Movement_mean` | never used; "movement-mean model" instead |
| FC-09 | `[0.116, 0.179]` is a percentile range, not a confidence interval | labelled as such in both Results and Methods |

`FC-02` and `FC-04` are defects in the frozen bundle and are flagged for upstream
repair. Neither changes a reported value.

### Guard added

`tests/testthat/test-figure1-behaviour-bridge.R`, 119 assertions. Before this,
**no test in the repository read anything under `manuscript/`** — the drafted
prose was mechanically unguarded. It pins bundle byte-integrity against recorded
SHA-256, the provenance chain, the absence of copied upstream code, the
`.gitattributes` rule, that §1 is written and §4 stays closed, that every
behavioural number in the draft is present in the frozen contract, the CombZ sign
contract, prohibited wording with the repository's standing denial exemption, and
that superseded phase-2B assertions are marked rather than left true.

### Not done, deliberately

Introduction, Discussion, Abstract. No Figure 1 graphic: none exists in this
repository by design, and Results §1 references panels `1a`–`1e` using the panel
assignments carried in the imported claim contract. No behavioural analysis was
re-run. No proteomics numerical output was touched.

## Phase 3 (HEAD 0a53d3d → this commit)

### Figure 1 frozen

Inspected at 183 mm publication width, all four panels rendered at their placed
size so native pt equals printed pt. No clipping, no overlapping annotation, no
unreadable text; the type scale is 50 elements at 6.50 pt and 12 at 7.00 pt,
which is the distribution Figures 2 and 3 use. Panel b reads as a definition
rather than a comparison: it carries the title *Later composite outcome and
phenotype definition*, the subtitle *the threshold defines these groups; it is
not a test of them*, no brackets and no stars.

**FIGURE_1_FROZEN = TRUE.** The renderer and its outputs are unchanged by this
phase.

One minor observation recorded and deliberately not acted on: in panel d the
group legend key for RES and the sex legend key for Female both render as grey
circles. The adjacent labels disambiguate them and the original design had the
same structure, so this did not meet the bar for reopening a frozen figure.

### Drafted

**Introduction** — 754 words, five paragraphs: adolescent plasticity and outcome
heterogeneity; the limits of endpoint behavioural sampling and the case for
continuous home-cage measurement; anatomical heterogeneity and the case for
spatial proteomics; the gap, stated as two separate questions rather than one
causal chain; and a study overview ending on aim rather than result. Six `[REF]`
markers, none invented.

**Discussion** — 2,509 words across ten subsections: principal findings; early
behaviour; sex; what the spatial measurement establishes; sparse protein-level
differences alongside coordinated program-level ones; spatial context; CA2-SLM;
co-abundance modules; candidate proteins; limitations; conclusion.

Three interpretive positions are worth recording because they were the hard
calls. The sparse-versus-coordinated section argues that protein-level testing
and ranked enrichment ask different questions and explicitly refuses the reading
that enrichment is more sensitive and therefore more true. The spatial section
states the claim as resolution rather than specificity, because establishing that
a program differs in one context and not another needs a heterogeneity test and
the omnibus tests did not survive correction. The CA2-SLM section presents the
quality-control audit as a strength rather than a caveat, on the grounds that an
uncorrected version of that result would have placed a striking anatomical claim
on the least reliable unit in the dataset.

### Provenance

- `manuscript/discussion_statement_provenance.csv` — new, 28 rows: 6 RESULT,
  7 INTERPRETATION, 9 LIMITATION, 1 LITERATURE_CONTEXT and 4 combinations. Every
  row resolves to an existing `C1-x`, `C2-x` or `C3-x` claim; no new scientific
  claim was created.
- `manuscript/citation_needs.csv` — new, 9 rows. Six required Introduction
  citations, two optional Discussion ones, and one already resolved by the
  reference the Results and Methods already carry.

### Semantic red-team

Twenty-five prohibited or risky terms scanned across the Introduction and
Discussion; 32 occurrences found and adjudicated individually. **Zero required
change.** Every occurrence is a denial (*not a molecular hotspot*, *rather than
specificity*, *not independent replication*, *constrains temporal order, not
mechanism*), a technical term (*data-independent-acquisition*, *cross-validation*,
*global structure of the dataset* in the PCA sense), or qualified in the same
clause (*sampled to favour microglia, not a purified population*). The single
bare use of "microglia" occurs inside its own disclaimer.

### Not done, deliberately

The Abstract, which should be written last. No analysis was rerun, no Results or
Methods value changed, Figures 2 and 3 were untouched, and the behaviour–proteomics
integration section remains closed on BH-006.

## Phase 4 (HEAD 59b496b → this commit)

### Scope

Close the submission-readiness gaps: the Methods TODO register, the behavioural
Extended Data, the figure export stage, the title, the Abstract, and an audit of
what is left. No analysis was rerun, no Results or Methods value was changed, and
no behavioural statistic was computed here.

### The finding that matters most

**Results §2 and §3 cite a figure generation that was never promoted.** The draft
was written against `final_truth_v9`, in which Figure 2 has panels a–h and
Figure 3 has panels a–i. What `figures/figure_02.R` and `figures/figure_03.R`
render, and what the export stage ships, is `manuscript_figures_v2` — Figure 2
with panels a–f and a Figure 3 that is a different figure, WGCNA-centred rather
than program-exemplar-centred. The v9 contract's own `status` field reads
`candidate_only_not_promoted`.

This is not a numbering drift and it is not repaired here. The v2 Figure 2 has no
CON-only spatial fingerprint panel and no bilateral reproducibility panel, and
the v2 Figure 3 has none of the three program-exemplar GSEA curves or protein
heatmaps that Results §3 points at. Renumbering the references would silently
substitute different evidence for the evidence the text describes. Both figure
sets exist and are rendered; what does not exist is a decision about which is the
manuscript's. Recorded as SR-01 to SR-03 in `manuscript/submission_readiness.csv`
and left for the author.

### Closed

- **Methods TODOs.** MT-01 software, seeds and ontology release, as a new
  *Software and reproducibility* subsection; MT-05 CA2-SLM prespecified
  thresholds; MT-12 network section retitled with its edge definition stated;
  MT-03 external reference cited. Every value came from a per-stage
  `sessionInfo.txt`, a committed config, or installed package metadata. The one
  value that could not be recovered — the acquisition and search settings — is
  marked in place and explained, not filled from convention.
- **Figure export.** `copy_export_targets` called `file.copy` without creating
  parent directories, so all 21 curated Figure 1–3 copies returned `FALSE` after
  ~5,200 extended-data files had been overwritten and before the manifests were
  rewritten — a silent partial export. Fixed and guarded; the stage now writes
  5,235 manifest rows including all 23 curated copies.
- **Behavioural Extended Data**, as Extended Data Figure 9: the complete a-priori
  model registry, the formal feature-by-sex interaction tests, and the
  sex-stratified correlations as descriptive context. Two tables imported
  byte-identically from the same upstream freeze the Figure 1 bridge rests on.
  The renderer computes nothing and refuses to draw panel c unless every
  interaction is still classified `FORMAL_INTERACTION_NOT_SUPPORTED`.
- **A dangling `Fig. 1e`.** Panel e was merged into panel d during the Figure 1
  rebuild, but Results §1 and two provenance rows still pointed at it. They now
  point at Extended Data Fig. 9b,c, which is where that evidence is.
- **Title**, selected from four candidates with each rejection recorded, and the
  **Abstract**, 243 words, with all eighteen quantities resolved to the Results
  statement they are quoted from.

### Not done, deliberately

SR-01 to SR-03, above. The acquisition settings and the cage-change schedule,
neither of which is in this repository. The PRIDE accession and the front-matter
sections. Discussion compression, which depends on the target journal.

## Phase 5 (HEAD 146b036 → this commit) — figure-generation adjudication

### Decision

**PROMOTE_V9 = NO, at this HEAD.** Not because v9 is wrong — on the evidence it is
the right generation — but because three defects would ship with it, one of which
is a rendering fault that makes Figure 3 assert a false anatomical claim.

### What the adjudication established

v9 **strictly dominates** v2 on evidence. Of 13 Results §2/§3 claims: 7
`SUPPORTED_BY_BOTH`, 6 `SUPPORTED_BY_V9`, **0 by v2 alone**, 0 by neither. v2
cannot support the bilateral-concordance claim (no such panel exists), the
CA2-SLM QC weakening (its DAP source has no robustness column at all), the
seven-program atlas, or the three exemplar curves.

The most instructive near-miss: manuscript Fig. 3b claims *"FDR-supported
constituent terms occur in 132 of the 378 theme × unit × contrast cells."* The v2
panel 3b source also has exactly 378 rows — and **zero** FDR-supported cells
(`tier_specific_fdr` minimum 0.2493). The row count is a coincidence and the
evidence is the opposite. Repointing the reference would have looked clean and
been badly wrong. This is why the Phase 4 instruction not to repoint blindly was
correct.

v9 is also current where it was most suspect: registry `manuscript_go_themes_v3`,
seven exact programs, 253 constituent GO terms, mitochondrial theme exactly 16
terms with glycolysis excluded and PDH/TCA retained; DAP arithmetic 37/28/6/9/15
exact; the three exemplars exact; leading edge 63 values with no invented
confidence interval and an explicit statement that none is individually
FDR-supported.

### Why it is still blocked

**PB-01 is the serious one.** `NF_RGT <- 7.6` is defined at
`R/nature_final_v7_figure3_panels.R:48` as *"mm reserved at the right for the
atlas legend, in a AND b"*. It is applied in all three DAP-track renderers (v7,
v8, v9) and **in no atlas renderer** — `f9_gsea_atlas` sets `legend.position =
"right"` with no matching `plot.margin`. So the "a AND b" contract was
implemented on one side only, in three successive generations.

Measured in the assembled SVG: panel a pitch 18.720 pt, panel b 17.540 pt. The
headline `28` sits at x = 251.89, which is panel b's **CA3-SO** column (254.46),
not CA2-SLM (236.92). Read as drawn, Figure 3 states that CA3 stratum oriens has
28 FDR-supported and 6 robustness-qualified proteins. Its true values are 0 and
0, and the CA2-SLM result the panel exists to carry is erased. A re-render does
not fix it: the defect is deterministic code.

PB-02 is structural — v9 Figure 2 panels b and c read their inputs from the *v2*
export namespace, so promoting v9 and retiring v2 removes their only producer.
PB-03 is a hard-coded stale legend constant that re-rendering cannot correct.

### Also found

The **currently shipped** Figure 3 is scientifically superseded (SR-21).
`results/manuscript/figure_3/panels/figure_03a.svg` draws CA2-SLM as a bar 25
units long against a maximum of 2 elsewhere, with no QC row — precisely what the
standing requirement forbids. Neither generation is shippable today.

`docs/REPOSITORY_ARCHITECTURE.md` already records `final_truth_v9` as the CURRENT
layer while the v9 contract's own `status` field still reads
`candidate_only_not_promoted`. That self-contradiction is the root of SR-01–03.

### Not done, deliberately

No promotion, no contract switch, no re-render, no export. Per the standing
instruction to stop before promotion when a critical panel fails. No analysis was
rerun and no Results value changed.

## Phase 5B (HEAD ebe3194 → this commit) — repair, promote, ship

**PROMOTE_V9 = YES.** Figures 2 and 3 are now `final_truth_v9`: Figure 2 a–h,
Figure 3 a–i, the structure Results §2 and §3 were written against. No Results
prose was rewritten, because the promoted panels *are* the evidence the prose
describes.

### The three repairs

**PB-01 — column misregistration.** Root cause was inter-panel layout geometry,
not data: both panels compute `xpos` from the same `sg_blocks()` order, but the
atlas let a right-hand ggplot legend take 18.07 mm of layout width while the DAP
track reserved the shared `NF_RGT` constant of 7.60 mm. Pitches were 17.540 and
18.720 pt, so the headline `28` sat on **CA3 SO** and the figure asserted a value
of 28 for a unit whose true value is 0. The atlas now declares the coupling
(`shares_column_geometry_with`) and moves its colour bar below the axis, pinning
the right gutter to `NF_RGT`. Measured after: both panels left inset 137.60 pt,
right inset 21.54 pt, pitch 18.720 pt — **18 of 18 columns misregistered before,
0 of 18 after, maximum offset 0.00 pt.**

**PB-02 — provenance.** Figure 2 b and c read from the *v2* export namespace, so
retiring v2 would have removed the only producer of their inputs. They now read
the canonical acquisition workbook and the canonical joint-QC PCA scores. Plotted
values identical before and after (323×5 and 323×4, `all.equal` TRUE).

**PB-03 — hard-coded constant.** The eps-floor disclosure said "90 of the 851"
and could not be corrected by re-running, because it was a literal. It is now
derived by `f9_eps_floor_disclosure()`, which also asserts that the smallest
positive raw p in the canonical output equals the declared floor. The legend now
reads **94 of 953**.

### Promotion

The manuscript assembler gained an opt-in `layout_mode: absolute`, because the
promoted figures cannot be expressed on an equal-cell grid without distorting
them. All 17 promoted panels are byte-identical copies of the v9 renders, so the
repaired geometry is preserved exactly.

### What shipped

Export ships Figure 1 a–d unchanged, Figure 2 a–h and Figure 3 a–i.
**SR-21 is closed on the artefact**: `figure_03a.svg` is byte-identical to the
post-QC v9 DAP track, draws both rows (37 total / 15 total), contains 28 and 6,
and contains no "hotspot". The pre-QC bar panel is no longer reachable from the
manuscript contract.

### Integrity

24 of 24 validity checks PASS. Changed-output manifest: **0
UNEXPECTED_SCIENTIFIC_CHANGE**. Figure 1 byte-identical. No analysis rerun, no
frozen source data changed, no Results value changed.

### Still open, deliberately

The four external-record items only: acquisition and search settings, the
cage-change schedule, the PRIDE accession and the front matter. Discussion
compression remains deferred as SR-18.

## Phase 6A (HEAD d47114b → this commit) — canonicalise and freeze Extended Data

The last publication-generation ambiguity is gone. Every numbered figure now has
exactly one canonical identity, recorded in
`manuscript/canonical_publication_registry.csv`: **9 canonical** and **3
withheld**, with the withheld identities reserved so nothing can quietly take
them.

### Promoted (6)

`extended_data_01` (bilateral), `extended_data_03` (CA2-SLM),
`extended_data_06` (atlas and exemplar curves), `extended_data_08` (networks),
plus the two behavioural figures below. Each was audited by reading its source
tables, not by trusting the generation name.

### Not promoted (3), and why

- **ED2** — panel c is presented on the figure, in its legend, in the contract
  narrative and in the supplementary table as the *complete* canonical GO
  inventory, and Figure 2h sends the reader to ED2 for "the complete evidence".
  It draws **14 rows across 7 contrasts** against **2,826** FDR-supported
  positive GO-BP terms across **11** contrasts — 0.5%. Verified independently.
  The drawn values are fine; the completeness claim is not.
- **ED4** (the WGCNA figure) — its phenotype panel annotates *"0 of 45 module ×
  contrast cells"* with **no compartment scope**, immediately beside an
  all-compartment *"0 of 35"*, in a figure whose panel a is explicitly
  three-compartment. 45 is neuropil-only; the all-compartment figure is 105.
  It also labels all 15 modules `peak <unit>` while 8 of 15 are classified
  `has_spatial_identity = FALSE` by the canonical atlas.
- **ED7** — scientifically current but cited nowhere. There is nothing to
  promote it into.

Nothing was substituted for any of them.

### The behavioural Extended Data was split

One figure was carrying two arguments. `extended_data_05` is now how the early
window was measured (111 animals, 50 contributing every expected slot, 61 missing
leading slots only, mean 98.6%, minimum 94.4%) plus the complete a-priori model
ladder. `extended_data_09` is the two secondary features and everything the study
can say about sex.

**One thing could not be built as specified.** Panel 9a was to show "Movement
RMSSD vs CombZ" and "Entropy ACF1 vs CombZ". The frozen bundle exports per-animal
values for `Movement_mean` **only** — `figure1c_movement_combz_source.csv` has no
RMSSD or entropy column and no other bridge table is per-animal. A scatter would
have required inventing the points. The panel draws effect sizes instead (ρ with
its bootstrap interval and BH q), which is exactly the quantity the manuscript
claims, and a test asserts that no per-animal source exists.

### Pre-restructure freeze

`manuscript/prerestructure_freeze_manifest.csv` — **230 objects, all present,
all hashed**: configuration contracts, manuscript provenance, the frozen upstream
bridge, every canonical panel and its source data, protected scientific state and
the guard tests. This is the equivalence oracle for the migration. See
`docs/PRERESTRUCTURE_FREEZE.md` for how to use it and for the two hash changes
that would be legitimate rather than defects.

### Not done, deliberately

No directory moved, no script renamed, no manuscript layer extracted, no import
path rewritten. This phase froze identity; it did not restructure.

## Phase 6A.5 (HEAD c6fc77f → this commit) — resolve Extended Data 2

**ED2 = CANONICAL**, with panels a and b. The last withheld identity the
manuscript depended on is gone, and the registry is now 10 canonical / 2
withheld.

### The proposal was wrong, and the check caught it

The obvious resolution — promote ED2 with the two sound panels, drop the
defective third, change nothing else — survived the structural attack and failed
the scientific one. Panel b, *the panel the manuscript actually cites*, carried
the same class of defect that withheld panel c: it marked significance from the
column named `p_adjust` and its legend and ST2 called that Benjamini-Hochberg
adjusted. In this analysis `p_adjust` is byte-identical to
`single_set_p_adjust` — BH over a family of one is a no-op — so a false
multiple-testing label sat on the study's **only externally anchored claim**.
Promoting it unfixed would have canonicalised exactly what the phase existed to
remove.

Repaired before promotion: the panel now uses `signature_FDR`, and its released
source data carries both statistics under their true names. **No mark moved** —
the two agree on all 30 rows, 28 significant either way.

### What else had to move with it

Dropping panel c does not remove a misrepresentation if its twin keeps shipping.
Corrected in the same pass:

- **ST3** was generated from the same 14-row sidecar under the title *"Every
  canonical GO term retained for every anatomical contrast"*. It now states its
  selection and names where the complete result actually is.
- **Figure 2h's pointer to ED2** existed in **four** places, one of them inside
  its own published source data (`legend_text_moved`). All four now name the
  released `control_anatomical_go_bp_gsea` supplementary table — 40,680 rows
  across all 11 contrasts.
- **ST2** relabelled the uncorrected p as BH-adjusted with a matching data
  dictionary. Both statistics are now named honestly.

### A detector that could not see the thing it was checking

Phase 5B closed SR-02 on "49 citation instances, 0 unresolved". That audit
matched only `Fig. Nx` and therefore examined **zero Extended Data citations**
while returning green. It is now `tools/audit_manuscript_references.R`, which
resolves each citation against the contract, the registry and the exported
artefact, and also reports the reverse defect. Current: **37 instances, 24 main
and 13 Extended Data, 0 unresolved, 0 canonical-but-uncited.**

That reverse check immediately found one: `extended_data_08` was canonical,
fully rendered and cited nowhere in the prose. The network paragraph states
exactly its evidence, so it now carries the citation.

### Still withheld, deliberately

ED4 and ED7, untouched. Their identities stay reserved.

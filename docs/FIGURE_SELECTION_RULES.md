# How the displayed proteins and pathways were selected

Every protein and every pathway drawn in Figures 2 and 3 is a **selection** from
a larger set. This document states, for each one, the rule that chose it, the
set it was chosen from, and whether that rule looked at the result the panel
displays.

It exists so that "why is this protein in the figure?" has a written answer that
does not depend on who is asking, and so that no selection in the package is
defended by memory.

> **The honest summary in one paragraph.** Everything in Figure 2 was selected
> by a rule fixed before any stress contrast was computed. Figure 3b is selected
> by the strictest rule in the package, which explicitly forbids selecting on
> effect size or phenotype direction. The three exemplar programs in Figure 3
> d–i are an **editorial choice made after the results were known** — one per
> compartment, to show distinct biology in each — and the proteins in 3 g/h/i
> are a top-7 display cut of whatever those programs contain. That last
> selection is the only unruled one, it is disclosed as such everywhere it
> appears, and the complete set it was drawn from is released.

---

## The three selection regimes

| Regime | Meaning | Panels |
|---|---|---|
| **A. Prespecified** | The rule and the set were fixed before any stress contrast existed. The displayed result cannot have influenced the choice. | 2d, 2e, 2g, 3b |
| **B. Result-ranked inside a prespecified set** | The *set* was fixed in advance; within it, the item shown is the extreme of the displayed quantity. | 2h |
| **C. Editorial** | Chosen by judgement after the results were known. | 3 d/e/f, and by inheritance 3 g/h/i |

A regime-C selection is legitimate for an **illustrative** panel whose complete
evidence is released elsewhere. It is not legitimate as the evidence itself.
That is exactly how Figure 3 is built: the atlas (3b) carries the evidence, the
exemplars illustrate it.

---

## Figure 2 — all selection is phenotype-blind

### 2d · spatial fingerprint — 19 proteins

**Rule.** For each of the **11 prespecified CON-only anatomical contrasts**, take
the **top 2 genes by BH-adjusted p** among those enriched on the target side
(logFC > 0). Ties break on −|logFC|, then on gene symbol, so the pick is
deterministic and reproducible. Rows are then ordered by each gene's baseline
CON peak spatial unit.

**Selected from.** The 11 contrasts, fitted on CON animals only. 22 picks yield
19 unique genes, because OXR1 and WFS1 each win two contrasts.

**Did the rule see the stress result?** No. The contrast set was fixed in
advance, the fit uses control animals only, and no stress group enters either
the selection or the row order.

> ADCY9 ANKRD63 ATP1A1 ATP8A1 CAMK1D CRACDL ELAVL2 GPC1 INA KCND2 NPTX1 OXR1
> PACSIN2 PDE1B PEX5L PITPNM2 SRGAP2 SYT12 WFS1

Implemented in `Exp9_manuscript/figures/spatial_v6_fingerprint_selection.R`
(`TOP_N <- 2L`).

### 2e · compartment markers — 10 proteins

**Rule.** A written list of canonical compartment markers, 3–4 per compartment,
each of which must then pass **primary and strict detection eligibility** in its
intended compartment. If any requested marker fails, the build stops rather than
silently substituting another.

**Selected from.** Canonical marker biology, gated by the detection-eligibility
audit.

**Did the rule see the result?** No. The selection provenance explicitly records
`marker_selection_used_observed_cross_compartment_direction = FALSE` and
`..._effect_magnitude = FALSE`, and treats the expected direction as
`post_selection_validation_only`. These markers were chosen to be
*recognisable*, not to perform well.

> nuclear / soma: Npm1, Ptbp2, Anp32a, Hdac5
> synaptic: Camk2a, Snap25, Syp
> microglia: P2ry12, C1qa, Ctss

Implemented in `R/qc/control_compartment_abundance_rendering.R`
(`ca_compact_recognizable_marker_config_v2()`).

### 2g · external validation — 10 pairings

**Rule.** Show every pairing for which an **expected anatomical
correspondence** was declared: CA1 to CA1, CA2/3 to CA2/3, a target stratum to
its own signature. Expectation is a property of anatomy and was written before
the results.

**Selected from.** All 30 tested internal-contrast × external-signature
pairings. All 30 are released in `ST2_external_signature_validation.csv`, so
the 20 off-target comparisons can be inspected rather than taken on trust.

**Did the rule see the result?** No.

Implemented in `R/spatial/control_spatial_identity_utils.R` (`expected_match`).

### 2h · internal anatomical GO programs — 7 terms

**Rule.** **One GO term per anatomical contrast: the maximum |NES|**, ties
broken on term description.

**Selected from.** The canonical GO BP GSEA of each anatomical contrast. Seven
contrasts, one term each.

**Did the rule see the result?** **Yes** — the term shown is the maximum of the
quantity shown. This is regime B: the contrast set was prespecified, but the
term within it was picked for being the largest. The panel is therefore
characterisation, not validation, and the complete term inventory is released
separately.

---

## Figure 3 — one strict rule, one editorial choice

### 3b · the seven-program atlas

**Rule.** Rows are the themes flagged **claim-eligible** in a version-controlled
registry, drawn in **registry display order**. Admission to the registry
requires all five of:

1. biological coherence,
2. recurrence across multiple canonical spatial contexts,
3. distinctness — not an ontology ancestor ladder, semantic duplicate or renamed
   subset,
4. non-QC status,
5. material representation — omitting it would underrepresent a recurrent
   concept.

Membership of each theme is generated **phenotype-independently** from explicit
GO anchors using only the `is_a` and `part_of` relationships.

**Selection must never use** SUS/RES direction, NES magnitude, how publishable an
example looks, or the recurrent leading-edge proteins. Rows are **never sorted**
by NES, FDR, direction or number of supported cells.

**Selected from.** 12,598 claim-eligible term-occurrences out of **203,073**
tested GO term × spatial unit × contrast combinations — 4,199 of them for the
SUS − RES contrast the panel draws, of which **348** are FDR-supported. A cell
carries a support dot when at least one of its constituent terms is
FDR-supported; the dot is not a measure of how broad that support is.

(12,598 counts term *occurrences*. The stored theme table has 12,778 rows here
because 180 terms legitimately sit in two themes and appear twice; the
inventory is keyed on the occurrence, so its count is the denominator to
quote.)

**Did the rule see the result?** No. This is the most defensible selection in
the paper, and the one to lead with.

Full contract: [ATLAS_PROGRAM_SELECTION_AND_NAMING_RULES.md](ATLAS_PROGRAM_SELECTION_AND_NAMING_RULES.md).

### 3 d/e/f · the three exemplar programs — EDITORIAL

This is the selection that previously had no written rule. It does now.

**The criterion, stated:**

> **One FDR-supported program per compartment, chosen to represent distinct
> biology in each of the three compartments.** Each exemplar must be (a)
> FDR-supported for SUS − RES in the spatial unit shown, (b) a member of a
> claim-eligible atlas theme, and (c) biologically distinct from the other two,
> so that the three together span the compartments rather than repeating one
> program. Within those constraints the choice is **editorial judgement, not an
> algorithmic maximum.**

| Program | Compartment · unit | GO term | Rank by \|NES\| in its own unit |
|---|---|---|---|
| synaptic signalling | neuropil · CA3-SR | GO:0099536 | **51 of 64** |
| mRNA processing | neuronal soma · CA2-SP | GO:0006397 | 6 of 24 |
| oxidative phosphorylation | microglia-enriched · CA1 | GO:0006119 | 3 of 22 |

Ranks are among FDR-supported claim-eligible terms in that spatial unit for
SUS − RES, deduplicated to one row per term occurrence.

**Did the rule see the result?** **Yes.** The three were chosen after the
phenotype-contrast results were known, and they are recorded as such in
`f3_program_example_selection.csv`
(`selected_after_seeing_phenotype_direction = TRUE`).

**Why this is defensible.** The three exemplars are *illustrations* of the atlas
in 3b, not independent evidence. The atlas is unselected and carries the claim;
the exemplars make one cell of it concrete and readable. Every place the
manuscript refers to them, it is bound by a semantic rule to say
*representative* or *illustrative* and is forbidden from saying *strongest*,
*top*, *dominant*, *most altered* or *major response*.

**Why they were not chosen as the strongest.** Choosing the maximum |NES| would
have made the figure less defensible, not more: at n = 3 per group the maximum
of a statistic is the value most inflated by sampling noise, selecting on the
maximum of a displayed quantity biases that quantity upward by construction, and
the atlas contract this project already committed to explicitly forbids
selecting on NES magnitude. It would also have collapsed the three-compartment
story, because RNA-processing terms carry larger NES than the synaptic exemplar
inside CA3-SR itself.

**The weak point, stated plainly.** The synaptic exemplar ranks 51st of 64 in its
own unit. Ranks of 6/24 and 3/22 are easy to describe as representative; 51/64
rests entirely on the "distinct biology per compartment" criterion above. If
that criterion is not accepted, that panel is the one to defend.

### 3 g/h/i · leading-edge proteins — 21 proteins

**Rule.** Take the **leading edge of the exemplar's own enrichment**, rank by
**|stored rank statistic|**, keep the **top 7**. The cut at 7 is a display
constraint chosen for legibility, not a statistical threshold.

**Selected from.** The complete leading edge of each exemplar term, for
SUS − RES: **255** genes (synaptic), **112** (mRNA processing), **48** (OXPHOS).

**What the released inventory does and does not let you check.**
`leading_edge_protein_inventory.csv` records **membership** — which proteins are
in each term's leading edge — and the **universe size** above. It does **not**
preserve the per-protein GSEA rank statistic or the original ranked-list order,
because the upstream theme table supplies `leading_edge_genes` as an
alphabetically sorted string, so neither the value nor its order reaches the
inventory. A reader can therefore confirm that the seven shown were eligible and
how many they were drawn from, but cannot reproduce *which* seven the rule
picked. Carrying the statistic would mean joining the GSEA ranked lists, which
is tracked separately as `RANKED_LIST_JOIN_PENDING`.

**Did the rule see the result?** **Yes, twice over** — both the parent program
and the ranking within it are outcome-dependent. This panel inherits the
regime-C status of the exemplar that produced it.

> synaptic: App Cnr1 Dbi Eif4ebp2 Ly6h Plppr4 Synpo
> mRNA processing: Cirbp Csdc2 Dcps Ddx23 Lsm3 Lsm8 Rbm8a
> OXPHOS: Cox5b Cox6b1 Iscu Ndufs8 Ndufv2 Ndufv3 Uqcrh

**The standing caveat.** **No protein in these panels is individually
FDR-supported.** The smallest BH FDR across all 63 displayed values is 0.53.
They are leading-edge contributors to an enrichment, and must never be called
validated, significant, key or driver proteins.

**An independent check that these proteins are not arbitrary.** They were
selected by the GSEA rank statistic. The atlas *recurrence* rule is a different
and unrelated criterion — a count of how many supported GO terms and how many
spatial contexts a protein appears in, which ignores effect size entirely.
Applying it after the fact:

| | Recurrent core | Intermediate | Single appearance |
|---|---|---|---|
| all leading-edge proteins in a claim-eligible theme (n = 1,992) | 43% | 43% | 14% |
| **the 21 shown in Figure 3 g/h/i** | **86%** | 14% | 0% |

18 of the 21 displayed proteins independently qualify as recurrent core, against
a 43% base rate. This does not make them individually significant and nothing
above is withdrawn — but it does mean the displayed proteins are the recurrent
ones, by a criterion that had no part in choosing them.

Implemented in `Exp9_manuscript/R/panels/final_truth_v9_panels.R`
(`f9_prot_values(..., top_n = 7L)`).

---

## The complete sets these were drawn from

Released by `analysis/integration/build_display_selection_inventories.R` into
`results/integration/build_display_selection_inventories/global/tables/`:

| Table | Rows | What it is |
|---|---|---|
| `pathway_enrichment_inventory.csv` | 203,073 | **Every** GO term tested in every spatial unit and contrast, with NES, raw p, BH FDR, theme assignment and claim-eligibility. The denominator for Figure 3b. |
| `pathway_enrichment_inventory_fdr_supported.csv` | 3,559 | The FDR-supported subset, for reading. |
| `leading_edge_protein_inventory.csv` | 282,296 | **Every** leading-edge protein of **every** FDR-supported term. The denominator for Figure 3 g/h/i. |
| `leading_edge_protein_recurrence.csv` | 12,227 | Each protein classified by the atlas recurrence rule (below). |
| `display_selection_disclosure.csv` | 7 | This document as a machine-readable table, one row per panel. |
| `inventory_data_dictionary.csv` | 72 | Column definitions and the four standing caveats. |

FDR-supported term-occurrences by contrast: **996** RES − CON, **1,361**
SUS − CON, **1,202** SUS − RES.

All six are in this repository's frozen export bundle, as publication identity
`supplementary_selection_inventories`:
`exports/publication_source_data/supplementary_selection_inventories/`, each
with a SHA-256 in the bundle manifest.

**Not yet imported by the manuscript.** The manuscript repository does not read
this directory at runtime — its trust boundary is a local, version-controlled
copy under `Exp9_manuscript/source_data/pRoteomics/`, refreshed by hand with
`tools/import_render_inputs.R`. That copy currently carries 54 rows across ten
publication identities and does **not** include these six. Releasing them to
readers therefore needs one further deliberate step: importing the identity into
the manuscript bundle. Until that happens the inventories are exported and
hash-manifested here, but are not part of what the manuscript ships.

### The figure

`analysis/integration/plot_display_selection_context.R` draws the same argument:
`results/integration/plot_display_selection_context/global/plots/display_selection_context.{svg,pdf}`.

- **a** — each exemplar inside the FDR-supported claim-eligible terms of its own
  spatial unit, so the 51/64, 6/24 and 3/22 positions are visible rather than
  asserted.
- **b** — the 21 displayed proteins against the recurrence rule they were not
  selected by.
- **c** — the funnel from 67,691 terms tested to the seven rows drawn.

This is a **diagnostic plot, not a manuscript panel.** It is deliberately built
in this repository rather than in the manuscript figure layer. Promoting it to
an Extended Data figure would touch the frozen figure contract, the legend
layer, the figure index and the promotion tests, and is a separate decision
that has not been taken.

### The recurrence classification

The atlas naming rules already define a protein-level rule that had never been
computed into an artefact:

> A protein counts as **recurrent core** only if it appears in **≥3 supported GO
> terms AND ≥3 spatial contexts**. Everything else is recorded separately as
> INTERMEDIATE or SINGLE_APPEARANCE.

`leading_edge_protein_recurrence.csv` applies it. For SUS − RES: **853**
recurrent-core, 855 intermediate, 284 single-appearance proteins.

This rule selects on **recurrence, not magnitude**, which is why it is a sounder
basis for any future protein-level display than a "largest fold change" cut
would be. The source document names three classes but fixes only the
recurrent-core boundary; this table's reading of the other two is stated in its
data dictionary rather than left implicit.

**It remains a count of enrichment memberships, not evidence about any
individual protein.** No protein in it is individually FDR-supported.

---

## What may and may not be said

| About | Say | Never say |
|---|---|---|
| the three exemplars | representative, illustrative, selected FDR-supported example | strongest, top, dominant, most altered, major response |
| leading-edge proteins | selected leading-edge proteins, contributors to the enrichment | validated, significant, key, driver proteins |
| atlas themes | theme-level summary, theme containing FDR-supported terms | FDR-significant theme, theme p-value |
| an unsupported result | did not survive correction, no detectable effect at this sample size | no difference, unchanged, equivalent |

Full list: the 20 semantic rules in
`Exp9_manuscript/figures/final_truth_v9_semantics.R`.

---

## Notes for anyone checking this document

- **The synaptic exemplar's rank has one defensible scope.** `PH-014` records
  that a hardcoded "51st of 60" in the semantics layer matches no computed
  scope. The value **51 of 64** given above is computed from the released
  inventory over FDR-supported claim-eligible terms within CA3-SR for SUS − RES,
  deduplicated to one row per term occurrence. Quote it with that scope
  attached, or not at all.
- **Nothing in the figures changed.** The inventories are additive releases. No
  panel, no source-data file and no manuscript claim was modified, and the
  frozen-object hashes checked by `audits/verify_scientific_contracts.R` are
  untouched.

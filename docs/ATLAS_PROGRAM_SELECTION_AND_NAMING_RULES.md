# Atlas program selection and naming rules

The governing contract for the curated spatial GO-program atlas. It exists so
that "why this row?" and "does this name fit?" have written answers that do not
depend on who is asking.

> **The atlas rows are curated descriptive biological program families and are
> not independently tested theme-level statistical endpoints.**

There is no theme-level p-value or FDR, and none is implied. All inference lives
at the level of the individual canonical GO term.

---

## Three contracts, never conflated

| | Question | Decided by |
|---|---|---|
| **A. Selection** | Why does this program deserve a row? | recurrence, coherence, distinctness |
| **B. Membership** | Which GO terms are in it? | ontology anchors and approved relationships |
| **C. Naming** | What label describes that membership? | dominant content + supported terms + protein core |

A program can pass A and still need a different name under C. That is exactly
what happened to Chromatin.

---

## A. Selection rule

A primary row is eligible only when **all five** hold:

1. **Biological coherence** — the terms form a coherent family or a defensible
   umbrella with clearly related substructure.
2. **Recurrence** — the biology recurs across multiple canonical spatial
   contexts, not one isolated result.
3. **Distinctness** — it is not an ontology ancestor ladder, a semantic
   duplicate, or a renamed subset of an existing row.
4. **Non-QC status** — not primarily technical, contamination-associated,
   epidermal/keratin, or immune-artefact biology.
5. **Material representation** — omitting it would materially underrepresent a
   recurrent biological concept.

**Selection must never use** desired SUS/RES direction, NES magnitude, how
publishable an example looks, or the recurrent leading-edge proteins.

> Leading-edge proteins are a **downstream naming and interpretation check**,
> never an admission criterion. A program is never admitted because its proteins
> are attractive.

### The ancestor-ladder test

A cluster can look large simply because GO tests a parent and its children
separately. For every candidate we therefore report what fraction of its
supported occurrences comes from terms that are **ancestors of other terms in the
same cluster**. Rejected clusters sat at 57–94%. This is the single most
discriminating check in the selection rule.

---

## B. Membership rule

Membership is generated phenotype-independently from:

- explicit GO-BP anchors in a version-controlled registry;
- the approved relationships `is_a` and `part_of` only;
- explicit ontology-based exclusions where scientifically necessary.

Membership is **never** defined from FDR-supported terms only, from attractive GO
names, from hand-picked individual terms, or from phenotype direction.

### Exclusion rule

A sub-DAG may be excluded only when it corresponds to a distinct biological
process, can be expressed as **one** ontology rule, does not depend on phenotype
statistics, and is documented.

The only exclusion in force: mitochondrial respiration / OXPHOS excludes the
`GO:0006096` glycolytic-process sub-DAG, because cytosolic glycolysis reaches the
theme only through `GO:0045333 cellular respiration`. The rule keeps pyruvate
decarboxylation and the TCA cycle, which a hand-built blacklist would plausibly
have removed. **Arbitrary individual-GO blacklists are forbidden.**

---

## C. Naming rule

A label is acceptable only if:

- **A** it describes the dominant semantic content of the full ontology-defined
  membership;
- **B** it accurately describes the supported constituent GO terms;
- **C** the recurrent leading-edge protein core is biologically compatible;
- **D** **blind recoverability** — with the label hidden, inspecting the supported
  GO terms and recurrent proteins would lead an informed biologist to
  approximately the same interpretation.

### Naming classifications

| Class | Meaning | Action |
|---|---|---|
| **STRONG** | blocks, supported terms and protein core all fit; no contradictory block | KEEP |
| **SUPPORTED_BUT_BROAD** | several subcomponents on one defensible axis, and the label names them | KEEP, state the umbrella |
| **MIXED** | still a real recurrent family, but a substantial component is not represented by the label | REFINE_WORDING before publication |
| **MISLEADING** | dominant evidence would naturally receive a materially different interpretation | RENAME or reconsider the row |

---

## Leading-edge protein rule

A protein counts as **recurrent core** only if it appears in **≥3 supported GO
terms AND ≥3 spatial contexts**. Everything else is recorded separately as
INTERMEDIATE or SINGLE_APPEARANCE.

Leading-edge membership is **not** individual protein significance. None of these
proteins is individually FDR-supported by this analysis, and they must never be
presented as if they were.

---

## Off-theme rule

A term is **not** off-theme merely because its Wang semantic similarity is low,
it sits in a minor cluster, or it is distant from one broad anchor.

Semantic metrics are **diagnostic evidence, not the verdict.** A term is
off-theme only when biological and ontology inspection shows its process is not
reasonably represented by the label.

Every flagged term is recorded with four fields: `semantic_flag`,
`biological_adjudication`, `final_off_theme_status`, `reason`.

> Applied to the current atlas: **13 terms were semantically flagged and 0 were
> biologically off-theme.** Without this rule, legitimate mitochondrial
> processes — pyruvate decarboxylation, proton-motive-force ATP synthesis — would
> have been mislabelled off-theme because of ontology geometry.

---

## Overlap rule

Primary themes need not be mutually exclusive; biologically legitimate overlap is
allowed. But overlap must be **quantified**, major shared branches **disclosed**,
and a unique-terms-only sensitivity check must show the overlap is not
manufacturing the pattern.

RNA processing and Translation / ribosome share six rRNA terms. Excluding every
shared term moves atlas cells by a median of **0.003 NES**, with **0** support-dot
changes. They are **related but non-identical programs** and both are retained.

---

## Completeness

The atlas is **curated, not exhaustive**: the seven rows carry 26.7% of
FDR-supported GO occurrences and 14.0% of unique supported GO identifiers. Those
fractions are low because GO stacks parents and children, not because recurrent
biology is missing — semantic review of everything outside the rows found **no
additional recurrent coherent program** meeting the criteria.

Distinguish **term coverage** from **biological-program coverage**. The atlas is
named for the second.

---

## Status

The seven-row set is **frozen**. Membership changes only to fix a factual
ontology defect. No row 8.

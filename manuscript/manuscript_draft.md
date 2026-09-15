# Manuscript draft

Phase 1: skeleton plus Results §2 and §3. Everything else is a placeholder.

Every quantitative statement in the drafted sections carries a row in
`manuscript/results_statement_provenance.csv`; every claim carries a row in
`manuscript/results_claim_provenance.csv`. Literature citations are `[REF]`
markers at this stage and are not to be invented.

---

## Title

_(candidates under discussion — see the drafting report)_

## Abstract

_[PLACEHOLDER — phase 3. Write last, once Results and Discussion are settled.]_

## Introduction

_[PLACEHOLDER — phase 3.]_

---

# Results

## 1. Adolescent social instability stress produces divergent later behavioural outcomes

_[PLACEHOLDER — phase 2, with Figure 1. Will define the CON / RES / SUS grouping
that every proteomic contrast in §2 and §3 depends on. No proteomic claim in §2
or §3 may characterise the behavioural phenotype beyond the grouping itself.]_

## 2. Spatially resolved hippocampal proteomics recovers reproducible anatomical molecular organisation

Interpreting outcome-associated proteomic differences requires first establishing
that the spatial measurement itself is anatomically faithful. We therefore
characterised the hippocampal proteome across three measurement compartments
before examining any stress variable. Laser-capture microdissection followed by
data-independent-acquisition mass spectrometry yielded 323 spatial acquisitions
from nine animals — 180 neuropil, 71 neuronal soma and 72 microglia-enriched ROI
— at a median of 4,668, 5,307 and 4,850 protein groups identified per
acquisition, respectively (Fig. 2a,b). Only the neuropil was sampled at region ×
layer resolution; the neuronal-soma and microglia-enriched compartments were
sampled at region level, and are treated as region-level throughout. After
filtering, the three analysis matrices comprised 5,054, 5,538 and 5,229 protein
groups, with 4,242 protein groups detected in all three compartments. Throughout,
the biological replicate is the animal (n = 3 per group); acquisitions are
repeated measurements within animals and are never treated as independent
replicates. [METHODS TODO: MT-01, MT-02]

Compartment, not stress group, dominates the global structure of the dataset. In
a joint principal-component analysis of the 4,242 shared protein groups across
all 323 acquisitions, the first two components captured 45.1% and 9.8% of the
variance, and compartment accounted for essentially all of the variance along PC1
(η² = 0.96) whereas experimental group accounted for effectively none (η² = 8.5 ×
10⁻⁵; Fig. 2c). This variance decomposition is computed across acquisitions and
is descriptive; because acquisitions are not independent replicates, no
hypothesis test is attached to it. The ordering is the one a spatial experiment
is designed to produce, with anatomy as the intended source of variation.

Within the three control animals, the proteome defines a coherent anatomical
fingerprint. Using eleven prespecified control-only anatomical contrasts and a
rank-based, phenotype-blind selection rule, 19 protein rows separated the 18
spatial units into anatomically coherent blocks (Fig. 2d), and ten prespecified
compartment markers were most abundant in their expected compartment (Fig. 2e).
Both panels are descriptive and no hypothesis test is applied to them. The
fingerprint rows were ranked from the same control-only contrasts the panel then
displays, so the panel shows the structure those contrasts capture rather than
testing it; no stress information enters either the selection or the row
ordering.

Bilateral concordance varies more by measurement compartment than by anatomical
scale. Comparing the per-protein anatomical contrast estimate obtained
independently from each hemisphere of the same three control animals, left–right
Pearson correlation was highest for neuronal-soma regional contrasts (median
0.80, range 0.72–0.92), intermediate and widest in the neuropil (median 0.74,
range 0.38–0.85) and lowest for the microglia-enriched ROI (median 0.56, range
0.50–0.76; Fig. 2f, Extended Data Fig. 1a). Anatomical scale as such is not the
discriminating factor: dentate-gyrus layer contrasts reproduced as well as the
better regional contrasts (r = 0.81), and CA1 stratum lacunosum-moleculare (r =
0.74) sat at the inventory median. The two lowest values in the inventory were
both CA1 laminar — stratum oriens (r = 0.38) and stratum radiatum (r = 0.49) —
followed by the four microglia-enriched regional contrasts (r = 0.50–0.76). We
therefore interpret CA1 stratum oriens and stratum radiatum distinctions, and
microglia-enriched regional distinctions, more conservatively than neuronal-soma
and dentate-gyrus contrasts, and this asymmetry is left visible rather than
filtered. Consistent with hemispheres being repeated tissue rather than
independent samples, averaging the two hemispheres raised the intraclass
correlation of endpoint scores from a median of 0.33 to 0.50 for compartment
scores and from 0.32 to 0.48 for reference-marker scores (Extended Data Fig. 1b),
which is why all downstream inference is performed on bilaterally aggregated,
animal-level values. [METHODS TODO: MT-13]

Finally, the neuropil and neuronal-soma spatial assignments recover independently
published hippocampal anatomy. Testing nine control-only internal anatomical
contrasts (five neuropil, four neuronal soma; three control animals) against
seven hippocampal subregion and synaptic signatures reported by an independent
study [REF], all ten a-priori expected contrast–signature pairings were recovered
in the expected direction and all ten were supported under the stored
signature-family correction (normalised enrichment score 1.95–3.69 for
positive-direction pairings; Fig. 2g, Extended Data Fig. 2b). No external
anatomical signature was available for the microglia-enriched compartment, which
is therefore not externally anchored. This is the only externally anchored
validation in the study; the gene-set annotation of the same contrasts shown
alongside it (Fig. 2h) uses the same proteomic data and is functional
characterisation, not independent validation. [METHODS TODO: MT-03, MT-04]

Together these results show that the workflow resolves hippocampal molecular
organisation that is reproducible between the two hemispheres of the same animal
— internal consistency rather than independent replication — at the anatomical
scale each compartment supports: region × layer in the neuropil and region level
in the neuronal-soma and microglia-enriched compartments, with externally
anchored support for the neuropil and neuronal-soma assignments, and with CA1
stratum oriens and stratum radiatum the least reproducible distinctions measured.

## 3. Later resilient and susceptible outcomes are associated with sparse protein-level and coordinated, spatially resolved molecular-program differences

We next asked how later resilient and susceptible outcomes are represented in
this spatially resolved proteome. At the level of individual proteins, the answer
is that they are represented sparsely. Across all 18 spatial units, 37 proteins
reached FDR support for the susceptible-versus-resilient contrast, and 12 of the
18 units contained none at all (Fig. 3a). Twenty-eight of the 37 were in a single
neuropil unit, CA2 stratum lacunosum-moleculare.

That apparent concentration weakens substantially under quality-control scrutiny
— from 28 proteins to 6 — and we report it as a qualified rather than a headline
result. CA2-SLM carried the highest pre-imputation missingness of the ten
neuropil units (6.2–24.8% across acquisitions) and that missingness was itself
unequal between groups. Because
per-sample median centring is sensitive to differential missingness, this
produced a systematic normalisation displacement between groups (−0.146 in
susceptible relative to resilient animals), and across acquisitions the
displacement tracked missingness almost exactly (Pearson r = 0.93). Applying the
prespecified robustness criteria, 6 of the 28 CA2-SLM proteins qualified; the 9
FDR-supported proteins outside CA2-SLM were never exposed to this artefact and
enter unchanged, giving 15 robustness-qualified proteins in total (Fig. 3a,
Extended Data Fig. 3). "Robustness-qualified" therefore means not excluded by the
CA2-SLM missingness audit, not passed an additional test. [METHODS TODO: MT-05]

Coordinated differences at the level of molecular programs are detectable,
including in units where no individual protein reaches FDR support. To
summarise them we mapped the canonical gene-set enrichment results onto a curated
atlas of seven ontology-defined program families — RNA processing / splicing /
RNP organisation, translation / ribosome biogenesis, chromatin organisation /
epigenetic regulation, mitochondrial respiration / oxidative phosphorylation,
synaptic signalling / vesicle-mediated transport, neuron projection development,
and autophagy / endolysosomal trafficking (253 constituent GO biological-process
terms; Fig. 3b, Extended Data Fig. 6a,b). FDR-supported constituent terms occur
in 132 of the 378 theme × unit × contrast cells, and as observed those supported
cells are distributed unevenly across spatial contexts and across the three
contrasts; no test of that unevenness was performed. The atlas is a
descriptive summary, not a test: cell colour is the median normalised enrichment
score of a family's constituent terms, no theme-level P value or FDR is computed
or implied, and the seven families deliberately capture a curated subset — 26.7%
of FDR-supported GO occurrences and 14.0% of unique supported GO identifiers —
with the complete results provided as source data. [METHODS TODO: MT-06, MT-07]

Three representative programs, one per measurement compartment, are shown (Fig.
3c–f). Each was chosen editorially from the terms already FDR-supported for the
susceptible-versus-resilient contrast in that compartment, so they illustrate the
form these program differences take rather than constituting an unbiased sample.
Throughout, gene-set FDR is conditional on the ranked per-gene contrast statistic
and is a statement about gene ranks rather than about the three animals per
group; enrichment P values are floored at the method tolerance of 1 × 10⁻¹⁰, and
each of the three terms below sits at that floor in at least one displayed cell.
A synaptic-signalling program in CA3 stratum radiatum neuropil (GO:0099536) was
lower in susceptible animals relative to both resilient animals and controls (NES
−1.72, FDR 4.1 × 10⁻⁸; NES −1.72, FDR 7.4 × 10⁻⁶), with no support for the
resilient-versus-control arm (FDR 0.19). An mRNA-processing program in CA2
neuronal soma (GO:0006397) showed the same asymmetry in the opposite direction
(NES +2.13, FDR 1.9 × 10⁻⁷ and NES +1.62, FDR 0.028; resilient-versus-control
FDR 0.78). Both are therefore susceptibility-associated. By contrast, a reduced
oxidative-phosphorylation program in the CA1 microglia-enriched ROI (GO:0006119)
was supported in all three contrasts in the same direction (NES −2.11, −2.76 and
−1.78; FDR 5.8 × 10⁻⁵, 1.5 × 10⁻⁸ and 3.1 × 10⁻³) and is better described as a
graded, stress-associated direction than as an outcome-specific one. These three
terms were selected as one illustrative example per compartment and are not the
strongest result in their own units. The microglia-enriched ROI is an enriched
measurement context rather than a purified population, so this result cannot
establish a cell-intrinsic microglial property. The three displayed pairwise
contrasts are algebraically related and are not independent replications.
[METHODS TODO: MT-08, MT-09]

Three further analyses bound how strongly these program-level results should be
read. Under a correlation-aware competitive gene-set sensitivity analysis applied
to the identical ranked statistics, directional structure was highly concordant
with the primary analysis (median per-comparison Spearman ρ = 0.93; all 3,559
FDR-supported terms direction-concordant) but inferential strength was attenuated,
with 22.1% retaining FDR support, and, for the susceptible-versus-resilient arm
displayed in Fig. 3, all three representative terms reproduced their direction
without retaining FDR support there (CAMERA FDR 0.10, 0.11 and 0.14); this is a
sensitivity analysis on the same data and not independent validation. Separately,
among the module × contrast cells reported in Extended Data, none survived
correction (0 of 45; smallest FDR 0.25) and no stress × spatial-unit interaction
omnibus test was FDR-supported (0 of 35 across the three compartments; smallest
FDR 0.27), so we do not claim a tested spatial heterogeneity of outcome effects.
At the network level, no whole-network group difference was detectable at three
animals per group (exact permutation P = 0.58, 0.80 and 0.74 for the three
compartments, against an attainable floor of 0.004), and among the eight neuropil
spatial-unit pairs tested against behavioural readouts no edge–behaviour
association survived correction (0 of 48 tests; smallest BH-adjusted P = 0.43);
the neuronal-soma and microglia-enriched compartments were not tested in that
analysis. Selected leading-edge proteins are
shown to expose which proteins carry each enrichment signal (Fig. 3g–i); none of
the 63 displayed values is individually FDR-supported at the protein level
(smallest BH FDR 0.53), and they are descriptive rather than independent
confirmation. [METHODS TODO: MT-10, MT-11, MT-12]

Taken together, later resilient and susceptible outcomes are associated with
sparse individual-protein differences and with coordinated molecular-program
differences that are resolved across hippocampal spatial contexts. The two are
assessed in separate multiple-testing families and were never placed on a common
scale, so neither is claimed to be the stronger.

## 4. Integration of behavioural outcome with the spatial proteome

_[PLACEHOLDER — phase 2 or later, only if the evidence warrants a separate
section. The edge–behaviour coupling analysis did not survive correction, so this
section may reduce to a limitation stated in the Discussion rather than a
standalone Results section. Do not invent a positive integration result.]_

---

# Discussion

_[PLACEHOLDER — phase 3.]_

# Methods

_[PLACEHOLDER — phase 2. See `manuscript/methods_todo.csv` for the structured
list of method topics the drafted Results already depend on.]_

# Data availability

_[PLACEHOLDER — PRIDE accession pending; see `docs/PRIDE_EXPORT.md`.]_

# Code availability

_[PLACEHOLDER — repository reference and the canonical entrypoints listed in
`docs/CANONICAL_ANALYSIS_ENTRYPOINTS.md`.]_

# Author contributions

_[PLACEHOLDER]_

# Acknowledgements

_[PLACEHOLDER]_

# References

_[PLACEHOLDER — all citations are `[REF]` markers at this stage.]_

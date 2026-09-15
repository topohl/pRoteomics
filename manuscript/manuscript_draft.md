# Manuscript draft

Phases 1–2C: Results §1, §2 and §3 and the Methods they depend on are drafted.
Results §4 is closed rather than pending. The Abstract, Introduction and
Discussion remain placeholders for phase 3.

Every quantitative statement in the drafted sections carries a row in
`manuscript/results_statement_provenance.csv`; every claim carries a row in
`manuscript/results_claim_provenance.csv`. Literature citations are `[REF]`
markers at this stage and are not to be invented.

Results §1 and the behavioural Methods are quoted from a frozen evidence bundle
imported from the upstream behavioural repository and mirrored byte-identically
in `manuscript/figure1_bridge_mmmsociability/`. No behavioural statistic is
computed in this repository. Provenance is in
`manuscript/figure1_bridge_provenance.csv` and
`manuscript/figure1_bridge_import_manifest.csv`; every point at which the bundle,
the upstream code or an older statement here disagreed is recorded in
`manuscript/figure1_bridge_conflicts.csv`.

---

## Title

_(candidates under discussion — see the drafting report)_

## Abstract

_[PLACEHOLDER — phase 3. Write last, once Results and Discussion are settled.]_

## Introduction

_[PLACEHOLDER — phase 3.]_

---

# Results

## 1. Early spontaneous home-cage activity predicts later composite stress outcome

Adolescent social instability produces outcomes that differ markedly between
individuals, and the question that motivates this work is whether that later
divergence is foreshadowed by behaviour recorded before the divergence exists. We
therefore separated the measurement timeline into a single early observation
window and a set of later outcome measures, with no overlap between them.
Radio-frequency identification tracking of undisturbed home-cage activity began at
postnatal day 25, in the first active phase following the first cage change, and
covered a fixed 12-h window from 18:30 to 06:30 in 10-min bins (72 slots;
Fig. 1b). Every component of the later outcome — novel-object recognition,
sucrose preference, weight deviation, delta corticosterone, adrenal weight and
spleen weight — was collected after this window had closed, as were the composite
score derived from those components and the resilient/susceptible labels derived
from that score. The predictor therefore precedes the outcome and the
classification by construction: neither existed at the time of recording.

Later outcome was summarised as a composite z-score (CombZ), the unweighted mean
of those six components, each z-scored against same-sex control animals and with
delta corticosterone, adrenal weight and spleen weight sign-inverted so that all
six point in the same direction. Higher CombZ indicates a more resilient-like
outcome; lower CombZ corresponds to a greater later stress burden. Stress-exposed
animals were classified as susceptible when CombZ fell more than one control
standard deviation below the same-sex control mean, and resilient otherwise;
control animals were never relabelled. Because these six measures define the
composite and the composite defines the classification, differences between the
resulting groups in those same measures are guaranteed by construction, and we do
not present them as independent confirmation of the phenotype (Fig. 1a).

Within this design, mean movement over the early window was negatively associated
with later CombZ (Spearman ρ = −0.39, 95% CI [−0.55, −0.21], q = 6.9 × 10⁻⁵;
n = 111 animals, 58 female and 53 male, comprising 24 control, 49 resilient and 38
susceptible; Fig. 1c). Given the orientation of the score, the negative sign means
that animals that were more active during the first undisturbed night after the
cage change tended towards a lower later CombZ — that is, towards a less
resilient-like outcome, corresponding to a greater later stress burden. The
short-timescale variability of the same signal, movement RMSSD, was associated in
the same direction but more weakly (ρ = −0.23, q = 0.026). A third prespecified
feature, the lag-one autocorrelation of binned activity entropy, did not reach
FDR support (ρ = −0.18, q = 0.067, with a bootstrap interval including zero) and
is not interpreted further.

To ask whether this association carries prospective information about individual
animals rather than only about the group, we used a model registry fixed before
fitting, in which the primary model takes mean movement as its sole predictor.
Held-out performance was estimated by leave-one-animal-out cross-validation,
refitting the model completely for each of the 111 animals and evaluating it on
the animal withheld. This movement-mean model explained approximately 16% of the
variation in later CombZ in held-out animals (R² = 0.159), against −0.018 for an
intercept-only baseline. A repeated grouped five-fold scheme with the animal as
the grouping unit (k = 5, 100 repeats) gave a closely matching estimate (mean
R² = 0.156; 2.5th–97.5th percentile range across repeats 0.116–0.179), and a
permutation test that repeated the entire fitting and cross-validation procedure
under permuted outcomes placed the observed value beyond every one of 1,000 draws
(p = 1/1001; Fig. 1d). No feature selection was performed, no model was chosen on
observed performance, and no outcome-derived label entered any model. Adding the
two remaining features, or sex, did not improve on mean movement alone. This is
internal validation: performance was estimated by withholding animals within a
single cohort, not by testing in an independent cohort.

The relationship did not differ detectably by sex. Formal feature-by-sex
interaction tests were unsupported for all three features (all q = 0.90), and the
sex-stratified estimates for mean movement, which are descriptive rather than a
test of difference, were near-identical (ρ = −0.41 in 58 females and −0.42 in 53
males; Fig. 1e).

Two limitations bound this result. Cage identity is not represented in the
analysis design, so cage-level dependence could be neither modelled nor assessed
retrospectively; this concerns how far the estimate generalises beyond the cages
sampled rather than offering any route by which outcome information could have
reached the predictor. And because validation is internal, the held-out estimate
may be optimistic with respect to structure shared within the cohort.

Early spontaneous behaviour therefore carries prospective information about where
an animal's later composite outcome will fall. What that later outcome
corresponds to in the brain is a separate question, and it is the one the
remainder of this work addresses. The proteomic analyses that follow use the
resilient and susceptible assignments only as group labels, in the nine animals
that contributed hippocampal tissue; they do not model `CombZ` as a continuous
variable and they inherit no part of the prediction analysis above. We first
establish what the spatially resolved proteome measures, and then ask how later
resilient and susceptible outcomes are represented within it.

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
no co-abundance module–outcome contrast survived correction — 0 of 45 module ×
contrast cells in the neuropil (smallest FDR 0.25), and 0 of 105 across all three
compartments (smallest FDR 0.16) — and no stress × spatial-unit interaction
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

_[CLOSED — no section will be written. Direct integration of the behavioural
outcome with the spatial proteome reaches no result at FDR < 0.05 in either
repository: the edge–behaviour coupling analysis here did not survive
correction, and the upstream behavioural analysis records the same conclusion as
finding BH-006. This is recorded as a limitation in the Discussion rather than
as a Results section. Do not fill this placeholder with an unsupported
integration result; the absence is the finding.]_

---

# Discussion

_[PLACEHOLDER — phase 3.]_

# Methods

Every parameter below traces to a row of
`manuscript/methods_statement_provenance.csv`. Items that could not be recovered
from this repository are marked `[METHOD DETAIL UNRESOLVED]` rather than filled
in from convention.

## Animals and experimental design

Nine animals contributed spatial proteomics: three control, three later
resilient and three later susceptible. **The biological replicate is the animal
throughout** (n = 3 per group). Spatial acquisitions, hemispheres and
animal × dataset network instances are repeated measurements and are never
treated as independent replicates. [M-15]

## Provenance of the behavioural analysis

The behavioural analysis was performed in a separate repository and is not
reproduced here. Every behavioural quantity reported in this work is quoted from
a frozen evidence bundle imported into `manuscript/figure1_bridge_mmmsociability/`
and recorded, with per-file hashes, in
`manuscript/figure1_bridge_import_manifest.csv` and
`manuscript/figure1_bridge_provenance.csv`. No behavioural statistic was
recomputed in this repository and no behavioural analysis code was copied into
it. The bundle derives from analysis commit `4b0f90f`, was frozen at commit
`53bc7e9`, and is byte-identical to the state at the verified source head
`a53d73f`. [M-19]

## Adolescent social instability stress

Animals underwent an adolescent social-instability paradigm consisting of
repeated changes of cage composition, beginning with the first cage change at
postnatal day 25. All later outcome assessment followed the paradigm.
_[PLACEHOLDER — phase 3: the number, spacing and composition rule of the cage
changes and the duration of the paradigm are design parameters and are not
specified in the frozen behavioural bundle, which fixes the analysis rather than
the husbandry protocol. To be supplied from the experimental record.]_ [M-20]

## Experimental timeline and early home-cage recording window

Home-cage activity was recorded continuously by radio-frequency identification
tracking. The early window used as the predictor is the first active-phase block
following the first cage change at postnatal day 25: a fixed clock window from
18:30 inclusive to 06:30 exclusive, 12 h in total, binned at 10 min to give an
expected 72 slots per animal. The window is defined on clock time, and
inactive-phase bins were never included. Seventy-two slots is the design
expectation rather than the realised coverage: 50 of 111 animals contributed all
72, and the remaining 61 were missing only leading slots, with no interior or
trailing gaps in any animal (mean coverage 98.6%, minimum 94.4%). Every component
of the later outcome, the composite score and the resilient/susceptible
classification derive from measurements taken after this window closed, so the
predictor precedes both the outcome and the group labels. [M-21]

## Composite outcome score

Later outcome was summarised as a composite z-score, `CombZ`, defined as the
unweighted mean of six components: novel-object recognition, sucrose preference,
weight deviation, delta corticosterone, adrenal weight and spleen weight. Each
component was z-scored against the control animals of the same sex, using twelve
control animals per sex and the population standard deviation. Delta
corticosterone, adrenal weight and spleen weight were sign-inverted before
averaging so that all six components share an orientation. Components contribute
with equal weight, one sixth each, and the mean is taken over the components
available for a given animal (`na.rm = TRUE`), so an animal missing a component is
scored on the remainder rather than dropped. No batch term enters the
construction. Higher `CombZ` denotes a more resilient-like outcome. [M-22]

## Resilient and susceptible classification

Among stress-exposed animals, an animal was classified susceptible when its
`CombZ` fell below the mean of same-sex control animals minus one same-sex
control population standard deviation, and resilient otherwise. The resulting
thresholds are −0.436641698 for males and −0.222390844 for females. Control
animals were never relabelled. Stored labels and labels reconstructed from the
rule agree for every stress-exposed animal. [M-23]

## Early behavioural features

Three features, fixed in advance, summarise the early window: `Movement_mean`,
the mean of the binned movement signal; `Movement_rmssd`, the root mean square of
successive differences of the same binned signal, a short-timescale variability
measure; and `Entropy_acf1`, the lag-one autocorrelation of the binned activity
entropy. All three are raw summaries computed directly from the window. No
scaling or transformation was applied, and no feature derived from a generalised
additive mixed model entered the analysis. [M-24]

## Association between early behaviour and later outcome

Each feature was related to `CombZ` by Spearman rank correlation, with 95%
confidence intervals from 5,000 non-parametric bootstrap resamples over animals
(percentile method, seed 123) and Benjamini–Hochberg correction across the three
prespecified features. The analysis population is 111 animals: 58 female and 53
male; 24 control, 49 resilient and 38 susceptible. [M-25]

## Out-of-sample prediction

The set of candidate models was fixed in a registry before any model was fitted,
and no model was selected on observed performance. The primary model, `movement_mean`,
takes `Movement_mean` as its sole predictor. Out-of-sample performance was
estimated by leave-one-animal-out cross-validation: for each of the 111 animals
the model was refitted completely on the remaining 110 and evaluated on the
withheld animal. A repeated grouped five-fold cross-validation was run as a
companion, with five folds, 100 repeats, the animal as the grouping unit and seed
521 for fold assignment; fold integrity was asserted rather than assumed, with
exactly one fold per animal per repeat. Performance is reported as R² against the
continuous `CombZ` target, with an intercept-only model as baseline. For the
repeated scheme, each repeat yields one out-of-fold R², and the reported spread is
the 2.5th–97.5th percentile range across the 100 repeat-level values; it is a
resampling range and is not a confidence interval. Significance was assessed by
permuting the outcome and repeating the complete fitting and cross-validation
procedure for each of 1,000 draws (seed 20260811), giving p = 1/1001 with no null
draw reaching the observed value. No feature selection, centring or scaling was
applied to the canonical models, and no outcome-derived group label was used as a
predictor in any of them; missing predictor values were imputed with the
training-fold median inside each split, so no information crosses a
cross-validation boundary. [M-26]

## Sex

Sex was examined by a formal feature-by-sex interaction term, fitted as a linear
model of `CombZ` on the feature, sex and their product, one model per feature,
with Benjamini–Hochberg correction across the three tests. Sex-stratified
correlations are reported as descriptive summaries and were not used to infer a
difference between sexes; where a within-sex correlation is quoted it is
uncorrected. [M-27]

## Interpretational constraints on the behavioural analysis

Three constraints are stated explicitly because they bound what the behavioural
result can support. First, the six outcome components define `CombZ` and `CombZ`
defines the resilient/susceptible classification, so differences between those
groups in those same components are guaranteed by construction; they are used to
describe the classification and are never treated as independent validation of
it. Second, cage identity is not represented in the analysis design, so
cage-level dependence was neither modelled nor assessable retrospectively from
this analysis; this is a limitation on generalisation rather than evidence of
information leaking from outcome to predictor, and no cage random effect was
introduced after the fact. Third, validation is internal — animals were withheld
within a single cohort — and the estimate may therefore be optimistic with
respect to structure shared within that cohort. A fourth point is recorded for
reproducibility rather than interpretation: the individual component z-scores
cannot all be regenerated under a single uniform derivation rule because the
source workbook stores them in positional per-sex blocks, whereas the downstream
chain from components to `CombZ` to group labels reproduces to numerical
precision (maximum absolute deviation 4.4 × 10⁻¹⁶). [M-28]

## Spatial proteomics sample collection and acquisition

Laser-capture microdissection followed by data-independent-acquisition mass
spectrometry yielded 323 spatial acquisitions (180 neuropil, 71 neuronal soma, 72
microglia-enriched ROI). Ten neuropil units carry region × layer resolution; the
neuronal-soma and microglia-enriched compartments are region-level only, giving
18 spatial units in total. _[METHODS TODO MT-01: instrument, gradient, DIA window
scheme and search settings to be added from the acquisition records.]_

## Preprocessing and animal-level bilateral aggregation

Protein groups were filtered at 70% missingness and imputed upstream, giving
analysis matrices of 5,054 (neuropil), 5,538 (soma) and 5,229
(microglia-enriched) protein groups, with 4,242 detected in all three
compartments. Technical rows (protein group × animal × hemisphere × spatial unit)
were aggregated to one animal-level value per protein group and spatial unit
before any group comparison, by equal-weight mean of the left and right
hemisphere where both were available, retaining the single available hemisphere
where one was missing. No hemisphere was imputed. Hemispheres are repeated tissue
within an animal, not independent biological replicates. [M-01]

## Differential protein abundance

Protein-level differential abundance was assessed with a moderated linear model
on animal-level values within each spatial unit (contract
`animal_level_protigy_da_v1`), with Benjamini–Hochberg correction **within each
comparison**. The primary contrast is susceptible versus resilient; resilient
versus control and susceptible versus control are reported alongside it and are
algebraically related to it (SUS−RES = SUS−CON − RES−CON). [M-02]

Proteins reaching FDR support in CA2 stratum lacunosum-moleculare were
additionally examined against prespecified robustness criteria following a
missingness and normalisation audit of that unit. **This qualification is not a
second FDR family and not an additional hypothesis test**: "robustness-qualified"
means not excluded by that audit. Proteins outside CA2-SLM were not exposed to
the artefact and enter unchanged. _[METHODS TODO MT-05: state the prespecified
thresholds.]_

## Ranked GO enrichment

Genes were ranked by the median moderated *t* statistic per official gene symbol
and tested against GO biological process with `clusterProfiler::gseGO`
(`minGSSize` 10, `maxGSSize` 800, Benjamini–Hochberg **within each comparison**),
across 54 comparisons. Enrichment *P* values are floored at the method tolerance
`eps = 1 × 10⁻¹⁰`, which is the `clusterProfiler` default and was never
overridden; a term reported at that value has a true *P* somewhere below it that
the method does not resolve, so its FDR bounds the evidence rather than measuring
it. Gene-set FDR is conditional on the ranked per-gene statistic and is a
statement about gene ranks, not about the three animals per group. [M-03, M-04,
M-05]

## Curated GO-program atlas

Canonical GO-BP results were mapped onto seven ontology-defined program families
using an anchor-based registry (`manuscript_go_themes_v3`), in which each family
is defined by one or more GO anchors matched by `anchor_and_descendants`,
`exact_go_id` or `exclude_anchor_and_descendants`. **Membership is defined by
ontology structure and is independent of phenotype.** The mitochondrial
respiration / oxidative phosphorylation family comprises 16 GO terms; the
glycolysis sub-DAG (`GO:0006096`) is excluded by an explicit exclusion rule,
while pyruvate decarboxylation to acetyl-CoA (`GO:0006086`) and the
tricarboxylic acid cycle (`GO:0006099`) are retained. Family colour is the
**median normalised enrichment score** of the constituent canonical terms; a
support marker indicates that at least one constituent term passed its own
prespecified FDR threshold. **The aggregation constitutes no additional
multiple-testing family and no theme-level *P* value or FDR is computed or
implied.** [M-06, M-07]

## CAMERA sensitivity analysis

As a sensitivity analysis, the identical canonical ranked statistics were tested
with `limma::cameraPR` (preranked competitive gene-set test) with the inter-gene
correlation fixed at the prespecified value of 0.01. This is a concordance check
on the same data and **is not independent validation, replication or
confirmation**. [M-08]

## WGCNA co-abundance modules

Module eigengenes were related to outcome group with a linear mixed model,
`eigengene ~ StressGroup + SpatialUnit + (1 | AnimalID)` (`lmerTest`), fitted
separately per compartment. Each compartment carries its own Benjamini–Hochberg
families: a primary family (modules × susceptible-versus-resilient), a secondary
family (modules × the two control contrasts) and an interaction omnibus family
(one StressGroup × SpatialUnit test per module). In the neuropil, which has 15
modules, the primary and secondary families together comprise **45 module ×
contrast cells, of which none is FDR-supported (smallest FDR 0.245)**. Pooling
the equivalent cells across all three compartments (35 modules) gives **105
tests, none supported, smallest FDR 0.156**. The interaction omnibus family
comprises **35 tests across the three compartments, none supported, smallest FDR
0.274**. [M-09, M-10]

## Spatial molecular-identity analyses

Bilateral concordance was computed in control animals only by fitting each
prespecified anatomical contrast separately in each hemisphere and correlating
the per-protein contrast estimates. Eleven contrasts are manuscript-locked (seven
neuropil, four neuronal soma); four microglia-enriched regional contrasts are
reported as context and are not part of the locked set. The two dentate-gyrus
layer contrasts are algebraic mirror images, because dentate-gyrus neuropil has
two layers, and are not independent. No hypothesis test is applied to these
concordance metrics. [M-14]

## External-reference validation

Nine control-only internal anatomical contrasts (five neuropil, four neuronal
soma; no microglia-enriched contrast has an external counterpart) were tested
against seven published hippocampal subregion and synaptic signatures [REF],
giving 30 contrast–signature pairings: ten designated a priori as expected
anatomical correspondences and twenty as off-target comparisons. **Each pairing
was run as a separate gene-set enrichment test against a single-signature
collection.** Because that collection contains one term, the per-pairing
Benjamini–Hochberg adjustment is a no-op and the stored `p_adjust` field is
numerically identical to the raw *P* value; **the multiple-testing correction
that Methods and Results rely on is the signature-family FDR**, applied within
three families of 12 (soma tissue), 6 (neuropil subregion) and 12 (CA1 laminar)
tests. Off-target pairings are reported for completeness; **no formal
expected-versus-off-target discrimination test was performed**, and 18 of the 20
off-target pairings also clear the threshold, so these comparisons should not be
read as establishing anatomical specificity. [M-11, M-12, M-13]

## Protein co-abundance network analysis

Animal-level spatial networks were compared with a leave-one-control-animal-out
consensus using exact permutation over group labels. Edge–behaviour associations
were tested for eight neuropil spatial-unit pairs against behavioural readouts by
Pearson correlation across the nine animals. _[METHODS TODO MT-12: state the
edge-definition rule and the behavioural readouts used.]_

## Statistics and reproducibility

All inference is on bilaterally aggregated, animal-level values with n = 3 per
group. Benjamini–Hochberg families are declared per analysis above and are never
pooled across analyses; protein-level and program-level results sit in separate
families that were never placed on a common scale. Null results are reported as
failure to survive correction at this sample size, never as evidence of absence.
_[METHODS TODO MT-01: software and package versions, seeds and the ontology
release used.]_

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

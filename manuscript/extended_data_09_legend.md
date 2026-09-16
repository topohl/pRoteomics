# Extended Data Figure 9 legend

Draft legend. Every number below is a cell of
`manuscript/figure1_bridge_mmmsociability/behavior_prediction_model_ladder.csv`
or
`manuscript/figure1_bridge_mmmsociability/behavior_sex_effect_contract.csv`,
both frozen upstream in `topohl/MMMSociability` at commit `53bc7e9` and imported
byte-identically with SHA-256. All statistics were computed in the upstream
behavioural repository; nothing on this figure is calculated in the proteomics
repository.

This figure exists to support two negative or non-selective statements made in
Results §1, neither of which Figure 1 can carry on its own: that the model shown
in Fig. 1d was fixed before fitting rather than chosen for its performance, and
that the association in Fig. 1c did not differ detectably by sex.

---

**Extended Data Figure 9 | The complete a-priori model registry, and the
feature-by-sex interaction tests.**

**(a)** All five models in the registry, in registry order. The registry was
fixed before any model was fitted; the rows are not ordered by performance and
no test compares one model with another. Filled point, leave-one-animal-out
estimate, which is the primary quantity; bar, the 2.5th–97.5th percentile range
across repeated grouped five-fold splits (k = 5, 100 repeats, animal as the
grouping unit), with the repeated-cross-validation mean marked as a tick. The bar
belongs to the repeated-cross-validation mean and is not a confidence interval
for the leave-one-animal-out point, which is why the two are drawn in different
ink. All models are fitted on the same n = 111 animals. Leave-one-animal-out
R² was −0.018 for the intercept-only baseline, 0.159 for mean movement alone,
0.152 for mean movement with movement RMSSD and entropy ACF1, 0.152 with sex
added to mean movement, and 0.142 with sex added to the three-feature model. The
four behavioural models differ from one another by less than the width of any one
of their intervals, so the figure states that adding features or sex did not
improve on mean movement alone and makes no stronger claim than that. A
full-refit outcome permutation was run for the two behaviour-only models, each
giving P = 1/1001 against 1,000 draws; it was not run for the baseline or for the
two sex-adjusted sensitivity models, and those rows are printed as "not run"
rather than left blank.

**(b)** The formal test of whether the early-behaviour association differs by
sex: one feature-by-sex interaction per early feature, with Benjamini–Hochberg
correction across the three. Point, interaction estimate; bar, 95% confidence
interval; q, the corrected value. Mean movement −0.015 (95% CI −0.238 to 0.208,
q = 0.90), movement RMSSD 0.052 (−0.162 to 0.266, q = 0.90), entropy ACF1 0.994
(−1.176 to 3.165, q = 0.90). Zero lies inside every interval and no interaction
is supported.

**(c)** The sex-stratified correlations of each early feature with later `CombZ`,
shown only because panel b is a negative result and a reader is entitled to see
the estimates the negative statement is about. Filled circle, females (n = 58);
open triangle, males (n = 53); the connecting rule spans the two. Mean movement
ρ = −0.41 in females and −0.42 in males, movement RMSSD −0.25 and −0.24, entropy
ACF1 −0.27 and −0.05. These are descriptive: they are not tests, they carry no
multiplicity correction, and the difference between any pair of them has not been
tested other than by the interaction tests in panel b, which are unsupported.

---

## Wording constraints applied

The upstream contract carries an explicit prohibited-wording list for the sex
panels: *female-specific*, *sex-specific*, *stronger in females*, *the effect was
driven by females*. None appears in this legend, in Results §1, or on the figure.
The permitted statement is that the association did not differ detectably by sex
and that the stratified estimates are descriptive.

The producer script refuses to render panel c unless every row of the sex
contract is still classified `FORMAL_INTERACTION_NOT_SUPPORTED`. If an
interaction ever becomes supported upstream, the render fails rather than
quietly redrawing a panel whose whole justification has changed.

Panel a states no winning model. "Improvement" is used only in the negative, and
the panel prints no ranking, no delta and no comparison test, because the
registry contains none.

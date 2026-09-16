# Figure 1 legend

Draft legend. Every number below is recoverable from
`manuscript/figure1_bridge_mmmsociability/source_data/figure1_panel_statistics.csv`,
which is the contract for what this figure may state. All statistics were
computed in the upstream behavioural repository and frozen; nothing on this
figure is calculated in the proteomics repository.

---

**Figure 1 | Early spontaneous home-cage activity predicts later composite stress
outcome.**

**(a)** Experimental timeline. Radio-frequency identification tracking of
undisturbed home-cage activity was recorded in the first active-phase block
following the first cage change at postnatal day 25 (P25; CC1), across a fixed
clock window from 18:30 inclusive to 06:30 exclusive — 12 h, binned at 10 min, giving
72 expected slots per animal. Seventy-two is the design expectation rather than
the realised coverage: 50 of 111 animals contributed all 72 slots and the
remaining 61 were missing leading slots only, with no interior or trailing gaps
in any animal (mean coverage 98.6%). Every component of the later outcome, the
composite score derived from them and the resilient/susceptible labels derived
from that score were obtained after this window closed, so the classification did
not exist at the time of recording.

**(b)** Later composite outcome and the rule that defines the phenotype. Each
point is one animal (n = 117), plotted by sex because the classification is
referenced within sex. `CombZ` is the unweighted mean of six components — novel
object recognition, sucrose preference, weight deviation, delta corticosterone,
adrenal weight and spleen weight — each z-scored against same-sex control animals
using the population standard deviation, with the last three sign-inverted so
that all six share an orientation. Higher `CombZ` indicates a more resilient-like
outcome. Grey line, same-sex control mean; dashed line, the susceptibility
threshold at one control population standard deviation below that mean
(−0.437 male, −0.222 female). Stress-exposed animals below the threshold were
classified susceptible and the remainder resilient; control animals were never
relabelled. This panel shows how the groups were defined. Because these six
components construct `CombZ` and `CombZ` constructs the classification,
differences between the resulting groups in those components are guaranteed and
are not presented here as independent confirmation of the phenotype.

**(c)** Early mean movement over the window in (a) against later `CombZ`; one
point per animal, n = 111 (58 female, 53 male; 24 control, 49 resilient, 38
susceptible). Spearman ρ = −0.39, 95% confidence interval [−0.55, −0.21] from
5,000 percentile bootstrap resamples, q = 6.9 × 10⁻⁵ after Benjamini–Hochberg
correction across three prespecified features. Given the orientation of the
score, the negative sign means that animals more active during the first
undisturbed night tended towards a lower later `CombZ`, that is, towards a less
resilient-like outcome. No model is fitted in this panel and no line is drawn
through the points; the reported statistic is a rank correlation. The panel is
deliberately not stratified by sex, because the formal feature-by-sex interaction
is unsupported (all q = 0.90).

**(d)** Held-out prediction of continuous `CombZ`. Each point is one animal, with
its predicted value obtained from a model refitted on the other 110 animals and
evaluated on that animal alone (leave-one-animal-out, n = 111). The model is the
prespecified movement-mean model, whose sole predictor is early mean movement;
it was fixed in a registry before fitting and carries no outcome-derived term.
Grey line, identity, not a fit. Leave-one-animal-out R² = 0.159 against −0.018
for an intercept-only baseline. The prediction target is the continuous score
throughout; no classifier was fitted, so no accuracy or area under the curve
exists and none is reported.

**(e)** Null distribution for the held-out result. The complete fitting and
cross-validation procedure was repeated under 1,000 permutations of the outcome,
refitting every model in full for each draw. Histogram, the 1,000 permuted
leave-one-animal-out R² values; vertical line, the observed value of 0.159, which
no permuted draw reached (p = 1/1001). A repeated grouped five-fold scheme with
the animal as the grouping unit gives a closely matching estimate (R² = 0.156;
2.5th–97.5th percentile range across 100 repeats 0.116–0.179 — a resampling range
across repeats, not a confidence interval).

The biological unit is the animal throughout. Validation is internal: animals
were withheld within a single cohort, not tested in an independent cohort. Cage
identity is not represented in the analysis design, so cage-level dependence
could be neither modelled nor assessed retrospectively; this bounds how far the
estimate generalises beyond the cages sampled and is not a route by which outcome
information could have reached the predictor.

---

## Wording constraints applied

| Constraint | Applied |
|---|---|
| prediction target is continuous `CombZ` | stated in (d) and (e) |
| never "predicts susceptibility" or "predicts resilience" | absent |
| never "movement-only" | the model is named the movement-mean model |
| never "independent" or "external" validation | stated as internal in the closing paragraph |
| no female-specific or sex-specific claim | (c) states the interaction is unsupported |
| Entropy ACF1 not presented as supported | absent from the figure entirely |
| RMSSD not promoted to equal status | absent from the figure entirely |
| 72 slots is the design expectation | (a) gives design and realised coverage separately |
| classification components not shown as validation | stated explicitly in (b) |
| repeated-CV range is not a confidence interval | stated explicitly in (e) |

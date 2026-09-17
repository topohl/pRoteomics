# Spatial systems foundations

Two contracts that the later spatial/cell atlas and animal-level network stages
depend on: a genuinely hemisphere-resolved data layer, and a phenotype-blind
reusable external cell-type annotation layer. Everything here is either
phenotype-blind or explicitly labelled as a phenotype overlay.

## 1. The biological replicate is AnimalID

There are 9 animals: 3 CON, 3 RES, 3 SUS. Every inferential statement about a
group is a statement about three animals.

One consequence is worth stating numerically, because it bounds what any later
edge-level or network-level test can show. An exact label permutation of SUS vs
RES has `C(6,3) = 20` assignments, so the smallest attainable two-sided p-value
is `2/20 = 0.10`. **No exact SUS-vs-RES permutation test can reach p < 0.05.**
The three-group CON/RES/SUS layout has 1,680 assignments and a floor of 0.0006.

## 2. A hemisphere is a repeated tissue sample, not a replicate

The metadata column recording the side is called `ReplicateGroup`, with values
`Left` and `Right`. That name is a hazard: a hemisphere is **not** a replicate
in any biological sense. The spatial systems layer therefore normalises it once
into a variable called `Hemisphere` with values `L` / `R`, and records
`source_hemisphere_field = "ReplicateGroup"` as provenance on every level.

Bilateral coverage is near complete:

| dataset | AnimalID x SpatialUnit cells | bilateral | one-sided |
|---|---|---|---|
| neuron_neuropil | 90 | 90 | 0 |
| microglia | 36 | 36 | 0 |
| neuron_soma | 36 | 35 | 1 (left only) |

## 3. The aggregation hierarchy

`R/data_contracts/spatial_systems_data_utils.R` builds four explicit levels:

```
LEVEL 0   canonical source samples
LEVEL 1   ProteinGroupID x AnimalID x Hemisphere x SpatialUnit
          technical rows collapsed WITHIN a side; sides kept apart
LEVEL 2   ProteinGroupID x AnimalID x SpatialUnit
          equal-weight mean of the available sides
LEVEL 3   explicit left-only / right-only / bilateral matrices
```

The one-sided rule is **not invented here**. Both canonical animal-level paths
already agree, and this layer reuses their vocabulary:

- `R/statistics/wgcna_group_effects_utils.R` — `equal_weight_mean_available_LR_after_within_hemisphere_mean`, `one_sided_observed_no_imputation`
- `R/data_contracts/protigy_input_utils.R` — `single_observed_hemisphere_no_imputation`

A one-sided cell contributes its observed side unchanged. Nothing is imputed.

`SpatialUnit` is `Region x Layer` for neuropil and `Region` for soma and
microglia.

## 4. Why both hemispheres are retained

Averaging the sides at the point of import throws away the only repeated
measurement in the design. Keeping them allows two things that are otherwise
impossible: a direct measure of within-animal reproducibility, and a variance
decomposition that separates between-animal signal from within-animal
hemispheric variation.

### The defect this replaces

The empirical ROI pipeline returned a field named `hemisphere_mat`. It is keyed
`AnimalID x spatial unit` with **no side component** — Left and Right were
already averaged — and for soma and microglia it is literally the same object as
`region_mat`. A bilateral analysis built on it returns an exactly-zero
left-right difference *by construction* and looks like a real result. It is now
documented as a hemisphere-AVERAGED deprecated alias of `animal_spatial_mat`,
with a genuinely side-resolved `side_mat` beside it.

## 5. What bilateral validation means

It asks whether an effect **reproduces in magnitude and direction** on the
opposite side of the same brain.

It deliberately does **not** use "significant on both sides" as the primary
metric. Each one-sided fit uses half the samples, so requiring independent
significance in both hemispheres measures statistical power, not reproducibility.
That fraction is emitted only as a sensitivity descriptor.

For WGCNA modules two different questions are kept apart, because they can
disagree:

1. **Absolute value reproducibility** — does the left value equal the right value?
2. **Spatial pattern reproducibility** — within one animal, does the module's
   profile across spatial units have the same shape on both sides?

A module with a constant side offset scores poorly on (1) while being perfectly
preserved on (2). Poor agreement is also **not** automatically technical failure:
it may be real hemispheric asymmetry, so the classes describe the pattern rather
than issuing a quality verdict.

## 6. Why it is not independent biological replication

Both hemispheres come from the same animal. Agreement between sides is evidence
about measurement reliability and within-animal symmetry. It does not increase
`n`, it does not generalise across animals, and it must never be reported as
replication. The correct term throughout is **cross-hemisphere validation**.

## 7. External cell-type enrichment vs empirical compartment affinity

These are different kinds of evidence and are kept in separate streams.

| | empirical compartment affinity | external cell-type enrichment |
|---|---|---|
| source | this experiment's own ROIs | an outside reference dataset |
| asks | is this protein enriched in the microglia-enriched ROI relative to a neuronal compartment? | are this module's genes preferentially expressed in a reference cell type? |
| independence class | phenotype-independent context (same animals) | external annotation |
| must not be read as | a cell-proportion estimate | a statement that the tissue contains that cell type |

A microglia-enriched ROI is **not** purified microglia.

### The EWCE defects this fixes

1. **The background could be silently dropped.** A retry without `bg` substitutes
   the full reference transcriptome, while the caller still recorded the intended
   background size and cached under the intended key — a wrong result was
   indistinguishable from a correct one, including on disk. The background is now
   mandatory and failure is loud.
2. **One global BH family spanned both arms.** `q_global <- p.adjust(p)` was
   applied across Baseline and Differential rows together, so a phenotype-blind
   result's FDR depended on how many phenotype tests happened to be significant
   in the same run — layer D was determining layer A's inference. Families are
   now explicit:
   - `ewce_differential_<dataset>_<level>`
   - `ewce_module_annotation_<dataset>_<scope>_<level>`, BH across ModuleID x CellType,
     with `scope` in {`all`, `core_kME06`, `top25`}
3. **The "Baseline" arm was group-specific**, so it was never phenotype-blind.
   The new API takes no group argument at all; a module's gene set is its
   membership. Historical group-specific outputs are more accurately described as
   `group_specific_expression_annotation`.

`run_ewce_gene_set_annotation()` is the callable API. It does not reimplement
EWCE — the analysis script delegates to the same `ewce_bootstrap_once()`, so one
code path performs the test.

One practical note: the WGCNA membership tables store gene symbols uppercased,
which are not valid `org.Mm.eg.db` keys and match nothing in the mouse
specificity reference (EWCE then reports "Only 0 provided").
`ewce_to_mouse_symbols()` resolves the canonical mouse casing, recovering ~99.6%
of symbols. The module annotation also passes `output_species = "mouse"` to stay
in the reference's native gene space; the canonical analysis keeps EWCE's own
default so its historical behaviour is unchanged.

## 8. Prerequisites for later stages

These outputs must exist and validate before the atlas and network passes:

| output | needed by |
|---|---|
| `11_spatial_systems/data_contract/spatial_systems_hemisphere_inventory.csv` | everything side-resolved |
| `.../spatial_systems_aggregation_validation.csv` | proves level-2 reproduces the canonical aggregation |
| `.../spatial_systems_evidence_dependence.csv` | prevents evidence double-counting in the atlas |
| `bilateral/bilateral_spatial_identity_*.csv` | baseline spatial identity layer |
| `bilateral/bilateral_empirical_compartment_*.csv` | compartment affinity layer |
| `bilateral/WGCNA_module_bilateral_*.csv` | module reproducibility annotation |
| `precision/bilateral_*.csv` | whether bilateral averaging may be claimed to improve precision |
| `celltype_annotation/WGCNA_module_external_celltype_affinity_*.csv` | the external annotation column of the atlas |
| `spatial_systems_foundation_validation.csv` | gates the whole layer; a critical FAIL exits non-zero |

Explicitly **not** built in this pass: the module spatial/cell atlas, the
SUS-RES effect-at-identity-peak classification, per-animal spatial networks,
differential networks, the integrated workbook and any manuscript figure.

# Microglia ROI Interpretation

The `microglia` dataset represents microglia-enriched regions of interest and their local microenvironment. It must not be described as purified microglia, isolated microglial cells, or cell-intrinsic microglial regulation from these data alone.

## Rationale

The microglia ROI workflow enriches for microglia-associated sampling regions, but the measured protein abundance can include nearby neuropil, synaptic, vascular, extracellular, and other local microenvironment signal. That makes the dataset useful for region-level microglia-enriched ROI interpretation, not for purified-cell attribution.

The active dataset contract therefore sets:

```text
region = true
layer = false
celltype_roi = true
purified_celltype = false
```

Layer-resolved scripts must fail clearly for `microglia` unless future metadata explicitly extend this contract.

## Neuropil Reference Annotation

Neuropil is used as a reference annotation layer because microglia-enriched ROIs can share or neighbor neuropil/synaptic abundance patterns, and because dataset families may have been normalized and imputed separately. The reference annotation identifies overlap, biological context, and interpretation risk.

This is not:

- subtraction of neuropil abundance values
- decontamination of microglia ROI tables
- correction of log fold changes
- proof that a signal is or is not microglial

## Allowed Claims

- Region-level microglia-enriched ROI abundance signals.
- Local microenvironment-associated enrichment patterns.
- Microglia-marker-supported program enrichment when targeted signatures support the statement.
- Neuropil/synaptic reference overlap as an annotation, caveat, or source of biological context.
- Differential abundance within the sampled ROI contract.

## Unsafe Claims

- Purified microglial cell-intrinsic regulation from these data alone.
- Layer-level microglia conclusions under the current metadata contract.
- Neuropil-corrected or decontaminated microglia abundance unless a separate validated method is added.
- Causal mechanism without orthogonal validation.
- Cell-type origin of a protein abundance change without independent cell-type-resolved evidence.

## Recommended Wording

Use language such as "microglia-enriched ROI", "microglia-associated local microenvironment", "ROI-level differential abundance", and "neuropil reference overlap".

Avoid language such as "purified microglia", "microglia-intrinsic protein regulation", "decontaminated microglia signal", or "layer-specific microglia effect" unless new metadata and validation explicitly support those claims.

# Microglia ROI Interpretation

The microglia dataset represents microglia-enriched regions of interest and their local microenvironment. It must not be described as purified microglia.

Allowed claims:

- Region-level microglia-enriched ROI signals.
- Local microenvironment-associated enrichment patterns.
- Microglia-marker-supported program enrichment when supported by targeted signatures.
- Neuropil/synaptic reference overlap as an annotation or caveat.

Disallowed claims:

- Purified microglial cell-intrinsic regulation from these data alone.
- Layer-level microglia conclusions unless new metadata explicitly support them.
- Subtractive correction against neuropil intensities or logFC values.
- Causal mechanism without orthogonal validation.

Neuropil is used as a reference annotation layer because microglia ROIs can contain neuropil/synaptic signal and because dataset families may be separately normalized and imputed. The reference annotation identifies overlap and interpretation risk; it is not a subtraction or decontamination procedure.

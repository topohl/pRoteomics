# Metadata provenance

`TPE9_sample_metadata_males.xlsx` is the canonical sample-metadata workbook for the E9 proteomics analyses. Where structured metadata fields and the raw acquisition filename disagree, the structured metadata should be treated as authoritative unless a later documented correction supersedes it.

## A0003 CA2 neuron sample reassignment

For animal `A0003`, four CA2 neuron samples were deliberately reassigned after QC/UMAP review indicated that the `sp` and `sr` sample identities were likely switched. The raw `sample_id` strings retain the original acquisition filenames, but the structured metadata fields contain the corrected biological identities.

Corrected assignments:

| Sample | Hemisphere | Raw filename layer | Corrected metadata layer | Corrected compartment |
| --- | --- | --- | --- | --- |
| `S129` | Left | `sp` | `sr` | `neuron_neuropil` |
| `S130` | Right | `sp` | `sr` | `neuron_neuropil` |
| `S131` | Left | `sr` | `sp` | `neuron_soma` |
| `S132` | Right | `sr` | `sp` | `neuron_soma` |

The correction is internally consistent across `layer`, `group`, `group2`, and `celltype_layer`, and occurs as reciprocal `sp <-> sr` swaps in both hemispheres. Therefore downstream code should match samples by exact `sample_id` but use the structured metadata columns (`region`, `layer`, `celltype_layer`, `ReplicateGroup`, `AnimalID`) for biological annotation rather than reparsing anatomy from the filename.

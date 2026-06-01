# Datasets

The active dataset families are defined in `R/dataset_config.R` and mirrored in `pipeline.yml`.

| dataset | region | layer | celltype_roi | purified_celltype | interpretation |
|---|---:|---:|---:|---:|---|
| `neuron_neuropil` | true | true | false | false | Region/layer-resolved neuron neuropil proteomics. |
| `neuron_soma` | true | true | false | false | Region/layer-resolved neuronal soma-enriched proteomics where metadata support layer analyses. |
| `microglia` | true | false | true | false | Region-resolved microglia-enriched ROI/local microenvironment proteomics; not purified microglia. |

## Allowed Analyses

| dataset | allowed |
|---|---|
| `neuron_neuropil` | region summaries; layer summaries; differential abundance contrasts; enrichment/GSEA; WGCNA; module activity; layer-resolved spatial networks; network-behavior coupling when behavior inputs are present |
| `neuron_soma` | region summaries; layer summaries where metadata are complete; differential abundance contrasts; enrichment/GSEA; WGCNA; module activity; layer-resolved spatial networks; network-behavior coupling when behavior inputs are present |
| `microglia` | region summaries; differential abundance contrasts; enrichment/GSEA; WGCNA/module summaries when input contracts are met; targeted microglia ROI signature enrichment; neuropil reference annotation |

## Disallowed Analyses

| dataset | disallowed without new metadata or validation |
|---|---|
| `neuron_neuropil` | purified neuron cell-intrinsic claims from bulk ROI data alone |
| `neuron_soma` | layer-level claims when layer metadata are missing or ambiguous; purified neuron cell-intrinsic claims from enriched samples alone |
| `microglia` | layer-level analysis; spatial network scripts requiring region/layer units; purified microglia cell-intrinsic regulation; neuropil subtraction/decontamination claims |

## Examples

```r
source("R/dataset_config.R")

assert_dataset_capability("neuron_neuropil", "layer", analysis = "spatial network analysis")
assert_dataset_capability("microglia", "region", analysis = "regional ROI summary")
dataset_has_capability("microglia", "layer")
```

The final call returns `FALSE`. Active scripts that require layer-resolved inputs should call `assert_dataset_capability(dataset, "layer", ...)` and fail clearly for `microglia`.

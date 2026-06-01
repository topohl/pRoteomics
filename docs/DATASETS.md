# Datasets

The active dataset families are defined in `R/dataset_config.R` and mirrored in `pipeline.yml`.

| dataset | region | layer | celltype_roi | purified_celltype | interpretation |
|---|---:|---:|---:|---:|---|
| `neuron_neuropil` | true | true | false | false | Region/layer-resolved neuron neuropil proteomics. |
| `neuron_soma` | true | true | false | false | Region/layer-resolved neuronal soma-enriched proteomics where metadata support layer analyses. |
| `microglia` | true | false | true | false | Region-resolved microglia-enriched ROI/local microenvironment proteomics; not purified microglia. |

Use `dataset_has_capability(dataset, "layer")` or `assert_dataset_capability(dataset, "layer")` before running layer-level analyses. Microglia analyses must remain region-only unless independent metadata justify a layer-resolved interpretation.

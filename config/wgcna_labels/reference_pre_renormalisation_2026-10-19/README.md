Reviewed WGCNA label registries for the network built before the 2026-10-19 re-normalisation (pRoteomics b478a13).
Kept as the versioned reference for the old-to-new lineage (tools/map_wgcna_module_lineage.R). They are not read by the pipeline:
module IDs of the rebuilt network need not denote the same modules. Labels for the new network are promoted into
config/wgcna_labels/<dataset>.csv by hand after review. gsea_wgcna_theme_module_mapping.csv is that network's theme-to-module
mapping (same commit), the lineage tool's old mapping; config/gsea_wgcna_theme_module_mapping.csv follows the current network.

# Allen Brain Cell Atlas, whole mouse brain 10x: mean expression per source cell type

Used by `analysis/spatial_validation/map_afferent_sources.R` (Figure 2k). It groups proteins by the cell type whose mRNA they are enriched in: entorhinal or other afferent neurons, CA3, granule, CA1 or CA2 cells, interneurons, glia, or a control class of non-projecting neurons. The thresholds are in that script's header. mRNA from other animals is a proxy for where a protein is made, not proof.

- **Source.** Yao et al. 2023, whole-mouse-brain taxonomy CCN20230722 (AIT21.0).
  - Expression is the official cluster-mean table (mean log2(CPM+1) per cluster) from the MapMyCells file `precomputed_stats_ABC_revision_230821.h5`.
  - Taxonomy and cell counts come from metadata/WMB-taxonomy/20231215.
  - Every file, URL, size and MD5 is listed in `download_manifest.csv`. The raw files are not kept here; they total 2.6 GB.
- **`abc_wmb_source_type_mean_log2cpm.csv`**: genes x 35 source types, the cell-weighted mean of the cluster means. SHA-256 `778e6a2dc82f55cbbbdd74060e56b579600e874b4be5add9227d4fc3e4e7dd88`.
- **`source_type_map.csv`**: the ABC clusters, supertypes or subclasses behind each source type, with the evidence and a confidence rating. SHA-256 `7c21a651c7be90420862c0f738931b1d80c2eb460170af068e85656716b9c4d7`.
- **`sanity_check_markers.csv`**: known markers checked in their expected types.
- **`build_info.txt`**: how the means were computed.
- **`scripts/`**: download, means and map, rerunnable in that order from `scripts/`: the raw files go to `raw/` here (or `ABC_WMB_RAW`), the tables to this folder. The released tables were built in a staging folder (paths in `download_manifest.csv`) with the same code. The scripts, this README, the manifests, the map and the small tables are tracked in git; the means table and the raw files are not, and `map_afferent_sources.R` refuses a means table whose SHA-256 differs from the one above.

Caveats:
- The cluster means are effectively 10Xv3-only.
- The values are means of log2(CPM+1), not the log of the mean CPM, so compare one gene across types.
- Small types: locus coeruleus 22 cells, CA2 143 cells, mossy cells 318 cells.
- The supramammillary types include posterior hypothalamus.
- 40 gene symbols are duplicated; the first is used.

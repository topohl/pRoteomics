# Published hippocampal stress signatures

Used by `analysis/spatial_validation/compare_external_stress_studies.R` (Extended Data). It asks whether published stress datasets show the same SUS / RES / CON differences as our spatial proteome. Each signature is one statistic per mouse gene and contrast, rebuilt from the deposited data. Positive values mean higher in the first group: `SUS_vs_CON`, `RES_vs_CON`, `SUS_vs_RES`, `stressed_vs_CON`, where stressed = (SUS + RES)/2, or the single stressed group.

## Layout

- **`scripts/download_stress_signatures.R`** fetches the 76 raw files (434 MB) into `raw/<ACCESSION>/`, or into the folder named by `STRESS_SIG_RAW`.
  - Files already present with the manifest MD5 are skipped.
  - `--check` checks presence and MD5 only, without network access.
  - Two inputs are GEO record pages (`GSE139684_gsm.txt`, `GPL13912.txt`), so their MD5 changes when GEO edits the record. The script then keeps the new file, reports it and exits with status 1; rebuild and compare.
- **`scripts/build_stress_signatures.R`** builds `raw/` into `signatures/<id>.csv` (gene, then one column per contrast) and `signature_info.csv`. Pass ids to build only some signatures. Runtime is about 3 min.
- **`signature_info.csv`** has one row per signature: source, design, n per group, pooling, contrasts, `match_units` (our units it is compared with, or `all`), `set`, `corrections`, notes, and the SHA-256 of `signatures/<id>.csv`. The analysis refuses a signature whose checksum differs.
- **`download_manifest.csv`** lists URL, size and MD5 for every raw file.
- **Not tracked:** `raw/` and `signatures/`. The 2026-10-29 raw downloads are kept in `C:/Users/topohl/Documents/exp9_fig2_rework_2026-10-27/stress_signatures_raw` (set `STRESS_SIG_RAW` to that folder).

## Methods shared by all signatures

- **Count data:** limma-voom. Genes are kept if CPM > 1 in at least the smallest group. The voom library size is a median-of-ratios size factor times the geometric-mean depth. eBayes is robust.
- **Log-scale data:** limma-trend with robust eBayes.
- **Gene mapping:** the key is the official MGI symbol (org.Mm.eg.db). Non-official symbols map through a unique alias. Rat genes (Henningsen) map rat to human to mouse with babelgene, keeping 1:1 mappings only.
- **Duplicate symbols:** when several rows map to one symbol, the row with the highest mean expression is kept.
- **Packages:** limma, data.table, readxl, babelgene, AnnotationDbi, org.Mm.eg.db, digest. Henningsen also needs org.Rn.eg.db and org.Hs.eg.db.
- **GSE81672 workbook:** this `.xls` cannot be opened by readxl. The build includes a small BIFF8 reader whose output equals the file's own CSV export cell for cell.
- **Checks:** every signature was built and then checked by an independent agent on 2026-10-29. The checks covered group mapping from the deposited annotation, contrast signs against the published results, and gene mapping. The rebuild reproduces the checked versions exactly.

## Signatures

| id | source | tissue, model, timing | n (CON / RES / SUS) | contrasts | set |
|---|---|---|---|---|---|
| Schmidt2010_dCA1, Schmidt2010_dDG | Schmidt et al. 2010, PMID 21159965, GSE11211 | dorsal CA1 / DG laser microdissection; adolescent social instability; 5 weeks; SUS/RES by basal corticosterone | one pool per group (4 technical arrays) | SUS_vs_RES (log2 ratio) | primary |
| Lehmann2018_microglia | Lehmann et al. 2018, PMID 30050134, GSE107941 | sorted whole-brain microglia; social defeat | 6 / 9 / 6 | SUS_vs_CON | primary |
| Henningsen2012_CA | Henningsen et al. 2012, PMID 22311638, PRD000546 | rat CA1 + CA3 proteome (iTRAQ); chronic mild stress | 3 pools of 3 per group | all four | primary |
| Bagot2016_vHIP_48h, Bagot2016_vHIP_28d | Bagot et al. 2016, PMID 27181059, GSE72343 (counts: Flati et al. 2020) | ventral hippocampus; social defeat; 48 h, 28 d | 3/3/3, 4/4/4 pools | all four | primary |
| Bagot2017_vHIP | Bagot et al. 2017, PMID 27569543, GSE81672 | ventral hippocampus; social defeat, saline arm | 4 / 4 / 3 | all four | primary |
| Laine2018_vHPC | Laine et al. 2018, PMID 30073192, GSE109315 | ventral hippocampus; social defeat; 6-8 d | 6 / 8 / 3 | all four | primary |
| Ayash2023_vHC | Ayash et al. 2023, PMID 37068238, GSE161726 | ventral hippocampus; social defeat; ~1 d | 6 / 7 / 12 | all four | primary |
| Butto2024_vHip_nuclei | Butto et al. 2024, PMID 38649752, GSE240573 | activated (Arc-tagged) nuclei, ventral hippocampus; social defeat | 6 / 7 / 6 pools of 2 | all four | primary |
| Bhatti2022_PV_TRAP | Bhatti et al. 2022, PMID 35813065, GSE184027 | parvalbumin-interneuron translatome, whole hippocampus; social defeat | 4 / 4 / 3 | all four | primary |
| Lyons2025_CA1, Lyons2025_CA3 | Lyons et al. 2025, PMID 39528642, GSE278620 | GeoMx CA1 / CA3; chronic subordination stress | 4 vs 6 mice | stressed_vs_CON | primary |
| Weber2019_microglia | Weber et al. 2019, PMID 30527629, GSE124169 | sorted whole-brain microglia; repeated social defeat; 24 d | 6 vs 7 | stressed_vs_CON | primary |
| Herve2017_DG | Herve et al. 2017, PMID 28848385, GSE84183 | dentate gyrus; unpredictable chronic mild stress | 8 vs 8 | stressed_vs_CON | primary |
| Sturman2021_vHC | Sturman et al. 2021, PMID 34527792, GSE172451 | ventral hippocampus; adolescent social instability; adult | 24 vs 24 | stressed_vs_CON | primary |
| Russo_astro_TRAP | unpublished deposit GSE139684 | astrocyte translatome, whole hippocampus; social defeat; 2 d | 3 / 3 / 5 pools of 2 | all four | exploratory (contamination dominates) |
| SENS_Schmidt2010_dCA1_asGEO | as Schmidt2010_dCA1 | GEO CA1 labels as deposited | | SUS_vs_RES | exploratory |

## Deposit corrections

- **GSE11211:** the CA1 labels on GEO are inverted relative to the published CA1 results. With the deposited labels, Gria1, Gria2 and the 1.9-fold Slc6a15 decrease all have the wrong sign. Schmidt2010_dCA1 follows the paper; the deposited orientation is kept as SENS_Schmidt2010_dCA1_asGEO. The DG labels agree with the paper.
- **GSE72343 (28 d):** GSM1860619/20 are titled and filed as susceptible, but their phenotype characteristic says resilient. The titles are used, giving the paper's n = 4 per group.
- **GSE84183:** the data columns of the series matrix are permuted relative to the sample records. The 16 per-sample raw files are used instead.
- **GSE107941:** all resilient arrays were scanned on a different day, so only SUS_vs_CON is provided.

## Caveats

- Most datasets are ventral or whole hippocampus (the microglia sets are whole brain); ours are dorsal subregions.
- Models, sampling times and definitions of susceptibility differ between studies.
- Several datasets use pooled samples (see `pooled_samples`).
- **Laine2018:** sequencing set is partly confounded with group.
- **Ayash2023:** SUS combines the two avoider groups (fixed before analysis).
- **Herve2017:** group is confounded with array position.
- Only the main statistic of each signature is packaged.

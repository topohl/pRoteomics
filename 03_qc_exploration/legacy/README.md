# Legacy QC Scripts

These scripts are retained for provenance only and are not part of the canonical
QC/exploration run order.

- `06_pcaPlot_Neha.r`: collaborator-specific PCA workflow for an external Neha
  input file. Use `../05_pca_confounding_qc.r` for canonical dataset-aware PCA.
- `08_boxplotBonanza.r`: old single-protein `TKNK_MOUSE` plotting script. It
  depended on hard-coded collaborator inputs and an undefined `available_groups`
  loop variable, so it was archived instead of being left in the active QC run.


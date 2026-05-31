#!/usr/bin/env bash
set -euo pipefail

# Non-destructive repository restructure for pRoteomics.
# This script uses git mv, preserving file history as much as Git can infer.
# Run from the repository root on a clean working tree:
#   bash tools/restructure_pipeline_folders.sh
# Then inspect with:
#   git status
#   git diff --stat

if [ ! -d .git ]; then
  echo "ERROR: run this from the repository root." >&2
  exit 1
fi

if [ -n "$(git status --porcelain)" ]; then
  echo "ERROR: working tree is not clean. Commit or stash changes first." >&2
  exit 1
fi

mkdir -p \
  01_preprocessing \
  02_id_mapping \
  03_qc_exploration \
  04_differential_expression_enrichment \
  05_celltype_enrichment_EWCE \
  06_modules_WGCNA \
  07_spatial_networks \
  08_behavior_physio_coupling \
  90_testing \
  99_deprecated \
  tools

move_if_exists() {
  local src="$1"
  local dst="$2"
  if [ -e "$src" ]; then
    mkdir -p "$(dirname "$dst")"
    git mv "$src" "$dst"
    echo "moved: $src -> $dst"
  else
    echo "skip missing: $src"
  fi
}

# 01 preprocessing / metadata preparation
move_if_exists "Formatting/gct_extractR.r" "01_preprocessing/01_gct_extractR.r"
move_if_exists "Formatting/excel_convert.r" "01_preprocessing/02_excel_convert.r"
move_if_exists "Formatting/format_metadata.r" "01_preprocessing/03_format_metadata.r"
move_if_exists "Analysis/metadata_create.r" "01_preprocessing/04_metadata_create.r"
move_if_exists "Formatting/impute.r" "01_preprocessing/05_impute.r"
move_if_exists "Formatting/merged_metadata_module_score.r" "01_preprocessing/06_merged_metadata_module_score.r"

# 02 ID mapping
move_if_exists "Mapping/MapThatProt.r" "02_id_mapping/01_MapThatProt.r"
move_if_exists "Mapping/MapThatProt_batch.r" "02_id_mapping/02_MapThatProt_batch.r"
move_if_exists "Mapping/MapThatProt_batch_WGCNAtoClusterProfiler.r" "02_id_mapping/03_MapThatProt_batch_WGCNAtoClusterProfiler.r"

# 03 QC / exploratory structure
move_if_exists "Analysis/qc_protein_peptide_plot.r" "03_qc_exploration/01_sample_qc_quicksearch.r"
move_if_exists "Analysis/rank_abundance_plot_E9.r" "03_qc_exploration/04_marker_rank_abundance_qc.r"
move_if_exists "Analysis/pcaPlot.r" "03_qc_exploration/05_pca_confounding_qc.r"
move_if_exists "Analysis/pcaPlot_v2.r" "03_qc_exploration/04_pcaPlot_v2.r"
move_if_exists "Analysis/pcaPlot_v3.r" "03_qc_exploration/05_pcaPlot_v3.r"
move_if_exists "Analysis/pcaPlot_Neha.r" "03_qc_exploration/06_pcaPlot_Neha.r"
move_if_exists "Analysis/varPart.r" "03_qc_exploration/06_variance_partitioning.r"
move_if_exists "Analysis/boxplotBonanza.r" "03_qc_exploration/08_boxplotBonanza.r"

# 04 differential expression / enrichment
move_if_exists "Analysis/clusterProfiler.r" "04_differential_expression_enrichment/01_clusterProfiler.r"
move_if_exists "Analysis/compareGO.r" "04_differential_expression_enrichment/02_compareGO.r"
move_if_exists "Analysis/compare_pathways.r" "04_differential_expression_enrichment/04_compare_pathways.r"
move_if_exists "Analysis/compare_sig_expr.r" "04_differential_expression_enrichment/05_compare_sig_expr.r"
move_if_exists "Analysis/control_strata_enrichment_figures.r" "04_differential_expression_enrichment/07_control_strata_enrichment_figures.r"

# 05 cell-type enrichment / EWCE
move_if_exists "Analysis/EWCE_E9.r" "05_celltype_enrichment_EWCE/01_EWCE_E9.r"
move_if_exists "Analysis/EWCE.r" "05_celltype_enrichment_EWCE/90_EWCE_legacy.r"

# 06 WGCNA / module-level analyses
move_if_exists "Analysis/WGCNA.r" "06_modules_WGCNA/01_WGCNA.r"
move_if_exists "Analysis/WGCNAtraitpreservation.r" "06_modules_WGCNA/02_WGCNAtraitpreservation.r"
move_if_exists "Analysis/module_spatial_networks.r" "06_modules_WGCNA/03_module_spatial_networks.r"
move_if_exists "testing/overlap_modules.r" "06_modules_WGCNA/04_overlap_modules.r"
move_if_exists "testing/module_score v.0.0.1.r" "06_modules_WGCNA/90_module_score_v0.0.1.r"
move_if_exists "testing/module_score v.0.0.2.r" "06_modules_WGCNA/91_module_score_v0.0.2.r"

# 07 spatial network analyses
move_if_exists "Analysis/network_spatial_relations.r" "07_spatial_networks/01_network_spatial_relations.r"
move_if_exists "Analysis/differential_networks.r" "07_spatial_networks/02_differential_networks.r"
move_if_exists "Analysis/bootstrap_network_stability.r" "07_spatial_networks/03_bootstrap_network_stability.r"
move_if_exists "Analysis/bootstrap_differential_network_stability.r" "07_spatial_networks/04_bootstrap_differential_network_stability.r"
move_if_exists "Analysis/bootstrap_differential_network_figures.r" "07_spatial_networks/05_bootstrap_differential_network_figures.r"
move_if_exists "Analysis/chord_diagram.r" "07_spatial_networks/06_chord_diagram.r"

# 08 behavior / physiology coupling
move_if_exists "Analysis/correlate_proteomics_with_behavior.r" "08_behavior_physio_coupling/01_correlate_proteomics_with_behavior.r"
move_if_exists "Analysis/network_behavior_coupling.r" "08_behavior_physio_coupling/02_network_behavior_coupling.r"

# testing and deprecated archives
move_if_exists "testing/WGCNAtraitpreservation.r" "90_testing/WGCNAtraitpreservation.r"
move_if_exists "testing/WGCNAtraitpreservation_2.r" "90_testing/WGCNAtraitpreservation_2.r"
move_if_exists "testing/chord_diagram.r" "90_testing/chord_diagram.r"
move_if_exists "testing/clusterProfiler v2.r" "90_testing/clusterProfiler_v2.r"
move_if_exists "testing/clusterProfiler v3 test 2.r" "90_testing/clusterProfiler_v3_test_2.r"
move_if_exists "testing/clusterProfiler v3 test.r" "90_testing/clusterProfiler_v3_test.r"
move_if_exists "testing/clusterProfiler v3.r" "90_testing/clusterProfiler_v3.r"
move_if_exists "testing/clusterProfilerBatch.r" "90_testing/clusterProfilerBatch.r"
move_if_exists "testing/clusterProfiler_test.r" "90_testing/clusterProfiler_test.r"
move_if_exists "testing/clusterProfiler_v0.0.1.r" "90_testing/clusterProfiler_v0.0.1.r"
move_if_exists "testing/compareGO reduced.r" "90_testing/compareGO_reduced.r"
move_if_exists "testing/compareGOBatch.r" "90_testing/compareGOBatch.r"
move_if_exists "testing/compare_GO_overlap v.0.0.1.r" "90_testing/compare_GO_overlap_v0.0.1.r"

move_if_exists "deprecated/clusterProfiler newest Sep16.r" "99_deprecated/clusterProfiler_newest_Sep16.r"
move_if_exists "deprecated/clusterProfiler newest Sep16_CONSOLIDATION + SUMMARY.r" "99_deprecated/clusterProfiler_newest_Sep16_CONSOLIDATION_SUMMARY.r"
move_if_exists "deprecated/clusterProfiler_jan2026.r" "99_deprecated/clusterProfiler_jan2026.r"
move_if_exists "deprecated/compareGO v2.r" "99_deprecated/compareGO_v2.r"

# Remove empty old buckets if they are empty.
rmdir Analysis Formatting Mapping testing deprecated 2>/dev/null || true

echo "Done. Inspect changes with: git status && git diff --stat"
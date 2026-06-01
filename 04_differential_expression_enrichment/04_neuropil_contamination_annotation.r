#!/usr/bin/env Rscript

# Neuropil contamination annotation workflow
# Purpose:
# Annotate microglia enrichment results using matched neuropil analyses.
# This script intentionally performs annotation and sensitivity assessment,
# not direct subtraction of neuropil effects.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(openxlsx)
})

message('Neuropil annotation workflow')
message('TODO: connect to clusterProfiler manifests')

classify_term <- function(microglia_fraction,
                          neuropil_fraction,
                          overlap_fraction) {

  if (!is.na(microglia_fraction) &&
      microglia_fraction >= 0.25 &&
      overlap_fraction < 0.20) {
    return('microglia_robust')
  }

  if (!is.na(neuropil_fraction) &&
      neuropil_fraction >= 0.25 &&
      overlap_fraction >= 0.50) {
    return('neuropil_sensitive')
  }

  if (overlap_fraction >= 0.20) {
    return('mixed_microenvironment')
  }

  'ambiguous'
}

marker_sets <- list(
  microglia = c('P2RY12','TMEM119','CX3CR1','CSF1R','TYROBP','HEXB','C1QA'),
  neuropil = c('SYN1','SYP','SNAP25','DLG4','CAMK2A','MAP2')
)

message('Marker sets loaded')
message('Future outputs:')
message('- GO term overlap annotation')
message('- compareGO overlap annotation')
message('- marker composition scoring')
message('- microglia_robust / mixed_microenvironment / neuropil_sensitive classification')

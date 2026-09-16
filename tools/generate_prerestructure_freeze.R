#!/usr/bin/env Rscript

# Pre-restructure publication and science freeze manifest.
#
# This exists to be an EQUIVALENCE ORACLE for a purely structural repository
# migration. The migration that follows moves files; it must not change science.
# The only way to prove that afterwards is to know, beforehand, the exact content
# of everything whose meaning must survive the move.
#
# What goes in is therefore chosen by one rule: if this object's CONTENT changed
# during a migration, would a scientific claim in the manuscript change with it?
# Rendered figures are included as artefacts, the tables behind them as content,
# the contracts that bind them as configuration, and the provenance tables that
# tie claims to evidence. Paths WILL change in the migration - that is the point
# of the migration - so the manifest records content hashes and a stable
# publication identity alongside the current path, never the path alone.
#
# Deliberately NOT a re-run of anything. It reads and hashes.

source(file.path("R", "paths.R"))
source(repo_path("R", "null_coalescing.R"))

OUT_CSV <- repo_path("manuscript", "prerestructure_freeze_manifest.csv")
OUT_MD <- repo_path("docs", "PRERESTRUCTURE_FREEZE.md")

sha <- function(p) if (file.exists(p) && !dir.exists(p)) {
  digest::digest(file = p, algo = "sha256")
} else NA_character_

rows <- list()
add <- function(class, publication_id, rel, role, note = "") {
  abs <- repo_path(rel)
  rows[[length(rows) + 1L]] <<- data.frame(
    object_class = class,
    publication_id = publication_id,
    repository_relative_path = rel,
    role = role,
    exists = file.exists(abs),
    bytes = if (file.exists(abs) && !dir.exists(abs)) file.size(abs) else NA_real_,
    sha256 = sha(abs),
    note = note,
    stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------- 1. contracts
for (f in c("figures/figure_contract.yml",
            "figures/figure_final_truth_v9_contract.yml",
            "config/manuscript_spatial_order.yml",
            "config/manuscript_palette.yml",
            "config/clusterProfiler_config.yml",
            "config/output_namespaces.yml",
            "pipeline.yml")) {
  add("configuration_contract", NA_character_, f, "binds panels to evidence")
}

# ------------------------------------------- 2. manuscript claim and provenance
for (f in list.files(repo_path("manuscript"), pattern = "[.]csv$")) {
  add("manuscript_provenance", NA_character_, file.path("manuscript", f),
      "claim / statement / audit provenance")
}
add("manuscript_text", NA_character_, "manuscript/manuscript_draft.md", "the manuscript")
for (f in c("manuscript/figure1_legend.md", "manuscript/extended_data_09_legend.md")) {
  add("manuscript_text", NA_character_, f, "figure legend")
}

# --------------------------------------------- 3. frozen behavioural bridge
br <- "manuscript/figure1_bridge_mmmsociability"
for (f in list.files(repo_path(br), recursive = TRUE)) {
  add("frozen_upstream_import", NA_character_, file.path(br, f),
      "byte-exact import from topohl/MMMSociability",
      "line-ending normalisation is disabled for this tree in .gitattributes")
}

# ------------------------------- 4. canonical figure artefacts and source data
contract <- yaml::read_yaml(repo_path("figures", "figure_contract.yml"))
for (key in names(contract$figures)) {
  fig <- contract$figures[[key]]
  if (!isTRUE(fig$is_numbered_manuscript_figure)) next
  pid <- as.character(fig$canonical_publication_id %||%
    if (!is.null(fig$extended_data_number))
      sprintf("extended_data_%02d", as.integer(fig$extended_data_number))
    else sprintf("figure_%02d", as.integer(key)))
  for (panel in fig$panels) {
    if (!is.null(panel$figure_source))
      add("canonical_figure_panel", pid, as.character(panel$figure_source),
          paste("panel", panel$id))
    if (!is.null(panel$primary_source))
      add("publication_source_data", pid, as.character(panel$primary_source),
          paste("source data behind panel", panel$id))
    for (d in as.character(unlist(panel$input_dependencies %||% character()))) {
      add("canonical_analysis_table", pid, d, paste("input to panel", panel$id))
    }
  }
  stub <- sub("^figure_0", "figure_", pid)
  for (ext in c("svg", "png", "pdf")) {
    add("canonical_assembled_figure", pid,
        file.path("results", "figures", "manuscript", pid, "assembled",
                  paste0(pid, ".", ext)),
        "assembled artefact")
  }
}

# ------------------------------------------------------ 5. protected state
for (f in c("results/tables/11_spatial_systems/atlas/protein_spatial_cell_affinity.csv",
            "results/tables/10_biological_integration/gsea_wgcna_concordance/global/ontology_aware_gsea_theme_assignments_all_contrasts.csv",
            "results/tables/11_spatial_systems/bilateral/bilateral_spatial_identity_summary.csv",
            "docs/publication_freeze_manifest.yml",
            "docs/ATLAS_PROGRAM_SELECTION_AND_NAMING_RULES.md",
            "docs/MANUSCRIPT_STATISTICAL_CONTRACT.md")) {
  add("protected_scientific_state", NA_character_, f, "must be invariant through migration")
}

# ------------------------------------------------------------- 6. guard tests
for (f in c("test-figure-promotion-v9.R", "test-figure-generation-adjudication.R",
            "test-extended-data-behaviour.R", "test-figure1-behaviour-bridge.R",
            "test-figure-01-renderer.R", "test-manuscript-figure-entrypoints.R",
            "test-publication-freeze-manifest.R", "test-output-namespace-contract.R",
            "test-pipeline-registry.R", "test-candidate-figure-layer.R")) {
  add("guard_test", NA_character_, file.path("tests", "testthat", f),
      "encodes a contract the migration must not break")
}

m <- unique(do.call(rbind, rows))
m <- m[order(m$object_class, m$publication_id, m$repository_relative_path), ]

con <- file(OUT_CSV, open = "wb")
write.csv(m, con, row.names = FALSE, na = "", eol = "\n")
close(con)

present <- sum(m$exists)
cat("pre-restructure freeze manifest written:", OUT_CSV, "\n")
cat("  objects:", nrow(m), " present:", present, " missing:", nrow(m) - present, "\n")
print(table(m$object_class))

md <- c(
  "# Pre-restructure publication and science freeze",
  "",
  "Generated by `tools/generate_prerestructure_freeze.R`. The machine-readable",
  "manifest is `manuscript/prerestructure_freeze_manifest.csv`.",
  "",
  "## What this is for",
  "",
  "A purely structural repository migration moves files. This manifest is the",
  "oracle that proves it moved nothing else. Every object listed here is one",
  "whose *content* must be identical before and after: if it changed, a claim in",
  "the manuscript would change with it.",
  "",
  "Paths are expected to change - that is what the migration is. The manifest",
  "therefore records a content hash and, where one exists, a stable publication",
  "identity, so an object can be re-found after it moves.",
  "",
  "## How to use it after the migration",
  "",
  "For every row, locate the object at its new path and compare `sha256`.",
  "A differing hash is a migration defect until proven otherwise, with two",
  "expected exceptions that must be argued explicitly rather than assumed:",
  "",
  "- a file whose content legitimately encodes its own path;",
  "- a file under `manuscript/figure1_bridge_mmmsociability/`, where",
  "  `.gitattributes` disables line-ending normalisation precisely so that the",
  "  bytes stay comparable. If those hashes change, normalisation was applied and",
  "  the byte-exactness of the upstream import has been lost.",
  "",
  sprintf("## Contents at generation: %d objects, %d present", nrow(m), present),
  "",
  "| class | objects |",
  "|---|---|",
  paste0("| ", names(table(m$object_class)), " | ", as.integer(table(m$object_class)), " |"),
  "",
  "Objects recorded as absent are declared dependencies that are not materialised",
  "in this checkout - typically regenerable `results/` artefacts. They are listed",
  "because their identity matters to the migration even when their bytes are not",
  "currently present.")
writeLines(md, OUT_MD, useBytes = TRUE)
cat("  doc written:", OUT_MD, "\n")

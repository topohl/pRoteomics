#!/usr/bin/env Rscript
# ================================================================
# Script: 06_modules_WGCNA/08_module_complex_architecture.r
# Stage: modules_downstream
# Scope: dataset_specific
# Consumes: existing WGCNA module definitions and biological annotations.
# Produces: protein-complex/organelle architecture evidence for modules.
# ================================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "wgcna_downstream_utils.R"))

SCRIPT_ID <- "06_modules_WGCNA/08_module_complex_architecture.r"
run <- integration_cli(allow_all = TRUE)

complex_sets <- list(
  mitochondrial_oxphos = c("Ndufs1", "Ndufs2", "Ndufv1", "Ndufa9", "Ndufb8", "Sdha", "Sdhb", "Uqcrc1", "Uqcrc2", "Cyc1", "Cox4i1", "Cox5a", "Atp5f1a", "Atp5f1b"),
  ribosome_translation = c("Rpl3", "Rpl4", "Rpl5", "Rpl7", "Rpl10", "Rps3", "Rps6", "Rps8", "Eef1a1", "Eef2"),
  synaptic_vesicle_psd = c("Snap25", "Syn1", "Syn2", "Syp", "Syt1", "Vamp2", "Dlg4", "Shank1", "Homer1"),
  lysosome_phagosome = c("Ctsb", "Ctsd", "Ctsz", "Laptm5", "Lamp1", "Tyrobp", "Trem2", "Apoe"),
  ecm_basement_membrane = c("Col4a1", "Col4a2", "Lama2", "Lamb1", "Lamb2", "Nid1", "Hspg2", "Itgb1")
)

make_dataset <- function(ds) {
  paths <- create_module_dirs("06_modules_WGCNA", file.path("module_complex_architecture", ds))
  files <- resolve_wgcna_files(ds)
  inputs <- list(definitions = files$definitions, module_annotation = path_results("tables", "06_modules_WGCNA", "module_annotation", ds, "WGCNA_module_biological_annotation.csv"))
  if (run$dry_run) {
    dry_run_inputs(paste(SCRIPT_ID, ds), inputs)
    return(NULL)
  }
  defs_loaded <- read_csv_optional(inputs$definitions, ds, "wgcna", "module_definitions", required = FALSE)
  ann_loaded <- read_csv_optional(inputs$module_annotation, ds, "wgcna", "module_annotation", required = FALSE)
  status <- rbind(defs_loaded$status, ann_loaded$status)
  defs <- defs_loaded$data
  if (is.null(defs) || !nrow(defs)) {
    out <- availability_evidence(ds, "module_complex_architecture", inputs$definitions, "WGCNA definitions unavailable.")
  } else {
    mod_col <- first_col(defs, c("ModuleID", "module_id", "ModuleColor"))
    gene_col <- first_col(defs, c("GeneSymbol", "gene_symbol", "Gene"))
    if (is.na(mod_col) || is.na(gene_col)) {
      out <- availability_evidence(ds, "module_complex_architecture", inputs$definitions, "WGCNA definitions lack module/gene columns.")
    } else {
      defs$gene_token <- normalize_token(defs[[gene_col]])
      rows <- list()
      for (module in unique(as.character(defs[[mod_col]]))) {
        module_genes <- unique(defs$gene_token[as.character(defs[[mod_col]]) == module])
        for (set_name in names(complex_sets)) {
          set_tokens <- normalize_token(complex_sets[[set_name]])
          overlap <- intersect(module_genes, set_tokens)
          rows[[length(rows) + 1L]] <- data.frame(
            dataset = ds,
            evidence_domain = "module_complex_architecture",
            evidence_id = paste(ds, module, set_name, sep = "::"),
            program_label = set_name,
            entity_type = "module",
            entity_id = module,
            support_count = length(overlap),
            evidence_status = ifelse(length(overlap) >= 3, "complex_or_organelle_supported", "not_supported"),
            interpretation_note = paste0("Matched proteins: ", paste(overlap, collapse = ";")),
            source_file = inputs$definitions,
            qc_flag = ifelse(length(overlap) >= 3, "PASS", "INFO"),
            stringsAsFactors = FALSE
          )
        }
      }
      out <- standardize_evidence(do.call(rbind, rows))
      if (!is.null(ann_loaded$data) && nrow(ann_loaded$data)) {
        ann <- ann_loaded$data
        ann_key <- first_col(ann, c("ModuleID", "module_id"))
        cls_col <- first_col(ann, c("microenvironment_class", "microenvironment_label"))
        if (!is.na(ann_key) && !is.na(cls_col)) {
          out$interpretation_note <- paste(out$interpretation_note, ann[[cls_col]][match(out$entity_id, ann[[ann_key]])], sep = "; microenvironment=")
        }
      }
    }
  }
  write_integration_table(out, paths, "module_complex_architecture.csv")
  write_csv_safe(status, file.path(paths$reports, "input_status.csv"))
  write_csv_safe(status, file.path(paths$source_data, "module_complex_architecture_input_status.csv"))
  write_integration_manifest(paths, inputs, list(tables = paths$tables, source_data = paths$source_data), list(dataset = ds), "Module protein-complex/organelle architecture from existing WGCNA definitions; no WGCNA recomputation.")
  out
}

res <- lapply(integration_datasets(run$dataset), make_dataset)
if (run$dry_run) quit(status = 0, save = "no")
message("Module complex architecture complete for: ", paste(integration_datasets(run$dataset), collapse = ", "))

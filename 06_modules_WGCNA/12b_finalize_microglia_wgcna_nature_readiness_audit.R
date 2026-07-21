#!/usr/bin/env Rscript

# Finalizes only the additive audit manifest after portable HTML packaging.

options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default) {
  hit <- match(flag, args)
  if (!is.na(hit) && hit < length(args)) args[[hit + 1L]] else default
}
repo_root <- normalizePath(arg_value("--repo-root", getwd()), winslash = "/", mustWork = TRUE)
output_dir <- file.path(
  repo_root,
  arg_value("--output-dir", "results/reviewer_audit/microglia_wgcna_nature_readiness")
)
if (!dir.exists(output_dir)) stop("Audit output directory does not exist: ", output_dir)

atomic_write_csv <- function(x, path) {
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  utils::write.csv(x, tmp, row.names = FALSE, na = "")
  if (!file.exists(tmp) || file.info(tmp)$size <= 0) stop("Incomplete write: ", path)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Could not publish: ", path)
}

required <- c(
  "module_eigengene_variance_partition.csv",
  "supermodule_eigengene_variance_partition.csv",
  "variance_partition_heatmap.svg",
  "variance_partition_heatmap.pdf",
  "model_estimability_audit.csv",
  "module_sensitivity_robustness.csv",
  "module_leave_one_animal_out.csv",
  "higher_order_block_robustness.csv",
  "module_biological_context_validation.csv",
  "module_to_gene_bridge.csv",
  "transcriptomic_matching_contract_audit.csv",
  "validation_table.csv",
  "protected_output_hash_audit.csv",
  "WGCNA_nature_readiness_report_artifact.json",
  "WGCNA_nature_readiness_report.html",
  "WGCNA_nature_readiness_report.md",
  "session_info.txt",
  "input_hashes.csv",
  "README.md"
)
required_paths <- file.path(output_dir, required)

robustness <- utils::read.csv(file.path(output_dir, "module_sensitivity_robustness.csv"), check.names = FALSE)
loao <- utils::read.csv(file.path(output_dir, "module_leave_one_animal_out.csv"), check.names = FALSE)
blocks <- utils::read.csv(file.path(output_dir, "higher_order_block_robustness.csv"), check.names = FALSE)
context <- utils::read.csv(file.path(output_dir, "module_biological_context_validation.csv"), check.names = FALSE)
protected <- utils::read.csv(file.path(output_dir, "protected_output_hash_audit.csv"), check.names = FALSE)
artifact <- jsonlite::read_json(file.path(output_dir, "WGCNA_nature_readiness_report_artifact.json"), simplifyVector = FALSE)

checks <- data.frame(
  check_id = c(
    "required_files_exist", "required_files_nonzero", "no_interrupted_temp_files",
    "module_sensitivity_unique", "loao_unique", "block_rows_unique",
    "context_rows_unique", "portable_artifact_contract", "portable_html_nonzero",
    "protected_hashes_unchanged"
  ),
  passed = c(
    all(file.exists(required_paths)),
    all(file.exists(required_paths)) && all(file.info(required_paths)$size > 0),
    !any(grepl("\\.tmp-|\\.part$|\\.incomplete$", list.files(output_dir, all.files = TRUE))),
    nrow(robustness) == 39L && !anyDuplicated(robustness[c("sensitivity_id", "ModuleID")]),
    nrow(loao) == 117L && !anyDuplicated(loao[c("omitted_AnimalID", "ModuleID")]),
    nrow(blocks) == 12L && !anyDuplicated(blocks[c("sensitivity_id", "SupermoduleID")]),
    nrow(context) == 13L && !anyDuplicated(context$ModuleID),
    identical(artifact$surface, "report") && identical(artifact$manifest$version, 1L) && identical(artifact$snapshot$version, 1L),
    file.exists(file.path(output_dir, "WGCNA_nature_readiness_report.html")) && file.info(file.path(output_dir, "WGCNA_nature_readiness_report.html"))$size > 100000,
    all(protected$unchanged)
  ),
  detail = c(
    paste(required, collapse = ";"),
    "all required audit artifacts have nonzero size",
    "no .tmp, .part or .incomplete files",
    "3 sensitivity matrices x 13 stable modules",
    "9 animals x 13 stable modules",
    "4 matrices x 3 multi-module blocks",
    "one context row per stable module",
    "surface=report; manifest.version=1; snapshot.version=1",
    "canonical portable builder output exceeds 100 KB",
    "primary state, Stage 05, registry and current publication figures retain pre-run MD5 hashes"
  ),
  stringsAsFactors = FALSE
)
atomic_write_csv(checks, file.path(output_dir, "finalization_validation.csv"))
if (!all(checks$passed)) stop("Finalization checks failed: ", paste(checks$check_id[!checks$passed], collapse = ", "))

files <- list.files(output_dir, recursive = TRUE, full.names = TRUE, all.files = TRUE)
files <- files[file.exists(files) & !dir.exists(files)]
files <- files[basename(files) != "audit_manifest.csv"]
relative <- substring(normalizePath(files, winslash = "/"), nchar(repo_root) + 2L)
manifest <- data.frame(
  relative_path = relative,
  bytes = unname(file.info(files)$size),
  md5 = unname(tools::md5sum(files)),
  generated_at = format(Sys.time(), tz = "Europe/Berlin", usetz = TRUE),
  audit_scope = "additive_fixed_membership_microglia_WGCNA",
  stringsAsFactors = FALSE
)
atomic_write_csv(manifest, file.path(output_dir, "audit_manifest.csv"))
message("Finalized audit manifest with ", nrow(manifest), " files")

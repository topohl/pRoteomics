#!/usr/bin/env Rscript
# Frozen, downstream-only source for the WGCNA module wheel and the theme x
# module composition panel (manuscript candidate). Copies the tables written by
# analysis/integration/test_theme_module_composition.R byte for byte after
# checking their schema, their internal consistency and that the run that wrote
# them used the producer code as committed at HEAD. Nothing is recomputed.
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(root, "pipeline.yml")))
  stop("Run from the pRoteomics repository root.", call. = FALSE)
rel <- function(...) file.path(root, ...)
git <- function(...) {
  out <- suppressWarnings(system2("git", c("-C", shQuote(root), ...),
                                  stdout = TRUE, stderr = TRUE))
  list(out = out, status = if (is.null(attr(out, "status"))) 0L else attr(out, "status"))
}
publication_id <- "wgcna_theme_composition_candidate"
exporter <- "tools/export_wgcna_theme_composition_candidate.R"
producers <- c("analysis/integration/test_theme_module_composition.R",
               "R/enrichment/theme_module_composition_utils.R",
               "config/manuscript_go_theme_registry.tsv",
               "config/gsea_wgcna_theme_module_mapping.csv")
src <- rel("results", "integration", "test_theme_module_composition", "global")
out <- rel("exports", "publication_source_data", publication_id)
manifest_path <- rel("exports", "publication_source_data", "manifest.csv")
if (file.exists(out)) stop("Candidate export already exists: ", out, call. = FALSE)

## Provenance: the exporter must be committed (so source_commit contains it),
## and the run must come from a commit whose producer code equals HEAD's.
head <- git("rev-parse", "HEAD")
if (head$status != 0L) stop("git rev-parse HEAD failed.", call. = FALSE)
head <- head$out[[1]]
dirty <- git("status", "--porcelain", "--", exporter, producers)$out
if (length(dirty)) stop("Exporter or producer files differ from HEAD:\n",
                        paste(dirty, collapse = "\n"), call. = FALSE)
run_manifest <- readLines(file.path(src, "manifests", "run_manifest.yml"), warn = FALSE)
run_commit <- sub("^git_commit:\\s*", "", grep("^git_commit:", run_manifest, value = TRUE))
if (length(run_commit) != 1L || !grepl("^[0-9a-f]{40}$", run_commit))
  stop("The composition run manifest has no git_commit.", call. = FALSE)
if (git("merge-base", "--is-ancestor", run_commit, head)$status != 0L)
  stop("Run commit ", run_commit, " is not an ancestor of HEAD.", call. = FALSE)
if (git("diff", "--quiet", run_commit, head, "--", producers)$status != 0L)
  stop("Producer code changed between the run (", substr(run_commit, 1, 7),
       ") and HEAD; re-run the composition step first.", call. = FALSE)

tables <- c("theme_entity_composition_long", "theme_composition_summary",
            "theme_composition_summary_all_variants", "theme_partition_agreement",
            "curated_mapping_composition_audit", "wgcna_module_wheel_structure",
            "composition_universe_audit")
paths <- setNames(file.path(src, "tables", paste0(tables, ".csv")), tables)
mirrors <- file.path(src, "tables", "source_data", paste0(tables, ".csv"))
if (!all(file.exists(paths))) stop("Missing composition tables: ",
                                   paste(paths[!file.exists(paths)], collapse = ", "), call. = FALSE)
for (i in seq_along(tables)) {
  if (file.exists(mirrors[[i]]) &&
      !identical(unname(tools::sha256sum(mirrors[[i]])), unname(tools::sha256sum(paths[[i]]))))
    stop("tables/ and tables/source_data/ disagree for ", tables[[i]], call. = FALSE)
}
x <- lapply(paths, utils::read.csv, stringsAsFactors = FALSE, check.names = FALSE)

need <- list(
  theme_entity_composition_long = c("dataset", "level", "entity_id", "member_modules", "theme_id",
    "theme_set_definition", "universe_definition", "entity_definition", "n_universe", "n_theme",
    "n_entity", "n_overlap", "frac_theme_in_entity", "fold_enrichment", "odds_ratio", "BH_FDR",
    "composition_supported", "evidence_role", "is_primary_variant"),
  theme_composition_summary = c("dataset", "theme_id", "n_theme", "n_carriers", "carriers",
    "share_captured_by_carriers", "top_carrier", "composition_class", "theme_label"),
  theme_partition_agreement = c("dataset", "level", "cramers_v", "cramers_v_perm_max",
    "cramers_v_perm_p", "nmi", "nmi_perm_max", "nmi_perm_p", "n_permutations", "permutation_seed"),
  wgcna_module_wheel_structure = c("dataset", "supermodule_id", "module_id", "n_proteins",
    "peak_unit", "spatial_tau", "bilateral_reproducibility_class", "n_member_modules",
    "is_multi_module_supermodule", "carried_themes", "top_theme", "top_theme_share"),
  curated_mapping_composition_audit = c("dataset", "theme_id", "entity_level", "entity_id",
    "curated_mapping_composition_status"),
  composition_universe_audit = c("dataset", "n_primary_universe", "n_modules",
    "theme_registry_version", "org_Mm_eg_db_version"))
for (nm in names(need)) {
  miss <- setdiff(need[[nm]], names(x[[nm]]))
  if (length(miss)) stop("Schema mismatch in ", nm, ": ", paste(miss, collapse = ", "), call. = FALSE)
}

datasets <- c("neuron_neuropil", "neuron_soma", "microglia")
long <- x$theme_entity_composition_long
wheel <- x$wgcna_module_wheel_structure
summ <- x$theme_composition_summary
if (!setequal(unique(long$dataset), datasets) || !setequal(unique(wheel$dataset), datasets) ||
    !setequal(unique(summ$dataset), datasets))
  stop("Every table must cover exactly the three compartments.", call. = FALSE)
if (!all(long$evidence_role == "structural_composition_not_phenotype"))
  stop("Unexpected evidence_role in the composition table.", call. = FALSE)
if (anyDuplicated(wheel[c("dataset", "module_id")]) || anyNA(wheel$supermodule_id))
  stop("The wheel table must list each module once, with a supermodule.", call. = FALSE)
prim <- long[long$is_primary_variant %in% TRUE & long$level == "module", , drop = FALSE]
if (!setequal(paste(prim$dataset, prim$entity_id), paste(wheel$dataset, wheel$module_id)))
  stop("Primary module rows and the wheel table list different modules.", call. = FALSE)
if (anyDuplicated(summ[c("dataset", "theme_id")]))
  stop("The primary summary has duplicate dataset x theme rows.", call. = FALSE)
for (i in seq_len(nrow(summ))) {
  cur <- prim[prim$dataset == summ$dataset[[i]] & prim$theme_id == summ$theme_id[[i]] &
                prim$composition_supported %in% TRUE, "entity_id"]
  stated <- if (nzchar(summ$carriers[[i]])) strsplit(summ$carriers[[i]], ";", fixed = TRUE)[[1]] else character()
  if (!setequal(cur, stated) || length(cur) != summ$n_carriers[[i]])
    stop("Summary carriers disagree with the supported primary rows: ",
         summ$dataset[[i]], " ", summ$theme_id[[i]], call. = FALSE)
}
for (i in seq_len(nrow(wheel))) {
  cur <- prim[prim$dataset == wheel$dataset[[i]] & prim$entity_id == wheel$module_id[[i]] &
                prim$composition_supported %in% TRUE, "theme_id"]
  stated <- if (nzchar(wheel$carried_themes[[i]])) strsplit(wheel$carried_themes[[i]], ";", fixed = TRUE)[[1]] else character()
  if (!setequal(cur, stated))
    stop("Wheel carried_themes disagree with the supported primary rows: ",
         wheel$dataset[[i]], " ", wheel$module_id[[i]], call. = FALSE)
}

## Theme order and registry labels, so the renderer never hard-codes them.
registry_path <- rel("config", "manuscript_go_theme_registry.tsv")
registry <- utils::read.delim(registry_path, stringsAsFactors = FALSE, check.names = FALSE)
theme_index <- unique(registry[registry$theme_role == "primary",
                               c("theme_id", "display_label", "display_order", "registry_version")])
theme_index <- theme_index[order(theme_index$display_order), , drop = FALSE]
if (anyDuplicated(theme_index$theme_id) || !setequal(theme_index$theme_id, unique(long$theme_id)) ||
    !identical(unique(theme_index$registry_version), unique(x$composition_universe_audit$theme_registry_version)))
  stop("Theme registry and composition tables disagree.", call. = FALSE)

dir.create(out, recursive = TRUE, showWarnings = FALSE)
for (nm in tables) {
  if (!file.copy(paths[[nm]], file.path(out, paste0(nm, ".csv"))))
    stop("Copy failed: ", nm, call. = FALSE)
}
utils::write.csv(theme_index, file.path(out, "theme_index.csv"), row.names = FALSE)
inputs <- data.frame(
  role = c(paste0("composition_table:", tables), "composition_run_manifest", "theme_registry"),
  path = c(paths, file.path(src, "manifests", "run_manifest.yml"), registry_path),
  stringsAsFactors = FALSE)
inputs$path <- substring(gsub("\\\\", "/", inputs$path), nchar(root) + 2L)
inputs$sha256 <- unname(tools::sha256sum(file.path(root, inputs$path)))
inputs$annotation_package <- unique(x$composition_universe_audit$org_Mm_eg_db_version)
inputs$source_run_commit <- run_commit
utils::write.csv(inputs, file.path(out, "input_manifest.csv"), row.names = FALSE)

## Publication manifest rows, stamped at HEAD (which contains this exporter).
man <- utils::read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
if (publication_id %in% man$publication_id)
  stop("manifest.csv already lists ", publication_id, call. = FALSE)
files <- c(paste0(tables, ".csv"), "theme_index.csv", "input_manifest.csv")
add <- do.call(rbind, lapply(files, function(f) {
  p <- file.path(out, f)
  d <- utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
  r <- file.path("exports", "publication_source_data", publication_id, f)
  data.frame(publication_id = publication_id, source_repo = "topohl/pRoteomics",
             source_commit = head, source_analysis = exporter, source_table = r,
             exported_file = r, rows = nrow(d), columns = ncol(d),
             sha256 = unname(tools::sha256sum(p)),
             contract_version = "publication_source_data_v1", stringsAsFactors = FALSE)
}))
raw <- readBin(manifest_path, "raw", file.info(manifest_path)$size)
eol <- if (any(raw == as.raw(13L))) "\r\n" else "\n"
q <- function(v) paste0("\"", gsub("\"", "\"\"", v, fixed = TRUE), "\"")
lines <- vapply(seq_len(nrow(add)), function(i) with(add[i, ], paste(
  q(publication_id), q(source_repo), q(source_commit), q(source_analysis), q(source_table),
  q(exported_file), rows, columns, q(sha256), q(contract_version), sep = ",")), character(1))
con <- file(manifest_path, open = "ab")
writeBin(charToRaw(paste0(paste(lines, collapse = eol), eol)), con)
close(con)
cat("Exported ", length(files), " files for ", publication_id, " at ", substr(head, 1, 7),
    " (run ", substr(run_commit, 1, 7), ") to ", out, "\n", sep = "")

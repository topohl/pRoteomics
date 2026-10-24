#!/usr/bin/env Rscript
# Does the WGCNA gene mapping change any module's GO evidence or label? (normalisation review, 2026-10-23). Read-only.
#
# The WGCNA build maps each protein's representative UniProt accession to Entrez with clusterProfiler::bitr (org.Mm.eg.db)
# and uses an older gene annotation; accessions missing from org.Mm.eg.db's UniProt cross-references get no Entrez ID and
# drop out of every module GO test. The differential-abundance tables carry the current official mapping
# (official_entrez_id, annotation status 'concordant_official_gene'). This check
#   1. counts, per network, proteins without an Entrez ID in the build but with an official one, and any disagreement;
#   2. recomputes every module's GO over-representation (BP, CC, MF; module 'all' set) exactly as the build does
#      (clusterProfiler::enrichGO, universe = the network's gene-level proteins, BH) with the official Entrez IDs added,
#      after confirming that the same code reproduces the stored results with the build's own mapping;
#   3. writes model-state copies carrying the corrected mapping and GO tables, for the rule-based labeller.
# Writes results/reviewer_audit/normalisation_choice_2026-10-21/entrez_mapping_check/ (tables) and the state copies to the
# directory given as the first argument.
# USAGE  Rscript tools/normalisation_review/check_entrez_mapping.R <state copy dir>      (working directory: repo root)

suppressPackageStartupMessages({ library(dplyr); library(readr); library(clusterProfiler); library(org.Mm.eg.db) })
args <- commandArgs(TRUE); ALT <- args[1]; dir.create(ALT, recursive = TRUE, showWarnings = FALSE)
ARCHIVE <- if (dir.exists("Q:/data")) "Q:" else
  "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/_archive_proteomics_pre_renormalisation_2026-10-19"
OUT <- file.path("results", "reviewer_audit", "normalisation_choice_2026-10-21", "entrez_mapping_check"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
STATE <- list(original = function(ds) file.path(ARCHIVE, "data", "processed", "06_modules_WGCNA", "01_WGCNA", ds, "wgcna_final_model_state.rds"),
              new = function(ds) file.path("results", "wgcna", "build_wgcna_modules", ds, "models", "wgcna_final_model_state.rds"))
DATASETS <- c("neuron_neuropil", "neuron_soma", "microglia")
is_module <- function(m) !is.na(m) & !grepl("grey|unassigned|m00", m)
# official mapping from the differential-abundance tables (same protein groups in every contrast file of a compartment)
official <- function(ds) {
  fs <- list.files(file.path("data", "processed", "02_id_mapping_animal_level", "mapped", ds, "forward", "per_file"), full.names = TRUE)
  bind_rows(lapply(fs, function(f) read_csv(f, show_col_types = FALSE, col_types = cols(.default = "c")) |>
                     transmute(original_identifier, official_entrez_id, official_gene_symbol, protein_group_gene_annotation_status))) |>
    distinct() |> filter(protein_group_gene_annotation_status == "concordant_official_gene", !is.na(official_entrez_id))
}
enrich <- function(genes, universe, ont) {
  if (length(genes) < 5) return(NULL)
  ego <- suppressWarnings(enrichGO(gene = genes, universe = universe, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = ont,
                                   pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = FALSE))
  df <- as.data.frame(ego); if (!nrow(df)) return(NULL)
  as_tibble(df) |> mutate(Ontology = ont)
}
counts <- list(); repro <- list()
for (v in names(STATE)) for (ds in DATASETS) {
  s <- readRDS(STATE[[v]](ds)); long <- as.data.frame(s$WGCNA_modules_long)
  off <- official(ds); stopifnot(!anyDuplicated(off$original_identifier))
  long$official <- off$official_entrez_id[match(long$original_identifier, off$original_identifier)]
  g <- long$gene_level_claim_allowed %in% TRUE
  counts[[length(counts) + 1L]] <- tibble(version = v, dataset = ds, proteins = nrow(long), gene_level = sum(g),
    entrez_in_build = sum(g & !is.na(long$EntrezID)), official_entrez = sum(g & !is.na(long$official)),
    added_by_official = sum(g & is.na(long$EntrezID) & !is.na(long$official)),
    disagree = sum(g & !is.na(long$EntrezID) & !is.na(long$official) & long$EntrezID != long$official),
    build_only = sum(g & !is.na(long$EntrezID) & is.na(long$official)), neither = sum(g & is.na(long$EntrezID) & is.na(long$official)))
  # corrected mapping: the official Entrez ID where the build has none (the build's own IDs are kept where they exist)
  fixed <- ifelse(g & is.na(long$EntrezID) & !is.na(long$official), long$official, as.character(long$EntrezID))
  uni_old <- unique(na.omit(as.character(long$EntrezID[g]))); uni_new <- unique(na.omit(fixed[g]))
  mods <- sort(unique(long$ModuleID[is_module(long$ModuleID)]))
  go_new <- list(); go_old <- list()
  done_f <- file.path(ALT, sprintf("%s_%s.rds", v, ds)); reuse <- file.exists(done_f) && ds != "neuron_soma"   # resume: reuse a finished network
  if (reuse) go_new <- list(readRDS(done_f)$GO_enrichment |> dplyr::select(-ModuleProteinSetType))
  if (!reuse) for (m in mods) for (ont in c("BP", "CC", "MF")) {
    sel <- g & long$ModuleID == m
    r <- enrich(unique(na.omit(fixed[sel])), uni_new, ont); if (!is.null(r)) go_new[[length(go_new) + 1L]] <- r |> mutate(ModuleID = m)
    if (ds == "neuron_soma") { r0 <- enrich(unique(na.omit(as.character(long$EntrezID[sel]))), uni_old, ont); if (!is.null(r0)) go_old[[length(go_old) + 1L]] <- r0 |> mutate(ModuleID = m) }
  }
  go_new <- bind_rows(go_new) |> mutate(ModuleProteinSetType = "all")
  if (ds == "neuron_soma") {   # the same code with the build's mapping must reproduce the stored GO tables
    st <- as.data.frame(s$GO_enrichment) |> filter(ModuleProteinSetType == "all") |> dplyr::select(ModuleID, Ontology, ID, p_stored = p.adjust, n_stored = Count)
    cmp <- bind_rows(go_old) |> dplyr::select(ModuleID, Ontology, ID, p_new = p.adjust, n_new = Count) |> full_join(st, by = c("ModuleID", "Ontology", "ID"))
    repro[[v]] <- tibble(version = v, dataset = ds, terms_stored = nrow(st), terms_recomputed = sum(!is.na(cmp$p_new)), unmatched = sum(is.na(cmp$p_new) | is.na(cmp$p_stored)),
                         max_abs_diff_q = max(abs(cmp$p_new - cmp$p_stored), na.rm = TRUE), count_mismatch = sum(cmp$n_new != cmp$n_stored, na.rm = TRUE))
    print(repro[[v]])
  }
  long$EntrezID <- fixed
  long$GeneSymbol <- ifelse(is.na(long$GeneSymbol) & !is.na(long$official), off$official_gene_symbol[match(long$original_identifier, off$original_identifier)], long$GeneSymbol)
  saveRDS(list(WGCNA_modules_long = long, GO_enrichment = go_new), file.path(ALT, sprintf("%s_%s.rds", v, ds)))
  cat(v, ds, ": added", counts[[length(counts)]]$added_by_official, "| significant terms before", sum(as.data.frame(s$GO_enrichment)$ModuleProteinSetType == "all" & as.data.frame(s$GO_enrichment)$p.adjust < 0.05),
      "after", sum(go_new$p.adjust < 0.05), "\n")
}
write_csv(bind_rows(counts), file.path(OUT, "mapping_counts.csv"))
write_csv(bind_rows(repro), file.path(OUT, "go_reproduction_check.csv"))
cat("done\n")

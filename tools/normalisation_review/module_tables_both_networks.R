# Side-by-side module and supermodule tables for both networks, and the Figure 3 exemplar check for each network
# (normalisation review 2026-10-22). Read-only. Functional annotation (scripted class, confidence, label), network
# stability (leave-one-animal-out) and biological group/phenotype statistics are kept in separate column blocks.
#
#   annotation   annotation_review/rule_based/{module_rule_classes,supermodule_rule_classes}.csv and labels.tsv
#   correspondence map_wgcna_module_lineage (shared proteins, share of each module, Jaccard)
#   group effects  each network's own module_group_effects.csv / supermodule effects (original: archive, original code;
#                  complete-case: working tree, corrected code): smallest p and FDR (primary global family), and the
#                  smallest local (unit-level, exploratory) FDR
#   phenotype      Behavior_coupling/stats/associations.csv (original network only): tests, smallest p and FDR, CombZ rho
#   Figure 3     per atlas exemplar term (neuropil synaptic signaling GO:0099536, soma mRNA processing GO:0006397,
#                microglia oxidative phosphorylation GO:0006119): significant BP terms in the exemplar's theme per
#                module, the exemplar term's enrichment per module, and how the term's GSEA leading edge (significant
#                contrasts of the exemplar unit; that network's own GSEA) distributes over the modules
# Writes original_assessment/modules_both_networks.csv, supermodules_both_networks.csv, figure3_exemplar_check.csv.
# USAGE  Rscript tools/normalisation_review/module_tables_both_networks.R      (working directory: repo root)

suppressPackageStartupMessages({ library(dplyr); library(readr); library(tidyr); library(org.Mm.eg.db) })
source(file.path("R", "paths.R")); source(repo_path("R", "manuscript_go_theme_utils.R"))
ARCHIVE <- if (dir.exists("Q:/data")) "Q:" else
  "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/_archive_proteomics_pre_renormalisation_2026-10-19"
# The original network: the original-centring regeneration on the corrected code (2026-10-24/25; worktree
# Analysis/proteomics_original, subst R:). Its membership equals the archive; module GO uses the same annotation
# and term-size limits as the complete-case network.
ORIG <- Sys.getenv("PROTEOMICS_ORIGINAL_ROOT", unset = if (dir.exists("R:/results")) "R:" else
  "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics_original")
REV <- file.path("results", "reviewer_audit", "normalisation_choice_2026-10-21")
RB <- file.path(REV, "annotation_review", "rule_based"); OUT <- file.path(REV, "original_assessment")
DATASETS <- c("neuron_neuropil", "neuron_soma", "microglia")
rd <- function(f) if (file.exists(f)) read_csv(f, show_col_types = FALSE, guess_max = 1e5) else NULL
mods <- rd(file.path(RB, "module_rule_classes.csv")); sms <- rd(file.path(RB, "supermodule_rule_classes.csv"))
lab <- read_tsv(file.path(RB, "labels.tsv"), show_col_types = FALSE, quote = "")
PATHS <- list(
  new = list(state = function(ds) file.path("results", "wgcna", "build_wgcna_modules", ds, "models", "wgcna_final_model_state.rds"),
             mge = function(ds) file.path("results", "wgcna", "test_module_phenotypes", ds, "tables", "module_group_effects.csv"),
             sge = function(ds) file.path("results", "wgcna", "test_module_phenotypes", ds, "tables", "supermodule_group_effects.csv"),
             cpl = function(ds) file.path("results", "integration", "test_module_behaviour_coupling", ds, "tables", "module_behavior_coupling.csv"),
             gsea = function(ds) file.path("results", "differential_abundance", "run_clusterprofiler_enrichment", ds, "models", "clusterProfiler_manifest.csv"),
             root = "P:/"),
  original = list(state = function(ds) file.path(ORIG, "results", "wgcna", "build_wgcna_modules", ds, "models", "wgcna_final_model_state.rds"),
                  mge = function(ds) file.path(ORIG, "results", "wgcna", "test_module_phenotypes", ds, "tables", "module_group_effects.csv"),
                  sge = function(ds) file.path(ORIG, "results", "wgcna", "test_module_phenotypes", ds, "tables", "supermodule_group_effects.csv"),
                  cpl = function(ds) file.path(ORIG, "results", "integration", "test_module_behaviour_coupling", ds, "tables", "module_behavior_coupling.csv"),
                  gsea = function(ds) file.path(ORIG, "results", "differential_abundance", "run_clusterprofiler_enrichment", ds, "models", "clusterProfiler_manifest.csv"),
                  root = paste0(ORIG, "/")))
fix_path <- function(p, root) sub("^[A-Za-z]:/+", root, gsub("\\\\", "/", p))
min_or_na <- function(x) if (length(x) && any(!is.na(x))) min(x, na.rm = TRUE) else NA_real_
# gene annotation and module GO limits of a WGCNA state (states built before 2026-10-23 used the legacy annotation and 10-500)
state_evidence <- function(s) paste(s$parameters$gene_annotation_contract_version %||% "legacy_uniprot_gene_name",
                                    s$parameters$module_go_min_gs_size %||% 10L, s$parameters$module_go_max_gs_size %||% 500L, sep = " / ")
evidence_seen <- character()

group_summary <- function(f, id_col) {
  x <- rd(f); if (is.null(x)) return(NULL)
  id <- intersect(c(id_col, "endpoint_id", "module_id", "supermodule_id"), names(x))[1]
  x$id <- x[[id]]
  glob <- if ("effect_scope" %in% names(x)) grepl("global", x$effect_scope) else rep(TRUE, nrow(x))
  fdr_g <- intersect(c("FDR_primary_global", "FDR_within_dataset_level", "FDR_global", "BH_q"), names(x))[1]
  fdr_l <- intersect(c("FDR_local_exploratory", "FDR_within_dataset_level"), names(x))[1]
  x |> mutate(glob = glob) |> group_by(id) |> summarise(
    group_min_p_global = min_or_na(p_value[glob]), group_min_FDR_global = min_or_na(.data[[fdr_g]][glob]),
    group_best_contrast = { k <- which(glob)[which.min(p_value[glob])]; if (length(k)) paste(contrast[k], sprintf("%+.2f", estimate[k])) else "" },
    group_min_FDR_local = min_or_na(.data[[fdr_l]][!glob]),
    group_best_local = { k <- which(!glob)[which.min(p_value[!glob])]; if (length(k)) paste(spatial_unit[k], contrast[k]) else "" }, .groups = "drop")
}
# module-phenotype associations (CombZ and independent traits, animal level): Behavior_coupling/ (built 2026-10-19 on the
# original network; not rerun on the complete-case network, whose pipeline coupling table carries no statistics)
coupling_summary <- function(v, ds, level = "module") {
  if (v != "original") return(NULL)
  x <- rd(file.path("Behavior_coupling", "stats", "associations.csv")); if (is.null(x)) return(NULL)
  # tests: CombZ (outcome) and the independent traits; components of CombZ are built in and not counted (Behavior_coupling/README)
  x |> filter(dataset == ds, .data$level == !!level, class %in% c("outcome", "independent")) |> mutate(trait = id, id = sub("^.*[|]", "", key)) |> group_by(id) |>
    summarise(phenotype_n_tests = n(), phenotype_n_p05 = sum(p_value < 0.05, na.rm = TRUE), phenotype_min_p = min_or_na(p_value),
              phenotype_min_FDR_family = min_or_na(p_adjusted), phenotype_min_FDR_all = min_or_na(p_adjusted_all),
              CombZ_rho = rho[trait == "combz"][1], CombZ_p = p_value[trait == "combz"][1], CombZ_p_perm = p_perm[trait == "combz"][1], .groups = "drop")
}

THEME_REG <- read_manuscript_go_theme_registry("config/manuscript_go_theme_registry.tsv")
EX <- tibble(dataset = DATASETS, unit = c("CA3 sr", "CA2 sp", "CA1"), term = c("GO:0099536", "GO:0006397", "GO:0006119"),
             term_name = c("synaptic signaling", "mRNA processing", "oxidative phosphorylation"),
             theme = c("synaptic_signaling_vesicle", "rna_processing_splicing_rnp", "mitochondrial_respiration_oxphos"))
frac <- function(r) as.numeric(sub("/.*", "", r)) / as.numeric(sub(".*/", "", r))

mt <- list(); st <- list(); fig <- list()
for (v in names(PATHS)) for (ds in DATASETS) {
  P <- PATHS[[v]]
  gm <- group_summary(P$mge(ds), "endpoint_id"); gsm <- group_summary(P$sge(ds), "supermodule_id"); cp <- coupling_summary(v, ds); cps <- coupling_summary(v, ds, "supermodule")
  m <- mods |> filter(version == v, dataset == ds) |>
    left_join(lab |> filter(version == v, level == "module") |> dplyr::select(dataset, ModuleID = entity_id, label), by = c("dataset", "ModuleID"))
  if (!is.null(gm)) m <- m |> left_join(gm, by = c("ModuleID" = "id"))
  if (!is.null(cp)) m <- m |> left_join(cp, by = c("ModuleID" = "id"))
  mt[[length(mt) + 1L]] <- m
  s <- sms |> filter(version == v, dataset == ds) |>
    left_join(lab |> filter(version == v, level == "supermodule") |> dplyr::select(dataset, SupermoduleID = entity_id, label), by = c("dataset", "SupermoduleID"))
  if (!is.null(gsm)) s <- s |> left_join(gsm, by = c("SupermoduleID" = "id"))
  if (!is.null(cps)) s <- s |> left_join(cps, by = c("SupermoduleID" = "id"))
  st[[length(st) + 1L]] <- s
  # ---- Figure 3 exemplar check on this network
  e <- EX[EX$dataset == ds, ]
  sobj <- readRDS(P$state(ds)); long <- as.data.frame(sobj$WGCNA_modules_long)
  evidence_seen <- unique(c(evidence_seen, state_evidence(sobj)))
  if (length(evidence_seen) > 1L) stop("The two networks carry different gene annotation or GO limits (", paste(evidence_seen, collapse = " vs "),
                                       "); regenerate both with the same code before comparing them.", call. = FALSE)
  go <- as.data.frame(sobj$GO_enrichment) |> filter(ModuleProteinSetType == "all")
  bp_sig <- go |> filter(Ontology == "BP", p.adjust < 0.05) |> distinct(ModuleID, ID, Description)
  asg <- map_go_terms_to_manuscript_themes(distinct(bp_sig, ID, Description), THEME_REG)$assignments |> distinct(ID = GO_ID, theme_id)
  in_theme <- bp_sig |> inner_join(asg, by = "ID", relationship = "many-to-many") |> filter(theme_id == e$theme) |> distinct(ModuleID, ID) |> count(ModuleID, name = "n_sig_BP_in_theme")
  term <- go |> filter(ID == e$term, Ontology == "BP") |> transmute(ModuleID, exemplar_padj = p.adjust, exemplar_overlap = GeneRatio, exemplar_fold = frac(GeneRatio) / frac(BgRatio))
  man <- rd(P$gsea(ds))
  le_share <- NULL
  if (!is.null(man)) {
    man <- man |> filter(result_type == "GSEA_GO", ontology == "BP", !as.logical(simplified))
    uk <- gsub(" ", "", e$unit)
    man <- man[grepl(paste0("^", uk), man$contrast) | grepl(paste0("^", sub(" .*", "", e$unit), "(microglia)?(sus|res)"), man$contrast), ]
    le <- bind_rows(lapply(seq_len(nrow(man)), function(i) {
      f <- fix_path(man$output_table[i], P$root); if (!file.exists(f)) return(NULL)
      x <- read_csv(f, show_col_types = FALSE) |> filter(ID == e$term); if (!nrow(x) || x$p.adjust[1] >= 0.05) return(NULL)
      genes <- strsplit(x$core_enrichment[1], "/", fixed = TRUE)[[1]]
      ent <- suppressMessages(AnnotationDbi::mapIds(org.Mm.eg.db, genes, "ENTREZID", "SYMBOL", multiVals = "first"))
      mod <- coalesce(long$ModuleID[match(as.character(ent), as.character(long$EntrezID))], long$ModuleID[match(toupper(genes), toupper(long$GeneSymbol))], "not in network")
      tibble(contrast = man$contrast[i], module = mod)
    }))
    if (nrow(le)) le_share <- le |> count(contrast, module) |> group_by(contrast) |> mutate(s = n / sum(n)) |> ungroup() |>
      complete(contrast, module, fill = list(n = 0, s = 0)) |> group_by(ModuleID = module) |>
      summarise(leading_edge_share = sprintf("%.0f-%.0f%% (%d contrasts)", 100 * min(s), 100 * max(s), n()), .groups = "drop")
  }
  f <- m |> dplyr::select(ModuleID, label, signal_class, rule_confidence, stability_median_jaccard, stability_n_recovered, lineage_module, lineage_share) |>
    left_join(in_theme, by = "ModuleID") |> left_join(term, by = "ModuleID")
  if (!is.null(le_share)) f <- f |> left_join(le_share, by = "ModuleID") else f$leading_edge_share <- "no significant contrast"
  fig[[length(fig) + 1L]] <- f |> mutate(version = v, dataset = ds, exemplar = paste(e$term_name, e$term, e$unit), .before = 1) |>
    filter(!is.na(n_sig_BP_in_theme) | !is.na(exemplar_padj)) |> arrange(desc(coalesce(n_sig_BP_in_theme, 0L)), exemplar_padj)
}
write_csv(bind_rows(mt), file.path(OUT, "modules_both_networks.csv"), na = "")
write_csv(bind_rows(st), file.path(OUT, "supermodules_both_networks.csv"), na = "")
write_csv(bind_rows(fig), file.path(OUT, "figure3_exemplar_check.csv"), na = "")
cat("done\n")

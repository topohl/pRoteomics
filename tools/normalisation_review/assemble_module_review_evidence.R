# Evidence for the scientific review of the WGCNA module and supermodule annotations of the re-normalised network
# (normalisation review 2026-10-21). Read-only with respect to the pipeline.
#
# Per module: size, hubs (kME), GO BP/CC/MF over-representation from the network build (adjusted P, overlap, term
# coverage), the automatic label, the lineage to the previous network (main source module, share, class, Jaccard) and
# the previous label. Per supermodule: constituent modules and the DIRECT enrichment of the reconstructed protein set
# (union of the member modules) against the same tested-proteome universe the module enrichment used, plus agreement
# between the members' own enrichments. Also the rule-based zoom check: for each atlas exemplar GO term, which module
# carries it and how much of the term's GSEA leading edge falls into each module (matched by Entrez ID), the fold
# enrichment of the term per module, and the number of significant BP terms per module in the exemplar's theme
# (manuscript theme registry). Every GO row keeps its gene symbols so gene-level statements can be checked.
#
# Writes to results/reviewer_audit/normalisation_choice_2026-10-21/annotation_review/.
# USAGE  Rscript tools/normalisation_review/assemble_module_review_evidence.R      (working directory: repo root)

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(clusterProfiler); library(org.Mm.eg.db)
})
ARCHIVE <- if (dir.exists("Q:/data")) "Q:" else
  "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/_archive_proteomics_pre_renormalisation_2026-10-19"
OUT <- file.path("results", "reviewer_audit", "normalisation_choice_2026-10-21", "annotation_review")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
w <- function(x, name) write_csv(x, file.path(OUT, name), na = "")
DATASETS <- c("neuron_neuropil", "neuron_soma", "microglia")
SHORT <- c(neuron_neuropil = "neuropil", neuron_soma = "soma", microglia = "microglia")
TOP <- 8L
EXEMPLARS <- tibble(dataset = c("neuron_neuropil", "neuron_soma", "microglia"),
                    exemplar_unit = c("CA3 sr", "CA2 sp", "CA1"),
                    exemplar_term = c("GO:0099536", "GO:0006397", "GO:0006119"),
                    exemplar_description = c("synaptic signaling", "mRNA processing", "oxidative phosphorylation"),
                    zoom_module_chosen = c("WGCNA_m04", "WGCNA_m01", "WGCNA_m09"))

EXEMPLAR_THEME <- c(neuron_neuropil = "synaptic_signaling_vesicle", neuron_soma = "rna_processing_splicing_rnp",
                    microglia = "mitochondrial_respiration_oxphos")
source(file.path("R", "paths.R")); source(repo_path("R", "manuscript_go_theme_utils.R"))
THEME_REG <- read_manuscript_go_theme_registry("config/manuscript_go_theme_registry.tsv")
fmt_term <- function(d, n = TOP) {
  if (!nrow(d)) return("")
  d <- head(arrange(d, p.adjust), n)
  paste(sprintf("%s (q=%.1e; %s; cov %s/%s)", d$Description, d$p.adjust, d$GeneRatio, sub("/.*$", "", d$GeneRatio),
                sub("/.*$", "", d$BgRatio)), collapse = " | ")
}
old_label_table <- function(ds) {
  f <- file.path("config", "wgcna_labels", "reference_pre_renormalisation_2026-10-19", paste0(ds, ".csv"))
  so <- readRDS(file.path(ARCHIVE, "data", "processed", "06_modules_WGCNA", "01_WGCNA", ds, "wgcna_final_model_state.rds"))
  auto <- as.data.frame(so$WGCNA_modules_long) |> distinct(ModuleID, previous_automatic_label = ModuleLabel_Final)
  if (file.exists(f)) {
    rv <- read_csv(f, show_col_types = FALSE) |> filter(level == "module")
    auto <- left_join(auto, transmute(rv, ModuleID = entity_id, previous_reviewed_label = reviewed_biological_label), by = "ModuleID")
  }
  if (!"previous_reviewed_label" %in% names(auto)) auto$previous_reviewed_label <- NA_character_
  auto
}

mod_rows <- list(); go_long <- list(); sm_rows <- list(); sm_go_long <- list(); zoom <- list()
for (ds in DATASETS) {
  cat("==", ds, "\n")
  s <- readRDS(file.path("results", "wgcna", "build_wgcna_modules", ds, "models", "wgcna_final_model_state.rds"))
  long <- as.data.frame(s$WGCNA_modules_long)
  universe <- unique(na.omit(as.character(long$EntrezID[long$gene_level_claim_allowed %in% TRUE])))
  ent2sym <- setNames(as.character(long$GeneSymbol), as.character(long$EntrezID))
  go <- as.data.frame(s$GO_enrichment) |> filter(ModuleProteinSetType == "all")
  stopifnot(length(unique(go$UniverseSize)) == 1L)
  universe_check <- c(reconstructed = length(universe), build = unique(go$UniverseSize))
  lin <- file.path("results", "wgcna", "map_wgcna_module_lineage", ds, "tables")
  origin <- read_csv(file.path(lin, "new_module_origin.csv"), show_col_types = FALSE)
  pairs <- read_csv(file.path(lin, "module_lineage_pairs.csv"), show_col_types = FALSE)
  carry <- if (file.exists(file.path(lin, "label_carry_decisions.csv"))) read_csv(file.path(lin, "label_carry_decisions.csv"), show_col_types = FALSE) else NULL
  oldlab <- old_label_table(ds)
  sm <- read_csv(file.path("results", "wgcna", "annotate_module_microenvironment", ds, "tables", "WGCNA_supermodule_biological_annotation.csv"),
                 show_col_types = FALSE) |> dplyr::select(SupermoduleID, member_modules, Supermodule_FinalLabel, safe_display_label)
  mod_to_sm <- sm |> separate_rows(member_modules, sep = ";") |> transmute(ModuleID = member_modules, SupermoduleID)
  modules <- sort(setdiff(unique(long$ModuleID), c(NA, "WGCNA_m00", "grey", "WGCNA_grey")))
  modules <- modules[!grepl("grey|unassigned|m00", modules)]

  for (m in modules) {
    mm <- long[long$ModuleID == m, , drop = FALSE]
    hubs <- mm |> arrange(desc(abs_kME)) |> head(10)
    g <- go[go$ModuleID == m, , drop = FALSE]
    sig <- g[g$p.adjust < 0.05, , drop = FALSE]
    best_bp <- sig[sig$Ontology == "BP", , drop = FALSE] |> arrange(p.adjust) |> head(1)
    o <- origin[origin$new_module == m, , drop = FALSE]
    pj <- pairs |> filter(.data[[intersect(c("new_module", "new"), names(pairs))[1]]] == m)
    jcol <- intersect(c("jaccard", "Jaccard"), names(pairs))[1]
    ocol <- intersect(c("old_module", "old"), names(pairs))[1]
    best_pair <- if (nrow(pj) && !is.na(jcol)) pj[which.max(pj[[jcol]]), , drop = FALSE] else NULL
    main_old <- if (nrow(o)) o$main_old[1] else NA_character_
    cd <- if (!is.null(carry)) carry[carry$level == "module" & carry$new_entity_id %in% m, , drop = FALSE] else NULL
    mod_rows[[length(mod_rows) + 1L]] <- tibble(
      dataset = ds, compartment = SHORT[[ds]], ModuleID = m,
      SupermoduleID = paste(mod_to_sm$SupermoduleID[mod_to_sm$ModuleID == m], collapse = ";"),
      n_proteins = nrow(mm), n_gene_mapped = sum(!is.na(mm$EntrezID)), median_abs_kME = median(mm$abs_kME, na.rm = TRUE),
      n_core_kME_0.6 = sum(mm$abs_kME >= 0.6, na.rm = TRUE),
      hub_proteins_top10 = paste(sprintf("%s (%.2f)", coalesce(hubs$GeneSymbol, hubs$FeatureDisplayLabel), hubs$kME), collapse = ", "),
      automatic_label = mm$ModuleLabel_Final[1], automatic_label_source = mm$ModuleLabel_Source[1],
      n_sig_BP = sum(sig$Ontology == "BP"), n_sig_CC = sum(sig$Ontology == "CC"), n_sig_MF = sum(sig$Ontology == "MF"),
      top_BP = fmt_term(sig[sig$Ontology == "BP", ]), top_CC = fmt_term(sig[sig$Ontology == "CC", ]), top_MF = fmt_term(sig[sig$Ontology == "MF", ]),
      best_BP_fraction_of_module = if (nrow(best_bp)) as.numeric(sub("/.*", "", best_bp$GeneRatio)) / as.numeric(sub(".*/", "", best_bp$GeneRatio)) else NA_real_,
      lineage_class = if (nrow(o)) o$class[1] else NA_character_, main_previous_module = main_old,
      share_from_main_previous = if (nrow(o)) o$share_from_main_old[1] else NA_real_,
      share_from_previous_unassigned = if (nrow(o)) o$share_from_unassigned[1] else NA_real_,
      lineage_sources = if (nrow(o)) o$sources[1] else NA_character_,
      best_jaccard_previous = if (!is.null(best_pair)) best_pair[[jcol]] else NA_real_,
      best_jaccard_previous_module = if (!is.null(best_pair)) best_pair[[ocol]] else NA_character_,
      previous_automatic_label_of_main = oldlab$previous_automatic_label[match(main_old, oldlab$ModuleID)],
      previous_reviewed_label_of_main = oldlab$previous_reviewed_label[match(main_old, oldlab$ModuleID)],
      label_carry_decision = if (!is.null(cd) && nrow(cd)) paste(sprintf("%s from %s: %s", cd$decision, cd$old_entity_id, cd$reason), collapse = " | ") else NA_character_)
    go_long[[length(go_long) + 1L]] <- sig |> group_by(Ontology) |> arrange(p.adjust, .by_group = TRUE) |> slice_head(n = 15) |> ungroup() |>
      transmute(dataset = ds, ModuleID = m, Ontology, ID, Description, GeneRatio, BgRatio, p.adjust, Count,
                genes = vapply(strsplit(as.character(geneID), "/", fixed = TRUE), function(e) paste(sort(unique(na.omit(ent2sym[e]))), collapse = "/"), ""))
  }

  # ---- supermodules: direct enrichment of the union of member proteins
  for (i in seq_len(nrow(sm))) {
    mem <- strsplit(sm$member_modules[i], ";", fixed = TRUE)[[1]]
    genes <- unique(na.omit(as.character(long$EntrezID[long$ModuleID %in% mem & long$gene_level_claim_allowed %in% TRUE])))
    res <- bind_rows(lapply(c("BP", "CC", "MF"), function(ont) {
      e <- tryCatch(suppressWarnings(enrichGO(genes, universe = universe, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = ont,
                                              pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                                              minGSSize = s$parameters$module_go_min_gs_size %||% 10L,   # the module GO limits of this build
                                              maxGSSize = s$parameters$module_go_max_gs_size %||% 500L, readable = TRUE)),
                    error = function(e) NULL)
      if (is.null(e)) return(NULL)
      as.data.frame(e) |> mutate(Ontology = ont)
    }))
    sig <- filter(res, p.adjust < 0.05)
    sm_go_long[[length(sm_go_long) + 1L]] <- sig |> group_by(Ontology) |> arrange(p.adjust, .by_group = TRUE) |> slice_head(n = 15) |> ungroup() |>
      transmute(dataset = ds, SupermoduleID = sm$SupermoduleID[i], Ontology, ID, Description, GeneRatio, BgRatio, p.adjust, Count, geneID)
    top_sm <- head(arrange(filter(sig, Ontology == "BP"), p.adjust), 10)$ID
    sig_bp_ids <- filter(sig, Ontology == "BP")$ID
    top5 <- head(arrange(filter(sig, Ontology == "BP"), p.adjust), 5)
    member_info <- lapply(mem, function(m) {
      gm <- go[go$ModuleID == m & go$p.adjust < 0.05 & go$Ontology == "BP", ]
      top_m <- head(arrange(gm, p.adjust), 10)
      mg <- unique(na.omit(as.character(long$EntrezID[long$ModuleID == m])))
      hit_share <- if (nrow(top5)) mean(vapply(strsplit(top5$geneID, "/", fixed = TRUE), function(gs) {
        sy <- toupper(gs); mean(sy %in% toupper(ent2sym[mg])) }, numeric(1))) else NA_real_
      list(id = m, n = sum(long$ModuleID == m), top = top_m$Description[1] %||% NA_character_,
           shares = length(intersect(top_m$ID, top_sm)), sig_in_sm_top = any(top_m$ID %in% top_sm),
           own_top_sig_in_union = sum(top_m$ID %in% sig_bp_ids), n_own_top = nrow(top_m), hit_share = hit_share)
    })
    sm_rows[[length(sm_rows) + 1L]] <- tibble(
      dataset = ds, compartment = SHORT[[ds]], SupermoduleID = sm$SupermoduleID[i], member_modules = sm$member_modules[i],
      n_member_modules = length(mem),
      member_sizes = paste(vapply(member_info, function(x) sprintf("%s:%d", sub("WGCNA_", "", x$id), x$n), ""), collapse = ", "),
      n_proteins = sum(long$ModuleID %in% mem), n_gene_mapped = length(genes),
      current_label = sm$Supermodule_FinalLabel[i], current_safe_display_label = sm$safe_display_label[i],
      direct_n_sig_BP = sum(sig$Ontology == "BP"), direct_n_sig_CC = sum(sig$Ontology == "CC"), direct_n_sig_MF = sum(sig$Ontology == "MF"),
      direct_top_BP = fmt_term(filter(sig, Ontology == "BP")), direct_top_CC = fmt_term(filter(sig, Ontology == "CC")),
      direct_top_MF = fmt_term(filter(sig, Ontology == "MF")),
      member_top_BP = paste(vapply(member_info, function(x) sprintf("%s: %s", sub("WGCNA_", "", x$id), coalesce(x$top, "none")), ""), collapse = " | "),
      member_agreement = paste(vapply(member_info, function(x) sprintf("%s: %d of own top-10 BP in union top-10; %d/%d of own top-10 significant in union; %.0f%% of union top-5 BP hit genes",
                                                                    sub("WGCNA_", "", x$id), x$shares, x$own_top_sig_in_union, x$n_own_top, 100 * x$hit_share), ""), collapse = " | "))
  }

  # ---- zoom check: exemplar term in each module, and leading-edge distribution
  ex <- EXEMPLARS[EXEMPLARS$dataset == ds, ]
  frac <- function(r) as.numeric(sub("/.*", "", r)) / as.numeric(sub(".*/", "", r))
  term_by_mod <- go |> filter(ID == ex$exemplar_term, Ontology == "BP") |>
    transmute(ModuleID, exemplar_padj = p.adjust, exemplar_overlap = GeneRatio, exemplar_bg = BgRatio,
              exemplar_fold = frac(GeneRatio) / frac(BgRatio))
  theme <- EXEMPLAR_THEME[[ds]]
  bp_sig <- go |> filter(Ontology == "BP", p.adjust < 0.05) |> distinct(ModuleID, ID, Description)
  asg <- map_go_terms_to_manuscript_themes(distinct(bp_sig, ID, Description), THEME_REG)$assignments |> distinct(ID = GO_ID, theme_id)
  in_theme <- bp_sig |> inner_join(asg, by = "ID", relationship = "many-to-many") |> filter(theme_id == theme) |>
    distinct(ModuleID, ID) |> count(ModuleID, name = "n_sig_BP_terms_in_theme")
  term_by_mod <- full_join(term_by_mod, in_theme, by = "ModuleID") |> mutate(theme_id = theme)
  man <- read_csv(file.path("results", "differential_abundance", "run_clusterprofiler_enrichment", ds, "models", "clusterProfiler_manifest.csv"), show_col_types = FALSE) |>
    filter(result_type == "GSEA_GO", ontology == "BP", !as.logical(simplified))
  unit_key <- gsub(" ", "", ex$exemplar_unit)
  man <- man[grepl(paste0("^", unit_key), man$contrast) | grepl(paste0("^", sub(" .*", "", ex$exemplar_unit), "(microglia)?(sus|res)"), man$contrast), ]
  le <- bind_rows(lapply(seq_len(nrow(man)), function(i) {
    x <- read_csv(sub("^P:/+", "P:/", man$output_table[i]), show_col_types = FALSE) |> filter(ID == ex$exemplar_term)
    if (!nrow(x)) return(NULL)
    genes <- strsplit(x$core_enrichment, "/", fixed = TRUE)[[1]]
    ent <- suppressMessages(AnnotationDbi::mapIds(org.Mm.eg.db, genes, "ENTREZID", "SYMBOL", multiVals = "first"))
    by_ent <- long$ModuleID[match(as.character(ent), as.character(long$EntrezID))]
    by_sym <- long$ModuleID[match(toupper(genes), toupper(long$GeneSymbol))]
    tibble(contrast = man$contrast[i], NES = x$NES, p.adjust = x$p.adjust, n_leading_edge = length(genes),
           module = coalesce(by_ent, by_sym, "not in network"))
  }))
  if (nrow(le)) {
    zoom[[ds]] <- le |> count(contrast, NES, p.adjust, n_leading_edge, module, name = "n_leading_edge_in_module") |>
      group_by(contrast) |> mutate(share_of_leading_edge = n_leading_edge_in_module / n_leading_edge) |> ungroup() |>
      left_join(term_by_mod, by = c("module" = "ModuleID")) |>
      mutate(dataset = ds, exemplar_term = ex$exemplar_term, exemplar_description = ex$exemplar_description,
             exemplar_unit = ex$exemplar_unit, zoom_module_chosen = ex$zoom_module_chosen, .before = 1) |>
      arrange(contrast, desc(n_leading_edge_in_module))
  }
  w(term_by_mod |> mutate(dataset = ds, exemplar_term = ex$exemplar_term, .before = 1) |> arrange(exemplar_padj),
    sprintf("zoom_exemplar_term_by_module_%s.csv", ds))
  cat("  universe reconstructed/build:", universe_check, "\n")
}
disp <- read_csv(file.path("results", "wgcna", "build_module_label_registry", "global", "tables", "WGCNA_canonical_display_label_resolution.csv"),
                 show_col_types = FALSE) |> transmute(dataset, level, entity_id, displayed_label = canonical_display_label, displayed_label_source = canonical_label_source)
w(disp, "displayed_labels.csv")
w(bind_rows(mod_rows), "module_evidence.csv")
w(bind_rows(go_long), "module_go_top_terms.csv")
w(bind_rows(sm_rows), "supermodule_evidence.csv")
w(bind_rows(sm_go_long), "supermodule_direct_go_top_terms.csv")
w(bind_rows(zoom), "zoom_exemplar_leading_edge_by_module.csv")
cat("done:", OUT, "\n")

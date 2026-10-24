#!/usr/bin/env Rscript
# Old-to-new WGCNA module lineage after the 2026-10-19 re-normalisation.
#
# Old network: the versioned reference in the pre-renormalisation archive (the hash-verified copy
# of results/ taken before the rerun). New network: this checkout's rebuilt Stage 01 outputs.
# Proteins are matched by ProteinGroupID; proteins in a network's feature universe but in no
# module count as "unassigned". Read-only on both networks and on config/.
#
# Writes, per dataset, to results/wgcna/map_wgcna_module_lineage/<dataset>/tables/:
#   module_lineage_pairs.csv        every old x new module pair sharing proteins
#   old_module_correspondence.csv   one row per old module: class and main new module(s)
#   new_module_origin.csv           one row per new module: class and main old module(s)
#   module_go_comparison.csv        GO BP evidence of each old module re-tested in its best new module
#   supermodule_lineage.csv         old x new supermodules on the proteins of their member modules
#   label_carry_decisions.csv       one row per reviewed old label: carried or not, and why
#   theme_mapping_carry.csv         one row per old theme-to-entity mapping: carried or not, and why
#   provisional_label_registry.csv  (microglia) reviewed-registry schema, all provisional_unreviewed
# and, across datasets, provisional_theme_module_mapping.csv (gww schema) in .../global/tables/.
#
# Classes (old module o, best new module n by Jaccard):
#   retained   o and n are each other's best match and Jaccard >= 0.5
#   merged     >= 50% of o is in n, and n also holds >= 50% of another old module
#   changed    >= 50% of o is in n, but not reciprocal or Jaccard < 0.5
#   split      no new module holds 50% of o; two or more hold >= 20% each
#   dissolved  otherwise (scattered or mostly unassigned)
# A label or theme is carried only for a retained entity AND when the new entity's own GO
# enrichment supports it: at least half of the old entity's ten most significant GO terms (BP, CC
# and MF; reviewed labels also rest on compartment evidence) are significant again (BH < 0.05) in
# the new entity, or, for a theme, the new entity has a significant BP term inside that theme. Overlap alone never carries anything; every other entity
# stays provisionally unannotated for review.
#
# USAGE  Rscript tools/map_wgcna_module_lineage.R [--dataset all|<dataset>] [--old-root <dir>]

source(if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R"))
source(repo_path("R", "wgcna_paths.R"))
source(repo_path("R", "wgcna_reviewed_label_registry.R"))
source(repo_path("R", "manuscript_go_theme_utils.R"))
suppressPackageStartupMessages({ library(dplyr); library(readr) })

args <- commandArgs(trailingOnly = TRUE)
arg <- function(flag, default) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[[i + 1L]] }
DATASETS <- { d <- arg("--dataset", "all"); if (identical(d, "all")) c("neuron_neuropil", "neuron_soma", "microglia") else d }
OLD_ROOT <- arg("--old-root", "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/_archive_proteomics_pre_renormalisation_2026-10-19")
REVIEWED_REFERENCE <- repo_path("config", "wgcna_labels", "reference_pre_renormalisation_2026-10-19")
OLD_THEME_MAPPING <- repo_path("config", "gsea_wgcna_theme_module_mapping.csv")
THEME_REGISTRY <- repo_path("config", "manuscript_go_theme_registry.tsv")
FDR <- 0.05; TOP_N <- 10L

old_path <- function(ds, ...) file.path(OLD_ROOT, "results", "tables", "06_modules_WGCNA", ...)
new_path <- function(ds, ...) repo_path("results", "wgcna", "build_wgcna_modules", ds, "tables", ...)
need <- function(p) { if (!file.exists(p)) stop("Missing input: ", p, call. = FALSE); p }
rd <- function(p) readr::read_csv(need(p), show_col_types = FALSE, progress = FALSE, guess_max = 1e5)

read_network <- function(ds, side) {
  if (side == "old") {
    mods <- rd(old_path(ds, "01_WGCNA", ds, "modules", "WGCNA_modules_long.csv"))
    uni <- rd(old_path(ds, "01_WGCNA", ds, "modules", "WGCNA_feature_universe.csv"))
    go <- rd(old_path(ds, "01_WGCNA", ds, "modules", "WGCNA_module_GO_enrichment_long.csv"))
    smap <- rd(old_path(ds, "identity_contract", ds, "WGCNA_module_supermodule_membership_contract.csv")) |>
      transmute(ModuleID = module_id, SupermoduleID = supermodule_id)
  } else {
    mods <- rd(new_path(ds, "modules", "WGCNA_modules_long.csv"))
    uni <- rd(new_path(ds, "modules", "WGCNA_feature_universe.csv"))
    go <- rd(new_path(ds, "modules", "WGCNA_module_GO_enrichment_long.csv"))
    smap <- rd(new_path(ds, "supermodules", "wgcna_module_supermodule_annotation.csv")) |>
      transmute(ModuleID = as.character(ModuleID), SupermoduleID = as.character(SupermoduleID))
  }
  smap <- distinct(smap) |> filter(!is.na(ModuleID), !is.na(SupermoduleID))
  member <- distinct(transmute(mods, ProteinGroupID = as.character(ProteinGroupID), ModuleID = as.character(ModuleID)))
  if (anyDuplicated(member$ProteinGroupID)) stop(side, " ", ds, ": a protein is in two modules.", call. = FALSE)
  universe <- unique(as.character(uni$ProteinGroupID))
  assign <- setNames(rep("unassigned", length(universe)), universe)
  assign[member$ProteinGroupID] <- member$ModuleID
  sig <- go |>
    filter(.data$ModuleProteinSetType == "all", .data$Ontology %in% c("BP", "CC", "MF"), .data$p.adjust < FDR) |>
    transmute(ModuleID = as.character(ModuleID), Ontology = as.character(Ontology), ID = as.character(ID),
              Description = as.character(Description), p.adjust = as.numeric(p.adjust))
  list(assign = assign, smap = smap, sig = sig, sig_bp = filter(sig, Ontology == "BP"))
}

overlap_pairs <- function(a, b, a_name, b_name) {
  ids <- intersect(names(a), names(b))
  tab <- as.data.frame(table(old = a[ids], new = b[ids]), stringsAsFactors = FALSE) |> filter(Freq > 0)
  n_old <- table(a); n_new <- table(b)
  tab |>
    transmute(old = .data$old, new = .data$new, n_shared = as.integer(Freq),
              n_old = as.integer(n_old[.data$old]), n_new = as.integer(n_new[.data$new]),
              jaccard = n_shared / (n_old + n_new - n_shared),
              share_of_old = n_shared / n_old, share_of_new = n_shared / n_new) |>
    setNames(c(a_name, b_name, "n_shared", paste0("n_", a_name), paste0("n_", b_name), "jaccard", "share_of_old", "share_of_new"))
}

classify_old <- function(pairs, old_col, new_col) {
  p <- filter(pairs, .data[[old_col]] != "unassigned")
  real <- filter(p, .data[[new_col]] != "unassigned")
  best_new <- real |> group_by(.data[[old_col]]) |> slice_max(jaccard, n = 1, with_ties = FALSE) |> ungroup()
  best_old_for_new <- real |> group_by(.data[[new_col]]) |> slice_max(jaccard, n = 1, with_ties = FALSE) |> ungroup() |>
    transmute(new_id = .data[[new_col]], best_old_for_new = .data[[old_col]])
  holds_other <- real |> filter(share_of_old >= 0.5) |> count(new_id = .data[[new_col]], name = "n_old_mostly_inside")
  per_old <- p |> group_by(o = .data[[old_col]]) |>
    summarise(n_targets_ge20 = sum(share_of_old >= 0.2 & .data[[new_col]] != "unassigned"),
              share_unassigned = sum(share_of_old[.data[[new_col]] == "unassigned"]),
              targets = paste(sprintf("%s (%.0f%%)", .data[[new_col]][order(-share_of_old)], 100 * sort(share_of_old, decreasing = TRUE))[seq_len(min(4, n()))], collapse = "; "),
              .groups = "drop")
  best_new |>
    transmute(o = .data[[old_col]], best_new = .data[[new_col]], best_jaccard = jaccard, share_of_old_in_best = share_of_old) |>
    left_join(best_old_for_new, by = c("best_new" = "new_id")) |>
    left_join(holds_other, by = c("best_new" = "new_id")) |>
    right_join(per_old, by = "o") |>
    mutate(
      reciprocal_best = !is.na(best_old_for_new) & best_old_for_new == o,
      n_old_mostly_inside = coalesce(n_old_mostly_inside, 0L),
      class = case_when(
        reciprocal_best & best_jaccard >= 0.5 ~ "retained",
        share_of_old_in_best >= 0.5 & n_old_mostly_inside >= 2 ~ "merged",
        share_of_old_in_best >= 0.5 ~ "changed",
        n_targets_ge20 >= 2 ~ "split",
        TRUE ~ "dissolved"
      )
    ) |>
    select(o, class, best_new, best_jaccard, share_of_old_in_best, reciprocal_best, n_targets_ge20, share_unassigned, targets)
}

classify_new <- function(pairs, old_col, new_col) {
  p <- filter(pairs, .data[[new_col]] != "unassigned")
  p |> group_by(n = .data[[new_col]]) |>
    summarise(
      main_old = .data[[old_col]][which.max(share_of_new)],
      share_from_main_old = max(share_of_new),
      n_old_contributing_ge20 = sum(share_of_new >= 0.2 & .data[[old_col]] != "unassigned"),
      share_from_unassigned = sum(share_of_new[.data[[old_col]] == "unassigned"]),
      sources = paste(sprintf("%s (%.0f%%)", .data[[old_col]][order(-share_of_new)], 100 * sort(share_of_new, decreasing = TRUE))[seq_len(min(4, n()))], collapse = "; "),
      .groups = "drop") |>
    mutate(class = case_when(
      n_old_contributing_ge20 >= 2 ~ "merged_from_several",
      share_from_main_old >= 0.5 & main_old != "unassigned" ~ "continues_one_old",
      TRUE ~ "new_or_reassembled"
    ))
}

top_terms <- function(sig, ids) sig |> filter(ModuleID %in% ids) |> group_by(ID) |>
  summarise(Description = paste0(first(Description), " (", first(Ontology), ")"), p.adjust = min(p.adjust), .groups = "drop") |> arrange(p.adjust)

go_support <- function(old_sig, new_sig, old_ids, new_ids) {
  ot <- top_terms(old_sig, old_ids); nt <- top_terms(new_sig, new_ids)
  top <- head(ot, TOP_N)
  kept <- top$ID %in% nt$ID
  tibble(
    n_old_sig_terms = nrow(ot), n_new_sig_terms = nrow(nt),
    go_set_jaccard = if (nrow(ot) + nrow(nt)) length(intersect(ot$ID, nt$ID)) / length(union(ot$ID, nt$ID)) else NA_real_,
    old_top_terms = paste(head(top$Description, 5), collapse = "; "),
    new_top_terms = paste(head(nt$Description, 5), collapse = "; "),
    n_old_top_terms = nrow(top), n_old_top_still_significant = sum(kept),
    go_supports_old = nrow(top) > 0 && mean(kept) >= 0.5
  )
}

theme_terms <- function(sig, registry) {
  if (!nrow(sig)) return(tibble(ModuleID = character(), theme_id = character()))
  m <- map_go_terms_to_manuscript_themes(distinct(sig, ID, Description), registry)$assignments
  if (is.null(m) || !nrow(m)) return(tibble(ModuleID = character(), theme_id = character()))
  sig |> inner_join(distinct(m, ID = GO_ID, theme_id), by = "ID", relationship = "many-to-many") |> distinct(ModuleID, theme_id)
}

registry <- read_manuscript_go_theme_registry(THEME_REGISTRY)
old_mapping <- if (file.exists(OLD_THEME_MAPPING)) rd(OLD_THEME_MAPPING) else NULL
provisional_mapping <- list()
overview <- list()

for (ds in DATASETS) {
  cat("== ", ds, "\n", sep = "")
  out <- wgcna_dirs("map_wgcna_module_lineage", ds, create = TRUE)$tables
  old <- read_network(ds, "old"); new <- read_network(ds, "new")

  pairs <- overlap_pairs(old$assign, new$assign, "old_module", "new_module") |> arrange(old_module, desc(jaccard))
  write_csv(mutate(pairs, dataset = ds, .before = 1), file.path(out, "module_lineage_pairs.csv"))
  oc <- classify_old(pairs, "old_module", "new_module") |> rename(old_module = o)
  nc <- classify_new(pairs, "old_module", "new_module") |> rename(new_module = n)

  go_cmp <- bind_rows(lapply(seq_len(nrow(oc)), function(i) {
    bind_cols(tibble(old_module = oc$old_module[i], best_new = oc$best_new[i], class = oc$class[i]),
              go_support(old$sig, new$sig, oc$old_module[i], oc$best_new[i]))
  }))
  oc <- left_join(oc, select(go_cmp, old_module, go_set_jaccard, n_old_top_still_significant, n_old_top_terms, go_supports_old, old_top_terms, new_top_terms), by = "old_module")
  write_csv(mutate(oc, dataset = ds, .before = 1), file.path(out, "old_module_correspondence.csv"))
  write_csv(mutate(nc, dataset = ds, .before = 1), file.path(out, "new_module_origin.csv"))
  write_csv(mutate(go_cmp, dataset = ds, .before = 1), file.path(out, "module_go_comparison.csv"))

  # Supermodules on the proteins of their member modules.
  sm_assign <- function(net) { s <- setNames(net$smap$SupermoduleID, net$smap$ModuleID); x <- unname(s[net$assign]); x[is.na(x)] <- "unassigned"; setNames(x, names(net$assign)) }
  sm_pairs <- overlap_pairs(sm_assign(old), sm_assign(new), "old_supermodule", "new_supermodule")
  sc <- classify_old(sm_pairs, "old_supermodule", "new_supermodule") |> rename(old_supermodule = o)
  sm_members <- function(net, sm) net$smap$ModuleID[net$smap$SupermoduleID == sm]
  sm_go <- bind_rows(lapply(seq_len(nrow(sc)), function(i) {
    go_support(old$sig, new$sig, sm_members(old, sc$old_supermodule[i]), sm_members(new, sc$best_new[i]))
  }))
  sc <- bind_cols(sc, select(sm_go, go_set_jaccard, n_old_top_still_significant, n_old_top_terms, go_supports_old, old_top_terms, new_top_terms))
  write_csv(mutate(sc, dataset = ds, .before = 1), file.path(out, "supermodule_lineage.csv"))

  # Reviewed old labels: carried only for a retained, GO-supported entity.
  reviewed_file <- file.path(REVIEWED_REFERENCE, paste0(ds, ".csv"))
  reviewed <- if (file.exists(reviewed_file)) rd(reviewed_file) else NULL
  carry <- NULL
  if (!is.null(reviewed) && nrow(reviewed)) {
    look <- bind_rows(
      transmute(oc, level = "module", entity_id = old_module, class, best_new, best_jaccard, go_supports_old, n_old_top_still_significant, n_old_top_terms),
      transmute(sc, level = "supermodule", entity_id = old_supermodule, class, best_new, best_jaccard, go_supports_old, n_old_top_still_significant, n_old_top_terms))
    carry <- reviewed |>
      transmute(dataset, level, entity_id = as.character(entity_id), reviewed_biological_label, reviewed_short_label, label_confidence) |>
      left_join(look, by = c("level", "entity_id")) |>
      mutate(decision = if_else(class %in% "retained" & go_supports_old %in% TRUE, "carry_provisional", "provisionally_unannotated"),
             reason = case_when(
               decision == "carry_provisional" ~ sprintf("retained (Jaccard %.2f) and %d of %d top GO terms significant again", best_jaccard, n_old_top_still_significant, n_old_top_terms),
               is.na(class) ~ "old entity not found in the reference network",
               class != "retained" ~ paste0("composition ", class, "; needs review"),
               TRUE ~ sprintf("retained but only %d of %d top GO terms significant again; needs review", n_old_top_still_significant, n_old_top_terms))) |>
      rename(old_entity_id = entity_id, new_entity_id = best_new)
    write_csv(carry, file.path(out, "label_carry_decisions.csv"))
  }

  # Microglia: provisional registry for the new network from lineage-carried labels only.
  if (identical(ds, "microglia")) {
    member_map <- wgcna_normalize_current_member_map(mutate(new$smap, dataset = ds), ds)
    carried <- if (is.null(carry)) tibble() else filter(carry, decision == "carry_provisional")
    entities <- bind_rows(tibble(level = "module", entity_id = sort(unique(member_map$ModuleID))),
                          tibble(level = "supermodule", entity_id = sort(unique(member_map$SupermoduleID))))
    lk <- entities |>
      left_join(transmute(carried, level, entity_id = new_entity_id, bio = reviewed_biological_label, short = reviewed_short_label,
                          conf = label_confidence, why = paste0("Carried from ", old_entity_id, ": ", reason)),
                by = c("level", "entity_id"), relationship = "many-to-many") |>
      group_by(level, entity_id) |>
      summarise(bio = if (n_distinct(bio, na.rm = TRUE) == 1) first(na.omit(bio)) else NA_character_,
                short = if (n_distinct(short, na.rm = TRUE) == 1) first(na.omit(short)) else NA_character_,
                conf = first(na.omit(conf)), why = paste(na.omit(why), collapse = " | "), .groups = "drop") |>
      transmute(dataset = ds, level, entity_id,
                best_data_driven_label = coalesce(bio, "Provisionally unannotated"),
                raw_top_GO_label = NA_character_, microenvironment_context = NA_character_,
                final_plot_label = paste0(entity_id, " \u00b7 ", coalesce(short, "unannotated (for review)")),
                label_confidence = if_else(is.na(bio), "low", coalesce(conf, "low")),
                label_basis = if_else(is.na(bio), "no retained, GO-supported predecessor", "lineage carry with GO support"),
                label_rationale = if_else(nzchar(why), why, "Flagged for review: no retained, GO-supported predecessor label."))
    reg <- wgcna_provisional_registry_from_lookup(lk, member_map, ds,
             prepared_by = "tools/map_wgcna_module_lineage.R (lineage carry with GO support)")
    write_csv(reg, file.path(out, "provisional_label_registry.csv"), na = "")
  }

  # Theme-to-entity mappings: carried only for a retained entity with in-theme GO support.
  if (!is.null(old_mapping)) {
    om <- filter(old_mapping, dataset == ds)
    if (nrow(om)) {
      nt_mod <- theme_terms(new$sig_bp, registry)
      sm_size <- count(new$smap, SupermoduleID, name = "n_members")
      look <- bind_rows(transmute(oc, entity_level = "module", entity_id = old_module, class, best_new),
                        transmute(sc, entity_level = "supermodule", entity_id = old_supermodule, class, best_new))
      tm <- om |> mutate(entity_id = as.character(entity_id)) |> left_join(look, by = c("entity_level", "entity_id")) |>
        rowwise() |>
        mutate(new_members = list(if (entity_level == "module") best_new else new$smap$ModuleID[new$smap$SupermoduleID %in% best_new]),
               theme_supported_in_new = any(nt_mod$ModuleID %in% new_members & nt_mod$theme_id == theme_id),
               new_is_singleton_supermodule = entity_level == "supermodule" && isTRUE(sm_size$n_members[match(best_new, sm_size$SupermoduleID)] == 1L)) |>
        ungroup() |>
        mutate(carry = mapping_status == "approved" & class %in% "retained" & theme_supported_in_new & !new_is_singleton_supermodule,
               reason = case_when(carry ~ "retained and the new entity has significant GO BP terms in this theme",
                                  mapping_status != "approved" ~ paste0("old status ", mapping_status, "; only approved mappings are carried"),
                                  is.na(class) ~ "old entity not found",
                                  class != "retained" ~ paste0("composition ", class),
                                  new_is_singleton_supermodule ~ "new supermodule is a singleton alias",
                                  TRUE ~ "no significant GO BP term in this theme in the new entity")) |>
        select(-new_members)
      write_csv(rename(tm, old_entity_id = entity_id, new_entity_id = best_new), file.path(out, "theme_mapping_carry.csv"))
      provisional_mapping[[ds]] <- tm |> filter(carry) |>
        transmute(dataset, theme_id, entity_level,
                  mapping_rationale = paste0("Provisional, not reviewed: carried from ", entity_id, " (approved for the previous network); retained entity with in-theme GO support in the re-normalised network."),
                  entity_id = best_new,
                  mapping_status = "provisional_supported",
                  annotation_source = "tools/map_wgcna_module_lineage.R", confidence, mapping_role,
                  approved_for_manuscript_interpretation = FALSE) |>
        distinct(dataset, theme_id, entity_level, entity_id, .keep_all = TRUE)
    }
  }

  overview[[ds]] <- tibble(dataset = ds, n_old_modules = nrow(oc), n_new_modules = nrow(nc),
                           old_retained = sum(oc$class == "retained"), old_changed = sum(oc$class == "changed"),
                           old_merged = sum(oc$class == "merged"), old_split = sum(oc$class == "split"),
                           old_dissolved = sum(oc$class == "dissolved"),
                           retained_with_go_support = sum(oc$class == "retained" & oc$go_supports_old %in% TRUE),
                           n_old_supermodules = nrow(sc), n_new_supermodules = n_distinct(new$smap$SupermoduleID),
                           supermodules_retained = sum(sc$class == "retained"))
  print(as.data.frame(select(oc, old_module, class, best_new, best_jaccard, go_supports_old)), row.names = FALSE)
}

g <- wgcna_dirs("map_wgcna_module_lineage", "global", create = TRUE)$tables
write_csv(bind_rows(overview), file.path(g, "lineage_overview.csv"))
if (length(provisional_mapping)) write_csv(bind_rows(provisional_mapping), file.path(g, "provisional_theme_module_mapping.csv"))
print(as.data.frame(bind_rows(overview)), row.names = FALSE)

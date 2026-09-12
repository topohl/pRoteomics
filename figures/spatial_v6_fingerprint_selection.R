#!/usr/bin/env Rscript

# Part-21 step 2: the PHENOTYPE-BLIND feature-selection rule for the Figure-2
# direct spatial molecular fingerprint, plus the two candidate display forms.
#
# THE RULE (one sentence, reproducible):
#   Every displayed row is defined by a source that was fixed before any stress
#   contrast was computed - either an EXTERNAL published hippocampal spatial
#   signature (Kaulich 2025), an EXTERNAL curated cell-class reference panel
#   (config/marker_panels/), or a PRESPECIFIED CON-ONLY anatomical contrast -
#   and the value shown is the within-protein standardised mean CON abundance
#   in each spatial unit. No stress group, contrast, effect size, p-value or
#   FDR enters selection or scoring at any point.
#
# Two display forms are built so the more legible one can be chosen (brief
# section 12 asks for both):
#   A. individual named proteins  - concrete, nameable, but only as spatial as
#      the proteins happen to be
#   B. external signature scores  - smoother and directly interpretable, and it
#      is the only form that gives the microglia-enriched compartment a
#      prespecified spatial row, because NO CON-only anatomical contrast exists
#      for microglia anywhere in the repository.
#
# CLAIM DISCIPLINE: descriptive only. Reshapes and averages canonical values.
# No model, no test, no p-value, no new statistic.

source(file.path("R", "paths.R"))
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "spatial_grammar_utils.R"))
suppressPackageStartupMessages({ library(readr) })
Sys.setenv(PROTEOMICS_SCRIPT_ID = "figures/spatial_v6_fingerprint_selection.R")

OUT <- function(...) {
  d <- path_results("tables", "manuscript_candidates", "spatial_v6")
  dir_create(d)
  file.path(d, ...)
}
rd <- function(p) as.data.frame(readr::read_csv(p, show_col_types = FALSE,
                                                progress = FALSE, guess_max = Inf))

BASE <- OUT("spatial_v6_con_baseline_profile_long.csv")
if (!file.exists(BASE)) {
  stop("missing_required_input: run figures/spatial_v6_baseline_profile.R first: ",
       BASE, call. = FALSE)
}
KAU <- repo_path("results", "tables", "04_differential_expression_enrichment",
                 "control_spatial_identity_validation", "global",
                 "kaulich_signature_mapping.csv")
ANA <- repo_path("results", "tables", "04_differential_expression_enrichment",
                 "control_spatial_identity_validation", "global",
                 "anatomical_protein_contrasts.csv")
REF <- repo_path("config", "marker_panels", "wgcna_reference_marker_sets.csv")

args <- commandArgs(trailingOnly = TRUE)
if ("--dry-run" %in% args || is_dry_run()) {
  for (p in c(BASE, KAU, ANA, REF)) {
    message("[DRY-RUN ", if (file.exists(p)) "PASS" else "WARN", "] ", relative_to(p))
  }
  quit(save = "no", status = 0L)
}

base <- rd(BASE)
base$gene <- toupper(trimws(as.character(base$GeneSymbol)))

# =====================================================================
# FORM B: external prespecified signature scores
# =====================================================================

# --- B1: external hippocampal spatial signatures (Kaulich 2025) ----------
kau <- rd(KAU)
kau$gene <- toupper(trimws(as.character(kau$canonical_official_gene_symbol)))
kau <- kau[!is.na(kau$gene) & nzchar(kau$gene), , drop = FALSE]
k_sets <- unique(kau[, c("external_signature", "gene")])
k_sets$signature <- paste0("Kaulich ", k_sets$external_signature)
k_sets$signature_family <- ifelse(
  k_sets$external_signature %in% c("CA1", "CA2/3", "DG"),
  "external hippocampal subregion signature",
  "external hippocampal strata signature")
k_sets$source <- "Kaulich 2025 external published hippocampal proteome"
k_sets$selection_rule <- paste0(
  "all genes of the external published signature that map to a canonical ",
  "official mouse symbol (kaulich_signature_mapping.csv, mapping_status = mapped)")
k_sets$declared_allowed_use <- "external spatial validation reference"

# --- B2: external curated cell-class reference panels --------------------
ref <- rd(REF)
ref$gene <- toupper(trimws(as.character(ref$gene_symbol)))
KEEP_SETS <- c("reference_hippocampal_excitatory_neuron",
               "reference_inhibitory_interneuron",
               "reference_astrocyte",
               "reference_oligodendrocyte",
               "reference_microglia_pvm",
               "reference_vascular",
               "canonical_neuronal_synaptic_neuropil")
ref <- ref[ref$marker_set %in% KEEP_SETS & nzchar(ref$gene), , drop = FALSE]
PRETTY <- c(reference_hippocampal_excitatory_neuron = "Excitatory neuron",
            reference_inhibitory_interneuron = "Inhibitory interneuron",
            reference_astrocyte = "Astrocyte",
            reference_oligodendrocyte = "Oligodendrocyte",
            reference_microglia_pvm = "Microglia / PVM",
            reference_vascular = "Vascular",
            canonical_neuronal_synaptic_neuropil = "Synaptic neuropil")
r_sets <- unique(ref[, c("marker_set", "gene", "source_name", "use_for")])
r_sets$signature <- unname(PRETTY[r_sets$marker_set])
r_sets$signature_family <- "external curated cell-class reference panel"
r_sets$source <- paste0("config/marker_panels/wgcna_reference_marker_sets.csv (",
                        r_sets$source_name, ")")
r_sets$selection_rule <- paste0(
  "every gene of the prespecified marker_set, fixed in a config file before ",
  "any stress contrast was computed")
r_sets$declared_allowed_use <- r_sets$use_for

sig <- rbind(
  k_sets[, c("signature", "signature_family", "gene", "source", "selection_rule",
             "declared_allowed_use")],
  r_sets[, c("signature", "signature_family", "gene", "source", "selection_rule",
             "declared_allowed_use")])
sig <- unique(sig)

# score = mean within-protein z of the signature's measured members, per
# (dataset, spatial_unit). A signature with too few measured members in a
# compartment is reported as NA rather than drawn from noise.
MIN_MEMBERS <- 5L
j <- merge(base[, c("dataset", "spatial_unit", "gene", "con_z", "ProteinGroupID")],
           sig, by = "gene", allow.cartesian = TRUE)
key <- paste(j$signature, j$dataset, j$spatial_unit, sep = "\r")
agg <- do.call(rbind, lapply(split(seq_len(nrow(j)), key), function(ix) {
  z <- j[ix, , drop = FALSE]
  data.frame(signature = z$signature[1], signature_family = z$signature_family[1],
             dataset = z$dataset[1], spatial_unit = z$spatial_unit[1],
             n_measured_members = length(unique(z$ProteinGroupID)),
             score = mean(z$con_z, na.rm = TRUE),
             stringsAsFactors = FALSE)
}))
rownames(agg) <- NULL
agg$score[agg$n_measured_members < MIN_MEMBERS] <- NA_real_
agg <- sg_annotate(agg, "spatial_unit", "dataset")
write_csv_safe(agg, OUT("figure2_spatial_fingerprint_scores.csv"))

# =====================================================================
# FORM A: individual proteins from PRESPECIFIED CON-ONLY contrasts
# =====================================================================
ana <- rd(ANA)
ana$gene <- toupper(trimws(as.character(ana$official_gene_symbol)))
ana <- ana[!is.na(ana$gene) & nzchar(ana$gene) & is.finite(ana$adj.P.Val), ,
           drop = FALSE]
# the ENRICHED direction of each prespecified anatomical contrast
ana <- ana[ana$logFC > 0, , drop = FALSE]
TOP_N <- 2L
pick <- do.call(rbind, lapply(split(seq_len(nrow(ana)), ana$contrast), function(ix) {
  z <- ana[ix, , drop = FALSE]
  # deterministic: order by FDR then by |logFC| then by symbol
  z <- z[order(z$adj.P.Val, -abs(z$logFC), z$gene), , drop = FALSE]
  z <- z[!duplicated(z$gene), , drop = FALSE]
  z <- utils::head(z, TOP_N)
  data.frame(gene = z$gene, contrast = z$contrast, dataset = z$dataset,
             logFC = z$logFC, adj_P_Val = z$adj.P.Val, stringsAsFactors = FALSE)
}))
rownames(pick) <- NULL
pick$signature_family <- "prespecified CON-only anatomical contrast"
pick$source <- "control_spatial_identity_validation/anatomical_protein_contrasts.csv"
pick$selection_rule <- paste0(
  "the top ", TOP_N, " genes by BH-adjusted p among the ENRICHED side ",
  "(logFC > 0) of each prespecified CON-only anatomical contrast; the ",
  "contrast set was fixed in advance and is fitted on CON animals only")
pick$declared_allowed_use <- "control spatial identity validation"

prot <- merge(base[, c("dataset", "spatial_unit", "gene", "con_z", "con_mean_log2",
                       "ProteinGroupID")],
              unique(pick[, c("gene", "contrast", "signature_family", "source",
                              "selection_rule", "declared_allowed_use")]),
              by = "gene")
prot <- sg_annotate(prot, "spatial_unit", "dataset")
write_csv_safe(prot, OUT("figure2_spatial_fingerprint_proteins.csv"))

# =====================================================================
# the selection provenance table required by the brief
# =====================================================================
selA <- unique(data.frame(
  display_form = "A_individual_proteins",
  row_id = pick$gene,
  source = pick$source,
  selection_rule = pick$selection_rule,
  phenotype_blind = "YES",
  phenotype_blind_evidence = paste0(
    "contrast fitted on StressGroup == CON only ",
    "(R/control_spatial_identity_utils.R); no stress group, contrast, effect ",
    "or FDR enters selection"),
  dataset = pick$dataset,
  biological_interpretation = pick$contrast,
  selected_from_external_or_prespecified_source = "prespecified CON-only contrast",
  stringsAsFactors = FALSE))
selB <- unique(data.frame(
  display_form = "B_signature_scores",
  row_id = sig$signature,
  source = sig$source,
  selection_rule = sig$selection_rule,
  phenotype_blind = "YES",
  phenotype_blind_evidence = ifelse(
    grepl("^Kaulich", sig$signature),
    "membership comes from an external published dataset with no contact with these animals",
    "membership fixed in a repository config file, external curated panel"),
  dataset = NA_character_,
  biological_interpretation = sig$signature_family,
  selected_from_external_or_prespecified_source = "external",
  stringsAsFactors = FALSE))
selB <- selB[!duplicated(selB$row_id), , drop = FALSE]
sel <- rbind(selA, selB)
sel$excluded_source_note <- paste0(
  "protein_baseline_spatial_profile.csv was NOT used for selection: its scope ",
  "is is_sus_res_fdr_supported | is_wgcna_candidate, and 494 of its 505 rows ",
  "carry a candidate_reason containing 'SUS - RES'. Its VALUES are CON-only ",
  "and are reproduced exactly here, but its ROW MEMBERSHIP is not phenotype-blind.")
write_csv_safe(sel, OUT("figure2_spatial_fingerprint_selection.csv"))

# ------------------------------------------------------------------ report
cat("\n===== phenotype-blind spatial fingerprint selection =====\n")
cat("form A - individual proteins from prespecified CON-only contrasts\n")
cat("  contrasts:", length(unique(pick$contrast)), " genes:",
    length(unique(pick$gene)), "\n")
print(table(pick$dataset))
cat("\n  selected genes by contrast:\n")
for (cc in sort(unique(pick$contrast))) {
  cat(sprintf("    %-36s %s\n", cc,
              paste(pick$gene[pick$contrast == cc], collapse = ", ")))
}
cat("\nform B - external signature scores\n")
cat("  signatures:", length(unique(agg$signature)), "\n")
cv <- tapply(agg$score, agg$signature, function(z) sum(!is.na(z)))
mm <- tapply(agg$n_measured_members, agg$signature, max)
for (s in names(cv)) {
  cat(sprintf("    %-32s units scored %2d   max measured members %4d\n",
              s, cv[[s]], mm[[s]]))
}
cat("\n  coverage by compartment (units with a usable score):\n")
print(table(agg$dataset[!is.na(agg$score)]))
cat("\nNOTE: no CON-only anatomical contrast exists for microglia anywhere in\n")
cat("the repo, so form A gives the microglia block no prespecified row of its\n")
cat("own; form B does, via the external cell-class panels.\n")
cat("\nwritten:\n  ", relative_to(OUT("figure2_spatial_fingerprint_selection.csv")),
    "\n  ", relative_to(OUT("figure2_spatial_fingerprint_scores.csv")),
    "\n  ", relative_to(OUT("figure2_spatial_fingerprint_proteins.csv")), "\n")

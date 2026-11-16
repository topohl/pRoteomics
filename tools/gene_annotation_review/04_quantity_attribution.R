#!/usr/bin/env Rscript
# Step 4: how much of a protein group's quantity comes from gene-unique peptides, and does including the shared
# ones bias the group comparisons?
#   Rscript tools/gene_annotation_review/04_quantity_attribution.R --out <dir of steps 1-3> --pr-matrix <...> \
#     --pg-matrix <...> --run-map <csv> [--accessions P23198,P23492]
# --run-map: one row per run with sample_number, run_number, AnimalID, ExpGroup (1 CON, 2 RES, 3 SUS),
#   celltype_layer (compartment), region, metadata_layer and action ("analysed" rows are used).
# Default accessions: those step 3 found with both gene-unique and shared peptides.
# For each accession and run: MaxLFQ over all precursors, over gene-unique precursors only and over shared ones,
# compared with the delivered protein-group value. At animal level (runs centred within spatial unit, then averaged
# per animal), the difference log2(unique-only) - log2(all) per group shows the bias of each group contrast.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "review_helpers.R"))
a <- review_args(c("out", "pr-matrix", "pg-matrix", "run-map"))
pep <- fread(file.path(a$out, "peptide_gene_matches.csv"), colClasses = "character")
summ <- fread(file.path(a$out, "peptide_specificity_summary.csv"), colClasses = "character")
targets <- if (!is.null(a$accessions)) strsplit(a$accessions, ",")[[1]] else
  summ[as.integer(gene_unique_peptides) > 0 & as.integer(shared_peptides) > 0]$accession
if (!length(targets)) { cat("no accession with both gene-unique and shared peptides\n"); quit(status = 0, save = "no") }
pr <- fread(a[["pr-matrix"]], sep = "\t", quote = ""); pg <- fread(a[["pg-matrix"]], sep = "\t", quote = "")
runs <- setdiff(names(pr), pr_meta_cols)
stopifnot(setequal(runs, setdiff(names(pg), pg_meta_cols)))
rm <- fread(a[["run-map"]])[action == "analysed"]
rm[, key := paste0(sample_number, "|", run_number)]
info <- rm[match(run_key(runs), key), .(AnimalID, ExpGroup, compartment = celltype_layer, unit = paste(region, metadata_layer))]
if (anyNA(info$AnimalID)) stop(sum(is.na(info$AnimalID)), " runs are not in the run map", call. = FALSE)
grp <- c(`1` = "CON", `2` = "RES", `3` = "SUS")

runs_out <- list(); comp_out <- list(); animal_out <- list()
for (x in targets) {
  sub <- pr[vapply(split_ids(Protein.Ids), function(i) x %in% i, logical(1))]
  cls <- unique(pep[accession == x, .(Stripped.Sequence = peptide, unique_to_gene = gene_unique %in% c("TRUE", "true", TRUE))])
  sub <- merge(sub, cls, by = "Stripped.Sequence", all.x = TRUE)
  X <- log2(as.matrix(sub[, ..runs])); X[!is.finite(X)] <- NA
  pgv <- log2(as.numeric(unlist(pg[Protein.Group == x, ..runs]))); if (!length(pgv)) pgv <- rep(NA_real_, length(runs)); pgv[!is.finite(pgv)] <- NA
  u <- sub$unique_to_gene %in% TRUE; s <- sub$unique_to_gene %in% FALSE
  tot <- colSums(2^X, na.rm = TRUE); utot <- colSums(2^X[u, , drop = FALSE], na.rm = TRUE)
  d <- data.table(accession = x, run = run_key(runs), info, pg_value = pgv, maxlfq_all = maxlfq(X),
                  maxlfq_unique = maxlfq(X[u, , drop = FALSE]), maxlfq_shared = maxlfq(X[s, , drop = FALSE]),
                  unique_share_of_intensity = ifelse(tot > 0, utot / tot, NA_real_))
  d[, logratio := maxlfq_unique - maxlfq_all]
  runs_out[[x]] <- d
  for (cp in unique(d$compartment)) {
    z <- d[compartment == cp]
    r <- function(p, q) suppressWarnings(stats::cor(p, q, use = "pairwise.complete.obs"))
    comp_out[[paste(x, cp)]] <- data.table(accession = x, compartment = cp, precursors_unique = sum(u), precursors_shared = sum(s),
      runs_with_value = sum(!is.na(z$pg_value)), r_pg_vs_all = round(r(z$pg_value, z$maxlfq_all), 3), r_pg_vs_unique = round(r(z$pg_value, z$maxlfq_unique), 3),
      r_unique_vs_shared = round(r(z$maxlfq_unique, z$maxlfq_shared), 3), unique_share_median = round(stats::median(z$unique_share_of_intensity, na.rm = TRUE), 2))
    z <- z[is.finite(logratio)]
    z[, lr_c := logratio - mean(logratio), by = unit]
    an <- z[, .(lr = mean(lr_c)), by = .(AnimalID, ExpGroup)][, G := grp[as.character(ExpGroup)]]
    m <- an[, .(v = mean(lr)), by = G]; dv <- function(p, q) round(m[G == p]$v - m[G == q]$v, 3)
    p <- tryCatch(stats::anova(stats::lm(lr ~ G, data = an))["G", "Pr(>F)"], error = function(e) NA_real_)
    animal_out[[paste(x, cp)]] <- data.table(accession = x, compartment = cp, animals = nrow(an), group_effect_P = signif(p, 2),
      `SUS-CON bias (log2)` = dv("SUS", "CON"), `RES-CON bias (log2)` = dv("RES", "CON"), `SUS-RES bias (log2)` = dv("SUS", "RES"))
  }
}
fwrite(rbindlist(runs_out), file.path(a$out, "quantity_attribution_runs.csv"))
C <- rbindlist(comp_out); A <- rbindlist(animal_out)
fwrite(C, file.path(a$out, "quantity_attribution_summary.csv")); fwrite(A, file.path(a$out, "quantity_attribution_animal.csv"))
options(width = 200); print(merge(C, A, by = c("accession", "compartment")))

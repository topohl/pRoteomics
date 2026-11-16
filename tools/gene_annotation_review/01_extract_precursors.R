#!/usr/bin/env Rscript
# Step 1: every DIA-NN precursor that maps to one of the reviewed accessions.
#   Rscript tools/gene_annotation_review/01_extract_precursors.R --pr-matrix <report.pr_matrix.tsv> \
#     --accessions <csv with column 'accession'> --out <dir>
# Writes precursors.csv (one row per accession x precursor: protein group, all protein IDs DIA-NN matched the
# peptide to, Proteotypic flag, sequence, runs detected), precursor_summary.csv and comembers.csv (every accession
# DIA-NN listed together with the reviewed ones; step 2 fetches their records).
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "review_helpers.R"))
a <- review_args(c("pr-matrix", "accessions", "out"))
dir.create(a$out, recursive = TRUE, showWarnings = FALSE)
acc <- unique(trimws(fread(a$accessions, colClasses = "character")$accession))
acc <- acc[nzchar(acc)]
pr <- fread(a[["pr-matrix"]], sep = "\t", quote = "")
stopifnot(all(pr_meta_cols %in% names(pr)))
runs <- setdiff(names(pr), pr_meta_cols)
v <- as.matrix(pr[, ..runs])
pr[, n_runs_detected := rowSums(!is.na(v) & v > 0)]
ids <- split_ids(pr$Protein.Ids); grp <- split_ids(pr$Protein.Group)
touch <- which(vapply(seq_len(nrow(pr)), function(i) any(ids[[i]] %in% acc) || any(grp[[i]] %in% acc), logical(1)))
sub <- pr[touch, c(pr_meta_cols, "n_runs_detected"), with = FALSE]
ids_s <- ids[touch]; grp_s <- grp[touch]
rows <- rbindlist(lapply(acc, function(x) {
  k <- which(vapply(seq_along(ids_s), function(i) x %in% ids_s[[i]] || x %in% grp_s[[i]], logical(1)))
  if (!length(k)) return(data.table(accession = x))
  cbind(data.table(accession = x), sub[k])
}), fill = TRUE)
fwrite(rows, file.path(a$out, "precursors.csv"))
summ <- rows[, .(n_precursors = sum(!is.na(Precursor.Id)), n_peptides = uniqueN(Stripped.Sequence, na.rm = TRUE),
                 n_proteotypic_peptides = uniqueN(Stripped.Sequence[Proteotypic %in% c(1, "1", TRUE)], na.rm = TRUE),
                 protein_groups = paste(unique(na.omit(Protein.Group)), collapse = " | "),
                 other_ids = paste(setdiff(unique(unlist(split_ids(na.omit(Protein.Ids)))), accession[1]), collapse = ";"),
                 max_runs_detected = suppressWarnings(max(n_runs_detected, na.rm = TRUE)), runs_total = length(runs)), by = accession]
fwrite(summ, file.path(a$out, "precursor_summary.csv"))
com <- sort(unique(c(acc, unlist(split_ids(c(na.omit(rows$Protein.Ids), na.omit(rows$Protein.Group)))))))
fwrite(data.table(accession = com, reviewed = com %in% acc), file.path(a$out, "comembers.csv"))
cat(sprintf("%d accessions, %d precursors, %d runs; %d accessions without precursors; %d co-listed accessions\n",
            length(acc), sum(!is.na(rows$Precursor.Id)), length(runs), sum(summ$n_precursors == 0), sum(!com %in% acc)))

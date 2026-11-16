#!/usr/bin/env Rscript
# Step 3: can the observed peptides tell the candidate genes of an accession apart?
#   Rscript tools/gene_annotation_review/03_peptide_gene_specificity.R --out <dir of steps 1-2> [--fasta <search FASTA>[,<more>]]
# Every observed peptide (Stripped.Sequence, I read as L) is searched in every RefSeq protein of every candidate
# gene, and with --fasta in every entry of a whole-proteome FASTA (the search FASTA if available, otherwise the
# complete UniProt mouse proteome). A peptide is gene-unique if it occurs in exactly one candidate gene's proteins
# and in no FASTA entry of another gene (entries of the same gene under older names or aliases do not count);
# every other peptide is shared.
# Writes peptide_gene_matches.csv, candidate_gene_products.csv and peptide_specificity_summary.csv.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "review_helpers.R"))
a <- review_args("out")
pre <- fread(file.path(a$out, "precursors.csv"), colClasses = "character")
cand <- fread(file.path(a$out, "candidate_genes.csv"), colClasses = "character")
prot <- fread(file.path(a$out, "refseq_proteins.csv"), colClasses = "character")
up <- fread(file.path(a$out, "uniprot_entries.csv"), colClasses = "character")
fasta <- if (!is.null(a$fasta)) rbindlist(lapply(strsplit(a$fasta, ",")[[1]], read_fasta)) else NULL
if (!is.null(fasta)) fasta[, seq_il := il(sequence)]
prot[, seq_il := il(sequence)]

res <- rbindlist(lapply(unique(pre$accession), function(x) {
  peps <- unique(na.omit(pre[accession == x & nzchar(Stripped.Sequence)]$Stripped.Sequence))
  cg <- cand[accession == x]
  if (!length(peps)) return(data.table(accession = x))
  # names under which a FASTA entry still belongs to a candidate gene (old Riken names, aliases, UniProt names)
  u <- up[accession == x]
  same_gene <- tolower(unique(na.omit(c(cg$symbol, cg$mgi_symbol, unlist(split_ids(cg$synonyms)),
                                        u$gene_name, unlist(split_ids(u$gene_synonyms))))))
  rbindlist(lapply(peps, function(p) {
    q <- il(p)
    hit <- vapply(cg$gene_id, function(id) { s <- prot[gene_id == id]$seq_il; if (!length(s)) NA else any(grepl(q, s, fixed = TRUE)) }, logical(1))
    row <- data.table(accession = x, peptide = p,
      in_uniprot_sequence = grepl(q, il(up[accession == x]$sequence %||% ""), fixed = TRUE),
      genes_containing = paste(cg$symbol[which(hit)], collapse = ";"), n_genes_containing = sum(hit, na.rm = TRUE),
      genes_without_refseq_protein = paste(cg$symbol[is.na(hit)], collapse = ";"),
      max_runs_detected = max(as.integer(pre[accession == x & Stripped.Sequence == p]$n_runs_detected), na.rm = TRUE),
      diann_protein_ids = paste(unique(pre[accession == x & Stripped.Sequence == p]$Protein.Ids), collapse = " | "))
    if (!is.null(fasta)) {
      f <- fasta[grepl(q, seq_il, fixed = TRUE)]
      other <- f[!is.na(gene) & !tolower(gene) %in% same_gene]
      row[, `:=`(fasta_entries = paste(f$accession, collapse = ";"), fasta_genes = paste(sort(unique(na.omit(f$gene))), collapse = ";"),
                 fasta_other_gene_entries = nrow(other), fasta_other_genes = paste(sort(unique(other$gene)), collapse = ";"))]
    }
    # gene-unique: in exactly one candidate gene's proteins and, when a FASTA is given, in no other gene's entry
    row[, gene_unique := n_genes_containing == 1 & (is.null(fasta) || fasta_other_gene_entries == 0)]
    row
  }))
}), fill = TRUE)
fwrite(res, file.path(a$out, "peptide_gene_matches.csv"))

# what the candidate genes' proteins are: curated (NP_) or predicted (XP_), and identical to the UniProt sequence?
prod <- rbindlist(lapply(unique(cand$accession), function(x) {
  s <- up[accession == x]$sequence %||% ""
  cand[accession == x, {
    ps <- prot[gene_id == .BY$gene_id]
    .(refseq_proteins = paste(unique(ps$refseq), collapse = ";"), curated_NP = any(grepl("^NP_", ps$refseq)),
      predicted_XP_only = nrow(ps) > 0 && all(grepl("^XP_", ps$refseq)), identical_to_uniprot = nzchar(s) && any(ps$sequence == s))
  }, by = .(accession, gene_id, symbol, type)]
}), fill = TRUE)
fwrite(prod, file.path(a$out, "candidate_gene_products.csv"))

summ <- res[!is.na(peptide), {
  u <- gene_unique %in% TRUE
  out <- .(n_peptides = .N, n_candidate_genes = nrow(cand[accession == .BY$accession]),
    gene_unique_peptides = sum(u), genes_with_unique_peptides = paste(sort(unique(genes_containing[u])), collapse = ";"),
    shared_peptides = sum(!u & n_genes_containing >= 1), peptides_in_no_candidate_protein = sum(n_genes_containing == 0))
  if ("fasta_other_gene_entries" %in% names(.SD)) out <- c(out, .(peptides_in_other_gene_fasta_entries = sum(fasta_other_gene_entries > 0)))
  out
}, by = accession]
summ[, verdict := fcase(
  n_candidate_genes == 0, "no candidate gene (no NCBI Gene ID)",
  n_candidate_genes == 1 & peptides_in_no_candidate_protein == 0, "one candidate gene; every peptide in its proteins",
  gene_unique_peptides > 0 & shared_peptides > 0, "gene-unique and shared peptides: identity supported, quantity mixed",
  gene_unique_peptides > 0, "gene-unique peptides only",
  default = "all peptides shared: gene of origin not resolvable")]
if ("peptides_in_other_gene_fasta_entries" %in% names(summ))
  summ[n_candidate_genes == 1 & peptides_in_other_gene_fasta_entries > 0,
       verdict := "one candidate gene, but some peptides also occur in another gene's FASTA entry"]
fwrite(summ, file.path(a$out, "peptide_specificity_summary.csv"))
print(summ[, .N, by = verdict])

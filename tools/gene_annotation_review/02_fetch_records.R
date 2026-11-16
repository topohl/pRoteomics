#!/usr/bin/env Rscript
# Step 2: public records for the accessions and every candidate gene (public identifiers only are sent).
#   Rscript tools/gene_annotation_review/02_fetch_records.R --out <dir of step 1> [--cache <dir>] [--extra-gene-ids <csv>]
# Candidate genes of an accession = the Gene IDs UniProt links to it, plus the Entrez IDs the installed
# org.Mm.eg.db links to it (the source the pipeline's gene annotation uses), plus optional extra IDs
# (--extra-gene-ids: csv with columns accession, gene_id). For each candidate gene: NCBI Gene (Datasets v2),
# the MGI marker (Alliance API) and all RefSeq proteins (E-utilities). Responses are cached in --cache (default
# <out>/cache) and every URL is logged with its retrieval time in <out>/fetch_log.csv; delete the cache to refresh.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "review_helpers.R"))
suppressPackageStartupMessages({ library(AnnotationDbi); library(org.Mm.eg.db) })
a <- review_args("out")
cache <- a$cache %||% file.path(a$out, "cache")
logf <- file.path(a$out, "fetch_log.csv")
com <- fread(file.path(a$out, "comembers.csv"), colClasses = "character")
rev_acc <- com[reviewed == "TRUE"]$accession

# UniProt entries
for (x in com$accession) fetch_cached(sprintf("https://rest.uniprot.org/uniprotkb/%s.json", x), file.path(cache, "uniprot", paste0(x, ".json")), logf)
xr <- function(e, db, prop = NULL) unlist(lapply(e$uniProtKBCrossReferences, function(r) if (identical(r$database, db)) {
  if (is.null(prop)) r$id else paste0(r$id, "=", unlist(lapply(r$properties, function(p) if (p$key == prop) p$value)))
}))
up <- rbindlist(lapply(com$accession, function(x) {
  f <- file.path(cache, "uniprot", paste0(x, ".json"))
  if (!file.exists(f) || !file.size(f)) return(data.table(accession = x, status = "not retrieved"))
  e <- fromJSON(f, simplifyVector = FALSE); gn <- e$genes[[1]] %||% list()
  data.table(accession = x, status = "ok", primary_accession = e$primaryAccession, entry_name = e$uniProtkbId, entry_type = e$entryType,
    secondary_accessions = paste(unlist(e$secondaryAccessions), collapse = ";"), taxon = as.character(e$organism$taxonId),
    gene_name = gn$geneName$value %||% NA_character_, gene_synonyms = paste(unlist(lapply(gn$synonyms, `[[`, "value")), collapse = ";"),
    geneid_xrefs = paste(xr(e, "GeneID"), collapse = ";"), mgi_xrefs = paste(xr(e, "MGI", "GeneName"), collapse = ";"),
    protein_name = e$proteinDescription$recommendedName$fullName$value %||% e$proteinDescription$submissionNames[[1]]$fullName$value %||% NA_character_,
    seq_version = as.character(e$entryAudit$sequenceVersion), sequence = e$sequence$value)
}), fill = TRUE)
fwrite(up, file.path(a$out, "uniprot_entries.csv"))

# candidate genes
orgdb <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db, keys = intersect(rev_acc, keys(org.Mm.eg.db, "UNIPROT")),
                                                keytype = "UNIPROT", columns = "ENTREZID"))
extra <- if (!is.null(a[["extra-gene-ids"]])) fread(a[["extra-gene-ids"]], colClasses = "character") else data.table(accession = character(), gene_id = character())
cand <- unique(rbind(
  up[accession %in% rev_acc & nzchar(geneid_xrefs), .(gene_id = unlist(split_ids(geneid_xrefs)), source = "uniprot"), by = accession],
  as.data.table(orgdb)[!is.na(ENTREZID), .(accession = UNIPROT, gene_id = ENTREZID, source = "org.Mm.eg.db")],
  extra[, .(accession, gene_id, source = "extra")]))
cand <- cand[, .(sources = paste(sort(unique(source)), collapse = ";")), by = .(accession, gene_id)]
gids <- sort(unique(cand$gene_id))
for (ch in split(gids, ceiling(seq_along(gids) / 20)))
  fetch_cached(sprintf("https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/%s", paste(ch, collapse = ",")), file.path(cache, "ncbi", sprintf("genes_%s.json", paste(ch, collapse = "_") |> substr(1, 60))), logf)
ncbi <- rbindlist(lapply(list.files(file.path(cache, "ncbi"), pattern = "^genes_.*json$", full.names = TRUE), function(f)
  rbindlist(lapply(fromJSON(f, simplifyVector = FALSE)$reports, function(r) { g <- r$gene; data.table(gene_id = g$gene_id, symbol = g$symbol,
    type = g$type %||% NA_character_, mgi = g$nomenclature_authority$identifier %||% NA_character_, synonyms = paste(unlist(g$synonyms), collapse = ";")) }), fill = TRUE)), fill = TRUE)
ncbi <- unique(ncbi[gene_id %in% gids])
for (m in sort(unique(na.omit(ncbi$mgi)))) fetch_cached(sprintf("https://www.alliancegenome.org/api/gene/%s", m), file.path(cache, "mgi", paste0(sub(":", "_", m), ".json")), logf, pause = 0.2)
mgi <- rbindlist(lapply(sort(unique(na.omit(ncbi$mgi))), function(m) {
  f <- file.path(cache, "mgi", paste0(sub(":", "_", m), ".json")); if (!file.exists(f) || !file.size(f)) return(NULL)
  g <- fromJSON(f, simplifyVector = FALSE)$gene
  data.table(mgi = m, mgi_symbol = g$geneSymbol$displayText %||% NA_character_, mgi_type = g$geneType$name %||% NA_character_)
}), fill = TRUE)
cand <- merge(merge(cand, ncbi, by = "gene_id", all.x = TRUE), mgi, by = "mgi", all.x = TRUE)
setcolorder(cand, c("accession", "gene_id", "symbol", "type", "mgi", "mgi_symbol", "mgi_type", "sources"))
fwrite(cand[order(accession, symbol)], file.path(a$out, "candidate_genes.csv"))

# RefSeq proteins of every candidate gene
for (id in gids) fetch_cached(sprintf("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=gene&db=protein&linkname=gene_protein_refseq&id=%s&retmode=json", id),
                              file.path(cache, "refseq", sprintf("elink_%s.json", id)), logf)
links <- rbindlist(lapply(gids, function(id) {
  j <- fromJSON(file.path(cache, "refseq", sprintf("elink_%s.json", id)), simplifyVector = FALSE)
  ls <- j$linksets[[1]]$linksetdbs
  data.table(gene_id = id, protein_uid = as.character(if (length(ls)) unlist(ls[[1]]$links) else character()))
}), fill = TRUE)
puids <- unique(na.omit(links$protein_uid))
for (ch in split(puids, ceiling(seq_along(puids) / 100))) {
  key <- digest::digest(paste(ch, collapse = ","), algo = "md5")
  fetch_cached(sprintf("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=%s&rettype=fasta&retmode=text", paste(ch, collapse = ",")), file.path(cache, "refseq", sprintf("fasta_%s.fa", key)), logf)
  fetch_cached(sprintf("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=protein&id=%s&retmode=json", paste(ch, collapse = ",")), file.path(cache, "refseq", sprintf("esum_%s.json", key)), logf)
}
esum <- rbindlist(lapply(list.files(file.path(cache, "refseq"), pattern = "^esum_.*json$", full.names = TRUE), function(f) {
  r <- fromJSON(f, simplifyVector = FALSE)$result
  rbindlist(lapply(unlist(r$uids), function(i) data.table(protein_uid = i, refseq = r[[i]]$accessionversion %||% NA_character_)))
}))
fa <- rbindlist(lapply(list.files(file.path(cache, "refseq"), pattern = "^fasta_.*fa$", full.names = TRUE), read_fasta))
fa[, refseq := sub("^(\\S+).*", "\\1", header)]
prot <- unique(merge(merge(links, unique(esum), by = "protein_uid"), unique(fa[, .(refseq, sequence)]), by = "refseq"))
fwrite(prot[, .(gene_id, refseq, length = nchar(sequence), sequence)], file.path(a$out, "refseq_proteins.csv"))
cat(sprintf("%d UniProt entries; %d candidate genes for %d accessions; %d RefSeq proteins\n",
            sum(up$status == "ok"), length(gids), uniqueN(cand$accession), uniqueN(prot$refseq)))

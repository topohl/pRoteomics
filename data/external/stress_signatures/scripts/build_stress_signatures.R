#!/usr/bin/env Rscript
# Build the published hippocampal stress signatures from the raw downloads.
#
#   Rscript scripts/build_stress_signatures.R            # all signatures
#   Rscript scripts/build_stress_signatures.R Weber2019_microglia Laine2018_vHPC   # selected ids
#
# Input : raw files in <package>/raw/<ACCESSION>/ (or the folder in env var STRESS_SIG_RAW),
#         fetched by scripts/download_stress_signatures.R.
# Output: <package>/signatures/<id>.csv  (column 1 = gene, then one column per contrast; positive = higher in
#         the first group), <package>/signature_info.csv (one row per signature).
# Every processing choice (labels, covariates, filters, weights, symbol mapping) is the one fixed by the checked
# loaders of 2026-10-29; do not change it without re-checking the signature.
# Packages: limma, data.table, readxl, babelgene, AnnotationDbi, org.Mm.eg.db, digest;
#           Henningsen2012_CA (rat) also needs org.Rn.eg.db and org.Hs.eg.db.

suppressPackageStartupMessages({
  library(limma); library(data.table); library(AnnotationDbi); library(org.Mm.eg.db)
})

args_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (!length(args_file)) stop("run with Rscript scripts/build_stress_signatures.R")
script_dir <- dirname(normalizePath(args_file))
root <- normalizePath(file.path(script_dir, ".."))
RAW <- Sys.getenv("STRESS_SIG_RAW", file.path(root, "raw"))
OUT <- file.path(root, "signatures")
dir.create(OUT, showWarnings = FALSE)

raw_file <- function(acc, file) {
  p <- file.path(RAW, acc, file)
  if (!file.exists(p)) stop("missing raw file ", p, " - run scripts/download_stress_signatures.R first")
  p
}

####################################################################################################
## Shared helpers
####################################################################################################

CONTRASTS <- c(SUS_vs_CON = "SUS-CON", RES_vs_CON = "RES-CON", SUS_vs_RES = "SUS-RES",
               stressed_vs_CON = "(SUS+RES)/2-CON")

# median-of-ratios size factors (DESeq2 style) on the genes counted in every sample
mor_sf <- function(counts) {
  pos <- counts[apply(counts > 0, 1, all), , drop = FALSE]
  lg <- rowMeans(log(pos))
  sf <- apply(pos, 2, function(x) exp(stats::median(log(x) - lg)))
  sf / exp(mean(log(sf)))
}

# limma-voom on counts: keep genes with CPM > min_cpm in >= min_samples samples; the voom library size is
# sf x geometric-mean depth (median-of-ratios factors already contain depth); robust eBayes.
voom_mor <- function(counts, design, contrasts, min_samples, min_cpm = 1) {
  counts <- as.matrix(counts); storage.mode(counts) <- "double"
  cpm <- t(t(counts) / colSums(counts)) * 1e6
  counts <- counts[rowSums(cpm > min_cpm) >= min_samples, , drop = FALSE]
  sf <- mor_sf(counts)
  v <- voom(counts, design, lib.size = sf * exp(mean(log(colSums(counts)))))
  cm <- makeContrasts(contrasts = contrasts, levels = design); colnames(cm) <- names(contrasts)
  eb <- eBayes(contrasts.fit(lmFit(v, design), cm), robust = TRUE)
  list(t = eb$t, ave = rowMeans(v$E))
}

# limma-trend on already normalised log-scale data, robust eBayes
trend_t <- function(logx, design, contrasts) {
  cm <- makeContrasts(contrasts = contrasts, levels = design); colnames(cm) <- names(contrasts)
  eBayes(contrasts.fit(lmFit(logx, design), cm), trend = TRUE, robust = TRUE)$t
}

# one row per symbol: the row with the highest mean expression ('ave')
collapse_by_symbol <- function(stat, sym, ave) {
  ok <- !is.na(sym) & nzchar(sym)
  o <- order(sym[ok], -ave[ok])
  st <- stat[ok, , drop = FALSE][o, , drop = FALSE]; s <- sym[ok][o]
  st <- st[!duplicated(s), , drop = FALSE]; rownames(st) <- s[!duplicated(s)]
  st
}

# Ensembl (version stripped) -> MGI symbol, first org.Mm.eg.db hit
ens_symbol <- function(ids) {
  unname(suppressMessages(mapIds(org.Mm.eg.db, keys = sub("[.].*$", "", ids), column = "SYMBOL",
                                 keytype = "ENSEMBL", multiVals = "first")))
}

# Ensembl -> MGI symbol with a deterministic choice for one-to-many IDs (prefer the gene for which the ID is its only
# Ensembl ID, then non-placeholder names). 'annot' (named by Ensembl ID) = symbols of the counting annotation: an ID
# reassigned since then (org.Mm.eg.db symbol differs, is not an alias of it, and the annotation symbol is a current
# MGI symbol not otherwise used) keeps the annotation symbol.
ens_symbol_unique <- function(ids, annot = NULL) {
  ids0 <- sub("[.].*$", "", ids)
  e2s <- unique(suppressMessages(AnnotationDbi::select(org.Mm.eg.db, keys = keys(org.Mm.eg.db, "ENSEMBL"),
                                                       keytype = "ENSEMBL", columns = "SYMBOL")))
  e2s <- e2s[!is.na(e2s$SYMBOL), ]
  ne <- table(e2s$SYMBOL); cand <- split(e2s$SYMBOL, e2s$ENSEMBL)
  pick <- function(c0) {
    if (is.null(c0)) return(NA_character_)
    if (length(c0) == 1) return(c0)
    c0[order(as.integer(ne[c0]), grepl("^Gm[0-9]+$|Rik$|^LOC[0-9]+$", c0))][1]
  }
  out <- vapply(ids0, function(e) pick(cand[[e]]), "")
  if (!is.null(annot)) {
    a <- unname(annot[ids0])
    al <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db, keys = unique(out[!is.na(out)]), keytype = "SYMBOL",
                                                 columns = "ALIAS"))
    alias_ok <- !is.na(out) & !is.na(a) & paste(out, a, sep = "\r") %in% paste(al$SYMBOL, al$ALIAS, sep = "\r")
    re <- !is.na(out) & !is.na(a) & a != out & !alias_ok & a %in% keys(org.Mm.eg.db, "SYMBOL") & !(a %in% out)
    out[re] <- a[re]
  }
  unname(out)
}

# Entrez -> MGI symbol, first org.Mm.eg.db hit
entrez_symbol <- function(eg) {
  eg <- as.character(eg); ok <- !is.na(eg) & nzchar(eg)
  out <- rep(NA_character_, length(eg))
  if (any(ok)) out[ok] <- unname(suppressMessages(mapIds(org.Mm.eg.db, keys = eg[ok], keytype = "ENTREZID", column = "SYMBOL",
                                                         multiVals = "first")))
  out
}

# Symbols from an older annotation -> current MGI symbol: official symbols kept; others replaced when the alias
# maps to exactly one official symbol. drop_clash: an alias may not take the name of an official symbol already in
# the input. keep_unmapped: leave unmapped non-official symbols as they are (else NA = dropped).
mgi_symbol <- function(sym, keep_unmapped = FALSE, drop_clash = TRUE) {
  sym <- as.character(sym)
  off <- keys(org.Mm.eg.db, keytype = "SYMBOL")
  isoff <- sym %in% off
  cand <- intersect(unique(sym[!isoff & !is.na(sym) & nzchar(sym)]), keys(org.Mm.eg.db, keytype = "ALIAS"))
  m <- character()
  if (length(cand)) {
    al <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db, keys = cand, keytype = "ALIAS", columns = "SYMBOL"))
    al <- al[!is.na(al$SYMBOL), ]
    n <- table(al$ALIAS); al <- al[al$ALIAS %in% names(n)[n == 1], ]
    m <- setNames(al$SYMBOL, al$ALIAS)
    if (drop_clash) m <- m[!(m %in% sym[isoff])]
  }
  out <- if (keep_unmapped) sym else ifelse(isoff, sym, NA_character_)
  hit <- !isoff & !is.na(sym) & sym %in% names(m)
  out[hit] <- unname(m[sym[hit]])
  out
}

# rat symbols -> mouse symbols through human, babelgene orthologs, 1:1 at both steps, pseudogene hits dropped
rat_to_mouse <- function(rat) {
  r2h <- babelgene::orthologs(unique(rat[!is.na(rat)]), species = "rat", human = FALSE, top = FALSE)
  gt_h <- suppressMessages(mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = unique(as.character(r2h$human_entrez)),
                                  column = "GENETYPE", keytype = "ENTREZID"))
  r2h <- r2h[!(gt_h[as.character(r2h$human_entrez)] %in% "pseudo"), ]
  h_per_r <- tapply(r2h$human_symbol, r2h$symbol, function(z) length(unique(z)))
  r_per_h <- tapply(r2h$symbol, r2h$human_symbol, function(z) length(unique(z)))
  r2h1 <- unique(r2h[r2h$symbol %in% names(h_per_r)[h_per_r == 1] &
                     r2h$human_symbol %in% names(r_per_h)[r_per_h == 1], c("symbol", "human_symbol")])
  h2m <- babelgene::orthologs(unique(r2h1$human_symbol), species = "mouse", human = TRUE, top = FALSE)
  gt_m <- suppressMessages(mapIds(org.Mm.eg.db, keys = unique(as.character(h2m$entrez)), column = "GENETYPE",
                                  keytype = "ENTREZID"))
  h2m <- h2m[!(gt_m[as.character(h2m$entrez)] %in% "pseudo") & !grepl("-ps[0-9]*$", h2m$symbol), ]
  m_per_h <- tapply(h2m$symbol, h2m$human_symbol, function(z) length(unique(z)))
  h_per_m <- tapply(h2m$human_symbol, h2m$symbol, function(z) length(unique(z)))
  h2m1 <- unique(h2m[h2m$human_symbol %in% names(m_per_h)[m_per_h == 1] &
                     h2m$symbol %in% names(h_per_m)[h_per_m == 1], c("human_symbol", "symbol")])
  h2m1$symbol[match(r2h1$human_symbol[match(rat, r2h1$symbol)], h2m1$human_symbol)]
}

# GEO series matrix: sample annotation (one column per characteristics line, named by its key in the first sample;
# a repeated key gets a suffix), the raw header lines, and optionally the data table
read_series_matrix <- function(path, data = FALSE) {
  L <- readLines(gzfile(path))
  get <- function(tag) lapply(L[startsWith(L, tag)], function(z) gsub('"', "", strsplit(z, "\t")[[1]][-1]))
  ann <- data.frame(gsm = get("!Sample_geo_accession")[[1]], title = get("!Sample_title")[[1]],
                    source = get("!Sample_source_name_ch1")[[1]], stringsAsFactors = FALSE)
  for (v in get("!Sample_characteristics_ch1")) {
    k <- make.unique(c(names(ann), sub(":.*", "", v[1])))[ncol(ann) + 1]
    ann[[k]] <- sub("^[^:]*: ", "", v)
  }
  dsc <- get("!Sample_description"); if (length(dsc)) ann$description <- dsc[[1]]
  sup <- get("!Sample_supplementary_file_1"); if (length(sup)) ann$suppl <- sup[[1]]
  for (v in get("!Sample_relation")) if (startsWith(v[1], "SRA:")) ann$sra <- v
  x <- NULL
  if (data) x <- read.delim(gzfile(path), comment.char = "!", check.names = FALSE)
  list(ann = ann, lines = L, data = x)
}

# Minimal BIFF8 (.xls) reader for a worksheet of NUMBER and shared-string cells (used for GSE81672, which readxl's
# libxls cannot open; the result equals Excel's own CSV export cell for cell).
read_xls_biff8 <- function(path) {
  con <- if (grepl("[.]gz$", path)) gzfile(path, "rb") else file(path, "rb")
  chunks <- list(); repeat { b <- readBin(con, "raw", 2^24); if (!length(b)) break; chunks[[length(chunks) + 1]] <- b }
  close(con); ole <- unlist(chunks, use.names = FALSE); rm(chunks)
  le32 <- function(r) readBin(r, "integer", n = length(r) %/% 4, size = 4, endian = "little")
  stopifnot(identical(ole[1:8], as.raw(c(0xD0, 0xCF, 0x11, 0xE0, 0xA1, 0xB1, 0x1A, 0xE1))))   # OLE2 compound file
  hd <- le32(ole[1:512]); ss <- 2L^(as.integer(ole[31]) + 256L * as.integer(ole[32]))
  sec <- function(s) {                                    # bytes of sectors s, joined over runs of adjacent sectors
    brk <- c(0L, which(diff(s) != 1L), length(s))
    unlist(lapply(seq_len(length(brk) - 1L), function(i) {
      a <- s[brk[i] + 1L]; z <- s[brk[i + 1L]]; ole[(512L + a * ss + 1L):(512L + (z + 1L) * ss)] }), use.names = FALSE)
  }
  difat <- hd[20:128]; nxt <- hd[18]
  for (i in seq_len(hd[19])) { v <- le32(ole[512L + nxt * ss + seq_len(ss)]); difat <- c(difat, v[-length(v)]); nxt <- v[length(v)] }
  fat <- le32(sec(difat[difat >= 0][seq_len(hd[12])]))
  chain <- function(s) { out <- integer(); while (s >= 0) { out <- c(out, s); s <- fat[s + 1L] }; out }
  dir <- sec(chain(hd[13]))
  ents <- lapply(seq_len(length(dir) %/% 128) - 1L, function(i) dir[i * 128 + 1:128])
  nm <- vapply(ents, function(e) { n <- as.integer(e[65]) + 256L * as.integer(e[66])
    if (n < 2) "" else rawToChar(e[seq(1, n - 2, by = 2)]) }, "")
  e <- ents[[which(nm %in% c("Workbook", "Book"))[1]]]
  len <- le32(e[121:124]); stopifnot(len >= 4096)        # large stream: stored in regular sectors
  wb <- sec(chain(le32(e[117:120])))[seq_len(len)]; rm(ole)
  b <- as.integer(wb); N <- length(b)                       # walk the BIFF records
  ty <- integer(N %/% 8); ln <- ty; of <- ty; k <- 0L; p <- 1L
  while (p + 3L <= N) {
    k <- k + 1L; ty[k] <- b[p] + 256L * b[p + 1L]; ln[k] <- b[p + 2L] + 256L * b[p + 3L]; of[k] <- p + 4L
    p <- p + 4L + ln[k]
  }
  ty <- ty[1:k]; ln <- ln[1:k]; of <- of[1:k]
  u16 <- function(o) b[o] + 256L * b[o + 1L]
  i_sst <- which(ty == 0x00FC)[1]; j <- i_sst + 1L       # shared string table + CONTINUE records
  while (ty[j] == 0x003C) j <- j + 1L
  segs <- lapply(i_sst:(j - 1L), function(i) wb[of[i] - 1L + seq_len(ln[i])])
  n_str <- le32(segs[[1]][5:8]); segs[[1]] <- segs[[1]][-(1:8)]
  si <- 1L; pos <- 1L
  take <- function(n) {
    out <- raw(0)
    while (n > 0) {
      av <- length(segs[[si]]) - pos + 1L
      if (av <= 0) { si <<- si + 1L; pos <<- 1L; next }
      m <- min(n, av); out <- c(out, segs[[si]][pos:(pos + m - 1L)]); pos <<- pos + m; n <- n - m
    }
    out
  }
  sst <- character(n_str)
  for (s in seq_len(n_str)) {
    h <- take(3); cch <- as.integer(h[1]) + 256L * as.integer(h[2]); fl <- as.integer(h[3])
    high <- bitwAnd(fl, 1L) > 0
    nrun <- if (bitwAnd(fl, 8L) > 0) { r <- take(2); as.integer(r[1]) + 256L * as.integer(r[2]) } else 0L
    next_ext <- if (bitwAnd(fl, 4L) > 0) le32(take(4)) else 0L
    txt <- ""; left <- cch
    while (left > 0) {                                       # characters may continue in the next record
      av <- length(segs[[si]]) - pos + 1L
      if (av <= 0) { si <- si + 1L; high <- bitwAnd(as.integer(segs[[si]][1]), 1L) > 0; pos <- 2L; next }
      bpc <- if (high) 2L else 1L; m <- min(left, av %/% bpc)
      by <- segs[[si]][pos:(pos + m * bpc - 1L)]
      txt <- paste0(txt, if (high) intToUtf8(readBin(by, "integer", n = m, size = 2, signed = FALSE, endian = "little"))
                         else rawToChar(by))
      pos <- pos + m * bpc; left <- left - m
    }
    if (nrun + next_ext > 0) take(4L * nrun + next_ext)
    sst[s] <- txt
  }
  bof <- which(ty == 0x0809); eof <- which(ty == 0x000A)   # first worksheet = second BOF .. next EOF
  ws <- bof[2]:eof[eof > bof[2]][1]
  bad <- intersect(ty[ws], c(0x027E, 0x00BD, 0x0006, 0x0204, 0x0205, 0x00D6))   # RK, MULRK, FORMULA, LABEL, BOOLERR, RSTRING
  if (length(bad)) stop("unsupported BIFF cell records: ", paste(sprintf("0x%04X", bad), collapse = " "))
  o <- of[ws][ty[ws] == 0x0203]                             # NUMBER
  num <- data.frame(r = u16(o), c = u16(o + 2L),
                    v = readBin(wb[as.vector(outer(0:7, o + 6L, "+"))], "double", n = length(o), size = 8, endian = "little"))
  o <- of[ws][ty[ws] == 0x00FD]                             # LABELSST
  lab <- data.frame(r = u16(o), c = u16(o + 2L), v = sst[le32(wb[as.vector(outer(0:3, o + 6L, "+"))]) + 1L])
  list(num = num, lab = lab)
}

####################################################################################################
## Datasets (one function each; each returns a list of signatures: id, stat matrix, info)
####################################################################################################

sig <- function(id, stat, ...) list(id = id, stat = stat, info = list(id = id, ...))

## --- GSE11211: Schmidt 2010, dorsal CA1 and DG, SUS vs RES pools --------------------------------
build_schmidt2010 <- function() {
  sm <- raw_file("GSE11211", "GSE11211_series_matrix.txt.gz")
  hdr <- readLines(gzfile(sm))
  fld <- function(key) { l <- grep(paste0("^!", key, "\t"), hdr, value = TRUE)[1]; gsub('"', "", strsplit(l, "\t")[[1]][-1]) }
  ann <- data.frame(gsm = fld("Sample_geo_accession"), title = fld("Sample_title"), source = fld("Sample_source_name_ch1"),
                    stringsAsFactors = FALSE)
  ann$group <- ifelse(grepl("^Susceptible phenotype", ann$title) & grepl("very high basal corticosterone", ann$source), "SUS",
               ifelse(grepl("^Resistant phenotype", ann$title) & grepl("very low basal corticosterone", ann$source), "RES", NA))
  ann$region <- ifelse(grepl(", CA1,", ann$title) & grepl("subregion CA1", ann$source), "CA1",
                ifelse(grepl(", DG,", ann$title) & grepl("subregion dentate gyrus", ann$source), "DG", NA))
  stopifnot(!anyNA(ann$group), !anyNA(ann$region), all(table(ann$region, ann$group) == 4))
  x <- read.delim(gzfile(sm), comment.char = "!", check.names = FALSE)
  rownames(x) <- as.character(x$ID_REF); x <- as.matrix(x[, ann$gsm])
  # noise floor: per-array mode and left-half SD of the background peak; shift peak to 0, floor at 3 SD, log2
  bg <- t(apply(x, 2, function(v) { w <- v[v > -400 & v < 400]; d <- density(w); m <- d$x[which.max(d$y)]
    l <- v[v < m]; c(mode = m, sd = sqrt(mean((l - m)^2))) }))
  xs <- sweep(x, 2, bg[, "mode"]); thr <- 3 * median(bg[, "sd"])
  # probe -> MGI symbol: unique RefSeq hit, else the annotated symbol if current, else a unique alias
  a <- read.csv(gzfile(raw_file("GSE11211", "GPL6753_Mouse_Ref-8_V1.csv.gz")), colClasses = "character")
  a$pid <- as.character(as.integer(a$ProbeId)); a <- a[match(rownames(x), a$pid), ]
  acc <- sub("[.].*$", "", a$Accession)
  rs <- suppressMessages(mapIds(org.Mm.eg.db, keys = unique(acc[!is.na(acc) & nzchar(acc)]), column = "SYMBOL",
                                keytype = "REFSEQ", multiVals = "list"))
  rs1 <- vapply(rs, function(s) if (length(s) == 1 && !is.na(s)) s else NA_character_, "")
  sym_rs <- unname(rs1[acc])
  valid <- keys(org.Mm.eg.db, keytype = "SYMBOL")
  s_ann <- a$Symbol
  sym_cur <- ifelse(!is.na(s_ann) & s_ann %in% valid, s_ann, NA)
  al <- suppressMessages(mapIds(org.Mm.eg.db, keys = unique(s_ann[!is.na(s_ann) & nzchar(s_ann) & !(s_ann %in% valid)]),
                                column = "SYMBOL", keytype = "ALIAS", multiVals = "list"))
  al1 <- vapply(al, function(s) if (length(unique(s)) == 1 && !is.na(s[1])) s[1] else NA_character_, "")
  sym_al <- unname(al1[s_ann])
  sym <- ifelse(!is.na(sym_rs), sym_rs, ifelse(!is.na(sym_cur), sym_cur, sym_al))
  # log2 ratio SUS - RES over the 4 technical replicate arrays; probes expressed in either group; flip = swap labels
  log_ratio <- function(reg, flip) {
    j <- ann$region == reg; g <- ann$group[j]; if (flip) g <- ifelse(g == "SUS", "RES", "SUS")
    xr <- xs[, j]
    expressed <- pmax(rowMeans(xr[, g == "SUS"]), rowMeans(xr[, g == "RES"])) > thr
    lx <- log2(pmax(xr, thr))
    lr <- rowMeans(lx[, g == "SUS"]) - rowMeans(lx[, g == "RES"])
    keep <- expressed & !is.na(sym)
    collapse_by_symbol(cbind(SUS_vs_RES = lr[keep]), sym[keep], rowMeans(lx)[keep])
  }
  base <- list(study = "Schmidt et al. 2010", pmid_or_doi = "PMID:21159965", accession = "GSE11211",
    species = "Mus musculus (outbred CD1)", sex = "male",
    assay = "Illumina MouseRef-8 v1.0 BeadChip (GPL6753), BeadStudio-normalised values from the series matrix",
    model = "chronic social instability in adolescence (cage composition changed twice weekly for 7 weeks); SUS = top 20% basal corticosterone, RES = lowest 20%",
    timing = "5 weeks after the end of stress", stat_type = "log_ratio", n_CON = NA, n_RES = 1L, n_SUS = 1L, n_stressed = 2L,
    pooled_samples = "yes: one RNA pool per group (6 mice per paper, 5 per GEO), 4 technical replicate arrays per pool; no unstressed control")
  note <- "Single pool per group: the log ratio mixes phenotype with pool genotype (outbred CD1) and dissection composition; no biological replication. Values floored at 3 SD above the background-noise peak before log2; probes expressed in either group; highest-mean probe per gene."
  ca1_units <- c("CA1so", "CA1sr", "CA1slm", "CA1sp", "CA1microglia")
  list(
    do.call(sig, c(list(id = "Schmidt2010_dCA1", stat = log_ratio("CA1", flip = TRUE)), base,
      list(region = "dorsal hippocampus CA1 (laser microdissected)", match_units = ca1_units, set = "primary",
           corrections = "GSE11211 CA1 labels inverted to match the paper: GEO CA1 arrays titled 'Susceptible' (GSM282626-29) are used as RES and 'Resistant' (GSM282630-33) as SUS, because the deposited orientation contradicts the published CA1 results (Gria1, Gria2, Slc6a15 1.9-fold down in SUS).",
           notes = paste(note, "The oriented SUS CA1 pool carries more DG (C1ql2, Prox1) and choroid plexus (Ttr) signal, i.e. a dissection-composition difference.")))),
    do.call(sig, c(list(id = "Schmidt2010_dDG", stat = log_ratio("DG", flip = FALSE)), base,
      list(region = "dorsal hippocampus dentate gyrus (laser microdissected)", match_units = c("DGmo", "DGpo", "DGsg", "DGmicroglia"),
           set = "primary", corrections = "none (DG labels as deposited agree with the paper)", notes = note))),
    do.call(sig, c(list(id = "SENS_Schmidt2010_dCA1_asGEO", stat = log_ratio("CA1", flip = FALSE)), base,
      list(region = "dorsal hippocampus CA1 (laser microdissected)", match_units = ca1_units, set = "exploratory",
           corrections = "none: sensitivity copy of Schmidt2010_dCA1 with the CA1 labels as deposited on GEO (which contradict the paper)",
           notes = paste("SENSITIVITY ONLY.", note))))
  )
}

## --- GSE107941: Lehmann 2018, whole-brain microglia (arrays) ------------------------------------
build_lehmann2018 <- function() {
  sm <- raw_file("GSE107941", "GSE107941_series_matrix.txt.gz")
  hdr <- readLines(gzfile(sm))
  fld <- function(pat) { l <- grep(pat, hdr, value = TRUE)[1]; gsub('"', "", strsplit(l, "\t")[[1]][-1]) }
  ann <- data.frame(gsm = fld("^!Sample_geo_accession\t"), title = fld("^!Sample_title\t"),
                    stress = fld("^!Sample_characteristics_ch1\t\"stress"), stringsAsFactors = FALSE)
  map <- c("stress: Control home caged" = "CON", "stress: Susceptible" = "SUS", "stress: Resilient" = "RES")
  ann$group <- unname(map[ann$stress]); stopifnot(!anyNA(ann$group))
  tl <- sub("^ML[0-9]+_([A-Z]+)[0-9].*$", "\\1", ann$title)
  stopifnot(all(c(HC = "CON", S = "SUS", R = "RES")[tl] == ann$group))
  x <- read.delim(gzfile(sm), comment.char = "!", check.names = FALSE)
  rownames(x) <- as.character(x$ID_REF); x <- as.matrix(x[, ann$gsm])
  # platform table (gene_assignment, category) from the family SOFT file
  L <- readLines(gzfile(raw_file("GSE107941", "GSE107941_family.soft.gz")), n = 42000)
  b <- grep("^!platform_table_begin", L); e <- grep("^!platform_table_end", L)
  p <- read.delim(text = L[(b + 1):(e - 1)], quote = "", check.names = FALSE, colClasses = "character")
  p <- p[match(rownames(x), p$ID), ]
  first <- sub(" /// .*$", "", p$gene_assignment)       # first gene assignment, field 5 = Entrez ID
  ez <- vapply(strsplit(first, " // ", fixed = TRUE), function(z)
    if (length(z) >= 5 && grepl("^[0-9]+$", trimws(z[5]))) trimws(z[5]) else NA_character_, "")
  ok <- !is.na(ez) & !grepl("^control|^normgene|^reporter", p$category)
  sym <- rep(NA_character_, nrow(x)); sym[ok] <- entrez_symbol(ez[ok])
  # expression filter: max group mean above the 95th percentile of the antigenomic background probe sets
  thr <- quantile(x[p$category %in% "control->bgp->antigenomic", , drop = FALSE], 0.95)
  gm <- sapply(c("CON", "SUS", "RES"), function(g) rowMeans(x[, ann$group == g, drop = FALSE]))
  keep <- apply(gm, 1, max) > thr & !is.na(sym)
  xe <- x[keep, ]
  grp <- factor(ann$group, levels = c("CON", "SUS", "RES"))
  des <- model.matrix(~ 0 + grp); colnames(des) <- levels(grp)
  tt <- trend_t(xe, des, CONTRASTS)                        # RES arrays stay in the model (variance estimate)
  st <- collapse_by_symbol(tt, sym[keep], rowMeans(xe))[, "SUS_vs_CON", drop = FALSE]
  list(sig("Lehmann2018_microglia", st, study = "Lehmann et al. 2018", pmid_or_doi = "PMID:30050134",
    accession = "GSE107941", species = "Mus musculus (C57BL/6N)", sex = "male",
    region = "whole brain minus cerebellum; CD11b-high magnetically sorted microglia",
    assay = "Affymetrix MoGene 2.0 ST (GPL16570), RMA log2 from the series matrix",
    model = "chronic social defeat, 14 d continuous sensory contact with a CD-1 aggressor; SUS = SI ratio < 1.2, RES = SI >= 1.2 and urine-marking preference > 20%",
    timing = "2 h after the last (day 14) defeat", stat_type = "moderated_t", n_CON = 6L, n_RES = 9L, n_SUS = 6L, n_stressed = 15L,
    pooled_samples = "no: one mouse per array", match_units = c("CA1microglia", "CA2microglia", "CA3microglia", "DGmicroglia"),
    set = "primary",
    corrections = "GSE107941: all RES arrays were scanned on another day (2014-11-06; CON and SUS 2014-10-30, CEL headers), so every RES contrast is confounded with scan batch and only SUS_vs_CON is provided.",
    notes = "limma-trend, robust eBayes; RES arrays kept in the model for the variance estimate. Whole-brain microglia, not hippocampus; tissue 2 h after the last defeat, so acute defeat responses are included. Probe sets: first Entrez ID in gene_assignment; filter on antigenomic background; highest-mean probe set per gene."))
}

## --- PRD000546: Henningsen 2012, rat CA1+CA3 iTRAQ proteome ----------------------------------------
build_henningsen2012 <- function() {
  td <- tempfile("prd000546"); dir.create(td)
  xls <- file.path(td, sprintf("Suppl_Table_%d.xls", 1:3))
  utils::unzip(raw_file("PRD000546", "mmc1.zip"), files = basename(xls), exdir = td)
  rd <- function(i) {
    d <- as.data.frame(readxl::read_excel(xls[i], col_names = FALSE, .name_repair = "minimal"))
    h <- paste(d[9, ], d[10, ]); b <- d[-(1:10), ]
    list(hdr = h, tab = data.frame(frac = b[, 1], ipi = b[, 2], desc = b[, 3], r1 = as.numeric(b[, 9]), r2 = as.numeric(b[, 10]),
                                   r3 = as.numeric(b[, 11]), meann = as.numeric(b[, 13]), stringsAsFactors = FALSE))
  }
  T1 <- rd(1); T2 <- rd(2); T3 <- rd(3)
  # ratio directions from the replicate column headers: T1 = C_k/A_k, T2 = R_k/A_k, T3 = C_k/R_k (A = anhedonic = SUS)
  stopifnot(identical(T1$hdr[9:11], c("iTRAQ C1/A1", "iTRAQ C2/A2", "iTRAQ C3/A3")),
            identical(T2$hdr[9:11], c("iTRAQ R1/A1", "iTRAQ R2/A2", "iTRAQ R3/A3")),
            identical(T3$hdr[9:11], c("iTRAQ C1/R1", "iTRAQ C2/R2", "iTRAQ C3/R3")))
  key <- function(t) paste(t$tab$frac, t$tab$ipi)
  k1 <- key(T1); stopifnot(!anyDuplicated(k1), setequal(k1, key(T2)), setequal(k1, key(T3)))
  a1 <- T1$tab; a2 <- T2$tab[match(k1, key(T2)), ]; a3 <- T3$tab[match(k1, key(T3)), ]
  stopifnot(identical(a1$desc, a2$desc), identical(a1$desc, a3$desc))
  lCA <- log2(as.matrix(a1[, c("r1", "r2", "r3")])); lRA <- log2(as.matrix(a2[, c("r1", "r2", "r3")]))
  lCR <- log2(as.matrix(a3[, c("r1", "r2", "r3")]))
  fr <- a1$frac
  ctr <- function(m) {                                     # centre each column on its per-run (study x fraction) median
    med <- sapply(c("M", "S"), function(f) apply(m[fr == f, , drop = FALSE], 2, median))
    for (f in c("M", "S")) m[fr == f, ] <- sweep(m[fr == f, , drop = FALSE], 2, med[, f])
    m }
  lCA <- ctr(lCA); lRA <- ctr(lRA); lCR <- ctr(lCR)
  M <- list(SUS_vs_CON = -lCA, RES_vs_CON = -lCR, SUS_vs_RES = -lRA, stressed_vs_CON = -(lCA + lCR) / 2)
  # rat symbol from the IPI description (first token, org.Rn.eg.db symbol or unique alias), then 1:1 mouse ortholog
  tok <- sub(" .*$", "", a1$desc)
  rs <- keys(org.Rn.eg.db::org.Rn.eg.db, keytype = "SYMBOL")
  fix1 <- function(s) {
    parts <- unlist(strsplit(s, ";", fixed = TRUE)); parts <- parts[!grepl("^LOC|-ps[0-9]*$", parts)]
    if (!length(parts)) return(NA_character_)
    cand <- unique(vapply(parts, function(p) {
      if (p %in% rs) return(p)
      p2 <- paste0(toupper(substr(p, 1, 1)), tolower(substring(p, 2))); if (p2 %in% rs) return(p2)
      a <- tryCatch(suppressMessages(mapIds(org.Rn.eg.db::org.Rn.eg.db, keys = p, column = "SYMBOL", keytype = "ALIAS",
                                            multiVals = "list"))[[1]], error = function(e) NA)
      a <- unique(a[!is.na(a)]); if (length(a) == 1) a else NA_character_ }, ""))
    cand <- cand[!is.na(cand)]
    if (length(cand) == 1) cand else NA_character_
  }
  ok_tok <- grepl("^[A-Za-z][A-Za-z0-9.;-]*$", tok) &
    !(tok %in% c("Similar", "similar", "Protein", "Isoform", "Putative", "Hypothetical", "Uncharacterized"))
  rat <- rep(NA_character_, length(tok)); u <- unique(tok[ok_tok]); um <- vapply(u, fix1, ""); rat[ok_tok] <- um[tok[ok_tok]]
  mouse <- rat_to_mouse(rat)
  # one-sample limma per contrast on the 3 per-study log2 ratios
  tt <- sapply(M, function(m) { rownames(m) <- k1; eBayes(lmFit(m, matrix(1, 3, 1, dimnames = list(NULL, "mean"))), robust = TRUE)$t[, 1] })
  # one row per mouse gene: the protein x fraction row with most quantified spectra
  scans <- a1$meann + a2$meann + a3$meann
  keep <- !is.na(mouse)
  o <- order(mouse[keep], -scans[keep]); idx <- which(keep)[o]; idx <- idx[!duplicated(mouse[idx])]
  st <- tt[idx, , drop = FALSE]; rownames(st) <- mouse[idx]
  unlink(td, recursive = TRUE)
  list(sig("Henningsen2012_CA", st, study = "Henningsen, Palmfeldt et al. 2012",
    pmid_or_doi = "PMID:22311638", accession = "PRD000546 (ratios from the MCP supplement mmc1.zip, Suppl_Table_1-3)",
    species = "Rattus norvegicus (Wistar); 1:1 mouse orthologs via human (babelgene)", sex = "male",
    region = "hippocampus CA1 + CA3 (gross subdissection, pooled), membrane and soluble fractions",
    assay = "iTRAQ LC-MS/MS proteomics, 3 independent studies x 2 fractions = 6 iTRAQ runs",
    model = "chronic mild stress 8 weeks; SUS = anhedonic-like (> 30% sucrose-intake drop), RES = resilient (< 10% drop), CON = unchallenged",
    timing = "end of 8-week chronic mild stress", stat_type = "moderated_t", n_CON = 3L, n_RES = 3L, n_SUS = 3L, n_stressed = 6L,
    pooled_samples = "yes: each replicate is a pool of 3 rats (9 rats per group); replicates paired by iTRAQ run",
    match_units = c("CA1so", "CA1sr", "CA1slm", "CA1sp", "CA1microglia", "CA2so", "CA2sr", "CA2slm", "CA2sp", "CA2microglia",
                    "CA3so", "CA3sr", "CA3sp", "CA3microglia"),
    set = "primary", corrections = "none",
    notes = "One-sample limma (robust eBayes, no trend) on the 3 per-study log2 ratios, each centred on its per-run median; residual df = 2, so the moderated t leans on the prior. iTRAQ ratio compression: effect sizes are underestimated, ranks are the usable information. Protein rows without a gene symbol in the IPI description are lost; membrane/soluble rows of one gene: most spectra kept. CA2 is included in match_units because a gross CA1 + CA3 dissection contains it."))
}

## --- GSE72343: Bagot 2016, ventral hippocampus, 48 h and 28 d arms ----------------------------------
build_bagot2016 <- function() {
  sm <- read_series_matrix(raw_file("GSE72343", "GSE72343_series_matrix.txt.gz"))$ann
  ann <- sm[sm$`brain region` == "adult ventral hippocampus" & !grepl("DKKL1|GFP", sm$title), ]
  ann$srx <- sub(".*term=", "", ann$sra)
  ann$tp <- c("48h post CSDS" = "48h", "28d post CSDS" = "28d", "28d +1h stress post CSDS" = "28d_1h")[ann$timepoint]
  ann$grp_char <- c(control = "CON", resilient = "RES", susceptible = "SUS")[ann$phenotype]
  ann$grp_title <- ifelse(grepl("_Control_", ann$title), "CON", ifelse(grepl("_Resilient_", ann$title), "RES",
                          ifelse(grepl("_Susceptible_", ann$title), "SUS", NA)))
  stopifnot(!anyNA(ann$tp), !anyNA(ann$grp_char), !anyNA(ann$grp_title))
  x <- read.csv(raw_file("GSE72343", "PRJNA293822_counts.csv"), check.names = FALSE, row.names = 1)
  cnt_all <- as.matrix(x[, ann$srx]); storage.mode(cnt_all) <- "double"
  stopifnot(all(cnt_all == round(cnt_all)), !anyNA(cnt_all))
  mgi <- mgi_symbol(rownames(cnt_all), drop_clash = FALSE)
  fit_arm <- function(tp, lab) {
    a <- ann[ann$tp == tp, ]
    g <- factor(a[[lab]], levels = c("CON", "RES", "SUS"))
    des <- model.matrix(~ 0 + g); colnames(des) <- levels(g)
    r <- voom_mor(cnt_all[, a$srx], des, CONTRASTS, min_samples = 3)
    collapse_by_symbol(r$t, mgi[match(rownames(r$t), rownames(cnt_all))], r$ave)
  }
  a28 <- ann[ann$tp == "28d", ]
  stopifnot(identical(sort(a28$gsm[a28$grp_char != a28$grp_title]), c("GSM1860619", "GSM1860620")))
  base <- list(study = "Bagot et al. 2016", pmid_or_doi = "PMID:27181059",
    accession = "GSE72343 (counts: Flati 2020 Stress Mice Portal re-quantification of PRJNA293822, figshare file 21808590)",
    species = "Mus musculus (C57BL/6J)", sex = "male", region = "ventral hippocampus (bulk punch)",
    assay = "bulk RNA-seq (total RNA), StringTie counts",
    model = "10-day chronic social defeat stress; SUS/RES by social interaction ratio", stat_type = "moderated_t",
    pooled_samples = "yes: each sample is a pool of 3-5 mice", match_units = "all", set = "primary")
  arm_note <- "voom with median-of-ratios size factors, ~0 + group, robust eBayes; each time point fitted separately as in the paper; no batch covariate available."
  list(
    do.call(sig, c(list(id = "Bagot2016_vHIP_28d", stat = fit_arm("28d", "grp_title")), base,
      list(timing = "28 d after the last defeat (home cage, no re-exposure)", n_CON = 4L, n_RES = 4L, n_SUS = 4L, n_stressed = 8L,
           corrections = "GSE72343 28d uses the GEO sample TITLE labels: GSM1860619/GSM1860620 are titled VHIP_Susceptible_28d_1/_2 and assigned to the susceptible processed file, but their 'phenotype' characteristic says resilient; titles match the paper (n = 4 per group) and the authors' own group means.",
           notes = paste(arm_note, "Only the 28 d arm without aggressor re-exposure is used (the stress-primed arm is excluded).")))),
    do.call(sig, c(list(id = "Bagot2016_vHIP_48h", stat = fit_arm("48h", "grp_char")), base,
      list(timing = "48 h after the last defeat (home cage)", n_CON = 3L, n_RES = 3L, n_SUS = 3L, n_stressed = 6L,
           corrections = "none", notes = arm_note)))
  )
}

## --- GSE81672: Bagot 2017, ventral hippocampus, saline arm ------------------------------------------
build_bagot2017 <- function() {
  cells <- read_xls_biff8(raw_file("GSE81672", "GSE81672_HTSeq_counts_matrix.xls.gz"))
  num <- cells$num; lab <- cells$lab
  hdr <- num[num$r == 0, ]; hdr <- hdr[order(hdr$c), ]               # header row = sample numbers
  ids <- lab[lab$c == 0 & lab$r > 0, ]; ids <- ids[order(ids$r), ]   # first column = Ensembl gene ID
  stopifnot(!any(lab$r == 0 & lab$c > 0), !anyDuplicated(ids$v), identical(ids$r, seq_len(nrow(ids))))
  body <- num[num$r > 0 & num$c > 0, ]
  x <- matrix(NA_real_, nrow(ids), nrow(hdr), dimnames = list(ids$v, format(hdr$v, scientific = FALSE, trim = TRUE)))
  x[cbind(body$r, match(body$c, hdr$c))] <- body$v
  stopifnot(!anyNA(x))
  # saline-treated vHIP samples, groups from the GEO sample titles
  smp <- c(`12` = "CON", `61` = "CON", `80` = "CON", `103` = "CON", `36` = "RES", `58` = "RES", `81` = "RES", `104` = "RES",
           `13` = "SUS", `105` = "SUS", `108` = "SUS")
  cnt <- x[, names(smp)]
  g <- factor(smp[colnames(cnt)]); des <- model.matrix(~ 0 + g); colnames(des) <- levels(g)
  r <- voom_mor(cnt, des, CONTRASTS, min_samples = min(table(g)))
  st <- collapse_by_symbol(r$t, ens_symbol(rownames(r$t)), rowMeans(log1p(cnt[rownames(r$t), ])))
  list(sig("Bagot2017_vHIP", st, study = "Bagot et al. 2017", pmid_or_doi = "PMID:27569543",
    accession = "GSE81672", species = "Mus musculus", sex = "male", region = "ventral hippocampus (bulk)",
    assay = "bulk RNA-seq, HTSeq counts (GSE81672_HTSeq_counts_matrix.xls)",
    model = "10-day chronic social defeat stress; saline-treated arm of the ketamine/imipramine study",
    timing = "after the saline treatment arm (days to weeks after defeat; see paper)", stat_type = "moderated_t",
    n_CON = 4L, n_RES = 4L, n_SUS = 3L, n_stressed = 7L, pooled_samples = "not recorded in the checked loader (see paper)",
    match_units = "all", set = "primary", corrections = "none",
    notes = "Saline-treated vHIP only: samples 12/61/80/103 CON, 36/58/81/104 RES, 13/105/108 SUS (GEO sample titles). voom with median-of-ratios size factors, robust eBayes; duplicate symbols -> highest mean log1p count. The .xls is read by a built-in BIFF8 reader (readxl cannot open it); identical to Excel's CSV export."))
}

## --- GSE109315: Laine 2018, ventral hippocampus, C57BL/6 only -----------------------------------
build_laine2018 <- function() {
  ann <- read_series_matrix(raw_file("GSE109315", "GSE109315-GPL13112_series_matrix.txt.gz"))$ann
  b6 <- ann[ann$source == "ventral hippocampus" & ann$strain == "C57BL/6Crl", ]
  b6$sample <- sub("_.*$", "", b6$title)
  b6$grp <- c(control = "CON", resilient = "RES", susceptible = "SUS")[b6$`csds phenotype`]
  stopifnot(!anyNA(b6$grp), all(b6$Sex == "male"))
  b6$cohort <- sub("[.].*$", "", b6$sample)
  x <- read.delim(gzfile(raw_file("GSE109315", "GSE109315_vHPC_counts_logCPM.txt.gz")), check.names = FALSE, row.names = 1)
  lx <- as.matrix(x[, b6$sample])
  lx <- lx[rowSums(lx > 0) >= 3, ]                          # CPM > 1 in >= 3 B6 samples
  g <- factor(b6$grp, levels = c("CON", "RES", "SUS"))
  b6$C2 <- as.numeric(b6$cohort == "C2")                    # defeat cohort (C1/C2 sample prefix)
  des <- model.matrix(~ 0 + g + C2, data = b6); colnames(des)[1:3] <- levels(g)
  st <- collapse_by_symbol(trend_t(lx, des, CONTRASTS), ens_symbol(rownames(lx)), rowMeans(lx))
  list(sig("Laine2018_vHPC", st, study = "Laine et al. 2018", pmid_or_doi = "PMID:30073192",
    accession = "GSE109315 (GSE109315_vHPC_counts_logCPM.txt.gz)", species = "Mus musculus (C57BL/6NCrl)", sex = "male",
    region = "ventral hippocampus (bulk dissection)", assay = "bulk RNA-seq (ribo-depleted total RNA), deposited logCPM",
    model = "10-day chronic social defeat stress; SUS/RES by social interaction ratio", timing = "6-8 d after the last defeat",
    stat_type = "moderated_t", n_CON = 6L, n_RES = 8L, n_SUS = 3L, n_stressed = 11L, pooled_samples = "no: individual mice",
    match_units = "all", set = "primary", corrections = "none",
    notes = "limma-trend (robust) on the deposited logCPM, ~0 + group + defeat cohort; C57BL/6 samples only (DBA/2 excluded). Sequencing set is partly confounded with group (all CON in set 1, 7/8 RES in set 2) and not modelled. Very strong CON-vs-defeated separation incl. a chr4 copy-number-variable cluster (Fam205a1/Fam205c/Gm3198) suggests part of it may be sub-colony background."))
}

## --- GSE161726: Ayash 2023, ventral hippocampus punches ----------------------------------------------
build_ayash2023 <- function() {
  ann <- read_series_matrix(raw_file("GSE161726", "GSE161726_series_matrix.txt.gz"))$ann
  ann <- ann[ann$tissue == "Ventral hippocampus (vHC)", ]
  stopifnot(nrow(ann) == 25)
  ann$sub <- c(Control = "Control", `Non-avoiders` = "NonAvoider", `Indiscriminate-avoiders` = "IndiscAvoider",
               `Discriminating-avoiders` = "DiscAvoider")[ann$`behavioural phenotype`]
  # a-priori mapping: Control = CON; Non-avoiders = RES; Indiscriminate + Discriminating avoiders = SUS
  ann$grp <- c(Control = "CON", NonAvoider = "RES", IndiscAvoider = "SUS", DiscAvoider = "SUS")[ann$sub]
  stopifnot(!anyNA(ann$grp))
  rd <- function(f) {                                       # featureCounts file; BAM path (header) gives the pool
    con <- gzfile(f); hdr <- readLines(con, n = 2); close(con)
    d <- read.delim(gzfile(f), comment.char = "#", check.names = FALSE)
    list(id = d$Geneid, cnt = d[[7]], path = strsplit(hdr[2], "\t")[[1]][7])
  }
  L <- lapply(basename(ann$suppl), function(f) rd(raw_file("GSE161726", f)))
  stopifnot(all(sapply(L, function(z) identical(z$id, L[[1]]$id))))
  cnt <- sapply(L, `[[`, "cnt"); rownames(cnt) <- L[[1]]$id; colnames(cnt) <- ann$gsm
  bam_dir <- sapply(L, function(z) sub("/transcript_quantification.*$", "", sub("^.*/Mueller_UM/", "", z$path)))
  ann$pool <- sub("_.*$", "", bam_dir)
  stopifnot(all(sub("^Pool[0-9]+_", "", bam_dir) == ann$description))
  g <- factor(ann$grp, levels = c("CON", "RES", "SUS"))
  des <- model.matrix(~ 0 + g + pool, data = ann); colnames(des)[1:3] <- levels(g)
  r <- voom_mor(cnt, des, CONTRASTS, min_samples = 5)
  st <- collapse_by_symbol(r$t, ens_symbol(rownames(r$t)), r$ave)
  list(sig("Ayash2023_vHC", st, study = "Ayash et al. 2023", pmid_or_doi = "PMID:37068238",
    accession = "GSE161726 (25 vHC featureCounts files)", species = "Mus musculus (C57BL/6)", sex = "male",
    region = "ventral hippocampus (0.75-1 mm punch)", assay = "low-input bulk RNA-seq (SMART-Seq v4), featureCounts GENCODE M17",
    model = "10-day chronic social defeat; subgroups from a social threat-safety test",
    timing = "~1 d after the last defeat (after the first social interaction test)", stat_type = "moderated_t",
    n_CON = 6L, n_RES = 7L, n_SUS = 12L, n_stressed = 19L, pooled_samples = "no: individual mice", match_units = "all",
    set = "primary", corrections = "none",
    notes = "SUS = Indiscriminate-avoiders (7) + Discriminating-avoiders (5), RES = Non-avoiders (7), fixed a priori; the authors call the Discriminating-avoiders resilient and place Non-avoiders outside the resilience-susceptibility axis. voom with median-of-ratios size factors, ~0 + group + sequencing pool, robust eBayes; all 6 CON were sequenced in Pool1."))
}

## --- GSE240573: Butto 2024, ventral hippocampus Arc-tagged nuclei ----------------------------------
build_butto2024 <- function() {
  rd <- function(f) {
    x <- read.delim(gzfile(raw_file("GSE240573", f)), check.names = FALSE, stringsAsFactors = FALSE)
    x[, nzchar(names(x)) & !grepl("^V[0-9]+$|^\\.\\.\\.", names(x)), drop = FALSE]
  }
  L <- lapply(c(CON = "GSE240573_vHip_Control_Counts.txt.gz", RES = "GSE240573_vHip_Resilience_Counts.txt.gz",
                SUS = "GSE240573_vHip_Susceptible_Counts.txt.gz"), rd)
  stopifnot(identical(L$CON$Gene, L$RES$Gene), identical(L$CON$Gene, L$SUS$Gene))
  cnt <- as.matrix(cbind(L$CON[, -1], L$RES[, -1], L$SUS[, -1])); rownames(cnt) <- L$CON$Gene
  grp <- rep(names(L), sapply(L, function(x) ncol(x) - 1))
  stopifnot(all(c(Control = "CON", Resilience = "RES", Susceptible = "SUS")[sub("_[Rr]ep[0-9]+$", "", colnames(cnt))] == grp))
  g <- factor(grp, levels = c("CON", "RES", "SUS")); des <- model.matrix(~ 0 + g); colnames(des) <- levels(g)
  r <- voom_mor(cnt, des, CONTRASTS, min_samples = min(table(grp)))
  st <- collapse_by_symbol(r$t, mgi_symbol(rownames(r$t), keep_unmapped = TRUE), r$ave)
  list(sig("Butto2024_vHip_nuclei", st, study = "Butto et al. 2024", pmid_or_doi = "PMID:38649752",
    accession = "GSE240573", species = "Mus musculus (Arc-CreERT2 x Sun1-GFP)", sex = "male", region = "ventral hippocampus",
    assay = "FANS-sorted Arc-GFP+ (tamoxifen-tagged activated) nuclei, nuclear RNA-seq, HTSeq counts",
    model = "10-day chronic social defeat; SI index split (RES > 100, SUS < 100); handled controls",
    timing = "tagged around the SI test ~8 d after the last defeat; tissue ~11 d after the end of stress",
    stat_type = "moderated_t", n_CON = 6L, n_RES = 7L, n_SUS = 6L, n_stressed = 13L,
    pooled_samples = "yes: nuclei pooled from 2 mice per sample", match_units = "all", set = "primary", corrections = "none",
    notes = "Activated (Arc-tagged), mostly glutamatergic neuronal nuclei; GFP+ identity of the deposited files inferred from the paper's supplementary tables. voom with median-of-ratios size factors, robust eBayes, no covariates; Resilience_Rep7 (lowest depth) is a PCA outlier, kept as by the authors."))
}

## --- GSE184027: Bhatti 2022, hippocampal PV-interneuron TRAP -----------------------------------------
build_bhatti2022 <- function() {
  x <- read.csv(gzfile(raw_file("GSE184027", "GSE184027_rawcount.txt.gz")), check.names = FALSE, stringsAsFactors = FALSE)
  cnt_all <- as.matrix(x[, -1]); rownames(cnt_all) <- x$Gene
  # GEO: CTRL = CON (C1, C6, C9, C11); Resilient = RES (R3, R5, R7, R8); 'Susceptiable (STR)' = SUS (S2, S4, S10); inputs dropped
  map <- c(C1 = "CON", C6 = "CON", C9 = "CON", C11 = "CON", R3 = "RES", R5 = "RES", R7 = "RES", R8 = "RES",
           S2 = "SUS", S4 = "SUS", S10 = "SUS", CI12 = "INPUT", CI13 = "INPUT", CI14 = "INPUT")
  stopifnot(setequal(names(map), colnames(cnt_all)))
  o <- as.data.frame(suppressMessages(readxl::read_excel(raw_file("GSE184027", "Table_2.xlsx"), sheet = "Original")))
  sym_all <- ens_symbol_unique(rownames(cnt_all), annot = setNames(o$GeneSymbol.x, o$Row.names))
  keep_s <- names(map)[map != "INPUT"]
  grp <- map[keep_s]
  g <- factor(grp, levels = c("CON", "RES", "SUS")); des <- model.matrix(~ 0 + g); colnames(des) <- levels(g)
  r <- voom_mor(cnt_all[, keep_s], des, CONTRASTS, min_samples = min(table(grp)))
  st <- collapse_by_symbol(r$t, sym_all[match(rownames(r$t), rownames(cnt_all))], r$ave)
  list(sig("Bhatti2022_PV_TRAP", st, study = "Bhatti et al. 2022", pmid_or_doi = "PMID:35813065",
    accession = "GSE184027", species = "Mus musculus (PV-Cre x Rosa26 EGFP-L10a)", sex = "male",
    region = "whole hippocampus (both hippocampi of one mouse per sample)",
    assay = "PV-interneuron TRAP (anti-EGFP polysome IP), SMART-Seq v4, HTSeq counts",
    model = "10-day chronic social defeat; SI ratio split (RES >= 1, SUS < 1); non-defeated controls",
    timing = "dissection >= ~2 d after the last defeat", stat_type = "moderated_t", n_CON = 4L, n_RES = 4L, n_SUS = 3L,
    n_stressed = 7L, pooled_samples = "no: one mouse per sample", match_units = "all", set = "primary", corrections = "none",
    notes = "PV+ interneuron translatome. Input samples CI12-CI14 excluded. voom with median-of-ratios size factors, robust eBayes, no covariates; a non-group PC1 (35%) inflates residual variance and no gene passes BH 5%. Ensembl IDs: deterministic choice for one-to-many IDs; IDs reassigned since the counting annotation keep the authors' symbol (Supp Table 2 'GeneSymbol.x')."))
}

## --- GSE139684: Milosevic & Russo (unpublished), hippocampal astrocyte TRAP --------------------------
build_russo_astro <- function() {
  ff <- sort(c("GSM4146188_1690_H_C_1_TAAGGCGAAGAGTAGA_L001_paired.sortedByCoord.bam.counts.txt.gz",
               "GSM4146189_1691_HIP_C1_CGTACTAGAGAGTAGA_L001_paired.sortedByCoord.bam.counts.txt.gz",
               "GSM4146190_1692_HIP_C3_AGGCAGAAAGAGTAGA_L001_paired.sortedByCoord.bam.counts.txt.gz",
               "GSM4146191_1693_H_RS1_TCCTGAGCGTAAGGAG_L001_paired.sortedByCoord.bam.counts.txt.gz",
               "GSM4146192_1694_H_RS_2_GGACTCCTGTAAGGAG_L001_paired.sortedByCoord.bam.counts.txt.gz",
               "GSM4146193_1695_HIP_R2_TAGGCATGGTAAGGAG_L001_paired.sortedByCoord.bam.counts.txt.gz",
               "GSM4146194_1696_HIP_SS1_CTCTCTACGTAAGGAG_L001_paired.sortedByCoord.bam.counts.txt.gz",
               "GSM4146195_1697_HIP_SS2_CAGAGAGGGTAAGGAG_L001_paired.sortedByCoord.bam.counts.txt.gz",
               "GSM4146196_1698_HIP_SS3_GCTACGCTACTGCATA_L001_paired.sortedByCoord.bam.counts.txt.gz",
               "GSM4146197_1699_H_SS_2_CGAGGCTGACTGCATA_L001_paired.sortedByCoord.bam.counts.txt.gz",
               "GSM4146198_1700_H_SS_3_AAGAGGCAACTGCATA_L001_paired.sortedByCoord.bam.counts.txt.gz"))
  # group from the GEO sample titles ('... hippocampus, control/resilient/stressed'), checked against the file tag
  gsm_txt <- readLines(raw_file("GSE139684", "GSE139684_gsm.txt"), warn = FALSE)
  acc <- sub("^\\^SAMPLE = ", "", gsm_txt[grepl("^\\^SAMPLE", gsm_txt)])
  tit <- sub("^!Sample_title = ", "", gsm_txt[grepl("^!Sample_title", gsm_txt)])
  stopifnot(length(acc) == length(tit))
  title <- tit[match(sub("_.*$", "", ff), acc)]
  grp <- ifelse(grepl("hippocampus, control", title), "CON", ifelse(grepl("hippocampus, resilient", title), "RES",
                ifelse(grepl("hippocampus, stressed", title), "SUS", NA)))
  tag <- sub("^GSM[0-9]+_[0-9]+_(.*)_[ACGT]{16}_L001.*$", "\\1", ff)
  grp_tag <- ifelse(grepl("_C_?[0-9]$|^H_C|^HIP_C", tag), "CON", ifelse(grepl("RS|_R[0-9]$", tag), "RES", ifelse(grepl("SS", tag), "SUS", NA)))
  batch <- ifelse(grepl("^HIP_", tag), "HIP", ifelse(grepl("^H_", tag), "H", NA))   # file-naming batch
  stopifnot(!anyNA(grp), all(grp == grp_tag), !anyNA(batch))
  L <- lapply(ff, function(f) { x <- read.delim(gzfile(raw_file("GSE139684", f)), header = FALSE, stringsAsFactors = FALSE); setNames(x$V2, x$V1) })
  stopifnot(all(sapply(L, function(v) identical(names(v), names(L[[1]])))))
  cnt <- do.call(cbind, L); colnames(cnt) <- paste(grp, tag, sep = "|")
  cnt <- cnt[!grepl("^__", rownames(cnt)), ]                 # drop HTSeq summary rows
  # a-priori contamination scores (mean centred log2-CPM): neuronal background and choroid-plexus carry-over
  lcpm <- log2(t(t(cnt) / colSums(cnt)) * 1e6 + 1)
  neu_mk <- c("Snap25", "Syt1", "Slc17a7", "Gad1", "Gad2", "Rbfox3", "Stmn2", "Syp", "Nrgn", "Camk2a")
  stopifnot(all(neu_mk %in% rownames(lcpm)))
  bg <- colMeans(lcpm[neu_mk, ] - rowMeans(lcpm[neu_mk, ]))
  cp_mk <- intersect(c("Ttr", "Kl", "Folr1", "Clic6", "Kcnj13", "Sostdc1"), rownames(lcpm)); stopifnot(length(cp_mk) >= 4)
  cp <- colMeans(lcpm[cp_mk, ] - rowMeans(lcpm[cp_mk, ]))
  g <- factor(grp, levels = c("CON", "RES", "SUS")); b <- factor(batch, levels = c("HIP", "H"))
  des <- model.matrix(~ 0 + g + b); colnames(des) <- c(levels(g), "batchH")
  des <- cbind(des, neuro_bg = bg, choroid = cp)
  r <- voom_mor(cnt, des, CONTRASTS, min_samples = min(table(grp)))
  st <- collapse_by_symbol(r$t, mgi_symbol(rownames(r$t), keep_unmapped = TRUE), r$ave)
  list(sig("Russo_astro_TRAP", st, study = "Milosevic & Russo, unpublished GEO deposit", pmid_or_doi = "none (unpublished; GEO GSE139684)",
    accession = "GSE139684", species = "Mus musculus (Tg(Aldh1l1-EGFP/Rpl10a)JD130Htz)", sex = "male (inferred: Xist absent, Ddx3y expressed)",
    region = "hippocampus (whole)", assay = "astrocyte TRAP (Aldh1l1-driven EGFP-L10a polysome IP), RNA-seq, STAR + HTSeq counts",
    model = "10-day chronic social defeat; SI split (RES SI > 1, SUS SI < 1); group-housed controls",
    timing = "2 d after the last defeat", stat_type = "moderated_t", n_CON = 3L, n_RES = 3L, n_SUS = 5L, n_stressed = 8L,
    pooled_samples = "yes: hippocampi of 2 mice per sample", match_units = "all", set = "exploratory", corrections = "none",
    notes = "LOW CONFIDENCE, exploratory only. Design ~0 + group + file-naming batch (H_/HIP_) + neuronal-background and choroid-plexus scores, voom with median-of-ratios size factors, robust eBayes, 5 residual df. IP purity and choroid-plexus carry-over vary strongly between pools; the adjustment over-corrects RES_vs_CON and cannot remove thalamic carry-over in the SUS pools. Unpublished: no external sign check."))
}

## --- GSE278620: Lyons 2025, GeoMx CA1 and CA3 AOIs, chronic subordination stress ------------------
build_lyons2025 <- function() {
  x <- read.csv(gzfile(raw_file("GSE278620", "GSE278620_Lyons_et_al_geomix_median_norm_q_datasets2.0.csv.gz")), check.names = FALSE)
  m <- as.matrix(x[, -(1:2)]); rownames(m) <- x$NAME
  neg <- m["NegProbe-WTX", ]; m <- m[rownames(m) != "NegProbe-WTX", ]
  ann <- read_series_matrix(raw_file("GSE278620", "GSE278620_series_matrix.txt.gz"))$ann
  norm_id <- function(s) gsub("[ _]", "", s)
  ann$col <- colnames(m)[match(norm_id(ann$title), norm_id(colnames(m)))]
  stopifnot(!anyNA(ann$col), !anyDuplicated(ann$col), setequal(ann$col, colnames(m)),
            all((substring(norm_id(ann$title), nchar(norm_id(ann$title))) == "+") == (ann$p16 == "high")),
            all(tapply(ann$group, ann$`mouse id`, function(g) length(unique(g))) == 1))
  one <- function(sf, units) {
    a <- ann[ann$subregion == sf, ]
    a$group <- factor(a$group, levels = c("ctrl", "CSS")); a$p16 <- factor(a$p16, levels = c("low", "high"))
    keep <- rowMeans(sweep(m[, a$col], 2, neg[a$col], ">")) >= 0.2   # above the AOI NegProbe in >= 20% of AOIs
    y <- log2(m[keep, a$col])
    design <- model.matrix(~ group + p16, data = a)
    dc <- duplicateCorrelation(y, design, block = a$`mouse id`)
    eb <- eBayes(lmFit(y, design, block = a$`mouse id`, correlation = dc$consensus.correlation), trend = TRUE, robust = TRUE)
    tt <- eb$t[, "groupCSS", drop = FALSE]; colnames(tt) <- "stressed_vs_CON"
    st <- collapse_by_symbol(tt, mgi_symbol(rownames(tt)), eb$Amean)
    nm <- tapply(a$`mouse id`, a$group, function(z) length(unique(z))); na <- table(a$group)
    sig(paste0("Lyons2025_", sf), st, study = "Lyons et al. 2025", pmid_or_doi = "PMID:39528642",
      accession = "GSE278620", species = "Mus musculus (p16-3MR transgenic, presumed)", sex = "male",
      region = paste0("hippocampus ", sf, " (GeoMx AOIs of ~230 nuclei on p16-high and p16-low areas; dorsoventral level not stated)"),
      assay = "Nanostring GeoMx DSP mouse Whole Transcriptome Atlas, FFPE; deposited median-normalised, LOQ-filtered counts",
      model = "chronic subordination stress (CSS), 4 weeks (CD-1 resident pairing); singly housed controls",
      timing = "end of the 4-week stress", stat_type = "moderated_t", n_CON = unname(nm[["ctrl"]]), n_RES = NA, n_SUS = NA,
      n_stressed = unname(nm[["CSS"]]),
      pooled_samples = sprintf("no: %d CON and %d CSS AOIs from the mice (mouse = random block)", na[["ctrl"]], na[["CSS"]]),
      match_units = units, set = "primary", corrections = "none",
      notes = "limma ~ group + p16 status on log2 counts, mouse as random block (duplicateCorrelation), eBayes(trend, robust); genes above the AOI NegProbe level in >= 20% of the subfield's AOIs. AOIs were placed on p16-high and p16-low areas, so the CSS effect is averaged over p16 status; layer composition is not annotated. Many genes are near the NegProbe background.")
  }
  list(one("CA1", c("CA1so", "CA1sr", "CA1slm", "CA1sp", "CA1microglia")),
       one("CA3", c("CA3so", "CA3sr", "CA3sp", "CA3microglia")))
}

## --- GSE124169: Weber 2019, sorted brain microglia 24 d after repeated social defeat ---------------
build_weber2019 <- function() {
  x0 <- read.csv(gzfile(raw_file("GSE124169", "GSE124169_countmatrix.csv.gz")), check.names = FALSE)
  # Excel date-mangled symbols: '15-Sep' = Sep15, 'N-Sep' = SeptN, 'N-Mar' = MarchN; the duplicated '2-Mar' is dropped
  g <- x0[[1]]
  dup_date <- g %in% g[duplicated(g)]
  g <- ifelse(g == "15-Sep", "Sep15", g)
  g <- sub("^([0-9]+)-Sep$", "Sept\\1", g)
  g <- sub("^([0-9]+)-Mar$", "March\\1", g)
  x <- as.matrix(x0[!dup_date, -1]); rownames(x) <- g[!dup_date]
  stopifnot(!anyDuplicated(rownames(x)))
  a <- read_series_matrix(raw_file("GSE124169", "GSE124169_series_matrix.txt.gz"))$ann
  a$num <- sub("^(CON|RSD) ", "", a$title)
  stopifnot(all(sub(" .*", "", a$title) == a$stress), setequal(a$num, colnames(x)), !anyDuplicated(a$num))
  a <- a[match(colnames(x), a$num), ]
  grp <- factor(ifelse(a$stress == "RSD", "RSD", "CON"), levels = c("CON", "RSD"))
  stopifnot(sum(grp == "CON") == 6, sum(grp == "RSD") == 7)
  y <- log2(x + 1)[rowSums(x >= 10) >= 6, ]                 # normalised count >= 10 in >= 6 samples
  design <- model.matrix(~ grp)
  eb <- eBayes(lmFit(y, design, weights = arrayWeights(y, design)), trend = TRUE, robust = TRUE)
  tt <- eb$t[, "grpRSD", drop = FALSE]; colnames(tt) <- "stressed_vs_CON"
  st <- collapse_by_symbol(tt, mgi_symbol(rownames(tt)), eb$Amean)
  list(sig("Weber2019_microglia", st, study = "Weber et al. 2019", pmid_or_doi = "PMID:30527629",
    accession = "GSE124169", species = "Mus musculus (C57BL/6)", sex = "male",
    region = "whole brain (Percoll-enriched, FACS-sorted CD11b+/CD45int microglia)",
    assay = "bulk RNA-seq of sorted microglia; deposited DESeq2-normalised count matrix",
    model = "repeated social defeat (6 cycles, CD-1 aggressor); stress sensitisation", timing = "24 d after the last defeat cycle",
    stat_type = "moderated_t", n_CON = 6L, n_RES = NA, n_SUS = NA, n_stressed = 7L,
    pooled_samples = "no: one sorted-microglia sample per mouse", match_units = c("CA1microglia", "CA2microglia", "CA3microglia", "DGmicroglia"),
    set = "primary", corrections = "none",
    notes = "limma on log2(normalised count + 1), RSD - CON, eBayes(trend, robust) with arrayWeights: two sorts contaminated by border macrophages/monocytes (CON 18) and neutrophils (RSD 8) are down-weighted. Excel date-mangled symbols restored. Whole-brain microglia, long recovery: small effects; the authors' RUVseq factors are not deposited."))
}

## --- GSE84183: Herve 2017, dentate gyrus, UCMS vs control (per-GSM raw arrays) ----------------------
build_herve2017 <- function() {
  sm <- read_series_matrix(raw_file("GSE84183", "GSE84183_series_matrix.txt.gz"))
  a <- sm$ann
  sup <- gsub('"', "", strsplit(sm$lines[startsWith(sm$lines, "!Sample_supplementary_file")], "\t")[[1]][-1])
  a$array <- sub(".*_(2528[0-9]+)_S01_GE1_105_Dec08_([0-9])_([0-9])[.]txt[.]gz$", "\\1_\\2_\\3", sup)
  a$campaign <- factor(a$`hybridization campaign`)
  untr <- a[a$flx == "No", ]
  stopifnot(all(untr$ucms %in% c("Stressed", "Unstressed")),
            all(grepl("-S-T_", untr$title) == (untr$ucms == "Stressed")), all(grepl("-NS-T_", untr$title) == (untr$ucms == "Unstressed")))
  s <- untr[untr$tissue == "Dentate Gyrus", ]
  stopifnot(nrow(s) == 16)
  GL <- readLines(raw_file("GSE84183", "GPL13912.txt"))
  gpl <- read.delim(text = GL[(which(GL == "!platform_table_begin") + 1):(which(GL == "!platform_table_end") - 1)],
                    quote = "", check.names = FALSE, colClasses = "character")
  gpl <- gpl[gpl$CONTROL_TYPE == "FALSE" & nzchar(gpl$NAME), ]
  pmap <- unique(gpl[, c("NAME", "GENE_ID", "GENE_SYMBOL")]); pmap <- pmap[!duplicated(pmap$NAME), ]
  read_fe <- function(f) {                                 # Agilent Feature Extraction file
    L <- readLines(gzfile(f))
    fp <- which(startsWith(L, "FEPARAMS")); fpv <- strsplit(L[fp + 1], "\t")[[1]]; names(fpv) <- strsplit(L[fp], "\t")[[1]]
    d <- read.delim(text = L[which(startsWith(L, "FEATURES")):length(L)], check.names = FALSE, quote = "", colClasses = "character")
    d <- d[d[[1]] == "DATA", ]
    list(barcode = unname(fpv["FeatureExtractor_Barcode"]),
         d = data.frame(ProbeName = d$ProbeName, ControlType = as.integer(d$ControlType), sig = as.numeric(d$gProcessedSignal),
                        wab = as.integer(d$gIsWellAboveBG), stringsAsFactors = FALSE))
  }
  fe <- lapply(s$gsm, function(g) read_fe(raw_file("GSE84183", paste0(g, ".txt.gz"))))
  stopifnot(all(sapply(fe, `[[`, "barcode") == s$array))  # raw file = array named in the GSM record
  d0 <- fe[[1]]$d; stopifnot(all(sapply(fe, function(z) identical(z$d$ProbeName, d0$ProbeName))))
  nc <- d0$ControlType == 0
  E <- normalizeBetweenArrays(sapply(fe, function(z) log2(pmax(z$d$sig[nc], 1))), method = "quantile")
  W <- sapply(fe, function(z) z$d$wab[nc])
  keep <- rowSums(W == 1) >= 8                              # well above background in >= 8 of 16 arrays
  E <- avereps(E[keep, ], ID = d0$ProbeName[nc][keep])
  grp <- factor(ifelse(s$ucms == "Stressed", "UCMS", "CON"), levels = c("CON", "UCMS"))
  design <- model.matrix(~ campaign + grp, data = data.frame(grp = grp, campaign = s$campaign))
  eb <- eBayes(lmFit(E, design), trend = TRUE, robust = TRUE)
  pm <- pmap[match(rownames(E), pmap$NAME), ]
  s1 <- entrez_symbol(pm$GENE_ID); s2 <- mgi_symbol(pm$GENE_SYMBOL)
  st <- collapse_by_symbol(cbind(stressed_vs_CON = eb$t[, "grpUCMS"]), ifelse(!is.na(s1), s1, s2), rowMeans(E))
  list(sig("Herve2017_DG", st, study = "Herve et al. 2017", pmid_or_doi = "PMID:28848385",
    accession = "GSE84183", species = "Mus musculus (BALB/c)", sex = "male",
    region = "dentate gyrus (microdissected; variable admixture of neighbouring tissue)",
    assay = "Agilent SurePrint G3 Mouse GE 8x60K (GPL13912), per-GSM Feature Extraction gProcessedSignal, log2, quantile-normalised over the 16 DG arrays",
    model = "unpredictable chronic mild stress (UCMS) 8-9 weeks, single-housed; group-housed controls",
    timing = "end of UCMS", stat_type = "moderated_t", n_CON = 8L, n_RES = NA, n_SUS = NA, n_stressed = 8L,
    pooled_samples = "no: one array per mouse", match_units = c("DGmo", "DGpo", "DGsg", "DGmicroglia"), set = "primary",
    corrections = "GSE84183: built from the per-GSM raw Feature Extraction files because the series-matrix data columns are permuted relative to the sample records (series column k holds the k-th array in barcode order).",
    notes = "limma ~ hybridization campaign + group, eBayes(trend, robust); fluoxetine arms dropped. Housing differs between groups (part of the contrast); group is confounded with array position on every slide (CON 1_1/1_2, UCMS 1_3/1_4); DG dissection purity varies."))
}

## --- GSE172451: Sturman 2021, ventral hippocampus, adolescent chronic social instability --------------
build_sturman2021 <- function() {
  x <- read.delim(gzfile(raw_file("GSE172451", "GSE172451_CSI2_genecounts.txt.gz")), check.names = FALSE, row.names = 1)
  cnt <- round(x[, grepl("Homecage", colnames(x))])         # home-cage arm only
  grp <- ifelse(grepl("\\.CSI\\.", colnames(cnt)), "CSI", "CON")
  g <- factor(grp); des <- model.matrix(~ 0 + g); colnames(des) <- levels(g)
  st <- voom_mor(cnt, des, c(stressed_vs_CON = "CSI-CON"), min_samples = min(table(grp)))$t
  list(sig("Sturman2021_vHC", st, study = "Sturman et al. 2021", pmid_or_doi = "PMID:34527792", accession = "GSE172451",
    species = "Mus musculus", sex = "male", region = "ventral hippocampus (bulk)", assay = "bulk RNA-seq, gene counts",
    model = "adolescent chronic social instability, 7 weeks", timing = "adulthood, weeks after the end of stress (home-cage arm)",
    stat_type = "moderated_t", n_CON = 24L, n_RES = NA, n_SUS = NA, n_stressed = 24L, pooled_samples = "no: individual mice",
    match_units = "all", set = "primary", corrections = "none",
    notes = "Home-cage samples only (acute-restraint arm excluded); no SUS/RES split. voom with median-of-ratios size factors, robust eBayes; gene names as deposited (no symbol remapping)."))
}

####################################################################################################
## Run
####################################################################################################

BUILDERS <- list(
  Schmidt2010 = list(f = build_schmidt2010, ids = c("Schmidt2010_dCA1", "Schmidt2010_dDG", "SENS_Schmidt2010_dCA1_asGEO")),
  Lehmann2018 = list(f = build_lehmann2018, ids = "Lehmann2018_microglia"),
  Henningsen2012 = list(f = build_henningsen2012, ids = "Henningsen2012_CA"),
  Bagot2016 = list(f = build_bagot2016, ids = c("Bagot2016_vHIP_28d", "Bagot2016_vHIP_48h")),
  Bagot2017 = list(f = build_bagot2017, ids = "Bagot2017_vHIP"),
  Laine2018 = list(f = build_laine2018, ids = "Laine2018_vHPC"),
  Ayash2023 = list(f = build_ayash2023, ids = "Ayash2023_vHC"),
  Butto2024 = list(f = build_butto2024, ids = "Butto2024_vHip_nuclei"),
  Bhatti2022 = list(f = build_bhatti2022, ids = "Bhatti2022_PV_TRAP"),
  Russo = list(f = build_russo_astro, ids = "Russo_astro_TRAP"),
  Lyons2025 = list(f = build_lyons2025, ids = c("Lyons2025_CA1", "Lyons2025_CA3")),
  Weber2019 = list(f = build_weber2019, ids = "Weber2019_microglia"),
  Herve2017 = list(f = build_herve2017, ids = "Herve2017_DG"),
  Sturman2021 = list(f = build_sturman2021, ids = "Sturman2021_vHC"))

INFO_COLS <- c("id", "study", "pmid_or_doi", "accession", "species", "sex", "region", "assay", "model", "timing", "stat_type",
               "n_CON", "n_RES", "n_SUS", "n_stressed", "pooled_samples", "contrasts", "match_units", "set", "corrections",
               "notes", "sha256")

write_signature <- function(s) {
  st <- s$stat
  stopifnot(is.matrix(st), is.numeric(st), !anyNA(st), !anyNA(rownames(st)), !anyDuplicated(rownames(st)),
            all(colnames(st) %in% names(CONTRASTS)))
  st <- st[order(rownames(st), method = "radix"), names(CONTRASTS)[names(CONTRASTS) %in% colnames(st)], drop = FALSE]
  f <- file.path(OUT, paste0(s$id, ".csv"))
  fwrite(data.table(gene = rownames(st), st), f, eol = "\n")
  i <- s$info
  i$contrasts <- paste(colnames(st), collapse = " ")
  i$match_units <- paste(i$match_units, collapse = " ")
  i$sha256 <- digest::digest(file = f, algo = "sha256")
  as.data.table(i[INFO_COLS])
}

want <- commandArgs(trailingOnly = TRUE)
all_ids <- unlist(lapply(BUILDERS, `[[`, "ids"))
if (length(bad <- setdiff(want, all_ids))) stop("unknown id(s): ", paste(bad, collapse = ", "))
info_file <- file.path(root, "signature_info.csv")
info <- if (file.exists(info_file)) fread(info_file, colClasses = "character") else NULL
t_all <- Sys.time()
for (b in names(BUILDERS)) {
  if (length(want) && !any(BUILDERS[[b]]$ids %in% want)) next
  t0 <- Sys.time()
  sigs <- BUILDERS[[b]]$f()
  stopifnot(setequal(vapply(sigs, `[[`, "", "id"), BUILDERS[[b]]$ids))
  rows <- rbindlist(lapply(sigs, write_signature))
  rows <- rows[, lapply(.SD, as.character)]
  info <- rbindlist(list(if (!is.null(info)) info[!info$id %in% rows$id], rows), fill = TRUE)
  cat(sprintf("%-15s %s (%.0f s)\n", b, paste(sprintf("%s %d genes", rows$id, sapply(sigs, function(s) nrow(s$stat))), collapse = "; "),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
info <- info[match(intersect(all_ids, info$id), info$id)]
fwrite(info, info_file, eol = "\n")
cat(sprintf("done in %.1f min; %d signatures in signature_info.csv\n", as.numeric(difftime(Sys.time(), t_all, units = "mins")), nrow(info)))

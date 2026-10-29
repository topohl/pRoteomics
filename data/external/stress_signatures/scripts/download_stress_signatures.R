#!/usr/bin/env Rscript
# Download every raw file that scripts/build_stress_signatures.R needs.
#
#   Rscript scripts/download_stress_signatures.R           # fetch missing files and files whose MD5 differs
#   Rscript scripts/download_stress_signatures.R --check   # no network: report presence and MD5 only
#
# Target: <package>/raw/<ACCESSION>/<file>, or the folder in env var STRESS_SIG_RAW.
# A file already present with the MD5 recorded in download_manifest.csv is skipped (with no manifest yet, any present
# file is kept). Only https URLs are used. download_manifest.csv (dataset_id, accession, file, url, size_bytes, md5,
# note) is rewritten with the size and MD5 of the files on disk; with --check it is only written if it does not exist.
# A download whose MD5 differs from the manifest is kept but reported (the source changed; rebuild and re-check), and
# the script then exits with status 1.

args_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (!length(args_file)) stop("run with Rscript scripts/download_stress_signatures.R")
script_dir <- dirname(normalizePath(args_file))
root <- normalizePath(file.path(script_dir, ".."))
RAW <- Sys.getenv("STRESS_SIG_RAW", file.path(root, "raw"))
MANIFEST <- file.path(root, "download_manifest.csv")
CHECK_ONLY <- "--check" %in% commandArgs(trailingOnly = TRUE)
UA <- "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"   # the Elsevier CDN refuses requests without a browser User-Agent
options(timeout = max(3600, getOption("timeout")))

geo_series <- function(acc, sub, f) sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%snnn/%s/%s/%s", substr(acc, 1, nchar(acc) - 3), acc, sub, f)
geo_sample <- function(f) { g <- sub("_.*$", "", f)
  sprintf("https://ftp.ncbi.nlm.nih.gov/geo/samples/%snnn/%s/suppl/%s", substr(g, 1, nchar(g) - 3), g, f) }
row <- function(dataset_id, accession, file, url, note) data.frame(dataset_id, accession, file, url, note)

## --- the files ---------------------------------------------------------------------------------------
ayash <- c("GSM4913335_C1_HC_S45", "GSM4913336_C5_HC_S48", "GSM4913337_C6_HC_S51", "GSM4913338_C3_HC_S54", "GSM4913339_C4_HC_S102",
           "GSM4913340_C5_HC_S105", "GSM4913351_D9_HC_S36", "GSM4913352_D15_HC_S39", "GSM4913353_D23_HC_S42", "GSM4913354_D9_HC_S87",
           "GSM4913355_D10_HC_S90", "GSM4913356_D18_HC_S93", "GSM4913357_D27_HC_S96", "GSM4913370_D4_HC_S15", "GSM4913371_D12_HC_S18",
           "GSM4913372_D22_HC_S24", "GSM4913373_D1_HC_S66", "GSM4913374_D35_HC_S72", "GSM4913384_D16_HC_S27", "GSM4913385_D18_HC_S30",
           "GSM4913386_D14_HC_S33", "GSM4913387_D2_HC_S75", "GSM4913388_D31_HC_S78", "GSM4913389_D34_HC_S81", "GSM4913390_D33_HC_S84")
russo <- c("GSM4146188_1690_H_C_1_TAAGGCGAAGAGTAGA", "GSM4146189_1691_HIP_C1_CGTACTAGAGAGTAGA", "GSM4146190_1692_HIP_C3_AGGCAGAAAGAGTAGA",
           "GSM4146191_1693_H_RS1_TCCTGAGCGTAAGGAG", "GSM4146192_1694_H_RS_2_GGACTCCTGTAAGGAG", "GSM4146193_1695_HIP_R2_TAGGCATGGTAAGGAG",
           "GSM4146194_1696_HIP_SS1_CTCTCTACGTAAGGAG", "GSM4146195_1697_HIP_SS2_CAGAGAGGGTAAGGAG", "GSM4146196_1698_HIP_SS3_GCTACGCTACTGCATA",
           "GSM4146197_1699_H_SS_2_CGAGGCTGACTGCATA", "GSM4146198_1700_H_SS_3_AAGAGGCAACTGCATA")
russo <- paste0(russo, "_L001_paired.sortedByCoord.bam.counts.txt.gz")
# GSE84183 untreated dentate-gyrus arrays: local name <GSM>.txt.gz, remote name = Feature Extraction file
herve <- c(GSM2228639 = "Dent_Gyr-NS-T_39 252800513928_1_1", GSM2228640 = "Dent_Gyr-NS-T_40 252800513928_1_2",
           GSM2228641 = "Dent_Gyr-NS-T_52 252800520267_1_1", GSM2228642 = "Dent_Gyr-NS-T_53 252800520267_1_2",
           GSM2228643 = "Dent_Gyr-NS-T_54 252800520293_1_1", GSM2228644 = "Dent_Gyr-NS-T_64 252800520293_1_2",
           GSM2228645 = "Dent_Gyr-NS-T_65 252800520299_1_1", GSM2228646 = "Dent_Gyr-NS-T_66 252800520299_1_2",
           GSM2228647 = "Dent_Gyr-S-T_1 252800513928_1_3", GSM2228648 = "Dent_Gyr-S-T_11 252800513928_1_4",
           GSM2228649 = "Dent_Gyr-S-T_13 252800520267_1_3", GSM2228650 = "Dent_Gyr-S-T_15 252800520267_1_4",
           GSM2228651 = "Dent_Gyr-S-T_19 252800520293_1_3", GSM2228652 = "Dent_Gyr-S-T_25 252800520293_1_4",
           GSM2228653 = "Dent_Gyr-S-T_32 252800520299_1_3", GSM2228654 = "Dent_Gyr-S-T_35 252800520299_1_4")
herve_remote <- sprintf("%s_US83700202_%s_S01_GE1_105_Dec08_%s.txt.gz", names(herve), sub("^.* ([0-9]+)_.*$", "\\1", herve),
                        sub("^.* [0-9]+_", "", herve))

S10 <- "Schmidt2010_dCA1 Schmidt2010_dDG SENS_Schmidt2010_dCA1_asGEO"; B16 <- "Bagot2016_vHIP_28d Bagot2016_vHIP_48h"
L25 <- "Lyons2025_CA1 Lyons2025_CA3"
FILES <- rbind(
  row(S10, "GSE11211", "GSE11211_series_matrix.txt.gz", geo_series("GSE11211", "matrix", "GSE11211_series_matrix.txt.gz"), "sample annotation + BeadStudio values"),
  row(S10, "GSE11211", "GPL6753_Mouse_Ref-8_V1.csv.gz", "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6753/suppl/GPL6753_Mouse_Ref-8_V1.csv.gz", "Illumina probe annotation (RefSeq, symbol)"),
  row("Lehmann2018_microglia", "GSE107941", "GSE107941_series_matrix.txt.gz", geo_series("GSE107941", "matrix", "GSE107941_series_matrix.txt.gz"), "sample annotation + RMA values"),
  row("Lehmann2018_microglia", "GSE107941", "GSE107941_family.soft.gz", geo_series("GSE107941", "soft", "GSE107941_family.soft.gz"), "platform table (gene_assignment)"),
  row("Henningsen2012_CA", "PRD000546", "mmc1.zip", "https://ars.els-cdn.com/content/image/1-s2.0-S1535947620330218-mmc1.zip", "Mol Cell Proteomics supplement; Suppl_Table_1-3.xls = per-study iTRAQ ratios (PRIDE PRD000546 holds the raw spectra)"),
  row(B16, "GSE72343", "GSE72343_series_matrix.txt.gz", geo_series("GSE72343", "matrix", "GSE72343_series_matrix.txt.gz"), "sample annotation (titles, characteristics, SRX)"),
  row(B16, "GSE72343", "PRJNA293822_counts.csv", "https://ndownloader.figshare.com/files/21808590", "Flati 2020 Stress Mice Portal counts for PRJNA293822 (figshare collection 10.6084/m9.figshare.c.4860843); integer counts despite the 'normalized' file name"),
  row("Bagot2017_vHIP", "GSE81672", "GSE81672_HTSeq_counts_matrix.xls.gz", geo_series("GSE81672", "suppl", "GSE81672_HTSeq_counts_matrix.xls.gz"), "HTSeq counts, BIFF8 .xls (read by the build's own reader)"),
  row("Laine2018_vHPC", "GSE109315", "GSE109315-GPL13112_series_matrix.txt.gz", geo_series("GSE109315", "matrix", "GSE109315-GPL13112_series_matrix.txt.gz"), "sample annotation (vHPC samples)"),
  row("Laine2018_vHPC", "GSE109315", "GSE109315_vHPC_counts_logCPM.txt.gz", geo_series("GSE109315", "suppl", "GSE109315_vHPC_counts_logCPM.txt.gz"), "deposited logCPM"),
  row("Ayash2023_vHC", "GSE161726", "GSE161726_series_matrix.txt.gz", geo_series("GSE161726", "matrix", "GSE161726_series_matrix.txt.gz"), "sample annotation + per-GSM file names"),
  row("Ayash2023_vHC", "GSE161726", paste0(ayash, ".tsv.gz"), geo_sample(paste0(ayash, ".tsv.gz")), "featureCounts, one vHC sample"),
  row("Butto2024_vHip_nuclei", "GSE240573", sprintf("GSE240573_vHip_%s_Counts.txt.gz", c("Control", "Resilience", "Susceptible")),
      geo_series("GSE240573", "suppl", sprintf("GSE240573_vHip_%s_Counts.txt.gz", c("Control", "Resilience", "Susceptible"))), "HTSeq counts, one group"),
  row("Bhatti2022_PV_TRAP", "GSE184027", "GSE184027_rawcount.txt.gz", geo_series("GSE184027", "suppl", "GSE184027_rawcount.txt.gz"), "HTSeq counts (CSV)"),
  row("Bhatti2022_PV_TRAP", "GSE184027", "Table_2.xlsx", "https://www.frontiersin.org/api/v4/articles/898851/file/Table_2.xlsx/898851_supplementary-materials_tables_2_xlsx/1", "paper Supplementary Table 2 (sheet 'Original': counting-annotation symbols for reassigned Ensembl IDs)"),
  row("Russo_astro_TRAP", "GSE139684", "GSE139684_gsm.txt", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139684&targ=gsm&view=full&form=text", "GEO sample records (titles); dynamic GEO page, MD5 changes if GEO edits the records"),
  row("Russo_astro_TRAP", "GSE139684", russo, geo_sample(russo), "HTSeq counts, one hippocampus pool"),
  row(L25, "GSE278620", "GSE278620_Lyons_et_al_geomix_median_norm_q_datasets2.0.csv.gz", geo_series("GSE278620", "suppl", "GSE278620_Lyons_et_al_geomix_median_norm_q_datasets2.0.csv.gz"), "GeoMx median-normalised counts"),
  row(L25, "GSE278620", "GSE278620_series_matrix.txt.gz", geo_series("GSE278620", "matrix", "GSE278620_series_matrix.txt.gz"), "AOI annotation"),
  row("Weber2019_microglia", "GSE124169", "GSE124169_countmatrix.csv.gz", geo_series("GSE124169", "suppl", "GSE124169_countmatrix.csv.gz"), "DESeq2-normalised counts"),
  row("Weber2019_microglia", "GSE124169", "GSE124169_series_matrix.txt.gz", geo_series("GSE124169", "matrix", "GSE124169_series_matrix.txt.gz"), "sample annotation"),
  row("Herve2017_DG", "GSE84183", "GSE84183_series_matrix.txt.gz", geo_series("GSE84183", "matrix", "GSE84183_series_matrix.txt.gz"), "sample annotation only (its data columns are permuted)"),
  row("Herve2017_DG", "GSE84183", "GPL13912.txt", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13912&targ=self&form=text&view=data", "Agilent platform table (probe -> Entrez/symbol); dynamic GEO page"),
  row("Herve2017_DG", "GSE84183", paste0(names(herve), ".txt.gz"), geo_sample(herve_remote), paste("Agilent FE file,", herve)),
  row("Sturman2021_vHC", "GSE172451", "GSE172451_CSI2_genecounts.txt.gz", geo_series("GSE172451", "suppl", "GSE172451_CSI2_genecounts.txt.gz"), "gene counts"))
stopifnot(all(startsWith(FILES$url, "https://")), !anyDuplicated(paste(FILES$accession, FILES$file)))

## --- download / check ----------------------------------------------------------------------------------
md5 <- function(p) unname(tools::md5sum(p))
old <- if (file.exists(MANIFEST)) read.csv(MANIFEST, colClasses = "character") else NULL
expected <- if (is.null(old)) character() else setNames(old$md5, paste(old$accession, old$file))
FILES$size_bytes <- NA_character_; FILES$md5 <- NA_character_; status <- character(nrow(FILES))
for (i in seq_len(nrow(FILES))) {
  key <- paste(FILES$accession[i], FILES$file[i]); exp <- unname(expected[key])
  dest <- file.path(RAW, FILES$accession[i], FILES$file[i])
  have <- file.exists(dest)
  if (have && (is.na(exp) || md5(dest) == exp)) { status[i] <- "ok"
  } else if (CHECK_ONLY) { status[i] <- if (have) "MD5 MISMATCH" else "MISSING"
  } else {
    dir.create(dirname(dest), showWarnings = FALSE, recursive = TRUE)
    tmp <- paste0(dest, ".part")
    r <- tryCatch(utils::download.file(FILES$url[i], tmp, mode = "wb", method = "libcurl", quiet = TRUE,
                                       headers = c(`User-Agent` = UA)), error = function(e) e)
    if (inherits(r, "error") || r != 0 || !file.exists(tmp) || file.size(tmp) == 0) {
      unlink(tmp); status[i] <- "DOWNLOAD FAILED"
      message("download failed: ", FILES$url[i], if (inherits(r, "error")) paste(":", conditionMessage(r)))
    } else {
      file.rename(tmp, dest)
      status[i] <- if (is.na(exp)) "downloaded" else if (md5(dest) == exp) "downloaded (MD5 ok)" else "downloaded, MD5 CHANGED"
    }
  }
  if (file.exists(dest)) { FILES$size_bytes[i] <- sprintf("%.0f", file.size(dest)); FILES$md5[i] <- md5(dest) }
  if (grepl("CHANGED|MISMATCH", status[i])) FILES$note[i] <- paste0(FILES$note[i], "; MD5 differs from the previous manifest (", exp, ")")
  cat(sprintf("%-22s %-10s %s\n", status[i], FILES$accession[i], FILES$file[i]))
}
if (!CHECK_ONLY || is.null(old)) {
  con <- file(MANIFEST, "wb")                               # LF line endings on every platform
  write.csv(FILES[, c("dataset_id", "accession", "file", "url", "size_bytes", "md5", "note")], con, row.names = FALSE)
  close(con)
}
bad <- !status %in% c("ok", "downloaded", "downloaded (MD5 ok)")
cat(sprintf("%d files: %d ok, %d problems; raw folder %s\n", nrow(FILES), sum(!bad), sum(bad), RAW))
if (any(bad)) quit(status = 1)

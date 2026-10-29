#!/usr/bin/env Rscript
# Build the per-animal behaviour outcomes the proteomic analyses relate to, from the frozen behaviour bundles.
#
#   Rscript data/external/behaviour_outcomes/scripts/build_behaviour_outcomes.R        (from the repository root)
#
# Input : three frozen bundles on the lab share (roots below, or env vars), each verified against its 00_manifest.csv
#         and the manifest's pinned SHA-256:
#           MMMSociability behaviour bundle ebb_v104_20261026_4b651ac (outcome version ov4): A4_combz_animals.csv
#             (CombZ = the integrated stress score, ISS; group; sex; batch) and A2_prediction_animals.csv (Stage 09 early
#             movement, the first active phase after the first cage change);
#           MMMSociability figure-support bundle fsb_v10_20261030_b74873b: S2_animal_series.csv (first-night movement after
#             each of the four cage changes, post hoc run 29f);
#           SISanalyzer bundle ctp_v14_20261026_f146fe4: B_values.csv (corticosterone rise cort_delta, sucrose preference
#             suc_pref).
# Output: <package>/behaviour_outcomes_ov4.csv, one row per animal with an ISS (117), and <package>/outcome_info.csv
#         (the bundles, their manifest SHA-256, the output SHA-256). Values are copied; the only arithmetic is the mean
#         of the four first-night movements of an animal.
# Replaces, for new analyses, the April 2026 workbook copy (E9_Behavior_Data.xlsx, zScore sheet) and the retired GAMM
# movement scores (auc_individual_animals_*.csv) in data/external/behavior/.

suppressPackageStartupMessages(library(data.table))
args_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (!length(args_file)) stop("run with Rscript data/external/behaviour_outcomes/scripts/build_behaviour_outcomes.R")
root <- normalizePath(file.path(dirname(normalizePath(args_file)), ".."), winslash = "/")

SHARE <- Sys.getenv("EXP9_ANALYSIS_ROOT", "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis")
BUNDLES <- list(
  ebb = list(id = "ebb_v104_20261026_4b651ac", manifest_sha256 = "9fe3c1a68db9a1d569791d684d710cb25b7098bc08f6ce03c53409bb87195c42",
             dir = file.path(SHARE, "Behavior/RFID/analysis_ready/canonical/behavior_bundle")),
  fsb = list(id = "fsb_v10_20261030_b74873b", manifest_sha256 = "3bc9db48a2a75e79d3cbc1aa7d7a939c80240e74edd53c817fef63db30685007",
             dir = file.path(SHARE, "Behavior/RFID/analysis_ready/canonical/figure_support_bundle")),
  ctp = list(id = "ctp_v14_20261026_f146fe4", manifest_sha256 = "a773bbefea9431a175f1ef2dd7c70a6c2f70ef1cd9e905229ca1dc2aa2290afd",
             dir = file.path(SHARE, "SIS_Analysis/bundles")))
sha <- function(p) unname(tools::sha256sum(p))
read_verified <- function(b, file) {
  d <- file.path(b$dir, b$id); man <- file.path(d, "00_manifest.csv")
  if (!identical(sha(man), b$manifest_sha256)) stop(b$id, ": 00_manifest.csv SHA-256 differs from the pinned value")
  m <- fread(man, colClasses = "character")
  key <- intersect(c("file", "path", "relative_path"), names(m))[1]; hcol <- intersect(c("sha256", "SHA256"), names(m))[1]
  if (is.na(key) || is.na(hcol)) stop(b$id, ": manifest without file / sha256 columns")
  want <- m[[hcol]][m[[key]] == file]
  if (length(want) != 1L || !identical(sha(file.path(d, file)), want)) stop(b$id, "/", file, ": not the file its manifest lists")
  fread(file.path(d, file), colClasses = list(character = intersect(c("AnimalNum", "Animal"), names(fread(file.path(d, file), nrows = 0L)))))
}

a4 <- read_verified(BUNDLES$ebb, "A4_combz_animals.csv")
a2 <- read_verified(BUNDLES$ebb, "A2_prediction_animals.csv")
s2 <- read_verified(BUNDLES$fsb, "S2_animal_series.csv")
bv <- read_verified(BUNDLES$ctp, "B_values.csv")

if (nrow(a4) != 117L || anyDuplicated(a4$AnimalNum)) stop("A4 must hold 117 animals once each")
mv4 <- s2[, .(n_cage_changes = .N, movement_cc1_cc4_mean = mean(as.numeric(movement))), by = AnimalNum]
if (any(mv4$n_cage_changes != 4L)) stop("S2: an animal lacks a cage change")
phys <- dcast(bv[id %in% c("cort_delta", "suc_pref", "CombZ")], Animal ~ id, value.var = "value")
setnames(phys, "Animal", "AnimalNum")
# the SISanalyzer bundle zero-pads purely numeric ids (0003, 00654); the behaviour bundle does not (3, 654)
phys[, AnimalNum := ifelse(grepl("^[0-9]+$", AnimalNum), sub("^0+", "", AnimalNum), AnimalNum)]
out <- merge(a4[, .(AnimalNum, Sex, Batch, Group, CombZ, rfid_tracked)], a2[, .(AnimalNum, early_movement = Movement_mean)],
             by = "AnimalNum", all.x = TRUE)
out <- merge(out, mv4[, .(AnimalNum, movement_cc1_cc4_mean)], by = "AnimalNum", all.x = TRUE)
out <- merge(out, phys[, .(AnimalNum, cort_delta, suc_pref, CombZ_ctp = CombZ)], by = "AnimalNum", all.x = TRUE)
if (nrow(out) != 117L) stop("join changed the animal count")
d <- max(abs(out$CombZ - out$CombZ_ctp), na.rm = TRUE)
if (anyNA(out$CombZ_ctp) || d > 1e-9) stop("CombZ of the SISanalyzer bundle differs from the behaviour bundle (", d, ")")
if (sum(!is.na(out$early_movement)) != 111L) stop("early movement must exist for the 111 tracked animals")
out[, CombZ_ctp := NULL]
# the numeric suffix of the behaviour id is the animal id the proteomics uses (config/animal_id_aliases.csv)
out[, animal_id := sub("^0+", "", sub("^[A-Za-z]+", "", AnimalNum))]
if (anyDuplicated(out$animal_id)) stop("animal number suffixes are not unique")
setcolorder(out, c("animal_id", "AnimalNum", "Sex", "Batch", "Group", "CombZ", "cort_delta", "suc_pref", "early_movement",
                   "movement_cc1_cc4_mean", "rfid_tracked"))
setorder(out, AnimalNum)
f <- file.path(root, "behaviour_outcomes_ov4.csv")
fwrite(out, f)
info <- data.table(item = c(paste0("bundle_", names(BUNDLES)), "output"),
                   id = c(vapply(BUNDLES, `[[`, "", "id"), "behaviour_outcomes_ov4.csv"),
                   sha256 = c(vapply(BUNDLES, `[[`, "", "manifest_sha256"), sha(f)),
                   note = c("00_manifest.csv of the MMMSociability behaviour bundle (A4, A2)",
                            "00_manifest.csv of the MMMSociability figure-support bundle (S2)",
                            "00_manifest.csv of the SISanalyzer bundle (B_values)",
                            "117 animals; CombZ = ISS; early_movement = Stage 09 Movement_mean; movement_cc1_cc4_mean = mean of S2"))
fwrite(info, file.path(root, "outcome_info.csv"))
cat("wrote", f, "(", nrow(out), "animals ) sha256", sha(f), "\n")

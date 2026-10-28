# Map hippocampal-input source types (and negative controls / glia) onto the ABC Atlas WMB
# taxonomy, compute one cell-weighted mean profile per source type, and run marker sanity checks.
#
# Requires: ../raw (download_abc_wmb.R) and the outputs of build_abc_wmb_means.R in this folder.
# Writes (this folder):
#   source_type_map.csv                     source_type, level, abc_ids, abc_names, evidence, confidence
#   abc_wmb_source_type_mean_log2cpm.csv    genes x source types (cell-weighted over the mapped clusters)
#   abc_wmb_source_type_info.csv            n_clusters / n_cells per source type
#   sanity_check_markers.csv                marker means in target vs comparator source types
#
# Level: "subclass" or "supertype" wherever the taxonomy unit is anatomically clean. "cluster" is
# used only where the supertype mixes the source population with other cells (EC layer II, where
# supertypes mix MEC, LEC and piriform L2; raphe, where serotonergic and Vglut3-only clusters share
# supertypes; locus coeruleus, whose supertype is 73% non-noradrenergic cells).
#
# Usage: Rscript build_source_type_map.R

suppressPackageStartupMessages({
  library(rhdf5)
  library(data.table)
  library(jsonlite)
})
setDTthreads(4)

args_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
here <- if (length(args_file)) dirname(normalizePath(args_file)) else getwd()
out_dir <- if (basename(here) == "scripts") dirname(here) else here   # the reference folder (data/external/allen_abc_wmb_10x)
raw <- normalizePath(Sys.getenv("ABC_WMB_RAW", file.path(out_dir, "raw")), winslash = "/")   # where download_abc_wmb.R puts them
h5_path <- file.path(raw, "mapmycells/WMB-taxonomy/20240831/precomputed_stats_ABC_revision_230821.h5")

tax <- fread(file.path(out_dir, "abc_wmb_taxonomy.csv"))
genes <- fread(file.path(out_dir, "abc_wmb_genes.csv"))

# ---- the map -------------------------------------------------------------------------------
# ids are the leading numbers of the official names (subclass "037", supertype "0100",
# cluster "0399"); `why` is the curated part of the evidence, `markers` are reported with values.
m <- function(source_type, level, ids, markers, confidence, why) {
  data.table(source_type, level, ids = list(ids), markers = list(markers), confidence, why)
}
map <- rbindlist(list(
  m("EC layer II - MEC stellate (ocean) cells", "cluster", c("0187", "0191"),
    c("Reln", "Calb1", "Wfs1", "Slc17a7"), "medium",
    "Reln-high/Calb1-low clusters of L2 IT ENT-po Glut located in ENTm; Yao 2021 cluster matches 149/150_L2 IT ENTm. Reln+ Calb1- is the stellate (ocean cell) signature. Supertypes 0051/0054 also contain LEC (0186) cells, hence cluster level."),
  m("EC layer II - MEC pyramidal island cells", "supertype", "0084",
    c("Wfs1", "Calb1", "Reln", "Slc17a7"), "medium",
    "Wfs1-high/Calb1-high/Reln-negative L2 cells (island-cell signature, Wfs1/Calb1) in ENTm2 and adjacent PAR by MERFISH; dissected from ENT; Yao 2021 cluster match '145_L2 IT PAR'. Taxonomy names it L2/3 IT PPP (parasubiculum-like), so the MEC-island vs PAR L2 split is not resolved."),
  m("EC layer II - LEC fan cells", "cluster", c("0155", "0158", "0186"),
    c("Reln", "Calb1", "Wfs1", "Slc17a7"), "medium",
    "Reln-high/Calb1-low clusters located in ENTl2 (MERFISH) and dissected from ENT; Yao 2021 cluster matches 151_L2 IT ENTl / 146_L2 IT ENTm. Supertype 0042 also holds piriform L2 clusters (0156, 0157, 0160), hence cluster level."),
  m("EC layer II - all MEC+LEC L2 excitatory", "cluster",
    c("0186", "0187", "0188", "0189", "0190", "0191", "0329", "0330", "0331", "0155", "0158", "0136", "0137"),
    c("Reln", "Calb1", "Wfs1", "Slc17a7"), "medium",
    "Union of the stellate, island and fan clusters plus the mixed Reln+/Calb1+ ENT-po clusters (0188-0190) and L2/3 ENTl clusters 0136/0137 (supertype 0038, MERFISH ENTl2). Use when the layer II input is treated as one population."),
  m("EC layer III - MEC", "supertype", "0036",
    c("Reln", "Calb1", "Pcp4", "Cux2", "Trps1", "Slc17a7"), "medium",
    "L2/3 IT ENT Glut_4: dissected from ENT, MERFISH ENTm3/ENTm5/ENTm2, curated 'L3 ENT po', Yao 2021 cluster matches 135/138_L3 IT ENTm (0131: 139_L3 IT ENTl). Reln- and Calb1-low, unlike layer II."),
  m("EC layer III - LEC", "supertype", c("0033", "0034", "0035", "0037"),
    c("Reln", "Calb1", "Pcp4", "Cux2", "Trps1", "Slc17a7"), "medium",
    "L2/3 IT ENT Glut_1/2/3/5: dissected from ENT, MERFISH ENTl3 (with ENTl5/ENTl2/TR), curated 'L3 ENTl' / 'L3 TR-ENTl' / 'L2/3 ENTl', Yao 2021 'L2/3 IT ENTl' / 'L3 IT ENT'. 0035 also has ENTl2 cells (23%)."),
  m("CA1 pyramidal cells", "supertype", c("0069", "0072"),
    c("Fibcd1", "Wfs1", "Rgs14", "Pcp4", "Amigo2", "Prox1"), "high",
    "CA1-ProS Glut_1 (curated 'CA1 do') and _4: >80% dissected from HIP, MERFISH CA1so/sr/sp, Yao 2021 'CA1'. Fibcd1/Wfs1-high. Excluded: 0070/0071/0074 (prosubiculum-rich, RHP-dissected) and 0073 (ventral CA1 bordering PA/HATA)."),
  m("CA2 pyramidal cells", "supertype", "0100",
    c("Pcp4", "Rgs14", "Amigo2", "Map3k15", "Fibcd1", "Prox1"), "high",
    "CA2-FC-IG Glut_1 (single cluster 0399, n=143): curated 'CA2', 100% HIP-dissected, Yao 2021 cluster match 360_CA2-IG-FC_2, MERFISH registration spreads it over the CA1/CA2/CA3 border strata (CCF CA2 is thin); Pcp4/Rgs14/Amigo2/Map3k15 all high (canonical CA2 markers). Sister supertype 0101 is indusium griseum/fasciola cinerea and is excluded."),
  m("CA3 pyramidal cells", "supertype", c("0075", "0077", "0078"),
    c("Bok", "Iyd", "Npy2r", "Calb2", "Drd2", "Prox1"), "high",
    "CA3 Glut_1/3/4: >95% HIP-dissected, MERFISH CA3sp/so/sr/slu, curated 'CA3' / 'CA3 ve', Yao 2021 'CA3'. Bok-high, Prox1-low. Excluded: 0076 (ventral-most CA3 mixed with posterior amygdala, MERFISH PA) and 0079 (mossy-cell candidate)."),
  m("Dentate granule cells", "subclass", "037",
    c("Prox1", "C1ql2", "Dock10", "Calb2", "Slc17a7"), "high",
    "DG Glut: 97% HIP-dissected, MERFISH DG-mo/sg/po, curated 'DG do' / 'DG ve' / 'DG po', Yao 2021 'DG'; Prox1/C1ql2-high. Immature granule cells are a separate subclass (038 DG-PIR Ex IMN, not included)."),
  m("Hilar mossy cells", "supertype", "0079",
    c("Calb2", "Drd2", "Prox1", "C1ql2", "Slc17a7"), "medium",
    "CA3 Glut_5 (clusters 0316/0317, n=318): Calb2-high and Drd2+ (established mossy-cell markers; Drd2 is ~0 in all other CA3/DG supertypes), Prox1/C1ql2-low (not granule cells), MERFISH CA3sr/sp plus DG-po (hilus) and DG-sg, curated 'CA3 ve'. Calb2+ fits ventral mossy cells; dorsal (Calb2-low) mossy cells are not resolved as a separate type in this taxonomy."),
  m("Hippocampal Sst (O-LM) interneurons", "supertype", c("0216", "0219"),
    c("Sst", "Chrna2", "Reln", "Grm1", "Elfn1", "Npy", "Pvalb"), "medium",
    "Sst Gaba supertypes with >=20% of cells HIP-dissected (0219: 82%, 0216: 26%; together 87% of the HIP-dissected Sst Gaba cells). Sst/Reln/Grm1/Elfn1-high (O-LM-like); MERFISH CA1so (0216) and CA3/DG (0219, includes hilar HIPP-like cells). Chrna2 (O-LM marker) is detected only weakly (<=2.6 log2CPM in any cluster). 0216 also holds subiculum/prosubiculum clusters."),
  m("Hippocampal Pvalb interneurons", "subclass", "052",
    c("Pvalb", "Lhx6", "Tac1", "Satb1", "Sst", "Slc32a1"), "medium",
    "Pvalb Gaba (basket/bistratified). Only 252 of its cells are HIP-dissected and they are spread over cortical-like supertypes (0206 53%, 0212 19%, 0208 10%), so the whole subclass is used; the only HIP-enriched supertype, 0206, is 65% amygdala/piriform cells (cluster 0738). Axo-axonic cells are the separate subclass 051 Pvalb chandelier Gaba (not included)."),
  m("Hippocampal Lamp5/neurogliaform interneurons", "supertype", c("0203", "0195"),
    c("Lamp5", "Ndnf", "Lhx6", "Nos1", "Npy", "Reln"), "medium",
    "0203 Lamp5 Lhx6 Gaba_1 (= subclass 050; 54% HIP-dissected, MERFISH DG-mo/CA strata; Lamp5/Lhx6/Nos1/Npy-high = MGE-derived ivy/neurogliaform) plus 0195 RHP-COA Ndnf Gaba_3 (Ndnf/Reln-high CGE neurogliaform-like, 44% HIP, MERFISH CA1sr/slm). Weighted mean is dominated by 0203 (8913 vs 695 cells). Cortical Lamp5 Gaba (049) excluded (4% HIP)."),
  m("Hippocampal Vip interneurons", "supertype", c("0179", "0181"),
    c("Vip", "Calb2", "Cck", "Cnr1", "Nos1", "Penk"), "medium",
    "Vip Gaba_7 (83% HIP) and _9 (56% HIP): the HIP-enriched Vip supertypes (42% of HIP-dissected Vip cells); Vip/Calb2-high = interneuron-selective (IS) cell profile. The remaining HIP Vip cells fall in cortical-like supertypes (0177, 0174, 0173); subclass 046 is the broader alternative."),
  m("Hippocampal Sncg (CCK) interneurons", "supertype", c("0189", "0193", "0196", "0197"),
    c("Sncg", "Cck", "Cnr1", "Vip", "Ndnf", "Slc32a1"), "low",
    "CCK/Cnr1-high supertypes with >=25% HIP-dissected cells: Sncg Gaba_5 (0189, 27% HIP) and RHP-COA Ndnf Gaba_1/4/5 (85-94% HIP), whose clusters Yao 2021 matched to 'Sncg'. Sncg itself is lower in the 048 supertypes (0.6-2.3) than in 047 (4.8-6.1). Alternative: subclass 047 Sncg Gaba (cortex-dominated)."),
  m("Nucleus reuniens glutamatergic", "supertype", "0670",
    c("Slc17a6", "Slc17a7", "Nox4", "Calb2", "Calb1", "Gad2"), "high",
    "RE-Xi Nox4 Glut_3: curated 'RE' / 'RE an', MERFISH RE 0.66 (PT 0.14), 87% TH-dissected; Slc17a6-high. Sister supertypes 0668/0671 (curated xiphoid, MERFISH RE+Xi) and 0669 (PR/SMT) excluded; subclass 152 is the broader alternative."),
  m("Supramammillary glutamatergic", "supertype", c("0622", "0625", "0626", "0685"),
    c("Slc17a6", "Slc32a1", "Calb2", "Nos1", "Pitx2", "Lmx1a"), "medium",
    "PH-SUM Foxa1 Glut_1/_4/_5 (MERFISH SUM 0.22-0.34, curated 'SUM'/'PH') and IF-RL-CLI-PAG Foxa1 Glut_1 (MERFISH SUM 0.43, curated 'right above SUM'); Slc17a6/Calb2/Nos1-high. 0626 contains Slc32a1+ Glut-GABA clusters (2555, 2556) matching the SuM->DG glutamate/GABA co-releasing neurons. All mix with posterior hypothalamus/mammillary cells."),
  m("MS/DB cholinergic", "supertype", "0259",
    c("Chat", "Slc5a7", "Slc18a3", "Slc32a1", "Lhx8", "Isl1"), "high",
    "PAL-STR Gaba-Chol_1: curated 'MS NDB', MERFISH NDB 0.33 / MS 0.15, PAL-dissected, Chat/Slc5a7/Slc18a3-high (co-express Gad2/Slc32a1). Sister 0260 = striatal cholinergic interneurons (CP/ACB/OT) and 0261 = Chat-negative; both excluded."),
  m("MS/DB GABAergic", "supertype", c("0245", "0246"),
    c("Gad2", "Slc32a1", "Lhx6", "Lhx8", "Pvalb", "Kcnc1"), "medium",
    "NDB-SI-MA-STRv Lhx8 Gaba_2 (curated 'MS', MERFISH MS 0.33) and _3 (curated 'NDB MS', MERFISH NDB 0.53 / MS 0.29); GABAergic, Lhx6/Lhx8+, Kcnc1-high. Pvalb is only moderately detected (cluster 0884: 4.2), so the PV-pacemaker septohippocampal subset is not cleanly separated."),
  m("MS/DB glutamatergic", "supertype", "0518",
    c("Slc17a6", "Slc32a1", "Gad2", "Isl1", "Calb2", "Pvalb"), "low",
    "MS-SF Bsx Glut_1 (= subclass 115): named for MS, MERFISH LSr 0.32 / MS 0.22 / NDB 0.13, Slc17a6-high, Slc32a1-low (Gad2+ without Vgat). Adjacent AVPV-MEPO-SFO Tbr1 Glut_1 (0519, MERFISH MS 0.39 but curated MEPO/MPO) is a possible alternative; the hippocampus-projecting Vglut2 septal population has no dedicated type."),
  m("Dorsal raphe serotonergic", "cluster", c("3880", "3881", "3882", "3883", "3884", "3885", "3896"),
    c("Tph2", "Slc6a4", "Fev", "Slc17a8", "En1", "Gad2"), "medium",
    "Tph2/Slc6a4-high clusters of MB-MY Tph2 Glut-Sero, midbrain-dissected, MERFISH DR/PAG (curated 'RAmb' / 'DR'). Supertype level is unsuitable: e.g. 0890 ('DR') and 0889 are mostly Tph2-negative Vglut3 (Slc17a8) clusters."),
  m("Median raphe serotonergic", "cluster", c("3886", "3903", "3904", "3911"),
    c("Tph2", "Slc6a4", "Fev", "Slc17a8", "En1", "Gad2"), "low",
    "Tph2/Slc6a4-high clusters dissected from pons with MERFISH CS (superior central nucleus = MR) / RPO / medial lemniscus and curated 'CS ve RPO'. 3910 (MERFISH DTN/PCG/CS) and 3902 (MRN, 'RAmb') are ambiguous between DR and MR and are excluded."),
  m("Locus coeruleus noradrenergic", "cluster", "4548",
    c("Dbh", "Th", "Slc6a2", "Slc18a2", "Phox2a", "Gal"), "medium",
    "Single cluster (n=22) of NTS Dbh Glut_4 curated 'LC SLC', pons-dissected, MERFISH PB/MV/SLC/LDT; Dbh/Th/Slc6a2/Slc18a2-high. The supertype (1007) adds cluster 4547 (n=59, Dbh 3.3, Th 0.7), so it is not used. LC is heavily under-sampled: the profile rests on 22 cells."),
  m("Negative control - striatal D1 MSN", "subclass", "061",
    c("Drd1", "Drd2", "Tac1", "Penk", "Ppp1r1b"), "high",
    "STR D1 Gaba: MERFISH CP/ACB/OT, STR-dissected, Drd1/Tac1-high."),
  m("Negative control - striatal D2 MSN", "subclass", "062",
    c("Drd1", "Drd2", "Tac1", "Penk", "Ppp1r1b"), "high",
    "STR D2 Gaba: MERFISH CP/ACB, STR-dissected, Drd2/Penk-high."),
  m("Negative control - cerebellar granule cells", "subclass", "314",
    c("Gabra6", "Cbln3", "Neurod1", "Slc17a7"), "high",
    "CB Granule Glut: CB-dissected, MERFISH cerebellar lobules, Gabra6-high."),
  m("Negative control - olfactory bulb granule cells", "subclass", "041",
    c("Gad1", "Slc32a1", "Frmd7", "Dlx1", "Prokr2"), "medium",
    "OB-in Frmd7 Gaba (inner granule cells): OLF-dissected, MERFISH MOBgr/MOBopl. Sister subclasses 042 OB-out (superficial/periglomerular-like) and 043 OB-mi Frmd7 Gaba are not included."),
  m("Astrocytes", "subclass", "319",
    c("Aqp4", "Gfap", "Aldh1l1", "Slc1a3", "Slc1a2"), "high",
    "Astro-TE NN (telencephalic astrocytes, which include hippocampal astrocytes; 11% of its cells are HIP-dissected). Non-telencephalic astrocytes are 318 Astro-NT NN."),
  m("Oligodendrocytes", "subclass", "327", c("Mog", "Plp1", "Mbp", "Olig2"), "high", "Oligo NN."),
  m("OPC", "subclass", "326", c("Pdgfra", "Cspg4", "Olig2"), "high", "OPC NN."),
  m("Microglia", "subclass", "334", c("P2ry12", "Tmem119", "Cx3cr1", "Csf1r"), "high",
    "Microglia NN (border-associated macrophages are 335 BAM NN, not included)."),
  m("Endothelial cells", "subclass", "333", c("Cldn5", "Flt1", "Pecam1"), "high", "Endo NN."),
  m("Pericytes", "subclass", "331", c("Kcnj8", "Vtn", "Pdgfrb", "Abcc9"), "high", "Peri NN."),
  m("VLMC", "subclass", "330", c("Dcn", "Col1a1", "Slc6a13", "Ptgds"), "high", "VLMC NN.")
))

# ---- resolve ids to clusters -------------------------------------------------------------------
name_col <- c(subclass = "subclass_name", supertype = "supertype_name", cluster = "cluster_name")
label_col <- c(subclass = "subclass_label", supertype = "supertype_label", cluster = "cluster_label")
lead_id <- function(x) sub(" .*", "", x)
resolve <- function(level, ids) {
  hit <- tax[lead_id(get(name_col[[level]])) %in% ids]
  found <- unique(lead_id(hit[[name_col[[level]]]]))
  if (!setequal(found, ids)) stop("unresolved ", level, " ids: ", paste(setdiff(ids, found), collapse = ","))
  hit
}

# ---- cluster means for every mapped cluster -----------------------------------------------------
c2r <- unlist(fromJSON(h5read(h5_path, "cluster_to_row")))
all_cl <- unique(unlist(lapply(seq_len(nrow(map)), function(i)
  resolve(map$level[i], map$ids[[i]])$cluster_label)))
rows <- sort(unname(c2r[all_cl]) + 1)
S <- h5read(h5_path, "sum", index = list(NULL, rows))
h5closeAll()
colnames(S) <- names(c2r)[match(rows - 1, c2r)]
stopifnot(setequal(colnames(S), all_cl))

# dissection region per cluster (for evidence strings)
cells <- fread(file.path(raw, "metadata/WMB-10X/20241115/cell_metadata.csv"),
               select = c("cluster_alias", "region_of_interest_acronym"))
roi <- cells[, .N, by = .(cluster_alias, roi = region_of_interest_acronym)]
rm(cells); invisible(gc())

fmt_freq <- function(names, n, cover = 0.8, min_frac = 0.03, max_n = 6) {
  ok <- !is.na(n) & n > 0
  if (!any(ok)) return(NA_character_)
  names <- names[ok]; n <- n[ok]
  o <- order(-n); fr <- n[o] / sum(n); nm <- names[o]
  keep <- seq_len(min(max_n, min(which(cumsum(fr) >= cover - 1e-9))))
  keep <- keep[fr[keep] >= min_frac | keep == 1]
  paste0(nm[keep], ":", sprintf("%.2f", fr[keep]), collapse = ",")
}
parse_freq <- function(x) {
  if (is.na(x) || !nzchar(x)) return(data.table(region = character(), frac = numeric()))
  parts <- strsplit(strsplit(x, ",(?! )", perl = TRUE)[[1]], ":", fixed = TRUE)
  data.table(region = vapply(parts, `[`, "", 1), frac = as.numeric(vapply(parts, `[`, "", 2)))
}
gidx <- function(sym) match(sym, genes$gene_symbol)   # first match for duplicated symbols

weighted_profile <- function(cl) {
  w <- cl$n_cells
  drop(S[, cl$cluster_label, drop = FALSE] %*% w) / sum(w)
}

# ---- evidence, map table, source-type profiles ---------------------------------------------------
out_map <- vector("list", nrow(map))
prof <- matrix(NA_real_, nrow(S), nrow(map), dimnames = list(NULL, map$source_type))
info <- vector("list", nrow(map))
for (i in seq_len(nrow(map))) {
  lev <- map$level[i]
  cl <- resolve(lev, map$ids[[i]])
  units <- unique(cl[, c(label_col[[lev]], name_col[[lev]]), with = FALSE])
  units <- units[order(get(name_col[[lev]]))]
  p <- weighted_profile(cl)
  prof[, i] <- p
  r <- roi[cluster_alias %in% cl$cluster_alias, .(N = sum(N)), by = roi]
  mer <- cl[, parse_freq(merfish_ccf_fine_freq), by = .(cluster_alias, n_cells)][
    , .(N = sum(frac * n_cells)), by = region]
  mk <- map$markers[[i]]
  mk_txt <- paste(sprintf("%s=%.1f", mk, p[gidx(mk)]), collapse = ", ")
  ann <- paste(unique(na.omit(cl$anatomical_annotation)), collapse = "; ")
  # Yao 2021 (cortex + hippocampus taxonomy) matches are only meaningful for pallial
  # neuron classes and non-neuronal cells; elsewhere they are forced nearest matches.
  y21_ok <- all(substr(cl$class_name, 1, 2) %in% c("01", "02", "04", "06", "07", "30", "31", "33", "34"))
  y21 <- if (y21_ok) paste(unique(na.omit(cl$yao2021_ctx_hpf_subclass[nzchar(cl$yao2021_ctx_hpf_subclass)])),
                           collapse = "; ") else ""
  nt <- unique(ifelse(is.na(cl$neurotransmitter) | !nzchar(cl$neurotransmitter), "NA", cl$neurotransmitter))
  parents <- if (lev == "cluster") paste(unique(cl$supertype_name), collapse = "; ") else
    if (lev == "supertype") paste(unique(cl$subclass_name), collapse = "; ") else
      paste(unique(cl$class_name), collapse = "; ")
  evidence <- paste0(
    map$why[i],
    " | parent: ", parents,
    " | n_clusters=", nrow(cl), ", n_cells=", sum(cl$n_cells),
    " (10Xv3 ", sprintf("%.0f%%", 100 * sum(cl$n_cells_10Xv3) / sum(cl$n_cells)), ")",
    " | neurotransmitter: ", paste(nt, collapse = "/"),
    " | dissection: ", fmt_freq(r$roi, r$N),
    " | MERFISH CCF (approx.): ", fmt_freq(mer$region, mer$N),
    if (nzchar(ann)) paste0(" | curated anatomy: ", ann) else "",
    if (nzchar(y21)) paste0(" | Yao 2021 subclass match: ", y21) else "",
    " | markers (mean log2CPM+1): ", mk_txt)
  out_map[[i]] <- data.table(
    source_type = map$source_type[i], level = lev,
    abc_ids = paste(units[[1]], collapse = ";"),
    abc_names = paste(units[[2]], collapse = ";"),
    evidence = evidence, confidence = map$confidence[i])
  info[[i]] <- data.table(source_type = map$source_type[i], level = lev,
                          n_units = nrow(units), n_clusters = nrow(cl), n_cells = sum(cl$n_cells),
                          n_cells_10Xv3 = sum(cl$n_cells_10Xv3), n_cells_10Xv2 = sum(cl$n_cells_10Xv2))
}
fwrite(rbindlist(out_map), file.path(out_dir, "source_type_map.csv"))
fwrite(rbindlist(info), file.path(out_dir, "abc_wmb_source_type_info.csv"))
fwrite(cbind(genes[, .(gene_symbol, gene_identifier)], as.data.table(round(prof, 3))),
       file.path(out_dir, "abc_wmb_source_type_mean_log2cpm.csv"))

# ---- sanity checks ------------------------------------------------------------------------------
sub_means <- fread(file.path(out_dir, "abc_wmb_subclass_mean_log2cpm.csv"))
sub_mat <- as.matrix(sub_means[, -(1:2)])
checks <- list(
  list("Dentate granule cells", c("Prox1", "C1ql2"),
       c("CA3 pyramidal cells", "CA1 pyramidal cells", "Hilar mossy cells")),
  list("CA2 pyramidal cells", c("Pcp4", "Rgs14", "Amigo2"),
       c("CA1 pyramidal cells", "CA3 pyramidal cells", "Dentate granule cells")),
  list("Hilar mossy cells", "Calb2",
       c("Dentate granule cells", "CA3 pyramidal cells", "CA1 pyramidal cells")),
  list("EC layer II - MEC stellate (ocean) cells", c("Reln", "Calb1"),
       c("EC layer II - MEC pyramidal island cells", "EC layer II - LEC fan cells", "EC layer III - MEC")),
  list("EC layer II - MEC pyramidal island cells", c("Reln", "Calb1"),
       c("EC layer II - MEC stellate (ocean) cells", "EC layer III - MEC")),
  list("EC layer II - LEC fan cells", c("Reln", "Calb1"),
       c("EC layer III - LEC", "EC layer II - MEC pyramidal island cells")),
  list("EC layer II - all MEC+LEC L2 excitatory", c("Reln", "Calb1"),
       c("EC layer III - MEC", "EC layer III - LEC")),
  list("Nucleus reuniens glutamatergic", "Slc17a6",
       c("CA1 pyramidal cells", "Negative control - striatal D1 MSN")),
  list("Supramammillary glutamatergic", "Slc17a6",
       c("CA1 pyramidal cells", "Negative control - striatal D1 MSN")),
  list("MS/DB cholinergic", "Chat", c("MS/DB GABAergic", "MS/DB glutamatergic")),
  list("Dorsal raphe serotonergic", "Tph2", c("Median raphe serotonergic", "Locus coeruleus noradrenergic")),
  list("Median raphe serotonergic", "Tph2", c("Dorsal raphe serotonergic", "Locus coeruleus noradrenergic")),
  list("Locus coeruleus noradrenergic", "Dbh", c("Dorsal raphe serotonergic", "Median raphe serotonergic")),
  list("Negative control - striatal D1 MSN", c("Drd1", "Drd2"), "Negative control - striatal D2 MSN"),
  list("Negative control - striatal D2 MSN", c("Drd1", "Drd2"), "Negative control - striatal D1 MSN"),
  list("Negative control - cerebellar granule cells", "Gabra6",
       c("Dentate granule cells", "Negative control - olfactory bulb granule cells"))
)
san <- rbindlist(lapply(checks, function(ck) rbindlist(lapply(ck[[2]], function(g) {
  gi <- gidx(g)
  v <- prof[gi, ck[[1]]]
  sv <- sub_mat[gi, ]
  data.table(
    target = ck[[1]], gene = g, target_mean = round(v, 2),
    comparators = paste(sprintf("%s=%.2f", ck[[3]], prof[gi, ck[[3]]]), collapse = "; "),
    subclass_percentile_of_target = round(mean(sv <= v + 5e-4), 3),
    brainwide_max_subclass = sprintf("%s=%.2f", colnames(sub_mat)[which.max(sv)], max(sv)),
    brainwide_median_subclass = round(median(sv), 2))
}))))
fwrite(san, file.path(out_dir, "sanity_check_markers.csv"))
message("source types: ", nrow(map), "; sanity rows: ", nrow(san))

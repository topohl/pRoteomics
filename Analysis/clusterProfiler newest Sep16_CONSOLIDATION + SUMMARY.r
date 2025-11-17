# =======================
# CONSOLIDATION + SUMMARY
# =======================
# Primary: GSEA (gseGO/gseKEGG/Reactome) using NES/FDR
# Secondary: ORA (separate tables/plots)
# Includes class-labeled heatmaps (SR/RC/SC:n) and comparison-level matrices. [web:22]

# 0) Setup --------------------------------------------------------------------
pkgs <- c("dplyr","tidyr","readr","stringr","purrr","ggplot2","igraph","ggraph", "patchwork", "tibble", "clusterProfiler", "org.Mm.eg.db", "DOSE", "enrichplot", "ReactomePA", "ComplexUpset")
missing_pkgs <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(missing_pkgs)) {
  message("Installing missing packages: ", paste(missing_pkgs, collapse=", "))
  tryCatch({
    install.packages(missing_pkgs, repos = c("https://cloud.r-project.org"))
  }, error = function(e) {
    warning("Primary CRAN failed: ", conditionMessage(e), " — retrying ggraph via r-universe")
    if ("ggraph" %in% missing_pkgs) install.packages("ggraph", repos = c("https://thomasp85.r-universe.dev","https://cloud.r-project.org"))
  })
}
suppressPackageStartupMessages({
  for (p in pkgs) {
    tryCatch(library(p, character.only = TRUE),
             error = function(e) warning(sprintf("Failed to load %s: %s", p, conditionMessage(e))))
  }
})

# Root paths (fixed) ----------------------------------------------------------
root_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Results/neuron-phenotypeWithinUnit"
out_dir  <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Results/ZZ_consolidated"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# New: subfolder layout for tables/plots --------------------------------------
dirs <- list(
  tables        = file.path(out_dir, "Tables"),
  tables_gsea   = file.path(out_dir, "Tables", "GSEA"),
  tables_ora    = file.path(out_dir, "Tables", "ORA"),
  tables_class  = file.path(out_dir, "Tables", "Class"),
  tables_cond   = file.path(out_dir, "Tables", "Conditions"),
  tables_qa     = file.path(out_dir, "Tables", "QA"),
  plots         = file.path(out_dir, "Plots"),
  plots_anat    = file.path(out_dir, "Plots", "Anatomy"),
  plots_class   = file.path(out_dir, "Plots", "Class"),
  plots_cond    = file.path(out_dir, "Plots", "Conditions"),
  plots_cm      = file.path(out_dir, "Plots", "ComparisonMatrices"),
  plots_ora     = file.path(out_dir, "Plots", "ORA")
)
invisible(lapply(dirs, function(d) dir.create(d, showWarnings = FALSE, recursive = TRUE)))  # recursive creation [web:8][web:12]

# 1) Helpers ------------------------------------------------------------------
read_if <- function(fp) {
  if (!file.exists(fp)) return(NULL)
  tryCatch(readr::read_csv(fp, guess_max = 10000, progress = FALSE, show_col_types = FALSE),
           error = function(e) { warning("Failed reading: ", fp, " — ", conditionMessage(e)); NULL })
}
read_context <- function(comp_dir) {
  readme <- file.path(comp_dir, "README.txt")
  if (file.exists(readme)) {
    rl <- readLines(readme, warn = FALSE)
    ctx <- rl[grepl("^Context: ", rl)]
    if (length(ctx)) return(sub("^Context:\\s*", "", ctx[1]))
  }
  basename(comp_dir)
}
`%||%` <- function(a,b) if (!is.null(a) && length(a)>0 && !is.na(a)) a else b
tokenize <- function(x){ x <- tolower(x); x <- gsub("[^a-z0-9]+"," ",x); unique(unlist(strsplit(x,"\\s+"))) }
jaccard  <- function(a,b){ ia<-tokenize(a); ib<-tokenize(b); if(!length(ia)||!length(ib)) return(0); length(intersect(ia,ib))/length(union(ia,ib)) }

# 2) GSEA collectors -----------------------------------------------------------
parse_gsea_tbl <- function(fp, db) {
  df <- read_if(fp); if (is.null(df) || !nrow(df)) return(NULL)
  nm <- names(df); has <- function(x) any(x %in% nm)
  if (!has("ID") || !has("Description")) return(NULL)
  if (!has("NES")) return(NULL)
  fdr_col <- c("p.adjust","qvalue","padj","FDR","adj.P.Val")[c("p.adjust","qvalue","padj","FDR","adj.P.Val") %in% nm][1]
  if (is.na(fdr_col)) return(NULL)
  out <- tibble::tibble(
    db = rep(db, nrow(df)),
    ID = df[["ID"]],
    Description = df[["Description"]],
    NES = suppressWarnings(as.numeric(df[["NES"]])),
    FDR = suppressWarnings(as.numeric(df[[fdr_col]]))
  ) %>% dplyr::filter(!is.na(ID), !is.na(NES), !is.na(FDR))
  if (!nrow(out)) return(NULL)
  out
}
empty_gsea <- function(){ tibble::tibble(db=character(), ID=character(), Description=character(), NES=double(), FDR=double()) }

collect_gsea_one <- function(comp_dir) {
  paths <- list(
    GO_BP    = file.path(comp_dir, "10_global","GO_BP","Tables"),
    KEGG     = file.path(comp_dir, "10_global","KEGG","Tables"),
    Reactome = file.path(comp_dir, "10_global","Reactome","Tables")
  )
  if (!any(dir.exists(unlist(paths)))) return(empty_gsea())

  gsego   <- if (dir.exists(paths$GO_BP))    list.files(paths$GO_BP,    pattern="^gseGO_BP_.*\\.csv$",      full.names=TRUE) else character(0)
  gsekegg <- if (dir.exists(paths$KEGG))     list.files(paths$KEGG,     pattern="^gseKEGG_.*\\.csv$",       full.names=TRUE) else character(0)
  rs      <- if (dir.exists(paths$Reactome)) list.files(paths$Reactome, pattern="^ReactomeGSEA_.*\\.csv$",  full.names=TRUE) else character(0)

  if (length(gsego)+length(gsekegg)+length(rs) == 0) return(empty_gsea())

  parts <- list(
    dplyr::bind_rows(lapply(gsego,   parse_gsea_tbl, db="GO_BP")),
    dplyr::bind_rows(lapply(gsekegg, parse_gsea_tbl, db="KEGG")),
    dplyr::bind_rows(lapply(rs,      parse_gsea_tbl, db="Reactome"))
  )
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (!length(parts)) return(empty_gsea())
  dplyr::bind_rows(parts)
}

# 3) ORA collector -------------------------------------------------------------
parse_ora_tbl <- function(fp, db) {
  df <- read_if(fp); if (is.null(df) || !nrow(df)) return(NULL)
  nm <- names(df); if (!("ID" %in% nm && "Description" %in% nm)) return(NULL)
  fdr_col <- c("p.adjust","qvalue","padj","FDR","adj.P.Val")[c("p.adjust","qvalue","padj","FDR","adj.P.Val") %in% nm][1]
  if (is.na(fdr_col)) return(NULL)
  out <- tibble::tibble(db=rep(db, nrow(df)), ID=df[["ID"]], Description=df[["Description"]], Padj=suppressWarnings(as.numeric(df[[fdr_col]]))) %>%
    dplyr::filter(!is.na(ID), !is.na(Padj))
  if (!nrow(out)) return(NULL)
  out
}
empty_ora <- function(){ tibble::tibble(db=character(), ID=character(), Description=character(), Padj=double()) }

collect_ora_one <- function(comp_dir) {
  paths <- list(
    GO_BP    = file.path(comp_dir, "10_global","GO_BP","Tables"),
    KEGG     = file.path(comp_dir, "10_global","KEGG","Tables"),
    Reactome = file.path(comp_dir, "10_global","Reactome","Tables")
  )
  if (!any(dir.exists(unlist(paths)))) return(empty_ora())

  files <- c(
    if (dir.exists(paths$GO_BP))    list.files(paths$GO_BP,    pattern="^ORA_.*\\.csv$",         full.names=TRUE) else character(0),
    if (dir.exists(paths$KEGG))     list.files(paths$KEGG,     pattern="^ORA_.*\\.csv$",         full.names=TRUE) else character(0),
    if (dir.exists(paths$Reactome)) list.files(paths$Reactome, pattern="^ReactomeORA_.*\\.csv$", full.names=TRUE) else character(0)
  )
  if (!length(files)) return(empty_ora())

  binders <- lapply(files, function(fp) {
    db <- if (grepl("/GO_BP/", fp)) "GO_BP" else if (grepl("/KEGG/", fp)) "KEGG" else if (grepl("/Reactome/", fp)) "Reactome" else "Unknown"
    parse_ora_tbl(fp, db = db)
  })
  binders <- binders[!vapply(binders, is.null, logical(1))]
  if (!length(binders)) return(empty_ora())
  dplyr::bind_rows(binders)
}

# 4) Scan comparison folders ---------------------------------------------------
comp_dirs <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
comp_dirs <- comp_dirs[dir.exists(file.path(comp_dirs, "10_global"))]
message(sprintf("Found %d comparison folders under neuron-phenotypeWithinUnit.", length(comp_dirs)))
if (!length(comp_dirs)) stop("No comparison folders with '10_global' found under the fixed root_dir.", call. = FALSE)

# 5) GSEA assembly -------------------------------------------------------------
gsea_list <- lapply(comp_dirs, function(cd) tibble::tibble(comp_dir = cd, context = read_context(cd), data = list(collect_gsea_one(cd))))
gsea_df <- purrr::map_dfr(gsea_list, function(x) { dx <- x$data[[1]]; if (!is.null(dx) && nrow(dx) > 0) cbind(x[1:2], dx) else NULL })
if (!nrow(gsea_df)) stop("No valid GSEA rows found across comparisons.", call. = FALSE)

req <- c("comp_dir","context","db","ID","Description","NES","FDR")
missing_cols <- setdiff(req, names(gsea_df))
if (length(missing_cols)) stop(sprintf("gsea_df missing columns: %s", paste(missing_cols, collapse=", ")), call. = FALSE)

gsea_df <- gsea_df %>% dplyr::mutate(comp_dir=as.character(comp_dir), context=as.character(context), db=as.character(db), ID=as.character(ID), Description=as.character(Description), NES=as.numeric(NES), FDR=as.numeric(FDR))
readr::write_csv(gsea_df, file.path(dirs$tables_gsea, "all_GSEA_long.csv"))  # moved [web:23]
message(sprintf("GSEA rows: %d", nrow(gsea_df)))

# 6) ORA assembly --------------------------------------------------------------
ora_list <- lapply(comp_dirs, function(cd) tibble::tibble(comp_dir = cd, context = read_context(cd), data = list(collect_ora_one(cd))))
ora_df <- purrr::map_dfr(ora_list, function(x) { dx <- x$data[[1]]; if (!is.null(dx) && nrow(dx) > 0) cbind(x[1:2], dx) else NULL })
if (nrow(ora_df)) readr::write_csv(ora_df, file.path(dirs$tables_ora, "all_ORA_long.csv"))

# 7) Theme building (GSEA only) -----------------------------------------------
build_themes_for_db <- function(df_db, top_n_terms = 200, sim_cut = 0.25) {
  if (is.null(df_db) || !nrow(df_db)) return(NULL)
  df_db <- dplyr::mutate(df_db, sig = FDR <= 0.05)

  agg <- df_db %>%
    dplyr::group_by(ID, Description, db) %>%
    dplyr::summarize(n_sig = sum(sig, na.rm = TRUE), n = dplyr::n(), meanNES = mean(NES, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(n_sig), dplyr::desc(abs(meanNES)))
  sel <- head(agg, top_n_terms)
  if (nrow(sel) < 2) return(dplyr::mutate(sel, theme = paste0(unique(df_db$db), "_", dplyr::row_number())))

  ids <- sel$ID; descs <- sel$Description; n <- length(ids)
  sim <- matrix(0, n, n); rownames(sim) <- colnames(sim) <- ids
  for (i in seq_len(n)) for (j in i:n) { s <- if (i==j) 1 else jaccard(descs[i], descs[j]); sim[i,j] <- s; sim[j,i] <- s }

  edge_df <- as.data.frame(as.table(sim), stringsAsFactors = FALSE) |>
    dplyr::rename(ID1 = Var1, ID2 = Var2, sim = Freq) |>
    dplyr::filter(ID1 != ID2, sim >= sim_cut)

  g <- igraph::make_empty_graph(directed = FALSE) |> igraph::add_vertices(n, name = ids)
  if (nrow(edge_df) > 0) {
    edges_mat <- unique(t(apply(edge_df[,c("ID1","ID2")], 1, function(x) sort(as.character(x)))))
    g <- igraph::add_edges(g, as.vector(t(edges_mat)))
  }
  if (igraph::vcount(g) == 0) { sel$theme <- paste0(unique(df_db$db), "_", seq_len(nrow(sel))); return(sel) }

  cl <- igraph::cluster_louvain(g)
  memb <- igraph::membership(cl)
  theme_map <- tibble::tibble(ID = names(memb), theme = paste0(unique(df_db$db), "_C", as.integer(memb)))
  out <- dplyr::left_join(sel, theme_map, by = "ID")
  if (any(is.na(out$theme))) out$theme[is.na(out$theme)] <- paste0(unique(df_db$db), "_", seq_len(sum(is.na(out$theme))))
  out
}

gsea_core <- gsea_df %>% dplyr::select(db, ID, Description, NES, FDR) %>% dplyr::mutate(db=as.character(db), ID=as.character(ID), Description=as.character(Description), NES=as.numeric(NES), FDR=as.numeric(FDR))
req_theme <- c("db","ID","Description","NES","FDR"); miss_theme <- setdiff(req_theme, names(gsea_core)); if (length(miss_theme)) stop(sprintf("gsea_core missing: %s", paste(miss_theme, collapse=", ")), call.=FALSE)
split_dbs <- split(gsea_core, gsea_core$db)
themes_db <- lapply(split_dbs, build_themes_for_db)
themes_db <- themes_db[!vapply(themes_db, is.null, logical(1))]
themes_tbl <- if (length(themes_db)) dplyr::bind_rows(themes_db) else tibble::tibble()
readr::write_csv(themes_tbl, file.path(dirs$tables_gsea, "themes_table.csv"))

# BP membership per theme (cluster), per DB
bp_members <- gsea_df %>%
  dplyr::filter(db == "GO_BP") %>%                                    # only BP here
  dplyr::inner_join(dplyr::select(themes_tbl, ID, theme), by="ID") %>%# map ID -> theme
  dplyr::group_by(db, theme, ID, Description) %>%
  dplyr::summarize(
    n_hits = dplyr::n(),
    meanNES = mean(NES, na.rm=TRUE),
    bestFDR = suppressWarnings(min(FDR, na.rm=TRUE)),
    .groups="drop"
  ) %>%
  dplyr::arrange(db, theme, dplyr::desc(n_hits), bestFDR)

# Save a single long CSV
readr::write_csv(bp_members, file.path(dirs$tables_gsea, "GO_BP_terms_by_theme_long.csv"))

# Optionally export one CSV per theme
bp_split <- split(bp_members, interaction(bp_members$db, bp_members$theme, drop=TRUE))
invisible(purrr::iwalk(bp_split, function(df, nm){
  parts <- strsplit(nm, "\\.")[[1]]
  dbn <- parts[1]; th <- parts[2]
  fp <- file.path(dirs$tables_gsea, paste0("GO_BP_terms_", th, ".csv"))
  readr::write_csv(df, fp)
}))

# 8) Aggregation per comparison x theme (GSEA) --------------------------------
assign_themes <- if (nrow(themes_tbl)) dplyr::inner_join(gsea_df, dplyr::select(themes_tbl, ID, theme), by = "ID") else tibble::tibble()
if (!nrow(assign_themes)) stop("No themes constructed (themes_tbl empty). Consider lowering sim_cut or increasing top_n_terms.", call. = FALSE)

theme_comp <- assign_themes %>%
  dplyr::group_by(db, theme, comp_dir, context) %>%
  dplyr::summarize(n_terms=dplyr::n(), meanNES=mean(NES, na.rm=TRUE), minFDR=suppressWarnings(min(FDR, na.rm=TRUE)), frac_sig=mean(FDR<=0.05, na.rm=TRUE), .groups="drop")
readr::write_csv(theme_comp, file.path(dirs$tables_gsea, "theme_by_comparison.csv"))

# Masking for plots
theme_comp$NES_masked <- ifelse(theme_comp$minFDR <= 0.05, theme_comp$meanNES, NA_real_)

# 9) Robust parsing — both sides (fixed) --------------------------------------
parse_both <- function(ctx) {
  parts <- unlist(strsplit(ctx, "(?:\\s+vs\\s+|_vs_)", perl=TRUE))
  parts <- trimws(parts)
  left  <- if (length(parts) >= 1) parts[1] else NA_character_
  right <- if (length(parts) >= 2) parts[2] else NA_character_

  region_pat <- "(DG|CA1|CA2|CA3)"
  base_layer <- "(mo|po|sr|sp|so|slm|sg)"
  cond_pat   <- "(res|sus|con)"

  parse_side <- function(side) {
    if (is.na(side) || side == "") return(list(region=NA_character_, layer=NA_character_, condition=NA_character_))
    m <- regexec(paste0("^", region_pat, "[_-]?", base_layer, "(", cond_pat, ")?$"), side, perl=TRUE)
    r <- regmatches(side, m)[[1]]
    if (length(r) == 0) {
      m2 <- regexec(paste0("^", region_pat, "[_-]?", base_layer, "(", cond_pat, ")?"), side, perl=TRUE)
      r  <- regmatches(side, m2)[[1]]
    }
    region <- NA_character_; layer <- NA_character_; condition <- NA_character_
    if (length(r) >= 3) {
      region    <- r[2]
      layer_raw <- tolower(r[3]); layer <- layer_raw
      if (length(r) >= 5 && nzchar(r[5])) condition <- tolower(r[5])
    } else {
      toks <- unlist(strsplit(side, "[_-]+"))
      region <- if (length(toks) >= 1 && grepl(paste0("^", region_pat, "$"), toks[1])) toks[1] else NA_character_
      layer  <- if (length(toks) >= 2 && grepl(paste0("^", base_layer, "$"), tolower(toks[2]))) tolower(toks[2]) else NA_character_
      condition <- if (length(toks) >= 3 && grepl(paste0("^", cond_pat, "$"), tolower(toks[3]))) tolower(toks[3]) else NA_character_
    }
    list(region=region, layer=layer, condition=condition)
  }

  L <- parse_side(left)
  R <- parse_side(right)
  list(L=L, R=R)
}

rl_both <- unique(theme_comp$context) %>%
  purrr::map_df(~{
    pr <- parse_both(.x)
    tibble::tibble(
      context   = .x,
      region_L  = pr$L$region, layer_L  = pr$L$layer, condition_L  = pr$L$condition,
      region_R  = pr$R$region, layer_R  = pr$R$layer, condition_R  = pr$R$condition
    )
  })
readr::write_csv(rl_both, file.path(dirs$tables, "context_region_layer_condition_both_sides.csv"))

# 10) Comparison classes (from parsed conditions) -----------------------------
comp_class_map <- rl_both %>%
  dplyr::transmute(
    context,
    cond_left  = condition_L,
    cond_right = condition_R,
    comparison_class = dplyr::case_when(
      cond_left == "sus" & cond_right == "res" ~ "sus vs res",
      cond_left == "res" & cond_right == "con" ~ "res vs con",
      cond_left == "sus" & cond_right == "con" ~ "sus vs con",
      TRUE ~ NA_character_
    )
  )
readr::write_csv(comp_class_map, file.path(dirs$tables_class, "context_comparison_class_map.csv"))

# 11) Localization tables ------------------------------------------------------
theme_loc <- theme_comp %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
  dplyr::filter(!is.na(region), !is.na(layer)) %>%
  dplyr::group_by(db, theme, region, layer) %>%
  dplyr::summarize(meanNES=mean(meanNES, na.rm=TRUE), frac_sig=mean(minFDR<=0.05, na.rm=TRUE), .groups="drop")
readr::write_csv(theme_loc, file.path(dirs$tables_gsea, "localization_table.csv"))

theme_loc_by_class <- theme_comp %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
  dplyr::left_join(comp_class_map %>% dplyr::select(context, comparison_class), by="context") %>%
  dplyr::filter(!is.na(region), !is.na(layer), !is.na(comparison_class)) %>%
  dplyr::group_by(db, theme, comparison_class, region, layer) %>%
  dplyr::summarize(meanNES = mean(meanNES, na.rm=TRUE), n_comp = dplyr::n_distinct(context), .groups="drop")
readr::write_csv(theme_loc_by_class, file.path(dirs$tables_class, "localization_by_comparison_class_table.csv"))

class_labels_df <- theme_comp %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
  dplyr::left_join(comp_class_map, by="context") %>%
  dplyr::filter(!is.na(region), !is.na(layer), !is.na(comparison_class)) %>%
  dplyr::group_by(db, theme, region, layer, comparison_class) %>%
  dplyr::summarize(n_comp = dplyr::n_distinct(context), .groups="drop") %>%
  dplyr::mutate(code = dplyr::case_when(
    comparison_class=="sus vs res" ~ "SR",
    comparison_class=="res vs con" ~ "RC",
    comparison_class=="sus vs con" ~ "SC",
    TRUE ~ ""
  ),
  label = paste0(code, ":", n_comp))
theme_loc_by_class_lab <- theme_loc_by_class %>%
  dplyr::left_join(class_labels_df %>% dplyr::select(db, theme, region, layer, comparison_class, label),
                   by=c("db","theme","region","layer","comparison_class"))
readr::write_csv(theme_loc_by_class_lab, file.path(dirs$tables_class, "localization_by_comparison_class_labeled.csv"))

# Condition-centric tables
theme_loc_left_cond <- theme_comp %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L, condition_L), by="context") %>%
  dplyr::filter(!is.na(region), !is.na(layer), !is.na(condition_L)) %>%
  dplyr::group_by(db, theme, region, layer, condition_L) %>%
  dplyr::summarize(meanNES = mean(meanNES, na.rm=TRUE), n_comp = dplyr::n_distinct(context), .groups="drop") %>%
  dplyr::rename(condition = condition_L)
readr::write_csv(theme_loc_left_cond, file.path(dirs$tables_cond, "localization_left_condition_table.csv"))

theme_loc_right_cond <- theme_comp %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_R, layer=layer_R, condition_R), by="context") %>%
  dplyr::filter(!is.na(region), !is.na(layer), !is.na(condition_R)) %>%
  dplyr::group_by(db, theme, region, layer, condition_R) %>%
  dplyr::summarize(meanNES = mean(meanNES, na.rm=TRUE), n_comp = dplyr::n_distinct(context), .groups="drop") %>%
  dplyr::rename(condition = condition_R)
readr::write_csv(theme_loc_right_cond, file.path(dirs$tables_cond, "localization_right_condition_table.csv"))

# Comparison-level long
comp_level <- theme_comp %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
  dplyr::left_join(comp_class_map %>% dplyr::select(context, comparison_class), by="context") %>%
  dplyr::filter(!is.na(region), !is.na(layer), !is.na(comparison_class)) %>%
  dplyr::mutate(ctx_label = paste0(comparison_class, "::", context))
readr::write_csv(comp_level, file.path(dirs$tables_class, "comparison_level_long.csv"))

# 11b) Directionality and bias summaries (tables)
dom_dir_tbl <- theme_comp %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
  dplyr::filter(!is.na(region), !is.na(layer), minFDR<=0.05) %>%
  dplyr::group_by(db, theme, region, layer) %>%
  dplyr::summarize(medNES = stats::median(meanNES, na.rm=TRUE),
                   n_sig  = dplyr::n(), .groups="drop") %>%
  dplyr::mutate(dir = dplyr::case_when(medNES > 0 ~ "up", medNES < 0 ~ "down", TRUE ~ "flat"))
readr::write_csv(dom_dir_tbl, file.path(dirs$tables_gsea, "dominant_direction_table.csv"))

symmetry_test <- theme_comp %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
  dplyr::filter(!is.na(region), !is.na(layer)) %>%
  dplyr::group_by(db, region, layer) %>%
  dplyr::summarize(p_sym = tryCatch({
      x <- meanNES; x <- x[is.finite(x)]
      if (length(x) >= 6) stats::wilcox.test(x, mu=0, exact=FALSE)$p.value else NA_real_
    }, error=function(e) NA_real_), .groups="drop")
readr::write_csv(symmetry_test, file.path(dirs$tables_qa, "directional_bias_pvalues.csv"))

# 11c) Susceptible drivers across classes (safe column names)
sr_sc_strength <- theme_comp %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region = region_L, layer = layer_L), by = "context") %>%
  dplyr::left_join(comp_class_map %>% dplyr::select(context, comparison_class), by = "context") %>%
  dplyr::filter(!is.na(region), !is.na(layer), !is.na(comparison_class)) %>%
  dplyr::group_by(db, theme, region, layer, comparison_class) %>%
  dplyr::summarize(meanNES_class = mean(meanNES, na.rm = TRUE),
                   n_ctx = dplyr::n_distinct(context), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = comparison_class,
    values_from = c(meanNES_class, n_ctx),
    values_fill = 0
  ) %>%
  # Rename wide columns with spaces to safe snake_case
  dplyr::rename(
    meanNES_class_sus_vs_res = `meanNES_class_sus vs res`,
    meanNES_class_sus_vs_con = `meanNES_class_sus vs con`,
    n_ctx_sus_vs_res         = `n_ctx_sus vs res`,
    n_ctx_sus_vs_con         = `n_ctx_sus vs con`
  ) %>%
  dplyr::mutate(
    SR_strength = meanNES_class_sus_vs_res,
    SC_strength = meanNES_class_sus_vs_con,
    driver_class = dplyr::case_when(
      abs(SR_strength) > abs(SC_strength) ~ "sus vs res",
      abs(SC_strength) > abs(SR_strength) ~ "sus vs con",
      TRUE ~ "tie"
    ),
    driver_gap = abs(SR_strength) - abs(SC_strength)
  )
readr::write_csv(sr_sc_strength, file.path(dirs$tables_class, "susceptibility_driver_strength.csv"))

# 11d) Plots: theme localization
# Driver map: which class (SR vs SC) dominates per theme and region-layer
plot_driver_map <- function(dbn) {
  df <- sr_sc_strength %>%
    dplyr::filter(db == dbn, theme %in% (top_themes %>% dplyr::filter(db==dbn) %>% dplyr::pull(theme)))
  if (!nrow(df)) return(NULL)

  # Compact label: nSR/nSC contributing contexts (uses safe column names)
  df <- df %>%
    dplyr::mutate(n_lab = paste0(n_ctx_sus_vs_res, "/", n_ctx_sus_vs_con))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = layer, y = region, fill = driver_class, alpha = abs(driver_gap))) +
    ggplot2::geom_tile(color="white") +
    ggplot2::geom_text(ggplot2::aes(label = n_lab), size=3) +
    ggplot2::scale_fill_manual(values = c("sus vs res" = "#377eb8", "sus vs con" = "#e41a1c", "tie" = "grey80")) +
    ggplot2::scale_alpha(range = c(0.3, 1), guide = "none") +
    ggplot2::facet_wrap(~ theme, ncol = 4) +
    ggplot2::labs(title = paste0(dbn, " — drivers of susceptibility (SR vs SC)"),
                  x = "Layer", y = "Region", fill = "Class") +
    ggplot2::theme_minimal(10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))

  fp <- file.path(dirs$plots_class, paste0("drivers_susceptibility_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=22, height=12); message("[drivers-map] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(sr_sc_strength$db), plot_driver_map))

# 11e) Plots: theme driver volcano
theme_driver_volcano <- sr_sc_strength %>%
  dplyr::group_by(db, theme) %>%
  dplyr::summarize(delta = mean(SR_strength - SC_strength, na.rm=TRUE),
                   abs_gap = mean(abs(SR_strength) - abs(SC_strength), na.rm=TRUE),
                   dom_class = dplyr::case_when(
                     abs_gap > 0 ~ "sus vs res",
                     abs_gap < 0 ~ "sus vs con",
                     TRUE ~ "tie"
                   ),
                   .groups="drop") %>%
  dplyr::left_join(theme_repl %>% dplyr::select(db, theme, repl_score, dom_dir), by=c("db","theme"))

plot_theme_driver_volcano <- function(dbn) {
  df <- theme_driver_volcano %>% dplyr::filter(db==dbn)
  if (!nrow(df)) return(NULL)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = delta, y = repl_score, color = dom_class, shape = dom_dir)) +
    ggplot2::geom_hline(yintercept = median(df$repl_score, na.rm=TRUE), linetype="dashed", color="grey60") +
    ggplot2::geom_vline(xintercept = 0, linetype="dotted", color="grey60") +
    ggplot2::geom_point(size=3, alpha=0.9) +
    ggplot2::scale_color_manual(values=c("sus vs res"="#377eb8","sus vs con"="#e41a1c","tie"="grey60")) +
    ggplot2::labs(title=paste0(dbn, " — theme-level SR vs SC dominance"),
                  x="SR_strength − SC_strength (mean NES)", y="Replication score",
                  color="Dominant class", shape="Dir") +
    ggplot2::theme_minimal(11)
  fp <- file.path(dirs$plots_class, paste0("theme_driver_volcano_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=10, height=7); message("[driver-volcano] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(theme_driver_volcano$db), plot_theme_driver_volcano))

# 12) NA diagnostics -----------------------------------------------------------
cov_tbl <- assign_themes %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by = "context") %>%
  dplyr::filter(!is.na(region), !is.na(layer)) %>%
  dplyr::group_by(db, theme, region, layer) %>%
  dplyr::summarize(n_comp_cover=dplyr::n_distinct(comp_dir), n_sig_cover=sum(FDR<=0.05, na.rm=TRUE), minFDR_all=suppressWarnings(min(FDR, na.rm=TRUE)), .groups="drop")
na_diag <- theme_loc %>%
  dplyr::full_join(cov_tbl, by = c("db","theme","region","layer")) %>%
  dplyr::mutate(cause = dplyr::case_when(
    is.na(meanNES) & (is.na(n_comp_cover) | n_comp_cover == 0) ~ "no_coverage",
    is.na(meanNES) & (n_comp_cover > 0) & (is.na(n_sig_cover) | n_sig_cover == 0) ~ "masked_nonsignificant",
    is.na(meanNES) ~ "other_missing",
    TRUE ~ "has_value"
  ))
readr::write_csv(na_diag, file.path(dirs$tables_qa, "localization_NA_diagnostics.csv"))

# 13) Replication summary ------------------------------------------------------
theme_repl <- theme_comp %>% dplyr::mutate(sign=sign(meanNES), sig=minFDR<=0.05) %>%
  dplyr::group_by(db, theme) %>% dplyr::summarize(n_comp=dplyr::n(), n_sig=sum(sig, na.rm=TRUE), n_pos=sum(sig & sign>0, na.rm=TRUE), n_neg=sum(sig & sign<0, na.rm=TRUE), dom_dir=ifelse(n_pos>=n_neg,"up","down"), repl_score=pmax(n_pos,n_neg), .groups="drop") %>%
  dplyr::arrange(dplyr::desc(repl_score))
readr::write_csv(theme_repl, file.path(dirs$tables_gsea, "theme_replication_scores.csv"))

top_themes <- theme_repl %>% dplyr::group_by(db) %>% dplyr::slice_max(order_by=repl_score, n=12) %>% dplyr::ungroup()

# 14) Plots -------------------------------------------------------------------
plot_theme_localization <- function(dbn) {
  df <- theme_loc %>% dplyr::filter(db==dbn, theme %in% (top_themes %>% dplyr::filter(db==dbn) %>% dplyr::pull(theme)))
  if (!nrow(df)) { message("[anatomy] No data for ", dbn); return(NULL) }
  p <- ggplot2::ggplot(df, ggplot2::aes(x=layer, y=region, fill=meanNES)) +
    ggplot2::geom_tile(color="grey90") +
    ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b") +
    ggplot2::facet_wrap(~ theme, scales="free", ncol=4) +
    ggplot2::labs(title=paste0("Localization by region-layer — ", dbn), x="Layer", y="Region") +
    ggplot2::theme_minimal(base_size=10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))
  fp <- file.path(dirs$plots_anat, paste0("localization_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=22, height=14); message("[anatomy] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(theme_loc$db), plot_theme_localization))

plot_theme_localization_by_class_labeled <- function(dbn) {
  df_all <- theme_loc_by_class_lab %>% dplyr::filter(db==dbn)
  tops <- top_themes %>% dplyr::filter(db==dbn) %>% dplyr::pull(theme)
  if (!length(tops) || !any(df_all$theme %in% tops)) {
    tops <- df_all %>% dplyr::count(theme, sort=TRUE) %>% dplyr::slice_head(n=12) %>% dplyr::pull(theme)
  }
  df <- df_all %>% dplyr::filter(theme %in% tops)
  if (!nrow(df)) { message("[class-labeled] No data for ", dbn, " after fallback."); return(NULL) }
  class_levels <- c("sus vs res","res vs con","sus vs con")
  df$comparison_class <- factor(df$comparison_class, levels = class_levels)
  p <- ggplot2::ggplot(df, ggplot2::aes(x=layer, y=region, fill=meanNES)) +
    ggplot2::geom_tile(color="grey90") +
    ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b") +
    ggplot2::facet_grid(rows = ggplot2::vars(theme), cols = ggplot2::vars(comparison_class), scales="free", drop=FALSE) +
    ggplot2::labs(title=paste0("Localization — comparison class (labeled) — ", dbn), x="Layer", y="Region") +
    ggplot2::theme_minimal(base_size=10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1)) +
    ggplot2::geom_text(ggplot2::aes(label = dplyr::coalesce(label, as.character(n_comp))), color="black", size=3, na.rm=TRUE)
  fp <- file.path(dirs$plots_class, paste0("localization_", dbn, "_by_comparison_class_labeled.svg"))
  ggplot2::ggsave(fp, p, width=30, height=20); message("[class-labeled] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(theme_loc_by_class_lab$db), plot_theme_localization_by_class_labeled))

plot_theme_localization_left_cond <- function(dbn) {
  df <- theme_loc_left_cond %>% dplyr::filter(db==dbn, theme %in% (top_themes %>% dplyr::filter(db==dbn) %>% dplyr::pull(theme)))
  if (!nrow(df)) { message("[left-cond] No data for ", dbn); return(NULL) }
  p <- ggplot2::ggplot(df, ggplot2::aes(x=layer, y=region, fill=meanNES)) +
    ggplot2::geom_tile(color="grey90") +
    ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b") +
    ggplot2::facet_grid(rows = ggplot2::vars(theme), cols = ggplot2::vars(condition), scales="free") +
    ggplot2::labs(title=paste0("Localization — tested condition (sus/res) — ", dbn), x="Layer", y="Region") +
    ggplot2::theme_minimal(base_size=10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))
  fp <- file.path(dirs$plots_cond, paste0("localization_", dbn, "_tested_condition.svg"))
  ggplot2::ggsave(fp, p, width=26, height=20); message("[left-cond] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(theme_loc_left_cond$db), plot_theme_localization_left_cond))

plot_theme_localization_right_cond <- function(dbn) {
  df <- theme_loc_right_cond %>% dplyr::filter(db==dbn, theme %in% (top_themes %>% dplyr::filter(db==dbn) %>% dplyr::pull(theme)))
  if (!nrow(df)) { message("[right-cond] No data for ", dbn); return(NULL) }
  p <- ggplot2::ggplot(df, ggplot2::aes(x=layer, y=region, fill=meanNES)) +
    ggplot2::geom_tile(color="grey90") +
    ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b") +
    ggplot2::facet_grid(rows = ggplot2::vars(theme), cols = ggplot2::vars(condition), scales="free") +
    ggplot2::labs(title=paste0("Localization — baseline condition (res/con) — ", dbn), x="Layer", y="Region") +
    ggplot2::theme_minimal(base_size=10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))
  fp <- file.path(dirs$plots_cond, paste0("localization_", dbn, "_baseline_condition.svg"))
  ggplot2::ggsave(fp, p, width=26, height=20); message("[right-cond] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(theme_loc_right_cond$db), plot_theme_localization_right_cond))

plot_theme_by_comparisons <- function(dbn, max_cols=60) {
  df_db <- comp_level %>% dplyr::filter(db==dbn)
  themes_here <- df_db %>% dplyr::count(theme, sort=TRUE) %>% dplyr::slice_head(n=12) %>% dplyr::pull(theme)
  if (!length(themes_here)) { message("[comp-matrix] No themes with class data for ", dbn); return(NULL) }
  lapply(themes_here, function(th) {
    df <- df_db %>% dplyr::filter(theme==th)
    ctx_ord <- df %>% dplyr::count(ctx_label, sort=TRUE) %>% dplyr::slice_head(n=max_cols) %>% dplyr::pull(ctx_label)
    df <- df %>% dplyr::filter(ctx_label %in% ctx_ord)
    df$ctx_label <- factor(df$ctx_label, levels=ctx_ord)
    p <- ggplot2::ggplot(df, ggplot2::aes(x=ctx_label, y=paste(region, layer, sep="_"), fill=meanNES)) +
      ggplot2::geom_tile(color="grey90") +
      ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b") +
      ggplot2::labs(title=paste0("Theme ", th, " — comparison-level NES — ", dbn), x="Comparison (class::context)", y="Region_Layer") +
      ggplot2::theme_minimal(base_size=10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle=60, hjust=1, vjust=1))
    fp <- file.path(dirs$plots_cm, paste0("theme_", gsub("[^A-Za-z0-9]+","_", th), "_", dbn, "_comparison_matrix.svg"))
    ggplot2::ggsave(fp, p, width=28, height=10); message("[comp-matrix] Saved: ", fp)
    invisible(p)
  })
}
invisible(lapply(unique(comp_level$db), plot_theme_by_comparisons))

# 14b) Compact panel (Anatomy + Class + Replication) per DB
plot_compact_panel <- function(dbn, n_themes = 6) {
  the <- top_themes %>% dplyr::filter(db==dbn) %>% dplyr::slice_head(n=n_themes) %>% dplyr::pull(theme)

  p1 <- ggplot2::ggplot(theme_loc %>% dplyr::filter(db==dbn, theme %in% the),
                        ggplot2::aes(x=layer,y=region,fill=meanNES)) +
    ggplot2::geom_tile(color="grey90") +
    ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b") +
    ggplot2::facet_wrap(~ theme, ncol=3) +
    ggplot2::labs(title=paste0(dbn, " — localization (NES)")) + ggplot2::theme_minimal(9)

  p2 <- ggplot2::ggplot(theme_loc_by_class_lab %>% dplyr::filter(db==dbn, theme %in% the),
                        ggplot2::aes(x=layer,y=region,fill=meanNES,label=label)) +
    ggplot2::geom_tile(color="grey95") +
    ggplot2::geom_text(size=2, color="black", na.rm=TRUE) +
    ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b") +
    ggplot2::facet_grid(rows=ggplot2::vars(theme), cols=ggplot2::vars(comparison_class), drop=FALSE) +
    ggplot2::labs(title="comparison class (labels=SR/RC/SC:n)") + ggplot2::theme_minimal(9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45,hjust=1))

  p3 <- ggplot2::ggplot(theme_repl %>% dplyr::filter(db==dbn, theme %in% the),
                        ggplot2::aes(x=reorder(theme, repl_score), y=repl_score, fill=dom_dir)) +
    ggplot2::geom_col() + ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values=c(down="#2166ac", up="#b2182b")) +
    ggplot2::labs(title="replication score", x="", y="max(up, down) significant") +
    ggplot2::theme_minimal(9)

  if (requireNamespace("patchwork", quietly=TRUE)) {
    panel <- p1 / p2 / p3 + patchwork::plot_layout(heights=c(1,1.1,0.6))
    fp <- file.path(dirs$plots, "Panels", paste0("panel_", dbn, "_compact.svg"))
    dir.create(dirname(fp), recursive=TRUE, showWarnings=FALSE)
    ggplot2::ggsave(fp, panel, width=14, height=16); message("[panel] Saved: ", fp)
  } else {
    message("[panel] patchwork not installed; skipping combined panel for ", dbn)
  }
}
invisible(lapply(unique(theme_loc$db), plot_compact_panel))

# 14c) Dominant direction lattice per DB (median NES, n_sig annotation)
plot_dir_lattice <- function(dbn) {
  df <- dom_dir_tbl %>% dplyr::filter(db==dbn, theme %in% (top_themes %>% dplyr::filter(db==dbn) %>% dplyr::pull(theme)))
  if (!nrow(df)) return(NULL)
  p <- ggplot2::ggplot(df, ggplot2::aes(x=layer, y=region, fill=dir, label=n_sig)) +
    ggplot2::geom_tile(color="white") + ggplot2::geom_text(size=3, color="black") +
    ggplot2::scale_fill_manual(values=c(down="#2166ac", flat="grey85", up="#b2182b")) +
    ggplot2::facet_wrap(~ theme, ncol=4) +
    ggplot2::labs(title=paste0(dbn, " — dominant direction (median NES, n_sig)"),
                  x="Layer", y="Region") + ggplot2::theme_minimal(9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))
  fp <- file.path(dirs$plots_anat, paste0("dominant_direction_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=22, height=12); message("[dir-lattice] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(theme_loc$db), plot_dir_lattice))

# 14d) Directional bias tiles per DB (−log10 p from Wilcoxon against 0)
plot_symmetry <- function(dbn) {
  df <- symmetry_test %>% dplyr::filter(db==dbn)
  p <- ggplot2::ggplot(df, ggplot2::aes(x=layer, y=region, fill=-log10(p_sym))) +
    ggplot2::geom_tile(color="white") +
    ggplot2::scale_fill_viridis_c(option="C", na.value="grey90") +
    ggplot2::labs(title=paste0(dbn, " — directional bias (−log10 p)"), x="Layer", y="Region") +
    ggplot2::theme_minimal(9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))
  fp <- file.path(dirs$plots_anat, paste0("directional_bias_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=8, height=6); message("[bias] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(theme_loc$db), plot_symmetry))

# 14e) Comparison-class balance per theme (stacked fractions across region-layer)
class_balance <- theme_comp %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
  dplyr::left_join(comp_class_map, by="context") %>%
  dplyr::filter(!is.na(region), !is.na(layer), !is.na(comparison_class)) %>%
  dplyr::group_by(db, theme, region, layer, comparison_class) %>%
  dplyr::summarize(n_cmp=dplyr::n_distinct(context), .groups="drop") %>%
  dplyr::group_by(db, theme, region, layer) %>%
  dplyr::mutate(frac = n_cmp/sum(n_cmp)) %>% dplyr::ungroup()
readr::write_csv(class_balance, file.path(dirs$tables_class, "comparison_class_balance.csv"))

plot_class_balance <- function(dbn, theme_id) {
  df <- class_balance %>% dplyr::filter(db==dbn, theme==theme_id)
  if (!nrow(df)) return(NULL)
  p <- ggplot2::ggplot(df, ggplot2::aes(x=paste(region,layer,sep="_"), y=frac, fill=comparison_class)) +
    ggplot2::geom_col(width=0.8) +
    ggplot2::scale_fill_manual(values=c("sus vs res"="#377eb8", "res vs con"="#4daf4a", "sus vs con"="#e41a1c")) +
    ggplot2::labs(title=paste0(theme_id, " — comparison balance (", dbn, ")"), x="Region_Layer", y="Fraction of contexts") +
    ggplot2::theme_minimal(9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=60, hjust=1))
  fp <- file.path(dirs$plots_class, paste0("class_balance_", gsub("[^A-Za-z0-9]+","_", theme_id), "_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=12, height=6); message("[class-balance] Saved: ", fp)
  invisible(p)
}
invisible(lapply(split(top_themes, top_themes$db), function(dfdb) {
  dbn <- unique(dfdb$db); ths <- head(dfdb$theme, 8)
  lapply(ths, function(th) plot_class_balance(dbn, th))
}))

# 14f) Contrast strip (compact comparison matrix subset) for top themes
plot_theme_strip <- function(dbn, theme_id, max_ctx=40) {
  df <- comp_level %>% dplyr::filter(db==dbn, theme==theme_id)
  if (!nrow(df)) return(NULL)
  ctx_ord <- df %>% dplyr::count(ctx_label, sort=TRUE) %>% dplyr::slice_head(n=max_ctx) %>% dplyr::pull(ctx_label)
  df <- df %>% dplyr::filter(ctx_label %in% ctx_ord)
  df$ctx_label <- factor(df$ctx_label, levels=ctx_ord)

  p <- ggplot2::ggplot(df, ggplot2::aes(x=ctx_label, y=paste(region,layer,sep="_"), fill=meanNES)) +
    ggplot2::geom_tile(color="grey95") +
    ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b") +
    ggplot2::facet_wrap(~ comparison_class, ncol=1, scales="free_x") +
    ggplot2::labs(title=paste0(theme_id, " — contrast strip (", dbn, ")"), x="Comparison", y="Region_Layer") +
    ggplot2::theme_minimal(9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=60, hjust=1))
  fp <- file.path(dirs$plots_cm, paste0("theme_", gsub("[^A-Za-z0-9]+","_", theme_id), "_", dbn, "_strip.svg"))
  ggplot2::ggsave(fp, p, width=14, height=8); message("[strip] Saved: ", fp)
  invisible(p)
}
invisible(lapply(split(top_themes, top_themes$db), function(dfdb) {
  dbn <- unique(dfdb$db); ths <- head(dfdb$theme, 6)
  lapply(ths, function(th) plot_theme_strip(dbn, th, max_ctx=40))
}))

# 14g) Simple localization graph per DB (themes ↔ region-layer nodes)
t_thr <- 0.6
edges_df <- theme_loc %>%
  dplyr::left_join(
    theme_comp %>%
      dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L),
                       by="context") %>%
      dplyr::filter(!is.na(region), !is.na(layer), minFDR<=0.05) %>%
      dplyr::group_by(db, theme, region, layer) %>%
      dplyr::summarize(n_sig = dplyr::n(), sign_dir = sign(mean(meanNES, na.rm=TRUE)), .groups="drop"),
    by=c("db","theme","region","layer")
  ) %>%
  # guard: keep rows with numeric NES and threshold on absolute value
  dplyr::filter(!is.na(meanNES), is.finite(meanNES), abs(meanNES) >= t_thr) %>%
  # vectorized fallback for missing n_sig
  dplyr::mutate(n_sig = ifelse(is.na(n_sig), 1L, n_sig),
                src = theme,
                dst = paste(region, layer, sep="_"),
                edge_col = ifelse(meanNES > 0, "#b2182b", "#2166ac"),
                edge_w = scales::rescale(pmax(n_sig, 1L), to=c(0.4, 2)))

plot_localization_graph <- function(dbn) {
  df <- edges_df %>% dplyr::filter(db==dbn)
  if (!nrow(df)) return(NULL)

  # Build graph
  g <- igraph::graph_from_data_frame(df %>% dplyr::select(src, dst), directed=FALSE)

  # Align edge attributes to graph edges
  e_df <- igraph::as_data_frame(g, what="edges")
  key_g  <- paste(e_df$from, e_df$to, sep="||")
  key_df <- paste(df$src, df$dst, sep="||")
  idx <- match(key_g, key_df)
  mis <- which(is.na(idx))
  if (length(mis)) {
    key_rev <- paste(e_df$to[mis], e_df$from[mis], sep="||")
    idx_rev <- match(key_rev, key_df)
    idx[mis] <- idx_rev
  }
  idx[is.na(idx)] <- 1L

  igraph::E(g)$edge_col <- df$edge_col[idx]
  igraph::E(g)$edge_w   <- df$edge_w[idx]

  set.seed(1)
  p <- ggraph::ggraph(g, layout="fr") +
    ggraph::geom_edge_link(ggplot2::aes(edge_colour = edge_col,
                                        edge_width  = edge_w), alpha=0.8) +
    ggraph::geom_node_point(ggplot2::aes(shape = ifelse(grepl("_", name), "RL", "Theme")), size=3) +
    ggraph::geom_node_text(ggplot2::aes(label=name), repel=TRUE, size=3) +
    ggraph::scale_edge_colour_identity() +
    ggraph::scale_edge_width(range=c(0.4,2)) +
    ggplot2::scale_shape_manual(values=c(Theme=16, RL=15)) +
    ggplot2::theme_void() +
    ggplot2::labs(title=paste0(dbn, " — localization graph (|NES|≥", t_thr, ")"))

  fp <- file.path(dirs$plots, "Graphs", paste0("localization_graph_", dbn, ".svg"))
  dir.create(dirname(fp), recursive=TRUE, showWarnings=FALSE)
  ggplot2::ggsave(fp, p, width=10, height=8); message("[graph] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(theme_loc$db), plot_localization_graph))


# 14h) Terms network per DB (Jaccard similarity of GSEA terms)
# --- Terms network (per DB) ---
# Recompute similarity within DB to match your themes more closely
# Term network colored by theme cluster
tokenize <- function(x){
  x <- tolower(x); x <- gsub("[^a-z0-9 ]+"," ",x)
  unique(unlist(strsplit(x,"\\s+")))
}
jaccard <- function(a,b){
  ia <- tokenize(a); ib <- tokenize(b)
  if (!length(ia) || !length(ib)) return(0)
  length(intersect(ia,ib)) / length(union(ia,ib))
}

rep_terms <- gsea_df %>%
  dplyr::inner_join(dplyr::select(themes_tbl, ID, theme), by="ID") %>%
  dplyr::group_by(db, theme, ID, Description) %>%
  dplyr::summarize(freq = dplyr::n(), .groups="drop") %>%
  dplyr::group_by(db, theme) %>%
  dplyr::slice_max(order_by=freq, n=10, with_ties=FALSE) %>%
  dplyr::ungroup()

plot_terms_network <- function(dbn, sim_cut=0.25, max_nodes=250) {
  df <- rep_terms %>% dplyr::filter(db==dbn)
  if (!nrow(df)) return(NULL)
  # control size: top 12 themes
  top_themes <- df %>% dplyr::count(theme, sort=TRUE) %>% dplyr::slice_head(n=12) %>% dplyr::pull(theme)
  df <- df %>% dplyr::filter(theme %in% top_themes)
  if (nrow(df) > max_nodes) {
    per <- max(5L, floor(max_nodes/length(top_themes)))
    df <- df %>% dplyr::group_by(theme) %>% dplyr::slice_head(n=per) %>% dplyr::ungroup()
  }

  terms <- df$Description; n <- length(terms)
  if (n < 2) return(NULL)
  edges <- list()
  for (i in seq_len(n-1)) for (j in (i+1):n) {
    s <- jaccard(terms[i], terms[j]); if (s >= sim_cut) edges[[length(edges)+1]] <- c(i,j,s)
  }
  if (!length(edges)) return(NULL)
  edges <- as.data.frame(do.call(rbind, edges)); names(edges) <- c("i","j","sim")
  edges$i <- as.integer(edges$i); edges$j <- as.integer(edges$j); edges$sim <- as.numeric(edges$sim)

  verts <- tibble::tibble(name = terms, theme = df$theme, term = df$Description)
  g <- igraph::graph_from_data_frame(
    d = tibble::tibble(from = verts$name[edges$i], to = verts$name[edges$j], sim = edges$sim),
    directed = FALSE, vertices = verts
  )

  set.seed(1)
  p <- ggraph::ggraph(g, layout="fr") +
    ggraph::geom_edge_link(ggplot2::aes(width=sim), colour="grey70", alpha=0.35) +
    ggraph::geom_node_point(ggplot2::aes(color=theme), size=2) +
    ggraph::geom_node_text(ggplot2::aes(label=term), size=2.4, repel=TRUE) +
    ggplot2::scale_edge_width(range=c(0.2,1.4), guide="none") +
    ggplot2::theme_void() +
    ggplot2::labs(title=paste0(dbn, " — terms network (sim≥", sim_cut, ")"), color="Theme")
  fp <- file.path(dirs$plots_anat, paste0("terms_network_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=16, height=12); message("[terms-net] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(rep_terms$db), plot_terms_network))

suppressWarnings({
  if (!requireNamespace("uwot", quietly=TRUE)) message("Optional: install.packages('uwot') for UMAP")
})

build_theme_kw_matrix <- function(dbn, top_kw=800) {
  df <- gsea_df %>%
    dplyr::filter(db==dbn) %>%
    dplyr::inner_join(dplyr::select(themes_tbl, ID, theme), by="ID") %>%
    dplyr::distinct(theme, Description)
  if (!nrow(df)) return(NULL)

  toks <- strsplit(tolower(gsub("[^a-z0-9 ]+"," ", df$Description)), "\\s+")
  toks <- lapply(toks, function(v) v[nzchar(v) & !v %in% c("and","or","of","the","to","in","by","for","process","regulation","cellular","activity","protein","pathway","signaling")])

  env <- new.env(parent=emptyenv())
  for (i in seq_len(nrow(df))) {
    th <- df$theme[i]
    for (k in toks[[i]]) {
      key <- paste(th,k,sep="||"); env[[key]] <- (env[[key]] %||% 0L) + 1L
    }
  }
  keys <- ls(env); if (!length(keys)) return(NULL)
  sp <- strsplit(keys,"\\|\\|"); ths <- vapply(sp, `[`, character(1), 1); kws <- vapply(sp, `[`, character(1), 2)
  val <- as.integer(mget(keys, env, ifnotfound=0L))
  mat <- tibble::tibble(theme=ths, kw=kws, val=val) %>%
    dplyr::group_by(kw) %>% dplyr::summarize(dfreq=dplyr::n(), .groups="drop") %>%
    dplyr::right_join(tibble::tibble(theme=ths, kw=kws, val=val), by="kw") %>%
    dplyr::filter(dfreq>=2) %>%
    dplyr::group_by(kw) %>% dplyr::mutate(tf = val/sum(val)) %>% dplyr::ungroup() %>%
    dplyr::group_by(kw) %>% dplyr::mutate(idf = log(n_distinct(theme)/n_distinct(theme[val>0])+1e-6)) %>% dplyr::ungroup() %>%
    dplyr::mutate(tfidf = tf*idf)

  vocab <- mat %>% dplyr::group_by(kw) %>% dplyr::summarize(s=sum(tfidf,na.rm=TRUE), .groups="drop") %>%
    dplyr::arrange(dplyr::desc(s)) %>% dplyr::slice_head(n=top_kw) %>% dplyr::pull(kw)
  mat <- mat %>% dplyr::filter(kw %in% vocab)
  if (!nrow(mat)) return(NULL)

  M <- tidyr::pivot_wider(mat, names_from=kw, values_from=tfidf, values_fill=0)
  M
}

plot_theme_embedding <- function(dbn) {
  M <- build_theme_kw_matrix(dbn)
  if (is.null(M) || nrow(M) < 2) return(NULL)
  meta <- M %>% dplyr::select(theme)
  X <- as.matrix(M %>% dplyr::select(-theme))
  rownames(X) <- meta$theme

  if (requireNamespace("uwot", quietly=TRUE)) {
    set.seed(1); emb <- uwot::umap(X, n_neighbors=10, min_dist=0.2, metric="cosine")
    emb <- as.data.frame(emb); colnames(emb) <- c("U1","U2")
  } else {
    pr <- stats::prcomp(X, scale.=TRUE); emb <- as.data.frame(pr$x[,1:2]); colnames(emb) <- c("U1","U2")
  }
  emb$theme <- rownames(X)

  p <- ggplot2::ggplot(emb, ggplot2::aes(x=U1, y=U2, label=theme)) +
    ggplot2::geom_point(color="#2c7fb8", size=2, alpha=0.85) +
    ggplot2::geom_text(size=2.8, nudge_y=0.02) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title=paste0(dbn, " — theme embedding (", if (requireNamespace("uwot", quietly=TRUE)) "UMAP" else "PCA", ")"),
                  x="Dim 1", y="Dim 2")
  fp <- file.path(dirs$plots_anat, paste0("theme_embedding_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=10, height=8); message("[theme-embed] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(themes_tbl$db), plot_theme_embedding))

plot_theme_dendrogram <- function(dbn) {
  M <- build_theme_kw_matrix(dbn)
  if (is.null(M) || nrow(M) < 2) return(NULL)
  X <- as.matrix(M %>% dplyr::select(-theme))
  rownames(X) <- M$theme
  # cosine distance
  A <- X + 1e-12
  nrm <- sqrt(rowSums(A*A))
  S <- tcrossprod(A / nrm)
  D <- 1 - S
  hc <- stats::hclust(stats::as.dist(D), method="average")

  if (!requireNamespace("ggdendro", quietly=TRUE)) {
    message("Optional: install.packages('ggdendro') for dendrogram plot"); return(NULL)
  }
  dd <- ggdendro::dendro_data(hc)
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data=dd$segments, ggplot2::aes(x=x, y=y, xend=xend, yend=yend)) +
    ggplot2::geom_text(data=dd$labels, ggplot2::aes(x=x, y=y, label=label), angle=90, hjust=1, vjust=0.5, size=2.6) +
    ggplot2::theme_void() + ggplot2::labs(title=paste0(dbn, " — theme dendrogram"))
  fp <- file.path(dirs$plots_anat, paste0("theme_dendrogram_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=14, height=10); message("[theme-dend] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(themes_tbl$db), plot_theme_dendrogram))

# 14i) GO BP terms per theme (bar, top by FDR/frequency)
plot_bp_terms_per_theme <- function(top_k = 8) {
  df <- gsea_df %>%
    dplyr::filter(db == "GO_BP") %>%
    dplyr::inner_join(dplyr::select(themes_tbl, ID, theme), by="ID") %>%
    dplyr::group_by(theme, Description) %>%
    dplyr::summarize(n_hits = dplyr::n(), meanNES = mean(NES, na.rm=TRUE),
                     bestFDR = suppressWarnings(min(FDR, na.rm=TRUE)), .groups="drop") %>%
    dplyr::group_by(theme) %>%
    dplyr::arrange(bestFDR, dplyr::desc(n_hits), .by_group=TRUE) %>%
    dplyr::slice_head(n = top_k) %>%
    dplyr::ungroup()

  if (!nrow(df)) { message("[bp-terms] No GO_BP data"); return(NULL) }

  # Reorder within each facet
  df <- df %>%
    dplyr::group_by(theme) %>%
    dplyr::mutate(term_order = reorder(Description, -n_hits)) %>%
    dplyr::ungroup()

  p <- ggplot2::ggplot(df, ggplot2::aes(x = n_hits, y = term_order, fill = -log10(bestFDR))) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_viridis_c(option="C", name = "-log10(FDR)") +
    ggplot2::facet_wrap(~ theme, scales = "free_y", ncol = 3) +
    ggplot2::labs(title = "GO BP terms per cluster (top by FDR/frequency)",
                  x = "Term count within cluster", y = "GO BP term") +
    ggplot2::theme_minimal(10)
  fp <- file.path(dirs$plots_anat, "GO_BP_terms_per_theme.svg")
  ggplot2::ggsave(fp, p, width=18, height=14); message("[bp-terms] Saved: ", fp)
  invisible(p)
}
plot_bp_terms_per_theme(top_k = 8)

bp_theme_summaries <- bp_members %>%
  dplyr::group_by(theme) %>%
  dplyr::arrange(bestFDR, dplyr::desc(n_hits), .by_group=TRUE) %>%
  dplyr::summarize(
    n_terms = dplyr::n(),
    top_terms = paste(head(Description, 6), collapse=" | "),
    medianNES = stats::median(meanNES, na.rm=TRUE),
    .groups="drop"
  )
readr::write_csv(bp_theme_summaries, file.path(dirs$tables_gsea, "GO_BP_theme_summaries.csv"))

# 15) ORA summaries (separate) ------------------------------------------------
if (nrow(ora_df)) {
  ora_repl <- ora_df %>% dplyr::group_by(db, ID, Description) %>% dplyr::summarize(n_comp=dplyr::n(), n_sig=sum(Padj<=0.05, na.rm=TRUE), meanPadj=mean(Padj, na.rm=TRUE), .groups="drop") %>% dplyr::arrange(dplyr::desc(n_sig), meanPadj)
  readr::write_csv(ora_repl, file.path(dirs$tables_ora, "ORA_replication_scores.csv"))
  top_ora <- ora_repl %>% dplyr::group_by(db) %>% dplyr::slice_max(n_sig, n=25) %>% dplyr::ungroup()
  p_ora <- ggplot2::ggplot(top_ora, ggplot2::aes(x=reorder(paste(db, Description, sep="::"), n_sig), y=n_sig, fill=db)) +
    ggplot2::geom_col() + ggplot2::coord_flip() +
    ggplot2::labs(title="ORA term replication (count of significant comparisons)", x="DB::Term", y="# significant comps") +
    ggplot2::theme_minimal(base_size=10)
  ggplot2::ggsave(file.path(dirs$plots_ora, "ORA_replication_bar.svg"), p_ora, width=18, height=12)
}

# ORA summary for sus vs res
if (nrow(ora_df)) {
  ora_susres <- ora_df %>%
    dplyr::left_join(comp_class_map, by="context") %>%
    dplyr::filter(comparison_class == "sus vs res") %>%
    dplyr::group_by(db, ID, Description) %>%
    dplyr::summarize(n_comp=dplyr::n(), n_sig=sum(Padj<=0.05, na.rm=TRUE), meanPadj=mean(Padj, na.rm=TRUE), .groups="drop") %>%
    dplyr::arrange(dplyr::desc(n_sig), meanPadj)
  readr::write_csv(ora_susres, file.path(dirs$tables_ora, "ORA_replication_scores_sus_vs_res.csv"))
  top_ora_sr <- ora_susres %>% dplyr::group_by(db) %>% dplyr::slice_max(n_sig, n=25) %>% dplyr::ungroup()
  p_ora_sr <- ggplot2::ggplot(top_ora_sr, ggplot2::aes(x=reorder(paste(db, Description, sep="::"), n_sig), y=n_sig, fill=db)) +
    ggplot2::geom_col() + ggplot2::coord_flip() +
    ggplot2::labs(title="ORA term replication in sus vs res (count of significant comparisons)", x="DB::Term", y="# significant comps") +
    ggplot2::theme_minimal(base_size=10)
  ggplot2::ggsave(file.path(dirs$plots_ora, "ORA_replication_bar_sus_vs_res.svg"), p_ora_sr, width=18, height=12)
}

# also list region / layer with most frequent ORA hits
if (nrow(ora_df)) {
  ora_rl_summary <- ora_df %>%
    dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
    dplyr::group_by(db, ID, Description, region, layer) %>%
    dplyr::summarize(n_comp=dplyr::n(), n_sig=sum(Padj<=0.05, na.rm=TRUE), meanPadj=mean(Padj, na.rm=TRUE), .groups="drop") %>%
    dplyr::arrange(dplyr::desc(n_sig), meanPadj)
  readr::write_csv(ora_rl_summary, file.path(dirs$tables_ora, "ORA_region_layer_summary.csv"))
}

# ORA summary for sus vs con
if (nrow(ora_df)) {
  ora_suscon <- ora_df %>%
    dplyr::left_join(comp_class_map, by="context") %>%
    dplyr::filter(comparison_class == "sus vs con") %>%
    dplyr::group_by(db, ID, Description) %>%
    dplyr::summarize(n_comp=dplyr::n(), n_sig=sum(Padj<=0.05, na.rm=TRUE), meanPadj=mean(Padj, na.rm=TRUE), .groups="drop") %>%
    dplyr::arrange(dplyr::desc(n_sig), meanPadj)
  readr::write_csv(ora_suscon, file.path(dirs$tables_ora, "ORA_replication_scores_sus_vs_con.csv"))
  top_ora_sc <- ora_suscon %>% dplyr::group_by(db) %>% dplyr::slice_max(n_sig, n=25) %>% dplyr::ungroup()
  p_ora_sc <- ggplot2::ggplot(top_ora_sc, ggplot2::aes(x=reorder(paste(db, Description, sep="::"), n_sig), y=n_sig, fill=db)) +
    ggplot2::geom_col() + ggplot2::coord_flip() +
    ggplot2::labs(title="ORA term replication in sus vs con (count of significant comparisons)", x="DB::Term", y="# significant comps") +
    ggplot2::theme_minimal(base_size=10)
  ggplot2::ggsave(file.path(dirs$plots_ora, "ORA_replication_bar_sus_vs_con.svg"), p_ora_sc, width=18, height=12)
}

# ORA summary for res vs con
if (nrow(ora_df)) {
  ora_rescon <- ora_df %>%
    dplyr::left_join(comp_class_map, by="context") %>%
    dplyr::filter(comparison_class == "res vs con") %>%
    dplyr::group_by(db, ID, Description) %>%
    dplyr::summarize(n_comp=dplyr::n(), n_sig=sum(Padj<=0.05, na.rm=TRUE), meanPadj=mean(Padj, na.rm=TRUE), .groups="drop") %>%
    dplyr::arrange(dplyr::desc(n_sig), meanPadj)
  readr::write_csv(ora_rescon, file.path(dirs$tables_ora, "ORA_replication_scores_res_vs_con.csv"))
  top_ora_rc <- ora_rescon %>% dplyr::group_by(db) %>% dplyr::slice_max(n_sig, n=25) %>% dplyr::ungroup()
  p_ora_rc <- ggplot2::ggplot(top_ora_rc, ggplot2::aes(x=reorder(paste(db, Description, sep="::"), n_sig), y=n_sig, fill=db)) +
    ggplot2::geom_col() + ggplot2::coord_flip() +
    ggplot2::labs(title="ORA term replication in res vs con (count of significant comparisons)", x="DB::Term", y="# significant comps") +
    ggplot2::theme_minimal(base_size=10)
  ggplot2::ggsave(file.path(dirs$plots_ora, "ORA_replication_bar_res_vs_con.svg"), p_ora_rc, width=18, height=12)
}

# 15b) ORA UpSet plots for sus vs res ------------------------------------------
# Filter sus vs res and significant terms with region/layer resolved
# 1) Derive sus/res indicator from `context`
# Adjust patterns to match actual naming conventions across contexts
# Examples handled:
#   "CA1_slmres vs CA1_slmcon"  -> "sus"
#   "CA1_slmcon vs CA1_slmres"  -> "res"
# If contexts contain other tokens, refine the regexes accordingly.
ora_df2 <- ora_df %>%
  mutate(
    # normalize spaces around 'vs'
    context_norm = str_squish(context),
    # detect which side of 'vs' carries the *_res token
    sus_res = case_when(
      str_detect(context_norm, "_res\\s+vs\\s+.*_con") ~ "sus",
      str_detect(context_norm, "_con\\s+vs\\s+.*_res") ~ "res",
      TRUE ~ NA_character_
    )
  )

# 2) Join mappings and filter significant terms with resolved region/layer and sus_res
sr_rl <- ora_df2 %>%
  left_join(comp_class_map, by = "context") %>%
  filter(comparison_class == "sus vs res") %>%
  left_join(
    rl_both %>% transmute(context, region = region_L, layer = layer_L),
    by = "context"
  ) %>%
  filter(Padj <= 0.05, !is.na(region), !is.na(layer), !is.na(sus_res))

if (nrow(sr_rl) > 0) {

  # 3) REGION UpSet faceted by sus/res
  region_sets <- sr_rl %>%
    distinct(ID, region, sus_res) %>%
    mutate(value = 1L) %>%
    pivot_wider(names_from = region, values_from = value, values_fill = 0L)

  region_cols <- setdiff(colnames(region_sets), c("ID", "sus_res"))

  if (length(region_cols) >= 2) {
    p_region <- upset(
      data = region_sets,
      intersect = region_cols,
      name = "Regions",
      base_annotations = list(
        "Intersection size" = intersection_size(
          mode = "exclusive_intersection",
          text = list(check_overlap = TRUE)
        )
      ),
      set_sizes = upset_set_size(),  # note: no 'text' arg here
      mode = "exclusive_intersection",
      wrap_by = "sus_res"  # facet by contrast direction
    ) +
      ggtitle("ORA sus vs res — Region intersections by contrast") +
      theme(plot.title = element_text(hjust = 0.5))

    fp_region <- file.path(dirs$plots_ora, "ORA_sus_vs_res_region_upset_by_contrast.png")
    ggsave(fp_region, p_region, width = 11, height = 6.5, dpi = 150)
    message("[ora-upset] Saved: ", fp_region)
  }

  # 4) LAYER UpSet with fill by sus/res for intersection bars
  layer_sets <- sr_rl %>%
    distinct(ID, layer, sus_res) %>%
    mutate(value = 1L) %>%
    pivot_wider(names_from = layer, values_from = value, values_fill = 0L)

  layer_cols <- setdiff(colnames(layer_sets), c("ID", "sus_res"))

  if (length(layer_cols) >= 2) {
    p_layer <- upset(
      data = layer_sets,
      intersect = layer_cols,
      name = "Layers",
      base_annotations = list(
        "Intersection size" = intersection_size(
          mode = "exclusive_intersection",
          mapping = aes(fill = sus_res),
          text = list(check_overlap = TRUE)
        )
      ),
      set_sizes = upset_set_size(),
      mode = "exclusive_intersection",
      queries = list(
        # keep legend for fill mapping and avoid over-highlighting the matrix
        upset_query(only_components = c("intersections_matrix", "intersections"),
                    color = "grey60", fill = NA)
      ),
      guides = "over"  # collect legends above set sizes for compactness
    ) +
      scale_fill_manual(values = c(sus = "#D55E00", res = "#0072B2"), na.translate = FALSE) +
      ggtitle("ORA sus vs res — Layer intersections (bar fill = contrast)") +
      theme(plot.title = element_text(hjust = 0.5))

    fp_layer <- file.path(dirs$plots_ora, "ORA_sus_vs_res_layer_upset_colored.png")
    ggsave(fp_layer, p_layer, width = 11, height = 6.5, dpi = 150)
    message("[ora-upset] Saved: ", fp_layer)
  }

  # 5) Region::Layer combined sets, faceted by sus/res
  region_layer_sets <- sr_rl %>%
    mutate(region_layer = paste(region, layer, sep = "::")) %>%
    distinct(ID, region_layer, sus_res) %>%
    mutate(value = 1L) %>%
    pivot_wider(names_from = region_layer, values_from = value, values_fill = 0L)

  rl_cols <- setdiff(colnames(region_layer_sets), c("ID", "sus_res"))

  if (length(rl_cols) >= 2) {
    p_rl <- upset(
      data = region_layer_sets,
      intersect = rl_cols,
      name = "Region::Layer",
      base_annotations = list(
        "Intersection size" = intersection_size(
          mode = "exclusive_intersection",
          mapping = aes(fill = sus_res),
          text = list(check_overlap = TRUE)
        )
      ),
      set_sizes = upset_set_size(),
      mode = "exclusive_intersection",
      wrap_by = "sus_res"
    ) +
      scale_fill_manual(values = c(sus = "#D55E00", res = "#0072B2"), na.translate = FALSE) +
      ggtitle("ORA sus vs res — Region::Layer intersections by contrast") +
      theme(plot.title = element_text(hjust = 0.5))

    fp_rl <- file.path(dirs$plots_ora, "ORA_sus_vs_res_region_layer_upset_by_contrast.png")
    ggsave(fp_rl, p_rl, width = 12, height = 7, dpi = 150)
    message("[ora-upset] Saved: ", fp_rl)
  }

} else {
  message("[ora-upset] No significant ORA terms with region/layer and sus/res resolved.")
}


# 15c) ORA UpSet plots highlighting major recurring terms -----------------------
# Focus on sus vs res, highlight major recurring terms across regions and layers
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(ComplexUpset)
})

# 1) Filter sus vs res significant terms and map region/layer
sr_rl <- ora_df %>%
  left_join(comp_class_map, by = "context") %>%
  filter(comparison_class == "sus vs res") %>%
  left_join(rl_both %>% transmute(context, region = region_L, layer = layer_L),
            by = "context") %>%
  filter(Padj <= 0.05, !is.na(region), !is.na(layer))

stopifnot(nrow(sr_rl) > 0)

# 2) Compute recurrence per term across regions and across layers
term_recur_region <- sr_rl %>%
  distinct(ID, region) %>%
  count(ID, name = "n_regions") %>%
  arrange(desc(n_regions))

term_recur_layer <- sr_rl %>%
  distinct(ID, layer) %>%
  count(ID, name = "n_layers") %>%
  arrange(desc(n_layers))

# Keep descriptions for labeling
term_descr <- sr_rl %>% distinct(ID, Description)

# Helper to pick top terms for annotation (regions)
pick_top_terms_region <- function(ids, recur_tbl, descr_tbl, k = 2) {
  tibble(ID = ids) %>%
    left_join(recur_tbl, by = "ID") %>%
    left_join(descr_tbl, by = "ID") %>%
    arrange(desc(n_regions)) %>%
    slice_head(n = k) %>%
    pull(Description) %>%
    unique() %>%
    { if (length(.) == 0) NA_character_ else paste(., collapse = " | ") }
}

# Helper to pick top terms for annotation (layers)
pick_top_terms_layer <- function(ids, recur_tbl, descr_tbl, k = 2) {
  tibble(ID = ids) %>%
    left_join(recur_tbl, by = "ID") %>%
    left_join(descr_tbl, by = "ID") %>%
    arrange(desc(n_layers)) %>%
    slice_head(n = k) %>%
    pull(Description) %>%
    unique() %>%
    { if (length(.) == 0) NA_character_ else paste(., collapse = " | ") }
}

# 3) REGION UPSET: Reshape for ComplexUpset
region_sets <- sr_rl %>%
  distinct(ID, region) %>%
  mutate(value = 1L) %>%
  pivot_wider(names_from = region, values_from = value, values_fill = 0L)

region_cols <- setdiff(names(region_sets), "ID")

# Compute intersection labels
int_labels_region <- region_sets %>%
  unite("pattern", all_of(region_cols), sep = "", remove = FALSE) %>%
  group_by(pattern) %>%
  summarize(
    ids = list(ID),
    size = n(),
    top_terms = pick_top_terms_region(ID, term_recur_region, term_descr, k = 2),
    .groups = "drop"
  ) %>%
  filter(size >= 3)

p_region <- upset(
  data = region_sets,
  intersect = region_cols,
  mode = "exclusive_intersection",
  min_size = 3,
  sort_intersections_by = c("cardinality", "degree"),
  sort_intersections = "descending",
  base_annotations = list(
    "Intersection size" = intersection_size(
      mode = "exclusive_intersection",
      text = list(check_overlap = TRUE)
    )
  ),
  set_sizes = upset_set_size()
) +
  ggtitle("sus vs res — GO/KEGG overlaps across regions") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(dirs$plots_ora, "ORA_sus_vs_res_region_upset_major_terms.svg"),
       plot = p_region, width = 12, height = 7, bg = "white")

# 4) LAYER UPSET: Reshape for ComplexUpset
layer_sets <- sr_rl %>%
  distinct(ID, layer) %>%
  mutate(value = 1L) %>%
  pivot_wider(names_from = layer, values_from = value, values_fill = 0L)

layer_cols <- setdiff(names(layer_sets), "ID")

# Compute intersection labels
int_labels_layer <- layer_sets %>%
  unite("pattern", all_of(layer_cols), sep = "", remove = FALSE) %>%
  group_by(pattern) %>%
  summarize(
    ids = list(ID),
    size = n(),
    top_terms = pick_top_terms_layer(ID, term_recur_layer, term_descr, k = 2),
    .groups = "drop"
  ) %>%
  filter(size >= 3)

p_layer <- upset(
  data = layer_sets,
  intersect = layer_cols,
  mode = "exclusive_intersection",
  min_size = 3,
  sort_intersections_by = c("cardinality", "degree"),
  sort_intersections = "descending",
  base_annotations = list(
    "Intersection size" = intersection_size(
      mode = "exclusive_intersection",
      text = list(check_overlap = TRUE)
    )
  ),
  set_sizes = upset_set_size()
) +
  ggtitle("sus vs res — GO/KEGG overlaps across layers") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(dirs$plots_ora, "ORA_sus_vs_res_layer_upset_major_terms.svg"),
       plot = p_layer, width = 12, height = 7, bg = "white")

# 5) Export ranked major terms (recurrence tables)
readr::write_csv(
  term_recur_region %>% left_join(term_descr, by = "ID"),
  file.path(dirs$tables_ora, "ORA_sus_vs_res_major_terms_by_region.csv")
)
readr::write_csv(
  term_recur_layer %>% left_join(term_descr, by = "ID"),
  file.path(dirs$tables_ora, "ORA_sus_vs_res_major_terms_by_layer.csv")
)

# 6) Export intersection summaries with top terms
readr::write_csv(
  int_labels_region,
  file.path(dirs$tables_ora, "ORA_sus_vs_res_region_intersections.csv")
)
readr::write_csv(
  int_labels_layer,
  file.path(dirs$tables_ora, "ORA_sus_vs_res_layer_intersections.csv")
)

# 7) Region-specific and layer-specific terms (unique to one location)
# Region-specific terms - FIXED with dplyr::select()
region_specific <- sr_rl %>%
  group_by(ID, Description) %>%
  summarize(regions = list(unique(region)),
            n_regions = n_distinct(region),
            layers = list(unique(layer)),
            .groups = "drop") %>%
  filter(n_regions == 1) %>%
  mutate(region = sapply(regions, `[`, 1)) %>%
  dplyr::select(ID, Description, region, layers)

readr::write_csv(
  region_specific,
  file.path(dirs$tables_ora, "ORA_sus_vs_res_region_specific_terms.csv")
)

# Layer-specific terms - FIXED with dplyr::select()
layer_specific <- sr_rl %>%
  group_by(ID, Description) %>%
  summarize(layers = list(unique(layer)),
            n_layers = n_distinct(layer),
            regions = list(unique(region)),
            .groups = "drop") %>%
  filter(n_layers == 1) %>%
  mutate(layer = sapply(layers, `[`, 1)) %>%
  dplyr::select(ID, Description, layer, regions)

readr::write_csv(
  layer_specific,
  file.path(dirs$tables_ora, "ORA_sus_vs_res_layer_specific_terms.csv")
)

# 8) Bar plots for region/layer-specific terms (top 15 each)
if (nrow(region_specific) > 0) {
  top_region_spec <- region_specific %>%
    group_by(region) %>%
    slice_head(n = 15) %>%
    ungroup()
  
  p_reg_spec <- ggplot(top_region_spec, aes(x = reorder(Description, region), y = 1, fill = region)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~ region, scales = "free_y", ncol = 2) +
    labs(title = "Region-specific ORA terms (sus vs res)", x = "Term", y = "") +
    theme_minimal(10) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  ggsave(file.path(dirs$plots_ora, "ORA_sus_vs_res_region_specific_bar.svg"),
         plot = p_reg_spec, width = 14, height = 10, bg = "white")
}

if (nrow(layer_specific) > 0) {
  top_layer_spec <- layer_specific %>%
    group_by(layer) %>%
    slice_head(n = 15) %>%
    ungroup()
  
  p_lay_spec <- ggplot(top_layer_spec, aes(x = reorder(Description, layer), y = 1, fill = layer)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~ layer, scales = "free_y", ncol = 2) +
    labs(title = "Layer-specific ORA terms (sus vs res)", x = "Term", y = "") +
    theme_minimal(10) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  ggsave(file.path(dirs$plots_ora, "ORA_sus_vs_res_layer_specific_bar.svg"),
         plot = p_lay_spec, width = 14, height = 10, bg = "white")
}


# create GO overrepresentation analysis plots for con 

# 16) Executive bundle + artifact check ---------------------------------------
readr::write_csv(top_themes,  file.path(dirs$tables_gsea, "top_themes_per_db.csv"))
readr::write_csv(theme_comp %>% dplyr::arrange(db, theme, dplyr::desc(abs(meanNES))), file.path(dirs$tables_gsea, "theme_by_comparison_sorted.csv"))

md <- c(
  "# Consolidated Summary",
  paste0("- Comparisons scanned: ", length(comp_dirs)),
  paste0("- GSEA rows: ", nrow(gsea_df)),
  "",
  "Artifacts:",
  "* Plots/Anatomy/localization_<DB>.svg",
  "* Plots/Class/localization_<DB>_by_comparison_class_labeled.svg",
  "* Plots/Conditions/localization_<DB>_tested_condition.svg",
  "* Plots/Conditions/localization_<DB>_baseline_condition.svg",
  "* Plots/ComparisonMatrices/theme_<THEME>_<DB>_comparison_matrix.svg",
  "* Tables/context_region_layer_condition_both_sides.csv",
  "* Tables/Class/context_comparison_class_map.csv",
  "* Tables/GSEA/localization_table.csv",
  "* Tables/Class/localization_by_comparison_class_table.csv",
  "* Tables/Class/localization_by_comparison_class_labeled.csv",
  "* Tables/Conditions/localization_left_condition_table.csv",
  "* Tables/Conditions/localization_right_condition_table.csv",
  "* Tables/Class/comparison_level_long.csv",
  "* Tables/QA/localization_NA_diagnostics.csv",
  "* Tables/GSEA/all_GSEA_long.csv / themes_table.csv / theme_by_comparison.csv / theme_replication_scores.csv",
  if (nrow(ora_df)) "* Tables/ORA/all_ORA_long.csv / ORA_replication_scores.csv; Plots/ORA/ORA_replication_bar.svg" else "* (no ORA detected)"
)
writeLines(md, con = file.path(out_dir, "CONSOLIDATED_SUMMARY.md"))

# Artifact existence check across subfolders
expected <- c(
  file.path(dirs$plots_anat,  paste0("localization_", unique(theme_loc$db), ".svg")),
  file.path(dirs$plots_class, paste0("localization_", unique(theme_loc_by_class_lab$db), "_by_comparison_class_labeled.svg")),
  file.path(dirs$plots_cond,  paste0("localization_", unique(theme_loc_left_cond$db), "_tested_condition.svg")),
  file.path(dirs$plots_cond,  paste0("localization_", unique(theme_loc_right_cond$db), "_baseline_condition.svg"))
)
exists_tbl <- tibble::tibble(file = unique(expected), exists = file.exists(unique(expected)))
readr::write_csv(exists_tbl, file.path(out_dir, "artifact_existence_check.csv"))
print(exists_tbl)

message("Consolidation complete. See ZZ_consolidated for outputs.")

























































# =======================
# CONSOLIDATION + SUMMARY
# =======================
# Primary: GSEA (gseGO/gseKEGG/Reactome) using NES/FDR
# Secondary: ORA (separate tables/plots)
# Includes class-labeled heatmaps (SR/RC/SC:n) and comparison-level matrices. [web:22]

# 0) Setup --------------------------------------------------------------------
pkgs <- c("dplyr","tidyr","readr","stringr","purrr","ggplot2","igraph","ggraph", "patchwork", "tibble", "clusterProfiler", "org.Mm.eg.db", "DOSE", "enrichplot", "ReactomePA")
missing_pkgs <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(missing_pkgs)) {
  message("Installing missing packages: ", paste(missing_pkgs, collapse=", "))
  tryCatch({
    install.packages(missing_pkgs, repos = c("https://cloud.r-project.org"))
  }, error = function(e) {
    warning("Primary CRAN failed: ", conditionMessage(e), " — retrying ggraph via r-universe")
    if ("ggraph" %in% missing_pkgs) install.packages("ggraph", repos = c("https://thomasp85.r-universe.dev","https://cloud.r-project.org"))
  })
}
suppressPackageStartupMessages({
  for (p in pkgs) {
    tryCatch(library(p, character.only = TRUE),
             error = function(e) warning(sprintf("Failed to load %s: %s", p, conditionMessage(e))))
  }
})

# Root paths (fixed) ----------------------------------------------------------
root_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Results/neuron-phenotypeWithinUnit"
out_dir  <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Results/ZZ_consolidated"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# New: subfolder layout for tables/plots --------------------------------------
dirs <- list(
  tables        = file.path(out_dir, "Tables"),
  tables_gsea   = file.path(out_dir, "Tables", "GSEA"),
  tables_ora    = file.path(out_dir, "Tables", "ORA"),
  tables_class  = file.path(out_dir, "Tables", "Class"),
  tables_cond   = file.path(out_dir, "Tables", "Conditions"),
  tables_qa     = file.path(out_dir, "Tables", "QA"),
  plots         = file.path(out_dir, "Plots"),
  plots_anat    = file.path(out_dir, "Plots", "Anatomy"),
  plots_class   = file.path(out_dir, "Plots", "Class"),
  plots_cond    = file.path(out_dir, "Plots", "Conditions"),
  plots_cm      = file.path(out_dir, "Plots", "ComparisonMatrices"),
  plots_ora     = file.path(out_dir, "Plots", "ORA")
)
invisible(lapply(dirs, function(d) dir.create(d, showWarnings = FALSE, recursive = TRUE)))  # recursive creation

# --- MODULE BRANCH: paths + filename parser ---
dirs$tables_modules <- file.path(out_dir, "Tables", "Modules")
dirs$plots_modules  <- file.path(out_dir, "Plots", "Modules")
invisible(lapply(c(dirs$tables_modules, dirs$plots_modules), dir.create, recursive=TRUE, showWarnings=FALSE))

modules_root <- file.path(dirname(root_dir), "20_modules")
if (!dir.exists(modules_root)) warning("modules_root not found: ", modules_root)

parse_module_fname <- function(fp) {
  bn <- basename(fp); noext <- sub("\\.csv$", "", bn, ignore.case = TRUE)
  m1 <- regexec("^gse([A-Za-z0-9]+)_([A-Za-z0-9]+)_([^_]+)_(.+)$", noext, perl=TRUE); r1 <- regmatches(noext, m1)[[1]]
  if (length(r1)) {
    db1 <- r1[2]; db2 <- r1[3]; db <- if (toupper(db1)=="GO" && toupper(db2)=="BP") "GO_BP" else paste0(db1,"_",db2)
    mod <- r1[4]; rem <- r1[5]
    gpos <- gregexpr("_(CA1|CA2|CA3|DG)", rem, perl=TRUE)[[1]]; split_pos <- if (length(gpos) && gpos[1] != -1) tail(gpos,1) else -1
    if (split_pos == -1) { toks <- strsplit(rem, "_")[[1]]; left <- paste(toks[1:(length(toks)-1)], collapse="_"); right <- toks[length(toks)] } else {
      left <- substr(rem, 1, split_pos-1); right <- substr(rem, split_pos+1, nchar(rem))
    }
    return(list(source="gsea", db=db, module=mod, context=paste0(left," vs ", right)))
  }
  m2 <- regexec("^ORA_([A-Za-z0-9]+)_([^_]+)_(.+)$", noext, perl=TRUE); r2 <- regmatches(noext, m2)[[1]]
  if (length(r2)) {
    db2 <- r2[2]; db <- if (toupper(db2)=="BP") "GO_BP" else db2; mod <- r2[3]; rem <- r2[4]
    gpos <- gregexpr("_(CA1|CA2|CA3|DG)", rem, perl=TRUE)[[1]]; split_pos <- if (length(gpos) && gpos[1] != -1) tail(gpos,1) else -1
    if (split_pos == -1) { toks <- strsplit(rem, "_")[[1]]; left <- paste(toks[1:(length(toks)-1)], collapse="_"); right <- toks[length(toks)] } else {
      left <- substr(rem, 1, split_pos-1); right <- substr(rem, split_pos+1, nchar(rem))
    }
    return(list(source="ora_deg", db=db, module=mod, context=paste0(left," vs ", right)))
  }
  m3 <- regexec("^ORA_fullModule_([A-Za-z0-9]+)_([^_]+)_(.+)$", noext, perl=TRUE); r3 <- regmatches(noext, m3)[[1]]
  if (length(r3)) {
    db2 <- r3[2]; db <- if (toupper(db2)=="BP") "GO_BP" else db2; mod <- r3[3]; rem <- r3[4]
    gpos <- gregexpr("_(CA1|CA2|CA3|DG)", rem, perl=TRUE)[[1]]; split_pos <- if (length(gpos) && gpos[1] != -1) tail(gpos,1) else -1
    if (split_pos == -1) { toks <- strsplit(rem, "_")[[1]]; left <- paste(toks[1:(length(toks)-1)], collapse="_"); right <- toks[length(toks)] } else {
      left <- substr(rem, 1, split_pos-1); right <- substr(rem, split_pos+1, nchar(rem))
    }
    return(list(source="ora_full", db=db, module=mod, context=paste0(left," vs ", right)))
  }
  NULL
}

# 1) Helpers ------------------------------------------------------------------
read_if <- function(fp) {
  if (!file.exists(fp)) return(NULL)
  tryCatch(readr::read_csv(fp, guess_max = 10000, progress = FALSE, show_col_types = FALSE),
           error = function(e) { warning("Failed reading: ", fp, " — ", conditionMessage(e)); NULL })
}
read_context <- function(comp_dir) {
  readme <- file.path(comp_dir, "README.txt")
  if (file.exists(readme)) {
    rl <- readLines(readme, warn = FALSE)
    ctx <- rl[grepl("^Context: ", rl)]
    if (length(ctx)) return(sub("^Context:\\s*", "", ctx[1]))
  }
  basename(comp_dir)
}
`%||%` <- function(a,b) if (!is.null(a) && length(a)>0 && !is.na(a)) a else b
tokenize <- function(x){ x <- tolower(x); x <- gsub("[^a-z0-9]+"," ",x); unique(unlist(strsplit(x,"\\s+"))) }
jaccard  <- function(a,b){ ia<-tokenize(a); ib<-tokenize(b); if(!length(ia)||!length(ib)) return(0); length(intersect(ia,ib))/length(union(ia,ib)) }

# 2) GSEA collectors -----------------------------------------------------------
parse_gsea_tbl <- function(fp, db) {
  df <- read_if(fp); if (is.null(df) || !nrow(df)) return(NULL)
  nm <- names(df); has <- function(x) any(x %in% nm)
  if (!has("ID") || !has("Description")) return(NULL)
  if (!has("NES")) return(NULL)
  fdr_col <- c("p.adjust","qvalue","padj","FDR","adj.P.Val")[c("p.adjust","qvalue","padj","FDR","adj.P.Val") %in% nm][1]
  if (is.na(fdr_col)) return(NULL)
  out <- tibble::tibble(
    db = rep(db, nrow(df)),
    ID = df[["ID"]],
    Description = df[["Description"]],
    NES = suppressWarnings(as.numeric(df[["NES"]])),
    FDR = suppressWarnings(as.numeric(df[[fdr_col]]))
  ) %>% dplyr::filter(!is.na(ID), !is.na(NES), !is.na(FDR))
  if (!nrow(out)) return(NULL)
  out
}
empty_gsea <- function(){ tibble::tibble(db=character(), ID=character(), Description=character(), NES=double(), FDR=double()) }

collect_gsea_one <- function(comp_dir) {
  paths <- list(
    GO_BP    = file.path(comp_dir, "10_global","GO_BP","Tables"),
    KEGG     = file.path(comp_dir, "10_global","KEGG","Tables"),
    Reactome = file.path(comp_dir, "10_global","Reactome","Tables")
  )
  if (!any(dir.exists(unlist(paths)))) return(empty_gsea())

  gsego   <- if (dir.exists(paths$GO_BP))    list.files(paths$GO_BP,    pattern="^gseGO_BP_.*\\.csv$",      full.names=TRUE) else character(0)
  gsekegg <- if (dir.exists(paths$KEGG))     list.files(paths$KEGG,     pattern="^gseKEGG_.*\\.csv$",       full.names=TRUE) else character(0)
  rs      <- if (dir.exists(paths$Reactome)) list.files(paths$Reactome, pattern="^ReactomeGSEA_.*\\.csv$",  full.names=TRUE) else character(0)

  if (length(gsego)+length(gsekegg)+length(rs) == 0) return(empty_gsea())

  parts <- list(
    dplyr::bind_rows(lapply(gsego,   parse_gsea_tbl, db="GO_BP")),
    dplyr::bind_rows(lapply(gsekegg, parse_gsea_tbl, db="KEGG")),
    dplyr::bind_rows(lapply(rs,      parse_gsea_tbl, db="Reactome"))
  )
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (!length(parts)) return(empty_gsea())
  dplyr::bind_rows(parts)
}

# 3) ORA collector -------------------------------------------------------------
parse_ora_tbl <- function(fp, db) {
  df <- read_if(fp); if (is.null(df) || !nrow(df)) return(NULL)
  nm <- names(df); if (!("ID" %in% nm && "Description" %in% nm)) return(NULL)
  fdr_col <- c("p.adjust","qvalue","padj","FDR","adj.P.Val")[c("p.adjust","qvalue","padj","FDR","adj.P.Val") %in% nm][1]
  if (is.na(fdr_col)) return(NULL)
  out <- tibble::tibble(db=rep(db, nrow(df)), ID=df[["ID"]], Description=df[["Description"]], Padj=suppressWarnings(as.numeric(df[[fdr_col]]))) %>%
    dplyr::filter(!is.na(ID), !is.na(Padj))
  if (!nrow(out)) return(NULL)
  out
}
empty_ora <- function(){ tibble::tibble(db=character(), ID=character(), Description=character(), Padj=double()) }

collect_ora_one <- function(comp_dir) {
  paths <- list(
    GO_BP    = file.path(comp_dir, "10_global","GO_BP","Tables"),
    KEGG     = file.path(comp_dir, "10_global","KEGG","Tables"),
    Reactome = file.path(comp_dir, "10_global","Reactome","Tables")
  )
  if (!any(dir.exists(unlist(paths)))) return(empty_ora())

  files <- c(
    if (dir.exists(paths$GO_BP))    list.files(paths$GO_BP,    pattern="^ORA_.*\\.csv$",         full.names=TRUE) else character(0),
    if (dir.exists(paths$KEGG))     list.files(paths$KEGG,     pattern="^ORA_.*\\.csv$",         full.names=TRUE) else character(0),
    if (dir.exists(paths$Reactome)) list.files(paths$Reactome, pattern="^ReactomeORA_.*\\.csv$", full.names=TRUE) else character(0)
  )
  if (!length(files)) return(empty_ora())

  binders <- lapply(files, function(fp) {
    db <- if (grepl("/GO_BP/", fp)) "GO_BP" else if (grepl("/KEGG/", fp)) "KEGG" else if (grepl("/Reactome/", fp)) "Reactome" else "Unknown"
    parse_ora_tbl(fp, db = db)
  })
  binders <- binders[!vapply(binders, is.null, logical(1))]
  if (!length(binders)) return(empty_ora())
  dplyr::bind_rows(binders)
}

# 4) Scan comparison folders ---------------------------------------------------
comp_dirs <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
comp_dirs <- comp_dirs[dir.exists(file.path(comp_dirs, "10_global"))]
message(sprintf("Found %d comparison folders under neuron-phenotypeWithinUnit.", length(comp_dirs)))
if (!length(comp_dirs)) stop("No comparison folders with '10_global' found under the fixed root_dir.", call. = FALSE)

# 5) GSEA assembly -------------------------------------------------------------
gsea_list <- lapply(comp_dirs, function(cd) tibble::tibble(comp_dir = cd, context = read_context(cd), data = list(collect_gsea_one(cd))))
gsea_df <- purrr::map_dfr(gsea_list, function(x) { dx <- x$data[[1]]; if (!is.null(dx) && nrow(dx) > 0) cbind(x[1:2], dx) else NULL })
if (!nrow(gsea_df)) stop("No valid GSEA rows found across comparisons.", call. = FALSE)

req <- c("comp_dir","context","db","ID","Description","NES","FDR")
missing_cols <- setdiff(req, names(gsea_df))
if (length(missing_cols)) stop(sprintf("gsea_df missing columns: %s", paste(missing_cols, collapse=", ")), call. = FALSE)

gsea_df <- gsea_df %>% dplyr::mutate(comp_dir=as.character(comp_dir), context=as.character(context), db=as.character(db), ID=as.character(ID), Description=as.character(Description), NES=as.numeric(NES), FDR=as.numeric(FDR))
readr::write_csv(gsea_df, file.path(dirs$tables_gsea, "all_GSEA_long.csv"))  # moved [web:23]
message(sprintf("GSEA rows: %d", nrow(gsea_df)))

# 6) ORA assembly --------------------------------------------------------------
ora_list <- lapply(comp_dirs, function(cd) tibble::tibble(comp_dir = cd, context = read_context(cd), data = list(collect_ora_one(cd))))
ora_df <- purrr::map_dfr(ora_list, function(x) { dx <- x$data[[1]]; if (!is.null(dx) && nrow(dx) > 0) cbind(x[1:2], dx) else NULL })
if (nrow(ora_df)) readr::write_csv(ora_df, file.path(dirs$tables_ora, "all_ORA_long.csv"))

# 6b) Module-level collectors and assembly --------------------------------
module_dirs <- if (dir.exists(modules_root)) list.dirs(modules_root, full.names=TRUE, recursive=FALSE) else character(0)
module_dirs <- module_dirs[grepl("/Module_", gsub("\\\\","/", module_dirs))]
message(sprintf("Found %d Module_* folders.", length(module_dirs)))

list_module_files <- function(module_dir) {
  tbl_dir <- file.path(module_dir, "Tables")
  if (!dir.exists(tbl_dir)) return(character(0))
  subdirs <- c("gsea","ora_deg","ora_full")
  unlist(lapply(subdirs, function(s) {
    d <- file.path(tbl_dir, s)
    if (dir.exists(d)) list.files(d, pattern="\\.csv$", full.names=TRUE) else character(0)
  }))
}

parse_module_csv <- function(fp) {
  meta <- parse_module_fname(fp); if (is.null(meta)) return(NULL)
  df <- read_if(fp); if (is.null(df) || !nrow(df)) return(NULL)
  nm <- names(df)
  if (identical(meta$source, "gsea")) {
    fdr_col <- c("p.adjust","qvalue","padj","FDR","adj.P.Val")[c("p.adjust","qvalue","padj","FDR","adj.P.Val") %in% nm][1]
    if (!all(c("ID","Description","NES") %in% nm) || is.na(fdr_col)) return(NULL)
    return(
      tibble::tibble(
        module = meta$module, db = meta$db, context = meta$context,
        ID = as.character(df$ID), Description = as.character(df$Description),
        NES = suppressWarnings(as.numeric(df$NES)), FDR = suppressWarnings(as.numeric(df[[fdr_col]]))
      ) %>% dplyr::filter(!is.na(ID), is.finite(NES), is.finite(FDR))
    )
  } else {
    padj_col <- c("p.adjust","qvalue","padj","FDR","adj.P.Val")[c("p.adjust","qvalue","padj","FDR","adj.P.Val") %in% nm][1]
    if (!all(c("ID","Description") %in% nm) || is.na(padj_col)) return(NULL)
    return(
      tibble::tibble(
        module = meta$module, db = meta$db, context = meta$context,
        ID = as.character(df$ID), Description = as.character(df$Description),
        Padj = suppressWarnings(as.numeric(df[[padj_col]])), source = meta$source
      ) %>% dplyr::filter(!is.na(ID), is.finite(Padj))
    )
  }
}

module_files <- unlist(lapply(module_dirs, list_module_files))
message(sprintf("Module CSV files detected: %d", length(module_files)))

mod_parsed <- lapply(module_files, parse_module_csv)
mod_parsed <- mod_parsed[!vapply(mod_parsed, is.null, logical(1))]

module_gsea_df <- purrr::map_dfr(mod_parsed, ~{ if (all(c("NES","FDR") %in% names(.x))) .x else NULL })
module_ora_df  <- purrr::map_dfr(mod_parsed, ~{ if ("Padj" %in% names(.x)) .x else NULL })

if (nrow(module_gsea_df)) {
  readr::write_csv(module_gsea_df, file.path(dirs$tables_modules, "module_GSEA_long.csv"))
  message(sprintf("Module GSEA rows: %d", nrow(module_gsea_df)))
} else message("No module GSEA rows.")

if (nrow(module_ora_df)) {
  readr::write_csv(module_ora_df, file.path(dirs$tables_modules, "module_ORA_long.csv"))
  message(sprintf("Module ORA rows: %d", nrow(module_ora_df)))
} else message("No module ORA rows.")

# 7) Theme building (GSEA only) -----------------------------------------------
build_themes_for_db <- function(df_db, top_n_terms = 200, sim_cut = 0.25) {
  if (is.null(df_db) || !nrow(df_db)) return(NULL)
  df_db <- dplyr::mutate(df_db, sig = FDR <= 0.05)

  agg <- df_db %>%
    dplyr::group_by(ID, Description, db) %>%
    dplyr::summarize(n_sig = sum(sig, na.rm = TRUE), n = dplyr::n(), meanNES = mean(NES, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(n_sig), dplyr::desc(abs(meanNES)))
  sel <- head(agg, top_n_terms)
  if (nrow(sel) < 2) return(dplyr::mutate(sel, theme = paste0(unique(df_db$db), "_", dplyr::row_number())))

  ids <- sel$ID; descs <- sel$Description; n <- length(ids)
  sim <- matrix(0, n, n); rownames(sim) <- colnames(sim) <- ids
  for (i in seq_len(n)) for (j in i:n) { s <- if (i==j) 1 else jaccard(descs[i], descs[j]); sim[i,j] <- s; sim[j,i] <- s }

  edge_df <- as.data.frame(as.table(sim), stringsAsFactors = FALSE) |>
    dplyr::rename(ID1 = Var1, ID2 = Var2, sim = Freq) |>
    dplyr::filter(ID1 != ID2, sim >= sim_cut)

  g <- igraph::make_empty_graph(directed = FALSE) |> igraph::add_vertices(n, name = ids)
  if (nrow(edge_df) > 0) {
    edges_mat <- unique(t(apply(edge_df[,c("ID1","ID2")], 1, function(x) sort(as.character(x)))))
    g <- igraph::add_edges(g, as.vector(t(edges_mat)))
  }
  if (igraph::vcount(g) == 0) { sel$theme <- paste0(unique(df_db$db), "_", seq_len(nrow(sel))); return(sel) }

  cl <- igraph::cluster_louvain(g)
  memb <- igraph::membership(cl)
  theme_map <- tibble::tibble(ID = names(memb), theme = paste0(unique(df_db$db), "_C", as.integer(memb)))
  out <- dplyr::left_join(sel, theme_map, by = "ID")
  if (any(is.na(out$theme))) out$theme[is.na(out$theme)] <- paste0(unique(df_db$db), "_", seq_len(sum(is.na(out$theme))))
  out
}

gsea_core <- gsea_df %>% dplyr::select(db, ID, Description, NES, FDR) %>% dplyr::mutate(db=as.character(db), ID=as.character(ID), Description=as.character(Description), NES=as.numeric(NES), FDR=as.numeric(FDR))
req_theme <- c("db","ID","Description","NES","FDR"); miss_theme <- setdiff(req_theme, names(gsea_core)); if (length(miss_theme)) stop(sprintf("gsea_core missing: %s", paste(miss_theme, collapse=", ")), call.=FALSE)
split_dbs <- split(gsea_core, gsea_core$db)
themes_db <- lapply(split_dbs, build_themes_for_db)
themes_db <- themes_db[!vapply(themes_db, is.null, logical(1))]
themes_tbl <- if (length(themes_db)) dplyr::bind_rows(themes_db) else tibble::tibble()
readr::write_csv(themes_tbl, file.path(dirs$tables_gsea, "themes_table.csv"))

# BP membership per theme (cluster), per DB
bp_members <- gsea_df %>%
  dplyr::filter(db == "GO_BP") %>%                                    # only BP here
  dplyr::inner_join(dplyr::select(themes_tbl, ID, theme), by="ID") %>%# map ID -> theme
  dplyr::group_by(db, theme, ID, Description) %>%
  dplyr::summarize(
    n_hits = dplyr::n(),
    meanNES = mean(NES, na.rm=TRUE),
    bestFDR = suppressWarnings(min(FDR, na.rm=TRUE)),
    .groups="drop"
  ) %>%
  dplyr::arrange(db, theme, dplyr::desc(n_hits), bestFDR)

# Save a single long CSV
readr::write_csv(bp_members, file.path(dirs$tables_gsea, "GO_BP_terms_by_theme_long.csv"))

# Optionally export one CSV per theme
bp_split <- split(bp_members, interaction(bp_members$db, bp_members$theme, drop=TRUE))
invisible(purrr::iwalk(bp_split, function(df, nm){
  parts <- strsplit(nm, "\\.")[[1]]
  dbn <- parts[1]; th <- parts[2]
  fp <- file.path(dirs$tables_gsea, paste0("GO_BP_terms_", th, ".csv"))
  readr::write_csv(df, fp)
}))

# 8) Aggregation per comparison x theme (GSEA) --------------------------------
assign_themes <- if (nrow(themes_tbl)) dplyr::inner_join(gsea_df, dplyr::select(themes_tbl, ID, theme), by = "ID") else tibble::tibble()
if (!nrow(assign_themes)) stop("No themes constructed (themes_tbl empty). Consider lowering sim_cut or increasing top_n_terms.", call. = FALSE)

theme_comp <- assign_themes %>%
  dplyr::group_by(db, theme, comp_dir, context) %>%
  dplyr::summarize(n_terms=dplyr::n(), meanNES=mean(NES, na.rm=TRUE), minFDR=suppressWarnings(min(FDR, na.rm=TRUE)), frac_sig=mean(FDR<=0.05, na.rm=TRUE), .groups="drop")
readr::write_csv(theme_comp, file.path(dirs$tables_gsea, "theme_by_comparison.csv"))

# Masking for plots
theme_comp$NES_masked <- ifelse(theme_comp$minFDR <= 0.05, theme_comp$meanNES, NA_real_)

# 9) Robust parsing — both sides (fixed) --------------------------------------
parse_both <- function(ctx) {
  parts <- unlist(strsplit(ctx, "(?:\\s+vs\\s+|_vs_)", perl=TRUE))
  parts <- trimws(parts)
  left  <- if (length(parts) >= 1) parts[1] else NA_character_
  right <- if (length(parts) >= 2) parts[2] else NA_character_

  region_pat <- "(DG|CA1|CA2|CA3)"
  base_layer <- "(mo|po|sr|sp|so|slm|sg)"
  cond_pat   <- "(res|sus|con)"

  parse_side <- function(side) {
    if (is.na(side) || side == "") return(list(region=NA_character_, layer=NA_character_, condition=NA_character_))
    m <- regexec(paste0("^", region_pat, "[_-]?", base_layer, "(", cond_pat, ")?$"), side, perl=TRUE)
    r <- regmatches(side, m)[[1]]
    if (length(r) == 0) {
      m2 <- regexec(paste0("^", region_pat, "[_-]?", base_layer, "(", cond_pat, ")?"), side, perl=TRUE)
      r  <- regmatches(side, m2)[[1]]
    }
    region <- NA_character_; layer <- NA_character_; condition <- NA_character_
    if (length(r) >= 3) {
      region    <- r[2]
      layer_raw <- tolower(r[3]); layer <- layer_raw
      if (length(r) >= 5 && nzchar(r[5])) condition <- tolower(r[5])
    } else {
      toks <- unlist(strsplit(side, "[_-]+"))
      region <- if (length(toks) >= 1 && grepl(paste0("^", region_pat, "$"), toks[1])) toks[1] else NA_character_
      layer  <- if (length(toks) >= 2 && grepl(paste0("^", base_layer, "$"), tolower(toks[2]))) tolower(toks[2]) else NA_character_
      condition <- if (length(toks) >= 3 && grepl(paste0("^", cond_pat, "$"), tolower(toks[3]))) tolower(toks[3]) else NA_character_
    }
    list(region=region, layer=layer, condition=condition)
  }

  L <- parse_side(left)
  R <- parse_side(right)
  list(L=L, R=R)
}

rl_both <- unique(theme_comp$context) %>%
  purrr::map_df(~{
    pr <- parse_both(.x)
    tibble::tibble(
      context   = .x,
      region_L  = pr$L$region, layer_L  = pr$L$layer, condition_L  = pr$L$condition,
      region_R  = pr$R$region, layer_R  = pr$R$layer, condition_R  = pr$R$condition
    )
  })
readr::write_csv(rl_both, file.path(dirs$tables, "context_region_layer_condition_both_sides.csv"))

# 10) Comparison classes (from parsed conditions) -----------------------------
comp_class_map <- rl_both %>%
  dplyr::transmute(
    context,
    cond_left  = condition_L,
    cond_right = condition_R,
    comparison_class = dplyr::case_when(
      cond_left == "sus" & cond_right == "res" ~ "sus vs res",
      cond_left == "res" & cond_right == "con" ~ "res vs con",
      cond_left == "sus" & cond_right == "con" ~ "sus vs con",
      TRUE ~ NA_character_
    )
  )
readr::write_csv(comp_class_map, file.path(dirs$tables_class, "context_comparison_class_map.csv"))

# 10b) Module context + class maps ----------------------------------------
if (exists("module_gsea_df") && nrow(module_gsea_df)) {
  mod_ctx <- unique(module_gsea_df$context)
} else if (exists("module_ora_df") && nrow(module_ora_df)) {
  mod_ctx <- unique(module_ora_df$context)
} else mod_ctx <- character(0)

if (length(mod_ctx)) {
  mod_rl_both <- tibble::tibble(context = mod_ctx) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(parsed = list(parse_both(context))) %>%
    dplyr::mutate(
      region_L  = parsed$L$region, layer_L  = parsed$L$layer, condition_L  = parsed$L$condition,
      region_R  = parsed$R$region, layer_R  = parsed$R$layer, condition_R  = parsed$R$condition
    ) %>% dplyr::ungroup() %>% dplyr::select(-parsed)
  readr::write_csv(mod_rl_both, file.path(dirs$tables_modules, "module_context_region_layer_condition_both_sides.csv"))

  mod_comp_class_map <- mod_rl_both %>%
    dplyr::transmute(
      context,
      cond_left  = condition_L,
      cond_right = condition_R,
      comparison_class = dplyr::case_when(
        cond_left == "sus" & cond_right == "res" ~ "sus vs res",
        cond_left == "res" & cond_right == "con" ~ "res vs con",
        cond_left == "sus" & cond_right == "con" ~ "sus vs con",
        TRUE ~ NA_character_
      )
    )
  readr::write_csv(mod_comp_class_map, file.path(dirs$tables_modules, "module_context_comparison_class_map.csv"))
}

# 11) Localization tables ------------------------------------------------------
theme_loc <- theme_comp %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
  dplyr::filter(!is.na(region), !is.na(layer)) %>%
  dplyr::group_by(db, theme, region, layer) %>%
  dplyr::summarize(meanNES=mean(meanNES, na.rm=TRUE), frac_sig=mean(minFDR<=0.05, na.rm=TRUE), .groups="drop")
readr::write_csv(theme_loc, file.path(dirs$tables_gsea, "localization_table.csv"))

theme_loc_by_class <- theme_comp %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
  dplyr::left_join(comp_class_map %>% dplyr::select(context, comparison_class), by="context") %>%
  dplyr::filter(!is.na(region), !is.na(layer), !is.na(comparison_class)) %>%
  dplyr::group_by(db, theme, comparison_class, region, layer) %>%
  dplyr::summarize(meanNES = mean(meanNES, na.rm=TRUE), n_comp = dplyr::n_distinct(context), .groups="drop")
readr::write_csv(theme_loc_by_class, file.path(dirs$tables_class, "localization_by_comparison_class_table.csv"))

class_labels_df <- theme_comp %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
  dplyr::left_join(comp_class_map, by="context") %>%
  dplyr::filter(!is.na(region), !is.na(layer), !is.na(comparison_class)) %>%
  dplyr::group_by(db, theme, region, layer, comparison_class) %>%
  dplyr::summarize(n_comp = dplyr::n_distinct(context), .groups="drop") %>%
  dplyr::mutate(code = dplyr::case_when(
    comparison_class=="sus vs res" ~ "SR",
    comparison_class=="res vs con" ~ "RC",
    comparison_class=="sus vs con" ~ "SC",
    TRUE ~ ""
  ),
  label = paste0(code, ":", n_comp))
theme_loc_by_class_lab <- theme_loc_by_class %>%
  dplyr::left_join(class_labels_df %>% dplyr::select(db, theme, region, layer, comparison_class, label),
                   by=c("db","theme","region","layer","comparison_class"))
readr::write_csv(theme_loc_by_class_lab, file.path(dirs$tables_class, "localization_by_comparison_class_labeled.csv"))

# Condition-centric tables
theme_loc_left_cond <- theme_comp %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L, condition_L), by="context") %>%
  dplyr::filter(!is.na(region), !is.na(layer), !is.na(condition_L)) %>%
  dplyr::group_by(db, theme, region, layer, condition_L) %>%
  dplyr::summarize(meanNES = mean(meanNES, na.rm=TRUE), n_comp = dplyr::n_distinct(context), .groups="drop") %>%
  dplyr::rename(condition = condition_L)
readr::write_csv(theme_loc_left_cond, file.path(dirs$tables_cond, "localization_left_condition_table.csv"))

theme_loc_right_cond <- theme_comp %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_R, layer=layer_R, condition_R), by="context") %>%
  dplyr::filter(!is.na(region), !is.na(layer), !is.na(condition_R)) %>%
  dplyr::group_by(db, theme, region, layer, condition_R) %>%
  dplyr::summarize(meanNES = mean(meanNES, na.rm=TRUE), n_comp = dplyr::n_distinct(context), .groups="drop") %>%
  dplyr::rename(condition = condition_R)
readr::write_csv(theme_loc_right_cond, file.path(dirs$tables_cond, "localization_right_condition_table.csv"))

# Comparison-level long
comp_level <- theme_comp %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
  dplyr::left_join(comp_class_map %>% dplyr::select(context, comparison_class), by="context") %>%
  dplyr::filter(!is.na(region), !is.na(layer), !is.na(comparison_class)) %>%
  dplyr::mutate(ctx_label = paste0(comparison_class, "::", context))
readr::write_csv(comp_level, file.path(dirs$tables_class, "comparison_level_long.csv"))

# 11b) Directionality and bias summaries (tables)
dom_dir_tbl <- theme_comp %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
  dplyr::filter(!is.na(region), !is.na(layer), minFDR<=0.05) %>%
  dplyr::group_by(db, theme, region, layer) %>%
  dplyr::summarize(medNES = stats::median(meanNES, na.rm=TRUE),
                   n_sig  = dplyr::n(), .groups="drop") %>%
  dplyr::mutate(dir = dplyr::case_when(medNES > 0 ~ "up", medNES < 0 ~ "down", TRUE ~ "flat"))
readr::write_csv(dom_dir_tbl, file.path(dirs$tables_gsea, "dominant_direction_table.csv"))

symmetry_test <- theme_comp %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
  dplyr::filter(!is.na(region), !is.na(layer)) %>%
  dplyr::group_by(db, region, layer) %>%
  dplyr::summarize(p_sym = tryCatch({
      x <- meanNES; x <- x[is.finite(x)]
      if (length(x) >= 6) stats::wilcox.test(x, mu=0, exact=FALSE)$p.value else NA_real_
    }, error=function(e) NA_real_), .groups="drop")
readr::write_csv(symmetry_test, file.path(dirs$tables_qa, "directional_bias_pvalues.csv"))

# 11c) Susceptible drivers across classes (safe column names)
sr_sc_strength <- theme_comp %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region = region_L, layer = layer_L), by = "context") %>%
  dplyr::left_join(comp_class_map %>% dplyr::select(context, comparison_class), by = "context") %>%
  dplyr::filter(!is.na(region), !is.na(layer), !is.na(comparison_class)) %>%
  dplyr::group_by(db, theme, region, layer, comparison_class) %>%
  dplyr::summarize(meanNES_class = mean(meanNES, na.rm = TRUE),
                   n_ctx = dplyr::n_distinct(context), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = comparison_class,
    values_from = c(meanNES_class, n_ctx),
    values_fill = 0
  ) %>%
  # Rename wide columns with spaces to safe snake_case
  dplyr::rename(
    meanNES_class_sus_vs_res = `meanNES_class_sus vs res`,
    meanNES_class_sus_vs_con = `meanNES_class_sus vs con`,
    n_ctx_sus_vs_res         = `n_ctx_sus vs res`,
    n_ctx_sus_vs_con         = `n_ctx_sus vs con`
  ) %>%
  dplyr::mutate(
    SR_strength = meanNES_class_sus_vs_res,
    SC_strength = meanNES_class_sus_vs_con,
    driver_class = dplyr::case_when(
      abs(SR_strength) > abs(SC_strength) ~ "sus vs res",
      abs(SC_strength) > abs(SR_strength) ~ "sus vs con",
      TRUE ~ "tie"
    ),
    driver_gap = abs(SR_strength) - abs(SC_strength)
  )
readr::write_csv(sr_sc_strength, file.path(dirs$tables_class, "susceptibility_driver_strength.csv"))

# 11d) Plots: theme localization
# Driver map: which class (SR vs SC) dominates per theme and region-layer
plot_driver_map <- function(dbn) {
  df <- sr_sc_strength %>%
    dplyr::filter(db == dbn, theme %in% (top_themes %>% dplyr::filter(db==dbn) %>% dplyr::pull(theme)))
  if (!nrow(df)) return(NULL)

  # Compact label: nSR/nSC contributing contexts (uses safe column names)
  df <- df %>%
    dplyr::mutate(n_lab = paste0(n_ctx_sus_vs_res, "/", n_ctx_sus_vs_con))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = layer, y = region, fill = driver_class, alpha = abs(driver_gap))) +
    ggplot2::geom_tile(color="white") +
    ggplot2::geom_text(ggplot2::aes(label = n_lab), size=3) +
    ggplot2::scale_fill_manual(values = c("sus vs res" = "#377eb8", "sus vs con" = "#e41a1c", "tie" = "grey80")) +
    ggplot2::scale_alpha(range = c(0.3, 1), guide = "none") +
    ggplot2::facet_wrap(~ theme, ncol = 4) +
    ggplot2::labs(title = paste0(dbn, " — drivers of susceptibility (SR vs SC)"),
                  x = "Layer", y = "Region", fill = "Class") +
    ggplot2::theme_minimal(10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))

  fp <- file.path(dirs$plots_class, paste0("drivers_susceptibility_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=22, height=12); message("[drivers-map] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(sr_sc_strength$db), plot_driver_map))

# 11e) Plots: theme driver volcano
theme_driver_volcano <- sr_sc_strength %>%
  dplyr::group_by(db, theme) %>%
  dplyr::summarize(delta = mean(SR_strength - SC_strength, na.rm=TRUE),
                   abs_gap = mean(abs(SR_strength) - abs(SC_strength), na.rm=TRUE),
                   dom_class = dplyr::case_when(
                     abs_gap > 0 ~ "sus vs res",
                     abs_gap < 0 ~ "sus vs con",
                     TRUE ~ "tie"
                   ),
                   .groups="drop") %>%
  dplyr::left_join(theme_repl %>% dplyr::select(db, theme, repl_score, dom_dir), by=c("db","theme"))

plot_theme_driver_volcano <- function(dbn) {
  df <- theme_driver_volcano %>% dplyr::filter(db==dbn)
  if (!nrow(df)) return(NULL)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = delta, y = repl_score, color = dom_class, shape = dom_dir)) +
    ggplot2::geom_hline(yintercept = median(df$repl_score, na.rm=TRUE), linetype="dashed", color="grey60") +
    ggplot2::geom_vline(xintercept = 0, linetype="dotted", color="grey60") +
    ggplot2::geom_point(size=3, alpha=0.9) +
    ggplot2::scale_color_manual(values=c("sus vs res"="#377eb8","sus vs con"="#e41a1c","tie"="grey60")) +
    ggplot2::labs(title=paste0(dbn, " — theme-level SR vs SC dominance"),
                  x="SR_strength − SC_strength (mean NES)", y="Replication score",
                  color="Dominant class", shape="Dir") +
    ggplot2::theme_minimal(11)
  fp <- file.path(dirs$plots_class, paste0("theme_driver_volcano_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=10, height=7); message("[driver-volcano] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(theme_driver_volcano$db), plot_theme_driver_volcano))

# --- MODULE ANALYSES (parallel to global; place after global sections 11–13 or right after 10b) ---
if (exists("module_gsea_df") && nrow(module_gsea_df) && nrow(mod_rl_both)) {
  # Aggregation per module x comparison (region-layer on left side)
  mod_theme_comp <- module_gsea_df %>%
    dplyr::left_join(mod_rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
    dplyr::filter(!is.na(region), !is.na(layer)) %>%
    dplyr::group_by(db, module, context, region, layer) %>%
    dplyr::summarize(meanNES = mean(NES, na.rm=TRUE),
                     minFDR  = suppressWarnings(min(FDR, na.rm=TRUE)),
                     .groups="drop")
  readr::write_csv(mod_theme_comp, file.path(dirs$tables_modules, "module_theme_by_comparison.csv"))

  # Localization table per module
  module_loc <- mod_theme_comp %>%
    dplyr::group_by(db, module, region, layer) %>%
    dplyr::summarize(meanNES = mean(meanNES, na.rm=TRUE),
                     frac_sig = mean(minFDR<=0.05, na.rm=TRUE),
                     .groups="drop")
  readr::write_csv(module_loc, file.path(dirs$tables_modules, "module_localization_table.csv"))

  # Replication per module
  module_repl <- mod_theme_comp %>%
    dplyr::mutate(sign=sign(meanNES), sig=minFDR<=0.05) %>%
    dplyr::group_by(db, module) %>%
    dplyr::summarize(n_comp=dplyr::n(),
                     n_sig=sum(sig, na.rm=TRUE),
                     n_pos=sum(sig & sign>0, na.rm=TRUE),
                     n_neg=sum(sig & sign<0, na.rm=TRUE),
                     dom_dir=ifelse(n_pos>=n_neg,"up","down"),
                     repl_score=pmax(n_pos,n_neg),
                     .groups="drop") %>%
    dplyr::arrange(dplyr::desc(repl_score))
  readr::write_csv(module_repl, file.path(dirs$tables_modules, "module_replication_scores.csv"))

  # Comparison class labeling
  module_loc_by_class <- mod_theme_comp %>%
    dplyr::left_join(mod_comp_class_map %>% dplyr::select(context, comparison_class), by="context") %>%
    dplyr::filter(!is.na(comparison_class)) %>%
    dplyr::group_by(db, module, comparison_class, region, layer) %>%
    dplyr::summarize(meanNES = mean(meanNES, na.rm=TRUE),
                     n_comp  = dplyr::n_distinct(context),
                     .groups="drop")
  readr::write_csv(module_loc_by_class, file.path(dirs$tables_modules, "module_localization_by_comparison_class.csv"))

  module_class_labels <- mod_theme_comp %>%
    dplyr::left_join(mod_comp_class_map, by="context") %>%
    dplyr::filter(!is.na(comparison_class)) %>%
    dplyr::group_by(db, module, region, layer, comparison_class) %>%
    dplyr::summarize(n_comp = dplyr::n_distinct(context), .groups="drop") %>%
    dplyr::mutate(code = dplyr::case_when(
      comparison_class=="sus vs res" ~ "SR",
      comparison_class=="res vs con" ~ "RC",
      comparison_class=="sus vs con" ~ "SC",
      TRUE ~ ""
    ),
    label = paste0(code, ":", n_comp))
  module_loc_by_class_lab <- module_loc_by_class %>%
    dplyr::left_join(module_class_labels %>% dplyr::select(db,module,region,layer,comparison_class,label),
                     by=c("db","module","region","layer","comparison_class"))
  readr::write_csv(module_loc_by_class_lab, file.path(dirs$tables_modules, "module_localization_by_comparison_class_labeled.csv"))

  # Driver strengths SR vs SC
  mod_sr_sc <- mod_theme_comp %>%
    dplyr::left_join(mod_comp_class_map %>% dplyr::select(context, comparison_class), by="context") %>%
    dplyr::filter(!is.na(comparison_class)) %>%
    dplyr::group_by(db, module, region, layer, comparison_class) %>%
    dplyr::summarize(meanNES_class = mean(meanNES, na.rm=TRUE),
                     n_ctx = dplyr::n_distinct(context), .groups="drop") %>%
    tidyr::pivot_wider(names_from = comparison_class, values_from = c(meanNES_class, n_ctx), values_fill = 0) %>%
    dplyr::rename(
      meanNES_class_sus_vs_res = `meanNES_class_sus vs res`,
      meanNES_class_sus_vs_con = `meanNES_class_sus vs con`,
      n_ctx_sus_vs_res         = `n_ctx_sus vs res`,
      n_ctx_sus_vs_con         = `n_ctx_sus vs con`
    ) %>%
    dplyr::mutate(
      SR_strength = meanNES_class_sus_vs_res,
      SC_strength = meanNES_class_sus_vs_con,
      driver_class = dplyr::case_when(
        abs(SR_strength) > abs(SC_strength) ~ "sus vs res",
        abs(SC_strength) > abs(SR_strength) ~ "sus vs con",
        TRUE ~ "tie"
      ),
      driver_gap = abs(SR_strength) - abs(SC_strength)
    )
  readr::write_csv(mod_sr_sc, file.path(dirs$tables_modules, "module_driver_strength.csv"))

  # Comparison-level long with class label
  mod_comp_level <- mod_theme_comp %>%
    dplyr::left_join(mod_comp_class_map %>% dplyr::select(context, comparison_class), by="context") %>%
    dplyr::filter(!is.na(comparison_class)) %>%
    dplyr::mutate(ctx_label = paste0(comparison_class, "::", context))
  readr::write_csv(mod_comp_level, file.path(dirs$tables_modules, "module_comparison_level_long.csv"))

  # Direction summaries
  mod_dom_dir <- mod_theme_comp %>%
    dplyr::left_join(mod_rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
    dplyr::filter(!is.na(region), !is.na(layer), minFDR<=0.05) %>%
    dplyr::group_by(db, module, region, layer) %>%
    dplyr::summarize(medNES = stats::median(meanNES, na.rm=TRUE),
                     n_sig  = dplyr::n(), .groups="drop") %>%
    dplyr::mutate(dir = dplyr::case_when(medNES > 0 ~ "up", medNES < 0 ~ "down", TRUE ~ "flat"))
  readr::write_csv(mod_dom_dir, file.path(dirs$tables_modules, "module_dominant_direction.csv"))

  mod_symmetry <- mod_theme_comp %>%
    dplyr::left_join(mod_rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
    dplyr::filter(!is.na(region), !is.na(layer)) %>%
    dplyr::group_by(db, region, layer) %>%
    dplyr::summarize(p_sym = tryCatch({
        x <- meanNES; x <- x[is.finite(x)]
        if (length(x) >= 6) stats::wilcox.test(x, mu=0, exact=FALSE)$p.value else NA_real_
      }, error=function(e) NA_real_), .groups="drop")
  readr::write_csv(mod_symmetry, file.path(dirs$tables_modules, "module_directional_bias_pvalues.csv"))

  # Plots: module localization heatmaps for top modules
  plot_module_localization <- function(dbn, n_mod=12) {
    tops <- module_repl %>% dplyr::filter(db==dbn) %>% dplyr::slice_max(repl_score, n=n_mod) %>% dplyr::pull(module)
    df <- module_loc %>% dplyr::filter(db==dbn, module %in% tops)
    if (!nrow(df)) return(NULL)
    p <- ggplot2::ggplot(df, ggplot2::aes(x=layer,y=region,fill=meanNES)) +
      ggplot2::geom_tile(color="grey90") +
      ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b") +
      ggplot2::facet_wrap(~ module, ncol=4) +
      ggplot2::labs(title=paste0("Module localization — ", dbn), x="Layer", y="Region") +
      ggplot2::theme_minimal(10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))
    fp <- file.path(dirs$plots_modules, paste0("module_localization_", dbn, ".svg"))
    ggplot2::ggsave(fp, p, width=22, height=14); message("[module-loc] Saved: ", fp)
    invisible(p)
  }
  invisible(lapply(unique(module_loc$db), plot_module_localization))

  # Plots: class-labeled module localization
plot_module_localization_by_class_labeled <- function(dbn) {
  df_all <- module_loc_by_class_lab %>% dplyr::filter(db==dbn)
  tops <- module_repl %>% dplyr::filter(db==dbn) %>% dplyr::slice_max(repl_score, n=12) %>% dplyr::pull(module)
  if (!length(tops) || !any(df_all$module %in% tops)) {
    tops <- df_all %>% dplyr::count(module, sort=TRUE) %>% dplyr::slice_head(n=12) %>% dplyr::pull(module)
  }
  df <- df_all %>% dplyr::filter(module %in% tops)
  if (!nrow(df)) { message("[module-class-labeled] No data for ", dbn, " after fallback."); return(NULL) }

  class_levels <- c("sus vs res","res vs con","sus vs con")
  df$comparison_class <- factor(df$comparison_class, levels = class_levels)

  p <- ggplot2::ggplot(df, ggplot2::aes(x=layer, y=region, fill=meanNES)) +
    ggplot2::geom_tile(color="grey90") +
    ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b", midpoint=0) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(module),
      cols = ggplot2::vars(comparison_class),
      scales = "free",
      drop = FALSE,
      labeller = ggplot2::labeller(module = ggplot2::label_wrap_gen(width = 18))
    ) +
    ggplot2::labs(title=paste0("Module localization — comparison class (labels=SR/RC/SC:n) — ", dbn),
                  x="Layer", y="Region") +
    ggplot2::theme_minimal(base_size=10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1)) +
    ggplot2::geom_text(ggplot2::aes(label = dplyr::coalesce(label, as.character(n_comp))),
                       color="black", size=3, na.rm=TRUE)

  fp <- file.path(dirs$plots_modules, paste0("module_localization_", dbn, "_by_comparison_class_labeled.svg"))
  ggplot2::ggsave(fp, p, width=30, height=20); message("[module-class-labeled] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(module_loc_by_class_lab$db), plot_module_localization_by_class_labeled))

# Plots: module comparison matrices (by module)
plot_module_by_comparisons <- function(dbn, max_cols=60) {
  df_db <- mod_comp_level %>% dplyr::filter(db==dbn)
  top_mods <- module_repl %>% dplyr::filter(db==dbn) %>% dplyr::slice_max(repl_score, n=8) %>% dplyr::pull(module)
  if (!length(top_mods)) { message("[mod-comp-matrix] No top modules for ", dbn); return(NULL) }
  lapply(top_mods, function(md) {
    df <- df_db %>% dplyr::filter(module==md)
    ctx_ord <- df %>% dplyr::count(ctx_label, sort=TRUE) %>% dplyr::slice_head(n=max_cols) %>% dplyr::pull(ctx_label)
    df <- df %>% dplyr::filter(ctx_label %in% ctx_ord)
    df$ctx_label <- factor(df$ctx_label, levels=ctx_ord)
    p <- ggplot2::ggplot(df, ggplot2::aes(x=ctx_label, y=paste(region, layer, sep="_"), fill=meanNES)) +
      ggplot2::geom_tile(color="grey90") +
      ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b", midpoint=0) +
      ggplot2::labs(title=paste0("Module ", md, " — comparison-level NES — ", dbn),
                    x="Comparison (class::context)", y="Region_Layer") +
      ggplot2::theme_minimal(base_size=10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle=60, hjust=1, vjust=1))
    fp <- file.path(dirs$plots_modules, paste0("module_", gsub("[^A-Za-z0-9]+","_", md), "_", dbn, "_comparison_matrix.svg"))
    ggplot2::ggsave(fp, p, width=28, height=10); message("[mod-comp-matrix] Saved: ", fp)
    invisible(p)
  })
}
invisible(lapply(unique(mod_comp_level$db), plot_module_by_comparisons))

# Plots: module driver map (SR vs SC)
plot_module_driver_map <- function(dbn) {
  top_mods <- module_repl %>% dplyr::filter(db==dbn) %>% dplyr::slice_max(repl_score, n=12) %>% dplyr::pull(module)
  df <- mod_sr_sc %>% dplyr::filter(db==dbn, module %in% top_mods)
  if (!nrow(df)) return(NULL)
  df <- df %>% dplyr::mutate(n_lab = paste0(n_ctx_sus_vs_res, "/", n_ctx_sus_vs_con))
  p <- ggplot2::ggplot(df, ggplot2::aes(x=layer, y=region, fill=driver_class, alpha=abs(driver_gap))) +
    ggplot2::geom_tile(color="white") +
    ggplot2::geom_text(ggplot2::aes(label = n_lab), size=3) +
    ggplot2::scale_fill_manual(values = c("sus vs res" = "#377eb8", "sus vs con" = "#e41a1c", "tie" = "grey80")) +
    ggplot2::scale_alpha(range = c(0.3, 1), guide = "none") +
    ggplot2::facet_wrap(~ module, ncol = 4, labeller = ggplot2::labeller(module = ggplot2::label_wrap_gen(18))) +
    ggplot2::labs(title = paste0(dbn, " — module drivers (SR vs SC)"),
                  x = "Layer", y = "Region", fill = "Class") +
    ggplot2::theme_minimal(10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))
  fp <- file.path(dirs$plots_modules, paste0("module_drivers_susceptibility_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=22, height=12); message("[mod-drivers] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(mod_sr_sc$db), plot_module_driver_map))

# 12) NA diagnostics -----------------------------------------------------------
cov_tbl <- assign_themes %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by = "context") %>%
  dplyr::filter(!is.na(region), !is.na(layer)) %>%
  dplyr::group_by(db, theme, region, layer) %>%
  dplyr::summarize(n_comp_cover=dplyr::n_distinct(comp_dir), n_sig_cover=sum(FDR<=0.05, na.rm=TRUE), minFDR_all=suppressWarnings(min(FDR, na.rm=TRUE)), .groups="drop")
na_diag <- theme_loc %>%
  dplyr::full_join(cov_tbl, by = c("db","theme","region","layer")) %>%
  dplyr::mutate(cause = dplyr::case_when(
    is.na(meanNES) & (is.na(n_comp_cover) | n_comp_cover == 0) ~ "no_coverage",
    is.na(meanNES) & (n_comp_cover > 0) & (is.na(n_sig_cover) | n_sig_cover == 0) ~ "masked_nonsignificant",
    is.na(meanNES) ~ "other_missing",
    TRUE ~ "has_value"
  ))
readr::write_csv(na_diag, file.path(dirs$tables_qa, "localization_NA_diagnostics.csv"))

# 13) Replication summary ------------------------------------------------------
theme_repl <- theme_comp %>% dplyr::mutate(sign=sign(meanNES), sig=minFDR<=0.05) %>%
  dplyr::group_by(db, theme) %>% dplyr::summarize(n_comp=dplyr::n(), n_sig=sum(sig, na.rm=TRUE), n_pos=sum(sig & sign>0, na.rm=TRUE), n_neg=sum(sig & sign<0, na.rm=TRUE), dom_dir=ifelse(n_pos>=n_neg,"up","down"), repl_score=pmax(n_pos,n_neg), .groups="drop") %>%
  dplyr::arrange(dplyr::desc(repl_score))
readr::write_csv(theme_repl, file.path(dirs$tables_gsea, "theme_replication_scores.csv"))

top_themes <- theme_repl %>% dplyr::group_by(db) %>% dplyr::slice_max(order_by=repl_score, n=12) %>% dplyr::ungroup()

# 14) Plots -------------------------------------------------------------------
plot_theme_localization <- function(dbn) {
  df <- theme_loc %>% dplyr::filter(db==dbn, theme %in% (top_themes %>% dplyr::filter(db==dbn) %>% dplyr::pull(theme)))
  if (!nrow(df)) { message("[anatomy] No data for ", dbn); return(NULL) }
  p <- ggplot2::ggplot(df, ggplot2::aes(x=layer, y=region, fill=meanNES)) +
    ggplot2::geom_tile(color="grey90") +
    ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b") +
    ggplot2::facet_wrap(~ theme, scales="free", ncol=4) +
    ggplot2::labs(title=paste0("Localization by region-layer — ", dbn), x="Layer", y="Region") +
    ggplot2::theme_minimal(base_size=10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))
  fp <- file.path(dirs$plots_anat, paste0("localization_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=22, height=14); message("[anatomy] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(theme_loc$db), plot_theme_localization))

plot_theme_localization_by_class_labeled <- function(dbn) {
  df_all <- theme_loc_by_class_lab %>% dplyr::filter(db==dbn)
  tops <- top_themes %>% dplyr::filter(db==dbn) %>% dplyr::pull(theme)
  if (!length(tops) || !any(df_all$theme %in% tops)) {
    tops <- df_all %>% dplyr::count(theme, sort=TRUE) %>% dplyr::slice_head(n=12) %>% dplyr::pull(theme)
  }
  df <- df_all %>% dplyr::filter(theme %in% tops)
  if (!nrow(df)) { message("[class-labeled] No data for ", dbn, " after fallback."); return(NULL) }
  class_levels <- c("sus vs res","res vs con","sus vs con")
  df$comparison_class <- factor(df$comparison_class, levels = class_levels)
  p <- ggplot2::ggplot(df, ggplot2::aes(x=layer, y=region, fill=meanNES)) +
    ggplot2::geom_tile(color="grey90") +
    ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b") +
    ggplot2::facet_grid(rows = ggplot2::vars(theme), cols = ggplot2::vars(comparison_class), scales="free", drop=FALSE) +
    ggplot2::labs(title=paste0("Localization — comparison class (labeled) — ", dbn), x="Layer", y="Region") +
    ggplot2::theme_minimal(base_size=10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1)) +
    ggplot2::geom_text(ggplot2::aes(label = dplyr::coalesce(label, as.character(n_comp))), color="black", size=3, na.rm=TRUE)
  fp <- file.path(dirs$plots_class, paste0("localization_", dbn, "_by_comparison_class_labeled.svg"))
  ggplot2::ggsave(fp, p, width=30, height=20); message("[class-labeled] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(theme_loc_by_class_lab$db), plot_theme_localization_by_class_labeled))

plot_theme_localization_left_cond <- function(dbn) {
  df <- theme_loc_left_cond %>% dplyr::filter(db==dbn, theme %in% (top_themes %>% dplyr::filter(db==dbn) %>% dplyr::pull(theme)))
  if (!nrow(df)) { message("[left-cond] No data for ", dbn); return(NULL) }
  p <- ggplot2::ggplot(df, ggplot2::aes(x=layer, y=region, fill=meanNES)) +
    ggplot2::geom_tile(color="grey90") +
    ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b") +
    ggplot2::facet_grid(rows = ggplot2::vars(theme), cols = ggplot2::vars(condition), scales="free") +
    ggplot2::labs(title=paste0("Localization — tested condition (sus/res) — ", dbn), x="Layer", y="Region") +
    ggplot2::theme_minimal(base_size=10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))
  fp <- file.path(dirs$plots_cond, paste0("localization_", dbn, "_tested_condition.svg"))
  ggplot2::ggsave(fp, p, width=26, height=20); message("[left-cond] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(theme_loc_left_cond$db), plot_theme_localization_left_cond))

plot_theme_localization_right_cond <- function(dbn) {
  df <- theme_loc_right_cond %>% dplyr::filter(db==dbn, theme %in% (top_themes %>% dplyr::filter(db==dbn) %>% dplyr::pull(theme)))
  if (!nrow(df)) { message("[right-cond] No data for ", dbn); return(NULL) }
  p <- ggplot2::ggplot(df, ggplot2::aes(x=layer, y=region, fill=meanNES)) +
    ggplot2::geom_tile(color="grey90") +
    ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b") +
    ggplot2::facet_grid(rows = ggplot2::vars(theme), cols = ggplot2::vars(condition), scales="free") +
    ggplot2::labs(title=paste0("Localization — baseline condition (res/con) — ", dbn), x="Layer", y="Region") +
    ggplot2::theme_minimal(base_size=10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))
  fp <- file.path(dirs$plots_cond, paste0("localization_", dbn, "_baseline_condition.svg"))
  ggplot2::ggsave(fp, p, width=26, height=20); message("[right-cond] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(theme_loc_right_cond$db), plot_theme_localization_right_cond))

plot_theme_by_comparisons <- function(dbn, max_cols=60) {
  df_db <- comp_level %>% dplyr::filter(db==dbn)
  themes_here <- df_db %>% dplyr::count(theme, sort=TRUE) %>% dplyr::slice_head(n=12) %>% dplyr::pull(theme)
  if (!length(themes_here)) { message("[comp-matrix] No themes with class data for ", dbn); return(NULL) }
  lapply(themes_here, function(th) {
    df <- df_db %>% dplyr::filter(theme==th)
    ctx_ord <- df %>% dplyr::count(ctx_label, sort=TRUE) %>% dplyr::slice_head(n=max_cols) %>% dplyr::pull(ctx_label)
    df <- df %>% dplyr::filter(ctx_label %in% ctx_ord)
    df$ctx_label <- factor(df$ctx_label, levels=ctx_ord)
    p <- ggplot2::ggplot(df, ggplot2::aes(x=ctx_label, y=paste(region, layer, sep="_"), fill=meanNES)) +
      ggplot2::geom_tile(color="grey90") +
      ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b") +
      ggplot2::labs(title=paste0("Theme ", th, " — comparison-level NES — ", dbn), x="Comparison (class::context)", y="Region_Layer") +
      ggplot2::theme_minimal(base_size=10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle=60, hjust=1, vjust=1))
    fp <- file.path(dirs$plots_cm, paste0("theme_", gsub("[^A-Za-z0-9]+","_", th), "_", dbn, "_comparison_matrix.svg"))
    ggplot2::ggsave(fp, p, width=28, height=10); message("[comp-matrix] Saved: ", fp)
    invisible(p)
  })
}
invisible(lapply(unique(comp_level$db), plot_theme_by_comparisons))

# 14b) Compact panel (Anatomy + Class + Replication) per DB
plot_compact_panel <- function(dbn, n_themes = 6) {
  the <- top_themes %>% dplyr::filter(db==dbn) %>% dplyr::slice_head(n=n_themes) %>% dplyr::pull(theme)

  p1 <- ggplot2::ggplot(theme_loc %>% dplyr::filter(db==dbn, theme %in% the),
                        ggplot2::aes(x=layer,y=region,fill=meanNES)) +
    ggplot2::geom_tile(color="grey90") +
    ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b") +
    ggplot2::facet_wrap(~ theme, ncol=3) +
    ggplot2::labs(title=paste0(dbn, " — localization (NES)")) + ggplot2::theme_minimal(9)

  p2 <- ggplot2::ggplot(theme_loc_by_class_lab %>% dplyr::filter(db==dbn, theme %in% the),
                        ggplot2::aes(x=layer,y=region,fill=meanNES,label=label)) +
    ggplot2::geom_tile(color="grey95") +
    ggplot2::geom_text(size=2, color="black", na.rm=TRUE) +
    ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b") +
    ggplot2::facet_grid(rows=ggplot2::vars(theme), cols=ggplot2::vars(comparison_class), drop=FALSE) +
    ggplot2::labs(title="comparison class (labels=SR/RC/SC:n)") + ggplot2::theme_minimal(9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45,hjust=1))

  p3 <- ggplot2::ggplot(theme_repl %>% dplyr::filter(db==dbn, theme %in% the),
                        ggplot2::aes(x=reorder(theme, repl_score), y=repl_score, fill=dom_dir)) +
    ggplot2::geom_col() + ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values=c(down="#2166ac", up="#b2182b")) +
    ggplot2::labs(title="replication score", x="", y="max(up, down) significant") +
    ggplot2::theme_minimal(9)

  if (requireNamespace("patchwork", quietly=TRUE)) {
    panel <- p1 / p2 / p3 + patchwork::plot_layout(heights=c(1,1.1,0.6))
    fp <- file.path(dirs$plots, "Panels", paste0("panel_", dbn, "_compact.svg"))
    dir.create(dirname(fp), recursive=TRUE, showWarnings=FALSE)
    ggplot2::ggsave(fp, panel, width=14, height=16); message("[panel] Saved: ", fp)
  } else {
    message("[panel] patchwork not installed; skipping combined panel for ", dbn)
  }
}
invisible(lapply(unique(theme_loc$db), plot_compact_panel))

# 14c) Dominant direction lattice per DB (median NES, n_sig annotation)
plot_dir_lattice <- function(dbn) {
  df <- dom_dir_tbl %>% dplyr::filter(db==dbn, theme %in% (top_themes %>% dplyr::filter(db==dbn) %>% dplyr::pull(theme)))
  if (!nrow(df)) return(NULL)
  p <- ggplot2::ggplot(df, ggplot2::aes(x=layer, y=region, fill=dir, label=n_sig)) +
    ggplot2::geom_tile(color="white") + ggplot2::geom_text(size=3, color="black") +
    ggplot2::scale_fill_manual(values=c(down="#2166ac", flat="grey85", up="#b2182b")) +
    ggplot2::facet_wrap(~ theme, ncol=4) +
    ggplot2::labs(title=paste0(dbn, " — dominant direction (median NES, n_sig)"),
                  x="Layer", y="Region") + ggplot2::theme_minimal(9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))
  fp <- file.path(dirs$plots_anat, paste0("dominant_direction_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=22, height=12); message("[dir-lattice] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(theme_loc$db), plot_dir_lattice))

# 14d) Directional bias tiles per DB (−log10 p from Wilcoxon against 0)
plot_symmetry <- function(dbn) {
  df <- symmetry_test %>% dplyr::filter(db==dbn)
  p <- ggplot2::ggplot(df, ggplot2::aes(x=layer, y=region, fill=-log10(p_sym))) +
    ggplot2::geom_tile(color="white") +
    ggplot2::scale_fill_viridis_c(option="C", na.value="grey90") +
    ggplot2::labs(title=paste0(dbn, " — directional bias (−log10 p)"), x="Layer", y="Region") +
    ggplot2::theme_minimal(9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))
  fp <- file.path(dirs$plots_anat, paste0("directional_bias_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=8, height=6); message("[bias] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(theme_loc$db), plot_symmetry))

# 14e) Comparison-class balance per theme (stacked fractions across region-layer)
class_balance <- theme_comp %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
  dplyr::left_join(comp_class_map, by="context") %>%
  dplyr::filter(!is.na(region), !is.na(layer), !is.na(comparison_class)) %>%
  dplyr::group_by(db, theme, region, layer, comparison_class) %>%
  dplyr::summarize(n_cmp=dplyr::n_distinct(context), .groups="drop") %>%
  dplyr::group_by(db, theme, region, layer) %>%
  dplyr::mutate(frac = n_cmp/sum(n_cmp)) %>% dplyr::ungroup()
readr::write_csv(class_balance, file.path(dirs$tables_class, "comparison_class_balance.csv"))

plot_class_balance <- function(dbn, theme_id) {
  df <- class_balance %>% dplyr::filter(db==dbn, theme==theme_id)
  if (!nrow(df)) return(NULL)
  p <- ggplot2::ggplot(df, ggplot2::aes(x=paste(region,layer,sep="_"), y=frac, fill=comparison_class)) +
    ggplot2::geom_col(width=0.8) +
    ggplot2::scale_fill_manual(values=c("sus vs res"="#377eb8", "res vs con"="#4daf4a", "sus vs con"="#e41a1c")) +
    ggplot2::labs(title=paste0(theme_id, " — comparison balance (", dbn, ")"), x="Region_Layer", y="Fraction of contexts") +
    ggplot2::theme_minimal(9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=60, hjust=1))
  fp <- file.path(dirs$plots_class, paste0("class_balance_", gsub("[^A-Za-z0-9]+","_", theme_id), "_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=12, height=6); message("[class-balance] Saved: ", fp)
  invisible(p)
}
invisible(lapply(split(top_themes, top_themes$db), function(dfdb) {
  dbn <- unique(dfdb$db); ths <- head(dfdb$theme, 8)
  lapply(ths, function(th) plot_class_balance(dbn, th))
}))

# 14f) Contrast strip (compact comparison matrix subset) for top themes
plot_theme_strip <- function(dbn, theme_id, max_ctx=40) {
  df <- comp_level %>% dplyr::filter(db==dbn, theme==theme_id)
  if (!nrow(df)) return(NULL)
  ctx_ord <- df %>% dplyr::count(ctx_label, sort=TRUE) %>% dplyr::slice_head(n=max_ctx) %>% dplyr::pull(ctx_label)
  df <- df %>% dplyr::filter(ctx_label %in% ctx_ord)
  df$ctx_label <- factor(df$ctx_label, levels=ctx_ord)

  p <- ggplot2::ggplot(df, ggplot2::aes(x=ctx_label, y=paste(region,layer,sep="_"), fill=meanNES)) +
    ggplot2::geom_tile(color="grey95") +
    ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b") +
    ggplot2::facet_wrap(~ comparison_class, ncol=1, scales="free_x") +
    ggplot2::labs(title=paste0(theme_id, " — contrast strip (", dbn, ")"), x="Comparison", y="Region_Layer") +
    ggplot2::theme_minimal(9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=60, hjust=1))
  fp <- file.path(dirs$plots_cm, paste0("theme_", gsub("[^A-Za-z0-9]+","_", theme_id), "_", dbn, "_strip.svg"))
  ggplot2::ggsave(fp, p, width=14, height=8); message("[strip] Saved: ", fp)
  invisible(p)
}
invisible(lapply(split(top_themes, top_themes$db), function(dfdb) {
  dbn <- unique(dfdb$db); ths <- head(dfdb$theme, 6)
  lapply(ths, function(th) plot_theme_strip(dbn, th, max_ctx=40))
}))

# 14g) Simple localization graph per DB (themes ↔ region-layer nodes)
t_thr <- 0.6
edges_df <- theme_loc %>%
  dplyr::left_join(
    theme_comp %>%
      dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L),
                       by="context") %>%
      dplyr::filter(!is.na(region), !is.na(layer), minFDR<=0.05) %>%
      dplyr::group_by(db, theme, region, layer) %>%
      dplyr::summarize(n_sig = dplyr::n(), sign_dir = sign(mean(meanNES, na.rm=TRUE)), .groups="drop"),
    by=c("db","theme","region","layer")
  ) %>%
  # guard: keep rows with numeric NES and threshold on absolute value
  dplyr::filter(!is.na(meanNES), is.finite(meanNES), abs(meanNES) >= t_thr) %>%
  # vectorized fallback for missing n_sig
  dplyr::mutate(n_sig = ifelse(is.na(n_sig), 1L, n_sig),
                src = theme,
                dst = paste(region, layer, sep="_"),
                edge_col = ifelse(meanNES > 0, "#b2182b", "#2166ac"),
                edge_w = scales::rescale(pmax(n_sig, 1L), to=c(0.4, 2)))

plot_localization_graph <- function(dbn) {
  df <- edges_df %>% dplyr::filter(db==dbn)
  if (!nrow(df)) return(NULL)

  # Build graph
  g <- igraph::graph_from_data_frame(df %>% dplyr::select(src, dst), directed=FALSE)

  # Align edge attributes to graph edges
  e_df <- igraph::as_data_frame(g, what="edges")
  key_g  <- paste(e_df$from, e_df$to, sep="||")
  key_df <- paste(df$src, df$dst, sep="||")
  idx <- match(key_g, key_df)
  mis <- which(is.na(idx))
  if (length(mis)) {
    key_rev <- paste(e_df$to[mis], e_df$from[mis], sep="||")
    idx_rev <- match(key_rev, key_df)
    idx[mis] <- idx_rev
  }
  idx[is.na(idx)] <- 1L

  igraph::E(g)$edge_col <- df$edge_col[idx]
  igraph::E(g)$edge_w   <- df$edge_w[idx]

  set.seed(1)
  p <- ggraph::ggraph(g, layout="fr") +
    ggraph::geom_edge_link(ggplot2::aes(edge_colour = edge_col,
                                        edge_width  = edge_w), alpha=0.8) +
    ggraph::geom_node_point(ggplot2::aes(shape = ifelse(grepl("_", name), "RL", "Theme")), size=3) +
    ggraph::geom_node_text(ggplot2::aes(label=name), repel=TRUE, size=3) +
    ggraph::scale_edge_colour_identity() +
    ggraph::scale_edge_width(range=c(0.4,2)) +
    ggplot2::scale_shape_manual(values=c(Theme=16, RL=15)) +
    ggplot2::theme_void() +
    ggplot2::labs(title=paste0(dbn, " — localization graph (|NES|≥", t_thr, ")"))

  fp <- file.path(dirs$plots, "Graphs", paste0("localization_graph_", dbn, ".svg"))
  dir.create(dirname(fp), recursive=TRUE, showWarnings=FALSE)
  ggplot2::ggsave(fp, p, width=10, height=8); message("[graph] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(theme_loc$db), plot_localization_graph))


# 14h) Terms network per DB (Jaccard similarity of GSEA terms)
# --- Terms network (per DB) ---
# Recompute similarity within DB to match your themes more closely
# Term network colored by theme cluster
tokenize <- function(x){
  x <- tolower(x); x <- gsub("[^a-z0-9 ]+"," ",x)
  unique(unlist(strsplit(x,"\\s+")))
}
jaccard <- function(a,b){
  ia <- tokenize(a); ib <- tokenize(b)
  if (!length(ia) || !length(ib)) return(0)
  length(intersect(ia,ib)) / length(union(ia,ib))
}

rep_terms <- gsea_df %>%
  dplyr::inner_join(dplyr::select(themes_tbl, ID, theme), by="ID") %>%
  dplyr::group_by(db, theme, ID, Description) %>%
  dplyr::summarize(freq = dplyr::n(), .groups="drop") %>%
  dplyr::group_by(db, theme) %>%
  dplyr::slice_max(order_by=freq, n=10, with_ties=FALSE) %>%
  dplyr::ungroup()

plot_terms_network <- function(dbn, sim_cut=0.25, max_nodes=250) {
  df <- rep_terms %>% dplyr::filter(db==dbn)
  if (!nrow(df)) return(NULL)
  # control size: top 12 themes
  top_themes <- df %>% dplyr::count(theme, sort=TRUE) %>% dplyr::slice_head(n=12) %>% dplyr::pull(theme)
  df <- df %>% dplyr::filter(theme %in% top_themes)
  if (nrow(df) > max_nodes) {
    per <- max(5L, floor(max_nodes/length(top_themes)))
    df <- df %>% dplyr::group_by(theme) %>% dplyr::slice_head(n=per) %>% dplyr::ungroup()
  }

  terms <- df$Description; n <- length(terms)
  if (n < 2) return(NULL)
  edges <- list()
  for (i in seq_len(n-1)) for (j in (i+1):n) {
    s <- jaccard(terms[i], terms[j]); if (s >= sim_cut) edges[[length(edges)+1]] <- c(i,j,s)
  }
  if (!length(edges)) return(NULL)
  edges <- as.data.frame(do.call(rbind, edges)); names(edges) <- c("i","j","sim")
  edges$i <- as.integer(edges$i); edges$j <- as.integer(edges$j); edges$sim <- as.numeric(edges$sim)

  verts <- tibble::tibble(name = terms, theme = df$theme, term = df$Description)
  g <- igraph::graph_from_data_frame(
    d = tibble::tibble(from = verts$name[edges$i], to = verts$name[edges$j], sim = edges$sim),
    directed = FALSE, vertices = verts
  )

  set.seed(1)
  p <- ggraph::ggraph(g, layout="fr") +
    ggraph::geom_edge_link(ggplot2::aes(width=sim), colour="grey70", alpha=0.35) +
    ggraph::geom_node_point(ggplot2::aes(color=theme), size=2) +
    ggraph::geom_node_text(ggplot2::aes(label=term), size=2.4, repel=TRUE) +
    ggplot2::scale_edge_width(range=c(0.2,1.4), guide="none") +
    ggplot2::theme_void() +
    ggplot2::labs(title=paste0(dbn, " — terms network (sim≥", sim_cut, ")"), color="Theme")
  fp <- file.path(dirs$plots_anat, paste0("terms_network_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=16, height=12); message("[terms-net] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(rep_terms$db), plot_terms_network))

suppressWarnings({
  if (!requireNamespace("uwot", quietly=TRUE)) message("Optional: install.packages('uwot') for UMAP")
})

build_theme_kw_matrix <- function(dbn, top_kw=800) {
  df <- gsea_df %>%
    dplyr::filter(db==dbn) %>%
    dplyr::inner_join(dplyr::select(themes_tbl, ID, theme), by="ID") %>%
    dplyr::distinct(theme, Description)
  if (!nrow(df)) return(NULL)

  toks <- strsplit(tolower(gsub("[^a-z0-9 ]+"," ", df$Description)), "\\s+")
  toks <- lapply(toks, function(v) v[nzchar(v) & !v %in% c("and","or","of","the","to","in","by","for","process","regulation","cellular","activity","protein","pathway","signaling")])

  env <- new.env(parent=emptyenv())
  for (i in seq_len(nrow(df))) {
    th <- df$theme[i]
    for (k in toks[[i]]) {
      key <- paste(th,k,sep="||"); env[[key]] <- (env[[key]] %||% 0L) + 1L
    }
  }
  keys <- ls(env); if (!length(keys)) return(NULL)
  sp <- strsplit(keys,"\\|\\|"); ths <- vapply(sp, `[`, character(1), 1); kws <- vapply(sp, `[`, character(1), 2)
  val <- as.integer(mget(keys, env, ifnotfound=0L))
  mat <- tibble::tibble(theme=ths, kw=kws, val=val) %>%
    dplyr::group_by(kw) %>% dplyr::summarize(dfreq=dplyr::n(), .groups="drop") %>%
    dplyr::right_join(tibble::tibble(theme=ths, kw=kws, val=val), by="kw") %>%
    dplyr::filter(dfreq>=2) %>%
    dplyr::group_by(kw) %>% dplyr::mutate(tf = val/sum(val)) %>% dplyr::ungroup() %>%
    dplyr::group_by(kw) %>% dplyr::mutate(idf = log(n_distinct(theme)/n_distinct(theme[val>0])+1e-6)) %>% dplyr::ungroup() %>%
    dplyr::mutate(tfidf = tf*idf)

  vocab <- mat %>% dplyr::group_by(kw) %>% dplyr::summarize(s=sum(tfidf,na.rm=TRUE), .groups="drop") %>%
    dplyr::arrange(dplyr::desc(s)) %>% dplyr::slice_head(n=top_kw) %>% dplyr::pull(kw)
  mat <- mat %>% dplyr::filter(kw %in% vocab)
  if (!nrow(mat)) return(NULL)

  M <- tidyr::pivot_wider(mat, names_from=kw, values_from=tfidf, values_fill=0)
  M
}

plot_theme_embedding <- function(dbn) {
  M <- build_theme_kw_matrix(dbn)
  if (is.null(M) || nrow(M) < 2) return(NULL)
  meta <- M %>% dplyr::select(theme)
  X <- as.matrix(M %>% dplyr::select(-theme))
  rownames(X) <- meta$theme

  if (requireNamespace("uwot", quietly=TRUE)) {
    set.seed(1); emb <- uwot::umap(X, n_neighbors=10, min_dist=0.2, metric="cosine")
    emb <- as.data.frame(emb); colnames(emb) <- c("U1","U2")
  } else {
    pr <- stats::prcomp(X, scale.=TRUE); emb <- as.data.frame(pr$x[,1:2]); colnames(emb) <- c("U1","U2")
  }
  emb$theme <- rownames(X)

  p <- ggplot2::ggplot(emb, ggplot2::aes(x=U1, y=U2, label=theme)) +
    ggplot2::geom_point(color="#2c7fb8", size=2, alpha=0.85) +
    ggplot2::geom_text(size=2.8, nudge_y=0.02) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title=paste0(dbn, " — theme embedding (", if (requireNamespace("uwot", quietly=TRUE)) "UMAP" else "PCA", ")"),
                  x="Dim 1", y="Dim 2")
  fp <- file.path(dirs$plots_anat, paste0("theme_embedding_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=10, height=8); message("[theme-embed] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(themes_tbl$db), plot_theme_embedding))

plot_theme_dendrogram <- function(dbn) {
  M <- build_theme_kw_matrix(dbn)
  if (is.null(M) || nrow(M) < 2) return(NULL)
  X <- as.matrix(M %>% dplyr::select(-theme))
  rownames(X) <- M$theme
  # cosine distance
  A <- X + 1e-12
  nrm <- sqrt(rowSums(A*A))
  S <- tcrossprod(A / nrm)
  D <- 1 - S
  hc <- stats::hclust(stats::as.dist(D), method="average")

  if (!requireNamespace("ggdendro", quietly=TRUE)) {
    message("Optional: install.packages('ggdendro') for dendrogram plot"); return(NULL)
  }
  dd <- ggdendro::dendro_data(hc)
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data=dd$segments, ggplot2::aes(x=x, y=y, xend=xend, yend=yend)) +
    ggplot2::geom_text(data=dd$labels, ggplot2::aes(x=x, y=y, label=label), angle=90, hjust=1, vjust=0.5, size=2.6) +
    ggplot2::theme_void() + ggplot2::labs(title=paste0(dbn, " — theme dendrogram"))
  fp <- file.path(dirs$plots_anat, paste0("theme_dendrogram_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=14, height=10); message("[theme-dend] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(themes_tbl$db), plot_theme_dendrogram))

# 14i) GO BP terms per theme (bar, top by FDR/frequency)
plot_bp_terms_per_theme <- function(top_k = 8) {
  df <- gsea_df %>%
    dplyr::filter(db == "GO_BP") %>%
    dplyr::inner_join(dplyr::select(themes_tbl, ID, theme), by="ID") %>%
    dplyr::group_by(theme, Description) %>%
    dplyr::summarize(n_hits = dplyr::n(), meanNES = mean(NES, na.rm=TRUE),
                     bestFDR = suppressWarnings(min(FDR, na.rm=TRUE)), .groups="drop") %>%
    dplyr::group_by(theme) %>%
    dplyr::arrange(bestFDR, dplyr::desc(n_hits), .by_group=TRUE) %>%
    dplyr::slice_head(n = top_k) %>%
    dplyr::ungroup()

  if (!nrow(df)) { message("[bp-terms] No GO_BP data"); return(NULL) }

  # Reorder within each facet
  df <- df %>%
    dplyr::group_by(theme) %>%
    dplyr::mutate(term_order = reorder(Description, -n_hits)) %>%
    dplyr::ungroup()

  p <- ggplot2::ggplot(df, ggplot2::aes(x = n_hits, y = term_order, fill = -log10(bestFDR))) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_viridis_c(option="C", name = "-log10(FDR)") +
    ggplot2::facet_wrap(~ theme, scales = "free_y", ncol = 3) +
    ggplot2::labs(title = "GO BP terms per cluster (top by FDR/frequency)",
                  x = "Term count within cluster", y = "GO BP term") +
    ggplot2::theme_minimal(10)
  fp <- file.path(dirs$plots_anat, "GO_BP_terms_per_theme.svg")
  ggplot2::ggsave(fp, p, width=18, height=14); message("[bp-terms] Saved: ", fp)
  invisible(p)
}
plot_bp_terms_per_theme(top_k = 8)

bp_theme_summaries <- bp_members %>%
  dplyr::group_by(theme) %>%
  dplyr::arrange(bestFDR, dplyr::desc(n_hits), .by_group=TRUE) %>%
  dplyr::summarize(
    n_terms = dplyr::n(),
    top_terms = paste(head(Description, 6), collapse=" | "),
    medianNES = stats::median(meanNES, na.rm=TRUE),
    .groups="drop"
  )
readr::write_csv(bp_theme_summaries, file.path(dirs$tables_gsea, "GO_BP_theme_summaries.csv"))

# 15) ORA summaries (separate) ------------------------------------------------
if (nrow(ora_df)) {
  ora_repl <- ora_df %>% dplyr::group_by(db, ID, Description) %>% dplyr::summarize(n_comp=dplyr::n(), n_sig=sum(Padj<=0.05, na.rm=TRUE), meanPadj=mean(Padj, na.rm=TRUE), .groups="drop") %>% dplyr::arrange(dplyr::desc(n_sig), meanPadj)
  readr::write_csv(ora_repl, file.path(dirs$tables_ora, "ORA_replication_scores.csv"))
  top_ora <- ora_repl %>% dplyr::group_by(db) %>% dplyr::slice_max(n_sig, n=25) %>% dplyr::ungroup()
  p_ora <- ggplot2::ggplot(top_ora, ggplot2::aes(x=reorder(paste(db, Description, sep="::"), n_sig), y=n_sig, fill=db)) +
    ggplot2::geom_col() + ggplot2::coord_flip() +
    ggplot2::labs(title="ORA term replication (count of significant comparisons)", x="DB::Term", y="# significant comps") +
    ggplot2::theme_minimal(base_size=10)
  ggplot2::ggsave(file.path(dirs$plots_ora, "ORA_replication_bar.svg"), p_ora, width=18, height=12)
}

# 16) Executive bundle + artifact check ---------------------------------------
readr::write_csv(top_themes,  file.path(dirs$tables_gsea, "top_themes_per_db.csv"))
readr::write_csv(theme_comp %>% dplyr::arrange(db, theme, dplyr::desc(abs(meanNES))), file.path(dirs$tables_gsea, "theme_by_comparison_sorted.csv"))

md <- c(
  "# Consolidated Summary",
  paste0("- Comparisons scanned: ", length(comp_dirs)),
  paste0("- GSEA rows: ", nrow(gsea_df)),
  "",
  "Artifacts:",
  "* Plots/Anatomy/localization_<DB>.svg",
  "* Plots/Class/localization_<DB>_by_comparison_class_labeled.svg",
  "* Plots/Conditions/localization_<DB>_tested_condition.svg",
  "* Plots/Conditions/localization_<DB>_baseline_condition.svg",
  "* Plots/ComparisonMatrices/theme_<THEME>_<DB>_comparison_matrix.svg",
  "* Tables/context_region_layer_condition_both_sides.csv",
  "* Tables/Class/context_comparison_class_map.csv",
  "* Tables/GSEA/localization_table.csv",
  "* Tables/Class/localization_by_comparison_class_table.csv",
  "* Tables/Class/localization_by_comparison_class_labeled.csv",
  "* Tables/Conditions/localization_left_condition_table.csv",
  "* Tables/Conditions/localization_right_condition_table.csv",
  "* Tables/Class/comparison_level_long.csv",
  "* Tables/QA/localization_NA_diagnostics.csv",
  "* Tables/GSEA/all_GSEA_long.csv / themes_table.csv / theme_by_comparison.csv / theme_replication_scores.csv",
  if (nrow(ora_df)) "* Tables/ORA/all_ORA_long.csv / ORA_replication_scores.csv; Plots/ORA/ORA_replication_bar.svg" else "* (no ORA detected)"
)
writeLines(md, con = file.path(out_dir, "CONSOLIDATED_SUMMARY.md"))

# Artifact existence check across subfolders
expected <- c(
  file.path(dirs$plots_anat,  paste0("localization_", unique(theme_loc$db), ".svg")),
  file.path(dirs$plots_class, paste0("localization_", unique(theme_loc_by_class_lab$db), "_by_comparison_class_labeled.svg")),
  file.path(dirs$plots_cond,  paste0("localization_", unique(theme_loc_left_cond$db), "_tested_condition.svg")),
  file.path(dirs$plots_cond,  paste0("localization_", unique(theme_loc_right_cond$db), "_baseline_condition.svg"))
)
exists_tbl <- tibble::tibble(file = unique(expected), exists = file.exists(unique(expected)))
readr::write_csv(exists_tbl, file.path(out_dir, "artifact_existence_check.csv"))
print(exists_tbl)

message("Consolidation complete. See ZZ_consolidated for outputs.")









































# =======================
# CONSOLIDATION + SUMMARY
# =======================
# Primary: GSEA (gseGO/gseKEGG/Reactome) using NES/FDR
# Secondary: ORA (separate tables/plots)
# Includes class-labeled heatmaps (SR/RC/SC:n) and comparison-level matrices. [web:41]

# 0) Setup --------------------------------------------------------------------
pkgs <- c("dplyr","tidyr","readr","stringr","purrr","ggplot2","igraph","ggraph", "patchwork", "tibble", "clusterProfiler", "org.Mm.eg.db", "DOSE", "enrichplot", "ReactomePA")
missing_pkgs <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(missing_pkgs)) {
  message("Installing missing packages: ", paste(missing_pkgs, collapse=", "))
  tryCatch({
    install.packages(missing_pkgs, repos = c("https://cloud.r-project.org"))
  }, error = function(e) {
    warning("Primary CRAN failed: ", conditionMessage(e), " — retrying ggraph via r-universe")
    if ("ggraph" %in% missing_pkgs) install.packages("ggraph", repos = c("https://thomasp85.r-universe.dev","https://cloud.r-project.org"))
  })
}
suppressPackageStartupMessages({
  for (p in pkgs) {
    tryCatch(library(p, character.only = TRUE),
             error = function(e) warning(sprintf("Failed to load %s: %s", p, conditionMessage(e))))
  }
})

# Root paths (fixed) ----------------------------------------------------------
root_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Results/neuron-phenotypeWithinUnit"
out_dir  <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Results/ZZ_consolidated"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# New: subfolder layout for tables/plots --------------------------------------
dirs <- list(
  tables        = file.path(out_dir, "Tables"),
  tables_gsea   = file.path(out_dir, "Tables", "GSEA"),
  tables_ora    = file.path(out_dir, "Tables", "ORA"),
  tables_class  = file.path(out_dir, "Tables", "Class"),
  tables_cond   = file.path(out_dir, "Tables", "Conditions"),
  tables_qa     = file.path(out_dir, "Tables", "QA"),
  plots         = file.path(out_dir, "Plots"),
  plots_anat    = file.path(out_dir, "Plots", "Anatomy"),
  plots_class   = file.path(out_dir, "Plots", "Class"),
  plots_cond    = file.path(out_dir, "Plots", "Conditions"),
  plots_cm      = file.path(out_dir, "Plots", "ComparisonMatrices"),
  plots_ora     = file.path(out_dir, "Plots", "ORA")
)
invisible(lapply(dirs, function(d) dir.create(d, showWarnings = FALSE, recursive = TRUE)))  # recursive creation [web:73][web:68]

# --- MODULE BRANCH: paths + filename parser (after dirs creation) ------------
dirs$tables_modules <- file.path(out_dir, "Tables", "Modules")
dirs$plots_modules  <- file.path(out_dir, "Plots",  "Modules")
invisible(lapply(c(dirs$tables_modules, dirs$plots_modules), dir.create, recursive=TRUE, showWarnings=FALSE))  # create module dirs [web:73][web:68]

modules_root <- file.path(dirname(root_dir), "20_modules")
if (!dir.exists(modules_root)) warning("modules_root not found: ", modules_root)  # diag [web:73]

parse_module_fname <- function(fp) {
  bn <- basename(fp); noext <- sub("\\.csv$", "", bn, ignore.case = TRUE)
  # gsea: gseGO_BP_<module>_<LEFT>_<RIGHT>
  m1 <- regexec("^gse([A-Za-z0-9]+)_([A-Za-z0-9]+)_([^_]+)_(.+)$", noext, perl=TRUE)
  r1 <- regmatches(noext, m1)[[1]]
  if (length(r1)) {
    db1 <- r1[2]; db2 <- r1[3]; db <- if (toupper(db1)=="GO" && toupper(db2)=="BP") "GO_BP" else paste0(db1,"_",db2)
    mod <- r1[4]; rem <- r1[5]
    gpos <- gregexpr("_(CA1|CA2|CA3|DG)", rem, perl=TRUE)[[1]]
    split_pos <- if (length(gpos) && gpos[1] != -1) tail(gpos,1) else -1
    if (split_pos == -1) { toks <- strsplit(rem, "_")[[1]]; left <- paste(toks[1:(length(toks)-1)], collapse="_"); right <- toks[length(toks)] } else {
      left <- substr(rem, 1, split_pos-1); right <- substr(rem, split_pos+1, nchar(rem))
    }
    return(list(source="gsea", db=db, module=mod, context=paste0(left," vs ", right)))
  }
  # ORA_BP_<module>_<LEFT>_<RIGHT>
  m2 <- regexec("^ORA_([A-Za-z0-9]+)_([^_]+)_(.+)$", noext, perl=TRUE)
  r2 <- regmatches(noext, m2)[[1]]
  if (length(r2)) {
    db2 <- r2[2]; db <- if (toupper(db2)=="BP") "GO_BP" else db2
    mod <- r2[3]; rem <- r2[4]
    gpos <- gregexpr("_(CA1|CA2|CA3|DG)", rem, perl=TRUE)[[1]]
    split_pos <- if (length(gpos) && gpos[1] != -1) tail(gpos,1) else -1
    if (split_pos == -1) { toks <- strsplit(rem, "_")[[1]]; left <- paste(toks[1:(length(toks)-1)], collapse="_"); right <- toks[length(toks)] } else {
      left <- substr(rem, 1, split_pos-1); right <- substr(rem, split_pos+1, nchar(rem))
    }
    return(list(source="ora_deg", db=db, module=mod, context=paste0(left," vs ", right)))
  }
  # ORA_fullModule_BP_<module>_<LEFT>_<RIGHT>
  m3 <- regexec("^ORA_fullModule_([A-Za-z0-9]+)_([^_]+)_(.+)$", noext, perl=TRUE)
  r3 <- regmatches(noext, m3)[[1]]
  if (length(r3)) {
    db2 <- r3[2]; db <- if (toupper(db2)=="BP") "GO_BP" else db2
    mod <- r3[3]; rem <- r3[4]
    gpos <- gregexpr("_(CA1|CA2|CA3|DG)", rem, perl=TRUE)[[1]]
    split_pos <- if (length(gpos) && gpos[1] != -1) tail(gpos,1) else -1
    if (split_pos == -1) { toks <- strsplit(rem, "_")[[1]]; left <- paste(toks[1:(length(toks)-1)], collapse="_"); right <- toks[length(toks)] } else {
      left <- substr(rem, 1, split_pos-1); right <- substr(rem, split_pos+1, nchar(rem))
    }
    return(list(source="ora_full", db=db, module=mod, context=paste0(left," vs ", right)))
  }
  NULL
}  # regex parser for provided file patterns [web:60][web:59]

# 1) Helpers ------------------------------------------------------------------
read_if <- function(fp) {
  if (!file.exists(fp)) return(NULL)
  tryCatch(readr::read_csv(fp, guess_max = 10000, progress = FALSE, show_col_types = FALSE),
           error = function(e) { warning("Failed reading: ", fp, " — ", conditionMessage(e)); NULL })
}
read_context <- function(comp_dir) {
  readme <- file.path(comp_dir, "README.txt")
  if (file.exists(readme)) {
    rl <- readLines(readme, warn = FALSE)
    ctx <- rl[grepl("^Context: ", rl)]
    if (length(ctx)) return(sub("^Context:\\s*", "", ctx[1]))
  }
  basename(comp_dir)
}
`%||%` <- function(a,b) if (!is.null(a) && length(a)>0 && !is.na(a)) a else b
tokenize <- function(x){ x <- tolower(x); x <- gsub("[^a-z0-9]+"," ",x); unique(unlist(strsplit(x,"\\s+"))) }
jaccard  <- function(a,b){ ia<-tokenize(a); ib<-tokenize(b); if(!length(ia)||!length(ib)) return(0); length(intersect(ia,ib))/length(union(ia,ib)) }

# 2) GSEA collectors -----------------------------------------------------------
parse_gsea_tbl <- function(fp, db) {
  df <- read_if(fp); if (is.null(df) || !nrow(df)) return(NULL)
  nm <- names(df); has <- function(x) any(x %in% nm)
  if (!has("ID") || !has("Description")) return(NULL)
  if (!has("NES")) return(NULL)
  fdr_col <- c("p.adjust","qvalue","padj","FDR","adj.P.Val")[c("p.adjust","qvalue","padj","FDR","adj.P.Val") %in% nm][1]
  if (is.na(fdr_col)) return(NULL)
  out <- tibble::tibble(
    db = rep(db, nrow(df)),
    ID = df[["ID"]],
    Description = df[["Description"]],
    NES = suppressWarnings(as.numeric(df[["NES"]])),
    FDR = suppressWarnings(as.numeric(df[[fdr_col]]))
  ) %>% dplyr::filter(!is.na(ID), !is.na(NES), !is.na(FDR))
  if (!nrow(out)) return(NULL)
  out
}
empty_gsea <- function(){ tibble::tibble(db=character(), ID=character(), Description=character(), NES=double(), FDR=double()) }

collect_gsea_one <- function(comp_dir) {
  paths <- list(
    GO_BP    = file.path(comp_dir, "10_global","GO_BP","Tables"),
    KEGG     = file.path(comp_dir, "10_global","KEGG","Tables"),
    Reactome = file.path(comp_dir, "10_global","Reactome","Tables")
  )
  if (!any(dir.exists(unlist(paths)))) return(empty_gsea())
  gsego   <- if (dir.exists(paths$GO_BP))    list.files(paths$GO_BP,    pattern="^gseGO_BP_.*\\.csv$",      full.names=TRUE) else character(0)
  gsekegg <- if (dir.exists(paths$KEGG))     list.files(paths$KEGG,     pattern="^gseKEGG_.*\\.csv$",       full.names=TRUE) else character(0)
  rs      <- if (dir.exists(paths$Reactome)) list.files(paths$Reactome, pattern="^ReactomeGSEA_.*\\.csv$",  full.names=TRUE) else character(0)
  if (length(gsego)+length(gsekegg)+length(rs) == 0) return(empty_gsea())
  parts <- list(
    dplyr::bind_rows(lapply(gsego,   parse_gsea_tbl, db="GO_BP")),
    dplyr::bind_rows(lapply(gsekegg, parse_gsea_tbl, db="KEGG")),
    dplyr::bind_rows(lapply(rs,      parse_gsea_tbl, db="Reactome"))
  )
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (!length(parts)) return(empty_gsea())
  dplyr::bind_rows(parts)
}

# 3) ORA collector -------------------------------------------------------------
parse_ora_tbl <- function(fp, db) {
  df <- read_if(fp); if (is.null(df) || !nrow(df)) return(NULL)
  nm <- names(df); if (!("ID" %in% nm && "Description" %in% nm)) return(NULL)
  fdr_col <- c("p.adjust","qvalue","padj","FDR","adj.P.Val")[c("p.adjust","qvalue","padj","FDR","adj.P.Val") %in% nm][1]
  if (is.na(fdr_col)) return(NULL)
  out <- tibble::tibble(db=rep(db, nrow(df)), ID=df[["ID"]], Description=df[["Description"]], Padj=suppressWarnings(as.numeric(df[[fdr_col]]))) %>%
    dplyr::filter(!is.na(ID), !is.na(Padj))
  if (!nrow(out)) return(NULL)
  out
}
empty_ora <- function(){ tibble::tibble(db=character(), ID=character(), Description=character(), Padj=double()) }

collect_ora_one <- function(comp_dir) {
  paths <- list(
    GO_BP    = file.path(comp_dir, "10_global","GO_BP","Tables"),
    KEGG     = file.path(comp_dir, "10_global","KEGG","Tables"),
    Reactome = file.path(comp_dir, "10_global","Reactome","Tables")
  )
  if (!any(dir.exists(unlist(paths)))) return(empty_ora())
  files <- c(
    if (dir.exists(paths$GO_BP))    list.files(paths$GO_BP,    pattern="^ORA_.*\\.csv$",         full.names=TRUE) else character(0),
    if (dir.exists(paths$KEGG))     list.files(paths$KEGG,     pattern="^ORA_.*\\.csv$",         full.names=TRUE) else character(0),
    if (dir.exists(paths$Reactome)) list.files(paths$Reactome, pattern="^ReactomeORA_.*\\.csv$", full.names=TRUE) else character(0)
  )
  if (!length(files)) return(empty_ora())
  binders <- lapply(files, function(fp) {
    db <- if (grepl("/GO_BP/", fp)) "GO_BP" else if (grepl("/KEGG/", fp)) "KEGG" else if (grepl("/Reactome/", fp)) "Reactome" else "Unknown"
    parse_ora_tbl(fp, db = db)
  })
  binders <- binders[!vapply(binders, is.null, logical(1))]
  if (!length(binders)) return(empty_ora())
  dplyr::bind_rows(binders)
}

# 4) Scan comparison folders ---------------------------------------------------
comp_dirs <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
comp_dirs <- comp_dirs[dir.exists(file.path(comp_dirs, "10_global"))]
message(sprintf("Found %d comparison folders under neuron-phenotypeWithinUnit.", length(comp_dirs)))
if (!length(comp_dirs)) stop("No comparison folders with '10_global' found under the fixed root_dir.", call. = FALSE)

# 5) GSEA assembly -------------------------------------------------------------
gsea_list <- lapply(comp_dirs, function(cd) tibble::tibble(comp_dir = cd, context = read_context(cd), data = list(collect_gsea_one(cd))))
gsea_df <- purrr::map_dfr(gsea_list, function(x) { dx <- x$data[[1]]; if (!is.null(dx) && nrow(dx) > 0) cbind(x[1:2], dx) else NULL })
if (!nrow(gsea_df)) stop("No valid GSEA rows found across comparisons.", call. = FALSE)
req <- c("comp_dir","context","db","ID","Description","NES","FDR")
missing_cols <- setdiff(req, names(gsea_df))
if (length(missing_cols)) stop(sprintf("gsea_df missing columns: %s", paste(missing_cols, collapse=", ")), call. = FALSE)
gsea_df <- gsea_df %>% dplyr::mutate(comp_dir=as.character(comp_dir), context=as.character(context), db=as.character(db), ID=as.character(ID), Description=as.character(Description), NES=as.numeric(NES), FDR=as.numeric(FDR))
readr::write_csv(gsea_df, file.path(dirs$tables_gsea, "all_GSEA_long.csv"))
message(sprintf("GSEA rows: %d", nrow(gsea_df)))

# 6) ORA assembly --------------------------------------------------------------
ora_list <- lapply(comp_dirs, function(cd) tibble::tibble(comp_dir = cd, context = read_context(cd), data = list(collect_ora_one(cd))))
ora_df <- purrr::map_dfr(ora_list, function(x) { dx <- x$data[[1]]; if (!is.null(dx) && nrow(dx) > 0) cbind(x[1:2], dx) else NULL })
if (nrow(ora_df)) readr::write_csv(ora_df, file.path(dirs$tables_ora, "all_ORA_long.csv"))

# --- 6b) Module-level collectors and assembly --------------------------------
module_dirs <- if (dir.exists(modules_root)) list.dirs(modules_root, full.names=TRUE, recursive=FALSE) else character(0)
module_dirs <- module_dirs[grepl("/Module_", gsub("\\\\","/", module_dirs))]
message(sprintf("Found %d Module_* folders.", length(module_dirs)))

list_module_files <- function(module_dir) {
  tbl_dir <- file.path(module_dir, "Tables")
  if (!dir.exists(tbl_dir)) return(character(0))
  subdirs <- c("gsea","ora_deg","ora_full")
  unlist(lapply(subdirs, function(s) {
    d <- file.path(tbl_dir, s)
    if (dir.exists(d)) list.files(d, pattern="\\.csv$", full.names=TRUE) else character(0)
  }))
}

parse_module_csv <- function(fp) {
  meta <- parse_module_fname(fp); if (is.null(meta)) return(NULL)
  df <- read_if(fp); if (is.null(df) || !nrow(df)) return(NULL)
  nm <- names(df)
  if (identical(meta$source, "gsea")) {
    fdr_col <- c("p.adjust","qvalue","padj","FDR","adj.P.Val")[c("p.adjust","qvalue","padj","FDR","adj.P.Val") %in% nm][1]
    if (!all(c("ID","Description","NES") %in% nm) || is.na(fdr_col)) return(NULL)
    tibble::tibble(
      module = meta$module, db = meta$db, context = meta$context,
      ID = as.character(df$ID), Description = as.character(df$Description),
      NES = suppressWarnings(as.numeric(df$NES)), FDR = suppressWarnings(as.numeric(df[[fdr_col]]))
    ) %>% dplyr::filter(!is.na(ID), is.finite(NES), is.finite(FDR))
  } else {
    padj_col <- c("p.adjust","qvalue","padj","FDR","adj.P.Val")[c("p.adjust","qvalue","padj","FDR","adj.P.Val") %in% nm][1]
    if (!all(c("ID","Description") %in% nm) || is.na(padj_col)) return(NULL)
    tibble::tibble(
      module = meta$module, db = meta$db, context = meta$context,
      ID = as.character(df$ID), Description = as.character(df$Description),
      Padj = suppressWarnings(as.numeric(df[[padj_col]])), source = meta$source
    ) %>% dplyr::filter(!is.na(ID), is.finite(Padj))
  }
}

module_files <- unlist(lapply(module_dirs, list_module_files))
message(sprintf("Module CSV files detected: %d", length(module_files)))
mod_parsed <- lapply(module_files, parse_module_csv)
mod_parsed <- mod_parsed[!vapply(mod_parsed, is.null, logical(1))]

module_gsea_df <- purrr::map_dfr(mod_parsed, ~{ if (all(c("NES","FDR") %in% names(.x))) .x else NULL })
module_ora_df  <- purrr::map_dfr(mod_parsed, ~{ if ("Padj" %in% names(.x)) .x else NULL })

if (nrow(module_gsea_df)) {
  readr::write_csv(module_gsea_df, file.path(dirs$tables_modules, "module_GSEA_long.csv"))
  message(sprintf("Module GSEA rows: %d", nrow(module_gsea_df)))
} else message("No module GSEA rows.")
if (nrow(module_ora_df)) {
  readr::write_csv(module_ora_df, file.path(dirs$tables_modules, "module_ORA_long.csv"))
  message(sprintf("Module ORA rows: %d", nrow(module_ora_df)))
} else message("No module ORA rows.")

# 7) Theme building (GSEA only) -----------------------------------------------
build_themes_for_db <- function(df_db, top_n_terms = 200, sim_cut = 0.25) {
  if (is.null(df_db) || !nrow(df_db)) return(NULL)
  df_db <- dplyr::mutate(df_db, sig = FDR <= 0.05)
  agg <- df_db %>%
    dplyr::group_by(ID, Description, db) %>%
    dplyr::summarize(n_sig = sum(sig, na.rm = TRUE), n = dplyr::n(), meanNES = mean(NES, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(n_sig), dplyr::desc(abs(meanNES)))
  sel <- head(agg, top_n_terms)
  if (nrow(sel) < 2) return(dplyr::mutate(sel, theme = paste0(unique(df_db$db), "_", dplyr::row_number())))
  ids <- sel$ID; descs <- sel$Description; n <- length(ids)
  sim <- matrix(0, n, n); rownames(sim) <- colnames(sim) <- ids
  for (i in seq_len(n)) for (j in i:n) { s <- if (i==j) 1 else jaccard(descs[i], descs[j]); sim[i,j] <- s; sim[j,i] <- s }
  edge_df <- as.data.frame(as.table(sim), stringsAsFactors = FALSE) |>
    dplyr::rename(ID1 = Var1, ID2 = Var2, sim = Freq) |>
    dplyr::filter(ID1 != ID2, sim >= sim_cut)
  g <- igraph::make_empty_graph(directed = FALSE) |> igraph::add_vertices(n, name = ids)
  if (nrow(edge_df) > 0) {
    edges_mat <- unique(t(apply(edge_df[,c("ID1","ID2")], 1, function(x) sort(as.character(x)))))
    g <- igraph::add_edges(g, as.vector(t(edges_mat)))
  }
  if (igraph::vcount(g) == 0) { sel$theme <- paste0(unique(df_db$db), "_", seq_len(nrow(sel))); return(sel) }
  cl <- igraph::cluster_louvain(g)
  memb <- igraph::membership(cl)
  theme_map <- tibble::tibble(ID = names(memb), theme = paste0(unique(df_db$db), "_C", as.integer(memb)))
  out <- dplyr::left_join(sel, theme_map, by = "ID")
  if (any(is.na(out$theme))) out$theme[is.na(out$theme)] <- paste0(unique(df_db$db), "_", seq_len(sum(is.na(out$theme))))
  out
}

gsea_core <- gsea_df %>% dplyr::select(db, ID, Description, NES, FDR) %>% dplyr::mutate(db=as.character(db), ID=as.character(ID), Description=as.character(Description), NES=as.numeric(NES), FDR=as.numeric(FDR))
req_theme <- c("db","ID","Description","NES","FDR"); miss_theme <- setdiff(req_theme, names(gsea_core)); if (length(miss_theme)) stop(sprintf("gsea_core missing: %s", paste(miss_theme, collapse=", ")), call.=FALSE)
split_dbs <- split(gsea_core, gsea_core$db)
themes_db <- lapply(split_dbs, build_themes_for_db)
themes_db <- themes_db[!vapply(themes_db, is.null, logical(1))]
themes_tbl <- if (length(themes_db)) dplyr::bind_rows(themes_db) else tibble::tibble()
readr::write_csv(themes_tbl, file.path(dirs$tables_gsea, "themes_table.csv"))

# BP membership per theme (cluster), per DB
bp_members <- gsea_df %>%
  dplyr::filter(db == "GO_BP") %>%
  dplyr::inner_join(dplyr::select(themes_tbl, ID, theme), by="ID") %>%
  dplyr::group_by(db, theme, ID, Description) %>%
  dplyr::summarize(
    n_hits = dplyr::n(),
    meanNES = mean(NES, na.rm=TRUE),
    bestFDR = suppressWarnings(min(FDR, na.rm=TRUE)),
    .groups="drop"
  ) %>%
  dplyr::arrange(db, theme, dplyr::desc(n_hits), bestFDR)
readr::write_csv(bp_members, file.path(dirs$tables_gsea, "GO_BP_terms_by_theme_long.csv"))

bp_split <- split(bp_members, interaction(bp_members$db, bp_members$theme, drop=TRUE))
invisible(purrr::iwalk(bp_split, function(df, nm){
  parts <- strsplit(nm, "\\.")[[1]]
  dbn <- parts[1]; th <- parts[2]
  fp <- file.path(dirs$tables_gsea, paste0("GO_BP_terms_", th, ".csv"))
  readr::write_csv(df, fp)
}))

# 8) Aggregation per comparison x theme (GSEA) --------------------------------
assign_themes <- if (nrow(themes_tbl)) dplyr::inner_join(gsea_df, dplyr::select(themes_tbl, ID, theme), by = "ID") else tibble::tibble()
if (!nrow(assign_themes)) stop("No themes constructed (themes_tbl empty). Consider lowering sim_cut or increasing top_n_terms.", call. = FALSE)

theme_comp <- assign_themes %>%
  dplyr::group_by(db, theme, comp_dir, context) %>%
  dplyr::summarize(n_terms=dplyr::n(), meanNES=mean(NES, na.rm=TRUE), minFDR=suppressWarnings(min(FDR, na.rm=TRUE)), frac_sig=mean(FDR<=0.05, na.rm=TRUE), .groups="drop")
readr::write_csv(theme_comp, file.path(dirs$tables_gsea, "theme_by_comparison.csv"))

theme_comp$NES_masked <- ifelse(theme_comp$minFDR <= 0.05, theme_comp$meanNES, NA_real_)

# 9) Robust parsing — both sides (fixed) --------------------------------------
parse_both <- function(ctx) {
  parts <- unlist(strsplit(ctx, "(?:\\s+vs\\s+|_vs_)", perl=TRUE))
  parts <- trimws(parts)
  left  <- if (length(parts) >= 1) parts[1] else NA_character_
  right <- if (length(parts) >= 2) parts[2] else NA_character_
  region_pat <- "(DG|CA1|CA2|CA3)"
  base_layer <- "(mo|po|sr|sp|so|slm|sg)"
  cond_pat   <- "(res|sus|con)"
  parse_side <- function(side) {
    if (is.na(side) || side == "") return(list(region=NA_character_, layer=NA_character_, condition=NA_character_))
    m <- regexec(paste0("^", region_pat, "[_-]?", base_layer, "(", cond_pat, ")?$"), side, perl=TRUE)
    r <- regmatches(side, m)[[1]]
    if (length(r) == 0) {
      m2 <- regexec(paste0("^", region_pat, "[_-]?", base_layer, "(", cond_pat, ")?"), side, perl=TRUE)
      r  <- regmatches(side, m2)[[1]]
    }
    region <- NA_character_; layer <- NA_character_; condition <- NA_character_
    if (length(r) >= 3) {
      region    <- r[2]
      layer_raw <- tolower(r[3]); layer <- layer_raw
      if (length(r) >= 5 && nzchar(r[5])) condition <- tolower(r[5])
    } else {
      toks <- unlist(strsplit(side, "[_-]+"))
      region <- if (length(toks) >= 1 && grepl(paste0("^", region_pat, "$"), toks[1])) toks[1] else NA_character_
      layer  <- if (length(toks) >= 2 && grepl(paste0("^", base_layer, "$"), tolower(toks[2]))) tolower(toks[2]) else NA_character_
      condition <- if (length(toks) >= 3 && grepl(paste0("^", cond_pat, "$"), tolower(toks[3]))) tolower(toks[3]) else NA_character_
    }
    list(region=region, layer=layer, condition=condition)
  }
  L <- parse_side(left)
  R <- parse_side(right)
  list(L=L, R=R)
}

rl_both <- unique(theme_comp$context) %>%
  purrr::map_df(~{
    pr <- parse_both(.x)
    tibble::tibble(
      context   = .x,
      region_L  = pr$L$region, layer_L  = pr$L$layer, condition_L  = pr$L$condition,
      region_R  = pr$R$region, layer_R  = pr$R$layer, condition_R  = pr$R$condition
    )
  })
readr::write_csv(rl_both, file.path(dirs$tables, "context_region_layer_condition_both_sides.csv"))

# 10) Comparison classes (from parsed conditions) -----------------------------
comp_class_map <- rl_both %>%
  dplyr::transmute(
    context,
    cond_left  = condition_L,
    cond_right = condition_R,
    comparison_class = dplyr::case_when(
      cond_left == "sus" & cond_right == "res" ~ "sus vs res",
      cond_left == "res" & cond_right == "con" ~ "res vs con",
      cond_left == "sus" & cond_right == "con" ~ "sus vs con",
      TRUE ~ NA_character_
    )
  )
readr::write_csv(comp_class_map, file.path(dirs$tables_class, "context_comparison_class_map.csv"))

# --- 10b) Module context + class maps ----------------------------------------
if (exists("module_gsea_df") && nrow(module_gsea_df)) {
  mod_ctx <- unique(module_gsea_df$context)
} else if (exists("module_ora_df") && nrow(module_ora_df)) {
  mod_ctx <- unique(module_ora_df$context)
} else mod_ctx <- character(0)

if (length(mod_ctx)) {
  mod_rl_both <- tibble::tibble(context = mod_ctx) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(parsed = list(parse_both(context))) %>%
    dplyr::mutate(
      region_L  = parsed$L$region, layer_L  = parsed$L$layer, condition_L  = parsed$L$condition,
      region_R  = parsed$R$region, layer_R  = parsed$R$layer, condition_R  = parsed$R$condition
    ) %>% dplyr::ungroup() %>% dplyr::select(-parsed)
  readr::write_csv(mod_rl_both, file.path(dirs$tables_modules, "module_context_region_layer_condition_both_sides.csv"))
  mod_comp_class_map <- mod_rl_both %>%
    dplyr::transmute(
      context,
      cond_left  = condition_L,
      cond_right = condition_R,
      comparison_class = dplyr::case_when(
        cond_left == "sus" & cond_right == "res" ~ "sus vs res",
        cond_left == "res" & cond_right == "con" ~ "res vs con",
        cond_left == "sus" & cond_right == "con" ~ "sus vs con",
        TRUE ~ NA_character_
      )
    )
  readr::write_csv(mod_comp_class_map, file.path(dirs$tables_modules, "module_context_comparison_class_map.csv"))
} else {
  mod_rl_both <- tibble::tibble()
  mod_comp_class_map <- tibble::tibble()
}

# 11) Localization tables (global) ... and all your existing global sections 11–16 remain unchanged

# ======================
# MODULE ANALYSES PARALLEL
# ======================
if (exists("module_gsea_df") && nrow(module_gsea_df) && nrow(mod_rl_both)) {
  # Aggregation per module x comparison (region-layer on left side)
  mod_theme_comp <- module_gsea_df %>%
    dplyr::left_join(mod_rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
    dplyr::filter(!is.na(region), !is.na(layer)) %>%
    dplyr::group_by(db, module, context, region, layer) %>%
    dplyr::summarize(meanNES = mean(NES, na.rm=TRUE),
                     minFDR  = suppressWarnings(min(FDR, na.rm=TRUE)),
                     .groups="drop")
  readr::write_csv(mod_theme_comp, file.path(dirs$tables_modules, "module_theme_by_comparison.csv"))

  # Localization table per module
  module_loc <- mod_theme_comp %>%
    dplyr::group_by(db, module, region, layer) %>%
    dplyr::summarize(meanNES = mean(meanNES, na.rm=TRUE),
                     frac_sig = mean(minFDR<=0.05, na.rm=TRUE),
                     .groups="drop")
  readr::write_csv(module_loc, file.path(dirs$tables_modules, "module_localization_table.csv"))

  # Replication per module
  module_repl <- mod_theme_comp %>%
    dplyr::mutate(sign=sign(meanNES), sig=minFDR<=0.05) %>%
    dplyr::group_by(db, module) %>%
    dplyr::summarize(n_comp=dplyr::n(),
                     n_sig=sum(sig, na.rm=TRUE),
                     n_pos=sum(sig & sign>0, na.rm=TRUE),
                     n_neg=sum(sig & sign<0, na.rm=TRUE),
                     dom_dir=ifelse(n_pos>=n_neg,"up","down"),
                     repl_score=pmax(n_pos,n_neg),
                     .groups="drop") %>%
    dplyr::arrange(dplyr::desc(repl_score))
  readr::write_csv(module_repl, file.path(dirs$tables_modules, "module_replication_scores.csv"))

  # Class-labeled localization
  module_loc_by_class <- mod_theme_comp %>%
    dplyr::left_join(mod_comp_class_map %>% dplyr::select(context, comparison_class), by="context") %>%
    dplyr::filter(!is.na(comparison_class)) %>%
    dplyr::group_by(db, module, comparison_class, region, layer) %>%
    dplyr::summarize(meanNES = mean(meanNES, na.rm=TRUE),
                     n_comp  = dplyr::n_distinct(context),
                     .groups="drop")
  readr::write_csv(module_loc_by_class, file.path(dirs$tables_modules, "module_localization_by_comparison_class.csv"))

  module_class_labels <- mod_theme_comp %>%
    dplyr::left_join(mod_comp_class_map, by="context") %>%
    dplyr::filter(!is.na(comparison_class)) %>%
    dplyr::group_by(db, module, region, layer, comparison_class) %>%
    dplyr::summarize(n_comp = dplyr::n_distinct(context), .groups="drop") %>%
    dplyr::mutate(code = dplyr::case_when(
      comparison_class=="sus vs res" ~ "SR",
      comparison_class=="res vs con" ~ "RC",
      comparison_class=="sus vs con" ~ "SC",
      TRUE ~ ""
    ),
    label = paste0(code, ":", n_comp))
  module_loc_by_class_lab <- module_loc_by_class %>%
    dplyr::left_join(module_class_labels %>% dplyr::select(db,module,region,layer,comparison_class,label),
                     by=c("db","module","region","layer","comparison_class"))
  readr::write_csv(module_loc_by_class_lab, file.path(dirs$tables_modules, "module_localization_by_comparison_class_labeled.csv"))

  # Driver strengths SR vs SC
  mod_sr_sc <- mod_theme_comp %>%
    dplyr::left_join(mod_comp_class_map %>% dplyr::select(context, comparison_class), by="context") %>%
    dplyr::filter(!is.na(comparison_class)) %>%
    dplyr::group_by(db, module, region, layer, comparison_class) %>%
    dplyr::summarize(meanNES_class = mean(meanNES, na.rm=TRUE),
                     n_ctx = dplyr::n_distinct(context), .groups="drop") %>%
    tidyr::pivot_wider(names_from = comparison_class, values_from = c(meanNES_class, n_ctx), values_fill = 0) %>%
    dplyr::rename(
      meanNES_class_sus_vs_res = `meanNES_class_sus vs res`,
      meanNES_class_sus_vs_con = `meanNES_class_sus vs con`,
      n_ctx_sus_vs_res         = `n_ctx_sus vs res`,
      n_ctx_sus_vs_con         = `n_ctx_sus vs con`
    ) %>%
    dplyr::mutate(
      SR_strength = meanNES_class_sus_vs_res,
      SC_strength = meanNES_class_sus_vs_con,
      driver_class = dplyr::case_when(
        abs(SR_strength) > abs(SC_strength) ~ "sus vs res",
        abs(SC_strength) > abs(SR_strength) ~ "sus vs con",
        TRUE ~ "tie"
      ),
      driver_gap = abs(SR_strength) - abs(SC_strength)
    )
  readr::write_csv(mod_sr_sc, file.path(dirs$tables_modules, "module_driver_strength.csv"))

  # Comparison-level long with class label
  mod_comp_level <- mod_theme_comp %>%
    dplyr::left_join(mod_comp_class_map %>% dplyr::select(context, comparison_class), by="context") %>%
    dplyr::filter(!is.na(comparison_class)) %>%
    dplyr::mutate(ctx_label = paste0(comparison_class, "::", context))
  readr::write_csv(mod_comp_level, file.path(dirs$tables_modules, "module_comparison_level_long.csv"))

  # Direction summaries (dominant direction and symmetry)
  mod_dom_dir <- mod_theme_comp %>%
    dplyr::left_join(mod_rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
    dplyr::filter(!is.na(region), !is.na(layer), minFDR<=0.05) %>%
    dplyr::group_by(db, module, region, layer) %>%
    dplyr::summarize(medNES = stats::median(meanNES, na.rm=TRUE),
                     n_sig  = dplyr::n(), .groups="drop") %>%
    dplyr::mutate(dir = dplyr::case_when(medNES > 0 ~ "up", medNES < 0 ~ "down", TRUE ~ "flat"))
  readr::write_csv(mod_dom_dir, file.path(dirs$tables_modules, "module_dominant_direction.csv"))

  mod_symmetry <- mod_theme_comp %>%
    dplyr::left_join(mod_rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
    dplyr::filter(!is.na(region), !is.na(layer)) %>%
    dplyr::group_by(db, region, layer) %>%
    dplyr::summarize(p_sym = tryCatch({
        x <- meanNES; x <- x[is.finite(x)]
        if (length(x) >= 6) stats::wilcox.test(x, mu=0, exact=FALSE)$p.value else NA_real_
      }, error=function(e) NA_real_), .groups="drop")
  readr::write_csv(mod_symmetry, file.path(dirs$tables_modules, "module_directional_bias_pvalues.csv"))

  # ----------------
  # Module plots
  # ----------------
  # Localization heatmaps for top modules
  plot_module_localization <- function(dbn, n_mod=12) {
    tops <- module_repl %>% dplyr::filter(db==dbn) %>% dplyr::slice_max(repl_score, n=n_mod) %>% dplyr::pull(module)
    df <- module_loc %>% dplyr::filter(db==dbn, module %in% tops)
    if (!nrow(df)) return(NULL)
    p <- ggplot2::ggplot(df, ggplot2::aes(x=layer,y=region,fill=meanNES)) +
      ggplot2::geom_tile(color="grey90") +
      ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b", midpoint=0) +
      ggplot2::facet_wrap(~ module, ncol=4, labeller = ggplot2::labeller(module = ggplot2::label_wrap_gen(18))) +
      ggplot2::labs(title=paste0("Module localization — ", dbn), x="Layer", y="Region") +
      ggplot2::theme_minimal(10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))
    fp <- file.path(dirs$plots_modules, paste0("module_localization_", dbn, ".svg"))
    ggplot2::ggsave(fp, p, width=22, height=14); message("[module-loc] Saved: ", fp)
    invisible(p)
  }
  invisible(lapply(unique(module_loc$db), plot_module_localization))

  # Plots: class-labeled module localization
  plot_module_localization_by_class_labeled <- function(dbn) {
    df_all <- module_loc_by_class_lab %>% dplyr::filter(db==dbn)
    tops <- module_repl %>% dplyr::filter(db==dbn) %>% dplyr::slice_max(repl_score, n=12) %>% dplyr::pull(module)
    if (!length(tops) || !any(df_all$module %in% tops)) {
      tops <- df_all %>% dplyr::count(module, sort=TRUE) %>% dplyr::slice_head(n=12) %>% dplyr::pull(module)
    }
    df <- df_all %>% dplyr::filter(module %in% tops)
    if (!nrow(df)) { message("[module-class-labeled] No data for ", dbn, " after fallback."); return(NULL) }

    class_levels <- c("sus vs res","res vs con","sus vs con")
    df$comparison_class <- factor(df$comparison_class, levels = class_levels)

    p <- ggplot2::ggplot(df, ggplot2::aes(x=layer, y=region, fill=meanNES)) +
      ggplot2::geom_tile(color="grey90") +
      ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b", midpoint=0) +
      ggplot2::facet_grid(
        rows = ggplot2::vars(module),
        cols = ggplot2::vars(comparison_class),
        scales = "free",
        drop = FALSE,
        labeller = ggplot2::labeller(module = ggplot2::label_wrap_gen(width = 18))
      ) +
      ggplot2::labs(title=paste0("Module localization — comparison class (labels=SR/RC/SC:n) — ", dbn),
                    x="Layer", y="Region") +
      ggplot2::theme_minimal(base_size=10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1)) +
      ggplot2::geom_text(ggplot2::aes(label = dplyr::coalesce(label, as.character(n_comp))),
                         color="black", size=3, na.rm=TRUE)

    fp <- file.path(dirs$plots_modules, paste0("module_localization_", dbn, "_by_comparison_class_labeled.svg"))
    ggplot2::ggsave(fp, p, width=30, height=20); message("[module-class-labeled] Saved: ", fp)
    invisible(p)
  }
  invisible(lapply(unique(module_loc_by_class_lab$db), plot_module_localization_by_class_labeled))

  # Plots: module comparison matrices (by module)
  plot_module_by_comparisons <- function(dbn, max_cols=60) {
    df_db <- mod_comp_level %>% dplyr::filter(db==dbn)
    top_mods <- module_repl %>% dplyr::filter(db==dbn) %>% dplyr::slice_max(repl_score, n=8) %>% dplyr::pull(module)
    if (!length(top_mods)) { message("[mod-comp-matrix] No top modules for ", dbn); return(NULL) }
    lapply(top_mods, function(md) {
      df <- df_db %>% dplyr::filter(module==md)
      ctx_ord <- df %>% dplyr::count(ctx_label, sort=TRUE) %>% dplyr::slice_head(n=max_cols) %>% dplyr::pull(ctx_label)
      df <- df %>% dplyr::filter(ctx_label %in% ctx_ord)
      df$ctx_label <- factor(df$ctx_label, levels=ctx_ord)
      p <- ggplot2::ggplot(df, ggplot2::aes(x=ctx_label, y=paste(region, layer, sep="_"), fill=meanNES)) +
        ggplot2::geom_tile(color="grey90") +
        ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b", midpoint=0) +
        ggplot2::labs(title=paste0("Module ", md, " — comparison-level NES — ", dbn),
                      x="Comparison (class::context)", y="Region_Layer") +
        ggplot2::theme_minimal(base_size=10) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle=60, hjust=1, vjust=1))
      fp <- file.path(dirs$plots_modules, paste0("module_", gsub("[^A-Za-z0-9]+","_", md), "_", dbn, "_comparison_matrix.svg"))
      ggplot2::ggsave(fp, p, width=28, height=10); message("[mod-comp-matrix] Saved: ", fp)
      invisible(p)
    })
  }
  invisible(lapply(unique(mod_comp_level$db), plot_module_by_comparisons))

  # Plots: module driver map (SR vs SC)
  plot_module_driver_map <- function(dbn) {
    top_mods <- module_repl %>% dplyr::filter(db==dbn) %>% dplyr::slice_max(repl_score, n=12) %>% dplyr::pull(module)
    df <- mod_sr_sc %>% dplyr::filter(db==dbn, module %in% top_mods)
    if (!nrow(df)) return(NULL)
    df <- df %>% dplyr::mutate(n_lab = paste0(n_ctx_sus_vs_res, "/", n_ctx_sus_vs_con))
    p <- ggplot2::ggplot(df, ggplot2::aes(x=layer, y=region, fill=driver_class, alpha=abs(driver_gap))) +
      ggplot2::geom_tile(color="white") +
      ggplot2::geom_text(ggplot2::aes(label = n_lab), size=3) +
      ggplot2::scale_fill_manual(values = c("sus vs res" = "#377eb8", "sus vs con" = "#e41a1c", "tie" = "grey80")) +
      ggplot2::scale_alpha(range = c(0.3, 1), guide = "none") +
      ggplot2::facet_wrap(~ module, ncol = 4, labeller = ggplot2::labeller(module = ggplot2::label_wrap_gen(18))) +
      ggplot2::labs(title = paste0(dbn, " — module drivers (SR vs SC)"),
                    x = "Layer", y = "Region", fill = "Class") +
      ggplot2::theme_minimal(10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))
    fp <- file.path(dirs$plots_modules, paste0("module_drivers_susceptibility_", dbn, ".svg"))
    ggplot2::ggsave(fp, p, width=22, height=12); message("[mod-drivers] Saved: ", fp)
    invisible(p)
  }
  invisible(lapply(unique(mod_sr_sc$db), plot_module_driver_map))
}

# 12) NA diagnostics -----------------------------------------------------------
cov_tbl <- assign_themes %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by = "context") %>%
  dplyr::filter(!is.na(region), !is.na(layer)) %>%
  dplyr::group_by(db, theme, region, layer) %>%
  dplyr::summarize(n_comp_cover=dplyr::n_distinct(comp_dir), n_sig_cover=sum(FDR<=0.05, na.rm=TRUE), minFDR_all=suppressWarnings(min(FDR, na.rm=TRUE)), .groups="drop")
na_diag <- theme_loc %>%
  dplyr::full_join(cov_tbl, by = c("db","theme","region","layer")) %>%
  dplyr::mutate(cause = dplyr::case_when(
    is.na(meanNES) & (is.na(n_comp_cover) | n_comp_cover == 0) ~ "no_coverage",
    is.na(meanNES) & (n_comp_cover > 0) & (is.na(n_sig_cover) | n_sig_cover == 0) ~ "masked_nonsignificant",
    is.na(meanNES) ~ "other_missing",
    TRUE ~ "has_value"
  ))
readr::write_csv(na_diag, file.path(dirs$tables_qa, "localization_NA_diagnostics.csv"))

# 13) Replication summary ------------------------------------------------------
theme_repl <- theme_comp %>% dplyr::mutate(sign=sign(meanNES), sig=minFDR<=0.05) %>%
  dplyr::group_by(db, theme) %>% dplyr::summarize(n_comp=dplyr::n(), n_sig=sum(sig, na.rm=TRUE), n_pos=sum(sig & sign>0, na.rm=TRUE), n_neg=sum(sig & sign<0, na.rm=TRUE), dom_dir=ifelse(n_pos>=n_neg,"up","down"), repl_score=pmax(n_pos,n_neg), .groups="drop") %>%
  dplyr::arrange(dplyr::desc(repl_score))
readr::write_csv(theme_repl, file.path(dirs$tables_gsea, "theme_replication_scores.csv"))

top_themes <- theme_repl %>% dplyr::group_by(db) %>% dplyr::slice_max(order_by=repl_score, n=12) %>% dplyr::ungroup()

# 14) Plots -------------------------------------------------------------------
plot_theme_localization <- function(dbn) {
  df <- theme_loc %>% dplyr::filter(db==dbn, theme %in% (top_themes %>% dplyr::filter(db==dbn) %>% dplyr::pull(theme)))
  if (!nrow(df)) { message("[anatomy] No data for ", dbn); return(NULL) }
  p <- ggplot2::ggplot(df, ggplot2::aes(x=layer, y=region, fill=meanNES)) +
    ggplot2::geom_tile(color="grey90") +
    ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b") +
    ggplot2::facet_wrap(~ theme, scales="free", ncol=4) +
    ggplot2::labs(title=paste0("Localization by region-layer — ", dbn), x="Layer", y="Region") +
    ggplot2::theme_minimal(base_size=10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))
  fp <- file.path(dirs$plots_anat, paste0("localization_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=22, height=14); message("[anatomy] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(theme_loc$db), plot_theme_localization))

plot_theme_localization_by_class_labeled <- function(dbn) {
  df_all <- theme_loc_by_class_lab %>% dplyr::filter(db==dbn)
  tops <- top_themes %>% dplyr::filter(db==dbn) %>% dplyr::pull(theme)
  if (!length(tops) || !any(df_all$theme %in% tops)) {
    tops <- df_all %>% dplyr::count(theme, sort=TRUE) %>% dplyr::slice_head(n=12) %>% dplyr::pull(theme)
  }
  df <- df_all %>% dplyr::filter(theme %in% tops)
  if (!nrow(df)) { message("[class-labeled] No data for ", dbn, " after fallback."); return(NULL) }
  class_levels <- c("sus vs res","res vs con","sus vs con")
  df$comparison_class <- factor(df$comparison_class, levels = class_levels)
  p <- ggplot2::ggplot(df, ggplot2::aes(x=layer, y=region, fill=meanNES)) +
    ggplot2::geom_tile(color="grey90") +
    ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b") +
    ggplot2::facet_grid(rows = ggplot2::vars(theme), cols = ggplot2::vars(comparison_class), scales="free", drop=FALSE) +
    ggplot2::labs(title=paste0("Localization — comparison class (labeled) — ", dbn), x="Layer", y="Region") +
    ggplot2::theme_minimal(base_size=10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1)) +
    ggplot2::geom_text(ggplot2::aes(label = dplyr::coalesce(label, as.character(n_comp))), color="black", size=3, na.rm=TRUE)
  fp <- file.path(dirs$plots_class, paste0("localization_", dbn, "_by_comparison_class_labeled.svg"))
  ggplot2::ggsave(fp, p, width=30, height=20); message("[class-labeled] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(theme_loc_by_class_lab$db), plot_theme_localization_by_class_labeled))

plot_theme_localization_left_cond <- function(dbn) {
  df <- theme_loc_left_cond %>% dplyr::filter(db==dbn, theme %in% (top_themes %>% dplyr::filter(db==dbn) %>% dplyr::pull(theme)))
  if (!nrow(df)) { message("[left-cond] No data for ", dbn); return(NULL) }
  p <- ggplot2::ggplot(df, ggplot2::aes(x=layer, y=region, fill=meanNES)) +
    ggplot2::geom_tile(color="grey90") +
    ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b") +
    ggplot2::facet_grid(rows = ggplot2::vars(theme), cols = ggplot2::vars(condition), scales="free") +
    ggplot2::labs(title=paste0("Localization — tested condition (sus/res) — ", dbn), x="Layer", y="Region") +
    ggplot2::theme_minimal(base_size=10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))
  fp <- file.path(dirs$plots_cond, paste0("localization_", dbn, "_tested_condition.svg"))
  ggplot2::ggsave(fp, p, width=26, height=20); message("[left-cond] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(theme_loc_left_cond$db), plot_theme_localization_left_cond))

plot_theme_localization_right_cond <- function(dbn) {
  df <- theme_loc_right_cond %>% dplyr::filter(db==dbn, theme %in% (top_themes %>% dplyr::filter(db==dbn) %>% dplyr::pull(theme)))
  if (!nrow(df)) { message("[right-cond] No data for ", dbn); return(NULL) }
  p <- ggplot2::ggplot(df, ggplot2::aes(x=layer, y=region, fill=meanNES)) +
    ggplot2::geom_tile(color="grey90") +
    ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b") +
    ggplot2::facet_grid(rows = ggplot2::vars(theme), cols = ggplot2::vars(condition), scales="free") +
    ggplot2::labs(title=paste0("Localization — baseline condition (res/con) — ", dbn), x="Layer", y="Region") +
    ggplot2::theme_minimal(base_size=10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))
  fp <- file.path(dirs$plots_cond, paste0("localization_", dbn, "_baseline_condition.svg"))
  ggplot2::ggsave(fp, p, width=26, height=20); message("[right-cond] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(theme_loc_right_cond$db), plot_theme_localization_right_cond))

plot_theme_by_comparisons <- function(dbn, max_cols=60) {
  df_db <- comp_level %>% dplyr::filter(db==dbn)
  themes_here <- df_db %>% dplyr::count(theme, sort=TRUE) %>% dplyr::slice_head(n=12) %>% dplyr::pull(theme)
  if (!length(themes_here)) { message("[comp-matrix] No themes with class data for ", dbn); return(NULL) }
  lapply(themes_here, function(th) {
    df <- df_db %>% dplyr::filter(theme==th)
    ctx_ord <- df %>% dplyr::count(ctx_label, sort=TRUE) %>% dplyr::slice_head(n=max_cols) %>% dplyr::pull(ctx_label)
    df <- df %>% dplyr::filter(ctx_label %in% ctx_ord)
    df$ctx_label <- factor(df$ctx_label, levels=ctx_ord)
    p <- ggplot2::ggplot(df, ggplot2::aes(x=ctx_label, y=paste(region, layer, sep="_"), fill=meanNES)) +
      ggplot2::geom_tile(color="grey90") +
      ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b") +
      ggplot2::labs(title=paste0("Theme ", th, " — comparison-level NES — ", dbn), x="Comparison (class::context)", y="Region_Layer") +
      ggplot2::theme_minimal(base_size=10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle=60, hjust=1, vjust=1))
    fp <- file.path(dirs$plots_cm, paste0("theme_", gsub("[^A-Za-z0-9]+","_", th), "_", dbn, "_comparison_matrix.svg"))
    ggplot2::ggsave(fp, p, width=28, height=10); message("[comp-matrix] Saved: ", fp)
    invisible(p)
  })
}
invisible(lapply(unique(comp_level$db), plot_theme_by_comparisons))

# 14b) Compact panel (Anatomy + Class + Replication) per DB
plot_compact_panel <- function(dbn, n_themes = 6) {
  the <- top_themes %>% dplyr::filter(db==dbn) %>% dplyr::slice_head(n=n_themes) %>% dplyr::pull(theme)

  p1 <- ggplot2::ggplot(theme_loc %>% dplyr::filter(db==dbn, theme %in% the),
                        ggplot2::aes(x=layer,y=region,fill=meanNES)) +
    ggplot2::geom_tile(color="grey90") +
    ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b") +
    ggplot2::facet_wrap(~ theme, ncol=3) +
    ggplot2::labs(title=paste0(dbn, " — localization (NES)")) + ggplot2::theme_minimal(9)

  p2 <- ggplot2::ggplot(theme_loc_by_class_lab %>% dplyr::filter(db==dbn, theme %in% the),
                        ggplot2::aes(x=layer,y=region,fill=meanNES,label=label)) +
    ggplot2::geom_tile(color="grey95") +
    ggplot2::geom_text(size=2, color="black", na.rm=TRUE) +
    ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b") +
    ggplot2::facet_grid(rows=ggplot2::vars(theme), cols=ggplot2::vars(comparison_class), drop=FALSE) +
    ggplot2::labs(title="comparison class (labels=SR/RC/SC:n)") + ggplot2::theme_minimal(9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45,hjust=1))

  p3 <- ggplot2::ggplot(theme_repl %>% dplyr::filter(db==dbn, theme %in% the),
                        ggplot2::aes(x=reorder(theme, repl_score), y=repl_score, fill=dom_dir)) +
    ggplot2::geom_col() + ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values=c(down="#2166ac", up="#b2182b")) +
    ggplot2::labs(title="replication score", x="", y="max(up, down) significant") +
    ggplot2::theme_minimal(9)

  if (requireNamespace("patchwork", quietly=TRUE)) {
    panel <- p1 / p2 / p3 + patchwork::plot_layout(heights=c(1,1.1,0.6))
    fp <- file.path(dirs$plots, "Panels", paste0("panel_", dbn, "_compact.svg"))
    dir.create(dirname(fp), recursive=TRUE, showWarnings=FALSE)
    ggplot2::ggsave(fp, panel, width=14, height=16); message("[panel] Saved: ", fp)
  } else {
    message("[panel] patchwork not installed; skipping combined panel for ", dbn)
  }
}
invisible(lapply(unique(theme_loc$db), plot_compact_panel))

# 14c) Dominant direction lattice per DB (median NES, n_sig annotation)
plot_dir_lattice <- function(dbn) {
  df <- dom_dir_tbl %>% dplyr::filter(db==dbn, theme %in% (top_themes %>% dplyr::filter(db==dbn) %>% dplyr::pull(theme)))
  if (!nrow(df)) return(NULL)
  p <- ggplot2::ggplot(df, ggplot2::aes(x=layer, y=region, fill=dir, label=n_sig)) +
    ggplot2::geom_tile(color="white") + ggplot2::geom_text(size=3, color="black") +
    ggplot2::scale_fill_manual(values=c(down="#2166ac", flat="grey85", up="#b2182b")) +
    ggplot2::facet_wrap(~ theme, ncol=4) +
    ggplot2::labs(title=paste0(dbn, " — dominant direction (median NES, n_sig)"),
                  x="Layer", y="Region") + ggplot2::theme_minimal(9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))
  fp <- file.path(dirs$plots_anat, paste0("dominant_direction_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=22, height=12); message("[dir-lattice] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(theme_loc$db), plot_dir_lattice))

# 14d) Directional bias tiles per DB (−log10 p from Wilcoxon against 0)
plot_symmetry <- function(dbn) {
  df <- symmetry_test %>% dplyr::filter(db==dbn)
  p <- ggplot2::ggplot(df, ggplot2::aes(x=layer, y=region, fill=-log10(p_sym))) +
    ggplot2::geom_tile(color="white") +
    ggplot2::scale_fill_viridis_c(option="C", na.value="grey90") +
    ggplot2::labs(title=paste0(dbn, " — directional bias (−log10 p)"), x="Layer", y="Region") +
    ggplot2::theme_minimal(9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))
  fp <- file.path(dirs$plots_anat, paste0("directional_bias_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=8, height=6); message("[bias] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(theme_loc$db), plot_symmetry))

# 14e) Comparison-class balance per theme (stacked fractions across region-layer)
class_balance <- theme_comp %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
  dplyr::left_join(comp_class_map, by="context") %>%
  dplyr::filter(!is.na(region), !is.na(layer), !is.na(comparison_class)) %>%
  dplyr::group_by(db, theme, region, layer, comparison_class) %>%
  dplyr::summarize(n_cmp=dplyr::n_distinct(context), .groups="drop") %>%
  dplyr::group_by(db, theme, region, layer) %>%
  dplyr::mutate(frac = n_cmp/sum(n_cmp)) %>% dplyr::ungroup()
readr::write_csv(class_balance, file.path(dirs$tables_class, "comparison_class_balance.csv"))

plot_class_balance <- function(dbn, theme_id) {
  df <- class_balance %>% dplyr::filter(db==dbn, theme==theme_id)
  if (!nrow(df)) return(NULL)
  p <- ggplot2::ggplot(df, ggplot2::aes(x=paste(region,layer,sep="_"), y=frac, fill=comparison_class)) +
    ggplot2::geom_col(width=0.8) +
    ggplot2::scale_fill_manual(values=c("sus vs res"="#377eb8", "res vs con"="#4daf4a", "sus vs con"="#e41a1c")) +
    ggplot2::labs(title=paste0(theme_id, " — comparison balance (", dbn, ")"), x="Region_Layer", y="Fraction of contexts") +
    ggplot2::theme_minimal(9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=60, hjust=1))
  fp <- file.path(dirs$plots_class, paste0("class_balance_", gsub("[^A-Za-z0-9]+","_", theme_id), "_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=12, height=6); message("[class-balance] Saved: ", fp)
  invisible(p)
}
invisible(lapply(split(top_themes, top_themes$db), function(dfdb) {
  dbn <- unique(dfdb$db); ths <- head(dfdb$theme, 8)
  lapply(ths, function(th) plot_class_balance(dbn, th))
}))

# 14f) Contrast strip (compact comparison matrix subset) for top themes
plot_theme_strip <- function(dbn, theme_id, max_ctx=40) {
  df <- comp_level %>% dplyr::filter(db==dbn, theme==theme_id)
  if (!nrow(df)) return(NULL)
  ctx_ord <- df %>% dplyr::count(ctx_label, sort=TRUE) %>% dplyr::slice_head(n=max_ctx) %>% dplyr::pull(ctx_label)
  df <- df %>% dplyr::filter(ctx_label %in% ctx_ord)
  df$ctx_label <- factor(df$ctx_label, levels=ctx_ord)

  p <- ggplot2::ggplot(df, ggplot2::aes(x=ctx_label, y=paste(region,layer,sep="_"), fill=meanNES)) +
    ggplot2::geom_tile(color="grey95") +
    ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b") +
    ggplot2::facet_wrap(~ comparison_class, ncol=1, scales="free_x") +
    ggplot2::labs(title=paste0(theme_id, " — contrast strip (", dbn, ")"), x="Comparison", y="Region_Layer") +
    ggplot2::theme_minimal(9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=60, hjust=1))
  fp <- file.path(dirs$plots_cm, paste0("theme_", gsub("[^A-Za-z0-9]+","_", theme_id), "_", dbn, "_strip.svg"))
  ggplot2::ggsave(fp, p, width=14, height=8); message("[strip] Saved: ", fp)
  invisible(p)
}
invisible(lapply(split(top_themes, top_themes$db), function(dfdb) {
  dbn <- unique(dfdb$db); ths <- head(dfdb$theme, 6)
  lapply(ths, function(th) plot_theme_strip(dbn, th, max_ctx=40))
}))

# 14g) Simple localization graph per DB (themes ↔ region-layer nodes)
t_thr <- 0.6
edges_df <- theme_loc %>%
  dplyr::left_join(
    theme_comp %>%
      dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L),
                       by="context") %>%
      dplyr::filter(!is.na(region), !is.na(layer), minFDR<=0.05) %>%
      dplyr::group_by(db, theme, region, layer) %>%
      dplyr::summarize(n_sig = dplyr::n(), sign_dir = sign(mean(meanNES, na.rm=TRUE)), .groups="drop"),
    by=c("db","theme","region","layer")
  ) %>%
  # guard: keep rows with numeric NES and threshold on absolute value
  dplyr::filter(!is.na(meanNES), is.finite(meanNES), abs(meanNES) >= t_thr) %>%
  # vectorized fallback for missing n_sig
  dplyr::mutate(n_sig = ifelse(is.na(n_sig), 1L, n_sig),
                src = theme,
                dst = paste(region, layer, sep="_"),
                edge_col = ifelse(meanNES > 0, "#b2182b", "#2166ac"),
                edge_w = scales::rescale(pmax(n_sig, 1L), to=c(0.4, 2)))

plot_localization_graph <- function(dbn) {
  df <- edges_df %>% dplyr::filter(db==dbn)
  if (!nrow(df)) return(NULL)

  # Build graph
  g <- igraph::graph_from_data_frame(df %>% dplyr::select(src, dst), directed=FALSE)

  # Align edge attributes to graph edges
  e_df <- igraph::as_data_frame(g, what="edges")
  key_g  <- paste(e_df$from, e_df$to, sep="||")
  key_df <- paste(df$src, df$dst, sep="||")
  idx <- match(key_g, key_df)
  mis <- which(is.na(idx))
  if (length(mis)) {
    key_rev <- paste(e_df$to[mis], e_df$from[mis], sep="||")
    idx_rev <- match(key_rev, key_df)
    idx[mis] <- idx_rev
  }
  idx[is.na(idx)] <- 1L

  igraph::E(g)$edge_col <- df$edge_col[idx]
  igraph::E(g)$edge_w   <- df$edge_w[idx]

  set.seed(1)
  p <- ggraph::ggraph(g, layout="fr") +
    ggraph::geom_edge_link(ggplot2::aes(edge_colour = edge_col,
                                        edge_width  = edge_w), alpha=0.8) +
    ggraph::geom_node_point(ggplot2::aes(shape = ifelse(grepl("_", name), "RL", "Theme")), size=3) +
    ggraph::geom_node_text(ggplot2::aes(label=name), repel=TRUE, size=3) +
    ggraph::scale_edge_colour_identity() +
    ggraph::scale_edge_width(range=c(0.4,2)) +
    ggplot2::scale_shape_manual(values=c(Theme=16, RL=15)) +
    ggplot2::theme_void() +
    ggplot2::labs(title=paste0(dbn, " — localization graph (|NES|≥", t_thr, ")"))

  fp <- file.path(dirs$plots, "Graphs", paste0("localization_graph_", dbn, ".svg"))
  dir.create(dirname(fp), recursive=TRUE, showWarnings=FALSE)
  ggplot2::ggsave(fp, p, width=10, height=8); message("[graph] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(theme_loc$db), plot_localization_graph))


# 14h) Terms network per DB (Jaccard similarity of GSEA terms)
# --- Terms network (per DB) ---
# Recompute similarity within DB to match your themes more closely
# Term network colored by theme cluster
tokenize <- function(x){
  x <- tolower(x); x <- gsub("[^a-z0-9 ]+"," ",x)
  unique(unlist(strsplit(x,"\\s+")))
}
jaccard <- function(a,b){
  ia <- tokenize(a); ib <- tokenize(b)
  if (!length(ia) || !length(ib)) return(0)
  length(intersect(ia,ib)) / length(union(ia,ib))
}

rep_terms <- gsea_df %>%
  dplyr::inner_join(dplyr::select(themes_tbl, ID, theme), by="ID") %>%
  dplyr::group_by(db, theme, ID, Description) %>%
  dplyr::summarize(freq = dplyr::n(), .groups="drop") %>%
  dplyr::group_by(db, theme) %>%
  dplyr::slice_max(order_by=freq, n=10, with_ties=FALSE) %>%
  dplyr::ungroup()

plot_terms_network <- function(dbn, sim_cut=0.25, max_nodes=250) {
  df <- rep_terms %>% dplyr::filter(db==dbn)
  if (!nrow(df)) return(NULL)
  # control size: top 12 themes
  top_themes <- df %>% dplyr::count(theme, sort=TRUE) %>% dplyr::slice_head(n=12) %>% dplyr::pull(theme)
  df <- df %>% dplyr::filter(theme %in% top_themes)
  if (nrow(df) > max_nodes) {
    per <- max(5L, floor(max_nodes/length(top_themes)))
    df <- df %>% dplyr::group_by(theme) %>% dplyr::slice_head(n=per) %>% dplyr::ungroup()
  }

  terms <- df$Description; n <- length(terms)
  if (n < 2) return(NULL)
  edges <- list()
  for (i in seq_len(n-1)) for (j in (i+1):n) {
    s <- jaccard(terms[i], terms[j]); if (s >= sim_cut) edges[[length(edges)+1]] <- c(i,j,s)
  }
  if (!length(edges)) return(NULL)
  edges <- as.data.frame(do.call(rbind, edges)); names(edges) <- c("i","j","sim")
  edges$i <- as.integer(edges$i); edges$j <- as.integer(edges$j); edges$sim <- as.numeric(edges$sim)

  verts <- tibble::tibble(name = terms, theme = df$theme, term = df$Description)
  g <- igraph::graph_from_data_frame(
    d = tibble::tibble(from = verts$name[edges$i], to = verts$name[edges$j], sim = edges$sim),
    directed = FALSE, vertices = verts
  )

  set.seed(1)
  p <- ggraph::ggraph(g, layout="fr") +
    ggraph::geom_edge_link(ggplot2::aes(width=sim), colour="grey70", alpha=0.35) +
    ggraph::geom_node_point(ggplot2::aes(color=theme), size=2) +
    ggraph::geom_node_text(ggplot2::aes(label=term), size=2.4, repel=TRUE) +
    ggplot2::scale_edge_width(range=c(0.2,1.4), guide="none") +
    ggplot2::theme_void() +
    ggplot2::labs(title=paste0(dbn, " — terms network (sim≥", sim_cut, ")"), color="Theme")
  fp <- file.path(dirs$plots_anat, paste0("terms_network_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=16, height=12); message("[terms-net] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(rep_terms$db), plot_terms_network))

suppressWarnings({
  if (!requireNamespace("uwot", quietly=TRUE)) message("Optional: install.packages('uwot') for UMAP")
})

build_theme_kw_matrix <- function(dbn, top_kw=800) {
  df <- gsea_df %>%
    dplyr::filter(db==dbn) %>%
    dplyr::inner_join(dplyr::select(themes_tbl, ID, theme), by="ID") %>%
    dplyr::distinct(theme, Description)
  if (!nrow(df)) return(NULL)

  toks <- strsplit(tolower(gsub("[^a-z0-9 ]+"," ", df$Description)), "\\s+")
  toks <- lapply(toks, function(v) v[nzchar(v) & !v %in% c("and","or","of","the","to","in","by","for","process","regulation","cellular","activity","protein","pathway","signaling")])

  env <- new.env(parent=emptyenv())
  for (i in seq_len(nrow(df))) {
    th <- df$theme[i]
    for (k in toks[[i]]) {
      key <- paste(th,k,sep="||"); env[[key]] <- (env[[key]] %||% 0L) + 1L
    }
  }
  keys <- ls(env); if (!length(keys)) return(NULL)
  sp <- strsplit(keys,"\\|\\|"); ths <- vapply(sp, `[`, character(1), 1); kws <- vapply(sp, `[`, character(1), 2)
  val <- as.integer(mget(keys, env, ifnotfound=0L))
  mat <- tibble::tibble(theme=ths, kw=kws, val=val) %>%
    dplyr::group_by(kw) %>% dplyr::summarize(dfreq=dplyr::n(), .groups="drop") %>%
    dplyr::right_join(tibble::tibble(theme=ths, kw=kws, val=val), by="kw") %>%
    dplyr::filter(dfreq>=2) %>%
    dplyr::group_by(kw) %>% dplyr::mutate(tf = val/sum(val)) %>% dplyr::ungroup() %>%
    dplyr::group_by(kw) %>% dplyr::mutate(idf = log(n_distinct(theme)/n_distinct(theme[val>0])+1e-6)) %>% dplyr::ungroup() %>%
    dplyr::mutate(tfidf = tf*idf)

  vocab <- mat %>% dplyr::group_by(kw) %>% dplyr::summarize(s=sum(tfidf,na.rm=TRUE), .groups="drop") %>%
    dplyr::arrange(dplyr::desc(s)) %>% dplyr::slice_head(n=top_kw) %>% dplyr::pull(kw)
  mat <- mat %>% dplyr::filter(kw %in% vocab)
  if (!nrow(mat)) return(NULL)

  M <- tidyr::pivot_wider(mat, names_from=kw, values_from=tfidf, values_fill=0)
  M
}

plot_theme_embedding <- function(dbn) {
  M <- build_theme_kw_matrix(dbn)
  if (is.null(M) || nrow(M) < 2) return(NULL)
  meta <- M %>% dplyr::select(theme)
  X <- as.matrix(M %>% dplyr::select(-theme))
  rownames(X) <- meta$theme

  if (requireNamespace("uwot", quietly=TRUE)) {
    set.seed(1); emb <- uwot::umap(X, n_neighbors=10, min_dist=0.2, metric="cosine")
    emb <- as.data.frame(emb); colnames(emb) <- c("U1","U2")
  } else {
    pr <- stats::prcomp(X, scale.=TRUE); emb <- as.data.frame(pr$x[,1:2]); colnames(emb) <- c("U1","U2")
  }
  emb$theme <- rownames(X)

  p <- ggplot2::ggplot(emb, ggplot2::aes(x=U1, y=U2, label=theme)) +
    ggplot2::geom_point(color="#2c7fb8", size=2, alpha=0.85) +
    ggplot2::geom_text(size=2.8, nudge_y=0.02) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title=paste0(dbn, " — theme embedding (", if (requireNamespace("uwot", quietly=TRUE)) "UMAP" else "PCA", ")"),
                  x="Dim 1", y="Dim 2")
  fp <- file.path(dirs$plots_anat, paste0("theme_embedding_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=10, height=8); message("[theme-embed] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(themes_tbl$db), plot_theme_embedding))

plot_theme_dendrogram <- function(dbn) {
  M <- build_theme_kw_matrix(dbn)
  if (is.null(M) || nrow(M) < 2) return(NULL)
  X <- as.matrix(M %>% dplyr::select(-theme))
  rownames(X) <- M$theme
  # cosine distance
  A <- X + 1e-12
  nrm <- sqrt(rowSums(A*A))
  S <- tcrossprod(A / nrm)
  D <- 1 - S
  hc <- stats::hclust(stats::as.dist(D), method="average")

  if (!requireNamespace("ggdendro", quietly=TRUE)) {
    message("Optional: install.packages('ggdendro') for dendrogram plot"); return(NULL)
  }
  dd <- ggdendro::dendro_data(hc)
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data=dd$segments, ggplot2::aes(x=x, y=y, xend=xend, yend=yend)) +
    ggplot2::geom_text(data=dd$labels, ggplot2::aes(x=x, y=y, label=label), angle=90, hjust=1, vjust=0.5, size=2.6) +
    ggplot2::theme_void() + ggplot2::labs(title=paste0(dbn, " — theme dendrogram"))
  fp <- file.path(dirs$plots_anat, paste0("theme_dendrogram_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=14, height=10); message("[theme-dend] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(themes_tbl$db), plot_theme_dendrogram))

# 14i) GO BP terms per theme (bar, top by FDR/frequency)
plot_bp_terms_per_theme <- function(top_k = 8) {
  df <- gsea_df %>%
    dplyr::filter(db == "GO_BP") %>%
    dplyr::inner_join(dplyr::select(themes_tbl, ID, theme), by="ID") %>%
    dplyr::group_by(theme, Description) %>%
    dplyr::summarize(n_hits = dplyr::n(), meanNES = mean(NES, na.rm=TRUE),
                     bestFDR = suppressWarnings(min(FDR, na.rm=TRUE)), .groups="drop") %>%
    dplyr::group_by(theme) %>%
    dplyr::arrange(bestFDR, dplyr::desc(n_hits), .by_group=TRUE) %>%
    dplyr::slice_head(n = top_k) %>%
    dplyr::ungroup()

  if (!nrow(df)) { message("[bp-terms] No GO_BP data"); return(NULL) }

  # Reorder within each facet
  df <- df %>%
    dplyr::group_by(theme) %>%
    dplyr::mutate(term_order = reorder(Description, -n_hits)) %>%
    dplyr::ungroup()

  p <- ggplot2::ggplot(df, ggplot2::aes(x = n_hits, y = term_order, fill = -log10(bestFDR))) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_viridis_c(option="C", name = "-log10(FDR)") +
    ggplot2::facet_wrap(~ theme, scales = "free_y", ncol = 3) +
    ggplot2::labs(title = "GO BP terms per cluster (top by FDR/frequency)",
                  x = "Term count within cluster", y = "GO BP term") +
    ggplot2::theme_minimal(10)
  fp <- file.path(dirs$plots_anat, "GO_BP_terms_per_theme.svg")
  ggplot2::ggsave(fp, p, width=18, height=14); message("[bp-terms] Saved: ", fp)
  invisible(p)
}
plot_bp_terms_per_theme(top_k = 8)

bp_theme_summaries <- bp_members %>%
  dplyr::group_by(theme) %>%
  dplyr::arrange(bestFDR, dplyr::desc(n_hits), .by_group=TRUE) %>%
  dplyr::summarize(
    n_terms = dplyr::n(),
    top_terms = paste(head(Description, 6), collapse=" | "),
    medianNES = stats::median(meanNES, na.rm=TRUE),
    .groups="drop"
  )
readr::write_csv(bp_theme_summaries, file.path(dirs$tables_gsea, "GO_BP_theme_summaries.csv"))

# 15) ORA summaries (separate) ------------------------------------------------
if (nrow(ora_df)) {
  ora_repl <- ora_df %>% dplyr::group_by(db, ID, Description) %>% dplyr::summarize(n_comp=dplyr::n(), n_sig=sum(Padj<=0.05, na.rm=TRUE), meanPadj=mean(Padj, na.rm=TRUE), .groups="drop") %>% dplyr::arrange(dplyr::desc(n_sig), meanPadj)
  readr::write_csv(ora_repl, file.path(dirs$tables_ora, "ORA_replication_scores.csv"))
  top_ora <- ora_repl %>% dplyr::group_by(db) %>% dplyr::slice_max(n_sig, n=25) %>% dplyr::ungroup()
  p_ora <- ggplot2::ggplot(top_ora, ggplot2::aes(x=reorder(paste(db, Description, sep="::"), n_sig), y=n_sig, fill=db)) +
    ggplot2::geom_col() + ggplot2::coord_flip() +
    ggplot2::labs(title="ORA term replication (count of significant comparisons)", x="DB::Term", y="# significant comps") +
    ggplot2::theme_minimal(base_size=10)
  ggplot2::ggsave(file.path(dirs$plots_ora, "ORA_replication_bar.svg"), p_ora, width=18, height=12)
}

# 16) Executive bundle + artifact check ---------------------------------------
readr::write_csv(top_themes,  file.path(dirs$tables_gsea, "top_themes_per_db.csv"))
readr::write_csv(theme_comp %>% dplyr::arrange(db, theme, dplyr::desc(abs(meanNES))), file.path(dirs$tables_gsea, "theme_by_comparison_sorted.csv"))

md <- c(
  "# Consolidated Summary",
  paste0("- Comparisons scanned: ", length(comp_dirs)),
  paste0("- GSEA rows: ", nrow(gsea_df)),
  "",
  "Artifacts:",
  "* Plots/Anatomy/localization_<DB>.svg",
  "* Plots/Class/localization_<DB>_by_comparison_class_labeled.svg",
  "* Plots/Conditions/localization_<DB>_tested_condition.svg",
  "* Plots/Conditions/localization_<DB>_baseline_condition.svg",
  "* Plots/ComparisonMatrices/theme_<THEME>_<DB>_comparison_matrix.svg",
  "* Tables/context_region_layer_condition_both_sides.csv",
  "* Tables/Class/context_comparison_class_map.csv",
  "* Tables/GSEA/localization_table.csv",
  "* Tables/Class/localization_by_comparison_class_table.csv",
  "* Tables/Class/localization_by_comparison_class_labeled.csv",
  "* Tables/Conditions/localization_left_condition_table.csv",
  "* Tables/Conditions/localization_right_condition_table.csv",
  "* Tables/Class/comparison_level_long.csv",
  "* Tables/QA/localization_NA_diagnostics.csv",
  "* Tables/GSEA/all_GSEA_long.csv / themes_table.csv / theme_by_comparison.csv / theme_replication_scores.csv",
  if (nrow(ora_df)) "* Tables/ORA/all_ORA_long.csv / ORA_replication_scores.csv; Plots/ORA/ORA_replication_bar.svg" else "* (no ORA detected)"
)
writeLines(md, con = file.path(out_dir, "CONSOLIDATED_SUMMARY.md"))

# Artifact existence check across subfolders
expected <- c(
  file.path(dirs$plots_anat,  paste0("localization_", unique(theme_loc$db), ".svg")),
  file.path(dirs$plots_class, paste0("localization_", unique(theme_loc_by_class_lab$db), "_by_comparison_class_labeled.svg")),
  file.path(dirs$plots_cond,  paste0("localization_", unique(theme_loc_left_cond$db), "_tested_condition.svg")),
  file.path(dirs$plots_cond,  paste0("localization_", unique(theme_loc_right_cond$db), "_baseline_condition.svg"))
)
exists_tbl <- tibble::tibble(file = unique(expected), exists = file.exists(unique(expected)))
readr::write_csv(exists_tbl, file.path(out_dir, "artifact_existence_check.csv"))
print(exists_tbl)

message("Consolidation complete. See ZZ_consolidated for outputs.")

# Optionally: append module artifacts to CONSOLIDATED_SUMMARY.md
try({
  md_mod <- c(
    "",
    "Module Artifacts:",
    "* Tables/Modules/module_GSEA_long.csv / module_ORA_long.csv",
    "* Tables/Modules/module_theme_by_comparison.csv / module_localization_table.csv",
    "* Tables/Modules/module_replication_scores.csv",
    "* Tables/Modules/module_localization_by_comparison_class.csv",
    "* Tables/Modules/module_localization_by_comparison_class_labeled.csv",
    "* Tables/Modules/module_driver_strength.csv",
    "* Tables/Modules/module_comparison_level_long.csv",
    "* Tables/Modules/module_dominant_direction.csv / module_directional_bias_pvalues.csv",
    "* Plots/Modules/module_localization_<DB>.svg",
    "* Plots/Modules/module_localization_<DB>_by_comparison_class_labeled.svg",
    "* Plots/Modules/module_<MODULE>_<DB>_comparison_matrix.svg",
    "* Plots/Modules/module_drivers_susceptibility_<DB>.svg"
  )
  cat(paste0(md_mod, collapse="\n"), file = file.path(out_dir, "CONSOLIDATED_SUMMARY.md"), append = TRUE)
}, silent = TRUE)








































# =======================
# CONSOLIDATION + SUMMARY
# =======================
# Primary: GSEA (gseGO/gseKEGG/Reactome) using NES/FDR
# Secondary: ORA (separate tables/plots)
# Includes class-labeled heatmaps (SR/RC/SC:n) and comparison-level matrices. [web:41]

# 0) Setup --------------------------------------------------------------------
pkgs <- c("dplyr","tidyr","readr","stringr","purrr","ggplot2","igraph","ggraph", "patchwork", "tibble", "clusterProfiler", "org.Mm.eg.db", "DOSE", "enrichplot", "ReactomePA")
missing_pkgs <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(missing_pkgs)) {
  message("Installing missing packages: ", paste(missing_pkgs, collapse=", "))
  tryCatch({
    install.packages(missing_pkgs, repos = c("https://cloud.r-project.org"))
  }, error = function(e) {
    warning("Primary CRAN failed: ", conditionMessage(e), " — retrying ggraph via r-universe")
    if ("ggraph" %in% missing_pkgs) install.packages("ggraph", repos = c("https://thomasp85.r-universe.dev","https://cloud.r-project.org"))
  })
}
suppressPackageStartupMessages({
  for (p in pkgs) {
    tryCatch(library(p, character.only = TRUE),
             error = function(e) warning(sprintf("Failed to load %s: %s", p, conditionMessage(e))))
  }
})

# Root paths (fixed) ----------------------------------------------------------
root_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Results/neuron-phenotypeWithinUnit"
out_dir  <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Results/ZZ_consolidated"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# New: subfolder layout for tables/plots --------------------------------------
dirs <- list(
  tables        = file.path(out_dir, "Tables"),
  tables_gsea   = file.path(out_dir, "Tables", "GSEA"),
  tables_ora    = file.path(out_dir, "Tables", "ORA"),
  tables_class  = file.path(out_dir, "Tables", "Class"),
  tables_cond   = file.path(out_dir, "Tables", "Conditions"),
  tables_qa     = file.path(out_dir, "Tables", "QA"),
  plots         = file.path(out_dir, "Plots"),
  plots_anat    = file.path(out_dir, "Plots", "Anatomy"),
  plots_class   = file.path(out_dir, "Plots", "Class"),
  plots_cond    = file.path(out_dir, "Plots", "Conditions"),
  plots_cm      = file.path(out_dir, "Plots", "ComparisonMatrices"),
  plots_ora     = file.path(out_dir, "Plots", "ORA")
)
invisible(lapply(dirs, function(d) dir.create(d, showWarnings = FALSE, recursive = TRUE)))  # recursive creation [web:73][web:68]

# --- MODULE BRANCH: paths + filename parser (after dirs creation) ------------
dirs$tables_modules <- file.path(out_dir, "Tables", "Modules")
dirs$plots_modules  <- file.path(out_dir, "Plots",  "Modules")
invisible(lapply(c(dirs$tables_modules, dirs$plots_modules), dir.create, recursive=TRUE, showWarnings=FALSE))  # create module dirs [web:73][web:68]

modules_root <- file.path(dirname(root_dir), "20_modules")
if (!dir.exists(modules_root)) warning("modules_root not found: ", modules_root)  # diag [web:73]

parse_module_fname <- function(fp) {
  bn <- basename(fp); noext <- sub("\\.csv$", "", bn, ignore.case = TRUE)
  # gsea: gseGO_BP_<module>_<LEFT>_<RIGHT>
  m1 <- regexec("^gse([A-Za-z0-9]+)_([A-Za-z0-9]+)_([^_]+)_(.+)$", noext, perl=TRUE)
  r1 <- regmatches(noext, m1)[[1]]
  if (length(r1)) {
    db1 <- r1[2]; db2 <- r1[3]; db <- if (toupper(db1)=="GO" && toupper(db2)=="BP") "GO_BP" else paste0(db1,"_",db2)
    mod <- r1[4]; rem <- r1[5]
    gpos <- gregexpr("_(CA1|CA2|CA3|DG)", rem, perl=TRUE)[[1]]
    split_pos <- if (length(gpos) && gpos[1] != -1) tail(gpos,1) else -1
    if (split_pos == -1) { toks <- strsplit(rem, "_")[[1]]; left <- paste(toks[1:(length(toks)-1)], collapse="_"); right <- toks[length(toks)] } else {
      left <- substr(rem, 1, split_pos-1); right <- substr(rem, split_pos+1, nchar(rem))
    }
    return(list(source="gsea", db=db, module=mod, context=paste0(left," vs ", right)))
  }
  # ORA_BP_<module>_<LEFT>_<RIGHT>
  m2 <- regexec("^ORA_([A-Za-z0-9]+)_([^_]+)_(.+)$", noext, perl=TRUE)
  r2 <- regmatches(noext, m2)[[1]]
  if (length(r2)) {
    db2 <- r2[2]; db <- if (toupper(db2)=="BP") "GO_BP" else db2
    mod <- r2[3]; rem <- r2[4]
    gpos <- gregexpr("_(CA1|CA2|CA3|DG)", rem, perl=TRUE)[[1]]
    split_pos <- if (length(gpos) && gpos[1] != -1) tail(gpos,1) else -1
    if (split_pos == -1) { toks <- strsplit(rem, "_")[[1]]; left <- paste(toks[1:(length(toks)-1)], collapse="_"); right <- toks[length(toks)] } else {
      left <- substr(rem, 1, split_pos-1); right <- substr(rem, split_pos+1, nchar(rem))
    }
    return(list(source="ora_deg", db=db, module=mod, context=paste0(left," vs ", right)))
  }
  # ORA_fullModule_BP_<module>_<LEFT>_<RIGHT>
  m3 <- regexec("^ORA_fullModule_([A-Za-z0-9]+)_([^_]+)_(.+)$", noext, perl=TRUE)
  r3 <- regmatches(noext, m3)[[1]]
  if (length(r3)) {
    db2 <- r3[2]; db <- if (toupper(db2)=="BP") "GO_BP" else db2
    mod <- r3[3]; rem <- r3[4]
    gpos <- gregexpr("_(CA1|CA2|CA3|DG)", rem, perl=TRUE)[[1]]
    split_pos <- if (length(gpos) && gpos[1] != -1) tail(gpos,1) else -1
    if (split_pos == -1) { toks <- strsplit(rem, "_")[[1]]; left <- paste(toks[1:(length(toks)-1)], collapse="_"); right <- toks[length(toks)] } else {
      left <- substr(rem, 1, split_pos-1); right <- substr(rem, split_pos+1, nchar(rem))
    }
    return(list(source="ora_full", db=db, module=mod, context=paste0(left," vs ", right)))
  }
  NULL
}  # regex parser for provided file patterns [web:60][web:59]

# 1) Helpers ------------------------------------------------------------------
read_if <- function(fp) {
  if (!file.exists(fp)) return(NULL)
  tryCatch(readr::read_csv(fp, guess_max = 10000, progress = FALSE, show_col_types = FALSE),
           error = function(e) { warning("Failed reading: ", fp, " — ", conditionMessage(e)); NULL })
}
read_context <- function(comp_dir) {
  readme <- file.path(comp_dir, "README.txt")
  if (file.exists(readme)) {
    rl <- readLines(readme, warn = FALSE)
    ctx <- rl[grepl("^Context: ", rl)]
    if (length(ctx)) return(sub("^Context:\\s*", "", ctx[1]))
  }
  basename(comp_dir)
}
`%||%` <- function(a,b) if (!is.null(a) && length(a)>0 && !is.na(a)) a else b
tokenize <- function(x){ x <- tolower(x); x <- gsub("[^a-z0-9]+"," ",x); unique(unlist(strsplit(x,"\\s+"))) }
jaccard  <- function(a,b){ ia<-tokenize(a); ib<-tokenize(b); if(!length(ia)||!length(ib)) return(0); length(intersect(ia,ib))/length(union(ia,ib)) }

# 2) GSEA collectors -----------------------------------------------------------
parse_gsea_tbl <- function(fp, db) {
  df <- read_if(fp); if (is.null(df) || !nrow(df)) return(NULL)
  nm <- names(df); has <- function(x) any(x %in% nm)
  if (!has("ID") || !has("Description")) return(NULL)
  if (!has("NES")) return(NULL)
  fdr_col <- c("p.adjust","qvalue","padj","FDR","adj.P.Val")[c("p.adjust","qvalue","padj","FDR","adj.P.Val") %in% nm][1]
  if (is.na(fdr_col)) return(NULL)
  out <- tibble::tibble(
    db = rep(db, nrow(df)),
    ID = df[["ID"]],
    Description = df[["Description"]],
    NES = suppressWarnings(as.numeric(df[["NES"]])),
    FDR = suppressWarnings(as.numeric(df[[fdr_col]]))
  ) %>% dplyr::filter(!is.na(ID), !is.na(NES), !is.na(FDR))
  if (!nrow(out)) return(NULL)
  out
}
empty_gsea <- function(){ tibble::tibble(db=character(), ID=character(), Description=character(), NES=double(), FDR=double()) }

collect_gsea_one <- function(comp_dir) {
  paths <- list(
    GO_BP    = file.path(comp_dir, "10_global","GO_BP","Tables"),
    KEGG     = file.path(comp_dir, "10_global","KEGG","Tables"),
    Reactome = file.path(comp_dir, "10_global","Reactome","Tables")
  )
  if (!any(dir.exists(unlist(paths)))) return(empty_gsea())
  gsego   <- if (dir.exists(paths$GO_BP))    list.files(paths$GO_BP,    pattern="^gseGO_BP_.*\\.csv$",      full.names=TRUE) else character(0)
  gsekegg <- if (dir.exists(paths$KEGG))     list.files(paths$KEGG,     pattern="^gseKEGG_.*\\.csv$",       full.names=TRUE) else character(0)
  rs      <- if (dir.exists(paths$Reactome)) list.files(paths$Reactome, pattern="^ReactomeGSEA_.*\\.csv$",  full.names=TRUE) else character(0)
  if (length(gsego)+length(gsekegg)+length(rs) == 0) return(empty_gsea())
  parts <- list(
    dplyr::bind_rows(lapply(gsego,   parse_gsea_tbl, db="GO_BP")),
    dplyr::bind_rows(lapply(gsekegg, parse_gsea_tbl, db="KEGG")),
    dplyr::bind_rows(lapply(rs,      parse_gsea_tbl, db="Reactome"))
  )
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (!length(parts)) return(empty_gsea())
  dplyr::bind_rows(parts)
}

# 3) ORA collector -------------------------------------------------------------
parse_ora_tbl <- function(fp, db) {
  df <- read_if(fp); if (is.null(df) || !nrow(df)) return(NULL)
  nm <- names(df); if (!("ID" %in% nm && "Description" %in% nm)) return(NULL)
  fdr_col <- c("p.adjust","qvalue","padj","FDR","adj.P.Val")[c("p.adjust","qvalue","padj","FDR","adj.P.Val") %in% nm][1]
  if (is.na(fdr_col)) return(NULL)
  out <- tibble::tibble(db=rep(db, nrow(df)), ID=df[["ID"]], Description=df[["Description"]], Padj=suppressWarnings(as.numeric(df[[fdr_col]]))) %>%
    dplyr::filter(!is.na(ID), !is.na(Padj))
  if (!nrow(out)) return(NULL)
  out
}
empty_ora <- function(){ tibble::tibble(db=character(), ID=character(), Description=character(), Padj=double()) }

collect_ora_one <- function(comp_dir) {
  paths <- list(
    GO_BP    = file.path(comp_dir, "10_global","GO_BP","Tables"),
    KEGG     = file.path(comp_dir, "10_global","KEGG","Tables"),
    Reactome = file.path(comp_dir, "10_global","Reactome","Tables")
  )
  if (!any(dir.exists(unlist(paths)))) return(empty_ora())
  files <- c(
    if (dir.exists(paths$GO_BP))    list.files(paths$GO_BP,    pattern="^ORA_.*\\.csv$",         full.names=TRUE) else character(0),
    if (dir.exists(paths$KEGG))     list.files(paths$KEGG,     pattern="^ORA_.*\\.csv$",         full.names=TRUE) else character(0),
    if (dir.exists(paths$Reactome)) list.files(paths$Reactome, pattern="^ReactomeORA_.*\\.csv$", full.names=TRUE) else character(0)
  )
  if (!length(files)) return(empty_ora())
  binders <- lapply(files, function(fp) {
    db <- if (grepl("/GO_BP/", fp)) "GO_BP" else if (grepl("/KEGG/", fp)) "KEGG" else if (grepl("/Reactome/", fp)) "Reactome" else "Unknown"
    parse_ora_tbl(fp, db = db)
  })
  binders <- binders[!vapply(binders, is.null, logical(1))]
  if (!length(binders)) return(empty_ora())
  dplyr::bind_rows(binders)
}

# 4) Scan comparison folders ---------------------------------------------------
comp_dirs <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
comp_dirs <- comp_dirs[dir.exists(file.path(comp_dirs, "10_global"))]
message(sprintf("Found %d comparison folders under neuron-phenotypeWithinUnit.", length(comp_dirs)))
if (!length(comp_dirs)) stop("No comparison folders with '10_global' found under the fixed root_dir.", call. = FALSE)

# 5) GSEA assembly -------------------------------------------------------------
gsea_list <- lapply(comp_dirs, function(cd) tibble::tibble(comp_dir = cd, context = read_context(cd), data = list(collect_gsea_one(cd))))
gsea_df <- purrr::map_dfr(gsea_list, function(x) { dx <- x$data[[1]]; if (!is.null(dx) && nrow(dx) > 0) cbind(x[1:2], dx) else NULL })
if (!nrow(gsea_df)) stop("No valid GSEA rows found across comparisons.", call. = FALSE)
req <- c("comp_dir","context","db","ID","Description","NES","FDR")
missing_cols <- setdiff(req, names(gsea_df))
if (length(missing_cols)) stop(sprintf("gsea_df missing columns: %s", paste(missing_cols, collapse=", ")), call. = FALSE)
gsea_df <- gsea_df %>% dplyr::mutate(comp_dir=as.character(comp_dir), context=as.character(context), db=as.character(db), ID=as.character(ID), Description=as.character(Description), NES=as.numeric(NES), FDR=as.numeric(FDR))
readr::write_csv(gsea_df, file.path(dirs$tables_gsea, "all_GSEA_long.csv"))
message(sprintf("GSEA rows: %d", nrow(gsea_df)))

# 6) ORA assembly --------------------------------------------------------------
ora_list <- lapply(comp_dirs, function(cd) tibble::tibble(comp_dir = cd, context = read_context(cd), data = list(collect_ora_one(cd))))
ora_df <- purrr::map_dfr(ora_list, function(x) { dx <- x$data[[1]]; if (!is.null(dx) && nrow(dx) > 0) cbind(x[1:2], dx) else NULL })
if (nrow(ora_df)) readr::write_csv(ora_df, file.path(dirs$tables_ora, "all_ORA_long.csv"))

# --- 6b) Module-level collectors and assembly --------------------------------
module_dirs <- if (dir.exists(modules_root)) list.dirs(modules_root, full.names=TRUE, recursive=FALSE) else character(0)
module_dirs <- module_dirs[grepl("/Module_", gsub("\\\\","/", module_dirs))]
message(sprintf("Found %d Module_* folders.", length(module_dirs)))

list_module_files <- function(module_dir) {
  tbl_dir <- file.path(module_dir, "Tables")
  if (!dir.exists(tbl_dir)) return(character(0))
  subdirs <- c("gsea","ora_deg","ora_full")
  unlist(lapply(subdirs, function(s) {
    d <- file.path(tbl_dir, s)
    if (dir.exists(d)) list.files(d, pattern="\\.csv$", full.names=TRUE) else character(0)
  }))
}

parse_module_csv <- function(fp) {
  meta <- parse_module_fname(fp); if (is.null(meta)) return(NULL)
  df <- read_if(fp); if (is.null(df) || !nrow(df)) return(NULL)
  nm <- names(df)
  if (identical(meta$source, "gsea")) {
    fdr_col <- c("p.adjust","qvalue","padj","FDR","adj.P.Val")[c("p.adjust","qvalue","padj","FDR","adj.P.Val") %in% nm][1]
    if (!all(c("ID","Description","NES") %in% nm) || is.na(fdr_col)) return(NULL)
    tibble::tibble(
      module = meta$module, db = meta$db, context = meta$context,
      ID = as.character(df$ID), Description = as.character(df$Description),
      NES = suppressWarnings(as.numeric(df$NES)), FDR = suppressWarnings(as.numeric(df[[fdr_col]]))
    ) %>% dplyr::filter(!is.na(ID), is.finite(NES), is.finite(FDR))
  } else {
    padj_col <- c("p.adjust","qvalue","padj","FDR","adj.P.Val")[c("p.adjust","qvalue","padj","FDR","adj.P.Val") %in% nm][1]
    if (!all(c("ID","Description") %in% nm) || is.na(padj_col)) return(NULL)
    tibble::tibble(
      module = meta$module, db = meta$db, context = meta$context,
      ID = as.character(df$ID), Description = as.character(df$Description),
      Padj = suppressWarnings(as.numeric(df[[padj_col]])), source = meta$source
    ) %>% dplyr::filter(!is.na(ID), is.finite(Padj))
  }
}

module_files <- unlist(lapply(module_dirs, list_module_files))
message(sprintf("Module CSV files detected: %d", length(module_files)))
mod_parsed <- lapply(module_files, parse_module_csv)
mod_parsed <- mod_parsed[!vapply(mod_parsed, is.null, logical(1))]

module_gsea_df <- purrr::map_dfr(mod_parsed, ~{ if (all(c("NES","FDR") %in% names(.x))) .x else NULL })
module_ora_df  <- purrr::map_dfr(mod_parsed, ~{ if ("Padj" %in% names(.x)) .x else NULL })

if (nrow(module_gsea_df)) {
  readr::write_csv(module_gsea_df, file.path(dirs$tables_modules, "module_GSEA_long.csv"))
  message(sprintf("Module GSEA rows: %d", nrow(module_gsea_df)))
} else message("No module GSEA rows.")
if (nrow(module_ora_df)) {
  readr::write_csv(module_ora_df, file.path(dirs$tables_modules, "module_ORA_long.csv"))
  message(sprintf("Module ORA rows: %d", nrow(module_ora_df)))
} else message("No module ORA rows.")

# 7) Theme building (GSEA only) -----------------------------------------------
build_themes_for_db <- function(df_db, top_n_terms = 200, sim_cut = 0.25) {
  if (is.null(df_db) || !nrow(df_db)) return(NULL)
  df_db <- dplyr::mutate(df_db, sig = FDR <= 0.05)
  agg <- df_db %>%
    dplyr::group_by(ID, Description, db) %>%
    dplyr::summarize(n_sig = sum(sig, na.rm = TRUE), n = dplyr::n(), meanNES = mean(NES, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(n_sig), dplyr::desc(abs(meanNES)))
  sel <- head(agg, top_n_terms)
  if (nrow(sel) < 2) return(dplyr::mutate(sel, theme = paste0(unique(df_db$db), "_", dplyr::row_number())))
  ids <- sel$ID; descs <- sel$Description; n <- length(ids)
  sim <- matrix(0, n, n); rownames(sim) <- colnames(sim) <- ids
  for (i in seq_len(n)) for (j in i:n) { s <- if (i==j) 1 else jaccard(descs[i], descs[j]); sim[i,j] <- s; sim[j,i] <- s }
  edge_df <- as.data.frame(as.table(sim), stringsAsFactors = FALSE) |>
    dplyr::rename(ID1 = Var1, ID2 = Var2, sim = Freq) |>
    dplyr::filter(ID1 != ID2, sim >= sim_cut)
  g <- igraph::make_empty_graph(directed = FALSE) |> igraph::add_vertices(n, name = ids)
  if (nrow(edge_df) > 0) {
    edges_mat <- unique(t(apply(edge_df[,c("ID1","ID2")], 1, function(x) sort(as.character(x)))))
    g <- igraph::add_edges(g, as.vector(t(edges_mat)))
  }
  if (igraph::vcount(g) == 0) { sel$theme <- paste0(unique(df_db$db), "_", seq_len(nrow(sel))); return(sel) }
  cl <- igraph::cluster_louvain(g)
  memb <- igraph::membership(cl)
  theme_map <- tibble::tibble(ID = names(memb), theme = paste0(unique(df_db$db), "_C", as.integer(memb)))
  out <- dplyr::left_join(sel, theme_map, by = "ID")
  if (any(is.na(out$theme))) out$theme[is.na(out$theme)] <- paste0(unique(df_db$db), "_", seq_len(sum(is.na(out$theme))))
  out
}

gsea_core <- gsea_df %>% dplyr::select(db, ID, Description, NES, FDR) %>% dplyr::mutate(db=as.character(db), ID=as.character(ID), Description=as.character(Description), NES=as.numeric(NES), FDR=as.numeric(FDR))
req_theme <- c("db","ID","Description","NES","FDR"); miss_theme <- setdiff(req_theme, names(gsea_core)); if (length(miss_theme)) stop(sprintf("gsea_core missing: %s", paste(miss_theme, collapse=", ")), call.=FALSE)
split_dbs <- split(gsea_core, gsea_core$db)
themes_db <- lapply(split_dbs, build_themes_for_db)
themes_db <- themes_db[!vapply(themes_db, is.null, logical(1))]
themes_tbl <- if (length(themes_db)) dplyr::bind_rows(themes_db) else tibble::tibble()
readr::write_csv(themes_tbl, file.path(dirs$tables_gsea, "themes_table.csv"))

# BP membership per theme (cluster), per DB
bp_members <- gsea_df %>%
  dplyr::filter(db == "GO_BP") %>%
  dplyr::inner_join(dplyr::select(themes_tbl, ID, theme), by="ID") %>%
  dplyr::group_by(db, theme, ID, Description) %>%
  dplyr::summarize(
    n_hits = dplyr::n(),
    meanNES = mean(NES, na.rm=TRUE),
    bestFDR = suppressWarnings(min(FDR, na.rm=TRUE)),
    .groups="drop"
  ) %>%
  dplyr::arrange(db, theme, dplyr::desc(n_hits), bestFDR)
readr::write_csv(bp_members, file.path(dirs$tables_gsea, "GO_BP_terms_by_theme_long.csv"))

bp_split <- split(bp_members, interaction(bp_members$db, bp_members$theme, drop=TRUE))
invisible(purrr::iwalk(bp_split, function(df, nm){
  parts <- strsplit(nm, "\\.")[[1]]
  dbn <- parts[1]; th <- parts[2]
  fp <- file.path(dirs$tables_gsea, paste0("GO_BP_terms_", th, ".csv"))
  readr::write_csv(df, fp)
}))

# 8) Aggregation per comparison x theme (GSEA) --------------------------------
assign_themes <- if (nrow(themes_tbl)) dplyr::inner_join(gsea_df, dplyr::select(themes_tbl, ID, theme), by = "ID") else tibble::tibble()
if (!nrow(assign_themes)) stop("No themes constructed (themes_tbl empty). Consider lowering sim_cut or increasing top_n_terms.", call. = FALSE)

theme_comp <- assign_themes %>%
  dplyr::group_by(db, theme, comp_dir, context) %>%
  dplyr::summarize(n_terms=dplyr::n(), meanNES=mean(NES, na.rm=TRUE), minFDR=suppressWarnings(min(FDR, na.rm=TRUE)), frac_sig=mean(FDR<=0.05, na.rm=TRUE), .groups="drop")
readr::write_csv(theme_comp, file.path(dirs$tables_gsea, "theme_by_comparison.csv"))

theme_comp$NES_masked <- ifelse(theme_comp$minFDR <= 0.05, theme_comp$meanNES, NA_real_)

# 9) Robust parsing — both sides (fixed) --------------------------------------
parse_both <- function(ctx) {
  parts <- unlist(strsplit(ctx, "(?:\\s+vs\\s+|_vs_)", perl=TRUE))
  parts <- trimws(parts)
  left  <- if (length(parts) >= 1) parts[1] else NA_character_
  right <- if (length(parts) >= 2) parts[2] else NA_character_
  region_pat <- "(DG|CA1|CA2|CA3)"
  base_layer <- "(mo|po|sr|sp|so|slm|sg)"
  cond_pat   <- "(res|sus|con)"
  parse_side <- function(side) {
    if (is.na(side) || side == "") return(list(region=NA_character_, layer=NA_character_, condition=NA_character_))
    m <- regexec(paste0("^", region_pat, "[_-]?", base_layer, "(", cond_pat, ")?$"), side, perl=TRUE)
    r <- regmatches(side, m)[[1]]
    if (length(r) == 0) {
      m2 <- regexec(paste0("^", region_pat, "[_-]?", base_layer, "(", cond_pat, ")?"), side, perl=TRUE)
      r  <- regmatches(side, m2)[[1]]
    }
    region <- NA_character_; layer <- NA_character_; condition <- NA_character_
    if (length(r) >= 3) {
      region    <- r[2]
      layer_raw <- tolower(r[3]); layer <- layer_raw
      if (length(r) >= 5 && nzchar(r[5])) condition <- tolower(r[5])
    } else {
      toks <- unlist(strsplit(side, "[_-]+"))
      region <- if (length(toks) >= 1 && grepl(paste0("^", region_pat, "$"), toks[1])) toks[1] else NA_character_
      layer  <- if (length(toks) >= 2 && grepl(paste0("^", base_layer, "$"), tolower(toks[2]))) tolower(toks[2]) else NA_character_
      condition <- if (length(toks) >= 3 && grepl(paste0("^", cond_pat, "$"), tolower(toks[3]))) tolower(toks[3]) else NA_character_
    }
    list(region=region, layer=layer, condition=condition)
  }
  L <- parse_side(left)
  R <- parse_side(right)
  list(L=L, R=R)
}

rl_both <- unique(theme_comp$context) %>%
  purrr::map_df(~{
    pr <- parse_both(.x)
    tibble::tibble(
      context   = .x,
      region_L  = pr$L$region, layer_L  = pr$L$layer, condition_L  = pr$L$condition,
      region_R  = pr$R$region, layer_R  = pr$R$layer, condition_R  = pr$R$condition
    )
  })
readr::write_csv(rl_both, file.path(dirs$tables, "context_region_layer_condition_both_sides.csv"))

# 10) Comparison classes (from parsed conditions) -----------------------------
comp_class_map <- rl_both %>%
  dplyr::transmute(
    context,
    cond_left  = condition_L,
    cond_right = condition_R,
    comparison_class = dplyr::case_when(
      cond_left == "sus" & cond_right == "res" ~ "sus vs res",
      cond_left == "res" & cond_right == "con" ~ "res vs con",
      cond_left == "sus" & cond_right == "con" ~ "sus vs con",
      TRUE ~ NA_character_
    )
  )
readr::write_csv(comp_class_map, file.path(dirs$tables_class, "context_comparison_class_map.csv"))

# --- 10b) Module context + class maps ----------------------------------------
if (exists("module_gsea_df") && nrow(module_gsea_df)) {
  mod_ctx <- unique(module_gsea_df$context)
} else if (exists("module_ora_df") && nrow(module_ora_df)) {
  mod_ctx <- unique(module_ora_df$context)
} else mod_ctx <- character(0)

if (length(mod_ctx)) {
  mod_rl_both <- tibble::tibble(context = mod_ctx) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(parsed = list(parse_both(context))) %>%
    dplyr::mutate(
      region_L  = parsed$L$region, layer_L  = parsed$L$layer, condition_L  = parsed$L$condition,
      region_R  = parsed$R$region, layer_R  = parsed$R$layer, condition_R  = parsed$R$condition
    ) %>% dplyr::ungroup() %>% dplyr::select(-parsed)
  readr::write_csv(mod_rl_both, file.path(dirs$tables_modules, "module_context_region_layer_condition_both_sides.csv"))
  mod_comp_class_map <- mod_rl_both %>%
    dplyr::transmute(
      context,
      cond_left  = condition_L,
      cond_right = condition_R,
      comparison_class = dplyr::case_when(
        cond_left == "sus" & cond_right == "res" ~ "sus vs res",
        cond_left == "res" & cond_right == "con" ~ "res vs con",
        cond_left == "sus" & cond_right == "con" ~ "sus vs con",
        TRUE ~ NA_character_
      )
    )
  readr::write_csv(mod_comp_class_map, file.path(dirs$tables_modules, "module_context_comparison_class_map.csv"))
} else {
  mod_rl_both <- tibble::tibble()
  mod_comp_class_map <- tibble::tibble()
}

# 11) Localization tables ------------------------------------------------------
theme_loc <- theme_comp %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
  dplyr::filter(!is.na(region), !is.na(layer)) %>%
  dplyr::group_by(db, theme, region, layer) %>%
  dplyr::summarize(meanNES=mean(meanNES, na.rm=TRUE), frac_sig=mean(minFDR<=0.05, na.rm=TRUE), .groups="drop")
readr::write_csv(theme_loc, file.path(dirs$tables_gsea, "localization_table.csv"))

theme_loc_by_class <- theme_comp %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
  dplyr::left_join(comp_class_map %>% dplyr::select(context, comparison_class), by="context") %>%
  dplyr::filter(!is.na(region), !is.na(layer), !is.na(comparison_class)) %>%
  dplyr::group_by(db, theme, comparison_class, region, layer) %>%
  dplyr::summarize(meanNES = mean(meanNES, na.rm=TRUE), n_comp = dplyr::n_distinct(context), .groups="drop")
readr::write_csv(theme_loc_by_class, file.path(dirs$tables_class, "localization_by_comparison_class_table.csv"))

class_labels_df <- theme_comp %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
  dplyr::left_join(comp_class_map, by="context") %>%
  dplyr::filter(!is.na(region), !is.na(layer), !is.na(comparison_class)) %>%
  dplyr::group_by(db, theme, region, layer, comparison_class) %>%
  dplyr::summarize(n_comp = dplyr::n_distinct(context), .groups="drop") %>%
  dplyr::mutate(code = dplyr::case_when(
    comparison_class=="sus vs res" ~ "SR",
    comparison_class=="res vs con" ~ "RC",
    comparison_class=="sus vs con" ~ "SC",
    TRUE ~ ""
  ),
  label = paste0(code, ":", n_comp))
theme_loc_by_class_lab <- theme_loc_by_class %>%
  dplyr::left_join(class_labels_df %>% dplyr::select(db, theme, region, layer, comparison_class, label),
                   by=c("db","theme","region","layer","comparison_class"))
readr::write_csv(theme_loc_by_class_lab, file.path(dirs$tables_class, "localization_by_comparison_class_labeled.csv"))

# Condition-centric tables
theme_loc_left_cond <- theme_comp %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L, condition_L), by="context") %>%
  dplyr::filter(!is.na(region), !is.na(layer), !is.na(condition_L)) %>%
  dplyr::group_by(db, theme, region, layer, condition_L) %>%
  dplyr::summarize(meanNES = mean(meanNES, na.rm=TRUE), n_comp = dplyr::n_distinct(context), .groups="drop") %>%
  dplyr::rename(condition = condition_L)
readr::write_csv(theme_loc_left_cond, file.path(dirs$tables_cond, "localization_left_condition_table.csv"))

theme_loc_right_cond <- theme_comp %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_R, layer=layer_R, condition_R), by="context") %>%
  dplyr::filter(!is.na(region), !is.na(layer), !is.na(condition_R)) %>%
  dplyr::group_by(db, theme, region, layer, condition_R) %>%
  dplyr::summarize(meanNES = mean(meanNES, na.rm=TRUE), n_comp = dplyr::n_distinct(context), .groups="drop") %>%
  dplyr::rename(condition = condition_R)
readr::write_csv(theme_loc_right_cond, file.path(dirs$tables_cond, "localization_right_condition_table.csv"))

# Comparison-level long
comp_level <- theme_comp %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
  dplyr::left_join(comp_class_map %>% dplyr::select(context, comparison_class), by="context") %>%
  dplyr::filter(!is.na(region), !is.na(layer), !is.na(comparison_class)) %>%
  dplyr::mutate(ctx_label = paste0(comparison_class, "::", context))
readr::write_csv(comp_level, file.path(dirs$tables_class, "comparison_level_long.csv"))

# 11b) Directionality and bias summaries (tables)
dom_dir_tbl <- theme_comp %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
  dplyr::filter(!is.na(region), !is.na(layer), minFDR<=0.05) %>%
  dplyr::group_by(db, theme, region, layer) %>%
  dplyr::summarize(medNES = stats::median(meanNES, na.rm=TRUE),
                   n_sig  = dplyr::n(), .groups="drop") %>%
  dplyr::mutate(dir = dplyr::case_when(medNES > 0 ~ "up", medNES < 0 ~ "down", TRUE ~ "flat"))
readr::write_csv(dom_dir_tbl, file.path(dirs$tables_gsea, "dominant_direction_table.csv"))

symmetry_test <- theme_comp %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
  dplyr::filter(!is.na(region), !is.na(layer)) %>%
  dplyr::group_by(db, region, layer) %>%
  dplyr::summarize(p_sym = tryCatch({
      x <- meanNES; x <- x[is.finite(x)]
      if (length(x) >= 6) stats::wilcox.test(x, mu=0, exact=FALSE)$p.value else NA_real_
    }, error=function(e) NA_real_), .groups="drop")
readr::write_csv(symmetry_test, file.path(dirs$tables_qa, "directional_bias_pvalues.csv"))

# 11c) Susceptible drivers across classes (safe column names)
sr_sc_strength <- theme_comp %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region = region_L, layer = layer_L), by = "context") %>%
  dplyr::left_join(comp_class_map %>% dplyr::select(context, comparison_class), by = "context") %>%
  dplyr::filter(!is.na(region), !is.na(layer), !is.na(comparison_class)) %>%
  dplyr::group_by(db, theme, region, layer, comparison_class) %>%
  dplyr::summarize(meanNES_class = mean(meanNES, na.rm = TRUE),
                   n_ctx = dplyr::n_distinct(context), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = comparison_class,
    values_from = c(meanNES_class, n_ctx),
    values_fill = 0
  ) %>%
  dplyr::rename(
    meanNES_class_sus_vs_res = `meanNES_class_sus vs res`,
    meanNES_class_sus_vs_con = `meanNES_class_sus vs con`,
    n_ctx_sus_vs_res         = `n_ctx_sus vs res`,
    n_ctx_sus_vs_con         = `n_ctx_sus vs con`
  ) %>%
  dplyr::mutate(
    SR_strength = meanNES_class_sus_vs_res,
    SC_strength = meanNES_class_sus_vs_con,
    driver_class = dplyr::case_when(
      abs(SR_strength) > abs(SC_strength) ~ "sus vs res",
      abs(SC_strength) > abs(SR_strength) ~ "sus vs con",
      TRUE ~ "tie"
    ),
    driver_gap = abs(SR_strength) - abs(SC_strength)
  )
readr::write_csv(sr_sc_strength, file.path(dirs$tables_class, "susceptibility_driver_strength.csv"))

# 12) NA diagnostics -----------------------------------------------------------
cov_tbl <- assign_themes %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by = "context") %>%
  dplyr::filter(!is.na(region), !is.na(layer)) %>%
  dplyr::group_by(db, theme, region, layer) %>%
  dplyr::summarize(n_comp_cover=dplyr::n_distinct(comp_dir), n_sig_cover=sum(FDR<=0.05, na.rm=TRUE), minFDR_all=suppressWarnings(min(FDR, na.rm=TRUE)), .groups="drop")
na_diag <- theme_loc %>%
  dplyr::full_join(cov_tbl, by = c("db","theme","region","layer")) %>%
  dplyr::mutate(cause = dplyr::case_when(
    is.na(meanNES) & (is.na(n_comp_cover) | n_comp_cover == 0) ~ "no_coverage",
    is.na(meanNES) & (n_comp_cover > 0) & (is.na(n_sig_cover) | n_sig_cover == 0) ~ "masked_nonsignificant",
    is.na(meanNES) ~ "other_missing",
    TRUE ~ "has_value"
  ))
readr::write_csv(na_diag, file.path(dirs$tables_qa, "localization_NA_diagnostics.csv"))

# 13) Replication summary ------------------------------------------------------
theme_repl <- theme_comp %>% dplyr::mutate(sign=sign(meanNES), sig=minFDR<=0.05) %>%
  dplyr::group_by(db, theme) %>% dplyr::summarize(n_comp=dplyr::n(), n_sig=sum(sig, na.rm=TRUE), n_pos=sum(sig & sign>0, na.rm=TRUE), n_neg=sum(sig & sign<0, na.rm=TRUE), dom_dir=ifelse(n_pos>=n_neg,"up","down"), repl_score=pmax(n_pos,n_neg), .groups="drop") %>%
  dplyr::arrange(dplyr::desc(repl_score))
readr::write_csv(theme_repl, file.path(dirs$tables_gsea, "theme_replication_scores.csv"))

top_themes <- theme_repl %>% dplyr::group_by(db) %>% dplyr::slice_max(order_by=repl_score, n=12) %>% dplyr::ungroup()

# 11d) Plots: theme localization
# Driver map: which class (SR vs SC) dominates per theme and region-layer
plot_driver_map <- function(dbn) {
  df <- sr_sc_strength %>%
    dplyr::filter(db == dbn, theme %in% (top_themes %>% dplyr::filter(db==dbn) %>% dplyr::pull(theme)))
  if (!nrow(df)) return(NULL)

  # Compact label: nSR/nSC contributing contexts (uses safe column names)
  df <- df %>%
    dplyr::mutate(n_lab = paste0(n_ctx_sus_vs_res, "/", n_ctx_sus_vs_con))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = layer, y = region, fill = driver_class, alpha = abs(driver_gap))) +
    ggplot2::geom_tile(color="white") +
    ggplot2::geom_text(ggplot2::aes(label = n_lab), size=3) +
    ggplot2::scale_fill_manual(values = c("sus vs res" = "#377eb8", "sus vs con" = "#e41a1c", "tie" = "grey80")) +
    ggplot2::scale_alpha(range = c(0.3, 1), guide = "none") +
    ggplot2::facet_wrap(~ theme, ncol = 4) +
    ggplot2::labs(title = paste0(dbn, " — drivers of susceptibility (SR vs SC)"),
                  x = "Layer", y = "Region", fill = "Class") +
    ggplot2::theme_minimal(10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))

  fp <- file.path(dirs$plots_class, paste0("drivers_susceptibility_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=22, height=12); message("[drivers-map] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(sr_sc_strength$db), plot_driver_map))

# 11e) Plots: theme driver volcano
theme_driver_volcano <- sr_sc_strength %>%
  dplyr::group_by(db, theme) %>%
  dplyr::summarize(delta = mean(SR_strength - SC_strength, na.rm=TRUE),
                   abs_gap = mean(abs(SR_strength) - abs(SC_strength), na.rm=TRUE),
                   dom_class = dplyr::case_when(
                     abs_gap > 0 ~ "sus vs res",
                     abs_gap < 0 ~ "sus vs con",
                     TRUE ~ "tie"
                   ),
                   .groups="drop") %>%
  dplyr::left_join(theme_repl %>% dplyr::select(db, theme, repl_score, dom_dir), by=c("db","theme"))

plot_theme_driver_volcano <- function(dbn) {
  df <- theme_driver_volcano %>% dplyr::filter(db==dbn)
  if (!nrow(df)) return(NULL)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = delta, y = repl_score, color = dom_class, shape = dom_dir)) +
    ggplot2::geom_hline(yintercept = median(df$repl_score, na.rm=TRUE), linetype="dashed", color="grey60") +
    ggplot2::geom_vline(xintercept = 0, linetype="dotted", color="grey60") +
    ggplot2::geom_point(size=3, alpha=0.9) +
    ggplot2::scale_color_manual(values=c("sus vs res"="#377eb8","sus vs con"="#e41a1c","tie"="grey60")) +
    ggplot2::labs(title=paste0(dbn, " — theme-level SR vs SC dominance"),
                  x="SR_strength − SC_strength (mean NES)", y="Replication score",
                  color="Dominant class", shape="Dir") +
    ggplot2::theme_minimal(11)
  fp <- file.path(dirs$plots_class, paste0("theme_driver_volcano_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=10, height=7); message("[driver-volcano] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(theme_driver_volcano$db), plot_theme_driver_volcano))

# --- MODULE ANALYSES (parallel to global; place after global sections 11–13 or right after 10b) ---
if (exists("module_gsea_df") && nrow(module_gsea_df) && nrow(mod_rl_both)) {
  # Aggregation per module x comparison (region-layer on left side)
  mod_theme_comp <- module_gsea_df %>%
    dplyr::left_join(mod_rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
    dplyr::filter(!is.na(region), !is.na(layer)) %>%
    dplyr::group_by(db, module, context, region, layer) %>%
    dplyr::summarize(meanNES = mean(NES, na.rm=TRUE),
                     minFDR  = suppressWarnings(min(FDR, na.rm=TRUE)),
                     .groups="drop")
  readr::write_csv(mod_theme_comp, file.path(dirs$tables_modules, "module_theme_by_comparison.csv"))

  # Localization table per module
  module_loc <- mod_theme_comp %>%
    dplyr::group_by(db, module, region, layer) %>%
    dplyr::summarize(meanNES = mean(meanNES, na.rm=TRUE),
                     frac_sig = mean(minFDR<=0.05, na.rm=TRUE),
                     .groups="drop")
  readr::write_csv(module_loc, file.path(dirs$tables_modules, "module_localization_table.csv"))

  # Replication per module
  module_repl <- mod_theme_comp %>%
    dplyr::mutate(sign=sign(meanNES), sig=minFDR<=0.05) %>%
    dplyr::group_by(db, module) %>%
    dplyr::summarize(n_comp=dplyr::n(),
                     n_sig=sum(sig, na.rm=TRUE),
                     n_pos=sum(sig & sign>0, na.rm=TRUE),
                     n_neg=sum(sig & sign<0, na.rm=TRUE),
                     dom_dir=ifelse(n_pos>=n_neg,"up","down"),
                     repl_score=pmax(n_pos,n_neg),
                     .groups="drop") %>%
    dplyr::arrange(dplyr::desc(repl_score))
  readr::write_csv(module_repl, file.path(dirs$tables_modules, "module_replication_scores.csv"))

  # Comparison class labeling
  module_loc_by_class <- mod_theme_comp %>%
    dplyr::left_join(mod_comp_class_map %>% dplyr::select(context, comparison_class), by="context") %>%
    dplyr::filter(!is.na(comparison_class)) %>%
    dplyr::group_by(db, module, comparison_class, region, layer) %>%
    dplyr::summarize(meanNES = mean(meanNES, na.rm=TRUE),
                     n_comp  = dplyr::n_distinct(context),
                     .groups="drop")
  readr::write_csv(module_loc_by_class, file.path(dirs$tables_modules, "module_localization_by_comparison_class.csv"))

  module_class_labels <- mod_theme_comp %>%
    dplyr::left_join(mod_comp_class_map, by="context") %>%
    dplyr::filter(!is.na(comparison_class)) %>%
    dplyr::group_by(db, module, region, layer, comparison_class) %>%
    dplyr::summarize(n_comp = dplyr::n_distinct(context), .groups="drop") %>%
    dplyr::mutate(code = dplyr::case_when(
      comparison_class=="sus vs res" ~ "SR",
      comparison_class=="res vs con" ~ "RC",
      comparison_class=="sus vs con" ~ "SC",
      TRUE ~ ""
    ),
    label = paste0(code, ":", n_comp))
  module_loc_by_class_lab <- module_loc_by_class %>%
    dplyr::left_join(module_class_labels %>% dplyr::select(db,module,region,layer,comparison_class,label),
                     by=c("db","module","region","layer","comparison_class"))
  readr::write_csv(module_loc_by_class_lab, file.path(dirs$tables_modules, "module_localization_by_comparison_class_labeled.csv"))

  # Driver strengths SR vs SC
  mod_sr_sc <- mod_theme_comp %>%
    dplyr::left_join(mod_comp_class_map %>% dplyr::select(context, comparison_class), by="context") %>%
    dplyr::filter(!is.na(comparison_class)) %>%
    dplyr::group_by(db, module, region, layer, comparison_class) %>%
    dplyr::summarize(meanNES_class = mean(meanNES, na.rm=TRUE),
                     n_ctx = dplyr::n_distinct(context), .groups="drop") %>%
    tidyr::pivot_wider(names_from = comparison_class, values_from = c(meanNES_class, n_ctx), values_fill = 0) %>%
    dplyr::rename(
      meanNES_class_sus_vs_res = `meanNES_class_sus vs res`,
      meanNES_class_sus_vs_con = `meanNES_class_sus vs con`,
      n_ctx_sus_vs_res         = `n_ctx_sus vs res`,
      n_ctx_sus_vs_con         = `n_ctx_sus vs con`
    ) %>%
    dplyr::mutate(
      SR_strength = meanNES_class_sus_vs_res,
      SC_strength = meanNES_class_sus_vs_con,
      driver_class = dplyr::case_when(
        abs(SR_strength) > abs(SC_strength) ~ "sus vs res",
        abs(SC_strength) > abs(SR_strength) ~ "sus vs con",
        TRUE ~ "tie"
      ),
      driver_gap = abs(SR_strength) - abs(SC_strength)
    )
  readr::write_csv(mod_sr_sc, file.path(dirs$tables_modules, "module_driver_strength.csv"))

  # Comparison-level long with class label
  mod_comp_level <- mod_theme_comp %>%
    dplyr::left_join(mod_comp_class_map %>% dplyr::select(context, comparison_class), by="context") %>%
    dplyr::filter(!is.na(comparison_class)) %>%
    dplyr::mutate(ctx_label = paste0(comparison_class, "::", context))
  readr::write_csv(mod_comp_level, file.path(dirs$tables_modules, "module_comparison_level_long.csv"))

  # Direction summaries
  mod_dom_dir <- mod_theme_comp %>%
    dplyr::left_join(mod_rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
    dplyr::filter(!is.na(region), !is.na(layer), minFDR<=0.05) %>%
    dplyr::group_by(db, module, region, layer) %>%
    dplyr::summarize(medNES = stats::median(meanNES, na.rm=TRUE),
                     n_sig  = dplyr::n(), .groups="drop") %>%
    dplyr::mutate(dir = dplyr::case_when(medNES > 0 ~ "up", medNES < 0 ~ "down", TRUE ~ "flat"))
  readr::write_csv(mod_dom_dir, file.path(dirs$tables_modules, "module_dominant_direction.csv"))

  mod_symmetry <- mod_theme_comp %>%
    dplyr::left_join(mod_rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
    dplyr::filter(!is.na(region), !is.na(layer)) %>%
    dplyr::group_by(db, region, layer) %>%
    dplyr::summarize(p_sym = tryCatch({
        x <- meanNES; x <- x[is.finite(x)]
        if (length(x) >= 6) stats::wilcox.test(x, mu=0, exact=FALSE)$p.value else NA_real_
      }, error=function(e) NA_real_), .groups="drop")
  readr::write_csv(mod_symmetry, file.path(dirs$tables_modules, "module_directional_bias_pvalues.csv"))

  # Plots: module localization heatmaps for top modules
  plot_module_localization <- function(dbn, n_mod=12) {
    tops <- module_repl %>% dplyr::filter(db==dbn) %>% dplyr::slice_max(repl_score, n=n_mod) %>% dplyr::pull(module)
    df <- module_loc %>% dplyr::filter(db==dbn, module %in% tops)
    if (!nrow(df)) return(NULL)
    p <- ggplot2::ggplot(df, ggplot2::aes(x=layer,y=region,fill=meanNES)) +
      ggplot2::geom_tile(color="grey90") +
      ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b") +
      ggplot2::facet_wrap(~ module, ncol=4) +
      ggplot2::labs(title=paste0("Module localization — ", dbn), x="Layer", y="Region") +
      ggplot2::theme_minimal(10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))
    fp <- file.path(dirs$plots_modules, paste0("module_localization_", dbn, ".svg"))
    ggplot2::ggsave(fp, p, width=22, height=14); message("[module-loc] Saved: ", fp)
    invisible(p)
  }
  invisible(lapply(unique(module_loc$db), plot_module_localization))

  # Plots: class-labeled module localization
plot_module_localization_by_class_labeled <- function(dbn) {
  df_all <- module_loc_by_class_lab %>% dplyr::filter(db==dbn)
  tops <- module_repl %>% dplyr::filter(db==dbn) %>% dplyr::slice_max(repl_score, n=12) %>% dplyr::pull(module)
  if (!length(tops) || !any(df_all$module %in% tops)) {
    tops <- df_all %>% dplyr::count(module, sort=TRUE) %>% dplyr::slice_head(n=12) %>% dplyr::pull(module)
  }
  df <- df_all %>% dplyr::filter(module %in% tops)
  if (!nrow(df)) { message("[module-class-labeled] No data for ", dbn, " after fallback."); return(NULL) }

  class_levels <- c("sus vs res","res vs con","sus vs con")
  df$comparison_class <- factor(df$comparison_class, levels = class_levels)

  p <- ggplot2::ggplot(df, ggplot2::aes(x=layer, y=region, fill=meanNES)) +
    ggplot2::geom_tile(color="grey90") +
    ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b", midpoint=0) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(module),
      cols = ggplot2::vars(comparison_class),
      scales = "free",
      drop = FALSE,
      labeller = ggplot2::labeller(module = ggplot2::label_wrap_gen(width = 18))
    ) +
    ggplot2::labs(title=paste0("Module localization — comparison class (labels=SR/RC/SC:n) — ", dbn),
                  x="Layer", y="Region") +
    ggplot2::theme_minimal(base_size=10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1)) +
    ggplot2::geom_text(ggplot2::aes(label = dplyr::coalesce(label, as.character(n_comp))),
                       color="black", size=3, na.rm=TRUE)

  fp <- file.path(dirs$plots_modules, paste0("module_localization_", dbn, "_by_comparison_class_labeled.svg"))
  ggplot2::ggsave(fp, p, width=30, height=20); message("[module-class-labeled] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(module_loc_by_class_lab$db), plot_module_localization_by_class_labeled))

# Plots: module comparison matrices (by module)
plot_module_by_comparisons <- function(dbn, max_cols=60) {
  df_db <- mod_comp_level %>% dplyr::filter(db==dbn)
  top_mods <- module_repl %>% dplyr::filter(db==dbn) %>% dplyr::slice_max(repl_score, n=8) %>% dplyr::pull(module)
  if (!length(top_mods)) { message("[mod-comp-matrix] No top modules for ", dbn); return(NULL) }
  lapply(top_mods, function(md) {
    df <- df_db %>% dplyr::filter(module==md)
    ctx_ord <- df %>% dplyr::count(ctx_label, sort=TRUE) %>% dplyr::slice_head(n=max_cols) %>% dplyr::pull(ctx_label)
    df <- df %>% dplyr::filter(ctx_label %in% ctx_ord)
    df$ctx_label <- factor(df$ctx_label, levels=ctx_ord)
    p <- ggplot2::ggplot(df, ggplot2::aes(x=ctx_label, y=paste(region, layer, sep="_"), fill=meanNES)) +
      ggplot2::geom_tile(color="grey90") +
      ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b", midpoint=0) +
      ggplot2::labs(title=paste0("Module ", md, " — comparison-level NES — ", dbn),
                    x="Comparison (class::context)", y="Region_Layer") +
      ggplot2::theme_minimal(base_size=10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle=60, hjust=1, vjust=1))
    fp <- file.path(dirs$plots_modules, paste0("module_", gsub("[^A-Za-z0-9]+","_", md), "_", dbn, "_comparison_matrix.svg"))
    ggplot2::ggsave(fp, p, width=28, height=10); message("[mod-comp-matrix] Saved: ", fp)
    invisible(p)
  })
}
invisible(lapply(unique(mod_comp_level$db), plot_module_by_comparisons))

# Plots: module driver map (SR vs SC)
plot_module_driver_map <- function(dbn) {
  top_mods <- module_repl %>% dplyr::filter(db==dbn) %>% dplyr::slice_max(repl_score, n=12) %>% dplyr::pull(module)
  df <- mod_sr_sc %>% dplyr::filter(db==dbn, module %in% top_mods)
  if (!nrow(df)) return(NULL)
  df <- df %>% dplyr::mutate(n_lab = paste0(n_ctx_sus_vs_res, "/", n_ctx_sus_vs_con))
  p <- ggplot2::ggplot(df, ggplot2::aes(x=layer, y=region, fill=driver_class, alpha=abs(driver_gap))) +
    ggplot2::geom_tile(color="white") +
    ggplot2::geom_text(ggplot2::aes(label = n_lab), size=3) +
    ggplot2::scale_fill_manual(values = c("sus vs res" = "#377eb8", "sus vs con" = "#e41a1c", "tie" = "grey80")) +
    ggplot2::scale_alpha(range = c(0.3, 1), guide = "none") +
    ggplot2::facet_wrap(~ module, ncol = 4, labeller = ggplot2::labeller(module = ggplot2::label_wrap_gen(18))) +
    ggplot2::labs(title = paste0(dbn, " — module drivers (SR vs SC)"),
                  x = "Layer", y = "Region", fill = "Class") +
    ggplot2::theme_minimal(10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))
  fp <- file.path(dirs$plots_modules, paste0("module_drivers_susceptibility_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=22, height=12); message("[mod-drivers] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(mod_sr_sc$db), plot_module_driver_map))

# ======================
# MODULE ANALYSES PARALLEL
# ======================
if (exists("module_gsea_df") && nrow(module_gsea_df) && nrow(mod_rl_both)) {
  # Aggregation per module x comparison (region-layer on left side)
  mod_theme_comp <- module_gsea_df %>%
    dplyr::left_join(mod_rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
    dplyr::filter(!is.na(region), !is.na(layer)) %>%
    dplyr::group_by(db, module, context, region, layer) %>%
    dplyr::summarize(meanNES = mean(NES, na.rm=TRUE),
                     minFDR  = suppressWarnings(min(FDR, na.rm=TRUE)),
                     .groups="drop")
  readr::write_csv(mod_theme_comp, file.path(dirs$tables_modules, "module_theme_by_comparison.csv"))

  # Localization table per module
  module_loc <- mod_theme_comp %>%
    dplyr::group_by(db, module, region, layer) %>%
    dplyr::summarize(meanNES = mean(meanNES, na.rm=TRUE),
                     frac_sig = mean(minFDR<=0.05, na.rm=TRUE),
                     .groups="drop")
  readr::write_csv(module_loc, file.path(dirs$tables_modules, "module_localization_table.csv"))

  # Replication per module
  module_repl <- mod_theme_comp %>%
    dplyr::mutate(sign=sign(meanNES), sig=minFDR<=0.05) %>%
    dplyr::group_by(db, module) %>%
    dplyr::summarize(n_comp=dplyr::n(),
                     n_sig=sum(sig, na.rm=TRUE),
                     n_pos=sum(sig & sign>0, na.rm=TRUE),
                     n_neg=sum(sig & sign<0, na.rm=TRUE),
                     dom_dir=ifelse(n_pos>=n_neg,"up","down"),
                     repl_score=pmax(n_pos,n_neg),
                     .groups="drop") %>%
    dplyr::arrange(dplyr::desc(repl_score))
  readr::write_csv(module_repl, file.path(dirs$tables_modules, "module_replication_scores.csv"))

  # Class-labeled localization
  module_loc_by_class <- mod_theme_comp %>%
    dplyr::left_join(mod_comp_class_map %>% dplyr::select(context, comparison_class), by="context") %>%
    dplyr::filter(!is.na(comparison_class)) %>%
    dplyr::group_by(db, module, comparison_class, region, layer) %>%
    dplyr::summarize(meanNES = mean(meanNES, na.rm=TRUE),
                     n_comp  = dplyr::n_distinct(context),
                     .groups="drop")
  readr::write_csv(module_loc_by_class, file.path(dirs$tables_modules, "module_localization_by_comparison_class.csv"))

  module_class_labels <- mod_theme_comp %>%
    dplyr::left_join(mod_comp_class_map, by="context") %>%
    dplyr::filter(!is.na(comparison_class)) %>%
    dplyr::group_by(db, module, region, layer, comparison_class) %>%
    dplyr::summarize(n_comp = dplyr::n_distinct(context), .groups="drop") %>%
    dplyr::mutate(code = dplyr::case_when(
      comparison_class=="sus vs res" ~ "SR",
      comparison_class=="res vs con" ~ "RC",
      comparison_class=="sus vs con" ~ "SC",
      TRUE ~ ""
    ),
    label = paste0(code, ":", n_comp))
  module_loc_by_class_lab <- module_loc_by_class %>%
    dplyr::left_join(module_class_labels %>% dplyr::select(db,module,region,layer,comparison_class,label),
                     by=c("db","module","region","layer","comparison_class"))
  readr::write_csv(module_loc_by_class_lab, file.path(dirs$tables_modules, "module_localization_by_comparison_class_labeled.csv"))

  # Driver strengths SR vs SC
  mod_sr_sc <- mod_theme_comp %>%
    dplyr::left_join(mod_comp_class_map %>% dplyr::select(context, comparison_class), by="context") %>%
    dplyr::filter(!is.na(comparison_class)) %>%
    dplyr::group_by(db, module, region, layer, comparison_class) %>%
    dplyr::summarize(meanNES_class = mean(meanNES, na.rm=TRUE),
                     n_ctx = dplyr::n_distinct(context), .groups="drop") %>%
    tidyr::pivot_wider(names_from = comparison_class, values_from = c(meanNES_class, n_ctx), values_fill = 0) %>%
    dplyr::rename(
      meanNES_class_sus_vs_res = `meanNES_class_sus vs res`,
      meanNES_class_sus_vs_con = `meanNES_class_sus vs con`,
      n_ctx_sus_vs_res         = `n_ctx_sus vs res`,
      n_ctx_sus_vs_con         = `n_ctx_sus vs con`
    ) %>%
    dplyr::mutate(
      SR_strength = meanNES_class_sus_vs_res,
      SC_strength = meanNES_class_sus_vs_con,
      driver_class = dplyr::case_when(
        abs(SR_strength) > abs(SC_strength) ~ "sus vs res",
        abs(SC_strength) > abs(SR_strength) ~ "sus vs con",
        TRUE ~ "tie"
      ),
      driver_gap = abs(SR_strength) - abs(SC_strength)
    )
  readr::write_csv(mod_sr_sc, file.path(dirs$tables_modules, "module_driver_strength.csv"))

  # Comparison-level long with class label
  mod_comp_level <- mod_theme_comp %>%
    dplyr::left_join(mod_comp_class_map %>% dplyr::select(context, comparison_class), by="context") %>%
    dplyr::filter(!is.na(comparison_class)) %>%
    dplyr::mutate(ctx_label = paste0(comparison_class, "::", context))
  readr::write_csv(mod_comp_level, file.path(dirs$tables_modules, "module_comparison_level_long.csv"))

  # Direction summaries (dominant direction and symmetry)
  mod_dom_dir <- mod_theme_comp %>%
    dplyr::left_join(mod_rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
    dplyr::filter(!is.na(region), !is.na(layer), minFDR<=0.05) %>%
    dplyr::group_by(db, module, region, layer) %>%
    dplyr::summarize(medNES = stats::median(meanNES, na.rm=TRUE),
                     n_sig  = dplyr::n(), .groups="drop") %>%
    dplyr::mutate(dir = dplyr::case_when(medNES > 0 ~ "up", medNES < 0 ~ "down", TRUE ~ "flat"))
  readr::write_csv(mod_dom_dir, file.path(dirs$tables_modules, "module_dominant_direction.csv"))

  mod_symmetry <- mod_theme_comp %>%
    dplyr::left_join(mod_rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
    dplyr::filter(!is.na(region), !is.na(layer)) %>%
    dplyr::group_by(db, region, layer) %>%
    dplyr::summarize(p_sym = tryCatch({
        x <- meanNES; x <- x[is.finite(x)]
        if (length(x) >= 6) stats::wilcox.test(x, mu=0, exact=FALSE)$p.value else NA_real_
      }, error=function(e) NA_real_), .groups="drop")
  readr::write_csv(mod_symmetry, file.path(dirs$tables_modules, "module_directional_bias_pvalues.csv"))

  # ----------------
  # Module plots
  # ----------------
  # Localization heatmaps for top modules
  plot_module_localization <- function(dbn, n_mod=12) {
    tops <- module_repl %>% dplyr::filter(db==dbn) %>% dplyr::slice_max(repl_score, n=n_mod) %>% dplyr::pull(module)
    df <- module_loc %>% dplyr::filter(db==dbn, module %in% tops)
    if (!nrow(df)) return(NULL)
    p <- ggplot2::ggplot(df, ggplot2::aes(x=layer,y=region,fill=meanNES)) +
      ggplot2::geom_tile(color="grey90") +
      ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b", midpoint=0) +
      ggplot2::facet_wrap(~ module, ncol=4, labeller = ggplot2::labeller(module = ggplot2::label_wrap_gen(18))) +
      ggplot2::labs(title=paste0("Module localization — ", dbn), x="Layer", y="Region") +
      ggplot2::theme_minimal(10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))
    fp <- file.path(dirs$plots_modules, paste0("module_localization_", dbn, ".svg"))
    ggplot2::ggsave(fp, p, width=22, height=14); message("[module-loc] Saved: ", fp)
    invisible(p)
  }
  invisible(lapply(unique(module_loc$db), plot_module_localization))

  # Plots: class-labeled module localization
  plot_module_localization_by_class_labeled <- function(dbn) {
    df_all <- module_loc_by_class_lab %>% dplyr::filter(db==dbn)
    tops <- module_repl %>% dplyr::filter(db==dbn) %>% dplyr::slice_max(repl_score, n=12) %>% dplyr::pull(module)
    if (!length(tops) || !any(df_all$module %in% tops)) {
      tops <- df_all %>% dplyr::count(module, sort=TRUE) %>% dplyr::slice_head(n=12) %>% dplyr::pull(module)
    }
    df <- df_all %>% dplyr::filter(module %in% tops)
    if (!nrow(df)) { message("[module-class-labeled] No data for ", dbn, " after fallback."); return(NULL) }

    class_levels <- c("sus vs res","res vs con","sus vs con")
    df$comparison_class <- factor(df$comparison_class, levels = class_levels)

    p <- ggplot2::ggplot(df, ggplot2::aes(x=layer, y=region, fill=meanNES)) +
      ggplot2::geom_tile(color="grey90") +
      ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b", midpoint=0) +
      ggplot2::facet_grid(
        rows = ggplot2::vars(module),
        cols = ggplot2::vars(comparison_class),
        scales = "free",
        drop = FALSE,
        labeller = ggplot2::labeller(module = ggplot2::label_wrap_gen(width = 18))
      ) +
      ggplot2::labs(title=paste0("Module localization — comparison class (labels=SR/RC/SC:n) — ", dbn),
                    x="Layer", y="Region") +
      ggplot2::theme_minimal(base_size=10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1)) +
      ggplot2::geom_text(ggplot2::aes(label = dplyr::coalesce(label, as.character(n_comp))),
                         color="black", size=3, na.rm=TRUE)

    fp <- file.path(dirs$plots_modules, paste0("module_localization_", dbn, "_by_comparison_class_labeled.svg"))
    ggplot2::ggsave(fp, p, width=30, height=20); message("[module-class-labeled] Saved: ", fp)
    invisible(p)
  }
  invisible(lapply(unique(module_loc_by_class_lab$db), plot_module_localization_by_class_labeled))

  # Plots: module comparison matrices (by module)
  plot_module_by_comparisons <- function(dbn, max_cols=60) {
    df_db <- mod_comp_level %>% dplyr::filter(db==dbn)
    top_mods <- module_repl %>% dplyr::filter(db==dbn) %>% dplyr::slice_max(repl_score, n=8) %>% dplyr::pull(module)
    if (!length(top_mods)) { message("[mod-comp-matrix] No top modules for ", dbn); return(NULL) }
    lapply(top_mods, function(md) {
      df <- df_db %>% dplyr::filter(module==md)
      ctx_ord <- df %>% dplyr::count(ctx_label, sort=TRUE) %>% dplyr::slice_head(n=max_cols) %>% dplyr::pull(ctx_label)
      df <- df %>% dplyr::filter(ctx_label %in% ctx_ord)
      df$ctx_label <- factor(df$ctx_label, levels=ctx_ord)
      p <- ggplot2::ggplot(df, ggplot2::aes(x=ctx_label, y=paste(region, layer, sep="_"), fill=meanNES)) +
        ggplot2::geom_tile(color="grey90") +
        ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b", midpoint=0) +
        ggplot2::labs(title=paste0("Module ", md, " — comparison-level NES — ", dbn),
                      x="Comparison (class::context)", y="Region_Layer") +
        ggplot2::theme_minimal(base_size=10) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle=60, hjust=1, vjust=1))
      fp <- file.path(dirs$plots_modules, paste0("module_", gsub("[^A-Za-z0-9]+","_", md), "_", dbn, "_comparison_matrix.svg"))
      ggplot2::ggsave(fp, p, width=28, height=10); message("[mod-comp-matrix] Saved: ", fp)
      invisible(p)
    })
  }
  invisible(lapply(unique(mod_comp_level$db), plot_module_by_comparisons))

  # Plots: module driver map (SR vs SC)
  plot_module_driver_map <- function(dbn) {
    top_mods <- module_repl %>% dplyr::filter(db==dbn) %>% dplyr::slice_max(repl_score, n=12) %>% dplyr::pull(module)
    df <- mod_sr_sc %>% dplyr::filter(db==dbn, module %in% top_mods)
    if (!nrow(df)) return(NULL)
    df <- df %>% dplyr::mutate(n_lab = paste0(n_ctx_sus_vs_res, "/", n_ctx_sus_vs_con))
    p <- ggplot2::ggplot(df, ggplot2::aes(x=layer, y=region, fill=driver_class, alpha=abs(driver_gap))) +
      ggplot2::geom_tile(color="white") +
      ggplot2::geom_text(ggplot2::aes(label = n_lab), size=3) +
      ggplot2::scale_fill_manual(values = c("sus vs res" = "#377eb8", "sus vs con" = "#e41a1c", "tie" = "grey80")) +
      ggplot2::scale_alpha(range = c(0.3, 1), guide = "none") +
      ggplot2::facet_wrap(~ module, ncol = 4, labeller = ggplot2::labeller(module = ggplot2::label_wrap_gen(18))) +
      ggplot2::labs(title = paste0(dbn, " — module drivers (SR vs SC)"),
                    x = "Layer", y = "Region", fill = "Class") +
      ggplot2::theme_minimal(10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))
    fp <- file.path(dirs$plots_modules, paste0("module_drivers_susceptibility_", dbn, ".svg"))
    ggplot2::ggsave(fp, p, width=22, height=12); message("[mod-drivers] Saved: ", fp)
    invisible(p)
  }
  invisible(lapply(unique(mod_sr_sc$db), plot_module_driver_map))
}

# 12) NA diagnostics -----------------------------------------------------------
cov_tbl <- assign_themes %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by = "context") %>%
  dplyr::filter(!is.na(region), !is.na(layer)) %>%
  dplyr::group_by(db, theme, region, layer) %>%
  dplyr::summarize(n_comp_cover=dplyr::n_distinct(comp_dir), n_sig_cover=sum(FDR<=0.05, na.rm=TRUE), minFDR_all=suppressWarnings(min(FDR, na.rm=TRUE)), .groups="drop")
na_diag <- theme_loc %>%
  dplyr::full_join(cov_tbl, by = c("db","theme","region","layer")) %>%
  dplyr::mutate(cause = dplyr::case_when(
    is.na(meanNES) & (is.na(n_comp_cover) | n_comp_cover == 0) ~ "no_coverage",
    is.na(meanNES) & (n_comp_cover > 0) & (is.na(n_sig_cover) | n_sig_cover == 0) ~ "masked_nonsignificant",
    is.na(meanNES) ~ "other_missing",
    TRUE ~ "has_value"
  ))
readr::write_csv(na_diag, file.path(dirs$tables_qa, "localization_NA_diagnostics.csv"))

# 13) Replication summary ------------------------------------------------------
theme_repl <- theme_comp %>% dplyr::mutate(sign=sign(meanNES), sig=minFDR<=0.05) %>%
  dplyr::group_by(db, theme) %>% dplyr::summarize(n_comp=dplyr::n(), n_sig=sum(sig, na.rm=TRUE), n_pos=sum(sig & sign>0, na.rm=TRUE), n_neg=sum(sig & sign<0, na.rm=TRUE), dom_dir=ifelse(n_pos>=n_neg,"up","down"), repl_score=pmax(n_pos,n_neg), .groups="drop") %>%
  dplyr::arrange(dplyr::desc(repl_score))
readr::write_csv(theme_repl, file.path(dirs$tables_gsea, "theme_replication_scores.csv"))

top_themes <- theme_repl %>% dplyr::group_by(db) %>% dplyr::slice_max(order_by=repl_score, n=12) %>% dplyr::ungroup()

# 14) Plots -------------------------------------------------------------------
plot_theme_localization <- function(dbn) {
  df <- theme_loc %>% dplyr::filter(db==dbn, theme %in% (top_themes %>% dplyr::filter(db==dbn) %>% dplyr::pull(theme)))
  if (!nrow(df)) { message("[anatomy] No data for ", dbn); return(NULL) }
  p <- ggplot2::ggplot(df, ggplot2::aes(x=layer, y=region, fill=meanNES)) +
    ggplot2::geom_tile(color="grey90") +
    ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b") +
    ggplot2::facet_wrap(~ theme, scales="free", ncol=4) +
    ggplot2::labs(title=paste0("Localization by region-layer — ", dbn), x="Layer", y="Region") +
    ggplot2::theme_minimal(base_size=10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))
  fp <- file.path(dirs$plots_anat, paste0("localization_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=22, height=14); message("[anatomy] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(theme_loc$db), plot_theme_localization))

plot_theme_localization_by_class_labeled <- function(dbn) {
  df_all <- theme_loc_by_class_lab %>% dplyr::filter(db==dbn)
  tops <- top_themes %>% dplyr::filter(db==dbn) %>% dplyr::pull(theme)
  if (!length(tops) || !any(df_all$theme %in% tops)) {
    tops <- df_all %>% dplyr::count(theme, sort=TRUE) %>% dplyr::slice_head(n=12) %>% dplyr::pull(theme)
  }
  df <- df_all %>% dplyr::filter(theme %in% tops)
  if (!nrow(df)) { message("[class-labeled] No data for ", dbn, " after fallback."); return(NULL) }
  class_levels <- c("sus vs res","res vs con","sus vs con")
  df$comparison_class <- factor(df$comparison_class, levels = class_levels)
  p <- ggplot2::ggplot(df, ggplot2::aes(x=layer, y=region, fill=meanNES)) +
    ggplot2::geom_tile(color="grey90") +
    ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b") +
    ggplot2::facet_grid(rows = ggplot2::vars(theme), cols = ggplot2::vars(comparison_class), scales="free", drop=FALSE) +
    ggplot2::labs(title=paste0("Localization — comparison class (labeled) — ", dbn), x="Layer", y="Region") +
    ggplot2::theme_minimal(base_size=10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1)) +
    ggplot2::geom_text(ggplot2::aes(label = dplyr::coalesce(label, as.character(n_comp))), color="black", size=3, na.rm=TRUE)
  fp <- file.path(dirs$plots_class, paste0("localization_", dbn, "_by_comparison_class_labeled.svg"))
  ggplot2::ggsave(fp, p, width=30, height=20); message("[class-labeled] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(theme_loc_by_class_lab$db), plot_theme_localization_by_class_labeled))

plot_theme_localization_left_cond <- function(dbn) {
  df <- theme_loc_left_cond %>% dplyr::filter(db==dbn, theme %in% (top_themes %>% dplyr::filter(db==dbn) %>% dplyr::pull(theme)))
  if (!nrow(df)) { message("[left-cond] No data for ", dbn); return(NULL) }
  p <- ggplot2::ggplot(df, ggplot2::aes(x=layer, y=region, fill=meanNES)) +
    ggplot2::geom_tile(color="grey90") +
    ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b") +
    ggplot2::facet_grid(rows = ggplot2::vars(theme), cols = ggplot2::vars(condition), scales="free") +
    ggplot2::labs(title=paste0("Localization — tested condition (sus/res) — ", dbn), x="Layer", y="Region") +
    ggplot2::theme_minimal(base_size=10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))
  fp <- file.path(dirs$plots_cond, paste0("localization_", dbn, "_tested_condition.svg"))
  ggplot2::ggsave(fp, p, width=26, height=20); message("[left-cond] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(theme_loc_left_cond$db), plot_theme_localization_left_cond))

plot_theme_localization_right_cond <- function(dbn) {
  df <- theme_loc_right_cond %>% dplyr::filter(db==dbn, theme %in% (top_themes %>% dplyr::filter(db==dbn) %>% dplyr::pull(theme)))
  if (!nrow(df)) { message("[right-cond] No data for ", dbn); return(NULL) }
  p <- ggplot2::ggplot(df, ggplot2::aes(x=layer, y=region, fill=meanNES)) +
    ggplot2::geom_tile(color="grey90") +
    ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b") +
    ggplot2::facet_grid(rows = ggplot2::vars(theme), cols = ggplot2::vars(condition), scales="free") +
    ggplot2::labs(title=paste0("Localization — baseline condition (res/con) — ", dbn), x="Layer", y="Region") +
    ggplot2::theme_minimal(base_size=10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))
  fp <- file.path(dirs$plots_cond, paste0("localization_", dbn, "_baseline_condition.svg"))
  ggplot2::ggsave(fp, p, width=26, height=20); message("[right-cond] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(theme_loc_right_cond$db), plot_theme_localization_right_cond))

plot_theme_by_comparisons <- function(dbn, max_cols=60) {
  df_db <- comp_level %>% dplyr::filter(db==dbn)
  themes_here <- df_db %>% dplyr::count(theme, sort=TRUE) %>% dplyr::slice_head(n=12) %>% dplyr::pull(theme)
  if (!length(themes_here)) { message("[comp-matrix] No themes with class data for ", dbn); return(NULL) }
  lapply(themes_here, function(th) {
    df <- df_db %>% dplyr::filter(theme==th)
    ctx_ord <- df %>% dplyr::count(ctx_label, sort=TRUE) %>% dplyr::slice_head(n=max_cols) %>% dplyr::pull(ctx_label)
    df <- df %>% dplyr::filter(ctx_label %in% ctx_ord)
    df$ctx_label <- factor(df$ctx_label, levels=ctx_ord)
    p <- ggplot2::ggplot(df, ggplot2::aes(x=ctx_label, y=paste(region, layer, sep="_"), fill=meanNES)) +
      ggplot2::geom_tile(color="grey90") +
      ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b") +
      ggplot2::labs(title=paste0("Theme ", th, " — comparison-level NES — ", dbn), x="Comparison (class::context)", y="Region_Layer") +
      ggplot2::theme_minimal(base_size=10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle=60, hjust=1, vjust=1))
    fp <- file.path(dirs$plots_cm, paste0("theme_", gsub("[^A-Za-z0-9]+","_", th), "_", dbn, "_comparison_matrix.svg"))
    ggplot2::ggsave(fp, p, width=28, height=10); message("[comp-matrix] Saved: ", fp)
    invisible(p)
  })
}
invisible(lapply(unique(comp_level$db), plot_theme_by_comparisons))

# 14b) Compact panel (Anatomy + Class + Replication) per DB
plot_compact_panel <- function(dbn, n_themes = 6) {
  the <- top_themes %>% dplyr::filter(db==dbn) %>% dplyr::slice_head(n=n_themes) %>% dplyr::pull(theme)

  p1 <- ggplot2::ggplot(theme_loc %>% dplyr::filter(db==dbn, theme %in% the),
                        ggplot2::aes(x=layer,y=region,fill=meanNES)) +
    ggplot2::geom_tile(color="grey90") +
    ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b") +
    ggplot2::facet_wrap(~ theme, ncol=3) +
    ggplot2::labs(title=paste0(dbn, " — localization (NES)")) + ggplot2::theme_minimal(9)

  p2 <- ggplot2::ggplot(theme_loc_by_class_lab %>% dplyr::filter(db==dbn, theme %in% the),
                        ggplot2::aes(x=layer,y=region,fill=meanNES,label=label)) +
    ggplot2::geom_tile(color="grey95") +
    ggplot2::geom_text(size=2, color="black", na.rm=TRUE) +
    ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b") +
    ggplot2::facet_grid(rows=ggplot2::vars(theme), cols=ggplot2::vars(comparison_class), drop=FALSE) +
    ggplot2::labs(title="comparison class (labels=SR/RC/SC:n)") + ggplot2::theme_minimal(9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45,hjust=1))

  p3 <- ggplot2::ggplot(theme_repl %>% dplyr::filter(db==dbn, theme %in% the),
                        ggplot2::aes(x=reorder(theme, repl_score), y=repl_score, fill=dom_dir)) +
    ggplot2::geom_col() + ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values=c(down="#2166ac", up="#b2182b")) +
    ggplot2::labs(title="replication score", x="", y="max(up, down) significant") +
    ggplot2::theme_minimal(9)

  if (requireNamespace("patchwork", quietly=TRUE)) {
    panel <- p1 / p2 / p3 + patchwork::plot_layout(heights=c(1,1.1,0.6))
    fp <- file.path(dirs$plots, "Panels", paste0("panel_", dbn, "_compact.svg"))
    dir.create(dirname(fp), recursive=TRUE, showWarnings=FALSE)
    ggplot2::ggsave(fp, panel, width=14, height=16); message("[panel] Saved: ", fp)
  } else {
    message("[panel] patchwork not installed; skipping combined panel for ", dbn)
  }
}
invisible(lapply(unique(theme_loc$db), plot_compact_panel))

# 14c) Dominant direction lattice per DB (median NES, n_sig annotation)
plot_dir_lattice <- function(dbn) {
  df <- dom_dir_tbl %>% dplyr::filter(db==dbn, theme %in% (top_themes %>% dplyr::filter(db==dbn) %>% dplyr::pull(theme)))
  if (!nrow(df)) return(NULL)
  p <- ggplot2::ggplot(df, ggplot2::aes(x=layer, y=region, fill=dir, label=n_sig)) +
    ggplot2::geom_tile(color="white") + ggplot2::geom_text(size=3, color="black") +
    ggplot2::scale_fill_manual(values=c(down="#2166ac", flat="grey85", up="#b2182b")) +
    ggplot2::facet_wrap(~ theme, ncol=4) +
    ggplot2::labs(title=paste0(dbn, " — dominant direction (median NES, n_sig)"),
                  x="Layer", y="Region") + ggplot2::theme_minimal(9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))
  fp <- file.path(dirs$plots_anat, paste0("dominant_direction_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=22, height=12); message("[dir-lattice] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(theme_loc$db), plot_dir_lattice))

# 14d) Directional bias tiles per DB (−log10 p from Wilcoxon against 0)
plot_symmetry <- function(dbn) {
  df <- symmetry_test %>% dplyr::filter(db==dbn)
  p <- ggplot2::ggplot(df, ggplot2::aes(x=layer, y=region, fill=-log10(p_sym))) +
    ggplot2::geom_tile(color="white") +
    ggplot2::scale_fill_viridis_c(option="C", na.value="grey90") +
    ggplot2::labs(title=paste0(dbn, " — directional bias (−log10 p)"), x="Layer", y="Region") +
    ggplot2::theme_minimal(9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))
  fp <- file.path(dirs$plots_anat, paste0("directional_bias_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=8, height=6); message("[bias] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(theme_loc$db), plot_symmetry))

# 14e) Comparison-class balance per theme (stacked fractions across region-layer)
class_balance <- theme_comp %>%
  dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L), by="context") %>%
  dplyr::left_join(comp_class_map, by="context") %>%
  dplyr::filter(!is.na(region), !is.na(layer), !is.na(comparison_class)) %>%
  dplyr::group_by(db, theme, region, layer, comparison_class) %>%
  dplyr::summarize(n_cmp=dplyr::n_distinct(context), .groups="drop") %>%
  dplyr::group_by(db, theme, region, layer) %>%
  dplyr::mutate(frac = n_cmp/sum(n_cmp)) %>% dplyr::ungroup()
readr::write_csv(class_balance, file.path(dirs$tables_class, "comparison_class_balance.csv"))

plot_class_balance <- function(dbn, theme_id) {
  df <- class_balance %>% dplyr::filter(db==dbn, theme==theme_id)
  if (!nrow(df)) return(NULL)
  p <- ggplot2::ggplot(df, ggplot2::aes(x=paste(region,layer,sep="_"), y=frac, fill=comparison_class)) +
    ggplot2::geom_col(width=0.8) +
    ggplot2::scale_fill_manual(values=c("sus vs res"="#377eb8", "res vs con"="#4daf4a", "sus vs con"="#e41a1c")) +
    ggplot2::labs(title=paste0(theme_id, " — comparison balance (", dbn, ")"), x="Region_Layer", y="Fraction of contexts") +
    ggplot2::theme_minimal(9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=60, hjust=1))
  fp <- file.path(dirs$plots_class, paste0("class_balance_", gsub("[^A-Za-z0-9]+","_", theme_id), "_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=12, height=6); message("[class-balance] Saved: ", fp)
  invisible(p)
}
invisible(lapply(split(top_themes, top_themes$db), function(dfdb) {
  dbn <- unique(dfdb$db); ths <- head(dfdb$theme, 8)
  lapply(ths, function(th) plot_class_balance(dbn, th))
}))

# 14f) Contrast strip (compact comparison matrix subset) for top themes
plot_theme_strip <- function(dbn, theme_id, max_ctx=40) {
  df <- comp_level %>% dplyr::filter(db==dbn, theme==theme_id)
  if (!nrow(df)) return(NULL)
  ctx_ord <- df %>% dplyr::count(ctx_label, sort=TRUE) %>% dplyr::slice_head(n=max_ctx) %>% dplyr::pull(ctx_label)
  df <- df %>% dplyr::filter(ctx_label %in% ctx_ord)
  df$ctx_label <- factor(df$ctx_label, levels=ctx_ord)

  p <- ggplot2::ggplot(df, ggplot2::aes(x=ctx_label, y=paste(region,layer,sep="_"), fill=meanNES)) +
    ggplot2::geom_tile(color="grey95") +
    ggplot2::scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b") +
    ggplot2::facet_wrap(~ comparison_class, ncol=1, scales="free_x") +
    ggplot2::labs(title=paste0(theme_id, " — contrast strip (", dbn, ")"), x="Comparison", y="Region_Layer") +
    ggplot2::theme_minimal(9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=60, hjust=1))
  fp <- file.path(dirs$plots_cm, paste0("theme_", gsub("[^A-Za-z0-9]+","_", theme_id), "_", dbn, "_strip.svg"))
  ggplot2::ggsave(fp, p, width=14, height=8); message("[strip] Saved: ", fp)
  invisible(p)
}
invisible(lapply(split(top_themes, top_themes$db), function(dfdb) {
  dbn <- unique(dfdb$db); ths <- head(dfdb$theme, 6)
  lapply(ths, function(th) plot_theme_strip(dbn, th, max_ctx=40))
}))

# 14g) Simple localization graph per DB (themes ↔ region-layer nodes)
t_thr <- 0.6
edges_df <- theme_loc %>%
  dplyr::left_join(
    theme_comp %>%
      dplyr::left_join(rl_both %>% dplyr::transmute(context, region=region_L, layer=layer_L),
                       by="context") %>%
      dplyr::filter(!is.na(region), !is.na(layer), minFDR<=0.05) %>%
      dplyr::group_by(db, theme, region, layer) %>%
      dplyr::summarize(n_sig = dplyr::n(), sign_dir = sign(mean(meanNES, na.rm=TRUE)), .groups="drop"),
    by=c("db","theme","region","layer")
  ) %>%
  # guard: keep rows with numeric NES and threshold on absolute value
  dplyr::filter(!is.na(meanNES), is.finite(meanNES), abs(meanNES) >= t_thr) %>%
  # vectorized fallback for missing n_sig
  dplyr::mutate(n_sig = ifelse(is.na(n_sig), 1L, n_sig),
                src = theme,
                dst = paste(region, layer, sep="_"),
                edge_col = ifelse(meanNES > 0, "#b2182b", "#2166ac"),
                edge_w = scales::rescale(pmax(n_sig, 1L), to=c(0.4, 2)))

plot_localization_graph <- function(dbn) {
  df <- edges_df %>% dplyr::filter(db==dbn)
  if (!nrow(df)) return(NULL)

  # Build graph
  g <- igraph::graph_from_data_frame(df %>% dplyr::select(src, dst), directed=FALSE)

  # Align edge attributes to graph edges
  e_df <- igraph::as_data_frame(g, what="edges")
  key_g  <- paste(e_df$from, e_df$to, sep="||")
  key_df <- paste(df$src, df$dst, sep="||")
  idx <- match(key_g, key_df)
  mis <- which(is.na(idx))
  if (length(mis)) {
    key_rev <- paste(e_df$to[mis], e_df$from[mis], sep="||")
    idx_rev <- match(key_rev, key_df)
    idx[mis] <- idx_rev
  }
  idx[is.na(idx)] <- 1L

  igraph::E(g)$edge_col <- df$edge_col[idx]
  igraph::E(g)$edge_w   <- df$edge_w[idx]

  set.seed(1)
  p <- ggraph::ggraph(g, layout="fr") +
    ggraph::geom_edge_link(ggplot2::aes(edge_colour = edge_col,
                                        edge_width  = edge_w), alpha=0.8) +
    ggraph::geom_node_point(ggplot2::aes(shape = ifelse(grepl("_", name), "RL", "Theme")), size=3) +
    ggraph::geom_node_text(ggplot2::aes(label=name), repel=TRUE, size=3) +
    ggraph::scale_edge_colour_identity() +
    ggraph::scale_edge_width(range=c(0.4,2)) +
    ggplot2::scale_shape_manual(values=c(Theme=16, RL=15)) +
    ggplot2::theme_void() +
    ggplot2::labs(title=paste0(dbn, " — localization graph (|NES|≥", t_thr, ")"))

  fp <- file.path(dirs$plots, "Graphs", paste0("localization_graph_", dbn, ".svg"))
  dir.create(dirname(fp), recursive=TRUE, showWarnings=FALSE)
  ggplot2::ggsave(fp, p, width=10, height=8); message("[graph] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(theme_loc$db), plot_localization_graph))


# 14h) Terms network per DB (Jaccard similarity of GSEA terms)
# --- Terms network (per DB) ---
# Recompute similarity within DB to match your themes more closely
# Term network colored by theme cluster
tokenize <- function(x){
  x <- tolower(x); x <- gsub("[^a-z0-9 ]+"," ",x)
  unique(unlist(strsplit(x,"\\s+")))
}
jaccard <- function(a,b){
  ia <- tokenize(a); ib <- tokenize(b)
  if (!length(ia) || !length(ib)) return(0)
  length(intersect(ia,ib)) / length(union(ia,ib))
}

rep_terms <- gsea_df %>%
  dplyr::inner_join(dplyr::select(themes_tbl, ID, theme), by="ID") %>%
  dplyr::group_by(db, theme, ID, Description) %>%
  dplyr::summarize(freq = dplyr::n(), .groups="drop") %>%
  dplyr::group_by(db, theme) %>%
  dplyr::slice_max(order_by=freq, n=10, with_ties=FALSE) %>%
  dplyr::ungroup()

plot_terms_network <- function(dbn, sim_cut=0.25, max_nodes=250) {
  df <- rep_terms %>% dplyr::filter(db==dbn)
  if (!nrow(df)) return(NULL)
  # control size: top 12 themes
  top_themes <- df %>% dplyr::count(theme, sort=TRUE) %>% dplyr::slice_head(n=12) %>% dplyr::pull(theme)
  df <- df %>% dplyr::filter(theme %in% top_themes)
  if (nrow(df) > max_nodes) {
    per <- max(5L, floor(max_nodes/length(top_themes)))
    df <- df %>% dplyr::group_by(theme) %>% dplyr::slice_head(n=per) %>% dplyr::ungroup()
  }

  terms <- df$Description; n <- length(terms)
  if (n < 2) return(NULL)
  edges <- list()
  for (i in seq_len(n-1)) for (j in (i+1):n) {
    s <- jaccard(terms[i], terms[j]); if (s >= sim_cut) edges[[length(edges)+1]] <- c(i,j,s)
  }
  if (!length(edges)) return(NULL)
  edges <- as.data.frame(do.call(rbind, edges)); names(edges) <- c("i","j","sim")
  edges$i <- as.integer(edges$i); edges$j <- as.integer(edges$j); edges$sim <- as.numeric(edges$sim)

  verts <- tibble::tibble(name = terms, theme = df$theme, term = df$Description)
  g <- igraph::graph_from_data_frame(
    d = tibble::tibble(from = verts$name[edges$i], to = verts$name[edges$j], sim = edges$sim),
    directed = FALSE, vertices = verts
  )

  set.seed(1)
  p <- ggraph::ggraph(g, layout="fr") +
    ggraph::geom_edge_link(ggplot2::aes(width=sim), colour="grey70", alpha=0.35) +
    ggraph::geom_node_point(ggplot2::aes(color=theme), size=2) +
    ggraph::geom_node_text(ggplot2::aes(label=term), size=2.4, repel=TRUE) +
    ggplot2::scale_edge_width(range=c(0.2,1.4), guide="none") +
    ggplot2::theme_void() +
    ggplot2::labs(title=paste0(dbn, " — terms network (sim≥", sim_cut, ")"), color="Theme")
  fp <- file.path(dirs$plots_anat, paste0("terms_network_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=16, height=12); message("[terms-net] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(rep_terms$db), plot_terms_network))

suppressWarnings({
  if (!requireNamespace("uwot", quietly=TRUE)) message("Optional: install.packages('uwot') for UMAP")
})

build_theme_kw_matrix <- function(dbn, top_kw=800) {
  df <- gsea_df %>%
    dplyr::filter(db==dbn) %>%
    dplyr::inner_join(dplyr::select(themes_tbl, ID, theme), by="ID") %>%
    dplyr::distinct(theme, Description)
  if (!nrow(df)) return(NULL)

  toks <- strsplit(tolower(gsub("[^a-z0-9 ]+"," ", df$Description)), "\\s+")
  toks <- lapply(toks, function(v) v[nzchar(v) & !v %in% c("and","or","of","the","to","in","by","for","process","regulation","cellular","activity","protein","pathway","signaling")])

  env <- new.env(parent=emptyenv())
  for (i in seq_len(nrow(df))) {
    th <- df$theme[i]
    for (k in toks[[i]]) {
      key <- paste(th,k,sep="||"); env[[key]] <- (env[[key]] %||% 0L) + 1L
    }
  }
  keys <- ls(env); if (!length(keys)) return(NULL)
  sp <- strsplit(keys,"\\|\\|"); ths <- vapply(sp, `[`, character(1), 1); kws <- vapply(sp, `[`, character(1), 2)
  val <- as.integer(mget(keys, env, ifnotfound=0L))
  mat <- tibble::tibble(theme=ths, kw=kws, val=val) %>%
    dplyr::group_by(kw) %>% dplyr::summarize(dfreq=dplyr::n(), .groups="drop") %>%
    dplyr::right_join(tibble::tibble(theme=ths, kw=kws, val=val), by="kw") %>%
    dplyr::filter(dfreq>=2) %>%
    dplyr::group_by(kw) %>% dplyr::mutate(tf = val/sum(val)) %>% dplyr::ungroup() %>%
    dplyr::group_by(kw) %>% dplyr::mutate(idf = log(n_distinct(theme)/n_distinct(theme[val>0])+1e-6)) %>% dplyr::ungroup() %>%
    dplyr::mutate(tfidf = tf*idf)

  vocab <- mat %>% dplyr::group_by(kw) %>% dplyr::summarize(s=sum(tfidf,na.rm=TRUE), .groups="drop") %>%
    dplyr::arrange(dplyr::desc(s)) %>% dplyr::slice_head(n=top_kw) %>% dplyr::pull(kw)
  mat <- mat %>% dplyr::filter(kw %in% vocab)
  if (!nrow(mat)) return(NULL)

  M <- tidyr::pivot_wider(mat, names_from=kw, values_from=tfidf, values_fill=0)
  M
}

plot_theme_embedding <- function(dbn) {
  M <- build_theme_kw_matrix(dbn)
  if (is.null(M) || nrow(M) < 2) return(NULL)
  meta <- M %>% dplyr::select(theme)
  X <- as.matrix(M %>% dplyr::select(-theme))
  rownames(X) <- meta$theme

  if (requireNamespace("uwot", quietly=TRUE)) {
    set.seed(1); emb <- uwot::umap(X, n_neighbors=10, min_dist=0.2, metric="cosine")
    emb <- as.data.frame(emb); colnames(emb) <- c("U1","U2")
  } else {
    pr <- stats::prcomp(X, scale.=TRUE); emb <- as.data.frame(pr$x[,1:2]); colnames(emb) <- c("U1","U2")
  }
  emb$theme <- rownames(X)

  p <- ggplot2::ggplot(emb, ggplot2::aes(x=U1, y=U2, label=theme)) +
    ggplot2::geom_point(color="#2c7fb8", size=2, alpha=0.85) +
    ggplot2::geom_text(size=2.8, nudge_y=0.02) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title=paste0(dbn, " — theme embedding (", if (requireNamespace("uwot", quietly=TRUE)) "UMAP" else "PCA", ")"),
                  x="Dim 1", y="Dim 2")
  fp <- file.path(dirs$plots_anat, paste0("theme_embedding_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=10, height=8); message("[theme-embed] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(themes_tbl$db), plot_theme_embedding))

plot_theme_dendrogram <- function(dbn) {
  M <- build_theme_kw_matrix(dbn)
  if (is.null(M) || nrow(M) < 2) return(NULL)
  X <- as.matrix(M %>% dplyr::select(-theme))
  rownames(X) <- M$theme
  # cosine distance
  A <- X + 1e-12
  nrm <- sqrt(rowSums(A*A))
  S <- tcrossprod(A / nrm)
  D <- 1 - S
  hc <- stats::hclust(stats::as.dist(D), method="average")

  if (!requireNamespace("ggdendro", quietly=TRUE)) {
    message("Optional: install.packages('ggdendro') for dendrogram plot"); return(NULL)
  }
  dd <- ggdendro::dendro_data(hc)
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data=dd$segments, ggplot2::aes(x=x, y=y, xend=xend, yend=yend)) +
    ggplot2::geom_text(data=dd$labels, ggplot2::aes(x=x, y=y, label=label), angle=90, hjust=1, vjust=0.5, size=2.6) +
    ggplot2::theme_void() + ggplot2::labs(title=paste0(dbn, " — theme dendrogram"))
  fp <- file.path(dirs$plots_anat, paste0("theme_dendrogram_", dbn, ".svg"))
  ggplot2::ggsave(fp, p, width=14, height=10); message("[theme-dend] Saved: ", fp)
  invisible(p)
}
invisible(lapply(unique(themes_tbl$db), plot_theme_dendrogram))

# 14i) GO BP terms per theme (bar, top by FDR/frequency)
plot_bp_terms_per_theme <- function(top_k = 8) {
  df <- gsea_df %>%
    dplyr::filter(db == "GO_BP") %>%
    dplyr::inner_join(dplyr::select(themes_tbl, ID, theme), by="ID") %>%
    dplyr::group_by(theme, Description) %>%
    dplyr::summarize(n_hits = dplyr::n(), meanNES = mean(NES, na.rm=TRUE),
                     bestFDR = suppressWarnings(min(FDR, na.rm=TRUE)), .groups="drop") %>%
    dplyr::group_by(theme) %>%
    dplyr::arrange(bestFDR, dplyr::desc(n_hits), .by_group=TRUE) %>%
    dplyr::slice_head(n = top_k) %>%
    dplyr::ungroup()

  if (!nrow(df)) { message("[bp-terms] No GO_BP data"); return(NULL) }

  # Reorder within each facet
  df <- df %>%
    dplyr::group_by(theme) %>%
    dplyr::mutate(term_order = reorder(Description, -n_hits)) %>%
    dplyr::ungroup()

  p <- ggplot2::ggplot(df, ggplot2::aes(x = n_hits, y = term_order, fill = -log10(bestFDR))) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_viridis_c(option="C", name = "-log10(FDR)") +
    ggplot2::facet_wrap(~ theme, scales = "free_y", ncol = 3) +
    ggplot2::labs(title = "GO BP terms per cluster (top by FDR/frequency)",
                  x = "Term count within cluster", y = "GO BP term") +
    ggplot2::theme_minimal(10)
  fp <- file.path(dirs$plots_anat, "GO_BP_terms_per_theme.svg")
  ggplot2::ggsave(fp, p, width=18, height=14); message("[bp-terms] Saved: ", fp)
  invisible(p)
}
plot_bp_terms_per_theme(top_k = 8)

bp_theme_summaries <- bp_members %>%
  dplyr::group_by(theme) %>%
  dplyr::arrange(bestFDR, dplyr::desc(n_hits), .by_group=TRUE) %>%
  dplyr::summarize(
    n_terms = dplyr::n(),
    top_terms = paste(head(Description, 6), collapse=" | "),
    medianNES = stats::median(meanNES, na.rm=TRUE),
    .groups="drop"
  )
readr::write_csv(bp_theme_summaries, file.path(dirs$tables_gsea, "GO_BP_theme_summaries.csv"))

# 15) ORA summaries (separate) ------------------------------------------------
if (nrow(ora_df)) {
  ora_repl <- ora_df %>% dplyr::group_by(db, ID, Description) %>% dplyr::summarize(n_comp=dplyr::n(), n_sig=sum(Padj<=0.05, na.rm=TRUE), meanPadj=mean(Padj, na.rm=TRUE), .groups="drop") %>% dplyr::arrange(dplyr::desc(n_sig), meanPadj)
  readr::write_csv(ora_repl, file.path(dirs$tables_ora, "ORA_replication_scores.csv"))
  top_ora <- ora_repl %>% dplyr::group_by(db) %>% dplyr::slice_max(n_sig, n=25) %>% dplyr::ungroup()
  p_ora <- ggplot2::ggplot(top_ora, ggplot2::aes(x=reorder(paste(db, Description, sep="::"), n_sig), y=n_sig, fill=db)) +
    ggplot2::geom_col() + ggplot2::coord_flip() +
    ggplot2::labs(title="ORA term replication (count of significant comparisons)", x="DB::Term", y="# significant comps") +
    ggplot2::theme_minimal(base_size=10)
  ggplot2::ggsave(file.path(dirs$plots_ora, "ORA_replication_bar.svg"), p_ora, width=18, height=12)
}

# 16) Executive bundle + artifact check ---------------------------------------
# 16) Executive bundle + artifact check (existing)
readr::write_csv(top_themes,  file.path(dirs$tables_gsea, "top_themes_per_db.csv"))
readr::write_csv(theme_comp %>% dplyr::arrange(db, theme, dplyr::desc(abs(meanNES))),
                 file.path(dirs$tables_gsea, "theme_by_comparison_sorted.csv"))

md <- c(
  "# Consolidated Summary",
  paste0("- Comparisons scanned: ", length(comp_dirs)),
  paste0("- GSEA rows: ", nrow(gsea_df)),
  "",
  "Artifacts:",
  "* Plots/Anatomy/localization_<DB>.svg",
  "* Plots/Class/localization_<DB>_by_comparison_class_labeled.svg",
  "* Plots/Conditions/localization_<DB>_tested_condition.svg",
  "* Plots/Conditions/localization_<DB>_baseline_condition.svg",
  "* Plots/ComparisonMatrices/theme_<THEME>_<DB>_comparison_matrix.svg",
  "* Tables/context_region_layer_condition_both_sides.csv",
  "* Tables/Class/context_comparison_class_map.csv",
  "* Tables/GSEA/localization_table.csv",
  "* Tables/Class/localization_by_comparison_class_table.csv",
  "* Tables/Class/localization_by_comparison_class_labeled.csv",
  "* Tables/Conditions/localization_left_condition_table.csv",
  "* Tables/Conditions/localization_right_condition_table.csv",
  "* Tables/Class/comparison_level_long.csv",
  "* Tables/QA/localization_NA_diagnostics.csv",
  "* Tables/GSEA/all_GSEA_long.csv / themes_table.csv / theme_by_comparison.csv / theme_replication_scores.csv",
  if (nrow(ora_df)) "* Tables/ORA/all_ORA_long.csv / ORA_replication_scores.csv; Plots/ORA/ORA_replication_bar.svg" else "* (no ORA detected)"
)
writeLines(md, con = file.path(out_dir, "CONSOLIDATED_SUMMARY.md"))

# Artifact existence check across subfolders (robust)
exp_list <- list()
if (exists("theme_loc") && nrow(theme_loc)) {
  exp_list[["anat"]] <- file.path(
    dirs$plots_anat, paste0("localization_", unique(theme_loc$db), ".svg")
  )
}
if (exists("theme_loc_by_class_lab") && nrow(theme_loc_by_class_lab)) {
  exp_list[["class"]] <- file.path(
    dirs$plots_class, paste0("localization_", unique(theme_loc_by_class_lab$db), "_by_comparison_class_labeled.svg")
  )
}
if (exists("theme_loc_left_cond") && nrow(theme_loc_left_cond)) {
  exp_list[["cond_left"]] <- file.path(
    dirs$plots_cond, paste0("localization_", unique(theme_loc_left_cond$db), "_tested_condition.svg")
  )
}
if (exists("theme_loc_right_cond") && nrow(theme_loc_right_cond)) {
  exp_list[["cond_right"]] <- file.path(
    dirs$plots_cond, paste0("localization_", unique(theme_loc_right_cond$db), "_baseline_condition.svg")
  )
}
expected <- unique(unlist(exp_list, use.names = FALSE))
exists_tbl <- tibble::tibble(file = expected, exists = file.exists(expected))
readr::write_csv(exists_tbl, file.path(out_dir, "artifact_existence_check.csv"))
print(exists_tbl)

message("Consolidation complete. See ZZ_consolidated for outputs.")

# Optionally: append module artifacts to CONSOLIDATED_SUMMARY.md (existing)
try({
  md_mod <- c(
    "",
    "Module Artifacts:",
    "* Tables/Modules/module_GSEA_long.csv / module_ORA_long.csv",
    "* Tables/Modules/module_theme_by_comparison.csv / module_localization_table.csv",
    "* Tables/Modules/module_replication_scores.csv",
    "* Tables/Modules/module_localization_by_comparison_class.csv",
    "* Tables/Modules/module_localization_by_comparison_class_labeled.csv",
    "* Tables/Modules/module_driver_strength.csv",
    "* Tables/Modules/module_comparison_level_long.csv",
    "* Tables/Modules/module_dominant_direction.csv / module_directional_bias_pvalues.csv",
    "* Plots/Modules/module_localization_<DB>.svg",
    "* Plots/Modules/module_localization_<DB>_by_comparison_class_labeled.svg",
    "* Plots/Modules/module_<MODULE>_<DB>_comparison_matrix.svg",
    "* Plots/Modules/module_drivers_susceptibility_<DB>.svg"
  )
  cat(paste0(md_mod, collapse="\n"),
      file = file.path(out_dir, "CONSOLIDATED_SUMMARY.md"),
      append = TRUE)
}, silent = TRUE)

# =========================
# 16x) Enriched summary
# =========================

# 16a) GO_BP exemplar terms per theme (top by best FDR, tiebreak by frequency)
exemplars <- gsea_df %>%
  dplyr::filter(db == "GO_BP") %>%
  dplyr::inner_join(dplyr::select(themes_tbl, ID, theme), by = "ID") %>%
  dplyr::group_by(theme, Description) %>%
  dplyr::summarize(n_hits = dplyr::n(),
                   bestFDR = suppressWarnings(min(FDR, na.rm = TRUE)),
                   meanNES = mean(NES, na.rm = TRUE),
                   .groups = "drop") %>%
  dplyr::arrange(bestFDR, dplyr::desc(n_hits)) %>%
  dplyr::group_by(theme) %>%
  dplyr::slice_head(n = 5) %>%
  dplyr::ungroup()
readr::write_csv(exemplars, file.path(dirs$tables_gsea, "GO_BP_theme_exemplar_terms.csv"))

# 16b) Top themes per DB with exemplars
top_themes_ex <- theme_repl %>%
  dplyr::group_by(db) %>%
  dplyr::slice_max(order_by = repl_score, n = 8, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(
    exemplars %>%
      dplyr::group_by(theme) %>%
      dplyr::summarize(exemplar = paste(head(Description, 3), collapse = " | "),
                       .groups = "drop"),
    by = "theme"
  )
readr::write_csv(top_themes_ex, file.path(dirs$tables_gsea, "top_themes_with_exemplars.csv"))

# 16c) Theme replication by comparison class and class-specific exemplars
theme_repl_class <- theme_comp %>%
  dplyr::left_join(comp_class_map, by = "context") %>%
  dplyr::filter(!is.na(comparison_class)) %>%
  dplyr::mutate(sign = sign(meanNES), sig = minFDR <= 0.05) %>%
  dplyr::group_by(db, comparison_class, theme) %>%
  dplyr::summarize(n_comp = dplyr::n(),
                   n_sig  = sum(sig, na.rm = TRUE),
                   n_pos  = sum(sig & sign > 0, na.rm = TRUE),
                   n_neg  = sum(sig & sign < 0, na.rm = TRUE),
                   dom_dir = ifelse(n_pos >= n_neg, "up", "down"),
                   repl_score = pmax(n_pos, n_neg),
                   .groups = "drop") %>%
  dplyr::arrange(db, comparison_class, dplyr::desc(repl_score))
readr::write_csv(theme_repl_class, file.path(dirs$tables_class, "theme_replication_by_comparison_class.csv"))

top_themes_by_class <- theme_repl_class %>%
  dplyr::group_by(db, comparison_class) %>%
  dplyr::slice_max(order_by = repl_score, n = 8, with_ties = FALSE) %>%
  dplyr::ungroup()

exemplars_by_class <- gsea_df %>%
  dplyr::inner_join(dplyr::select(themes_tbl, ID, theme), by = "ID") %>%
  dplyr::inner_join(comp_class_map, by = "context") %>%
  dplyr::filter(db == "GO_BP", !is.na(comparison_class)) %>%
  dplyr::group_by(db, comparison_class, theme, Description) %>%
  dplyr::summarize(n_hits = dplyr::n(),
                   bestFDR = suppressWarnings(min(FDR, na.rm = TRUE)),
                   meanNES = mean(NES, na.rm = TRUE),
                   .groups = "drop") %>%
  dplyr::arrange(bestFDR, dplyr::desc(n_hits)) %>%
  dplyr::group_by(db, comparison_class, theme) %>%
  dplyr::slice_head(n = 5) %>%
  dplyr::ungroup()

top_themes_by_class_ex <- top_themes_by_class %>%
  dplyr::left_join(
    exemplars_by_class %>%
      dplyr::group_by(db, comparison_class, theme) %>%
      dplyr::summarize(exemplar = paste(head(Description, 3), collapse = " | "),
                       .groups = "drop"),
    by = c("db", "comparison_class", "theme")
  )
readr::write_csv(top_themes_by_class_ex, file.path(dirs$tables_class, "top_themes_by_class_with_exemplars.csv"))

# 16d) ORA top replicated terms overall and by class (if present)
if (nrow(ora_df)) {
  ora_top <- ora_df %>%
    dplyr::group_by(db, ID, Description) %>%
    dplyr::summarize(n_comp = dplyr::n(),
                     n_sig = sum(Padj <= 0.05, na.rm = TRUE),
                     meanPadj = mean(Padj, na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(n_sig), meanPadj) %>%
    dplyr::group_by(db) %>%
    dplyr::slice_head(n = 10) %>%
    dplyr::ungroup()
  readr::write_csv(ora_top, file.path(dirs$tables_ora, "ORA_top_terms_per_db.csv"))

  ora_top_by_class <- ora_df %>%
    dplyr::inner_join(comp_class_map, by = "context") %>%
    dplyr::filter(!is.na(comparison_class)) %>%
    dplyr::group_by(db, comparison_class, ID, Description) %>%
    dplyr::summarize(n_comp = dplyr::n(),
                     n_sig = sum(Padj <= 0.05, na.rm = TRUE),
                     meanPadj = mean(Padj, na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(n_sig), meanPadj) %>%
    dplyr::group_by(db, comparison_class) %>%
    dplyr::slice_head(n = 10) %>%
    dplyr::ungroup()
  readr::write_csv(ora_top_by_class, file.path(dirs$tables_class, "ORA_top_terms_per_db_by_class.csv"))
}

# 16e) Coverage / QA counters
qa_counts <- list(
  gsea_by_db = as.data.frame(table(gsea_df$db)),
  n_themes   = length(unique(themes_tbl$theme)),
  n_comp     = length(unique(theme_comp$context))
)
if (exists("na_diag")) {
  qa_counts$na_causes <- as.data.frame(table(na_diag$cause, useNA = "ifany"))
}
saveRDS(qa_counts, file.path(dirs$tables_qa, "coverage_counts.rds"))

# 16f) Append enriched Markdown to CONSOLIDATED_SUMMARY.md
classes <- c("sus vs res","sus vs con","res vs con")
md_more <- c(
  "",
  "## Highlights (expanded)",
  paste0("* Databases: ", paste(unique(gsea_df$db), collapse = ", ")),
  paste0("* Comparisons: ", qa_counts$n_comp),
  paste0("* Themes: ", qa_counts$n_themes),
  "",
  "### Top themes per DB"
)

# Per-DB top themes with exemplars
for (dbn in unique(top_themes_ex$db)) {
  d <- top_themes_ex %>% dplyr::filter(db == dbn)
  md_more <- c(md_more, paste0("- ", dbn, ":"))
  if (nrow(d)) {
    lines <- paste0("  - ", d$theme,
                    " (repl=", d$repl_score,
                    ", dir=", d$dom_dir, "): ",
                    dplyr::coalesce(d$exemplar, "—"))
    md_more <- c(md_more, lines)
  }
}

# Per-class top themes with exemplars
md_more <- c(md_more, "", "### Top themes by comparison class")
if (nrow(top_themes_by_class_ex)) {
  for (dbn in unique(top_themes_by_class_ex$db)) {
    for (cl in classes) {
      d <- top_themes_by_class_ex %>% dplyr::filter(db == dbn, comparison_class == cl)
      if (!nrow(d)) next
      md_more <- c(md_more, paste0("- ", dbn, " • ", cl, ":"))
      lines <- paste0("  - ", d$theme,
                      " (repl=", d$repl_score,
                      ", dir=", d$dom_dir, "): ",
                      dplyr::coalesce(d$exemplar, "—"))
      md_more <- c(md_more, lines)
    }
  }
}

# ORA sections if present
if (exists("ora_top") && nrow(ora_df)) {
  md_more <- c(md_more, "", "### ORA replicated terms (top 10 per DB)")
  for (dbn in unique(ora_top$db)) {
    d <- ora_top %>% dplyr::filter(db == dbn)
    md_more <- c(md_more, paste0("- ", dbn, ":"))
    if (nrow(d)) {
      lines <- paste0("  - ", d$Description, " (n_sig=", d$n_sig,
                      ", meanPadj=", sprintf("%.3g", d$meanPadj), ")")
      md_more <- c(md_more, lines)
    }
  }
}

if (exists("ora_top_by_class")) {
  md_more <- c(md_more, "", "### ORA by comparison class (top 10)")
  for (dbn in unique(ora_top_by_class$db)) {
    for (cl in classes) {
      d <- ora_top_by_class %>% dplyr::filter(db == dbn, comparison_class == cl)
      if (!nrow(d)) next
      md_more <- c(md_more, paste0("- ", dbn, " • ", cl, ":"))
      lines <- paste0("  - ", d$Description, " (n_sig=", d$n_sig,
                      ", meanPadj=", sprintf("%.3g", d$meanPadj), ")")
      md_more <- c(md_more, lines)
    }
  }
}

cat(paste0(md_more, collapse = "\n"),
    file = file.path(out_dir, "CONSOLIDATED_SUMMARY.md"),
    append = TRUE)

# 16g) Optional JSON bundle for dashboards
if (requireNamespace("jsonlite", quietly = TRUE)) {
  bundle <- list(
    top_themes       = top_themes_ex,
    exemplars        = exemplars,
    top_by_class     = top_themes_by_class_ex,
    exemplars_by_cls = exemplars_by_class,
    ora_top          = if (exists("ora_top")) ora_top else NULL,
    ora_top_by_class = if (exists("ora_top_by_class")) ora_top_by_class else NULL,
    coverage         = qa_counts
  )
  jsonlite::write_json(bundle,
    path = file.path(out_dir, "bundle_summary.json"),
    pretty = TRUE, auto_unbox = TRUE)
}

# Close any remaining open braces from incomplete blocks above
}

message("Consolidation complete. Enhanced summary with exemplars and class-specific analysis.")
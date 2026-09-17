source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "wgcna_label_adjudication_utils.R"))
source(repo_path("R", "wgcna_label_activation_utils.R"))

# ---------------------------------------------------------------- fixtures

wal_fx_entities <- function() {
  data.frame(
    dataset = "d", level = c("module", "module", "supermodule"),
    entity_id = c("WGCNA_m01", "WGCNA_m02", "SM01"), stringsAsFactors = FALSE)
}
wal_fx_stage07 <- function() {
  data.frame(dataset = "d", level = c("module", "module", "supermodule"),
             entity_id = c("WGCNA_m01", "WGCNA_m02", "SM01"),
             canonical_biological_label = c("s7 one", "s7 two", "s7 sm"),
             stringsAsFactors = FALSE)
}
wal_fx_stage01 <- function() {
  data.frame(dataset = "d", level = "module",
             entity_id = c("WGCNA_m01", "WGCNA_m02"),
             ModuleLabel_Final = c("s1 one", "s1 two"), stringsAsFactors = FALSE)
}
wal_fx_reviewed <- function(status = "reviewed", reviewer = "A Human") {
  data.frame(dataset = "d", level = c("module", "module", "supermodule"),
             entity_id = c("WGCNA_m01", "WGCNA_m02", "SM01"),
             reviewed_biological_label = c("rev one", "rev two", "rev sm"),
             confidence = "high", adjudication_status = status,
             reviewer = reviewer, stringsAsFactors = FALSE)
}

# =====================================================================
# ACTIVE vs PROPOSED
# =====================================================================

testthat::test_that("a proposal is never treated as an active reviewed registry", {
  proposed <- wal_fx_reviewed(status = "proposed", reviewer = NA_character_)
  testthat::expect_false(wal_is_active_reviewed(proposed))
  # status alone is not enough
  testthat::expect_false(wal_is_active_reviewed(
    wal_fx_reviewed(status = "reviewed", reviewer = NA_character_)))
  # a reviewer alone is not enough
  testthat::expect_false(wal_is_active_reviewed(
    wal_fx_reviewed(status = "proposed", reviewer = "A Human")))
  # placeholder reviewer strings do not count
  for (r in c("", "  ", "NA", "none")) {
    testthat::expect_false(wal_is_active_reviewed(
      wal_fx_reviewed(status = "reviewed", reviewer = r)), info = r)
  }
  # both together, with a real label, do
  testthat::expect_true(wal_is_active_reviewed(wal_fx_reviewed()))
  # an empty label blocks activation even when status and reviewer are fine
  blank <- wal_fx_reviewed(); blank$reviewed_biological_label[2] <- ""
  testthat::expect_false(wal_is_active_reviewed(blank))
})

testthat::test_that("the generated proposal files on disk are inactive", {
  for (ds in c("neuron_neuropil", "neuron_soma", "microglia")) {
    p <- path_results("reviewer_audit", "wgcna_label_adjudication",
                      paste0("proposed_", ds, "_reviewed_labels.csv"))
    if (!file.exists(p)) next
    r <- utils::read.csv(p, stringsAsFactors = FALSE)
    testthat::expect_false(wal_is_active_reviewed(r), info = ds)
  }
  # Part-29 finalization: neuron_neuropil now has a registry, but it must
  # contain ONLY the single adjudicated module and must never assert a
  # cell-type identity. neuron_soma still has none. This keeps the original
  # guard - no bulk auto-activation - while allowing the reviewed entry.
  testthat::expect_false(
    file.exists(repo_path("config", "wgcna_labels", "neuron_soma.csv")))
  np <- repo_path("config", "wgcna_labels", "neuron_neuropil.csv")
  if (file.exists(np)) {
    r <- utils::read.csv(np, stringsAsFactors = FALSE)
    testthat::expect_identical(nrow(r), 1L)
    testthat::expect_identical(r$entity_id, "WGCNA_m11")
    testthat::expect_identical(r$adjudication_status, "reviewed")
    testthat::expect_true(grepl("myelin-associated",
                                r$reviewed_biological_label, fixed = TRUE))
    # a co-abundance module may be named for its protein composition, never
    # for a cell type it cannot establish
    testthat::expect_false(grepl("oligodendrocyte",
                                 r$reviewed_biological_label, ignore.case = TRUE))
    testthat::expect_false(grepl("oligodendrocyte",
                                 r$reviewed_short_label, ignore.case = TRUE))
  }
})

# =====================================================================
# PRECEDENCE
# =====================================================================

testthat::test_that("precedence is active-reviewed > Stage07 > Stage01 > id", {
  ent <- wal_fx_entities()

  # 1. active reviewed wins over everything
  r <- resolve_wgcna_display_label(ent, wal_fx_reviewed(), wal_fx_stage07(), wal_fx_stage01())
  testthat::expect_identical(r$canonical_display_label, c("rev one", "rev two", "rev sm"))
  testthat::expect_true(all(r$canonical_label_source == "active_reviewed_registry"))

  # 2. a PROPOSED registry is ignored, so Stage07 wins
  r2 <- resolve_wgcna_display_label(
    ent, wal_fx_reviewed(status = "proposed", reviewer = NA_character_),
    wal_fx_stage07(), wal_fx_stage01())
  testthat::expect_identical(r2$canonical_display_label, c("s7 one", "s7 two", "s7 sm"))
  testthat::expect_true(all(r2$canonical_label_source == "stage07_canonical"))

  # 3. no Stage07 -> Stage01 fallback for modules, id fallback for the supermodule
  r3 <- resolve_wgcna_display_label(ent, NULL, NULL, wal_fx_stage01())
  testthat::expect_identical(r3$canonical_display_label, c("s1 one", "s1 two", "SM01"))
  testthat::expect_identical(r3$canonical_label_source,
                             c("stage01_fallback", "stage01_fallback", "entity_id_fallback"))

  # 4. nothing at all -> never empty
  r4 <- resolve_wgcna_display_label(ent, NULL, NULL, NULL)
  testthat::expect_identical(r4$canonical_display_label, ent$entity_id)
  testthat::expect_true(all(nzchar(r4$canonical_display_label)))
})

testthat::test_that("provenance labels are retained, never discarded", {
  r <- resolve_wgcna_display_label(wal_fx_entities(), wal_fx_reviewed(),
                                   wal_fx_stage07(), wal_fx_stage01())
  testthat::expect_true(all(c("Stage01_ModuleLabel_Final", "Stage07_label",
                              "reviewed_label") %in% names(r)))
  # the superseded labels are still there
  testthat::expect_identical(r$Stage07_label, c("s7 one", "s7 two", "s7 sm"))
  testthat::expect_identical(r$Stage01_ModuleLabel_Final,
                             c("s1 one", "s1 two", NA_character_))
})

testthat::test_that("one canonical display label per dataset and entity", {
  r <- resolve_wgcna_display_label(wal_fx_entities(), NULL, wal_fx_stage07(), NULL)
  testthat::expect_silent(wal_assert_one_label_per_entity(r))
  # the same entity carrying two labels must be rejected
  bad <- rbind(r, r[1, ])
  bad$canonical_display_label[nrow(bad)] <- "a different name"
  testthat::expect_error(wal_assert_one_label_per_entity(bad),
                         "more than one display label")
  # the same label repeated across rows is fine
  ok <- rbind(r, r[1, ])
  testthat::expect_silent(wal_assert_one_label_per_entity(ok))
})

# =====================================================================
# HIERARCHY REDUNDANCY
# =====================================================================

testthat::test_that("hierarchy redundancy is detected and deterministic", {
  mods <- data.frame(
    dataset = "d", entity_id = c("m1", "m2", "m3"),
    parent_entity_id = "SM01",
    proposed_final_label = c("shared name", "shared name", "distinct name"),
    stringsAsFactors = FALSE)
  sups <- data.frame(dataset = "d", entity_id = "SM01",
                     proposed_final_label = "shared name", stringsAsFactors = FALSE)
  w <- wal_hierarchy_redundancy(mods, sups)
  testthat::expect_false(is.na(w))
  testthat::expect_match(w, "matches 2 of 3 member modules", fixed = TRUE)
  # deterministic
  testthat::expect_identical(w, wal_hierarchy_redundancy(mods, sups))

  # distinct member names clear the parent-match warning
  mods2 <- mods
  mods2$proposed_final_label <- c("a", "b", "c")
  testthat::expect_true(is.na(wal_hierarchy_redundancy(mods2, sups)))

  # a singleton is never flagged
  single <- data.frame(dataset = "d", entity_id = "m9", parent_entity_id = "SM09",
                       proposed_final_label = "same", stringsAsFactors = FALSE)
  sup9 <- data.frame(dataset = "d", entity_id = "SM09",
                     proposed_final_label = "same", stringsAsFactors = FALSE)
  testthat::expect_true(is.na(wal_hierarchy_redundancy(single, sup9)))

  # members sharing a name with each other are flagged even if the parent differs
  sups_diff <- sups; sups_diff$proposed_final_label <- "umbrella"
  w2 <- wal_hierarchy_redundancy(mods, sups_diff)
  testthat::expect_match(w2, "shared by multiple member modules", fixed = TRUE)
})

testthat::test_that("label matching ignores case and punctuation only", {
  mods <- data.frame(dataset = "d", entity_id = c("m1", "m2"),
                     parent_entity_id = "SM01",
                     proposed_final_label = c("Synaptic Organization / Signalling",
                                              "synaptic organization - signalling"),
                     stringsAsFactors = FALSE)
  sups <- data.frame(dataset = "d", entity_id = "SM01",
                     proposed_final_label = "synaptic organization signalling",
                     stringsAsFactors = FALSE)
  testthat::expect_false(is.na(wal_hierarchy_redundancy(mods, sups)))
})

# =====================================================================
# ACTIVATION RECOMMENDATION
# =====================================================================

testthat::test_that("low confidence and mixed classes are never activation-ready", {
  testthat::expect_true(wal_recommend_activation("high", "high_confidence_coherent", "KEEP"))
  testthat::expect_true(wal_recommend_activation("moderate", "coherent_but_wording_poor", "REFINE_WORDING"))
  testthat::expect_false(wal_recommend_activation("low", "high_confidence_coherent", "KEEP"))
  testthat::expect_false(wal_recommend_activation("high", "mixed_biology", "REFINE_WORDING"))
  testthat::expect_false(wal_recommend_activation("high", "unresolved", "REFINE_WORDING"))
  testthat::expect_false(wal_recommend_activation("high", "peripheral_enrichment_only", "RENAME"))
  testthat::expect_false(wal_recommend_activation("high", "high_confidence_coherent", "MIXED"))
  testthat::expect_false(wal_recommend_activation("high", "high_confidence_coherent", "UNRESOLVED"))
})

testthat::test_that("the approval table never fabricates a human decision", {
  ent <- data.frame(dataset = "d", level = "module", entity_id = "m1",
                    confidence = "high", evidence_class = "high_confidence_coherent",
                    adjudication_action = "KEEP", stringsAsFactors = FALSE)
  tab <- wal_build_approval_table(ent)
  testthat::expect_true(all(wal_approval_columns() %in% names(tab)))
  testthat::expect_true(all(is.na(tab$human_decision)))
  testthat::expect_true(tab$recommended_for_activation)
  # mixed/unresolved rows accepted as legitimate values
  ent2 <- ent; ent2$evidence_class <- "mixed_biology"
  ent2$adjudication_action <- "MIXED"; ent2$confidence <- "low"
  tab2 <- wal_build_approval_table(ent2)
  testthat::expect_false(tab2$recommended_for_activation)
  testthat::expect_true(is.na(tab2$human_decision))
})

# =====================================================================
# CURATED ADJUDICATIONS (structure, not biology)
# =====================================================================

testthat::test_that("curated adjudications are well formed and phenotype-blind", {
  cur <- wal_curated_adjudications()
  testthat::expect_true(all(c("dataset", "level", "entity_id",
                              "proposed_final_label", "confidence", "rationale")
                            %in% names(cur)))
  testthat::expect_identical(anyDuplicated(paste(cur$dataset, cur$entity_id)), 0L)
  testthat::expect_true(all(cur$confidence %in% c("high", "moderate", "low")))
  testthat::expect_true(all(nzchar(cur$rationale)))
  testthat::expect_true(all(grepl("^WGCNA_m[0-9]{2}$", cur$entity_id)))

  # no rationale text may cite phenotype evidence
  for (pat in c("\\bSUS\\b", "\\bRES\\b", "susceptib", "resilien",
                "differential abundance", "candidate tier", "group effect",
                "log2FC", "\\bDAP", "fold change")) {
    # note: a bare "FDR" is allowed - GO enrichment FDR is phenotype-blind
    # evidence, unlike a differential-abundance FDR.
    testthat::expect_false(any(grepl(pat, cur$rationale)), info = pat)
  }
  testthat::expect_silent(wla_assert_phenotype_blind(cur, "curated adjudications"))
  sup <- wal_curated_supermodule_adjudications()
  testthat::expect_true(all(sup$level == "supermodule"))
  testthat::expect_identical(anyDuplicated(paste(sup$dataset, sup$entity_id)), 0L)
})

# =====================================================================
# CLAIMS-TABLE CONTRACT
# =====================================================================

testthat::test_that("the claims-table resolver is a no-op without entity columns", {
  x <- data.frame(dataset = "d", other = 1, stringsAsFactors = FALSE)
  testthat::expect_identical(attach_canonical_wgcna_display_label(x), x)
  testthat::expect_silent(assert_one_canonical_label_per_wgcna_entity(x))
})

testthat::test_that("the claims-table resolver accepts both entity column spellings", {
  a <- data.frame(wgcna_entity_id = "m1", wgcna_level = "module", stringsAsFactors = FALSE)
  b <- data.frame(module_or_supermodule_id = "m1", entity_level = "module", stringsAsFactors = FALSE)
  testthat::expect_identical(wal_claims_entity_columns(a)$id, "wgcna_entity_id")
  testthat::expect_identical(wal_claims_entity_columns(b)$id, "module_or_supermodule_id")
})

testthat::test_that("the real claims table carries one label per entity with provenance", {
  p <- path_results("tables", "biological_claims_table.csv")
  testthat::skip_if_not(file.exists(p), "claims table has not been generated")
  x <- utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)

  testthat::expect_true(all(c("canonical_display_label", "canonical_label_source",
                              "Stage01_ModuleLabel_Final", "Stage07_label",
                              "Stage06_label") %in% names(x)))
  testthat::expect_silent(assert_one_canonical_label_per_wgcna_entity(x))

  cols <- wal_claims_entity_columns(x)
  e <- x[!is.na(x$canonical_display_label) & nzchar(x$canonical_display_label), , drop = FALSE]
  testthat::skip_if(nrow(e) == 0L)
  # every source value is from the declared vocabulary
  testthat::expect_true(all(e$canonical_label_source %in% wal_label_sources()))
  # and the invariant holds directly
  k <- paste(e$dataset, e[[cols$level]], e[[cols$id]])
  testthat::expect_equal(
    max(tapply(e$canonical_display_label, k, function(v) length(unique(v)))), 1)
})

testthat::test_that("raw Stage-01 failures never surface as canonical displays", {
  p <- path_results("tables", "biological_claims_table.csv")
  testthat::skip_if_not(file.exists(p))
  x <- utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
  # these are known Stage-01 microglia labels that the reviewed registry replaces
  for (bad in c("Zona Pellucida", "Binding Sperm", "Skin Development")) {
    testthat::expect_false(any(grepl(bad, x$canonical_display_label, ignore.case = TRUE)),
                           info = bad)
  }
})

# =====================================================================
# REAL APPROVAL TABLE
# =====================================================================

wal_approval_path <- function(f) {
  path_results("reviewer_audit", "wgcna_label_approval", f)
}

testthat::test_that("the approval table covers every entity exactly once", {
  p <- wal_approval_path("WGCNA_final_label_approval_table.csv")
  testthat::skip_if_not(file.exists(p), "approval table has not been generated")
  a <- utils::read.csv(p, stringsAsFactors = FALSE)

  testthat::expect_identical(anyDuplicated(paste(a$dataset, a$level, a$entity_id)), 0L)
  testthat::expect_true(all(is.na(a$human_decision) | a$human_decision == ""))

  for (ds in unique(a$dataset)) {
    ct <- path_results("tables", "06_modules_WGCNA", "identity_contract", ds,
                       "WGCNA_module_supermodule_membership_contract.csv")
    if (!file.exists(ct)) next
    m <- utils::read.csv(ct, stringsAsFactors = FALSE)
    sub <- a[a$dataset == ds, , drop = FALSE]
    # IDs unchanged and complete
    testthat::expect_setequal(sub$entity_id[sub$level == "module"], unique(m$module_id))
    testthat::expect_setequal(sub$entity_id[sub$level == "supermodule"],
                              unique(m$supermodule_id))
  }
  # every proposed label is non-empty, and mixed/unresolved is a legal value
  testthat::expect_true(all(nzchar(a$proposed_final_label)))
  testthat::expect_true(any(grepl("unresolved", a$proposed_final_label)))
})

testthat::test_that("nothing low-confidence or mixed is recommended for activation", {
  p <- wal_approval_path("WGCNA_final_label_approval_table.csv")
  testthat::skip_if_not(file.exists(p))
  a <- utils::read.csv(p, stringsAsFactors = FALSE)
  rec <- a[a$recommended_for_activation %in% TRUE, , drop = FALSE]
  testthat::expect_false(any(rec$confidence == "low"))
  testthat::expect_false(any(grepl("mixed|unresolved|peripheral", rec$evidence_class)))
  testthat::expect_false(any(rec$adjudication_action %in% c("MIXED", "UNRESOLVED")))
  testthat::expect_false(any(grepl("unresolved", rec$proposed_final_label,
                                   ignore.case = TRUE)))
})

testthat::test_that("no phenotype field reaches the approval table", {
  p <- wal_approval_path("WGCNA_final_label_approval_table.csv")
  testthat::skip_if_not(file.exists(p))
  a <- utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
  testthat::expect_silent(wla_assert_phenotype_blind(a, "approval table"))
})

# =====================================================================
# ACTIVATION RULES: defined, never executed
# =====================================================================

testthat::test_that("activation rules are declared and nothing executes them", {
  r <- wal_activation_rules()
  testthat::expect_gt(nrow(r), 0L)
  testthat::expect_true(all(c("step", "rule_id", "rule",
                              "executed_by_this_script") %in% names(r)))
  testthat::expect_identical(r$step, seq_len(nrow(r)))
  testthat::expect_false(any(r$executed_by_this_script))
  testthat::expect_identical(anyDuplicated(r$rule_id), 0L)
  testthat::expect_true(all(nzchar(r$rule)))

  # the rules must name the two things that make a registry active
  live <- paste(r$rule, collapse = " ")
  testthat::expect_match(live, "config/wgcna_labels", fixed = TRUE)
  testthat::expect_match(live, "reviewed", fixed = TRUE)
})

testthat::test_that("the approval script writes the rules but activates nothing", {
  p <- wal_approval_path("WGCNA_label_activation_rules.csv")
  testthat::skip_if_not(file.exists(p), "approval package has not been generated")
  r <- utils::read.csv(p, stringsAsFactors = FALSE)
  testthat::expect_false(any(as.logical(r$executed_by_this_script)))

  # the script must not contain a write into the live registry directory
  src <- readLines(repo_path("analysis/wgcna",
                             "build_module_label_registry.R"), warn = FALSE)
  writes <- grep("write_csv_safe|write[.]csv|writeLines|saveRDS|file[.]copy|file[.]rename",
                 src, value = TRUE)
  testthat::expect_false(any(grepl("config", writes)))
  testthat::expect_false(any(grepl("wgcna_labels", writes)))
})

# =====================================================================
# ACTION / LABEL / CLASS CONSISTENCY
# =====================================================================

testthat::test_that("the action always agrees with the adjudicated label", {
  p <- wal_approval_path("WGCNA_final_label_approval_table.csv")
  testthat::skip_if_not(file.exists(p))
  a <- utils::read.csv(p, stringsAsFactors = FALSE)
  lab <- tolower(trimws(a$proposed_final_label))
  act <- toupper(trimws(a$adjudication_action))

  # a mixed/unresolved label must carry the matching action
  testthat::expect_true(all(act[lab %in% c("mixed / unresolved", "mixed/unresolved")] == "MIXED"))
  testthat::expect_true(all(act[lab %in% "unresolved"] == "UNRESOLVED"))
  # and the inverse: a MIXED/UNRESOLVED action must not sit on a named label
  named <- !lab %in% c("mixed / unresolved", "mixed/unresolved", "unresolved")
  testthat::expect_false(any(act[named] %in% c("MIXED", "UNRESOLVED")))
})

testthat::test_that("a manual override keeps the algorithmic class as provenance", {
  p <- wal_approval_path("WGCNA_final_label_approval_table.csv")
  testthat::skip_if_not(file.exists(p))
  a <- utils::read.csv(p, stringsAsFactors = FALSE)
  testthat::expect_true(all(c("algorithmic_evidence_class",
                              "evidence_class_source") %in% names(a)))
  testthat::expect_true(all(a$evidence_class_source %in%
                              c("manual_adjudication", "algorithmic")))
  # where a human overrode the algorithm, the original class is still recorded
  ov <- a[a$evidence_class == "manually_adjudicated_coherent", , drop = FALSE]
  if (nrow(ov)) {
    testthat::expect_true(all(nzchar(ov$algorithmic_evidence_class)))
    testthat::expect_true(all(ov$evidence_class_source == "manual_adjudication"))
    # an override must never be a silent upgrade: it has to carry a rationale
    testthat::expect_true(all(nchar(ov$rationale) > 80))
  }
})

testthat::test_that("every curated entity carries a substantive rationale", {
  cur <- wal_curated_adjudications()
  testthat::expect_true(all(nchar(cur$rationale) > 80))
  # a label that differs from mixed/unresolved must cite at least one p-value
  named <- !tolower(cur$proposed_final_label) %in%
    c("mixed / unresolved", "unresolved")
  testthat::expect_true(all(grepl("p=", cur$rationale[named])))
})

# =====================================================================
# AN ACTIVE REVIEWED LABEL IS NEVER PROPOSED FOR REPLACEMENT
# =====================================================================

testthat::test_that("KEEP_ACTIVE_REVIEWED can never be recommended", {
  testthat::expect_false(
    wal_recommend_activation("high", "high_confidence_coherent", "KEEP_ACTIVE_REVIEWED"))
  testthat::expect_false(
    wal_recommend_activation("moderate", "manually_adjudicated_coherent", "KEEP_ACTIVE_REVIEWED"))
})

testthat::test_that("a dataset with an active registry proposes no replacements", {
  p <- wal_approval_path("WGCNA_final_label_approval_table.csv")
  testthat::skip_if_not(file.exists(p))
  a <- utils::read.csv(p, stringsAsFactors = FALSE)
  rev <- a[a$adjudication_action %in% "KEEP_ACTIVE_REVIEWED", , drop = FALSE]
  testthat::skip_if(nrow(rev) == 0L)

  # the proposal must BE the reviewed label, not a replacement for it
  testthat::expect_identical(rev$proposed_final_label, rev$current_active_label)
  testthat::expect_false(any(rev$recommended_for_activation %in% TRUE))
  # and the superseded automatic proposal is kept, not discarded
  testthat::expect_true("automatic_proposal_label" %in% names(a))
  testthat::expect_true(any(nzchar(rev$automatic_proposal_label)))
})

# =====================================================================
# BP BLIND-SPOT DETECTOR
# =====================================================================

testthat::test_that("the BP blind-spot detector is advisory and deterministic", {
  rows <- data.frame(
    dataset = "d", entity_id = c("m1", "m2", "m3"),
    proposed_final_label = c("mixed / unresolved", "a real name", "unresolved"),
    algorithmic_evidence_class = c("mixed_biology", "coherent", "unresolved"),
    stringsAsFactors = FALSE)
  go <- data.frame(
    dataset = "d", ModuleID = c("m1", "m2", "m3"),
    Ontology = "CC", Description = c("some complex", "other complex", "third"),
    p.adjust = c(1e-20, 1e-20, 0.5), stringsAsFactors = FALSE)
  w <- wal_bp_blind_spot(rows, go)
  # flagged: mixed with strong CC
  testthat::expect_false(is.na(w[1]))
  testthat::expect_match(w[1], "BP-only theme layer", fixed = TRUE)
  # not flagged: already has a real name
  testthat::expect_true(is.na(w[2]))
  # not flagged: unresolved but the CC evidence is not significant
  testthat::expect_true(is.na(w[3]))
  # deterministic, and never alters a label
  testthat::expect_identical(w, wal_bp_blind_spot(rows, go))
  testthat::expect_identical(rows$proposed_final_label,
                             c("mixed / unresolved", "a real name", "unresolved"))
  # no GO table at all is a safe no-op
  testthat::expect_true(all(is.na(wal_bp_blind_spot(rows, NULL))))
})

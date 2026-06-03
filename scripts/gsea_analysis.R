#!/usr/bin/env Rscript
# ==============================================================================
# title: gsea_analysis.R
# project: Ketamine Proteomics Analysis Project
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-05-28
# last modified: 2026-05-28
#
# purpose:
#   Run Gene Set Enrichment Analysis (GSEA) on the full ranked list of
#   quantitative proteins (~846), testing whether members of GO Biological
#   Process, KEGG, and Reactome gene sets cluster toward the extremes of the
#   ranking. Complements ORA by capturing subthreshold pathway-level signal
#   that the ORA significance cutoff discards (Section 4.1, decision 2026-03-17).
#
#   Two ranking metrics are run in parallel as a sensitivity analysis:
#     (1) log2 fold change  (primary; matches project notes framing)
#     (2) signed -log10(p)  (combines effect size + confidence; supplementary)
#
#   Two min_size values are run in parallel (default 10, 15) so the synaptic
#   priming / SNARE-assembly gene sets (typical size 10-12 in mouse) appear
#   in the primary (min_size=10) run and the stricter MSigDB convention is
#   available as a sensitivity check (min_size=15). max_size defaults to 250
#   because a 500-gene set would cover >50% of the 846-protein ranked list and
#   admit very broad terms whose enrichment is driven by ranked-list coverage
#   rather than biological specificity.
#
#   Presence/absence proteins (PD placeholders FC=100 and FC=0.01) are
#   excluded by design via a positive whitelist on the input CSV's `category`
#   column (must equal "quantitative"). A defensive |log2FC| assertion
#   guarantees no placeholders survive the category filter; see [[ketamine-
#   data-limitations]] for full rationale.
#
# inputs:
#   - results/quantitative/all_quantitative_proteins.csv (default)
#     Required columns: 'Gene Symbol', 'log2_fold_change',
#                       'Abundance Ratio Adj. P-Value: (ketamine) / (control)',
#                       'category' (must be "quantitative").
#
# outputs:
#   results/gsea/
#   |-- ranked_list_log2fc.rnk             (gene symbol -> log2FC, descending)
#   |-- ranked_list_signed_logp.rnk        (gene symbol -> signed -log10(p))
#   |-- unmapped_symbols.txt               (symbols not mapped to ENTREZ)
#   |-- kegg_version.txt                   (KEGG release captured at run time)
#   |-- min10/                             (primary; recovers small synaptic sets)
#   |   |-- log2fc/
#   |   |   |-- go_bp_gsea_results.csv     (with convergence_key column)
#   |   |   |-- go_bp_dotplot.{pdf,png}
#   |   |   |-- go_bp_ridgeplot.{pdf,png}
#   |   |   |-- go_bp_top_terms_running_score.{pdf,png}
#   |   |   |-- kegg_*.{...}
#   |   |   `-- reactome_*.{...}
#   |   `-- signed_logp/                   (same layout)
#   |-- min15/                             (sensitivity; stricter MSigDB cutoff)
#   |   `-- ... (same layout)
#   |-- ranking_metric_concordance.csv     (computed on primary min_size only)
#   |-- gsea_analysis_report.txt
#   `-- sessionInfo.txt
#
# usage:
#   Rscript gsea_analysis.R \
#       --input ../results/quantitative/all_quantitative_proteins.csv \
#       --output_dir ../results/gsea \
#       --log_dir ../logs \
#       --fdr 0.05 \
#       --min_sizes 10,15 \
#       --max_size 250 \
#       --n_perm 10000 \
#       --seed 42
#
#   copy/paste: Rscript gsea_analysis.R --input ../results/quantitative/all_quantitative_proteins.csv --output_dir ../results/gsea --log_dir ../logs --fdr 0.05 --min_sizes 10,15 --max_size 250 --n_perm 10000 --seed 42
#
# dependencies (CRAN + Bioconductor):
#   install.packages(c("tidyverse", "optparse", "cowplot", "ggridges"))
#   BiocManager::install(c("clusterProfiler", "org.Mm.eg.db", "fgsea",
#                          "ReactomePA", "enrichplot", "DOSE",
#                          "AnnotationDbi", "KEGGREST", "reactome.db"))
#   ggridges must be loaded with `library(ggridges)`, not just installed,
#   because enrichplot's ridgeplot uses NSE that requires ggridges on the
#   search path.
#
# notes:
#   - gseKEGG queries the KEGG REST API; internet required. KEGGREST::keggInfo
#     captures the KEGG release identifier at run time to make any drift
#     between reruns auditable. Pinning is not possible: KEGG does not offer
#     versioned queries through its REST API.
#   - Convergence-key column: each *_gsea_results.csv includes a
#     `convergence_key` column normalized for joining against the existing
#     g:Profiler ORA outputs in results/pathway_analysis/. Format:
#       GO:BP    -> "GO:0001234" (unchanged from clusterProfiler default)
#       KEGG     -> "KEGG:04020"  (organism prefix stripped, "KEGG:" added)
#       Reactome -> "Reactome:R-MMU-12345" (prefix added)
#   - Duplicate gene symbols are collapsed once by max |log2FC| before either
#     ranking metric is computed, so the two rankings differ only in metric
#     and not in which protein row represents each gene. This is required
#     for the metric-concordance analysis to be interpretable.
#   - Reproducibility under parallelism is not guaranteed by set.seed() alone
#     if BiocParallel is configured with a multicore backend. If exact bitwise
#     reproducibility across hosts becomes needed, add
#     BiocParallel::register(SerialParam(RNGseed = opt$seed)) near the top.
# ==============================================================================


# --- Library Loading ---
suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(fgsea)
  library(ReactomePA)
  library(enrichplot)
  library(DOSE)
  library(cowplot)
  library(AnnotationDbi)
  library(KEGGREST)
  library(reactome.db)
  # ggridges must be ATTACHED (not just installed) because enrichplot's
  # ridgeplot relies on NSE that resolves a `selected` variable from
  # ggridges's search-path scope. With ggridges merely installed but not
  # attached, ridgeplot errors with `object 'selected' not found`.
  library(ggridges)
})


# --- Configuration & Constants ---

# Column names in input CSV
GENE_COL           <- "Gene Symbol"
LOG2FC_COL         <- "log2_fold_change"
PVAL_COL           <- "Abundance Ratio Adj. P-Value: (ketamine) / (control)"
CATEGORY_COL       <- "category"
CATEGORY_WHITELIST <- "quantitative"   # only rows with this label enter GSEA

# Species annotation
ORGANISM_DB       <- org.Mm.eg.db
KEGG_ORGANISM     <- "mmu"      # KEGG three-letter code for mouse
REACTOME_ORGANISM <- "mouse"    # ReactomePA expects 'mouse'

# Defensive assertion only. The primary placeholder safeguard is
# CATEGORY_WHITELIST above; this cutoff guarantees no FC=100 / FC=0.01
# placeholders survive should the upstream category labels ever change.
PLACEHOLDER_LOG2FC_CUTOFF <- 6.0

# Visual constants (used in figures elsewhere; kept for cross-script consistency)
COLOR_UP   <- "#D62728"
COLOR_DOWN <- "#1F77B4"

# Null/NA-coalescing helper for the report writer
`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a


# --- CLI Argument Parsing ---

option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              default = "../results/quantitative/all_quantitative_proteins.csv",
              help = "Input CSV (quantitative proteins ranked list source) [%default]"),
  make_option(c("-o", "--output_dir"), type = "character",
              default = "../results/gsea",
              help = "Output directory [%default]"),
  make_option(c("-l", "--log_dir"), type = "character",
              default = "../logs",
              help = "Log directory [%default]"),
  make_option("--fdr", type = "numeric", default = 0.05,
              help = "FDR threshold for declaring a gene set significant [%default]"),
  make_option("--min_sizes", type = "character", default = "10,15",
              help = "Comma-separated min_size values; first is primary [%default]"),
  make_option("--max_size", type = "integer", default = 250,
              help = "Maximum gene set size [%default]"),
  make_option("--n_perm", type = "integer", default = 10000,
              help = "Initial fgsea permutation count (adaptive) [%default]"),
  make_option("--seed", type = "integer", default = 42,
              help = "Random seed for reproducibility [%default]"),
  make_option("--top_n_plot", type = "integer", default = 20,
              help = "Number of top terms in dotplot/ridgeplot [%default]"),
  make_option("--top_n_running", type = "integer", default = 5,
              help = "Number of top terms per direction in gseaplot2 running score [%default]"),
  make_option("--dpi", type = "integer", default = 600,
              help = "DPI for PNG outputs (Word-embedding target) [%default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Parse --min_sizes from CSV string into integer vector. First value is primary.
opt$min_sizes <- as.integer(trimws(strsplit(opt$min_sizes, ",", fixed = TRUE)[[1]]))
if (any(is.na(opt$min_sizes)) || length(opt$min_sizes) == 0) {
  stop("--min_sizes must be a comma-separated list of integers (e.g., '10,15').")
}


# --- Logging Helpers ---

init_log <- function(log_dir) {
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  log_path <- file.path(log_dir, paste0("gsea_analysis_", ts, ".log"))
  cat(
    strrep("=", 70), "\n",
    "KETAMINE PROTEOMICS ANALYSIS PROJECT - GSEA LOG\n",
    strrep("=", 70), "\n\n",
    sprintf("Log file: %s\n", basename(log_path)),
    sprintf("Started: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sep = "", file = log_path
  )
  log_path
}

log_msg <- function(log_path, msg, console = TRUE) {
  ts   <- format(Sys.time(), "%H:%M:%S")
  line <- sprintf("[%s] %s", ts, msg)
  cat(line, "\n", file = log_path, append = TRUE)
  if (console) message(line)
}


# --- Convergence-Key Normalization ---

# Produce a canonical key per gene set for joining GSEA output against the
# existing g:Profiler ORA output (results/pathway_analysis/). See script
# header for the per-database convention.
normalize_convergence_key <- function(ids, source) {
  switch(source,
    "go_bp"    = ids,
    "kegg"     = sub(paste0("^", KEGG_ORGANISM), "KEGG:", ids),
    "reactome" = paste0("Reactome:", ids),
    ids
  )
}


# --- Data Preparation ---

load_ranked_list <- function(input_path, log_path) {
  # Step 1: load the categorized protein CSV
  log_msg(log_path, sprintf("Loading ranked-list source: %s", input_path))
  df <- read_csv(input_path, show_col_types = FALSE)
  log_msg(log_path, sprintf("  Rows in input: %d", nrow(df)))

  # Step 2: primary placeholder safeguard - whitelist by category column.
  # This rejects any row whose label is not the expected "quantitative",
  # protecting against upstream regenerations that introduce other categories.
  if (!CATEGORY_COL %in% colnames(df)) {
    stop(sprintf("Required column '%s' not found in input. ",  CATEGORY_COL),
         "Cannot apply category whitelist filter.")
  }
  n_pre_cat <- nrow(df)
  df <- df |> filter(.data[[CATEGORY_COL]] == CATEGORY_WHITELIST)
  log_msg(log_path, sprintf("  After category whitelist ('%s' only): %d (excluded %d)",
                            CATEGORY_WHITELIST, nrow(df), n_pre_cat - nrow(df)))

  # Step 3: standardize column names and drop incomplete rows.
  # rename(any_of(named_vec)) is the dplyr idiom for programmatic renaming
  # where the source name is held in a variable; names of the vector are
  # new names, values are old names.
  required <- c(GENE_COL, LOG2FC_COL, PVAL_COL)
  missing_cols <- setdiff(required, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Input is missing required column(s): ",
         paste(missing_cols, collapse = ", "))
  }
  # dplyr::rename is explicit because S4Vectors (a transitive dep of
  # reactome.db / ReactomePA) exports a `rename` generic that otherwise
  # shadows dplyr::rename and does not understand tidyselect helpers.
  df <- df |>
    dplyr::rename(any_of(c(gene   = GENE_COL,
                           log2fc = LOG2FC_COL,
                           pval   = PVAL_COL))) |>
    filter(!is.na(gene), gene != "",
           !is.na(log2fc), is.finite(log2fc),
           !is.na(pval),   is.finite(pval), pval > 0)
  log_msg(log_path, sprintf("  After dropping NA gene/log2FC/p-value: %d", nrow(df)))

  # Step 4: defensive assertion - no PD placeholder should survive the
  # category filter. If this fires, the upstream pipeline produced rows
  # mislabeled as "quantitative" that look like placeholders.
  n_placeholder <- sum(abs(df$log2fc) >= PLACEHOLDER_LOG2FC_CUTOFF, na.rm = TRUE)
  if (n_placeholder > 0) {
    stop(sprintf(
      "Assertion failed: %d rows have |log2FC| >= %.1f after category whitelist. ",
      n_placeholder, PLACEHOLDER_LOG2FC_CUTOFF),
      "This indicates the input contains presence/absence placeholders ",
      "that the category filter did not catch. Inspect the upstream CSV.")
  }

  # Step 5: collapse duplicated gene symbols ONCE by max |log2FC|, then
  # compute both ranking metrics from the surviving rows. This guarantees
  # the two rankings differ only in metric and not in which protein-accession
  # row represents each gene, which is required for metric-concordance to be
  # interpretable as ranking-metric sensitivity rather than row-selection noise.
  df <- df |>
    group_by(gene) |>
    slice_max(abs(log2fc), n = 1, with_ties = FALSE) |>
    ungroup() |>
    mutate(
      rank_log2fc      = log2fc,
      rank_signed_logp = -log10(pval) * sign(log2fc)
    ) |>
    # Guard against any Inf that might arise from p-values rounding to zero
    # downstream of the pval > 0 filter (defensive; should not occur).
    filter(is.finite(rank_log2fc), is.finite(rank_signed_logp))
  log_msg(log_path, sprintf("  Proteins after duplicate collapse: %d", nrow(df)))

  list(
    raw_n = nrow(df),
    rankings = list(
      log2fc      = setNames(df$rank_log2fc,      df$gene),
      signed_logp = setNames(df$rank_signed_logp, df$gene)
    )
  )
}

map_to_entrez <- function(gene_vec, log_path, label,
                          unmapped_out_path = NULL,
                          recovered_out_path = NULL) {
  # Map gene-symbol-keyed ranking values to ENTREZ-ID-keyed ranking values
  # in three steps, recovering signal that a naive bitr() lookup loses:
  #
  #   1. EXPAND multi-symbol aggregate rows. PD packs paralogous proteins
  #      that share peptides into a single row whose Gene Symbol value is
  #      a semicolon-separated list (e.g., "H3c1; H3c10; H3c11; H3c8").
  #      Each component inherits the aggregate row's score. This is the
  #      same value PD reports for each member anyway, but it does add
  #      tied ranks; see Section 8.1.5 of the project notes.
  #
  #   2. PRIMARY lookup SYMBOL -> ENTREZID via org.Mm.eg.db.
  #
  #   3. ALIAS fallback for symbols still unmapped after step 2. The MGI
  #      alias table catches most HGNC-vs-MGI naming differences such as
  #      Ca2 -> Car2, Fh -> Fh1, Gpi -> Gpi1, Txn -> Txn1, Znf706 -> Zfp706.
  #      Recovered mappings are written to a sidecar file for auditing.
  #
  # After mapping, the rankings are collapsed once more by ENTREZID to
  # handle the case where two source rows resolve to the same gene
  # (e.g., the alias and the canonical symbol both appearing).
  expanded <- tibble(orig = names(gene_vec), score = unname(gene_vec)) |>
    tidyr::separate_rows(orig, sep = "\\s*;\\s*") |>
    filter(orig != "")

  n_input  <- length(gene_vec)
  n_aggregate <- sum(grepl(";", names(gene_vec), fixed = TRUE))
  n_expanded  <- nrow(expanded) - n_input
  if (n_aggregate > 0) {
    log_msg(log_path,
            sprintf("  [%s] expanded %d aggregate row(s) into %d additional entries",
                     label, n_aggregate, n_expanded))
  }

  unique_syms <- unique(expanded$orig)

  # Step 2: primary SYMBOL -> ENTREZID
  m_sym <- suppressWarnings(
    bitr(unique_syms, fromType = "SYMBOL", toType = "ENTREZID",
         OrgDb = ORGANISM_DB)
  )
  m_sym <- if (nrow(m_sym) > 0) {
    m_sym |> dplyr::distinct(SYMBOL, .keep_all = TRUE)
  } else {
    tibble(SYMBOL = character(0), ENTREZID = character(0))
  }

  # Step 3: ALIAS fallback for symbols still unmapped
  still_unmapped <- setdiff(unique_syms, m_sym$SYMBOL)
  m_alias <- if (length(still_unmapped) > 0) {
    tmp <- suppressWarnings(
      bitr(still_unmapped, fromType = "ALIAS", toType = "ENTREZID",
           OrgDb = ORGANISM_DB)
    )
    if (nrow(tmp) > 0) {
      tmp |> dplyr::rename(SYMBOL = ALIAS) |>
        dplyr::distinct(SYMBOL, .keep_all = TRUE)
    } else {
      tibble(SYMBOL = character(0), ENTREZID = character(0))
    }
  } else {
    tibble(SYMBOL = character(0), ENTREZID = character(0))
  }

  # Step 4: mt- prefix retry for mitochondrially-encoded genes.
  # PD writes the human-style "Mtco2"; MGI uses "mt-Co2". For symbols still
  # unmapped after ALIAS, try paste0("mt-", first_letter_upper(rest_lower)).
  still_after_alias <- setdiff(still_unmapped, m_alias$SYMBOL)
  m_mt <- if (length(still_after_alias) > 0) {
    mt_candidates <- ifelse(
      grepl("^[Mm][Tt]", still_after_alias),
      paste0("mt-",
             toupper(substr(still_after_alias, 3, 3)),
             tolower(substr(still_after_alias, 4, nchar(still_after_alias)))),
      NA_character_
    )
    keep <- !is.na(mt_candidates)
    if (any(keep)) {
      tmp <- suppressWarnings(
        bitr(mt_candidates[keep], fromType = "SYMBOL", toType = "ENTREZID",
             OrgDb = ORGANISM_DB)
      )
      if (nrow(tmp) > 0) {
        # Map the recovered ENTREZID back to the ORIGINAL PD-style symbol
        # so the downstream join uses the original.
        lookup <- tibble(
          orig_sym = still_after_alias[keep],
          mt_sym   = mt_candidates[keep]
        ) |> inner_join(tmp, by = c("mt_sym" = "SYMBOL"))
        tibble(SYMBOL = lookup$orig_sym, ENTREZID = lookup$ENTREZID)
      } else {
        tibble(SYMBOL = character(0), ENTREZID = character(0))
      }
    } else {
      tibble(SYMBOL = character(0), ENTREZID = character(0))
    }
  } else {
    tibble(SYMBOL = character(0), ENTREZID = character(0))
  }

  log_msg(log_path,
          sprintf("  [%s] mapped via SYMBOL: %d; via ALIAS: %d; via mt- prefix: %d; total unique: %d",
                   label, nrow(m_sym), nrow(m_alias), nrow(m_mt),
                   nrow(m_sym) + nrow(m_alias) + nrow(m_mt)))

  # Audit sidecar: list all symbols recovered by the fallback paths.
  recovered_total <- nrow(m_alias) + nrow(m_mt)
  if (recovered_total > 0 && !is.null(recovered_out_path)) {
    lines <- c(
      "# Symbols recovered via fallback lookups (not via primary SYMBOL key)",
      "# Format: ORIGINAL_SYMBOL\\tRECOVERED_ENTREZID\\tVIA"
    )
    if (nrow(m_alias) > 0) {
      lines <- c(lines,
                 paste(m_alias$SYMBOL, m_alias$ENTREZID, "ALIAS", sep = "\t"))
    }
    if (nrow(m_mt) > 0) {
      lines <- c(lines,
                 paste(m_mt$SYMBOL, m_mt$ENTREZID, "MT_PREFIX", sep = "\t"))
    }
    writeLines(lines, recovered_out_path)
    log_msg(log_path,
            sprintf("    Wrote %d fallback-recovered mappings to %s",
                     recovered_total, basename(recovered_out_path)))
  }

  mapping <- bind_rows(m_sym, m_alias, m_mt)
  unmapped <- setdiff(unique_syms, mapping$SYMBOL)
  if (length(unmapped) > 0 && !is.null(unmapped_out_path)) {
    writeLines(unmapped, unmapped_out_path)
    log_msg(log_path,
            sprintf("    Wrote %d still-unmapped symbols to %s",
                     length(unmapped), basename(unmapped_out_path)))
  }

  # Join expanded ranking against the mapping; collapse to one row per
  # ENTREZID by max |score| in case the same gene is reached via two
  # source symbols (canonical + alias) or via expansion + a separate row.
  ranked_df <- expanded |>
    inner_join(mapping |> dplyr::select(SYMBOL, ENTREZID),
               by = c("orig" = "SYMBOL")) |>
    filter(!is.na(score), is.finite(score)) |>
    group_by(ENTREZID) |>
    slice_max(abs(score), n = 1, with_ties = FALSE) |>
    ungroup() |>
    arrange(desc(score))

  setNames(ranked_df$score, ranked_df$ENTREZID)
}


# --- KEGG Version Capture ---

capture_kegg_version <- function(out_path, log_path) {
  # One-time KEGG release capture. Forensic record only; does not prevent
  # drift between reruns because the KEGG REST API does not support
  # versioned queries.
  v <- tryCatch(KEGGREST::keggInfo("kegg"),
                error = function(e) sprintf("ERROR: %s", conditionMessage(e)))
  writeLines(c(
    sprintf("Captured: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    "KEGGREST::keggInfo('kegg') output:",
    v
  ), out_path)
  log_msg(log_path, sprintf("KEGG release version captured to %s", basename(out_path)))
}


# --- Gene-Set Size Survey ---
# For each database we count gene-set sizes after intersection with the
# ranked-list universe, so the report can show how many sets were excluded
# by min_size and max_size at the chosen cutoffs.

survey_go_bp_sizes <- function(entrez_universe) {
  m <- suppressMessages(AnnotationDbi::select(
    ORGANISM_DB, keys = entrez_universe,
    columns = c("GOALL", "ONTOLOGYALL"), keytype = "ENTREZID"
  ))
  m <- m |> filter(ONTOLOGYALL == "BP", !is.na(GOALL))
  if (nrow(m) == 0) return(integer(0))
  m |> distinct(GOALL, ENTREZID) |> count(GOALL, name = "size") |> pull(size)
}

survey_kegg_sizes <- function(entrez_universe, organism = KEGG_ORGANISM) {
  link <- tryCatch(KEGGREST::keggLink("pathway", organism),
                   error = function(e) NULL)
  if (is.null(link) || length(link) == 0) return(integer(0))
  entrez_part <- sub(paste0("^", organism, ":"), "", names(link))
  in_universe <- entrez_part %in% entrez_universe
  if (!any(in_universe)) return(integer(0))
  as.integer(table(link[in_universe]))
}

survey_reactome_sizes <- function(entrez_universe) {
  r2g <- tryCatch(as.list(reactome.db::reactomePATHID2EXTID),
                  error = function(e) NULL)
  if (is.null(r2g)) return(integer(0))
  r2g_mouse <- r2g[grepl("^R-MMU-", names(r2g))]
  if (length(r2g_mouse) == 0) return(integer(0))
  sapply(r2g_mouse, function(g) sum(g %in% entrez_universe))
}

summarize_size_exclusions <- function(sizes, min_size, max_size) {
  if (length(sizes) == 0) return(tibble(
    n_potential = NA_integer_, n_below_min = NA_integer_,
    n_above_max = NA_integer_, n_within = NA_integer_
  ))
  tibble(
    n_potential = length(sizes),
    n_below_min = sum(sizes < min_size),
    n_above_max = sum(sizes > max_size),
    n_within    = sum(sizes >= min_size & sizes <= max_size)
  )
}


# --- GSEA Runners ---

run_go_bp <- function(ranking, min_size, max_size, n_perm, seed) {
  set.seed(seed)
  gseGO(
    geneList      = ranking,
    OrgDb         = ORGANISM_DB,
    ont           = "BP",
    keyType       = "ENTREZID",
    minGSSize     = min_size,
    maxGSSize     = max_size,
    pvalueCutoff  = 1,
    pAdjustMethod = "BH",
    eps           = 0,
    nPermSimple   = n_perm,
    verbose       = FALSE
  )
}

run_kegg <- function(ranking, min_size, max_size, n_perm, seed) {
  set.seed(seed)
  gseKEGG(
    geneList      = ranking,
    organism      = KEGG_ORGANISM,
    keyType       = "ncbi-geneid",   # explicit; gseKEGG accepts NCBI Gene IDs
    minGSSize     = min_size,
    maxGSSize     = max_size,
    pvalueCutoff  = 1,
    pAdjustMethod = "BH",
    eps           = 0,
    nPermSimple   = n_perm,
    verbose       = FALSE
  )
}

run_reactome <- function(ranking, min_size, max_size, n_perm, seed) {
  set.seed(seed)
  gsePathway(
    geneList      = ranking,
    organism      = REACTOME_ORGANISM,
    minGSSize     = min_size,
    maxGSSize     = max_size,
    pvalueCutoff  = 1,
    pAdjustMethod = "BH",
    eps           = 0,
    nPermSimple   = n_perm,
    verbose       = FALSE
  )
}


# --- Saving & Plotting Helpers ---

save_dual <- function(plot, path_no_ext, w = 9, h = 7, dpi = 600, log_path = NULL) {
  # Save PDF (archival) and high-DPI PNG (Word embed) per project convention.
  # On systems where the X11/Cairo runtime libs are missing, cairo_pdf emits
  # a *warning* (not an error) and does not produce a file. The post-write
  # file.size check distinguishes both failure modes from real success.
  pdf_path <- paste0(path_no_ext, ".pdf")
  if (file.exists(pdf_path)) file.remove(pdf_path)

  cairo_ok <- tryCatch({
    suppressWarnings(
      ggsave(pdf_path, plot, width = w, height = h, device = cairo_pdf)
    )
    file.exists(pdf_path) && file.size(pdf_path) > 1000L
  }, error = function(e) FALSE)

  if (!cairo_ok) {
    if (!is.null(log_path)) {
      log_msg(log_path,
              sprintf("    NOTE: cairo_pdf unavailable; using default pdf device for %s",
                       basename(pdf_path)))
    }
    ggsave(pdf_path, plot, width = w, height = h)
  }

  ggsave(paste0(path_no_ext, ".png"), plot, width = w, height = h, dpi = dpi)
}

save_results_table <- function(gsea_obj, csv_path, fdr, source_label) {
  if (is.null(gsea_obj) || nrow(as.data.frame(gsea_obj)) == 0) return(NULL)
  df <- as.data.frame(gsea_obj) |>
    arrange(p.adjust) |>
    mutate(
      direction       = ifelse(NES > 0, "up_in_ketamine", "down_in_ketamine"),
      passes_fdr      = p.adjust <= fdr,
      convergence_key = normalize_convergence_key(ID, source_label)
    ) |>
    relocate(convergence_key, .after = ID)
  write_csv(df, csv_path)
  invisible(df)
}

make_dotplot <- function(gsea_obj, title, n_terms, fdr, log_path = NULL, db_label = "") {
  # S4 copy-on-modify: assigning to a local `sig` makes a copy; the original
  # gsea_obj@result is untouched. Filtering @result on the copy lets
  # enrichplot's fortify/color scales reflect only significant rows.
  res <- as.data.frame(gsea_obj)
  if (sum(res$p.adjust <= fdr, na.rm = TRUE) == 0) return(NULL)
  sig <- gsea_obj
  sig@result <- sig@result |> filter(p.adjust <= fdr) |> arrange(p.adjust)
  tryCatch(
    dotplot(sig, showCategory = n_terms, split = ".sign",
            color = "p.adjust") +
      facet_grid(. ~ .sign) +
      ggtitle(title) +
      theme(plot.title = element_text(face = "bold")),
    error = function(e) {
      if (!is.null(log_path))
        log_msg(log_path, sprintf("    dotplot failed (%s): %s",
                                   db_label, conditionMessage(e)))
      NULL
    }
  )
}

make_ridgeplot <- function(gsea_obj, title, n_terms, fdr, log_path = NULL, db_label = "") {
  res <- as.data.frame(gsea_obj)
  if (sum(res$p.adjust <= fdr, na.rm = TRUE) == 0) return(NULL)
  sig <- gsea_obj
  sig@result <- sig@result |> filter(p.adjust <= fdr) |> arrange(p.adjust)
  # enrichplot::ridgeplot.gseaResult checks `inherits(showCategory, "numeric")`,
  # which returns FALSE for integers (optparse delivers --top_n_plot as int).
  # The check is is.numeric-incompatible; cast to double to take the intended
  # integer-count branch.
  tryCatch(
    ridgeplot(sig, showCategory = as.double(n_terms)) +
      ggtitle(title) +
      theme(plot.title = element_text(face = "bold")),
    error = function(e) {
      if (!is.null(log_path))
        log_msg(log_path, sprintf("    ridgeplot failed (%s): %s",
                                   db_label, conditionMessage(e)))
      NULL
    }
  )
}

make_top_running_score <- function(gsea_obj, title, top_n, fdr, log_path = NULL, db_label = "") {
  sig <- as.data.frame(gsea_obj) |> filter(p.adjust <= fdr)
  if (nrow(sig) == 0) return(NULL)
  top_up <- sig |> filter(NES > 0) |> slice_min(p.adjust, n = top_n, with_ties = FALSE) |> pull(ID)
  top_dn <- sig |> filter(NES < 0) |> slice_min(p.adjust, n = top_n, with_ties = FALSE) |> pull(ID)
  ids <- c(top_up, top_dn)
  if (length(ids) == 0) return(NULL)
  tryCatch(
    gseaplot2(gsea_obj, geneSetID = ids, title = title, pvalue_table = TRUE),
    error = function(e) {
      if (!is.null(log_path))
        log_msg(log_path, sprintf("    running-score failed (%s): %s",
                                   db_label, conditionMessage(e)))
      NULL
    }
  )
}


# --- Per-Database Pipeline ---

run_db <- function(db_label, runner, ranking, opt, min_size, out_dir, log_path,
                   size_survey_fn = NULL) {
  log_msg(log_path, sprintf("  Running GSEA: %s (min_size=%d, max_size=%d)",
                            db_label, min_size, opt$max_size))

  # Set-size survey (independent of GSEA; counts what is excluded by cutoffs)
  size_summary <- if (!is.null(size_survey_fn)) {
    sizes <- tryCatch(size_survey_fn(unique(names(ranking))),
                      error = function(e) {
                        log_msg(log_path, sprintf("    Size-survey error (%s): %s",
                                                   db_label, conditionMessage(e)))
                        integer(0)
                      })
    s <- summarize_size_exclusions(sizes, min_size, opt$max_size)
    log_msg(log_path,
            sprintf("    %s set-size survey: total=%s, below_min=%s, above_max=%s, within=%s",
                    db_label,
                    format(s$n_potential %||% "NA"),
                    format(s$n_below_min %||% "NA"),
                    format(s$n_above_max %||% "NA"),
                    format(s$n_within    %||% "NA")))
    s
  } else NULL

  gsea <- tryCatch(
    runner(ranking, min_size, opt$max_size, opt$n_perm, opt$seed),
    error = function(e) {
      log_msg(log_path, sprintf("    ERROR (%s): %s", db_label, conditionMessage(e)))
      NULL
    }
  )

  if (is.null(gsea) || nrow(as.data.frame(gsea)) == 0) {
    log_msg(log_path, sprintf("    %s: no enrichment results", db_label))
    return(list(gsea = NULL, df = NULL, size_summary = size_summary))
  }

  # Convert ENTREZ IDs in leading-edge sets to symbols for readability
  gsea <- tryCatch(
    setReadable(gsea, OrgDb = ORGANISM_DB, keyType = "ENTREZID"),
    error = function(e) {
      log_msg(log_path,
              sprintf("    NOTE: setReadable failed for %s (%s); leaving IDs in source format",
                       db_label, conditionMessage(e)))
      gsea
    }
  )

  res_df <- save_results_table(
    gsea,
    file.path(out_dir, sprintf("%s_gsea_results.csv", db_label)),
    opt$fdr, db_label
  )
  n_sig <- sum(res_df$passes_fdr, na.rm = TRUE)
  log_msg(log_path,
          sprintf("    %s: %d gene sets tested, %d significant at FDR %.3f",
                  db_label, nrow(res_df), n_sig, opt$fdr))

  if (n_sig > 0) {
    p_dot <- make_dotplot(gsea, sprintf("GSEA: %s", db_label),
                          opt$top_n_plot, opt$fdr, log_path, db_label)
    if (!is.null(p_dot)) save_dual(p_dot,
                                   file.path(out_dir, sprintf("%s_dotplot", db_label)),
                                   w = 10, h = 8, dpi = opt$dpi, log_path = log_path)

    p_ridge <- make_ridgeplot(gsea, sprintf("GSEA: %s", db_label),
                              opt$top_n_plot, opt$fdr, log_path, db_label)
    if (!is.null(p_ridge)) save_dual(p_ridge,
                                     file.path(out_dir, sprintf("%s_ridgeplot", db_label)),
                                     w = 10, h = 10, dpi = opt$dpi, log_path = log_path)

    p_run <- make_top_running_score(gsea,
                                    sprintf("GSEA: %s top terms", db_label),
                                    top_n = opt$top_n_running, opt$fdr,
                                    log_path = log_path, db_label = db_label)
    if (!is.null(p_run)) save_dual(p_run,
                                   file.path(out_dir,
                                             sprintf("%s_top_terms_running_score", db_label)),
                                   w = 12, h = 9, dpi = opt$dpi, log_path = log_path)
  }

  list(gsea = gsea, df = res_df, size_summary = size_summary)
}


# --- Concordance Between Ranking Metrics ---

compute_concordance <- function(results, runner_names) {
  map_dfr(runner_names, function(db) {
    a <- results$log2fc[[db]]$df
    b <- results$signed_logp[[db]]$df
    if (is.null(a) || is.null(b)) {
      return(tibble(
        database          = db,
        n_sig_log2fc      = NA_integer_,
        n_sig_signed_logp = NA_integer_,
        n_overlap         = NA_integer_,
        jaccard           = NA_real_,
        nes_spearman      = NA_real_,
        comment           = "no_results_in_one_or_both_rankings"
      ))
    }
    a_sig <- a |> filter(passes_fdr) |> pull(ID)
    b_sig <- b |> filter(passes_fdr) |> pull(ID)
    overlap <- intersect(a_sig, b_sig)
    uni     <- union(a_sig, b_sig)

    if (length(a_sig) == 0 && length(b_sig) == 0) {
      jacc <- NA_real_; comment_j <- "no_significant_sets_in_either_ranking"
    } else if (length(a_sig) == 0 || length(b_sig) == 0) {
      jacc <- NA_real_; comment_j <- "no_significant_sets_in_one_ranking"
    } else {
      jacc <- length(overlap) / length(uni); comment_j <- ""
    }

    # dplyr::select is explicit because AnnotationDbi (loaded transitively
    # via reactome.db / ReactomePA) exports a `select` generic that shadows
    # dplyr::select.
    shared <- inner_join(
      a |> dplyr::select(ID, NES_a = NES),
      b |> dplyr::select(ID, NES_b = NES),
      by = "ID"
    )
    rho <- if (nrow(shared) > 5) {
      suppressWarnings(cor(shared$NES_a, shared$NES_b, method = "spearman"))
    } else NA_real_

    tibble(
      database          = db,
      n_sig_log2fc      = length(a_sig),
      n_sig_signed_logp = length(b_sig),
      n_overlap         = length(overlap),
      jaccard           = jacc,
      nes_spearman      = rho,
      comment           = comment_j
    )
  })
}


# --- Report Writer ---

write_report <- function(all_results, concordance, opt, prep,
                          kegg_version_path, path) {
  con <- file(path, "w")
  on.exit(close(con))
  writeLines(c(
    strrep("=", 70),
    "KETAMINE PROTEOMICS - GSEA ANALYSIS REPORT",
    strrep("=", 70),
    sprintf("Run:                 %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("Input file:          %s", opt$input),
    sprintf("Proteins ranked:     %d (post-category-filter, post-collapse)", prep$raw_n),
    sprintf("FDR threshold:       %.3f", opt$fdr),
    sprintf("min_size values:     %s (primary = %d)",
            paste(opt$min_sizes, collapse = ", "), opt$min_sizes[1]),
    sprintf("max_size:            %d", opt$max_size),
    sprintf("Permutations:        %d (fgsea adaptive; actual count varies per set)",
            opt$n_perm),
    sprintf("Random seed:         %d", opt$seed),
    sprintf("KEGG version:        captured in %s", basename(kegg_version_path)),
    ""
  ), con)

  for (ms_label in names(all_results)) {
    writeLines(c(
      strrep("=", 70),
      sprintf("MIN_SIZE RUN: %s", ms_label),
      strrep("=", 70)
    ), con)

    results <- all_results[[ms_label]]
    for (rname in names(results)) {
      writeLines(c(
        strrep("-", 70),
        sprintf("Ranking metric: %s", rname),
        strrep("-", 70)
      ), con)

      for (db in names(results[[rname]])) {
        r <- results[[rname]][[db]]

        if (!is.null(r$size_summary)) {
          s <- r$size_summary
          writeLines(sprintf(
            "  %-9s set-size survey: %s potential / %s below min / %s above max / %s within",
            db,
            format(s$n_potential %||% "NA"),
            format(s$n_below_min %||% "NA"),
            format(s$n_above_max %||% "NA"),
            format(s$n_within    %||% "NA")
          ), con)
        }

        if (is.null(r$df) || nrow(r$df) == 0) {
          writeLines(sprintf("  %-9s: no enrichment results", db), con)
          writeLines("", con)
          next
        }
        df    <- r$df
        n_sig <- sum(df$passes_fdr, na.rm = TRUE)

        writeLines(c(
          sprintf("  %s:", db),
          sprintf("    Gene sets tested:        %d", nrow(df)),
          sprintf("    Significant (FDR<=%.3f): %d", opt$fdr, n_sig)
        ), con)

        if (n_sig > 0) {
          sig    <- df |> filter(passes_fdr) |> arrange(p.adjust)
          top_up <- sig |> filter(NES > 0) |> slice_head(n = 10)
          top_dn <- sig |> filter(NES < 0) |> slice_head(n = 10)

          if (nrow(top_up) > 0) {
            writeLines("    Top up-enriched (positive NES; ranking high in ketamine):", con)
            for (i in seq_len(nrow(top_up))) {
              writeLines(sprintf("      [%2d] NES=%+.2f  FDR=%.2e  %s",
                                 i, top_up$NES[i], top_up$p.adjust[i],
                                 substr(top_up$Description[i], 1, 70)), con)
            }
          }
          if (nrow(top_dn) > 0) {
            writeLines("    Top down-enriched (negative NES; ranking low in ketamine):", con)
            for (i in seq_len(nrow(top_dn))) {
              writeLines(sprintf("      [%2d] NES=%+.2f  FDR=%.2e  %s",
                                 i, top_dn$NES[i], top_dn$p.adjust[i],
                                 substr(top_dn$Description[i], 1, 70)), con)
            }
          }
        }
        writeLines("", con)
      }
    }
  }

  writeLines(c(
    strrep("=", 70),
    sprintf("RANKING METRIC CONCORDANCE (log2FC vs signed -log10(p); min_size = %d)",
            opt$min_sizes[1]),
    strrep("=", 70),
    "Per-database overlap of significant gene sets and Spearman correlation of",
    "NES. Spearman is computed on the intersection of gene sets tested in both",
    "rankings, not the union; very low overlap renders the metric uninformative.",
    "Concordance is computed on the primary min_size run only.",
    ""
  ), con)
  for (i in seq_len(nrow(concordance))) {
    cmt <- if (nzchar(concordance$comment[i])) sprintf("  [%s]",
                                                       concordance$comment[i]) else ""
    writeLines(sprintf(
      "  %-9s  log2FC_sig=%3s  signedlogp_sig=%3s  overlap=%3s  Jaccard=%s  NES_Spearman=%s%s",
      concordance$database[i],
      format(concordance$n_sig_log2fc[i]      %||% "NA"),
      format(concordance$n_sig_signed_logp[i] %||% "NA"),
      format(concordance$n_overlap[i]         %||% "NA"),
      if (is.na(concordance$jaccard[i])) "NA"
        else sprintf("%.2f", concordance$jaccard[i]),
      if (is.na(concordance$nes_spearman[i])) "NA"
        else sprintf("%.2f", concordance$nes_spearman[i]),
      cmt
    ), con)
  }
  writeLines(c(
    "",
    "Interpretation hints (defer setting hard thresholds until first run):",
    "  - Jaccard quantifies overlap of which gene sets pass FDR. Very low",
    "    overlap suggests metric choice is load-bearing and both rankings",
    "    should appear in Results.",
    "  - NES_Spearman quantifies effect-size agreement over sets shared by both",
    "    rankings. Computed only on the intersection; not a global metric.",
    "  - When one ranking returns zero significant sets, Jaccard is reported",
    "    as NA rather than 0 to avoid misleading interpretation."
  ), con)
}


# --- Main ---

main <- function(opt) {
  log_path <- init_log(opt$log_dir)
  log_msg(log_path, "GSEA pipeline starting")
  log_msg(log_path, sprintf(
    "Parameters: FDR=%.3f, min_sizes=[%s], max_size=%d, n_perm=%d, seed=%d, DPI=%d",
    opt$fdr, paste(opt$min_sizes, collapse = ","),
    opt$max_size, opt$n_perm, opt$seed, opt$dpi
  ))

  # Step 1: load and prepare the ranking source. Both ranking metrics share
  # the same protein-level row selection by design (see load_ranked_list).
  prep     <- load_ranked_list(opt$input, log_path)
  rankings <- prep$rankings

  # Step 2: ensure output dir exists; save .rnk files (shared across min_size runs)
  if (!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)
  for (rname in names(rankings)) {
    vec <- sort(rankings[[rname]], decreasing = TRUE)
    write_tsv(
      tibble(gene = names(vec), score = vec),
      file.path(opt$output_dir, sprintf("ranked_list_%s.rnk", rname)),
      col_names = FALSE   # standard GSEA .rnk format is headerless
    )
  }
  log_msg(log_path, sprintf("Saved .rnk files to %s", opt$output_dir))

  # Step 3: capture KEGG release version (forensic record)
  capture_kegg_version(file.path(opt$output_dir, "kegg_version.txt"), log_path)

  # Step 4: map SYMBOL -> ENTREZ once per ranking; both rankings share the
  # same input gene universe so unmapped and alias-recovered symbols are
  # identical between them. Write the unmapped + recovered sidecars from
  # the log2fc pass only.
  unmapped_path  <- file.path(opt$output_dir, "unmapped_symbols.txt")
  recovered_path <- file.path(opt$output_dir, "alias_recovered_symbols.txt")
  log_msg(log_path, "Mapping gene symbols to ENTREZ IDs")
  rankings_entrez <- list(
    log2fc      = map_to_entrez(rankings$log2fc,      log_path, "log2fc",
                                 unmapped_out_path  = unmapped_path,
                                 recovered_out_path = recovered_path),
    signed_logp = map_to_entrez(rankings$signed_logp, log_path, "signed_logp",
                                 unmapped_out_path  = NULL,
                                 recovered_out_path = NULL)
  )

  runners <- list(go_bp = run_go_bp, kegg = run_kegg, reactome = run_reactome)
  size_surveys <- list(
    go_bp    = survey_go_bp_sizes,
    kegg     = survey_kegg_sizes,
    reactome = survey_reactome_sizes
  )

  # Step 5: outer loop over min_size values, inner loops over ranking metric
  # and database.
  all_results <- list()
  for (ms in opt$min_sizes) {
    ms_label <- sprintf("min%d", ms)
    log_msg(log_path, sprintf("=== min_size = %d (%s) ===", ms, ms_label))
    ms_dir <- file.path(opt$output_dir, ms_label)
    if (!dir.exists(ms_dir)) dir.create(ms_dir, recursive = TRUE)

    all_results[[ms_label]] <- list()
    for (rname in names(rankings_entrez)) {
      log_msg(log_path, sprintf("Ranking: %s", rname))
      out_dir <- file.path(ms_dir, rname)
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

      all_results[[ms_label]][[rname]] <- list()
      for (db in names(runners)) {
        all_results[[ms_label]][[rname]][[db]] <- run_db(
          db_label       = db,
          runner         = runners[[db]],
          ranking        = rankings_entrez[[rname]],
          opt            = opt,
          min_size       = ms,
          out_dir        = out_dir,
          log_path       = log_path,
          size_survey_fn = size_surveys[[db]]
        )
      }
    }
  }

  # Step 6: ranking-metric concordance computed on the primary min_size only.
  # Concordance is about metric robustness, not size-cutoff robustness;
  # running it on both min_sizes would dilute the message.
  primary_ms <- sprintf("min%d", opt$min_sizes[1])
  log_msg(log_path,
          sprintf("Computing ranking-metric concordance on primary min_size (%s)",
                   primary_ms))
  concordance <- compute_concordance(all_results[[primary_ms]], names(runners))
  write_csv(concordance, file.path(opt$output_dir, "ranking_metric_concordance.csv"))

  # Step 7: summary report
  log_msg(log_path, "Writing summary report")
  write_report(
    all_results       = all_results,
    concordance       = concordance,
    opt               = opt,
    prep              = prep,
    kegg_version_path = file.path(opt$output_dir, "kegg_version.txt"),
    path              = file.path(opt$output_dir, "gsea_analysis_report.txt")
  )

  # Step 8: capture session info for reproducibility
  capture.output(sessionInfo(),
                 file = file.path(opt$output_dir, "sessionInfo.txt"))

  log_msg(log_path, "GSEA pipeline complete")
}


# --- Execute ---
main(opt)

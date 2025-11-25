# title: 'create_differential_abundance_table.R'
# project: 'Astrocyte Proteomic Responses to Ketamine in Mice'
# author: 'Reina Hastings'
# contact: 'reinahastings13@gmail.com'
# date created: 11/14/2025
# last modified: 11/24/2025
# purpose:
#   Use vendor-computed statistics from Proteome Discoverer to build a
#   differential abundance table for ketamine vs control.
#   Classify proteins as upregulated, downregulated, ketamine-exclusive,
#   control-exclusive, or not significant based on:
#     - Abundance Ratio (ketamine/control)
#     - Adjusted p-value (FDR q-value)
#     - Presence/absence patterns (file membership)
# inputs:
#   - CSV with proteins detected in both groups
#   - CSV with proteins detected only in control
#   - CSV with proteins detected only in ketamine
# outputs:
#   - <outdir>/differential_abundance_table.tsv
#   - <outdir>/significant_proteins.tsv
#   - <outdir>/class_counts.tsv
#   - <logfile> (script progress and summary)
# required packages:
#   - tidyverse (readr, dplyr, tibble, stringr, purrr)
# usage:
#   Rscript create_differential_abundance_table.R \
#     -b ../data/proteins_both.csv \
#     -c ../data/proteins_control.csv \
#     -k ../data/proteins_ketamine.csv \
#     -o ../results/diff_abundance \
#     -l ../logs/create_diff_abundance_table.log

suppressPackageStartupMessages({
  library(tidyverse)
})

# --------------------------- helper functions --------------------------- #

log_msg <- function(...) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg <- paste0("[", ts, "] ", paste0(..., collapse = ""))
  cat(msg, "\n")
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  out <- list(
    both_file     = "../data/proteins_both.csv",
    control_file  = "../data/proteins_control.csv",
    ketamine_file = "../data/proteins_ketamine.csv",
    outdir        = "../results/diff_abundance",
    logfile       = NULL
  )

  if (length(args) == 0) return(out)

  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    val <- if (i + 1 <= length(args)) args[i + 1] else NA_character_

    if (key %in% c("-b", "--both")) {
      out$both_file <- val
    } else if (key %in% c("-c", "--control")) {
      out$control_file <- val
    } else if (key %in% c("-k", "--ketamine")) {
      out$ketamine_file <- val
    } else if (key %in% c("-o", "--outdir")) {
      out$outdir <- val
    } else if (key %in% c("-l", "--logfile")) {
      out$logfile <- val
    } else {
      warning("Unknown argument: ", key)
    }
    i <- i + 2
  }

  if (is.null(out$logfile)) {
    if (!dir.exists("logs")) dir.create("logs", recursive = TRUE)
    out$logfile <- file.path("logs", "create_diff_abundance_table.log")
  }

  out
}

# read CSV with logging
read_csv_logged <- function(path) {
  log_msg("Reading file: ", path)

  if (!file.exists(path)) {
    stop("Input file does not exist: ", path)
  }

  df <- readr::read_csv(path, show_col_types = FALSE)
  log_msg("  -> Read ", nrow(df), " rows and ", ncol(df), " columns.")
  df
}

# --------------------------- main script --------------------------- #

main <- function() {

  args <- parse_args()

  if (!dir.exists(args$outdir)) {
    dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
  }

  log_con <- file(args$logfile, open = "wt")
  sink(log_con, split = TRUE)
  on.exit({
    sink()
    close(log_con)
  }, add = TRUE)

  log_msg("------------------------------------------------------------")
  log_msg("Starting create_diff_abundance_table.R")
  log_msg("both_file     = ", args$both_file)
  log_msg("control_file  = ", args$control_file)
  log_msg("ketamine_file = ", args$ketamine_file)
  log_msg("outdir        = ", args$outdir)
  log_msg("logfile       = ", args$logfile)

  q_threshold     <- 0.05
  log2fc_moderate <- 0.58
  log2fc_strong   <- 1.0

  log_msg("Using thresholds:")
  log_msg("  q_threshold     = ", q_threshold)
  log_msg("  log2fc_moderate = ", log2fc_moderate)
  log_msg("  log2fc_strong   = ", log2fc_strong)

  # 1) read inputs ------------------------------------------------------- #
  both_df     <- read_csv_logged(args$both_file)
  control_df  <- read_csv_logged(args$control_file)
  ketamine_df <- read_csv_logged(args$ketamine_file)

  # 2) standardize colnames ---------------------------------------------- #
  standardize_cols <- function(df, presence_pattern) {

    df2 <- df %>%
      rename(
        ProteinID = any_of(c("Accession")),
        GeneSymbol = any_of(c("Gene Symbol", "GeneSymbol")),
        AbundanceRatio = any_of(c(
          "Abundance Ratio: (ketamine) / (control)",
          "Abundance Ratio (ketamine)/(control)"
        )),
        AbundanceRatioAdjP = any_of(c(
          "Abundance Ratio Adj. P-Value: (ketamine) / (control)",
          "Abundance Ratio Adj. P-Value (ketamine)/(control)"
        )),
        AbundanceRatioVarPct = any_of(c(
          "Abundance Ratio Variability [%]: (ketamine) / (control)",
          "Abundance Ratio Variability [%] (ketamine)/(control)"
        ))
      ) %>%
      mutate(presence_pattern = presence_pattern)

    df2
  }

  log_msg("Standardizing columns and tagging presence_pattern...")

  both_df_std     <- standardize_cols(both_df, "both_groups")
  control_df_std  <- standardize_cols(control_df, "control_only")
  ketamine_df_std <- standardize_cols(ketamine_df, "ketamine_only")

  log_msg("Rows after standardization:")
  log_msg("  both_groups   : ", nrow(both_df_std))
  log_msg("  control_only  : ", nrow(control_df_std))
  log_msg("  ketamine_only : ", nrow(ketamine_df_std))

  # 3) compute log2FC ---------------------------------------------------- #
  log_msg("Computing log2 fold-change for proteins detected in both groups...")

  both_df_std <- both_df_std %>%
    mutate(
      log2FC = if_else(
        !is.na(AbundanceRatio) & AbundanceRatio > 0,
        log2(AbundanceRatio),
        NA_real_
      )
    )

  # convert all columns to characters where necessary to avoid type conflicts
  harmonize_types <- function(df) {
    df %>% mutate(across(everything(), as.character))
  }

  master_df <- bind_rows(
    harmonize_types(both_df_std),
    harmonize_types(control_df_std),
    harmonize_types(ketamine_df_std)
  )

  log_msg("Total rows in master_df: ", nrow(master_df))

  # 4) classification ---------------------------------------------------- #
  log_msg("Applying classification rules...")

  master_df <- master_df %>%
    mutate(
      AbundanceRatio = as.numeric(AbundanceRatio),
      AbundanceRatioAdjP = as.numeric(AbundanceRatioAdjP),
      log2FC = as.numeric(log2FC),

      is_significant = !is.na(AbundanceRatioAdjP) &
                       AbundanceRatioAdjP <= q_threshold,

      change_class = case_when(
        presence_pattern == "ketamine_only" ~ "exclusive_ketamine",
        presence_pattern == "control_only"  ~ "exclusive_control",
        presence_pattern == "both_groups" &
          is_significant & !is.na(log2FC) & log2FC >= log2fc_moderate ~ "upregulated",
        presence_pattern == "both_groups" &
          is_significant & !is.na(log2FC) & log2FC <= -log2fc_moderate ~ "downregulated",
        TRUE ~ "not_significant"
      ),

      strong_effect = case_when(
        !is.na(log2FC) & abs(log2FC) >= log2fc_strong ~ TRUE,
        TRUE ~ FALSE
      )
    )

  class_counts <- master_df %>%
    count(change_class, name = "n_proteins") %>%
    arrange(desc(n_proteins))

  log_msg("Classification counts:")
  for (i in seq_len(nrow(class_counts))) {
    log_msg("  ", class_counts$change_class[i], ": ",
            class_counts$n_proteins[i])
  }

  sig_df <- master_df %>%
    filter(change_class %in% c(
      "upregulated", "downregulated",
      "exclusive_ketamine", "exclusive_control"
    ))

  log_msg("Total significant/biologically flagged proteins: ", nrow(sig_df))

  # 5) write outputs ---------------------------------------------------- #
  out_table_path  <- file.path(args$outdir, "diff_abundance_table.tsv")
  out_sig_path    <- file.path(args$outdir, "significant_proteins.tsv")
  out_counts_path <- file.path(args$outdir, "class_counts.tsv")

  log_msg("Writing outputs:")
  log_msg("  master table      -> ", out_table_path)
  log_msg("  significant table -> ", out_sig_path)
  log_msg("  class counts      -> ", out_counts_path)

  readr::write_tsv(master_df, out_table_path)
  readr::write_tsv(sig_df, out_sig_path)
  readr::write_tsv(class_counts, out_counts_path)

  log_msg("Done. create_diff_abundance_table.R finished successfully.")
  log_msg("------------------------------------------------------------")
}

# run main() only when script is executed directly, not when sourced
if (identical(environment(), globalenv()) && !length(sys.calls())) {
  main()
}
# title: 'plot_diff_abundance.R'
# project: 'Astrocyte Proteomic Responses to Ketamine in Mice'
# author: 'Reina Hastings'
# contact: 'reinahastings13@gmail.com'
# date created: 11/24/2025
# last modified: 11/24/2025
# purpose:
#   Create core visualizations from the differential abundance results:
#     - Volcano plot of log2 fold-change vs -log10(q-value)
#     - Bar plot of protein class counts (up, down, exclusive, not_significant)
#     - Optional heatmap of top |log2FC| proteins using group-level abundances
#   Inputs are the outputs of create_differential_abundance_table.R.
#
# required packages:
#   - tidyverse
#   - viridis
#
# inputs:
#   - differential abundance table:
#       ../results/diff_abundance/diff_abundance_table.tsv
#   - (optional) class counts table:
#       ../results/diff_abundance/class_counts.tsv
#
# outputs:
#   - <outdir>/volcano_plot.png
#   - <outdir>/class_counts_bar.png
#   - <outdir>/top_N_heatmap.png
#
# usage:
#   Rscript plot_diff_abundance.R \
#     -i ../results/diff_abundance/diff_abundance_table.tsv \
#     -o ../results/figs
# (single line version): Rscript plot_diff_abundance.R -i ../results/diff_abundance/diff_abundance_table.tsv -o ../results/figs

suppressPackageStartupMessages({
  library(tidyverse)
  library(viridis)
})

# --------------------------- helper functions --------------------------- #

log_msg <- function(...) {
  ts  <- format(Sys.time(), '%Y-%m-%d %H:%M:%S')
  msg <- paste0('[', ts, '] ', paste0(..., collapse = ''))
  cat(msg, '\n')
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  out <- list(
    da_file   = '../results/diff_abundance/differential_abundance_table.tsv',
    outdir    = '../results/figs',
    top_n     = 50
  )

  if (length(args) == 0) return(out)

  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    val <- if (i + 1 <= length(args)) args[i + 1] else NA_character_

    if (key %in% c('-i', '--input')) {
      out$da_file <- val
    } else if (key %in% c('-o', '--outdir')) {
      out$outdir <- val
    } else if (key %in% c('-n', '--top_N')) {
      out$top_n <- as.integer(val)
    } else {
      warning('Unknown argument: ', key)
    }
    i <- i + 2
  }

  out
}

read_tsv_logged <- function(path) {
  log_msg('Reading file: ', path)
  if (!file.exists(path)) {
    stop('Input file does not exist: ', path)
  }
  df <- readr::read_tsv(path, show_col_types = FALSE)
  log_msg('  -> Read ', nrow(df), ' rows and ', ncol(df), ' columns.')
  df
}

# --------------------------- plotting functions ------------------------- #

make_volcano_plot <- function(df, out_path,
                              q_threshold = 0.05,
                              log2fc_cut = 0.58) {

  log_msg('Preparing data for volcano plot...')

  # restrict to proteins with log2FC (both_groups)
  vol_df <- df %>%
    filter(!is.na(log2FC)) %>%
    mutate(
      neg_log10_q = case_when(
        !is.na(AbundanceRatioAdjP) & AbundanceRatioAdjP > 0 ~ -log10(AbundanceRatioAdjP),
        TRUE ~ NA_real_
      ),
      change_class = factor(
        change_class,
        levels = c('upregulated',
                   'downregulated',
                   'exclusive_ketamine',
                   'exclusive_control',
                   'not_significant')
      )
    )

  log_msg('Creating volcano plot...')
  p <- ggplot(vol_df, aes(x = log2FC, y = neg_log10_q, colour = change_class)) +
    geom_jitter(width = 0.05, height = 0.05, alpha = 0.8, size = 1.8, na.rm = TRUE) +
    geom_vline(xintercept = c(-log2fc_cut, log2fc_cut),
               linetype = 'dashed', linewidth = 0.4) +
    geom_hline(yintercept = -log10(q_threshold),
               linetype = 'dashed', linewidth = 0.4) +
    scale_colour_manual(
      values = c(
        'upregulated'        = '#E64B35',
        'downregulated'      = '#4DBBD5',
        'exclusive_ketamine' = '#F39B7F',
        'exclusive_control'  = '#91D1C2',
        'not_significant'    = 'grey70'
      ),
      na.translate = FALSE
    ) +
    labs(
      title    = 'Volcano plot: Ketamine vs Control',
      x        = 'log2 Fold-change (ketamine / control)',
      y        = expression(-log[10]('Adjusted p-value')),
      colour   = 'Class'
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title   = element_text(hjust = 0.5, face = 'bold'),
      legend.title = element_text(face = 'bold')
    )

  ggsave(out_path, p, width = 7, height = 5, dpi = 300)
  log_msg('  -> Volcano plot saved to: ', out_path)
}

make_class_counts_bar <- function(df, out_path) {
  log_msg('Creating class counts bar plot...')

  counts <- df %>%
    count(change_class, name = 'n_proteins') %>%
    mutate(
      change_class = factor(
        change_class,
        levels = c('upregulated',
                   'downregulated',
                   'exclusive_ketamine',
                   'exclusive_control',
                   'not_significant')
      )
    ) %>%
    arrange(change_class)

  p <- ggplot(counts, aes(x = change_class, y = n_proteins, fill = change_class)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = n_proteins),
              vjust = -0.3, size = 3) +
    scale_fill_manual(
      values = c(
        'upregulated'        = '#E64B35',
        'downregulated'      = '#4DBBD5',
        'exclusive_ketamine' = '#F39B7F',
        'exclusive_control'  = '#91D1C2',
        'not_significant'    = 'grey70'
      ),
      guide = 'none'
    ) +
    labs(
      title = 'Protein class counts',
      x     = 'Class',
      y     = 'Number of proteins'
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = 'bold'),
      axis.text.x = element_text(angle = 30, hjust = 1)
    )

  ggsave(out_path, p, width = 6, height = 4, dpi = 300)
  log_msg('  -> Class counts bar plot saved to: ', out_path)
}

make_top_N_heatmap <- function(df, out_path, top_n = 50) {
  log_msg('Preparing data for top-N heatmap...')

  # rename group abundance columns (vendor names) if present
  df2 <- df %>%
    rename(
      Abund_control = any_of('Abundances (Grouped): control'),
      Abund_ketamine = any_of('Abundances (Grouped): ketamine')
    )

  if (!all(c('Abund_control', 'Abund_ketamine') %in% colnames(df2))) {
    log_msg('WARNING: Grouped abundance columns not found; skipping heatmap.')
    return(invisible(NULL))
  }

  # numeric conversion for means
  df2 <- df2 %>%
    mutate(
      Abund_control  = as.numeric(Abund_control),
      Abund_ketamine = as.numeric(Abund_ketamine)
    )

  # choose a label column for y-axis: GeneSymbol if present, else ProteinID
  label_col <- if ('GeneSymbol' %in% colnames(df2)) 'GeneSymbol' else 'ProteinID'

  # select top N by |log2FC| among proteins with defined log2FC
  top_df <- df2 %>%
    filter(!is.na(log2FC)) %>%
    arrange(desc(abs(log2FC))) %>%
    slice_head(n = top_n)

  if (nrow(top_df) == 0) {
    log_msg('WARNING: No proteins with log2FC available for heatmap; skipping.')
    return(invisible(NULL))
  }

  # wide -> long: two columns (control vs ketamine) per protein
  long_df <- top_df %>%
    select(any_of(c(label_col, 'log2FC', 'Abund_control', 'Abund_ketamine'))) %>%
    pivot_longer(
      cols = c('Abund_control', 'Abund_ketamine'),
      names_to = 'group',
      values_to = 'abundance'
    ) %>%
    mutate(
      group = recode(group,
                     'Abund_control'  = 'Control',
                     'Abund_ketamine' = 'Ketamine')
    )

  # z-score per protein to highlight patterns rather than absolute scale
  long_df <- long_df %>%
    group_by(.data[[label_col]]) %>%
    mutate(
      abundance_z = if (all(is.na(abundance))) NA_real_
                    else as.numeric(scale(abundance))
    ) %>%
    ungroup()

  # to keep ordering by effect size
  long_df <- long_df %>%
    mutate(
      !!label_col := factor(
        .data[[label_col]],
        levels = unique(top_df[[label_col]][order(top_df$log2FC, decreasing = TRUE)])
      )
    )

  log_msg('Creating top-N heatmap (N = ', top_n, ')...')

  p <- ggplot(long_df, aes(x = group, y = .data[[label_col]], fill = abundance_z)) +
    geom_tile(color = 'grey20', linewidth = 0.1) +
    scale_fill_viridis(option = 'magma', na.value = 'grey90') +
    labs(
      title = paste0('Top ', top_n, ' proteins by |log2FC|'),
      x     = 'Group',
      y     = label_col,
      fill  = 'Z-score\nwithin protein'
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = 'bold'),
      axis.text.y = element_text(size = 6)
    )

  ggsave(out_path, p, width = 5, height = 0.18 * nrow(top_df) + 1.5, dpi = 300)
  log_msg('  -> Heatmap saved to: ', out_path)
}

# --------------------------- main script --------------------------- #

main <- function() {
  args <- parse_args()

  if (!dir.exists(args$outdir)) {
    dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
  }

  log_msg('------------------------------------------------------------')
  log_msg('Starting plot_differential_abundance.R')
  log_msg('da_file = ', args$da_file)
  log_msg('outdir  = ', args$outdir)
  log_msg('top_n   = ', args$top_n)

  da_df <- read_tsv_logged(args$da_file)

  # ensure numeric types for key columns
  da_df <- da_df %>%
    mutate(
      log2FC             = as.numeric(log2FC),
      AbundanceRatioAdjP = suppressWarnings(as.numeric(AbundanceRatioAdjP))
    )

  # volcano plot
  volcano_path <- file.path(args$outdir, 'volcano_plot.png')
  make_volcano_plot(da_df, volcano_path)

  # class counts bar plot
  class_bar_path <- file.path(args$outdir, 'class_counts_bar.png')
  make_class_counts_bar(da_df, class_bar_path)

  # top-N heatmap (optional; will quietly skip if required cols missing)
  heatmap_path <- file.path(args$outdir, 'top_N_heatmap.png')
  make_top_N_heatmap(da_df, heatmap_path, top_n = args$top_n)

  log_msg('All plots generated.')
  log_msg('------------------------------------------------------------')
}

if (identical(environment(), globalenv()) && !length(sys.calls())) {
  main()
}
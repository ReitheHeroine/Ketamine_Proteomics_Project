# title: pathway_analysis.py
# project: Ketamine Proteomics Analysis Project
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2025-12-09
# last modified: 2026-06-01
#
# purpose:
#   Perform pathway enrichment analysis on differential abundance results using
#   over-representation analysis (ORA). Queries GO Biological Process, KEGG, and
#   Reactome (REAC) databases via gProfiler. Analyzes upregulated and
#   downregulated protein sets separately, including presence/absence proteins.
#
#   Databases are selectable with --sources (default: GO:BP KEGG; pass REAC to
#   add Reactome). Each source also emits a <set>_<token>_full.csv carrying
#   effective_domain_size and intersection gene symbols, which the GSEA-style
#   ORA figure scripts (plot_ora_lollipops.py, plot_ora_publication_figure.py)
#   consume offline. --enrich-only / --full-only run targeted, non-clobbering
#   source additions (see usage).
#
#   Also supports REVIGO filtering mode to reduce GO term redundancy.
#
# inputs:
#   Normal mode:
#     - all_proteins_categorized.csv: Master file from diff_abundance_analysis.py
#   
#   REVIGO mode:
#     - REVIGO output table (TSV from http://revigo.irb.hr/)
#     - Original pathway_analysis results directory
#
# outputs:
#   results/pathway_analysis/
#   ├── upregulated/
#   │   ├── upregulated_GO_BP.csv
#   │   ├── upregulated_KEGG.csv
#   │   ├── upregulated_GO_BP_plot.pdf/html
#   │   ├── upregulated_GO_BP_dotplot.pdf/html
#   │   ├── upregulated_GO_BP_table.pdf/html
#   │   ├── upregulated_GO_BP_network.pdf/html
#   │   ├── upregulated_GO_BP_clustered.pdf/html
#   │   ├── upregulated_GO_BP_publication_figure.pdf/html
#   │   ├── upregulated_KEGG_plot.pdf/html
#   │   ├── upregulated_KEGG_dotplot.pdf/html
#   │   ├── upregulated_KEGG_table.pdf/html
#   │   ├── upregulated_KEGG_network.pdf/html
#   │   ├── upregulated_KEGG_publication_figure.pdf/html
#   │   └── upregulated_summary.txt
#   ├── downregulated/
#   │   └── (similar structure, if significant terms found)
#   ├── revigo/  (if REVIGO mode used)
#   │   ├── upregulated/
#   │   │   ├── upregulated_GO_BP_revigo.csv
#   │   │   └── (all visualizations with _revigo suffix)
#   │   └── revigo_filtering_summary.txt
#   └── pathway_analysis_report.txt
#
# usage:
#   # Normal mode - run full pathway analysis
#   python pathway_analysis.py \
#       --input ../results/combined/all_proteins_categorized.csv \
#       --output_dir ../results/pathway_analysis \
#       --log_dir ../logs \
#       --pval_threshold 0.05 \
#       --fdr_threshold 0.05
#
#   copy/paste: python pathway_analysis.py --input ../results/combined/all_proteins_categorized.csv --output_dir ../results/pathway_analysis --log_dir ../logs --pval_threshold 0.05 --fdr_threshold 0.05
#
#   # Add Reactome only (clean + _full CSVs; no figures/summary/report)
#   python pathway_analysis.py --input data/all_proteins_categorized.csv \
#       --output_dir results/pathway_analysis --log_dir logs \
#       --sources REAC --enrich-only
#
#   # Backfill _full CSVs for GO:BP/KEGG without touching their clean CSVs
#   python pathway_analysis.py --input data/all_proteins_categorized.csv \
#       --output_dir results/pathway_analysis --log_dir logs \
#       --sources "GO:BP" KEGG --full-only
#
#   # REVIGO mode - filter GO terms using REVIGO output
#   python pathway_analysis.py \
#       --revigo ../results/pathway_analysis/revigo/Revigo_BP_Table.tsv \
#       --original_results ../results/pathway_analysis \
#       --output_dir ../results/pathway_analysis \
#       --log_dir ../logs \
#       --input ../results/combined/all_proteins_categorized.csv
#
#   copy/paste: python pathway_analysis.py --revigo ../results/pathway_analysis/revigo/Revigo_BP_Table.tsv --original_results ../results/pathway_analysis --output_dir ../results/pathway_analysis --log_dir ../logs --input ../results/combined/all_proteins_categorized.csv
#
# dependencies:
#   pip install gprofiler-official pandas numpy plotly kaleido networkx scipy
#
# notes:
#   - Normal mode requires internet connection to query gProfiler API
#   - Background set: all 903 detected proteins
#   - Upregulated set: significant up (p ≤ 0.05) + ketamine-specific
#   - Downregulated set: significant down (p ≤ 0.05) + control-specific
#   - REVIGO (Supek et al., 2011) reduces GO term redundancy by clustering
#     semantically similar terms

import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import argparse
import os
from datetime import datetime
from gprofiler import GProfiler
import networkx as nx
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist

# =============================================================================
# CONFIGURATION AND CONSTANTS
# =============================================================================

# --- g:Profiler version pin ---
# Pinned to the archived release that was live on 2026-01-06, when this
# analysis was first run for the thesis. Without pinning, g:Profiler serves
# whichever GO + Ensembl bundle is currently live; that bundle was refreshed
# on 2026-01-23 (GO release classes/2026-01-23) and the term sizes, gene
# annotations, and FDR rankings drifted noticeably between runs. Pinning
# the base_url to gprofiler_archive3/e113_eg59_p19 (GO classes/2025-03-16,
# Ensembl 113, Ensembl Genomes 59) reproduces the 2026-01-06 enrichment
# table exactly (verified against documented top-6 terms in
# results/tables/GO_BP_enrichment_lollipop_for_llm.md).
# Archive index: https://biit.cs.ut.ee/gprofiler/api/util/list_archives
GPROFILER_ARCHIVE_URL = 'https://biit.cs.ut.ee/gprofiler_archive3/e113_eg59_p19/'

# Column names from input data
GENE_COL = 'Gene Symbol'
CATEGORY_COL = 'category'
DIRECTION_COL = 'direction'
PVAL_COL = 'Abundance Ratio Adj. P-Value: (ketamine) / (control)'

# gProfiler settings
ORGANISM = 'mmusculus'  # Mus musculus (mouse)
SOURCES = ['GO:BP', 'KEGG']  # Default sources; overridable via --sources
#   GO:BP = GO Biological Process, KEGG = KEGG pathways, REAC = Reactome.

# --- Source naming helpers ---------------------------------------------------
# Map a gProfiler source token to (a) a human-readable display name used in
# figure titles/reports and (b) a filesystem-safe token used in output
# filenames. Centralizing these avoids the scattered
# "'GO Biological Process' if source == 'GO:BP' else 'KEGG Pathways'" ternaries
# and lets Reactome (REAC) be labelled/filed correctly. REAC has no colon, so
# the legacy source.replace(':', '_') idiom would have produced 'REAC'; the
# token map yields 'Reactome' to mirror the GSEA naming scheme.
SOURCE_DISPLAY_NAMES = {
    'GO:BP': 'GO Biological Process',
    'KEGG': 'KEGG Pathways',
    'REAC': 'Reactome Pathways',
}
SOURCE_FILE_TOKENS = {
    'GO:BP': 'GO_BP',
    'KEGG': 'KEGG',
    'REAC': 'Reactome',
}


def source_display_name(source):
    '''Human-readable database name for figure titles and reports.'''
    return SOURCE_DISPLAY_NAMES.get(source, source)


def source_file_token(source):
    '''Filesystem-safe database token for output filenames (REAC -> Reactome).'''
    return SOURCE_FILE_TOKENS.get(source, source.replace(':', '_'))


# Visual constants
COLOR_UP = '#D62728'       # Red for upregulated
COLOR_DOWN = '#1F77B4'     # Blue for downregulated

# =============================================================================
# LOGGING FUNCTIONS
# =============================================================================

def setup_logging(log_dir):
    '''Initialize timestamped log file.'''
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_filename = f'pathway_analysis_{timestamp}.log'
    log_path = os.path.join(log_dir, log_filename)
    
    with open(log_path, 'w') as f:
        f.write('=' * 70 + '\n')
        f.write('KETAMINE PROTEOMICS ANALYSIS PROJECT - PATHWAY ANALYSIS LOG\n')
        f.write('=' * 70 + '\n\n')
        f.write(f'Log file: {log_filename}\n')
        f.write(f'Started: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n\n')
    
    return log_path


def log_message(log_path, message, print_to_console=True):
    '''Write message to log file and optionally to console.'''
    timestamp = datetime.now().strftime('%H:%M:%S')
    formatted_msg = f'[{timestamp}] {message}'
    
    with open(log_path, 'a') as f:
        f.write(formatted_msg + '\n')
    
    if print_to_console:
        print(formatted_msg)


# =============================================================================
# DATA PREPARATION FUNCTIONS
# =============================================================================

def load_and_prepare_gene_sets(input_path, pval_threshold, log_path):
    '''
    Load data and prepare gene sets for pathway analysis.
    
    Returns:
        dict: Contains 'upregulated', 'downregulated', and 'background' gene lists
    '''
    log_message(log_path, f'Loading data from {input_path}')
    
    df = pd.read_csv(input_path)
    log_message(log_path, f'  Loaded {len(df)} proteins')
    
    # ---------------------------------------------------------------------
    # Extract gene symbols (remove any missing values)
    # ---------------------------------------------------------------------
    df = df[df[GENE_COL].notna() & (df[GENE_COL] != '')].copy()
    log_message(log_path, f'  Proteins with valid gene symbols: {len(df)}')
    
    # ---------------------------------------------------------------------
    # Define background set (all detected proteins)
    # ---------------------------------------------------------------------
    background = df[GENE_COL].unique().tolist()
    log_message(log_path, f'  Background set size: {len(background)}')
    
    # ---------------------------------------------------------------------
    # Define upregulated set:
    #   - Significant upregulated (quantitative, p <= threshold, up direction)
    #   - Ketamine-specific (presence/absence)
    # ---------------------------------------------------------------------
    sig_up = df[
        (df[CATEGORY_COL] == 'quantitative') &
        (df[PVAL_COL] <= pval_threshold) &
        (df[DIRECTION_COL] == 'up_in_ketamine')
    ][GENE_COL].tolist()
    
    ket_specific = df[
        df[CATEGORY_COL] == 'presence_absence_ketamine_specific'
    ][GENE_COL].tolist()
    
    upregulated = list(set(sig_up + ket_specific))
    
    log_message(log_path, f'  Upregulated set:')
    log_message(log_path, f'    Significant up (p <= {pval_threshold}): {len(sig_up)}')
    log_message(log_path, f'    Ketamine-specific: {len(ket_specific)}')
    log_message(log_path, f'    Total (unique): {len(upregulated)}')
    
    # ---------------------------------------------------------------------
    # Define downregulated set:
    #   - Significant downregulated (quantitative, p <= threshold, down direction)
    #   - Control-specific (presence/absence)
    # ---------------------------------------------------------------------
    sig_down = df[
        (df[CATEGORY_COL] == 'quantitative') &
        (df[PVAL_COL] <= pval_threshold) &
        (df[DIRECTION_COL] == 'down_in_ketamine')
    ][GENE_COL].tolist()
    
    ctrl_specific = df[
        df[CATEGORY_COL] == 'presence_absence_control_specific'
    ][GENE_COL].tolist()
    
    downregulated = list(set(sig_down + ctrl_specific))
    
    log_message(log_path, f'  Downregulated set:')
    log_message(log_path, f'    Significant down (p <= {pval_threshold}): {len(sig_down)}')
    log_message(log_path, f'    Control-specific: {len(ctrl_specific)}')
    log_message(log_path, f'    Total (unique): {len(downregulated)}')
    
    return {
        'upregulated': upregulated,
        'downregulated': downregulated,
        'background': background,
        'protein_data': df,  # Return full protein data for fold change info
        'counts': {
            'sig_up': len(sig_up),
            'ket_specific': len(ket_specific),
            'sig_down': len(sig_down),
            'ctrl_specific': len(ctrl_specific)
        }
    }


# =============================================================================
# PATHWAY ANALYSIS FUNCTIONS
# =============================================================================

def run_gprofiler_enrichment(gene_list, background, source, log_path):
    '''
    Run gProfiler enrichment analysis for a single source (GO:BP or KEGG).
    
    Parameters:
        gene_list (list): List of gene symbols to analyze
        background (list): Background gene list
        source (str): Database source ('GO:BP' or 'KEGG')
        log_path (str): Path to log file
    
    Returns:
        pd.DataFrame: Enrichment results
    '''
    log_message(log_path, f'    Querying {source}...')

    gp = GProfiler(return_dataframe=True, base_url=GPROFILER_ARCHIVE_URL)
    
    try:
        results = gp.profile(
            organism=ORGANISM,
            query=gene_list,
            background=background,
            sources=[source],
            user_threshold=0.05,  # Initial threshold (we'll filter by FDR later)
            significance_threshold_method='fdr',  # Benjamini-Hochberg
            no_evidences=False  # Include evidence codes
        )
        
        if results is not None and len(results) > 0:
            log_message(log_path, f'      Found {len(results)} enriched terms')
            return results
        else:
            log_message(log_path, f'      No enriched terms found')
            return pd.DataFrame()
            
    except Exception as e:
        log_message(log_path, f'      ERROR: {str(e)}')
        return pd.DataFrame()


def run_pathway_analysis(gene_list, background, set_name, output_dir, log_path,
                         fdr_threshold, write_clean=True):
    '''
    Run full pathway analysis for a gene set across the configured SOURCES.

    For every source two artifacts can be written:
      - the "clean" CSV (<set>_<token>.csv): the curated, renamed columns used
        by the legacy in-script figures and the thesis tables.
      - the "_full" CSV (<set>_<token>_full.csv): always written. Carries the
        extra fields the GSEA-style ORA figure scripts need offline, namely
        effective_domain_size (for fold enrichment) and the semicolon-joined
        intersection gene symbols (the ORA analog of GSEA core_enrichment).

    Parameters:
        gene_list (list): List of gene symbols
        background (list): Background gene list
        set_name (str): Name of the set ('upregulated' or 'downregulated')
        output_dir (str): Directory for outputs
        log_path (str): Path to log file
        fdr_threshold (float): FDR threshold for significance
        write_clean (bool): If False, skip writing the clean CSV and only emit
            the _full CSV. Used to backfill _full files for GO:BP/KEGG without
            overwriting their already-committed clean CSVs.

    Returns:
        dict: Results for each source
    '''
    log_message(log_path, f'  Running pathway analysis for {set_name} set ({len(gene_list)} genes)...')
    
    # Create output subdirectory
    set_output_dir = os.path.join(output_dir, set_name)
    os.makedirs(set_output_dir, exist_ok=True)
    
    results = {}
    
    for source in SOURCES:
        # ---------------------------------------------------------------------
        # Run enrichment
        # ---------------------------------------------------------------------
        df = run_gprofiler_enrichment(gene_list, background, source, log_path)
        
        if len(df) == 0:
            results[source] = {'clean': pd.DataFrame(), 'full': pd.DataFrame()}
            continue
        
        # Store full results (including intersections for network/clustering)
        full_results = df.copy()
        
        # ---------------------------------------------------------------------
        # Filter by FDR threshold
        # ---------------------------------------------------------------------
        df_sig = df[df['p_value'] <= fdr_threshold].copy()
        log_message(log_path, f'      Significant at FDR <= {fdr_threshold}: {len(df_sig)}')
        
        # ---------------------------------------------------------------------
        # Clean up and format results
        # ---------------------------------------------------------------------
        if len(df_sig) > 0:
            # Select and rename relevant columns
            cols_to_keep = [
                'native', 'name', 'p_value', 'term_size', 
                'query_size', 'intersection_size', 'precision', 'recall'
            ]
            cols_available = [c for c in cols_to_keep if c in df_sig.columns]
            df_clean = df_sig[cols_available].copy()
            
            # Rename columns for clarity
            df_clean = df_clean.rename(columns={
                'native': 'term_id',
                'name': 'term_name',
                'p_value': 'fdr_pvalue',
                'intersection_size': 'gene_count',
                'precision': 'gene_ratio'
            })
            
            # Sort by p-value
            df_clean = df_clean.sort_values('fdr_pvalue')
            
            # Add gene ratio as fraction string
            if 'gene_count' in df_clean.columns and 'query_size' in df_clean.columns:
                df_clean['ratio'] = df_clean.apply(
                    lambda r: f"{r['gene_count']}/{r['query_size']}", axis=1
                )
            
            results[source] = {'clean': df_clean, 'full': full_results}

            source_name = source_file_token(source)

            # --- Clean CSV (legacy figures + thesis tables) -------------------
            if write_clean:
                csv_path = os.path.join(set_output_dir,
                                        f'{set_name}_{source_name}.csv')
                df_clean.to_csv(csv_path, index=False)
                log_message(log_path, f'      Saved: {csv_path}')

            # --- Full CSV (offline inputs for GSEA-style ORA figures) ---------
            # Always written. Adds effective_domain_size (denominator for fold
            # enrichment) and the intersection gene symbols, neither of which
            # survive into the clean CSV.
            full_csv_path = os.path.join(set_output_dir,
                                         f'{set_name}_{source_name}_full.csv')
            write_full_results_csv(df_sig, full_csv_path, log_path)
        else:
            results[source] = {'clean': pd.DataFrame(), 'full': full_results}

    return results


def write_full_results_csv(df_sig, out_path, log_path):
    '''
    Persist the significant gProfiler rows with the fields the GSEA-style ORA
    figure scripts consume offline.

    Columns (mirrors the clean CSV where they overlap, plus the two extras):
        term_id, term_name, fdr_pvalue, term_size, query_size,
        intersection_size, effective_domain_size, precision, recall,
        intersection_genes (semicolon-joined gene symbols)

    The gProfiler 'intersections' column is a Python list of gene symbols when
    the query was run with no_evidences=False; it is serialized here as a
    ';'-joined string so the CSV round-trips cleanly.
    '''
    def _join_intersections(value):
        # gProfiler returns a list of gene symbols; guard against the rare
        # NaN / string fallbacks seen elsewhere in this module.
        if isinstance(value, list):
            return ';'.join(str(g) for g in value)
        if isinstance(value, str):
            return ';'.join(g.strip() for g in value.split(',') if g.strip())
        return ''

    cols = ['native', 'name', 'p_value', 'term_size', 'query_size',
            'intersection_size', 'effective_domain_size', 'precision',
            'recall']
    available = [c for c in cols if c in df_sig.columns]
    df_full = df_sig[available].copy()
    df_full = df_full.rename(columns={
        'native': 'term_id',
        'name': 'term_name',
        'p_value': 'fdr_pvalue',
    })
    if 'intersections' in df_sig.columns:
        df_full['intersection_genes'] = df_sig['intersections'].apply(
            _join_intersections
        ).values
    else:
        df_full['intersection_genes'] = ''
    df_full = df_full.sort_values('fdr_pvalue')
    df_full.to_csv(out_path, index=False)
    log_message(log_path, f'      Saved: {out_path}')


# =============================================================================
# VISUALIZATION FUNCTIONS
# =============================================================================

def create_enrichment_bar_plot(df, source, set_name, output_dir, log_path, n_terms=15, filename_suffix=''):
    '''
    Create horizontal bar plot of top enriched terms.
    
    Parameters:
        df (pd.DataFrame): Enrichment results
        source (str): Database source name
        set_name (str): 'upregulated' or 'downregulated'
        output_dir (str): Output directory
        log_path (str): Path to log file
        n_terms (int): Number of top terms to display
        filename_suffix (str): Optional suffix for output filenames (e.g., '_revigo')
    '''
    if len(df) == 0:
        log_message(log_path, f'    Skipping plot for {source} (no significant terms)')
        return
    
    # Select top terms
    plot_df = df.head(n_terms).copy()
    
    # Reverse order for horizontal bar plot (top term at top)
    plot_df = plot_df.iloc[::-1]
    
    # Set color based on direction
    color = COLOR_UP if set_name == 'upregulated' else COLOR_DOWN
    
    # Create figure
    fig = go.Figure()
    
    fig.add_trace(go.Bar(
        x=-np.log10(plot_df['fdr_pvalue']),
        y=plot_df['term_name'],
        orientation='h',
        marker_color=color,
        text=[f"{row['gene_count']} genes" for _, row in plot_df.iterrows()],
        textposition='outside',
        textfont=dict(size=9),
        hovertemplate='<b>%{y}</b><br>' +
                      '-log₁₀(FDR): %{x:.2f}<br>' +
                      '<extra></extra>'
    ))
    
    # Format source name for title
    source_display = source_display_name(source)
    direction_display = 'Upregulated' if set_name == 'upregulated' else 'Downregulated'
    revigo_label = ' (REVIGO filtered)' if filename_suffix else ''
    
    # Layout
    fig.update_layout(
        title=dict(
            text=f'{source_display}: {direction_display} Proteins{revigo_label}',
            font=dict(size=16, family='Arial Black')
        ),
        xaxis=dict(
            title=dict(text='-log₁₀(FDR p-value)', font=dict(size=12)),
            gridcolor='rgba(0,0,0,0.1)'
        ),
        yaxis=dict(
            title=dict(text=''),
            tickfont=dict(size=10)
        ),
        plot_bgcolor='white',
        width=800,
        height=max(400, 50 + len(plot_df) * 25),
        margin=dict(l=300, r=80, t=60, b=60)
    )
    
    # Save outputs
    source_name = source_file_token(source)
    set_output_dir = os.path.join(output_dir, set_name)
    
    pdf_path = os.path.join(set_output_dir, f'{set_name}_{source_name}_plot{filename_suffix}.pdf')
    html_path = os.path.join(set_output_dir, f'{set_name}_{source_name}_plot{filename_suffix}.html')
    
    fig.write_image(pdf_path, scale=2)
    fig.write_html(html_path)
    
    log_message(log_path, f'    Saved: {pdf_path}')
    log_message(log_path, f'    Saved: {html_path}')


def create_dot_plot(df, source, set_name, output_dir, log_path, n_terms=15, filename_suffix=''):
    '''
    Create dot plot showing gene ratio, gene count, and p-value.
    
    Parameters:
        df (pd.DataFrame): Enrichment results
        source (str): Database source name
        set_name (str): 'upregulated' or 'downregulated'
        output_dir (str): Output directory
        log_path (str): Path to log file
        n_terms (int): Number of top terms to display
        filename_suffix (str): Optional suffix for output filenames (e.g., '_revigo')
    '''
    if len(df) == 0 or 'gene_ratio' not in df.columns:
        return
    
    # Select top terms
    plot_df = df.head(n_terms).copy()
    
    # Reverse order for plot (top term at top)
    plot_df = plot_df.iloc[::-1]
    
    # Create figure
    fig = go.Figure()
    
    fig.add_trace(go.Scatter(
        x=plot_df['gene_ratio'],
        y=plot_df['term_name'],
        mode='markers',
        marker=dict(
            size=plot_df['gene_count'] * 2 + 5,  # Scale by gene count
            color=-np.log10(plot_df['fdr_pvalue']),
            colorscale='Reds' if set_name == 'upregulated' else 'Blues',
            showscale=True,
            colorbar=dict(title='-log₁₀(FDR)')
        ),
        text=[f"<b>{row['term_name']}</b><br>" +
              f"Genes: {row['gene_count']}<br>" +
              f"Gene ratio: {row['gene_ratio']:.3f}<br>" +
              f"FDR: {row['fdr_pvalue']:.2e}" 
              for _, row in plot_df.iterrows()],
        hoverinfo='text'
    ))
    
    # Format source name for title
    source_display = source_display_name(source)
    direction_display = 'Upregulated' if set_name == 'upregulated' else 'Downregulated'
    revigo_label = ' (REVIGO filtered)' if filename_suffix else ''
    
    # Layout
    fig.update_layout(
        title=dict(
            text=f'{source_display}: {direction_display} Proteins (Dot Plot){revigo_label}',
            font=dict(size=16, family='Arial Black')
        ),
        xaxis=dict(
            title=dict(text='Gene Ratio', font=dict(size=12)),
            gridcolor='rgba(0,0,0,0.1)'
        ),
        yaxis=dict(
            title=dict(text=''),
            tickfont=dict(size=10)
        ),
        plot_bgcolor='white',
        width=850,
        height=max(400, 50 + len(plot_df) * 25),
        margin=dict(l=300, r=100, t=60, b=60)
    )
    
    # Save outputs
    source_name = source_file_token(source)
    set_output_dir = os.path.join(output_dir, set_name)
    
    pdf_path = os.path.join(set_output_dir, f'{set_name}_{source_name}_dotplot{filename_suffix}.pdf')
    html_path = os.path.join(set_output_dir, f'{set_name}_{source_name}_dotplot{filename_suffix}.html')
    
    fig.write_image(pdf_path, scale=2)
    fig.write_html(html_path)
    
    log_message(log_path, f'    Saved: {pdf_path}')
    log_message(log_path, f'    Saved: {html_path}')


def create_results_table(df, source, set_name, output_dir, log_path, filename_suffix=''):
    '''
    Create formatted HTML and PDF tables of enrichment results.
    
    Parameters:
        df (pd.DataFrame): Enrichment results
        source (str): Database source name
        set_name (str): 'upregulated' or 'downregulated'
        output_dir (str): Output directory
        log_path (str): Path to log file
        filename_suffix (str): Optional suffix for output filenames (e.g., '_revigo')
    '''
    if len(df) == 0:
        return
    
    # Format source name for display
    source_display = source_display_name(source)
    direction_display = 'Upregulated' if set_name == 'upregulated' else 'Downregulated'
    header_color = COLOR_UP if set_name == 'upregulated' else COLOR_DOWN
    
    # ---------------------------------------------------------------------
    # Prepare table data
    # ---------------------------------------------------------------------
    table_df = df.copy()
    table_df['fdr_pvalue_fmt'] = table_df['fdr_pvalue'].apply(lambda x: f'{x:.2e}')
    table_df['gene_ratio_fmt'] = table_df['gene_ratio'].apply(lambda x: f'{x:.3f}')
    
    # Add rank column
    table_df = table_df.reset_index(drop=True)
    table_df.index = table_df.index + 1
    table_df['Rank'] = table_df.index
    
    # ---------------------------------------------------------------------
    # Create HTML table
    # ---------------------------------------------------------------------
    html_content = f'''<!DOCTYPE html>
<html>
<head>
    <title>{source_display}: {direction_display} Proteins</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 40px;
            background-color: #f9f9f9;
        }}
        h1 {{
            color: #333;
            border-bottom: 3px solid {header_color};
            padding-bottom: 10px;
        }}
        table {{
            border-collapse: collapse;
            width: 100%;
            background-color: white;
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
        }}
        th {{
            background-color: {header_color};
            color: white;
            padding: 12px 15px;
            text-align: left;
            font-weight: bold;
            position: sticky;
            top: 0;
        }}
        td {{
            padding: 10px 15px;
            border-bottom: 1px solid #ddd;
        }}
        tr:nth-child(even) {{
            background-color: #f8f9fa;
        }}
        tr:hover {{
            background-color: #e8f4f8;
        }}
        .term-name {{
            max-width: 400px;
        }}
        .footer {{
            margin-top: 20px;
            font-size: 12px;
            color: #666;
        }}
        .count {{
            font-size: 14px;
            color: #555;
            margin-bottom: 15px;
        }}
    </style>
</head>
<body>
    <h1>{source_display}: {direction_display} Proteins</h1>
    <p class="count">Showing {len(table_df)} enriched terms (FDR ≤ 0.05)</p>
    <table>
        <thead>
            <tr>
                <th>Rank</th>
                <th>Term ID</th>
                <th>Term Name</th>
                <th>Gene Count</th>
                <th>Gene Ratio</th>
                <th>FDR p-value</th>
            </tr>
        </thead>
        <tbody>
'''
    
    for _, row in table_df.iterrows():
        html_content += f'''            <tr>
                <td>{row['Rank']}</td>
                <td>{row['term_id']}</td>
                <td class="term-name">{row['term_name']}</td>
                <td>{row['gene_count']}</td>
                <td>{row['gene_ratio_fmt']}</td>
                <td>{row['fdr_pvalue_fmt']}</td>
            </tr>
'''
    
    html_content += f'''        </tbody>
    </table>
    <p class="footer">Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")} | 
    Ketamine Proteomics Analysis Project</p>
</body>
</html>
'''
    
    # Save HTML
    source_name = source_file_token(source)
    set_output_dir = os.path.join(output_dir, set_name)
    html_path = os.path.join(set_output_dir, f'{set_name}_{source_name}_table{filename_suffix}.html')
    
    with open(html_path, 'w') as f:
        f.write(html_content)
    
    log_message(log_path, f'    Saved: {html_path}')
    
    # ---------------------------------------------------------------------
    # Create PDF table using plotly
    # ---------------------------------------------------------------------
    fig = go.Figure(data=[go.Table(
        header=dict(
            values=['<b>Rank</b>', '<b>Term ID</b>', '<b>Term Name</b>', 
                    '<b>Genes</b>', '<b>Gene Ratio</b>', '<b>FDR p-value</b>'],
            fill_color=header_color,
            font=dict(color='white', size=10),
            align='left',
            height=30
        ),
        cells=dict(
            values=[
                table_df['Rank'],
                table_df['term_id'],
                table_df['term_name'],
                table_df['gene_count'],
                table_df['gene_ratio_fmt'],
                table_df['fdr_pvalue_fmt']
            ],
            fill_color=[['white', '#f8f9fa'] * (len(table_df) // 2 + 1)],
            font=dict(size=9),
            align='left',
            height=22
        )
    )])
    
    revigo_label = ' (REVIGO filtered)' if filename_suffix else ''
    fig.update_layout(
        title=dict(
            text=f'{source_display}: {direction_display} Proteins{revigo_label}<br><sup>Showing {len(table_df)} enriched terms</sup>',
            font=dict(size=14, family='Arial Black')
        ),
        width=900,
        height=max(400, 80 + len(table_df) * 24),
        margin=dict(t=80, b=20, l=20, r=20)
    )
    
    pdf_path = os.path.join(set_output_dir, f'{set_name}_{source_name}_table{filename_suffix}.pdf')
    fig.write_image(pdf_path, scale=2)
    
    log_message(log_path, f'    Saved: {pdf_path}')


def create_gene_term_network(df, source, set_name, output_dir, log_path, enrichment_results_full, filename_suffix=''):
    '''
    Create gene-term network visualization showing relationships between
    enriched terms and the genes that drive them.
    
    Parameters:
        df (pd.DataFrame): Enrichment results (filtered)
        source (str): Database source name
        set_name (str): 'upregulated' or 'downregulated'
        output_dir (str): Output directory
        log_path (str): Path to log file
        enrichment_results_full (pd.DataFrame): Full gProfiler results with intersections
        filename_suffix (str): Optional suffix for output filenames (e.g., '_revigo')
    '''
    if len(df) == 0:
        return
    
    log_message(log_path, f'    Creating gene-term network for {source}...')
    
    # ---------------------------------------------------------------------
    # Build network from term-gene relationships
    # ---------------------------------------------------------------------
    G = nx.Graph()
    
    # Get term IDs from our filtered results
    term_ids = set(df['term_id'].tolist())
    
    # Check if we have intersection data
    if 'intersections' not in enrichment_results_full.columns:
        log_message(log_path, f'      No intersection data available, skipping network')
        return
    
    # Filter full results to our significant terms
    full_df = enrichment_results_full[enrichment_results_full['native'].isin(term_ids)].copy()
    
    if len(full_df) == 0:
        log_message(log_path, f'      No matching terms found in full results')
        return
    
    # Add term nodes and gene nodes
    term_nodes = []
    gene_nodes = set()
    edges = []
    
    for _, row in full_df.iterrows():
        term_id = row['native']
        term_name = row['name']
        genes = row['intersections']
        
        if genes is None or (isinstance(genes, float) and np.isnan(genes)):
            continue
        
        # Handle different formats of intersection data
        if isinstance(genes, str):
            gene_list = [g.strip() for g in genes.split(',')]
        elif isinstance(genes, list):
            gene_list = genes
        else:
            continue
        
        # Add term node
        term_nodes.append({
            'id': term_id,
            'name': term_name,
            'type': 'term',
            'size': len(gene_list),
            'pvalue': df[df['term_id'] == term_id]['fdr_pvalue'].values[0] if term_id in df['term_id'].values else 0.05
        })
        
        # Add gene nodes and edges
        for gene in gene_list:
            gene_nodes.add(gene)
            edges.append((term_id, gene))
    
    if len(term_nodes) == 0 or len(gene_nodes) == 0:
        log_message(log_path, f'      No valid term-gene relationships found')
        return
    
    # Build networkx graph
    for term in term_nodes:
        G.add_node(term['id'], **term)
    
    for gene in gene_nodes:
        G.add_node(gene, type='gene', name=gene, size=1)
    
    for edge in edges:
        G.add_edge(edge[0], edge[1])
    
    # ---------------------------------------------------------------------
    # Calculate layout using spring layout
    # ---------------------------------------------------------------------
    pos = nx.spring_layout(G, k=2, iterations=50, seed=42)
    
    # ---------------------------------------------------------------------
    # Create plotly figure
    # ---------------------------------------------------------------------
    # Edge traces
    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])
    
    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines'
    )
    
    # Term node trace
    term_x = [pos[t['id']][0] for t in term_nodes]
    term_y = [pos[t['id']][1] for t in term_nodes]
    term_sizes = [max(15, min(50, t['size'] * 3)) for t in term_nodes]
    term_colors = [-np.log10(t['pvalue']) for t in term_nodes]
    term_text = [f"<b>{t['name']}</b><br>ID: {t['id']}<br>Genes: {t['size']}<br>FDR: {t['pvalue']:.2e}" 
                 for t in term_nodes]
    
    term_trace = go.Scatter(
        x=term_x, y=term_y,
        mode='markers',
        hoverinfo='text',
        text=term_text,
        marker=dict(
            size=term_sizes,
            color=term_colors,
            colorscale='Reds' if set_name == 'upregulated' else 'Blues',
            showscale=True,
            colorbar=dict(title='-log₁₀(FDR)', x=1.02),
            line=dict(width=1, color='white')
        ),
        name='Enriched Terms'
    )
    
    # Gene node trace
    gene_list = list(gene_nodes)
    gene_x = [pos[g][0] for g in gene_list]
    gene_y = [pos[g][1] for g in gene_list]
    gene_degrees = [G.degree(g) for g in gene_list]
    gene_text = [f"<b>{g}</b><br>Connected to {G.degree(g)} terms" for g in gene_list]
    
    gene_trace = go.Scatter(
        x=gene_x, y=gene_y,
        mode='markers',
        hoverinfo='text',
        text=gene_text,
        marker=dict(
            size=[max(6, min(20, d * 2)) for d in gene_degrees],
            color='#2CA02C',  # Green for genes
            line=dict(width=0.5, color='white')
        ),
        name='Genes'
    )
    
    # Create figure
    source_display = source_display_name(source)
    direction_display = 'Upregulated' if set_name == 'upregulated' else 'Downregulated'
    
    fig = go.Figure(data=[edge_trace, gene_trace, term_trace])
    
    fig.update_layout(
        title=dict(
            text=f'Gene-Term Network: {source_display} ({direction_display})<br>' +
                 f'<sup>{len(term_nodes)} terms, {len(gene_nodes)} genes</sup>',
            font=dict(size=16, family='Arial Black')
        ),
        showlegend=True,
        legend=dict(x=0, y=1, bgcolor='rgba(255,255,255,0.8)'),
        hovermode='closest',
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        plot_bgcolor='white',
        width=1000,
        height=800,
        margin=dict(t=80, b=20, l=20, r=100)
    )
    
    # Save outputs
    source_name = source_file_token(source)
    set_output_dir = os.path.join(output_dir, set_name)
    
    pdf_path = os.path.join(set_output_dir, f'{set_name}_{source_name}_network{filename_suffix}.pdf')
    html_path = os.path.join(set_output_dir, f'{set_name}_{source_name}_network{filename_suffix}.html')
    
    fig.write_image(pdf_path, scale=2)
    fig.write_html(html_path)
    
    log_message(log_path, f'    Saved: {pdf_path}')
    log_message(log_path, f'    Saved: {html_path}')


def create_term_clustering(df, source, set_name, output_dir, log_path, enrichment_results_full):
    '''
    Create clustered visualization of enriched terms based on gene overlap.
    Terms with similar gene sets are grouped together to reduce redundancy.
    
    Parameters:
        df (pd.DataFrame): Enrichment results (filtered)
        source (str): Database source name
        set_name (str): 'upregulated' or 'downregulated'
        output_dir (str): Output directory
        log_path (str): Path to log file
        enrichment_results_full (pd.DataFrame): Full gProfiler results with intersections
    '''
    if len(df) < 3:  # Need at least 3 terms to cluster
        log_message(log_path, f'    Skipping clustering for {source} (need at least 3 terms)')
        return
    
    log_message(log_path, f'    Creating term clustering for {source}...')
    
    # ---------------------------------------------------------------------
    # Build term-gene matrix for clustering
    # ---------------------------------------------------------------------
    term_ids = df['term_id'].tolist()
    
    # Check if we have intersection data
    if 'intersections' not in enrichment_results_full.columns:
        log_message(log_path, f'      No intersection data available, skipping clustering')
        return
    
    # Get gene sets for each term
    term_genes = {}
    all_genes = set()
    
    for term_id in term_ids:
        row = enrichment_results_full[enrichment_results_full['native'] == term_id]
        if len(row) == 0:
            continue
        
        genes = row['intersections'].values[0]
        
        if genes is None or (isinstance(genes, float) and np.isnan(genes)):
            continue
        
        if isinstance(genes, str):
            gene_set = set(g.strip() for g in genes.split(','))
        elif isinstance(genes, list):
            gene_set = set(genes)
        else:
            continue
        
        term_genes[term_id] = gene_set
        all_genes.update(gene_set)
    
    if len(term_genes) < 3:
        log_message(log_path, f'      Not enough terms with gene data for clustering')
        return
    
    # ---------------------------------------------------------------------
    # Calculate Jaccard similarity between terms
    # ---------------------------------------------------------------------
    term_list = list(term_genes.keys())
    n_terms = len(term_list)
    
    # Create distance matrix (1 - Jaccard similarity)
    distance_matrix = np.zeros((n_terms, n_terms))
    
    for i in range(n_terms):
        for j in range(i + 1, n_terms):
            set_i = term_genes[term_list[i]]
            set_j = term_genes[term_list[j]]
            
            intersection = len(set_i & set_j)
            union = len(set_i | set_j)
            
            jaccard = intersection / union if union > 0 else 0
            distance = 1 - jaccard
            
            distance_matrix[i, j] = distance
            distance_matrix[j, i] = distance
    
    # ---------------------------------------------------------------------
    # Perform hierarchical clustering
    # ---------------------------------------------------------------------
    # Convert to condensed distance matrix
    condensed_dist = pdist(distance_matrix)
    
    # Perform clustering
    linkage_matrix = linkage(condensed_dist, method='average')
    
    # Assign clusters (using distance threshold of 0.7)
    clusters = fcluster(linkage_matrix, t=0.7, criterion='distance')
    
    # ---------------------------------------------------------------------
    # Prepare data for visualization
    # ---------------------------------------------------------------------
    # Map term IDs to their info
    term_info = []
    for i, term_id in enumerate(term_list):
        row = df[df['term_id'] == term_id].iloc[0]
        term_info.append({
            'term_id': term_id,
            'term_name': row['term_name'],
            'fdr_pvalue': row['fdr_pvalue'],
            'gene_count': row['gene_count'],
            'cluster': clusters[i]
        })
    
    term_info_df = pd.DataFrame(term_info)
    term_info_df = term_info_df.sort_values(['cluster', 'fdr_pvalue'])
    
    # Assign cluster labels (most significant term in each cluster)
    cluster_labels = {}
    for cluster_id in term_info_df['cluster'].unique():
        cluster_terms = term_info_df[term_info_df['cluster'] == cluster_id]
        # Use the most significant term as the cluster representative
        representative = cluster_terms.iloc[0]['term_name']
        # Truncate if too long
        if len(representative) > 40:
            representative = representative[:37] + '...'
        cluster_labels[cluster_id] = f"Cluster {cluster_id}: {representative}"
    
    term_info_df['cluster_label'] = term_info_df['cluster'].map(cluster_labels)
    
    # ---------------------------------------------------------------------
    # Create clustered bar plot
    # ---------------------------------------------------------------------
    source_display = source_display_name(source)
    direction_display = 'Upregulated' if set_name == 'upregulated' else 'Downregulated'
    
    # Generate colors for clusters
    n_clusters = len(cluster_labels)
    if set_name == 'upregulated':
        base_colors = ['#D62728', '#FF7F0E', '#E377C2', '#8C564B', '#9467BD', '#17BECF', '#BCBD22']
    else:
        base_colors = ['#1F77B4', '#2CA02C', '#9467BD', '#8C564B', '#E377C2', '#17BECF', '#BCBD22']
    
    cluster_colors = {c: base_colors[i % len(base_colors)] for i, c in enumerate(sorted(cluster_labels.keys()))}
    
    fig = go.Figure()
    
    # Reverse for horizontal bar plot
    plot_df = term_info_df.iloc[::-1]
    
    fig.add_trace(go.Bar(
        x=-np.log10(plot_df['fdr_pvalue']),
        y=plot_df['term_name'],
        orientation='h',
        marker_color=[cluster_colors[c] for c in plot_df['cluster']],
        text=[f"Cluster {c}" for c in plot_df['cluster']],
        textposition='inside',
        textfont=dict(size=8, color='white'),
        hovertemplate='<b>%{y}</b><br>' +
                      '-log₁₀(FDR): %{x:.2f}<br>' +
                      '%{text}<br>' +
                      '<extra></extra>'
    ))
    
    fig.update_layout(
        title=dict(
            text=f'{source_display}: {direction_display} Proteins (Clustered)<br>' +
                 f'<sup>{len(term_info_df)} terms in {n_clusters} clusters (by gene overlap)</sup>',
            font=dict(size=16, family='Arial Black')
        ),
        xaxis=dict(
            title=dict(text='-log₁₀(FDR p-value)', font=dict(size=12)),
            gridcolor='rgba(0,0,0,0.1)'
        ),
        yaxis=dict(
            title=dict(text=''),
            tickfont=dict(size=9)
        ),
        plot_bgcolor='white',
        width=900,
        height=max(500, 50 + len(plot_df) * 22),
        margin=dict(l=350, r=50, t=80, b=60)
    )
    
    # Save outputs
    source_name = source_file_token(source)
    set_output_dir = os.path.join(output_dir, set_name)
    
    pdf_path = os.path.join(set_output_dir, f'{set_name}_{source_name}_clustered.pdf')
    html_path = os.path.join(set_output_dir, f'{set_name}_{source_name}_clustered.html')
    
    fig.write_image(pdf_path, scale=2)
    fig.write_html(html_path)
    
    log_message(log_path, f'    Saved: {pdf_path}')
    log_message(log_path, f'    Saved: {html_path}')
    
    # ---------------------------------------------------------------------
    # Save cluster assignments to CSV
    # ---------------------------------------------------------------------
    cluster_csv_path = os.path.join(set_output_dir, f'{set_name}_{source_name}_clusters.csv')
    term_info_df.to_csv(cluster_csv_path, index=False)
    log_message(log_path, f'    Saved: {cluster_csv_path}')


def wrap_label(label, max_line_len=30):
    '''
    Wrap long labels to two lines for better readability.
    Splits at a space near the middle of the label.
    
    Parameters:
        label (str): The label text to wrap
        max_line_len (int): Maximum characters per line before wrapping
    
    Returns:
        str: Label with <br> tag inserted for line break if needed
    '''
    if len(label) <= max_line_len:
        return label
    
    # Find a good split point (space nearest to middle)
    mid = len(label) // 2
    
    # Search for space near middle, expanding outward
    best_split = None
    for offset in range(mid):
        # Check position after middle
        if mid + offset < len(label) and label[mid + offset] == ' ':
            best_split = mid + offset
            break
        # Check position before middle
        if mid - offset >= 0 and label[mid - offset] == ' ':
            best_split = mid - offset
            break
    
    # If no space found, just split at max_line_len
    if best_split is None:
        best_split = max_line_len
    
    # Insert HTML line break
    return label[:best_split] + '<br>' + label[best_split:].strip()


def create_publication_figure(df, source, set_name, output_dir, log_path, 
                               enrichment_results_full, protein_data, n_top_terms=6):
    '''
    Create publication-style figure inspired by Herzog et al. 2021 Figure 4.
    
    Uses an improved 3-row layout to prevent label overlap:
    - Row 1: Panel A (main heatmap) spanning full width
    - Row 2: Panels B, C, D (gene heatmaps)
    - Row 3: Panels E, F (gene heatmaps)
    
    Parameters:
        df (pd.DataFrame): Enrichment results (filtered)
        source (str): Database source name
        set_name (str): 'upregulated' or 'downregulated'
        output_dir (str): Output directory
        log_path (str): Path to log file
        enrichment_results_full (pd.DataFrame): Full gProfiler results with intersections
        protein_data (pd.DataFrame): Original protein data with fold changes
        n_top_terms (int): Number of top terms to display
    '''
    if len(df) < 2:
        log_message(log_path, f'    Skipping publication figure for {source} (need at least 2 terms)')
        return
    
    log_message(log_path, f'    Creating publication-style figure for {source}...')
    
    # ---------------------------------------------------------------------
    # Get top terms and their genes
    # ---------------------------------------------------------------------
    top_terms = df.head(n_top_terms).copy()
    term_ids = top_terms['term_id'].tolist()
    
    # Check if we have intersection data
    if 'intersections' not in enrichment_results_full.columns:
        log_message(log_path, f'      No intersection data available')
        return
    
    # Get gene sets for each term
    term_genes = {}
    all_genes = set()
    
    for term_id in term_ids:
        row = enrichment_results_full[enrichment_results_full['native'] == term_id]
        if len(row) == 0:
            continue
        
        genes = row['intersections'].values[0]
        
        if genes is None or (isinstance(genes, float) and np.isnan(genes)):
            continue
        
        if isinstance(genes, str):
            gene_list = [g.strip() for g in genes.split(',')]
        elif isinstance(genes, list):
            gene_list = genes
        else:
            continue
        
        term_genes[term_id] = gene_list
        all_genes.update(gene_list)
    
    if len(term_genes) < 2:
        log_message(log_path, f'      Not enough terms with gene data')
        return
    
    # ---------------------------------------------------------------------
    # Order terms by significance (FDR p-value) - NO clustering
    # This ensures Panel A and Panels B-F are in the same order
    # ---------------------------------------------------------------------
    # Terms from top_terms are already sorted by fdr_pvalue (most significant first)
    ordered_term_ids = top_terms['term_id'].tolist()
    ordered_term_names = top_terms['term_name'].tolist()
    ordered_pvalues = top_terms['fdr_pvalue'].tolist()
    
    # Filter to only terms that have gene data available
    filtered_ids = []
    filtered_names = []
    filtered_pvalues = []
    
    for tid, tname, pval in zip(ordered_term_ids, ordered_term_names, ordered_pvalues):
        if tid in term_genes and len(term_genes[tid]) > 0:
            filtered_ids.append(tid)
            filtered_names.append(tname)
            filtered_pvalues.append(pval)
    
    ordered_term_ids = filtered_ids
    ordered_term_names = filtered_names
    ordered_pvalues = filtered_pvalues
    
    log_message(log_path, f'      Terms ordered by significance:')
    for i, (tname, pval) in enumerate(zip(ordered_term_names, ordered_pvalues)):
        log_message(log_path, f'        {i+1}. {tname} (FDR={pval:.2e})')
    
    # ---------------------------------------------------------------------
    # Create figure with improved 3-row layout
    # Row 1: Panel A (main heatmap) - full width
    # Row 2: Panels B, C, D (gene heatmaps)
    # Row 3: Panels E, F (gene heatmaps)
    # ---------------------------------------------------------------------
    n_pathway_panels = min(5, len(ordered_term_ids))  # Max 5 pathway panels
    
    # Build subplot specs for 3-row layout
    specs = [
        [{'colspan': 3, 'type': 'heatmap'}, None, None],
        [{'type': 'heatmap'}] * 3,
        [{'type': 'heatmap'}] * 3
    ]
    
    # --- Subplot titles: descriptive subheading only ---
    # Per style guide section 3, in-figure subheadings are 12 pt bold sentence case.
    # Panel letters (A, B, C, ...) are placed separately at the upper-left corner
    # of each panel via add_annotation() at 15 pt bold (see section 3 + 11.1).
    def _sentence_case(s):
        return s[0].upper() + s[1:] if s else s

    subplot_titles = [f'<b>Top biological functions</b>']

    # Row 2 panels: B, C, D
    for i in range(min(3, n_pathway_panels)):
        name = wrap_label(_sentence_case(ordered_term_names[i]), 25)
        subplot_titles.append(f'<b>{name}</b>')
    # Pad row 2 if fewer than 3 panels
    while len(subplot_titles) < 4:
        subplot_titles.append('')

    # Row 3 panels: E, F
    for i in range(3, n_pathway_panels):
        name = wrap_label(_sentence_case(ordered_term_names[i]), 25)
        subplot_titles.append(f'<b>{name}</b>')
    # Pad row 3
    while len(subplot_titles) < 7:
        subplot_titles.append('')
    
    fig = make_subplots(
        rows=3, cols=3,
        row_heights=[0.32, 0.34, 0.34],
        specs=specs,
        subplot_titles=subplot_titles,
        vertical_spacing=0.14,
        horizontal_spacing=0.10
    )
    
    # ---------------------------------------------------------------------
    # Panel A: Top biological functions heatmap
    # Reverse order so most significant appears at TOP (Plotly y-axis is bottom-to-top)
    # ---------------------------------------------------------------------
    # Wrap term names for display (allows two lines)
    display_term_names = [wrap_label(name, 35) for name in reversed(ordered_term_names)]
    heatmap_values = [[-np.log10(p)] for p in reversed(ordered_pvalues)]

    # Align the -log10(FDR) colorbar's top edge and full span to Panel A's
    # plotting domain (the central "Top biological functions" heatmap).
    # make_subplots() has already set yaxis.domain for subplot 1, so we read it here.
    panel_a_y_domain = fig.layout.yaxis.domain
    panel_a_top = panel_a_y_domain[1]
    panel_a_height = panel_a_y_domain[1] - panel_a_y_domain[0]

    fig.add_trace(
        go.Heatmap(
            z=heatmap_values,
            y=display_term_names,
            x=['Ketamine vs Control'],
            colorscale='Purples',
            showscale=True,
            colorbar=dict(
                title=dict(
                    text='-log<sub>10</sub>(FDR)',
                    font=dict(size=14, color='#000000', family='Arial')
                ),
                tickfont=dict(size=12, color='#000000', family='Arial'),
                x=1.02,
                len=panel_a_height,
                y=panel_a_top,
                yanchor='top',
                thickness=15
            ),
            hovertemplate='<b>%{y}</b><br>-log<sub>10</sub>(FDR): %{z:.2f}<extra></extra>'
        ),
        row=1, col=1
    )
    
    # ---------------------------------------------------------------------
    # Panels B-F: Gene-level heatmaps for each pathway
    # ---------------------------------------------------------------------
    # Get log2FC for all genes from protein data
    gene_fc = {}
    if 'Gene Symbol' in protein_data.columns and 'log2_fold_change' in protein_data.columns:
        for _, row in protein_data.iterrows():
            gene = row['Gene Symbol']
            if pd.notna(gene) and gene != '':
                gene_fc[gene] = row['log2_fold_change']
    
    for panel_idx in range(n_pathway_panels):
        term_id = ordered_term_ids[panel_idx]
        genes = term_genes.get(term_id, [])
        
        if len(genes) == 0:
            continue
        
        # Get fold changes for genes in this pathway
        gene_values = []
        gene_names = []
        
        for gene in genes[:15]:  # Limit to 15 genes per panel
            fc = gene_fc.get(gene, 0)
            if not np.isnan(fc):
                gene_values.append([fc])
                gene_names.append(gene)
        
        if len(gene_values) == 0:
            continue
        
        # Sort by fold change: ascending order so HIGHEST appears at TOP of heatmap
        # (Plotly heatmaps display y-axis from bottom to top)
        sorted_data = sorted(zip(gene_names, gene_values), key=lambda x: x[1][0], reverse=False)
        gene_names = [x[0] for x in sorted_data]
        gene_values = [x[1] for x in sorted_data]

        # In a proteomics fold-change figure, the symbol labels the protein product
        # detected in mass spec (style guide section 4.3): set in uppercase roman.
        gene_names = [g.upper() for g in gene_names]

        # Determine row and column in 3-row layout
        if panel_idx < 3:
            plot_row = 2
            plot_col = panel_idx + 1
        else:
            plot_row = 3
            plot_col = panel_idx - 3 + 1

        # Only show colorbar on the last panel
        show_colorbar = (panel_idx == n_pathway_panels - 1)

        # Project palette (style guide section 9.1.1):
        #   down-in-ketamine -> control color #7FB3D8 (light blue)
        #   up-in-ketamine   -> ketamine color #E8735A (coral)
        # Replaces the previous green/red diverging scale, which violated
        # the section 9.2 ban on red/green as the only distinguishing channel.
        fig.add_trace(
            go.Heatmap(
                z=gene_values,
                y=gene_names,
                x=['log<sub>2</sub>FC'],
                colorscale=[[0, '#7FB3D8'], [0.5, '#FFFFFF'], [1, '#E8735A']],
                zmid=0,
                zmin=-4,
                zmax=4,
                showscale=show_colorbar,
                colorbar=dict(
                    title=dict(
                        text='log<sub>2</sub>FC',
                        font=dict(size=14, color='#000000', family='Arial')
                    ),
                    tickfont=dict(size=12, color='#000000', family='Arial'),
                    x=1.02,
                    len=0.28,
                    y=0.17,
                    thickness=15
                ) if show_colorbar else None,
                hovertemplate='<b>%{y}</b><br>log<sub>2</sub>FC: %{z:.2f}<extra></extra>'
            ),
            row=plot_row, col=plot_col
        )
    
    # ---------------------------------------------------------------------
    # Update layout
    # ---------------------------------------------------------------------
    # Per style guide section 11.2, no on-figure title is set. The figure
    # number, title sentence, and "Upregulated proteins in Ketamine vs Control"
    # statement all live in the Word document caption block below the figure.
    fig.update_layout(
        height=1300,
        width=1100,
        plot_bgcolor='white',
        paper_bgcolor='white',
        font=dict(family='Arial', size=12, color='#000000'),
        margin=dict(l=180, r=80, t=60, b=50)
    )

    # Tick labels: 12 pt black Arial across all panels (style guide section 3)
    fig.update_yaxes(tickfont=dict(size=12, color='#000000', family='Arial'))
    fig.update_xaxes(tickfont=dict(size=12, color='#000000', family='Arial'))

    # Subplot-title annotations: 12 pt bold sentence case in black
    # (Section 3 figure-interior subheading; bold is supplied by the <b> wrapper
    # in subplot_titles, since plotly's annotation font dict has no weight key.)
    for annotation in fig['layout']['annotations']:
        annotation['font'] = dict(size=12, color='#000000', family='Arial')

    # Panel letters (A, B, C, ...): 15 pt bold black at upper-left of each
    # panel's plotting domain (style guide section 3).
    panel_axis_suffixes = ['', '2', '3', '4', '5', '6']
    panel_letters = ['A', 'B', 'C', 'D', 'E', 'F']
    n_total_panels = 1 + n_pathway_panels  # Panel A + pathway panels
    for i in range(min(6, n_total_panels)):
        ax_suffix = panel_axis_suffixes[i]
        x_domain = fig.layout[f'xaxis{ax_suffix}'].domain
        y_domain = fig.layout[f'yaxis{ax_suffix}'].domain
        fig.add_annotation(
            x=x_domain[0] - 0.005,
            y=y_domain[1] + 0.015,
            xref='paper',
            yref='paper',
            text=f'<b>{panel_letters[i]}</b>',
            showarrow=False,
            font=dict(size=15, color='#000000', family='Arial'),
            xanchor='right',
            yanchor='bottom'
        )
    
    # ---------------------------------------------------------------------
    # Save outputs
    # ---------------------------------------------------------------------
    source_name = source_file_token(source)
    set_output_dir = os.path.join(output_dir, set_name)
    
    pdf_path = os.path.join(set_output_dir, f'{set_name}_{source_name}_publication_figure.pdf')
    html_path = os.path.join(set_output_dir, f'{set_name}_{source_name}_publication_figure.html')
    
    fig.write_image(pdf_path, scale=2)
    fig.write_html(html_path)
    
    log_message(log_path, f'    Saved: {pdf_path}')
    log_message(log_path, f'    Saved: {html_path}')


# =============================================================================
# REPORT GENERATION
# =============================================================================

def generate_set_summary(results, gene_list, set_name, output_dir, fdr_threshold, log_path):
    '''
    Generate summary report for a gene set.
    '''
    set_output_dir = os.path.join(output_dir, set_name)
    summary_path = os.path.join(set_output_dir, f'{set_name}_summary.txt')
    
    with open(summary_path, 'w') as f:
        f.write('=' * 70 + '\n')
        f.write(f'PATHWAY ANALYSIS SUMMARY: {set_name.upper()}\n')
        f.write('=' * 70 + '\n\n')
        f.write(f'Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n')
        f.write(f'Input genes: {len(gene_list)}\n')
        f.write(f'FDR threshold: {fdr_threshold}\n\n')
        
        for source in SOURCES:
            source_display = source_display_name(source)
            f.write(f'\n{"-" * 50}\n')
            f.write(f'{source_display}\n')
            f.write(f'{"-" * 50}\n')
            
            # Handle new results structure
            result_data = results.get(source, {'clean': pd.DataFrame()})
            df = result_data.get('clean', pd.DataFrame()) if isinstance(result_data, dict) else result_data
            
            if len(df) == 0:
                f.write('No significantly enriched terms found.\n')
            else:
                f.write(f'Significant terms: {len(df)}\n\n')
                f.write('Top 10 terms:\n')
                
                for i, (_, row) in enumerate(df.head(10).iterrows(), 1):
                    f.write(f"\n{i}. {row['term_name']}\n")
                    f.write(f"   ID: {row['term_id']}\n")
                    f.write(f"   FDR p-value: {row['fdr_pvalue']:.2e}\n")
                    f.write(f"   Genes: {row['gene_count']}/{row.get('query_size', 'N/A')}\n")
    
    log_message(log_path, f'  Saved: {summary_path}')


def generate_final_report(gene_sets, all_results, output_dir, pval_threshold, fdr_threshold, log_path):
    '''
    Generate final comprehensive report.
    '''
    report_path = os.path.join(output_dir, 'pathway_analysis_report.txt')
    
    with open(report_path, 'w') as f:
        f.write('=' * 70 + '\n')
        f.write('KETAMINE PROTEOMICS ANALYSIS PROJECT\n')
        f.write('PATHWAY ENRICHMENT ANALYSIS REPORT\n')
        f.write('=' * 70 + '\n\n')
        
        f.write(f'Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n\n')
        
        # ---------------------------------------------------------------------
        # Analysis parameters
        # ---------------------------------------------------------------------
        f.write('-' * 70 + '\n')
        f.write('ANALYSIS PARAMETERS\n')
        f.write('-' * 70 + '\n')
        f.write(f'Method: Over-representation analysis (ORA)\n')
        f.write(f'Tool: gProfiler\n')
        f.write(f'Organism: Mus musculus (mouse)\n')
        f.write(f'Databases: GO Biological Process, KEGG\n')
        f.write(f'Significance threshold (differential abundance): p <= {pval_threshold}\n')
        f.write(f'Significance threshold (pathway enrichment): FDR <= {fdr_threshold}\n')
        f.write(f'Multiple testing correction: Benjamini-Hochberg\n')
        f.write(f'Background set: All detected proteins ({len(gene_sets["background"])})\n\n')
        
        # ---------------------------------------------------------------------
        # Gene set composition
        # ---------------------------------------------------------------------
        f.write('-' * 70 + '\n')
        f.write('GENE SET COMPOSITION\n')
        f.write('-' * 70 + '\n\n')
        
        counts = gene_sets['counts']
        
        f.write('UPREGULATED SET:\n')
        f.write(f'  Significant upregulated (p <= {pval_threshold}): {counts["sig_up"]}\n')
        f.write(f'  Ketamine-specific (presence/absence): {counts["ket_specific"]}\n')
        f.write(f'  Total: {len(gene_sets["upregulated"])}\n\n')
        
        f.write('DOWNREGULATED SET:\n')
        f.write(f'  Significant downregulated (p <= {pval_threshold}): {counts["sig_down"]}\n')
        f.write(f'  Control-specific (presence/absence): {counts["ctrl_specific"]}\n')
        f.write(f'  Total: {len(gene_sets["downregulated"])}\n\n')
        
        # ---------------------------------------------------------------------
        # Results summary
        # ---------------------------------------------------------------------
        f.write('-' * 70 + '\n')
        f.write('RESULTS SUMMARY\n')
        f.write('-' * 70 + '\n\n')
        
        for set_name in ['upregulated', 'downregulated']:
            f.write(f'{set_name.upper()}:\n')
            results = all_results.get(set_name, {})
            
            for source in SOURCES:
                source_display = 'GO:BP' if source == 'GO:BP' else 'KEGG'
                result_data = results.get(source, {'clean': pd.DataFrame()})
                df = result_data.get('clean', pd.DataFrame()) if isinstance(result_data, dict) else result_data
                n_terms = len(df) if len(df) > 0 else 0
                f.write(f'  {source_display}: {n_terms} enriched terms\n')
            f.write('\n')
        
        # ---------------------------------------------------------------------
        # Top findings - GO:BP
        # ---------------------------------------------------------------------
        f.write('-' * 70 + '\n')
        f.write('TOP FINDINGS - GO BIOLOGICAL PROCESS\n')
        f.write('-' * 70 + '\n\n')

        for set_name in ['upregulated', 'downregulated']:
            f.write(f'{set_name.upper()} - Top 10 GO:BP terms:\n')
            results = all_results.get(set_name, {})
            result_data = results.get('GO:BP', {'clean': pd.DataFrame(), 'full': pd.DataFrame()})
            # Use 'full' results to get intersection data (gene lists)
            df_full = result_data.get('full', pd.DataFrame()) if isinstance(result_data, dict) else result_data
            df_clean = result_data.get('clean', pd.DataFrame()) if isinstance(result_data, dict) else result_data

            if len(df_clean) > 0 and len(df_full) > 0:
                # Get top 10 term IDs from clean, then look up full data
                top_term_ids = df_clean.head(10)['term_id'].tolist() if 'term_id' in df_clean.columns else []

                for i, (_, row) in enumerate(df_clean.head(10).iterrows(), 1):
                    term_name = row['term_name']
                    term_id = row.get('term_id', '')
                    fdr = row['fdr_pvalue']
                    gene_count = row.get('intersection_size', row.get('gene_count', 'N/A'))

                    # Get genes from full results
                    genes_str = ''
                    if term_id and len(df_full) > 0:
                        full_row = df_full[df_full['native'] == term_id]
                        if len(full_row) > 0:
                            genes = full_row.iloc[0].get('intersections', '')
                            if isinstance(genes, list):
                                genes_str = ', '.join(sorted(genes))
                            elif isinstance(genes, str) and genes:
                                genes_str = ', '.join(sorted(genes.split(',')))

                    f.write(f"  {i:2d}. {term_name}\n")
                    f.write(f"      FDR = {fdr:.2e}, Gene count = {gene_count}\n")
                    if genes_str:
                        f.write(f"      Genes: {genes_str}\n")
            else:
                f.write('  No significant terms found.\n')
            f.write('\n')

        # ---------------------------------------------------------------------
        # Top findings - KEGG
        # ---------------------------------------------------------------------
        f.write('-' * 70 + '\n')
        f.write('TOP FINDINGS - KEGG PATHWAYS\n')
        f.write('-' * 70 + '\n\n')

        for set_name in ['upregulated', 'downregulated']:
            f.write(f'{set_name.upper()} - All KEGG pathways:\n')
            results = all_results.get(set_name, {})
            result_data = results.get('KEGG', {'clean': pd.DataFrame(), 'full': pd.DataFrame()})
            # Use 'full' results to get intersection data (gene lists)
            df_full = result_data.get('full', pd.DataFrame()) if isinstance(result_data, dict) else result_data
            df_clean = result_data.get('clean', pd.DataFrame()) if isinstance(result_data, dict) else result_data

            if len(df_clean) > 0 and len(df_full) > 0:
                for i, (_, row) in enumerate(df_clean.iterrows(), 1):
                    term_name = row['term_name']
                    term_id = row.get('term_id', row.get('native', 'N/A'))
                    fdr = row['fdr_pvalue']
                    gene_count = row.get('intersection_size', row.get('gene_count', 'N/A'))

                    # Get genes from full results
                    genes_str = ''
                    if term_id and len(df_full) > 0:
                        full_row = df_full[df_full['native'] == term_id]
                        if len(full_row) > 0:
                            genes = full_row.iloc[0].get('intersections', '')
                            if isinstance(genes, list):
                                genes_str = ', '.join(sorted(genes))
                            elif isinstance(genes, str) and genes:
                                genes_str = ', '.join(sorted(genes.split(',')))

                    f.write(f"  {i}. {term_name} ({term_id})\n")
                    f.write(f"     FDR = {fdr:.2e}, Gene count = {gene_count}\n")
                    if genes_str:
                        f.write(f"     Genes: {genes_str}\n")
            else:
                f.write('  No significant pathways found.\n')
            f.write('\n')

        # ---------------------------------------------------------------------
        # Detailed gene lists
        # ---------------------------------------------------------------------
        f.write('-' * 70 + '\n')
        f.write('DETAILED GENE INFORMATION\n')
        f.write('-' * 70 + '\n\n')

        for set_name in ['upregulated', 'downregulated']:
            gene_list = gene_sets.get(set_name, [])
            f.write(f'{set_name.upper()} GENES ({len(gene_list)} total):\n')
            if gene_list:
                # Format genes in columns
                genes_sorted = sorted(gene_list)
                line = '  '
                for i, gene in enumerate(genes_sorted):
                    if i > 0 and i % 8 == 0:
                        f.write(line.rstrip(', ') + '\n')
                        line = '  '
                    line += f'{gene}, '
                if line.strip():
                    f.write(line.rstrip(', ') + '\n')
            else:
                f.write('  None\n')
            f.write('\n')

        # ---------------------------------------------------------------------
        # Biological interpretation notes
        # ---------------------------------------------------------------------
        f.write('-' * 70 + '\n')
        f.write('INTERPRETATION NOTES\n')
        f.write('-' * 70 + '\n\n')

        # Check for synaptic/vesicular enrichment
        results = all_results.get('upregulated', {})
        gobp_data = results.get('GO:BP', {'clean': pd.DataFrame()})
        gobp_df = gobp_data.get('clean', pd.DataFrame()) if isinstance(gobp_data, dict) else gobp_data
        kegg_data = results.get('KEGG', {'clean': pd.DataFrame()})
        kegg_df = kegg_data.get('clean', pd.DataFrame()) if isinstance(kegg_data, dict) else kegg_data

        if len(gobp_df) > 0 or len(kegg_df) > 0:
            f.write('Key biological themes identified:\n\n')

            # Check for synaptic terms
            synaptic_terms = []
            if len(gobp_df) > 0:
                synaptic_mask = gobp_df['term_name'].str.contains('synap|vesicle|exocyt', case=False, na=False)
                synaptic_terms = gobp_df[synaptic_mask]['term_name'].head(5).tolist()

            if synaptic_terms:
                f.write('  1. SYNAPTIC/VESICULAR FUNCTION:\n')
                f.write('     Strong enrichment for synaptic vesicle and exocytosis pathways.\n')
                f.write('     This suggests ketamine treatment enhances vesicular release\n')
                f.write('     machinery in astrocytes.\n\n')

            # Check for metabolic terms
            if len(kegg_df) > 0:
                metabolic_mask = kegg_df['term_name'].str.contains('metabol|oxidative|phosphoryl', case=False, na=False)
                if metabolic_mask.any():
                    f.write('  2. METABOLIC CHANGES:\n')
                    f.write('     Enrichment in metabolic pathways indicates potential changes\n')
                    f.write('     in cellular energy production or utilization.\n\n')

            # Note about spurious pathways
            if len(kegg_df) > 0:
                spurious = ['Rheumatoid arthritis', 'Collecting duct']
                spurious_found = [t for t in kegg_df['term_name'].tolist()
                                  if any(s in t for s in spurious)]
                if spurious_found:
                    f.write('  NOTE ON SPURIOUS PATHWAYS:\n')
                    f.write(f'     The following may be false positives due to gene promiscuity:\n')
                    for sp in spurious_found:
                        f.write(f'     - {sp}\n')
                    f.write('     These pathways share genes with synaptic/secretory processes\n')
                    f.write('     but are not biologically relevant in astrocyte context.\n\n')
        else:
            f.write('No significant enrichment detected for interpretation.\n\n')

        # ---------------------------------------------------------------------
        # File locations
        # ---------------------------------------------------------------------
        f.write('-' * 70 + '\n')
        f.write('OUTPUT FILE LOCATIONS\n')
        f.write('-' * 70 + '\n\n')
        f.write(f'Results directory: {output_dir}\n\n')
        f.write('Upregulated results:\n')
        f.write('  - upregulated/upregulated_GO_BP.csv\n')
        f.write('  - upregulated/upregulated_KEGG.csv\n')
        f.write('  - upregulated/upregulated_*_publication_figure.pdf\n\n')
        f.write('Downregulated results:\n')
        f.write('  - downregulated/downregulated_GO_BP.csv\n')
        f.write('  - downregulated/downregulated_KEGG.csv\n\n')
        f.write('Visualizations (PDF and interactive HTML):\n')
        f.write('  - *_plot.pdf/html      : Bar plots\n')
        f.write('  - *_dotplot.pdf/html   : Dot plots\n')
        f.write('  - *_table.pdf/html     : Interactive tables\n')
        f.write('  - *_network.pdf/html   : Gene-term networks\n')
        f.write('  - *_clustered.pdf/html : Clustered heatmaps\n')

    log_message(log_path, f'Saved final report: {report_path}')


# =============================================================================
# REVIGO INTEGRATION FUNCTIONS
# =============================================================================

def parse_revigo_output(revigo_path, log_path):
    '''
    Parse REVIGO output table to get list of non-redundant GO term IDs.
    
    REVIGO output format:
    - TermID: GO term ID (with quotes)
    - Name: Term name
    - Value: -log10(p-value) as negative number
    - Representative: NaN for representative terms, GO ID (as number) for redundant terms
    
    Parameters:
        revigo_path (str): Path to REVIGO TSV output file
        log_path (str): Path to log file
    
    Returns:
        set: Set of non-redundant GO term IDs to keep
    '''
    log_message(log_path, f'Parsing REVIGO output: {revigo_path}')
    
    df = pd.read_csv(revigo_path, sep='\t')
    
    # Clean up quotes in TermID column
    df['TermID'] = df['TermID'].str.replace('"', '')
    
    # Representatives have NaN in the Representative column
    representatives = df[df['Representative'].isna()]
    redundant = df[df['Representative'].notna()]
    
    # Get set of non-redundant term IDs
    keep_terms = set(representatives['TermID'].tolist())
    
    log_message(log_path, f'  Total terms in REVIGO output: {len(df)}')
    log_message(log_path, f'  Non-redundant representatives: {len(representatives)}')
    log_message(log_path, f'  Redundant terms removed: {len(redundant)}')
    log_message(log_path, f'  Reduction: {len(df)} → {len(keep_terms)} ({100*len(redundant)/len(df):.1f}% removed)')
    
    return keep_terms


def run_revigo_mode(args, log_path):
    '''
    Run pathway analysis in REVIGO filtering mode.
    
    This mode:
    1. Reads existing pathway analysis results
    2. Filters GO:BP terms to keep only REVIGO non-redundant representatives
    3. Regenerates all visualizations with filtered terms
    4. Saves results to a 'revigo' subdirectory
    
    Parameters:
        args: Parsed command-line arguments
        log_path (str): Path to log file
    '''
    log_message(log_path, '')
    log_message(log_path, '=' * 50)
    log_message(log_path, 'REVIGO FILTERING MODE')
    log_message(log_path, '=' * 50)
    
    # -------------------------------------------------------------------------
    # Validate inputs
    # -------------------------------------------------------------------------
    if not args.original_results:
        raise ValueError('--original_results is required for REVIGO mode')
    
    if not os.path.exists(args.revigo):
        raise FileNotFoundError(f'REVIGO file not found: {args.revigo}')
    
    if not os.path.exists(args.original_results):
        raise FileNotFoundError(f'Original results directory not found: {args.original_results}')
    
    # -------------------------------------------------------------------------
    # Parse REVIGO output to get non-redundant terms
    # -------------------------------------------------------------------------
    keep_terms = parse_revigo_output(args.revigo, log_path)
    
    # -------------------------------------------------------------------------
    # Setup output directory
    # -------------------------------------------------------------------------
    revigo_output_dir = os.path.join(args.output_dir, 'revigo')
    os.makedirs(revigo_output_dir, exist_ok=True)
    log_message(log_path, f'  Output directory: {revigo_output_dir}')
    
    # -------------------------------------------------------------------------
    # Process each gene set (upregulated/downregulated)
    # -------------------------------------------------------------------------
    for set_name in ['upregulated', 'downregulated']:
        set_input_dir = os.path.join(args.original_results, set_name)
        
        if not os.path.exists(set_input_dir):
            log_message(log_path, f'\n  Skipping {set_name} (directory not found)')
            continue
        
        log_message(log_path, f'\n  Processing {set_name}...')
        
        # Create output subdirectory
        set_output_dir = os.path.join(revigo_output_dir, set_name)
        os.makedirs(set_output_dir, exist_ok=True)
        
        # ---------------------------------------------------------------------
        # Load original GO:BP results
        # ---------------------------------------------------------------------
        go_bp_path = os.path.join(set_input_dir, f'{set_name}_GO_BP.csv')
        
        if not os.path.exists(go_bp_path):
            log_message(log_path, f'    GO:BP results not found: {go_bp_path}')
            continue
        
        df_original = pd.read_csv(go_bp_path)
        log_message(log_path, f'    Loaded {len(df_original)} original GO:BP terms')
        
        # ---------------------------------------------------------------------
        # Filter to keep only REVIGO non-redundant terms
        # ---------------------------------------------------------------------
        df_filtered = df_original[df_original['term_id'].isin(keep_terms)].copy()
        log_message(log_path, f'    After REVIGO filtering: {len(df_filtered)} terms')
        
        if len(df_filtered) == 0:
            log_message(log_path, f'    No terms remaining after filtering!')
            continue
        
        # Save filtered results
        filtered_csv_path = os.path.join(set_output_dir, f'{set_name}_GO_BP_revigo.csv')
        df_filtered.to_csv(filtered_csv_path, index=False)
        log_message(log_path, f'    Saved: {filtered_csv_path}')
        
        # ---------------------------------------------------------------------
        # Load protein data for visualizations (needed for publication figure)
        # ---------------------------------------------------------------------
        protein_data = None
        if args.input and os.path.exists(args.input):
            protein_data = pd.read_csv(args.input)
            log_message(log_path, f'    Loaded protein data: {len(protein_data)} proteins')
        
        # ---------------------------------------------------------------------
        # Regenerate visualizations with filtered terms
        # ---------------------------------------------------------------------
        log_message(log_path, f'    Creating visualizations...')
        
        source = 'GO:BP'
        
        # Bar plot
        create_enrichment_bar_plot(
            df_filtered, source, set_name, revigo_output_dir, log_path,
            filename_suffix='_revigo'
        )
        
        # Dot plot
        create_dot_plot(
            df_filtered, source, set_name, revigo_output_dir, log_path,
            filename_suffix='_revigo'
        )
        
        # Results table
        create_results_table(
            df_filtered, source, set_name, revigo_output_dir, log_path,
            filename_suffix='_revigo'
        )
        
        # -----------------------------------------------------------------
        # Gene-term network, clustering, and publication figure require
        # the full g:Profiler results (with 'intersections' column).
        # The full results were not saved to disk in normal mode, so
        # re-query g:Profiler to obtain them.
        # -----------------------------------------------------------------
        df_full_filtered = None

        go_bp_full_path = os.path.join(set_input_dir, f'{set_name}_GO_BP_full.csv')
        if os.path.exists(go_bp_full_path):
            # If a full CSV was saved, use it
            df_full = pd.read_csv(go_bp_full_path)
            df_full_filtered = df_full[df_full['native'].isin(keep_terms)].copy()
        elif protein_data is not None:
            # Re-query g:Profiler to get intersection data
            log_message(log_path, f'    Full results CSV not found; re-querying g:Profiler for intersection data...')
            gene_sets = load_and_prepare_gene_sets(args.input, 0.05, log_path)
            gene_list = gene_sets.get(set_name, [])
            background = gene_sets.get('background', [])
            if len(gene_list) > 0:
                df_full_query = run_gprofiler_enrichment(gene_list, background, source, log_path)
                if len(df_full_query) > 0:
                    df_full_filtered = df_full_query[df_full_query['native'].isin(keep_terms)].copy()
                    log_message(log_path, f'    Re-query successful: {len(df_full_filtered)} REVIGO terms with intersection data')

        if df_full_filtered is not None and len(df_full_filtered) > 0:
            # Gene-term network
            create_gene_term_network(
                df_filtered, source, set_name, revigo_output_dir, log_path,
                df_full_filtered, filename_suffix='_revigo'
            )

            # Term clustering
            create_term_clustering(
                df_filtered, source, set_name, revigo_output_dir, log_path,
                df_full_filtered
            )

            # Publication figure
            if protein_data is not None:
                create_publication_figure(
                    df_filtered, source, set_name, revigo_output_dir, log_path,
                    df_full_filtered, protein_data
                )
        else:
            log_message(log_path, f'    Skipping network/clustering/publication figures (no intersection data available)')
        
        # ---------------------------------------------------------------------
        # Also copy KEGG results (unchanged by REVIGO)
        # ---------------------------------------------------------------------
        kegg_path = os.path.join(set_input_dir, f'{set_name}_KEGG.csv')
        if os.path.exists(kegg_path):
            import shutil
            kegg_dest = os.path.join(set_output_dir, f'{set_name}_KEGG.csv')
            shutil.copy(kegg_path, kegg_dest)
            log_message(log_path, f'    Copied KEGG results (unchanged by REVIGO)')
    
    # -------------------------------------------------------------------------
    # Generate summary report
    # -------------------------------------------------------------------------
    summary_path = os.path.join(revigo_output_dir, 'revigo_filtering_summary.txt')
    with open(summary_path, 'w') as f:
        f.write('=' * 70 + '\n')
        f.write('REVIGO FILTERING SUMMARY\n')
        f.write('=' * 70 + '\n\n')
        f.write(f'Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n')
        f.write(f'REVIGO input: {args.revigo}\n')
        f.write(f'Original results: {args.original_results}\n')
        f.write(f'Output directory: {revigo_output_dir}\n\n')
        f.write(f'Non-redundant terms retained: {len(keep_terms)}\n\n')
        f.write('This analysis uses REVIGO (Supek et al., 2011) to reduce\n')
        f.write('redundancy among GO terms by clustering semantically similar\n')
        f.write('terms and retaining representative terms.\n\n')
        f.write('Reference:\n')
        f.write('Supek F, Bošnjak M, Škunca N, Šmuc T. REVIGO summarizes and\n')
        f.write('visualizes long lists of Gene Ontology terms. PLoS ONE 2011;6(7):e21800.\n')
    
    log_message(log_path, f'\n  Saved summary: {summary_path}')
    log_message(log_path, '')
    log_message(log_path, '=' * 50)
    log_message(log_path, 'REVIGO filtering complete!')
    log_message(log_path, '=' * 50)


# =============================================================================
# MAIN FUNCTION
# =============================================================================

def main():
    '''Main entry point for pathway analysis.'''
    
    # ---------------------------------------------------------------------
    # Parse arguments
    # ---------------------------------------------------------------------
    parser = argparse.ArgumentParser(
        description='Pathway enrichment analysis for differential abundance results'
    )
    parser.add_argument(
        '--input', '-i',
        type=str,
        required=False,
        help='Path to all_proteins_categorized.csv (required for normal mode)'
    )
    parser.add_argument(
        '--output_dir', '-o',
        type=str,
        required=True,
        help='Output directory for pathway analysis results'
    )
    parser.add_argument(
        '--log_dir', '-l',
        type=str,
        required=True,
        help='Directory for log files'
    )
    parser.add_argument(
        '--pval_threshold', '-p',
        type=float,
        default=0.05,
        help='P-value threshold for differential abundance (default: 0.05)'
    )
    parser.add_argument(
        '--fdr_threshold', '-f',
        type=float,
        default=0.05,
        help='FDR threshold for pathway enrichment (default: 0.05)'
    )
    
    # Source selection / targeted-run arguments
    parser.add_argument(
        '--sources',
        nargs='+',
        default=['GO:BP', 'KEGG'],
        help=("gProfiler source tokens to query (default: GO:BP KEGG). Pass "
              "REAC to add Reactome. Example: --sources REAC")
    )
    parser.add_argument(
        '--enrich-only',
        action='store_true',
        help=("Run enrichment and write CSVs only; skip the legacy in-script "
              "figures, per-set summaries, and the final report. Use this to "
              "add a source (e.g. REAC) without clobbering the shared report "
              "and summary files.")
    )
    parser.add_argument(
        '--full-only',
        action='store_true',
        help=("Write only the <set>_<token>_full.csv files (skip the clean "
              "CSVs). Implies --enrich-only. Use to backfill _full files for "
              "GO:BP/KEGG without overwriting their committed clean CSVs.")
    )

    # REVIGO mode arguments
    parser.add_argument(
        '--revigo', '-r',
        type=str,
        default=None,
        help='Path to REVIGO output table (TSV). Enables REVIGO filtering mode.'
    )
    parser.add_argument(
        '--original_results',
        type=str,
        default=None,
        help='Path to original pathway_analysis results directory (required for REVIGO mode)'
    )
    
    args = parser.parse_args()

    # ---------------------------------------------------------------------
    # Apply source selection + targeted-run semantics
    #   --full-only implies --enrich-only (data-only backfill).
    #   enrich_only => skip legacy figures, per-set summaries, final report so
    #   targeted runs (e.g. adding REAC) do not clobber shared GO:BP/KEGG
    #   outputs.
    # ---------------------------------------------------------------------
    global SOURCES
    SOURCES = args.sources
    enrich_only = args.enrich_only or args.full_only
    write_clean = not args.full_only

    # ---------------------------------------------------------------------
    # Setup
    # ---------------------------------------------------------------------
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(args.log_dir, exist_ok=True)
    log_path = setup_logging(args.log_dir)
    
    # ---------------------------------------------------------------------
    # Check for REVIGO mode
    # ---------------------------------------------------------------------
    if args.revigo:
        # Run in REVIGO filtering mode
        log_message(log_path, '=' * 50)
        log_message(log_path, 'KETAMINE PROTEOMICS ANALYSIS PROJECT')
        log_message(log_path, 'Pathway Enrichment Analysis - REVIGO MODE')
        log_message(log_path, '=' * 50)
        log_message(log_path, f'REVIGO input:      {args.revigo}')
        log_message(log_path, f'Original results:  {args.original_results}')
        log_message(log_path, f'Output directory:  {args.output_dir}')
        
        run_revigo_mode(args, log_path)
        return
    
    # ---------------------------------------------------------------------
    # Normal mode - validate inputs
    # ---------------------------------------------------------------------
    if not args.input:
        raise ValueError('--input is required for normal mode (use --revigo for REVIGO filtering mode)')
    
    log_message(log_path, '=' * 50)
    log_message(log_path, 'KETAMINE PROTEOMICS ANALYSIS PROJECT')
    log_message(log_path, 'Pathway Enrichment Analysis')
    log_message(log_path, '=' * 50)
    log_message(log_path, f'Input file:        {args.input}')
    log_message(log_path, f'Output directory:  {args.output_dir}')
    log_message(log_path, f'P-value threshold: {args.pval_threshold}')
    log_message(log_path, f'FDR threshold:     {args.fdr_threshold}')
    log_message(log_path, f'Sources:           {", ".join(SOURCES)}')
    if enrich_only:
        log_message(log_path, f'Mode:              enrich-only '
                              f'(write_clean={write_clean}; figures/summary/'
                              f'report skipped)')
    log_message(log_path, '')

    # ---------------------------------------------------------------------
    # Load and prepare gene sets
    # ---------------------------------------------------------------------
    gene_sets = load_and_prepare_gene_sets(args.input, args.pval_threshold, log_path)
    
    # ---------------------------------------------------------------------
    # Run pathway analysis for each set
    # ---------------------------------------------------------------------
    all_results = {}
    
    for set_name in ['upregulated', 'downregulated']:
        log_message(log_path, '')
        gene_list = gene_sets[set_name]
        
        if len(gene_list) < 5:
            log_message(log_path, f'  Skipping {set_name} (only {len(gene_list)} genes, need at least 5)')
            all_results[set_name] = {}
            continue
        
        # Run enrichment
        results = run_pathway_analysis(
            gene_list=gene_list,
            background=gene_sets['background'],
            set_name=set_name,
            output_dir=args.output_dir,
            log_path=log_path,
            fdr_threshold=args.fdr_threshold,
            write_clean=write_clean
        )

        all_results[set_name] = results

        # In enrich-only mode, stop after writing CSVs: skip the legacy
        # in-script figures, the per-set summary, and (below) the final
        # report so a targeted source run does not overwrite shared outputs.
        if enrich_only:
            continue

        # Create visualizations
        log_message(log_path, f'  Creating visualizations...')
        for source in SOURCES:
            result_data = results.get(source, {'clean': pd.DataFrame(), 'full': pd.DataFrame()})
            df_clean = result_data.get('clean', pd.DataFrame())
            df_full = result_data.get('full', pd.DataFrame())
            
            if len(df_clean) > 0:
                # Basic visualizations
                create_enrichment_bar_plot(
                    df_clean, source, set_name, args.output_dir, log_path
                )
                create_dot_plot(
                    df_clean, source, set_name, args.output_dir, log_path
                )
                
                # Results table
                create_results_table(
                    df_clean, source, set_name, args.output_dir, log_path
                )
                
                # Gene-term network
                create_gene_term_network(
                    df_clean, source, set_name, args.output_dir, log_path, df_full
                )
                
                # Term clustering (only for GO:BP where we expect redundancy)
                if source == 'GO:BP':
                    create_term_clustering(
                        df_clean, source, set_name, args.output_dir, log_path, df_full
                    )
                
                # Publication-style figure (only for upregulated since downregulated had no results)
                if set_name == 'upregulated':
                    create_publication_figure(
                        df_clean, source, set_name, args.output_dir, log_path,
                        df_full, gene_sets['protein_data']
                    )
        
        # Generate set summary
        generate_set_summary(
            results, gene_list, set_name, args.output_dir, args.fdr_threshold, log_path
        )
    
    # ---------------------------------------------------------------------
    # Generate final report (skipped in enrich-only mode to preserve the
    # existing GO:BP/KEGG report)
    # ---------------------------------------------------------------------
    if not enrich_only:
        log_message(log_path, '')
        generate_final_report(
            gene_sets, all_results, args.output_dir,
            args.pval_threshold, args.fdr_threshold, log_path
        )
    
    # ---------------------------------------------------------------------
    # Finalize
    # ---------------------------------------------------------------------
    log_message(log_path, '')
    log_message(log_path, '=' * 50)
    log_message(log_path, 'Pathway analysis complete!')
    log_message(log_path, '=' * 50)


if __name__ == '__main__':
    main()
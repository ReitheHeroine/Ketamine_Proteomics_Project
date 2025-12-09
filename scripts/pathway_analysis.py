# title: pathway_analysis.py
# project: Ketamine Proteomics Analysis Project
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2025-12-09
# last modified: 2025-12-09
#
# purpose:
#   Perform pathway enrichment analysis on differential abundance results using
#   over-representation analysis (ORA). Queries GO Biological Process and KEGG
#   databases via gProfiler. Analyzes upregulated and downregulated protein sets
#   separately, including presence/absence proteins.
#
# inputs:
#   - all_proteins_categorized.csv: Master file from diff_abundance_analysis.py
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
#   │   ├── upregulated_KEGG_plot.pdf/html
#   │   ├── upregulated_KEGG_dotplot.pdf/html
#   │   ├── upregulated_KEGG_table.pdf/html
#   │   ├── upregulated_KEGG_network.pdf/html
#   │   └── upregulated_summary.txt
#   ├── downregulated/
#   │   └── (similar structure, if significant terms found)
#   └── pathway_analysis_report.txt
#
# usage:
#   python pathway_analysis.py \
#       --input ../results/combined/all_proteins_categorized.csv \
#       --output_dir ../results/pathway_analysis \
#       --log_dir ../logs \
#       --pval_threshold 0.05 \
#       --fdr_threshold 0.05
#
# dependencies:
#   pip install gprofiler-official pandas numpy plotly kaleido networkx scipy
#
# notes:
#   - Requires internet connection to query gProfiler API
#   - Background set: all 903 detected proteins
#   - Upregulated set: significant up (p ≤ 0.05) + ketamine-specific
#   - Downregulated set: significant down (p ≤ 0.05) + control-specific

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

# Column names from input data
GENE_COL = 'Gene Symbol'
CATEGORY_COL = 'category'
DIRECTION_COL = 'direction'
PVAL_COL = 'Abundance Ratio Adj. P-Value: (ketamine) / (control)'

# gProfiler settings
ORGANISM = 'mmusculus'  # Mus musculus (mouse)
SOURCES = ['GO:BP', 'KEGG']  # GO Biological Process and KEGG

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
    
    gp = GProfiler(return_dataframe=True)
    
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


def run_pathway_analysis(gene_list, background, set_name, output_dir, log_path, fdr_threshold):
    '''
    Run full pathway analysis for a gene set (GO:BP and KEGG).
    
    Parameters:
        gene_list (list): List of gene symbols
        background (list): Background gene list
        set_name (str): Name of the set ('upregulated' or 'downregulated')
        output_dir (str): Directory for outputs
        log_path (str): Path to log file
        fdr_threshold (float): FDR threshold for significance
    
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
            
            # Save to CSV
            source_name = source.replace(':', '_')
            csv_path = os.path.join(set_output_dir, f'{set_name}_{source_name}.csv')
            df_clean.to_csv(csv_path, index=False)
            log_message(log_path, f'      Saved: {csv_path}')
        else:
            results[source] = {'clean': pd.DataFrame(), 'full': full_results}
    
    return results


# =============================================================================
# VISUALIZATION FUNCTIONS
# =============================================================================

def create_enrichment_bar_plot(df, source, set_name, output_dir, log_path, n_terms=15):
    '''
    Create horizontal bar plot of top enriched terms.
    
    Parameters:
        df (pd.DataFrame): Enrichment results
        source (str): Database source name
        set_name (str): 'upregulated' or 'downregulated'
        output_dir (str): Output directory
        log_path (str): Path to log file
        n_terms (int): Number of top terms to display
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
    source_display = 'GO Biological Process' if source == 'GO:BP' else 'KEGG Pathways'
    direction_display = 'Upregulated' if set_name == 'upregulated' else 'Downregulated'
    
    # Layout
    fig.update_layout(
        title=dict(
            text=f'{source_display}: {direction_display} Proteins',
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
    source_name = source.replace(':', '_')
    set_output_dir = os.path.join(output_dir, set_name)
    
    pdf_path = os.path.join(set_output_dir, f'{set_name}_{source_name}_plot.pdf')
    html_path = os.path.join(set_output_dir, f'{set_name}_{source_name}_plot.html')
    
    fig.write_image(pdf_path, scale=2)
    fig.write_html(html_path)
    
    log_message(log_path, f'    Saved: {pdf_path}')
    log_message(log_path, f'    Saved: {html_path}')


def create_dot_plot(df, source, set_name, output_dir, log_path, n_terms=15):
    '''
    Create dot plot showing gene ratio, gene count, and p-value.
    
    Parameters:
        df (pd.DataFrame): Enrichment results
        source (str): Database source name
        set_name (str): 'upregulated' or 'downregulated'
        output_dir (str): Output directory
        log_path (str): Path to log file
        n_terms (int): Number of top terms to display
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
    source_display = 'GO Biological Process' if source == 'GO:BP' else 'KEGG Pathways'
    direction_display = 'Upregulated' if set_name == 'upregulated' else 'Downregulated'
    
    # Layout
    fig.update_layout(
        title=dict(
            text=f'{source_display}: {direction_display} Proteins (Dot Plot)',
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
    source_name = source.replace(':', '_')
    set_output_dir = os.path.join(output_dir, set_name)
    
    pdf_path = os.path.join(set_output_dir, f'{set_name}_{source_name}_dotplot.pdf')
    html_path = os.path.join(set_output_dir, f'{set_name}_{source_name}_dotplot.html')
    
    fig.write_image(pdf_path, scale=2)
    fig.write_html(html_path)
    
    log_message(log_path, f'    Saved: {pdf_path}')
    log_message(log_path, f'    Saved: {html_path}')


def create_results_table(df, source, set_name, output_dir, log_path):
    '''
    Create formatted HTML and PDF tables of enrichment results.
    
    Parameters:
        df (pd.DataFrame): Enrichment results
        source (str): Database source name
        set_name (str): 'upregulated' or 'downregulated'
        output_dir (str): Output directory
        log_path (str): Path to log file
    '''
    if len(df) == 0:
        return
    
    # Format source name for display
    source_display = 'GO Biological Process' if source == 'GO:BP' else 'KEGG Pathways'
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
    source_name = source.replace(':', '_')
    set_output_dir = os.path.join(output_dir, set_name)
    html_path = os.path.join(set_output_dir, f'{set_name}_{source_name}_table.html')
    
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
    
    fig.update_layout(
        title=dict(
            text=f'{source_display}: {direction_display} Proteins<br><sup>Showing {len(table_df)} enriched terms</sup>',
            font=dict(size=14, family='Arial Black')
        ),
        width=900,
        height=max(400, 80 + len(table_df) * 24),
        margin=dict(t=80, b=20, l=20, r=20)
    )
    
    pdf_path = os.path.join(set_output_dir, f'{set_name}_{source_name}_table.pdf')
    fig.write_image(pdf_path, scale=2)
    
    log_message(log_path, f'    Saved: {pdf_path}')


def create_gene_term_network(df, source, set_name, output_dir, log_path, enrichment_results_full):
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
    source_display = 'GO Biological Process' if source == 'GO:BP' else 'KEGG Pathways'
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
    source_name = source.replace(':', '_')
    set_output_dir = os.path.join(output_dir, set_name)
    
    pdf_path = os.path.join(set_output_dir, f'{set_name}_{source_name}_network.pdf')
    html_path = os.path.join(set_output_dir, f'{set_name}_{source_name}_network.html')
    
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
    source_display = 'GO Biological Process' if source == 'GO:BP' else 'KEGG Pathways'
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
    source_name = source.replace(':', '_')
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


def create_publication_figure(df, source, set_name, output_dir, log_path, 
                               enrichment_results_full, protein_data, n_top_terms=6):
    '''
    Create publication-style figure inspired by Herzog et al. 2021 Figure 4.
    
    Panel A: Top biological functions with hierarchical clustering and significance heatmap
    Panels B-F: Gene-level heatmaps for each top pathway showing log2 fold changes
    
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
    # Calculate similarity matrix for hierarchical clustering
    # ---------------------------------------------------------------------
    term_list = list(term_genes.keys())
    n_terms = len(term_list)
    
    similarity_matrix = np.zeros((n_terms, n_terms))
    for i in range(n_terms):
        for j in range(n_terms):
            set_i = set(term_genes[term_list[i]])
            set_j = set(term_genes[term_list[j]])
            intersection = len(set_i & set_j)
            union = len(set_i | set_j)
            similarity_matrix[i, j] = intersection / union if union > 0 else 0
    
    # Perform hierarchical clustering
    distance_matrix = 1 - similarity_matrix
    np.fill_diagonal(distance_matrix, 0)
    
    condensed_dist = pdist(distance_matrix)
    if len(condensed_dist) > 0:
        linkage_matrix = linkage(condensed_dist, method='average')
        from scipy.cluster.hierarchy import dendrogram
        dendro = dendrogram(linkage_matrix, no_plot=True)
        order = dendro['leaves']
    else:
        order = list(range(n_terms))
    
    # Reorder terms based on clustering
    ordered_term_ids = [term_list[i] for i in order]
    ordered_term_names = []
    ordered_pvalues = []
    
    for tid in ordered_term_ids:
        row = top_terms[top_terms['term_id'] == tid]
        if len(row) > 0:
            ordered_term_names.append(row['term_name'].values[0])
            ordered_pvalues.append(row['fdr_pvalue'].values[0])
    
    # ---------------------------------------------------------------------
    # Create figure with subplots
    # Panel A: Main heatmap with clustering
    # Panels B onwards: Individual pathway gene heatmaps
    # ---------------------------------------------------------------------
    n_pathway_panels = min(5, len(ordered_term_ids))  # Max 5 pathway panels
    
    # Calculate subplot layout
    fig = make_subplots(
        rows=2, cols=max(3, n_pathway_panels),
        row_heights=[0.5, 0.5],
        column_widths=[0.4] + [0.6 / max(2, n_pathway_panels - 1)] * (max(2, n_pathway_panels - 1)),
        specs=[[{'colspan': max(3, n_pathway_panels), 'type': 'heatmap'}, None, None] + [None] * (max(0, n_pathway_panels - 3)),
               [{'type': 'heatmap'}] * max(3, n_pathway_panels)],
        subplot_titles=['A. Top Biological Functions'] + 
                       [f'{chr(66 + i)}. {ordered_term_names[i][:30]}...' if len(ordered_term_names[i]) > 30 
                        else f'{chr(66 + i)}. {ordered_term_names[i]}' 
                        for i in range(n_pathway_panels)],
        vertical_spacing=0.15,
        horizontal_spacing=0.05
    )
    
    # ---------------------------------------------------------------------
    # Panel A: Top biological functions heatmap
    # ---------------------------------------------------------------------
    # Create heatmap data
    heatmap_values = [[-np.log10(p)] for p in ordered_pvalues]
    
    fig.add_trace(
        go.Heatmap(
            z=heatmap_values,
            y=ordered_term_names,
            x=['Ketamine vs Control'],
            colorscale='Purples',
            showscale=True,
            colorbar=dict(
                title='-log₁₀(FDR)',
                x=0.35,
                len=0.4,
                y=0.75
            ),
            hovertemplate='<b>%{y}</b><br>-log₁₀(FDR): %{z:.2f}<extra></extra>'
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
        
        # Add heatmap for this pathway
        fig.add_trace(
            go.Heatmap(
                z=gene_values,
                y=gene_names,
                x=['log₂FC'],
                colorscale=[[0, '#2CA02C'], [0.5, '#FFFFFF'], [1, '#D62728']],  # Green-White-Red
                zmid=0,
                zmin=-4,
                zmax=4,
                showscale=(panel_idx == n_pathway_panels - 1),  # Only show scale on last panel
                colorbar=dict(
                    title='log₂FC',
                    x=1.02,
                    len=0.4,
                    y=0.25
                ) if panel_idx == n_pathway_panels - 1 else None,
                hovertemplate='<b>%{y}</b><br>log₂FC: %{z:.2f}<extra></extra>'
            ),
            row=2, col=panel_idx + 1
        )
    
    # ---------------------------------------------------------------------
    # Update layout
    # ---------------------------------------------------------------------
    source_display = 'GO Biological Process' if source == 'GO:BP' else 'KEGG Pathways'
    
    fig.update_layout(
        title=dict(
            text=f'Pathway Enrichment Analysis: {source_display}<br>' +
                 f'<sup>Upregulated proteins in Ketamine vs Control</sup>',
            font=dict(size=16, family='Arial Black'),
            x=0.5
        ),
        height=900,
        width=1400,
        plot_bgcolor='white',
        paper_bgcolor='white',
        font=dict(size=10)
    )
    
    # Update y-axes for gene heatmaps
    for i in range(n_pathway_panels):
        fig.update_yaxes(tickfont=dict(size=9, color='#1F77B4'), row=2, col=i + 1)
    
    # ---------------------------------------------------------------------
    # Save outputs
    # ---------------------------------------------------------------------
    source_name = source.replace(':', '_')
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
            source_display = 'GO Biological Process' if source == 'GO:BP' else 'KEGG Pathways'
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
        # Top findings
        # ---------------------------------------------------------------------
        f.write('-' * 70 + '\n')
        f.write('TOP FINDINGS\n')
        f.write('-' * 70 + '\n\n')
        
        for set_name in ['upregulated', 'downregulated']:
            f.write(f'{set_name.upper()} - Top 5 GO Biological Process terms:\n')
            results = all_results.get(set_name, {})
            result_data = results.get('GO:BP', {'clean': pd.DataFrame()})
            df = result_data.get('clean', pd.DataFrame()) if isinstance(result_data, dict) else result_data
            
            if len(df) > 0:
                for i, (_, row) in enumerate(df.head(5).iterrows(), 1):
                    f.write(f"  {i}. {row['term_name']} (FDR={row['fdr_pvalue']:.2e})\n")
            else:
                f.write('  No significant terms found.\n')
            f.write('\n')
    
    log_message(log_path, f'Saved final report: {report_path}')


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
        required=True,
        help='Path to all_proteins_categorized.csv'
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
    
    args = parser.parse_args()
    
    # ---------------------------------------------------------------------
    # Setup
    # ---------------------------------------------------------------------
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(args.log_dir, exist_ok=True)
    log_path = setup_logging(args.log_dir)
    
    log_message(log_path, '=' * 50)
    log_message(log_path, 'KETAMINE PROTEOMICS ANALYSIS PROJECT')
    log_message(log_path, 'Pathway Enrichment Analysis')
    log_message(log_path, '=' * 50)
    log_message(log_path, f'Input file:        {args.input}')
    log_message(log_path, f'Output directory:  {args.output_dir}')
    log_message(log_path, f'P-value threshold: {args.pval_threshold}')
    log_message(log_path, f'FDR threshold:     {args.fdr_threshold}')
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
            fdr_threshold=args.fdr_threshold
        )
        
        all_results[set_name] = results
        
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
    # Generate final report
    # ---------------------------------------------------------------------
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
# title: visualize_diff_abundance.py
# project: Ketamine Proteomics Analysis Project
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2025-12-09
# last modified: 2025-12-09
#
# purpose:
#   Generate publication-quality visualizations for differential abundance
#   proteomics results. Creates volcano plots, summary bar charts, and top
#   proteins tables/charts. Outputs both static (PDF) and interactive (HTML).
#
# inputs:
#   - all_proteins_categorized.csv: Master file from diff_abundance_analysis.py
#
# outputs:
#   results/figures/
#   ├── volcano_plot.pdf
#   ├── volcano_plot.html
#   ├── ma_plot.pdf
#   ├── ma_plot.html
#   ├── summary_bar_chart.pdf
#   ├── top_proteins_bar_chart.pdf
#   ├── top_proteins_bar_chart.html
#   ├── top20_upregulated_table.pdf
#   ├── top20_upregulated_table.html
#   ├── top20_downregulated_table.pdf
#   └── top20_downregulated_table.html
#
# usage:
#   python visualize_diff_abundance.py \
#       --input ../results/combined/all_proteins_categorized.csv \
#       --output_dir ../results/figures \
#       --log_dir ../logs \
#       --pval_threshold 0.05
#
#   copy/paste: python visualize_diff_abundance.py --input ../results/combined/all_proteins_categorized.csv --output_dir ../results/figures --log_dir ../logs --pval_threshold 0.05
#
# notes:
#   - Volcano plot shows quantitative proteins only
#   - Interactive HTML files can be opened in any web browser

import pandas as pd
import numpy as np
import plotly.graph_objects as go
import argparse
import os
from datetime import datetime

# =============================================================================
# CONFIGURATION AND CONSTANTS
# =============================================================================

# Column names
RATIO_COL = 'Abundance Ratio: (ketamine) / (control)'
PVAL_COL = 'Abundance Ratio Adj. P-Value: (ketamine) / (control)'
LOG2FC_COL = 'log2_fold_change'
ACCESSION_COL = 'Accession'
GENE_COL = 'Gene Symbol'
CATEGORY_COL = 'category'
DIRECTION_COL = 'direction'
SIGNIFICANT_COL = 'significant'
ABUNDANCE_CONTROL_COL = 'Abundances (Grouped): control'
ABUNDANCE_KETAMINE_COL = 'Abundances (Grouped): ketamine'
SOURCE_FILE_COL = 'source_file'

# Visual constants
COLOR_UP = '#D62728'       # Red for upregulated
COLOR_DOWN = '#1F77B4'     # Blue for downregulated
COLOR_NS = '#7F7F7F'       # Gray for not significant
COLOR_PA_KET = '#FF7F0E'   # Orange for ketamine-specific
COLOR_PA_CTRL = '#9467BD'  # Purple for control-specific

# =============================================================================
# LOGGING FUNCTIONS
# =============================================================================

def setup_logging(log_dir):
    '''Initialize timestamped log file.'''
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_filename = f'visualize_diff_abundance_{timestamp}.log'
    log_path = os.path.join(log_dir, log_filename)
    
    with open(log_path, 'w') as f:
        f.write('=' * 70 + '\n')
        f.write('KETAMINE PROTEOMICS ANALYSIS PROJECT - VISUALIZATION LOG\n')
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

def load_and_prepare_data(input_path, pval_threshold, log_path):
    '''
    Load the categorized proteins file and prepare for visualization.
    '''
    log_message(log_path, f'Loading data from {input_path}')
    
    df = pd.read_csv(input_path)
    log_message(log_path, f'  Loaded {len(df)} proteins')
    
    # ---------------------------------------------------------------------
    # Calculate -log10(p-value) for volcano plot
    # ---------------------------------------------------------------------
    df['neg_log10_pval'] = -np.log10(df[PVAL_COL].replace(0, 1e-20))
    
    # ---------------------------------------------------------------------
    # Re-apply significance based on current threshold
    # ---------------------------------------------------------------------
    df['significant'] = (df[PVAL_COL] <= pval_threshold) & (df[CATEGORY_COL] == 'quantitative')
    
    # ---------------------------------------------------------------------
    # Assign colors and legend groups for plotting
    # ---------------------------------------------------------------------
    def assign_plot_properties(row):
        if row[CATEGORY_COL] == 'presence_absence_ketamine_specific':
            return COLOR_PA_KET, 'Ketamine-specific'
        elif row[CATEGORY_COL] == 'presence_absence_control_specific':
            return COLOR_PA_CTRL, 'Control-specific'
        elif row[CATEGORY_COL] == 'quantitative':
            if row['significant']:
                if row[DIRECTION_COL] == 'up_in_ketamine':
                    return COLOR_UP, 'Up in Ketamine (sig.)'
                else:
                    return COLOR_DOWN, 'Down in Ketamine (sig.)'
            else:
                return COLOR_NS, 'Not significant'
        else:
            return COLOR_NS, 'Other'
    
    df[['color', 'legend_group']] = df.apply(
        assign_plot_properties, axis=1, result_type='expand'
    )
    
    log_message(log_path, '  Data preparation complete')
    
    return df


# =============================================================================
# VISUALIZATION FUNCTIONS
# =============================================================================

def create_volcano_plot(df, pval_threshold, output_dir, log_path):
    '''
    Create volcano plot for quantitative proteins only.
    Labels top 10 proteins by p-value.
    '''
    log_message(log_path, 'Creating volcano plot...')
    
    # ---------------------------------------------------------------------
    # Filter to quantitative proteins only
    # ---------------------------------------------------------------------
    quant_df = df[df[CATEGORY_COL] == 'quantitative'].copy()
    
    fig = go.Figure()
    
    # ---------------------------------------------------------------------
    # Plot each group in order
    # ---------------------------------------------------------------------
    legend_order = [
        ('Up in Ketamine (sig.)', COLOR_UP),
        ('Down in Ketamine (sig.)', COLOR_DOWN),
        ('Not significant', COLOR_NS),
    ]
    
    for legend_name, color in legend_order:
        subset = quant_df[quant_df['legend_group'] == legend_name]
        if len(subset) == 0:
            continue
            
        fig.add_trace(go.Scatter(
            x=subset[LOG2FC_COL],
            y=subset['neg_log10_pval'],
            mode='markers',
            marker=dict(
                color=color,
                size=8,
                line=dict(width=0.5, color='white')
            ),
            name=f'{legend_name} (n={len(subset)})',
            text=subset.apply(
                lambda r: f"Gene: {r[GENE_COL]}<br>"
                          f"Accession: {r[ACCESSION_COL]}<br>"
                          f"log2FC: {r[LOG2FC_COL]:.2f}<br>"
                          f"p-value: {r[PVAL_COL]:.2e}",
                axis=1
            ),
            hoverinfo='text'
        ))
    
    # ---------------------------------------------------------------------
    # Add significance threshold line
    # ---------------------------------------------------------------------
    fig.add_hline(
        y=-np.log10(pval_threshold),
        line_dash='dash',
        line_color='black',
        line_width=1,
        annotation_text=f'p = {pval_threshold}',
        annotation_position='top right'
    )
    
    # ---------------------------------------------------------------------
    # Label top 10 proteins by p-value
    # ---------------------------------------------------------------------
    top10 = quant_df.nsmallest(10, PVAL_COL)
    
    for _, row in top10.iterrows():
        fig.add_annotation(
            x=row[LOG2FC_COL],
            y=row['neg_log10_pval'],
            text=row[GENE_COL],
            showarrow=True,
            arrowhead=0,
            arrowsize=0.5,
            arrowwidth=1,
            ax=20,
            ay=-20,
            font=dict(size=9, color='black'),
            bgcolor='rgba(255,255,255,0.7)',
            borderpad=2
        )
    
    # ---------------------------------------------------------------------
    # Layout
    # ---------------------------------------------------------------------
    fig.update_layout(
        title=dict(
            text='Volcano Plot: Ketamine vs Control',
            font=dict(size=18, family='Arial Black')
        ),
        xaxis=dict(
            title=dict(text='log₂(Fold Change)', font=dict(size=14)),
            zeroline=True,
            zerolinecolor='lightgray',
            zerolinewidth=1,
            gridcolor='rgba(0,0,0,0.1)'
        ),
        yaxis=dict(
            title=dict(text='-log₁₀(Adjusted p-value)', font=dict(size=14)),
            gridcolor='rgba(0,0,0,0.1)'
        ),
        legend=dict(
            title=dict(text='Category'),
            yanchor='top',
            y=0.99,
            xanchor='left',
            x=1.02,
            bgcolor='rgba(255,255,255,0.8)',
            bordercolor='lightgray',
            borderwidth=1
        ),
        plot_bgcolor='white',
        width=900,
        height=650,
        margin=dict(r=200)
    )
    
    # ---------------------------------------------------------------------
    # Save outputs
    # ---------------------------------------------------------------------
    html_path = os.path.join(output_dir, 'volcano_plot.html')
    pdf_path = os.path.join(output_dir, 'volcano_plot.pdf')
    
    fig.write_html(html_path)
    fig.write_image(pdf_path, scale=2)
    
    log_message(log_path, f'  Saved: {html_path}')
    log_message(log_path, f'  Saved: {pdf_path}')


def create_ma_plot(df, pval_threshold, output_dir, log_path):
    '''
    Create MA plot (Mean-Average plot) for quantitative proteins.
    X-axis: Average log2 abundance (A)
    Y-axis: log2 Fold Change (M)
    
    This plot helps identify abundance-dependent bias in fold change estimates.
    '''
    log_message(log_path, 'Creating MA plot...')
    
    # ---------------------------------------------------------------------
    # Filter to quantitative proteins with valid abundance values
    # ---------------------------------------------------------------------
    quant_df = df[df[CATEGORY_COL] == 'quantitative'].copy()
    
    # Remove rows with missing or zero abundance values
    quant_df = quant_df[
        (quant_df[ABUNDANCE_CONTROL_COL] > 0) & 
        (quant_df[ABUNDANCE_KETAMINE_COL] > 0)
    ].copy()
    
    # ---------------------------------------------------------------------
    # Calculate average log2 abundance (A value)
    # A = (log2(control) + log2(ketamine)) / 2 = log2(sqrt(control * ketamine))
    # ---------------------------------------------------------------------
    quant_df['log2_avg_abundance'] = (
        np.log2(quant_df[ABUNDANCE_CONTROL_COL]) + 
        np.log2(quant_df[ABUNDANCE_KETAMINE_COL])
    ) / 2
    
    fig = go.Figure()
    
    # ---------------------------------------------------------------------
    # Plot each group in order
    # ---------------------------------------------------------------------
    legend_order = [
        ('Up in Ketamine (sig.)', COLOR_UP),
        ('Down in Ketamine (sig.)', COLOR_DOWN),
        ('Not significant', COLOR_NS),
    ]
    
    for legend_name, color in legend_order:
        subset = quant_df[quant_df['legend_group'] == legend_name]
        if len(subset) == 0:
            continue
            
        fig.add_trace(go.Scatter(
            x=subset['log2_avg_abundance'],
            y=subset[LOG2FC_COL],
            mode='markers',
            marker=dict(
                color=color,
                size=8,
                line=dict(width=0.5, color='white')
            ),
            name=f'{legend_name} (n={len(subset)})',
            text=subset.apply(
                lambda r: f"Gene: {r[GENE_COL]}<br>"
                          f"Accession: {r[ACCESSION_COL]}<br>"
                          f"log2FC: {r[LOG2FC_COL]:.2f}<br>"
                          f"Avg log2 Abundance: {r['log2_avg_abundance']:.2f}<br>"
                          f"p-value: {r[PVAL_COL]:.2e}",
                axis=1
            ),
            hoverinfo='text'
        ))
    
    # ---------------------------------------------------------------------
    # Add horizontal line at y=0 (no change)
    # ---------------------------------------------------------------------
    fig.add_hline(
        y=0,
        line_dash='dash',
        line_color='black',
        line_width=1
    )
    
    # ---------------------------------------------------------------------
    # Layout
    # ---------------------------------------------------------------------
    fig.update_layout(
        title=dict(
            text='MA Plot: Ketamine vs Control',
            font=dict(size=18, family='Arial Black')
        ),
        xaxis=dict(
            title=dict(text='Average log₂(Abundance)', font=dict(size=14)),
            gridcolor='rgba(0,0,0,0.1)'
        ),
        yaxis=dict(
            title=dict(text='log₂(Fold Change)', font=dict(size=14)),
            zeroline=True,
            zerolinecolor='lightgray',
            zerolinewidth=1,
            gridcolor='rgba(0,0,0,0.1)'
        ),
        legend=dict(
            title=dict(text='Category'),
            yanchor='top',
            y=0.99,
            xanchor='left',
            x=1.02,
            bgcolor='rgba(255,255,255,0.8)',
            bordercolor='lightgray',
            borderwidth=1
        ),
        plot_bgcolor='white',
        width=900,
        height=650,
        margin=dict(r=200)
    )
    
    # ---------------------------------------------------------------------
    # Add annotation explaining the plot
    # ---------------------------------------------------------------------
    fig.add_annotation(
        text='Dashed line: no change (log₂FC = 0)',
        xref='paper', yref='paper',
        x=0, y=-0.1,
        showarrow=False,
        font=dict(size=10, color='gray'),
        align='left'
    )
    
    # ---------------------------------------------------------------------
    # Save outputs
    # ---------------------------------------------------------------------
    html_path = os.path.join(output_dir, 'ma_plot.html')
    pdf_path = os.path.join(output_dir, 'ma_plot.pdf')
    
    fig.write_html(html_path)
    fig.write_image(pdf_path, scale=2)
    
    log_message(log_path, f'  Saved: {html_path}')
    log_message(log_path, f'  Saved: {pdf_path}')


def create_summary_bar_chart(df, pval_threshold, output_dir, log_path):
    '''
    Create summary bar chart showing counts of differential proteins.
    Fixed formatting to prevent label cutoff and overlap.
    '''
    log_message(log_path, 'Creating summary bar chart...')
    
    # ---------------------------------------------------------------------
    # Calculate counts
    # ---------------------------------------------------------------------
    quant_df = df[df[CATEGORY_COL] == 'quantitative']
    sig_up = len(quant_df[(quant_df['significant']) & (quant_df[DIRECTION_COL] == 'up_in_ketamine')])
    sig_down = len(quant_df[(quant_df['significant']) & (quant_df[DIRECTION_COL] == 'down_in_ketamine')])
    
    ket_specific = len(df[df[CATEGORY_COL] == 'presence_absence_ketamine_specific'])
    ctrl_specific = len(df[df[CATEGORY_COL] == 'presence_absence_control_specific'])
    
    # ---------------------------------------------------------------------
    # Create figure
    # ---------------------------------------------------------------------
    categories = [
        'Upregulated<br>(Quantitative)',
        'Downregulated<br>(Quantitative)',
        'Ketamine-<br>specific',
        'Control-<br>specific'
    ]
    counts = [sig_up, sig_down, ket_specific, ctrl_specific]
    colors = [COLOR_UP, COLOR_DOWN, COLOR_PA_KET, COLOR_PA_CTRL]
    
    fig = go.Figure()
    
    fig.add_trace(go.Bar(
        x=categories,
        y=counts,
        marker_color=colors,
        text=counts,
        textposition='outside',
        textfont=dict(size=14, color='black')
    ))
    
    # ---------------------------------------------------------------------
    # Layout - fixed margins and y-axis range
    # ---------------------------------------------------------------------
    max_count = max(counts)
    
    fig.update_layout(
        title=dict(
            text=f'Summary of Differential Proteins (p ≤ {pval_threshold})',
            font=dict(size=16, family='Arial Black'),
            y=0.95
        ),
        xaxis=dict(
            title=dict(text='', font=dict(size=12)),
            tickfont=dict(size=11)
        ),
        yaxis=dict(
            title=dict(text='Number of Proteins', font=dict(size=14)),
            gridcolor='rgba(0,0,0,0.1)',
            range=[0, max_count * 1.2]  # Add 20% headroom for labels
        ),
        plot_bgcolor='white',
        width=650,
        height=550,
        showlegend=False,
        bargap=0.3,
        margin=dict(t=80, b=120, l=60, r=40)  # Increased margins
    )
    
    # ---------------------------------------------------------------------
    # Add footnote with proper positioning
    # ---------------------------------------------------------------------
    fig.add_annotation(
        text=f'Quantitative: adj. p ≤ {pval_threshold}<br>Presence/Absence: detected in one condition only',
        xref='paper', yref='paper',
        x=0.5, y=-0.22,
        showarrow=False,
        font=dict(size=10, color='gray'),
        align='center'
    )
    
    # ---------------------------------------------------------------------
    # Save output
    # ---------------------------------------------------------------------
    pdf_path = os.path.join(output_dir, 'summary_bar_chart.pdf')
    fig.write_image(pdf_path, scale=2)
    
    log_message(log_path, f'  Saved: {pdf_path}')


def create_top_proteins_bar_chart(df, output_dir, log_path, n_proteins=10):
    '''
    Create horizontal bar chart of top proteins by fold change.
    '''
    log_message(log_path, f'Creating top {n_proteins} proteins bar chart...')
    
    # ---------------------------------------------------------------------
    # Get top proteins (significant only, by absolute fold change)
    # ---------------------------------------------------------------------
    sig_df = df[(df[CATEGORY_COL] == 'quantitative') & (df['significant'])].copy()
    sig_df['abs_log2fc'] = sig_df[LOG2FC_COL].abs()
    top_df = sig_df.nlargest(n_proteins, 'abs_log2fc').sort_values(LOG2FC_COL)
    
    # ---------------------------------------------------------------------
    # Create figure
    # ---------------------------------------------------------------------
    colors = [COLOR_UP if x > 0 else COLOR_DOWN for x in top_df[LOG2FC_COL]]
    
    fig = go.Figure()
    
    fig.add_trace(go.Bar(
        y=top_df[GENE_COL],
        x=top_df[LOG2FC_COL],
        orientation='h',
        marker_color=colors,
        text=[f'{x:.2f}' for x in top_df[LOG2FC_COL]],
        textposition='outside',
        textfont=dict(size=10)
    ))
    
    # ---------------------------------------------------------------------
    # Layout
    # ---------------------------------------------------------------------
    fig.update_layout(
        title=dict(
            text=f'Top {n_proteins} Significant Proteins by Fold Change',
            font=dict(size=16, family='Arial Black')
        ),
        xaxis=dict(
            title=dict(text='log₂(Fold Change)', font=dict(size=14)),
            zeroline=True,
            zerolinecolor='black',
            zerolinewidth=1,
            gridcolor='rgba(0,0,0,0.1)'
        ),
        yaxis=dict(
            title=dict(text=''),
            tickfont=dict(size=11)
        ),
        plot_bgcolor='white',
        width=700,
        height=500,
        margin=dict(l=120, r=80)
    )
    
    # ---------------------------------------------------------------------
    # Save outputs
    # ---------------------------------------------------------------------
    html_path = os.path.join(output_dir, 'top_proteins_bar_chart.html')
    pdf_path = os.path.join(output_dir, 'top_proteins_bar_chart.pdf')
    
    fig.write_html(html_path)
    fig.write_image(pdf_path, scale=2)
    
    log_message(log_path, f'  Saved: {html_path}')
    log_message(log_path, f'  Saved: {pdf_path}')


def create_top_proteins_table(df, output_dir, log_path, direction, n_proteins=20):
    '''
    Create table of top proteins (HTML and PDF).
    
    Parameters:
        direction (str): 'up' for upregulated, 'down' for downregulated
    '''
    # ---------------------------------------------------------------------
    # Determine direction-specific parameters
    # ---------------------------------------------------------------------
    if direction == 'up':
        direction_filter = 'up_in_ketamine'
        title_text = f'Top {n_proteins} Upregulated Proteins (Ketamine vs Control)'
        file_prefix = 'top20_upregulated'
        direction_symbol = '↑ Up'
        header_color = COLOR_UP
    else:
        direction_filter = 'down_in_ketamine'
        title_text = f'Top {n_proteins} Downregulated Proteins (Ketamine vs Control)'
        file_prefix = 'top20_downregulated'
        direction_symbol = '↓ Down'
        header_color = COLOR_DOWN
    
    log_message(log_path, f'Creating {direction}regulated proteins table...')
    
    # ---------------------------------------------------------------------
    # Get proteins by p-value for the specified direction
    # ---------------------------------------------------------------------
    sig_df = df[
        (df[CATEGORY_COL] == 'quantitative') & 
        (df['significant']) &
        (df[DIRECTION_COL] == direction_filter)
    ].copy()
    
    # Get top N (or all if fewer than N)
    actual_n = min(n_proteins, len(sig_df))
    top_df = sig_df.nsmallest(actual_n, PVAL_COL)
    
    if len(top_df) == 0:
        log_message(log_path, f'  No {direction}regulated proteins found, skipping table')
        return
    
    # ---------------------------------------------------------------------
    # Prepare table data
    # ---------------------------------------------------------------------
    table_df = top_df[[ACCESSION_COL, GENE_COL, LOG2FC_COL, PVAL_COL, SOURCE_FILE_COL]].copy()
    table_df.columns = ['Accession', 'Gene', 'log2FC', 'Adj. p-value', 'Source']
    table_df['log2FC'] = table_df['log2FC'].round(3)
    table_df['Adj. p-value'] = table_df['Adj. p-value'].apply(lambda x: f'{x:.2e}')
    table_df = table_df.reset_index(drop=True)
    table_df.index = table_df.index + 1
    table_df.index.name = 'Rank'
    table_df = table_df.reset_index()
    
    # ---------------------------------------------------------------------
    # Create HTML table
    # ---------------------------------------------------------------------
    html_content = f'''<!DOCTYPE html>
<html>
<head>
    <title>{title_text}</title>
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
            max-width: 800px;
            background-color: white;
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
        }}
        th {{
            background-color: {header_color};
            color: white;
            padding: 12px 15px;
            text-align: left;
            font-weight: bold;
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
        .gene {{
            font-weight: bold;
        }}
        .footer {{
            margin-top: 20px;
            font-size: 12px;
            color: #666;
        }}
        .count {{
            margin-bottom: 20px;
            font-size: 14px;
            color: #555;
        }}
    </style>
</head>
<body>
    <h1>{title_text}</h1>
    <p class="count">Showing {actual_n} proteins</p>
    <table>
        <thead>
            <tr>
                <th>Rank</th>
                <th>Accession</th>
                <th>Gene</th>
                <th>log₂FC</th>
                <th>Adj. p-value</th>
                <th>Source</th>
            </tr>
        </thead>
        <tbody>
'''
    
    for _, row in table_df.iterrows():
        html_content += f'''            <tr>
                <td>{row['Rank']}</td>
                <td>{row['Accession']}</td>
                <td class="gene">{row['Gene']}</td>
                <td>{row['log2FC']}</td>
                <td>{row['Adj. p-value']}</td>
                <td>{row['Source']}</td>
            </tr>
'''
    
    html_content += f'''        </tbody>
    </table>
    <p class="footer">Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")} | 
    Ketamine Proteomics Analysis Project<br>
    Source: ketamine = high confidence in ketamine only | control = high confidence in control only | both = high confidence in both</p>
</body>
</html>
'''
    
    html_path = os.path.join(output_dir, f'{file_prefix}_table.html')
    with open(html_path, 'w') as f:
        f.write(html_content)
    
    log_message(log_path, f'  Saved: {html_path}')
    
    # ---------------------------------------------------------------------
    # Create PDF table using plotly
    # ---------------------------------------------------------------------
    fig = go.Figure(data=[go.Table(
        header=dict(
            values=['<b>Rank</b>', '<b>Accession</b>', '<b>Gene</b>', 
                    '<b>log₂FC</b>', '<b>Adj. p-value</b>', '<b>Source</b>'],
            fill_color=header_color,
            font=dict(color='white', size=11),
            align='left',
            height=30
        ),
        cells=dict(
            values=[
                table_df['Rank'],
                table_df['Accession'],
                table_df['Gene'],
                table_df['log2FC'],
                table_df['Adj. p-value'],
                table_df['Source']
            ],
            fill_color=[['white', '#f8f9fa'] * (actual_n // 2 + 1)],
            font=dict(size=10),
            align='left',
            height=25
        )
    )])
    
    fig.update_layout(
        title=dict(
            text=f'{title_text}<br><sup>Showing {actual_n} proteins | Source: ketamine/control/both = high confidence detection file</sup>',
            font=dict(size=14, family='Arial Black')
        ),
        width=850,
        height=max(400, 50 + actual_n * 28),  # Dynamic height based on rows
        margin=dict(t=80, b=20, l=20, r=20)
    )
    
    pdf_path = os.path.join(output_dir, f'{file_prefix}_table.pdf')
    fig.write_image(pdf_path, scale=2)
    
    log_message(log_path, f'  Saved: {pdf_path}')


# =============================================================================
# MAIN FUNCTION
# =============================================================================

def main():
    '''Main entry point for visualization script.'''
    
    # ---------------------------------------------------------------------
    # Parse arguments
    # ---------------------------------------------------------------------
    parser = argparse.ArgumentParser(
        description='Generate visualizations for differential abundance results'
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
        help='Output directory for figures'
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
        help='P-value threshold for significance (default: 0.05)'
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
    log_message(log_path, 'Visualization Generation')
    log_message(log_path, '=' * 50)
    log_message(log_path, f'Input file:        {args.input}')
    log_message(log_path, f'Output directory:  {args.output_dir}')
    log_message(log_path, f'P-value threshold: {args.pval_threshold}')
    log_message(log_path, '')
    
    # ---------------------------------------------------------------------
    # Load and prepare data
    # ---------------------------------------------------------------------
    df = load_and_prepare_data(args.input, args.pval_threshold, log_path)
    
    # ---------------------------------------------------------------------
    # Generate visualizations
    # ---------------------------------------------------------------------
    create_volcano_plot(df, args.pval_threshold, args.output_dir, log_path)
    create_ma_plot(df, args.pval_threshold, args.output_dir, log_path)
    create_summary_bar_chart(df, args.pval_threshold, args.output_dir, log_path)
    create_top_proteins_bar_chart(df, args.output_dir, log_path, n_proteins=10)
    create_top_proteins_table(df, args.output_dir, log_path, direction='up', n_proteins=20)
    create_top_proteins_table(df, args.output_dir, log_path, direction='down', n_proteins=20)
    
    # ---------------------------------------------------------------------
    # Finalize
    # ---------------------------------------------------------------------
    log_message(log_path, '')
    log_message(log_path, '=' * 50)
    log_message(log_path, 'Visualization generation complete!')
    log_message(log_path, '=' * 50)


if __name__ == '__main__':
    main()
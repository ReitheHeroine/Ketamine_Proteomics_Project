# title: visualize_diff_abundance.py
# project: Ketamine Proteomics Analysis Project
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2025-12-09
# last modified: 2026-05-27 (restyled create_summary_bar_chart per style guide;
#                            volcano text bumped to ~1.5x slide variant and set
#                            to black for presentation readability)
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
#   ├── volcano_plot.pdf/.png/.html
#   ├── ma_plot.pdf/.png/.html
#   ├── summary_bar_chart.pdf/.png
#   ├── top_proteins_bar_chart.pdf/.png/.html
#   ├── top20_upregulated_table.pdf/.png/.html
#   ├── top20_downregulated_table.pdf/.png/.html
#   ├── variability_fc_pvalue_relationship.pdf/.png/.html
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
#   - Variability/FC/p-value plot shows statistical consistency of Proteome
#     Discoverer output (quantitative proteins only)

import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
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
VARIABILITY_COL = 'Abundance Ratio Variability [%]: (ketamine) / (control)'

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
    Labels: Ina (most upregulated), the second most upregulated protein
    by log2FC, and all significantly downregulated proteins.

    Styled per project_notes/figure_and_table_style_guide.Rmd (2026-05-27).
    Uses local color overrides rather than the module-level COLOR_UP/DOWN/NS
    constants because the other figures in this script have not yet been
    restyled (per the guide's "aspirational, not retroactive" rule).
    '''
    log_message(log_path, 'Creating volcano plot...')

    # ---------------------------------------------------------------------
    # Local style constants (figure_and_table_style_guide.Rmd, 2026-05-27)
    # ---------------------------------------------------------------------
    SG_UP = '#E8735A'         # up-in-ketamine = ketamine series color (coral)
    SG_DOWN = '#7FB3D8'       # down-in-ketamine = control series color (light blue)
    SG_NS = '#999999'         # not significant (gray)
    SG_THRESHOLD = '#666666'  # threshold / zero reference lines
    SG_TEXT = '#000000'       # all text (presentation override; style guide default is #2C3E50)
    SG_ARROW = '#888888'      # annotation leader lines (mid gray)
    SG_GRID = '#E0E0E0'       # faint reference grid lines
    FONT_FAMILY = 'Arial'

    # ---------------------------------------------------------------------
    # Filter to quantitative proteins only
    # ---------------------------------------------------------------------
    quant_df = df[df[CATEGORY_COL] == 'quantitative'].copy()

    fig = go.Figure()

    # ---------------------------------------------------------------------
    # Plot each group in order
    # ---------------------------------------------------------------------
    legend_order = [
        ('Up in Ketamine (sig.)', SG_UP),
        ('Down in Ketamine (sig.)', SG_DOWN),
        ('Not significant', SG_NS),
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
                size=14,
                line=dict(width=0.75, color='white')
            ),
            name=f'{legend_name} (n={len(subset)})',
            text=subset.apply(
                lambda r: f"Gene: {r[GENE_COL]}<br>"
                          f"Accession: {r[ACCESSION_COL]}<br>"
                          f"Log<sub>2</sub> FC: {r[LOG2FC_COL]:.2f}<br>"
                          f"Adj. <i>p</i>-value: {r[PVAL_COL]:.2e}",
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
        line_color=SG_THRESHOLD,
        line_width=1,
        annotation_text=f'<i>p</i> = {pval_threshold}',
        annotation_position='top right',
        annotation_font=dict(size=18, color=SG_TEXT, family=FONT_FAMILY)
    )

    # ---------------------------------------------------------------------
    # Label a curated subset of proteins:
    #   - Ina (most upregulated by log2FC)
    #   - Second most upregulated protein by log2FC
    #   - All significantly downregulated proteins
    # Per style guide section 4.3, protein labels in mass-spec figures are
    # set in all-caps roman (e.g., INA, PLP1), not the mouse-gene convention
    # (Ina, Plp1) used in the underlying data file.
    # ---------------------------------------------------------------------
    sig_df = quant_df[quant_df['significant']].copy()
    sig_up = sig_df[sig_df[DIRECTION_COL] == 'up_in_ketamine']
    sig_down = sig_df[sig_df[DIRECTION_COL] == 'down_in_ketamine']

    # Top 2 upregulated by log2FC (expected: Ina first, then second-highest)
    top_up = sig_up.nlargest(2, LOG2FC_COL)

    # Combine: top 2 upregulated + all downregulated
    labeled = pd.concat([top_up, sig_down])

    log_message(log_path, f'  Labeling {len(labeled)} proteins on volcano plot:')
    for _, row in labeled.iterrows():
        log_message(log_path,
                    f'    {row[GENE_COL]}: log2FC={row[LOG2FC_COL]:.2f}, '
                    f'p={row[PVAL_COL]:.2e}')

    for _, row in labeled.iterrows():
        fig.add_annotation(
            x=row[LOG2FC_COL],
            y=row['neg_log10_pval'],
            text=row[GENE_COL].upper(),
            showarrow=True,
            arrowhead=0,
            arrowsize=0.5,
            arrowwidth=1,
            arrowcolor=SG_ARROW,
            ax=30,
            ay=-30,
            font=dict(size=18, color=SG_TEXT, family=FONT_FAMILY),
            bgcolor='rgba(255,255,255,0.7)',
            borderpad=2
        )

    # ---------------------------------------------------------------------
    # Layout
    # ---------------------------------------------------------------------
    # No on-figure title per style guide section 11.2: the title sentence
    # belongs in the Word document caption block, not inside the exported
    # figure. Axis titles remain in-figure per section 11.1.
    fig.update_layout(
        font=dict(family=FONT_FAMILY, size=22, color=SG_TEXT),
        xaxis=dict(
            title=dict(
                text='Log<sub>2</sub> fold change',
                font=dict(family=FONT_FAMILY, size=22, color=SG_TEXT)
            ),
            tickfont=dict(family=FONT_FAMILY, size=20, color=SG_TEXT),
            zeroline=True,
            zerolinecolor=SG_GRID,
            zerolinewidth=1,
            gridcolor=SG_GRID,
            showline=True,
            linecolor='black',
            linewidth=2,
            mirror=True,
            ticks='outside',
            tickcolor='black',
            tickwidth=2,
            ticklen=6
        ),
        yaxis=dict(
            title=dict(
                text='-Log<sub>10</sub>(adjusted <i>p</i>-value)',
                font=dict(family=FONT_FAMILY, size=22, color=SG_TEXT)
            ),
            tickfont=dict(family=FONT_FAMILY, size=20, color=SG_TEXT),
            gridcolor=SG_GRID,
            showline=True,
            linecolor='black',
            linewidth=2,
            mirror=True,
            ticks='outside',
            tickcolor='black',
            tickwidth=2,
            ticklen=6
        ),
        # Legend placed inside the plot in the top-left empty quadrant
        # (negative log2FC, high -log10(p) region is sparsely populated).
        # x is offset enough to clear the y-axis line and tick labels.
        legend=dict(
            font=dict(family=FONT_FAMILY, size=20, color=SG_TEXT),
            yanchor='top',
            y=0.98,
            xanchor='left',
            x=0.10,
            bgcolor='rgba(255,255,255,0.95)',
            bordercolor='#D0D0D0',
            borderwidth=1
        ),
        plot_bgcolor='white',
        paper_bgcolor='white',
        width=1100,
        height=800,
        margin=dict(r=80, t=40, b=90, l=120)
    )
    
    # ---------------------------------------------------------------------
    # Save outputs
    # ---------------------------------------------------------------------
    html_path = os.path.join(output_dir, 'volcano_plot.html')
    pdf_path = os.path.join(output_dir, 'volcano_plot.pdf')
    png_path = os.path.join(output_dir, 'volcano_plot.png')

    fig.write_html(html_path)
    fig.write_image(pdf_path, scale=2)
    fig.write_image(png_path, scale=2)

    log_message(log_path, f'  Saved: {html_path}')
    log_message(log_path, f'  Saved: {pdf_path}')
    log_message(log_path, f'  Saved: {png_path}')


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
    png_path = os.path.join(output_dir, 'ma_plot.png')

    fig.write_html(html_path)
    fig.write_image(pdf_path, scale=2)
    fig.write_image(png_path, scale=2)

    log_message(log_path, f'  Saved: {html_path}')
    log_message(log_path, f'  Saved: {pdf_path}')
    log_message(log_path, f'  Saved: {png_path}')


def create_summary_bar_chart(df, pval_threshold, output_dir, log_path):
    '''
    Create summary bar chart showing counts of differentially abundant
    proteins by category (quantitative up/down + presence/absence).

    Styled per project_notes/figure_and_table_style_guide.Rmd (2026-05-27).
    Uses local color overrides rather than the module-level COLOR_UP/DOWN/
    PA_* constants because the other figures in this script have not yet
    been restyled (per the guide's "aspirational, not retroactive" rule).

    Palette:
      - Upregulated (quant)   -> #E8735A (style guide sec. 9.1.1: ketamine series)
      - Downregulated (quant) -> #7FB3D8 (style guide sec. 9.1.1: control series)
      - Ketamine-specific P/A -> #FF7F0E (style guide sec. 9.1.3: orange)
      - Control-specific P/A  -> #9467BD (style guide sec. 9.1.3: purple)

    Four distinct colors fully disambiguate the categories on their own, so
    no pattern fill is needed (style guide sec. 9.2). Per section 11.1-11.2,
    the figure interior carries only the data and minimum labels; the figure
    title, abbreviation definitions, p-value threshold, and statistical-
    method statement live in the Word document caption block (ordering per
    section 11.4).

    Type sizes follow the manuscript defaults from style guide section 3
    (as updated 2026-05-27): axis title 14 pt, tick labels 12 pt, data
    labels 10 pt. Figure exported at 650 x 400 px, scale = 2. For a slide-
    deck rendering, scale every size by ~1.25x to the slide-variant column
    of the section 3 table (17 / 15 / 12 pt) and bump the bar outline
    weight proportionally.
    '''
    log_message(log_path, 'Creating summary bar chart...')

    # ---------------------------------------------------------------------
    # Local style constants (figure_and_table_style_guide.Rmd, 2026-05-27)
    # ---------------------------------------------------------------------
    SG_UP = '#E8735A'         # up-in-ketamine (coral, sec. 9.1.1)
    SG_DOWN = '#7FB3D8'       # down-in-ketamine (light blue, sec. 9.1.1)
    SG_PA_KET = '#FF7F0E'     # ketamine-specific P/A (orange, sec. 9.1.3)
    SG_PA_CTRL = '#9467BD'    # control-specific P/A (purple, sec. 9.1.3)
    SG_TEXT = '#000000'       # all figure-interior text (black, sec. 9.1.4)
    SG_OUTLINE = '#2C3E50'    # bar / marker outlines (dark slate, sec. 9.1.4)
    SG_GRID = '#E0E0E0'       # faint reference grid lines (sec. 9.1.4)
    FONT_FAMILY = 'Arial'

    # ---------------------------------------------------------------------
    # Calculate counts
    # ---------------------------------------------------------------------
    quant_df = df[df[CATEGORY_COL] == 'quantitative']
    sig_up = len(quant_df[(quant_df['significant']) & (quant_df[DIRECTION_COL] == 'up_in_ketamine')])
    sig_down = len(quant_df[(quant_df['significant']) & (quant_df[DIRECTION_COL] == 'down_in_ketamine')])

    ket_specific = len(df[df[CATEGORY_COL] == 'presence_absence_ketamine_specific'])
    ctrl_specific = len(df[df[CATEGORY_COL] == 'presence_absence_control_specific'])

    # ---------------------------------------------------------------------
    # Build figure
    # Sentence-case category labels per style guide sec. 5.
    # ---------------------------------------------------------------------
    categories = [
        'Upregulated<br>(quantitative)',
        'Downregulated<br>(quantitative)',
        'Ketamine-<br>specific',
        'Control-<br>specific',
    ]
    counts = [sig_up, sig_down, ket_specific, ctrl_specific]
    colors = [SG_UP, SG_DOWN, SG_PA_KET, SG_PA_CTRL]

    fig = go.Figure()

    fig.add_trace(go.Bar(
        x=categories,
        y=counts,
        marker=dict(
            color=colors,
            line=dict(width=0.6, color=SG_OUTLINE)
        ),
        text=counts,
        textposition='outside',
        textfont=dict(family=FONT_FAMILY, size=10, color=SG_TEXT),
        cliponaxis=False
    ))

    # ---------------------------------------------------------------------
    # Layout (no in-figure title; per sec. 11.2, figure number and title
    # sentence live in the Word document caption block, not the raster).
    # Type sizes per style guide section 3 (manuscript defaults, updated
    # 2026-05-27).
    # ---------------------------------------------------------------------
    max_count = max(counts) if counts else 1

    fig.update_layout(
        font=dict(family=FONT_FAMILY, size=14, color=SG_TEXT),
        xaxis=dict(
            tickfont=dict(family=FONT_FAMILY, size=12, color=SG_TEXT),
            showgrid=False,
            showline=True,
            linecolor=SG_TEXT,
            linewidth=1,
            mirror=True  # draw opposing border at the top
        ),
        yaxis=dict(
            title=dict(
                text='Number of proteins',
                font=dict(family=FONT_FAMILY, size=14, color=SG_TEXT)
            ),
            tickfont=dict(family=FONT_FAMILY, size=12, color=SG_TEXT),
            gridcolor=SG_GRID,
            range=[0, max_count * 1.2],  # headroom for outside data labels
            showline=True,
            linecolor=SG_TEXT,
            linewidth=1,
            mirror=True  # draw opposing border on the right
        ),
        plot_bgcolor='white',
        paper_bgcolor='white',
        width=650,
        height=400,
        showlegend=False,
        bargap=0.3,
        margin=dict(t=30, b=90, l=80, r=40)
    )

    # ---------------------------------------------------------------------
    # Save outputs (scale=2 per style guide sec. 14)
    # ---------------------------------------------------------------------
    pdf_path = os.path.join(output_dir, 'summary_bar_chart.pdf')
    png_path = os.path.join(output_dir, 'summary_bar_chart.png')
    fig.write_image(pdf_path, scale=2)
    fig.write_image(png_path, scale=2)

    log_message(log_path, f'  Saved: {pdf_path}')
    log_message(log_path, f'  Saved: {png_path}')


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
    png_path = os.path.join(output_dir, 'top_proteins_bar_chart.png')

    fig.write_html(html_path)
    fig.write_image(pdf_path, scale=2)
    fig.write_image(png_path, scale=2)

    log_message(log_path, f'  Saved: {html_path}')
    log_message(log_path, f'  Saved: {pdf_path}')
    log_message(log_path, f'  Saved: {png_path}')


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
    png_path = os.path.join(output_dir, f'{file_prefix}_table.png')
    fig.write_image(pdf_path, scale=2)
    fig.write_image(png_path, scale=2)

    log_message(log_path, f'  Saved: {pdf_path}')
    log_message(log_path, f'  Saved: {png_path}')


# =============================================================================
# VARIABILITY vs FOLD CHANGE vs P-VALUE RELATIONSHIP PLOT
# =============================================================================

def create_variability_fc_pvalue_plot(df, pval_threshold, output_dir, log_path):
    '''
    Create scatter plot showing the three-way relationship between fold change,
    variability %, and p-value for quantitative proteins.

    Purpose:
      1. Validate that Proteome Discoverer statistics behave as expected
         (high-variability proteins should need larger FC for significance)
      2. Communicate n=3 statistical behavior to thesis committee

    Layout:
      - Central scatter: log2(FC) on X-axis, Variability % on Y-axis
      - Points colored by -log10(p-value) with continuous color scale
      - Significant proteins outlined with black border
      - Key proteins of interest labeled
      - Marginal histograms on top and right edges

    Red flags checked:
      - High variability (>60%) + low p-value (<0.01) at modest |FC| (<1)
    '''
    log_message(log_path, 'Creating variability vs FC vs p-value relationship plot...')

    # -----------------------------------------------------------------
    # Step 1: Filter to quantitative proteins with valid variability
    # -----------------------------------------------------------------
    quant_df = df[df[CATEGORY_COL] == 'quantitative'].copy()
    quant_df = quant_df[quant_df[VARIABILITY_COL].notna()].copy()

    log_message(log_path, f'  Quantitative proteins with variability data: {len(quant_df)}')

    # -----------------------------------------------------------------
    # Step 2: Compute derived columns
    # -----------------------------------------------------------------
    quant_df['abs_log2fc'] = quant_df[LOG2FC_COL].abs()
    quant_df['neg_log10_pval'] = -np.log10(quant_df[PVAL_COL].replace(0, 1e-20))

    # -----------------------------------------------------------------
    # Step 3: Red flag check - high variability + low p at modest FC
    # -----------------------------------------------------------------
    red_flag_mask = (
        (quant_df[VARIABILITY_COL] > 60) &
        (quant_df[PVAL_COL] < 0.01) &
        (quant_df['abs_log2fc'] < 1)
    )
    n_red_flags = red_flag_mask.sum()

    if n_red_flags > 0:
        log_message(log_path, f'  *** RED FLAG: {n_red_flags} protein(s) with '
                    f'variability >60%, p<0.01, and |log2FC|<1 ***')
        flagged = quant_df[red_flag_mask][[GENE_COL, LOG2FC_COL, VARIABILITY_COL, PVAL_COL]]
        for _, row in flagged.iterrows():
            log_message(log_path,
                        f'    {row[GENE_COL]}: log2FC={row[LOG2FC_COL]:.3f}, '
                        f'var={row[VARIABILITY_COL]:.1f}%, '
                        f'p={row[PVAL_COL]:.2e}')
    else:
        log_message(log_path, '  Red flag check PASSED: no high-variability + '
                    'low-p + modest-FC proteins found')

    # -----------------------------------------------------------------
    # Step 4: Create figure with marginal histograms
    # -----------------------------------------------------------------
    fig = make_subplots(
        rows=2, cols=2,
        column_widths=[0.82, 0.18],
        row_heights=[0.18, 0.82],
        shared_xaxes=True,
        shared_yaxes=True,
        horizontal_spacing=0.02,
        vertical_spacing=0.02
    )

    # -----------------------------------------------------------------
    # Step 5: Main scatter plot (row=2, col=1)
    # -----------------------------------------------------------------
    # Split into significant and non-significant for distinct marker styles
    sig_mask = quant_df['significant']
    ns_df = quant_df[~sig_mask]
    sig_df = quant_df[sig_mask]

    # 5a. Non-significant proteins (no border)
    fig.add_trace(go.Scatter(
        x=ns_df[LOG2FC_COL],
        y=ns_df[VARIABILITY_COL],
        mode='markers',
        marker=dict(
            color=ns_df['neg_log10_pval'],
            colorscale='Viridis',
            size=7,
            opacity=0.6,
            cmin=0,
            cmax=quant_df['neg_log10_pval'].quantile(0.95),
            line=dict(width=0),
            showscale=False
        ),
        name=f'Not significant (n={len(ns_df)})',
        text=ns_df.apply(
            lambda r: f"Gene: {r[GENE_COL]}<br>"
                      f"log2FC: {r[LOG2FC_COL]:.3f}<br>"
                      f"Variability: {r[VARIABILITY_COL]:.1f}%<br>"
                      f"p-value: {r[PVAL_COL]:.2e}<br>"
                      f"-log10(p): {r['neg_log10_pval']:.2f}",
            axis=1
        ),
        hoverinfo='text',
        showlegend=True
    ), row=2, col=1)

    # 5b. Significant proteins (black border)
    fig.add_trace(go.Scatter(
        x=sig_df[LOG2FC_COL],
        y=sig_df[VARIABILITY_COL],
        mode='markers',
        marker=dict(
            color=sig_df['neg_log10_pval'],
            colorscale='Viridis',
            size=9,
            opacity=0.9,
            cmin=0,
            cmax=quant_df['neg_log10_pval'].quantile(0.95),
            line=dict(width=1.5, color='black'),
            colorbar=dict(
                title=dict(text='-log₁₀(p)', font=dict(size=12)),
                x=1.02,
                len=0.6,
                y=0.35,
                thickness=15
            )
        ),
        name=f'Significant, p ≤ {pval_threshold} (n={len(sig_df)})',
        text=sig_df.apply(
            lambda r: f"Gene: {r[GENE_COL]}<br>"
                      f"log2FC: {r[LOG2FC_COL]:.3f}<br>"
                      f"Variability: {r[VARIABILITY_COL]:.1f}%<br>"
                      f"p-value: {r[PVAL_COL]:.2e}<br>"
                      f"-log10(p): {r['neg_log10_pval']:.2f}",
            axis=1
        ),
        hoverinfo='text',
        showlegend=True
    ), row=2, col=1)

    # -----------------------------------------------------------------
    # Step 6: Label key proteins of interest
    # -----------------------------------------------------------------
    # Tier 1 SNARE/vesicle + contamination examples
    label_genes = ['Snap25', 'Syt1', 'Stx1a', 'Stxbp1', 'Vamp2',
                   'Ina', 'Plp1', 'Sv2a']

    labeled = quant_df[quant_df[GENE_COL].isin(label_genes)]

    for _, row in labeled.iterrows():
        fig.add_annotation(
            x=row[LOG2FC_COL],
            y=row[VARIABILITY_COL],
            text=row[GENE_COL],
            showarrow=True,
            arrowhead=0,
            arrowsize=0.5,
            arrowwidth=1,
            ax=25,
            ay=-18,
            font=dict(size=9, color='black'),
            bgcolor='rgba(255,255,255,0.8)',
            borderpad=2,
            xref='x',
            yref='y',
            row=2, col=1
        )

    # -----------------------------------------------------------------
    # Step 7: Marginal histogram - top (log2FC distribution)
    # -----------------------------------------------------------------
    fig.add_trace(go.Histogram(
        x=quant_df[LOG2FC_COL],
        nbinsx=50,
        marker_color='rgba(100, 100, 100, 0.5)',
        showlegend=False,
        hoverinfo='skip'
    ), row=1, col=1)

    # -----------------------------------------------------------------
    # Step 8: Marginal histogram - right (variability distribution)
    # -----------------------------------------------------------------
    fig.add_trace(go.Histogram(
        y=quant_df[VARIABILITY_COL],
        nbinsy=40,
        marker_color='rgba(100, 100, 100, 0.5)',
        showlegend=False,
        hoverinfo='skip'
    ), row=2, col=2)

    # -----------------------------------------------------------------
    # Step 9: Layout and styling
    # -----------------------------------------------------------------
    fig.update_layout(
        title=dict(
            text='Fold Change vs. Variability vs. Statistical Significance',
            font=dict(size=16, family='Arial Black'),
            x=0.4,
            y=0.98
        ),
        plot_bgcolor='white',
        width=950,
        height=750,
        margin=dict(t=60, b=80, l=70, r=120),
        legend=dict(
            yanchor='top',
            y=0.99,
            xanchor='left',
            x=0.01,
            bgcolor='rgba(255,255,255,0.9)',
            bordercolor='lightgray',
            borderwidth=1,
            font=dict(size=10)
        )
    )

    # Main scatter axes
    fig.update_xaxes(
        title_text='log₂(Fold Change)',
        title_font=dict(size=13),
        gridcolor='rgba(0,0,0,0.08)',
        zeroline=True,
        zerolinecolor='rgba(0,0,0,0.2)',
        zerolinewidth=1,
        row=2, col=1
    )
    fig.update_yaxes(
        title_text='Variability [%]',
        title_font=dict(size=13),
        gridcolor='rgba(0,0,0,0.08)',
        row=2, col=1
    )

    # Hide marginal axes ticks/labels
    fig.update_xaxes(showticklabels=False, row=1, col=1)
    fig.update_yaxes(showticklabels=False, row=1, col=1)
    fig.update_xaxes(showticklabels=False, row=2, col=2)
    fig.update_yaxes(showticklabels=False, row=2, col=2)

    # White background for marginal panels
    fig.update_xaxes(gridcolor='rgba(0,0,0,0)', row=1, col=1)
    fig.update_yaxes(gridcolor='rgba(0,0,0,0)', row=1, col=1)
    fig.update_xaxes(gridcolor='rgba(0,0,0,0)', row=2, col=2)
    fig.update_yaxes(gridcolor='rgba(0,0,0,0)', row=2, col=2)

    # Hide empty subplot (top-right corner)
    fig.update_xaxes(visible=False, row=1, col=2)
    fig.update_yaxes(visible=False, row=1, col=2)

    # -----------------------------------------------------------------
    # Step 10: Add annotation with interpretation guidance
    # -----------------------------------------------------------------
    fig.add_annotation(
        text=('Color = -log₁₀(p-value). Black-bordered points = significant (p ≤ 0.05).<br>'
              'Quantitative proteins only (n=' + str(len(quant_df)) + '). '
              'Presence/absence proteins excluded (no variability data).'),
        xref='paper', yref='paper',
        x=0.4, y=-0.08,
        showarrow=False,
        font=dict(size=10, color='gray'),
        align='center'
    )

    # -----------------------------------------------------------------
    # Step 11: Compute and log summary statistics for QC
    # -----------------------------------------------------------------
    # Correlation between |FC| and -log10(p)
    corr_fc_p = np.corrcoef(quant_df['abs_log2fc'], quant_df['neg_log10_pval'])[0, 1]
    # Correlation between variability and -log10(p)
    corr_var_p = np.corrcoef(quant_df[VARIABILITY_COL], quant_df['neg_log10_pval'])[0, 1]
    # Correlation between variability and |FC|
    corr_var_fc = np.corrcoef(quant_df[VARIABILITY_COL], quant_df['abs_log2fc'])[0, 1]

    log_message(log_path, f'  Correlations (quantitative proteins, n={len(quant_df)}):')
    log_message(log_path, f'    |log2FC| vs -log10(p): r = {corr_fc_p:.3f}')
    log_message(log_path, f'    Variability vs -log10(p): r = {corr_var_p:.3f}')
    log_message(log_path, f'    Variability vs |log2FC|:  r = {corr_var_fc:.3f}')

    # Variability ranges among significant proteins
    if len(sig_df) > 0:
        log_message(log_path, f'  Significant proteins variability range: '
                    f'{sig_df[VARIABILITY_COL].min():.1f}% - '
                    f'{sig_df[VARIABILITY_COL].max():.1f}% '
                    f'(median: {sig_df[VARIABILITY_COL].median():.1f}%)')
        log_message(log_path, f'  Non-significant proteins variability range: '
                    f'{ns_df[VARIABILITY_COL].min():.1f}% - '
                    f'{ns_df[VARIABILITY_COL].max():.1f}% '
                    f'(median: {ns_df[VARIABILITY_COL].median():.1f}%)')

    # -----------------------------------------------------------------
    # Step 12: Save outputs
    # -----------------------------------------------------------------
    html_path = os.path.join(output_dir, 'variability_fc_pvalue_relationship.html')
    pdf_path = os.path.join(output_dir, 'variability_fc_pvalue_relationship.pdf')
    png_path = os.path.join(output_dir, 'variability_fc_pvalue_relationship.png')

    fig.write_html(html_path)
    fig.write_image(pdf_path, scale=2)
    fig.write_image(png_path, scale=2)

    log_message(log_path, f'  Saved: {html_path}')
    log_message(log_path, f'  Saved: {pdf_path}')
    log_message(log_path, f'  Saved: {png_path}')


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
    create_variability_fc_pvalue_plot(df, args.pval_threshold, args.output_dir, log_path)

    # ---------------------------------------------------------------------
    # Finalize
    # ---------------------------------------------------------------------
    log_message(log_path, '')
    log_message(log_path, '=' * 50)
    log_message(log_path, 'Visualization generation complete!')
    log_message(log_path, '=' * 50)


if __name__ == '__main__':
    main()
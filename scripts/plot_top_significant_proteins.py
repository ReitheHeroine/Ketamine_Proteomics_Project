# title: plot_top_significant_proteins.py
# project: Ketamine Astrocyte Proteomics
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-03-19
# last modified: 2026-03-19
#
# purpose:
#   Generate a horizontal bar chart of the top significant proteins ranked
#   by adjusted p-value, excluding suspected contaminant proteins (neuronal
#   and oligodendrocyte markers detected in the dataset).
#
# inputs:
#   - data/all_proteins_categorized.csv
#   - data/cell_type_markers.csv
#
# outputs:
#   results/figures/
#   ├── top_significant_proteins_no_contaminants.png
#   ├── top_significant_proteins_no_contaminants.pdf
#   └── top_significant_proteins_no_contaminants.html
#
# usage:
#   python scripts/plot_top_significant_proteins.py
#   python scripts/plot_top_significant_proteins.py --n 20
#   python scripts/plot_top_significant_proteins.py --include-contaminants

import pandas as pd
import numpy as np
import plotly.graph_objects as go
import argparse
import os

# ======================================================================
# 1. CONFIGURATION
# ======================================================================

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROTEINS_FILE = os.path.join(BASE_DIR, 'data', 'all_proteins_categorized.csv')
MARKERS_FILE = os.path.join(BASE_DIR, 'data', 'cell_type_markers.csv')
OUTPUT_DIR = os.path.join(BASE_DIR, 'results', 'figures')

COL_GENE = 'Gene Symbol'
COL_LOG2FC = 'log2_fold_change'
COL_PVAL = 'Abundance Ratio Adj. P-Value: (ketamine) / (control)'
COL_CATEGORY = 'category'
COL_SIGNIFICANT = 'significant'
COL_DIRECTION = 'direction'

# Cell types considered suspected contaminants
CONTAMINANT_CELL_TYPES = ['Neuron', 'Oligodendrocyte']

# ======================================================================
# 2. PARSE ARGUMENTS
# ======================================================================

parser = argparse.ArgumentParser(
    description='Bar chart of top significant proteins by adjusted p-value, '
                'excluding suspected contaminants.'
)
parser.add_argument('--n', type=int, default=25,
                    help='Number of top proteins to display (default: 25)')
parser.add_argument('--include-contaminants', action='store_true',
                    help='Include suspected contaminant markers in the figure')
args = parser.parse_args()

# ======================================================================
# 3. LOAD DATA AND FILTER
# ======================================================================

print('Loading data...')
proteins_df = pd.read_csv(PROTEINS_FILE)
markers_df = pd.read_csv(MARKERS_FILE)

# -- identify suspected contaminant gene symbols --
contaminant_markers = markers_df[
    markers_df['cell_type'].isin(CONTAMINANT_CELL_TYPES)
]
contaminant_genes = set(contaminant_markers['gene_symbol'].str.lower())

# -- filter to significant quantitative proteins --
sig_df = proteins_df[
    (proteins_df[COL_SIGNIFICANT] == True) &
    (proteins_df[COL_CATEGORY] == 'quantitative')
].copy()

print(f'Significant quantitative proteins: {len(sig_df)}')

# -- exclude contaminants unless flag is set --
if not args.include_contaminants:
    contam_mask = sig_df[COL_GENE].str.lower().isin(contaminant_genes)
    excluded = sig_df[contam_mask]
    sig_df = sig_df[~contam_mask].copy()
    print(f'Excluded {len(excluded)} suspected contaminant markers:')
    for _, row in excluded.iterrows():
        ct = contaminant_markers[
            contaminant_markers['gene_symbol'].str.lower() == row[COL_GENE].lower()
        ].iloc[0]['cell_type']
        print(f'  {row[COL_GENE]:10s}  ({ct})')
    print(f'Remaining significant proteins: {len(sig_df)}')

# -- rank by adjusted p-value (most significant first) --
sig_df = sig_df.sort_values(COL_PVAL, ascending=True)
top_df = sig_df.head(args.n).copy()

print(f'\nTop {len(top_df)} proteins by adjusted p-value:')

# -- reverse for horizontal bar chart (most significant at top) --
top_df = top_df.iloc[::-1].reset_index(drop=True)

# ======================================================================
# 4. BUILD FIGURE
# ======================================================================

print('\nCreating figure...')

# -- compute -log10(p-value) for display --
top_df['neg_log10_pval'] = -np.log10(top_df[COL_PVAL].replace(0, 1e-20))

# -- color by direction --
colors = []
for _, row in top_df.iterrows():
    if row[COL_DIRECTION] == 'up_in_ketamine':
        colors.append('#E8735A')   # red for upregulated
    elif row[COL_DIRECTION] == 'down_in_ketamine':
        colors.append('#7FB3D8')   # blue for downregulated
    else:
        colors.append('#999999')

fig = go.Figure()

fig.add_trace(go.Bar(
    y=top_df[COL_GENE],
    x=top_df['neg_log10_pval'],
    orientation='h',
    marker_color=colors,
    marker_line=dict(width=0.5, color='#555555'),
    hovertemplate=(
        '<b>%{y}</b><br>'
        '-log₁₀(adj. p): %{x:.2f}<br>'
        'log₂FC: %{customdata[0]:.3f}<br>'
        'adj. p-value: %{customdata[1]:.2e}'
        '<extra></extra>'
    ),
    customdata=list(zip(top_df[COL_LOG2FC], top_df[COL_PVAL])),
    showlegend=False,
))

# -- significance threshold line --
threshold_x = -np.log10(0.05)
fig.add_vline(
    x=threshold_x,
    line=dict(color='#888888', width=1.5, dash='dot'),
    annotation_text='p = 0.05',
    annotation_position='top',
    annotation_font=dict(size=9, color='#888888'),
)

# -- legend traces for direction --
fig.add_trace(go.Bar(
    y=[None], x=[None],
    marker_color='#E8735A',
    name='Upregulated in ketamine',
    showlegend=True,
))
fig.add_trace(go.Bar(
    y=[None], x=[None],
    marker_color='#7FB3D8',
    name='Downregulated in ketamine',
    showlegend=True,
))

# ======================================================================
# 5. LAYOUT
# ======================================================================

n_displayed = len(top_df)
contam_note = '' if args.include_contaminants else ' (excluding suspected contaminants)'
title_text = (
    f'Top {n_displayed} Significant Proteins by Adjusted P-Value'
    f'{contam_note}'
)

fig.update_layout(
    title=dict(
        text=title_text,
        font=dict(size=14, color='#2C3E50'),
        x=0.5,
        xanchor='center',
    ),
    xaxis=dict(
        title=dict(text='-log₁₀(adjusted p-value)', font=dict(size=12)),
        tickfont=dict(size=10),
    ),
    yaxis=dict(
        tickfont=dict(size=9),
    ),
    legend=dict(
        orientation='h',
        yanchor='bottom',
        y=1.02,
        xanchor='left',
        x=0.0,
        font=dict(size=11),
    ),
    plot_bgcolor='white',
    paper_bgcolor='white',
    width=800,
    height=max(500, n_displayed * 22 + 150),
    margin=dict(t=100, b=80, l=100, r=30),
    bargap=0.2,
    showlegend=True,
)

fig.update_xaxes(
    showgrid=True,
    gridcolor='#E8E8E8',
    gridwidth=1,
    zeroline=True,
    zerolinecolor='#CCCCCC',
)
fig.update_yaxes(showgrid=False)

# -- footnote --
if not args.include_contaminants:
    fig.add_annotation(
        xref='paper', yref='paper',
        x=0, y=-0.12,
        text='Excluded: neuronal and oligodendrocyte markers detected in dataset '
             '(see cell_type_markers.csv for full list).',
        showarrow=False,
        font=dict(size=8, color='#666666'),
        xanchor='left',
    )

# ======================================================================
# 6. SAVE OUTPUTS
# ======================================================================

os.makedirs(OUTPUT_DIR, exist_ok=True)

suffix = '' if not args.include_contaminants else '_with_contaminants'
png_path = os.path.join(OUTPUT_DIR, f'top_significant_proteins_no_contaminants{suffix}.png')
pdf_path = os.path.join(OUTPUT_DIR, f'top_significant_proteins_no_contaminants{suffix}.pdf')
html_path = os.path.join(OUTPUT_DIR, f'top_significant_proteins_no_contaminants{suffix}.html')

fig.write_image(png_path, scale=2)
fig.write_image(pdf_path, scale=2)
fig.write_html(html_path)

print(f'\nSaved: {png_path}')
print(f'Saved: {pdf_path}')
print(f'Saved: {html_path}')
print('\nDone.')

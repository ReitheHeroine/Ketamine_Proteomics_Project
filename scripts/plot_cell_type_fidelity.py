# title: plot_cell_type_fidelity.py
# project: Ketamine Astrocyte Proteomics
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-03-17
# last modified: 2026-03-17
#
# purpose:
#   Generate a grouped bar chart for the thesis proposal EDA/validation slide
#   that visually communicates both sample validation (astrocyte markers) and
#   contamination assessment in a single figure. Reads cell-type marker
#   definitions from data/cell_type_markers.csv and cross-references with
#   the proteomics dataset, including only markers that were detected.
#   Presence/absence proteins (NaN in one condition) are annotated rather
#   than plotted with misleading placeholder values.
#
# inputs:
#   - data/cell_type_markers.csv       (marker definitions by cell type)
#   - data/all_proteins_categorized.csv (proteomics dataset)
#
# outputs:
#   results/figures/
#   ├── cell_type_fidelity_barchart.png
#   ├── cell_type_fidelity_barchart.pdf
#   └── cell_type_fidelity_barchart.html
#
# usage example:
#   python scripts/plot_cell_type_fidelity.py              # full version (all detected markers)
#   python scripts/plot_cell_type_fidelity.py --simplified  # presentation version (curated subset)

import pandas as pd
import plotly.graph_objects as go
import numpy as np
import argparse
import os

# ======================================================================
# 1. CONFIGURATION
# ======================================================================

# -- file paths --
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MARKERS_FILE = os.path.join(BASE_DIR, 'data', 'cell_type_markers.csv')
PROTEINS_FILE = os.path.join(BASE_DIR, 'data', 'all_proteins_categorized.csv')
OUTPUT_DIR = os.path.join(BASE_DIR, 'results', 'figures')

# -- cell types to include (in display order) --
# Maps cell_type values from markers CSV to display labels
CELL_TYPE_ORDER = [
    ('Astrocyte', 'Astrocyte Markers'),
    ('Neuron', 'Neuronal Markers'),
    ('Oligodendrocyte', 'Oligodendrocyte Markers'),
]

# -- curated subset for presentation version --
# Most canonical, unambiguous markers per cell type
PRESENTATION_SUBSET = [
    'Gfap', 'Aqp4', 'Aldh1l1', 'Slc1a3', 'Glul', 'Gja1',  # astrocyte
    'Ina', 'Nefl', 'Nefm',                                    # neuronal
    'Plp1', 'Mbp', 'Mog',                                     # oligodendrocyte
]

# -- colors (matching existing figure style) --
COLOR_CONTROL = '#7FB3D8'       # light blue for control
COLOR_KETAMINE = '#E8735A'      # light red for ketamine

# -- column names in proteomics data --
COL_CONTROL = 'Abundances (Grouped): control'
COL_KETAMINE = 'Abundances (Grouped): ketamine'
COL_GENE = 'Gene Symbol'

# ======================================================================
# 2. LOAD AND CROSS-REFERENCE DATA
# ======================================================================

# ======================================================================
# PARSE ARGUMENTS
# ======================================================================

parser = argparse.ArgumentParser(description='Cell-type fidelity bar chart')
parser.add_argument('--simplified', action='store_true',
                    help='Generate simplified presentation version with curated marker subset')
args = parser.parse_args()

SIMPLIFIED = args.simplified

print('Loading marker definitions...')
markers_df = pd.read_csv(MARKERS_FILE)

print('Loading proteomics data...')
proteins_df = pd.read_csv(PROTEINS_FILE)

# -- build case-insensitive lookup from proteomics data --
proteins_df['gene_lower'] = proteins_df[COL_GENE].str.lower()

# -- filter markers to cell types we want --
cell_types_wanted = [ct for ct, _ in CELL_TYPE_ORDER]
markers_df = markers_df[markers_df['cell_type'].isin(cell_types_wanted)].copy()

# -- cross-reference: keep only markers detected in the dataset --
markers_df['gene_lower'] = markers_df['gene_symbol'].str.lower()
detected_mask = markers_df['gene_lower'].isin(proteins_df['gene_lower'])
detected_markers = markers_df[detected_mask].copy()
excluded_markers = markers_df[~detected_mask].copy()

print(f'\nMarker cross-reference:')
print(f'  Total markers in selected cell types: {len(markers_df)}')
print(f'  Detected in proteomics dataset: {len(detected_markers)}')
print(f'  Not detected (excluded): {len(excluded_markers)}')

# -- apply presentation subset filter if requested --
if SIMPLIFIED:
    subset_lower = [g.lower() for g in PRESENTATION_SUBSET]
    detected_markers = detected_markers[detected_markers['gene_lower'].isin(subset_lower)].copy()
    print(f'  Simplified mode: filtered to {len(detected_markers)} curated markers')

# -- merge with proteomics abundance data --
target_df = detected_markers.merge(
    proteins_df[[COL_GENE, COL_CONTROL, COL_KETAMINE, 'category', 'gene_lower']],
    on='gene_lower',
    how='left'
)
target_df[COL_CONTROL] = pd.to_numeric(target_df[COL_CONTROL], errors='coerce')
target_df[COL_KETAMINE] = pd.to_numeric(target_df[COL_KETAMINE], errors='coerce')

# -- display label: uppercase gene symbol --
target_df['display_label'] = target_df['gene_symbol'].str.upper()

# -- assign display group name --
ct_to_label = dict(CELL_TYPE_ORDER)
target_df['display_group'] = target_df['cell_type'].map(ct_to_label)

# -- order: by cell type group, then by order in the markers CSV --
ct_order = {ct: i for i, (ct, _) in enumerate(CELL_TYPE_ORDER)}
target_df['ct_sort'] = target_df['cell_type'].map(ct_order)
target_df = target_df.sort_values(['ct_sort', 'gene_lower']).reset_index(drop=True)

# -- identify presence/absence proteins (NaN in one condition) --
pa_mask = target_df[COL_CONTROL].isna() | target_df[COL_KETAMINE].isna()
pa_proteins = target_df[pa_mask]
quant_proteins = target_df[~pa_mask]

# ======================================================================
# 3. REPORT DETECTION RESULTS
# ======================================================================

print('\n=== Detected markers by cell type ===')
for ct, label in CELL_TYPE_ORDER:
    ct_all = markers_df[markers_df['cell_type'] == ct]
    ct_detected = target_df[target_df['cell_type'] == ct]
    ct_excluded = excluded_markers[excluded_markers['cell_type'] == ct]
    print(f'\n  {label}: {len(ct_detected)}/{len(ct_all)} detected')
    for _, row in ct_detected.iterrows():
        ctrl = row[COL_CONTROL]
        ket = row[COL_KETAMINE]
        ctrl_str = f'{ctrl:.1f}' if not np.isnan(ctrl) else 'N/D'
        ket_str = f'{ket:.1f}' if not np.isnan(ket) else 'N/D'
        pa_flag = ' [P/A]' if np.isnan(ctrl) or np.isnan(ket) else ''
        print(f'    {row["gene_symbol"]:10s}  ctrl={ctrl_str:>6s}  ket={ket_str:>6s}{pa_flag}')
    if len(ct_excluded) > 0:
        nd_list = ', '.join(ct_excluded['gene_symbol'].tolist())
        print(f'    Not detected: {nd_list}')

if len(pa_proteins) > 0:
    print(f'\nNote: {len(pa_proteins)} presence/absence protein(s) will be annotated:')
    for _, row in pa_proteins.iterrows():
        print(f'  {row["gene_symbol"]} ({row["cell_type"]}) - '
              f'{"ketamine" if np.isnan(row[COL_CONTROL]) else "control"}-specific')

# -- abundance summary for QC --
print('\nAbundance summary (quantitative proteins only):')
for ct, label in CELL_TYPE_ORDER:
    subset = quant_proteins[quant_proteins['cell_type'] == ct]
    if len(subset) > 0:
        ctrl_mean = subset[COL_CONTROL].mean()
        ket_mean = subset[COL_KETAMINE].mean()
        print(f'  {label}:')
        print(f'    Control mean: {ctrl_mean:.1f}  |  Ketamine mean: {ket_mean:.1f}')

# ======================================================================
# 4. BUILD FIGURE
# ======================================================================

print('\nCreating cell-type fidelity bar chart...')

fig = go.Figure()

# -- x-axis labels --
x_labels = target_df['display_label'].tolist()

# -- for presence/absence proteins, plot 0 instead of NaN --
ctrl_values = target_df[COL_CONTROL].fillna(0).tolist()
ket_values = target_df[COL_KETAMINE].fillna(0).tolist()

# -- control bars --
fig.add_trace(go.Bar(
    name='Control',
    x=x_labels,
    y=ctrl_values,
    marker_color=COLOR_CONTROL,
    marker_line=dict(width=0.5, color='#555555'),
    hovertemplate='%{x}<br>Control: %{y:.1f}<extra></extra>',
))

# -- ketamine bars --
fig.add_trace(go.Bar(
    name='Ketamine',
    x=x_labels,
    y=ket_values,
    marker_color=COLOR_KETAMINE,
    marker_line=dict(width=0.5, color='#555555'),
    hovertemplate='%{x}<br>Ketamine: %{y:.1f}<extra></extra>',
))

# -- annotate presence/absence proteins --
for _, row in pa_proteins.iterrows():
    gene_upper = row['gene_symbol'].upper()
    if np.isnan(row[COL_CONTROL]):
        annotation_text = 'N/D'
        annotation_x = gene_upper
        # place annotation where the control bar would be (slightly left of center)
        fig.add_annotation(
            x=gene_upper, y=5,
            text='N/D',
            showarrow=False,
            font=dict(size=8, color='#888888'),
            xshift=-12,
        )
    if np.isnan(row[COL_KETAMINE]):
        fig.add_annotation(
            x=gene_upper, y=5,
            text='N/D',
            showarrow=False,
            font=dict(size=8, color='#888888'),
            xshift=12,
        )

# ======================================================================
# 5. ADD CLUSTER SEPARATORS AND LABELS
# ======================================================================

# -- count proteins per group --
group_counts = []
for ct, label in CELL_TYPE_ORDER:
    n = len(target_df[target_df['cell_type'] == ct])
    if n > 0:
        group_counts.append((label, n, ct))

# -- vertical separator lines and group labels --
y_max = max(max(ctrl_values), max(ket_values)) * 1.05
label_y = y_max * 1.02

cumulative = 0
separator_positions = []
for i, (label, n, ct) in enumerate(group_counts):
    # add separator line before this group (except the first)
    if i > 0:
        sep_x = cumulative - 0.5
        separator_positions.append(sep_x)
        fig.add_shape(
            type='line',
            x0=sep_x, x1=sep_x,
            y0=0, y1=y_max,
            line=dict(color='#888888', width=1.5, dash='dot'),
        )

    # add group label centered over the cluster
    center_x = cumulative + (n - 1) / 2
    # color: dark blue for astrocyte, red for contaminants
    label_color = '#2C3E50' if ct == 'Astrocyte' else '#C0392B'
    fig.add_annotation(
        x=center_x,
        y=label_y,
        text=f'<b>{label}</b>',
        showarrow=False,
        font=dict(size=11, color=label_color),
        yanchor='bottom',
    )

    cumulative += n

# ======================================================================
# 6. LAYOUT
# ======================================================================

fig.update_layout(
    title=dict(
        text='Cell-Type Fidelity Assessment: Grouped Abundances by Condition',
        font=dict(size=16, color='#2C3E50'),
        x=0.5,
        xanchor='center',
    ),
    barmode='group',
    bargap=0.15,
    bargroupgap=0.05,
    xaxis=dict(
        title='',
        tickangle=-45,
        tickfont=dict(size=11 if SIMPLIFIED else 9),
    ),
    yaxis=dict(
        title=dict(text='Grouped Abundance', font=dict(size=13)),
        tickfont=dict(size=11),
    ),
    legend=dict(
        orientation='h',
        yanchor='bottom',
        y=1.08,
        xanchor='left',
        x=0.0,
        font=dict(size=12),
    ),
    plot_bgcolor='white',
    paper_bgcolor='white',
    width=900 if SIMPLIFIED else 1200,
    height=550 if SIMPLIFIED else 650,
    margin=dict(t=150 if SIMPLIFIED else 180, b=100, l=70, r=30),
)

# -- grid lines --
fig.update_yaxes(
    showgrid=True,
    gridcolor='#E8E8E8',
    gridwidth=1,
    zeroline=True,
    zerolinecolor='#CCCCCC',
    zerolinewidth=1,
)

fig.update_xaxes(showgrid=False)

# -- add footnote for presence/absence proteins --
if len(pa_proteins) > 0:
    pa_genes = ', '.join(pa_proteins['gene_symbol'].str.upper().tolist())
    fig.add_annotation(
        xref='paper', yref='paper',
        x=0, y=-0.22,
        text=f'N/D = not detected in that condition (presence/absence protein). '
             f'Applies to: {pa_genes}.',
        showarrow=False,
        font=dict(size=9, color='#666666'),
        xanchor='left',
    )

# ======================================================================
# 7. SAVE OUTPUTS
# ======================================================================

os.makedirs(OUTPUT_DIR, exist_ok=True)

suffix = '_simplified' if SIMPLIFIED else ''
png_path = os.path.join(OUTPUT_DIR, f'cell_type_fidelity_barchart{suffix}.png')
pdf_path = os.path.join(OUTPUT_DIR, f'cell_type_fidelity_barchart{suffix}.pdf')
html_path = os.path.join(OUTPUT_DIR, f'cell_type_fidelity_barchart{suffix}.html')

fig.write_image(png_path, scale=2)
fig.write_image(pdf_path, scale=2)
fig.write_html(html_path)

print(f'\nSaved: {png_path}')
print(f'Saved: {pdf_path}')
print(f'Saved: {html_path}')
print('\nDone.')

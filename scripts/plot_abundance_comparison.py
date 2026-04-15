# title: plot_abundance_comparison.py
# project: Ketamine Astrocyte Proteomics
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-02-20
# last modified: 2026-02-20
#
# purpose:
#   Compare grouped abundances between top pathway proteins (synaptic vesicle
#   cycle / SNARE machinery) and canonical astrocyte marker proteins. Generates:
#     - Figure set 1: All proteins on one chart, separated by condition
#       (control vs ketamine)
#     - Figure set 2: Control vs ketamine side-by-side, separated by protein
#       group (astrocyte markers vs pathway proteins)
#   Presence/absence proteins are annotated rather than plotted with placeholder
#   ratios.
#
# inputs:
#   - results/combined/all_proteins_categorized.csv
#
# outputs:
#   results/figures/
#   ├── abundance_comparison_control.png / .html
#   ├── abundance_comparison_ketamine.png / .html
#   ├── abundance_comparison_astrocyte_markers.png / .html
#   └── abundance_comparison_pathway_proteins.png / .html
#
# usage example:
#   python scripts/plot_abundance_comparison.py

import pandas as pd
import plotly.graph_objects as go
import os

# ======================================================================
# 1. CONFIGURATION
# ======================================================================

# -- file paths --
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INPUT_FILE = os.path.join(BASE_DIR, 'results', 'combined', 'all_proteins_categorized.csv')
OUTPUT_DIR = os.path.join(BASE_DIR, 'results', 'figures')

# -- protein groups --
PATHWAY_GENES = [
    'Stx1a', 'Vamp2', 'Snap25', 'Rab3a', 'Stx1b', 'Stxbp1',
    'Syt1', 'Camk2a', 'Sv2a', 'Syp', 'Prrt2', 'Snca',
    'Syn1', 'Syn2', 'Syngr1', 'Dnm1'
]

ASTROCYTE_MARKERS = [
    'Gfap', 'Vim', 'S100b', 'Slc1a3', 'Slc1a2', 'Glul',
    'Aldh1l1', 'Aldoc', 'Aqp4', 'Gja1', 'Fabp7', 'Atp1b2'
]

# -- display labels for astrocyte markers --
MARKER_LABELS = {
    'Gfap': 'GFAP',
    'Vim': 'VIM',
    'S100b': 'S100B',
    'Slc1a3': 'SLC1A3\n(GLAST)',
    'Slc1a2': 'SLC1A2\n(GLT-1)',
    'Glul': 'GLUL',
    'Aldh1l1': 'ALDH1L1',
    'Aldoc': 'ALDOC',
    'Aqp4': 'AQP4',
    'Gja1': 'GJA1\n(Cx43)',
    'Fabp7': 'FABP7\n(BLBP)',
    'Atp1b2': 'ATP1B2\n(ACSA-2)',
}

# -- display labels for pathway proteins --
PATHWAY_LABELS = {
    'Stx1a': 'STX1A',
    'Vamp2': 'VAMP2',
    'Snap25': 'SNAP25',
    'Rab3a': 'RAB3A',
    'Stx1b': 'STX1B',
    'Stxbp1': 'STXBP1',
    'Syt1': 'SYT1',
    'Camk2a': 'CAMK2A',
    'Sv2a': 'SV2A*',
    'Syp': 'SYP',
    'Prrt2': 'PRRT2',
    'Snca': 'SNCA*',
    'Syn1': 'SYN1',
    'Syn2': 'SYN2',
    'Syngr1': 'SYNGR1',
    'Dnm1': 'DNM1',
}

# -- colors --
COLOR_PATHWAY = '#E74C3C'           # red for pathway proteins
COLOR_ASTROCYTE = '#3498DB'         # blue for astrocyte markers
COLOR_CONTROL = '#7FB3D8'           # light blue for control
COLOR_KETAMINE = '#E8735A'          # light red for ketamine
COLOR_PRESENCE_ONLY = '#F5B041'     # amber for presence-only proteins

# ======================================================================
# 2. LOAD AND FILTER DATA
# ======================================================================

print('Loading data...')
df = pd.read_csv(INPUT_FILE)

# -- extract pathway proteins --
pathway_df = df[df['Gene Symbol'].isin(PATHWAY_GENES)].copy()
pathway_df['group'] = 'Synaptic Vesicle / SNARE'
pathway_df['display_label'] = pathway_df['Gene Symbol'].map(PATHWAY_LABELS)

# -- extract astrocyte markers --
marker_df = df[df['Gene Symbol'].isin(ASTROCYTE_MARKERS)].copy()
marker_df['group'] = 'Astrocyte Markers'
marker_df['display_label'] = marker_df['Gene Symbol'].map(MARKER_LABELS)

# -- report what was found --
print(f'Pathway proteins found: {len(pathway_df)}/{len(PATHWAY_GENES)}')
print(f'Astrocyte markers found: {len(marker_df)}/{len(ASTROCYTE_MARKERS)}')

missing_pathway = set(PATHWAY_GENES) - set(pathway_df['Gene Symbol'])
missing_markers = set(ASTROCYTE_MARKERS) - set(marker_df['Gene Symbol'])
if missing_pathway:
    print(f'  Missing pathway genes: {missing_pathway}')
if missing_markers:
    print(f'  Missing marker genes: {missing_markers}')

# -- identify presence/absence proteins --
pa_proteins = pathway_df[pathway_df['category'].str.contains('presence_absence', na=False)]
pa_genes = set(pa_proteins['Gene Symbol'].tolist())
if len(pa_proteins) > 0:
    print(f'\nNote: {len(pa_proteins)} pathway protein(s) are presence/absence '
          f'(detected only in ketamine):')
    for _, row in pa_proteins.iterrows():
        print(f'  {row["Gene Symbol"]} - ketamine-specific')

# -- order proteins by input list order --
pathway_order = [g for g in PATHWAY_GENES if g in pathway_df['Gene Symbol'].values]
marker_order = [g for g in ASTROCYTE_MARKERS if g in marker_df['Gene Symbol'].values]

combined = pd.concat([pathway_df, marker_df], ignore_index=True)
gene_order = pathway_order + marker_order
combined['sort_key'] = combined['Gene Symbol'].map({g: i for i, g in enumerate(gene_order)})
combined = combined.sort_values('sort_key').reset_index(drop=True)

# -- presence/absence annotation text --
PA_NOTE = ('* Presence/absence protein: detected only in ketamine samples. '
           'No quantitative abundance in control; ketamine value is a '
           'placeholder (not a real measurement).')


# ======================================================================
# 3. FIGURE SET 1: BY CONDITION (all proteins, one condition per figure)
# ======================================================================

def create_by_condition_figure(condition, abundance_col, title_suffix):
    '''
    Create a bar chart showing all proteins for one condition (control or
    ketamine), with pathway proteins and astrocyte markers colored differently.
    Presence/absence proteins are annotated instead of using placeholder values.
    '''

    fig = go.Figure()

    # ---- pathway proteins ----
    pw_data = combined[combined['group'] == 'Synaptic Vesicle / SNARE'].copy()
    pw_abundances = pd.to_numeric(pw_data[abundance_col], errors='coerce')
    pw_labels = pw_data['display_label'].tolist()
    pw_genes = pw_data['Gene Symbol'].tolist()
    pw_categories = pw_data['category'].tolist()

    pw_colors = []
    pw_text = []
    pw_y_vals = []
    for gene, val, cat in zip(pw_genes, pw_abundances, pw_categories):
        is_pa = 'presence_absence' in str(cat)
        if is_pa and condition == 'control':
            # -- no data in control for presence/absence proteins --
            pw_colors.append(COLOR_PRESENCE_ONLY)
            pw_text.append('N/D')
            pw_y_vals.append(0)
        elif is_pa and condition == 'ketamine':
            # -- ketamine value is a placeholder, annotate it --
            pw_colors.append(COLOR_PRESENCE_ONLY)
            pw_text.append('Ket. only')
            pw_y_vals.append(0)
        else:
            pw_colors.append(COLOR_PATHWAY)
            pw_text.append(f'{val:.1f}')
            pw_y_vals.append(val if not pd.isna(val) else 0)

    fig.add_trace(go.Bar(
        x=pw_labels,
        y=pw_y_vals,
        name='Synaptic Vesicle / SNARE',
        marker_color=pw_colors,
        text=pw_text,
        textposition='outside',
        textfont=dict(size=9),
        hovertemplate='<b>%{x}</b><br>Abundance: %{text}<extra>Pathway Protein</extra>',
    ))

    # ---- astrocyte markers ----
    am_data = combined[combined['group'] == 'Astrocyte Markers'].copy()
    am_abundances = pd.to_numeric(am_data[abundance_col], errors='coerce')
    am_labels = am_data['display_label'].tolist()

    fig.add_trace(go.Bar(
        x=am_labels,
        y=am_abundances,
        name='Astrocyte Markers',
        marker_color=COLOR_ASTROCYTE,
        text=[f'{v:.1f}' for v in am_abundances],
        textposition='outside',
        textfont=dict(size=9),
        hovertemplate='<b>%{x}</b><br>Abundance: %{y:.1f}<extra>Astrocyte Marker</extra>',
    ))

    # ---- divider line ----
    n_pathway = len(pw_labels)
    fig.add_vline(x=n_pathway - 0.5, line_dash='dash',
                  line_color='grey', line_width=1.5)

    # ---- group labels ----
    fig.add_annotation(
        x=(n_pathway - 1) / 2, y=1.08, xref='x', yref='paper',
        text='<b>Synaptic Vesicle / SNARE Proteins</b>',
        showarrow=False, font=dict(size=12, color=COLOR_PATHWAY),
    )
    n_markers = len(am_labels)
    fig.add_annotation(
        x=n_pathway + (n_markers - 1) / 2, y=1.08, xref='x', yref='paper',
        text='<b>Astrocyte Marker Proteins</b>',
        showarrow=False, font=dict(size=12, color=COLOR_ASTROCYTE),
    )

    # ---- presence/absence footnote ----
    fig.add_annotation(
        x=0.01, y=-0.22, xref='paper', yref='paper',
        text=PA_NOTE,
        showarrow=False, font=dict(size=9, color='grey'), align='left',
    )

    # ---- layout ----
    fig.update_layout(
        title=dict(
            text=f'Grouped Protein Abundance — {title_suffix}',
            font=dict(size=16), x=0.5,
        ),
        yaxis_title='Grouped Abundance (scaled)',
        xaxis=dict(tickangle=-45, tickfont=dict(size=10)),
        yaxis=dict(rangemode='tozero'),
        plot_bgcolor='white', paper_bgcolor='white',
        showlegend=False,
        margin=dict(b=140, t=100, l=80, r=40),
        width=1200, height=600,
    )
    fig.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='lightgrey')

    return fig


# ======================================================================
# 4. FIGURE SET 2: BY PROTEIN GROUP (control vs ketamine side-by-side)
# ======================================================================

def create_by_group_figure(group_name, gene_list, label_map, title, color_accent):
    '''
    Create a grouped bar chart showing control vs ketamine abundances
    side-by-side for one protein group. Presence/absence proteins are
    annotated instead of using placeholder values.
    '''

    group_df = combined[combined['Gene Symbol'].isin(gene_list)].copy()
    group_df['sort_key'] = group_df['Gene Symbol'].map(
        {g: i for i, g in enumerate(gene_list)})
    group_df = group_df.sort_values('sort_key').reset_index(drop=True)

    labels = group_df['display_label'].tolist()
    genes = group_df['Gene Symbol'].tolist()
    categories = group_df['category'].tolist()

    ctrl_vals = pd.to_numeric(
        group_df['Abundances (Grouped): control'], errors='coerce')
    ket_vals = pd.to_numeric(
        group_df['Abundances (Grouped): ketamine'], errors='coerce')

    # -- prepare control trace --
    ctrl_y = []
    ctrl_text = []
    ctrl_colors = []
    for gene, val, cat in zip(genes, ctrl_vals, categories):
        is_pa = 'presence_absence' in str(cat)
        if is_pa:
            ctrl_y.append(0)
            ctrl_text.append('N/D')
            ctrl_colors.append(COLOR_PRESENCE_ONLY)
        else:
            ctrl_y.append(val if not pd.isna(val) else 0)
            ctrl_text.append(f'{val:.1f}' if not pd.isna(val) else 'N/D')
            ctrl_colors.append(COLOR_CONTROL)

    # -- prepare ketamine trace --
    ket_y = []
    ket_text = []
    ket_colors = []
    for gene, val, cat in zip(genes, ket_vals, categories):
        is_pa = 'presence_absence' in str(cat)
        if is_pa:
            ket_y.append(0)
            ket_text.append('Ket. only*')
            ket_colors.append(COLOR_PRESENCE_ONLY)
        else:
            ket_y.append(val if not pd.isna(val) else 0)
            ket_text.append(f'{val:.1f}' if not pd.isna(val) else 'N/D')
            ket_colors.append(COLOR_KETAMINE)

    fig = go.Figure()

    # -- control bars --
    fig.add_trace(go.Bar(
        x=labels,
        y=ctrl_y,
        name='Control (Saline)',
        marker_color=ctrl_colors,
        text=ctrl_text,
        textposition='outside',
        textfont=dict(size=9),
        hovertemplate='<b>%{x}</b><br>Control: %{text}<extra></extra>',
    ))

    # -- ketamine bars --
    fig.add_trace(go.Bar(
        x=labels,
        y=ket_y,
        name='Ketamine (10 mg/kg)',
        marker_color=ket_colors,
        text=ket_text,
        textposition='outside',
        textfont=dict(size=9),
        hovertemplate='<b>%{x}</b><br>Ketamine: %{text}<extra></extra>',
    ))

    # -- presence/absence footnote if needed --
    has_pa = any('presence_absence' in str(c) for c in categories)
    if has_pa:
        fig.add_annotation(
            x=0.01, y=-0.22, xref='paper', yref='paper',
            text=PA_NOTE,
            showarrow=False, font=dict(size=9, color='grey'), align='left',
        )

    # ---- layout ----
    fig.update_layout(
        title=dict(
            text=title,
            font=dict(size=16), x=0.5,
        ),
        yaxis_title='Grouped Abundance (scaled)',
        xaxis=dict(tickangle=-45, tickfont=dict(size=10)),
        yaxis=dict(rangemode='tozero'),
        barmode='group',
        plot_bgcolor='white', paper_bgcolor='white',
        legend=dict(
            orientation='h',
            yanchor='bottom', y=1.02,
            xanchor='center', x=0.5,
            font=dict(size=12),
        ),
        margin=dict(b=140, t=100, l=80, r=40),
        width=1200, height=600,
    )
    fig.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='lightgrey')

    return fig


# ======================================================================
# 5. GENERATE ALL FIGURES
# ======================================================================

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ---- Figure Set 1: by condition ----
print('\n--- Figure Set 1: By Condition ---')

print('Generating control abundance figure...')
fig_ctrl = create_by_condition_figure(
    condition='control',
    abundance_col='Abundances (Grouped): control',
    title_suffix='Control (Saline)',
)
ctrl_png = os.path.join(OUTPUT_DIR, 'abundance_comparison_control.png')
ctrl_html = os.path.join(OUTPUT_DIR, 'abundance_comparison_control.html')
fig_ctrl.write_image(ctrl_png, scale=3)
fig_ctrl.write_html(ctrl_html)
print(f'  Saved: {ctrl_png}')
print(f'  Saved: {ctrl_html}')

print('Generating ketamine abundance figure...')
fig_ket = create_by_condition_figure(
    condition='ketamine',
    abundance_col='Abundances (Grouped): ketamine',
    title_suffix='Ketamine (10 mg/kg)',
)
ket_png = os.path.join(OUTPUT_DIR, 'abundance_comparison_ketamine.png')
ket_html = os.path.join(OUTPUT_DIR, 'abundance_comparison_ketamine.html')
fig_ket.write_image(ket_png, scale=3)
fig_ket.write_html(ket_html)
print(f'  Saved: {ket_png}')
print(f'  Saved: {ket_html}')

# ---- Figure Set 2: by protein group (control vs ketamine) ----
print('\n--- Figure Set 2: By Protein Group (Control vs Ketamine) ---')

print('Generating astrocyte markers comparison figure...')
fig_markers = create_by_group_figure(
    group_name='Astrocyte Markers',
    gene_list=marker_order,
    label_map=MARKER_LABELS,
    title='Astrocyte Marker Proteins — Control vs Ketamine',
    color_accent=COLOR_ASTROCYTE,
)
markers_png = os.path.join(OUTPUT_DIR, 'abundance_comparison_astrocyte_markers.png')
markers_html = os.path.join(OUTPUT_DIR, 'abundance_comparison_astrocyte_markers.html')
fig_markers.write_image(markers_png, scale=3)
fig_markers.write_html(markers_html)
print(f'  Saved: {markers_png}')
print(f'  Saved: {markers_html}')

print('Generating pathway proteins comparison figure...')
fig_pathway = create_by_group_figure(
    group_name='Synaptic Vesicle / SNARE',
    gene_list=pathway_order,
    label_map=PATHWAY_LABELS,
    title='Synaptic Vesicle / SNARE Proteins — Control vs Ketamine',
    color_accent=COLOR_PATHWAY,
)
pathway_png = os.path.join(OUTPUT_DIR, 'abundance_comparison_pathway_proteins.png')
pathway_html = os.path.join(OUTPUT_DIR, 'abundance_comparison_pathway_proteins.html')
fig_pathway.write_image(pathway_png, scale=3)
fig_pathway.write_html(pathway_html)
print(f'  Saved: {pathway_png}')
print(f'  Saved: {pathway_html}')

print('\nDone. All 8 figures generated.')

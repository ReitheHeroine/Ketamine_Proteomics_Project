# title: plot_cell_type_fidelity.py
# project: Ketamine Astrocyte Proteomics
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-05-26  (replaces prior compositional grouped-bar version
#                            from 2026-03-17; see commit history)
# last modified: 2026-05-27  (style guide compliance pass: Arial fonts,
#                             section-3 type sizes, log subscript notation,
#                             directional P/A colors, removal of on-figure
#                             title and footnote per section 11.2; see
#                             project_notes/figure_and_table_style_guide.Rmd.
#                             Same-day revisions:
#                             - slide-variant text sizes (~1.25x manuscript
#                               spec, section 3 table 2);
#                             - split into one figure per cluster (the
#                               combined version was too tall to read);
#                             - shared x-axis range across the three figures,
#                               excluding P/A markers from the range;
#                               significance asterisks removed (legend defines
#                               significance instead);
#                             - zero line changed to dotted; reference lines
#                               stop below the cluster header so they do not
#                               cross the header text;
#                             - P/A and N/D row labels placed next to the zero
#                               line on their row (label text replaces marker
#                               point for those rows) rather than at the
#                               axis extreme;
#                             - cluster labels prefixed with "Canonical".)
#
# purpose:
#   Forest plot of log2 fold change (ketamine / control) for cell-type marker
#   proteins, organized into Astrocyte / Neuronal / Oligodendrocyte clusters.
#   Replaces the previous "Cell-Type Fidelity Assessment" grouped-bar chart,
#   which plotted Proteome Discoverer's `Abundances (Grouped):` columns. Those
#   values are per-protein rescaled so that each protein's two condition means
#   sum to 200 (Scaling Mode: On All Average), which makes them compositional
#   and not interpretable as absolute or even consistent relative abundance.
#   This version instead reads PD's `Abundance Ratio`, `Abundance Ratio Adj.
#   P-Value`, and `Abundance Ratio Variability [%]` columns, which PD computes
#   from peptide-level intensities via pairwise ratios and which are NOT
#   constrained to a constant per-protein sum.
#
#   The figure communicates three things per marker:
#     - signed log2 fold change between conditions (point position)
#     - whether the shift is statistically significant at adj.P < 0.05
#       (filled vs open point; significance-level stars)
#     - PD's pairwise-ratio variability, mapped onto the log2 axis as
#       approximate multiplicative spread (horizontal segment around the
#       point). This is a visualization aid, not a confidence interval.
#
#   Markers in the reference panel that were not detected in the proteomics
#   dataset appear as N/D rows. Presence/absence markers (detected in only
#   one condition; PD assigns ratio = 100.0 or 0.01) are plotted at the PD
#   default position and annotated as P/A.
#
# inputs:
#   - data/cell_type_markers.csv        (marker definitions by cell type)
#   - data/all_proteins_categorized.csv (master proteomics table from
#                                        diff_abundance_analysis.py output)
#
# outputs:
#   results/figures/
#   |- cell_type_marker_forest_astrocyte.{png,pdf,html}
#   |- cell_type_marker_forest_neuron.{png,pdf,html}
#   `- cell_type_marker_forest_oligodendrocyte.{png,pdf,html}
#   (with `_simplified` suffix appended when --simplified is set, e.g.,
#    cell_type_marker_forest_astrocyte_simplified.png)
#
# usage example:
#   python scripts/plot_cell_type_fidelity.py              # full panel
#   python scripts/plot_cell_type_fidelity.py --simplified # curated 12-marker subset

import argparse
import os
import sys

import numpy as np
import pandas as pd
import plotly.graph_objects as go


# ============================================================================
# 1. CONFIGURATION
# ============================================================================

# --- file paths ---
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MARKERS_FILE = os.path.join(BASE_DIR, "data", "cell_type_markers.csv")
PROTEINS_FILE = os.path.join(BASE_DIR, "data", "all_proteins_categorized.csv")
OUTPUT_DIR = os.path.join(BASE_DIR, "results", "figures")

# --- column names in the proteomics master file ---
COL_GENE = "Gene Symbol"
COL_RATIO = "Abundance Ratio: (ketamine) / (control)"
COL_PVAL = "Abundance Ratio Adj. P-Value: (ketamine) / (control)"
COL_VAR = "Abundance Ratio Variability [%]: (ketamine) / (control)"
COL_CATEGORY = "category"

# --- cluster display order ---
# Each cluster is now rendered as its own figure (top-down: Astrocyte,
# Neuronal, Oligodendrocyte). Tuple = (cell_type value in markers CSV,
# display label, accent color, output-filename slug).
# Accent colors are used for both the cluster header text and the significant
# marker dots within that cluster, so each panel reads as a matched set.
# NOTE: the astrocyte accent here (#1E8449, Flat-UI "Nephritis" forest green)
# overrides the style guide section 9.1.2 default of #2C3E50 dark slate; the
# green pairs better with the neuronal red and oligodendrocyte purple as a
# three-color set. Dark slate remains the neutral text/outline color
# (COLOR_POINT below) and the open-circle outline for not-significant points.
CLUSTER_ORDER = [
    ("Astrocyte",       "Canonical astrocyte markers",       "#1E8449", "astrocyte"),
    ("Neuron",          "Canonical neuronal markers",        "#C0392B", "neuron"),
    ("Oligodendrocyte", "Canonical oligodendrocyte markers", "#7D3C98", "oligodendrocyte"),
]

# --- curated 12-marker subset for --simplified mode ---
PRESENTATION_SUBSET = [
    "Gfap", "Aqp4", "Aldh1l1", "Slc1a3", "Glul", "Gja1",  # astrocyte
    "Ina", "Nefl", "Nefm",                                  # neuron
    "Plp1", "Mbp", "Mog",                                   # oligo
]

# --- significance threshold (matches diff_abundance_analysis.py default) ---
PVAL_SIG = 0.05

# --- styling (palette per figure_and_table_style_guide.Rmd section 9) ---
BLACK = "#000000"              # pure black; all in-figure text (except the
                               # semantic-colored row labels listed below) and
                               # the x / y axis lines. Project-wide default
                               # was dark slate (#2C3E50, section 9.1.4); this
                               # figure overrides to pure black for higher
                               # print contrast on the per-cluster panels.
COLOR_POINT = "#2C3E50"        # dark slate; retained for the open-circle
                               # outline on not-significant markers (softer
                               # than pure black so the colored significant
                               # dots remain the visual focus).
COLOR_ERR_BAR = "#B0B0B0"      # mid gray; error bars and secondary annotations
COLOR_REF_ZERO = "#666666"     # darker gray; zero / threshold reference line
COLOR_REF_GRID = "#E0E0E0"     # light gray; faint integer-log2FC reference lines
COLOR_ND_TEXT = "#999999"      # gray; N/D row labels and "not detected" text
COLOR_SEPARATOR = "#D0D0D0"    # mid-light gray; cluster separator dotted lines

# P/A markers point at the extreme right (ketamine-only) or extreme left
# (control-only) of the figure. Reuse the directional treatment-group colors
# per section 9.1.1: "up-in-ketamine" -> ketamine color, "down-in-ketamine"
# -> control color. This is the same direction-as-color convention used in
# scripts/plot_top_significant_proteins.py and scripts/visualize_diff_abundance.py.
COLOR_PA_KET = "#E8735A"       # coral; P/A ketamine-only (up-in-ketamine)
COLOR_PA_CTRL = "#7FB3D8"      # light blue; P/A control-only (down-in-ketamine)

# PD assigns ratio = 100.0 (or 0.01) for presence/absence proteins. With
# log2(100) ~= 6.64, P/A markers sit at the far right of the plot; annotate.
PA_RATIO_HIGH = 100.0
PA_RATIO_LOW = 0.01


# ============================================================================
# 2. ARGUMENT PARSING
# ============================================================================

parser = argparse.ArgumentParser(
    description="Forest plot of cell-type marker log2 fold change "
                "(ketamine / control) using PD pairwise-ratio data."
)
parser.add_argument(
    "--simplified",
    action="store_true",
    help="Restrict to the curated 12-marker canonical subset.",
)
args = parser.parse_args()
SIMPLIFIED = args.simplified


# ============================================================================
# 3. LOAD AND CROSS-REFERENCE
# ============================================================================

print("Loading marker panel and proteomics data...")
markers = pd.read_csv(MARKERS_FILE)
proteins = pd.read_csv(PROTEINS_FILE)

# --- normalize gene-symbol case for matching ---
markers["gene_lower"] = markers["gene_symbol"].astype(str).str.lower()
proteins["gene_lower"] = proteins[COL_GENE].astype(str).str.lower()

# --- keep only the three clusters of interest, in the configured order ---
cluster_names = [c[0] for c in CLUSTER_ORDER]
markers = markers[markers["cell_type"].isin(cluster_names)].copy()

if SIMPLIFIED:
    subset_lower = [g.lower() for g in PRESENTATION_SUBSET]
    markers = markers[markers["gene_lower"].isin(subset_lower)].copy()
    print(f"  --simplified: panel reduced to {len(markers)} curated markers")

# --- merge marker panel with proteomics ratio columns ---
df = markers.merge(
    proteins[[COL_GENE, COL_RATIO, COL_PVAL, COL_VAR, COL_CATEGORY, "gene_lower"]],
    on="gene_lower",
    how="left",
)

# --- coerce ratio columns to numeric and derive log2FC ---
df[COL_RATIO] = pd.to_numeric(df[COL_RATIO], errors="coerce")
df[COL_PVAL] = pd.to_numeric(df[COL_PVAL], errors="coerce")
df[COL_VAR] = pd.to_numeric(df[COL_VAR], errors="coerce")
df["log2FC"] = np.log2(df[COL_RATIO])

# --- flags ---
df["detected"] = df[COL_RATIO].notna()
df["is_pa_ket"] = df[COL_RATIO] == PA_RATIO_HIGH
df["is_pa_ctrl"] = df[COL_RATIO] == PA_RATIO_LOW
df["is_pa"] = df["is_pa_ket"] | df["is_pa_ctrl"]
df["is_quantitative"] = df["detected"] & ~df["is_pa"]
df["sig"] = df["detected"] & (df[COL_PVAL] < PVAL_SIG)

# --- map PD pairwise variability % to approximate half-width on log2 scale ---
# Var% behaves like a coefficient of variation on the multiplicative ratio.
# log2(1 + Var/100) approximates one-sigma spread on the log2 axis.
df["err_log2"] = np.log2(1.0 + df[COL_VAR].fillna(0.0) / 100.0)

print(f"\nMarker panel size: {len(df)}")
print(f"  detected:        {df['detected'].sum()}")
print(f"  quantitative:    {df['is_quantitative'].sum()}")
print(f"  presence/abs:    {df['is_pa'].sum()}")
print(f"  not detected:    {(~df['detected']).sum()}")


# ============================================================================
# 4. ROW SORTING (within-cluster only; y positions are assigned per-figure)
# ============================================================================

# Within each cluster: sort detected markers by log2FC ascending (so the eye
# scans from "no shift" at top of cluster to "largest shift" at bottom). N/D
# markers go to the bottom of their cluster.
def within_cluster_sort_key(row: pd.Series) -> float:
    if not row["detected"]:
        return float("inf")
    return float(row["log2FC"])


df["_within_sort"] = df.apply(within_cluster_sort_key, axis=1)
df["_cluster_idx"] = df["cell_type"].map({c[0]: i for i, c in enumerate(CLUSTER_ORDER)})
df = df.sort_values(["_cluster_idx", "_within_sort"]).reset_index(drop=True)


# ============================================================================
# 5. SHARED X-AXIS RANGE (excludes P/A markers)
# ============================================================================

# Shared across all three per-cluster figures so effect sizes are directly
# comparable across panels. P/A markers (PD ratio = 100 or 0.01, log2FC = ±6.64)
# are excluded from the range because they are qualitative placeholders, not
# quantitative measurements; they are rendered as row labels next to the zero
# line instead of as points (see build_cluster_figure below).
quant_lfc = df.loc[df["is_quantitative"], "log2FC"]
if len(quant_lfc) > 0:
    XMIN = float(min(quant_lfc.min() - 1.0, -1.5))
    XMAX = float(max(quant_lfc.max() + 1.0, 1.5))
else:
    XMIN, XMAX = -3.0, 3.0


# ============================================================================
# 6. BUILD ONE FIGURE PER CLUSTER
# ============================================================================

# Layout constants shared across all per-cluster figures
HEADER_Y = 0.2         # y position of cluster header (data coord). Slightly
                       # above the first marker row (y=1.0) so the header has
                       # clear space above the protein-name column.
HEADER_X_PAPER = -0.10 # x position of cluster header (paper coord). Negative
                       # paper x places the header in the left margin area
                       # where the protein-name y-axis tick labels live, so
                       # the header reads as a section title for that column
                       # rather than as a label inside the plot area.
LINE_TOP_Y = 0.7       # reference lines start here (below the header)
ROW_HEIGHT = 1.0       # y spacing between consecutive marker rows
LABEL_X_OFFSET = 0.15  # horizontal offset of P/A and N/D labels from zero
ROW_FIRST_Y = 1.0      # y position of the first marker row


def build_cluster_figure(cluster_name: str,
                         cluster_label: str,
                         cluster_color: str) -> go.Figure | None:
    """Build a forest plot figure for a single cluster.

    The figure follows style guide sections 3, 9, and 11.2:
    - No on-figure title; caption content lives in the Word document.
    - Arial fonts at slide-variant sizes (section 3 table 2).
    - Reference lines start below the cluster header so they do not cross
      the header text (the "Canonical X markers" subheading).
    - Zero line is dotted.
    - P/A and N/D rows are labeled at the zero line in their row rather
      than at the axis extreme; no marker point is drawn for those rows.
    """
    sub = df[df["cell_type"] == cluster_name].copy().reset_index(drop=True)
    if len(sub) == 0:
        return None

    # --- y positions: header at top, marker rows below ---
    sub["y_pos"] = [ROW_FIRST_Y + ROW_HEIGHT * i for i in range(len(sub))]
    y_max = float(sub["y_pos"].max() + 0.5)

    fig = go.Figure()

    # --- faint reference lines at integer log2FC (skip zero) ---
    # Start at LINE_TOP_Y (below the header) so the lines do not cross the
    # "Canonical X markers" header text.
    for x in range(int(np.floor(XMIN)) - 1, int(np.ceil(XMAX)) + 2):
        if x == 0 or x < XMIN or x > XMAX:
            continue
        fig.add_shape(
            type="line",
            x0=x, x1=x, y0=LINE_TOP_Y, y1=y_max,
            line=dict(color=COLOR_REF_GRID, width=1, dash="dot"),
            layer="below",
        )

    # --- dotted zero line (also starts below the header) ---
    fig.add_shape(
        type="line",
        x0=0, x1=0, y0=LINE_TOP_Y, y1=y_max,
        line=dict(color=COLOR_REF_ZERO, width=1.5, dash="dot"),
        layer="below",
    )

    # Subsets used for the trace blocks below.
    quant_sub = sub[sub["is_quantitative"]].copy()
    sig_sub = quant_sub[quant_sub["sig"]].copy()
    ns_sub = quant_sub[~quant_sub["sig"]].copy()

    # --- variability segment (invisible markers, used only for error_x) ---
    if len(quant_sub) > 0:
        fig.add_trace(go.Scatter(
            x=quant_sub["log2FC"],
            y=quant_sub["y_pos"],
            error_x=dict(
                type="data",
                array=quant_sub["err_log2"],
                thickness=1.9,
                width=0,
                color=COLOR_ERR_BAR,
            ),
            mode="markers",
            marker=dict(size=0.01, color="rgba(0,0,0,0)"),
            showlegend=False,
            hoverinfo="skip",
        ))

    # Both legend entries (Significant + Not significant) are emitted in
    # every per-cluster figure so the three panels share an identical legend
    # regardless of how their markers happen to split. When a category is
    # empty for this cluster, the trace is added with a single sentinel
    # [None] coordinate: plotly registers the legend entry but does not
    # render a marker (None values are skipped at draw time).
    HOVER_TEMPLATE = (
        "<b>%{customdata[0]}</b><br>"
        "log<sub>2</sub> FC = %{x:.2f}<br>"
        "Fold change = %{customdata[2]:.1f}x<br>"
        "Adj. <i>p</i> = %{customdata[1]:.2e}<br>"
        "Var%% = %{customdata[3]:.1f}<extra></extra>"
    )

    def _trace_data(s):
        """Return (x, y, customdata) for a marker trace, with a sentinel
        [None] row when the subset is empty so the legend entry still
        appears."""
        if len(s) > 0:
            return (s["log2FC"], s["y_pos"],
                    s[["gene_symbol", COL_PVAL, COL_RATIO, COL_VAR]].values)
        return [None], [None], [[None, None, None, None]]

    sig_x, sig_y, sig_cd = _trace_data(sig_sub)
    ns_x, ns_y, ns_cd = _trace_data(ns_sub)

    # --- significant points (filled, in the cluster's accent color) ---
    fig.add_trace(go.Scatter(
        x=sig_x, y=sig_y,
        mode="markers",
        marker=dict(
            size=14,
            color=cluster_color,
            symbol="circle",
            line=dict(width=1, color=cluster_color),
        ),
        name="Significant",
        showlegend=True,
        customdata=sig_cd,
        hovertemplate=HOVER_TEMPLATE,
    ))

    # --- non-significant points (open, neutral dark-slate outline) ---
    fig.add_trace(go.Scatter(
        x=ns_x, y=ns_y,
        mode="markers",
        marker=dict(
            size=14,
            color="white",
            symbol="circle",
            line=dict(width=1.5, color=COLOR_POINT),
        ),
        name="Not significant",
        showlegend=True,
        customdata=ns_cd,
        hovertemplate=HOVER_TEMPLATE,
    ))

    # --- P/A ketamine-only labels: just right of the zero line, no marker drawn ---
    for _, row in sub[sub["is_pa_ket"]].iterrows():
        fig.add_annotation(
            x=LABEL_X_OFFSET, y=row["y_pos"],
            text="<b>Ketamine-only</b>",
            showarrow=False,
            font=dict(size=11, color=COLOR_PA_KET),
            xanchor="left",
            yanchor="middle",
        )

    # --- P/A control-only labels: just left of the zero line, no marker drawn ---
    for _, row in sub[sub["is_pa_ctrl"]].iterrows():
        fig.add_annotation(
            x=-LABEL_X_OFFSET, y=row["y_pos"],
            text="<b>Control-only</b>",
            showarrow=False,
            font=dict(size=11, color=COLOR_PA_CTRL),
            xanchor="right",
            yanchor="middle",
        )

    # --- "not detected" labels: just right of the zero line (so text does not ---
    # --- overlap the zero line, which is also dotted gray) ---
    for _, row in sub[~sub["detected"]].iterrows():
        fig.add_annotation(
            x=LABEL_X_OFFSET, y=row["y_pos"],
            text="<i>not detected</i>",
            showarrow=False,
            font=dict(size=11, color=COLOR_ND_TEXT),
            xanchor="left",
            yanchor="middle",
        )

    # --- cluster header above the y-axis tick-label column ---
    # x in paper coords (independent of XMIN); y in data coords so the header
    # tracks the first marker row at y=ROW_FIRST_Y regardless of x-axis range.
    # Header text is black so that the colored signal in the figure comes
    # only from the significant marker dots (in cluster_color) and the
    # row-status labels (P/A coral/blue, N/D gray).
    fig.add_annotation(
        xref="paper",
        yref="y",
        x=HEADER_X_PAPER, y=HEADER_Y,
        text=f"<b>{cluster_label}</b>",
        showarrow=False,
        font=dict(size=14, color=BLACK),
        xanchor="left",
        yanchor="middle",
    )

    # --- custom y-axis line (drawn as a shape, not via yaxis.showline) ---
    # The default plotly y-axis line spans the full plot area, which
    # visually competes with the cluster header at the top. This custom
    # line stops at LINE_TOP_Y (between the header and the first marker
    # row) and reaches down to the bottom of the y range so it still meets
    # the x-axis line at the corner.
    fig.add_shape(
        type="line",
        xref="paper", yref="y",
        x0=0, x1=0,
        y0=LINE_TOP_Y, y1=y_max + 0.3,
        line=dict(color=BLACK, width=2),
        layer="below",
    )

    # --- y-axis tick labels (uppercase gene symbols; roman per section 4.3) ---
    y_tick_vals = sub["y_pos"].tolist()
    y_tick_text = [g.upper() for g in sub["gene_symbol"]]

    # --- layout ---
    # Per style guide section 11.2, NO on-figure title is set; the title
    # sentence, significance ladder, and abbreviation definitions live in the
    # Word document caption block (and in cell_type_fidelity_for_llm.md), NOT
    # inside the figure.
    #
    # Legend y is computed in absolute pixels and converted to plotly's paper
    # coordinates (which span the plot area, not the full figure). This keeps
    # the legend at a fixed pixel distance below the x-axis title regardless
    # of how short the panel is - critical for the simplified variant where
    # only 3 markers / cluster gives a small plot area.
    TOP_MARGIN, BOTTOM_MARGIN = 30, 130
    ROW_PIXELS = 38
    total_h = max(380, ROW_PIXELS * len(sub) + 180)
    plot_h = total_h - TOP_MARGIN - BOTTOM_MARGIN
    legend_y_paper = -70.0 / plot_h   # ~70 px below plot area, clearing the axis title

    fig.update_layout(
        font=dict(family="Arial", size=14, color=BLACK),
        xaxis=dict(
            title=dict(
                text="Log<sub>2</sub> fold change (ketamine / control)",
                font=dict(size=14),
            ),
            tickfont=dict(size=12),
            range=[XMIN, XMAX],
            showgrid=False,
            zeroline=False,
            showline=True,
            linewidth=2,
            linecolor=BLACK,
        ),
        yaxis=dict(
            tickmode="array",
            tickvals=y_tick_vals,
            ticktext=y_tick_text,
            tickfont=dict(size=12),
            range=[y_max + 0.3, -0.1],   # reversed so header is at top
            showgrid=False,
            zeroline=False,
            # Built-in y-axis line disabled: the custom shape drawn above in
            # build_cluster_figure provides a y-axis line that stops below
            # the cluster header rather than spanning the full plot area.
            showline=False,
        ),
        plot_bgcolor="white",
        paper_bgcolor="white",
        width=900,
        height=total_h,
        margin=dict(t=TOP_MARGIN, b=BOTTOM_MARGIN, l=140, r=100),
        showlegend=True,
        legend=dict(
            orientation="h",
            yanchor="top",
            y=legend_y_paper,
            xanchor="center",
            x=0.5,
            font=dict(size=12),
        ),
    )

    return fig


# ============================================================================
# 7. SAVE (one PNG/PDF/HTML triplet per cluster)
# ============================================================================

os.makedirs(OUTPUT_DIR, exist_ok=True)
suffix = "_simplified" if SIMPLIFIED else ""

print("\nBuilding and writing per-cluster forest plots...")
for cluster_name, cluster_label, cluster_color, slug in CLUSTER_ORDER:
    fig = build_cluster_figure(cluster_name, cluster_label, cluster_color)
    if fig is None:
        print(f"  {cluster_name}: no markers in panel; skipping.")
        continue
    base = os.path.join(OUTPUT_DIR, f"cell_type_marker_forest_{slug}{suffix}")
    png_path, pdf_path, html_path = f"{base}.png", f"{base}.pdf", f"{base}.html"
    try:
        fig.write_image(png_path, scale=2)
        print(f"  {png_path}")
    except Exception as exc:
        print(f"  PNG render failed for {slug} ({exc}); install/upgrade 'kaleido' if needed.",
              file=sys.stderr)
    try:
        fig.write_image(pdf_path, scale=2)
        print(f"  {pdf_path}")
    except Exception as exc:
        print(f"  PDF render failed for {slug} ({exc}); install/upgrade 'kaleido' if needed.",
              file=sys.stderr)
    fig.write_html(html_path)
    print(f"  {html_path}")


# ============================================================================
# 9. CONSOLE SUMMARY
# ============================================================================

print("\nCluster-level summary (quantitative markers only):")
for cluster_name, cluster_label, _, _ in CLUSTER_ORDER:
    sub = df[(df["cell_type"] == cluster_name) & df["is_quantitative"]]
    if len(sub) == 0:
        continue
    med = sub["log2FC"].median()
    nsig = int(sub["sig"].sum())
    print(
        f"  {cluster_label:36s}  n={len(sub):2d}  "
        f"median log2FC = {med:+.2f}  (FC = {2**med:.2f}x)  "
        f"significant: {nsig}/{len(sub)} at adj.P < {PVAL_SIG}"
    )

print("\nDone.")

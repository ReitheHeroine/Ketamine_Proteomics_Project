#!/usr/bin/env python3
"""
================================================================================
plot_go_bp_clusters_thesis.py

Title:         GO Biological Process clustered enrichment figure (thesis)
Project:       Ketamine Astrocyte Proteomics
Author:        Reina Hastings (reinahastings13@gmail.com)
Date created:  2026-05-28
Last modified: 2026-05-28

Purpose:
    Re-plot the GO:BP cluster figure produced by pathway_analysis.py so that
    the y-axis ordering is interpretable. Clusters are ranked by the minimum
    FDR among their member terms (most significant cluster on top), and
    within each cluster terms are ranked by FDR (most significant term on
    top of that cluster's band). Bars are colored by cluster ID using a
    curated 15-color categorical palette. Style compliance follows
    project_notes/figure_and_table_style_guide.Rmd: Arial, black text, no
    on-figure title, sentence-case axis label, 14/12/10 pt hierarchy.

Inputs:
    --source full (default): results/pathway_analysis/upregulated/
                             upregulated_GO_BP_clusters.csv
    --source revigo:         results/pathway_analysis/revigo/upregulated/
                             upregulated_GO_BP_clusters.csv
    CSV must contain columns:
        term_id, term_name, fdr_pvalue, gene_count, cluster, cluster_label

Outputs (results/figures/ by default):
    Default run (--scope both, the default) writes BOTH:
        GO_BP_clusters_top15_thesis.{pdf,png,html}
        GO_BP_clusters_full_thesis.{pdf,png,html}
    --scope top15 or --scope full writes only the corresponding pair.

Usage examples:
    python scripts/plot_go_bp_clusters_thesis.py
    python scripts/plot_go_bp_clusters_thesis.py --scope top15
    python scripts/plot_go_bp_clusters_thesis.py --scope full
    python scripts/plot_go_bp_clusters_thesis.py --source revigo --top-n 10

Dependencies:
    pandas, numpy, plotly, kaleido (for write_image)

================================================================================
"""

from __future__ import annotations

# --- Imports ----------------------------------------------------------------
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go


# =============================================================================
# --- Project paths ----------------------------------------------------------
# =============================================================================
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_FULL_CSV = (
    PROJECT_ROOT
    / "results/pathway_analysis/upregulated/upregulated_GO_BP_clusters.csv"
)
DEFAULT_REVIGO_CSV = (
    PROJECT_ROOT
    / "results/pathway_analysis/revigo/upregulated/upregulated_GO_BP_clusters.csv"
)
DEFAULT_OUTDIR = PROJECT_ROOT / "results/figures/ORA"
DEFAULT_BASENAME_TOPN = "GO_BP_clusters_top15_thesis"
DEFAULT_BASENAME_FULL = "GO_BP_clusters_full_thesis"


# =============================================================================
# --- Style guide constants --------------------------------------------------
#   Anchored in project_notes/figure_and_table_style_guide.Rmd:
#     - Section 2:    Arial inside figures.
#     - Section 3:    manuscript type sizes (14/12/10 pt hierarchy).
#     - Section 9.1:  project palette.
#     - Section 9.1.4: ALL figure-interior text in pure black #000000.
#     - Section 11.2: no on-figure title; title sentence belongs in the
#                     Word caption block, not in the rendered file.
# =============================================================================
FONT_FAMILY = "Arial"
FONT_COLOR = "#000000"
AXIS_TITLE_SIZE = 14
TICK_LABEL_SIZE = 12
LEGEND_TEXT_SIZE = 12
DATA_LABEL_SIZE = 10
SUBHEADING_SIZE = 12

DARK_ACCENT = "#2C3E50"          # Bar outlines (Section 9.1.4).
SEPARATOR_COLOR = "#D0D0D0"      # Dashed cluster separator.
REFERENCE_LINE_COLOR = "#666666"
GRID_COLOR = "#E0E0E0"
WHITE = "#FFFFFF"

# --- Curated 15-color cluster palette ---------------------------------------
#   Anchored at the project ketamine coral (most-significant cluster gets
#   the treatment color), then extended with project blues / purples and a
#   set of accessible accent hues. With more than 15 clusters the palette
#   cycles. Cluster identity is reinforced by the in-bar "Cluster N" text
#   label, so colors function as a secondary cue (Section 9.2: never rely
#   on color alone).
# =============================================================================
CLUSTER_PALETTE = [
    "#E8735A",  # 1  Ketamine coral (project palette anchor)
    "#7FB3D8",  # 2  Control blue (project palette)
    "#9467BD",  # 3  Purple
    "#FF7F0E",  # 4  Orange
    "#2CA02C",  # 5  Green
    "#17BECF",  # 6  Cyan
    "#BCBD22",  # 7  Olive
    "#E377C2",  # 8  Pink
    "#8C564B",  # 9  Brown
    "#1F77B4",  # 10 Blue
    "#D62728",  # 11 Red
    "#7D3C98",  # 12 Deep purple
    "#27AE60",  # 13 Emerald
    "#F39C12",  # 14 Mustard
    "#34495E",  # 15 Navy
]


# =============================================================================
# --- Data loading -----------------------------------------------------------
# =============================================================================
def load_clusters(csv_path: Path) -> pd.DataFrame:
    """Load cluster-assignment CSV and add derived columns."""
    if not csv_path.exists():
        sys.exit(f"ERROR: cluster CSV not found: {csv_path}")
    df = pd.read_csv(csv_path)
    required = {"term_id", "term_name", "fdr_pvalue", "gene_count", "cluster"}
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"ERROR: required columns missing in {csv_path}: {missing}")

    # --- Derived columns ----------------------------------------------------
    df["neg_log10_fdr"] = -np.log10(df["fdr_pvalue"])
    # Sentence case (Section 5): GO term names are stored lowercase.
    df["term_name_display"] = df["term_name"].apply(
        lambda s: (s[0].upper() + s[1:]) if isinstance(s, str) and s else s
    )
    return df


def rank_clusters_and_terms(
    df: pd.DataFrame, top_n_clusters: int | None
) -> pd.DataFrame:
    """Rank clusters by min FDR and terms within each cluster by FDR.

    Parameters
    ----------
    df : DataFrame
        Term-level rows with a 'cluster' column.
    top_n_clusters : int or None
        If a positive integer, keep only the top-N clusters by min FDR.
        None or 0 keeps all clusters.

    Returns
    -------
    DataFrame
        Re-ranked frame with new columns:
            cluster_min_fdr : min FDR within cluster
            cluster_rank    : 1 = most significant cluster
            cluster_color   : palette color for this cluster
        The frame is sorted so the most significant cluster's most
        significant term is at the top (row 0).
    """
    # --- Compute cluster-level summary -------------------------------------
    cluster_summary = (
        df.groupby("cluster")
        .agg(cluster_min_fdr=("fdr_pvalue", "min"))
        .reset_index()
        .sort_values("cluster_min_fdr")
        .reset_index(drop=True)
    )
    cluster_summary["cluster_rank"] = np.arange(1, len(cluster_summary) + 1)

    # --- Top-N cluster filtering (optional) --------------------------------
    if top_n_clusters and top_n_clusters > 0:
        cluster_summary = cluster_summary.head(top_n_clusters).copy()

    # --- Assign palette colors keyed by cluster rank -----------------------
    cluster_summary["cluster_color"] = [
        CLUSTER_PALETTE[(r - 1) % len(CLUSTER_PALETTE)]
        for r in cluster_summary["cluster_rank"]
    ]

    # --- Merge and re-rank --------------------------------------------------
    out = df.merge(cluster_summary, on="cluster", how="inner")
    out = out.sort_values(
        ["cluster_rank", "fdr_pvalue"], ascending=[True, True]
    ).reset_index(drop=True)
    return out


# =============================================================================
# --- Figure construction ----------------------------------------------------
# =============================================================================
def build_figure(
    df: pd.DataFrame,
    show_cluster_headers: bool,
    row_height_px: int = 22,
    min_height_px: int = 500,
) -> go.Figure:
    """Build the clustered horizontal bar figure.

    df is assumed pre-ranked: most-significant cluster's most-significant
    term is row 0, and the y-axis is laid out so row 0 sits at the top.
    """
    if len(df) == 0:
        sys.exit("ERROR: input frame is empty; nothing to plot.")

    n_rows = len(df)

    # --- Y-coordinate layout (descending so row 0 sits at the top) ---------
    df = df.copy()
    df["__y"] = list(range(n_rows - 1, -1, -1))

    # --- Cluster boundary rows (transitions in cluster_rank) ----------------
    cluster_block_first_idx = (
        df.reset_index().groupby("cluster_rank", sort=False)["index"].first().tolist()
    )

    # --- Plot area horizontal range ----------------------------------------
    xmax = float(df["neg_log10_fdr"].max())
    xpad_right = max(1.5, xmax * 0.08)
    x_axis_max = xmax + xpad_right

    fig = go.Figure()

    # --- Bars: one trace per cluster so the legend reads cleanly -----------
    for cluster_id, group in df.groupby("cluster_rank", sort=True):
        bar_color = group["cluster_color"].iloc[0]
        cluster_native = int(group["cluster"].iloc[0])
        representative = group["term_name_display"].iloc[0]
        # Truncate representative for the legend entry only
        legend_repr = (
            representative if len(representative) <= 40
            else representative[:37] + "..."
        )
        fig.add_trace(
            go.Bar(
                x=group["neg_log10_fdr"],
                y=group["__y"],
                orientation="h",
                marker=dict(
                    color=bar_color,
                    line=dict(color=DARK_ACCENT, width=0.6),
                ),
                text=[f"C{cluster_native}"] * len(group),
                textposition="inside",
                insidetextanchor="start",
                textfont=dict(
                    family=FONT_FAMILY,
                    size=DATA_LABEL_SIZE,
                    color=WHITE,
                ),
                name=f"C{cluster_native}: {legend_repr}",
                customdata=group[
                    ["term_id", "term_name_display", "gene_count",
                     "fdr_pvalue", "cluster"]
                ].values,
                hovertemplate=(
                    "<b>%{customdata[1]}</b><br>"
                    "ID: %{customdata[0]}<br>"
                    "Cluster: %{customdata[4]}<br>"
                    "Gene count: %{customdata[2]}<br>"
                    "FDR: %{customdata[3]:.2e}<br>"
                    "−log10(FDR): %{x:.2f}<extra></extra>"
                ),
                showlegend=show_cluster_headers,
            )
        )

    # --- Dashed separators between cluster blocks --------------------------
    # Place a separator at y = (top of next block's y) + 0.5.
    for idx in cluster_block_first_idx[1:]:
        sep_y = df.loc[idx, "__y"] + 0.5
        fig.add_shape(
            type="line",
            x0=0, x1=x_axis_max,
            y0=sep_y, y1=sep_y,
            line=dict(color=SEPARATOR_COLOR, width=0.8, dash="dash"),
            layer="below",
        )

    # --- Y-axis tick labels (term names) ------------------------------------
    yticks = df["__y"].tolist()
    ylabels = df["term_name_display"].tolist()

    # --- Figure height (proportional to row count) -------------------------
    fig_height = max(min_height_px, 90 + n_rows * row_height_px)

    # --- Layout (NO on-figure title; Section 11.2) -------------------------
    fig.update_layout(
        font=dict(family=FONT_FAMILY, size=AXIS_TITLE_SIZE, color=FONT_COLOR),
        plot_bgcolor=WHITE,
        paper_bgcolor=WHITE,
        width=1100,
        height=fig_height,
        margin=dict(l=380, r=240 if show_cluster_headers else 60, t=40, b=70),
        barmode="overlay",
        bargap=0.15,
        xaxis=dict(
            title=dict(
                text="−log<sub>10</sub>(FDR)",
                font=dict(
                    family=FONT_FAMILY, size=AXIS_TITLE_SIZE, color=FONT_COLOR
                ),
            ),
            tickfont=dict(
                family=FONT_FAMILY, size=TICK_LABEL_SIZE, color=FONT_COLOR
            ),
            range=[0, x_axis_max],
            showgrid=True,
            gridcolor=GRID_COLOR,
            gridwidth=0.5,
            zeroline=True,
            zerolinecolor=REFERENCE_LINE_COLOR,
            zerolinewidth=0.8,
            showline=True,
            linecolor=FONT_COLOR,
            linewidth=1,
            ticks="outside",
            ticklen=4,
        ),
        yaxis=dict(
            tickmode="array",
            tickvals=yticks,
            ticktext=ylabels,
            tickfont=dict(
                family=FONT_FAMILY, size=TICK_LABEL_SIZE, color=FONT_COLOR
            ),
            ticklabelstandoff=10,
            range=[-0.7, n_rows - 1 + 0.7],
            showgrid=False,
            zeroline=False,
            showline=False,
            ticks="",
        ),
        legend=dict(
            font=dict(
                family=FONT_FAMILY, size=LEGEND_TEXT_SIZE, color=FONT_COLOR
            ),
            title=dict(
                text="<b>Cluster (representative term)</b>",
                font=dict(
                    family=FONT_FAMILY,
                    size=LEGEND_TEXT_SIZE,
                    color=FONT_COLOR,
                ),
            ),
            x=1.02, xanchor="left",
            y=1.0, yanchor="top",
            bgcolor="rgba(255,255,255,0.85)",
            bordercolor=SEPARATOR_COLOR,
            borderwidth=0.5,
            itemsizing="constant",
            traceorder="normal",
        ),
    )
    return fig


# =============================================================================
# --- Export -----------------------------------------------------------------
# =============================================================================
def export_figure(
    fig: go.Figure, outdir: Path, basename: str, png_scale: int
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    pdf_path = outdir / f"{basename}.pdf"
    png_path = outdir / f"{basename}.png"
    html_path = outdir / f"{basename}.html"

    fig.write_image(pdf_path)                      # vector master
    fig.write_image(png_path, scale=png_scale)     # ~600 DPI for Word
    fig.write_html(html_path, include_plotlyjs="cdn")

    print(f"Wrote: {pdf_path}")
    print(f"Wrote: {png_path}  (scale={png_scale})")
    print(f"Wrote: {html_path}")


# =============================================================================
# --- CLI --------------------------------------------------------------------
# =============================================================================
def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Re-plot the clustered GO:BP enrichment figure with clusters "
            "ranked by min FDR and terms within clusters ranked by FDR, "
            "styled per project_notes/figure_and_table_style_guide.Rmd."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--source",
        choices=["full", "revigo"],
        default="full",
        help=(
            "Cluster-assignment CSV to read. 'full' = the non-REVIGO "
            "clusters CSV in results/pathway_analysis/upregulated/; "
            "'revigo' = the REVIGO-reduced version."
        ),
    )
    parser.add_argument(
        "--clusters-csv",
        type=Path,
        default=None,
        help="Override clusters CSV path (otherwise inferred from --source).",
    )
    parser.add_argument(
        "--scope",
        choices=["both", "top15", "full"],
        default="both",
        help=(
            "Which figure(s) to render. 'top15' = top 15 clusters by min "
            "FDR; 'full' = all clusters; 'both' writes both with distinct "
            "basenames."
        ),
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=15,
        help="Cluster count used when scope includes 'top15'.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=DEFAULT_OUTDIR,
        help="Output directory.",
    )
    parser.add_argument(
        "--basename-top",
        default=DEFAULT_BASENAME_TOPN,
        help="Basename (no extension) for the top-N figure.",
    )
    parser.add_argument(
        "--basename-full",
        default=DEFAULT_BASENAME_FULL,
        help="Basename (no extension) for the all-clusters figure.",
    )
    parser.add_argument(
        "--png-scale",
        type=int,
        default=5,
        help=(
            "Plotly scale factor for PNG export. scale=5 on a 1100-wide "
            "base gives ~600 DPI when embedded at ~7.5 in. width."
        ),
    )
    return parser.parse_args(argv)


def main() -> None:
    args = parse_args()

    # --- Resolve input CSV --------------------------------------------------
    if args.clusters_csv is None:
        args.clusters_csv = (
            DEFAULT_FULL_CSV if args.source == "full" else DEFAULT_REVIGO_CSV
        )

    print(f"Reading: {args.clusters_csv}")
    df = load_clusters(args.clusters_csv)
    n_clusters = df["cluster"].nunique()
    print(f"  rows: {len(df)} terms across {n_clusters} clusters")

    # --- Render requested scopes -------------------------------------------
    targets: list[tuple[str, int | None]] = []
    if args.scope in ("both", "top15"):
        targets.append((args.basename_top, args.top_n))
    if args.scope in ("both", "full"):
        targets.append((args.basename_full, None))

    for basename, top_n_clusters in targets:
        scope_label = (
            f"top {top_n_clusters} clusters" if top_n_clusters else
            "all clusters"
        )
        print(f"\nBuilding figure: {basename} ({scope_label})")

        ranked = rank_clusters_and_terms(df, top_n_clusters)
        n_kept_clusters = ranked["cluster"].nunique()
        n_kept_terms = len(ranked)
        print(f"  rendering {n_kept_terms} terms in {n_kept_clusters} clusters")

        # Show the cluster legend only when the count is legend-friendly
        show_legend = n_kept_clusters <= 20

        fig = build_figure(ranked, show_cluster_headers=show_legend)
        export_figure(fig, args.outdir, basename, args.png_scale)


if __name__ == "__main__":
    main()

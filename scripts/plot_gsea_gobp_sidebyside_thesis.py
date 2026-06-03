#!/usr/bin/env python3
"""
================================================================================
plot_gsea_gobp_sidebyside_thesis.py

Title:         Side-by-side GSEA GO:BP lollipop figure (thesis) - up and down
               panels sharing one diverging NES colorbar and one size legend
Project:       Ketamine Astrocyte Proteomics
Author:        Reina Hastings (reinahastings13@gmail.com)
Date created:  2026-06-01
Last modified: 2026-06-01

Purpose:
    Place the two standalone GSEA GO:BP lollipops (up / down) next to each
    other as two side-by-side panels (up on the LEFT, down on the RIGHT),
    each retaining its own y-axis of term labels and its own |NES|-descending
    ranking (most pronounced pathway at the top), but sharing:
      - ONE diverging NES colorbar spanning the OBSERVED NES range, from the
        most positive (ketamine coral) to the most negative (control blue);
      - ONE leading-edge gene-count size legend.

    This is the side-by-side counterpart to plot_gsea_gobp_combined_thesis.py
    (which stacks the two directions in a single signed-NES-ordered column).
    Both share the project's style constants, the diverging colorscale, and
    the cross-figure marker-size ruler, so the three GSEA GO:BP figures (two
    standalone single-direction + these two combined variants) stay mutually
    comparable.

    Shared visual encoding (per panel):
      y position : |NES| rank within that direction (top = most pronounced)
      x position : -log10(FDR), shared fixed range, FDR=0.05 reference line
      marker size: leading-edge gene count (Count), area-scaled against a
                   shared maximum
      marker color: SIGNED NES on the one diverging coral<->blue scale,
                    mapped to the observed NES range across both panels

    REVIGO simplification (Supek et al., 2011) is applied per direction
    before top-N selection (GO:BP only), reusing the cached TSVs under
    results/gsea/revigo/. Disable with --no-revigo.

    Style follows project_notes/figure_and_table_style_guide.Rmd: Arial,
    14/12/10 pt hierarchy, sentence case, black interior text, project
    palette, no on-figure title.

Inputs:
    --input (optional): GSEA GO:BP results CSV; default inferred from
        --min-size and --metric as
        results/gsea/min{min_size}/{metric}/go_bp_gsea_results.csv
    Required CSV columns: ID, Description, NES, p.adjust, passes_fdr,
        direction, setSize, core_enrichment.

    REVIGO cache (reused for GO:BP):
        results/gsea/revigo/go_bp_up_revigo.tsv
        results/gsea/revigo/go_bp_down_revigo.tsv

Outputs (results/figures/thesis/ by default):
    GSEA_GO_BP_lollipop_sidebyside_thesis.{pdf,png,html}

Usage examples:
    python scripts/plot_gsea_gobp_sidebyside_thesis.py
    python scripts/plot_gsea_gobp_sidebyside_thesis.py --up-n 15 --down-n 15
    python scripts/plot_gsea_gobp_sidebyside_thesis.py --no-revigo
    python scripts/plot_gsea_gobp_sidebyside_thesis.py --refresh-revigo

    copy/paste (defaults, from project root):
        python scripts/plot_gsea_gobp_sidebyside_thesis.py

Dependencies:
    pandas, numpy, plotly, kaleido (PNG/PDF export). Shares helper functions
    and style constants with plot_gsea_gobp_combined_thesis.py (imported).
    Run with the project's conda env, e.g.
        /Users/reina/miniconda3/envs/ketamine_project/bin/python
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
from plotly.subplots import make_subplots

# Reuse the shared loaders, REVIGO handling, selection, style constants, and
# the diverging colorscale from the stacked-version module (same scripts/ dir).
sys.path.insert(0, str(Path(__file__).resolve().parent))
import plot_gsea_gobp_combined_thesis as v1   # noqa: E402


# =============================================================================
# --- Paths / reused constants -----------------------------------------------
# =============================================================================
PROJECT_ROOT = v1.PROJECT_ROOT
DEFAULT_OUTDIR = PROJECT_ROOT / "results/figures/thesis"
DEFAULT_REVIGO_DIR = v1.DEFAULT_REVIGO_DIR

FONT_FAMILY = v1.FONT_FAMILY
FONT_COLOR = v1.FONT_COLOR
AXIS_TITLE_SIZE = v1.AXIS_TITLE_SIZE
TICK_LABEL_SIZE = v1.TICK_LABEL_SIZE
LEGEND_TEXT_SIZE = v1.LEGEND_TEXT_SIZE
DATA_LABEL_SIZE = v1.DATA_LABEL_SIZE
SUBHEADING_SIZE = v1.SUBHEADING_SIZE

KETAMINE_COLOR = v1.KETAMINE_COLOR
CONTROL_COLOR = v1.CONTROL_COLOR
DARK_ACCENT = v1.DARK_ACCENT
SEPARATOR_COLOR = v1.SEPARATOR_COLOR
REFERENCE_LINE_COLOR = v1.REFERENCE_LINE_COLOR
GRID_COLOR = v1.GRID_COLOR
STEM_COLOR = v1.STEM_COLOR
WHITE = v1.WHITE

DEFAULT_X_FLOOR = v1.DEFAULT_X_FLOOR
DEFAULT_X_CEILING = v1.DEFAULT_X_CEILING
FDR_REFLINE_X = v1.FDR_REFLINE_X
MAX_MARKER_SIZE = v1.MAX_MARKER_SIZE
SHARED_MAX_COUNT = v1.SHARED_MAX_COUNT
SHARED_LEGEND_SIZES = v1.SHARED_LEGEND_SIZES
NES_DIVERGING = v1.NES_DIVERGING


# =============================================================================
# --- Helpers ----------------------------------------------------------------
# =============================================================================
def _observed_ticks(cmin: float, cmax: float) -> list[float]:
    """Integer NES ticks for the colorbar within the observed [cmin, cmax]."""
    lo = int(np.ceil(cmin))
    hi = int(np.floor(cmax))
    return [float(v) for v in range(lo, hi + 1)]


def _add_panel(
    fig: go.Figure,
    df_dir: pd.DataFrame,
    col: int,
    cmin: float,
    cmax: float,
    sizeref: float,
    x_floor: float,
    x_ceiling: float,
    show_colorbar: bool,
) -> int:
    """Add one direction's lollipop panel to subplot column `col`.

    df_dir is pre-sorted by |NES| descending (most pronounced first). Returns
    the number of terms plotted (panel height in y-units).
    """
    n = len(df_dir)
    if n == 0:
        return 0
    df_dir = df_dir.copy()
    df_dir["__y"] = list(range(n, 0, -1))      # top of panel = highest |NES|

    # --- Stems ------------------------------------------------------------
    for _, row in df_dir.iterrows():
        fig.add_shape(
            type="line",
            x0=x_floor, x1=row["neg_log10_padj"],
            y0=row["__y"], y1=row["__y"],
            line=dict(color=STEM_COLOR, width=1.2),
            layer="below",
            row=1, col=col,
        )

    # --- FDR=0.05 reference line ------------------------------------------
    if x_floor <= FDR_REFLINE_X <= x_ceiling:
        fig.add_shape(
            type="line",
            x0=FDR_REFLINE_X, x1=FDR_REFLINE_X,
            y0=0.3, y1=n + 1.3,
            line=dict(color=REFERENCE_LINE_COLOR, width=1.2, dash="dash"),
            layer="below",
            row=1, col=col,
        )

    # --- Markers (signed NES color on the shared diverging scale) --------
    colorbar = None
    if show_colorbar:
        colorbar = dict(
            title=dict(
                text="NES",
                font=dict(family=FONT_FAMILY, size=LEGEND_TEXT_SIZE,
                          color=FONT_COLOR),
                side="right",
            ),
            tickfont=dict(family=FONT_FAMILY, size=TICK_LABEL_SIZE,
                          color=FONT_COLOR),
            tickmode="array",
            tickvals=_observed_ticks(cmin, cmax),
            thickness=16,
            len=0.55,
            y=0.74,
            yanchor="middle",
            x=1.015,
            xanchor="left",
            outlinewidth=0,
        )

    fig.add_trace(
        go.Scatter(
            x=df_dir["neg_log10_padj"],
            y=df_dir["__y"],
            mode="markers",
            marker=dict(
                size=df_dir["count"].clip(lower=1).astype(float),
                sizemode="area",
                sizeref=sizeref,
                sizemin=4,
                color=df_dir["NES"].astype(float),
                colorscale=NES_DIVERGING,
                cmin=cmin,
                cmax=cmax,
                showscale=show_colorbar,
                line=dict(color=DARK_ACCENT, width=0.6),
                colorbar=colorbar,
            ),
            customdata=df_dir[
                ["ID", "term_name_display", "NES", "p.adjust",
                 "count", "setSize"]
            ].values,
            hovertemplate=(
                "<b>%{customdata[1]}</b><br>"
                "ID: %{customdata[0]}<br>"
                "NES: %{customdata[2]:+.2f}<br>"
                "FDR: %{customdata[3]:.2e}<br>"
                "Leading-edge count: %{customdata[4]}<br>"
                "Set size: %{customdata[5]}<extra></extra>"
            ),
            showlegend=False,
        ),
        row=1, col=col,
    )

    # --- Y tick labels for this panel ------------------------------------
    axis_key = "yaxis" if col == 1 else f"yaxis{col}"
    fig.update_layout(**{
        axis_key: dict(
            tickmode="array",
            tickvals=df_dir["__y"].tolist(),
            ticktext=df_dir["term_name_display"].tolist(),
            tickfont=dict(family=FONT_FAMILY, size=TICK_LABEL_SIZE,
                          color=FONT_COLOR),
            ticklabelstandoff=8,
            range=[0.3, n + 1.7],
            showgrid=False, zeroline=False, showline=False, ticks="",
        )
    })
    return n


# =============================================================================
# --- Figure construction ----------------------------------------------------
# =============================================================================
def build_sidebyside_figure(
    df_up: pd.DataFrame,
    df_down: pd.DataFrame,
    x_floor: float = DEFAULT_X_FLOOR,
    x_ceiling: float = DEFAULT_X_CEILING,
) -> go.Figure:
    """Build the side-by-side up|down GSEA GO:BP lollipop figure."""
    if len(df_up) + len(df_down) == 0:
        sys.exit("ERROR: no terms to plot.")

    # --- Shared diverging color range = observed NES range ----------------
    all_nes = pd.concat([df_up["NES"], df_down["NES"]]).astype(float)
    cmin = float(all_nes.min())
    cmax = float(all_nes.max())
    if cmax - cmin < 1e-6:
        cmax = cmin + 1e-3

    # --- Shared marker-size ruler ----------------------------------------
    sizeref = 2.0 * SHARED_MAX_COUNT / (MAX_MARKER_SIZE ** 2)

    # Wide horizontal spacing so the right panel's y-axis term labels have
    # room to sit between the two plot areas.
    fig = make_subplots(
        rows=1, cols=2,
        horizontal_spacing=0.30,
        subplot_titles=("", ""),
    )

    # Up on the left (col 1) carries the shared colorbar; down on the right
    # (col 2) reuses the same cmin/cmax so its colors are on the same scale.
    n_up = _add_panel(fig, df_up, col=1, cmin=cmin, cmax=cmax,
                      sizeref=sizeref, x_floor=x_floor, x_ceiling=x_ceiling,
                      show_colorbar=True)
    n_down = _add_panel(fig, df_down, col=2, cmin=cmin, cmax=cmax,
                        sizeref=sizeref, x_floor=x_floor, x_ceiling=x_ceiling,
                        show_colorbar=False)

    # --- Shared leading-edge count size legend (dummy traces) ------------
    for n in SHARED_LEGEND_SIZES:
        fig.add_trace(
            go.Scatter(
                x=[None], y=[None], mode="markers",
                marker=dict(
                    size=[float(n)], sizemode="area", sizeref=sizeref,
                    sizemin=4, color=WHITE,
                    line=dict(color=DARK_ACCENT, width=0.6),
                ),
                name=f"{n}", showlegend=True, hoverinfo="skip",
            ),
            row=1, col=1,
        )

    # --- Direction labels above each panel -------------------------------
    x_center = (x_floor + x_ceiling) / 2.0
    fig.add_annotation(
        x=x_center, y=n_up + 0.9, xref="x", yref="y",
        text="<b>Up in ketamine</b>", showarrow=False,
        font=dict(family=FONT_FAMILY, size=SUBHEADING_SIZE,
                  color=KETAMINE_COLOR),
        align="center",
    )
    fig.add_annotation(
        x=x_center, y=n_down + 0.9, xref="x2", yref="y2",
        text="<b>Down in ketamine</b>", showarrow=False,
        font=dict(family=FONT_FAMILY, size=SUBHEADING_SIZE,
                  color=CONTROL_COLOR),
        align="center",
    )

    # --- FDR labels (one per panel, in paper-relative y) -----------------
    for xref in ("x", "x2"):
        fig.add_annotation(
            x=FDR_REFLINE_X, y=0, xref=xref, yref="paper",
            text="FDR = 0.05", showarrow=False,
            font=dict(family=FONT_FAMILY, size=DATA_LABEL_SIZE,
                      color=REFERENCE_LINE_COLOR),
            xanchor="left", yanchor="bottom", xshift=4, yshift=4,
        )

    # --- Shared x-axis styling for both panels ---------------------------
    n_rows = max(n_up, n_down)
    x_axis_style = dict(
        title=dict(
            text="−log<sub>10</sub>(FDR)",
            font=dict(family=FONT_FAMILY, size=AXIS_TITLE_SIZE,
                      color=FONT_COLOR),
        ),
        tickfont=dict(family=FONT_FAMILY, size=TICK_LABEL_SIZE,
                      color=FONT_COLOR),
        range=[x_floor, x_ceiling],
        showgrid=True, gridcolor=GRID_COLOR, gridwidth=0.5,
        zeroline=True, zerolinecolor=REFERENCE_LINE_COLOR, zerolinewidth=0.8,
        showline=True, linecolor=FONT_COLOR, linewidth=1,
        ticks="outside", ticklen=4,
    )
    fig.update_layout(xaxis=x_axis_style, xaxis2=x_axis_style)

    # --- Global layout (NO on-figure title; Section 11.2) ----------------
    fig.update_layout(
        font=dict(family=FONT_FAMILY, size=AXIS_TITLE_SIZE, color=FONT_COLOR),
        plot_bgcolor=WHITE,
        paper_bgcolor=WHITE,
        width=1500,
        height=int(max(480, 36 * n_rows + 170)),
        margin=dict(l=300, r=180, t=60, b=70),
        legend=dict(
            font=dict(family=FONT_FAMILY, size=LEGEND_TEXT_SIZE,
                      color=FONT_COLOR),
            title=dict(
                text="<b>Leading-edge<br>gene count</b>",
                font=dict(family=FONT_FAMILY, size=LEGEND_TEXT_SIZE,
                          color=FONT_COLOR),
            ),
            x=1.015, xanchor="left",
            y=0.12, yanchor="bottom",
            bgcolor="rgba(255,255,255,0.85)",
            bordercolor=SEPARATOR_COLOR, borderwidth=0.5,
            itemsizing="trace",
        ),
    )
    return fig


# =============================================================================
# --- CLI --------------------------------------------------------------------
# =============================================================================
def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Place the up- and down-direction GSEA GO:BP lollipops side by "
            "side (up left, down right) sharing one diverging NES colorbar "
            "(observed range) and one leading-edge-count size legend. Styled "
            "per project_notes/figure_and_table_style_guide.Rmd."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input", type=Path, default=None,
                        help="Override GO:BP GSEA results CSV path.")
    parser.add_argument("--min-size", type=int, default=10, choices=[10, 15],
                        help="Which gsea_analysis.R min_size run to read.")
    parser.add_argument("--metric", choices=["log2fc", "signed_logp"],
                        default="log2fc",
                        help="Which ranking-metric subdirectory to read.")
    parser.add_argument("--up-n", type=int, default=15,
                        help="Number of top upregulated terms (left panel).")
    parser.add_argument("--down-n", type=int, default=15,
                        help="Number of top downregulated terms (right panel).")
    parser.add_argument("--x-floor", type=float, default=DEFAULT_X_FLOOR,
                        help="Left edge of the -log10(FDR) axis.")
    parser.add_argument("--x-ceiling", type=float, default=DEFAULT_X_CEILING,
                        help="Right edge of the -log10(FDR) axis (shared).")
    parser.add_argument("--no-revigo", action="store_true",
                        help="Disable REVIGO simplification.")
    parser.add_argument("--refresh-revigo", action="store_true",
                        help="Force re-fetch of REVIGO outputs.")
    parser.add_argument("--revigo-dir", type=Path, default=DEFAULT_REVIGO_DIR,
                        help="Directory for REVIGO input/output cache.")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR,
                        help="Figure output directory.")
    parser.add_argument("--basename", type=str,
                        default="GSEA_GO_BP_lollipop_sidebyside_thesis",
                        help="Output file basename (no extension).")
    parser.add_argument("--png-scale", type=int, default=5,
                        help="Plotly scale factor for PNG export.")
    return parser.parse_args(argv)


def main() -> None:
    args = parse_args()
    use_revigo = not args.no_revigo

    if args.input is None:
        args.input = v1.default_input_csv(args.min_size, args.metric)

    print(f"Reading: {args.input}")
    df = v1.load_gsea(args.input)
    n_sig = int(df["passes_fdr"].sum())
    print(f"  rows total: {len(df)} ; significant (FDR<=0.05): {n_sig}")

    revigo_reps = None
    if use_revigo:
        print("REVIGO simplification enabled (GO:BP).")
        revigo_reps = v1.revigo_filter_ids(
            df_sig=df[df["passes_fdr"]],
            cache_dir=args.revigo_dir,
            refresh=args.refresh_revigo,
        )
    else:
        print("REVIGO simplification disabled by --no-revigo.")

    df_up = v1.select_top_for_direction(df, "up", args.up_n, revigo_reps)
    df_down = v1.select_top_for_direction(df, "down", args.down_n, revigo_reps)

    if len(df_up) < args.up_n:
        print(f"  NOTE: only {len(df_up)} up terms available "
              f"(requested {args.up_n}).")
    if len(df_down) < args.down_n:
        print(f"  NOTE: only {len(df_down)} down terms available "
              f"(requested {args.down_n}).")

    print(f"\nLeft panel: top {len(df_up)} up terms (|NES| desc).")
    print(f"Right panel: top {len(df_down)} down terms (|NES| desc).")
    nes_all = pd.concat([df_up['NES'], df_down['NES']])
    print(f"Shared colorbar NES range (observed): "
          f"{nes_all.min():+.2f} to {nes_all.max():+.2f}")

    fig = build_sidebyside_figure(
        df_up, df_down, x_floor=args.x_floor, x_ceiling=args.x_ceiling,
    )
    v1.export_figure(fig, args.outdir, args.basename, args.png_scale)


if __name__ == "__main__":
    main()

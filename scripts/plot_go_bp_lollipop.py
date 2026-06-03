#!/usr/bin/env python3
"""
================================================================================
plot_go_bp_lollipop.py

Title:         GO Biological Process enrichment lollipop figure (thesis)
Project:       Ketamine Astrocyte Proteomics
Author:        Reina Hastings (reinahastings13@gmail.com)
Date created:  2026-05-27
Last modified: 2026-05-27

Purpose:
    Generate a publication-grade lollipop figure summarising GO:BP
    over-representation results from g:Profiler. By default reads the
    REVIGO-filtered upregulated set and splits the display into the
    top N most-enriched terms above a dashed separator and a curated
    list of 'additional terms of interest' below it. The figure is
    compliant with project_notes/figure_and_table_style_guide.Rmd
    (Arial, 14/12/10 pt, sentence case, black interior text, project
    palette, no on-figure title).

Inputs:
    --source revigo (default): results/pathway_analysis/revigo/upregulated/
                               upregulated_GO_BP_revigo.csv
    --source full:             results/pathway_analysis/upregulated/
                               upregulated_GO_BP.csv
    Either CSV must contain columns: term_id, term_name, fdr_pvalue,
    gene_count.

Outputs (written to results/figures/ by default):
    GO_BP_enrichment_lollipop_thesis.pdf   (vector, canonical)
    GO_BP_enrichment_lollipop_thesis.png   (raster, scale=5 -> ~600 DPI
                                            when embedded at ~7.5 in.
                                            width in Word)
    GO_BP_enrichment_lollipop_thesis.html  (interactive plotly export)

Usage examples:
    python scripts/plot_go_bp_lollipop.py
    python scripts/plot_go_bp_lollipop.py --colormap viridis
    python scripts/plot_go_bp_lollipop.py --source full --top-n 12
    python scripts/plot_go_bp_lollipop.py --additional-terms \
        "phospholipid homeostasis" "neuron projection regeneration"

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
DEFAULT_REVIGO_CSV = (
    PROJECT_ROOT
    / "results/pathway_analysis/revigo/upregulated/upregulated_GO_BP_revigo.csv"
)
DEFAULT_FULL_CSV = (
    PROJECT_ROOT / "results/pathway_analysis/upregulated/upregulated_GO_BP.csv"
)
# DEPRECATED: this GO:BP-only ORA lollipop has been superseded by
# plot_ora_lollipops.py (--db go_bp), which mirrors the GSEA lollipop style and
# naming (ORA_GO_BP_lollipop_up_thesis) and also covers KEGG and Reactome. The
# GO_BP_enrichment_lollipop_thesis output has been retired from results/figures/.
# Kept for provenance; outputs now land in results/figures/ORA/ if re-run.
DEFAULT_OUTDIR = PROJECT_ROOT / "results/figures/ORA"
DEFAULT_BASENAME = "GO_BP_enrichment_lollipop_thesis"


# =============================================================================
# --- Style guide constants --------------------------------------------------
#   Anchored in project_notes/figure_and_table_style_guide.Rmd:
#     - Section 2:    Arial inside figures.
#     - Section 3:    manuscript type sizes.
#     - Section 9.1:  project palette (ketamine coral, separators, grays).
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

KETAMINE_COLOR = "#E8735A"       # Treatment direction; ramp anchor.
DARK_ACCENT = "#2C3E50"          # Marker outlines (Section 9.1.4).
SEPARATOR_COLOR = "#D0D0D0"      # Dashed cluster separator.
REFERENCE_LINE_COLOR = "#666666"
GRID_COLOR = "#E0E0E0"
STEM_COLOR = "#B0B0B0"           # Lollipop stems (secondary annotation).
WHITE = "#FFFFFF"

# Single-hue sequential ramp from very light coral to project ketamine
# accent. Section 9.2 prefers viridis/cividis for sequential ramps; this
# figure uses a coral ramp because the encoded magnitude is semantically
# "up in ketamine," and the high end matches the project palette's
# ketamine direction color. Switchable via --colormap. Documented per
# Section 15 (deviations are documented inline).
KETAMINE_SEQUENTIAL = [
    [0.00, "#FCEDE8"],
    [0.25, "#F5B8A6"],
    [0.50, "#EF947D"],
    [0.75, "#EB825F"],
    [1.00, KETAMINE_COLOR],
]


# =============================================================================
# --- Curated default 'additional terms of interest' -------------------------
#   Matches the 2026-02-02 thesis figure
#   (results/figures/GO_BP_enrichment_lollipop_thesis.png). Override
#   with --additional-terms on the CLI.
# =============================================================================
DEFAULT_ADDITIONAL_TERMS = [
    "regulation of cell projection organization",
    "neuron projection regeneration",
    "regulation of fatty acid metabolic process",
    "regulation of small molecule metabolic process",
    "phospholipid homeostasis",
]


# =============================================================================
# --- Data loading -----------------------------------------------------------
# =============================================================================
def load_enrichment(csv_path: Path) -> pd.DataFrame:
    """Load g:Profiler enrichment CSV; add neg_log10_fdr and display name."""
    if not csv_path.exists():
        sys.exit(f"ERROR: enrichment CSV not found: {csv_path}")
    df = pd.read_csv(csv_path)
    required = {"term_id", "term_name", "fdr_pvalue", "gene_count"}
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


def select_top_and_additional(
    df: pd.DataFrame,
    top_n: int,
    additional_names: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Split enrichment table into (top-N by FDR, curated additional list)."""
    df_sorted = df.sort_values("fdr_pvalue").reset_index(drop=True)
    top = df_sorted.head(top_n).copy()

    # --- Case-insensitive lookup of curated additional terms ---------------
    lower_to_canonical = {n.lower(): n for n in df_sorted["term_name"]}
    matched = []
    missing = []
    for name in additional_names:
        canonical = lower_to_canonical.get(name.lower())
        if canonical is not None:
            matched.append(canonical)
        else:
            missing.append(name)
    if missing:
        print(
            f"WARN: additional terms not found in input (skipped): {missing}",
            file=sys.stderr,
        )

    additional = df_sorted[df_sorted["term_name"].isin(matched)].copy()
    additional = additional[~additional["term_id"].isin(top["term_id"])]

    # --- Preserve the user-supplied order for the curated set --------------
    order_lookup = {n.lower(): i for i, n in enumerate(additional_names)}
    additional["__order"] = additional["term_name"].str.lower().map(order_lookup)
    additional = (
        additional.sort_values("__order")
        .drop(columns="__order")
        .reset_index(drop=True)
    )
    return top, additional


# =============================================================================
# --- Figure construction ----------------------------------------------------
# =============================================================================
def build_figure(
    top: pd.DataFrame,
    additional: pd.DataFrame,
    colorscale,
    max_marker_size: int = 30,
) -> go.Figure:
    """Build the lollipop figure with two y-stacked groups and a separator."""
    n_top = len(top)
    n_add = len(additional)
    if n_top == 0:
        sys.exit("ERROR: top group is empty; nothing to plot.")

    # --- Y-coordinate layout ------------------------------------------------
    #   Higher y => higher on figure.
    #   Stack (top -> bottom of figure):
    #     top group:   y = n_top+n_add .. n_add+1   (most significant on top)
    #     separator:   y = n_add                     (dashed line + label)
    #     additional:  y = n_add-1 .. 0              (first curated just below sep)
    top = top.copy()
    additional = additional.copy()
    top["__y"] = list(range(n_top + n_add, n_add, -1))
    additional["__y"] = list(range(n_add - 1, -1, -1))
    sep_y = n_add

    combined = pd.concat([top, additional], ignore_index=True)

    # --- Color scale range across all displayed terms ----------------------
    cmin = float(combined["neg_log10_fdr"].min())
    cmax = float(combined["neg_log10_fdr"].max())

    # --- Marker sizing (sizeref per plotly bubble convention) --------------
    size_values = combined["gene_count"].clip(lower=1).astype(float)
    sizeref = 2.0 * size_values.max() / (max_marker_size ** 2)

    # --- Plot area horizontal range ----------------------------------------
    xmax = float(combined["neg_log10_fdr"].max())
    xpad_right = max(1.5, xmax * 0.10)
    x_axis_max = xmax + xpad_right

    fig = go.Figure()

    # --- Lollipop stems (thin gray; markers carry the FDR color) -----------
    for _, row in combined.iterrows():
        fig.add_shape(
            type="line",
            x0=0, x1=row["neg_log10_fdr"],
            y0=row["__y"], y1=row["__y"],
            line=dict(color=STEM_COLOR, width=1.2),
            layer="below",
        )

    # --- Markers ------------------------------------------------------------
    fig.add_trace(
        go.Scatter(
            x=combined["neg_log10_fdr"],
            y=combined["__y"],
            mode="markers",
            marker=dict(
                size=size_values,
                sizemode="area",
                sizeref=sizeref,
                sizemin=4,
                color=combined["neg_log10_fdr"],
                colorscale=colorscale,
                cmin=cmin,
                cmax=cmax,
                line=dict(color=DARK_ACCENT, width=0.6),
                colorbar=dict(
                    title=dict(
                        text="−log<sub>10</sub>(FDR)",
                        font=dict(
                            family=FONT_FAMILY,
                            size=LEGEND_TEXT_SIZE,
                            color=FONT_COLOR,
                        ),
                        side="right",
                    ),
                    tickfont=dict(
                        family=FONT_FAMILY,
                        size=TICK_LABEL_SIZE,
                        color=FONT_COLOR,
                    ),
                    thickness=14,
                    len=0.55,
                    y=0.72,
                    yanchor="middle",
                    x=1.02,
                    xanchor="left",
                    outlinewidth=0,
                ),
            ),
            customdata=combined[
                ["term_id", "term_name_display", "gene_count", "fdr_pvalue"]
            ].values,
            hovertemplate=(
                "<b>%{customdata[1]}</b><br>"
                "ID: %{customdata[0]}<br>"
                "Gene count: %{customdata[2]}<br>"
                "FDR: %{customdata[3]:.2e}<br>"
                "−log10(FDR): %{x:.2f}<extra></extra>"
            ),
            showlegend=False,
        )
    )

    # --- Gene-count size legend (3 dummy traces, visually scaled correctly)
    legend_counts = [5, 10, 20]
    for n in legend_counts:
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker=dict(
                    size=[float(n)],
                    sizemode="area",
                    sizeref=sizeref,
                    sizemin=4,
                    color=WHITE,
                    line=dict(color=DARK_ACCENT, width=0.6),
                ),
                name=f"{n} genes",
                showlegend=True,
                hoverinfo="skip",
            )
        )

    # --- Y-axis term labels -------------------------------------------------
    yticks = combined["__y"].tolist()
    ylabels = combined["term_name_display"].tolist()

    # --- Dashed separator line between the two groups ----------------------
    fig.add_shape(
        type="line",
        x0=0, x1=x_axis_max,
        y0=sep_y, y1=sep_y,
        line=dict(color=SEPARATOR_COLOR, width=1.2, dash="dash"),
        layer="below",
    )

    # --- In-figure section labels (Section 11.1: cluster identifiers OK) ---
    label_x = x_axis_max / 2
    fig.add_annotation(
        x=label_x, y=n_top + n_add + 1.0,
        xref="x", yref="y",
        text="<b>Top enriched terms</b>",
        showarrow=False,
        font=dict(family=FONT_FAMILY, size=SUBHEADING_SIZE, color=FONT_COLOR),
        align="center",
    )
    if n_add > 0:
        fig.add_annotation(
            x=label_x, y=sep_y,
            xref="x", yref="y",
            text="<b>Additional terms of interest</b>",
            showarrow=False,
            font=dict(
                family=FONT_FAMILY, size=SUBHEADING_SIZE, color=FONT_COLOR
            ),
            align="center",
            bgcolor=WHITE,         # mask the dashed line behind the text
            borderpad=3,
        )

    # --- Layout (NO on-figure title; Section 11.2) -------------------------
    fig.update_layout(
        font=dict(family=FONT_FAMILY, size=AXIS_TITLE_SIZE, color=FONT_COLOR),
        plot_bgcolor=WHITE,
        paper_bgcolor=WHITE,
        width=900,
        height=680,
        margin=dict(l=315, r=150, t=40, b=70),
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
            ticklabelstandoff=10,    # gap between labels and the y axis line
            range=[-0.7, n_top + n_add + 1.7],
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
                text="<b>Gene count</b>",
                font=dict(
                    family=FONT_FAMILY,
                    size=LEGEND_TEXT_SIZE,
                    color=FONT_COLOR,
                ),
            ),
            x=1.02, xanchor="left",
            y=0.05, yanchor="bottom",
            bgcolor="rgba(255,255,255,0.85)",
            bordercolor=SEPARATOR_COLOR,
            borderwidth=0.5,
            itemsizing="trace",
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
            "Generate the GO:BP enrichment lollipop thesis figure, "
            "styled per project_notes/figure_and_table_style_guide.Rmd."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--source",
        choices=["revigo", "full"],
        default="revigo",
        help="REVIGO-filtered or full upregulated GO:BP ORA CSV.",
    )
    parser.add_argument(
        "--enrichment-csv",
        type=Path,
        default=None,
        help="Override enrichment CSV path (otherwise inferred from --source).",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=10,
        help="Number of top-enriched terms shown above the separator.",
    )
    parser.add_argument(
        "--additional-terms",
        nargs="*",
        default=DEFAULT_ADDITIONAL_TERMS,
        help=(
            "Term names to show below the separator (case-insensitive "
            "match on term_name; pass an empty list to omit this group)."
        ),
    )
    parser.add_argument(
        "--colormap",
        choices=["ketamine", "viridis", "cividis"],
        default="ketamine",
        help=(
            "Sequential colormap for -log10(FDR). 'ketamine' is the "
            "single-hue coral ramp anchored at #E8735A; 'viridis' and "
            "'cividis' follow the Section 9.2 default preference."
        ),
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=DEFAULT_OUTDIR,
        help="Output directory.",
    )
    parser.add_argument(
        "--basename",
        default=DEFAULT_BASENAME,
        help="Output file basename (no extension).",
    )
    parser.add_argument(
        "--png-scale",
        type=int,
        default=5,
        help=(
            "Plotly scale factor for PNG export. scale=5 on a 900x680 "
            "base gives ~600 DPI when embedded at ~7.5 in. width."
        ),
    )
    return parser.parse_args(argv)


def main() -> None:
    args = parse_args()

    # --- Resolve input CSV --------------------------------------------------
    if args.enrichment_csv is None:
        args.enrichment_csv = (
            DEFAULT_REVIGO_CSV if args.source == "revigo" else DEFAULT_FULL_CSV
        )

    print(f"Reading: {args.enrichment_csv}")
    df = load_enrichment(args.enrichment_csv)
    print(f"  rows: {len(df)}")

    # --- Build the two displayed groups ------------------------------------
    top, additional = select_top_and_additional(
        df, args.top_n, args.additional_terms
    )
    print(f"Top {len(top)} enriched terms:")
    for _, r in top.iterrows():
        print(f"  - {r['term_name_display']} (FDR {r['fdr_pvalue']:.2e})")
    print(f"Additional terms of interest ({len(additional)}):")
    for _, r in additional.iterrows():
        print(f"  - {r['term_name_display']} (FDR {r['fdr_pvalue']:.2e})")

    # --- Resolve colormap ---------------------------------------------------
    if args.colormap == "ketamine":
        colorscale = KETAMINE_SEQUENTIAL
    else:
        colorscale = args.colormap.capitalize()   # plotly accepts string names

    # --- Render and write ---------------------------------------------------
    fig = build_figure(top, additional, colorscale)
    export_figure(fig, args.outdir, args.basename, args.png_scale)


if __name__ == "__main__":
    main()
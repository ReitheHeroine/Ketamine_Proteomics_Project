#!/usr/bin/env python3
"""
================================================================================
plot_ora_lollipops.py

Title:         ORA enrichment lollipop figures (thesis) - upregulated, one DB
               per figure, mirroring the GSEA lollipops
Project:       Ketamine Astrocyte Proteomics
Author:        Reina Hastings (reinahastings13@gmail.com)
Date created:  2026-06-01
Last modified: 2026-06-01

Purpose:
    Generate publication-grade single-direction lollipop figures summarising
    over-representation analysis (ORA, g:Profiler) results for one of three
    databases (GO:BP, KEGG, Reactome). These mirror the GSEA lollipops
    (plot_gsea_lollipops.py) in layout, style, and file naming so the ORA and
    GSEA panels can sit side by side in the thesis. ORA in this project is
    upregulated-only (the downregulated set yields no significantly enriched
    terms), so each invocation writes a single "up" figure per database.

    Visual encoding (per figure) - mirrors the GSEA lollipop style, with fold
    enrichment as the ORA effect-size analog of GSEA's NES:
      y-axis ordering:   terms ranked by FDR (most significant at the top).
                         Unlike GSEA (which orders by its color metric |NES|),
                         ORA orders by significance because fold enrichment
                         collapses to a single value for every fully-covered
                         term (recall = 1) and would otherwise front-load a
                         block of tied small terms.
      x-axis:            -log10(FDR), so the most statistically robust
                         enrichments extend furthest right.
      marker size:       intersection gene count (query genes annotated to the
                         term), area-scaled - the ORA analog of GSEA's
                         leading-edge count.
      marker color:      sequential ketamine-coral ramp by fold enrichment.
                         Fold enrichment = (intersection/query) /
                         (term_size/effective_domain_size): the observed/
                         expected over-representation ratio. It is the ORA
                         effect-size analog to GSEA's |NES| and here varies
                         independently of the FDR-driven y-order, adding an
                         orthogonal effect-size signal.

    For GO:BP, REVIGO simplification (Supek et al., 2011) is applied before
    top-N selection so semantically redundant terms do not dominate the figure.
    REVIGO representatives are read from the cached ORA REVIGO output produced by
    pathway_analysis.py. REVIGO is GO-only and is skipped for KEGG and Reactome.

    Style follows project_notes/figure_and_table_style_guide.Rmd: Arial,
    14/12/10 pt hierarchy, sentence case, black interior text, project palette,
    no on-figure title. Cross-figure scaling constants make all three ORA
    lollipops directly comparable (shared x-range, marker-size scale, legend
    reference dots), exactly as the GSEA set does.

Inputs:
    --db {go_bp, kegg, reactome}: which database to plot.
    --input (optional): ORA "_full" results CSV; default inferred from --db as
        results/pathway_analysis/upregulated/upregulated_{token}_full.csv
        Required columns: term_id, term_name, fdr_pvalue, term_size,
        query_size, intersection_size, effective_domain_size.
        (Produced by pathway_analysis.py; see write_full_results_csv there.)

    REVIGO representatives (GO:BP only):
        results/pathway_analysis/revigo/upregulated/upregulated_GO_BP_revigo.csv
        Required column: term_id.

Outputs (results/figures/ORA/ by default):
    ORA_GO_BP_lollipop_up_thesis.{pdf,png,html}
    ORA_KEGG_lollipop_up_thesis.{pdf,png,html}
    ORA_Reactome_lollipop_up_thesis.{pdf,png,html}

Usage examples:
    python scripts/plot_ora_lollipops.py --db go_bp
    python scripts/plot_ora_lollipops.py --db kegg --top-n 12
    python scripts/plot_ora_lollipops.py --db reactome
    python scripts/plot_ora_lollipops.py --db go_bp --no-revigo

    copy/paste (all three DBs):
        for db in go_bp kegg reactome; do \\
            python scripts/plot_ora_lollipops.py --db $db; \\
        done

Dependencies:
    pandas, numpy, plotly, kaleido (PNG export). Imports style constants and
    the figure exporter from plot_gsea_lollipops.py (must live in the same
    scripts/ directory).

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

# Shared style constants, color ramp, sentence-casing, and figure exporter.
# Reusing them guarantees the ORA lollipops are visually identical to the GSEA
# lollipops except for the encoded metric.
from plot_gsea_lollipops import (
    AXIS_TITLE_SIZE,
    DARK_ACCENT,
    DATA_LABEL_SIZE,
    FDR_REFLINE_X,
    FDR_THRESHOLD,
    FONT_COLOR,
    FONT_FAMILY,
    GRID_COLOR,
    KETAMINE_COLOR,
    KETAMINE_SEQUENTIAL,
    LEGEND_TEXT_SIZE,
    MAX_MARKER_SIZE,
    REFERENCE_LINE_COLOR,
    SEPARATOR_COLOR,
    STEM_COLOR,
    SUBHEADING_SIZE,
    TICK_LABEL_SIZE,
    WHITE,
    _sentence_case_term,
    export_figure,
)


# =============================================================================
# --- Project paths ----------------------------------------------------------
# =============================================================================
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUTDIR = PROJECT_ROOT / "results/figures/ORA"
DEFAULT_REVIGO_CSV = (
    PROJECT_ROOT
    / "results/pathway_analysis/revigo/upregulated/upregulated_GO_BP_revigo.csv"
)


# =============================================================================
# --- Database configuration -------------------------------------------------
# =============================================================================
DB_CONFIG = {
    "go_bp":    {"label": "GO Biological Process", "short": "GO_BP",
                  "token": "GO_BP", "supports_revigo": True},
    "kegg":     {"label": "KEGG",                  "short": "KEGG",
                  "token": "KEGG", "supports_revigo": False},
    "reactome": {"label": "Reactome",              "short": "Reactome",
                  "token": "Reactome", "supports_revigo": False},
}


# =============================================================================
# --- Cross-figure scaling constants (ORA-specific) --------------------------
#   Tuned to the ORA upregulated value ranges, which are smaller than GSEA's:
#     - intersection counts top out around 21 (GO:BP), 13 (KEGG), 8 (Reactome).
#     - -log10(FDR) tops out near 12.9 (GO:BP synaptic vesicle cycle).
#   Holding these constant across the three ORA lollipops gives them a shared
#   ruler, the same purpose the GSEA scaling constants serve there.
# =============================================================================
DEFAULT_X_FLOOR = 1.0
DEFAULT_X_CEILING = 13.5
SHARED_MAX_COUNT = 25          # marker-size reference (covers GO:BP max ~21)
SHARED_LEGEND_SIZES = [5, 10, 20]


# =============================================================================
# --- Data loading -----------------------------------------------------------
# =============================================================================
def load_ora_full(csv_path: Path) -> pd.DataFrame:
    """Load an ORA _full results CSV; add fold enrichment, -log10(FDR), names.

    Adds:
      neg_log10_fdr     : -log10(fdr_pvalue), clipped to avoid inf.
      fold_enrichment   : (intersection/query) / (term_size/effective_domain).
      count             : intersection gene count (size encoding).
      term_name_display : sentence-cased term name.
    """
    if not csv_path.exists():
        sys.exit(f"ERROR: ORA results CSV not found: {csv_path}\n"
                 f"       Generate it with pathway_analysis.py "
                 f"(writes <set>_<token>_full.csv).")
    df = pd.read_csv(csv_path)
    required = {"term_id", "term_name", "fdr_pvalue", "term_size",
                "query_size", "intersection_size", "effective_domain_size"}
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"ERROR: required columns missing in {csv_path}: {missing}")

    df["neg_log10_fdr"] = -np.log10(df["fdr_pvalue"].clip(lower=1e-300))

    # Fold enrichment: observed fraction of the query in the term divided by
    # the fraction expected by chance given the term's size relative to the
    # effective annotation universe. Effective_domain_size is constant within a
    # source/query, so the ranking is monotonic in (precision / term_size), but
    # the absolute value is the interpretable over-representation ratio.
    observed = df["intersection_size"] / df["query_size"]
    expected = df["term_size"] / df["effective_domain_size"]
    df["fold_enrichment"] = observed / expected.replace(0, np.nan)

    df["count"] = df["intersection_size"].astype(float)
    df["term_name_display"] = df["term_name"].apply(_sentence_case_term)
    return df


def default_input_csv(token: str) -> Path:
    """Construct default ORA _full CSV path from the database file token."""
    return (
        PROJECT_ROOT / "results/pathway_analysis/upregulated"
        / f"upregulated_{token}_full.csv"
    )


def load_revigo_representative_ids(csv_path: Path) -> set[str]:
    """Return the set of REVIGO representative GO term IDs (GO:BP only)."""
    if not csv_path.exists():
        sys.exit(f"ERROR: REVIGO representatives CSV not found: {csv_path}\n"
                 f"       Run pathway_analysis.py REVIGO mode, or pass "
                 f"--no-revigo.")
    df = pd.read_csv(csv_path)
    if "term_id" not in df.columns:
        sys.exit(f"ERROR: REVIGO CSV missing 'term_id' column: {csv_path}")
    return set(df["term_id"].astype(str))


# =============================================================================
# --- Top-N selection --------------------------------------------------------
# =============================================================================
def select_top(
    df: pd.DataFrame,
    top_n: int,
    revigo_ids: set[str] | None,
) -> pd.DataFrame:
    """Return top-N terms ordered by FDR (most significant first).

    Ordering is by significance rather than fold enrichment: fold enrichment
    collapses to effective_domain_size/query_size (a single value) for every
    term the query fully covers (recall = 1), so ordering by it would front-
    load a block of tied small terms and flatten the color encoding. Ordering
    by FDR keeps the headline pathways and lets fold enrichment (the color
    metric) vary as an orthogonal effect-size signal. Fold enrichment ties
    break the FDR order. Optionally restricts to REVIGO representatives first
    (GO:BP).
    """
    sub = df.copy()
    if revigo_ids is not None:
        sub = sub[sub["term_id"].astype(str).isin(revigo_ids)]
    sub = sub.sort_values(
        ["fdr_pvalue", "fold_enrichment"], ascending=[True, False]
    ).head(top_n).reset_index(drop=True)
    return sub


# =============================================================================
# --- Figure construction ----------------------------------------------------
# =============================================================================
def build_figure(
    df_dir: pd.DataFrame,
    db_label: str,
    x_floor: float = DEFAULT_X_FLOOR,
    x_ceiling: float = DEFAULT_X_CEILING,
) -> go.Figure:
    """Build the single-direction (upregulated) ORA lollipop figure.

    df_dir contains the selected top-N terms, pre-sorted by FDR ascending
    (most significant first). Encoding and cross-figure standardization mirror
    plot_gsea_lollipops.build_figure, with fold enrichment as the color metric
    and intersection count as the size metric.
    """
    n_terms = len(df_dir)
    if n_terms == 0:
        sys.exit("ERROR: no terms to plot.")

    df_dir = df_dir.copy()
    # Top of figure (highest y) = first row (highest fold enrichment).
    df_dir["__y"] = list(range(n_terms, 0, -1))

    colorscale = KETAMINE_SEQUENTIAL
    direction_label_text = "Up in ketamine"

    # Color values = fold enrichment; cmin/cmax span the panel-local range so
    # the full coral ramp is used regardless of absolute magnitude.
    color_values = df_dir["fold_enrichment"].astype(float)
    cmin = float(color_values.min())
    cmax = float(color_values.max())
    if cmax - cmin < 1e-6:
        cmax = cmin + 1e-3

    # Marker area sized to intersection count, calibrated against
    # SHARED_MAX_COUNT (not the per-figure max) so the same count renders at the
    # same diameter across all three ORA lollipops.
    size_values = df_dir["count"].clip(lower=1).astype(float)
    sizeref = 2.0 * SHARED_MAX_COUNT / (MAX_MARKER_SIZE ** 2)

    x_axis_min = x_floor
    x_axis_max = x_ceiling

    # Warn (don't auto-expand) if a term sits past the shared ceiling.
    x_max_observed = float(df_dir["neg_log10_fdr"].max())
    if x_max_observed > x_axis_max:
        n_clipped = int((df_dir["neg_log10_fdr"] > x_axis_max).sum())
        print(f"  WARNING: {n_clipped} term(s) have -log10(FDR) > "
              f"x_ceiling={x_axis_max:.2f} (max observed {x_max_observed:.2f}). "
              f"Consider --x-ceiling {np.ceil(x_max_observed) + 0.5:.1f} "
              f"(apply consistently across all ORA figures).")

    fig = go.Figure()

    # --- Stems --------------------------------------------------------------
    for _, row in df_dir.iterrows():
        fig.add_shape(
            type="line",
            x0=x_axis_min, x1=row["neg_log10_fdr"],
            y0=row["__y"], y1=row["__y"],
            line=dict(color=STEM_COLOR, width=1.2),
            layer="below",
        )

    # --- Markers ------------------------------------------------------------
    fig.add_trace(
        go.Scatter(
            x=df_dir["neg_log10_fdr"],
            y=df_dir["__y"],
            mode="markers",
            marker=dict(
                size=size_values,
                sizemode="area",
                sizeref=sizeref,
                sizemin=4,
                color=color_values,
                colorscale=colorscale,
                cmin=cmin,
                cmax=cmax,
                line=dict(color=DARK_ACCENT, width=0.6),
                colorbar=dict(
                    title=dict(
                        text="Fold<br>enrichment",
                        font=dict(family=FONT_FAMILY,
                                  size=LEGEND_TEXT_SIZE,
                                  color=FONT_COLOR),
                        side="right",
                    ),
                    tickfont=dict(family=FONT_FAMILY,
                                  size=TICK_LABEL_SIZE,
                                  color=FONT_COLOR),
                    thickness=14,
                    len=0.55,
                    y=0.72,
                    yanchor="middle",
                    x=1.02,
                    xanchor="left",
                    outlinewidth=0,
                ),
            ),
            customdata=df_dir[
                ["term_id", "term_name_display", "fold_enrichment",
                 "fdr_pvalue", "count", "term_size"]
            ].values,
            hovertemplate=(
                "<b>%{customdata[1]}</b><br>"
                "ID: %{customdata[0]}<br>"
                "Fold enrichment: %{customdata[2]:.2f}<br>"
                "FDR: %{customdata[3]:.2e}<br>"
                "Intersection count: %{customdata[4]}<br>"
                "Term size: %{customdata[5]}<extra></extra>"
            ),
            showlegend=False,
        )
    )

    # --- Intersection count legend (fixed reference sizes, shared) ---------
    for n in SHARED_LEGEND_SIZES:
        fig.add_trace(
            go.Scatter(
                x=[None], y=[None], mode="markers",
                marker=dict(
                    size=[float(n)],
                    sizemode="area",
                    sizeref=sizeref,
                    sizemin=4,
                    color=WHITE,
                    line=dict(color=DARK_ACCENT, width=0.6),
                ),
                name=f"{n}",
                showlegend=True,
                hoverinfo="skip",
            )
        )

    # --- Y-axis tick labels -------------------------------------------------
    yticks = df_dir["__y"].tolist()
    ylabels = df_dir["term_name_display"].tolist()

    # --- FDR=0.05 reference line (shared across all ORA lollipops) ---------
    if x_axis_min <= FDR_REFLINE_X <= x_axis_max:
        fig.add_shape(
            type="line",
            x0=FDR_REFLINE_X, x1=FDR_REFLINE_X,
            y0=0.3, y1=n_terms + 1.3,
            line=dict(color=REFERENCE_LINE_COLOR, width=1.2, dash="dash"),
            layer="below",
        )
        fig.add_annotation(
            x=FDR_REFLINE_X, y=0,
            xref="x", yref="paper",
            text="FDR = 0.05",
            showarrow=False,
            font=dict(family=FONT_FAMILY, size=DATA_LABEL_SIZE,
                      color=REFERENCE_LINE_COLOR),
            xanchor="left", yanchor="bottom",
            xshift=4, yshift=4,
        )

    # --- Direction subheading (top of figure, centered over plot area) -----
    plot_center_x = (x_axis_min + x_axis_max) / 2
    fig.add_annotation(
        x=plot_center_x,
        y=n_terms + 0.9,
        xref="x", yref="y",
        text=f"<b>{direction_label_text}</b>",
        showarrow=False,
        font=dict(family=FONT_FAMILY, size=SUBHEADING_SIZE,
                  color=FONT_COLOR),
        align="center",
    )

    # --- Layout (NO on-figure title; Section 11.2) --------------------------
    fig.update_layout(
        font=dict(family=FONT_FAMILY, size=AXIS_TITLE_SIZE, color=FONT_COLOR),
        plot_bgcolor=WHITE,
        paper_bgcolor=WHITE,
        width=1000,
        height=max(420, 36 * n_terms + 160),
        margin=dict(l=420, r=180, t=50, b=70),
        xaxis=dict(
            title=dict(
                text="−log<sub>10</sub>(FDR)",
                font=dict(family=FONT_FAMILY, size=AXIS_TITLE_SIZE,
                          color=FONT_COLOR),
            ),
            tickfont=dict(family=FONT_FAMILY, size=TICK_LABEL_SIZE,
                          color=FONT_COLOR),
            range=[x_axis_min, x_axis_max],
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
            tickfont=dict(family=FONT_FAMILY, size=TICK_LABEL_SIZE,
                          color=FONT_COLOR),
            ticklabelstandoff=10,
            range=[0.3, n_terms + 1.7],
            showgrid=False,
            zeroline=False,
            showline=False,
            ticks="",
        ),
        legend=dict(
            font=dict(family=FONT_FAMILY, size=LEGEND_TEXT_SIZE,
                      color=FONT_COLOR),
            title=dict(
                text="<b>Intersection<br>gene count</b>",
                font=dict(family=FONT_FAMILY, size=LEGEND_TEXT_SIZE,
                          color=FONT_COLOR),
            ),
            x=1.02, xanchor="left",
            y=0.10, yanchor="bottom",
            bgcolor="rgba(255,255,255,0.85)",
            bordercolor=SEPARATOR_COLOR,
            borderwidth=0.5,
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
            "Generate single-direction (upregulated) ORA lollipop thesis "
            "figures for one of GO:BP, KEGG, or Reactome, mirroring the GSEA "
            "lollipops (plot_gsea_lollipops.py) in style and naming."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--db", choices=list(DB_CONFIG.keys()), required=True,
        help="Which ORA database to plot.",
    )
    parser.add_argument(
        "--input", type=Path, default=None,
        help=("Override ORA _full results CSV path. Default inferred from "
              "--db as results/pathway_analysis/upregulated/"
              "upregulated_{token}_full.csv"),
    )
    parser.add_argument(
        "--top-n", type=int, default=15,
        help="Number of top terms shown in the figure.",
    )
    parser.add_argument(
        "--x-floor", type=float, default=DEFAULT_X_FLOOR,
        help="Left edge of the -log10(FDR) axis.",
    )
    parser.add_argument(
        "--x-ceiling", type=float, default=DEFAULT_X_CEILING,
        help=("Right edge of the -log10(FDR) axis. Shared across all three ORA "
              "lollipops so they are directly comparable."),
    )
    parser.add_argument(
        "--no-revigo", action="store_true",
        help=("Disable REVIGO simplification for GO:BP and use raw significant "
              "ORA terms instead."),
    )
    parser.add_argument(
        "--revigo-csv", type=Path, default=DEFAULT_REVIGO_CSV,
        help="REVIGO representatives CSV (GO:BP only).",
    )
    parser.add_argument(
        "--outdir", type=Path, default=DEFAULT_OUTDIR,
        help="Figure output directory.",
    )
    parser.add_argument(
        "--png-scale", type=int, default=5,
        help="Plotly scale factor for PNG export (~600 DPI at default size).",
    )
    return parser.parse_args(argv)


def main() -> None:
    args = parse_args()

    db_cfg = DB_CONFIG[args.db]
    use_revigo = db_cfg["supports_revigo"] and not args.no_revigo

    if args.input is None:
        args.input = default_input_csv(db_cfg["token"])

    print(f"Reading: {args.input}")
    df = load_ora_full(args.input)
    print(f"  rows total (significant ORA terms): {len(df)}")

    # --- REVIGO step (GO:BP only) ------------------------------------------
    revigo_ids = None
    if use_revigo:
        print(f"REVIGO simplification enabled (GO:BP): {args.revigo_csv}")
        revigo_ids = load_revigo_representative_ids(args.revigo_csv)
        print(f"  REVIGO representatives: {len(revigo_ids)}")
    elif db_cfg["supports_revigo"]:
        print("REVIGO simplification disabled by --no-revigo.")
    else:
        print(f"REVIGO not applicable for {db_cfg['label']} (GO-only tool).")

    # --- Top-N selection + figure ------------------------------------------
    df_top = select_top(df, args.top_n, revigo_ids)
    if len(df_top) == 0:
        sys.exit("ERROR: no terms selected; nothing to plot.")

    print(f"Top {len(df_top)} terms (by FDR, most significant first):")
    for _, r in df_top.iterrows():
        print(f"  FDR={r['fdr_pvalue']:.2e}  FE={r['fold_enrichment']:5.2f}  "
              f"count={int(r['count'])}  {r['term_name_display']}")

    fig = build_figure(
        df_top, db_cfg["label"],
        x_floor=args.x_floor, x_ceiling=args.x_ceiling,
    )
    basename = f"ORA_{db_cfg['short']}_lollipop_up_thesis"
    export_figure(fig, args.outdir, basename, args.png_scale)


if __name__ == "__main__":
    main()

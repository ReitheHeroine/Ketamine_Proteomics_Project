#!/usr/bin/env python3
"""
================================================================================
plot_gsea_lollipops.py

Title:         GSEA enrichment lollipop figures (thesis) - one direction per figure
Project:       Ketamine Astrocyte Proteomics
Author:        Reina Hastings (reinahastings13@gmail.com)
Date created:  2026-05-29
Last modified: 2026-05-29

Purpose:
    Generate publication-grade single-direction lollipop figures
    summarising GSEA results for one of three databases (GO:BP, KEGG,
    Reactome). Each invocation produces separate up- and down-direction
    figures (unless --direction restricts to one).

    Visual encoding (per figure):
      y-axis ordering:   terms ranked by |NES| descending, so the most
                         pronounced pathway is at the top.
      x-axis:            -log10(FDR), so the most statistically robust
                         enrichments extend furthest right.
      marker size:       setSize (number of annotated genes in the
                         pathway), area-scaled.
      marker color:      sequential ramp by |NES| - ketamine coral for
                         up-direction figures, control blue for down.
                         Reinforces the y-axis ranking and keeps NES
                         information visually accessible without the
                         narrow-band readability problem of a diverging
                         color encoding.

    For GO:BP, REVIGO simplification (Supek et al., 2011) is applied per
    direction before top-N selection so semantically redundant terms
    (e.g., 'establishment of localization in cell' vs 'cellular
    localization') don't dominate the figure. REVIGO is GO-only and is
    therefore skipped for KEGG and Reactome. REVIGO outputs are cached
    under results/gsea/revigo/ so re-runs of the plot do not re-fetch.

    Style follows project_notes/figure_and_table_style_guide.Rmd: Arial,
    14/12/10 pt hierarchy, sentence case, black interior text, project
    palette, no on-figure title.

Inputs:
    --db {go_bp, kegg, reactome}: which database to plot.
    --direction {up, down, both}: which direction(s) to generate.
                                  'both' (default) writes two files.
    --input (optional): GSEA results CSV; default inferred from
        --db, --min-size, --metric as
        results/gsea/{min_size}/{metric}/{db}_gsea_results.csv
    Required CSV columns: ID, Description, NES, p.adjust, passes_fdr,
        direction, setSize.

    REVIGO cache (auto-generated for GO:BP):
        results/gsea/revigo/go_bp_up_revigo.tsv
        results/gsea/revigo/go_bp_down_revigo.tsv

Outputs (results/figures/ by default):
    GSEA_GO_BP_lollipop_up_thesis.{pdf,png,html}
    GSEA_GO_BP_lollipop_down_thesis.{pdf,png,html}
    GSEA_KEGG_lollipop_up_thesis.{pdf,png,html}
    GSEA_KEGG_lollipop_down_thesis.{pdf,png,html}
    GSEA_Reactome_lollipop_up_thesis.{pdf,png,html}
    GSEA_Reactome_lollipop_down_thesis.{pdf,png,html}

Usage examples:
    python scripts/plot_gsea_lollipops.py --db go_bp
    python scripts/plot_gsea_lollipops.py --db kegg --top-n 12
    python scripts/plot_gsea_lollipops.py --db reactome --direction up
    python scripts/plot_gsea_lollipops.py --db go_bp --no-revigo
    python scripts/plot_gsea_lollipops.py --db go_bp --refresh-revigo

    copy/paste (all three DBs, both directions): \\
        for db in go_bp kegg reactome; do \\
            python scripts/plot_gsea_lollipops.py --db $db; \\
        done

Dependencies:
    pandas, numpy, plotly, kaleido (PNG export), requests (transitive via
    revigo_fetch.py).

================================================================================
"""

from __future__ import annotations

# --- Imports ----------------------------------------------------------------
import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go


# =============================================================================
# --- Project paths ----------------------------------------------------------
# =============================================================================
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUTDIR = PROJECT_ROOT / "results/figures/GSEA"
DEFAULT_REVIGO_DIR = PROJECT_ROOT / "results/gsea/revigo"
REVIGO_FETCH_SCRIPT = PROJECT_ROOT / "scripts/revigo_fetch.py"


# =============================================================================
# --- Database configuration -------------------------------------------------
# =============================================================================
DB_CONFIG = {
    "go_bp":    {"label": "GO Biological Process", "short": "GO_BP",
                  "supports_revigo": True},
    "kegg":     {"label": "KEGG",                  "short": "KEGG",
                  "supports_revigo": False},
    "reactome": {"label": "Reactome",              "short": "Reactome",
                  "supports_revigo": False},
}


# =============================================================================
# --- Style guide constants --------------------------------------------------
#   Anchored in project_notes/figure_and_table_style_guide.Rmd:
#     - Section 2:    Arial inside figures.
#     - Section 3:    manuscript type sizes.
#     - Section 9.1:  project palette.
#     - Section 9.1.4: ALL figure-interior text in pure black #000000.
#     - Section 11.2: no on-figure title.
# =============================================================================
FONT_FAMILY = "Arial"
FONT_COLOR = "#000000"
AXIS_TITLE_SIZE = 14
TICK_LABEL_SIZE = 12
LEGEND_TEXT_SIZE = 12
DATA_LABEL_SIZE = 10
SUBHEADING_SIZE = 12

KETAMINE_COLOR = "#E8735A"       # Up-direction accent.
CONTROL_COLOR = "#1F77B4"        # Down-direction accent (matches
                                 # gsea_analysis.R COLOR_DOWN).
DARK_ACCENT = "#2C3E50"          # Marker outlines.
SEPARATOR_COLOR = "#D0D0D0"
REFERENCE_LINE_COLOR = "#666666"
GRID_COLOR = "#E0E0E0"
STEM_COLOR = "#B0B0B0"
WHITE = "#FFFFFF"

# --- Cross-figure scaling constants -----------------------------------------
# These make all six GSEA lollipop figures (3 DBs x 2 directions) directly
# visually comparable: same x-axis range, same marker-size scale, same
# legend reference dots. A reader can lay any two figures side by side and
# the positions and dot sizes mean the same thing.
#
# X-axis (-log10(FDR)) defaults:
#   floor 1.0 sits just below -log10(0.05)=1.301 (the significance cutoff,
#   which is drawn as a dashed reference line). Ceiling 10.0 accommodates
#   the most-significant terms observed in this project (GO:BP up tops out
#   around 9.75). CLI flags --x-floor / --x-ceiling override per figure.
DEFAULT_X_FLOOR = 1.0
DEFAULT_X_CEILING = 10.0
FDR_THRESHOLD = 0.05
FDR_REFLINE_X = -np.log10(FDR_THRESHOLD)   # ~1.301

# Marker-size scaling: sizeref is computed against SHARED_MAX_COUNT
# rather than per-figure max, so e.g. a 30-leading-edge-gene pathway
# renders at the same diameter in every figure. Count (the leading-edge
# gene count from `core_enrichment`) is preferred over setSize for GSEA
# visualizations per the clusterProfiler/enrichplot convention because
# it reflects the data-driven enrichment signal rather than the
# database-driven pathway annotation. SHARED_MAX_COUNT = 100 covers the
# observed max in this project (GO:BP up: ~85-99 leading-edge genes for
# the largest pathways).
MAX_MARKER_SIZE = 32
SHARED_MAX_COUNT = 100
# Fixed legend reference sizes shown on every figure; chosen to span the
# observed leading-edge count range across all six lollipops
# (~3 to ~99) while remaining visually distinguishable.
SHARED_LEGEND_SIZES = [10, 30, 75]

# Single-hue sequential ramps. Up uses the project ketamine-coral anchor;
# down uses a desaturated-to-saturated blue ramp anchored at CONTROL_COLOR.
# Both ramps map low |NES| to a near-white tint and high |NES| to the
# saturated brand color so the most pronounced pathway reads visually
# darkest, matching the y-axis (top-of-figure) ranking.
KETAMINE_SEQUENTIAL = [
    [0.00, "#FCEDE8"],
    [0.25, "#F5B8A6"],
    [0.50, "#EF947D"],
    [0.75, "#EB825F"],
    [1.00, KETAMINE_COLOR],
]
CONTROL_SEQUENTIAL = [
    [0.00, "#EAF3FA"],
    [0.25, "#B5D2EA"],
    [0.50, "#7FB3D8"],
    [0.75, "#4A95C6"],
    [1.00, CONTROL_COLOR],
]


# =============================================================================
# --- Data loading -----------------------------------------------------------
# =============================================================================
def _sentence_case_term(s: str) -> str:
    """Sentence-case a term name while preserving biological prefixes.

    GO:BP / Reactome term strings start either lowercase (gseGO output) or
    with a lower-then-upper biological prefix (mRNA, miRNA, tRNA, snRNA,
    rRNA, snoRNA, pre-mRNA, ncRNA, etc.). Naive str.capitalize() would
    convert 'mRNA' to 'Mrna' and the s[0].upper()+s[1:] idiom would
    produce 'MRNA'. Both lose the biologically meaningful prefix casing.

    Rule: if the second character is already uppercase, the first
    character is part of a chemistry/biology prefix and is left as-is.
    Otherwise the first character is uppercased. Preserves internal
    capitalization (e.g. 'L1CAM interactions' is unchanged).
    """
    if not isinstance(s, str) or not s:
        return s
    if len(s) >= 2 and s[1].isupper():
        return s
    return s[0].upper() + s[1:]


def load_gsea(csv_path: Path) -> pd.DataFrame:
    """Load a GSEA results CSV; add neg_log10_padj, count, sentence-case names.

    Adds a `count` column = number of leading-edge genes (parsed from the
    slash-separated `core_enrichment` string). This is the GSEA "Count"
    metric used by clusterProfiler::dotplot for size encoding.
    """
    if not csv_path.exists():
        sys.exit(f"ERROR: GSEA results CSV not found: {csv_path}")
    df = pd.read_csv(csv_path)
    required = {"ID", "Description", "NES", "p.adjust", "passes_fdr",
                "direction", "setSize", "core_enrichment"}
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"ERROR: required columns missing in {csv_path}: {missing}")

    # passes_fdr is written as boolean by clusterProfiler->R but read by
    # pandas as object/string. Normalize to bool.
    df["passes_fdr"] = (
        df["passes_fdr"].astype(str).str.upper().isin(["TRUE", "T", "1"])
    )

    df["neg_log10_padj"] = -np.log10(df["p.adjust"].clip(lower=1e-300))
    df["term_name_display"] = df["Description"].apply(_sentence_case_term)

    # Leading-edge count: number of slash-separated gene symbols in
    # core_enrichment. Falls back to 0 for any (rare) missing rows so
    # downstream sizing math doesn't NaN out.
    df["count"] = (
        df["core_enrichment"].fillna("").apply(
            lambda s: len([g for g in s.split("/") if g.strip()])
        )
    )
    return df


def default_input_csv(db: str, min_size: int, metric: str) -> Path:
    """Construct default GSEA results CSV path from CLI parameters."""
    return (
        PROJECT_ROOT / "results/gsea" / f"min{min_size}" / metric
        / f"{db}_gsea_results.csv"
    )


# =============================================================================
# --- REVIGO integration (GO:BP only) ----------------------------------------
# =============================================================================
def _write_revigo_input(df_dir: pd.DataFrame, out_path: Path) -> None:
    """Write the (term_id, fdr_pvalue) CSV that revigo_fetch.py consumes."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df_dir.rename(
        columns={"ID": "term_id", "p.adjust": "fdr_pvalue"}
    )[["term_id", "fdr_pvalue"]].to_csv(out_path, index=False)


def _run_revigo_fetch(input_csv: Path, output_tsv: Path) -> None:
    """Invoke scripts/revigo_fetch.py as a subprocess."""
    cmd = [
        sys.executable,
        str(REVIGO_FETCH_SCRIPT),
        "--input", str(input_csv),
        "--output", str(output_tsv),
        # Mouse taxon, cutoff=0.7 ("medium"), SIMREL, BP namespace
        # are revigo_fetch.py defaults.
    ]
    print(f"  Running REVIGO fetch: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(
            f"ERROR: revigo_fetch.py failed (exit {result.returncode}).\n"
            f"  stdout: {result.stdout}\n"
            f"  stderr: {result.stderr}"
        )


def _parse_revigo_representatives(tsv_path: Path) -> set[str]:
    """Return the set of representative term IDs from a REVIGO TSV.

    Matches pathway_analysis.py::parse_revigo_output: representatives are
    rows whose Representative column is NaN (REVIGO leaves it empty for
    the chosen representative of each cluster).
    """
    if not tsv_path.exists():
        sys.exit(f"ERROR: REVIGO output TSV not found: {tsv_path}")
    df = pd.read_csv(tsv_path, sep="\t")
    if "Representative" not in df.columns or "TermID" not in df.columns:
        sys.exit(
            f"ERROR: REVIGO TSV missing required columns "
            f"(TermID, Representative): {tsv_path}"
        )
    return set(df[df["Representative"].isna()]["TermID"].astype(str))


def revigo_filter_ids(
    df_sig: pd.DataFrame,
    cache_dir: Path,
    db_short: str,
    refresh: bool,
) -> dict[str, set[str]]:
    """Run (or load cached) REVIGO per direction; return rep IDs per direction."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    out_by_direction: dict[str, set[str]] = {}

    for direction_label, direction_value in [
        ("up", "up_in_ketamine"),
        ("down", "down_in_ketamine"),
    ]:
        sub = df_sig[df_sig["direction"] == direction_value].copy()
        if sub.empty:
            print(f"  REVIGO: no significant {direction_label}-direction "
                  f"terms; skipping.")
            out_by_direction[direction_label] = set()
            continue

        base = f"{db_short.lower()}_{direction_label}"
        input_csv = cache_dir / f"{base}_input.csv"
        output_tsv = cache_dir / f"{base}_revigo.tsv"

        if refresh or not output_tsv.exists():
            _write_revigo_input(sub, input_csv)
            _run_revigo_fetch(input_csv, output_tsv)
        else:
            print(f"  REVIGO: using cached {output_tsv.name}")

        reps = _parse_revigo_representatives(output_tsv)
        print(f"  REVIGO {direction_label}: {len(reps)} representative "
              f"terms (from {len(sub)} significant)")
        out_by_direction[direction_label] = reps

    return out_by_direction


# =============================================================================
# --- Top-N selection (single direction) -------------------------------------
# =============================================================================
def select_top_for_direction(
    df: pd.DataFrame,
    direction: str,
    top_n: int,
    revigo_reps: dict[str, set[str]] | None,
) -> pd.DataFrame:
    """Return top-N significant terms for one direction, ordered by |NES| desc.

    direction is 'up' or 'down'. Filters to passes_fdr == True, restricts
    to the requested direction, optionally restricts to REVIGO
    representatives (GO:BP only), sorts by |NES| descending, then takes
    the top N.
    """
    dir_value = "up_in_ketamine" if direction == "up" else "down_in_ketamine"
    sub = df[df["passes_fdr"] & (df["direction"] == dir_value)].copy()
    if revigo_reps is not None and revigo_reps.get(direction):
        sub = sub[sub["ID"].isin(revigo_reps[direction])]
    sub["__abs_nes"] = sub["NES"].abs()
    # Tie-break by p.adjust ascending so two terms with identical |NES|
    # surface the more statistically robust one first.
    sub = sub.sort_values(
        ["__abs_nes", "p.adjust"], ascending=[False, True]
    ).head(top_n).reset_index(drop=True)
    return sub


# =============================================================================
# --- Figure construction ----------------------------------------------------
# =============================================================================
def build_figure(
    df_dir: pd.DataFrame,
    direction: str,
    db_label: str,
    x_floor: float = DEFAULT_X_FLOOR,
    x_ceiling: float = DEFAULT_X_CEILING,
) -> go.Figure:
    """Build a single-direction GSEA lollipop figure.

    df_dir contains the selected top-N terms for one direction, pre-sorted
    by |NES| descending (most pronounced first).

    Visual encoding:
      y position : row index in df_dir (top of figure = first row = highest |NES|)
      x position : -log10(FDR)
      size       : Count = leading-edge gene count (the clusterProfiler/
                   enrichplot convention for GSEA dotplots), area-scaled
      color      : sequential ramp by |NES|, coral for up / blue for down

    Cross-figure standardization:
      - x range is fixed at [x_floor, x_ceiling] (default [1.0, 10.0]) so
        all GSEA lollipops in this project share an x-axis ruler.
      - Marker sizeref is calibrated against SHARED_MAX_COUNT rather
        than the panel-local max, so a given leading-edge count renders
        at the same diameter in every figure.
      - A dashed vertical reference line at -log10(FDR_THRESHOLD)=1.301
        marks the FDR=0.05 significance cutoff.
    """
    n_terms = len(df_dir)
    if n_terms == 0:
        sys.exit(f"ERROR: no terms to plot for direction '{direction}'.")

    df_dir = df_dir.copy()
    # Top of figure (highest y) = first row (highest |NES|).
    df_dir["__y"] = list(range(n_terms, 0, -1))

    if direction == "up":
        colorscale = KETAMINE_SEQUENTIAL
        accent_color = KETAMINE_COLOR
        direction_label_text = "Up in ketamine"
    else:
        colorscale = CONTROL_SEQUENTIAL
        accent_color = CONTROL_COLOR
        direction_label_text = "Down in ketamine"

    # Color values (always positive); cmin/cmax set so the panel-local
    # |NES| range spans the full ramp regardless of absolute magnitude.
    color_values = df_dir["__abs_nes"].astype(float)
    cmin = float(color_values.min())
    cmax = float(color_values.max())
    if cmax - cmin < 1e-6:
        # All ties; nudge cmax so plotly doesn't divide by zero.
        cmax = cmin + 1e-3

    # Marker area sized to Count (leading-edge gene count). sizeref is
    # computed against SHARED_MAX_COUNT (not per-figure max) so the same
    # leading-edge count renders at the same pixel diameter in every
    # GSEA lollipop figure across the thesis.
    size_values = df_dir["count"].clip(lower=1).astype(float)
    sizeref = 2.0 * SHARED_MAX_COUNT / (MAX_MARKER_SIZE ** 2)

    # X-axis: fixed shared range. The right pad already lives in
    # x_ceiling; the floor sits just below FDR=0.05 so the reference
    # line at -log10(0.05) has a bit of breathing room.
    x_axis_min = x_floor
    x_axis_max = x_ceiling

    # Warn if any displayed term sits past the ceiling. Auto-expanding
    # would defeat the cross-figure-shared-scale purpose, so leave the
    # decision to the user.
    x_max_observed = float(df_dir["neg_log10_padj"].max())
    if x_max_observed > x_axis_max:
        n_clipped = int((df_dir["neg_log10_padj"] > x_axis_max).sum())
        print(f"  WARNING: {n_clipped} term(s) have -log10(FDR) > "
              f"x_ceiling={x_axis_max:.2f} (max observed "
              f"{x_max_observed:.2f}). Markers past the ceiling will "
              f"render at the right edge.")
        print(f"           Consider --x-ceiling "
              f"{np.ceil(x_max_observed) + 0.5:.1f} (apply consistently "
              f"across all figures for shared scaling).")

    fig = go.Figure()

    # --- Stems --------------------------------------------------------------
    # Stems begin at x_axis_min (the left edge of the visible plot area)
    # rather than at x=0 because the shared x-range starts above zero;
    # this avoids the stems being silently clipped by plotly.
    for _, row in df_dir.iterrows():
        fig.add_shape(
            type="line",
            x0=x_axis_min, x1=row["neg_log10_padj"],
            y0=row["__y"], y1=row["__y"],
            line=dict(color=STEM_COLOR, width=1.2),
            layer="below",
        )

    # --- Markers ------------------------------------------------------------
    fig.add_trace(
        go.Scatter(
            x=df_dir["neg_log10_padj"],
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
                        text="|NES|",
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
        )
    )

    # --- Leading-edge count legend (fixed reference sizes, shared) ---------
    # Using a fixed legend rather than per-panel sizes means a 30-gene
    # reference dot looks the same across every GSEA lollipop figure in
    # the thesis, which is the visual equivalent of a shared ruler.
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

    # --- FDR=0.05 reference line (shared across all GSEA lollipops) -------
    # All plotted terms are significant (passes_fdr) so they sit to the
    # right of this line. The line gives the reader a calibrated visual
    # anchor for "where significance starts" that is identical across
    # every figure.
    if x_axis_min <= FDR_REFLINE_X <= x_axis_max:
        fig.add_shape(
            type="line",
            x0=FDR_REFLINE_X, x1=FDR_REFLINE_X,
            y0=0.3, y1=n_terms + 1.3,
            line=dict(color=REFERENCE_LINE_COLOR, width=1.2, dash="dash"),
            layer="below",
        )
        # Label rendered in paper coordinates so it sits flush with the
        # bottom of the plot area regardless of the y-data range.
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
                text="<b>Leading-edge<br>gene count</b>",
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
# --- Export -----------------------------------------------------------------
# =============================================================================
def export_figure(fig: go.Figure, outdir: Path, basename: str,
                   png_scale: int) -> None:
    """Write PDF (vector), PNG (high DPI for Word), and HTML (interactive)."""
    outdir.mkdir(parents=True, exist_ok=True)
    pdf_path = outdir / f"{basename}.pdf"
    png_path = outdir / f"{basename}.png"
    html_path = outdir / f"{basename}.html"

    fig.write_image(pdf_path)
    fig.write_image(png_path, scale=png_scale)
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
            "Generate single-direction GSEA lollipop thesis figures "
            "(one figure per direction) for one of GO:BP, KEGG, or "
            "Reactome, styled per project_notes/"
            "figure_and_table_style_guide.Rmd."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--db", choices=list(DB_CONFIG.keys()), required=True,
        help="Which GSEA database to plot.",
    )
    parser.add_argument(
        "--direction", choices=["up", "down", "both"], default="both",
        help="Which direction figure(s) to generate.",
    )
    parser.add_argument(
        "--input", type=Path, default=None,
        help=("Override GSEA results CSV path. Default inferred from "
              "--db, --min-size, --metric."),
    )
    parser.add_argument(
        "--min-size", type=int, default=10, choices=[10, 15],
        help="Which gsea_analysis.R min_size run to read.",
    )
    parser.add_argument(
        "--metric", choices=["log2fc", "signed_logp"], default="log2fc",
        help="Which ranking-metric subdirectory to read.",
    )
    parser.add_argument(
        "--top-n", type=int, default=15,
        help="Number of top terms shown in each single-direction figure.",
    )
    parser.add_argument(
        "--x-floor", type=float, default=DEFAULT_X_FLOOR,
        help=("Left edge of the -log10(FDR) axis. Default sits just below "
              "the FDR=0.05 cutoff so the reference line has visual room."),
    )
    parser.add_argument(
        "--x-ceiling", type=float, default=DEFAULT_X_CEILING,
        help=("Right edge of the -log10(FDR) axis. Default accommodates "
              "the most-significant terms in this project (GO:BP up tops "
              "out near 9.75). Shared across all six lollipop figures so "
              "they are directly visually comparable."),
    )
    parser.add_argument(
        "--no-revigo", action="store_true",
        help=("Disable REVIGO simplification for GO:BP and use raw "
              "significant GSEA terms instead."),
    )
    parser.add_argument(
        "--refresh-revigo", action="store_true",
        help=("Force re-fetch of REVIGO outputs even if cached TSVs exist "
              "under results/gsea/revigo/."),
    )
    parser.add_argument(
        "--revigo-dir", type=Path, default=DEFAULT_REVIGO_DIR,
        help="Directory for REVIGO input/output cache.",
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
        args.input = default_input_csv(args.db, args.min_size, args.metric)

    print(f"Reading: {args.input}")
    df = load_gsea(args.input)
    n_sig = int(df["passes_fdr"].sum())
    print(f"  rows total: {len(df)} ; significant (FDR<=0.05): {n_sig}")

    # --- REVIGO step (GO:BP only) ------------------------------------------
    revigo_reps = None
    if use_revigo:
        print("REVIGO simplification enabled (GO:BP).")
        revigo_reps = revigo_filter_ids(
            df_sig=df[df["passes_fdr"]],
            cache_dir=args.revigo_dir,
            db_short=db_cfg["short"],
            refresh=args.refresh_revigo,
        )
    elif db_cfg["supports_revigo"]:
        print("REVIGO simplification disabled by --no-revigo.")
    else:
        print(f"REVIGO not applicable for {db_cfg['label']} (GO-only tool).")

    # --- Per-direction generation ------------------------------------------
    directions = (["up", "down"] if args.direction == "both"
                  else [args.direction])

    for direction in directions:
        print(f"\n--- {direction.upper()} direction ---")
        df_dir = select_top_for_direction(
            df, direction, args.top_n, revigo_reps
        )
        if len(df_dir) == 0:
            print(f"  No significant {direction}-direction terms; "
                  f"skipping figure.")
            continue

        print(f"Top {len(df_dir)} {direction}-direction terms "
              f"(by |NES| descending):")
        for _, r in df_dir.iterrows():
            sign = "+" if direction == "up" else "-"
            print(f"  {sign} NES={r['NES']:+.2f}  FDR={r['p.adjust']:.2e}  "
                  f"setSize={r['setSize']}  {r['term_name_display']}")

        fig = build_figure(
            df_dir, direction, db_cfg["label"],
            x_floor=args.x_floor, x_ceiling=args.x_ceiling,
        )
        basename = (f"GSEA_{db_cfg['short']}_lollipop_"
                    f"{direction}_thesis")
        export_figure(fig, args.outdir, basename, args.png_scale)


if __name__ == "__main__":
    main()
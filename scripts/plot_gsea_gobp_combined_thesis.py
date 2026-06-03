#!/usr/bin/env python3
"""
================================================================================
plot_gsea_gobp_combined_thesis.py

Title:         Combined GSEA GO:BP lollipop figure (thesis) - up and down in one
               panel with a single diverging NES colorbar
Project:       Ketamine Astrocyte Proteomics
Author:        Reina Hastings (reinahastings13@gmail.com)
Date created:  2026-06-01
Last modified: 2026-06-01

Purpose:
    Combine the two single-direction GSEA GO:BP lollipop figures
    (GSEA_GO_BP_lollipop_up_thesis / _down_thesis) into ONE thesis panel.
    The most upregulated pathway sits at the very top of the figure and the
    most downregulated at the very bottom, with terms ordered by SIGNED NES
    descending so the color gradient runs monotonically coral -> blue down
    the page.

    This differs from plot_gsea_lollipops.py in two deliberate ways:
      1. Color encodes SIGNED NES on a single diverging scale (saturated
         ketamine coral at the most positive NES, pale neutral at NES = 0,
         control blue at the most negative NES), shown as one long colorbar
         shared by both directions. The per-direction single-hue |NES| ramps
         of the standalone figures are replaced by this one continuous ruler.
      2. Up and down terms live in a single stacked column separated by a
         small gap and a faint horizontal rule, rather than in two files.

    Shared visual encoding (unchanged from plot_gsea_lollipops.py so the
    panel stays comparable to the other GSEA lollipops):
      y position : signed-NES rank (top = most positive, bottom = most negative)
      x position : -log10(FDR), shared fixed range, FDR=0.05 reference line
      marker size: leading-edge gene count (Count), area-scaled against a
                   shared maximum so a given count renders at the same
                   diameter across every GSEA lollipop in the thesis
      marker color: signed NES on the diverging coral<->blue scale

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

    REVIGO cache (auto-generated/reused for GO:BP):
        results/gsea/revigo/go_bp_up_revigo.tsv
        results/gsea/revigo/go_bp_down_revigo.tsv

Outputs (results/figures/thesis/ by default):
    GSEA_GO_BP_lollipop_combined_thesis.{pdf,png,html}

Usage examples:
    python scripts/plot_gsea_gobp_combined_thesis.py
    python scripts/plot_gsea_gobp_combined_thesis.py --up-n 15 --down-n 5
    python scripts/plot_gsea_gobp_combined_thesis.py --up-n 20 --down-n 10
    python scripts/plot_gsea_gobp_combined_thesis.py --no-revigo
    python scripts/plot_gsea_gobp_combined_thesis.py --refresh-revigo

    copy/paste (defaults, from project root):
        python scripts/plot_gsea_gobp_combined_thesis.py

Dependencies:
    pandas, numpy, plotly, kaleido (PNG/PDF export), requests (transitive via
    revigo_fetch.py). Run with the project's conda env, e.g.
        /Users/reina/miniconda3/envs/ketamine_project/bin/python
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
DEFAULT_OUTDIR = PROJECT_ROOT / "results/figures/thesis"
DEFAULT_REVIGO_DIR = PROJECT_ROOT / "results/gsea/revigo"
REVIGO_FETCH_SCRIPT = PROJECT_ROOT / "scripts/revigo_fetch.py"


# =============================================================================
# --- Style guide constants --------------------------------------------------
#   Anchored in project_notes/figure_and_table_style_guide.Rmd:
#     - Section 2:     Arial inside figures.
#     - Section 3:     manuscript type sizes.
#     - Section 9.1:   project palette.
#     - Section 9.1.4: ALL figure-interior text in pure black #000000.
#     - Section 11.2:  no on-figure title.
#   Mirrors plot_gsea_lollipops.py so this panel stays visually consistent.
# =============================================================================
FONT_FAMILY = "Arial"
FONT_COLOR = "#000000"
AXIS_TITLE_SIZE = 14
TICK_LABEL_SIZE = 12
LEGEND_TEXT_SIZE = 12
DATA_LABEL_SIZE = 10
SUBHEADING_SIZE = 12

KETAMINE_COLOR = "#E8735A"       # Positive-NES (up) anchor.
CONTROL_COLOR = "#1F77B4"        # Negative-NES (down) anchor (matches
                                 # gsea_analysis.R COLOR_DOWN).
DARK_ACCENT = "#2C3E50"          # Marker outlines.
SEPARATOR_COLOR = "#D0D0D0"
REFERENCE_LINE_COLOR = "#666666"
GRID_COLOR = "#E0E0E0"
STEM_COLOR = "#B0B0B0"
WHITE = "#FFFFFF"

# --- Cross-figure scaling constants -----------------------------------------
# Identical to plot_gsea_lollipops.py so this combined panel shares the same
# x-axis ruler and marker-size scale as the standalone lollipops.
DEFAULT_X_FLOOR = 1.0
DEFAULT_X_CEILING = 10.0
FDR_THRESHOLD = 0.05
FDR_REFLINE_X = -np.log10(FDR_THRESHOLD)   # ~1.301

MAX_MARKER_SIZE = 32
SHARED_MAX_COUNT = 100
SHARED_LEGEND_SIZES = [10, 30, 75]

# --- Diverging NES colorscale -----------------------------------------------
# One continuous ruler for SIGNED NES. With a symmetric cmin/cmax centered on
# zero, NES = 0 lands at the 0.5 stop (pale neutral), the most positive NES at
# 1.0 (ketamine coral) and the most negative NES at 0.0 (control blue). The
# intermediate stops reuse the tints from the two single-hue ramps in
# plot_gsea_lollipops.py so the coral and blue arms match the standalone
# figures.
NES_DIVERGING = [
    [0.00, CONTROL_COLOR],   # most negative NES -> saturated blue
    [0.25, "#7FB3D8"],
    [0.45, "#D6E6F2"],
    [0.50, "#F7F7F7"],       # NES = 0 -> pale neutral
    [0.55, "#FBE2D9"],
    [0.75, "#EF947D"],
    [1.00, KETAMINE_COLOR],  # most positive NES -> saturated coral
]


# =============================================================================
# --- Data loading -----------------------------------------------------------
# =============================================================================
def _sentence_case_term(s: str) -> str:
    """Sentence-case a term name while preserving biological prefixes.

    If the second character is already uppercase, the first character is part
    of a chemistry/biology prefix (mRNA, miRNA, tRNA, ...) and is left as-is;
    otherwise the first character is uppercased. Internal capitalization is
    preserved. (Same rule as plot_gsea_lollipops.py.)
    """
    if not isinstance(s, str) or not s:
        return s
    if len(s) >= 2 and s[1].isupper():
        return s
    return s[0].upper() + s[1:]


def load_gsea(csv_path: Path) -> pd.DataFrame:
    """Load a GSEA results CSV; add neg_log10_padj, count, sentence-case names.

    `count` = number of leading-edge genes parsed from the slash-separated
    `core_enrichment` string (the clusterProfiler/enrichplot size metric).
    """
    if not csv_path.exists():
        sys.exit(f"ERROR: GSEA results CSV not found: {csv_path}")
    df = pd.read_csv(csv_path)
    required = {"ID", "Description", "NES", "p.adjust", "passes_fdr",
                "direction", "setSize", "core_enrichment"}
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"ERROR: required columns missing in {csv_path}: {missing}")

    df["passes_fdr"] = (
        df["passes_fdr"].astype(str).str.upper().isin(["TRUE", "T", "1"])
    )
    df["neg_log10_padj"] = -np.log10(df["p.adjust"].clip(lower=1e-300))
    df["term_name_display"] = df["Description"].apply(_sentence_case_term)
    df["count"] = (
        df["core_enrichment"].fillna("").apply(
            lambda s: len([g for g in s.split("/") if g.strip()])
        )
    )
    return df


def default_input_csv(min_size: int, metric: str) -> Path:
    """Construct default GO:BP GSEA results CSV path from CLI parameters."""
    return (
        PROJECT_ROOT / "results/gsea" / f"min{min_size}" / metric
        / "go_bp_gsea_results.csv"
    )


# =============================================================================
# --- REVIGO integration (cached TSVs; GO:BP only) ---------------------------
#   Mirrors the REVIGO handling in plot_gsea_lollipops.py: reuse the cached
#   per-direction TSVs unless --refresh-revigo is passed.
# =============================================================================
def _write_revigo_input(df_dir: pd.DataFrame, out_path: Path) -> None:
    """Write the (term_id, fdr_pvalue) CSV that revigo_fetch.py consumes."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df_dir.rename(
        columns={"ID": "term_id", "p.adjust": "fdr_pvalue"}
    )[["term_id", "fdr_pvalue"]].to_csv(out_path, index=False)


def _run_revigo_fetch(input_csv: Path, output_tsv: Path) -> None:
    """Invoke scripts/revigo_fetch.py as a subprocess (defaults: mouse, 0.7)."""
    cmd = [
        sys.executable, str(REVIGO_FETCH_SCRIPT),
        "--input", str(input_csv),
        "--output", str(output_tsv),
    ]
    print(f"  Running REVIGO fetch: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(
            f"ERROR: revigo_fetch.py failed (exit {result.returncode}).\n"
            f"  stdout: {result.stdout}\n  stderr: {result.stderr}"
        )


def _parse_revigo_representatives(tsv_path: Path) -> set[str]:
    """Return the set of representative term IDs from a REVIGO TSV.

    Representatives are rows whose Representative column is NaN (REVIGO leaves
    it empty for the chosen representative of each cluster).
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

        base = f"go_bp_{direction_label}"
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
# --- Top-N selection --------------------------------------------------------
# =============================================================================
def select_top_for_direction(
    df: pd.DataFrame,
    direction: str,
    top_n: int,
    revigo_reps: dict[str, set[str]] | None,
) -> pd.DataFrame:
    """Return top-N significant terms for one direction, ordered by |NES| desc.

    Selection is by |NES| (the most pronounced enrichments per direction),
    matching plot_gsea_lollipops.py. The combined figure later re-orders the
    union by SIGNED NES for the stacked layout, but the per-direction TOP-N
    membership is identical to the standalone figures.
    """
    dir_value = "up_in_ketamine" if direction == "up" else "down_in_ketamine"
    sub = df[df["passes_fdr"] & (df["direction"] == dir_value)].copy()
    if revigo_reps is not None and revigo_reps.get(direction):
        sub = sub[sub["ID"].isin(revigo_reps[direction])]
    sub["__abs_nes"] = sub["NES"].abs()
    sub = sub.sort_values(
        ["__abs_nes", "p.adjust"], ascending=[False, True]
    ).head(top_n).reset_index(drop=True)
    return sub


# =============================================================================
# --- Figure construction ----------------------------------------------------
# =============================================================================
def build_combined_figure(
    df_up: pd.DataFrame,
    df_down: pd.DataFrame,
    x_floor: float = DEFAULT_X_FLOOR,
    x_ceiling: float = DEFAULT_X_CEILING,
    block_gap: float = 1.4,
) -> go.Figure:
    """Build the combined up+down GSEA GO:BP lollipop figure.

    Layout (top -> bottom): up terms (positive NES) then down terms (negative
    NES), each block internally sorted by signed NES descending, separated by
    `block_gap` blank y-units carrying a faint horizontal rule and the two
    direction labels. Color encodes signed NES on NES_DIVERGING with a
    symmetric range centered at 0.
    """
    n_up = len(df_up)
    n_down = len(df_down)
    if n_up + n_down == 0:
        sys.exit("ERROR: no terms to plot.")

    # --- Assemble the stacked y-coordinate column --------------------------
    #   Down block occupies y = 1 .. n_down (bottom). A gap follows. Up block
    #   occupies the top. Within each block, signed NES descending means the
    #   highest y carries the most positive / least negative NES, so the
    #   whole column is monotonic in signed NES from top to bottom.
    up = df_up.sort_values("NES", ascending=False).reset_index(drop=True)
    down = df_down.sort_values("NES", ascending=False).reset_index(drop=True)

    # Down block: top of block (highest y in block) = least negative NES.
    down = down.copy()
    down["__y"] = list(range(n_down, 0, -1))                  # n_down .. 1
    # Up block: sits above the gap; highest y = most positive NES.
    up = up.copy()
    up_base = n_down + block_gap
    up["__y"] = [up_base + (n_up - i) for i in range(n_up)]   # top..bottom

    plotted = pd.concat([up, down], ignore_index=True)
    plotted["__abs_nes"] = plotted["NES"].abs()

    y_top = up["__y"].max() if n_up else down["__y"].max()
    y_bottom = 1.0
    boundary_y = n_down + block_gap / 2.0                     # gap center

    # --- Symmetric diverging color range -----------------------------------
    #   Centered at 0 so coral and blue arms are balanced; the most extreme
    #   |NES| sets both ends.
    cabs = float(plotted["__abs_nes"].max())
    cabs = max(cabs, 1e-3)
    cmin, cmax = -cabs, cabs

    # --- Marker sizing (shared scale with the standalone lollipops) --------
    size_values = plotted["count"].clip(lower=1).astype(float)
    sizeref = 2.0 * SHARED_MAX_COUNT / (MAX_MARKER_SIZE ** 2)

    x_axis_min, x_axis_max = x_floor, x_ceiling

    # Warn (do not auto-expand) if a term sits past the shared ceiling.
    x_max_observed = float(plotted["neg_log10_padj"].max())
    if x_max_observed > x_axis_max:
        n_clipped = int((plotted["neg_log10_padj"] > x_axis_max).sum())
        print(f"  WARNING: {n_clipped} term(s) have -log10(FDR) > "
              f"x_ceiling={x_axis_max:.2f} (max observed "
              f"{x_max_observed:.2f}). Consider --x-ceiling "
              f"{np.ceil(x_max_observed) + 0.5:.1f}.")

    fig = go.Figure()

    # --- Stems --------------------------------------------------------------
    for _, row in plotted.iterrows():
        fig.add_shape(
            type="line",
            x0=x_axis_min, x1=row["neg_log10_padj"],
            y0=row["__y"], y1=row["__y"],
            line=dict(color=STEM_COLOR, width=1.2),
            layer="below",
        )

    # --- Markers (signed NES color) ----------------------------------------
    fig.add_trace(
        go.Scatter(
            x=plotted["neg_log10_padj"],
            y=plotted["__y"],
            mode="markers",
            marker=dict(
                size=size_values,
                sizemode="area",
                sizeref=sizeref,
                sizemin=4,
                color=plotted["NES"].astype(float),
                colorscale=NES_DIVERGING,
                cmin=cmin,
                cmax=cmax,
                line=dict(color=DARK_ACCENT, width=0.6),
                colorbar=dict(
                    title=dict(
                        text="NES",
                        font=dict(family=FONT_FAMILY, size=LEGEND_TEXT_SIZE,
                                  color=FONT_COLOR),
                        side="right",
                    ),
                    tickfont=dict(family=FONT_FAMILY, size=TICK_LABEL_SIZE,
                                  color=FONT_COLOR),
                    # Symmetric ticks through zero so the sign crossing reads
                    # clearly on the one long colorbar.
                    tickmode="array",
                    tickvals=_symmetric_ticks(cabs),
                    thickness=16,
                    len=0.62,
                    y=0.72,
                    yanchor="middle",
                    x=1.02,
                    xanchor="left",
                    outlinewidth=0,
                ),
            ),
            customdata=plotted[
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

    # --- FDR=0.05 reference line -------------------------------------------
    if x_axis_min <= FDR_REFLINE_X <= x_axis_max:
        fig.add_shape(
            type="line",
            x0=FDR_REFLINE_X, x1=FDR_REFLINE_X,
            y0=y_bottom - 0.7, y1=y_top + 0.7,
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

    # --- Up/down separator rule (in the gap) -------------------------------
    if n_up and n_down:
        fig.add_shape(
            type="line",
            x0=x_axis_min, x1=x_axis_max,
            y0=boundary_y, y1=boundary_y,
            line=dict(color=SEPARATOR_COLOR, width=1.0, dash="dot"),
            layer="below",
        )

    # --- Direction labels ---------------------------------------------------
    plot_center_x = (x_axis_min + x_axis_max) / 2.0
    if n_up:
        fig.add_annotation(
            x=plot_center_x, y=y_top + 0.9,
            xref="x", yref="y",
            text="<b>Up in ketamine</b>",
            showarrow=False,
            font=dict(family=FONT_FAMILY, size=SUBHEADING_SIZE,
                      color=KETAMINE_COLOR),
            align="center",
        )
    if n_down:
        # Sit just below the separator, introducing the down block.
        fig.add_annotation(
            x=plot_center_x, y=boundary_y - 0.15,
            xref="x", yref="y",
            text="<b>Down in ketamine</b>",
            showarrow=False,
            font=dict(family=FONT_FAMILY, size=SUBHEADING_SIZE,
                      color=CONTROL_COLOR),
            align="center",
            yanchor="top",
        )

    # --- Axes / layout (NO on-figure title; Section 11.2) ------------------
    yticks = plotted["__y"].tolist()
    ylabels = plotted["term_name_display"].tolist()
    n_rows_equiv = n_up + n_down + block_gap

    fig.update_layout(
        font=dict(family=FONT_FAMILY, size=AXIS_TITLE_SIZE, color=FONT_COLOR),
        plot_bgcolor=WHITE,
        paper_bgcolor=WHITE,
        width=1000,
        height=int(max(480, 36 * n_rows_equiv + 170)),
        margin=dict(l=440, r=190, t=60, b=70),
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
            range=[y_bottom - 1.2, y_top + 1.8],
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


def _symmetric_ticks(cabs: float) -> list[float]:
    """Return integer-ish symmetric colorbar ticks spanning [-cabs, +cabs].

    Picks a tick step so there are roughly 4-6 ticks each side including 0,
    keeping the long colorbar legible without crowding.
    """
    step = 1.0 if cabs <= 3.0 else 2.0
    ticks = [0.0]
    v = step
    while v <= cabs + 1e-9:
        ticks = [-v] + ticks + [v]
        v += step
    return ticks


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
            "Combine the up- and down-direction GSEA GO:BP lollipops into one "
            "thesis panel (most upregulated at top, most downregulated at "
            "bottom) sharing a single diverging NES colorbar and one "
            "leading-edge-count size legend. Styled per project_notes/"
            "figure_and_table_style_guide.Rmd."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input", type=Path, default=None,
        help=("Override GO:BP GSEA results CSV path. Default inferred from "
              "--min-size and --metric."),
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
        "--up-n", type=int, default=15,
        help="Number of top upregulated (positive NES) terms.",
    )
    parser.add_argument(
        "--down-n", type=int, default=5,
        help="Number of top downregulated (negative NES) terms.",
    )
    parser.add_argument(
        "--x-floor", type=float, default=DEFAULT_X_FLOOR,
        help="Left edge of the -log10(FDR) axis.",
    )
    parser.add_argument(
        "--x-ceiling", type=float, default=DEFAULT_X_CEILING,
        help=("Right edge of the -log10(FDR) axis. Shared with the standalone "
              "lollipops so the panels are directly comparable."),
    )
    parser.add_argument(
        "--block-gap", type=float, default=1.4,
        help="Blank y-units between the up and down blocks (carries the rule).",
    )
    parser.add_argument(
        "--no-revigo", action="store_true",
        help="Disable REVIGO simplification; use raw significant terms.",
    )
    parser.add_argument(
        "--refresh-revigo", action="store_true",
        help="Force re-fetch of REVIGO outputs even if cached TSVs exist.",
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
        "--basename", type=str, default="GSEA_GO_BP_lollipop_combined_thesis",
        help="Output file basename (no extension).",
    )
    parser.add_argument(
        "--png-scale", type=int, default=5,
        help="Plotly scale factor for PNG export (~600 DPI at default size).",
    )
    return parser.parse_args(argv)


def main() -> None:
    args = parse_args()
    use_revigo = not args.no_revigo

    if args.input is None:
        args.input = default_input_csv(args.min_size, args.metric)

    print(f"Reading: {args.input}")
    df = load_gsea(args.input)
    n_sig = int(df["passes_fdr"].sum())
    print(f"  rows total: {len(df)} ; significant (FDR<=0.05): {n_sig}")

    # --- REVIGO step --------------------------------------------------------
    revigo_reps = None
    if use_revigo:
        print("REVIGO simplification enabled (GO:BP).")
        revigo_reps = revigo_filter_ids(
            df_sig=df[df["passes_fdr"]],
            cache_dir=args.revigo_dir,
            refresh=args.refresh_revigo,
        )
    else:
        print("REVIGO simplification disabled by --no-revigo.")

    # --- Per-direction top-N selection -------------------------------------
    df_up = select_top_for_direction(df, "up", args.up_n, revigo_reps)
    df_down = select_top_for_direction(df, "down", args.down_n, revigo_reps)

    if len(df_up) < args.up_n:
        print(f"  NOTE: only {len(df_up)} up-direction terms available "
              f"(requested {args.up_n}).")
    if len(df_down) < args.down_n:
        print(f"  NOTE: only {len(df_down)} down-direction terms available "
              f"(requested {args.down_n}).")

    print(f"\nTop {len(df_up)} up-direction terms (signed-NES order, "
          f"top of figure):")
    for _, r in df_up.sort_values("NES", ascending=False).iterrows():
        print(f"  + NES={r['NES']:+.2f}  FDR={r['p.adjust']:.2e}  "
              f"count={r['count']}  {r['term_name_display']}")
    print(f"\nTop {len(df_down)} down-direction terms (signed-NES order, "
          f"bottom of figure):")
    for _, r in df_down.sort_values("NES", ascending=False).iterrows():
        print(f"  - NES={r['NES']:+.2f}  FDR={r['p.adjust']:.2e}  "
              f"count={r['count']}  {r['term_name_display']}")

    fig = build_combined_figure(
        df_up, df_down,
        x_floor=args.x_floor, x_ceiling=args.x_ceiling,
        block_gap=args.block_gap,
    )
    export_figure(fig, args.outdir, args.basename, args.png_scale)


if __name__ == "__main__":
    main()

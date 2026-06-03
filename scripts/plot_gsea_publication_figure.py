#!/usr/bin/env python3
"""
================================================================================
plot_gsea_publication_figure.py

Title:         GSEA GO:BP publication-style composite figure (thesis) - one
               direction per figure
Project:       Ketamine Astrocyte Proteomics
Author:        Reina Hastings (reinahastings13@gmail.com)
Date created:  2026-05-29
Last modified: 2026-05-29

Purpose:
    Build the GSEA GO:BP publication-style composite figure for a single
    direction (up- or down-in-ketamine). Mirrors the layout of the ORA
    pathway_analysis.py::create_publication_figure but uses the GSEA
    leading-edge subset (`core_enrichment`) for the per-pathway gene
    heatmaps.

    Layout (3-row plotly subplot grid, one direction per figure):
      Row 1, Panel A:   top 5 pathways for the chosen direction, as a
                        column heatmap colored by |NES| on a sequential
                        ramp (ketamine coral for up, control blue for
                        down). Y-ordering: |NES| descending so the most
                        pronounced pathway sits at the top.
      Row 2, Panels B-D: per-pathway leading-edge gene heatmaps for the
                        top-3 pathways (cells colored by log2FC on a
                        blue-white-coral diverging scale).
      Row 3, Panels E-F: leading-edge gene heatmaps for pathways 4 and 5.

    Terms are pre-filtered through REVIGO (Supek et al., 2011) so the top
    pathways are not semantically redundant, then further deduplicated by
    leading-edge Jaccard overlap so two pathways with near-identical gene
    sets do not both appear (e.g., 'protein-DNA complex assembly' vs
    'protein-DNA complex organization').

    Style follows project_notes/figure_and_table_style_guide.Rmd: Arial,
    14/12/10 pt, sentence case, black interior text, no on-figure title,
    project palette.

Inputs:
    --direction {up, down, both} (default: both, writes two files)
    --gsea-csv: GSEA GO:BP results CSV
        Default: results/gsea/min10/log2fc/go_bp_gsea_results.csv
        Required columns: ID, Description, NES, p.adjust, passes_fdr,
            direction, setSize, core_enrichment.
    --protein-csv: per-protein quantitative table providing log2FC
        Default: results/quantitative/all_quantitative_proteins.csv
        Required columns: Gene Symbol, log2_fold_change.

    REVIGO cache (shared with plot_gsea_lollipops.py):
        results/gsea/revigo/go_bp_{up,down}_revigo.tsv

Outputs (results/figures/ by default):
    GSEA_GO_BP_publication_figure_up_thesis.{pdf,png,html}
    GSEA_GO_BP_publication_figure_down_thesis.{pdf,png,html}

Usage examples:
    python scripts/plot_gsea_publication_figure.py
    python scripts/plot_gsea_publication_figure.py --direction up
    python scripts/plot_gsea_publication_figure.py --no-revigo
    python scripts/plot_gsea_publication_figure.py --refresh-revigo
    python scripts/plot_gsea_publication_figure.py --n-panels 5 \\
        --max-genes-per-panel 15

Dependencies:
    pandas, numpy, plotly, kaleido (PNG export); imports helpers from
    plot_gsea_lollipops.py (must live in the same scripts/ directory).

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

# Shared helpers from the lollipop script (REVIGO subprocess + parsing,
# loaders, colorscale ramps, style constants). Both scripts share the
# same REVIGO cache and color palette.
from plot_gsea_lollipops import (
    CONTROL_SEQUENTIAL,
    DEFAULT_REVIGO_DIR,
    KETAMINE_COLOR,
    KETAMINE_SEQUENTIAL,
    load_gsea,
    revigo_filter_ids,
)


# =============================================================================
# --- Project paths ----------------------------------------------------------
# =============================================================================
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_GSEA_CSV = (
    PROJECT_ROOT / "results/gsea/min10/log2fc/go_bp_gsea_results.csv"
)
DEFAULT_PROTEIN_CSV = (
    PROJECT_ROOT / "results/quantitative/all_quantitative_proteins.csv"
)
DEFAULT_OUTDIR = PROJECT_ROOT / "results/figures/GSEA"


# =============================================================================
# --- Style guide constants --------------------------------------------------
# =============================================================================
FONT_FAMILY = "Arial"
FONT_COLOR = "#000000"
AXIS_TITLE_SIZE = 14
TICK_LABEL_SIZE = 12
SUBHEADING_SIZE = 12
PANEL_LETTER_SIZE = 15

WHITE = "#FFFFFF"
CONTROL_COLOR_HEATMAP = "#7FB3D8"     # Light blue end of diverging gene
                                       # heatmap scale; matches the ORA
                                       # composite.

# Single-direction colorscales for per-pathway gene heatmaps (log2FC).
# Multi-stop ramps with extended endpoints and a tinted floor (not pure
# white) so:
#   1. Small-magnitude values (the bulk of leading-edge gene log2FCs in
#      this dataset) get a visible tint rather than blending into pure
#      white.
#   2. The dark end extends past the brand-anchor color, giving extreme
#      values more saturation headroom and pulling perceptual
#      differentiation across the entire range.
#   3. Multiple intermediate stops create more discriminable color steps
#      than a 2-stop linear RGB interpolation, which is not perceptually
#      uniform.
SINGLE_DIRECTION_UP = [    # 0 -> +max; pale tint -> deep coral
    [0.00, "#FFF5F0"],     # near-white with faint warm tint
    [0.20, "#F5C5B5"],     # light coral
    [0.40, "#EF947D"],
    [0.60, KETAMINE_COLOR],  # #E8735A (brand anchor)
    [0.80, "#C95642"],
    [1.00, "#9A3C2D"],     # deep brick red - extreme value saturation
]
SINGLE_DIRECTION_DOWN = [  # -max -> 0; deep navy -> pale tint
    [0.00, "#0E3A5F"],     # deep navy - extreme value saturation
    [0.20, "#1F77B4"],     # control blue (brand anchor; matches
                            # gsea_analysis.R COLOR_DOWN)
    [0.40, "#4A95C6"],
    [0.60, "#7FB3D8"],
    [0.80, "#B5D2EA"],
    [1.00, "#F0F6FB"],     # near-white with faint cool tint
]
# Reversed version of SINGLE_DIRECTION_DOWN: pale tint at low position,
# deep navy at high position. Used together with negated z values (so the
# heatmap stores |log2FC| instead of log2FC) to flip the down-composite
# colorbar's reading direction. Plotly's colorbar always renders zmax at
# the top, so to display most-negative at the colorbar top we make
# |log2FC| (which is largest for the most-negative log2FC) the z, then
# customize tickvals/ticktext to display the signed log2FC value.
SINGLE_DIRECTION_DOWN_FLIPPED = [
    [0.00, "#F0F6FB"],     # pale tint at z=0 (no effect)
    [0.20, "#B5D2EA"],
    [0.40, "#7FB3D8"],
    [0.60, "#4A95C6"],
    [0.80, "#1F77B4"],
    [1.00, "#0E3A5F"],     # deep navy at z=|max| (most negative log2FC)
]

# log2FC color range bounds. zmin/zmax are computed from the 95th
# percentile of |log2FC| across all displayed leading-edge genes - using
# the percentile (rather than the absolute max) prevents a single
# extreme outlier (e.g., INA at log2FC ~3 in the up direction) from
# stretching the color range so wide that typical genes (~0.5-1.0)
# become indistinguishable. The outlier saturates at the deepest color;
# the bulk of genes get full perceptual range.
FC_RANGE_FLOOR = 0.4        # absolute minimum cap (so a panel with
                             # only tiny effects still gets some range)
FC_RANGE_CEILING = 4.0       # generous absolute maximum cap so the data
                             # max (e.g., INA at log2FC ~3 in up) is not
                             # clipped under the default 100th-percentile
                             # ("use data max") cap.
DEFAULT_FC_CAP_SKIP_TOP_N = 1   # by default the single most-extreme
                                 # gene (e.g., INA at log2FC ~3 in up,
                                 # TEX15 at ~-1.2 in down) saturates at
                                 # the colorbar end. The cap is set to
                                 # the 2nd-most-extreme value so the
                                 # remaining ~70+ leading-edge genes get
                                 # the full color range. Increase via
                                 # --fc-cap-skip-top-n to let more
                                 # outliers saturate; set to 0 to use
                                 # the absolute max (no saturation).


# =============================================================================
# --- Data loading -----------------------------------------------------------
# =============================================================================
def load_protein_fc(csv_path: Path) -> dict[str, float]:
    """Return {Gene Symbol -> log2_fold_change} from the quantitative table."""
    if not csv_path.exists():
        sys.exit(f"ERROR: protein CSV not found: {csv_path}")
    df = pd.read_csv(csv_path, usecols=["Gene Symbol", "log2_fold_change"])
    df = df.dropna(subset=["Gene Symbol", "log2_fold_change"])
    # If a gene appears more than once (e.g., aggregate rows with
    # semicolon-joined symbols), keep the row with the largest |log2FC|.
    df["__abs"] = df["log2_fold_change"].abs()
    df = df.sort_values("__abs", ascending=False).drop_duplicates("Gene Symbol")
    return dict(zip(df["Gene Symbol"], df["log2_fold_change"]))


def parse_leading_edge(core_enrichment: str) -> list[str]:
    """Split clusterProfiler's slash-separated core_enrichment string."""
    if not isinstance(core_enrichment, str) or not core_enrichment:
        return []
    return [g.strip() for g in core_enrichment.split("/") if g.strip()]


# =============================================================================
# --- Term selection ---------------------------------------------------------
# =============================================================================
def _greedy_jaccard_dedupe(
    df_sorted: pd.DataFrame, n: int, threshold: float
) -> pd.DataFrame:
    """Greedy selection of top-N terms with leading-edge Jaccard <= threshold.

    REVIGO deduplicates by GO-ontology semantic distance, but for the
    publication composite we additionally want to avoid picking two terms
    whose leading-edge gene sets are near-identical (e.g., 'protein-DNA
    complex assembly' and 'protein-DNA complex organization' share the
    same 7-gene leading edge). Walks df_sorted top to bottom and keeps a
    term only if its leading edge has Jaccard <= threshold against every
    already-selected term.
    """
    selected_rows: list[pd.Series] = []
    selected_le: list[set[str]] = []
    for _, row in df_sorted.iterrows():
        le = set(parse_leading_edge(row.get("core_enrichment", "")))
        if not le:
            continue
        max_j = 0.0
        for prev in selected_le:
            inter = len(le & prev)
            union = len(le | prev) or 1
            max_j = max(max_j, inter / union)
        if max_j <= threshold:
            selected_rows.append(row)
            selected_le.append(le)
            if len(selected_rows) >= n:
                break
    if not selected_rows:
        return df_sorted.head(0).reset_index(drop=True)
    return pd.DataFrame(selected_rows).reset_index(drop=True)


def select_top_for_direction(
    df: pd.DataFrame,
    direction: str,
    n: int,
    revigo_reps: dict[str, set[str]] | None,
    jaccard_threshold: float,
) -> pd.DataFrame:
    """Top-N significant terms for one direction, ordered by |NES| desc.

    Pipeline:
      1. Filter to passes_fdr == True and direction == requested.
      2. (Optional) restrict to REVIGO representative IDs.
      3. Sort by |NES| descending; tie-break by p.adjust ascending.
      4. Greedy leading-edge Jaccard dedup to remove near-duplicate
         pathways (e.g., 'X assembly' vs 'X organization').
    """
    dir_value = "up_in_ketamine" if direction == "up" else "down_in_ketamine"
    sub = df[df["passes_fdr"] & (df["direction"] == dir_value)].copy()
    if revigo_reps is not None and revigo_reps.get(direction):
        sub = sub[sub["ID"].isin(revigo_reps[direction])]
    sub["__abs_nes"] = sub["NES"].abs()
    sub = sub.sort_values(
        ["__abs_nes", "p.adjust"], ascending=[False, True]
    )
    if jaccard_threshold >= 1.0:
        return sub.head(n).reset_index(drop=True)
    return _greedy_jaccard_dedupe(sub, n, jaccard_threshold)


# =============================================================================
# --- Helpers ----------------------------------------------------------------
# =============================================================================
def wrap_label(text: str, width: int = 25) -> str:
    """Soft-wrap a long term name onto two lines for subplot titles."""
    if not isinstance(text, str) or len(text) <= width:
        return text
    words = text.split()
    line1, line2 = "", ""
    for w in words:
        candidate = (line1 + " " + w).strip()
        if len(candidate) <= width and not line2:
            line1 = candidate
        else:
            line2 = (line2 + " " + w).strip()
    return f"{line1}<br>{line2}" if line2 else line1


# =============================================================================
# --- Figure construction ----------------------------------------------------
# =============================================================================
def build_publication_figure(
    df_dir: pd.DataFrame,
    direction: str,
    gene_fc: dict[str, float],
    n_panels: int,
    max_genes_per_panel: int,
    fc_cap_skip_top_n: int = DEFAULT_FC_CAP_SKIP_TOP_N,
) -> go.Figure:
    """Assemble the single-direction 3-row publication composite.

    df_dir : top pathways for the chosen direction, sorted by |NES| desc.
             At most n_panels are used for gene heatmaps; Panel A includes
             all rows in df_dir for context.
    """
    if df_dir.empty:
        sys.exit(f"ERROR: no terms passed filtering for direction "
                 f"'{direction}'; nothing to plot.")

    # --- Subplot specs ------------------------------------------------------
    # Row 1: Panel A spanning 3 columns.
    # Rows 2-3: up to 5 pathway gene heatmaps (B-F), 3 in row 2 and 2 in
    # row 3 (3rd slot of row 3 left empty so the layout matches the ORA
    # composite exactly).
    specs = [
        [{"colspan": 3, "type": "heatmap"}, None, None],
        [{"type": "heatmap"}] * 3,
        [{"type": "heatmap"}] * 3,
    ]

    pathway_rows = df_dir.head(n_panels)
    n_pathway_panels = len(pathway_rows)

    # --- Subplot titles -----------------------------------------------------
    direction_subhead = ("Up in ketamine" if direction == "up"
                         else "Down in ketamine")
    subplot_titles = [
        f"<b>Top biological functions ({direction_subhead.lower()})</b>"
    ]
    # Pad to 7 (Panel A + 3 in row 2 + 3 in row 3) so make_subplots gets
    # the expected count regardless of how many pathway panels are filled.
    for i in range(6):
        if i < n_pathway_panels:
            subplot_titles.append(
                f"<b>{wrap_label(pathway_rows.iloc[i]['term_name_display'], 25)}</b>"
            )
        else:
            subplot_titles.append("")

    fig = make_subplots(
        rows=3, cols=3,
        row_heights=[0.30, 0.35, 0.35],
        specs=specs,
        subplot_titles=subplot_titles,
        vertical_spacing=0.14,
        horizontal_spacing=0.10,
    )

    # --- Panel A: top biological functions ---------------------------------
    # |NES| on a single-direction sequential ramp. Y-axis ordering puts the
    # most pronounced pathway at the top of the rendered heatmap (Plotly
    # y-axis is bottom-to-top, so we reverse the lists for plotting).
    ordered_names = df_dir["term_name_display"].tolist()
    ordered_abs_nes = df_dir["__abs_nes"].tolist()
    ordered_padj = df_dir["p.adjust"].tolist()
    ordered_nes = df_dir["NES"].tolist()

    z_values = [[abs_nes] for abs_nes in reversed(ordered_abs_nes)]
    y_labels = list(reversed(ordered_names))
    customdata_a = list(reversed([
        [nes, padj] for nes, padj in zip(ordered_nes, ordered_padj)
    ]))

    abs_nes_lo = min(ordered_abs_nes)
    abs_nes_hi = max(ordered_abs_nes)
    if abs_nes_hi - abs_nes_lo < 1e-6:
        abs_nes_hi = abs_nes_lo + 1e-3

    colorscale_a = (KETAMINE_SEQUENTIAL if direction == "up"
                    else CONTROL_SEQUENTIAL)

    panel_a_y_domain = fig.layout.yaxis.domain
    panel_a_top = panel_a_y_domain[1]
    panel_a_height = panel_a_y_domain[1] - panel_a_y_domain[0]

    fig.add_trace(
        go.Heatmap(
            z=z_values,
            y=y_labels,
            x=["Ketamine vs Control"],
            colorscale=colorscale_a,
            zmin=abs_nes_lo,
            zmax=abs_nes_hi,
            customdata=customdata_a,
            showscale=True,
            colorbar=dict(
                title=dict(
                    text="|NES|",
                    font=dict(size=AXIS_TITLE_SIZE, color=FONT_COLOR,
                              family=FONT_FAMILY),
                ),
                tickfont=dict(size=TICK_LABEL_SIZE, color=FONT_COLOR,
                              family=FONT_FAMILY),
                x=1.02,
                len=panel_a_height,
                y=panel_a_top,
                yanchor="top",
                thickness=15,
            ),
            hovertemplate=(
                "<b>%{y}</b><br>NES: %{customdata[0]:+.2f}<br>"
                "FDR: %{customdata[1]:.2e}<extra></extra>"
            ),
        ),
        row=1, col=1,
    )

    # --- Shared (whole-figure) log2FC color range --------------------------
    # All gene heatmaps in this figure use the SAME zmin/zmax so the
    # color scale is directly comparable between panels and only ONE
    # colorbar appears on the figure. The cap is the requested
    # percentile of |log2FC| across all displayed leading-edge genes
    # (default 95th) so the most extreme few saturate at the deepest
    # color, freeing the rest of the color range for the bulk of genes.
    # Panels with smaller-effect leading edges will visibly blend more
    # than panels with larger effects - this is the honest signal.
    same_sign_fcs: list[float] = []
    for _, term_row in pathway_rows.iterrows():
        for g in parse_leading_edge(term_row["core_enrichment"])[
            :max_genes_per_panel
        ]:
            fc = gene_fc.get(g)
            if fc is None or not np.isfinite(fc):
                continue
            if (direction == "up" and fc > 0) or \
               (direction == "down" and fc < 0):
                same_sign_fcs.append(float(fc))
    if same_sign_fcs:
        abs_fcs = sorted((abs(v) for v in same_sign_fcs), reverse=True)
        # Skip the top N extreme values; the (N+1)-th most extreme
        # becomes the cap. Saturating outliers (those exceeding the cap)
        # are rendered at the deepest color but the cap itself sits
        # inside the bulk of the data, giving the remaining genes the
        # full color range. Clamp the skip index so we never overshoot
        # the list of available values.
        skip = min(fc_cap_skip_top_n, len(abs_fcs) - 1)
        rank_cap = abs_fcs[skip]
        fc_zlim = float(np.clip(rank_cap,
                                FC_RANGE_FLOOR, FC_RANGE_CEILING))
    else:
        fc_zlim = FC_RANGE_FLOOR

    if direction == "up":
        fc_zmin, fc_zmax = 0.0, fc_zlim
        gene_colorscale = SINGLE_DIRECTION_UP
    else:
        # For the down composite we render |log2FC| (not log2FC) on the
        # color scale and use the flipped colorscale; the actual signed
        # log2FC is preserved via customdata for hover and via the
        # custom colorbar tick labels. This is the workaround needed to
        # put the most-negative value at the TOP of the colorbar
        # (matching the gene order in each panel), because plotly's
        # colorbar always renders zmax at the top.
        fc_zmin, fc_zmax = 0.0, fc_zlim
        gene_colorscale = SINGLE_DIRECTION_DOWN_FLIPPED

    # --- Panels B-G: per-pathway leading-edge gene heatmaps ----------------
    # Layout:
    #   row 2 cols 1-3 = pathways 1-3 (panels B, C, D)
    #   row 3 cols 1-3 = pathways 4-6 (panels E, F, G)
    # 3+3 grid below Panel A so both rows are visually balanced.
    panel_grid = [
        (2, 1), (2, 2), (2, 3),
        (3, 1), (3, 2), (3, 3),
    ]

    for panel_idx in range(n_pathway_panels):
        if panel_idx >= len(panel_grid):
            break
        row, col = panel_grid[panel_idx]
        term_row = pathway_rows.iloc[panel_idx]
        leading_genes = parse_leading_edge(term_row["core_enrichment"])
        if not leading_genes:
            continue

        gene_values = []
        gene_names = []
        for g in leading_genes:
            fc = gene_fc.get(g)
            if fc is None or not np.isfinite(fc):
                continue
            gene_values.append([fc])
            gene_names.append(g)
            if len(gene_names) >= max_genes_per_panel:
                break

        if not gene_values:
            continue

        # Sort by |log2FC| ascending so the largest-magnitude effect ends
        # up at the top of the rendered cell stack (Plotly y is
        # bottom-to-top, so the last item in the list renders at the
        # top). This puts the most-upregulated gene at the top of up
        # panels and the most-downregulated gene at the top of down
        # panels - "most-affected at top" works symmetrically for both
        # directions and matches reader intuition.
        sorted_pairs = sorted(
            zip(gene_names, gene_values), key=lambda x: abs(x[1][0])
        )
        gene_names = [g for g, _ in sorted_pairs]
        gene_values = [v for _, v in sorted_pairs]

        # Gene symbols denote protein products in this proteomics figure;
        # style guide section 4.3 uses uppercase roman for protein-product
        # labels (matches the ORA composite).
        gene_names = [g.upper() for g in gene_names]

        is_last_panel = (panel_idx == n_pathway_panels - 1)

        # For down composite: render |log2FC| on the color scale (z is
        # forced positive) so the colorbar's top can display the most-
        # negative value. customdata carries the original signed
        # log2FC so the hover tooltip still shows real values. The
        # colorbar tick labels are overridden to display the signed
        # (negative) values matching the gene direction.
        if direction == "down":
            z_trace = [[abs(v[0])] for v in gene_values]
            customdata_trace = [[v[0]] for v in gene_values]
            hover_template = (
                "<b>%{y}</b><br>"
                "log<sub>2</sub>FC: %{customdata[0]:.2f}<extra></extra>"
            )
            # Build 5 evenly-spaced ticks on [0, fc_zmax] and label them
            # with the corresponding NEGATIVE log2FC values so the user
            # reads the colorbar in the natural "0 at bottom, most-
            # downregulated at top" direction.
            n_ticks = 5
            tick_positions = np.linspace(0, fc_zmax, n_ticks)
            tickvals = tick_positions.tolist()
            ticktext = [
                f"{-t:.2f}" if t > 0 else "0"
                for t in tick_positions
            ]
        else:
            z_trace = gene_values
            customdata_trace = None
            hover_template = (
                "<b>%{y}</b><br>"
                "log<sub>2</sub>FC: %{z:.2f}<extra></extra>"
            )
            tickvals = None
            ticktext = None

        colorbar_kwargs = dict(
            title=dict(
                text="log<sub>2</sub>FC",
                font=dict(size=AXIS_TITLE_SIZE, color=FONT_COLOR,
                          family=FONT_FAMILY),
            ),
            tickfont=dict(size=TICK_LABEL_SIZE, color=FONT_COLOR,
                          family=FONT_FAMILY),
            x=1.02,
            len=0.28,
            y=0.17,
            thickness=15,
        )
        if tickvals is not None:
            colorbar_kwargs["tickmode"] = "array"
            colorbar_kwargs["tickvals"] = tickvals
            colorbar_kwargs["ticktext"] = ticktext

        fig.add_trace(
            go.Heatmap(
                z=z_trace,
                y=gene_names,
                x=["log<sub>2</sub>FC"],
                colorscale=gene_colorscale,
                zmin=fc_zmin,
                zmax=fc_zmax,
                showscale=is_last_panel,
                customdata=customdata_trace,
                colorbar=(colorbar_kwargs if is_last_panel else None),
                hovertemplate=hover_template,
            ),
            row=row, col=col,
        )

    # --- Layout -------------------------------------------------------------
    fig.update_layout(
        height=1300,
        width=1100,
        plot_bgcolor=WHITE,
        paper_bgcolor=WHITE,
        font=dict(family=FONT_FAMILY, size=TICK_LABEL_SIZE, color=FONT_COLOR),
        margin=dict(l=180, r=80, t=70, b=50),
    )
    fig.update_yaxes(
        tickfont=dict(family=FONT_FAMILY, size=TICK_LABEL_SIZE,
                      color=FONT_COLOR)
    )
    fig.update_xaxes(
        tickfont=dict(family=FONT_FAMILY, size=TICK_LABEL_SIZE,
                      color=FONT_COLOR)
    )

    # Subplot-title annotation styling: 12 pt bold black sentence case.
    for annotation in fig["layout"]["annotations"]:
        annotation["font"] = dict(
            family=FONT_FAMILY, size=SUBHEADING_SIZE, color=FONT_COLOR
        )

    # Panel letters (A-G) at the upper-left corner of each subplot domain.
    # Axis suffix sequence follows the make_subplots default:
    #   (1,1) -> no suffix; (2,1)-(2,3) -> 2-4; (3,1)-(3,3) -> 5-7.
    panel_axis_suffixes = ["", "2", "3", "4", "5", "6", "7"]
    panel_letters = ["A", "B", "C", "D", "E", "F", "G"]
    n_total_panels = 1 + n_pathway_panels
    for i in range(min(7, n_total_panels)):
        ax_suffix = panel_axis_suffixes[i]
        x_domain = fig.layout[f"xaxis{ax_suffix}"].domain
        y_domain = fig.layout[f"yaxis{ax_suffix}"].domain
        fig.add_annotation(
            x=x_domain[0] - 0.005,
            y=y_domain[1] + 0.015,
            xref="paper",
            yref="paper",
            text=f"<b>{panel_letters[i]}</b>",
            showarrow=False,
            font=dict(family=FONT_FAMILY, size=PANEL_LETTER_SIZE,
                      color=FONT_COLOR),
            xanchor="right",
            yanchor="bottom",
        )

    return fig


# =============================================================================
# --- Export -----------------------------------------------------------------
# =============================================================================
def export_figure(fig: go.Figure, outdir: Path, basename: str,
                   png_scale: int) -> None:
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
            "Generate the single-direction GSEA GO:BP publication-style "
            "composite figure (3-row layout, top-N |NES| Panel A + "
            "leading-edge gene heatmaps)."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--direction", choices=["up", "down", "both"], default="both",
        help="Which direction figure(s) to generate.",
    )
    parser.add_argument(
        "--gsea-csv", type=Path, default=DEFAULT_GSEA_CSV,
        help="GSEA GO:BP results CSV.",
    )
    parser.add_argument(
        "--protein-csv", type=Path, default=DEFAULT_PROTEIN_CSV,
        help="Per-protein quantitative table providing log2FC for genes.",
    )
    parser.add_argument(
        "--n-panels", type=int, default=6,
        help=("Number of pathway gene-heatmap panels per direction (max 6; "
              "the composite has 6 pathway slots = B,C,D,E,F,G arranged "
              "as 3+3 across rows 2 and 3)."),
    )
    parser.add_argument(
        "--max-genes-per-panel", type=int, default=15,
        help="Maximum leading-edge genes shown per pathway heatmap.",
    )
    parser.add_argument(
        "--jaccard-dedup", type=float, default=0.7,
        help=("Maximum leading-edge Jaccard between any two selected "
              "pathways per direction. Lower = stricter dedup. Set to "
              "1.0 to disable."),
    )
    parser.add_argument(
        "--fc-cap-skip-top-n", type=int,
        default=DEFAULT_FC_CAP_SKIP_TOP_N,
        help=("Number of top-extreme |log2FC| values to skip when "
              "setting the color-scale cap. Default 1 uses the "
              "2nd-most-extreme value as the cap so a single outlier "
              "(e.g., INA at log2FC ~3 in up, TEX15 at ~-1.2 in down) "
              "saturates at the deepest color and the bulk of genes "
              "get the full color range. 0 = use absolute max (no "
              "saturation). 2+ = let more outliers saturate."),
    )
    parser.add_argument(
        "--no-revigo", action="store_true",
        help="Disable REVIGO simplification; use raw GSEA terms.",
    )
    parser.add_argument(
        "--refresh-revigo", action="store_true",
        help="Force re-fetch REVIGO outputs even if cached.",
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
        "--basename-suffix", default="",
        help=("Optional suffix appended to the output basename (before "
              "the extension). Useful for keeping multiple variants of "
              "the same figure side by side, e.g., "
              "'--basename-suffix _fullrange' produces "
              "'GSEA_GO_BP_publication_figure_up_thesis_fullrange.pdf'."),
    )
    parser.add_argument(
        "--png-scale", type=int, default=3,
        help=("Plotly scale factor for PNG export. The composite is large "
              "(1100x1300); scale=3 keeps file size reasonable while still "
              "yielding ~300-400 DPI at Word-embed size."),
    )
    args = parser.parse_args(argv)
    if args.n_panels > 6:
        parser.error("--n-panels must be <= 6 (composite has 6 pathway "
                     "slots: B, C, D, E, F, G).")
    return args


def main() -> None:
    args = parse_args()

    print(f"Reading GSEA results: {args.gsea_csv}")
    df = load_gsea(args.gsea_csv)
    print(f"  rows total: {len(df)} ; significant (FDR<=0.05): "
          f"{int(df['passes_fdr'].sum())}")

    print(f"Reading protein log2FC: {args.protein_csv}")
    gene_fc = load_protein_fc(args.protein_csv)
    print(f"  gene -> log2FC entries: {len(gene_fc)}")

    # --- REVIGO step --------------------------------------------------------
    revigo_reps = None
    if not args.no_revigo:
        print("REVIGO simplification enabled (GO:BP).")
        revigo_reps = revigo_filter_ids(
            df_sig=df[df["passes_fdr"]],
            cache_dir=args.revigo_dir,
            db_short="GO_BP",
            refresh=args.refresh_revigo,
        )
    else:
        print("REVIGO simplification disabled by --no-revigo.")

    # --- Per-direction generation ------------------------------------------
    directions = (["up", "down"] if args.direction == "both"
                  else [args.direction])

    for direction in directions:
        print(f"\n--- {direction.upper()} direction ---")
        df_dir = select_top_for_direction(
            df, direction, args.n_panels, revigo_reps,
            jaccard_threshold=args.jaccard_dedup,
        )
        if df_dir.empty:
            print(f"  No pathways selected for {direction}; skipping.")
            continue

        print(f"Selected {len(df_dir)} {direction}-direction terms "
              f"(by |NES| desc, leading-edge Jaccard <= "
              f"{args.jaccard_dedup}):")
        for _, r in df_dir.iterrows():
            sign = "+" if direction == "up" else "-"
            print(f"  {sign} NES={r['NES']:+.2f}  FDR={r['p.adjust']:.2e}  "
                  f"setSize={r['setSize']}  {r['term_name_display']}")

        fig = build_publication_figure(
            df_dir, direction, gene_fc,
            args.n_panels, args.max_genes_per_panel,
            fc_cap_skip_top_n=args.fc_cap_skip_top_n,
        )
        basename = (f"GSEA_GO_BP_publication_figure_{direction}_thesis"
                    f"{args.basename_suffix}")
        export_figure(fig, args.outdir, basename, args.png_scale)


if __name__ == "__main__":
    main()
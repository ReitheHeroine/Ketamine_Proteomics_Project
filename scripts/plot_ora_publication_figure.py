#!/usr/bin/env python3
"""
================================================================================
plot_ora_publication_figure.py

Title:         ORA GO:BP publication-style composite figure (thesis) -
               upregulated, mirroring the GSEA composite
Project:       Ketamine Astrocyte Proteomics
Author:        Reina Hastings (reinahastings13@gmail.com)
Date created:  2026-06-01
Last modified: 2026-06-01

Purpose:
    Build the ORA GO:BP publication-style composite figure for the upregulated
    set, mirroring plot_gsea_publication_figure.py in layout, style, and file
    naming. Where the GSEA composite uses |NES| (Panel A) and the leading-edge
    gene subset (per-pathway gene heatmaps), this ORA composite uses fold
    enrichment (Panel A) and the g:Profiler intersection genes (the query
    proteins annotated to each term).

    Layout (3-row plotly subplot grid):
      Row 1, Panel A:   top pathways as a column heatmap colored by fold
                        enrichment on the ketamine-coral sequential ramp.
                        Y-ordering: FDR (most significant pathway at the top).
                        Fold enrichment drives color but not order - it
                        collapses to one value for fully-covered terms, so
                        FDR drives selection/ordering and fold enrichment
                        varies as an orthogonal color gradient.
      Row 2, Panels B-D: per-pathway intersection-gene heatmaps for the top-3
                        pathways (cells colored by log2FC on the single-
                        direction up ramp).
      Row 3, Panels E-G: intersection-gene heatmaps for pathways 4-6.

    Terms are pre-filtered through REVIGO (Supek et al., 2011) representatives,
    then further deduplicated by intersection-gene Jaccard overlap so two terms
    with near-identical gene sets do not both appear.

    Intersection genes that lack a quantitative log2FC (presence/absence
    ketamine-specific proteins) are omitted from the gene heatmaps, exactly as
    the GSEA composite omits genes missing from its log2FC map. Such genes are
    still counted in the ORA enrichment itself; they simply cannot be placed on
    a fold-change color scale. This caveat belongs in the thesis caption.

    Style follows project_notes/figure_and_table_style_guide.Rmd: Arial,
    14/12/10 pt, sentence case, black interior text, no on-figure title.

Inputs:
    --ora-csv: ORA GO:BP "_full" results CSV
        Default: results/pathway_analysis/upregulated/upregulated_GO_BP_full.csv
        Required columns: term_id, term_name, fdr_pvalue, term_size,
            query_size, intersection_size, effective_domain_size,
            intersection_genes (semicolon-joined).
    --protein-csv: per-protein table providing log2FC
        Default: results/pathway_analysis/all_proteins_categorized_expanded.csv
        Required columns: Gene Symbol, log2_fold_change, category.
        (This is the same expanded categorized table pathway_analysis.py used,
        so the gene universe matches the committed ORA results.) Only
        category == 'quantitative' rows are used for log2FC, so presence/absence
        ketamine-specific proteins - whose log2_fold_change is the sentinel
        log2(100) = 6.644 rather than a measurement - are excluded from the
        gene heatmaps (see the caveat below).
    --revigo-csv: REVIGO representatives CSV (GO:BP)
        Default: results/pathway_analysis/revigo/upregulated/
            upregulated_GO_BP_revigo.csv
        Required column: term_id.

Outputs (results/figures/ORA/ by default):
    ORA_GO_BP_publication_figure_up_thesis.{pdf,png,html}
    ORA_GO_BP_publication_figure_up_thesis_fullrange.{pdf,png,html}  (via
        --fc-cap-skip-top-n 0 --basename-suffix _fullrange)

Usage examples:
    python scripts/plot_ora_publication_figure.py
    python scripts/plot_ora_publication_figure.py --no-revigo
    python scripts/plot_ora_publication_figure.py --fc-cap-skip-top-n 0 \\
        --basename-suffix _fullrange

Dependencies:
    pandas, numpy, plotly, kaleido (PNG export); imports helpers from
    plot_gsea_publication_figure.py and plot_gsea_lollipops.py (must live in
    the same scripts/ directory).

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

# Shared helpers/constants from the GSEA composite + lollipop scripts. Reusing
# them keeps the ORA composite visually identical to the GSEA composite except
# for the encoded metric (fold enrichment vs |NES|) and the gene subset
# (intersection vs leading edge).
from plot_gsea_publication_figure import (
    AXIS_TITLE_SIZE,
    DEFAULT_FC_CAP_SKIP_TOP_N,
    FC_RANGE_CEILING,
    FC_RANGE_FLOOR,
    FONT_COLOR,
    FONT_FAMILY,
    PANEL_LETTER_SIZE,
    SINGLE_DIRECTION_UP,
    SUBHEADING_SIZE,
    TICK_LABEL_SIZE,
    WHITE,
    _greedy_jaccard_dedupe,
    export_figure,
    parse_leading_edge,
    wrap_label,
)
from plot_gsea_lollipops import KETAMINE_SEQUENTIAL, _sentence_case_term


# =============================================================================
# --- Project paths ----------------------------------------------------------
# =============================================================================
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_ORA_CSV = (
    PROJECT_ROOT
    / "results/pathway_analysis/upregulated/upregulated_GO_BP_full.csv"
)
DEFAULT_PROTEIN_CSV = (
    PROJECT_ROOT
    / "results/pathway_analysis/all_proteins_categorized_expanded.csv"
)
DEFAULT_REVIGO_CSV = (
    PROJECT_ROOT
    / "results/pathway_analysis/revigo/upregulated/upregulated_GO_BP_revigo.csv"
)
DEFAULT_OUTDIR = PROJECT_ROOT / "results/figures/ORA"


# =============================================================================
# --- Data loading -----------------------------------------------------------
# =============================================================================
def load_ora_full(csv_path: Path) -> pd.DataFrame:
    """Load the ORA _full GO:BP CSV; add fold enrichment, names, core_enrichment.

    The intersection_genes column (semicolon-joined) is re-exposed as a
    slash-joined `core_enrichment` column so the GSEA helpers
    (parse_leading_edge, _greedy_jaccard_dedupe) can be reused unchanged.
    """
    if not csv_path.exists():
        sys.exit(f"ERROR: ORA results CSV not found: {csv_path}")
    df = pd.read_csv(csv_path)
    required = {"term_id", "term_name", "fdr_pvalue", "term_size",
                "query_size", "intersection_size", "effective_domain_size",
                "intersection_genes"}
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"ERROR: required columns missing in {csv_path}: {missing}")

    observed = df["intersection_size"] / df["query_size"]
    expected = df["term_size"] / df["effective_domain_size"]
    df["fold_enrichment"] = observed / expected.replace(0, np.nan)
    df["term_name_display"] = df["term_name"].apply(_sentence_case_term)
    # Slash-join so parse_leading_edge ('/' splitter) works on intersection genes.
    df["core_enrichment"] = (
        df["intersection_genes"].fillna("").astype(str)
        .apply(lambda s: "/".join(g.strip() for g in s.split(";") if g.strip()))
    )
    return df


def load_quant_protein_fc(csv_path: Path) -> dict[str, float]:
    """Return {Gene Symbol -> log2_fold_change} for QUANTITATIVE proteins only.

    The GSEA composite draws log2FC from the quantitative-only table, so every
    gene there has a measured fold change. The ORA intersection, by contrast,
    can include presence/absence ketamine-specific proteins whose
    log2_fold_change is a sentinel (log2(100) = 6.644), not a measurement.
    Restricting the map to category == 'quantitative' keeps the gene heatmaps
    honest (measured changes only) and matches the GSEA composite. The omitted
    presence/absence genes still count toward the enrichment itself; they are
    simply not placeable on a measured-fold-change color scale.
    """
    if not csv_path.exists():
        sys.exit(f"ERROR: protein CSV not found: {csv_path}")
    df = pd.read_csv(
        csv_path, usecols=["Gene Symbol", "log2_fold_change", "category"]
    )
    df = df[df["category"] == "quantitative"]
    df = df.dropna(subset=["Gene Symbol", "log2_fold_change"])
    # If a gene appears more than once, keep the largest-|log2FC| row.
    df["__abs"] = df["log2_fold_change"].abs()
    df = df.sort_values("__abs", ascending=False).drop_duplicates("Gene Symbol")
    return dict(zip(df["Gene Symbol"], df["log2_fold_change"]))


def load_revigo_representative_ids(csv_path: Path) -> set[str]:
    """Return the set of REVIGO representative GO term IDs."""
    if not csv_path.exists():
        sys.exit(f"ERROR: REVIGO representatives CSV not found: {csv_path}")
    df = pd.read_csv(csv_path)
    if "term_id" not in df.columns:
        sys.exit(f"ERROR: REVIGO CSV missing 'term_id' column: {csv_path}")
    return set(df["term_id"].astype(str))


# =============================================================================
# --- Term selection ---------------------------------------------------------
# =============================================================================
def select_top(
    df: pd.DataFrame,
    n: int,
    revigo_ids: set[str] | None,
    jaccard_threshold: float,
) -> pd.DataFrame:
    """Top-N terms ordered by FDR (most significant first), Jaccard deduped.

    Ordering is by significance, not fold enrichment: fold enrichment collapses
    to a single value (effective_domain_size/query_size) for every fully-
    covered term (recall = 1), so ordering/selecting by it yields a block of
    tied terms and a flat Panel A. Ordering by FDR keeps the headline pathways
    and lets the fold-enrichment color encoding (Panel A) vary as a real
    gradient. Dedup is on the intersection gene set (carried on the reused
    `core_enrichment` column), mirroring the GSEA composite's leading-edge
    Jaccard dedup.
    """
    sub = df.copy()
    if revigo_ids is not None:
        sub = sub[sub["term_id"].astype(str).isin(revigo_ids)]
    sub = sub.sort_values(
        ["fdr_pvalue", "fold_enrichment"], ascending=[True, False]
    )
    if jaccard_threshold >= 1.0:
        return sub.head(n).reset_index(drop=True)
    return _greedy_jaccard_dedupe(sub, n, jaccard_threshold)


# =============================================================================
# --- Figure construction ----------------------------------------------------
# =============================================================================
def build_publication_figure(
    df_top: pd.DataFrame,
    gene_fc: dict[str, float],
    n_panels: int,
    max_genes_per_panel: int,
    fc_cap_skip_top_n: int = DEFAULT_FC_CAP_SKIP_TOP_N,
) -> go.Figure:
    """Assemble the upregulated 3-row ORA composite (mirrors the GSEA composite).

    df_top : top pathways sorted by FDR (most significant first). At most
             n_panels are used for gene heatmaps; Panel A includes all rows for
             context.
    """
    if df_top.empty:
        sys.exit("ERROR: no terms passed filtering; nothing to plot.")

    # --- Subplot specs: Panel A spanning 3 cols (row 1), 3+3 gene heatmaps ---
    specs = [
        [{"colspan": 3, "type": "heatmap"}, None, None],
        [{"type": "heatmap"}] * 3,
        [{"type": "heatmap"}] * 3,
    ]

    pathway_rows = df_top.head(n_panels)
    n_pathway_panels = len(pathway_rows)

    # --- Subplot titles -----------------------------------------------------
    subplot_titles = [
        "<b>Top biological functions (up in ketamine)</b>"
    ]
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

    # --- Panel A: top biological functions, FDR-ordered, FE-colored --------
    # Rows arrive FDR-sorted (most significant first). Plotly y-axis renders
    # bottom-to-top, so reverse the lists to put the most significant pathway
    # at the top of the rendered heatmap. Cell color = fold enrichment.
    ordered_names = df_top["term_name_display"].tolist()
    ordered_fe = df_top["fold_enrichment"].tolist()
    ordered_fdr = df_top["fdr_pvalue"].tolist()

    z_values = [[fe] for fe in reversed(ordered_fe)]
    y_labels = list(reversed(ordered_names))
    customdata_a = list(reversed([
        [fe, fdr] for fe, fdr in zip(ordered_fe, ordered_fdr)
    ]))

    fe_lo = min(ordered_fe)
    fe_hi = max(ordered_fe)
    if fe_hi - fe_lo < 1e-6:
        fe_hi = fe_lo + 1e-3

    panel_a_y_domain = fig.layout.yaxis.domain
    panel_a_top = panel_a_y_domain[1]
    panel_a_height = panel_a_y_domain[1] - panel_a_y_domain[0]

    fig.add_trace(
        go.Heatmap(
            z=z_values,
            y=y_labels,
            x=["Ketamine vs Control"],
            colorscale=KETAMINE_SEQUENTIAL,
            zmin=fe_lo,
            zmax=fe_hi,
            customdata=customdata_a,
            showscale=True,
            colorbar=dict(
                title=dict(
                    text="Fold<br>enrichment",
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
                "<b>%{y}</b><br>Fold enrichment: %{customdata[0]:.2f}<br>"
                "FDR: %{customdata[1]:.2e}<extra></extra>"
            ),
        ),
        row=1, col=1,
    )

    # --- Shared (whole-figure) log2FC color range --------------------------
    # All gene heatmaps share one zmin/zmax (and one colorbar). The cap is the
    # (fc_cap_skip_top_n + 1)-th most extreme |log2FC| across displayed
    # intersection genes, so a few outliers saturate while the bulk get the
    # full color range. ORA is upregulated-only, so only positive log2FCs are
    # color-mapped.
    pos_fcs: list[float] = []
    for _, term_row in pathway_rows.iterrows():
        for g in parse_leading_edge(term_row["core_enrichment"])[
            :max_genes_per_panel
        ]:
            fc = gene_fc.get(g)
            if fc is None or not np.isfinite(fc):
                continue
            if fc > 0:
                pos_fcs.append(float(fc))
    if pos_fcs:
        abs_fcs = sorted((abs(v) for v in pos_fcs), reverse=True)
        skip = min(fc_cap_skip_top_n, len(abs_fcs) - 1)
        rank_cap = abs_fcs[skip]
        fc_zlim = float(np.clip(rank_cap, FC_RANGE_FLOOR, FC_RANGE_CEILING))
    else:
        fc_zlim = FC_RANGE_FLOOR
    fc_zmin, fc_zmax = 0.0, fc_zlim

    # --- Panels B-G: per-pathway intersection-gene log2FC heatmaps ----------
    panel_grid = [
        (2, 1), (2, 2), (2, 3),
        (3, 1), (3, 2), (3, 3),
    ]

    for panel_idx in range(n_pathway_panels):
        if panel_idx >= len(panel_grid):
            break
        row, col = panel_grid[panel_idx]
        term_row = pathway_rows.iloc[panel_idx]
        genes = parse_leading_edge(term_row["core_enrichment"])
        if not genes:
            continue

        gene_values = []
        gene_names = []
        for g in genes:
            fc = gene_fc.get(g)
            if fc is None or not np.isfinite(fc):
                continue   # presence/absence gene without a log2FC; omit
            gene_values.append([fc])
            gene_names.append(g)
            if len(gene_names) >= max_genes_per_panel:
                break

        if not gene_values:
            continue

        # Sort by |log2FC| ascending so the largest-magnitude effect renders at
        # the top of the cell stack (Plotly y is bottom-to-top).
        sorted_pairs = sorted(
            zip(gene_names, gene_values), key=lambda x: abs(x[1][0])
        )
        gene_names = [g for g, _ in sorted_pairs]
        gene_values = [v for _, v in sorted_pairs]

        # Protein-product labels: uppercase roman (style guide 4.3).
        gene_names = [g.upper() for g in gene_names]

        is_last_panel = (panel_idx == n_pathway_panels - 1)

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

        fig.add_trace(
            go.Heatmap(
                z=gene_values,
                y=gene_names,
                x=["log<sub>2</sub>FC"],
                colorscale=SINGLE_DIRECTION_UP,
                zmin=fc_zmin,
                zmax=fc_zmax,
                showscale=is_last_panel,
                colorbar=(colorbar_kwargs if is_last_panel else None),
                hovertemplate=(
                    "<b>%{y}</b><br>"
                    "log<sub>2</sub>FC: %{z:.2f}<extra></extra>"
                ),
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

    # Subplot-title annotation styling.
    for annotation in fig["layout"]["annotations"]:
        annotation["font"] = dict(
            family=FONT_FAMILY, size=SUBHEADING_SIZE, color=FONT_COLOR
        )

    # Panel letters (A-G) at the upper-left of each subplot domain.
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
# --- CLI --------------------------------------------------------------------
# =============================================================================
def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate the upregulated ORA GO:BP publication-style composite "
            "figure (3-row layout, top-N fold-enrichment Panel A + "
            "intersection-gene log2FC heatmaps), mirroring "
            "plot_gsea_publication_figure.py."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--ora-csv", type=Path, default=DEFAULT_ORA_CSV,
        help="ORA GO:BP _full results CSV.",
    )
    parser.add_argument(
        "--protein-csv", type=Path, default=DEFAULT_PROTEIN_CSV,
        help="Per-protein table providing log2FC for genes.",
    )
    parser.add_argument(
        "--revigo-csv", type=Path, default=DEFAULT_REVIGO_CSV,
        help="REVIGO representatives CSV (GO:BP).",
    )
    parser.add_argument(
        "--n-panels", type=int, default=6,
        help=("Number of pathway gene-heatmap panels (max 6; B-G arranged 3+3 "
              "across rows 2 and 3)."),
    )
    parser.add_argument(
        "--max-genes-per-panel", type=int, default=15,
        help="Maximum intersection genes shown per pathway heatmap.",
    )
    parser.add_argument(
        "--jaccard-dedup", type=float, default=0.7,
        help=("Maximum intersection-gene Jaccard between any two selected "
              "pathways. Lower = stricter dedup. Set to 1.0 to disable."),
    )
    parser.add_argument(
        "--fc-cap-skip-top-n", type=int,
        default=DEFAULT_FC_CAP_SKIP_TOP_N,
        help=("Number of top-extreme |log2FC| values to skip when setting the "
              "color-scale cap. Default 1 lets a single outlier saturate; 0 "
              "uses the absolute max (no saturation, the '_fullrange' variant)."),
    )
    parser.add_argument(
        "--no-revigo", action="store_true",
        help="Disable REVIGO simplification; use raw significant ORA terms.",
    )
    parser.add_argument(
        "--outdir", type=Path, default=DEFAULT_OUTDIR,
        help="Figure output directory.",
    )
    parser.add_argument(
        "--basename-suffix", default="",
        help=("Optional suffix appended to the output basename (before the "
              "extension), e.g. '--basename-suffix _fullrange'."),
    )
    parser.add_argument(
        "--png-scale", type=int, default=3,
        help="Plotly scale factor for PNG export (~300-400 DPI at this size).",
    )
    args = parser.parse_args(argv)
    if args.n_panels > 6:
        parser.error("--n-panels must be <= 6 (composite has 6 pathway slots).")
    return args


def main() -> None:
    args = parse_args()

    print(f"Reading ORA results: {args.ora_csv}")
    df = load_ora_full(args.ora_csv)
    print(f"  significant GO:BP terms: {len(df)}")

    print(f"Reading protein log2FC (quantitative only): {args.protein_csv}")
    gene_fc = load_quant_protein_fc(args.protein_csv)
    print(f"  quantitative gene -> log2FC entries: {len(gene_fc)}")

    revigo_ids = None
    if not args.no_revigo:
        print(f"REVIGO simplification enabled (GO:BP): {args.revigo_csv}")
        revigo_ids = load_revigo_representative_ids(args.revigo_csv)
        print(f"  REVIGO representatives: {len(revigo_ids)}")
    else:
        print("REVIGO simplification disabled by --no-revigo.")

    df_top = select_top(
        df, args.n_panels, revigo_ids, jaccard_threshold=args.jaccard_dedup
    )
    if df_top.empty:
        sys.exit("ERROR: no pathways selected; nothing to plot.")

    print(f"Selected {len(df_top)} terms (by FDR, most significant first; "
          f"intersection Jaccard <= {args.jaccard_dedup}):")
    for _, r in df_top.iterrows():
        print(f"  FDR={r['fdr_pvalue']:.2e}  FE={r['fold_enrichment']:5.2f}  "
              f"{r['term_name_display']}")

    fig = build_publication_figure(
        df_top, gene_fc,
        args.n_panels, args.max_genes_per_panel,
        fc_cap_skip_top_n=args.fc_cap_skip_top_n,
    )
    basename = f"ORA_GO_BP_publication_figure_up_thesis{args.basename_suffix}"
    export_figure(fig, args.outdir, basename, args.png_scale)


if __name__ == "__main__":
    main()

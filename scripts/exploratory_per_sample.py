# title: exploratory_per_sample.py
# project: Ketamine Proteomics Analysis Project
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-06-02
# last modified: 2026-06-02
# purpose:
#   Exploratory, NOT-for-publication look at the newly received per-sample
#   normalized abundance matrix (904 proteins x 8 runs). Assesses missingness,
#   runs PCA and a Welch t-test volcano on the complete-case 3v3 biological
#   design, and quantifies how the technical replicates (F5, F11) should be
#   handled by comparing "drop" vs "average-into-biological-partner".
# inputs:
#   results/exploratory/per_sample_normalized_abundance.csv
#     (tidy export of the 8 normalized-abundance columns + Accession/GeneSymbol)
# outputs:
#   results/exploratory/missingness_summary.txt
#   results/exploratory/pca_all8.png  pca_3v3.png
#   results/exploratory/volcano_3v3.png
#   results/exploratory/diffabundance_3v3_welch.csv
#   results/exploratory/techrep_handling_comparison.txt
# usage:
#   conda activate ketamine_project
#   python scripts/exploratory_per_sample.py
# NOTE: exploratory only. n=3 biological/group, proof-of-concept dataset.

from __future__ import annotations
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

# --- Config -----------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[1]
EXP = ROOT / "results" / "exploratory"
CSV = EXP / "per_sample_normalized_abundance.csv"

CTRL_BIO = ["F1_ctrl", "F4_ctrl", "F9_ctrl"]   # F5 = technical reinjection of F4
KET_BIO  = ["F2_ket",  "F6_ket",  "F10_ket"]   # F11 = technical reinjection of F10
ALL8 = ["F1_ctrl", "F4_ctrl", "F5_ctrl", "F9_ctrl",
        "F2_ket",  "F6_ket",  "F10_ket", "F11_ket"]
TECH_PAIRS = [("F4_ctrl", "F5_ctrl"), ("F10_ket", "F11_ket")]


def bh_fdr(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(p, float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def welch_volcano(df: pd.DataFrame, ctrl: list[str], ket: list[str],
                  min_per_group: int = 2) -> pd.DataFrame:
    """Per-protein Welch t-test on log2 abundance; log2FC = ket - ctrl."""
    rows = []
    for _, r in df.iterrows():
        c = np.log2(r[ctrl].dropna().astype(float).values)
        k = np.log2(r[ket].dropna().astype(float).values)
        if len(c) >= min_per_group and len(k) >= min_per_group:
            t, p = stats.ttest_ind(k, c, equal_var=False)
            rows.append((r["Accession"], r.get("GeneSymbol"),
                         len(c), len(k), k.mean() - c.mean(), p))
    res = pd.DataFrame(rows, columns=["Accession", "GeneSymbol",
                                      "n_ctrl", "n_ket", "log2FC", "p"])
    res["FDR"] = bh_fdr(res["p"].values)
    return res.sort_values("p").reset_index(drop=True)


def main() -> None:
    df = pd.read_csv(CSV)
    present = df[ALL8].notna()

    # --- 1. Missingness ----------------------------------------------------
    lines = ["MISSINGNESS SUMMARY (new per-sample matrix)", "=" * 50,
             f"Total proteins: {len(df)}", ""]
    lines.append("Per-run detection (# proteins with a value):")
    for s in ALL8:
        lines.append(f"  {s:10s} {int(present[s].sum()):4d}  "
                     f"({present[s].mean()*100:4.1f}%)")
    det8 = present.sum(axis=1)
    lines.append("\nDetection breadth (# of 8 runs a protein appears in):")
    for k in range(8, 0, -1):
        lines.append(f"  in {k}/8 runs: {(det8 == k).sum():4d}")
    cb = df[CTRL_BIO].notna().sum(axis=1)
    kb = df[KET_BIO].notna().sum(axis=1)
    lines += ["", "On the 3v3 biological design (F5,F11 dropped):",
              f"  complete (3 ctrl & 3 ket): {((cb == 3) & (kb == 3)).sum()}",
              f"  Welch-testable (>=2 & >=2): {((cb >= 2) & (kb >= 2)).sum()}",
              f"  detected both groups (>=1 & >=1): {((cb >= 1) & (kb >= 1)).sum()}",
              f"  control-bio only: {((cb >= 1) & (kb == 0)).sum()}",
              f"  ketamine-bio only: {((kb >= 1) & (cb == 0)).sum()}"]
    # MNAR gradient: abundance vs detection breadth
    lines.append("\nMNAR check (median abundance by detection breadth):")
    longvals = df[ALL8].values
    for k in range(1, 9):
        mask = det8 == k
        if mask.sum():
            med = np.nanmedian(longvals[mask.values])
            lines.append(f"  in {k}/8: n={mask.sum():4d}  median_abund={med:11.0f}")
    (EXP / "missingness_summary.txt").write_text("\n".join(lines))
    print("\n".join(lines))

    # --- 2. PCA (complete-case) -------------------------------------------
    def run_pca(cols: list[str], title: str, fname: str) -> None:
        sub = df[cols].dropna()
        X = np.log2(sub.values).T                      # samples x proteins
        X = X - X.mean(axis=0, keepdims=True)          # center per protein
        pcs = PCA(n_components=2).fit(X)
        Y = pcs.transform(X)
        ev = pcs.explained_variance_ratio_ * 100
        fig, ax = plt.subplots(figsize=(5, 4.2))
        for i, s in enumerate(cols):
            col = "#1f77b4" if "ctrl" in s else "#d62728"
            ax.scatter(Y[i, 0], Y[i, 1], c=col, s=90,
                       edgecolor="k", linewidth=0.5, zorder=3)
            ax.annotate(s.split("_")[0], (Y[i, 0], Y[i, 1]),
                        xytext=(5, 4), textcoords="offset points", fontsize=8)
        ax.set_xlabel(f"PC1 ({ev[0]:.1f}%)")
        ax.set_ylabel(f"PC2 ({ev[1]:.1f}%)")
        ax.set_title(f"{title}\n(complete-case n={sub.shape[0]} proteins)",
                     fontsize=9)
        ax.scatter([], [], c="#1f77b4", label="control")
        ax.scatter([], [], c="#d62728", label="ketamine")
        ax.legend(fontsize=8, frameon=False)
        fig.tight_layout()
        fig.savefig(EXP / fname, dpi=150)
        plt.close(fig)
        print(f"  PCA[{title}] PC1={ev[0]:.1f}% PC2={ev[1]:.1f}% "
              f"on {sub.shape[0]} proteins")

    print("\nPCA:")
    run_pca(ALL8, "All 8 runs (incl. technical reps)", "pca_all8.png")
    run_pca(CTRL_BIO + KET_BIO, "3v3 biological (F5,F11 dropped)", "pca_3v3.png")

    # --- 3. Volcano on the drop-F5/F11 3v3 design -------------------------
    res = welch_volcano(df, CTRL_BIO, KET_BIO, min_per_group=2)
    res.to_csv(EXP / "diffabundance_3v3_welch.csv", index=False)
    sig_p = (res["p"] < 0.05).sum()
    sig_fdr = (res["FDR"] < 0.05).sum()
    print(f"\nVolcano (3v3, drop F5/F11): {len(res)} tested, "
          f"{sig_p} raw p<0.05, {sig_fdr} FDR<0.05")

    fig, ax = plt.subplots(figsize=(5.5, 4.6))
    x = res["log2FC"].values
    y = -np.log10(res["p"].values)
    sig = (res["p"] < 0.05) & (res["log2FC"].abs() > 1)
    ax.scatter(x[~sig], y[~sig], s=10, c="#bbbbbb", linewidth=0)
    ax.scatter(x[sig], y[sig], s=14, c="#d62728", linewidth=0)
    ax.axhline(-np.log10(0.05), ls="--", c="k", lw=0.6)
    ax.axvline(1, ls="--", c="k", lw=0.6)
    ax.axvline(-1, ls="--", c="k", lw=0.6)
    for _, r in res[sig].head(12).iterrows():
        lbl = r["GeneSymbol"] if isinstance(r["GeneSymbol"], str) else r["Accession"]
        ax.annotate(lbl, (r["log2FC"], -np.log10(r["p"])),
                    fontsize=6, xytext=(2, 2), textcoords="offset points")
    ax.set_xlabel("log2 fold change (ketamine / control)")
    ax.set_ylabel("-log10 p (Welch)")
    ax.set_title("Exploratory volcano, 3v3 biological\n"
                 "n=3/group proof-of-concept, raw p", fontsize=9)
    fig.tight_layout()
    fig.savefig(EXP / "volcano_3v3.png", dpi=150)
    plt.close(fig)

    # --- 4. Technical-replicate handling: drop vs average -----------------
    avg = df.copy()
    avg["F4_ctrl"] = df[["F4_ctrl", "F5_ctrl"]].mean(axis=1, skipna=True)
    avg["F10_ket"] = df[["F10_ket", "F11_ket"]].mean(axis=1, skipna=True)
    res_avg = welch_volcano(avg, CTRL_BIO, KET_BIO, min_per_group=2)

    m = res.merge(res_avg, on="Accession", suffixes=("_drop", "_avg"))
    rfc = np.corrcoef(m["log2FC_drop"], m["log2FC_avg"])[0, 1]
    # technical-rep CV
    cvs = []
    for a, b in TECH_PAIRS:
        pair = df[[a, b]].dropna()
        cvs += list(pair.std(axis=1, ddof=0) / pair.mean(axis=1) * 100)
    cvs = np.array(cvs)
    tl = ["TECHNICAL REPLICATE HANDLING", "=" * 50,
          "Reinjection pairs: F4/F5 (control), F10/F11 (ketamine)",
          f"Technical CV: median {np.median(cvs):.1f}%, mean {cvs.mean():.1f}%, "
          f"{(cvs < 20).mean()*100:.0f}% under 20%", "",
          "Effect on the volcano:",
          f"  DROP F5/F11   : {len(res)} tested, "
          f"{(res['p']<0.05).sum()} raw p<0.05, {(res['FDR']<0.05).sum()} FDR<0.05",
          f"  AVERAGE F4/F5 & F10/F11: {len(res_avg)} tested, "
          f"{(res_avg['p']<0.05).sum()} raw p<0.05, {(res_avg['FDR']<0.05).sum()} FDR<0.05",
          f"  log2FC correlation drop-vs-average: r = {rfc:.4f}",
          f"  proteins p<0.05 in both: "
          f"{((m['p_drop']<0.05) & (m['p_avg']<0.05)).sum()}"]
    (EXP / "techrep_handling_comparison.txt").write_text("\n".join(tl))
    print("\n" + "\n".join(tl))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""plot_burden.py — Stage 14 burden QC plots (CallForge).

  burden_qq.{png,svg}        QQ plot of gene-burden p-values (+ genomic-inflation lambda)
  burden_manhattan.{png,svg} gene-level Manhattan (-log10 p per unit)

Reads burden_results.tsv + burden_summary.json. PNG+SVG, captions to $CALLFORGE_CAPTIONS.
Caption states the small-N validity caveat plainly.
"""
import argparse, csv, json, math, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def save(fig, outdir, stem, caption):
    os.makedirs(outdir, exist_ok=True)
    fig.text(0.5, 0.005, caption, ha="center", va="bottom", fontsize=8, style="italic", wrap=True)
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(outdir, f"{stem}.{ext}"), bbox_inches="tight", dpi=140)
    plt.close(fig)
    with open(os.path.join(outdir, os.environ.get("CALLFORGE_CAPTIONS", "captions.tsv")), "a") as fh:
        fh.write(f"{stem}\t{caption}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    summ = json.load(open(a.summary)) if os.path.isfile(a.summary) else {}
    lam = summ.get("genomic_inflation_lambda")

    # gene-level results with a numeric p
    pts = []
    if os.path.isfile(a.results):
        for r in csv.DictReader(open(a.results), delimiter="\t"):
            if r.get("kind") == "gene" and r.get("p") not in ("NA", "", None):
                try:
                    pts.append((r["unit"].split(":", 1)[1], float(r["p"])))
                except ValueError:
                    pass

    # ---- QQ ----
    fig, ax = plt.subplots(figsize=(5, 5))
    if pts:
        ps = sorted(p for _, p in pts)
        n = len(ps)
        obs = [-math.log10(max(p, 1e-12)) for p in ps]
        exp = [-math.log10((i + 0.5) / n) for i in range(n)]
        mx = max(obs + exp + [1])
        ax.plot([0, mx], [0, mx], color="grey", ls="--", lw=1)
        ax.scatter(exp, obs, c="#1f77b4")
        ax.set_xlabel("expected -log10(p)"); ax.set_ylabel("observed -log10(p)")
        if lam is not None:
            ax.text(0.05, 0.92, f"λ = {lam}", transform=ax.transAxes, fontsize=10)
    else:
        ax.text(0.5, 0.5, "no testable gene units\n(no qualifying rare+functional variants)",
                ha="center", va="center", transform=ax.transAxes)
    ax.set_title("Gene-burden QQ plot")
    save(fig, a.outdir, "burden_qq",
         "QQ plot of gene-burden p-values vs uniform expectation; lambda = genomic inflation. "
         "Small-N synthetic results prove the machinery only — not valid association.")

    # ---- Manhattan ----
    fig, ax = plt.subplots(figsize=(max(5, 0.6 * max(1, len(pts)) + 2), 4))
    if pts:
        pts.sort(key=lambda x: x[0])
        genes = [g for g, _ in pts]
        ax.scatter(range(len(genes)), [-math.log10(max(p, 1e-12)) for _, p in pts], c="#2ca02c")
        ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes, rotation=45, ha="right")
        ax.axhline(-math.log10(0.05 / max(1, len(genes))), color="#d62728", ls="--", lw=1,
                   label="Bonferroni 0.05")
        ax.set_ylabel("-log10(p)"); ax.legend(fontsize=7)
    else:
        ax.text(0.5, 0.5, "no testable gene units", ha="center", va="center", transform=ax.transAxes)
        ax.set_xticks([])
    ax.set_title("Gene-burden Manhattan")
    save(fig, a.outdir, "burden_manhattan",
         "Per-gene burden significance (-log10 p) with a Bonferroni line; collapses rare+functional "
         "variants per gene. Validity depends on N, phenotype definition, ancestry matching.")
    print(f"[plot_burden] {len(pts)} gene units plotted; lambda={lam}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""plot_coverage.py — Stage 4 capture/coverage QC plots (CallForge).

  on_target_enrichment.{png,svg}  % on-target + fold-enrichment per sample
  dup_rate.{png,svg}              duplication rate per sample
  coverage_uniformity.{png,svg}   FOLD_80_BASE_PENALTY per sample (lower = more uniform)
  coverage_cumulative.{png,svg}   % target bases >= depth (mosdepth region dist)
  per_gene_depth.{png,svg}        per-gene mean depth, genes x samples (CaptureForge track style)

PNG+SVG, captions embedded + appended to <outdir>/captions.tsv.
"""
import argparse, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def save(fig, outdir, stem, caption):
    os.makedirs(outdir, exist_ok=True)
    fig.text(0.5, 0.005, caption, ha="center", va="bottom", fontsize=8, style="italic", wrap=True)
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(outdir, f"{stem}.{ext}"), bbox_inches="tight", dpi=140)
    plt.close(fig)
    with open(os.path.join(outdir, os.environ.get("CALLFORGE_CAPTIONS", "captions.tsv")), "a") as fh:
        fh.write(f"{stem}\t{caption}\n")


def bars(qs, sids, key, ylabel, title, stem, caption, outdir, pct=False, hline=None):
    vals = [(q.get(key) or 0) * (100 if pct else 1) for q in qs]
    fig, ax = plt.subplots(figsize=(max(5, 0.6 * len(sids) + 2), 4))
    ax.bar(sids, vals, color="#1f77b4")
    if hline is not None:
        ax.axhline(hline, color="#d62728", ls="--", lw=1, label=f"threshold {hline}")
        ax.legend(fontsize=7)
    ax.set_ylabel(ylabel); ax.set_title(title); ax.tick_params(axis="x", rotation=45)
    save(fig, outdir, stem, caption)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qc-json", nargs="+", required=True)
    ap.add_argument("--region-dist", nargs="*", default=[], help="mosdepth *.region.dist.txt")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    qs = sorted((json.load(open(p)) for p in a.qc_json), key=lambda x: x["sample_id"])
    sids = [q["sample_id"] for q in qs]

    bars(qs, sids, "pct_on_target", "% on-target", "On-target rate per sample",
         "on_target_enrichment_ontarget", pct=True, hline=40,
         caption="Percent of bases on the capture target per sample (capture specificity).", outdir=a.outdir)
    bars(qs, sids, "fold_enrichment", "fold enrichment", "Capture fold-enrichment per sample",
         "on_target_enrichment_fold",
         caption="Fold-enrichment of on-target vs genome-average coverage (capture performance).", outdir=a.outdir)
    bars(qs, sids, "dup_rate", "duplication rate", "Duplication rate per sample",
         "dup_rate", hline=0.40,
         caption="Library duplication rate per sample (coordinate dedup; lower is better).", outdir=a.outdir)
    bars(qs, sids, "fold_80_base_penalty", "FOLD_80 penalty", "Coverage uniformity (FOLD_80)",
         "coverage_uniformity",
         caption="FOLD_80_BASE_PENALTY: fold extra sequencing for 80% of targets to reach the mean (1.0 = perfectly uniform).", outdir=a.outdir)

    # cumulative coverage from mosdepth region dist
    if a.region_dist:
        fig, ax = plt.subplots(figsize=(6, 4))
        plotted = False
        for p in a.region_dist:
            sid = os.path.basename(p).split(".")[0]
            xs, ys = [], []
            for ln in open(p):
                f = ln.rstrip("\n").split("\t")
                if len(f) == 3 and f[0] == "total":
                    xs.append(int(f[1])); ys.append(float(f[2]) * 100)
            if xs:
                order = np.argsort(xs)
                ax.plot(np.array(xs)[order], np.array(ys)[order], label=sid, lw=1.2); plotted = True
        if plotted:
            ax.set_xlabel("depth (x)"); ax.set_ylabel("% target bases >= depth")
            ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.legend(fontsize=7)
            ax.set_title("Cumulative target coverage")
            save(fig, a.outdir, "coverage_cumulative",
                 "Fraction of target bases covered at >= each depth, per sample (sensitivity for calling).")
        else:
            plt.close(fig)

    # per-gene depth heatmap (genes x samples) — CaptureForge track style
    genes = sorted({g for q in qs for g in (q.get("per_gene_depth") or {})})
    if genes:
        M = np.array([[ (q.get("per_gene_depth") or {}).get(g, np.nan) for q in qs] for g in genes])
        fig, ax = plt.subplots(figsize=(max(4, 0.6 * len(sids) + 2), max(3, 0.3 * len(genes) + 2)))
        im = ax.imshow(M, aspect="auto", cmap="viridis")
        ax.set_xticks(range(len(sids))); ax.set_xticklabels(sids, rotation=45, ha="right")
        ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes, fontsize=7)
        fig.colorbar(im, ax=ax, label="mean depth (x)")
        ax.set_title("Per-gene mean depth")
        save(fig, a.outdir, "per_gene_depth",
             "Mean coverage depth per target gene per sample; spot under-covered genes (e.g. CR1/CFH).")
    print(f"[plot_coverage] wrote coverage QC plots to {a.outdir}")


if __name__ == "__main__":
    main()

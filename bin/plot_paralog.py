#!/usr/bin/env python3
"""plot_paralog.py — Stage 10 paralog-aware QC plot (CallForge).

  paralog_confidence.{png,svg}  variants in each paralog gene region, split
                                confident vs flagged_low_confidence.

PNG+SVG, caption embedded + appended to $CALLFORGE_CAPTIONS.
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    s = json.load(open(a.summary))
    per_gene = s.get("per_gene", {})
    genes = sorted(per_gene) or s.get("paralog_genes", [])

    fig, ax = plt.subplots(figsize=(max(5, 0.7 * max(1, len(genes)) + 2), 4))
    if per_gene:
        conf = [per_gene[g].get("confident", 0) for g in genes]
        flag = [per_gene[g].get("flagged_low_confidence", 0) for g in genes]
        x = np.arange(len(genes))
        ax.bar(x, conf, label="confident (uniquely callable)", color="#2ca02c")
        ax.bar(x, flag, bottom=conf, label="flagged low-confidence", color="#d62728")
        ax.set_xticks(x); ax.set_xticklabels(genes, rotation=45, ha="right")
        ax.set_ylabel("variants"); ax.legend(fontsize=7)
    else:
        ax.text(0.5, 0.5, f"no variants in paralog regions\n(paralog genes: {', '.join(s.get('paralog_genes', []))})",
                ha="center", va="center", transform=ax.transAxes); ax.set_xticks([])
    ax.set_title("Paralog-region variants: confident vs flagged")
    save(fig, a.outdir, "paralog_confidence",
         "Variants in each paralog-ambiguous gene region, split into confident (uniquely "
         "callable, high MQ) vs flagged low-confidence (multimapping-ambiguous, low MQ).")
    print(f"[plot_paralog] wrote paralog plot to {a.outdir}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""plot_cnv.py — Stage 8 CNV QC plots (CallForge).

  cnv_copyratio.{png,svg}      per-target log2 copy-ratio, coloured by sample
  cnv_calls_per_sample.{png,svg} number of non-neutral CNV calls per sample
  cnv_callability.{png,svg}    per-gene CNV-callability view: callable vs
                               breakpoint-blind (CaptureForge labels) so
                               low-confidence loci (CR1/CFH) are visibly marked

PNG+SVG, captions embedded + appended to $CALLFORGE_CAPTIONS.
"""
import argparse, csv, json, os
from collections import defaultdict
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
    ap.add_argument("--copyratio", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    # per-target copy ratio
    by_sample = defaultdict(list)
    with open(a.copyratio) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                by_sample[r["sample"]].append((r["gene"], float(r["log2"])))
            except ValueError:
                pass
    fig, ax = plt.subplots(figsize=(7, 4))
    for sample, pts in sorted(by_sample.items()):
        ax.scatter(range(len(pts)), [p[1] for p in pts], s=8, alpha=0.6, label=sample)
    ax.axhline(0, color="grey", lw=0.8)
    ax.set_xlabel("target bin (ordered)"); ax.set_ylabel("log2 copy ratio")
    ax.set_title("Per-target copy ratio"); ax.legend(fontsize=7)
    save(fig, a.outdir, "cnv_copyratio",
         "Per-target log2 copy ratio per sample; deviations from 0 indicate copy-number change.")

    # calls per sample
    summary = json.load(open(a.summary))
    cps = summary.get("calls_per_sample", {})
    fig, ax = plt.subplots(figsize=(max(5, 0.6 * max(1, len(cps)) + 2), 4))
    if cps:
        ax.bar(list(cps), list(cps.values()), color="#1f77b4")
    else:
        ax.text(0.5, 0.5, "no non-neutral CNV calls", ha="center", va="center", transform=ax.transAxes)
        ax.set_xticks([])
    ax.set_ylabel("CNV calls"); ax.set_title("CNV calls per sample"); ax.tick_params(axis="x", rotation=45)
    save(fig, a.outdir, "cnv_calls_per_sample",
         "Number of non-neutral copy-number calls per sample.")

    # per-gene callability view
    genes, flags = [], []
    with open(a.metadata) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("is_cnv_target") == "yes":
                genes.append(r["gene"]); flags.append(r.get("cnv_callable", "unknown"))
    fig, ax = plt.subplots(figsize=(max(4, 0.7 * max(1, len(genes)) + 2), 3.5))
    colmap = {"yes": "#2ca02c", "no": "#d62728", "unknown": "#999999"}
    if genes:
        ax.bar(genes, [1] * len(genes), color=[colmap.get(f, "#999999") for f in flags])
        for i, f in enumerate(flags):
            lbl = {"yes": "callable", "no": "breakpoint-blind", "unknown": "unknown"}.get(f, f)
            ax.text(i, 0.5, lbl, ha="center", va="center", rotation=90, color="white", fontweight="bold", fontsize=8)
        ax.set_yticks([])
    else:
        ax.text(0.5, 0.5, "no CNV-target genes", ha="center", va="center", transform=ax.transAxes); ax.set_xticks([])
    ax.set_title("Per-gene CNV callability (CaptureForge)"); ax.tick_params(axis="x", rotation=45)
    save(fig, a.outdir, "cnv_callability",
         "CNV-target genes coloured by CaptureForge callability: green=depth-callable, "
         "red=breakpoint-blind (low-confidence, e.g. CR1/CFH).")
    print(f"[plot_cnv] wrote CNV plots to {a.outdir}")


if __name__ == "__main__":
    main()

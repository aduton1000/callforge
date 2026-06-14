#!/usr/bin/env python3
"""plot_filter_qc.py — Stage 7 hard-filter QC plots (CallForge).

  filter_before_after.{png,svg}  total vs PASS counts for SNP / indel
  filter_pass_rate.{png,svg}     PASS rate per type (+ overall)
  filter_reasons.{png,svg}       breakdown of hard-filter rejection reasons

PNG+SVG, captions embedded + appended to $CALLFORGE_CAPTIONS (default captions.tsv).
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
    ap.add_argument("--counts-json", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    c = json.load(open(a.counts_json))

    # before vs after
    types = ["snp", "indel"]
    totals = [c.get(f"{t}_total", 0) for t in types]
    passes = [c.get(f"{t}_pass", 0) for t in types]
    x = np.arange(len(types))
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.bar(x - 0.2, totals, width=0.4, label="called (pre-filter)", color="#999999")
    ax.bar(x + 0.2, passes, width=0.4, label="PASS (post-filter)", color="#2ca02c")
    ax.set_xticks(x); ax.set_xticklabels(["SNP", "indel"]); ax.set_ylabel("sites"); ax.legend(fontsize=8)
    ax.set_title("Hard-filter: called vs PASS")
    save(fig, a.outdir, "filter_before_after",
         "Variant counts before and after GATK hard-filtering, by type (sites retained).")

    # pass rate
    rates = [(c.get(f"{t}_pass_rate") or 0) * 100 for t in types] + [(c.get("pass_rate_overall") or 0) * 100]
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.bar(["SNP", "indel", "overall"], rates, color=["#1f77b4", "#1f77b4", "#ff7f0e"])
    ax.set_ylabel("% PASS"); ax.set_ylim(0, 100); ax.set_title("Hard-filter PASS rate")
    save(fig, a.outdir, "filter_pass_rate",
         "Fraction of called variants passing the GATK hard-filters, by type and overall.")

    # reasons
    reasons = c.get("filter_reasons", {})
    fig, ax = plt.subplots(figsize=(max(5, 0.6 * len(reasons) + 2), 4))
    if reasons:
        names = sorted(reasons, key=lambda k: -reasons[k])
        ax.bar(names, [reasons[n] for n in names], color="#d62728")
        ax.set_ylabel("sites rejected"); ax.tick_params(axis="x", rotation=45)
    else:
        ax.text(0.5, 0.5, "no variants rejected by hard-filters", ha="center", va="center",
                transform=ax.transAxes); ax.set_xticks([])
    ax.set_title("Hard-filter rejection reasons")
    save(fig, a.outdir, "filter_reasons",
         "Number of sites rejected by each hard-filter criterion (which filter removed what).")
    print(f"[plot_filter_qc] wrote filter QC plots to {a.outdir}")


if __name__ == "__main__":
    main()

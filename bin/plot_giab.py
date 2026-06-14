#!/usr/bin/env python3
"""plot_giab.py — Stage 13 GIAB benchmarking plot (CallForge).

Parses hap.py summary.csv (SNV + indel, PASS) into precision / recall / F1 and
renders:
  giab_precision_recall.{png,svg}  precision & recall per variant type
  giab_f1.{png,svg}                F1 per variant type

The hap.py comparison is RESTRICTED TO THE PANEL BED (on-target only) — scoring
genome-wide would distort a panel's measured sensitivity. On synthetic test data
this only proves the benchmarking machinery runs; real numbers need a real GIAB
control + truth v4.2.1.
"""
import argparse, csv, json, os
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


def num(x):
    try:
        return float(x)
    except (ValueError, TypeError):
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", required=True, help="hap.py summary.csv")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--restricted-to", default="panel BED", help="region restriction (for the caption)")
    a = ap.parse_args()

    metrics = {}
    with open(a.summary) as fh:
        for r in csv.DictReader(fh):
            if r.get("Filter") != "PASS":
                continue
            t = r.get("Type")  # SNP / INDEL
            metrics[t] = {"recall": num(r.get("METRIC.Recall")),
                          "precision": num(r.get("METRIC.Precision")),
                          "f1": num(r.get("METRIC.F1_Score")),
                          "truth_total": num(r.get("TRUTH.TOTAL")),
                          "tp": num(r.get("TRUTH.TP")), "fn": num(r.get("TRUTH.FN")),
                          "fp": num(r.get("QUERY.FP"))}
    with open(os.path.join(a.outdir, "giab_metrics.json"), "w") as fh:
        json.dump({"restricted_to": a.restricted_to, "metrics": metrics}, fh, indent=2)

    types = [t for t in ("SNP", "INDEL") if t in metrics]
    # precision / recall grouped bars
    fig, ax = plt.subplots(figsize=(5.5, 4))
    if types:
        import numpy as np
        x = np.arange(len(types))
        ax.bar(x - 0.2, [(metrics[t]["precision"] or 0) for t in types], 0.4, label="precision", color="#1f77b4")
        ax.bar(x + 0.2, [(metrics[t]["recall"] or 0) for t in types], 0.4, label="recall", color="#ff7f0e")
        ax.set_xticks(x); ax.set_xticklabels(types); ax.set_ylim(0, 1.05); ax.legend(fontsize=8)
    else:
        ax.text(0.5, 0.5, "no hap.py metrics", ha="center", va="center", transform=ax.transAxes)
    ax.set_title("GIAB precision / recall (on-target)")
    save(fig, a.outdir, "giab_precision_recall",
         f"hap.py precision & recall vs GIAB truth, restricted to the {a.restricted_to} (on-target). "
         "Synthetic test only proves the machinery; real numbers need a GIAB control + truth v4.2.1.")

    # F1
    fig, ax = plt.subplots(figsize=(5, 4))
    if types:
        ax.bar(types, [(metrics[t]["f1"] or 0) for t in types], color="#2ca02c"); ax.set_ylim(0, 1.05)
        for i, t in enumerate(types):
            ax.text(i, (metrics[t]["f1"] or 0) + 0.02, f"{metrics[t]['f1']:.3f}", ha="center", fontsize=8)
    else:
        ax.text(0.5, 0.5, "no hap.py metrics", ha="center", va="center", transform=ax.transAxes)
    ax.set_ylabel("F1"); ax.set_title("GIAB F1 (on-target)")
    save(fig, a.outdir, "giab_f1",
         f"hap.py F1 per variant type vs GIAB truth, restricted to the {a.restricted_to}.")
    print(f"[plot_giab] metrics: {metrics}")


if __name__ == "__main__":
    main()

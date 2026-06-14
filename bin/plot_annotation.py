#!/usr/bin/env python3
"""plot_annotation.py — Stage 11 annotation QC plots (CallForge).

  annotation_landing.{png,svg}   % of variants annotated per source (the landing proof)
  consequence_dist.{png,svg}     VEP consequence distribution
  novel_vs_known.{png,svg}       known (dbSNP rsID) vs novel
  gnomad_afr_af.{png,svg}        distribution of gnomAD-AFR AF (+ % with any AFR AF)
  phylop_dist.{png,svg}          PhyloP conservation-score distribution

Reads the flattened per-variant TSV (flatten_vcf.py) + annotation_landing.json.
PNG+SVG, captions to $CALLFORGE_CAPTIONS.
"""
import argparse, csv, json, math, os
from collections import Counter
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
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--landing", required=True)
    ap.add_argument("--af-field", default="AF_afr")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    afr = f"gnomAD_{a.af_field}"

    rows = list(csv.DictReader(open(a.tsv), delimiter="\t"))
    n = len(rows)

    # ---- landing (% annotated per source) ----
    landing = json.load(open(a.landing)).get("per_source", {})
    labels, pcts = [], []
    for k, v in landing.items():
        if v.get("status") == "ok":
            labels.append(k); pcts.append(v.get("pct_annotated", 0))
    fig, ax = plt.subplots(figsize=(max(5, 0.8 * max(1, len(labels)) + 2), 4))
    if labels:
        ax.bar(labels, pcts, color="#2ca02c"); ax.set_ylim(0, 100); ax.set_ylabel("% variants annotated")
        for i, p in enumerate(pcts):
            ax.text(i, p + 1, f"{p}%", ha="center", fontsize=8)
    else:
        ax.text(0.5, 0.5, "no sources applied", ha="center", va="center", transform=ax.transAxes)
    ax.set_title("Annotation landing: % variants annotated per source")
    save(fig, a.outdir, "annotation_landing",
         "Fraction of callset variants that actually received each annotation (proof the "
         "annotation landed, not silently missed on a contig-name mismatch).")

    # ---- consequence distribution ----
    cons = Counter(r["consequence"] for r in rows if r.get("consequence") not in (".", ""))
    fig, ax = plt.subplots(figsize=(7, 4))
    if cons:
        items = cons.most_common(15)
        ax.barh([k for k, _ in items][::-1], [v for _, v in items][::-1], color="#1f77b4")
        ax.set_xlabel("variants")
    else:
        ax.text(0.5, 0.5, "no VEP consequences (VEP not run?)", ha="center", va="center", transform=ax.transAxes)
    ax.set_title("VEP consequence distribution")
    save(fig, a.outdir, "consequence_dist",
         "Most-severe VEP consequence per variant across the callset.")

    # ---- novel vs known (dbSNP) ----
    known = sum(1 for r in rows if r.get("rsID", ".").startswith("rs"))
    novel = n - known
    fig, ax = plt.subplots(figsize=(4.5, 4))
    ax.bar(["known (dbSNP)", "novel"], [known, novel], color=["#1f77b4", "#ff7f0e"])
    ax.set_ylabel("variants"); ax.set_title("Novel vs known (dbSNP rsID)")
    save(fig, a.outdir, "novel_vs_known",
         "Variants with a dbSNP rsID (known) vs without (novel) — novelty rate.")

    # ---- gnomAD-AFR AF distribution ----
    afs = []
    for r in rows:
        v = r.get(afr, ".")
        try:
            f = float(v)
            if f > 0:
                afs.append(f)
        except (ValueError, TypeError):
            pass
    with_afr = sum(1 for r in rows if r.get(afr, ".") not in (".", "", None))
    fig, ax = plt.subplots(figsize=(6, 4))
    if afs:
        ax.hist([math.log10(x) for x in afs], bins=30, color="#9467bd")
        ax.set_xlabel("log10 gnomAD-AFR AF"); ax.set_ylabel("variants")
    else:
        ax.text(0.5, 0.5, "no gnomAD-AFR AF values", ha="center", va="center", transform=ax.transAxes)
    ax.set_title(f"gnomAD-AFR AF distribution ({round(100*with_afr/max(1,n),1)}% of variants have an AFR AF)")
    save(fig, a.outdir, "gnomad_afr_af",
         "Distribution of gnomAD African-ancestry allele frequencies; title shows the % of "
         "variants carrying an AFR AF (key for rare-variant filtering / burden).")

    # ---- PhyloP distribution ----
    phylop = []
    for r in rows:
        try:
            phylop.append(float(r.get("PhyloP", ".")))
        except (ValueError, TypeError):
            pass
    fig, ax = plt.subplots(figsize=(6, 4))
    if phylop:
        ax.hist(phylop, bins=30, color="#17becf")
        ax.set_xlabel("PhyloP score"); ax.set_ylabel("variants")
    else:
        ax.text(0.5, 0.5, "no PhyloP scores (bigWig absent)", ha="center", va="center", transform=ax.transAxes)
    ax.set_title("PhyloP conservation distribution")
    save(fig, a.outdir, "phylop_dist",
         "PhyloP conservation-score distribution; positive = conserved (likely functional).")
    print(f"[plot_annotation] wrote annotation plots to {a.outdir}")


if __name__ == "__main__":
    main()

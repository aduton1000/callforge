#!/usr/bin/env python3
"""plot_calling_qc.py — Stage 6 SNV/indel calling QC plots (CallForge).

  variants_per_sample.{png,svg}  variant count per sample (het + hom-alt)
  titv_per_sample.{png,svg}      transition/transversion ratio per sample (+ overall)
  het_hom_per_sample.{png,svg}   heterozygous/hom-alt ratio per sample
  qual_dist.{png,svg}            site QUAL distribution (PASS vs filtered)
  dp_dist.{png,svg}              site INFO/DP distribution

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
    ap.add_argument("--stats-json", required=True, help="parsed bcftools stats JSON (joint, pre-filter)")
    ap.add_argument("--qual-dp-tsv", help="QUAL<TAB>DP<TAB>TYPE<TAB>FILTER per site")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    st = json.load(open(a.stats_json))
    ps = st.get("per_sample", {})
    sids = sorted(ps)

    def bar(values, ylabel, title, stem, caption, hline=None, hlabel=None):
        fig, ax = plt.subplots(figsize=(max(5, 0.6 * len(sids) + 2), 4))
        ax.bar(sids, values, color="#1f77b4")
        if hline is not None:
            ax.axhline(hline, color="#d62728", ls="--", lw=1, label=hlabel)
            ax.legend(fontsize=7)
        ax.set_ylabel(ylabel); ax.set_title(title); ax.tick_params(axis="x", rotation=45)
        save(fig, a.outdir, stem, caption)

    bar([ps[s]["n_variants"] for s in sids], "variants (het + hom-alt)",
        "Variants per sample", "variants_per_sample",
        "Per-sample count of non-reference genotypes in the joint callset (calling yield).")
    bar([ps[s]["ts_tv"] or 0 for s in sids], "Ti/Tv", "Ti/Tv per sample",
        "titv_per_sample",
        "Transition/transversion ratio per sample; whole-exome ~3.0, smaller panels vary (calling quality).",
        hline=st.get("ts_tv_overall"), hlabel=f"overall {st.get('ts_tv_overall')}")
    bar([ps[s]["het_hom_ratio"] or 0 for s in sids], "het / hom-alt",
        "Het:Hom-alt ratio per sample", "het_hom_per_sample",
        "Heterozygous-to-homozygous-alt ratio per sample; extreme values flag sample-quality issues.")

    if a.qual_dp_tsv and os.path.isfile(a.qual_dp_tsv):
        quals_pass, quals_fail, dps = [], [], []
        for ln in open(a.qual_dp_tsv):
            f = ln.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            q, d, typ, filt = f[0], f[1], f[2], f[3]
            try:
                qv = float(q)
                (quals_pass if filt in ("PASS", ".") else quals_fail).append(qv)
            except ValueError:
                pass
            try:
                dps.append(float(d))
            except ValueError:
                pass
        fig, ax = plt.subplots(figsize=(6, 4))
        if quals_pass: ax.hist(quals_pass, bins=30, alpha=0.7, label="PASS", color="#2ca02c")
        if quals_fail: ax.hist(quals_fail, bins=30, alpha=0.7, label="filtered", color="#d62728")
        ax.set_xlabel("QUAL"); ax.set_ylabel("sites"); ax.legend(fontsize=7)
        ax.set_title("Site QUAL distribution")
        save(fig, a.outdir, "qual_dist",
             "Variant-site QUAL distribution, PASS vs hard-filtered (call confidence).")
        if dps:
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.hist(dps, bins=30, color="#1f77b4")
            ax.set_xlabel("INFO/DP (site depth)"); ax.set_ylabel("sites")
            ax.set_title("Site depth (DP) distribution")
            save(fig, a.outdir, "dp_dist",
                 "Per-site total depth (INFO/DP) distribution across the joint callset.")
    print(f"[plot_calling_qc] wrote calling QC plots to {a.outdir}")


if __name__ == "__main__":
    main()

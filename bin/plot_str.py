#!/usr/bin/env python3
"""plot_str.py — Stage 9 STR QC plots (CallForge).

  str_allele_sizes.{png,svg}   repeat-size (allele) distribution per locus
  str_call_rate.{png,svg}      per-locus call rate across the cohort
  str_read_support.{png,svg}   spanning-read support per sample x locus

PNG+SVG, captions embedded + appended to $CALLFORGE_CAPTIONS.
"""
import argparse, csv, json, os
from collections import defaultdict
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
    ap.add_argument("--calls", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    alleles = defaultdict(list)       # locus -> [repeat units]
    support = defaultdict(dict)       # (sample) -> {locus: reads}
    with open(a.calls) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            for k in ("allele1", "allele2"):
                if r.get(k) not in ("", None):
                    try:
                        alleles[r["locus"]].append(int(r[k]))
                    except ValueError:
                        pass
            try:
                support[r["sample"]][r["locus"]] = int(r["spanning_reads"])
            except (ValueError, TypeError):
                pass

    # allele-size distribution per locus
    loci = sorted(alleles)
    fig, ax = plt.subplots(figsize=(6, 4))
    if loci:
        for loc in loci:
            ax.hist(alleles[loc], bins=range(0, max(alleles[loc]) + 3), alpha=0.6, label=loc)
        ax.set_xlabel("repeat units (allele size)"); ax.set_ylabel("alleles"); ax.legend(fontsize=7)
    else:
        ax.text(0.5, 0.5, "no STR alleles genotyped", ha="center", va="center", transform=ax.transAxes)
    ax.set_title("STR allele/repeat-size distribution")
    save(fig, a.outdir, "str_allele_sizes",
         "Distribution of genotyped repeat-unit counts (allele sizes) per STR locus.")

    # call rate per locus
    summary = json.load(open(a.summary))
    cr = summary.get("call_rate_per_locus", {})
    fig, ax = plt.subplots(figsize=(max(5, 0.7 * max(1, len(cr)) + 2), 4))
    if cr:
        ax.bar(list(cr), [v * 100 for v in cr.values()], color="#1f77b4"); ax.set_ylim(0, 100)
        ax.set_ylabel("% samples genotyped")
    else:
        ax.text(0.5, 0.5, "no STR loci", ha="center", va="center", transform=ax.transAxes); ax.set_xticks([])
    ax.set_title("STR call rate per locus"); ax.tick_params(axis="x", rotation=45)
    save(fig, a.outdir, "str_call_rate",
         "Fraction of samples successfully genotyped at each STR locus (call rate).")

    # read support
    samples = sorted(support)
    fig, ax = plt.subplots(figsize=(max(5, 0.6 * max(1, len(samples)) + 2), 4))
    if samples and loci:
        for loc in loci:
            ax.plot(samples, [support[s].get(loc, 0) for s in samples], marker="o", label=loc)
        ax.set_ylabel("spanning reads"); ax.legend(fontsize=7)
    else:
        ax.text(0.5, 0.5, "no read-support data", ha="center", va="center", transform=ax.transAxes); ax.set_xticks([])
    ax.set_title("STR spanning-read support"); ax.tick_params(axis="x", rotation=45)
    save(fig, a.outdir, "str_read_support",
         "Spanning-read support per sample at each STR locus (genotype confidence).")
    print(f"[plot_str] wrote STR plots to {a.outdir}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""plot_stage0.py — Stage 0 QC plots (CallForge).

Renders, each as PNG+SVG with a one-line plain caption (embedded + sidecar):
  resource_availability.{png,svg}   resource manifest found/missing summary
  input_readcounts.{png,svg}        per-sample input read pairs (if counts given)

Captions are appended to <outdir>/captions.tsv  (figure<TAB>caption) so the
cohort dashboard can collate every figure's caption in one place.
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
    with open(os.path.join(outdir, "captions.tsv"), "a") as fh:
        fh.write(f"{stem}\t{caption}\n")


def plot_resources(manifest, outdir):
    res = manifest.get("resources", {})
    names, states = [], []
    statemap = {"ok": 2, "unindexed": 1, "build_mismatch": 1, "incomplete": 1,
                "missing": 0, "n/a": -1}
    for k, v in res.items():
        names.append(k)
        st = v.get("status", "missing")
        states.append(statemap.get(st, 0))
    colors = {2: "#2ca02c", 1: "#ff7f0e", 0: "#d62728", -1: "#999999"}
    fig, ax = plt.subplots(figsize=(7, 0.5 * len(names) + 1.5))
    ax.barh(names, [1] * len(names), color=[colors[s] for s in states])
    for i, (n, s) in enumerate(zip(names, states)):
        label = {2: "found", 1: "needs attn", 0: "missing", -1: "n/a"}[s]
        ax.text(0.5, i, label, ha="center", va="center", color="white", fontsize=9, fontweight="bold")
    ax.set_xlim(0, 1); ax.set_xticks([])
    ax.set_title(f"Resource availability — {manifest.get('species','')} {manifest.get('genome_build','')}")
    save(fig, outdir, "resource_availability",
         "Annotation/known-sites resources discovered locally; orange = present but "
         "needs indexing/build check, red = missing (stage will skip gracefully).")


def plot_readcounts(counts_tsv, outdir):
    samples, reads = [], []
    with open(counts_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            samples.append(r["sample_id"]); reads.append(int(r["read_pairs"]))
    fig, ax = plt.subplots(figsize=(max(5, 0.5 * len(samples) + 2), 4))
    ax.bar(samples, reads, color="#1f77b4")
    ax.set_ylabel("input read pairs"); ax.set_title("Per-sample input read pairs")
    ax.tick_params(axis="x", rotation=45)
    save(fig, outdir, "input_readcounts",
         "Raw paired-end read count per sample from the input FASTQs (sequencing yield).")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--readcounts", default=None)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    with open(a.manifest) as fh:
        manifest = json.load(fh)
    plot_resources(manifest, a.outdir)
    if a.readcounts and os.path.isfile(a.readcounts):
        plot_readcounts(a.readcounts, a.outdir)
    print(f"[plot_stage0] wrote plots to {a.outdir}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""plot_align_qc.py — Stage 3/4 alignment QC plots (CallForge).

  mapping_rate.{png,svg}    % mapped and % properly-paired per sample
  insert_size.{png,svg}     insert-size distribution per sample (from *.ishist.tsv)

PNG+SVG, one-line caption embedded + appended to <outdir>/captions.tsv.
"""
import argparse, glob, json, os
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
    ap.add_argument("--qc-json", nargs="+", required=True)
    ap.add_argument("--ishist", nargs="*", default=[])
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    qs = sorted((json.load(open(p)) for p in a.qc_json), key=lambda x: x["sample_id"])
    sids = [q["sample_id"] for q in qs]

    # mapping rate
    fig, ax = plt.subplots(figsize=(max(5, 0.6 * len(sids) + 2), 4))
    x = range(len(sids))
    ax.bar([i - 0.2 for i in x], [q.get("pct_mapped") or 0 for q in qs], width=0.4, label="% mapped")
    ax.bar([i + 0.2 for i in x], [q.get("pct_properly_paired") or 0 for q in qs], width=0.4,
           label="% properly paired")
    ax.set_xticks(list(x)); ax.set_xticklabels(sids, rotation=45, ha="right")
    ax.set_ylabel("% of reads"); ax.set_ylim(0, 100); ax.legend()
    ax.set_title("Alignment rate per sample")
    save(fig, a.outdir, "mapping_rate",
         "Percent of reads mapped and properly paired per sample (alignment success).")

    # insert size
    if a.ishist:
        fig, ax = plt.subplots(figsize=(6, 4))
        plotted = False
        for p in a.ishist:
            sid = os.path.basename(p).split(".")[0]
            xs, ys = [], []
            for ln in open(p):
                if ln.startswith("insert_size"):
                    continue
                f = ln.split("\t")
                if len(f) >= 2:
                    xs.append(int(f[0])); ys.append(int(f[1]))
            if xs:
                ax.plot(xs, ys, label=sid, lw=1); plotted = True
        if plotted:
            ax.set_xlabel("insert size (bp)"); ax.set_ylabel("read pairs")
            ax.set_xlim(0, 1000); ax.legend(fontsize=7)
            ax.set_title("Insert-size distribution")
            save(fig, a.outdir, "insert_size",
                 "Fragment insert-size distribution per sample (library size; capture efficiency).")
        else:
            plt.close(fig)
    print(f"[plot_align_qc] wrote alignment QC plots to {a.outdir}")


if __name__ == "__main__":
    main()

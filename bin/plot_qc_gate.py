#!/usr/bin/env python3
"""plot_qc_gate.py — Stage 5 QC gate plots (CallForge).

  qc_scorecard.{png,svg}        samples x checks grid, pass/FAIL/NA coloured + decision
  cohort_coverage_dist.{png,svg} distribution of mean target depth across the cohort
  contamination.{png,svg}        VerifyBamID2 FREEMIX per sample (or 'not assessed')
  sex_check.{png,svg}            declared vs inferred sex (or 'not assessed')

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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", required=True, help="cohort_qc_summary.json")
    ap.add_argument("--qc-json", nargs="+", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    summary = json.load(open(a.summary))
    qs = sorted((json.load(open(p)) for p in a.qc_json), key=lambda x: x["sample_id"])
    sids = [q["sample_id"] for q in qs]
    decisions = summary.get("decisions", {})

    # ── scorecard grid ──
    checks = ["depth", "dup", "on_target", "contamination", "sex"]
    code = {"pass": 2, "warn": 1, "NA": -1, "FAIL": 0}
    # recompute per-check status from thresholds applied earlier is not stored per sample;
    # derive a simple grid from values vs thresholds in the summary.
    thr = summary.get("thresholds", {})
    grid = np.full((len(sids), len(checks)), -1.0)
    for i, q in enumerate(qs):
        def st(ok, assessed=True):
            return 2 if (assessed and ok) else (0 if assessed else -1)
        d = q.get("mean_target_depth"); grid[i, 0] = st(d is not None and d >= thr.get("min_mean_target_depth", 30), d is not None)
        du = q.get("dup_rate");         grid[i, 1] = st(du is not None and du <= thr.get("max_dup_rate", 0.4), du is not None)
        ot = q.get("pct_on_target");    grid[i, 2] = st(ot is not None and ot >= thr.get("min_on_target", 0.4), ot is not None)
        if q.get("contamination_assessed"):
            con = q.get("contamination"); grid[i, 3] = st(con is not None and con <= thr.get("max_contamination", 0.03), con is not None)
        if q.get("sex_assessed"):
            grid[i, 4] = st(q.get("declared_sex") == q.get("sex_inferred"))
    cmap = plt.matplotlib.colors.ListedColormap(["#d62728", "#999999", "#ff7f0e", "#2ca02c"])  # FAIL,NA,warn,pass -> -1..2 remapped
    disp = np.where(grid < 0, 0.5, grid + 1)  # NA(-1)->0.5 bucket
    fig, ax = plt.subplots(figsize=(max(5, len(checks) + 2), max(3, 0.5 * len(sids) + 2)))
    ax.imshow(grid, aspect="auto", cmap=plt.matplotlib.colors.ListedColormap(
        ["#d62728", "#999999", "#2ca02c"]), vmin=-1, vmax=2)
    ax.set_xticks(range(len(checks))); ax.set_xticklabels(checks, rotation=30, ha="right")
    ax.set_yticks(range(len(sids)))
    ax.set_yticklabels([f"{s} [{decisions.get(s,'?')}]" for s in sids], fontsize=8)
    for i in range(len(sids)):
        for j in range(len(checks)):
            v = grid[i, j]
            ax.text(j, i, {2: "✓", 0: "✗", -1: "—"}.get(v, "?"), ha="center", va="center",
                    color="white", fontsize=9)
    ax.set_title(f"Per-sample QC scorecard  (PASS={summary['n_pass']}  QUARANTINE={summary['n_quarantine']})")
    save(fig, a.outdir, "qc_scorecard",
         "Per-sample QC checks (green=pass, red=fail, grey=not assessed); sample label shows the gate decision. Quarantined samples are kept and reported, not dropped.")

    # ── cohort coverage distribution ──
    depths = [q.get("mean_target_depth") for q in qs if q.get("mean_target_depth") is not None]
    fig, ax = plt.subplots(figsize=(6, 4))
    if depths:
        ax.hist(depths, bins=max(5, len(depths)), color="#1f77b4", edgecolor="white")
        ax.axvline(thr.get("min_mean_target_depth", 30), color="#d62728", ls="--", label="min depth")
        ax.legend(fontsize=7)
    ax.set_xlabel("mean target depth (x)"); ax.set_ylabel("samples")
    ax.set_title("Cohort mean target depth")
    save(fig, a.outdir, "cohort_coverage_dist",
         "Distribution of per-sample mean target depth across the cohort vs the gate threshold.")

    # ── contamination ──
    fig, ax = plt.subplots(figsize=(max(5, 0.6 * len(sids) + 2), 4))
    assessed = any(q.get("contamination_assessed") for q in qs)
    if assessed:
        ax.bar(sids, [q.get("contamination") or 0 for q in qs], color="#9467bd")
        ax.axhline(thr.get("max_contamination", 0.03), color="#d62728", ls="--", label="max")
        ax.set_ylabel("FREEMIX (contamination)"); ax.legend(fontsize=7)
    else:
        ax.text(0.5, 0.5, "contamination not assessed\n(no VerifyBamID2 SVD panel supplied)",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_xticks([])
    ax.set_title("Contamination per sample"); ax.tick_params(axis="x", rotation=45)
    save(fig, a.outdir, "contamination",
         "VerifyBamID2 FREEMIX contamination estimate per sample vs threshold (grey note if no SVD panel available).")

    # ── sex check ──
    fig, ax = plt.subplots(figsize=(max(5, 0.6 * len(sids) + 2), 3.5))
    if any(q.get("sex_assessed") for q in qs):
        ax.scatter(sids, [q.get("declared_sex") for q in qs], label="declared", marker="o")
        ax.scatter(sids, [q.get("sex_inferred") for q in qs], label="inferred", marker="x")
        ax.legend(fontsize=7)
    else:
        ax.text(0.5, 0.5, "sex not assessed\n(panel has no X/Y targets)",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_xticks([])
    ax.set_title("Sex concordance"); ax.tick_params(axis="x", rotation=45)
    save(fig, a.outdir, "sex_check",
         "Declared vs inferred sex per sample (grey note if the panel has no X/Y targets to infer from).")
    print(f"[plot_qc_gate] wrote QC-gate plots to {a.outdir}")


if __name__ == "__main__":
    main()

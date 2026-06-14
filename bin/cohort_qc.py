#!/usr/bin/env python3
"""cohort_qc.py — Stage 12 cohort QC: relatedness, ancestry PCA, sex, missingness.

Consumes the joint callset + somalier outputs and produces:
  relatedness_heatmap.{png,svg}  pairwise relatedness (somalier)
  ancestry_pca.{png,svg}         genotype PCA (numpy) of the cohort
  sex_check.{png,svg}            somalier-inferred vs declared sex (concordance)
  missingness.{png,svg}          per-sample genotype missingness
  cohort_qc.json                 tidy summary (relatedness pairs, sex calls, missingness, PCA)

Notes for a tight targeted panel are embedded in captions: relatedness/ancestry
rely on enough genome-wide (largely off-target) sites; on few on-target sites the
signal is weak. Real ancestry assignment uses somalier's 1000G-labelled reference.
"""
import argparse, csv, json, os, subprocess
import numpy as np
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


def samples_of(vcf):
    return subprocess.run(["bcftools", "query", "-l", vcf], check=True,
                          capture_output=True, text=True).stdout.split()


def genotype_matrix(vcf, samples):
    """variants x samples dosage matrix (0/1/2; nan=missing)."""
    out = subprocess.run(["bcftools", "query", "-f", "[%GT\t]\n", vcf],
                         check=True, capture_output=True, text=True).stdout.splitlines()
    rows = []
    for line in out:
        gts = line.rstrip("\t").split("\t")
        if len(gts) != len(samples):
            continue
        row = []
        for gt in gts:
            alleles = gt.replace("|", "/").split("/")
            if "." in alleles or gt in (".", "./.", ".|."):
                row.append(np.nan)
            else:
                row.append(sum(1 for x in alleles if x not in ("0", ".")))
        rows.append(row)
    return np.array(rows, dtype=float) if rows else np.zeros((0, len(samples)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--somalier-samples", default=None, help="somalier .samples.tsv")
    ap.add_argument("--somalier-pairs", default=None, help="somalier .pairs.tsv")
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    samples = samples_of(a.vcf)
    declared = {}
    with open(a.samplesheet) as fh:
        for r in csv.DictReader(fh):
            declared[r["sample_id"]] = (r.get("sex") or "U")

    summary = {"samples": samples, "n_samples": len(samples)}

    # ---- missingness ----
    M = genotype_matrix(a.vcf, samples)
    nvar = M.shape[0]
    miss = (np.isnan(M).sum(axis=0) / nvar) if nvar else np.zeros(len(samples))
    summary["n_variants"] = int(nvar)
    summary["missingness"] = {s: round(float(m), 4) for s, m in zip(samples, miss)}
    fig, ax = plt.subplots(figsize=(max(5, 0.6 * len(samples) + 2), 4))
    ax.bar(samples, miss * 100, color="#1f77b4"); ax.set_ylabel("% missing genotypes")
    ax.set_title("Per-sample missingness"); ax.tick_params(axis="x", rotation=45)
    save(fig, a.outdir, "missingness",
         "Per-sample fraction of no-call genotypes in the joint callset (sample quality).")

    # ---- ancestry PCA (genotype PCA via numpy SVD) ----
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    pcs = {}
    if nvar >= 2 and len(samples) >= 2:
        X = M.copy()
        colmean = np.nanmean(X, axis=0)
        inds = np.where(np.isnan(X))
        X[inds] = np.take(colmean, inds[1])           # mean-impute
        X = X - X.mean(axis=0)
        var = X.var(axis=0); X = X[:, :]              # keep all samples (columns)
        # PCA over samples: SVD of variants x samples
        try:
            U, S, Vt = np.linalg.svd(X, full_matrices=False)
            comp = Vt.T[:, :2] if Vt.shape[0] >= 2 else np.column_stack([Vt.T[:, 0], np.zeros(len(samples))])
            for i, s in enumerate(samples):
                pcs[s] = [float(comp[i, 0]), float(comp[i, 1])]
            ax.scatter(comp[:, 0], comp[:, 1], c="#1f77b4")
            for i, s in enumerate(samples):
                ax.annotate(s, (comp[i, 0], comp[i, 1]), fontsize=8)
            ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
        except np.linalg.LinAlgError:
            ax.text(0.5, 0.5, "PCA failed", ha="center", va="center", transform=ax.transAxes)
    else:
        ax.text(0.5, 0.5, f"too few variants for PCA (n={nvar})", ha="center", va="center", transform=ax.transAxes)
    summary["pca"] = pcs
    ax.set_title("Ancestry PCA (cohort genotypes)")
    save(fig, a.outdir, "ancestry_pca",
         "Genotype PCA of the cohort; on a tight panel few on-target sites give weak structure — "
         "real ancestry assignment uses somalier with the 1000G-labelled reference (off-target sites).")

    # ---- relatedness heatmap (somalier pairs) ----
    rel = {}
    if a.somalier_pairs and os.path.isfile(a.somalier_pairs):
        with open(a.somalier_pairs) as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            relcol = next((c for c in (rd.fieldnames or []) if "relatedness" in c.lower()), None)
            ida = next((c for c in rd.fieldnames if c in ("#sample_a", "sample_a")), rd.fieldnames[0])
            for r in rd:
                rel[(r[ida], r["sample_b"])] = float(r[relcol]) if relcol and r.get(relcol) else 0.0
    fig, ax = plt.subplots(figsize=(5, 4.5))
    n = len(samples)
    Mr = np.eye(n)
    for i, si in enumerate(samples):
        for j, sj in enumerate(samples):
            if i != j:
                Mr[i, j] = rel.get((si, sj), rel.get((sj, si), 0.0))
    im = ax.imshow(Mr, cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(range(n)); ax.set_xticklabels(samples, rotation=45, ha="right")
    ax.set_yticks(range(n)); ax.set_yticklabels(samples)
    fig.colorbar(im, ax=ax, label="relatedness")
    ax.set_title("Relatedness (somalier)" if rel else "Relatedness (somalier — no pairs)")
    summary["relatedness"] = {f"{a_}|{b_}": v for (a_, b_), v in rel.items()}
    save(fig, a.outdir, "relatedness_heatmap",
         "Pairwise somalier relatedness; synthetic test samples are near-identical so values run high "
         "(machinery proven) — a real cohort shows true relatedness structure.")

    # ---- sex inference from somalier X signal ----
    # somalier's `sex` column is the GIVEN sex (often -9/unknown without a ped); the
    # actual inference comes from X heterozygosity + X-depth-vs-autosomal-depth:
    #   Female = het calls on X + ~diploid X depth; Male = ~0 X het + ~half X depth.
    sex_calls = {}
    sex_evidence = {}
    if a.somalier_samples and os.path.isfile(a.somalier_samples):
        with open(a.somalier_samples) as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            idc = next((c for c in rd.fieldnames if c in ("#family_id", "#sample_id", "sample_id", "sample")), rd.fieldnames[0])
            def fnum(r, *keys):
                for k in keys:
                    if k in r and r[k] not in ("", ".", None):
                        try: return float(r[k])
                        except ValueError: pass
                return None
            for r in rd:
                sid = (r.get("sample_id") or r.get(idc) or "").strip()
                xhet = fnum(r, "X_het"); xdepth = fnum(r, "X_depth_mean")
                auto = fnum(r, "gt_depth_mean", "depth_mean")
                ratio = (xdepth / auto) if (xdepth is not None and auto) else None
                # X heterozygosity is the primary signal (two X copies -> het calls;
                # hemizygous male -> ~0 het). X/autosomal depth ratio is the tiebreaker.
                call = "U"
                if xhet is not None:
                    if xhet >= 2:    call = "F"
                    elif xhet == 0:  call = "M"
                    else:            call = "F" if (ratio is None or ratio > 0.7) else "M"
                elif ratio is not None:
                    call = "F" if ratio > 0.7 else "M"
                sex_calls[sid] = call
                sex_evidence[sid] = {"X_het": xhet, "X_depth": xdepth, "auto_depth": auto,
                                     "X_auto_ratio": round(ratio, 3) if ratio is not None else None}
    summary["sex_inferred"] = sex_calls
    summary["sex_evidence"] = sex_evidence
    summary["sex_declared"] = {s: declared.get(s, "U") for s in samples}
    fig, ax = plt.subplots(figsize=(max(5, 0.6 * n + 2), 3.8))
    if sex_calls:
        ax.scatter(samples, [declared.get(s, "U") for s in samples], marker="o", s=80, label="declared")
        ax.scatter(samples, [sex_calls.get(s, "U") for s in samples], marker="x", s=80, label="somalier")
        ax.legend(fontsize=7)
        concordant = sum(1 for s in samples if declared.get(s, "U") == sex_calls.get(s, "?"))
        summary["sex_concordant"] = concordant; summary["sex_total"] = len(samples)
    else:
        ax.text(0.5, 0.5, "no somalier sex calls", ha="center", va="center", transform=ax.transAxes); ax.set_xticks([])
    ax.set_title("Sex: declared vs somalier-inferred"); ax.tick_params(axis="x", rotation=45)
    save(fig, a.outdir, "sex_check",
         "Declared vs somalier-inferred sex per sample (somalier infers from X/Y signal incl. "
         "off-target reads — the backstop for an autosomal panel with no X/Y targets).")

    with open(os.path.join(a.outdir, "cohort_qc.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"[cohort_qc] {len(samples)} samples, {nvar} variants; "
          f"sex_inferred={sex_calls}; missingness={summary['missingness']}")


if __name__ == "__main__":
    main()

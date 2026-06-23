#!/usr/bin/env python3
"""simulate_cohort.py — emit the reviewable cohort PLAN + sample sheet (CallForge sim).

This is the Mac-buildable, dependency-free half of the harness. From cohort_design.py
it writes the artifacts that do NOT need a genome or a read simulator:

  * samplesheet.csv      the CallForge --input sheet (sample_id,fastq_1,fastq_2,sex,
                         phenotype,covariate_age,covariate_pc1,batch). FASTQ paths point
                         at <fixture-dir>/fastq/<sid>_R{1,2}.fastq.gz, which build_fixture.py
                         then populates.
  * cohort_plan.json     machine-readable plan: samples, the per-feature carrier matrix,
                         intended classes/consequences, and the checker thresholds.
  * cohort_plan.tsv      human-readable genotype matrix (samples x planted features).
  * gene_set.tsv         the chosen genes + their simulation roles + chromosome.

The PLAN is what a reviewer reads before any cluster run. Exact genomic coordinates and
ref/alt alleles are resolved later (they depend on the generated sequence) and written
into truth_manifest.json by build_fixture.py — the checker uses that resolved manifest.

stdlib-only. Deterministic.
"""
import argparse
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cohort_design as D


def write_samplesheet(path, fixture_dir):
    fqdir = os.path.join(fixture_dir, "fastq")
    cols = ["sample_id", "fastq_1", "fastq_2", "sex", "phenotype",
            "covariate_age", "covariate_pc1", "batch"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for s in D.SAMPLES:
            sid = s["sample_id"]
            w.writerow({
                "sample_id": sid,
                "fastq_1": os.path.join(fqdir, f"{sid}_R1.fastq.gz"),
                "fastq_2": os.path.join(fqdir, f"{sid}_R2.fastq.gz"),
                "sex": s["sex"],
                "phenotype": s["phenotype"],
                "covariate_age": s["covariate_age"],
                "covariate_pc1": s["covariate_pc1"],
                "batch": s["batch"],
            })


def feature_rows():
    """Flatten every planted feature into reviewable rows (class, gene, detail, carriers)."""
    rows = []
    for v in D.CODING_VARIANTS:
        rows.append({"feature": v["id"], "class": "coding", "gene": v["gene"],
                     "detail": D.EXPECTED_CSQ[v["kind"]], "carriers": v["carriers"]})
    for v in D.BURDEN_VARIANTS:
        rows.append({"feature": v["id"], "class": "burden_coding", "gene": v["gene"],
                     "detail": D.EXPECTED_CSQ[v["kind"]], "carriers": v["carriers"]})
    for v in D.PROMOTER_VARIANTS:
        rows.append({"feature": v["id"], "class": "promoter", "gene": v["gene"],
                     "detail": D.EXPECTED_CSQ["promoter"], "carriers": v["carriers"]})
    s = D.STR_LOCUS
    rows.append({"feature": f"str_{s['gene']}", "class": "str", "gene": s["gene"],
                 "detail": f"({s['motif']})n normal={s['normal_units']} expanded={s['expanded_units']}",
                 "carriers": s["carriers"]})
    for e in D.CNV_EVENTS:
        rows.append({"feature": e["id"], "class": "cnv", "gene": e["gene"],
                     "detail": f"{e['state']} fold={e['fold']}", "carriers": e["carriers"]})
    return rows


def write_plan_tsv(path, rows):
    with open(path, "w") as fh:
        fh.write("feature\tclass\tgene\tdetail\t" + "\t".join(D.SAMPLE_IDS) + "\n")
        for r in rows:
            carriers = set(r["carriers"])
            cells = ["1" if sid in carriers else "." for sid in D.SAMPLE_IDS]
            fh.write(f"{r['feature']}\t{r['class']}\t{r['gene']}\t{r['detail']}\t"
                     + "\t".join(cells) + "\n")


def write_gene_set_tsv(path):
    with open(path, "w") as fh:
        fh.write("gene\tchrom\troles\n")
        for g, meta in D.GENES.items():
            fh.write(f"{g}\t{meta['chrom']}\t{','.join(meta['roles'])}\n")


def write_plan_json(path, fixture_dir, rows):
    plan = {
        "seed": D.SEED,
        "fixture_dir": fixture_dir,
        "n_samples": len(D.SAMPLES),
        "cases": D.CASES,
        "controls": D.CONTROLS,
        "samples": D.SAMPLES,
        "genes": D.GENES,
        "burden_gene": D.BURDEN_GENE,
        "cnv_target_genes": D.CNV_TARGET_GENES,
        "str_locus": D.STR_LOCUS,
        "features": rows,
        "thresholds": D.THRESHOLDS,
        "note": ("PLAN only — exact coords/ref/alt are resolved into truth_manifest.json "
                 "by build_fixture.py; the checker uses that resolved manifest."),
    }
    with open(path, "w") as fh:
        json.dump(plan, fh, indent=2)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default=".", help="where to write the plan files (default: cwd)")
    ap.add_argument("--fixture-dir", default=None,
                    help="dir build_fixture.py will write the reference + FASTQs into "
                         "(used to compose FASTQ paths in the sample sheet; default: <outdir>)")
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    fixture_dir = os.path.abspath(a.fixture_dir or a.outdir)

    rows = feature_rows()
    write_samplesheet(os.path.join(a.outdir, "samplesheet.csv"), fixture_dir)
    write_plan_tsv(os.path.join(a.outdir, "cohort_plan.tsv"), rows)
    write_gene_set_tsv(os.path.join(a.outdir, "gene_set.tsv"))
    write_plan_json(os.path.join(a.outdir, "cohort_plan.json"), fixture_dir, rows)

    n_feat = len(rows)
    print(f"[simulate_cohort] wrote plan to {a.outdir}: "
          f"{len(D.SAMPLES)} samples, {n_feat} planted features "
          f"({len(D.CASES)} cases / {len(D.CONTROLS)} controls); "
          f"FASTQ paths -> {fixture_dir}/fastq/")


if __name__ == "__main__":
    main()

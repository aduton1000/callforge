#!/usr/bin/env python3
"""simulate_cohort.py — reviewable cohort PLAN + sample sheet (host, stdlib-only).

The host-buildable half of the REAL-reference harness. From cohort_design.py it writes
the artifacts that need NO genome and NO read simulator:

  * samplesheet.csv      CallForge --input sheet (sample_id,fastq_1,fastq_2,sex,phenotype,
                         covariate_age,covariate_pc1,batch). FASTQ paths point at
                         <fixture-dir>/fastq/<sid>_R{1,2}.fastq.gz, populated on the cluster
                         by generate_reads.sh (dwgsim against real GRCh38).
  * cohort_plan.json     machine-readable plan: samples, the per-feature carrier matrix,
                         intended classes/consequences, real-coordinate WINDOWS, the
                         per-stage real-reference dependency map, and the thresholds.
  * cohort_plan.tsv      human-readable genotype matrix (samples x planted features).
  * gene_set.tsv         chosen genes + roles + REAL GRCh38 chrom/window/transcript.
  * manifest_skeleton.json  the truth-manifest shape BEFORE cluster resolution (exact
                         coords/ref/alt/consequence are filled by resolve_variants.py).

This is what a reviewer reads before any cluster run. Exact planted positions are resolved
against the real FASTA + gene model on the cluster (resolve_variants.py) — NOT here.

Deterministic. stdlib-only.
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
            w.writerow({"sample_id": sid,
                        "fastq_1": os.path.join(fqdir, f"{sid}_R1.fastq.gz"),
                        "fastq_2": os.path.join(fqdir, f"{sid}_R2.fastq.gz"),
                        "sex": s["sex"], "phenotype": s["phenotype"],
                        "covariate_age": s["covariate_age"],
                        "covariate_pc1": s["covariate_pc1"], "batch": s["batch"]})


def feature_rows():
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
                 "detail": f"{s['chrom']}:{s['ref_region'][0]}-{s['ref_region'][1]} "
                           f"({s['motif']})n {s['normal_units']}->{s['expanded_units']}",
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
        fh.write("gene\tchrom\tstrand\twindow_start\twindow_end\ttranscript\troles\n")
        for g, m in D.GENES.items():
            ws, we = m["window"]
            fh.write(f"{g}\t{m['chrom']}\t{m['strand']}\t{ws}\t{we}\t{m['tx']}\t{','.join(m['roles'])}\n")


def write_plan_json(path, fixture_dir, rows):
    plan = {
        "seed": D.SEED, "genome_build": D.GENOME_BUILD, "fixture_dir": fixture_dir,
        "n_samples": len(D.SAMPLES), "cases": D.CASES, "controls": D.CONTROLS,
        "samples": D.SAMPLES, "genes": D.GENES, "str_locus": D.STR_LOCUS,
        "sex_regions": D.SEX_REGIONS, "burden_gene": D.BURDEN_GENE,
        "cnv_target_genes": D.CNV_TARGET_GENES, "features": rows,
        "stage_refs": D.STAGE_REFS, "thresholds": D.THRESHOLDS,
        "note": ("PLAN only — exact coords/ref/alt/consequence are RESOLVED against the "
                 "real GRCh38 FASTA + gene model by resolve_variants.py on the cluster and "
                 "written into truth_manifest.json; the checker uses that resolved manifest."),
    }
    with open(path, "w") as fh:
        json.dump(plan, fh, indent=2)


def write_manifest_skeleton(path):
    """The truth-manifest shape before cluster resolution (coords = null placeholders)."""
    skel = {
        "_status": "SKELETON — resolve on cluster (resolve_variants.py) before checking",
        "seed": D.SEED, "genome_build": D.GENOME_BUILD,
        "genome_fasta": None, "target_bed": None, "gene_metadata": None,
        "str_catalog": None, "somalier_sites": None,
        "samples": D.SAMPLES, "cases": D.CASES, "controls": D.CONTROLS,
        "coding_variants": [{"id": v["id"], "gene": v["gene"], "kind": v["kind"],
                             "consequence": D.EXPECTED_CSQ[v["kind"]],
                             "carriers": v["carriers"],
                             "contig": None, "pos": None, "ref": None, "alt": None,
                             "vartype": None} for v in D.ALL_CODING],
        "promoter_variants": [{"id": v["id"], "gene": v["gene"],
                               "consequence": D.EXPECTED_CSQ["promoter"],
                               "carriers": v["carriers"],
                               "contig": None, "pos": None, "ref": None, "alt": None}
                              for v in D.PROMOTER_VARIANTS],
        "str": {"gene": D.STR_LOCUS["gene"], "chrom": D.STR_LOCUS["chrom"],
                "motif": D.STR_LOCUS["motif"], "ref_region": D.STR_LOCUS["ref_region"],
                "normal_units": D.STR_LOCUS["normal_units"],
                "expanded_units": D.STR_LOCUS["expanded_units"],
                "carriers": D.STR_LOCUS["carriers"], "locus_id": None, "repeat_region": None},
        "cnv_events": D.CNV_EVENTS,
        "burden": {"gene": D.BURDEN_GENE,
                   "qualifying_feature_ids": [v["id"] for v in D.BURDEN_VARIANTS],
                   "case_carriers": sorted({c for v in D.BURDEN_VARIANTS
                                            for c in v["carriers"] if c in D.CASES}),
                   "control_carriers": sorted({c for v in D.BURDEN_VARIANTS
                                               for c in v["carriers"] if c in D.CONTROLS})},
        "cohortqc": {"sex_by_sample": {s["sample_id"]: s["sex"] for s in D.SAMPLES}},
        "stage_refs": D.STAGE_REFS, "thresholds": D.THRESHOLDS,
    }
    with open(path, "w") as fh:
        json.dump(skel, fh, indent=2)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default=".", help="where to write the plan files")
    ap.add_argument("--fixture-dir", default=None,
                    help="dir the cluster writes reference + FASTQs into (composes FASTQ "
                         "paths in the sample sheet; default: <outdir>)")
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    fixture_dir = os.path.abspath(a.fixture_dir or a.outdir)
    rows = feature_rows()
    write_samplesheet(os.path.join(a.outdir, "samplesheet.csv"), fixture_dir)
    write_plan_tsv(os.path.join(a.outdir, "cohort_plan.tsv"), rows)
    write_gene_set_tsv(os.path.join(a.outdir, "gene_set.tsv"))
    write_plan_json(os.path.join(a.outdir, "cohort_plan.json"), fixture_dir, rows)
    write_manifest_skeleton(os.path.join(a.outdir, "manifest_skeleton.json"))
    print(f"[simulate_cohort] wrote plan to {a.outdir}: {len(D.SAMPLES)} samples, "
          f"{len(rows)} planted features ({len(D.CASES)} cases / {len(D.CONTROLS)} controls); "
          f"FASTQ paths -> {fixture_dir}/fastq/ (populated on the cluster).")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""paralog_flag.py — Stage 10 paralog-aware flagging (CallForge).

For CaptureForge-flagged paralog genes (e.g. CD209/CASP1/CR1/HP/CFH), marks
variants falling in those paralog-ambiguous target regions and assigns a
multimapping-aware confidence from mapping quality (MQ): MQ>=min_mq -> 'confident'
(uniquely callable, recovered); else 'flagged_low_confidence'. The labels are
WRITTEN BACK into the VCF (INFO/PARALOG_GENE, INFO/PARALOG_CONF) so a variant in
a paralog region visibly carries its flag downstream — the "recover what's
uniquely callable, flag the rest" approach.

Emits: paralog.annotated.vcf.gz(.tbi), paralog_variants.tsv, paralog_summary.json
Requires bcftools/bgzip/tabix on PATH.
"""
import argparse, csv, json, os, subprocess, sys
from collections import defaultdict


def sh(cmd, **kw):
    return subprocess.run(cmd, check=True, capture_output=True, text=True, **kw).stdout


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--bed", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--min-mq", type=float, default=50.0)
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    od = a.outdir

    # paralog genes from metadata
    paralog = set()
    with open(a.metadata) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("is_paralog") == "yes":
                paralog.add(r["gene"])

    # paralog regions = target intervals of paralog genes
    regbed = os.path.join(od, "paralog_regions.bed")
    n_reg = 0
    with open(a.bed) as fh, open(regbed, "w") as out:
        rows = []
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            gene = (f[3] if len(f) > 3 else "NA|x").split("|")[0]
            if gene in paralog:
                rows.append((f[0], int(f[1]), int(f[2]), gene)); n_reg += 1
        for (c, s, e, g) in sorted(rows):
            out.write(f"{c}\t{s}\t{e}\t{g}\n")

    out_vcf = os.path.join(od, "paralog.annotated.vcf.gz")
    per_gene = defaultdict(lambda: {"confident": 0, "flagged_low_confidence": 0})
    var_rows = []

    if n_reg == 0:
        # no paralog regions on this panel -> pass VCF through unchanged
        sh(["bcftools", "view", a.vcf, "-Oz", "-o", out_vcf])
        sh(["tabix", "-f", "-p", "vcf", out_vcf])
    else:
        sh(["bgzip", "-f", regbed]); sh(["tabix", "-f", "-p", "bed", regbed + ".gz"])
        hdr = os.path.join(od, "paralog.hdr")
        with open(hdr, "w") as fh:
            fh.write('##INFO=<ID=PARALOG_GENE,Number=1,Type=String,Description="Paralog-ambiguous gene region (CaptureForge paralog flag)">\n')
            fh.write('##INFO=<ID=PARALOG_CONF,Number=1,Type=String,Description="Paralog-region confidence: confident|flagged_low_confidence (by MQ)">\n')
        # step 1: add PARALOG_GENE for variants overlapping paralog regions
        v1 = os.path.join(od, "_p1.vcf.gz")
        sh(["bcftools", "annotate", "-a", regbed + ".gz", "-h", hdr,
            "-c", "CHROM,FROM,TO,PARALOG_GENE", a.vcf, "-Oz", "-o", v1])
        sh(["tabix", "-f", "-p", "vcf", v1])
        # step 2: classify confidence by MQ; build a per-site annotation
        q = sh(["bcftools", "query", "-f", "%CHROM\t%POS\t%INFO/PARALOG_GENE\t%INFO/MQ\n", v1])
        annot = os.path.join(od, "_conf.tsv")
        with open(annot, "w") as fh:
            for line in q.splitlines():
                c, p, gene, mq = (line.split("\t") + ["", "", "", ""])[:4]
                if gene in (".", "", None):
                    continue
                try:
                    conf = "confident" if float(mq) >= a.min_mq else "flagged_low_confidence"
                except ValueError:
                    conf = "flagged_low_confidence"     # missing MQ -> conservative
                per_gene[gene][conf] += 1
                var_rows.append((c, p, gene, conf, mq))
                fh.write(f"{c}\t{p}\t{conf}\n")
        if var_rows:
            sh(["bgzip", "-f", annot]); sh(["tabix", "-f", "-s", "1", "-b", "2", "-e", "2", annot + ".gz"])
            hdr2 = os.path.join(od, "paralog_conf.hdr")
            with open(hdr2, "w") as fh:
                fh.write('##INFO=<ID=PARALOG_CONF,Number=1,Type=String,Description="Paralog-region confidence: confident|flagged_low_confidence (by MQ)">\n')
            sh(["bcftools", "annotate", "-a", annot + ".gz", "-h", hdr2,
                "-c", "CHROM,POS,PARALOG_CONF", v1, "-Oz", "-o", out_vcf])
        else:
            sh(["bcftools", "view", v1, "-Oz", "-o", out_vcf])
        sh(["tabix", "-f", "-p", "vcf", out_vcf])

    with open(os.path.join(od, "paralog_variants.tsv"), "w") as fh:
        fh.write("chrom\tpos\tparalog_gene\tconfidence\tMQ\n")
        for row in var_rows:
            fh.write("\t".join(str(x) for x in row) + "\n")

    summary = {"paralog_genes": sorted(paralog), "n_paralog_regions": n_reg,
               "n_paralog_variants": len(var_rows),
               "n_confident": sum(g["confident"] for g in per_gene.values()),
               "n_flagged": sum(g["flagged_low_confidence"] for g in per_gene.values()),
               "per_gene": {g: dict(v) for g, v in per_gene.items()},
               "min_mq": a.min_mq}
    with open(os.path.join(od, "paralog_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"[paralog_flag] {len(var_rows)} variants in {n_reg} paralog regions "
          f"({summary['n_confident']} confident, {summary['n_flagged']} flagged); genes={sorted(per_gene)}")


if __name__ == "__main__":
    main()

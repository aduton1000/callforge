#!/usr/bin/env python3
"""filter_counts.py — Stage 7 hard-filter before/after counts (CallForge).

Reads the FILTER-annotated joint VCF (via bcftools query) and tallies SNP/indel
totals vs PASS, plus the breakdown of filter reasons. Emits JSON for the report
and the filter QC plot. The applied thresholds are recorded verbatim.
"""
import argparse, json, subprocess, sys
from collections import Counter


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True, help="FILTER-annotated joint VCF")
    ap.add_argument("--snp-expr", default="")
    ap.add_argument("--indel-expr", default="")
    ap.add_argument("--out-json", required=True)
    a = ap.parse_args()

    try:
        out = subprocess.run(["bcftools", "query", "-f", "%TYPE\t%FILTER\n", a.vcf],
                             check=True, capture_output=True, text=True).stdout
    except (subprocess.SubprocessError, OSError) as e:
        sys.exit(f"[filter_counts] bcftools query failed: {e}")

    tot = Counter(); passed = Counter(); reasons = Counter()
    for line in out.splitlines():
        typ, _, filt = line.partition("\t")
        klass = "snp" if typ.upper() == "SNP" else "indel" if typ.upper() in ("INDEL", "MNP", "OTHER") else "other"
        tot[klass] += 1
        if filt in ("PASS", "."):
            passed[klass] += 1
        else:
            for r in filt.split(";"):
                reasons[r] += 1

    def rate(k):
        return round(passed[k] / tot[k], 4) if tot[k] else None
    res = {
        "vcf": a.vcf,
        "snp_total": tot["snp"], "snp_pass": passed["snp"], "snp_pass_rate": rate("snp"),
        "indel_total": tot["indel"], "indel_pass": passed["indel"], "indel_pass_rate": rate("indel"),
        "other_total": tot["other"], "other_pass": passed["other"],
        "total": sum(tot.values()), "pass_total": sum(passed.values()),
        "pass_rate_overall": round(sum(passed.values()) / sum(tot.values()), 4) if sum(tot.values()) else None,
        "filter_reasons": dict(reasons),
        "thresholds": {"snp_filter_expr": a.snp_expr, "indel_filter_expr": a.indel_expr},
    }
    with open(a.out_json, "w") as fh:
        json.dump(res, fh, indent=2)
    print(f"[filter_counts] SNP {res['snp_pass']}/{res['snp_total']} pass; "
          f"indel {res['indel_pass']}/{res['indel_total']} pass; "
          f"overall {res['pass_total']}/{res['total']}")


if __name__ == "__main__":
    main()

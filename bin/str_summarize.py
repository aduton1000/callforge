#!/usr/bin/env python3
"""str_summarize.py — Stage 9 STR genotyping summary (CallForge).

Parses ExpansionHunter per-sample JSON outputs into a tidy table + summary:
per locus per sample — genotype (repeat units), allele sizes, read support
(spanning reads); and per-locus call rate across the cohort.

Inputs: --eh-json <sample>.json ... (one per sample; name prefix = sample id)
Outputs: str_calls.tsv, str_summary.json
"""
import argparse, json, os, re
from collections import defaultdict


def alleles_from_genotype(gt):
    if not gt:
        return []
    out = []
    for x in str(gt).split("/"):
        x = x.strip()
        if x.isdigit():
            out.append(int(x))
    return out


def spanning_count(v):
    # CountsOfSpanningReads like "(20, 14), (21, 3)" -> total spanning reads
    s = v.get("CountsOfSpanningReads", "")
    return sum(int(m) for m in re.findall(r"\(\s*\d+\s*,\s*(\d+)\s*\)", s)) if s else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--eh-json", nargs="+", required=True)
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    rows = []
    loci = set()
    samples = set()
    per_locus_called = defaultdict(int)
    for path in a.eh_json:
        sample = os.path.basename(path).split(".")[0]
        samples.add(sample)
        try:
            data = json.load(open(path))
        except (json.JSONDecodeError, OSError):
            continue
        for locus_id, lr in (data.get("LocusResults") or {}).items():
            loci.add(locus_id)
            for vid, v in (lr.get("Variants") or {}).items():
                gt = v.get("Genotype")
                alleles = alleles_from_genotype(gt)
                called = bool(alleles)
                if called:
                    per_locus_called[locus_id] += 1
                rows.append({"sample": sample, "locus": locus_id, "variant": vid,
                             "genotype": gt or "./.", "allele1": alleles[0] if alleles else "",
                             "allele2": alleles[1] if len(alleles) > 1 else "",
                             "repeat_unit": v.get("RepeatUnit", ""),
                             "spanning_reads": spanning_count(v),
                             "called": "yes" if called else "no"})

    cols = ["sample", "locus", "variant", "genotype", "allele1", "allele2",
            "repeat_unit", "spanning_reads", "called"]
    with open(os.path.join(a.outdir, "str_calls.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    n = len(samples) or 1
    call_rate = {loc: round(per_locus_called[loc] / n, 3) for loc in loci}
    summary = {"n_samples": len(samples), "loci": sorted(loci),
               "call_rate_per_locus": call_rate,
               "n_calls": sum(1 for r in rows if r["called"] == "yes")}
    with open(os.path.join(a.outdir, "str_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"[str_summarize] {len(samples)} samples x {len(loci)} loci; "
          f"call_rate={call_rate}")


if __name__ == "__main__":
    main()

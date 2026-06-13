#!/usr/bin/env python3
"""ingest_captureforge.py — Stage 0 CaptureForge handoff ingestion (CallForge).

Derives the per-gene metadata table that drives the matching downstream steps
(CNV / STR / paralog-aware / burden) from the CaptureForge handoff:

  final_covered_targets.bed   name column = 'GENE|class' -> design classes per gene
  metrics.json (optional)     cnv_bin_spacing[gene].depth_callable -> CNV callability
                              qc_gates.low_coverage_genes          -> low-coverage flag
  baits.csv   (optional)      off_target_flag aggregated per gene  -> specificity signal
  --paralog-genes (param)     comma list of paralog-ambiguous genes (Stage 10)

Emits a gene_metadata.tsv (one row per gene) + a JSON summary. Works from the
BED alone (CNV callability -> 'unknown' when metrics.json is absent), so it is
organism-/panel-generic.
"""
import argparse, csv, json, os, sys
from collections import defaultdict


def parse_bed(path):
    genes = defaultdict(lambda: defaultdict(int))   # gene -> class -> n_intervals
    spans = defaultdict(lambda: defaultdict(lambda: [None, None]))  # gene -> class -> [chrom, [min,max]]
    loci = defaultdict(lambda: defaultdict(list))    # gene -> class -> [(chrom,start,end)]
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, start, end = f[0], int(f[1]), int(f[2])
            name = f[3] if len(f) > 3 else "NA|coding"
            gene, _, klass = name.partition("|")
            klass = klass or "coding"
            genes[gene][klass] += 1
            loci[gene][klass].append((chrom, start, end))
    return genes, loci


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed", required=True)
    ap.add_argument("--metrics", default=None, help="CaptureForge metrics.json")
    ap.add_argument("--baits", default=None, help="CaptureForge baits.csv")
    ap.add_argument("--paralog-genes", default="", help="comma list of paralog genes")
    ap.add_argument("--burden-group-map", default=None,
                    help="optional TSV: gene<TAB>burden_group (else burden_group=gene)")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-json", required=True)
    a = ap.parse_args()

    genes, loci = parse_bed(a.bed)
    if not genes:
        sys.stderr.write(f"[ingest_captureforge] ERROR: no intervals parsed from {a.bed}\n")
        sys.exit(1)

    # metrics.json — CNV callability + low-coverage flags
    cnv_callable, cnv_reason, low_cov = {}, {}, set()
    metrics = {}
    if a.metrics and os.path.isfile(a.metrics):
        with open(a.metrics) as fh:
            metrics = json.load(fh)
        for g, m in (metrics.get("cnv_bin_spacing") or {}).items():
            cnv_callable[g] = bool(m.get("depth_callable"))
            cnv_reason[g] = (f"median_gap={m.get('median_gap')} max_gap={m.get('max_gap')} "
                             f"cv_gap={m.get('cv_gap')}")
        for entry in (metrics.get("qc_gates", {}).get("low_coverage_genes") or []):
            low_cov.add(entry.get("gene"))

    # baits.csv — per-gene off-target (specificity) signal
    off_frac = {}
    if a.baits and os.path.isfile(a.baits):
        tot, off = defaultdict(int), defaultdict(int)
        with open(a.baits, newline="") as fh:
            for r in csv.DictReader(fh):
                g = r.get("gene")
                if not g:
                    continue
                tot[g] += 1
                try:
                    off[g] += 1 if int(r.get("off_target_flag", 0)) else 0
                except ValueError:
                    pass
        off_frac = {g: round(off[g] / tot[g], 4) for g in tot if tot[g]}

    paralogs = {g.strip() for g in a.paralog_genes.split(",") if g.strip()}

    burden_map = {}
    if a.burden_group_map and os.path.isfile(a.burden_group_map):
        with open(a.burden_group_map) as fh:
            for line in fh:
                if line.strip() and not line.startswith("#"):
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= 2:
                        burden_map[parts[0]] = parts[1]

    cols = ["gene", "classes", "n_intervals", "n_coding", "n_promoter", "n_anchor",
            "n_str", "n_cnv", "is_cnv_target", "cnv_callable", "cnv_callable_reason",
            "is_str_target", "str_loci", "is_paralog", "off_target_frac",
            "low_coverage", "burden_group"]
    rows = []
    for g in sorted(genes):
        klasses = genes[g]
        is_cnv = "cnv" in klasses
        callable_str = "n/a"
        if is_cnv:
            callable_str = ({True: "yes", False: "no"}.get(cnv_callable[g]) if g in cnv_callable
                            else "unknown")
        is_str = "str" in klasses
        str_loci = ";".join(f"{c}:{s}-{e}" for (c, s, e) in loci[g].get("str", []))
        rows.append({
            "gene": g,
            "classes": ",".join(sorted(klasses)),
            "n_intervals": sum(klasses.values()),
            "n_coding": klasses.get("coding", 0),
            "n_promoter": klasses.get("promoter", 0),
            "n_anchor": klasses.get("anchor", 0),
            "n_str": klasses.get("str", 0),
            "n_cnv": klasses.get("cnv", 0),
            "is_cnv_target": "yes" if is_cnv else "no",
            "cnv_callable": callable_str,
            "cnv_callable_reason": cnv_reason.get(g, ""),
            "is_str_target": "yes" if is_str else "no",
            "str_loci": str_loci,
            "is_paralog": "yes" if g in paralogs else "no",
            "off_target_frac": off_frac.get(g, ""),
            "low_coverage": "yes" if g in low_cov else "no",
            "burden_group": burden_map.get(g, g),
        })

    with open(a.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    summary = {
        "n_genes": len(rows),
        "classes_present": sorted({k for g in genes for k in genes[g]}),
        "cnv_targets": [r["gene"] for r in rows if r["is_cnv_target"] == "yes"],
        "cnv_callable_yes": [r["gene"] for r in rows if r["cnv_callable"] == "yes"],
        "cnv_callable_no":  [r["gene"] for r in rows if r["cnv_callable"] == "no"],
        "cnv_callable_unknown": [r["gene"] for r in rows if r["cnv_callable"] == "unknown"],
        "str_targets": [r["gene"] for r in rows if r["is_str_target"] == "yes"],
        "paralog_genes": sorted(paralogs),
        "low_coverage_genes": sorted(low_cov),
        "metrics_used": bool(metrics),
        "baits_used": bool(off_frac),
    }
    with open(a.out_json, "w") as fh:
        json.dump(summary, fh, indent=2)

    print(f"[ingest_captureforge] {len(rows)} genes; classes={summary['classes_present']}; "
          f"CNV targets={len(summary['cnv_targets'])} "
          f"(callable={len(summary['cnv_callable_yes'])}, not={len(summary['cnv_callable_no'])}, "
          f"unknown={len(summary['cnv_callable_unknown'])}); "
          f"STR={len(summary['str_targets'])}; paralog={len(paralogs)}")


if __name__ == "__main__":
    main()

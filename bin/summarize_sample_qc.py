#!/usr/bin/env python3
"""summarize_sample_qc.py — per-sample post-alignment QC summary (CallForge, Stage 4).

Parses the standard tool outputs into one tidy per-sample JSON that the QC gate
(Stage 5) and the plots consume. Every optional input degrades to null/NA rather
than failing, so the summary is always produced.

Inputs (any may be omitted):
  --flagstat     samtools flagstat
  --stats        samtools stats               (insert size, error rate)
  --markdup      Picard MarkDuplicates metrics (PERCENT_DUPLICATION)
  --hsmetrics    Picard CollectHsMetrics       (on-target, fold-enrichment, depth, %≥Nx)
  --mosdepth-summary  mosdepth *.summary.txt   (region mean depth)
  --mosdepth-regions  mosdepth *.regions.bed.gz (per-target -> per-gene depth)
  --verifybamid  VerifyBamID2 *.selfSM         (FREEMIX contamination)
  --sex          declared sex (M/F/U) from the sample sheet
"""
import argparse, gzip, json, os, re, sys
from collections import defaultdict


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def parse_flagstat(p):
    if not p or not os.path.isfile(p):
        return {}
    txt = open(p).read()
    def grab(pat):
        m = re.search(pat, txt)
        return int(m.group(1)) if m else None
    total = grab(r"(\d+) \+ \d+ in total")
    mapped = grab(r"(\d+) \+ \d+ mapped")
    pp = grab(r"(\d+) \+ \d+ properly paired")
    out = {"total_reads": total, "mapped_reads": mapped, "properly_paired_reads": pp}
    if total:
        if mapped is not None: out["pct_mapped"] = round(100 * mapped / total, 3)
        if pp is not None:     out["pct_properly_paired"] = round(100 * pp / total, 3)
    return out


def parse_stats(p):
    if not p or not os.path.isfile(p):
        return {}, []
    out, hist = {}, []
    for line in open(p):
        if line.startswith("SN\t"):
            f = line.rstrip("\n").split("\t")
            key = f[1].rstrip(":")
            if key == "insert size average": out["insert_size_mean"] = num(f[2])
            elif key == "insert size standard deviation": out["insert_size_sd"] = num(f[2])
            elif key == "error rate": out["error_rate"] = num(f[2])
            elif key == "reads mapped": out["sn_reads_mapped"] = num(f[2])
        elif line.startswith("IS\t"):
            f = line.split("\t")
            hist.append((int(f[1]), int(f[2])))   # insert_size, count(all pairs)
    return out, hist


def parse_picard_metrics(p):
    """Return the single METRICS CLASS data row as a dict (col -> value)."""
    if not p or not os.path.isfile(p):
        return {}
    lines = open(p).read().splitlines()
    for i, ln in enumerate(lines):
        if ln.startswith("## METRICS CLASS"):
            header = lines[i + 1].split("\t")
            values = lines[i + 2].split("\t") if i + 2 < len(lines) else []
            return dict(zip(header, values))
    return {}


def parse_mosdepth_summary(p):
    if not p or not os.path.isfile(p):
        return {}
    for line in open(p):
        f = line.rstrip("\n").split("\t")
        if f and f[0] == "total_region":
            return {"mean_region_depth": num(f[3])}   # chrom,len,bases,mean,min,max
    return {}


def parse_mosdepth_regions(p):
    """Per-gene mean depth from mosdepth regions BED (name col = GENE|class)."""
    if not p or not os.path.isfile(p):
        return {}
    opener = gzip.open if p.endswith(".gz") else open
    by_gene = defaultdict(lambda: [0.0, 0])  # gene -> [weighted_depth_sum, total_len]
    with opener(p, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            c, s, e, name, depth = f[0], int(f[1]), int(f[2]), f[3], num(f[4])
            gene = name.split("|")[0]
            L = e - s
            by_gene[gene][0] += (depth or 0) * L
            by_gene[gene][1] += L
    return {g: round(v[0] / v[1], 2) for g, v in by_gene.items() if v[1]}


def parse_verifybamid(p):
    if not p or not os.path.isfile(p):
        return {"contamination": None, "contamination_assessed": False}
    lines = open(p).read().splitlines()
    if len(lines) >= 2:
        header = lines[0].split("\t")
        row = dict(zip(header, lines[1].split("\t")))
        fm = num(row.get("FREEMIX"))
        return {"contamination": fm, "contamination_assessed": fm is not None}
    return {"contamination": None, "contamination_assessed": False}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--flagstat"); ap.add_argument("--stats")
    ap.add_argument("--markdup"); ap.add_argument("--hsmetrics")
    ap.add_argument("--mosdepth-summary"); ap.add_argument("--mosdepth-regions")
    ap.add_argument("--verifybamid"); ap.add_argument("--sex", default="U")
    ap.add_argument("--out-json", required=True)
    a = ap.parse_args()

    q = {"sample_id": a.sample_id, "declared_sex": a.sex}
    q.update(parse_flagstat(a.flagstat))
    stats, hist = parse_stats(a.stats)
    q.update(stats)

    md = parse_picard_metrics(a.markdup)
    if md:
        q["dup_rate"] = num(md.get("PERCENT_DUPLICATION"))

    hs = parse_picard_metrics(a.hsmetrics)
    if hs:
        q["pct_on_target"] = num(hs.get("PCT_SELECTED_BASES"))
        q["fold_enrichment"] = num(hs.get("FOLD_ENRICHMENT"))
        q["mean_target_depth"] = num(hs.get("MEAN_TARGET_COVERAGE"))
        q["pct_target_10x"] = num(hs.get("PCT_TARGET_BASES_10X"))
        q["pct_target_20x"] = num(hs.get("PCT_TARGET_BASES_20X"))
        q["pct_target_30x"] = num(hs.get("PCT_TARGET_BASES_30X"))
        q["fold_80_base_penalty"] = num(hs.get("FOLD_80_BASE_PENALTY"))  # uniformity

    q.update(parse_mosdepth_summary(a.mosdepth_summary))
    per_gene = parse_mosdepth_regions(a.mosdepth_regions)
    if per_gene:
        q["per_gene_depth"] = per_gene
    # Fall back to mosdepth region mean if HsMetrics absent.
    if q.get("mean_target_depth") is None and q.get("mean_region_depth") is not None:
        q["mean_target_depth"] = q["mean_region_depth"]

    q.update(parse_verifybamid(a.verifybamid))

    # Sex inference from X/Y depth is only possible if the panel has X/Y targets;
    # per-gene depths here are autosomal in many panels -> mark NA (gate handles it).
    xy = [g for g in per_gene if g.upper() in ("SRY", "ZFY")]
    q["sex_inferred"] = "NA"
    q["sex_assessed"] = False

    with open(a.out_json, "w") as fh:
        json.dump(q, fh, indent=2)
    # tiny insert-size histogram sidecar for the plot (downsampled)
    if hist:
        with open(a.out_json.replace(".json", ".ishist.tsv"), "w") as fh:
            fh.write("insert_size\tcount\n")
            for isz, c in hist:
                fh.write(f"{isz}\t{c}\n")
    print(f"[summarize_sample_qc] {a.sample_id}: depth={q.get('mean_target_depth')} "
          f"on_target={q.get('pct_on_target')} dup={q.get('dup_rate')} "
          f"mapped%={q.get('pct_mapped')}")


if __name__ == "__main__":
    main()

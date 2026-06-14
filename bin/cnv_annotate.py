#!/usr/bin/env python3
"""cnv_annotate.py — Stage 8 CNV call annotation (CallForge).

Maps CNVkit per-sample segment calls (.call.cns) and per-target copy ratios
(.cnr) to panel genes (via the GENE|class target BED) and LABELS each call with
its CaptureForge CNV-callability (from gene_metadata.tsv): depth-callable genes
(e.g. GYPC/HP) are 'callable'; breakpoint-blind genes (e.g. CR1/CFH) are
'breakpoint_blind/low_confidence'. The label is PROPAGATED onto every emitted
call so a CNV at CR1 visibly carries 'breakpoint_blind' — not computed then dropped.

Inputs:
  --cns        space/tab CNVkit .call.cns files (one per sample; name = <sample>.call.cns)
  --cnr        CNVkit .cnr files (per-target log2/copy-ratio; name = <sample>.cnr)
  --bed        target BED (GENE|class names) for region->gene mapping
  --metadata   gene_metadata.tsv (cnv_callable column)
Outputs:
  cnv_calls.tsv        one row per non-neutral CNV call, with gene + callability label
  cnv_copyratio.tsv    per-target copy ratio (long: sample, gene, chrom, start, end, log2, depth)
  cnv_summary.json     calls/sample, calls by callability, genes flagged low-confidence
"""
import argparse, csv, json, os
from collections import defaultdict


def load_bed(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            rows.append((f[0], int(f[1]), int(f[2]), (f[3] if len(f) > 3 else "NA|coding").split("|")[0]))
    return rows


def gene_for(bed, chrom, start, end):
    best, ov = None, 0
    for (c, s, e, g) in bed:
        if c != chrom:
            continue
        o = min(end, e) - max(start, s)
        if o > ov:
            ov, best = o, g
    return best


def load_callability(path):
    lab = {}
    if path and os.path.isfile(path):
        with open(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                lab[r["gene"]] = r.get("cnv_callable", "unknown")
    return lab


def confidence(callable_flag):
    return {"yes": "callable", "no": "breakpoint_blind_low_confidence",
            "unknown": "unknown_confidence", "n/a": "not_cnv_target"}.get(callable_flag, "unknown_confidence")


def read_cns(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cns", nargs="+", required=True)
    ap.add_argument("--cnr", nargs="*", default=[])
    ap.add_argument("--bed", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()

    bed = load_bed(a.bed)
    callab = load_callability(a.metadata)
    os.makedirs(a.outdir, exist_ok=True)

    calls = []
    per_sample = defaultdict(int)
    by_conf = defaultdict(int)
    for path in a.cns:
        sample = os.path.basename(path).split(".")[0]
        for r in read_cns(path):
            chrom = r.get("chromosome") or r.get("chrom")
            start, end = int(float(r.get("start", 0))), int(float(r.get("end", 0)))
            cn = r.get("cn")
            cn = int(cn) if (cn not in (None, "", ".")) else None
            log2 = float(r.get("log2", 0))
            # non-neutral call: integer CN != 2 (diploid) when available, else |log2|>0.3
            is_cnv = (cn is not None and cn != 2) or (cn is None and abs(log2) > 0.3)
            if not is_cnv:
                continue
            gene = gene_for(bed, chrom, start, end) or "intergenic"
            cflag = callab.get(gene, "unknown")
            conf = confidence(cflag)
            calls.append({"sample": sample, "gene": gene, "chrom": chrom, "start": start,
                          "end": end, "log2": round(log2, 4), "cn": cn,
                          "type": "DEL" if (cn is not None and cn < 2) or log2 < 0 else "DUP",
                          "cnv_callable": cflag, "confidence": conf})
            per_sample[sample] += 1
            by_conf[conf] += 1

    cols = ["sample", "gene", "chrom", "start", "end", "type", "cn", "log2", "cnv_callable", "confidence"]
    with open(os.path.join(a.outdir, "cnv_calls.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t"); w.writeheader()
        for c in calls:
            w.writerow(c)

    # per-target copy ratios (for the per-target copy-ratio plot)
    with open(os.path.join(a.outdir, "cnv_copyratio.tsv"), "w") as fh:
        fh.write("sample\tgene\tchrom\tstart\tend\tlog2\tdepth\n")
        for path in a.cnr:
            sample = os.path.basename(path).split(".")[0]
            with open(path) as cfh:
                for r in csv.DictReader(cfh, delimiter="\t"):
                    chrom = r.get("chromosome") or r.get("chrom")
                    s, e = int(float(r.get("start", 0))), int(float(r.get("end", 0)))
                    gene = (r.get("gene") or gene_for(bed, chrom, s, e) or "intergenic")
                    gene = gene.split(",")[0].split("|")[0]
                    fh.write(f"{sample}\t{gene}\t{chrom}\t{s}\t{e}\t{r.get('log2','0')}\t{r.get('depth','')}\n")

    # genes flagged low-confidence among CNV targets (propagated label visibility)
    low_conf_genes = sorted(g for g, f in callab.items() if f == "no")
    summary = {"n_calls": len(calls), "calls_per_sample": dict(per_sample),
               "calls_by_confidence": dict(by_conf),
               "cnv_target_genes": sorted(g for g, f in callab.items() if f in ("yes", "no")),
               "breakpoint_blind_genes": low_conf_genes,
               "callable_genes": sorted(g for g, f in callab.items() if f == "yes")}
    with open(os.path.join(a.outdir, "cnv_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"[cnv_annotate] {len(calls)} CNV calls; by confidence={dict(by_conf)}; "
          f"breakpoint-blind genes={low_conf_genes}; callable={summary['callable_genes']}")


if __name__ == "__main__":
    main()

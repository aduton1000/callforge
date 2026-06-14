#!/usr/bin/env python3
"""annotate_phylop.py — add INFO/PhyloP from a conservation bigWig (CallForge Stage 11).

pyBigWig is queried at each variant position; the bigWig's contig naming
(usually chr-prefixed) is reconciled to the query VCF's naming so scores are not
silently missed on a chr/no-chr mismatch. Scores are written back as INFO/PhyloP
via bcftools annotate (annotation keyed on the query VCF's own CHROM).
"""
import argparse, math, os, subprocess, sys
import pyBigWig


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--bigwig", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    bw = pyBigWig.open(a.bigwig)
    bw_chroms = set(bw.chroms().keys())

    def bw_contig(c):
        for cand in (c, "chr" + c, c[3:] if c.lower().startswith("chr") else c,
                     "chrM" if c == "MT" else c):
            if cand in bw_chroms:
                return cand
        return None

    q = subprocess.run(["bcftools", "query", "-f", "%CHROM\t%POS\n", a.vcf],
                       check=True, capture_output=True, text=True).stdout
    annot = a.out + ".phylop.tsv"
    n_written = 0
    with open(annot, "w") as fh:
        for line in q.splitlines():
            chrom, pos = line.split("\t")
            pos = int(pos)
            bc = bw_contig(chrom)
            if not bc:
                continue
            try:
                vals = bw.values(bc, pos - 1, pos)
                v = vals[0] if vals else None
            except (RuntimeError, IndexError):
                v = None
            if v is None or (isinstance(v, float) and math.isnan(v)):
                continue
            fh.write(f"{chrom}\t{pos}\t{round(float(v), 4)}\n")
            n_written += 1
    bw.close()

    if n_written == 0:
        sys.stderr.write("[annotate_phylop] no PhyloP scores matched (contig mismatch?); passing through\n")
        subprocess.run(["bcftools", "view", a.vcf, "-Oz", "-o", a.out], check=True)
        subprocess.run(["tabix", "-f", "-p", "vcf", a.out], check=True)
        return

    subprocess.run(["bgzip", "-f", annot], check=True)
    subprocess.run(["tabix", "-f", "-s", "1", "-b", "2", "-e", "2", annot + ".gz"], check=True)
    hdr = a.out + ".phylop.hdr"
    with open(hdr, "w") as fh:
        fh.write('##INFO=<ID=PhyloP,Number=1,Type=Float,Description="PhyloP conservation score (bigWig)">\n')
    subprocess.run(["bcftools", "annotate", "-a", annot + ".gz", "-h", hdr,
                    "-c", "CHROM,POS,PhyloP", a.vcf, "-Oz", "-o", a.out], check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", a.out], check=True)
    print(f"[annotate_phylop] PhyloP written for {n_written} variants")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""check_reference_invariant.py — Stage 0 reference invariant (CallForge).

Enforces, fail-loud, the invariant inherited from CaptureForge:
  1. The reference is a NO-ALT primary assembly — no alt/decoy/HLA contigs.
  2. The target BED's build matches the reference: every BED contig exists in
     the reference .fai, chr-prefix style is consistent, and no interval runs
     off the end of its contig.

Reads the reference .fai (call `samtools faidx` first if absent — the module
does). stdlib-only.

Exit 0 = invariant holds; exit 2 = violated (with a precise reason).
"""
import argparse, json, os, re, sys

# Contigs that must NOT appear in a no-alt primary assembly.
ALT_DECOY_PAT = re.compile(r"(_alt\b|_alt$|_decoy|hs38d1|chrEBV|\bEBV\b|^HLA-|_hap\d)", re.I)


def fail(report, path, reason):
    report["pass"] = False
    report["reason"] = reason
    if path:
        with open(path, "w") as fh:
            json.dump(report, fh, indent=2)
    sys.stderr.write(f"[reference_invariant] FAIL: {reason}\n")
    sys.exit(2)


def load_fai(fai):
    contigs = {}
    with open(fai) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                contigs[f[0]] = int(f[1])
    return contigs


def chr_style(names):
    # 'ucsc' if (almost) all contigs are chr-prefixed; 'ensembl' otherwise.
    chrpref = sum(1 for n in names if n.startswith("chr"))
    return "ucsc" if chrpref > len(names) / 2 else "ensembl"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fai", required=True, help="reference .fai")
    ap.add_argument("--bed", required=True, help="target BED")
    ap.add_argument("--genome-build", default="GRCh38")
    ap.add_argument("--out-json", default=None)
    a = ap.parse_args()

    report = {"pass": True, "genome_build": a.genome_build, "fai": a.fai, "bed": a.bed}

    if not os.path.isfile(a.fai):
        fail(report, a.out_json, f"reference .fai not found: {a.fai} (run samtools faidx)")
    if not os.path.isfile(a.bed):
        fail(report, a.out_json, f"target BED not found: {a.bed}")

    contigs = load_fai(a.fai)
    if not contigs:
        fail(report, a.out_json, f"empty/unreadable .fai: {a.fai}")

    # (1) no-alt assertion
    alts = [c for c in contigs if ALT_DECOY_PAT.search(c)]
    report["n_contigs"] = len(contigs)
    report["ref_chr_style"] = chr_style(list(contigs))
    if alts:
        fail(report, a.out_json,
             f"reference is NOT no-alt: {len(alts)} alt/decoy contigs present "
             f"(e.g. {alts[:5]}). CallForge requires the same no-alt primary "
             f"assembly CaptureForge designed against.")

    # (2) BED-vs-reference build match
    bed_contigs, max_end = {}, {}
    n_intervals = 0
    with open(a.bed) as fh:
        for ln, line in enumerate(fh, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                fail(report, a.out_json, f"BED line {ln}: fewer than 3 columns")
            c, start, end = f[0], int(f[1]), int(f[2])
            bed_contigs[c] = bed_contigs.get(c, 0) + 1
            max_end[c] = max(max_end.get(c, 0), end)
            n_intervals += 1
    report["n_bed_intervals"] = n_intervals
    report["bed_chr_style"] = chr_style(list(bed_contigs))

    if report["bed_chr_style"] != report["ref_chr_style"]:
        fail(report, a.out_json,
             f"chromosome-naming mismatch: BED is '{report['bed_chr_style']}' "
             f"but reference is '{report['ref_chr_style']}'. Re-export the BED "
             f"or reference with matching contig names (e.g. '1' vs 'chr1').")

    missing = sorted(c for c in bed_contigs if c not in contigs)
    if missing:
        fail(report, a.out_json,
             f"{len(missing)} BED contig(s) absent from the reference: {missing[:8]}")

    overruns = [(c, max_end[c], contigs[c]) for c in bed_contigs if max_end[c] > contigs[c]]
    if overruns:
        c, e, L = overruns[0]
        fail(report, a.out_json,
             f"BED interval runs past contig end on '{c}' ({e} > {L}); "
             f"BED and reference are different builds.")

    report["bed_contigs_used"] = len(bed_contigs)
    report["reason"] = "no-alt reference confirmed; BED build matches reference"
    if a.out_json:
        with open(a.out_json, "w") as fh:
            json.dump(report, fh, indent=2)
    print(f"[reference_invariant] OK: no-alt {a.genome_build} ({len(contigs)} contigs, "
          f"{report['ref_chr_style']} naming); BED uses {len(bed_contigs)} contigs / "
          f"{n_intervals} intervals, build matches.")


if __name__ == "__main__":
    main()

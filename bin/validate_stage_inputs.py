#!/usr/bin/env python3
"""validate_stage_inputs.py — fail-loud input validation for CallForge stage subcommands.

A thin, stdlib-only gate that the per-stage entry workflows run BEFORE invoking the
real stage modules. It does NOT reimplement stage logic — it only asserts the
externally-provided inputs are well-formed and mutually consistent, exiting non-zero
(2) on any mismatch so a standalone stage run fails as loudly as the full pipeline.

Checks (opt-in via flags):
  --vcf FILE ...          bgzip magic + companion .tbi present (index); contig names
                          (reconciled for an optional 'chr' prefix) are a subset of the
                          reference contigs when --fai is given.
  --bam FILE ...          companion .bai/.csi present; @SQ contigs subset of --fai
                          (uses samtools if on PATH; otherwise index-only check).
  --fai FILE              reference .fai supplying the allowed contig set + build check.
  --samplesheet FILE      with --require-phenotype: assert >=1 non-empty phenotype value
                          (burden needs it; the pipeline would otherwise write a skip
                          marker — standalone burden stops loudly instead).

Exit 0 = all requested checks pass; exit 2 = a hard mismatch (message on stderr).
"""
import argparse, csv, gzip, os, re, shutil, subprocess, sys

def die(msg):
    sys.stderr.write(f"[validate_stage_inputs] ERROR: {msg}\n")
    sys.exit(2)

def ok(msg):
    sys.stderr.write(f"[validate_stage_inputs] OK: {msg}\n")

def norm_contig(c):
    return c[3:] if c.lower().startswith("chr") else c

def fai_contigs(fai):
    out = set()
    with open(fai) as fh:
        for line in fh:
            if line.strip():
                out.add(norm_contig(line.split("\t")[0]))
    return out

def is_bgzip(path):
    with open(path, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"

def vcf_header_contigs(path):
    contigs = []
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", errors="replace") as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            m = re.match(r"##contig=<ID=([^,>]+)", line)
            if m:
                contigs.append(m.group(1))
    return contigs

def check_vcf(path, fai_set):
    if not os.path.isfile(path):
        die(f"VCF not found: {path}")
    if path.endswith(".gz") and not is_bgzip(path):
        die(f"VCF '{path}' is .gz but not bgzip-compressed (gzip magic absent)")
    tbi = path + ".tbi"
    csi = path + ".csi"
    if not (os.path.isfile(tbi) or os.path.isfile(csi)):
        die(f"VCF index missing for {path} (expected {os.path.basename(tbi)} or .csi) — "
            f"run `tabix -p vcf {os.path.basename(path)}`")
    if fai_set:
        vc = vcf_header_contigs(path)
        if vc:
            stray = sorted({c for c in vc if norm_contig(c) not in fai_set})
            if stray:
                die(f"VCF '{path}' has contig(s) absent from the reference "
                    f".fai (build/contig mismatch): {', '.join(stray[:8])}"
                    + (" …" if len(stray) > 8 else ""))
    ok(f"VCF {os.path.basename(path)} (bgzip + index"
       + (" + contigs" if fai_set else "") + ")")

def bam_header_contigs(path):
    if not shutil.which("samtools"):
        return None
    try:
        out = subprocess.run(["samtools", "view", "-H", path], capture_output=True,
                             text=True, check=True).stdout
    except subprocess.CalledProcessError as e:
        die(f"samtools could not read BAM header for {path}: {e}")
    return [l.split("\tSN:")[1].split("\t")[0] for l in out.splitlines()
            if l.startswith("@SQ") and "\tSN:" in l]

def check_bam(path, fai_set):
    if not os.path.isfile(path):
        die(f"BAM not found: {path}")
    if not (os.path.isfile(path + ".bai") or os.path.isfile(re.sub(r"\.bam$", ".bai", path))
            or os.path.isfile(path + ".csi")):
        die(f"BAM index missing for {path} (expected .bai/.csi) — run `samtools index {os.path.basename(path)}`")
    contigs = bam_header_contigs(path)
    if fai_set and contigs:
        stray = sorted({c for c in contigs if norm_contig(c) not in fai_set})
        if stray:
            die(f"BAM '{path}' has @SQ contig(s) absent from the reference .fai "
                f"(build/contig mismatch): {', '.join(stray[:8])}")
    ok(f"BAM {os.path.basename(path)} (index"
       + (" + contigs" if (fai_set and contigs) else "") + ")")

def check_phenotype(sheet):
    if not os.path.isfile(sheet):
        die(f"sample sheet not found: {sheet}")
    with open(sheet, newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows or "phenotype" not in rows[0]:
        die(f"sample sheet '{sheet}' has no 'phenotype' column — burden needs a "
            f"phenotype (case/control or quantitative). Add it (see §5.1) and re-run.")
    n = sum(1 for r in rows if (r.get("phenotype") or "").strip())
    if n == 0:
        die(f"sample sheet '{sheet}' has a phenotype column but every value is blank — "
            f"burden cannot run. Fill phenotype for the cohort and re-run.")
    ok(f"phenotype present for {n}/{len(rows)} sample(s)")

def main(argv=None):
    ap = argparse.ArgumentParser(prog="validate_stage_inputs.py")
    ap.add_argument("--vcf", nargs="*", default=[])
    ap.add_argument("--bam", nargs="*", default=[])
    ap.add_argument("--fai", default=None)
    ap.add_argument("--samplesheet", default=None)
    ap.add_argument("--require-phenotype", action="store_true")
    ap.add_argument("--stage", default="stage")
    a = ap.parse_args(argv)

    fai_set = fai_contigs(a.fai) if a.fai else set()
    if a.fai and not fai_set:
        die(f"reference .fai is empty: {a.fai}")
    for v in a.vcf:
        check_vcf(v, fai_set)
    for b in a.bam:
        check_bam(b, fai_set)
    if a.require_phenotype:
        if not a.samplesheet:
            die("--require-phenotype needs --samplesheet")
        check_phenotype(a.samplesheet)
    print(f"[validate_stage_inputs] {a.stage}: all input checks passed")
    return 0

if __name__ == "__main__":
    sys.exit(main())

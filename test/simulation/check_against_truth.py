#!/usr/bin/env python3
"""check_against_truth.py — assert CallForge results against the planted truth manifest.

Given truth_manifest.json (from build_fixture.py) and a CallForge results/ directory,
checks each stage for CORRECTNESS (did it recover what we planted?) and prints a
PASS/FAIL line per stage plus an overall verdict. Thresholds are NAMED constants carried
in the manifest (see cohort_design.THRESHOLDS) — printed so they are visible and tunable.

Tolerances are deliberately loose where the stage is probabilistic (calling/CNV/STR): we
assert the planted SIGNAL surfaces, not exact metrics. See README "Checker thresholds".

Result files consumed (relative to --results):
  stage7_filter/joint.filtered.vcf.gz      calling
  stage11_annotation/variants.flat.tsv     annotation consequence
  stage9_str/str_calls.tsv                 STR expansion
  stage8_cnv/cnv_calls.tsv                 CNV del/dup
  stage14_burden/burden_results.tsv        burden enrichment
  stage12_cohortqc/cohort_qc.json          sex concordance + relatedness

stdlib-only. Exit code 0 if every (non-skipped) stage passes, else 1.
"""
import argparse
import csv
import gzip
import json
import math
import os
import sys

PASS, FAIL, SKIP = "PASS", "FAIL", "SKIP"


# --------------------------------------------------------------- minimal VCF reader
def read_vcf(path):
    """Yield (chrom, pos, ref, [alts], samples, {sample: gt_string}) for a small VCF.gz."""
    op = gzip.open if path.endswith(".gz") else open
    samples = []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            chrom, pos, ref, alt = f[0], int(f[1]), f[3], f[4].split(",")
            fmt = f[8].split(":")
            gi = fmt.index("GT") if "GT" in fmt else 0
            gts = {}
            for s, col in zip(samples, f[9:]):
                parts = col.split(":")
                gts[s] = parts[gi] if gi < len(parts) else "./."
            yield chrom, pos, ref, alt, samples, gts


def gt_nonref(gt):
    for a in gt.replace("|", "/").split("/"):
        if a not in ("0", ".", ""):
            return True
    return False


def load_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def frac(n, d):
    return (n / d) if d else 0.0


# ------------------------------------------------------------------------- checks
def check_calling(man, results):
    vcf = os.path.join(results, "stage7_filter", "joint.filtered.vcf.gz")
    if not os.path.isfile(vcf):
        return SKIP, [f"calling VCF not found: {vcf}"]
    # index VCF records by (chrom,pos)
    recs = {}
    for chrom, pos, ref, alts, samples, gts in read_vcf(vcf):
        recs.setdefault((chrom, pos), []).append((ref, alts, gts))
    thr = man["thresholds"]["CALLING_MIN_RECALL"]
    total = recovered = 0
    misses = []
    for t in man["coding_variants"]:
        total += 1
        carriers = t["carriers"]
        need = math.ceil(0.5 * len(carriers))
        ok = False
        if t["vartype"] == "snv":
            for (ref, alts, gts) in recs.get((t["contig"], t["pos"]), []):
                if t["alt"] in alts and sum(gt_nonref(gts.get(c, "./.")) for c in carriers) >= need:
                    ok = True
                    break
        else:  # indel (frameshift) — match any indel call near the locus in carriers
            for (c2, p2) in recs:
                if c2 == t["contig"] and abs(p2 - t["pos"]) <= 30:
                    for (ref, alts, gts) in recs[(c2, p2)]:
                        if any(len(a) != len(ref) for a in alts) and \
                           sum(gt_nonref(gts.get(c, "./.")) for c in carriers) >= need:
                            ok = True
                            break
                if ok:
                    break
        recovered += ok
        if not ok:
            misses.append(t["id"])
    rec = frac(recovered, total)
    status = PASS if rec >= thr else FAIL
    lines = [f"planted coding variants recovered (genotype-aware): {recovered}/{total} "
             f"= {rec:.2f} (threshold >= {thr})"]
    if misses:
        lines.append(f"  missed: {', '.join(misses)}")
    return status, lines


def check_annotation(man, results):
    flat = os.path.join(results, "stage11_annotation", "variants.flat.tsv")
    if not os.path.isfile(flat):
        return SKIP, [f"variants.flat.tsv not found: {flat}"]
    rows = load_tsv(flat)
    by_pos = {(r.get("chrom"), r.get("pos")): r for r in rows}
    thr = man["thresholds"]["ANNOTATION_MIN_CSQ_MATCH"]
    checked = matched = 0
    proms_ok = proms_total = 0
    bad = []
    for t in man["coding_variants"] + man["promoter_variants"]:
        r = by_pos.get((t["contig"], str(t["pos"])))
        is_prom = t in man["promoter_variants"]
        if is_prom:
            proms_total += 1
        else:
            checked += 1
        if not r:
            bad.append(f"{t['id']}(not annotated)")
            continue
        csq = (r.get("consequence") or "").lower()
        want = t["consequence"]
        if is_prom:
            if "upstream" in csq or "regulatory" in csq:
                proms_ok += 1
            else:
                bad.append(f"{t['id']}({csq}!=upstream)")
        else:
            if want in csq:
                matched += 1
            else:
                bad.append(f"{t['id']}({csq}!={want})")
    m = frac(matched, checked)
    pm = frac(proms_ok, proms_total)
    status = PASS if (m >= thr and pm >= thr) else FAIL
    lines = [f"coding consequence exact match: {matched}/{checked} = {m:.2f} (>= {thr})",
             f"promoter -> upstream/regulatory: {proms_ok}/{proms_total} = {pm:.2f} (>= {thr})"]
    if bad:
        lines.append("  mismatches: " + ", ".join(bad[:8]) + (" …" if len(bad) > 8 else ""))
    return status, lines


def check_str(man, results):
    s = man.get("str")
    path = os.path.join(results, "stage9_str", "str_calls.tsv")
    if not s:
        return SKIP, ["no STR locus in manifest"]
    if not os.path.isfile(path):
        return SKIP, [f"str_calls.tsv not found: {path}"]
    rows = load_tsv(path)
    delta = man["thresholds"]["STR_EXPANSION_MIN_DELTA_UNITS"]
    tol = man["thresholds"]["STR_LEN_TOL_UNITS"]
    rec_thr = man["thresholds"]["STR_MIN_CARRIER_RECALL"]
    carriers = set(s["carriers"])
    exp, norm = s["expanded_units"], s["normal_units"]

    def alleles(sample):
        vals = []
        for r in rows:
            if r.get("sample") == sample and r.get("called") == "yes":
                for k in ("allele1", "allele2"):
                    v = r.get(k, "")
                    if str(v).strip() not in ("", "."):
                        vals.append(int(float(v)))
        return vals

    carrier_hits = 0
    detail = []
    for c in sorted(carriers):
        al = alleles(c)
        big = max(al) if al else 0
        expanded = bool(al) and (big - min(al) >= delta or big >= norm + delta) and abs(big - exp) <= tol
        carrier_hits += expanded
        if not expanded:
            detail.append(f"{c}:alleles={al}")
    # spot-check non-carriers stay near normal (informational, not a hard gate)
    noncarrier_expanded = []
    for samp in [x["sample_id"] for x in man["samples"] if x["sample_id"] not in carriers]:
        al = alleles(samp)
        if al and max(al) >= norm + delta:
            noncarrier_expanded.append(samp)
    rec = frac(carrier_hits, len(carriers))
    status = PASS if rec >= rec_thr else FAIL
    lines = [f"STR carriers flagged expanded (~{exp}±{tol} units): {carrier_hits}/{len(carriers)} "
             f"= {rec:.2f} (>= {rec_thr})"]
    if detail:
        lines.append("  carrier misses: " + ", ".join(detail))
    if noncarrier_expanded:
        lines.append(f"  WARN spurious expansion in non-carriers: {', '.join(noncarrier_expanded)}")
    return status, lines


def check_cnv(man, results):
    path = os.path.join(results, "stage8_cnv", "cnv_calls.tsv")
    if not os.path.isfile(path):
        return SKIP, [f"cnv_calls.tsv not found: {path}"]
    rows = load_tsv(path)
    thr = man["thresholds"]["CNV_MIN_RECALL"]
    total = recovered = 0
    misses = []
    for e in man["cnv_events"]:
        for c in e["carriers"]:
            total += 1
            hit = any(r.get("sample") == c and r.get("gene") == e["gene"]
                      and (r.get("type") or "").upper() == e["state"] for r in rows)
            recovered += hit
            if not hit:
                misses.append(f"{e['gene']}/{c}/{e['state']}")
    rec = frac(recovered, total)
    status = PASS if rec >= thr else FAIL
    lines = [f"planted CNV events called (sample+gene+type): {recovered}/{total} "
             f"= {rec:.2f} (>= {thr})"]
    if misses:
        lines.append("  missed: " + ", ".join(misses[:10]) + (" …" if len(misses) > 10 else ""))
    return status, lines


def check_burden(man, results):
    path = os.path.join(results, "stage14_burden", "burden_results.tsv")
    if not os.path.isfile(path):
        return SKIP, [f"burden_results.tsv not found: {path}"]
    rows = load_tsv(path)
    gene = man["burden"]["gene"]
    pmax = man["thresholds"]["BURDEN_P_MAX"]

    def pval(r):
        try:
            return float(r.get("p", "nan"))
        except ValueError:
            return float("nan")

    units = [(r.get("unit", ""), pval(r)) for r in rows if not math.isnan(pval(r))]
    if not units:
        return FAIL, ["no usable p-values in burden_results.tsv"]
    target = [(u, p) for (u, p) in units if u.endswith(gene) or u == gene or gene in u]
    if not target:
        return FAIL, [f"burden gene {gene} not present in results"]
    best_target_p = min(p for _, p in target)
    overall_best = min(p for _, p in units)
    is_top = abs(best_target_p - overall_best) < 1e-12
    status = PASS if (is_top or best_target_p <= pmax) else FAIL
    lines = [f"burden gene {gene}: p={best_target_p:.3g} "
             f"(rank#1={is_top}, threshold p<= {pmax}); cohort best p={overall_best:.3g}"]
    return status, lines


def check_cohortqc(man, results):
    path = os.path.join(results, "stage12_cohortqc", "cohort_qc.json")
    if not os.path.isfile(path):
        return SKIP, [f"cohort_qc.json not found: {path}"]
    with open(path) as fh:
        qc = json.load(fh)
    declared = man["cohortqc"]["sex_by_sample"]
    inferred = qc.get("sex_inferred", {})
    conc = sum(1 for s, sx in declared.items() if inferred.get(s) == sx)
    total = len(declared)
    sex_frac = frac(conc, total)
    sex_thr = man["thresholds"]["COHORTQC_SEX_MIN_CONCORDANCE"]
    rel = qc.get("relatedness", {})
    rel_thr = man["thresholds"]["COHORTQC_MAX_RELATEDNESS"]
    high = {k: v for k, v in rel.items() if isinstance(v, (int, float)) and v > rel_thr}
    status = PASS if (sex_frac >= sex_thr and not high) else FAIL
    lines = [f"somalier sex concordance: {conc}/{total} = {sex_frac:.2f} (>= {sex_thr})",
             f"max pairwise relatedness <= {rel_thr}: "
             f"{'ok' if not high else 'VIOLATED ' + str(high)}"]
    return status, lines


CHECKS = [
    ("calling (stage 6/7)", check_calling),
    ("annotation (stage 11)", check_annotation),
    ("STR (stage 9)", check_str),
    ("CNV (stage 8)", check_cnv),
    ("burden (stage 14)", check_burden),
    ("cohort QC (stage 12)", check_cohortqc),
]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True, help="truth_manifest.json from build_fixture.py")
    ap.add_argument("--results", required=True, help="CallForge results/ directory")
    ap.add_argument("--skip-missing", action="store_true",
                    help="treat absent stage outputs as SKIP (default) rather than FAIL")
    a = ap.parse_args()

    with open(a.manifest) as fh:
        man = json.load(fh)

    print(f"== CallForge simulation truth check ==")
    print(f"manifest: {a.manifest}")
    print(f"results:  {a.results}\n")

    summary = []
    any_fail = False
    for name, fn in CHECKS:
        try:
            status, lines = fn(man, a.results)
        except Exception as e:                      # never let one stage abort the report
            status, lines = FAIL, [f"checker error: {e}"]
        if status == SKIP and not a.skip_missing:
            status = FAIL
        summary.append((name, status))
        if status == FAIL:
            any_fail = True
        print(f"[{status}] {name}")
        for ln in lines:
            print(f"        {ln}")
    print("\n== summary ==")
    for name, status in summary:
        print(f"  {status:4s}  {name}")
    overall = "FAIL" if any_fail else "PASS"
    print(f"\nOVERALL: {overall}")
    return 1 if any_fail else 0


if __name__ == "__main__":
    sys.exit(main())

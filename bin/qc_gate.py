#!/usr/bin/env python3
"""qc_gate.py — Stage 5 per-sample QC gate & quarantine (CallForge).

Aggregates the per-sample QC JSONs (summarize_sample_qc.py) and applies the
documented thresholds. Failing samples are FLAGGED and QUARANTINED (excluded from
joint calling) — never silently dropped; they remain in every report with the
reason. Metrics that could not be assessed (e.g. contamination without a
VerifyBamID2 panel, sex without X/Y targets) are recorded NA and do not fail the
sample.

Emits:
  qc_scorecard.tsv        one row per sample, per-metric pass/fail + decision
  qc_pass.txt             sample_ids that enter joint calling
  qc_quarantine.txt       sample_ids quarantined (+ reasons)
  cohort_qc_summary.json  thresholds, counts, per-sample decisions
"""
import argparse, json, glob, os, sys


def load(paths):
    out = []
    for p in paths:
        with open(p) as fh:
            out.append(json.load(fh))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qc-json", nargs="+", required=True, help="per-sample QC JSONs")
    ap.add_argument("--min-mean-target-depth", type=float, default=30)
    ap.add_argument("--max-dup-rate", type=float, default=0.40)
    ap.add_argument("--min-on-target", type=float, default=0.40)
    ap.add_argument("--max-contamination", type=float, default=0.03)
    ap.add_argument("--enforce-sex-check", action="store_true")
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()

    samples = load(a.qc_json)
    os.makedirs(a.outdir, exist_ok=True)
    thr = {"min_mean_target_depth": a.min_mean_target_depth, "max_dup_rate": a.max_dup_rate,
           "min_on_target": a.min_on_target, "max_contamination": a.max_contamination,
           "enforce_sex_check": a.enforce_sex_check}

    rows, passed, quarantined = [], [], []
    for q in sorted(samples, key=lambda x: x["sample_id"]):
        sid = q["sample_id"]
        reasons, checks = [], {}

        def check(name, value, ok, fail_msg):
            if value is None:
                checks[name] = "NA"
                return
            checks[name] = "pass" if ok else "FAIL"
            if not ok:
                reasons.append(fail_msg)

        depth = q.get("mean_target_depth")
        check("depth", depth, depth is not None and depth >= a.min_mean_target_depth,
              f"mean_target_depth {depth} < {a.min_mean_target_depth}")
        dup = q.get("dup_rate")
        check("dup", dup, dup is not None and dup <= a.max_dup_rate,
              f"dup_rate {dup} > {a.max_dup_rate}")
        ont = q.get("pct_on_target")
        check("on_target", ont, ont is not None and ont >= a.min_on_target,
              f"pct_on_target {ont} < {a.min_on_target}")

        if q.get("contamination_assessed"):
            con = q.get("contamination")
            check("contamination", con, con is not None and con <= a.max_contamination,
                  f"contamination {con} > {a.max_contamination}")
        else:
            checks["contamination"] = "NA"

        if q.get("sex_assessed"):
            concord = q.get("declared_sex", "U") == q.get("sex_inferred")
            if a.enforce_sex_check:
                check("sex", concord, concord,
                      f"sex mismatch declared={q.get('declared_sex')} inferred={q.get('sex_inferred')}")
            else:
                checks["sex"] = "pass" if concord else "warn"
        else:
            checks["sex"] = "NA"

        decision = "QUARANTINE" if reasons else "PASS"
        (quarantined if reasons else passed).append(sid)
        rows.append({"sample_id": sid, "decision": decision, "checks": checks,
                     "reasons": "; ".join(reasons),
                     "mean_target_depth": depth, "dup_rate": dup, "pct_on_target": ont,
                     "fold_enrichment": q.get("fold_enrichment"),
                     "pct_target_20x": q.get("pct_target_20x"),
                     "contamination": q.get("contamination"),
                     "fold_80_base_penalty": q.get("fold_80_base_penalty")})

    # scorecard TSV
    cols = ["sample_id", "decision", "mean_target_depth", "dup_rate", "pct_on_target",
            "fold_enrichment", "pct_target_20x", "fold_80_base_penalty", "contamination",
            "check_depth", "check_dup", "check_on_target", "check_contamination", "check_sex",
            "reasons"]
    with open(os.path.join(a.outdir, "qc_scorecard.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            c = r["checks"]
            fh.write("\t".join(str(x) for x in [
                r["sample_id"], r["decision"], r["mean_target_depth"], r["dup_rate"],
                r["pct_on_target"], r["fold_enrichment"], r["pct_target_20x"],
                r["fold_80_base_penalty"], r["contamination"],
                c.get("depth"), c.get("dup"), c.get("on_target"),
                c.get("contamination"), c.get("sex"), r["reasons"]]) + "\n")

    with open(os.path.join(a.outdir, "qc_pass.txt"), "w") as fh:
        fh.write("\n".join(passed) + ("\n" if passed else ""))
    with open(os.path.join(a.outdir, "qc_quarantine.txt"), "w") as fh:
        for r in rows:
            if r["decision"] == "QUARANTINE":
                fh.write(f"{r['sample_id']}\t{r['reasons']}\n")

    summary = {"thresholds": thr, "n_samples": len(rows),
               "n_pass": len(passed), "n_quarantine": len(quarantined),
               "passed": passed, "quarantined": quarantined,
               "decisions": {r["sample_id"]: r["decision"] for r in rows}}
    with open(os.path.join(a.outdir, "cohort_qc_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    print(f"[qc_gate] {len(rows)} samples: {len(passed)} PASS, {len(quarantined)} QUARANTINE")
    if quarantined:
        print(f"[qc_gate] quarantined (kept + reported, excluded from joint calling): {quarantined}")
    if not passed:
        sys.stderr.write("[qc_gate] WARNING: no samples passed the QC gate — joint calling will be empty.\n")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""make_provenance.py — per-run provenance.json (CallForge, Stage 15).

Captures, for reproducibility: pipeline version + git commit, genome build, the
resource manifest (found/missing + build/md5), which samples were quarantined and
why, sample-sheet summary, burden ancestry-PC gate decision, and pointers to the
pinned env/*.yml for exact tool versions (each stage runs in its own env/container).
stdlib only.
"""
import argparse, json, os, subprocess


def jload(p):
    try:
        return json.load(open(p)) if p and os.path.isfile(p) else None
    except json.JSONDecodeError:
        return None


def git_commit():
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip() or None
    except OSError:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--version", default="0.1.0")
    ap.add_argument("--profile", default="")
    ap.add_argument("--nextflow-version", default="")
    ap.add_argument("--genome-build", default="GRCh38")
    ap.add_argument("--genome-fasta", default="")
    ap.add_argument("--manifest", default=None)
    ap.add_argument("--sheet-summary", default=None)
    ap.add_argument("--qc-summary", default=None)      # cohort_qc_summary.json (Stage 5 gate)
    ap.add_argument("--burden-meta", default=None)
    ap.add_argument("--quarantine", default=None)      # qc_quarantine.txt
    ap.add_argument("--out", default="provenance.json")
    a = ap.parse_args()

    manifest = jload(a.manifest) or {}
    res = {}
    for name, r in (manifest.get("resources") or {}).items():
        if isinstance(r, dict):
            res[name] = {k: r.get(k) for k in ("path", "status", "build", "release",
                                               "md5", "md5:first1MiB", "has_afr") if k in r}

    quarantined = []
    if a.quarantine and os.path.isfile(a.quarantine):
        for line in open(a.quarantine):
            if line.strip():
                p = line.rstrip("\n").split("\t")
                quarantined.append({"sample": p[0], "reason": p[1] if len(p) > 1 else ""})

    qc = jload(a.qc_summary) or {}
    burden = jload(a.burden_meta) or {}
    prov = {
        "pipeline": "CallForge", "version": a.version,
        "git_commit": git_commit(), "profile": a.profile,
        "nextflow_version": a.nextflow_version,
        "reference": {"genome_build": a.genome_build, "genome_fasta": a.genome_fasta,
                      "invariant": "no-alt primary assembly == design genome; BED build asserted"},
        "resources": res,
        "samplesheet_summary": jload(a.sheet_summary),
        "qc_gate": {"thresholds": qc.get("thresholds"), "passed": qc.get("passed"),
                    "quarantined": qc.get("quarantined")},
        "quarantined_samples": quarantined,
        "burden": {"status": burden.get("status"), "engine_note": "engine stamped in burden_results.tsv",
                   "ancestry_pcs": burden.get("ancestry_pcs"),
                   "ancestry_pcs_used": burden.get("ancestry_pcs_used"),
                   "n_cases": burden.get("n_cases"), "n_controls": burden.get("n_controls")},
        "tool_versions_note": "exact pinned versions in env/*.yml; VEP/hap.py via pinned containers "
                              "(vep_image / happy_image); each stage runs in its own env/container.",
    }
    with open(a.out, "w") as fh:
        json.dump(prov, fh, indent=2)
    print(f"[make_provenance] commit={prov['git_commit']}, "
          f"quarantined={[q['sample'] for q in quarantined]} -> {a.out}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""validate_samplesheet.py — Stage 0 sample-sheet validation (CallForge).

Validates the CallForge sample sheet and emits (a) a normalized sheet and
(b) a JSON summary that downstream stages and the report consume.

Sample sheet (CSV) columns:
  sample_id   (required, unique)
  fastq_1     (required, exists)
  fastq_2     (required, exists; must differ from fastq_1)
  sex         (optional: M/F/male/female/1/2/unknown)
  phenotype   (optional: case/control 1/2/0 or quantitative; presence enables burden)
  covariate_* (optional: any number of numeric covariate columns)
  batch       (optional)

Fail-loud on: missing required columns, duplicate sample_id, missing FASTQ,
R1==R2. Warn (not fail) on: odd sex tokens, mixed phenotype types.

stdlib-only so it runs without the conda env active.
"""
import argparse, csv, json, os, sys, re

REQUIRED = ["sample_id", "fastq_1", "fastq_2"]
SEX_MAP = {"m": "M", "male": "M", "1": "M",
           "f": "F", "female": "F", "2": "F",
           "u": "U", "unknown": "U", "": "U", "0": "U"}


def err(msg):
    sys.stderr.write(f"[validate_samplesheet] ERROR: {msg}\n")
    sys.exit(1)


def warn(msg):
    sys.stderr.write(f"[validate_samplesheet] WARN: {msg}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sheet", required=True)
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--check-files", action="store_true",
                    help="assert FASTQ paths exist (skip for dry config parsing)")
    a = ap.parse_args()

    if not os.path.isfile(a.sheet):
        err(f"sample sheet not found: {a.sheet}")

    with open(a.sheet, newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        err("sample sheet is empty")

    cols = list(rows[0].keys())
    for c in REQUIRED:
        if c not in cols:
            err(f"missing required column '{c}' (found: {cols})")

    cov_cols = sorted(c for c in cols if c.startswith("covariate_"))
    has_phenotype = "phenotype" in cols and any((r.get("phenotype") or "").strip() for r in rows)
    has_sex = "sex" in cols
    has_batch = "batch" in cols

    seen, samples = set(), []
    pheno_vals = []
    for i, r in enumerate(rows, 1):
        sid = (r.get("sample_id") or "").strip()
        if not sid:
            err(f"row {i}: empty sample_id")
        if not re.fullmatch(r"[A-Za-z0-9._-]+", sid):
            err(f"row {i}: sample_id '{sid}' has illegal characters (use [A-Za-z0-9._-])")
        if sid in seen:
            err(f"duplicate sample_id '{sid}'")
        seen.add(sid)

        f1 = (r.get("fastq_1") or "").strip()
        f2 = (r.get("fastq_2") or "").strip()
        if not f1 or not f2:
            err(f"sample {sid}: fastq_1 and fastq_2 are both required (paired-end)")
        if f1 == f2:
            err(f"sample {sid}: fastq_1 == fastq_2")
        if a.check_files:
            for f in (f1, f2):
                if not os.path.isfile(f):
                    err(f"sample {sid}: FASTQ not found: {f}")
                if not re.search(r"\.(fastq|fq)(\.gz)?$", f):
                    warn(f"sample {sid}: '{f}' does not look like a FASTQ")

        sex_raw = (r.get("sex") or "").strip().lower() if has_sex else ""
        sex = SEX_MAP.get(sex_raw)
        if sex is None:
            warn(f"sample {sid}: unrecognized sex '{r.get('sex')}' -> U")
            sex = "U"

        pheno = (r.get("phenotype") or "").strip() if has_phenotype else ""
        if pheno:
            pheno_vals.append(pheno)

        rec = {"sample_id": sid, "fastq_1": f1, "fastq_2": f2, "sex": sex,
               "phenotype": pheno, "batch": (r.get("batch") or "").strip() if has_batch else ""}
        for c in cov_cols:
            rec[c] = (r.get(c) or "").strip()
        samples.append(rec)

    # Phenotype typing (binary vs quantitative) — informational.
    pheno_type = "none"
    if pheno_vals:
        uniq = set(pheno_vals)
        if uniq <= {"0", "1", "2", "case", "control", "ctrl"}:
            pheno_type = "binary"
        else:
            try:
                [float(v) for v in pheno_vals]
                pheno_type = "quantitative"
            except ValueError:
                pheno_type = "mixed"
                warn("phenotype column mixes binary/quantitative/text tokens")

    # Normalized CSV
    out_cols = ["sample_id", "fastq_1", "fastq_2", "sex", "phenotype", "batch"] + cov_cols
    with open(a.out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_cols)
        w.writeheader()
        for s in samples:
            w.writerow({c: s.get(c, "") for c in out_cols})

    summary = {
        "n_samples": len(samples),
        "sample_ids": [s["sample_id"] for s in samples],
        "has_phenotype": has_phenotype,
        "phenotype_type": pheno_type,
        "n_with_phenotype": len(pheno_vals),
        "covariate_columns": cov_cols,
        "has_sex": has_sex,
        "has_batch": has_batch,
        "burden_eligible": bool(has_phenotype),
    }
    with open(a.out_json, "w") as fh:
        json.dump(summary, fh, indent=2)

    print(f"[validate_samplesheet] OK: {len(samples)} samples, "
          f"phenotype={'yes('+pheno_type+')' if has_phenotype else 'no -> burden will be skipped'}, "
          f"covariates={len(cov_cols)}")


if __name__ == "__main__":
    main()

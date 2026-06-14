#!/usr/bin/env python3
"""burden_collapse_test.py — Stage 14 collapsing burden test (CallForge).

The lightweight, dependency-free burden engine (CMC-style collapsing): per unit
(gene / gene-set) a 2x2 case/control x carrier/non-carrier chi-square test, p via
erfc (1 df), plus the genomic-inflation factor lambda. Used for small cohorts and
as the always-available default; regenie (Firth, covariate-adjusted) and SKAT-O/
STAAR are the production engines (separate processes) for real cohorts.

NB: a collapsing chi-square does not adjust for covariates/ancestry PCs the way
regenie/SKAT do — small-N results are NOT valid association (see the report caveat).
stdlib only.
"""
import argparse, csv, json, math, os, statistics


def chi2_2x2(a, b, c, d):
    """Pearson chi-square for [[a,b],[c,d]]; returns (chi2, p) with 1 df via erfc."""
    n = a + b + c + d
    r1, r2, c1, c2 = a + b, c + d, a + c, b + d
    if min(r1, r2, c1, c2) == 0 or n == 0:
        return 0.0, 1.0
    chi2 = 0.0
    for obs, (rr, cc) in zip((a, b, c, d), ((r1, c1), (r1, c2), (r2, c1), (r2, c2))):
        exp = rr * cc / n
        if exp > 0:
            chi2 += (obs - exp) ** 2 / exp
    p = math.erfc(math.sqrt(chi2 / 2.0)) if chi2 > 0 else 1.0
    return chi2, p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--phenotype", required=True)
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    pheno = {}
    for r in csv.DictReader(open(a.phenotype), delimiter="\t"):
        try:
            pheno[r["sample_id"]] = int(r["phenotype"])
        except (ValueError, KeyError):
            pass

    rows = list(csv.reader(open(a.matrix), delimiter="\t"))
    results, chi2s = [], []
    if rows:
        header = rows[0]
        samples = header[1:]
        for row in rows[1:]:
            if not row:
                continue
            unit = row[0]
            carr = {s: int(v) for s, v in zip(samples, row[1:])}
            # 2x2: case/control x carrier/non-carrier (only samples with a phenotype)
            cc = {(p, c): 0 for p in (1, 0) for c in (1, 0)}
            for s, ph in pheno.items():
                if s in carr:
                    cc[(ph, 1 if carr[s] else 0)] += 1
            a_ = cc[(1, 1)]; b_ = cc[(1, 0)]; c_ = cc[(0, 1)]; d_ = cc[(0, 0)]
            chi2, p = chi2_2x2(a_, b_, c_, d_)
            n_carr = a_ + c_
            if n_carr > 0:
                chi2s.append(chi2)
            results.append({"unit": unit, "kind": unit.split(":")[0], "engine": "collapse_chi2",
                            "n_carriers": n_carr, "case_carriers": a_, "control_carriers": c_,
                            "chi2": round(chi2, 4),
                            "p": ("NA" if n_carr == 0 else round(p, 6))})

    # genomic inflation lambda (median chi2 / median of chi2_1 = 0.4549)
    lam = round(statistics.median(chi2s) / 0.4549, 4) if chi2s else None

    results.sort(key=lambda r: (r["p"] == "NA", r["p"] if r["p"] != "NA" else 1))
    with open(os.path.join(a.outdir, "burden_results.tsv"), "w") as fh:
        # 'engine' column so a result is never mistaken for a regenie association
        # when the file is read without the dashboard.
        cols = ["unit", "kind", "engine", "n_carriers", "case_carriers", "control_carriers", "chi2", "p"]
        fh.write("\t".join(cols) + "\n")
        for r in results:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    summary = {"engine": "collapse_chi2", "n_units": len(results),
               "n_tested": len(chi2s), "genomic_inflation_lambda": lam,
               "top": results[:10],
               "caveat": ("Collapsing chi-square; NOT covariate/PC-adjusted. Small-N results are "
                          "machinery-only, not valid association. Use regenie/SKAT-O on a real cohort.")}
    with open(os.path.join(a.outdir, "burden_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"[burden_collapse_test] {len(results)} units, {len(chi2s)} tested, lambda={lam}")


if __name__ == "__main__":
    main()

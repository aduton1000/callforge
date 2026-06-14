#!/usr/bin/env python3
"""verify_annotation.py — Stage 11 annotation-landing proof (CallForge).

After VEP + vcfanno + PhyloP, PROVE that each source's values actually attached
to variants (not silently missed on a chr/no-chr mismatch). For each applied
source, count variants carrying its field and show concrete examples; emit
per-source % annotated + a per-variant landing table. Fails loud if an APPLIED
source landed on 0 variants (the silent-miss symptom).
"""
import argparse, json, os, subprocess, sys


def q(vcf, fmt):
    return subprocess.run(["bcftools", "query", "-f", fmt, vcf],
                          check=True, capture_output=True, text=True).stdout.splitlines()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--sources", required=True, help="annotation_sources.json from annotate_dbs.py")
    ap.add_argument("--gnomad-af-field", default="AF_afr")
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    applied = set(json.load(open(a.sources)).get("sources_applied", []))

    total = len(q(a.vcf, "%CHROM\n"))
    afr = f"gnomAD_{a.gnomad_af_field}"
    # probes: field -> bcftools query expr that is non-'.' when present
    probes = {
        "vep_csq":  ("CSQ",          "%INFO/CSQ\n"),
        "gnomad":   (afr,            f"%INFO/{afr}\n"),
        "clinvar":  ("ClinVar_CLNSIG","%INFO/ClinVar_CLNSIG\n"),
        "dbsnp":    ("ID(rsID)",     "%ID\n"),
        "phylop":   ("PhyloP",       "%INFO/PhyloP\n"),
    }
    # which to check: VEP/CSQ always; others only if applied
    check = {"vep_csq": True, "gnomad": "gnomad" in applied, "clinvar": "clinvar" in applied,
             "dbsnp": "dbsnp" in applied, "phylop": "phylop" in applied}

    report = {"n_variants": total, "per_source": {}, "failures": []}
    for key, do in check.items():
        field, fmt = probes[key]
        if not do:
            report["per_source"][key] = {"field": field, "status": "not_applied"}
            continue
        try:
            vals = q(a.vcf, fmt)
        except subprocess.SubprocessError:
            report["per_source"][key] = {"field": field, "status": "field_absent_in_header"}
            report["failures"].append(f"{key}: field {field} not in VCF header though source applied")
            continue
        if key == "dbsnp":
            present = sum(1 for v in vals if v.startswith("rs"))
        else:
            present = sum(1 for v in vals if v not in (".", "", None))
        pct = round(100 * present / total, 2) if total else 0.0
        report["per_source"][key] = {"field": field, "n_annotated": present,
                                     "pct_annotated": pct, "status": "ok" if present else "LANDED_ZERO"}
        if present == 0:
            report["failures"].append(f"{key}: APPLIED but landed on 0/{total} variants "
                                      f"(silent-miss / contig-name mismatch?)")

    # per-variant landing table (first 25 variants x sources)
    cols = ["CHROM", "POS", "REF", "ALT", "ID"]
    fmt = "%CHROM\t%POS\t%REF\t%ALT\t%ID"
    extra = []
    if check["gnomad"]:  fmt += f"\t%INFO/{afr}";          extra.append(afr)
    if check["clinvar"]: fmt += "\t%INFO/ClinVar_CLNSIG";  extra.append("ClinVar_CLNSIG")
    if check["phylop"]:  fmt += "\t%INFO/PhyloP";          extra.append("PhyloP")
    fmt += "\n"
    rows = q(a.vcf, fmt)[:25]
    with open(os.path.join(a.outdir, "annotation_landing_table.tsv"), "w") as fh:
        fh.write("\t".join(cols + extra + ["has_CSQ"]) + "\n")
        csq = q(a.vcf, "%CHROM:%POS\t%INFO/CSQ\n")
        csq_map = {l.split("\t")[0]: (l.split("\t")[1] not in (".", "")) for l in csq if "\t" in l}
        for r in rows:
            f = r.split("\t")
            key = f"{f[0]}:{f[1]}"
            fh.write(r + "\t" + ("yes" if csq_map.get(key) else "no") + "\n")

    with open(os.path.join(a.outdir, "annotation_landing.json"), "w") as fh:
        json.dump(report, fh, indent=2)

    print(f"[verify_annotation] {total} variants; "
          + "; ".join(f"{k}={v.get('pct_annotated','-')}%" for k, v in report["per_source"].items()
                      if v.get("status") == "ok"))
    if report["failures"]:
        sys.stderr.write("[verify_annotation] LANDING FAILURES:\n  " + "\n  ".join(report["failures"]) + "\n")
        sys.exit(3)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""build_str_catalog.py — Stage 9 fallback ExpansionHunter catalog (CallForge).

Used ONLY when no curated --str_catalog is supplied. Builds a best-effort
ExpansionHunter variant catalog from the STR-class targets in gene_metadata.tsv
(str_loci column). The repeat motif is ASSUMED (default GT) and flagged in a
warning — for production, supply a curated catalog (PIEZO1 E756del / HMOX1 (GT)n /
SLC11A1 (GT)n etc.) via --str_catalog, since motif/structure are locus-specific.
"""
import argparse, csv, json, sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--motif", default="GT")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    catalog = []
    with open(a.metadata) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("is_str_target") != "yes":
                continue
            for locus in (r.get("str_loci") or "").split(";"):
                locus = locus.strip()
                if not locus:
                    continue
                catalog.append({"LocusId": f"{r['gene']}_{locus.replace(':','_')}",
                                "LocusStructure": f"({a.motif})*",
                                "ReferenceRegion": locus, "VariantType": "Repeat"})
    with open(a.out, "w") as fh:
        json.dump(catalog, fh, indent=2)
    sys.stderr.write(f"[build_str_catalog] WARNING: motif assumed '{a.motif}' for {len(catalog)} "
                     f"loci — supply a curated --str_catalog for production.\n")
    print(f"[build_str_catalog] wrote {len(catalog)} loci to {a.out}")


if __name__ == "__main__":
    main()

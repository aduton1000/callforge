#!/usr/bin/env python3
"""parse_bcftools_stats.py — Stage 6/7 callset stats -> JSON (CallForge).

Parses `bcftools stats -s -` output into a tidy JSON used by the calling/filter
plots and the gate report:
  totals (records, SNPs, indels, multiallelic), overall Ti/Tv,
  per-sample: nNonRefHom, nHets, het/hom ratio, Ti/Tv, nIndels, nMissing,
              variant count (= nNonRefHom + nHets), avg depth.
"""
import argparse, json, sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stats", required=True, help="bcftools stats -s - output")
    ap.add_argument("--label", default="callset")
    ap.add_argument("--out-json", required=True)
    a = ap.parse_args()

    sn, tstv, per_sample = {}, {}, {}
    for line in open(a.stats):
        if line.startswith("SN\t"):
            f = line.rstrip("\n").split("\t")
            key = f[2].replace("number of ", "").rstrip(":")
            try:
                sn[key] = int(f[3])
            except ValueError:
                pass
        elif line.startswith("TSTV\t"):
            f = line.rstrip("\n").split("\t")
            # TSTV id ts tv ts/tv ts(1stALT) tv(1stALT) ts/tv(1stALT)
            tstv = {"ts": int(f[2]), "tv": int(f[3]), "ts_tv": float(f[4])}
        elif line.startswith("PSC\t"):
            f = line.rstrip("\n").split("\t")
            # PSC id sample nRefHom nNonRefHom nHets nTs nTv nIndels avgDepth nSingletons nHapRef nHapAlt nMissing
            sample = f[2]
            nnonrefhom = int(f[4]); nhets = int(f[5])
            nts = int(f[6]); ntv = int(f[7])
            per_sample[sample] = {
                "nNonRefHom": nnonrefhom, "nHets": nhets,
                "het_hom_ratio": round(nhets / nnonrefhom, 3) if nnonrefhom else None,
                "nTs": nts, "nTv": ntv,
                "ts_tv": round(nts / ntv, 3) if ntv else None,
                "nIndels": int(f[8]),
                "avgDepth": float(f[9]) if f[9] not in ("", ".") else None,
                "nMissing": int(f[13]) if len(f) > 13 else None,
                "n_variants": nnonrefhom + nhets,
            }

    out = {"label": a.label,
           "n_records": sn.get("records"),
           "n_snps": sn.get("SNPs"),
           "n_indels": sn.get("indels"),
           "n_multiallelic": sn.get("multiallelic sites"),
           "n_samples": len(per_sample),
           "ts_tv_overall": tstv.get("ts_tv"),
           "ts_overall": tstv.get("ts"), "tv_overall": tstv.get("tv"),
           "per_sample": per_sample,
           "samples": sorted(per_sample)}
    with open(a.out_json, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"[parse_bcftools_stats] {a.label}: records={out['n_records']} "
          f"SNPs={out['n_snps']} indels={out['n_indels']} Ti/Tv={out['ts_tv_overall']} "
          f"samples={out['n_samples']} ({','.join(out['samples'])})")


if __name__ == "__main__":
    main()

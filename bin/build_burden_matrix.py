#!/usr/bin/env python3
"""build_burden_matrix.py — Stage 14 burden input matrix (CallForge).

Builds the rare-variant burden analysis matrix from the annotated callset:
  1. keep RARE (gnomAD-AFR AF < af_max, or novel/absent) + FUNCTIONAL (VEP
     consequence in the qualifying set) variants;
  2. collapse to carrier indicators per gene AND per gene-set (CaptureForge
     burden_group from gene_metadata);
  3. assemble phenotype (case/control) + covariates (covariate_* + ancestry PCs).

REQUIRES a phenotype: if the sample sheet has none, writes a 'skipped_no_phenotype'
marker and empty outputs (the pipeline then stops at the annotated callset).
stdlib + bcftools.
"""
import argparse, csv, json, os, re, subprocess
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--metadata", required=True)       # gene_metadata.tsv (burden_group)
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--cohort-qc", default=None)        # cohort_qc.json (ancestry PCs)
    ap.add_argument("--af-field", default="AF_afr")
    ap.add_argument("--af-max", type=float, default=0.01)
    ap.add_argument("--csq", required=True)             # comma list of qualifying consequences
    ap.add_argument("--n-pcs", type=int, default=4)
    ap.add_argument("--covariates", default="")         # comma list of covariate_* cols ('' = all)
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    afr = f"gnomAD_{a.af_field}"
    csq_set = set(c.strip() for c in a.csq.split(",") if c.strip())

    # sample sheet: phenotype + covariates
    rows = list(csv.DictReader(open(a.samplesheet)))
    cov_cols = [c for c in (rows[0].keys() if rows else []) if c.startswith("covariate_")]
    if a.covariates:
        cov_cols = [c for c in a.covariates.split(",") if c in cov_cols]
    pheno = {}
    for r in rows:
        v = (r.get("phenotype") or "").strip().lower()
        pheno[r["sample_id"]] = {"case": 1, "2": 1, "control": 0, "ctrl": 0, "1": 0, "0": 0}.get(v)
    has_pheno = any(v is not None for v in pheno.values())
    meta_info = {"af_field": afr, "af_max": a.af_max, "qualifying_csq": sorted(csq_set)}

    if not has_pheno:
        meta_info["status"] = "skipped_no_phenotype"
        with open(os.path.join(a.outdir, "burden_meta.json"), "w") as fh:
            json.dump(meta_info, fh, indent=2)
        for f in ("burden_matrix.tsv", "phenotype.tsv", "covariates.tsv", "qualifying_variants.tsv"):
            open(os.path.join(a.outdir, f), "w").close()
        print("[build_burden_matrix] no phenotype in sample sheet -> burden SKIPPED")
        return

    # gene -> burden group
    group = {}
    for r in csv.DictReader(open(a.metadata), delimiter="\t"):
        group[r["gene"]] = r.get("burden_group") or r["gene"]

    # ancestry PCs
    pcs = {}
    if a.cohort_qc and os.path.isfile(a.cohort_qc):
        pcs = (json.load(open(a.cohort_qc)).get("pca") or {})

    samples = subprocess.run(["bcftools", "query", "-l", a.vcf], check=True,
                             capture_output=True, text=True).stdout.split()
    # CSQ format
    h = subprocess.run(["bcftools", "view", "-h", a.vcf], capture_output=True, text=True).stdout
    m = re.search(r"##INFO=<ID=CSQ.*?Format: ([^\">]+)", h)
    cfmt = m.group(1).split("|") if m else []
    ci = cfmt.index("Consequence") if "Consequence" in cfmt else 1
    gi = cfmt.index("SYMBOL") if "SYMBOL" in cfmt else 0

    q = subprocess.run(["bcftools", "query", "-f",
                        f"%CHROM\t%POS\t%INFO/{afr}\t%INFO/CSQ[\t%GT]\n", a.vcf],
                       capture_output=True, text=True).stdout.splitlines()

    # unit -> sample -> carrier(0/1); also collect qualifying variants
    units = defaultdict(lambda: defaultdict(int))   # "gene:NAME" / "group:NAME"
    qual_variants = []
    variant_geno = []    # variant-level dosage (for regenie/SKAT engines)
    for line in q:
        f = line.split("\t")
        if len(f) < 4 + len(samples):
            continue
        chrom, pos, af, csq = f[0], f[1], f[2], f[3]
        gts = f[4:4 + len(samples)]
        # rare?
        rare = True
        if af not in (".", "", None):
            try:
                rare = float(af) < a.af_max
            except ValueError:
                rare = True
        if not rare:
            continue
        # functional? (any transcript consequence in the set) + gene
        gene, cons_hit = None, None
        for tx in csq.split(","):
            parts = tx.split("|")
            if len(parts) <= max(ci, gi):
                continue
            for cons in parts[ci].split("&"):
                if cons in csq_set:
                    cons_hit = cons; gene = parts[gi] or gene
                    break
            if cons_hit:
                gene = parts[gi] or gene; break
        if not cons_hit or not gene:
            continue
        grp = group.get(gene, gene)
        qual_variants.append((chrom, pos, gene, grp, cons_hit, af))
        dos = []
        for s, gt in zip(samples, gts):
            al = gt.replace("|", "/").split("/")
            d = 0 if ("." in al or gt in (".", "./.")) else sum(1 for x in al if x not in ("0", "."))
            dos.append(d)
            if d > 0:
                units[f"gene:{gene}"][s] = 1
                units[f"group:{grp}"][s] = 1
        variant_geno.append((f"{chrom}:{pos}", gene, grp, dos))

    # write outputs
    with open(os.path.join(a.outdir, "qualifying_variants.tsv"), "w") as fh:
        fh.write("chrom\tpos\tgene\tburden_group\tconsequence\tgnomAD_AFR_AF\n")
        for v in qual_variants:
            fh.write("\t".join(map(str, v)) + "\n")

    # variant-level dosage matrix (for regenie/SKAT, which group variants per gene)
    with open(os.path.join(a.outdir, "genotypes.tsv"), "w") as fh:
        fh.write("variant\tgene\tburden_group\t" + "\t".join(samples) + "\n")
        for (vid, gene, grp, dos) in variant_geno:
            fh.write(f"{vid}\t{gene}\t{grp}\t" + "\t".join(map(str, dos)) + "\n")

    unit_names = sorted(units)
    with open(os.path.join(a.outdir, "burden_matrix.tsv"), "w") as fh:
        fh.write("unit\t" + "\t".join(samples) + "\n")
        for u in unit_names:
            fh.write(u + "\t" + "\t".join(str(units[u].get(s, 0)) for s in samples) + "\n")

    with open(os.path.join(a.outdir, "phenotype.tsv"), "w") as fh:
        fh.write("sample_id\tphenotype\n")
        for s in samples:
            fh.write(f"{s}\t{pheno.get(s, '')}\n")

    npc = min(a.n_pcs, max((len(v) for v in pcs.values()), default=0))
    with open(os.path.join(a.outdir, "covariates.tsv"), "w") as fh:
        hdr = ["sample_id"] + cov_cols + [f"PC{i+1}" for i in range(npc)]
        fh.write("\t".join(hdr) + "\n")
        covmap = {r["sample_id"]: r for r in rows}
        for s in samples:
            vals = [covmap.get(s, {}).get(c, "") for c in cov_cols]
            vals += [f"{pcs.get(s, [0]*npc)[i]:.5f}" if s in pcs else "0" for i in range(npc)]
            fh.write("\t".join([s] + [str(v) for v in vals]) + "\n")

    meta_info.update({"status": "ok", "n_samples": len(samples),
                      "n_cases": sum(1 for v in pheno.values() if v == 1),
                      "n_controls": sum(1 for v in pheno.values() if v == 0),
                      "n_qualifying_variants": len(qual_variants),
                      "n_units": len(unit_names),
                      "n_gene_units": sum(1 for u in unit_names if u.startswith("gene:")),
                      "n_groupset_units": sum(1 for u in unit_names if u.startswith("group:")),
                      "covariates": cov_cols, "n_pcs": npc})
    with open(os.path.join(a.outdir, "burden_meta.json"), "w") as fh:
        json.dump(meta_info, fh, indent=2)
    print(f"[build_burden_matrix] {len(qual_variants)} qualifying variants -> {len(unit_names)} units "
          f"(N={len(samples)}: {meta_info['n_cases']} case / {meta_info['n_controls']} control)")


if __name__ == "__main__":
    main()

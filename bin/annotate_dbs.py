#!/usr/bin/env python3
"""annotate_dbs.py — Stage 11 vcfanno/dbSNP/PhyloP annotation with chr-reconciliation.

Adds, to a (VEP-annotated) callset, from manifest-discovered resources:
  - gnomAD population AFs (namespaced gnomAD_*, so the cohort's own AF/AC/AN are
    NEVER overwritten — incl. the AFR field) via vcfanno
  - ClinVar clinical significance fields via vcfanno
  - dbSNP rsIDs into the ID column via bcftools annotate
  - PhyloP conservation score via pyBigWig (per-variant)

CHR-RECONCILIATION: the design reference / callset is Ensembl-named (1,2,..),
gnomAD/PhyloP are usually chr-prefixed, NCBI ClinVar is Ensembl-named. Naming is
detected PER resource; the query is renamed to a resource's style only for the
resources that differ, then renamed back — so an annotation is never silently
missed on a contig-name mismatch. Sources that are missing/unfit are skipped
(organism-generic: non-human resources simply absent).

Requires bcftools/vcfanno/tabix; pyBigWig only if a PhyloP bigWig is supplied.
"""
import argparse, json, os, re, shutil, subprocess, sys


def sh(cmd, **kw):
    return subprocess.run(cmd, check=True, capture_output=True, text=True, **kw).stdout


def vcf_contigs(vcf):
    out = sh(["bcftools", "view", "-h", vcf])
    cs = re.findall(r"##contig=<ID=([^,>]+)", out)
    if not cs:  # fall back to first record
        rec = subprocess.run(["bcftools", "view", "-H", vcf], capture_output=True, text=True).stdout.splitlines()
        cs = [rec[0].split("\t")[0]] if rec else []
    return cs


def style_of(contigs):
    chrpref = sum(1 for c in contigs if c.lower().startswith("chr"))
    return "ucsc" if chrpref > len(contigs) / 2 else "ensembl"


def db_style(path):
    """Contig style of a tabixed VCF (from the index)."""
    try:
        out = sh(["tabix", "-l", path])
        cs = [x.strip() for x in out.splitlines() if x.strip()]
        return style_of(cs) if cs else "unknown"
    except subprocess.SubprocessError:
        return "unknown"


def header_info_fields(path):
    out = sh(["bcftools", "view", "-h", path])
    return set(re.findall(r"##INFO=<ID=([^,]+)", out))


def build_rename_map(fai, to_style, path):
    """Write a bcftools --rename-chrs map (current->target) from the .fai contigs."""
    lines = []
    with open(fai) as fh:
        for ln in fh:
            name = ln.split("\t")[0]
            if to_style == "ucsc":
                tgt = name if name.lower().startswith("chr") else ("chrM" if name == "MT" else "chr" + name)
            else:  # ensembl
                tgt = name[3:] if name.lower().startswith("chr") else name
                if tgt == "M":
                    tgt = "MT"
            if tgt != name:
                lines.append(f"{name}\t{tgt}")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + ("\n" if lines else ""))
    return path


def rename(vcf, mapfile, out):
    sh(["bcftools", "annotate", "--rename-chrs", mapfile, vcf, "-Oz", "-o", out])
    sh(["tabix", "-f", "-p", "vcf", out])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--gnomad-af-field", default="AF_afr")
    ap.add_argument("--phylop", default=None, help="PhyloP bigWig (optional)")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--out", default="annotated.vcf.gz")
    a = ap.parse_args()
    od = a.outdir
    os.makedirs(od, exist_ok=True)

    M = json.load(open(a.manifest)).get("resources", {})
    def res(name):
        r = M.get(name, {})
        p = r.get("path")
        return p if (p and r.get("status") in ("ok", "unindexed")) else None
    gnomad = res("gnomad_vcf"); clinvar = res("clinvar_vcf"); dbsnp = res("dbsnp_vcf")
    phylop = a.phylop if (a.phylop and os.path.isfile(a.phylop)) else res("phylop_bw")

    qstyle = style_of(vcf_contigs(a.vcf))
    applied = {"sources_applied": [], "sources_skipped": [], "query_style": qstyle,
               "gnomad_af_field": a.gnomad_af_field}

    # group annotations by required contig style
    cur = a.vcf
    work = []  # cleanup

    def ensure_style(vcf_in, target_style, tag):
        """Return a vcf whose contigs are in target_style (rename if needed)."""
        if qstyle == target_style:
            return vcf_in, False
        mp = build_rename_map(a.fai, target_style, os.path.join(od, f".rn_{tag}.txt"))
        out = os.path.join(od, f".q_{tag}.vcf.gz")
        rename(vcf_in, mp, out); work.append(out)
        return out, True

    # ---- gnomAD + ClinVar via vcfanno, grouped by their style ----
    # Build per-style vcfanno configs.
    def vcfanno_block(path, prefix, want_fields, ops_default="self"):
        have = header_info_fields(path)
        fields = [f for f in want_fields if f in have]
        if not fields:
            return None, []
        names = [f"{prefix}{f}" for f in fields]
        block = ("[[annotation]]\n"
                 f'file="{os.path.abspath(path)}"\n'
                 f'fields=[{", ".join(chr(34)+f+chr(34) for f in fields)}]\n'
                 f'names=[{", ".join(chr(34)+n+chr(34) for n in names)}]\n'
                 f'ops=[{", ".join(chr(34)+ops_default+chr(34) for _ in fields)}]\n')
        return block, names

    # Determine each resource's style
    styles = {}
    for nm, p in (("gnomad", gnomad), ("clinvar", clinvar), ("dbsnp", dbsnp)):
        if p:
            styles[nm] = db_style(p)

    # Process resources grouped by style so we rename the query at most once per style.
    for target_style in sorted(set(styles.values())):
        members = [nm for nm, s in styles.items() if s == target_style]
        # build vcfanno config for gnomad/clinvar members
        cfg_blocks, applied_names = [], []
        gnomad_fields = ["AF", a.gnomad_af_field, "AF_grpmax", "AF_popmax", "nhomalt"]
        if "gnomad" in members:
            blk, names = vcfanno_block(gnomad, "gnomAD_", gnomad_fields)
            if blk: cfg_blocks.append(blk); applied_names += names; applied["sources_applied"].append("gnomad")
            else: applied["sources_skipped"].append("gnomad(no fields)")
        if "clinvar" in members:
            blk, names = vcfanno_block(clinvar, "ClinVar_",
                                       ["CLNSIG", "CLNSIGCONF", "CLNREVSTAT", "CLNDN", "CLNHGVS"])
            if blk: cfg_blocks.append(blk); applied_names += names; applied["sources_applied"].append("clinvar")
            else: applied["sources_skipped"].append("clinvar(no fields)")

        if not cfg_blocks and "dbsnp" not in members:
            continue
        qstyled, renamed = ensure_style(cur, target_style, target_style)

        # vcfanno (gnomad/clinvar)
        if cfg_blocks:
            toml = os.path.join(od, f".vcfanno_{target_style}.toml")
            with open(toml, "w") as fh:
                fh.write("\n".join(cfg_blocks))
            outv = os.path.join(od, f".anno_{target_style}.vcf.gz")
            with open(outv, "wb") as ofh:
                p1 = subprocess.Popen(["vcfanno", toml, qstyled], stdout=subprocess.PIPE)
                p2 = subprocess.Popen(["bgzip"], stdin=p1.stdout, stdout=ofh)
                p1.stdout.close(); p2.communicate()
            sh(["tabix", "-f", "-p", "vcf", outv]); work.append(outv)
            qstyled = outv

        # dbSNP rsID via bcftools annotate
        if "dbsnp" in members:
            outd = os.path.join(od, f".dbsnp_{target_style}.vcf.gz")
            sh(["bcftools", "annotate", "-a", dbsnp, "-c", "ID", qstyled, "-Oz", "-o", outd])
            sh(["tabix", "-f", "-p", "vcf", outd]); work.append(outd)
            qstyled = outd; applied["sources_applied"].append("dbsnp")

        # rename back to query style if we renamed
        if renamed:
            mp = build_rename_map(a.fai, qstyle, os.path.join(od, f".rnback_{target_style}.txt"))
            outb = os.path.join(od, f".back_{target_style}.vcf.gz")
            rename(qstyled, mp, outb); work.append(outb)
            cur = outb
        else:
            cur = qstyled

    for nm in ("gnomad", "clinvar", "dbsnp"):
        src = {"gnomad": gnomad, "clinvar": clinvar, "dbsnp": dbsnp}[nm]
        if not src:
            applied["sources_skipped"].append(f"{nm}(missing/unfit)")

    # ---- PhyloP via pyBigWig (bigWig is chr-prefixed; reconcile by query naming) ----
    if phylop:
        try:
            import pyBigWig  # noqa
            outp = os.path.join(od, ".phylop.vcf.gz")
            rc = subprocess.run([sys.executable, os.path.join(os.path.dirname(__file__), "annotate_phylop.py"),
                                 "--vcf", cur, "--bigwig", phylop, "--out", outp], capture_output=True, text=True)
            if rc.returncode == 0 and os.path.isfile(outp):
                cur = outp; work.append(outp); applied["sources_applied"].append("phylop")
            else:
                applied["sources_skipped"].append("phylop(annotation error)")
                sys.stderr.write(rc.stderr[-500:] + "\n")
        except ImportError:
            applied["sources_skipped"].append("phylop(no pyBigWig)")
    else:
        applied["sources_skipped"].append("phylop(missing)")

    # finalize
    final = os.path.join(od, a.out)
    sh(["bcftools", "view", cur, "-Oz", "-o", final])
    sh(["tabix", "-f", "-p", "vcf", final])
    applied["output"] = final
    with open(os.path.join(od, "annotation_sources.json"), "w") as fh:
        json.dump(applied, fh, indent=2)
    # cleanup intermediates
    for f in work:
        for ext in ("", ".tbi"):
            try: os.remove(f + ext)
            except OSError: pass
    print(f"[annotate_dbs] applied={applied['sources_applied']} skipped={applied['sources_skipped']}")


if __name__ == "__main__":
    main()

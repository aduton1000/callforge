#!/usr/bin/env python3
"""resolve_variants.py — place planted variants at REAL GRCh38 coords (cluster step).

Resolves the cohort_design.py PLAN against the REAL reference + gene model: for every
planted coding/promoter variant it finds a real position of the chosen transcript that
ACTUALLY produces the intended consequence (verified by translating the real CDS), is
mappable (not soft-masked/N/low-complexity), and is NOVEL (absent from an optional
known-variants VCF so it reads as rare/qualifying). It then writes the resolved truth
manifest, the real-coordinate target BED, the real HTT STR catalog, and gene_metadata.tsv.

Nothing is faked: positions, ref/alt and consequences are real-reference facts. The
read simulator (simulate_reads.py) and the pipeline both operate on these real coords.

Inputs:
  --fasta    real GRCh38 no-alt FASTA (.fai beside it)              [REQUIRED]
  --gtf      gene-model GTF (Ensembl, gzip ok) — authoritative CDS/exons  [REQUIRED]
  --known    optional gnomAD/dbSNP VCF(.gz, tabixed) for novelty/collision checks
  --outdir   fixture dir (writes ref/, truth_manifest.json, sim.params.yaml)
  --make-gene-metadata  path to bin/make_gene_metadata.py (default: repo bin/)

Requires pysam (in the CallForge .sif). Run inside the container (generate_reads.sh does).
"""
import argparse
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cohort_design as D

CODON = {  # standard genetic code
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L",
    "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V",
    "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P",
    "CCC": "P", "CCA": "P", "CCG": "P", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*",
    "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "TGT": "C",
    "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G",
    "GGG": "G",
}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
FLANK = 150            # exon padding in the target BED (so reads anchor)
PROMOTER_OFFSET = 900  # bp upstream of TSS for the promoter SNV (<5kb -> upstream_gene_variant)
MAX_EXONS = 8          # cap exons simulated per gene (CNVkit needs several bins; keeps it fast)


def comp(b):
    return b.translate(COMP)


def revcomp(s):
    return s.translate(COMP)[::-1]


def attr(s, key):
    i = s.find(key + ' "')
    if i < 0:
        return None
    j = s.find('"', i + len(key) + 2)
    return s[i + len(key) + 2:j]


def parse_gtf(path, genes):
    """Return {gene: {tx_id: {'strand':..,'chrom':..,'exon':[(s,e)],'cds':[(s,e)]}}} for wanted genes."""
    want = set(genes)
    want_tx = {m["tx"] for m in genes.values()}
    out = {}
    op = (lambda p: __import__("gzip").open(p, "rt")) if path.endswith(".gz") else (lambda p: open(p))
    with op(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("exon", "CDS"):
                continue
            gname = attr(f[8], "gene_name")
            if gname not in want:
                continue
            txid = (attr(f[8], "transcript_id") or "").split(".")[0]
            g = out.setdefault(gname, {})
            t = g.setdefault(txid, {"strand": f[6], "chrom": f[0], "exon": [], "cds": []})
            t[f[2].lower() if f[2] == "exon" else "cds"].append((int(f[3]), int(f[4])))
    return out, want_tx


def pick_tx(gene, txs, want_tx):
    pref = D.GENES[gene]["tx"].split(".")[0]
    if pref in txs and txs[pref]["cds"]:
        return pref
    # fallback: transcript with the longest CDS
    coding = [(tid, sum(e - s + 1 for s, e in t["cds"])) for tid, t in txs.items() if t["cds"]]
    if coding:
        return max(coding, key=lambda x: x[1])[0]
    # no CDS (e.g. HTT chosen tx may lack CDS in this GTF) -> longest-exon transcript
    return max(txs, key=lambda tid: sum(e - s + 1 for s, e in txs[tid]["exon"]))


def build_cds(fa, tx):
    """Return (cds_seq, coords) in CODING order. coords[i]=(genomic_pos_1based, fwd_ref_base)."""
    chrom, strand = tx["chrom"], tx["strand"]
    exons = sorted(tx["cds"])
    seq, coords = [], []
    order = exons if strand == "+" else list(reversed(exons))
    for (s, e) in order:
        block = fa.fetch(chrom, s - 1, e).upper()   # forward strand, 1-based inclusive
        positions = range(s, e + 1) if strand == "+" else range(e, s - 1, -1)
        for gp in positions:
            fwd = block[gp - s] if strand == "+" else block[gp - s]
            coords.append((gp, fwd))
        seq.append(block if strand == "+" else revcomp(block))
    return "".join(seq), coords


def coding_base_at(coords, cds_idx, strand):
    """The base as seen in CODING orientation at cds index."""
    gp, fwd = coords[cds_idx]
    return fwd if strand == "+" else comp(fwd)


def to_forward(cds_base, strand):
    return cds_base if strand == "+" else comp(cds_base)


def mappable(fa, chrom, pos, win=25):
    """Heuristic: uppercase ACGT, no N nearby, not low-complexity/homopolymer context."""
    s = fa.fetch(chrom, max(0, pos - 1 - win), pos - 1 + win).upper()
    if not s or "N" in s:
        return False
    base = fa.fetch(chrom, pos - 1, pos).upper()
    if base not in "ACGT":
        return False
    # reject if any base makes up >70% of the window (low complexity / homopolymer)
    return max(s.count(b) for b in "ACGT") <= 0.70 * len(s)


def is_known(known, chrom, pos):
    if known is None:
        return False
    try:
        for _ in known.fetch(chrom, pos - 1, pos):
            return True
    except (ValueError, OSError):
        return False
    return False


def find_coding(fa, tx, coords, cds_seq, kind, known, used, start_codon):
    """Find a real coding edit producing `kind`. Returns dict or None.
    used: set of genomic positions already taken in this gene (avoid collisions)."""
    strand = tx["strand"]
    ncod = len(cds_seq) // 3
    lo, hi = max(1, ncod // 6), ncod - ncod // 6        # middle region, skip start/end
    for k in range(max(start_codon, lo), hi):
        codon = cds_seq[3 * k:3 * k + 3]
        if codon not in CODON or CODON[codon] == "*":
            continue
        if kind == "frameshift":
            # 1-bp deletion at the codon's 2nd base, anchored on the previous base (VCF-norm)
            ci = 3 * k + 1
            gp, _ = coords[ci]
            anch_gp = gp - 1 if strand == "+" else gp + 1     # forward anchor = pos-1
            if gp in used or not mappable(fa, tx["chrom"], gp) or is_known(known, tx["chrom"], gp):
                continue
            apos = min(gp, anch_gp)
            ref = fa.fetch(tx["chrom"], apos - 1, apos + 1).upper()   # 2 fwd bases
            if len(ref) != 2 or "N" in ref:
                continue
            return {"pos": apos, "ref": ref, "alt": ref[0], "vartype": "del",
                    "consequence": "frameshift_variant", "genomic": gp}
        for ci in range(3):                                 # which codon base to substitute
            cds_pos = 3 * k + ci
            gp, _ = coords[cds_pos]
            if gp in used or not mappable(fa, tx["chrom"], gp) or is_known(known, tx["chrom"], gp):
                continue
            cur = coding_base_at(coords, cds_pos, strand)
            for nb in "ACGT":
                if nb == cur:
                    continue
                newcod = codon[:ci] + nb + codon[ci + 1:]
                aa0, aa1 = CODON[codon], CODON.get(newcod)
                if aa1 is None:
                    continue
                ok = (kind == "stop" and aa1 == "*") or \
                     (kind == "missense" and aa1 != "*" and aa1 != aa0)
                if not ok:
                    continue
                fref = to_forward(cur, strand)
                falt = to_forward(nb, strand)
                # verify against the real FASTA forward base
                if fa.fetch(tx["chrom"], gp - 1, gp).upper() != fref:
                    continue
                return {"pos": gp, "ref": fref, "alt": falt, "vartype": "snv",
                        "consequence": D.EXPECTED_CSQ[kind], "genomic": gp}
    return None


def find_promoter(fa, tx, used, known):
    chrom, strand = tx["chrom"], tx["strand"]
    exons = sorted(tx["exon"]) or sorted(tx["cds"])
    tss = exons[0][0] if strand == "+" else exons[-1][1]
    base = tss - PROMOTER_OFFSET if strand == "+" else tss + PROMOTER_OFFSET
    for d in range(0, 400):                                 # nudge until mappable/novel/unused
        for pos in (base + d, base - d):
            if pos <= 1 or pos in used:
                continue
            if not mappable(fa, chrom, pos) or is_known(known, chrom, pos):
                continue
            ref = fa.fetch(chrom, pos - 1, pos).upper()
            alt = {"A": "T", "T": "A", "C": "G", "G": "C"}.get(ref)   # transversion
            if alt:
                return {"pos": pos, "ref": ref, "alt": alt, "consequence": "upstream_gene_variant"}
    return None


def select_exons(tx, near_positions):
    """Exons to simulate: those overlapping planted positions + fill to MAX_EXONS."""
    exons = sorted(tx["exon"]) or sorted(tx["cds"])
    keep = [iv for iv in exons if any(iv[0] <= p <= iv[1] for p in near_positions)]
    for iv in exons:
        if len(keep) >= MAX_EXONS:
            break
        if iv not in keep:
            keep.append(iv)
    return sorted(set(keep))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--known", default=None, help="gnomAD/dbSNP VCF.gz (tabixed) for novelty checks")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--make-gene-metadata", default=None)
    a = ap.parse_args()

    try:
        import pysam
    except ImportError:
        sys.exit("[resolve_variants] ERROR: pysam required (run inside the CallForge .sif).")

    here = os.path.dirname(os.path.abspath(__file__))
    mgm = a.make_gene_metadata or os.path.normpath(os.path.join(here, "..", "..", "bin", "make_gene_metadata.py"))
    refdir = os.path.join(a.outdir, "ref")
    os.makedirs(refdir, exist_ok=True)
    fa = pysam.FastaFile(a.fasta)
    known = pysam.TabixFile(a.known) if a.known else None

    gtf, want_tx = parse_gtf(a.gtf, D.GENES)
    txsel = {}
    for g in D.GENES:
        if g not in gtf:
            sys.exit(f"[resolve_variants] ERROR: gene {g} not found in GTF {a.gtf}")
        txsel[g] = gtf[g][pick_tx(g, gtf[g], want_tx)]

    # resolve coding + promoter variants
    coding_truth, promoter_truth = [], []
    used_by_gene = {}
    cds_cache = {}
    # group planted coding by gene to spread codon offsets and avoid collisions
    by_gene = {}
    for v in D.ALL_CODING:
        by_gene.setdefault(v["gene"], []).append(v)
    for gene, variants in by_gene.items():
        tx = txsel[gene]
        if not tx["cds"]:
            sys.exit(f"[resolve_variants] ERROR: chosen transcript for {gene} has no CDS in GTF")
        cds_seq, coords = build_cds(fa, tx)
        cds_cache[gene] = (cds_seq, coords)
        used = used_by_gene.setdefault(gene, set())
        offset = max(1, (len(cds_seq) // 3) // (len(variants) + 4))
        for n, v in enumerate(variants):
            r = find_coding(fa, tx, coords, cds_seq, v["kind"], known, used, start_codon=(n + 1) * offset)
            if r is None:
                sys.exit(f"[resolve_variants] ERROR: could not place {v['id']} ({v['kind']}) in {gene}")
            used.add(r["genomic"])
            coding_truth.append({"id": v["id"], "gene": gene, "contig": tx["chrom"],
                                 "pos": r["pos"], "ref": r["ref"], "alt": r["alt"],
                                 "vartype": r["vartype"], "consequence": r["consequence"],
                                 "zygosity": "het", "carriers": v["carriers"]})
    for v in D.PROMOTER_VARIANTS:
        tx = txsel[v["gene"]]
        used = used_by_gene.setdefault(v["gene"], set())
        r = find_promoter(fa, tx, used, known)
        if r is None:
            sys.exit(f"[resolve_variants] ERROR: could not place promoter variant for {v['gene']}")
        used.add(r["pos"])
        promoter_truth.append({"id": v["id"], "gene": v["gene"], "contig": tx["chrom"],
                               "pos": r["pos"], "ref": r["ref"], "alt": r["alt"],
                               "consequence": "upstream_gene_variant", "zygosity": "het",
                               "carriers": v["carriers"]})

    # STR locus — validate the motif at the real HTT tract
    s = D.STR_LOCUS
    rs, re_ = s["ref_region"]
    tract = fa.fetch(s["chrom"], rs - 1, re_).upper()
    motif_ok = tract.startswith(s["motif"]) and set(tract) <= set("ACGT")
    str_region = f"{s['chrom']}:{rs}-{re_}"
    locus_id = f"{s['gene']}_{s['chrom']}_{rs}_{re_}"
    str_truth = {"gene": s["gene"], "contig": s["chrom"], "repeat_region": str_region,
                 "locus_id": locus_id, "motif": s["motif"], "motif_validated": bool(motif_ok),
                 "normal_units": s["normal_units"], "expanded_units": s["expanded_units"],
                 "carriers": s["carriers"]}
    if not motif_ok:
        sys.stderr.write(f"[resolve_variants] WARN: HTT tract at {str_region} did not start "
                         f"with motif {s['motif']} (got {tract[:9]}…) — verify the STR coords.\n")

    # ---- target BED (real coords) ----
    cnv_pos = {g: [] for g in D.CNV_TARGET_GENES}
    plant_pos = {}
    for t in coding_truth + promoter_truth:
        plant_pos.setdefault(t["gene"], []).append(t["pos"])
    bed_rows = []
    for gene, m in D.GENES.items():
        tx = txsel[gene]
        roles = m["roles"]
        if "str" in roles:
            bed_rows.append((m["chrom"], rs - FLANK, re_ + FLANK, f"{gene}|str"))
            continue
        klass = "cnv" if "cnv" in roles else "coding"
        exons = select_exons(tx, plant_pos.get(gene, []) or [tx["exon"][0][0] if tx["exon"] else tx["cds"][0][0]])
        for (es, ee) in exons:
            bed_rows.append((m["chrom"], es - FLANK, ee + FLANK, f"{gene}|{klass}"))
        if "promoter" in roles:
            for t in promoter_truth:
                if t["gene"] == gene:
                    bed_rows.append((m["chrom"], t["pos"] - FLANK, t["pos"] + FLANK, f"{gene}|promoter"))
    targets_bed = os.path.join(refdir, "targets.bed")
    with open(targets_bed, "w") as fh:
        for (c, st, en, name) in sorted(set(bed_rows), key=lambda r: (r[0], r[1])):
            fh.write(f"{c}\t{max(0, st)}\t{en}\t{name}\t0\t+\n")

    # ---- STR catalog (real HTT) ----
    catalog = [{"LocusId": locus_id, "LocusStructure": f"({s['motif']})*",
                "ReferenceRegion": str_region, "VariantType": "Repeat"}]
    catalog_path = os.path.join(refdir, "str_catalog.json")
    with open(catalog_path, "w") as fh:
        json.dump(catalog, fh, indent=2)

    # ---- gene_metadata.tsv via the REAL builder ----
    strloci = os.path.join(refdir, "_str_loci.tsv")
    with open(strloci, "w") as fh:
        fh.write(f"{s['gene']}\t{str_region}\n")
    burden = os.path.join(refdir, "_burden_groups.tsv")
    with open(burden, "w") as fh:
        fh.write(f"{D.BURDEN_GENE}\t{D.BURDEN_GENE}\n")
    meta_tsv = os.path.join(refdir, "gene_metadata.tsv")
    subprocess.run([sys.executable, mgm, "--bed", targets_bed,
                    "--cnv-genes", ",".join(D.CNV_TARGET_GENES),
                    "--str-loci", strloci, "--burden-groups", burden,
                    "--out-tsv", meta_tsv, "--out-json", os.path.join(refdir, "gene_metadata.json"),
                    "--overwrite"], check=True)

    # ---- resolved truth manifest ----
    manifest = {
        "seed": D.SEED, "genome_build": D.GENOME_BUILD, "_status": "RESOLVED (real GRCh38)",
        "genome_fasta": os.path.abspath(a.fasta), "gtf": os.path.abspath(a.gtf),
        "target_bed": os.path.abspath(targets_bed),
        "gene_metadata": os.path.abspath(meta_tsv),
        "str_catalog": os.path.abspath(catalog_path),
        "sex_regions": D.SEX_REGIONS,
        "samples": D.SAMPLES, "cases": D.CASES, "controls": D.CONTROLS,
        "coding_variants": coding_truth, "promoter_variants": promoter_truth,
        "str": str_truth, "cnv_events": D.CNV_EVENTS,
        "burden": {"gene": D.BURDEN_GENE,
                   "qualifying_feature_ids": [v["id"] for v in D.BURDEN_VARIANTS],
                   "case_carriers": sorted({c for v in D.BURDEN_VARIANTS
                                            for c in v["carriers"] if c in D.CASES}),
                   "control_carriers": sorted({c for v in D.BURDEN_VARIANTS
                                               for c in v["carriers"] if c in D.CONTROLS})},
        "cohortqc": {"sex_by_sample": {x["sample_id"]: x["sex"] for x in D.SAMPLES}},
        "stage_refs": D.STAGE_REFS, "thresholds": D.THRESHOLDS,
    }
    with open(os.path.join(a.outdir, "truth_manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)
    write_manifest_tsv(os.path.join(a.outdir, "truth_manifest.tsv"), manifest)

    print(f"[resolve_variants] resolved against real {D.GENOME_BUILD}: "
          f"{len(coding_truth)} coding + {len(promoter_truth)} promoter variants placed & "
          f"consequence-verified, STR motif_validated={motif_ok}, "
          f"{len(set(bed_rows))} BED intervals -> {a.outdir}/truth_manifest.json")
    for t in coding_truth + promoter_truth:
        print(f"  {t['id']:14s} {t['gene']:6s} {t['contig']}:{t['pos']} "
              f"{t['ref']}>{t['alt']} {t['consequence']}")


def write_manifest_tsv(path, m):
    with open(path, "w") as fh:
        fh.write("class\tfeature\tgene\tcontig\tcoord\tdetail\tcarriers\n")
        for t in m["coding_variants"]:
            fh.write(f"coding\t{t['id']}\t{t['gene']}\t{t['contig']}\t{t['pos']}\t"
                     f"{t['ref']}>{t['alt']}|{t['consequence']}\t{','.join(t['carriers'])}\n")
        for t in m["promoter_variants"]:
            fh.write(f"promoter\t{t['id']}\t{t['gene']}\t{t['contig']}\t{t['pos']}\t"
                     f"{t['ref']}>{t['alt']}|{t['consequence']}\t{','.join(t['carriers'])}\n")
        s = m["str"]
        fh.write(f"str\tstr_{s['gene']}\t{s['gene']}\t{s['contig']}\t{s['repeat_region']}\t"
                 f"({s['motif']})n {s['normal_units']}->{s['expanded_units']}\t{','.join(s['carriers'])}\n")
        for t in m["cnv_events"]:
            fh.write(f"cnv\t{t['id']}\t{t['gene']}\t{D.GENES[t['gene']]['chrom']}\t-\t"
                     f"{t['state']} fold={t['fold']}\t{','.join(t['carriers'])}\n")


if __name__ == "__main__":
    main()

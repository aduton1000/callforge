#!/usr/bin/env python3
"""build_fixture.py — generate the CallForge simulation fixture (reference + reads + truth).

Consumes cohort_design.py and builds, under --outdir, a fully self-contained, hermetic
fixture for the 16-sample cohort:

  ref/genome.fa(.fai)        synthetic mini-genome; contigs named with Ensembl labels
                             (1,2,…,X — NO `chr` prefix). Each gene occupies SYNTHETIC
                             local coordinates on its real chromosome's contig.
  ref/targets.bed            GENE|class target BED (coding/promoter/str/cnv).
  ref/genes.gtf(.gz/.tbi)    one protein_coding transcript per CDS gene -> VEP (gtf mode)
                             computes consequences offline.
  ref/gene_metadata.tsv      via the REAL bin/make_gene_metadata.py (schema-correct).
  ref/str_catalog.json       curated ExpansionHunter catalog for the HTT (CAG)n locus.
  ref/somalier_sites.vcf.gz  autosomal + X sites for somalier (relatedness + sex).
  ref/truth.vcf.gz / .bed    GIAB-style truth for one control sample.
  resources/gnomad_sim.vcf.gz faked gnomAD-AFR AF at planted coding sites (rare -> burden
                             qualifiers) + one site per contig (genomic-scope FIT).
  fastq/<sid>_R{1,2}.fastq.gz per-sample paired reads carrying the planted variants, with
                             per-sample DEPTH modulation over CNV regions and an EXPANDED
                             (CAG)n repeat for STR carriers.
  truth_manifest.json / .tsv RESOLVED truth (exact coords/ref/alt/consequence + the
                             sample x feature matrix) — the spec the checker asserts against.
  sim.params.yaml            turnkey params for the cluster run (fill container_image +
                             slurm_partition).

Reads are emitted by deterministic pure-Python fragment tiling of two haplotypes per
sample (hapA=ref, hapB=variant) so heterozygous planted variants land at ~0.5 allele
balance. No external read simulator is required; only samtools/bgzip/tabix (for indices).

Runs anywhere samtools/bgzip/tabix exist (the cluster .sif, or a Mac with htslib). It
does NOT run the pipeline — that is the separate cluster step.
"""
import argparse
import gzip
import json
import os
import subprocess
import sys
import zlib

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cohort_design as D

# ------------------------------------------------------------------- tunables
READLEN = 150
FRAG = 320
BASE_DEPTH = 30                 # total het depth; each haplotype emits ~half
PROMOTER_LEN = 1600             # upstream region (<5kb -> VEP upstream_gene_variant)
CDS_CODONS = 120               # protein-coding exon length in codons (incl. ATG + stop)
STR_FLANK = 1200               # unique flank each side of the repeat (> EH window)
CNV_REGION = 1600              # CNV target length (depth-modulated)
LEAD = 300                     # unique lead/trail flank around each gene
GENE_GAP = 500                 # spacer between genes on a contig
N_X_SITES = 12                 # informative X sites for somalier sex inference

# Non-stop codons (mapping-unique filler for CDS). No TAA/TAG/TGA.
SAFE_CODONS = ["GCT", "GGA", "ACC", "AAC", "TTC", "GTT", "CAC", "ATC",
               "CCA", "TGG", "CTG", "GAC", "TAC", "AGC", "GTG", "CGT"]
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def run(cmd):
    subprocess.run(cmd, check=True, capture_output=True, text=True)


class Rng:
    """Tiny deterministic LCG -> stable across machines/Python versions (no hash salt)."""
    def __init__(self, seed):
        self.x = (seed & 0xFFFFFFFF) or 1

    def next(self):
        self.x = (1103515245 * self.x + 12345) & 0x7FFFFFFF
        return self.x

    def base(self):
        return "ACGT"[self.next() % 4]

    def seq(self, n):
        return "".join(self.base() for _ in range(n))

    def codon(self):
        return SAFE_CODONS[self.next() % len(SAFE_CODONS)]


def gseed(gene):
    return (D.SEED ^ zlib.crc32(gene.encode())) & 0xFFFFFFFF


# --------------------------------------------------------------- gene construction
def has_cds(gene):
    roles = D.GENES[gene]["roles"]
    return any(r in roles for r in ("coding", "coding_negative", "promoter", "burden"))


def coding_sites_for(gene):
    """Planted coding variants (showcase + burden) for this gene, with codon indices."""
    plants = [v for v in D.ALL_CODING if v["gene"] == gene]
    out, idx = [], 10
    for v in plants:
        out.append({"id": v["id"], "kind": v["kind"], "carriers": v["carriers"], "codon_idx": idx})
        idx += 7
    return out


def build_gene(gene):
    """Return a dict describing the gene's local sequence, intervals, and planted sites."""
    rng = Rng(gseed(gene))
    roles = D.GENES[gene]["roles"]
    seq = []
    pos = 0

    def emit(s):
        nonlocal pos
        seq.append(s)
        start = pos
        pos += len(s)
        return start

    info = {"gene": gene, "chrom": D.GENES[gene]["chrom"], "bed": [],
            "coding": [], "promoter": None, "str": None, "cnv": None, "cds": None}

    emit(rng.seq(LEAD))

    is_promoter = "promoter" in roles
    if "cnv" in roles:
        reg_start = emit(rng.seq(CNV_REGION))
        info["cnv"] = {"start": reg_start, "end": pos}
        info["bed"].append((reg_start, pos, "cnv"))
        emit(rng.seq(LEAD))

    elif "str" in roles:
        f1 = emit(rng.seq(STR_FLANK))
        s = D.STR_LOCUS
        rep_start = emit(s["motif"] * s["normal_units"])
        rep_end = pos
        emit(rng.seq(STR_FLANK))
        emit(rng.seq(LEAD))
        info["str"] = {"rep_start": rep_start, "rep_end": rep_end,
                       "motif": s["motif"], "normal_units": s["normal_units"],
                       "expanded_units": s["expanded_units"], "carriers": s["carriers"]}
        info["bed"].append((f1 + 50, pos - LEAD - 50, "str"))   # span flanks for anchored reads

    else:  # CDS gene (coding / coding_negative / promoter-with-transcript / burden)
        if is_promoter:
            prom_start = emit(rng.seq(PROMOTER_LEN))
            ptgt = prom_start + 800
            rb = seq_char(seq, ptgt)
            ab = {"A": "T", "T": "A", "C": "G", "G": "C"}[rb]
            info["promoter"] = {"pos": ptgt, "ref": rb, "alt": ab}
            # promoter target interval (gets reads so the SNV is callable)
            info["bed"].append((prom_start + 700, prom_start + 900, "promoter"))

        # CDS: ATG + middle (filler/planted) + stop. Planted codons resolved to exact bases.
        plants = {p["codon_idx"]: p for p in coding_sites_for(gene)}
        codons = ["ATG"]
        sites = []
        for i in range(1, CDS_CODONS - 1):
            if i in plants:
                p = plants[i]
                if p["kind"] == "stop":
                    cod, voff, ref, alt, csq, vtype = "CAA", 0, "C", "T", "stop_gained", "snv"
                elif p["kind"] == "missense":
                    cod, voff, ref, alt, csq, vtype = "GCT", 1, "C", "T", "missense_variant", "snv"
                else:  # frameshift = 1-bp deletion of the codon's first base
                    cod = rng.codon()
                    voff, ref, alt, csq, vtype = 0, cod[0], "", "frameshift_variant", "del1"
                codons.append(cod)
                sites.append({"id": p["id"], "carriers": p["carriers"], "cds_off": 3 * i + voff,
                              "ref": ref, "alt": alt, "consequence": csq, "vartype": vtype})
            else:
                codons.append(rng.codon())
        codons.append("TAA")
        cds_seq = "".join(codons)
        cds_start = emit(cds_seq)
        cds_end = pos
        emit(rng.seq(LEAD))
        info["cds"] = {"start": cds_start, "end": cds_end}
        info["bed"].append((cds_start, cds_end, "coding"))
        for s in sites:
            s["pos"] = cds_start + s["cds_off"]      # gene-local 0-based
            info["coding"].append(s)

    info["seq"] = "".join(seq)
    info["length"] = len(info["seq"])
    return info


def seq_char(seq_parts, pos):
    """char at absolute index `pos` across the list-of-strings being built."""
    acc = 0
    for part in seq_parts:
        if pos < acc + len(part):
            return part[pos - acc]
        acc += len(part)
    return "A"


# ----------------------------------------------------------------- read emission
def sample_edits(sample, gi):
    """List of (kind, pos, ref, alt) gene-local edits this sample carries for gene gi."""
    edits = []
    for s in gi["coding"]:
        if sample in s["carriers"]:
            edits.append((s["vartype"], s["pos"], s["ref"], s["alt"]))
    if gi["promoter"] and sample in [c for v in D.PROMOTER_VARIANTS
                                     if v["gene"] == gi["gene"] for c in v["carriers"]]:
        p = gi["promoter"]
        edits.append(("snv", p["pos"], p["ref"], p["alt"]))
    return edits


def haplotype_b(sample, gi):
    """Build the variant haplotype sequence for (sample, gene). STR handled separately."""
    seq = list(gi["seq"])
    edits = sample_edits(sample, gi)
    # SNVs first (no shift), then single-base deletions right-to-left.
    for kind, pos, ref, alt in edits:
        if kind == "snv" and pos < len(seq):
            seq[pos] = alt
    for kind, pos, ref, alt in sorted([e for e in edits if e[0] == "del1"],
                                      key=lambda e: -e[1]):
        if pos < len(seq):
            del seq[pos]
    s = "".join(seq)
    # STR expansion (carriers are heterozygous: hapB carries the expanded allele).
    if gi["str"] and sample in gi["str"]["carriers"]:
        st = gi["str"]
        exp = st["motif"] * st["expanded_units"]
        s = s[:st["rep_start"]] + exp + s[st["rep_end"]:]
    return s


def tile_reads(o1, o2, seq, lo, hi, depth, tag):
    """Tile paired reads from seq[lo:hi] at ~depth; write to gz handles."""
    per = max(1, depth)
    step = max(1, (READLEN * 2) // per)
    q = "I" * READLEN
    rid = 0
    end = min(len(seq), hi) - FRAG
    for start in range(max(0, lo - 30), max(0, end), step):
        frag = seq[start:start + FRAG]
        if len(frag) < FRAG:
            continue
        m1 = frag[:READLEN]
        m2 = revcomp(frag[-READLEN:])
        o1.write(f"@{tag}_{rid}/1\n{m1}\n+\n{q}\n")
        o2.write(f"@{tag}_{rid}/2\n{m2}\n+\n{q}\n")
        rid += 1
    return rid


def cnv_fold(sample, gene):
    for e in D.CNV_EVENTS:
        if e["gene"] == gene and sample in e["carriers"]:
            return e["fold"]
    return 1.0


def emit_sample(sample, sex, genes, xseq, x_sites, fqdir):
    r1p = os.path.join(fqdir, f"{sample}_R1.fastq.gz")
    r2p = os.path.join(fqdir, f"{sample}_R2.fastq.gz")
    with gzip.open(r1p, "wt") as o1, gzip.open(r2p, "wt") as o2:
        # --- X off-target reads -> somalier sex (F: full depth + het; M: half depth hom-ref)
        xdepth = BASE_DEPTH if sex == "F" else max(1, BASE_DEPTH // 2)
        step = max(1, (READLEN * 2) // xdepth)
        q = "I" * READLEN
        xc = 0
        for start in range(0, len(xseq) - FRAG, step):
            frag = xseq[start:start + FRAG]
            if sex == "F" and xc % 2 == 1:               # alt haplotype -> het at X sites
                fl = list(frag)
                for (p, rb, ab) in x_sites:
                    if start <= p < start + FRAG:
                        fl[p - start] = ab
                frag = "".join(fl)
            o1.write(f"@{sample}_X{xc}/1\n{frag[:READLEN]}\n+\n{q}\n")
            o2.write(f"@{sample}_X{xc}/2\n{revcomp(frag[-READLEN:])}\n+\n{q}\n")
            xc += 1
        # --- per-gene reads from two haplotypes (het balance), CNV-depth-modulated
        for gi in genes:
            fold = cnv_fold(sample, gi["gene"])
            hap_depth = max(1, int(round(BASE_DEPTH * fold / 2)))
            hapA = gi["seq"]
            hapB = haplotype_b(sample, gi)
            for (s, e, klass) in gi["bed"]:
                tile_reads(o1, o2, hapA, s, e, hap_depth, f"{sample}_{gi['gene']}_{klass}_A")
                # hapB may be longer/shorter (indel/STR) -> tile its full equivalent span
                hb_hi = e + (len(hapB) - len(hapA)) if e >= (gi["str"]["rep_start"] if gi["str"] else 10**9) else e
                tile_reads(o1, o2, hapB, s, max(e, hb_hi), hap_depth, f"{sample}_{gi['gene']}_{klass}_B")
    return r1p, r2p


# --------------------------------------------------------------------- writers
def write_genome(path, contigs):
    with open(path, "w") as fh:
        for name, seq in contigs:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
    run(["samtools", "faidx", path])


def write_bed(path, rows):
    with open(path, "w") as fh:
        for (c, s, e, name) in sorted(rows):
            fh.write(f"{c}\t{s}\t{e}\t{name}\t0\t+\n")


def write_gtf(path, gtf_genes):
    with open(path, "w") as fh:
        for (c, g, cs, ce) in gtf_genes:
            at = (f'gene_id "{g}"; transcript_id "{g}_t1"; gene_name "{g}"; '
                  f'transcript_biotype "protein_coding";')
            fh.write(f'{c}\tcallforge_sim\tgene\t{cs+1}\t{ce}\t.\t+\t.\t'
                     f'gene_id "{g}"; gene_name "{g}"; gene_biotype "protein_coding";\n')
            fh.write(f"{c}\tcallforge_sim\ttranscript\t{cs+1}\t{ce}\t.\t+\t.\t{at}\n")
            fh.write(f'{c}\tcallforge_sim\texon\t{cs+1}\t{ce}\t.\t+\t.\t{at} exon_number "1";\n')
            fh.write(f'{c}\tcallforge_sim\tCDS\t{cs+1}\t{ce}\t.\t+\t0\t{at} exon_number "1";\n')
    run(["bash", "-c", f"sort -k1,1 -k4,4n '{path}' | bgzip > '{path}.gz'"])
    run(["tabix", "-f", "-p", "gff", path + ".gz"])


def write_vcf_bgzip(path, contigs, header_extra, records):
    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        for name, seq in contigs:
            fh.write(f"##contig=<ID={name},length={len(seq)}>\n")
        for h in header_extra:
            fh.write(h + "\n")
        fh.write(records)
    run(["bgzip", "-f", path])
    run(["tabix", "-p", "vcf", path + ".gz"])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--make-gene-metadata", default=None,
                    help="path to bin/make_gene_metadata.py (default: ../../bin relative to this script)")
    a = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    mgm = a.make_gene_metadata or os.path.normpath(os.path.join(here, "..", "..", "bin", "make_gene_metadata.py"))
    refdir = os.path.join(a.outdir, "ref")
    fqdir = os.path.join(a.outdir, "fastq")
    resdir = os.path.join(a.outdir, "resources")
    for d in (refdir, fqdir, resdir):
        os.makedirs(d, exist_ok=True)

    # ---- build each gene, then lay genes out per contig (Ensembl-named) ----
    genes = {g: build_gene(g) for g in D.GENES}
    by_chrom = {}
    for g, gi in genes.items():
        by_chrom.setdefault(gi["chrom"], []).append(g)

    contigs = []                 # (name, seq) in stable order
    bed_rows = []                # (contig, start, end, "GENE|class")
    gtf_genes = []               # (contig, gene, cds_start_abs, cds_end_abs)
    str_catalog = []
    str_truth = None
    cnv_truth = []
    coding_truth = []
    promoter_truth = []
    gnomad_recs = {}             # (contig,pos1)->(ref,alt,af)
    somalier_sites = {}          # (contig,pos1)->(ref,alt)

    for chrom in sorted(by_chrom, key=lambda c: (len(c), c)):
        seq_parts = []
        offset = 0
        for g in sorted(by_chrom[chrom]):
            gi = genes[g]
            gi["offset"] = offset
            gi["contig"] = chrom
            seq_parts.append(gi["seq"])
            # BED (absolute contig coords)
            for (s, e, klass) in gi["bed"]:
                bed_rows.append((chrom, offset + s, offset + e, f"{g}|{klass}"))
            # GTF
            if gi["cds"]:
                gtf_genes.append((chrom, g, offset + gi["cds"]["start"], offset + gi["cds"]["end"]))
            # coding truth + gnomAD + somalier site at each coding variant
            for s in gi["coding"]:
                pos1 = offset + s["pos"] + 1
                coding_truth.append({"id": s["id"], "gene": g, "contig": chrom, "pos": pos1,
                                     "ref": s["ref"], "alt": s["alt"], "vartype": s["vartype"],
                                     "consequence": s["consequence"], "zygosity": "het",
                                     "carriers": s["carriers"]})
                if s["vartype"] == "snv":
                    gnomad_recs[(chrom, pos1)] = (s["ref"], s["alt"], "0.002")
                    somalier_sites[(chrom, pos1)] = (s["ref"], s["alt"])
            # promoter truth
            if gi["promoter"]:
                pv = next(v for v in D.PROMOTER_VARIANTS if v["gene"] == g)
                pos1 = offset + gi["promoter"]["pos"] + 1
                promoter_truth.append({"id": pv["id"], "gene": g, "contig": chrom, "pos": pos1,
                                       "ref": gi["promoter"]["ref"], "alt": gi["promoter"]["alt"],
                                       "consequence": "upstream_gene_variant",
                                       "zygosity": "het", "carriers": pv["carriers"]})
            # STR catalog + truth
            if gi["str"]:
                st = gi["str"]
                rs, re_ = offset + st["rep_start"] + 1, offset + st["rep_end"]
                region = f"{chrom}:{rs}-{re_}"
                locus_id = f"{g}_{chrom}_{rs}_{re_}"
                str_catalog.append({"LocusId": locus_id, "LocusStructure": f"({st['motif']})*",
                                    "ReferenceRegion": region, "VariantType": "Repeat"})
                str_truth = {"gene": g, "contig": chrom, "repeat_region": region,
                             "locus_id": locus_id, "motif": st["motif"],
                             "normal_units": st["normal_units"],
                             "expanded_units": st["expanded_units"], "carriers": st["carriers"]}
            # CNV truth
            if gi["cnv"]:
                ev = next(e for e in D.CNV_EVENTS if e["gene"] == g)
                cnv_truth.append({"id": ev["id"], "gene": g, "contig": chrom,
                                  "start": offset + gi["cnv"]["start"], "end": offset + gi["cnv"]["end"],
                                  "state": ev["state"], "fold": ev["fold"], "carriers": ev["carriers"]})
            offset += gi["length"] + GENE_GAP
            seq_parts.append("A" * GENE_GAP)            # inter-gene spacer
        contigs.append((chrom, "".join(seq_parts)))

    # ---- X off-target contig for somalier sex ----
    xrng = Rng(gseed("chrX_offtarget"))
    xseq = xrng.seq(20000)
    x_sites = []
    for i in range(N_X_SITES):
        p = 1000 + i * 1400
        if p < len(xseq):
            rb = xseq[p]
            ab = {"A": "G", "G": "A", "C": "T", "T": "C"}[rb]
            x_sites.append((p, rb, ab))
            somalier_sites[("X", p + 1)] = (rb, ab)
    contigs.append(("X", xseq))

    # one gnomAD + somalier site per contig (genomic-scope FIT for discovery)
    for (cname, cseq) in contigs:
        if cname == "X":
            continue
        if not any(c == cname for (c, _) in gnomad_recs):
            rb = cseq[200]
            ab = {"A": "G", "G": "A", "C": "T", "T": "C"}.get(rb, "G")
            gnomad_recs[(cname, 201)] = (rb, ab, "0.40")

    # ---- write reference, BED, GTF ----
    genome_fa = os.path.join(refdir, "genome.fa")
    write_genome(genome_fa, contigs)
    targets_bed = os.path.join(refdir, "targets.bed")
    write_bed(targets_bed, bed_rows)
    gtf_path = os.path.join(refdir, "genes.gtf")
    write_gtf(gtf_path, gtf_genes)

    # ---- STR catalog (curated) ----
    catalog_path = os.path.join(refdir, "str_catalog.json")
    with open(catalog_path, "w") as fh:
        json.dump(str_catalog, fh, indent=2)

    # ---- str-loci + burden-groups TSVs, then the REAL gene-metadata builder ----
    strloci_tsv = os.path.join(refdir, "_str_loci.tsv")
    with open(strloci_tsv, "w") as fh:
        if str_truth:
            fh.write(f"{str_truth['gene']}\t{str_truth['repeat_region']}\n")
    burden_tsv = os.path.join(refdir, "_burden_groups.tsv")
    with open(burden_tsv, "w") as fh:
        fh.write(f"{D.BURDEN_GENE}\t{D.BURDEN_GENE}\n")
    meta_tsv = os.path.join(refdir, "gene_metadata.tsv")
    run([sys.executable, mgm, "--bed", targets_bed,
         "--cnv-genes", ",".join(D.CNV_TARGET_GENES),
         "--str-loci", strloci_tsv, "--burden-groups", burden_tsv,
         "--out-tsv", meta_tsv, "--out-json", os.path.join(refdir, "gene_metadata.json"),
         "--overwrite"])

    # ---- faked gnomAD-AFR resource ----
    gnv = os.path.join(resdir, "gnomad_sim.vcf")
    recs = "##INFO=<ID=AF_afr,Number=A,Type=Float,Description=\"AFR allele frequency\">\n##reference=GRCh38\n"
    recs += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    cidx = {c: i for i, (c, _) in enumerate(contigs)}
    for (c, p) in sorted(gnomad_recs, key=lambda k: (cidx.get(k[0], 99), k[1])):
        rb, ab, af = gnomad_recs[(c, p)]
        recs += f"{c}\t{p}\t.\t{rb}\t{ab}\t.\t.\tAF_afr={af}\n"
    write_vcf_bgzip(gnv, contigs, [], recs)

    # ---- somalier sites ----
    site_recs = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    for (c, p) in sorted(somalier_sites, key=lambda k: (cidx.get(k[0], 99), k[1])):
        rb, ab = somalier_sites[(c, p)]
        site_recs += f"{c}\t{p}\t{c}_{p}\t{rb}\t{ab}\t.\t.\t.\n"
    write_vcf_bgzip(os.path.join(refdir, "somalier_sites.vcf"), contigs, [], site_recs)

    # ---- per-sample reads ----
    gene_list = [genes[g] for g in D.GENES]
    sex_by_sample = {}
    for s in D.SAMPLES:
        emit_sample(s["sample_id"], s["sex"], gene_list, xseq, x_sites, fqdir)
        sex_by_sample[s["sample_id"]] = s["sex"]

    # ---- GIAB-style truth for one control (S09) over its SNV calls ----
    giab_id = "S09"
    truth_lines = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    truth_lines += f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{giab_id}\n"
    gt_rows = sorted([t for t in coding_truth
                      if t["vartype"] == "snv" and giab_id in t["carriers"]],
                     key=lambda t: (cidx.get(t["contig"], 99), t["pos"]))
    for t in gt_rows:
        truth_lines += f"{t['contig']}\t{t['pos']}\t.\t{t['ref']}\t{t['alt']}\t100\tPASS\t.\tGT\t0/1\n"
    write_vcf_bgzip(os.path.join(refdir, "truth.vcf"), contigs, [], truth_lines)
    with open(os.path.join(refdir, "truth.bed"), "w") as fh:
        for (c, s, e, _) in sorted(bed_rows):
            fh.write(f"{c}\t{s}\t{e}\n")

    # ---- resolved truth manifest ----
    manifest = {
        "seed": D.SEED,
        "genome_fasta": os.path.abspath(genome_fa),
        "target_bed": os.path.abspath(targets_bed),
        "gtf": os.path.abspath(gtf_path + ".gz"),
        "gene_metadata": os.path.abspath(meta_tsv),
        "str_catalog": os.path.abspath(catalog_path),
        "somalier_sites": os.path.abspath(os.path.join(refdir, "somalier_sites.vcf.gz")),
        "gnomad_vcf": os.path.abspath(gnv + ".gz"),
        "contigs": {name: len(seq) for name, seq in contigs},
        "samples": D.SAMPLES,
        "cases": D.CASES, "controls": D.CONTROLS,
        "coding_variants": coding_truth,
        "promoter_variants": promoter_truth,
        "str": str_truth,
        "cnv_events": cnv_truth,
        "burden": {"gene": D.BURDEN_GENE,
                   "qualifying_feature_ids": [v["id"] for v in D.BURDEN_VARIANTS],
                   "case_carriers": sorted({c for v in D.BURDEN_VARIANTS
                                            for c in v["carriers"] if c in D.CASES}),
                   "control_carriers": sorted({c for v in D.BURDEN_VARIANTS
                                               for c in v["carriers"] if c in D.CONTROLS})},
        "cohortqc": {"sex_by_sample": sex_by_sample},
        "giab": {"control_id": giab_id,
                 "truth_vcf": os.path.abspath(os.path.join(refdir, "truth.vcf.gz")),
                 "truth_bed": os.path.abspath(os.path.join(refdir, "truth.bed"))},
        "thresholds": D.THRESHOLDS,
    }
    with open(os.path.join(a.outdir, "truth_manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)
    write_manifest_tsv(os.path.join(a.outdir, "truth_manifest.tsv"), manifest)

    # ---- turnkey params for the cluster run ----
    write_params_yaml(os.path.join(a.outdir, "sim.params.yaml"), manifest, a.outdir)

    print(f"[build_fixture] fixture ready in {a.outdir}: "
          f"{len(contigs)} contigs, {len(bed_rows)} targets, "
          f"{len(coding_truth)} coding + {len(promoter_truth)} promoter + "
          f"{len(cnv_truth)} CNV planted, STR={str_truth['gene'] if str_truth else None}, "
          f"{len(D.SAMPLES)} samples.")


def write_manifest_tsv(path, m):
    with open(path, "w") as fh:
        fh.write("class\tfeature\tgene\tcontig\tcoord\tdetail\tcarriers\n")
        for t in m["coding_variants"]:
            fh.write(f"coding\t{t['id']}\t{t['gene']}\t{t['contig']}\t{t['pos']}\t"
                     f"{t['ref']}>{t['alt'] or '-'}|{t['consequence']}\t{','.join(t['carriers'])}\n")
        for t in m["promoter_variants"]:
            fh.write(f"promoter\t{t['id']}\t{t['gene']}\t{t['contig']}\t{t['pos']}\t"
                     f"{t['ref']}>{t['alt']}|{t['consequence']}\t{','.join(t['carriers'])}\n")
        if m["str"]:
            s = m["str"]
            fh.write(f"str\tstr_{s['gene']}\t{s['gene']}\t{s['contig']}\t{s['repeat_region']}\t"
                     f"({s['motif']})n {s['normal_units']}->{s['expanded_units']}\t{','.join(s['carriers'])}\n")
        for t in m["cnv_events"]:
            fh.write(f"cnv\t{t['id']}\t{t['gene']}\t{t['contig']}\t{t['start']}-{t['end']}\t"
                     f"{t['state']} fold={t['fold']}\t{','.join(t['carriers'])}\n")


def write_params_yaml(path, m, outdir):
    lines = [
        "# CallForge simulation-cohort params (generated by build_fixture.py).",
        "# Fill container_image + slurm_partition for the cluster, then:",
        "#   callforge-run -params-file sim.params.yaml -profile hpc_slurm,singularity",
        f"genome_build:    GRCh38",
        f"genome_fasta:    {m['genome_fasta']}",
        f"target_bed:      {m['target_bed']}",
        f"gene_metadata:   {m['gene_metadata']}",
        f"str_catalog:     {m['str_catalog']}",
        f"somalier_sites:  {m['somalier_sites']}",
        f"gnomad_vcf:      {m['gnomad_vcf']}",
        f"gnomad_af_field: AF_afr",
        f"gtf:             {m['gtf']}",
        f"vep_mode:        gtf",
        f"paralog_genes:   ''",
        f"cnv_caller:      cnvkit",
        f"burden_engine:   collapse   # tiny cohort -> collapse/SKAT; regenie needs more samples",
        f"burden_af_max:   0.01",
        f"giab_control_id: {m['giab']['control_id']}",
        f"giab_truth_vcf:  {m['giab']['truth_vcf']}",
        f"giab_truth_bed:  {m['giab']['truth_bed']}",
        f"outdir:          {os.path.abspath(os.path.join(outdir, 'results'))}",
        "container_image: /PATH/TO/callforge.sif   # FILL ME",
        "slurm_partition: FILL_ME                  # FILL ME",
        "bind_paths:      /hpc                      # adjust to where this fixture lives",
    ]
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""make_test_reference.py — build the CallForge `test` profile fixture.

From the local GRCh38 no-alt FASTA + the CaptureForge target BED, builds a tiny,
self-contained, hermetic test fixture:
  * test_genome.fa(.fai)   each selected gene's locus ± flank becomes a short
                           contig (named test_<GENE>); no alt/decoy contigs.
  * test_targets.bed       the genes' intervals remapped to local contig coords
                           (same naming style as the reference -> invariant holds).
  * synthetic paired FASTQs simulated by tiling each target at ~depthx (150 bp).
  * test_gene_metadata.tsv via ingest_captureforge.py on the subset BED.
  * test_str_catalog.json  minimal ExpansionHunter catalog for STR-class targets.
  * test_truth.vcf.gz(.tbi)+ test_truth.bed  tiny GIAB-style truth for the control.
  * resources/             a fake VEP cache + tiny gnomAD VCF so discovery FINDS
                           something (exercises Stage 0 discovery offline).
  * test_samplesheet.csv   CTRL (control/GIAB) + S1 (case) + S2 (control).

Uses `samtools faidx` (host) for extraction — no pysam dependency. Deterministic
(no RNG): reads are a fixed tiling, so re-runs are byte-stable.
"""
import argparse, gzip, json, os, subprocess, sys
from collections import defaultdict

FLANK = 400
MAX_INTERVALS_PER_GENE = 4
READLEN = 150
FRAG = 300


def run(cmd):
    return subprocess.run(cmd, check=True, capture_output=True, text=True).stdout


def revcomp(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def parse_bed(path):
    g = defaultdict(lambda: defaultdict(list))  # gene -> class -> [(chrom,start,end)]
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, s, e = f[0], int(f[1]), int(f[2])
            gene, _, klass = (f[3] if len(f) > 3 else "NA|coding").partition("|")
            g[gene][klass or "coding"].append((chrom, s, e))
    return g


def pick_genes(genes):
    """Pick a small, class-diverse set: one STR, one CNV, one paralog, one coding."""
    chosen = []
    def first_with(klass):
        for gene in sorted(genes):
            if klass in genes[gene] and gene not in chosen:
                return gene
        return None
    for klass in ("str", "cnv"):
        gname = first_with(klass)
        if gname:
            chosen.append(gname)
    for p in ("HP", "CR1", "CFH", "CD209", "CASP1"):
        if p in genes and p not in chosen:
            chosen.append(p); break
    for gene in sorted(genes):
        if gene not in chosen:
            chosen.append(gene); break
    return chosen[:4]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", required=True)
    ap.add_argument("--bed", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--ingest", required=True, help="path to ingest_captureforge.py")
    ap.add_argument("--depth", type=int, default=30)
    a = ap.parse_args()

    refdir = os.path.join(a.outdir, "ref")
    fqdir = os.path.join(a.outdir, "fastq")
    resdir = os.path.join(a.outdir, "resources")
    for d in (refdir, fqdir, resdir):
        os.makedirs(d, exist_ok=True)

    genes = parse_bed(a.bed)
    chosen = pick_genes(genes)
    sys.stderr.write(f"[make_test] selected genes: {chosen}\n")

    contigs = {}          # contig_name -> sequence
    bed_rows = []         # (contig, start, end, name)
    str_rows = []         # (contig, start, end, gene)
    for gene in chosen:
        # gather intervals across classes, cap count, keep on a single chrom
        ivs = []
        for klass, locs in genes[gene].items():
            for (c, s, e) in locs:
                ivs.append((c, s, e, klass))
        ivs.sort()
        chrom0 = ivs[0][0]
        ivs = [iv for iv in ivs if iv[0] == chrom0][:MAX_INTERVALS_PER_GENE]
        win_s = max(0, min(s for _, s, _, _ in ivs) - FLANK)
        win_e = max(e for _, _, e, _ in ivs) + FLANK
        contig = f"test_{gene}"
        seq = run(["samtools", "faidx", a.genome, f"{chrom0}:{win_s+1}-{win_e}"])
        seq = "".join(seq.splitlines()[1:]).upper()
        if not seq or set(seq) <= {"N"}:
            sys.stderr.write(f"[make_test] WARN: empty/N-only locus for {gene}; skipping\n")
            continue
        contigs[contig] = seq
        for (c, s, e, klass) in ivs:
            ls, le = s - win_s, e - win_s
            le = min(le, len(seq))
            if le <= ls:
                continue
            bed_rows.append((contig, ls, le, f"{gene}|{klass}"))
            if klass == "str":
                str_rows.append((contig, ls, le, gene))

    if not contigs:
        sys.exit("[make_test] ERROR: no usable loci extracted")

    # ---- synthetic STR locus (gene HMOX1, (GT)n) so ExpansionHunter has a real
    #      repeat to genotype; flanks are revcomp'd real sequence (unique mapping).
    #      Flanks must exceed EH's ~1000 bp extension window each side. ----
    src = max(contigs.values(), key=len)
    flank = 1100
    if len(src) < 2 * flank:
        src = (src * (2 * flank // max(1, len(src)) + 1))       # tile up if needed
    n_gt = 20
    str_seq = revcomp(src[:flank]) + ("GT" * n_gt) + revcomp(src[flank:2 * flank])
    contigs["test_HMOX1"] = str_seq
    rep_s, rep_e = flank, flank + 2 * n_gt                       # precise repeat coords
    bed_rows.append(("test_HMOX1", 50, len(str_seq) - 50, "HMOX1|str"))  # wide target -> reads
    str_rows.append(("test_HMOX1", rep_s, rep_e, "HMOX1"))               # catalog at the repeat

    # ---- write mini genome + faidx ----
    genome_fa = os.path.join(refdir, "test_genome.fa")
    with open(genome_fa, "w") as fh:
        for name, seq in contigs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + "\n")
    run(["samtools", "faidx", genome_fa])

    # ---- target BED (local coords) ----
    test_bed = os.path.join(refdir, "test_targets.bed")
    with open(test_bed, "w") as fh:
        for (c, s, e, name) in sorted(bed_rows):
            fh.write(f"{c}\t{s}\t{e}\t{name}\t0\t+\n")

    # ---- CaptureForge-style metrics.json: CNV callability so the test exercises
    #      callability labelling (CFH/CR1 breakpoint-blind; others depth-callable) ----
    cnv_genes = sorted({n.split("|")[0] for (c, s, e, n) in bed_rows if n.split("|")[1] == "cnv"})
    metrics = {"cnv_bin_spacing": {}, "qc_gates": {"low_coverage_genes": []}}
    for g in cnv_genes:
        callable_ = g not in ("CFH", "CR1")
        metrics["cnv_bin_spacing"][g] = {
            "n_bins": 200, "median_gap": 52,
            "max_gap": 18000 if not callable_ else 700,
            "cv_gap": 4.4 if not callable_ else 1.1, "depth_callable": callable_}
        if not callable_:
            metrics["qc_gates"]["low_coverage_genes"].append({"gene": g, "raw_cov": 0.42})
    test_metrics = os.path.join(refdir, "test_metrics.json")
    with open(test_metrics, "w") as fh:
        json.dump(metrics, fh, indent=2)

    # ---- gene metadata via the real ingest script (BED + metrics) ----
    subprocess.run([sys.executable, a.ingest, "--bed", test_bed, "--metrics", test_metrics,
                    "--paralog-genes", "HP,CR1,CFH,CD209,CASP1",
                    "--out-tsv", os.path.join(refdir, "test_gene_metadata.tsv"),
                    "--out-json", os.path.join(refdir, "test_gene_metadata.json")], check=True)

    # ---- ExpansionHunter catalog (curated) ----
    # Only the synthetic HMOX1 (GT)n is a genuine repeat locus; other str-class
    # targets are not EH loci (the real run supplies a curated catalog for
    # PIEZO1 E756del / HMOX1 (GT)n / SLC11A1 (GT)n via --str_catalog).
    catalog = [{"LocusId": g, "LocusStructure": "(GT)*", "ReferenceRegion": f"{c}:{s}-{e}",
                "VariantType": "Repeat"} for (c, s, e, g) in str_rows if c == "test_HMOX1"]
    with open(os.path.join(refdir, "test_str_catalog.json"), "w") as fh:
        json.dump(catalog, fh, indent=2)

    # ---- synthetic paired reads (deterministic tiling) ----
    # cnv_del: set of contigs to render at HALF coverage for this sample, so CNVkit
    # calls a deletion there (used to force a CFH breakpoint-blind CNV call).
    def simulate(sample, spike=None, depth=None, cnv_del=None):
        depth = depth or a.depth
        cnv_del = cnv_del or set()
        r1p = os.path.join(fqdir, f"{sample}_R1.fastq.gz")
        r2p = os.path.join(fqdir, f"{sample}_R2.fastq.gz")
        spike = spike or {}
        with gzip.open(r1p, "wt") as o1, gzip.open(r2p, "wt") as o2:
            rid = 0
            for (c, s, e, name) in bed_rows:
                seq = list(contigs[c])
                for pos, alt in spike.get(c, []):
                    if s <= pos < e and pos < len(seq):
                        seq[pos] = alt
                seq = "".join(seq)
                eff_depth = max(1, depth // 2) if c in cnv_del else depth
                step = max(1, READLEN // max(1, eff_depth) * 2)
                for start in range(max(0, s - 50), min(len(seq), e + 50) - FRAG, step):
                    frag = seq[start:start + FRAG]
                    if len(frag) < FRAG:
                        continue
                    m1 = frag[:READLEN]
                    m2 = revcomp(frag[-READLEN:])
                    q = "I" * READLEN
                    o1.write(f"@{sample}_{rid}/1\n{m1}\n+\n{q}\n")
                    o2.write(f"@{sample}_{rid}/2\n{m2}\n+\n{q}\n")
                    rid += 1
        return r1p, r2p

    # spiked SNVs (carried by CTRL + S1): one in a non-paralog gene (truth) and one
    # in a CFH paralog region so paralog-aware flagging has a variant to label.
    def make_spike(contig, pos):
        rb = contigs[contig][pos]
        ab = {"A": "G", "G": "A", "C": "T", "T": "C", "N": "A"}.get(rb, "A")
        return (contig, pos, rb, ab)
    spikes = [make_spike(bed_rows[0][0], bed_rows[0][1] + 20)]
    cfh_rows = [(c, s, e) for (c, s, e, n) in bed_rows if c == "test_CFH"]
    if cfh_rows:
        c, s, e = max(cfh_rows, key=lambda r: r[2] - r[1])     # longest CFH row (gets reads)
        spikes.append(make_spike(c, s + 50))
    spike_map = {}
    for (c, p, rb, ab) in spikes:
        spike_map.setdefault(c, []).append((p, ab))

    # (sample, spike, depth, cnv_del) — LOWQC is deliberately under-covered to
    # exercise quarantine; S1 carries a half-coverage CFH deletion so a CFH
    # (breakpoint-blind) CNV call appears and its low-confidence label is proven.
    samples = [("CTRL", spike_map, 30, set()), ("S1", spike_map, 30, {"test_CFH"}),
               ("S2", {}, 30, set()), ("LOWQC", {}, 2, set())]
    sheet = os.path.join(refdir, "..", "test_samplesheet.csv")
    sheet = os.path.normpath(sheet)
    with open(sheet, "w") as fh:
        fh.write("sample_id,fastq_1,fastq_2,sex,phenotype,covariate_age,batch\n")
        meta = {"CTRL": ("F", "control", "30", "b1"),
                "S1":   ("M", "case", "45", "b1"),
                "S2":   ("F", "control", "52", "b2"),
                "LOWQC":("F", "control", "39", "b2")}
        for sample, sp, dp, cd in samples:
            r1, r2 = simulate(sample, sp, dp, cd)
            sx, ph, ag, bt = meta[sample]
            fh.write(f"{sample},{os.path.abspath(r1)},{os.path.abspath(r2)},{sx},{ph},{ag},{bt}\n")

    # ---- GIAB-style truth (control spiked SNV) ----
    truth_vcf = os.path.join(refdir, "test_truth.vcf")
    with open(truth_vcf, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        for name, seq in contigs.items():
            fh.write(f"##contig=<ID={name},length={len(seq)}>\n")
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCTRL\n")
        for (c, p, rb, ab) in sorted(spikes):
            fh.write(f"{c}\t{p+1}\t.\t{rb}\t{ab}\t100\tPASS\t.\tGT\t0/1\n")
    run(["bgzip", "-f", truth_vcf])
    run(["tabix", "-p", "vcf", truth_vcf + ".gz"])
    with open(os.path.join(refdir, "test_truth.bed"), "w") as fh:
        for (c, s, e, _) in sorted(bed_rows):
            fh.write(f"{c}\t{s}\t{e}\n")

    # ---- fake resource fixture for offline discovery ----
    vep = os.path.join(resdir, "homo_sapiens", "110_GRCh38")
    os.makedirs(vep, exist_ok=True)
    with open(os.path.join(vep, "info.txt"), "w") as fh:
        fh.write("# fake VEP cache fixture for CallForge test discovery\nassembly\tGRCh38\n")
    gnv = os.path.join(resdir, "gnomad_test.vcf")
    # One record per target contig so the toy gnomAD is genomic-scope FIT (spans
    # the whole test panel) — exercises the scope gate's pass path.
    first_iv = {}
    for (c, s, e, _) in bed_rows:
        first_iv.setdefault(c, (s, e))
    with open(gnv, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n##INFO=<ID=AF_afr,Number=A,Type=Float,Description=\"AFR AF\">\n")
        for name, seq in contigs.items():
            fh.write(f"##contig=<ID={name},length={len(seq)}>\n")
        fh.write("##reference=GRCh38\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for c in contigs:
            pos = first_iv.get(c, (10, 20))[0] + 5
            rb = contigs[c][pos] if pos < len(contigs[c]) else "A"
            ab = {"A": "G", "G": "A", "C": "T", "T": "C", "N": "A"}.get(rb, "A")
            fh.write(f"{c}\t{pos+1}\t.\t{rb}\t{ab}\t.\t.\tAF_afr=0.01\n")
    run(["bgzip", "-f", gnv]); run(["tabix", "-p", "vcf", gnv + ".gz"])

    print(f"[make_test] fixture ready in {a.outdir} "
          f"({len(contigs)} contigs, {len(bed_rows)} targets, samples=CTRL,S1,S2)")


if __name__ == "__main__":
    main()

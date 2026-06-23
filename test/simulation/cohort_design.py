#!/usr/bin/env python3
"""cohort_design.py — SINGLE SOURCE OF TRUTH for the CallForge REAL-reference sim cohort.

Declares (with NO I/O) the 16-sample cohort used to validate CallForge against the REAL
GRCh38 no-alt reference + REAL annotation databases on the cluster: the generic gene set
with REAL GRCh38 coordinates, the per-class planted-variant matrix (which sample carries
what), the case/control + covariate design, the named checker thresholds, and the
per-stage real-reference dependency map.

Reads are simulated FROM the real reference at REAL coordinates (see simulate_reads.py /
dwgsim) and the pipeline runs against the REAL staged DBs — NOTHING is faked. Exact
planted positions/ref/alt and the verified consequence are RESOLVED on the cluster from
the real FASTA + gene model (resolve_variants.py) and written into truth_manifest.json;
this module carries only the cohort DESIGN and reviewable real-coordinate WINDOWS.

Scripts that import this:
  * simulate_cohort.py   -> reviewable plan + sample sheet + manifest skeleton (host, stdlib)
  * resolve_variants.py  -> real exonic placement + consequence verification (cluster)
  * simulate_reads.py    -> dwgsim reads from real reference (cluster)
  * check_against_truth.py -> staged PASS/FAIL vs the resolved manifest

stdlib-only so it loads without the conda env active.
"""

# Master seed — every deterministic choice (carrier assignment, dwgsim -z) derives from it.
SEED = 20260623
GENOME_BUILD = "GRCh38"

# --------------------------------------------------------------------------- genes
# Generic, well-characterized, NON-panel genes (no malaria-panel gene). Coordinates are
# REAL GRCh38 (Ensembl no-alt naming: "17","4",… — NO `chr` prefix). The `window` is a
# padded gene-locus span used to BOUND read simulation and as a BED fallback; the
# resolver treats the staged gene model (GTF) as AUTHORITATIVE for exon/CDS positions and
# validates these against it. `tx` is the MANE Select / canonical transcript (intent).
#   roles: which class(es) the gene carries.
GENES = {
    "TP53":  {"chrom": "17", "strand": "-", "window": (7668421, 7687550),  "tx": "ENST00000269305", "roles": ["coding"]},
    "BRCA1": {"chrom": "17", "strand": "-", "window": (43044295, 43125483),"tx": "ENST00000357654", "roles": ["coding"]},
    "NF1":   {"chrom": "17", "strand": "+", "window": (31094927, 31377677),"tx": "ENST00000358273", "roles": ["coding_negative"]},
    "EGFR":  {"chrom": "7",  "strand": "+", "window": (55019017, 55211628),"tx": "ENST00000275493", "roles": ["coding"]},
    "PTEN":  {"chrom": "10", "strand": "+", "window": (87863625, 87971930),"tx": "ENST00000371953", "roles": ["coding"]},
    "KRAS":  {"chrom": "12", "strand": "-", "window": (25205246, 25250936),"tx": "ENST00000256078", "roles": ["coding"]},
    "APC":   {"chrom": "5",  "strand": "+", "window": (112707498, 112846239),"tx": "ENST00000257430", "roles": ["coding", "promoter", "burden"]},
    "MLH1":  {"chrom": "3",  "strand": "+", "window": (36993332, 37050918),"tx": "ENST00000231790", "roles": ["promoter"]},
    "RB1":   {"chrom": "13", "strand": "+", "window": (48303748, 48481890),"tx": "ENST00000267163", "roles": ["cnv"]},
    "BRCA2": {"chrom": "13", "strand": "+", "window": (32315508, 32400268),"tx": "ENST00000380152", "roles": ["cnv"]},
    "MYC":   {"chrom": "8",  "strand": "+", "window": (127735434, 127742951),"tx": "ENST00000621592", "roles": ["cnv"]},
    "HTT":   {"chrom": "4",  "strand": "+", "window": (3074681, 3243960),  "tx": "ENST00000355072", "roles": ["str"]},
}

# --------------------------------------------------------------------- STR (real HTT)
# Real GRCh38 HTT exon-1 CAG tract (the standard ExpansionHunter HTT locus). The resolver
# VERIFIES the motif at these coordinates against the real FASTA before use.
#   normal_units : within the normal HTT range (median ~17-19 CAG).
#   expanded_units: a clearly pathogenic but SHORT-READ-DETECTABLE expansion. 45 CAG = 135 bp
#     repeat; with 150 bp reads ExpansionHunter detects the PRESENCE of the expansion via
#     flanking/in-repeat reads, but exact sizing of repeats near/over read length is noisy
#     -> the checker requires "detectably expanded vs normal", not an exact length (see
#     STR_LEN_TOL_UNITS / STR_EXPANSION_MIN_DELTA_UNITS).
STR_LOCUS = {
    "gene": "HTT", "chrom": "4", "motif": "CAG",
    "ref_region": (3074877, 3074933),     # 1-based inclusive CAG tract (GRCh38); resolver validates
    "normal_units": 18, "expanded_units": 45,
    "carriers": ["S02", "S04", "S06", "S08"],
}

# ------------------------------------------------------ sex regions (real X/Y, off-target)
# A real targeted panel infers sex from OFF-TARGET X/Y reads (depth ratio + Y presence).
# This sim mirrors that: small REAL X and Y windows (male-specific Y, outside PAR) are
# simulated at ploidy-scaled depth (F: 2x X, 0 Y; M: 1x X, 1x Y) so somalier infers sex
# from real reference coverage — NOT a synthetic contig. Resolver/simulator verify non-N.
SEX_REGIONS = {
    "X": {"chrom": "X", "window": (100000000, 100050000)},   # unique X (away from PAR/segdups)
    "Y": {"chrom": "Y", "window": (6900000, 6950000)},       # male-specific Y (outside PAR1)
}

# --------------------------------------------------------------------------- samples
# 16 samples. Cases = S01..S08, controls = S09..S16. Sex assigned independently of case
# status (8 M / 8 F). Covariates numeric.
def _build_samples():
    out = []
    for i in range(1, 17):
        sid = f"S{i:02d}"
        out.append({"sample_id": sid,
                    "phenotype": "case" if i <= 8 else "control",
                    "sex": "M" if i % 2 == 1 else "F",
                    "covariate_age": 30 + (i * 37) % 35,
                    "covariate_pc1": round(((i * 53) % 100) / 100.0 - 0.5, 3),
                    "batch": "b1" if i % 2 == 0 else "b2"})
    return out

SAMPLES = _build_samples()
SAMPLE_IDS = [s["sample_id"] for s in SAMPLES]
CASES = [s["sample_id"] for s in SAMPLES if s["phenotype"] == "case"]
CONTROLS = [s["sample_id"] for s in SAMPLES if s["phenotype"] == "control"]

# --------------------------------------------------------------- coding showcase set
# Heterozygous coding variants for CALLING recall (genotype-aware) + ANNOTATION
# consequence accuracy. `kind` -> the intended REAL VEP consequence; the resolver places
# each at a real exonic position of `tx` that ACTUALLY produces that consequence (verified
# by translating the real CDS) and is novel/mappable.
CODING_VARIANTS = [
    {"id": "cv_TP53_stop",  "gene": "TP53",  "kind": "stop",       "carriers": ["S01", "S05", "S12"]},
    {"id": "cv_EGFR_mis",   "gene": "EGFR",  "kind": "missense",   "carriers": ["S02", "S07", "S10"]},
    {"id": "cv_BRCA1_fs",   "gene": "BRCA1", "kind": "frameshift", "carriers": ["S03", "S09", "S14"]},
    {"id": "cv_PTEN_stop",  "gene": "PTEN",  "kind": "stop",       "carriers": ["S04", "S11"]},
    {"id": "cv_KRAS_mis",   "gene": "KRAS",  "kind": "missense",   "carriers": ["S06", "S13", "S16"]},
]
EXPECTED_CSQ = {"stop": "stop_gained", "missense": "missense_variant",
                "frameshift": "frameshift_variant", "promoter": "upstream_gene_variant"}

# ------------------------------------------------------------------ promoter variants
# SNVs a short distance upstream of the transcript TSS (real coords) -> VEP
# upstream_gene_variant (regulatory class).
PROMOTER_VARIANTS = [
    {"id": "pv_APC_prom",  "gene": "APC",  "carriers": ["S01", "S08"]},
    {"id": "pv_MLH1_prom", "gene": "MLH1", "carriers": ["S05", "S15"]},
]

# ------------------------------------------------------------------------ CNV events
# Per-sample read DEPTH is modulated over the gene's real exons so CNVkit (depth-based,
# panel-of-normals from the cohort) calls del/dup. fold<1 -> DEL, fold>1 -> DUP.
CNV_EVENTS = [
    {"id": "cnv_RB1_del",   "gene": "RB1",   "state": "DEL", "fold": 0.45, "carriers": ["S03", "S07"]},
    {"id": "cnv_BRCA2_del", "gene": "BRCA2", "state": "DEL", "fold": 0.45, "carriers": ["S10", "S13"]},
    {"id": "cnv_MYC_dup",   "gene": "MYC",   "state": "DUP", "fold": 1.8,  "carriers": ["S05", "S12"]},
]
CNV_TARGET_GENES = sorted({e["gene"] for e in CNV_EVENTS})

# --------------------------------------------------------------------- burden signal
# APC carries the planted burden signal: DISTINCT novel rare coding variants enriched in
# cases. Each is a real het coding SNV (so calling/annotation see it); against REAL
# gnomAD-AFR they read as rare/qualifying. Collectively they make APC the most
# case-enriched gene so the burden test surfaces APC. Mostly cases + one control.
BURDEN_GENE = "APC"
BURDEN_VARIANTS = [
    {"id": "bv_APC_1", "gene": "APC", "kind": "missense", "carriers": ["S01"]},
    {"id": "bv_APC_2", "gene": "APC", "kind": "missense", "carriers": ["S02"]},
    {"id": "bv_APC_3", "gene": "APC", "kind": "stop",     "carriers": ["S03"]},
    {"id": "bv_APC_4", "gene": "APC", "kind": "missense", "carriers": ["S04", "S05"]},
    {"id": "bv_APC_5", "gene": "APC", "kind": "missense", "carriers": ["S06"]},
    {"id": "bv_APC_6", "gene": "APC", "kind": "stop",     "carriers": ["S07", "S09"]},  # S09 = lone control
]

# Every planted variant that lands in coding sequence (showcase + burden).
ALL_CODING = CODING_VARIANTS + BURDEN_VARIANTS

# --------------------------------------------- per-stage real-reference dependency map
# What each stage's CHECK requires to be staged on the cluster. The SPINE
# (align->QC->call) needs ONLY the genome (+.dict/+bwa-mem2 index built on the cluster) so
# it can be validated BEFORE the heavy annotation DBs are downloaded. The checker SKIPs a
# stage (with a "needs <ref>" message) when its output is absent.
#   stage number is the CallForge stage the check targets (for --through).
STAGE_REFS = {
    "calling":  {"num": 7,  "needs": ["genome"]},
    "cnv":      {"num": 8,  "needs": ["genome"]},
    "str":      {"num": 9,  "needs": ["genome", "str_catalog (provided by harness)"]},
    "annotation": {"num": 11, "needs": ["vep_cache", "phylop_bigwig"]},
    "cohortqc": {"num": 12, "needs": ["somalier_sites"]},
    "burden":   {"num": 14, "needs": ["gnomad_afr (rarity)", "vep (qualifying consequence)"]},
}

# --------------------------------------------------------------------- checker knobs
# Named thresholds with rationale. Tolerances are loose where a stage is probabilistic
# (calling/CNV/STR) — we assert the planted signal SURFACES against the REAL reference,
# not exact metrics. See README "Checker thresholds".
THRESHOLDS = {
    # CALLING (6/7): fraction of planted het coding variants recovered, genotype-aware
    # (right sample, right ALT) by GATK against real GRCh38.
    "CALLING_MIN_RECALL": 0.80,
    # ANNOTATION (11): fraction of recovered coding variants whose REAL VEP consequence
    # matches the intended term; promoter -> upstream/regulatory. Deterministic given a
    # call -> set high.
    "ANNOTATION_MIN_CSQ_MATCH": 0.80,
    # STR (9): a carrier's larger allele must exceed normal by >= this many units to count
    # as a detected expansion (presence, not exact size — short reads near/over the repeat
    # length size noisily)...
    "STR_EXPANSION_MIN_DELTA_UNITS": 8,
    # ...and, when ExpansionHunter does size it, within this (wide) tolerance of planted.
    "STR_LEN_TOL_UNITS": 20,
    "STR_MIN_CARRIER_RECALL": 0.75,
    # CNV (8): fraction of planted del/dup events called in the right sample+gene with the
    # correct type. CNVkit on a tiny cohort PoN is noisy -> loose.
    "CNV_MIN_RECALL": 0.50,
    # BURDEN (14): the planted burden gene must rank #1 (smallest p) OR p below this.
    "BURDEN_P_MAX": 0.05,
    # COHORT QC (12): somalier sex (from real X/Y off-target coverage) vs the sheet; max
    # pairwise relatedness (independent samples).
    "COHORTQC_SEX_MIN_CONCORDANCE": 0.90,
    "COHORTQC_MAX_RELATEDNESS": 0.25,
}


def carrier_map():
    """feature_id -> set(carrier sample_ids), across every planted class."""
    m = {}
    for v in ALL_CODING + PROMOTER_VARIANTS + CNV_EVENTS:
        m[v["id"]] = set(v["carriers"])
    m[f"str_{STR_LOCUS['gene']}"] = set(STR_LOCUS["carriers"])
    return m


if __name__ == "__main__":
    print(f"build={GENOME_BUILD} samples={len(SAMPLES)} cases={len(CASES)} controls={len(CONTROLS)}")
    print(f"genes={len(GENES)} coding_variants={len(ALL_CODING)} "
          f"promoter={len(PROMOTER_VARIANTS)} cnv={len(CNV_EVENTS)} "
          f"str={STR_LOCUS['gene']}({STR_LOCUS['normal_units']}->{STR_LOCUS['expanded_units']}) "
          f"burden_gene={BURDEN_GENE}")

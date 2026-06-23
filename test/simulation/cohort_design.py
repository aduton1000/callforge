#!/usr/bin/env python3
"""cohort_design.py — the SINGLE SOURCE OF TRUTH for the CallForge simulation cohort.

This module declares (with NO I/O) the synthetic 16-sample cohort: the generic gene
set, the per-class planted variants, the sample x feature genotype matrix, the
case/control + covariate design, and the named checker thresholds.

Three scripts import it so the plan, the generated reads, and the checks can never
drift apart:
  * simulate_cohort.py     -> writes the human-reviewable plan + sample sheet (Mac, no deps)
  * build_fixture.py       -> resolves coords, plants variants, emits reads + truth manifest
  * check_against_truth.py -> asserts each stage's output against the manifest

Design intent (see README.md):
  * NEUTRAL gene set — well-known cancer/Mendelian illustrative genes (TP53, EGFR, …)
    plus HTT for the STR locus. NONE is a malaria-panel gene.
  * Ensembl no-alt contig naming (1,2,…,X — no `chr` prefix), matching the real pipeline.
  * Every class the pipeline calls is exercised: coding SNV/indel (stop/missense/
    frameshift), promoter/regulatory, STR expansion, CNV del/dup, and a case/control
    burden signal.
  * Deterministic / seeded so re-runs are byte-stable.

stdlib-only (no third-party imports) so it loads without the conda env active.
"""

# Master seed — every random choice downstream derives from this (see build_fixture.py).
SEED = 20260623

# --------------------------------------------------------------------------- genes
# Generic, well-characterized, NON-panel genes. `chrom` is the real Ensembl
# chromosome label (no `chr` prefix); the synthetic mini-genome lays each gene out on
# a contig of that name at SYNTHETIC local coordinates (NOT real genomic positions —
# see README honest-scope note). Multiple genes may share a chromosome contig.
#   role: which class(es) the gene carries in the simulation.
GENES = {
    "TP53":  {"chrom": "17", "roles": ["coding"]},                 # showcase stop_gained
    "BRCA1": {"chrom": "17", "roles": ["coding"]},                 # showcase frameshift
    "NF1":   {"chrom": "17", "roles": ["coding_negative"]},        # negative control (no planted variant)
    "EGFR":  {"chrom": "7",  "roles": ["coding"]},                 # showcase missense
    "PTEN":  {"chrom": "10", "roles": ["coding"]},                 # showcase stop_gained
    "KRAS":  {"chrom": "12", "roles": ["coding"]},                 # showcase missense
    "APC":   {"chrom": "5",  "roles": ["coding", "promoter", "burden"]},  # burden signal + promoter
    "MLH1":  {"chrom": "3",  "roles": ["promoter"]},               # promoter/regulatory
    "RB1":   {"chrom": "13", "roles": ["cnv"]},                    # CNV deletion target
    "BRCA2": {"chrom": "13", "roles": ["cnv"]},                    # CNV deletion target
    "MYC":   {"chrom": "8",  "roles": ["cnv"]},                    # CNV duplication target
    "HTT":   {"chrom": "4",  "roles": ["str"]},                    # STR (CAG)n expansion locus
}

# --------------------------------------------------------------------------- samples
# 16 samples. Cases = S01..S08, controls = S09..S16. Sex is assigned independently of
# case status (8 M / 8 F) and is encoded deterministically into X off-target reads so
# somalier (stage 12) can be checked for sex concordance. Covariates are numeric.
def _build_samples():
    out = []
    for i in range(1, 17):
        sid = f"S{i:02d}"
        pheno = "case" if i <= 8 else "control"
        sex = "M" if i % 2 == 1 else "F"
        age = 30 + (i * 37) % 35                 # deterministic spread 30..64
        pc1 = round(((i * 53) % 100) / 100.0 - 0.5, 3)   # deterministic ~[-0.5,0.5)
        batch = "b1" if i % 2 == 0 else "b2"
        out.append({"sample_id": sid, "phenotype": pheno, "sex": sex,
                    "covariate_age": age, "covariate_pc1": pc1, "batch": batch})
    return out

SAMPLES = _build_samples()
SAMPLE_IDS = [s["sample_id"] for s in SAMPLES]
CASES = [s["sample_id"] for s in SAMPLES if s["phenotype"] == "case"]
CONTROLS = [s["sample_id"] for s in SAMPLES if s["phenotype"] == "control"]

# --------------------------------------------------------------- coding showcase set
# Clearly-callable, heterozygous coding variants used to measure CALLING recall
# (genotype-aware) and ANNOTATION consequence accuracy. `kind` drives how the
# generator plants it and which exact VEP consequence the checker requires.
#   kind: "stop"=stop_gained, "missense"=missense_variant, "frameshift"=frameshift_variant
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
# SNVs upstream of the transcript start -> VEP upstream_gene_variant (regulatory class).
PROMOTER_VARIANTS = [
    {"id": "pv_APC_prom",  "gene": "APC",  "carriers": ["S01", "S08"]},
    {"id": "pv_MLH1_prom", "gene": "MLH1", "carriers": ["S05", "S15"]},
]

# ------------------------------------------------------------------------- STR locus
# HTT (CAG)n — the textbook tandem repeat (Huntingtin); neutral, non-panel, added
# explicitly as the STR test locus. Carriers are HETEROZYGOUS expanded (one normal,
# one expanded allele); non-carriers are normal/normal.
STR_LOCUS = {
    "gene": "HTT", "motif": "CAG", "normal_units": 18, "expanded_units": 45,
    "carriers": ["S02", "S04", "S06", "S08"],
}

# ------------------------------------------------------------------------ CNV events
# Read DEPTH is modulated per sample over the target region so CNVkit (depth-based,
# panel-of-normals from the cohort) calls del/dup. fold<1 -> DEL, fold>1 -> DUP.
CNV_EVENTS = [
    {"id": "cnv_RB1_del",   "gene": "RB1",   "state": "DEL", "fold": 0.45, "carriers": ["S03", "S07"]},
    {"id": "cnv_BRCA2_del", "gene": "BRCA2", "state": "DEL", "fold": 0.45, "carriers": ["S10", "S13"]},
    {"id": "cnv_MYC_dup",   "gene": "MYC",   "state": "DUP", "fold": 1.8,  "carriers": ["S05", "S12"]},
]
CNV_TARGET_GENES = sorted({e["gene"] for e in CNV_EVENTS})

# --------------------------------------------------------------------- burden signal
# APC carries the planted burden signal: several DISTINCT rare coding variants
# enriched in cases. Each is a real heterozygous coding SNV (so calling/annotation
# also see it); collectively they make APC the most case-enriched gene so the burden
# test (stage 14) surfaces APC. Mostly cases + a single control for realism.
BURDEN_GENE = "APC"
BURDEN_VARIANTS = [
    {"id": "bv_APC_1", "gene": "APC", "kind": "missense", "carriers": ["S01"]},
    {"id": "bv_APC_2", "gene": "APC", "kind": "missense", "carriers": ["S02"]},
    {"id": "bv_APC_3", "gene": "APC", "kind": "stop",     "carriers": ["S03"]},
    {"id": "bv_APC_4", "gene": "APC", "kind": "missense", "carriers": ["S04", "S05"]},
    {"id": "bv_APC_5", "gene": "APC", "kind": "missense", "carriers": ["S06"]},
    {"id": "bv_APC_6", "gene": "APC", "kind": "stop",     "carriers": ["S07", "S09"]},  # S09 = lone control
]

# All planted variants that land in coding sequence (showcase + burden) — the
# generator plants every one as a real read-level heterozygous variant.
ALL_CODING = CODING_VARIANTS + BURDEN_VARIANTS

# --------------------------------------------------------------------- checker knobs
# Named thresholds with rationale. Calling/CNV/STR are probabilistic on short
# synthetic reads, so tolerances are deliberately loose — we assert the planted
# signal SURFACES, not exact metrics. See README "Checker thresholds".
THRESHOLDS = {
    # CALLING (stage 6/7): fraction of planted het coding variants recovered in the
    # right sample with the right ALT. High-depth clearly-callable sites should mostly
    # be recovered; allow misses for indels/edge cases.
    "CALLING_MIN_RECALL": 0.80,
    # ANNOTATION (stage 11): of the recovered showcase coding variants, fraction
    # carrying the exact expected VEP consequence; promoter variants must be
    # upstream/regulatory. Annotation is deterministic given a call, so set high.
    "ANNOTATION_MIN_CSQ_MATCH": 0.80,
    # STR (stage 9): a carrier's larger allele must exceed its normal allele by at
    # least this many repeat units to count as a detected expansion...
    "STR_EXPANSION_MIN_DELTA_UNITS": 10,
    # ...and the called expanded allele should be within this tolerance of the planted
    # length (EH genotyping is noisy on short reads).
    "STR_LEN_TOL_UNITS": 12,
    "STR_MIN_CARRIER_RECALL": 0.75,   # fraction of STR carriers flagged expanded
    # CNV (stage 8): fraction of planted del/dup events called in the right sample &
    # gene with the correct type. CNVkit on a tiny cohort PoN is noisy -> loose.
    "CNV_MIN_RECALL": 0.50,
    # BURDEN (stage 14): the planted burden gene must rank #1 (smallest p) OR have
    # p below this ceiling. We assert the signal surfaces, not an exact p-value.
    "BURDEN_P_MAX": 0.05,
    # COHORT QC (stage 12): somalier sex inference vs the sheet, and max pairwise
    # relatedness (independent samples).
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
    # Quick self-summary (no side effects beyond printing).
    print(f"samples={len(SAMPLES)} cases={len(CASES)} controls={len(CONTROLS)}")
    print(f"genes={len(GENES)} coding_variants={len(ALL_CODING)} "
          f"promoter={len(PROMOTER_VARIANTS)} cnv={len(CNV_EVENTS)} "
          f"str_carriers={len(STR_LOCUS['carriers'])} burden_gene={BURDEN_GENE}")

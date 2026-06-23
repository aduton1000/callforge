# CallForge cohort simulation harness

A self-contained, deterministic **16-sample synthetic cohort** with **planted ground
truth across every class the pipeline calls**, used to validate the full 16-stage
CallForge pipeline END-TO-END on the cluster for **correctness** — did each stage detect
what we planted? — not merely "did it run."

It emits a synthetic reference + per-sample FASTQs + sample sheet + phenotype/covariates
+ a machine-readable **truth manifest**, and ships a **checker** that asserts each stage's
output against that manifest.

> ### Honest scope note — read this first
> Simulated reads validate the pipeline **plumbing** and **planted-signal recovery**: that
> reads flow through all 16 stages and that the variants we deliberately inserted are
> recovered with the right class/consequence/genotype. They are **NOT** a substitute for
> real-data calibration or formal GIAB precision/recall, which require a real GIAB sample
> and real genomic backgrounds. "16 stages green on the sim" means **the wiring is correct
> and planted signals surface** — it does **not** mean the panel is clinically validated.
> The mini-genome is synthetic: contigs carry Ensembl chromosome **labels** (`1,2,…,X`) but
> the sequence and coordinates are fabricated, not real genomic positions.

This fixture is a **separate, cluster-run** test asset. It is **not** wired into any
Mac/offline CI (it needs the `.sif` + a Nextflow run). It is distinct from the small
`conf/test.config` fixture (`bin/make_test_reference.py`), which is a 4-sample smoke test.

---

## Files

| File | Runs where | Purpose |
|------|-----------|---------|
| `cohort_design.py` | imported | **Single source of truth**: gene set, planted variants, sample×feature matrix, case/control + covariates, checker thresholds. |
| `simulate_cohort.py` | host (stdlib only) | Writes the reviewable **plan** + `samplesheet.csv` (no genome/simulator needed). |
| `build_fixture.py` | in-container (needs htslib) | Builds the mini-genome, plants variants, emits reads (het balance, CNV depth modulation, STR expansion), writes the **resolved `truth_manifest.json`**. |
| `generate_reads.sh` | cluster | Wrapper: runs the plan on the host + `build_fixture.py` inside the `.sif`. |
| `check_against_truth.py` | cluster (after the run) | Asserts each stage's results against the manifest; PASS/FAIL per stage + overall. |

---

## How to run on the cluster

```bash
# 0) clone/pull the repo on the cluster, cd into it. Choose a WRITABLE fixture dir
#    (scratch/project space — NOT the read-only install).
cd test/simulation

# 1) generate the fixture (reference + reads + truth manifest), in-container:
CALLFORGE_SIF=/path/to/callforge.sif ./generate_reads.sh /scratch/$USER/cf_sim singularity

# 2) fill the two FILL_ME params (container_image, slurm_partition) in the generated
#    /scratch/$USER/cf_sim/sim.params.yaml, then run the FULL pipeline from a writable dir:
cd /scratch/$USER/cf_sim
callforge-run -params-file sim.params.yaml -profile hpc_slurm,singularity

# 3) check the results against the planted truth:
python3 <repo>/test/simulation/check_against_truth.py \
    --manifest /scratch/$USER/cf_sim/truth_manifest.json \
    --results  /scratch/$USER/cf_sim/results
```

The generator is hermetic (pure-Python read emission; only htslib is required), so it can
also be built on a workstation with `samtools/bgzip/tabix` for inspection — only the
Nextflow run itself needs the cluster + `.sif`. Write the fixture **outside** the repo
(scratch) so generated reads/results are never committed.

---

## The gene set (generic, non-panel) and why each was chosen

Well-known cancer/Mendelian **illustrative** genes (plus `HTT` for the STR locus). **None
is a malaria-panel gene.** Each gene is assigned its real chromosome (for plausible
Ensembl naming) but laid out at synthetic local coordinates.

| Gene | Chrom | Role in the sim |
|------|-------|-----------------|
| `TP53`  | 17 | coding **stop-gain** (calling + annotation showcase) |
| `EGFR`  | 7  | coding **missense** |
| `BRCA1` | 17 | coding **frameshift** (1-bp deletion) |
| `PTEN`  | 10 | coding **stop-gain** |
| `KRAS`  | 12 | coding **missense** |
| `APC`   | 5  | **burden gene** (rare coding variants enriched in cases) **+ promoter** SNV |
| `MLH1`  | 3  | **promoter / regulatory** SNV (upstream of TSS) |
| `RB1`   | 13 | **CNV deletion** target |
| `BRCA2` | 13 | **CNV deletion** target |
| `MYC`   | 8  | **CNV duplication** target |
| `HTT`   | 4  | **STR** `(CAG)n` expansion locus (textbook tandem repeat; added explicitly as the STR test) |
| `NF1`   | 17 | **negative control** gene (no planted variant) |

The set deliberately spans all classes: coding SNV/indel (stop/missense/frameshift),
promoter/regulatory, STR expansion, CNV del/dup, and a case/control burden signal. `HTT`
is added as the STR locus because none of the cancer genes carries a canonical disease STR.

### Planted-variant matrix (summary)

16 samples `S01..S16`; **cases = S01–S08**, **controls = S09–S16**. All planted small
variants are **heterozygous**. Carriers (see `cohort_plan.tsv` / `truth_manifest.tsv` for
the full matrix):

- **Coding showcase:** TP53 stop (S01,S05,S12), EGFR missense (S02,S07,S10), BRCA1
  frameshift (S03,S09,S14), PTEN stop (S04,S11), KRAS missense (S06,S13,S16).
- **Promoter:** APC (S01,S08), MLH1 (S05,S15) → `upstream_gene_variant`.
- **STR:** HTT `(CAG)18 → (CAG)45` expanded in S02,S04,S06,S08; normal elsewhere.
- **CNV:** RB1 DEL (S03,S07), BRCA2 DEL (S10,S13), MYC DUP (S05,S12) — via per-sample
  read-**depth** modulation over the region (CNVkit is depth/PoN based).
- **Burden:** 6 distinct rare APC coding variants concentrated in cases (S01–S07) with a
  single control (S09) for realism → APC is the most case-enriched gene.

---

## Truth manifest schema (`truth_manifest.json`)

Resolved by `build_fixture.py` (exact coords/ref/alt depend on the generated sequence):

```
seed, genome_fasta, target_bed, gtf, gene_metadata, str_catalog, somalier_sites,
gnomad_vcf, contigs:{name:length}, samples:[{sample_id,phenotype,sex,covariate_*,batch}],
cases:[…], controls:[…],
coding_variants:[{id,gene,contig,pos(1-based),ref,alt,vartype,consequence,zygosity,carriers}],
promoter_variants:[{… consequence:"upstream_gene_variant" …}],
str:{gene,contig,repeat_region,locus_id,motif,normal_units,expanded_units,carriers},
cnv_events:[{id,gene,contig,start,end,state(DEL|DUP),fold,carriers}],
burden:{gene,qualifying_feature_ids,case_carriers,control_carriers},
cohortqc:{sex_by_sample:{sample:M|F}},
giab:{control_id,truth_vcf,truth_bed},
thresholds:{…}
```

A human-readable `truth_manifest.tsv` (class / feature / gene / contig / coord / detail /
carriers) is written alongside it.

---

## Checker thresholds (named, with rationale)

Defined as `THRESHOLDS` in `cohort_design.py` and carried in the manifest. Tolerances are
loose where a stage is probabilistic — we assert the planted signal **surfaces**, not exact
metrics.

| Threshold | Value | Rationale |
|-----------|-------|-----------|
| `CALLING_MIN_RECALL` | 0.80 | Fraction of planted het coding variants recovered **genotype-aware** (right sample, right ALT). Clearly-callable high-depth sites should mostly be recovered; allow indel/edge misses. |
| `ANNOTATION_MIN_CSQ_MATCH` | 0.80 | Of recovered coding variants, fraction carrying the **exact** expected VEP consequence (`stop_gained`/`missense_variant`/`frameshift_variant`); promoter → `upstream`/`regulatory`. Annotation is deterministic given a call → set high. |
| `STR_EXPANSION_MIN_DELTA_UNITS` | 10 | A carrier's larger allele must exceed normal by ≥ this many repeat units to count as a detected expansion. |
| `STR_LEN_TOL_UNITS` | 12 | Called expanded allele within ± this of the planted length (EH genotyping is noisy on short reads). |
| `STR_MIN_CARRIER_RECALL` | 0.75 | Fraction of STR carriers flagged expanded. |
| `CNV_MIN_RECALL` | 0.50 | Fraction of planted del/dup events called in the right sample+gene with correct type. CNVkit on a tiny cohort PoN is noisy → loose. |
| `BURDEN_P_MAX` | 0.05 | Planted burden gene (APC) must rank #1 (smallest p) **or** have p below this. We assert the signal surfaces, not an exact p. |
| `COHORTQC_SEX_MIN_CONCORDANCE` | 0.90 | somalier-inferred sex vs the sheet (sex is encoded deterministically into X off-target reads). |
| `COHORTQC_MAX_RELATEDNESS` | 0.25 | Samples are independent → no pair above this. |

What each stage's check asserts:

- **calling (6/7)** — `joint.filtered.vcf.gz`: ≥ `CALLING_MIN_RECALL` of planted coding
  variants recovered, genotype-aware (majority of carriers non-ref at the site; indels
  matched as an indel call within ±30 bp of the locus).
- **annotation (11)** — `variants.flat.tsv`: planted stop-gain/missense/frameshift carry the
  expected VEP `consequence`; promoter variants flagged `upstream`/`regulatory`.
- **STR (9)** — `str_calls.tsv`: carriers show an expanded allele (~45 units ± tol); spurious
  expansion in non-carriers is reported as a warning.
- **CNV (8)** — `cnv_calls.tsv`: planted DEL/DUP events called in the right sample+gene with
  the correct `type`.
- **burden (14)** — `burden_results.tsv`: the APC unit ranks #1 or p ≤ `BURDEN_P_MAX`.
- **cohort QC (12)** — `cohort_qc.json`: somalier `sex_inferred` matches the sheet; no
  pairwise relatedness above the ceiling.

---

## Determinism

Everything derives from `cohort_design.SEED` via a small fixed LCG (no `hash()` salt, no
`Math.random`), so re-runs are byte-stable across machines. Change the cohort by editing
`cohort_design.py` only — the plan, the generated reads, and the checks stay in lock-step.

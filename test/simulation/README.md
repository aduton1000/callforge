# CallForge cohort simulation harness — REAL-reference validation

A deterministic **16-sample cohort** with **planted ground truth across every class the
pipeline calls**, used to validate the full 16-stage CallForge pipeline END-TO-END on the
cluster **against the REAL GRCh38 no-alt reference and the REAL annotation databases** —
nothing faked. "Green" here means **CallForge works against real references on the
cluster**, which is what de-risks the real cohort.

Reads are simulated **from the real reference at the genes' real coordinates** with a real
Illumina-error-model simulator (**dwgsim**); planted variants sit at real exonic positions
whose consequence is **verified against the real gene model**; and the pipeline runs against
**real** VEP/gnomAD-AFR/dbSNP/ClinVar/PhyloP/somalier references.

> ### Honest scope note — read this first
> Simulated reads validate that CallForge runs correctly **against real references** and
> recovers **planted** signals at real coordinates. They are **NOT** a substitute for
> real-data calibration or formal GIAB precision/recall (which need a real GIAB sample and
> real genomic backgrounds). GIAB benchmarking applies to the **real** cohort, not this sim
> (the sim has no real GIAB sample) — so the sim checker does **not** include a GIAB stage,
> while `bin/fetch_annotation_refs.sh` still stages HG002 truth for the real run.

This is a **separate, cluster-run** asset (read-sim + pipeline run happen on the cluster).
It is **not** wired into any Mac/offline CI. It is distinct from the small `conf/test.config`
fixture (`bin/make_test_reference.py`), which is a 4-sample hermetic smoke test.

---

## What changed from the mini-genome version (commit 9565e12)

**Rebased onto real GRCh38** (kept the good parts): the cohort design, the planted-variant
matrix, the truth-manifest schema, the checkers + named thresholds, and the generic gene
set (TP53, EGFR, BRCA1, PTEN, KRAS, APC, MLH1, RB1, BRCA2, MYC, HTT, NF1 — no panel genes).
**Removed the fakes**: no synthetic mini-genome, no hand-rolled Python read emission, no
faked gnomAD/VEP/GIAB. Reads now come from the real FASTA via **dwgsim**; consequences are
verified against the real CDS; the pipeline runs against the real staged DBs.

---

## Files

| File | Runs where | Purpose |
|------|-----------|---------|
| `cohort_design.py` | imported | **Single source of truth**: gene set + REAL GRCh38 coords, planted matrix, samples/covariates, thresholds, per-stage ref map. |
| `simulate_cohort.py` | host (stdlib) | Reviewable plan + `samplesheet.csv` + manifest **skeleton** (no genome needed). |
| `resolve_variants.py` | in-`.sif` | Places each variant at a **real** exonic/upstream position, **verifies the consequence against the real CDS**, checks mappability + novelty; writes the **resolved** `truth_manifest.json` + real-coord `targets.bed` + real HTT `str_catalog.json` + `gene_metadata.tsv`. Needs FASTA + GTF (+ optional known-VCF). |
| `simulate_reads.py` | in-`.sif` | **dwgsim** reads from the real reference: het haplotypes, CNV depth modulation, expanded HTT repeat, real X/Y off-target for sex. |
| `generate_reads.sh` | cluster | Wrapper: host plan, then resolve + simulate inside the `.sif`; writes `sim.params.yaml`. |
| `check_against_truth.py` | cluster | **Staged** PASS/FAIL vs the manifest (`--stages` / `--through`; SKIP-if-ref-absent). |

`bin/fetch_annotation_refs.sh` (repo `bin/`) stages the real DBs (resumable + verified).

---

## Real-reference dependency map (what each stage's CHECK needs)

| Stage check | `--through` | Needs staged |
|-------------|:-----------:|--------------|
| calling (6/7) | 7 | **genome only** (+ .dict/.bwa-mem2 index built on the cluster) |
| CNV (8) | 8 | genome only |
| STR (9) | 9 | genome + the STR catalog this harness provides |
| annotation (11) | 11 | **VEP cache** + **PhyloP** |
| cohort QC (12) | 12 | **somalier sites** |
| burden (14) | 14 | **gnomAD-AFR** (rarity) + VEP (qualifying consequence) |

The **spine** (`--through 9`: calling/cnv/str) needs **only the genome**, so it is validated
**before** the heavy annotation DBs finish downloading. Read GENERATION additionally needs a
gene-model **GTF** (small, ~50 MB) to place variants at real coords with verified
consequences — this is separate from the heavy annotation DBs.

---

## How to run on the cluster

```bash
# 0) on the cluster, build the genome .dict + bwa-mem2 index once (standard CallForge prep).
#    Pick a WRITABLE fixture dir (scratch/project — NOT the read-only install).

# 1) generate the fixture: resolve real-coord variants + simulate reads (both in the .sif)
cd test/simulation
CALLFORGE_SIF=/path/callforge.sif ./generate_reads.sh \
    /scratch/$USER/cf_sim \
    /refs/GRCh38_noalt.clean.fa \
    /refs/Homo_sapiens.GRCh38.110.gtf.gz \
    singularity \
    /refs/gnomad.afr.sites.vcf.gz        # optional 5th arg -> novelty/collision checks

# 2) run + check the SPINE first (genome only — before annotation DBs are staged):
cd /scratch/$USER/cf_sim
callforge-run -params-file sim.params.yaml -profile hpc_slurm,singularity   # spine stages
python3 <repo>/test/simulation/check_against_truth.py \
    --manifest truth_manifest.json --results results --through 9

# 3) stage the annotation/burden/cohortqc DBs (resumable + verified), then fill the
#    commented params in sim.params.yaml and re-run those stages:
OUTDIR=/refs/callforge <repo>/bin/fetch_annotation_refs.sh
callforge-run -params-file sim.params.yaml -profile hpc_slurm,singularity   # later stages
python3 <repo>/test/simulation/check_against_truth.py \
    --manifest truth_manifest.json --results results            # all stages
```

`dwgsim` must be in the `.sif` — it was added to `env/callforge.yml`, so **rebuild the
image** (`apptainer build callforge.sif env/callforge.def`) before step 1.

---

## How planted variants are placed at real coordinates

`resolve_variants.py` (in-container, with the real FASTA + GTF):

1. Looks up each gene's chosen transcript in the **real gene model** (MANE/canonical from
   `cohort_design.GENES`, else the longest-CDS transcript) and builds the **spliced CDS**
   from the real FASTA, strand-aware.
2. For each planted coding variant, **scans real codons** for a single-base change that
   produces the intended consequence and **verifies it by translating the real CDS**: a
   `stop_gained` must actually create a stop at that real position; a `missense_variant`
   must change the amino acid (not synonymous, not stop); a `frameshift_variant` is a 1-bp
   exonic indel. Forward-strand ref/alt are read back from the FASTA and checked.
3. **Mappability**: rejects soft-masked/N/low-complexity (homopolymer) context.
   **Novelty/collision**: with an optional gnomAD/dbSNP VCF, rejects positions carrying a
   known variant, so planted variants read as **novel → rare → burden-qualifying** and
   don't confound the annotation/burden checks.
4. Promoter variants are placed ~900 bp upstream of the real TSS (strand-aware) →
   `upstream_gene_variant`. The **HTT** CAG tract motif is validated against the real FASTA
   at the standard GRCh38 locus.

Exact resolved positions/ref/alt/consequence are written into `truth_manifest.json` (the
checker's spec) and echoed at the end of the resolve step for review.

### STR detectability (short reads)

HTT is planted at `(CAG)18 → (CAG)45`. A 45-unit tract (135 bp) is near the 150 bp read
length, so ExpansionHunter reliably detects the **presence** of the expansion via
flanking/in-repeat reads but sizes it noisily. The checker therefore requires the carrier's
larger allele to be **detectably expanded vs normal** (`STR_EXPANSION_MIN_DELTA_UNITS`), with
a **wide** size tolerance (`STR_LEN_TOL_UNITS`) — not an exact length.

### Sex / cohort QC

The gene panel is autosomal, so — like a real targeted panel — sex is inferred from
**real X/Y off-target reads**: small real X and (male-specific) Y windows are simulated at
ploidy-scaled depth (F: 2× X, 0 Y; M: 1× X, 1× Y) so somalier infers sex from **real
reference coverage**, not a synthetic contig.

---

## Truth manifest schema (`truth_manifest.json`, resolved)

```
seed, genome_build, genome_fasta, gtf, target_bed, gene_metadata, str_catalog, sex_regions,
samples:[{sample_id,phenotype,sex,covariate_*,batch}], cases, controls,
coding_variants:[{id,gene,contig,pos,ref,alt,vartype,consequence,zygosity,carriers}],
promoter_variants:[{… consequence:"upstream_gene_variant" …}],
str:{gene,contig,repeat_region,locus_id,motif,motif_validated,normal_units,expanded_units,carriers},
cnv_events:[{id,gene,state(DEL|DUP),fold,carriers}],
burden:{gene,qualifying_feature_ids,case_carriers,control_carriers},
cohortqc:{sex_by_sample}, stage_refs, thresholds
```

A human-readable `truth_manifest.tsv` is written alongside.

---

## Checker thresholds (named, with rationale)

Defined as `THRESHOLDS` in `cohort_design.py` and carried in the manifest.

| Threshold | Value | Rationale |
|-----------|-------|-----------|
| `CALLING_MIN_RECALL` | 0.80 | Planted het coding variants recovered genotype-aware by GATK vs real GRCh38. |
| `ANNOTATION_MIN_CSQ_MATCH` | 0.80 | Recovered variants carrying the exact **real VEP** consequence; promoter → upstream/regulatory. |
| `STR_EXPANSION_MIN_DELTA_UNITS` | 8 | Carrier's larger allele must exceed normal by ≥ this (presence, not exact size). |
| `STR_LEN_TOL_UNITS` | 20 | Wide size tolerance — short reads size near/over-read-length repeats noisily. |
| `STR_MIN_CARRIER_RECALL` | 0.75 | Fraction of STR carriers flagged expanded. |
| `CNV_MIN_RECALL` | 0.50 | Planted del/dup called in the right sample+gene+type; CNVkit on a tiny PoN is noisy. |
| `BURDEN_P_MAX` | 0.05 | Planted burden gene (APC) ranks #1 or p below this. |
| `COHORTQC_SEX_MIN_CONCORDANCE` | 0.90 | somalier sex (from real X/Y coverage) vs the sheet. |
| `COHORTQC_MAX_RELATEDNESS` | 0.25 | Independent samples → no pair above this. |

Run a subset with `--stages calling,cnv,str` or `--through 9`; a stage whose output is
absent SKIPs with a `needs: <ref>` message (use `--strict` to turn SKIPs into failures).

---

## Cluster sequence (summary)

1. build genome `.dict` + bwa-mem2 index (cluster);
2. `generate_reads.sh` → resolve real-coord variants + dwgsim reads from real GRCh38;
3. run + check the **spine** (`--through 9`, genome only);
4. `fetch_annotation_refs.sh` (resumable + checksum-verified) → stage VEP/gnomAD-AFR/PhyloP/…;
5. run + check **annotation / burden / cohort QC**.

Determinism: everything derives from `cohort_design.SEED` (carrier matrix, dwgsim `-z`).
Edit `cohort_design.py` only — plan, reads, and checks stay in lock-step.

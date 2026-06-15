# CallForge

**CallForge — a reproducible Nextflow pipeline for targeted-capture germline analysis:
FASTQ → annotated joint callset (SNV/indel, CNV, STR, paralog-aware) with cohort QC, GIAB
benchmarking, and rare-variant burden testing. Pairs with CaptureForge; runs standalone.**

## What it is

CallForge takes demultiplexed paired-end FASTQs from a hybridization-capture panel through
alignment → per-sample QC (flag/quarantine) → SNV/indel joint calling → CNV / STR /
paralog-aware analysis → annotation → cohort QC → GIAB benchmarking → rare-variant burden
testing — disease- and organism-agnostic, driven by inputs and config.

**It runs standalone.** The minimum inputs are **FASTQs + a reference + a target BED**.
**[CaptureForge](https://github.com/aduton1000/captureforge) is optional enrichment, not a
prerequisite** — its per-gene metadata improves the interpretation layer (CNV callability,
STR loci, paralog flags, burden gene-sets), and `callforge metadata` generates an equivalent
table when you don't have CaptureForge. *Runs standalone; better with CaptureForge.*

## Key features

- **Variant classes:** SNV/indel (GATK HaplotypeCaller→GenomicsDB→GenotypeGVCFs, or
  DeepVariant+GLnexus), CNV (CNVkit), STR (ExpansionHunter), paralog-aware flagging
  (`INFO/PARALOG_GENE`/`PARALOG_CONF`).
- **Per-sample QC gate** that **flags & quarantines** failing samples (depth, dup, on-target,
  contamination, sex) — never silently drops data.
- **Manifest-driven annotation** (VEP + vcfanno: gnomAD-AFR / ClinVar / dbSNP + PhyloP) with a
  **genomic-scope gate** that fails loud on a region-limited DB, and chr-reconciliation.
- **GIAB benchmarking** (hap.py) restricted to the panel BED (on-target only).
- **Rare-variant burden** collapsed by gene & gene-set, with an **ancestry-PC stability gate**.
- **Disease/organism-agnostic**; resources are **discovered, not downloaded**.
- **The `callforge` CLI** front door, **ten per-stage subcommands** (run one stage standalone),
  and per-process **compute controls** (`--max_cpus`/`--max_memory`).

The authoritative guide is **[`docs/manual/callforge_manual.md`](docs/manual/callforge_manual.md)**
(rendered: [`callforge_manual.pdf`](docs/manual/callforge_manual.pdf)) — this README is the front door.

## Install

**Local — any single Linux or macOS machine (Windows via WSL2):**

```bash
git clone https://github.com/aduton1000/callforge.git
cd callforge
mamba env create -f env/callforge.yml   # core env; per-stage envs build on demand
pip install -e .                        # registers the `callforge` command on PATH
callforge --version
```

You also need **Nextflow (≥23.10, Java 17+)** and either conda/mamba or a container engine.
A run composes **two** profiles: an **execution** profile (`local` = this machine, `hpc_slurm`
= a SLURM cluster) and a **packaging** profile (`conda` / `docker` / `apptainer` = how
dependencies are provided) — e.g. `-profile local,conda`. `local` is **not** Mac-specific
(`mac_local` remains a deprecated alias). VEP and hap.py always run from their official
containers.

**HPC:** shared conda envs + Apptainer + an Lmod module put `callforge` on PATH cluster-wide —
see **[`docs/hpc_deployment.md`](docs/hpc_deployment.md)**.

> The raw forms still work (`nextflow run main.nf …`, `python3 bin/…`), but `callforge` is the
> recommended interface. Once the repo is public (or for users with access), you can also run it
> straight from GitHub: `nextflow run aduton1000/callforge -r <tag> -profile local,conda …`.

## Usage — the `callforge` CLI

```text
callforge run          run the full pipeline (wrapper over `nextflow run`)
callforge samplesheet  scaffold/assemble the sample sheet (pairs FASTQs, joins metadata)
callforge metadata     generate gene_metadata.tsv WITHOUT CaptureForge
callforge resources    annotation-resource discovery / genomic-scope checks
callforge --help | --version
```

**Ten per-stage subcommands** run a single stage standalone on its own inputs (e.g. re-annotate
after a DB update, re-burden with a new threshold), writing to a **non-destructive
`results/standalone/<stage>_<timestamp>/`** with their own report + provenance:

| Subcommand | Purpose |
|:-----------|:--------|
| `align`    | FASTQ → analysis-ready BAM(s) (bwa-mem2 + dedup + BQSR) + alignment QC |
| `coverage` | BAM(s) → coverage/enrichment metrics + per-sample QC gate (pass/quarantine) |
| `call`     | QC-pass BAM(s) → joint hard-filtered VCF (GATK, or `--caller deepvariant`) |
| `cnv`      | BAM(s) → CNVkit calls, labelled by CaptureForge callability |
| `str`      | BAM(s) → ExpansionHunter genotypes on the (CaptureForge) STR catalog |
| `paralog`  | VCF → paralog-aware flagging (`PARALOG_GENE`/`PARALOG_CONF`) |
| `anno`     | re-annotate a VCF (VEP + vcfanno DBs + PhyloP), e.g. after a DB update |
| `cohortqc` | somalier relatedness/sex + ancestry PCA + missingness |
| `giab`     | hap.py precision/recall/F1 vs GIAB truth, restricted to the panel BED |
| `burden`   | re-run rare-variant burden with new thresholds/engine on an annotated VCF |

Run `callforge <stage> --help` for that stage's required + optional input files.

## Examples

**Quickstart (test profile — a tiny synthetic fixture, runs in minutes):**

```bash
bash test/make_test_data.sh
callforge run -profile test,docker
```

**Build a sample sheet** (pair a FASTQ directory and join your clinical metadata):

```bash
callforge samplesheet \
    --fastq-dir /data/fastqs \
    --metadata clinical.csv \
    --out sample_sheet.csv
```

**Full cohort run** (`-params-file` supplies reference, target BED, resources, sample sheet):

```bash
callforge run \
    -profile local,conda \
    -params-file params.yaml
```

**Re-annotate after an annotation-database update** (only Stage 11, non-destructive):

```bash
callforge anno \
    --input_vcf results/stage10_paralog/paralog.annotated.vcf.gz \
    --genome_fasta /refs/GRCh38_noalt.fa \
    --target_bed /panel/targets.bed \
    --resource_dirs /refs/annotation_db_2025_06 \
    -profile local,conda
```

**Re-run burden with a new AF cutoff and engine** (reuse the annotated callset):

```bash
callforge burden \
    --input_vcf results/stage11_annotation/annotated.vcf.gz \
    --input sample_sheet.csv \
    --target_bed /panel/targets.bed \
    --burden_engine regenie \
    --burden_af_max 0.005 \
    -profile local,conda
```

**Without CaptureForge** — generate the panel metadata, then run:

```bash
callforge metadata \
    --bed /panel/targets.bed \
    --paralog-genes CD209,CFH \
    --burden-groups burden_groups.tsv \
    --out-tsv gene_metadata.tsv
# then: callforge run -profile local,conda -params-file params.yaml \
#         --gene_metadata gene_metadata.tsv
```

Commands are copy-pasteable; required pipeline params are `--input`, `--genome_fasta`,
`--target_bed` (see `params.example.yaml`).

## Inputs

- **Sample sheet** (CSV): `sample_id,fastq_1,fastq_2,sex?,phenotype?,covariate_*?,batch?` —
  `phenotype` is required only for burden. Build it with `callforge samplesheet`.
- **Target BED** — the panel regions (CaptureForge `final_covered_targets.bed`, or any BED).
- **Reference** — the no-alt GRCh38 primary assembly the panel was designed against.
- **CaptureForge metadata** (optional) — `--captureforge_dir` (auto-finds `metrics.json` +
  `baits.csv`), or generate an equivalent `gene_metadata.tsv` with `callforge metadata`.
- **Annotation databases** (optional) — discovered from `--resource_dirs`
  (`callforge resources`); fetch a panel slice with `bin/fetch_annotation_dbs.sh`.

## Outputs

Annotated joint VCF (`stage11_annotation/annotated.vcf.gz`) + flattened per-variant TSV
(`variants.flat.tsv`); CNV (`cnv_calls.tsv`), STR (`str_calls.tsv`), paralog
(`paralog.annotated.vcf.gz`); a self-contained **cohort QC dashboard**
(`cohort_qc_dashboard.html`); GIAB metrics (`happy.summary.csv`); burden results
(`burden_results.tsv`); and `provenance.json` (tool versions, resources, git commit).

## CaptureForge pairing

CallForge is the **analysis** half of a pair: **[CaptureForge](https://github.com/aduton1000/captureforge)**
does the capture **design** (probe/bait design, target definition, specificity, QC) and CallForge
does the downstream analysis. CaptureForge's handoff (`final_covered_targets.bed` + per-gene
metadata, including design-specific CNV callability) enriches CallForge — but CallForge runs
without it.

## Honesty note

The shipped test outputs and the bundled `docs/cohort_qc_dashboard.html` are generated from a
**synthetic test cohort** — they demonstrate the machinery and plots, **not** real-data quality
or results. Some heavy/alternative paths (DeepVariant+GLnexus, regenie/SKAT-O, the full
align→coverage→call spine) are verified by preview/compile + component tests and carry a
**validate-at-deployment checklist** (`docs/hpc_deployment.md`) for one live end-to-end run
before being relied on for real data.

## Docs · License · Citation

- **Manual:** [`docs/manual/callforge_manual.md`](docs/manual/callforge_manual.md) ·
  [PDF](docs/manual/callforge_manual.pdf) · HPC: [`docs/hpc_deployment.md`](docs/hpc_deployment.md)
  · design notes: `docs/decisions.md`, `docs/reuse.md`.
- **License:** MIT — see [`LICENSE`](LICENSE).
- **Citation:** see [`CITATION.cff`](CITATION.cff).

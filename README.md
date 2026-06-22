# CallForge

**A reproducible Nextflow pipeline for targeted-capture germline analysis — from FASTQ to an annotated joint callset, cohort QC, GIAB benchmarking, and rare-variant burden testing.**

![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)
![Nextflow](https://img.shields.io/badge/Nextflow-DSL2%20%E2%89%A523.10-23aa62.svg)
![CLI: Python](https://img.shields.io/badge/CLI-Python-3776ab.svg)

## Overview

CallForge takes demultiplexed paired-end FASTQs from a hybridization-capture panel through alignment, per-sample QC, SNV/indel joint calling, CNV / STR / paralog-aware analysis, annotation, cohort QC, GIAB benchmarking, and a rare-variant burden framework. It is disease- and organism-agnostic: change inputs and config, never code.

It **runs standalone** — the minimum inputs are FASTQs, a reference, and a target BED. **[CaptureForge](https://github.com/aduton1000/captureforge) metadata is optional enrichment, not a prerequisite**: its per-gene metadata sharpens the interpretation layer (CNV callability, STR loci, paralog flags, burden gene-sets), and `callforge metadata` generates an equivalent table when CaptureForge is not in use. CallForge is the **analysis** counterpart to CaptureForge's **design**.

## Features

- **Four variant classes**: SNV/indel (GATK HaplotypeCaller → GenomicsDB → GenotypeGVCFs, or DeepVariant + GLnexus), CNV (CNVkit), STR (ExpansionHunter), and paralog-aware flagging (`INFO/PARALOG_GENE` / `PARALOG_CONF`).
- **Per-sample QC gate** that flags and quarantines failing samples (depth, duplication, on-target, contamination, sex) rather than silently dropping data.
- **Manifest-driven annotation** (VEP + vcfanno: gnomAD-AFR / ClinVar / dbSNP + PhyloP) with a **genomic-scope gate** that fails loud on a region-limited database, plus chr-reconciliation.
- **GIAB benchmarking** (hap.py) restricted to the panel BED for on-target precision/recall/F1.
- **Rare-variant burden framework** — collapse, regenie, or SKAT-O — with an ancestry-PC stability gate.
- **Reproducible**: per-stage pinned conda environments or containers; `provenance.json` on every run; annotation resources are discovered and validated locally, not downloaded mid-run.
- **The `callforge` CLI** plus **ten per-stage subcommands** for re-running a single stage, and per-process compute controls (`--max_cpus` / `--max_memory`).

## Requirements

- **Nextflow** ≥ 23.10 (Java 17+).
- **conda / mamba** for the per-stage environments, **or** a container engine. VEP and hap.py always run from their official containers.
- **Platforms**: Linux and macOS (Windows via WSL2). On Apple Silicon, per-stage environments build under Rosetta where a tool lacks a native `osx-arm64` build; the report flags any emulated stage.

## Installation

```bash
git clone https://github.com/aduton1000/callforge.git
cd callforge
mamba env create -f env/callforge.yml   # core env; per-stage envs build on demand
pip install -e .                        # registers the `callforge` command
callforge --version
```

A run composes **two** profiles with a comma: an **execution** profile — `local` (this machine: any Linux/macOS, Windows via WSL2) or `hpc_slurm` (a SLURM cluster) — and a **packaging** profile — `conda`, `docker`, `apptainer`, or `singularity` (SingularityCE) — that selects how dependencies are provided. For example, `-profile local,conda` runs on this machine using conda environments; `-profile hpc_slurm,singularity` runs on a SingularityCE cluster. On a cluster, the `callforge-run` wrapper (with a sourced site env file) or an Lmod module reduces a run to a single command for all users (see [`docs/hpc_deployment.md`](docs/hpc_deployment.md)).

> **Advanced/developer.** The raw forms work directly: `nextflow run main.nf …` and `python3 bin/…`. With repository access, the pipeline also runs straight from GitHub once a tag is published — `nextflow run aduton1000/callforge -r <tag> -profile local,conda …` (the repository is currently private, so this form is available to users with access).

## Quickstart

Run the full pipeline end-to-end on the bundled test fixture:

```bash
bash test/make_test_data.sh
callforge run -profile test,docker
```

This builds a small synthetic fixture and runs all stages in minutes, writing the annotated callset, CNV/STR/paralog outputs, the cohort QC dashboard, and a GIAB report to `results/`. The `test` profile uses a **synthetic cohort**: it demonstrates the machinery and plots, not real-data quality or results.

## Usage

**Full pipeline** — `-params-file` supplies the reference, target BED, annotation resources, and sample sheet:

```bash
callforge run \
    -profile local,conda \
    -params-file params.yaml
```

Required parameters: `--input`, `--genome_fasta`, `--target_bed` (see [`params.example.yaml`](params.example.yaml)).

**Helpers:**

```bash
callforge samplesheet --fastq-dir /data/fastqs --metadata clinical.csv --out sample_sheet.csv
callforge metadata    --bed /panel/targets.bed --paralog-genes CD209,CFH --out-tsv gene_metadata.tsv
callforge resources   --dirs /refs --target-bed /panel/targets.bed --out-json resource_manifest.json
```

`samplesheet` pairs FASTQs and joins clinical metadata; `metadata` generates the per-gene table **without CaptureForge** (then run with `--gene_metadata gene_metadata.tsv`); `resources` discovers and scope-checks annotation databases.

**Per-stage subcommands** run one stage standalone on its own inputs, writing to a non-destructive `results/standalone/<stage>_<timestamp>/` with their own report and provenance:

| Subcommand | Purpose |
|:-----------|:--------|
| `align`    | FASTQ → analysis-ready BAM(s) (bwa-mem2 + dedup + BQSR) + alignment QC |
| `coverage` | BAM(s) → coverage/enrichment metrics + per-sample QC gate (pass/quarantine) |
| `call`     | QC-pass BAM(s) → joint hard-filtered VCF (GATK, or `--caller deepvariant`) |
| `cnv`      | BAM(s) → CNVkit calls, labelled by CaptureForge callability |
| `str`      | BAM(s) → ExpansionHunter genotypes on the STR catalog |
| `paralog`  | VCF → paralog-aware flagging (`PARALOG_GENE` / `PARALOG_CONF`) |
| `anno`     | re-annotate a VCF (VEP + vcfanno + PhyloP) |
| `cohortqc` | somalier relatedness/sex + ancestry PCA + missingness |
| `giab`     | hap.py precision/recall/F1 vs GIAB truth, restricted to the panel BED |
| `burden`   | rare-variant burden with new thresholds/engine on an annotated VCF |

For example, re-annotate an existing callset after an annotation-database update — without rerunning upstream stages:

```bash
callforge anno \
    --input_vcf results/stage10_paralog/paralog.annotated.vcf.gz \
    --genome_fasta /refs/GRCh38_noalt.fa \
    --target_bed /panel/targets.bed \
    --resource_dirs /refs/annotation_db_2025_06 \
    -profile local,conda
```

Run `callforge <stage> --help` for a stage's required and optional inputs.

## Inputs and outputs

**Inputs.** A sample sheet (`sample_id,fastq_1,fastq_2` plus optional `sex,phenotype,covariate_*,batch`; `phenotype` is required only for burden) — build it with `callforge samplesheet`. A target BED (a CaptureForge `final_covered_targets.bed`, or any panel BED) and a no-alt reference assembly. Optionally, CaptureForge metadata via `--captureforge_dir`, or an equivalent table from `callforge metadata`. Annotation databases are discovered from `--resource_dirs` (`callforge resources`); `bin/fetch_annotation_dbs.sh` fetches a panel slice.

**Outputs (`results/`).** Annotated joint VCF (`stage11_annotation/annotated.vcf.gz`) and a flattened per-variant TSV (`variants.flat.tsv`); CNV (`cnv_calls.tsv`), STR (`str_calls.tsv`), and paralog (`paralog.annotated.vcf.gz`) calls; a self-contained cohort QC dashboard (`cohort_qc_dashboard.html`); GIAB metrics (`happy.summary.csv`); burden results (`burden_results.tsv`); and `provenance.json` (tool versions, resources, git commit).

## Status and scope

Every metric CallForge currently produces comes from the **synthetic test fixture** and demonstrates the machinery, not real-data quality. GIAB precision/recall, burden statistics, ancestry PCs, relatedness, and coverage numbers from the test profile are **machinery checks on synthetic data** — real benchmarking and association require a real cohort. Burden testing is a framework the pipeline runs; its validity depends on real sample size, phenotype, and ancestry matching. The alternative and heavy engines — DeepVariant + GLnexus, regenie, and SKAT-O — are wired and verified by compile and component tests, and carry a validate-at-deployment checklist ([`docs/hpc_deployment.md`](docs/hpc_deployment.md)) for one live run before being relied on for real data.

## Documentation

The full user manual covers concepts, laptop and HPC installs, every stage and option, the annotation databases, QC interpretation, and worked scenarios:

- **Manual**: [`docs/manual/callforge_manual.md`](docs/manual/callforge_manual.md) · [PDF](docs/manual/callforge_manual.pdf)
- **HPC deployment**: [`docs/hpc_deployment.md`](docs/hpc_deployment.md)
- **Design notes**: [`docs/decisions.md`](docs/decisions.md) · [`docs/reuse.md`](docs/reuse.md)

## CallForge and CaptureForge

CallForge (analysis) and **[CaptureForge](https://github.com/aduton1000/captureforge)** (design) are a pair: CaptureForge defines, tiles, and QCs the capture baits and emits `final_covered_targets.bed` plus per-gene metadata; CallForge consumes that handoff to align, call, annotate, and run cohort QC and rare-variant burden testing. Each runs independently — CallForge needs only FASTQs, a reference, and a target BED.

## Citation, license, and issues

- **Cite**: see [`CITATION.cff`](CITATION.cff).
- **License**: MIT — see [`LICENSE`](LICENSE).
- **Issues / contributions**: open an issue or pull request on the [GitHub repository](https://github.com/aduton1000/callforge).

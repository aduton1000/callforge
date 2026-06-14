# CallForge

**A disease- and organism-agnostic Nextflow DSL2 pipeline that turns demultiplexed
paired-end FASTQs from a hybridization-capture panel into a high-quality,
joint-genotyped, annotated, QC'd callset — SNV/indel, CNV, STR and paralog-aware —
and on to a rare-variant burden/association layer.**

CallForge is the analysis counterpart to **[CaptureForge](../captureforge)** (probe
design) and consumes its handoff: `final_covered_targets.bed` as the target, plus
per-gene metadata (design class, CNV callability, paralog flags, STR loci, burden
groups) that drives the matching downstream steps.

Immediate use: the 59-gene human malaria-susceptibility panel on African-ancestry
DBS samples. It generalizes to any panel/cohort/organism by changing inputs + config.

> **Build status:** Phase 1 complete (scaffold, envs, profiles, resource discovery,
> reference invariant) and green end-to-end on the `test` profile. Phases 2–8
> (align→QC gate→calling→CNV/STR/paralog→annotation→cohort QC→GIAB→burden→reporting)
> are scaffolded with a marked extension point in `subworkflows/callforge.nf`.

---

## Engineering stance

- **Stand on validated foundations.** Subworkflows are adapted from the validated
  nf-core/sarek + nf-core/raredisease patterns (alignment, dedup/BQSR,
  HaplotypeCaller/DeepVariant, VEP, CNV, STR); the bespoke layers (CaptureForge
  handoff, paralog-aware handling, burden, resource auto-discovery, QC dashboard)
  are built here. `docs/reuse.md` tracks reused-vs-built per stage.
- **Reference invariant (inherited from CaptureForge):** align to the **same GRCh38
  no-alt primary assembly** the design used; assert no alt/decoy contigs and a
  matching reference/BED build — **fail loud** (`bin/check_reference_invariant.py`).
- **QC at every stage with plots**, each carrying a one-line plain caption; gates
  **flag & quarantine** failing samples — never silently drop data.
- **Discover, don't download.** Local annotation/known-sites resources are found and
  validated (`bin/discover_resources.py`); nothing is fetched unless missing **and**
  required **and** opted in (`--allow_download`).
- **Genomic-scope gate.** Discovery validates not just build/index but whether each
  annotation DB's data actually **spans the target panel**. A found-but-region-limited
  DB (e.g. a single-locus gnomAD) is flagged `unfit_scope` and **fails loud** — it
  would otherwise silently return "no AF / no ClinVar" outside its region and corrupt
  rare-variant filtering and burden testing. Waive with `--ignore_unfit_resources`
  (not advised); tune with `--scope_min`.

## Quick start

```bash
# 1. environments (Apple Silicon: prefix CONDA_SUBDIR=osx-64 for x86-only tools)
mamba env create -f env/callforge.yml          # core
# specialized envs created on demand: env/{annotate,cohortqc,str,cnv,happy,burden,deepvariant}.yml

# 2. tiny end-to-end sanity run (generates the fixture, runs in minutes)
bash test/make_test_data.sh
nextflow run main.nf -profile test,conda

# 3. full run (point --input at your cohort sample sheet)
nextflow run main.nf -profile mac_local,conda -params-file params.full.yaml   # authoring Mac
nextflow run main.nf -profile hpc_slurm,apptainer -params-file params.full.yaml # production cluster
```

Required params: `--input`, `--genome_fasta` (no-alt), `--target_bed`. See
`params.example.yaml` and `nextflow.config` for the full surface.

## Inputs

| Input | Notes |
|-------|-------|
| **Sample sheet** (CSV) | `sample_id,fastq_1,fastq_2,sex?,phenotype?,covariate_*?,batch?`. `phenotype` is required only for the burden layer; without it the pipeline stops at the annotated callset and skips burden with a clear message. |
| **Target BED** | CaptureForge `final_covered_targets.bed` (as-built coverage). |
| **CaptureForge metadata** | `metrics.json` + `baits.csv` (auto-found under `--captureforge_dir`), or supply `--gene_metadata` / `--paralog_genes` / `--str_catalog` directly. |
| **Reference** | The no-alt GRCh38 primary assembly the panel was designed against. |
| **GIAB** (optional) | A control sample in the sheet + GIAB truth VCF/BED (v4.2.1) for benchmarking restricted to the panel BED. |

## Pipeline stages (0–15)

0. Inputs & resource discovery · 1. Raw-read QC · 2. Trimming · 3. Alignment ·
4. Post-alignment (dedup/BQSR/HsMetrics/coverage) · 5. Per-sample QC gate &
quarantine · 6. SNV/indel joint calling · 7. Hard-filtering · 8. CNV ·
9. STR · 10. Paralog-aware handling · 11. Annotation (VEP + vcfanno) ·
12. Cohort QC (relatedness/ancestry) · 13. Callset QC & GIAB validation ·
14. Burden/association · 15. Reporting (MultiQC + cohort dashboard).

Each stage writes metrics (JSON/TSV), generates captioned plots (PNG+SVG), and
feeds MultiQC + a self-contained cohort QC dashboard.

## Profiles

- `test` — 3 samples, a few small genes, a subset reference; minutes on a Mac.
- `mac_local` — full pipeline on one macOS machine (heavy steps slow; see Apple-Silicon note).
- `hpc_slurm` — SLURM + Apptainer (or shared conda at `/hpc/opt/conda`); see `docs/hpc_deployment.md`.

Compose with a container/engine profile: `-profile test,conda` / `,docker` / `,apptainer`.

### Apple Silicon (arm64)

Some bioconda tools (gatk4, bwa-mem2, deepvariant, expansionhunter) lack
`osx-arm64` builds. Create envs under Rosetta so x86_64 builds resolve:

```bash
CONDA_SUBDIR=osx-64 mamba env create -f env/callforge.yml
```

The run report flags which stages executed x86-emulated. Linux/HPC is native.

## Repository layout

```
main.nf  nextflow.config  params.example.yaml  params.full.yaml
conf/{base,test,mac_local,hpc_slurm}.config
modules/        # one .nf per process group (stage0_inputs.nf, ...)
subworkflows/   # callforge.nf — top-level orchestration
bin/            # discover_resources.py, check_reference_invariant.py, ingest_captureforge.py, ...
env/            # pinned conda envs + Dockerfile + Apptainer def
assets/  docs/  test/   results/(gitignored)
```

## License

MIT — see `LICENSE`. Cite via `CITATION.cff`.

# CallForge — operations manual

End-to-end guide for running CallForge (FASTQ → annotated joint callset → burden).
Pairs with `README.md` (overview), `docs/decisions.md` (design decisions/caveats),
`docs/hpc_deployment.md` (cluster), and `docs/reuse.md` (reused vs built).

## 1. Environments
The pipeline runs each stage in its own pinned conda env or container (`env/*.yml`,
`vep_image`, `happy_image`). Create them via `-profile conda` (Nextflow builds each
`env/*.yml` on first use) or build the core container (`env/Dockerfile` / `env/callforge.def`).

- **Linux / HPC**: native; `-profile hpc_slurm,apptainer` (or `,conda`).
- **macOS / Apple Silicon**: `CONDA_SUBDIR=osx-64 mamba env create -f env/<x>.yml`
  (Rosetta) where a tool lacks an `osx-arm64` build; the report flags emulated stages.
  VEP and hap.py run from their official **containers** (robust across platforms).

## 2. Profiles
| profile | use |
|---|---|
| `test` | 3–4 synthetic samples, subset reference; minutes. Add `,docker` for VEP/hap.py. |
| `mac_local` | full pipeline on one macOS machine (heavy steps slow). |
| `hpc_slurm` | SLURM + Apptainer / shared conda; the production cohort run. |

Compose an engine: `-profile test,docker`, `mac_local,conda`, `hpc_slurm,apptainer`.

## 3. Inputs (see `params.example.yaml`)
- `--input` sample sheet CSV: `sample_id,fastq_1,fastq_2,sex?,phenotype?,covariate_*?,batch?`.
  `phenotype` (case/control) is required only for burden (else it stops at the annotated callset).
- `--genome_fasta` the **no-alt GRCh38 primary assembly the panel was designed against**.
- `--target_bed` CaptureForge `final_covered_targets.bed`; `--captureforge_dir` for `metrics.json`+`baits.csv`.
- `--resource_dirs` to scan for annotation/known-sites resources (discovered, not downloaded).

## 4. Run
```bash
bash test/make_test_data.sh && nextflow run main.nf -profile test,docker          # sanity
nextflow run main.nf -profile mac_local,conda     -params-file params.full.yaml    # full (Mac)
nextflow run main.nf -profile hpc_slurm,apptainer -params-file params.full.yaml    # full (HPC)
```

## 5. Outputs (`results/`)
`stage0_inputs/` (manifest, reference invariant, gene metadata, validated sheet) ·
`stage1_rawqc` `stage2_trim` `stage3_align` `stage4_postalign` (HsMetrics, mosdepth) ·
`stage5_qc_gate/` (scorecard, **qc_pass.txt / qc_quarantine.txt**) ·
`stage6_calling/` (`joint.vcf.gz`, stats) · `stage7_filter/` (`joint.filtered.vcf.gz`) ·
`stage8_cnv` (calls + callability) · `stage9_str` · `stage10_paralog/` (`paralog.annotated.vcf.gz`) ·
`stage11_annotation/` (**`annotated.vcf.gz`**, **`variants.flat.tsv`**, `annotation_landing.json`) ·
`stage12_cohortqc` · `stage13_giab` · `stage14_burden/` (`burden_results.tsv` — engine-stamped) ·
**`stage15_report/`** (`cohort_qc_dashboard.html`, `multiqc_report.html`, `provenance.json`).

Every stage emits captioned plots (PNG+SVG) collated into the dashboard.

## 6. Key behaviours
- **Reference invariant**: fails loud on alt/decoy contigs or a BED/reference build mismatch.
- **Resource scope gate**: an annotation DB that doesn't span the panel (or is a truncated
  `.bgz`) is flagged unfit and **fails loud** — never silently under-annotates.
- **QC gate**: failing samples are **quarantined** (excluded from joint calling AND burden
  counts), never dropped; reasons recorded.
- **Annotation landing** is verified (fail-loud) after chr-reconciliation.
- **Ancestry-PC gate**: unstable PCs are dropped and burden runs unadjusted (stated in provenance).

## 7. Real-run checklist (before the production cohort)
- [ ] Demultiplexed cohort FASTQs + sample sheet (with `phenotype` for burden).
- [ ] Genome-wide **gnomAD** (AFR), **ClinVar**, **dbSNP** (also a BQSR known-site), **PhyloP**
      — all must pass the scope gate. dbSNP not yet fetched locally.
- [ ] BQSR known-sites (dbSNP + Mills + 1000G indels) via `--known_sites`.
- [ ] **1000G GRCh38 somalier sites** via `--somalier_sites`.
- [ ] GIAB control in the sheet + truth v4.2.1 (`--giab_control_id`, `--giab_truth_vcf/bed`).
- [ ] Confirm **ancestry PCs are stable** before using as burden covariates (gate enforces this).
- [ ] **Validate before production**: DeepVariant+GLnexus, regenie, SKAT-O (all wired, not yet executed).
- [ ] Re-generate the genome-wide gnomAD panel slice on the HPC (fast network).

## 8. Troubleshooting
- *VEP "CacheDir" error in gtf mode* — don't pass `--offline` with `--gtf` (cache mode only).
- *mosdepth "exec format error" on arm64* — use mosdepth ≥0.3.11 or `CONDA_SUBDIR=osx-64`.
- *Scope gate fails on a region-subset DB* — supply a genome-wide replacement (the message says so).
- *No samples reach calling* — check `qc_quarantine.txt`; relax thresholds only with justification.

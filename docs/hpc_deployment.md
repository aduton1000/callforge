# CallForge — HPC deployment (Ubuntu 24.04 SLURM cluster)

CallForge is authored on macOS (Apple Silicon) but the full cohort run is intended
for the cluster, where the toolchain is **native linux-64** (no Rosetta).

## 1. Conda envs (shared)

The cluster's conda lives at `/hpc/opt/conda`:

```bash
source /hpc/opt/conda/etc/profile.d/conda.sh
# Build each env once into the shared env tree:
mamba env create -p /hpc/opt/conda/envs/callforge          -f env/callforge.yml
mamba env create -p /hpc/opt/conda/envs/callforge-annotate -f env/annotate.yml
# ... and cohortqc / str / cnv / happy / burden / deepvariant as needed.
```

`conf/hpc_slurm.config` sets `conda.cacheDir` so per-process envs are built once and
reused cluster-wide. Run with `-profile hpc_slurm,conda`.

## 2. Containers (Apptainer)

Docker is usually forbidden on shared clusters. Build the Apptainer image:

```bash
apptainer build callforge.sif env/callforge.def
```

Run with `-profile hpc_slurm,apptainer`. `autoMounts=true` binds the work dir;
**bind-mount the resource directories** so resource discovery and annotation can see
them inside the container, e.g.:

```bash
export APPTAINER_BINDPATH=/data/refs,/data/annotation_db,$HOME/.vep
```

## 3. SLURM

`conf/hpc_slurm.config` uses the SLURM executor (global, cluster-wide across nodes).
Override per-cluster knobs on the command line or in the params file:

```bash
callforge run -profile hpc_slurm,apptainer -params-file params.full.yaml \
  --slurm_partition compute --slurm_account mylab --scratch_dir /scratch/$USER
```

Heavy steps (genome `bwa-mem2 index`, joint genotyping, gCNV) carry the `index`/
`align`/`large` labels and get bigger slots (see the `withLabel` blocks).

## 4. Resources on HPC

Copy the discovered local resources (VEP cache, gnomAD-AFR, dbSNP, ClinVar, PhyloP,
BQSR known-sites) to a shared path and point `--resource_dirs` at it, **or** pass the
explicit `--vep_cache/--gnomad_vcf/...` overrides. The Stage-0 manifest records what
was found vs missing so a cluster run is reproducible. Lustre/GPFS scratch is
recommended for `work/` (`--scratch_dir`).

## 5. Reproducibility

Every run emits `provenance.json` (tool versions, reference + resource releases/md5,
git commit, quarantined samples + reasons). Pin Nextflow with `nextflowVersion` in
`manifest{}` and the conda/container versions in `env/*.yml`.

## 6. Validate-at-deployment checklist (stage subcommands)

The stage subcommands (`callforge <stage>`, see the manual) were built and verified by
`-preview` compile + component/CLI tests; the pure-Python paths were also run live
locally (`burden`, `paralog`). Before any heavy standalone path is relied on for real
data, run it **end-to-end once** here on the cluster (where the conda envs / Apptainer
images exist natively). Preview-compile + component tests are sufficient for development
gates; they are **not** sufficient to call a heavy path production-ready.

- [ ] `anno` — the DB-swap path (VEP container + vcfanno/`callforge-annotate`).
- [ ] `align` → `coverage` → `call` — the spine (bwa-mem2 / picard / mosdepth / GATK).
- [ ] `cnv` (CNVkit / `callforge-cnv`) and `str` (ExpansionHunter / `callforge-str`).
- [ ] Alternative engines: DeepVariant + GLnexus (`--caller deepvariant`),
      regenie and SKAT-O (`--burden_engine regenie|skat`).

Live-verified locally so far: full pipeline (`-profile test`), `callforge burden`
(collapse), `callforge paralog` (PARALOG_GENE/PARALOG_CONF written).

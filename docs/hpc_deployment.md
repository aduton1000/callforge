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
nextflow run main.nf -profile hpc_slurm,apptainer -params-file params.full.yaml \
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

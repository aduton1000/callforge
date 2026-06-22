# CallForge — HPC deployment (SLURM + Apptainer/Singularity)

A cohort-scale germline run is best suited to the HPC: alignment, joint genotyping, and
DeepVariant/GLnexus are heavy, and the toolchain is **native linux-64** (no emulation).
This is a single central install other users invoke. All `/hpc/...`, `/shared/...` paths
below are **examples** — substitute your site's layout.

## 0. One-time: central install

```bash
# clone into a shared location (read-execute for the group)
cd /shared/apps
git clone <your-private-remote>/callforge.git
cd callforge
git checkout v0.1.0     # pin a tagged release
```

> **Shared, read-only conda?** If you install the `callforge` CLI (or build an env)
> against a cluster's shared conda whose package cache is not group-writable, conda
> fails with *"cannot set permissions"* on the pkgs cache. Point the cache at a writable
> dir first:
> ```bash
> export CONDA_PKGS_DIRS="$HOME/.conda/pkgs"
> ```
> On a container-only deployment (below) you don't need a host conda env at all — the
> tools come from the `.sif` and the per-stage official images.

## 1. Build the container image (native linux-64 toolchain; Docker usually banned on HPC)

The same `env/callforge.def` builds one core `.sif` that runs under **either** Apptainer
**or** SingularityCE — use whichever engine your cluster provides:

```bash
cd /shared/apps/callforge

# Apptainer sites:
apptainer build callforge.sif env/callforge.def

# SingularityCE sites (no `apptainer` binary):
sudo singularity build callforge.sif env/callforge.def
# unprivileged alternative (needs fakeroot configured — a /etc/subuid + /etc/subgid
# entry for your user; ask your admin if `--fakeroot` errors):
#   singularity build --fakeroot callforge.sif env/callforge.def

# → contains bwa-mem2, samtools/bcftools/htslib, gatk4, picard, bedtools, mosdepth,
#   multiqc, the Python plotting/QC layer, AND procps (`ps`) + base shell utils —
#   Nextflow calls `ps` to poll tasks, so a stripped image makes every task fail.
#   Build from the provided `.def` (don't trim it).
```

> **Stage-specific images are pulled separately.** A few stages pin their own official
> images and are **not** part of `callforge.sif` — VEP (`ensemblorg/ensembl-vep`), hap.py
> (`jmcdani20/hap.py`), and DeepVariant (`google/deepvariant`, only with
> `--caller deepvariant`). With a container engine enabled, Nextflow pulls these on first
> use into the shared `container_cache` (see §4). The annotation/STR/CNV/cohort-QC/burden
> conda envs (`env/annotate.yml`, `env/str.yml`, `env/cnv.yml`, `env/cohortqc.yml`,
> `env/happy.yml`, `env/burden.yml`, `env/deepvariant.yml`) are used under
> `-profile …,conda`. See the **verify-at-deployment checklist** (§7) for what must be
> exercised live before a heavy path is relied on.

> **Verify a tool is in the image** with a NON-login shell — `singularity exec <img>
> <tool> --version`. A login shell (`bash -lc …`) re-sources profiles and can reset PATH,
> hiding the image's tools and giving a misleading "not found":
> ```bash
> singularity exec callforge.sif bwa-mem2 version
> singularity exec callforge.sif gatk --version
> singularity exec callforge.sif ps --version          # procps present (Nextflow needs it)
> singularity exec callforge.sif python -c 'import matplotlib; print(matplotlib.__version__)'
> ```

## 2. Centralized reference + resource data (stage once, shared read-only)

```bash
export CALLFORGE_REFS=/shared/refs/callforge/GRCh38
mkdir -p "$CALLFORGE_REFS"
```

Stage **once** (these do not change between runs):

- **Reference genome** — the **same no-alt primary assembly** CaptureForge designed
  against (`--genome_fasta`). The `.fai` and `.dict` are built on demand if absent.
- **Annotation databases** — VEP cache (`--vep_cache`), gnomAD (prefer one carrying AFR
  frequencies, `--gnomad_vcf`), dbSNP (`--dbsnp_vcf`), ClinVar (`--clinvar_vcf`), PhyloP
  bigWig (`--phylop_bw`). `bin/fetch_annotation_dbs.sh` fetches a panel slice; or point
  `--resource_dirs` at the directory holding them and let Stage 0 discover them.
- **BQSR known-sites** (`--known_sites`: dbSNP + Mills + 1000G indels) — else BQSR is
  skipped.
- **Optional validation refs** — somalier sites (`--somalier_sites`, cohort QC), GIAB
  truth VCF/BED (`--giab_truth_vcf` / `--giab_truth_bed`, benchmarking).

The **target BED** (`--target_bed`) is per-panel: a CaptureForge `final_covered_targets.bed`
or any panel BED. CaptureForge per-gene metadata (`--captureforge_dir`) is optional
enrichment; without it, build an equivalent table with `callforge metadata`.

> **Contig-naming consistency (lesson learned — verify at deployment).** CallForge
> enforces, fail-loud, that the reference is **no-alt** and that the target BED's build
> matches it: every BED contig exists in the reference `.fai`, the chr-prefix style
> (`chr1` vs `1`) is **consistent** between BED and reference, and no interval runs off a
> contig end (`bin/check_reference_invariant.py`, exit 2 on violation). Because CallForge
> also consumes BAM/VCF inputs (standalone stages) and external annotation databases, the
> **same chr-prefix convention must hold across the reference, the target BED, the input
> BAMs/VCFs, and every annotation DB** — a gnomAD VCF named `1..22` against a `chr1`
> reference silently under-annotates. Confirm the convention across all of these at
> deployment; the genomic-scope gate (`--scope_min`) additionally fails loud on a
> region-limited DB.

## 3. Resolve resources (discovery, or explicit paths)

```bash
callforge resources --dirs "$CALLFORGE_REFS" --target-bed "$CALLFORGE_REFS"/targets.bed \
    --out-json "$CALLFORGE_REFS"/resource_manifest.json
```

The Stage-0 manifest records what was found vs missing (and whether each DB spans the
target contigs) so a cluster run is reproducible. Pass `--resource_manifest` to skip
re-discovery, or give the explicit `--vep_cache/--gnomad_vcf/--dbsnp_vcf/...` overrides.

## 4. Run on SLURM

Compose the packaging profile that matches your engine: **`hpc_slurm,singularity`** on a
SingularityCE cluster, or **`hpc_slurm,apptainer`** on an Apptainer one (both run the same
`.sif`); or **`hpc_slurm,conda`** for the shared-conda route.

```bash
cd /scratch/$USER/run1
callforge run \
    --pipeline /shared/apps/callforge/main.nf \
    -profile hpc_slurm,singularity \
    -params-file my.params.yaml \
    --input         samplesheet.csv \
    --genome_fasta  $CALLFORGE_REFS/GRCh38_noalt_primary.fa \
    --target_bed    $CALLFORGE_REFS/final_covered_targets.bed \
    --resource_dirs $CALLFORGE_REFS/annotation_db \
    --slurm_partition <queue> --scratch_dir $TMPDIR \
    --container_image /shared/apps/callforge/callforge.sif \
    --outdir results
# refs outside $HOME/CWD (e.g. annotation DBs / VEP cache under /hpc) must be bound into
# task containers — the reliable route is an exported env var (see callout below):
export SINGULARITY_BIND=/hpc      # APPTAINER_BIND=/hpc on an Apptainer site; [cluster-specific]
```

> **`--container_image` is required when an engine is enabled.** It sets
> `process.container`, so every task runs via `singularity exec <image> …` (the `.sif`
> from §1). Without it, Nextflow turns the engine on but wires no image to processes —
> tasks then run on the **host** (silently, until a stage needs a tool the host lacks,
> e.g. `gatk` or `matplotlib`). The pipeline now **fails loudly** in that case. For
> singularity/apptainer use a local `.sif` path; for docker use the image tag. Put it in
> the params file as `container_image:` so it's not retyped each run. (VEP, hap.py, and
> DeepVariant pin their own images and are unaffected.)

> **Container bind path.** Singularity/Apptainer only auto-mount `$HOME`, `/tmp`, and the
> CWD. Annotation databases, the VEP cache, and known-sites are **referenced in place —
> not staged into work dirs** — so if your references live outside those (e.g. under
> `/hpc/...`) you must bind that root into the task containers, or they fail with "No such
> file or directory".
>
> - **Primary — `export SINGULARITY_BIND=/hpc`** (or `APPTAINER_BIND` on an Apptainer
>   site). Nextflow scrubs the task environment, but the `singularity`/`apptainer`
>   profiles **whitelist** these vars (`envWhitelist`) so they survive the scrub and the
>   engine honors them. This is the most robust route.
> - **Alternative — `--bind_paths /hpc`** (or `bind_paths:` in the params file). Applied
>   per task via `process.containerOptions` (a closure resolved after the params merge —
>   *not* engine `runOptions`, which would evaluate before the value merges and never
>   bind).
>
> Either works; you don't need both. Default is empty (no bind). The `callforge-run`
> wrapper sets both from `CALLFORGE_BIND`.

Tune in `conf/hpc_slurm.config`: `slurm_partition` (**required** for SLURM runs — there is
no safe default; discover yours with `sinfo -s` and set `--slurm_partition <queue>` or
`slurm_partition: <queue>` in the params file. If unset the run **fails loudly** before
any job is submitted, rather than silently submitting to a non-existent partition),
`slurm_account`, `slurm_queue_opts`, `scratch_dir` (Lustre/GPFS scratch recommended for
`work/`), `max_cpus/memory/time`, and `container_cache` (the shared image cache;
`apptainer_cache` remains a back-compat alias). Heavy steps carry the `index`/`align`/
`large` labels and get bigger slots (16 cpu / 48–64 GB).

## 5. Shared invocation: wrapper on PATH (no module system) or Lmod

`share/bin/callforge-run` is a thin launcher that forwards to
`nextflow run … -profile hpc_slurm,<engine>`. It defaults to the **singularity** engine
and resolves the install root, `.sif`, refs, and bind from environment variables (with
install-relative defaults) — so it needs **no module system**.

**No-module clusters — source the site env once, then run short.** Copy
`share/callforge-env.sh.example` to a shared path, fill in the values once, and have users
source it from `~/.bashrc`. After that, `callforge-run` is on PATH and the bind/image/refs
are configured — no per-run exports:

```bash
# one-time, per site (admin):
cp share/callforge-env.sh.example /shared/apps/callforge/callforge-env.sh
$EDITOR /shared/apps/callforge/callforge-env.sh      # set HOME/SIF/REFS/ENGINE/BIND

# one-time, per user:
echo 'source /shared/apps/callforge/callforge-env.sh' >> ~/.bashrc

# every run thereafter (queue + refs + image all come from the params file / env):
callforge-run -params-file site.params.yaml --input samplesheet.csv --outdir results
```

Start `site.params.yaml` from the annotated template
**`conf/cluster.params.example.yaml`** (copy and edit once): it carries `genome_fasta`,
`target_bed`, the annotation resources, `bind_paths`, `container_image`, and
`slurm_partition` — so they are not retyped per run. The wrapper still accepts overrides
(`--slurm_partition <queue>`, `--outdir …`, any `nextflow` flag).

**Lmod sites:** install the modulefile stub so users `module load callforge` (it sets the
same env vars), then run the wrapper:

```bash
module use /shared/apps/callforge/share/modulefiles
module load callforge/0.1.0
callforge-run -params-file site.params.yaml --input samplesheet.csv --slurm_partition <queue>
```

## 6. Reproducibility

Every run emits `provenance.json` (tool versions, reference + resource releases/md5, git
commit, quarantined samples + reasons). Pin Nextflow with `nextflowVersion` in
`manifest{}` and the conda/container versions in `env/*.yml` / `env/callforge.def`.

## 7. Validate-at-deployment checklist (heavy paths)

The pipeline spine and the stage subcommands compile-preview and pass component/CLI tests,
and the pure-Python paths (`burden` collapse, `paralog`) have been run live locally.
**Before any heavy path is relied on for real data, run it end-to-end once on the cluster**
(where the conda envs / container images exist natively). Preview-compile + component tests
gate development; they do **not** certify a heavy path production-ready.

- [ ] `align` → `coverage` → `call` — the spine (bwa-mem2 / picard / mosdepth / GATK).
- [ ] `anno` — VEP container + vcfanno/PhyloP, with the real annotation DBs and bind.
- [ ] `cnv` (CNVkit) and `str` (ExpansionHunter) on real BAMs.
- [ ] `cohortqc` (somalier) and `giab` (hap.py container) on a real control.
- [ ] Alternative/heavy engines: DeepVariant + GLnexus (`--caller deepvariant`), regenie
      and SKAT-O (`--burden_engine regenie|skat`) — container pull + resource sizing.
- [ ] Contig-naming consistency across reference / target BED / input BAMs+VCFs / every
      annotation DB (§2), and the genomic-scope gate against the real gnomAD release.

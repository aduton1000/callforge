#!/usr/bin/env python3
"""callforge — unified command-line front door for the CallForge pipeline.

Subcommands:
  callforge run          run the full pipeline (thin wrapper over `nextflow run`)
  callforge samplesheet  scaffold/assemble a sample sheet (bin/init_sample_sheet.py)
  callforge resources    resource discovery / scope checks (bin/discover_resources.py)
  callforge --version    print the CallForge version
  callforge --help       list subcommands

Design notes
------------
* This is a FRIENDLY FRONT DOOR, not a re-architecture. The raw forms keep working
  unchanged: `python3 bin/init_sample_sheet.py …` and `nextflow run main.nf …`.
* It does NOT touch the per-stage conda environments. CallForge deliberately uses
  multiple small, isolated per-stage envs (callforge, callforge-cnv, callforge-str,
  callforge-annotate, callforge-cohortqc, …). With `-profile conda`/`apptainer`,
  Nextflow activates the correct env/container per process automatically — the user
  never activates them by hand. `callforge run` only launches Nextflow.

Per-stage subcommands
---------------------
The ten per-stage subcommands (align, coverage, call, cnv, str, paralog, anno,
cohortqc, giab, burden) run ONE analytical stage standalone. Each maps to
`nextflow run main.nf --stage <name> …` (a `--stage` param, not `-entry`: portable
across Nextflow versions; the strict parser disallows `-entry` for named workflows),
reusing the SAME modules as the full pipeline. See the STAGES table below.
"""
import argparse
import importlib.util
import os
import shlex
import subprocess
import sys

from . import __version__

# repo root = parent of this package dir (works for `pip install -e .` / in-repo use)
PKG_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(PKG_DIR)
BIN = os.path.join(REPO_ROOT, "bin")


# ----------------------------------------------------------------- helpers
def _load_bin(script):
    """Import a bin/<script>.py file as a module object."""
    path = os.path.join(BIN, script)
    if not os.path.isfile(path):
        sys.stderr.write(f"[callforge] ERROR: cannot find {path}\n")
        sys.exit(2)
    spec = importlib.util.spec_from_file_location(script.replace(".py", ""), path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _run_bin_via_argv(script, rest):
    """Invoke a bin script's main() by temporarily setting sys.argv (no edits to it)."""
    mod = _load_bin(script)
    old = sys.argv
    sys.argv = [os.path.join(BIN, script), *rest]
    try:
        return mod.main() or 0
    finally:
        sys.argv = old


def _run_state(passthrough):
    """Pin Nextflow's run state to a USER-WRITABLE dir, never the read-only shared
    install. Nextflow anchors .nextflow/ (history, history.lock, cache) and the log to
    launchDir == the process CWD, and -work-dir moves ONLY work/ — not those. So we
    launch from the user's CWD (never cd into the install) and, belt-and-suspenders,
    point the log + work dir + the launcher's own scratch (NXF_HOME) at user space.

    Returns (rundir, core_flags, run_flags, env): core_flags go BEFORE `run` (-log is a
    core option), run_flags after it (-work-dir is a run option)."""
    rundir = os.environ.get("CALLFORGE_RUNDIR") or os.getcwd()
    core, run = [], []
    has_log = any(a == "-log" or a.startswith("-log=") for a in passthrough)
    has_work = any(a in ("-work-dir", "-w") or a.startswith("-work-dir=") for a in passthrough)
    if not has_log:
        core += ["-log", os.path.join(rundir, ".nextflow.log")]
    if not has_work and not os.environ.get("NXF_WORK"):
        run += ["-work-dir", os.path.join(rundir, "work")]
    env = os.environ.copy()
    env.setdefault("NXF_HOME", os.path.join(os.path.expanduser("~"), ".nextflow"))
    return rundir, core, run, env


def _require_writable(rundir):
    """Refuse to run from a non-writable dir with a clear message — instead of letting
    Nextflow emit the cryptic '.nextflow/history.lock (Permission denied)'."""
    if not os.access(rundir, os.W_OK):
        sys.stderr.write(
            f"[callforge] ERROR: run directory is not writable: {rundir}\n"
            "  CallForge writes run state (.nextflow/, logs, work/) to this directory.\n"
            "  Run from a writable directory (your home or a project dir), e.g.:\n"
            "    cd ~ && callforge run -params-file site.params.yaml --input samplesheet.csv --outdir results\n"
            "  (or set CALLFORGE_RUNDIR=/path/to/writable to override).\n")
        sys.exit(1)


def _warn_env_if_local(passthrough):
    """Backstop: warn (do NOT auto-activate) if the `callforge` conda env is not on
    PATH for a LOCAL (non-container) run. The user stays in control.

    Suppressed for ANY container engine (the tools come from the image / per-stage
    containers, not a host conda env) — including singularity, so an in-container
    cluster run never emits a misleading 'conda activate callforge'. Also suppressed
    for the test profile (its packaging profile is chosen explicitly)."""
    joined = " ".join(passthrough)
    engines = ("docker", "apptainer", "singularity", "podman", "shifter", "charliecloud", "sarus")
    if any(e in joined for e in engines) or "--container_image" in passthrough:
        return                                  # a container provides the tools
    env = os.environ.get("CONDA_DEFAULT_ENV", "")
    prefix = os.environ.get("CONDA_PREFIX", "")
    if env.startswith("callforge") or prefix.rstrip("/").rsplit("/", 1)[-1].startswith("callforge"):
        return
    sys.stderr.write(
        "[callforge] WARN: the `callforge` conda env was not detected. Local "
        "(non-container) runs — especially Apple-Silicon arm64 — need its tools on PATH:\n"
        "    conda activate callforge\n"
        "  (Nextflow still activates the correct per-stage env per process with -profile conda.)\n"
        "  Proceeding anyway (use a container engine, e.g. -profile <exec>,docker, to avoid this).\n")


def _resolve_pipeline(override):
    """Resolve the Nextflow entry. Default: repo's main.nf. Override may be a path,
    a project dir, or a GitHub `owner/repo` name — passed to Nextflow verbatim."""
    if override:
        return override
    cand = os.path.join(REPO_ROOT, "main.nf")
    return cand if os.path.isfile(cand) else "main.nf"


# Pipeline parameters are declared in Nextflow (nextflow.config / conf/), NOT in
# argparse — `callforge run` forwards them verbatim. argparse therefore can't own
# their --help, so we document them here (printed as the `run` epilog) WITH provenance
# so a reader knows, per flag, whether they set it or it is already wired.
#   [install]          once in the site env file (CALLFORGE_HOME/_SIF/_REFS/_BIND/_ENGINE)
#   [site]             once per cluster in the params file (-params-file)
#   [per-run]          you provide every run (the minimal set)
#   [optional]         has a sane default; set only to override (default shown)
#   [stage-dependent]  only needed when that stage/engine/branch runs
# Defaults below are the ACTUAL nextflow.config / conf defaults (no invented values).
PASSTHROUGH_HELP = """\
Nextflow pipeline parameters (forwarded to `nextflow run`; declared in
nextflow.config / conf/, so they are passthrough — not listed under "options" above).
Provenance: [install]=env file once · [site]=params file once per cluster ·
[per-run]=every run · [optional]=defaulted override · [stage-dependent]=only when that stage runs.

  -profile NAMES        [per-run] execution+packaging, composed: test|local|hpc_slurm
                        with conda|docker|apptainer|singularity (e.g. local,conda)
  -params-file FILE     [site] YAML carrying the site values (genome, target BED, DBs, …)

 inputs (REQUIRED)
  --input PATH          [per-run] REQUIRED sample sheet CSV
                        (sample_id,fastq_1,fastq_2[,sex,phenotype,covariate_*,batch])
  --genome_fasta PATH   [site] REQUIRED no-alt primary-assembly FASTA (same build/contig
                        naming as the target BED; .fai/.dict built if absent)
  --target_bed PATH     [site] REQUIRED panel BED (CaptureForge final_covered_targets.bed
                        or any panel BED)

 CaptureForge handoff (optional enrichment — CNV callability / STR / paralog / burden groups)
  --captureforge_dir PATH [optional] results/<run> dir (auto-finds metrics.json + baits.csv)
  --cf_metrics_json PATH  [optional] explicit CaptureForge metrics.json
  --gene_metadata PATH    [optional] pre-built per-gene metadata TSV (use WITHOUT CaptureForge;
                          build it with `callforge metadata`)
  --paralog_genes LIST    [optional] comma list (e.g. CD209,CASP1,CR1,HP,CFH)
  --str_catalog PATH      [optional] ExpansionHunter catalog JSON (else built from STR targets)

 reference / annotation resources (referenced in place → cover their root with --bind_paths)
  --resource_dirs LIST  [site] comma list scanned by discovery (default $HOME,/usr/local/share,/opt)
  --allow_download BOOL [optional] fetch a missing required DB (default false)
  --vep_cache PATH      [site] VEP cache dir (else discovered)        --vep_release N (110)
  --gnomad_vcf PATH     [site] gnomAD VCF (prefer AFR)               --gnomad_af_field X (AF_afr)
  --dbsnp_vcf PATH      [site] dbSNP VCF                              --clinvar_vcf PATH [site]
  --phylop_bw PATH      [site] PhyloP bigWig                         --known_sites LIST [site] (BQSR)
  --vcfanno_toml PATH   [site] bring-your-own vcfanno DBs (organism-generic)
  --scope_min F         [optional] min fraction of target contigs a DB must span (0.9);
                        a region-limited DB is flagged UNFIT and fails loud
  --species NAME        [optional] VEP species (default homo_sapiens)
  --genome_build STR    [optional] label (default GRCh38)

 container / HPC
  --container_image X   [site] image every task runs in (.sif path | docker tag).
                        REQUIRED when a container engine is enabled (else fails loudly);
                        the callforge-run wrapper supplies it from CALLFORGE_SIF
  --bind_paths PATHS    [site] host path(s) to bind into containers for in-place refs
                        (genome / annotation DBs / VEP cache). Default: none. Comma-separated
  --slurm_partition Q   [site] SLURM queue. REQUIRED under -profile hpc_slurm (no safe
                        default; fails loudly if unset). Find one with `sinfo -s`
  --slurm_account S     [optional] SLURM account (--account=...)
  --slurm_queue_opts S  [optional] extra sbatch/clusterOptions
  --scratch_dir PATH    [optional] per-task scratch under SLURM (default /tmp)

 analysis options (all [optional]/[stage-dependent]; defaults shown)
  --caller X (gatk: gatk|deepvariant)        --ploidy N (2)   --joint_method X (genomicsdb)
  --cnv_caller X (cnvkit)  --cnv_method X (hybrid)  --cnv_segment_method X (cbs)
  --cnv_min_samples_pon N (5)   --str_motif X (GT)   --paralog_min_mq N (50)
  --vep_mode X (cache|gtf)  --gtf PATH (for vep_mode=gtf)  --vep_image X (ensemblorg/ensembl-vep:release_112.0)
  --min_mean_target_depth N (30)  --max_dup_rate F (0.40)  --min_on_target F (0.40)
  --max_contamination F (0.03)  --enforce_sex_check BOOL (true)
  --run_burden BOOL (true; auto-skip if no phenotype)  --burden_engine X (regenie|skat|collapse)
  --burden_af_max F (0.01)  --burden_csq LIST  --n_ancestry_pcs N (4)  --covariates LIST
  --giab_control_id ID  --giab_truth_vcf PATH  --giab_truth_bed PATH  --happy_image X
  --somalier_sites PATH [stage-dependent] (cohortqc)
  --max_cpus N  --max_memory X  --max_time X  (per-profile compute ceilings)

 output
  --panel_name STR      [per-run/optional] label on outputs (default callforge_panel)
  --outdir PATH         [per-run] results directory (default results)

See `params.example.yaml`, `conf/cluster.params.example.yaml`, and the manual's
Parameter Reference for the full table. Per-stage conda envs/containers are activated
automatically by Nextflow; you never activate them by hand."""


# ----------------------------------------------------------------- subcommands
def cmd_run(rest):
    ap = argparse.ArgumentParser(
        prog="callforge run", allow_abbrev=False, add_help=True,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "Run the full CallForge pipeline (FASTQ -> annotated joint callset -> burden). "
            "Thin wrapper over `nextflow run`; -profile, -params-file, -resume and every "
            "--<param> pass straight through to Nextflow.\n\n"
            "Most options are set ONCE (env file + params file), so a typical run is short:\n"
            "    callforge run -profile local,conda -params-file params.yaml \\\n"
            "        --input samplesheet.csv --outdir results\n"
            "On a deployed cluster the wrapper makes it shorter still:\n"
            "    callforge-run -params-file site.params.yaml --input samplesheet.csv --outdir results\n"
            "The flags below OVERRIDE params-file values; you rarely type them all."),
        epilog=PASSTHROUGH_HELP)
    ap.add_argument("--pipeline", metavar="PATH|owner/repo",
                    help="[optional] Nextflow entry to run (default: the repo's main.nf). "
                         "May be a path, a project directory, or a GitHub owner/repo.")
    ap.add_argument("--print-cmd", action="store_true",
                    help="[optional] print the `nextflow run …` command and exit (do not "
                         "execute) — same as setting CALLFORGE_DRY_RUN")
    known, passthrough = ap.parse_known_args(rest)
    _warn_env_if_local(passthrough)
    return _launch_nextflow(known.pipeline, [], passthrough, known.print_cmd)


def _launch_nextflow(pipeline_override, prefix_args, passthrough, print_cmd):
    """Build + run (or print) `nextflow run <pipeline> <prefix_args> <passthrough>`.
    prefix_args are CLI-injected flags (e.g. --stage anno); passthrough is everything
    the user typed (so -profile/-params-file/-resume/--<param> reach Nextflow)."""
    pipeline = _resolve_pipeline(pipeline_override)
    # Run from the user's CWD (launchDir), NOT the install — so .nextflow/ + logs +
    # work/ land in user space, never the read-only shared repo.
    rundir, core, run_flags, env = _run_state(passthrough)
    cmd = ["nextflow", *core, "run", pipeline, *run_flags, *prefix_args, *passthrough]
    sys.stderr.write("[callforge] + " + " ".join(shlex.quote(c) for c in cmd) + "\n")
    if print_cmd or os.environ.get("CALLFORGE_DRY_RUN"):
        # transparency / testability: echo without launching Nextflow
        print(" ".join(shlex.quote(c) for c in cmd))
        return 0
    _require_writable(rundir)
    try:
        return subprocess.call(cmd, env=env)
    except FileNotFoundError:
        sys.stderr.write("[callforge] ERROR: `nextflow` not found on PATH. Install "
                         "Nextflow (Java 17+) — see the manual §4.\n")
        return 127


def cmd_samplesheet(rest):
    # init_sample_sheet.main(argv) is importable; raw `python3 bin/init_sample_sheet.py`
    # still works too.
    mod = _load_bin("init_sample_sheet.py")
    return mod.main(rest) or 0


def cmd_resources(rest):
    # discover_resources.py is unmodified; invoke its main() via an argv shim.
    return _run_bin_via_argv("discover_resources.py", rest)


def cmd_metadata(rest):
    # generate gene_metadata.tsv WITHOUT CaptureForge (same schema the stages consume).
    # raw form `python3 bin/make_gene_metadata.py …` still works.
    mod = _load_bin("make_gene_metadata.py")
    return mod.main(rest) or 0


# ----------------------------------------------------------------- stage subcommands
# Each stage runs ONE pipeline stage standalone via `nextflow run main.nf --stage <name>`,
# reusing the SAME modules as the full pipeline. The CLI injects `--stage <name>` and
# passes everything else through to Nextflow. `required`/`optional` document the input
# files the entry workflow consumes (kept in sync with subworkflows/entries.nf).
#   Provenance tags (same model as `callforge run`):
#     [per-run]          a data file you provide for THIS standalone run (the VCF/BAMs/sheet)
#     [site]             a cluster reference set once (genome, target BED, annotation DBs)
#     [optional]         defaulted; set only to override (default shown)
#     [stage-dependent]  only needed for that engine/branch of the stage
STAGES = {
    "align": {
        "desc": "FASTQ -> analysis-ready BAM(s) (bwa-mem2 + dedup + BQSR) + alignment QC",
        "required": [
            "[per-run] --input <CSV>             sample sheet (sample_id,fastq_1,fastq_2[,sex])",
            "          OR --fastq_1 <R1> --fastq_2 <R2> [--sample_id <id>]   (single sample)",
            "[site]    --genome_fasta <FASTA>    no-alt primary assembly (indices built if absent)",
            "[site]    --target_bed <BED>        panel BED (the QC metrics feed the alignment plots)",
        ],
        "optional": ["[site] --known_sites <vcf,...>   BQSR known-sites (else BQSR is skipped)"],
        "outputs": "stage4_postalign/<sample>.analysis.bam + align_report.md "
                   "(mapping_rate, insert_size, dup_rate plots)",
    },
    "coverage": {
        "desc": "BAM(s) -> coverage/enrichment metrics + per-sample QC gate (pass/quarantine)",
        "required": [
            "[per-run] --input_bams <glob|csv>   analysis BAM(s) (with .bai)",
            "[site]    --genome_fasta <FASTA>    reference (HsMetrics)",
            "[site]    --target_bed <BED>        panel BED",
        ],
        "optional": [
            "[per-run] --input <CSV>             sample sheet to supply sex (else inferred 'U')",
            "[optional] thresholds: --min_mean_target_depth/--min_on_target/--max_dup_rate/…",
        ],
        "outputs": "stage5_qc_gate/qc_scorecard.tsv + qc_pass.txt/qc_quarantine.txt + "
                   "coverage_report.md (on_target/coverage_uniformity/per_gene_depth/cumulative/scorecard)",
    },
    "call": {
        "desc": "QC-PASS BAM(s) -> joint hard-filtered VCF (GATK, or --caller deepvariant)",
        "required": [
            "[per-run] --input_bams <glob|csv>   QC-PASS BAM(s) with .bai (call does NOT re-run the",
            "                                    QC gate — run `callforge coverage` first, pass the pass BAMs)",
            "[site]    --genome_fasta <FASTA>    reference",
            "[site]    --target_bed <BED>        panel BED",
        ],
        "optional": [
            "[optional]        --caller deepvariant   DeepVariant + GLnexus (default: gatk)",
            "[optional]        --joint_method <m>     genomicsdb | combinegvcfs (default: genomicsdb)",
        ],
        "outputs": "stage7_filter/joint.filtered.vcf.gz + call_report.md "
                   "(variants_per_sample, titv, het_hom, qual_dist, filter_*)",
    },
    "cnv": {
        "desc": "BAM(s) -> CNVkit calls, labelled by CaptureForge callability",
        "required": [
            "[per-run] --input_bams <glob|csv>   BAM(s) with .bai (PoN pooled from these)",
            "[site]    --genome_fasta <FASTA>    reference",
            "[site]    --target_bed <BED>        panel BED",
        ],
        "optional": [
            "[optional] CaptureForge callability via --cf_metrics_json / --captureforge_dir",
            "[optional] --cnv_method / --cnv_segment_method   (only --cnv_caller cnvkit is implemented)",
        ],
        "outputs": "stage8_cnv/cnv_calls.tsv (each call stamped callable | "
                   "breakpoint_blind_low_confidence) + cnv_report.md",
    },
    "str": {
        "desc": "BAM(s) -> ExpansionHunter genotypes on the CaptureForge STR catalog",
        "required": [
            "[per-run] --input_bams <glob|csv>   BAM(s) with .bai",
            "[site]    --genome_fasta <FASTA>    reference",
            "[site]    --target_bed <BED>        panel BED",
        ],
        "optional": [
            "[optional] --str_catalog <JSON>     CaptureForge STR catalog (else built from the",
            "                                    CaptureForge STR loci in gene_metadata)",
        ],
        "outputs": "stage9_str/str_calls.tsv + str_report.md "
                   "(str_allele_sizes, str_call_rate, str_read_support)",
    },
    "paralog": {
        "desc": "VCF -> paralog-aware flagging (INFO/PARALOG_GENE, PARALOG_CONF)",
        "required": [
            "[per-run] --input_vcf <VCF.gz>      callset VCF (with .tbi); uses its own MQ, not BAMs",
            "[site]    --target_bed <BED>        panel BED (CaptureForge paralog regions via metadata)",
        ],
        "optional": ["[optional] --paralog_min_mq <int>   MQ below which a paralog-region call is flagged (default 50)"],
        "outputs": "stage10_paralog/paralog.annotated.vcf.gz (INFO/PARALOG_GENE + "
                   "PARALOG_CONF) + paralog_report.md",
    },
    "anno": {
        "desc": "re-annotate a VCF (VEP + vcfanno DBs + PhyloP), e.g. after a DB update",
        "required": [
            "[per-run] --input_vcf <VCF.gz>      paralog-flagged or filtered VCF (with .tbi)",
            "[site]    --genome_fasta <FASTA>    no-alt primary assembly (reference)",
            "[site]    --target_bed <BED>        panel BED (reference invariant + DB scope gate)",
        ],
        "optional": [
            "[site]     --vcfanno_toml <TOML>    bring-your-own vcfanno databases",
            "[optional] --species <name>         VEP species (default homo_sapiens)",
            "[site]     annotation DBs are discovered from --resource_dirs / manifest",
        ],
        "outputs": "stage11_annotation/annotated.vcf.gz + variants.flat.tsv + anno_report.md",
    },
    "cohortqc": {
        "desc": "somalier relatedness/sex (off-target X/Y backstop) + ancestry PCA + missingness",
        "required": [
            "[per-run] --input_bams <glob|csv>   BAM(s) with .bai (somalier extract / sex)",
            "[per-run] --input_vcf <VCF.gz>      joint VCF (missingness + ancestry PCA)",
            "[site]    --genome_fasta <FASTA>    reference",
            "[site]    --somalier_sites <VCF.gz> somalier sites (+ .tbi)",
            "[per-run] --input <CSV>             sample sheet (declared sex)",
        ],
        "optional": ["[optional] (ancestry PCs feed `callforge burden --cohort_qc_json cohort_qc.json`)"],
        "outputs": "stage12_cohortqc/cohort_qc.json + cohortqc_report.md "
                   "(relatedness_heatmap, ancestry_pca, cohort_sex_check, missingness)",
    },
    "giab": {
        "desc": "hap.py precision/recall/F1 vs GIAB truth, restricted to the panel BED (on-target)",
        "required": [
            "[per-run] --input_vcf <VCF.gz>      joint callset VCF (with .tbi)",
            "[per-run] --giab_control_id <id>    control sample_id in the VCF to benchmark",
            "[site]    --giab_truth_vcf <VCF.gz> GIAB truth (+ .tbi); --giab_truth_bed <BED>",
            "[site]    --target_bed <BED>        panel BED (on-target restriction, hap.py -T)",
            "[site]    --genome_fasta <FASTA>    reference",
        ],
        "optional": ["[optional] (on-target restriction is recorded in happy.runinfo.json)"],
        "outputs": "stage13_giab/happy.summary.csv + happy.runinfo.json + giab_report.md "
                   "(giab_precision_recall, giab_f1)",
    },
    "burden": {
        "desc": "re-run rare-variant burden with new thresholds/engine on an annotated VCF",
        "required": [
            "[per-run] --input_vcf <VCF.gz>      annotated VCF (with .tbi)",
            "[per-run] --input <CSV>             sample sheet with phenotype + covariate_* columns",
            "[site]    --target_bed <BED>        panel BED (CaptureForge burden groups via metadata)",
        ],
        "optional": [
            "[optional]        --burden_engine <e>   collapse | regenie | skat (default: regenie)",
            "[optional]        --burden_af_max <f>   rare-AF threshold (default 0.01)",
            "[optional]        --burden_csq <list>   qualifying consequences",
            "[stage-dependent] --cohort_qc_json <JSON>  ancestry PCs (else none; PC-stability gate applies)",
        ],
        "outputs": "stage14_burden/burden_results.tsv + burden_report.md",
    },
}


def _stage_help(stage):
    s = STAGES[stage]
    lines = [f"callforge {stage} — {s['desc']}", "",
             f"usage: callforge {stage} [--print-cmd] [--in-place] [--pipeline PATH] <inputs> "
             f"[-profile …] [-params-file …] [-resume]", "",
             "Runs `nextflow run main.nf --stage %s …`, reusing the SAME module as the "
             "full pipeline (no logic is reimplemented). Output goes to a fresh, "
             "non-destructive results/standalone/%s_<timestamp>/ unless you pass --outdir "
             "or --in-place." % (stage, stage), "",
             "Provenance: [per-run]=a data file for THIS run · [site]=a cluster reference "
             "set once (often in -params-file) · [optional]=defaulted override · "
             "[stage-dependent]=only for that engine/branch.", "",
             "required input files:"]
    lines += [f"  {r}" for r in s["required"]]
    lines += ["", "optional:"]
    lines += [f"  {o}" for o in s["optional"]]
    lines += ["", f"key outputs: {s['outputs']}", "",
              "Typical (refs/queue/image from the params file):",
              f"  callforge {stage} <inputs> -profile hpc_slurm,singularity -params-file site.params.yaml",
              "", "Any other flag (-profile, -params-file, -resume, --<param>) is passed "
              "straight to Nextflow. See `callforge run --help` for the full parameter "
              "reference with provenance."]
    return "\n".join(lines)


def cmd_stage(stage, rest):
    if any(h in rest for h in ("-h", "--help")):
        print(_stage_help(stage))
        return 0
    ap = argparse.ArgumentParser(prog=f"callforge {stage}", allow_abbrev=False, add_help=False)
    ap.add_argument("--pipeline", metavar="PATH|owner/repo")
    ap.add_argument("--print-cmd", action="store_true")
    ap.add_argument("--in-place", action="store_true",
                    help="update the pipeline's results/ in place (default: a fresh "
                         "results/standalone/<stage>_<timestamp>/ — non-destructive)")
    known, passthrough = ap.parse_known_args(rest)

    # Non-destructive by default: send outputs to a fresh standalone dir unless the user
    # picked an --outdir or opted into --in-place. (The re-annotate / re-burden use cases
    # are comparisons — never clobber the baseline results/.)
    prefix = ["--stage", stage]
    if known.in_place:
        prefix += ["--in_place", "true"]
    elif "--outdir" not in passthrough:
        import datetime
        ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        outdir = os.path.join("results", "standalone", f"{stage}_{ts}")
        prefix += ["--outdir", outdir]
        sys.stderr.write(f"[callforge] standalone outputs -> {outdir} "
                         f"(use --outdir <dir> to change, or --in-place to update results/)\n")
    return _launch_nextflow(known.pipeline, prefix, passthrough, known.print_cmd)


# ----------------------------------------------------------------- dispatch
SUBCOMMANDS = {
    "run": (cmd_run, "run the full CallForge pipeline (wrapper over `nextflow run`)"),
    "samplesheet": (cmd_samplesheet, "scaffold/assemble a sample sheet (modes A/B/C)"),
    "metadata": (cmd_metadata, "generate gene_metadata.tsv without CaptureForge"),
    "resources": (cmd_resources, "resource discovery / genomic-scope checks"),
}
# stage subcommands (run one stage standalone, reusing the pipeline's modules)
for _st in STAGES:
    SUBCOMMANDS[_st] = ((lambda s: (lambda rest: cmd_stage(s, rest)))(_st),
                        "stage: " + STAGES[_st]["desc"])


def _top_help():
    lines = ["callforge — unified front door for the CallForge pipeline",
             f"version {__version__}", "",
             "usage: callforge <command> [options]", "", "commands:"]
    for name, (_, desc) in SUBCOMMANDS.items():
        lines.append(f"  {name:<13} {desc}")
    lines += ["", "  --version      print version",
              "  --help         show this help", "",
              "Raw forms still work: `python3 bin/init_sample_sheet.py …`, "
              "`nextflow run main.nf …`.",
              "Run `callforge <command> --help` for command-specific options."]
    return "\n".join(lines)


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in ("-h", "--help"):
        print(_top_help())
        return 0
    if argv[0] in ("-V", "--version"):
        print(__version__)
        return 0
    cmd, rest = argv[0], argv[1:]
    if cmd not in SUBCOMMANDS:
        sys.stderr.write(f"[callforge] unknown command '{cmd}'\n\n{_top_help()}\n")
        return 2
    return SUBCOMMANDS[cmd][0](rest)


if __name__ == "__main__":
    sys.exit(main())

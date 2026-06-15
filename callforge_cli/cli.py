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

Future enhancement (NOT built here)
-----------------------------------
Per-stage subcommands such as `callforge cnv` / `callforge anno` would map to
Nextflow `-entry <SUBWORKFLOW>` entry points to run/resume a single stage. The
extension point is marked `# >>> STAGE-SUBCOMMAND EXTENSION POINT <<<` in dispatch().
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


def _resolve_pipeline(override):
    """Resolve the Nextflow entry. Default: repo's main.nf. Override may be a path,
    a project dir, or a GitHub `owner/repo` name — passed to Nextflow verbatim."""
    if override:
        return override
    cand = os.path.join(REPO_ROOT, "main.nf")
    return cand if os.path.isfile(cand) else "main.nf"


# ----------------------------------------------------------------- subcommands
def cmd_run(rest):
    ap = argparse.ArgumentParser(
        prog="callforge run", allow_abbrev=False, add_help=True,
        description="Run the full CallForge pipeline (FASTQ -> annotated joint "
                    "callset -> burden). Thin wrapper over `nextflow run`; all other "
                    "flags (-profile, -params-file, -resume, --<param> overrides) are "
                    "passed straight through to Nextflow.",
        epilog="`-params-file params.full.yaml` supplies the reference (--genome_fasta), "
               "target BED + CaptureForge handoff (--target_bed/--captureforge_dir), "
               "annotation resources, and the sample sheet (--input). Per-stage conda "
               "envs/containers are activated automatically by Nextflow; you never "
               "activate them by hand.")
    ap.add_argument("--pipeline", metavar="PATH|owner/repo",
                    help="Nextflow entry to run (default: the repo's main.nf). May be a "
                         "path, a project directory, or a GitHub owner/repo.")
    ap.add_argument("--print-cmd", action="store_true",
                    help="print the `nextflow run …` command and exit (do not execute)")
    known, passthrough = ap.parse_known_args(rest)
    return _launch_nextflow(known.pipeline, [], passthrough, known.print_cmd)


def _launch_nextflow(pipeline_override, prefix_args, passthrough, print_cmd):
    """Build + run (or print) `nextflow run <pipeline> <prefix_args> <passthrough>`.
    prefix_args are CLI-injected flags (e.g. --stage anno); passthrough is everything
    the user typed (so -profile/-params-file/-resume/--<param> reach Nextflow)."""
    pipeline = _resolve_pipeline(pipeline_override)
    cmd = ["nextflow", "run", pipeline, *prefix_args, *passthrough]
    sys.stderr.write("[callforge] + " + " ".join(shlex.quote(c) for c in cmd) + "\n")
    if print_cmd or os.environ.get("CALLFORGE_DRY_RUN"):
        # transparency / testability: echo without launching Nextflow
        print(" ".join(shlex.quote(c) for c in cmd))
        return 0
    try:
        return subprocess.call(cmd)
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


# ----------------------------------------------------------------- stage subcommands
# Each stage runs ONE pipeline stage standalone via `nextflow run main.nf --stage <name>`,
# reusing the SAME modules as the full pipeline. The CLI injects `--stage <name>` and
# passes everything else through to Nextflow. `required`/`optional` document the input
# files the entry workflow consumes (kept in sync with subworkflows/entries.nf).
STAGES = {
    "align": {
        "desc": "FASTQ -> analysis-ready BAM(s) (bwa-mem2 + dedup + BQSR) + alignment QC",
        "required": [
            "--input <CSV>             sample sheet (sample_id,fastq_1,fastq_2[,sex])",
            "   OR --fastq_1 <R1> --fastq_2 <R2> [--sample_id <id>]   (single sample)",
            "--genome_fasta <FASTA>    no-alt primary assembly (indices built if absent)",
            "--target_bed <BED>        panel BED (the QC metrics feed the alignment plots)",
        ],
        "optional": ["--known_sites <vcf,...>   BQSR known-sites (else BQSR is skipped)"],
        "outputs": "stage4_postalign/<sample>.analysis.bam + align_report.md "
                   "(mapping_rate, insert_size, dup_rate plots)",
    },
    "coverage": {
        "desc": "BAM(s) -> coverage/enrichment metrics + per-sample QC gate (pass/quarantine)",
        "required": [
            "--input_bams <glob|csv>   analysis BAM(s) (with .bai)",
            "--genome_fasta <FASTA>    reference (HsMetrics)",
            "--target_bed <BED>        panel BED",
        ],
        "optional": [
            "--input <CSV>             sample sheet to supply sex (else inferred 'U')",
            "thresholds: --min_mean_target_depth/--min_on_target/--max_dup_rate/…",
        ],
        "outputs": "stage5_qc_gate/qc_scorecard.tsv + qc_pass.txt/qc_quarantine.txt + "
                   "coverage_report.md (on_target/coverage_uniformity/per_gene_depth/cumulative/scorecard)",
    },
    "call": {
        "desc": "QC-PASS BAM(s) -> joint hard-filtered VCF (GATK, or --caller deepvariant)",
        "required": [
            "--input_bams <glob|csv>   QC-PASS BAM(s) with .bai (call does NOT re-run the",
            "                          QC gate — run `callforge coverage` first and pass the pass BAMs)",
            "--genome_fasta <FASTA>    reference",
            "--target_bed <BED>        panel BED",
        ],
        "optional": [
            "--caller deepvariant      DeepVariant + GLnexus (default: gatk)",
            "--joint_method <m>        genomicsdb | combinegvcfs",
        ],
        "outputs": "stage7_filter/joint.filtered.vcf.gz + call_report.md "
                   "(variants_per_sample, titv, het_hom, qual_dist, filter_*)",
    },
    "anno": {
        "desc": "re-annotate a VCF (VEP + vcfanno DBs + PhyloP), e.g. after a DB update",
        "required": [
            "--input_vcf <VCF.gz>      paralog-flagged or filtered VCF (with .tbi)",
            "--genome_fasta <FASTA>    no-alt primary assembly (reference)",
            "--target_bed <BED>        panel BED (reference invariant + DB scope gate)",
        ],
        "optional": [
            "--vcfanno_toml <TOML>     bring-your-own vcfanno databases",
            "--species <name>          VEP species (default homo_sapiens)",
            "annotation DBs are discovered from --resource_dirs / manifest",
        ],
        "outputs": "stage11_annotation/annotated.vcf.gz + variants.flat.tsv + anno_report.md",
    },
    "burden": {
        "desc": "re-run rare-variant burden with new thresholds/engine on an annotated VCF",
        "required": [
            "--input_vcf <VCF.gz>      annotated VCF (with .tbi)",
            "--input <CSV>             sample sheet with phenotype + covariate_* columns",
            "--target_bed <BED>        panel BED (CaptureForge burden groups via metadata)",
        ],
        "optional": [
            "--burden_engine <e>       collapse | regenie | skat",
            "--burden_af_max <f>       rare-AF threshold (default 0.01)",
            "--burden_csq <list>       qualifying consequences",
            "--cohort_qc_json <JSON>   ancestry PCs (else none; PC-stability gate applies)",
        ],
        "outputs": "stage14_burden/burden_results.tsv + burden_report.md",
    },
}


def _stage_help(stage):
    s = STAGES[stage]
    lines = [f"callforge {stage} — {s['desc']}", "",
             f"usage: callforge {stage} [--print-cmd] [--pipeline PATH] <inputs> "
             f"[-profile …] [-params-file …] [-resume]", "",
             "Runs `nextflow run main.nf --stage %s …`, reusing the same module as the "
             "full pipeline." % stage, "",
             "required input files:"]
    lines += [f"  {r}" for r in s["required"]]
    lines += ["", "optional:"]
    lines += [f"  {o}" for o in s["optional"]]
    lines += ["", f"key outputs: {s['outputs']}",
              "", "Any other flag is passed straight to Nextflow."]
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

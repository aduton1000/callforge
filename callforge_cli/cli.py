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

    pipeline = _resolve_pipeline(known.pipeline)
    cmd = ["nextflow", "run", pipeline, *passthrough]
    sys.stderr.write("[callforge] + " + " ".join(shlex.quote(c) for c in cmd) + "\n")
    if known.print_cmd or os.environ.get("CALLFORGE_DRY_RUN"):
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


# ----------------------------------------------------------------- dispatch
SUBCOMMANDS = {
    "run": (cmd_run, "run the full CallForge pipeline (wrapper over `nextflow run`)"),
    "samplesheet": (cmd_samplesheet, "scaffold/assemble a sample sheet (modes A/B/C)"),
    "resources": (cmd_resources, "resource discovery / genomic-scope checks"),
    # >>> STAGE-SUBCOMMAND EXTENSION POINT <<<
    # Future: per-stage subcommands (e.g. "cnv", "anno") dispatching to
    #   callforge run --pipeline <main.nf> -entry <SUBWORKFLOW> …
    # to run/resume a single stage. Not built yet — see module docstring.
}


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

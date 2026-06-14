#!/usr/bin/env python3
"""regen_figures.py — regenerate the manual's QC figures from existing results data.

Manual-only helper: it re-runs the EXISTING pipeline plot scripts (bin/plot_*.py,
unmodified) against the existing results/ data, under a monkeypatch that suppresses
the baked-in bottom caption (the manual provides its own caption beneath each image)
so rotated x-axis labels are not overlapped. It does NOT modify pipeline code and
does NOT run the Nextflow pipeline.

Run from the repo root:  python3 docs/manual/regen_figures.py
Requires: matplotlib + numpy (for the plot scripts) and bcftools on PATH (cohort_qc).
"""
import glob, os, runpy, sys, traceback

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
BIN = os.path.join(REPO, "bin")
RES = os.path.join(REPO, "results")
FIG = os.path.join(REPO, "docs", "manual", "figures")
os.makedirs(FIG, exist_ok=True)
os.environ["CALLFORGE_CAPTIONS"] = "/tmp/_regen_caps.tsv"   # don't clobber figures/captions.tsv

# ---- monkeypatch: drop the baked bottom caption (y < 0.06) so labels get room ----
import matplotlib
matplotlib.use("Agg")
import matplotlib.figure as mfig
_orig_text = mfig.Figure.text
def _patched_text(self, x, y, s, *a, **k):
    if isinstance(y, (int, float)) and y < 0.06:     # the save() caption call
        return None
    return _orig_text(self, x, y, s, *a, **k)
mfig.Figure.text = _patched_text

# ---- monkeypatch: pad every title so it clears top-of-axes value labels ----
# (e.g. annotation_landing draws "100%" bar labels at ylim-top; with the default
#  tiny title pad the title and those labels collide. Extra pad separates them.)
import matplotlib.axes as maxes
_orig_set_title = maxes.Axes.set_title
def _patched_set_title(self, label, *a, **k):
    k.setdefault("pad", 18)
    return _orig_set_title(self, label, *a, **k)
maxes.Axes.set_title = _patched_set_title

# give rotated tick labels breathing room and pad the tight bbox
matplotlib.rcParams.update({"savefig.pad_inches": 0.25, "figure.autolayout": False})

qc = sorted(glob.glob(f"{RES}/stage5_qc_gate/per_sample/*.qc.json"))
ishist = sorted(glob.glob(f"{RES}/stage5_qc_gate/per_sample/*.ishist.tsv"))
regiondist = sorted(glob.glob(f"{RES}/stage4_postalign/mosdepth/*.mosdepth.region.dist.txt"))

JOBS = [
    ("plot_stage0.py", ["--manifest", f"{RES}/stage0_inputs/resource_manifest.json",
                        "--readcounts", f"{RES}/stage0_inputs/plots/input_readcounts.tsv", "--outdir", FIG]),
    ("plot_align_qc.py", ["--qc-json", *qc, "--ishist", *ishist, "--outdir", FIG]),
    ("plot_coverage.py", ["--qc-json", *qc, "--region-dist", *regiondist, "--outdir", FIG]),
    ("plot_qc_gate.py", ["--summary", f"{RES}/stage5_qc_gate/cohort_qc_summary.json",
                         "--qc-json", *qc, "--outdir", FIG]),
    ("plot_calling_qc.py", ["--stats-json", f"{RES}/stage6_calling/calling_stats.json",
                            "--qual-dp-tsv", f"{RES}/stage7_filter/qual_dp.tsv", "--outdir", FIG]),
    ("plot_filter_qc.py", ["--counts-json", f"{RES}/stage7_filter/filter_counts.json", "--outdir", FIG]),
    ("plot_cnv.py", ["--copyratio", f"{RES}/stage8_cnv/cnv_copyratio.tsv",
                     "--summary", f"{RES}/stage8_cnv/cnv_summary.json",
                     "--metadata", f"{RES}/stage0_inputs/gene_metadata.tsv", "--outdir", FIG]),
    ("plot_str.py", ["--calls", f"{RES}/stage9_str/str_calls.tsv",
                     "--summary", f"{RES}/stage9_str/str_summary.json", "--outdir", FIG]),
    ("plot_paralog.py", ["--summary", f"{RES}/stage10_paralog/paralog_summary.json", "--outdir", FIG]),
    ("plot_annotation.py", ["--tsv", f"{RES}/stage11_annotation/variants.flat.tsv",
                            "--landing", f"{RES}/stage11_annotation/annotation_landing.json",
                            "--af-field", "AF_afr", "--outdir", FIG]),
    ("cohort_qc.py", ["--vcf", f"{RES}/stage7_filter/joint.filtered.vcf.gz",
                      "--somalier-samples", f"{RES}/stage12_cohortqc/somalier/cohort.samples.tsv",
                      "--somalier-pairs", f"{RES}/stage12_cohortqc/somalier/cohort.pairs.tsv",
                      "--samplesheet", f"{RES}/stage0_inputs/samplesheet.valid.csv", "--outdir", FIG]),
    ("plot_giab.py", ["--summary", f"{RES}/stage13_giab/happy.summary.csv",
                      "--restricted-to", "panel BED (on-target)", "--outdir", FIG]),
    ("plot_burden.py", ["--results", f"{RES}/stage14_burden/burden_results.tsv",
                        "--summary", f"{RES}/stage14_burden/burden_summary.json", "--outdir", FIG]),
]

ok, fail = 0, []
for script, argv in JOBS:
    path = os.path.join(BIN, script)
    if not os.path.isfile(path):
        fail.append(f"{script} (missing)"); continue
    sys.argv = [path, *argv]
    try:
        runpy.run_path(path, run_name="__main__")
        ok += 1
    except SystemExit:
        ok += 1
    except Exception:
        fail.append(script); sys.stderr.write(f"[regen] {script} FAILED:\n{traceback.format_exc()}\n")

# cohort_qc.py also writes cohort_qc.json into FIG — harmless; remove stray non-figure files
for stray in ("cohort_qc.json",):
    p = os.path.join(FIG, stray)
    if os.path.exists(p):
        os.remove(p)
print(f"[regen] regenerated figures from {ok}/{len(JOBS)} plotters; failures={fail}")

#!/usr/bin/env python3
"""stage_report.py — focused per-stage report for a CallForge stage subcommand.

Stdlib-only. Writes a Markdown report for ONE stage that embeds that stage's QC
plots with their existing one-line captions (read from the `captions_*.tsv` the
pipeline's plot generators already emit) and renders that stage's key summary
metrics. It does NOT recompute anything — it presents the real stage outputs.

Inputs are staged into the cwd by Nextflow; the script globs `*.png` there, reads
the caption file(s), and renders any summary JSON/TSV passed via --summary.

  stage_report.py --stage anno --title "Annotation (Stage 11)" \
      --captions captions_annotation.tsv \
      --summary annotation_landing.json --summary variants.flat.tsv \
      --outdir . --out anno_report.md [--synthetic]
"""
import argparse, csv, glob, json, os, sys, datetime

def read_captions(paths):
    cap = {}
    for p in paths:
        if p and os.path.isfile(p):
            with open(p) as fh:
                for line in fh:
                    if "\t" in line:
                        stem, text = line.rstrip("\n").split("\t", 1)
                        cap[stem] = text
    return cap

def render_summary(path, max_rows=12):
    """Render a JSON (compact) or TSV (markdown table) summary file."""
    base = os.path.basename(path)
    if not os.path.isfile(path):
        return [f"- `{base}` (not found)"]
    if path.endswith(".json"):
        try:
            obj = json.load(open(path))
        except Exception as e:
            return [f"- `{base}`: unreadable JSON ({e})"]
        lines = [f"**`{base}`**", "", "```json", json.dumps(obj, indent=2)[:1800], "```"]
        return lines
    # TSV -> markdown table (first max_rows rows)
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    if not rows:
        return [f"- `{base}` (empty)"]
    out = [f"**`{base}`** (first {min(max_rows, len(rows)-1)} of {len(rows)-1} row(s))", ""]
    hdr = rows[0]
    out.append("| " + " | ".join(hdr) + " |")
    out.append("|" + "|".join("---" for _ in hdr) + "|")
    for r in rows[1:max_rows + 1]:
        out.append("| " + " | ".join(r) + " |")
    return out

def main(argv=None):
    ap = argparse.ArgumentParser(prog="stage_report.py")
    ap.add_argument("--stage", required=True)
    ap.add_argument("--title", required=True)
    ap.add_argument("--captions", nargs="*", default=[])
    ap.add_argument("--summary", action="append", default=[])
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--out", required=True)
    ap.add_argument("--synthetic", action="store_true",
                    help="note that inputs are synthetic/test (machinery, not data quality)")
    ap.add_argument("--input-note", default="")
    a = ap.parse_args(argv)

    caps = read_captions(a.captions)
    pngs = sorted(os.path.basename(p) for p in glob.glob(os.path.join(a.outdir, "*.png")))

    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    md = [f"# CallForge — {a.title}", "",
          f"*Standalone stage run · `callforge {a.stage}` · generated {ts}*", ""]
    if a.input_note:
        md += [f"**Inputs:** {a.input_note}", ""]
    if a.synthetic:
        md += ["> **Note:** inputs are synthetic/test fixtures — the figures demonstrate the "
               "machinery, not real data quality.", ""]

    md += ["## Key metrics", ""]
    if a.summary:
        for s in a.summary:
            md += render_summary(s) + [""]
    else:
        md += ["_(no summary files provided)_", ""]

    md += ["## QC plots", ""]
    if pngs:
        for png in pngs:
            stem = png[:-4]
            cap = caps.get(stem, stem.replace("_", " "))
            md += [f"### {stem}", "", f"![{cap}]({png})", "", f"*{cap}*", ""]
    else:
        md += ["_(no plots found)_", ""]

    out_path = os.path.join(a.outdir, a.out)
    with open(out_path, "w") as fh:
        fh.write("\n".join(md) + "\n")
    print(f"[stage_report] wrote {out_path} ({len(pngs)} plot(s), {len(a.summary)} summary file(s))")
    return 0

if __name__ == "__main__":
    sys.exit(main())

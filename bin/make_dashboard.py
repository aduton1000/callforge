#!/usr/bin/env python3
"""make_dashboard.py — self-contained cohort QC dashboard (CallForge, Stage 15).

Scans the working directory for captioned plots (captions_<stage>.tsv + matching
PNGs, staged from every stage) and key summary JSON/TSVs, and emits ONE
self-contained HTML (PNGs base64-embedded, every figure with its one-line caption,
grouped by stage) plus the per-sample QC scorecard and headline tables. stdlib only.
"""
import argparse, base64, csv, glob, html, json, os

STAGE_TITLES = {
    "stage0": "Stage 0 — Inputs & resource discovery",
    "align": "Stages 1-4 — Alignment & capture QC",
    "coverage": "Stages 1-4 — Coverage",
    "qc": "Stage 5 — Per-sample QC gate",
    "calling": "Stage 6 — SNV/indel calling",
    "filter": "Stage 7 — Hard-filtering",
    "cnv": "Stage 8 — CNV (CaptureForge callability)",
    "str": "Stage 9 — STR genotyping",
    "paralog": "Stage 10 — Paralog-aware flagging",
    "annotation": "Stage 11 — Annotation (VEP + vcfanno)",
    "cohortqc": "Stage 12 — Cohort QC (relatedness/ancestry/sex)",
    "giab": "Stage 13 — GIAB benchmarking (on-target)",
    "burden": "Stage 14 — Rare-variant burden",
}


def b64img(path):
    with open(path, "rb") as fh:
        return base64.b64encode(fh.read()).decode()


def load_captions():
    """stem -> (caption, stage_key) from every captions_*.tsv in cwd."""
    caps = {}
    for cf in sorted(glob.glob("captions_*.tsv")) + (["captions.tsv"] if os.path.exists("captions.tsv") else []):
        stage = cf.replace("captions_", "").replace(".tsv", "")
        for line in open(cf):
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                caps[parts[0]] = (parts[1], stage)
    return caps


def table_html(tsv, max_rows=30):
    if not os.path.isfile(tsv):
        return "<p><em>not available</em></p>"
    rows = list(csv.reader(open(tsv), delimiter="\t"))
    if not rows:
        return "<p><em>empty</em></p>"
    out = ["<table>"]
    out.append("<tr>" + "".join(f"<th>{html.escape(c)}</th>" for c in rows[0]) + "</tr>")
    for r in rows[1:max_rows + 1]:
        out.append("<tr>" + "".join(f"<td>{html.escape(c)}</td>" for c in r) + "</tr>")
    out.append("</table>")
    if len(rows) - 1 > max_rows:
        out.append(f"<p><em>… {len(rows)-1-max_rows} more rows</em></p>")
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", default="callforge")
    ap.add_argument("--out", default="cohort_qc_dashboard.html")
    a = ap.parse_args()

    caps = load_captions()
    # group plots by stage
    by_stage = {}
    for png in sorted(glob.glob("*.png")):
        stem = os.path.basename(png)[:-4]
        caption, stage = caps.get(stem, ("", "other"))
        by_stage.setdefault(stage, []).append((stem, png, caption))

    order = list(STAGE_TITLES.keys()) + sorted(set(by_stage) - set(STAGE_TITLES))
    parts = ["""<!DOCTYPE html><html><head><meta charset="utf-8"><title>CallForge cohort QC</title>
<style>body{font-family:-apple-system,Helvetica,Arial,sans-serif;margin:0;background:#f7f7f9;color:#222}
header{background:#1f2a44;color:#fff;padding:18px 28px}h1{margin:0;font-size:22px}
h2{border-bottom:2px solid #1f2a44;padding-bottom:4px;margin-top:34px;color:#1f2a44}
.wrap{max-width:1100px;margin:0 auto;padding:0 28px 60px}.fig{display:inline-block;vertical-align:top;width:340px;margin:10px;background:#fff;border:1px solid #ddd;border-radius:6px;padding:8px}
.fig img{width:100%;height:auto}.cap{font-size:11px;color:#555;font-style:italic;margin-top:4px}
table{border-collapse:collapse;background:#fff;font-size:12px;margin:8px 0}th,td{border:1px solid #ccc;padding:3px 7px;text-align:left}
th{background:#eef}nav a{color:#9cf;margin-right:12px;font-size:13px}</style></head><body>"""]
    parts.append(f'<header><h1>CallForge — cohort QC dashboard</h1><div>panel: {html.escape(a.panel)}</div>')
    parts.append('<nav>' + " ".join(f'<a href="#{s}">{STAGE_TITLES.get(s,s).split(" — ")[0]}</a>'
                                     for s in order if s in by_stage) + '</nav></header>')
    parts.append('<div class="wrap">')

    # headline tables (if present)
    parts.append('<h2 id="summary">Run summary</h2>')
    for label, f in [("Resource manifest", "resource_manifest.json"),
                     ("Per-sample QC scorecard", "qc_scorecard.tsv"),
                     ("Annotation landing", "annotation_landing.json"),
                     ("GIAB (on-target)", "giab_metrics.json"),
                     ("Burden results", "burden_results.tsv")]:
        if f.endswith(".tsv") and os.path.isfile(f):
            parts.append(f"<h3>{label}</h3>" + table_html(f))
        elif f.endswith(".json") and os.path.isfile(f):
            parts.append(f"<h3>{label}</h3><pre style='background:#fff;border:1px solid #ccc;padding:8px;font-size:11px;max-height:240px;overflow:auto'>"
                         + html.escape(json.dumps(json.load(open(f)), indent=2)[:4000]) + "</pre>")

    # plots by stage
    for s in order:
        if s not in by_stage:
            continue
        parts.append(f'<h2 id="{s}">{html.escape(STAGE_TITLES.get(s, s))}</h2>')
        for stem, png, caption in by_stage[s]:
            parts.append(f'<div class="fig"><img src="data:image/png;base64,{b64img(png)}">'
                         f'<div class="cap"><b>{html.escape(stem)}</b><br>{html.escape(caption)}</div></div>')
    parts.append("</div></body></html>")

    with open(a.out, "w") as fh:
        fh.write("\n".join(parts))
    n_fig = sum(len(v) for v in by_stage.values())
    print(f"[make_dashboard] {n_fig} figures across {len(by_stage)} stages -> {a.out}")


if __name__ == "__main__":
    main()

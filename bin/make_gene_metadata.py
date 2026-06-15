#!/usr/bin/env python3
"""make_gene_metadata.py — generate gene_metadata.tsv WITHOUT CaptureForge (CallForge).

The standalone on-ramp for the CNV / STR / paralog / burden interpretation layer, for
users who do not have a CaptureForge handoff. It emits the SAME schema
ingest_captureforge.py produces (so the stages consume it identically), built from a
gene list and/or a target BED, plus the user's own panel/study declarations.

It NEVER invents the design-specific fields a generic generator cannot know:
  * cnv_callable / cnv_callable_reason — capture-design breakpoint-blindness: set to
    'unknown' / 'user_declared_not_design_derived' for declared CNV targets (never a
    fabricated yes/no callability).
  * off_target_frac / low_coverage — design/run-derived: left blank (not assessed).
  * burden_group — irreducibly the user's (like phenotype): blank unless declared.

Stdlib-only. Mirrors the init_sample_sheet.py pattern (importable main(argv),
reconciliation report, fail-loud exit codes).

  make_gene_metadata.py --bed panel.bed \
      --cnv-genes CFH,CR1 --paralog-genes CD209,CFH \
      --str-loci str_loci.tsv --burden-groups burden_groups.tsv \
      --out-tsv gene_metadata.tsv --out-json gene_metadata.json
"""
import argparse, csv, difflib, json, os, re, sys
from collections import defaultdict

# identical column order to ingest_captureforge.py
COLS = ["gene", "classes", "n_intervals", "n_coding", "n_promoter", "n_anchor",
        "n_str", "n_cnv", "is_cnv_target", "cnv_callable", "cnv_callable_reason",
        "is_str_target", "str_loci", "is_paralog", "off_target_frac",
        "low_coverage", "burden_group"]
LOCUS_RE = re.compile(r"^[^:]+:\d+-\d+$")

def die(msg):
    sys.stderr.write(f"[make_gene_metadata] ERROR: {msg}\n")
    sys.exit(2)

def warn(msg):
    sys.stderr.write(f"[make_gene_metadata] WARN: {msg}\n")

def parse_bed(path):
    """gene -> {class: n_intervals} and gene -> {class: [(chrom,start,end)]} (GENE|class names)."""
    genes = defaultdict(lambda: defaultdict(int))
    loci = defaultdict(lambda: defaultdict(list))
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            chrom, start, end = f[0], int(f[1]), int(f[2])
            gene, _, klass = (f[3] if len(f) > 3 else "NA|coding").partition("|")
            klass = klass or "coding"
            genes[gene][klass] += 1
            loci[gene][klass].append((chrom, start, end))
    return genes, loci

def gtf_genes(path):
    out = set()
    op = __import__("gzip").open if path.endswith(".gz") else open
    with op(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or "\t" not in line:
                continue
            fields = line.split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            m = re.search(r'gene_name "([^"]+)"', fields[8]) or re.search(r'gene_id "([^"]+)"', fields[8])
            if m:
                out.add(m.group(1))
    return out

def read_list(spec):
    return [g.strip() for g in (spec or "").split(",") if g.strip()]

def read_2col(path, sep=None):
    """gene -> value (TSV/CSV: gene<sep>value). Returns dict + duplicate list."""
    out, dups = {}, []
    with open(path, newline="") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = re.split(r"[\t,]" if sep is None else sep, line.rstrip("\n"), maxsplit=1)
            if len(parts) < 2:
                continue
            g, v = parts[0].strip(), parts[1].strip()
            if g in out:
                dups.append(g)
            out[g] = v
    return out, dups


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="callforge metadata",
        description="Generate gene_metadata.tsv without CaptureForge (same schema the "
                    "stages consume). Never invents design-specific callability or "
                    "burden groups.")
    src = ap.add_argument_group("gene universe (need --bed and/or --genes/--gene-list)")
    src.add_argument("--bed", help="target BED (col4 'GENE|class'; classes/intervals/loci)")
    src.add_argument("--genes", help="comma list of genes")
    src.add_argument("--gene-list", help="file with one gene per line")
    src.add_argument("--gtf", help="gene annotation GTF (resolve/validate gene names)")
    dec = ap.add_argument_group("panel/study declarations (user-supplied; never inferred)")
    dec.add_argument("--cnv-genes", default="", help="comma list of CNV-target genes")
    dec.add_argument("--paralog-genes", default="", help="comma list of paralog-ambiguous genes")
    dec.add_argument("--str-loci", help="TSV/CSV: gene<sep>chr:start-end[;chr:start-end]")
    dec.add_argument("--burden-groups", help="TSV/CSV: gene<sep>burden_group (else blank)")
    ap.add_argument("--out-tsv", default="gene_metadata.tsv")
    ap.add_argument("--out-json", default=None)
    ap.add_argument("--report", default="", help="report path (default <out-tsv>_report.txt)")
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args(argv)

    if not (a.bed or a.genes or a.gene_list):
        die("need a gene universe: --bed and/or --genes / --gene-list")
    if os.path.exists(a.out_tsv) and not a.overwrite:
        die(f"output exists: {a.out_tsv} (use --overwrite)")

    # ---- gene universe ----
    bed_genes, bed_loci = (parse_bed(a.bed) if a.bed and os.path.isfile(a.bed) else ({}, {}))
    if a.bed and not bed_genes:
        die(f"no intervals parsed from --bed {a.bed}")
    listed, dup_listed = [], []
    if a.gene_list:
        if not os.path.isfile(a.gene_list):
            die(f"--gene-list not found: {a.gene_list}")
        seen = set()
        for line in open(a.gene_list):
            g = line.strip()
            if g and not g.startswith("#"):
                (dup_listed.append(g) if g in seen else None); seen.add(g); listed.append(g)
    listed += read_list(a.genes)
    gtf_set = gtf_genes(a.gtf) if a.gtf and os.path.isfile(a.gtf) else set()

    # ---- declarations ----
    cnv_genes = set(read_list(a.cnv_genes))
    paralog_genes = set(read_list(a.paralog_genes))
    str_map, dup_str = (read_2col(a.str_loci) if a.str_loci and os.path.isfile(a.str_loci) else ({}, []))
    burden_map, dup_burden = (read_2col(a.burden_groups) if a.burden_groups and os.path.isfile(a.burden_groups) else ({}, []))

    # validate str loci format (hard fail on malformed — would silently break the catalog)
    for g, loci_str in str_map.items():
        for locus in loci_str.split(";"):
            if locus and not LOCUS_RE.match(locus.strip()):
                die(f"--str-loci: malformed locus for {g}: '{locus}' (expect chr:start-end)")

    universe = set(bed_genes) | set(listed) | cnv_genes | paralog_genes | set(str_map) | set(burden_map)
    if not universe:
        die("no genes resolved from inputs")

    # ---- assemble rows (identical schema to ingest_captureforge.py) ----
    rows = []
    for g in sorted(universe):
        klasses = dict(bed_genes.get(g, {}))
        if g in cnv_genes:
            klasses.setdefault("cnv", klasses.get("cnv", 0))
        if g in str_map:
            klasses.setdefault("str", klasses.get("str", 0))
        if not klasses:
            klasses = {"coding": 0}          # default class when nothing else is known
        is_cnv = ("cnv" in klasses) or (g in cnv_genes)
        is_str = ("str" in klasses) or (g in str_map)
        # str_loci: user declaration wins; else BED str intervals
        sl = str_map.get(g) or ";".join(f"{c}:{s}-{e}" for (c, s, e) in bed_loci.get(g, {}).get("str", []))
        rows.append({
            "gene": g,
            "classes": ",".join(sorted(klasses)),
            "n_intervals": sum(klasses.values()),
            "n_coding": klasses.get("coding", 0),
            "n_promoter": klasses.get("promoter", 0),
            "n_anchor": klasses.get("anchor", 0),
            "n_str": klasses.get("str", 0),
            "n_cnv": klasses.get("cnv", 0),
            "is_cnv_target": "yes" if is_cnv else "no",
            # NEVER fabricate design callability:
            "cnv_callable": "unknown" if is_cnv else "n/a",
            "cnv_callable_reason": "user_declared_not_design_derived" if is_cnv else "",
            "is_str_target": "yes" if is_str else "no",
            "str_loci": sl,
            "is_paralog": "yes" if g in paralog_genes else "no",
            "off_target_frac": "",      # not design-derived
            "low_coverage": "",         # not assessed
            "burden_group": burden_map.get(g, ""),   # never auto-assign
        })

    if not a.overwrite and os.path.exists(a.out_tsv):
        die(f"output exists: {a.out_tsv} (use --overwrite)")

    # ---- reconciliation ----
    hard, dups, decl_not_in_bed, bed_not_listed = [], [], [], []
    for label, dl in (("--gene-list", dup_listed), ("--str-loci", dup_str), ("--burden-groups", dup_burden)):
        if dl:
            dups.append(f"{label}: {', '.join(sorted(set(dl)))}")
    if dups:
        hard.append("duplicate gene(s): " + "; ".join(dups))
    if bed_genes:
        bedset = set(bed_genes)
        for g in sorted((cnv_genes | paralog_genes | set(str_map))):
            if g not in bedset:
                near = difflib.get_close_matches(g, list(bedset), n=1, cutoff=0.6)
                decl_not_in_bed.append(f"{g}" + (f"  [closest BED gene: {near[0]}]" if near else "  [no near-miss]"))
        if listed:
            for g in sorted(bedset - set(listed)):
                bed_not_listed.append(g)

    # ---- write outputs (only if no hard error) ----
    if not hard:
        with open(a.out_tsv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t")
            w.writeheader(); w.writerows(rows)
        if a.out_json:
            summary = {
                "n_genes": len(rows), "source": ("bed" if a.bed else "gene-list"),
                "cnv_targets": [r["gene"] for r in rows if r["is_cnv_target"] == "yes"],
                "str_targets": [r["gene"] for r in rows if r["is_str_target"] == "yes"],
                "paralog_genes": [r["gene"] for r in rows if r["is_paralog"] == "yes"],
                "burden_grouped": [r["gene"] for r in rows if r["burden_group"]],
                "cnv_callable_note": "user_declared_not_design_derived (use CaptureForge for design callability)",
                "generator": "make_gene_metadata.py (standalone; not CaptureForge)",
            }
            json.dump(summary, open(a.out_json, "w"), indent=2)

    # ---- report ----
    report_path = a.report or (os.path.splitext(a.out_tsv)[0] + "_report.txt")
    n_cnv = sum(1 for r in rows if r["is_cnv_target"] == "yes")
    n_str = sum(1 for r in rows if r["is_str_target"] == "yes")
    n_par = sum(1 for r in rows if r["is_paralog"] == "yes")
    n_bg = sum(1 for r in rows if r["burden_group"])
    lines = ["CallForge gene_metadata generator — reconciliation report", "=" * 58,
             f"summary: {len(rows)} gene(s); {n_cnv} CNV target(s), {n_str} STR target(s), "
             f"{n_par} paralog gene(s), {n_bg} with a burden_group", ""]

    def block(title, items):
        lines.append(f"[{title}] ({len(items)})")
        lines.extend([f"  - {it}" for it in items] or ["  (none)"])
        lines.append("")

    block("DECLARED GENE NOT IN TARGET BED (typo? — with near-miss hint)", decl_not_in_bed)
    block("BED GENE WITH NO EXPLICIT DECLARATION (defaults applied)", bed_not_listed)
    block("CNV CALLABILITY NOT DESIGN-DERIVED (CaptureForge value-add)",
          [f"{r['gene']}: cnv_callable=unknown" for r in rows if r["is_cnv_target"] == "yes"])
    block("BURDEN GROUP NOT ASSIGNED (user must declare for gene-set burden)",
          [] if n_bg == len(rows) else
          [f"{len(rows) - n_bg} of {len(rows)} gene(s) have a blank burden_group "
           f"-> burden collapses per-gene only (supply --burden-groups for gene-set burden)"])
    if hard:
        lines.append("HARD ERRORS (no metadata written):")
        lines += [f"  ! {h}" for h in hard] + [""]
    report_txt = "\n".join(lines)
    open(report_path, "w").write(report_txt + "\n")
    sys.stderr.write(report_txt + "\n")

    for it in decl_not_in_bed:
        warn(f"declared gene not in target BED: {it}")
    if hard:
        for h in hard:
            sys.stderr.write(f"[make_gene_metadata] ERROR: {h}\n")
        sys.stderr.write(f"[make_gene_metadata] no metadata written; see {report_path}\n")
        sys.exit(1)

    print(f"[make_gene_metadata] wrote {len(rows)} genes -> {a.out_tsv} "
          f"(cnv={n_cnv}, str={n_str}, paralog={n_par}, burden_groups={n_bg}); report -> {report_path}")
    if n_bg < len(rows):
        print("[make_gene_metadata] NOTE: genes without a burden_group collapse per-gene only; "
              "supply --burden-groups for gene-set burden. CNV callability is user-declared "
              "(not design-derived) — CaptureForge provides breakpoint-blindness labels.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

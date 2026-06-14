#!/usr/bin/env python3
"""init_sample_sheet.py — scaffold / assemble a CallForge sample sheet (stdlib-only).

CallForge's sample sheet (see bin/validate_samplesheet.py) has columns:

    sample_id   (required, unique, [A-Za-z0-9._-])
    fastq_1     (required, paired-end R1)
    fastq_2     (required, paired-end R2; must differ from fastq_1)
    sex         (optional: M/F/male/female/1/2/unknown)
    phenotype   (optional: case/control 1/2/0 or quantitative; presence enables burden)
    covariate_* (optional: any number of numeric covariate columns)
    batch       (optional)

This helper ONLY pairs FASTQs and joins user-supplied metadata. It NEVER invents,
guesses, or defaults a phenotype (or any clinical field) — the phenotype is
irreducibly the user's. Anything it cannot resolve it REPORTS; it does not silently
fill. It reads the FASTQ/metadata files only (never moves or modifies them).

Three modes (all emit the same schema):

  A  scaffold from a FASTQ directory
       init_sample_sheet.py --fastq-dir DIR --out sample_sheet.csv
     -> sample_id/fastq_1/fastq_2 filled (absolute paths), metadata columns blank.

  B  scaffold FASTQs AND join a metadata/clinical CSV
       init_sample_sheet.py --fastq-dir DIR --metadata clinical.csv --out sample_sheet.csv
     -> complete sheet (sequencing + phenotype) in one step.

  C  join metadata onto an existing sample_id,fastq_1,fastq_2 pairing (e.g. a LIMS export)
       init_sample_sheet.py --fastq-csv manifest.csv --metadata clinical.csv --out sample_sheet.csv

Alongside the CSV it writes a reconciliation report (<out>_report.txt) listing orphan
FASTQs, sequence-without-metadata, metadata-without-sequence (with near-miss ID
candidates), duplicate IDs, rows missing a required field, and blank phenotypes. It
EXITS NON-ZERO on hard problems (duplicate sample_id, a written row missing a required
field, fastq_1 == fastq_2) so a broken sheet can never slip through; soft problems
(orphan FASTQs, unmatched IDs, blank phenotype) are warned but do not fail.
"""
import argparse, csv, difflib, os, re, sys

REQUIRED = ["sample_id", "fastq_1", "fastq_2"]
# canonical optional order mirrors validate_samplesheet.py's normalized output
META_ORDER = ["sex", "phenotype", "batch"]
ID_RE = re.compile(r"[A-Za-z0-9._-]+")
FASTQ_EXT = re.compile(r"\.(fastq|fq)(\.gz)?$", re.I)

# read-marker defaults: match the R1/R2 token (+ optional _001 lane chunk) sitting
# right before the FASTQ extension. Covers _R1_/_R1./.R1./_1./.1. conventions.
DEFAULT_R1 = r"(?i)(?:_R1|\.R1|_1|\.1)(?:_\d{3})?(?=\.(?:fastq|fq)(?:\.gz)?$)"
DEFAULT_R2 = r"(?i)(?:_R2|\.R2|_2|\.2)(?:_\d{3})?(?=\.(?:fastq|fq)(?:\.gz)?$)"


def warn(msg):
    sys.stderr.write(f"[init_sample_sheet] WARN: {msg}\n")


def die(msg):
    sys.stderr.write(f"[init_sample_sheet] ERROR: {msg}\n")
    sys.exit(2)


# ---------------------------------------------------------------- FASTQ pairing
def split_read(name, r1_re, r2_re):
    """Return (stem, read) where read is 1 or 2, or (None, None) if no marker."""
    for read, rx in ((1, r1_re), (2, r2_re)):
        m = rx.search(name)
        if m:
            stem = name[:m.start()] + name[m.end():]
            return FASTQ_EXT.sub("", stem).strip("._-"), read
    return None, None


def infer_sample_id(stem, id_regex):
    if id_regex:
        m = re.search(id_regex, stem)
        if m:
            return (m.group(1) if m.groups() else m.group(0)).strip("._-")
    s = re.sub(r"(?i)_S\d+(?=$|_)", "", stem)   # bcl2fastq sample-number token
    s = re.sub(r"(?i)_L\d{3}(?=$|_)", "", s)     # lane token
    return s.strip("._-")


def scan_fastq_dir(fastq_dir, recursive, r1_re, r2_re, id_regex):
    """Pair FASTQs in a directory. Returns (pairs, orphans).
    pairs: list of {sample_id, fastq_1, fastq_2, _stem}; orphans: list of paths."""
    files = []
    if recursive:
        for root, _, names in os.walk(fastq_dir):
            files += [os.path.join(root, n) for n in names if FASTQ_EXT.search(n)]
    else:
        files = [os.path.join(fastq_dir, n) for n in sorted(os.listdir(fastq_dir))
                 if FASTQ_EXT.search(n)]
    by_stem = {}   # stem -> {1: path, 2: path}
    orphans = []
    for path in sorted(files):
        stem, read = split_read(os.path.basename(path), r1_re, r2_re)
        if read is None:
            orphans.append(os.path.abspath(path))   # no R1/R2 marker at all
            continue
        by_stem.setdefault(stem, {})[read] = os.path.abspath(path)
    pairs = []
    for stem, reads in sorted(by_stem.items()):
        if 1 in reads and 2 in reads:
            pairs.append({"sample_id": infer_sample_id(stem, id_regex),
                          "fastq_1": reads[1], "fastq_2": reads[2], "_stem": stem})
        else:
            orphans += list(reads.values())          # only one mate present
    return pairs, orphans


def read_fastq_csv(path):
    """Mode C: read an existing sample_id,fastq_1,fastq_2 manifest."""
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        die(f"--fastq-csv is empty: {path}")
    for c in REQUIRED:
        if c not in rows[0]:
            die(f"--fastq-csv missing required column '{c}' (found: {list(rows[0])})")
    pairs = []
    for r in rows:
        pairs.append({"sample_id": (r.get("sample_id") or "").strip(),
                      "fastq_1": (r.get("fastq_1") or "").strip(),
                      "fastq_2": (r.get("fastq_2") or "").strip(), "_stem": ""})
    return pairs


# ------------------------------------------------------------------- metadata
def read_metadata(path, id_col):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        die(f"--metadata is empty: {path}")
    if id_col not in rows[0]:
        die(f"--metadata missing id column '{id_col}' (found: {list(rows[0])})")
    return rows, list(rows[0].keys())


def norm_key(s, strip_suffixes, lowercase):
    s = (s or "").strip()
    for suf in strip_suffixes:
        if suf and s.endswith(suf):
            s = s[:-len(suf)]
    return s.lower() if lowercase else s


# ---------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(
        description="Scaffold/assemble a CallForge sample sheet (pairs FASTQs, joins "
                    "user metadata; never invents phenotype).")
    src = ap.add_argument_group("inputs (choose a mode)")
    src.add_argument("--fastq-dir", help="Mode A/B: directory of paired FASTQs")
    src.add_argument("--fastq-csv", help="Mode C: existing sample_id,fastq_1,fastq_2 CSV")
    src.add_argument("--metadata", help="Mode B/C: clinical/study metadata CSV to join")
    ap.add_argument("--recursive", action="store_true", help="recurse into --fastq-dir")
    ap.add_argument("--r1-pattern", default=DEFAULT_R1, help="regex for the R1 read marker")
    ap.add_argument("--r2-pattern", default=DEFAULT_R2, help="regex for the R2 read marker")
    ap.add_argument("--id-regex", default="",
                    help="regex (capture group 1 = sample_id) applied to the paired stem")
    ap.add_argument("--metadata-id-col", default="sample_id",
                    help="metadata column to join on (default sample_id)")
    ap.add_argument("--keep-extra-cols", action="store_true",
                    help="pass through metadata columns beyond the known schema")
    ap.add_argument("--strip-suffixes", default="",
                    help="comma list of suffixes to strip from IDs before matching")
    ap.add_argument("--lowercase-ids", action="store_true",
                    help="case-insensitive ID matching")
    ap.add_argument("--out", default="sample_sheet.csv")
    ap.add_argument("--report", default="", help="report path (default: <out>_report.txt)")
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args()

    if not a.fastq_dir and not a.fastq_csv:
        die("need --fastq-dir (mode A/B) or --fastq-csv (mode C)")
    if a.fastq_dir and a.fastq_csv:
        die("use either --fastq-dir or --fastq-csv, not both")
    if a.fastq_csv and not a.metadata:
        die("--fastq-csv requires --metadata (mode C joins metadata onto a pairing)")
    if os.path.exists(a.out) and not a.overwrite:
        die(f"output exists: {a.out} (use --overwrite)")

    mode = "C" if a.fastq_csv else ("B" if a.metadata else "A")
    strip_suffixes = [s for s in a.strip_suffixes.split(",") if s] if a.strip_suffixes else []

    # ---- gather FASTQ pairings ----
    if a.fastq_csv:
        if not os.path.isfile(a.fastq_csv):
            die(f"--fastq-csv not found: {a.fastq_csv}")
        pairs, orphans = read_fastq_csv(a.fastq_csv), []
    else:
        if not os.path.isdir(a.fastq_dir):
            die(f"--fastq-dir not found: {a.fastq_dir}")
        try:
            r1_re, r2_re = re.compile(a.r1_pattern), re.compile(a.r2_pattern)
        except re.error as e:
            die(f"bad --r1/--r2 pattern: {e}")
        pairs, orphans = scan_fastq_dir(a.fastq_dir, a.recursive, r1_re, r2_re, a.id_regex)

    # ---- detect duplicate sample_ids among the pairings ----
    seen, dup_fastq = set(), []
    for p in pairs:
        if p["sample_id"] in seen:
            dup_fastq.append(p["sample_id"])
        seen.add(p["sample_id"])

    # ---- load + join metadata ----
    meta_rows, meta_cols, meta_extra_cols = [], [], []
    if a.metadata:
        if not os.path.isfile(a.metadata):
            die(f"--metadata not found: {a.metadata}")
        meta_rows, meta_cols = read_metadata(a.metadata, a.metadata_id_col)
        known = set(REQUIRED) | set(META_ORDER) | {a.metadata_id_col}
        cov_cols = sorted(c for c in meta_cols if c.startswith("covariate_"))
        meta_extra_cols = [c for c in meta_cols
                           if c not in known and not c.startswith("covariate_")]
    else:
        cov_cols = []

    # index metadata by normalized key; flag duplicate metadata ids
    meta_index, dup_meta, meta_norm_keys = {}, [], {}
    for r in meta_rows:
        raw = (r.get(a.metadata_id_col) or "").strip()
        k = norm_key(raw, strip_suffixes, a.lowercase_ids)
        if k in meta_index:
            dup_meta.append(raw)
        meta_index[k] = r
        meta_norm_keys[k] = raw

    # ---- assemble rows ----
    out_meta_cols = META_ORDER + cov_cols + (meta_extra_cols if a.keep_extra_cols else [])
    out_cols = REQUIRED + out_meta_cols
    matched_meta_keys, seq_no_meta, rows_out = set(), [], []
    for p in pairs:
        rec = {"sample_id": p["sample_id"], "fastq_1": p["fastq_1"], "fastq_2": p["fastq_2"]}
        for c in out_meta_cols:
            rec[c] = ""
        if a.metadata:
            k = norm_key(p["sample_id"], strip_suffixes, a.lowercase_ids)
            m = meta_index.get(k)
            if m:
                matched_meta_keys.add(k)
                for c in out_meta_cols:
                    rec[c] = (m.get(c) or "").strip()
            else:
                seq_no_meta.append(p["sample_id"])
        rows_out.append(rec)

    # metadata rows that never matched a FASTQ pairing (with near-miss candidates)
    fastq_keys = [norm_key(p["sample_id"], strip_suffixes, a.lowercase_ids) for p in pairs]
    meta_no_seq = []
    for k, raw in meta_norm_keys.items():
        if k not in matched_meta_keys and a.metadata:
            near = difflib.get_close_matches(k, fastq_keys, n=1, cutoff=0.6)
            cand = next((p["sample_id"] for p in pairs
                         if norm_key(p["sample_id"], strip_suffixes, a.lowercase_ids) == near[0]),
                        None) if near else None
            meta_no_seq.append((raw, cand))

    # ---- validation of the assembled rows ----
    hard, missing_required, bad_id, blank_pheno = [], [], [], []
    for rec in rows_out:
        for c in REQUIRED:
            if not rec.get(c):
                missing_required.append(rec.get("sample_id") or "<blank>")
                break
        sid = rec.get("sample_id", "")
        if sid and not ID_RE.fullmatch(sid):
            bad_id.append(sid)
        if rec["fastq_1"] and rec["fastq_1"] == rec["fastq_2"]:
            hard.append(f"fastq_1 == fastq_2 for '{sid}'")
        if "phenotype" in out_cols and not rec.get("phenotype"):
            blank_pheno.append(sid)

    dup_all = sorted(set(dup_fastq) | set(dup_meta))
    if dup_all:
        hard.append(f"duplicate sample_id(s): {', '.join(dup_all)}")
    if missing_required:
        hard.append(f"row(s) missing a required field: {', '.join(missing_required)}")
    if bad_id:
        hard.append(f"sample_id(s) with illegal characters: {', '.join(bad_id)}")

    # ---- write CSV (only if no hard error would emit a broken sheet) ----
    if not hard:
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=out_cols)
            w.writeheader()
            for rec in rows_out:
                w.writerow({c: rec.get(c, "") for c in out_cols})

    # ---- reconciliation report ----
    report_path = a.report or (os.path.splitext(a.out)[0] + "_report.txt")
    n_total = len(rows_out)
    n_joined = sum(1 for r in rows_out if r.get("phenotype")) if a.metadata else 0
    n_seqonly = len(seq_no_meta)
    n_metaonly = len(meta_no_seq)
    lines = []
    lines.append(f"CallForge sample-sheet reconciliation report (mode {mode})")
    lines.append("=" * 60)
    lines.append(f"summary: {n_total} sample(s) written, {n_joined} fully-joined, "
                 f"{n_seqonly} sequence-only, {n_metaonly} metadata-only, "
                 f"{len(orphans)} orphan FASTQ(s)")
    lines.append("")

    def block(title, items):
        lines.append(f"[{title}] ({len(items)})")
        for it in items:
            lines.append(f"  - {it}")
        if not items:
            lines.append("  (none)")
        lines.append("")

    block("ORPHAN FASTQ (found but not paired — no mate)", orphans)
    block("SEQUENCE WITHOUT METADATA (FASTQ paired, no clinical row)", seq_no_meta)
    block("METADATA WITHOUT SEQUENCE (clinical row, no FASTQ)",
          [f"{raw}" + (f"  [closest FASTQ id: {cand}]" if cand else "  [no near-miss]")
           for raw, cand in meta_no_seq])
    block("DUPLICATE sample_id", dup_all)
    block("ROWS MISSING A REQUIRED FIELD", missing_required)
    block("BLANK PHENOTYPE (burden needs it — fill manually, never auto-filled)",
          blank_pheno if a.metadata else [])
    if hard:
        lines.append("HARD ERRORS (no sheet written):")
        for h in hard:
            lines.append(f"  ! {h}")
        lines.append("")

    report_txt = "\n".join(lines)
    with open(report_path, "w") as fh:
        fh.write(report_txt + "\n")
    sys.stderr.write(report_txt + "\n")

    # ---- warnings + exit ----
    for o in orphans:
        warn(f"orphan FASTQ (no mate): {o}")
    for s in seq_no_meta:
        warn(f"sequence-only (no metadata) for sample_id: {s}")
    for raw, cand in meta_no_seq:
        warn(f"metadata-only (no FASTQ) for '{raw}'" + (f" (closest: {cand})" if cand else ""))
    for s in blank_pheno:
        warn(f"blank phenotype for '{s}' — burden stage needs it (not auto-filled)")

    if hard:
        for h in hard:
            sys.stderr.write(f"[init_sample_sheet] ERROR: {h}\n")
        sys.stderr.write(f"[init_sample_sheet] no sheet written; see {report_path}\n")
        sys.exit(1)

    print(f"[init_sample_sheet] OK (mode {mode}): wrote {n_total} sample(s) -> {a.out}")
    print(f"[init_sample_sheet] report -> {report_path}")
    if a.metadata and blank_pheno:
        print(f"[init_sample_sheet] NOTE: {len(blank_pheno)} row(s) have a blank phenotype; "
              "fill them before the burden stage (the pipeline otherwise stops at the "
              "annotated callset).")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""discover_resources.py — Stage 0 reference/resource auto-discovery (CallForge, §5).

Scans the given directories for annotation / known-sites resources, VALIDATES
each (build = the design build, indexed/tabixed, non-empty), and emits a
resource MANIFEST (path, type, build, fingerprint, status, notes).

Resources located (organism-generic):
  vep_cache    VEP cache for <species>+<assembly>   (dir: <species>/<release>_<assembly>)
  gnomad_vcf   gnomAD sites VCF (prefers one carrying AFR frequencies)
  dbsnp_vcf    dbSNP VCF
  clinvar_vcf  ClinVar VCF
  phylop_bw    PhyloP conservation bigWig
  known_sites  BQSR known-sites (dbSNP + Mills + 1000G indels)

Behaviour: USE WHAT'S FOUND; DOWNLOAD NOTHING. Human-only resources gracefully
report "missing" when absent (e.g. a non-human panel). Explicit --override paths
are validated and take precedence over a scan hit.

stdlib-only. Fingerprint is fast (size + md5 of the first 1 MiB) unless --full-md5.
"""
import argparse, glob, gzip, hashlib, json, os, re, shutil, subprocess, sys
from collections import defaultdict

SKIP_DIRS = {".git", "work", ".nextflow", "node_modules", "__pycache__",
             ".cache", "Library", ".Trash", "results", "results_full"}

HAVE_TABIX = shutil.which("tabix") is not None


def bed_spans(path):
    """Per-contig (min_start, max_end) across the target BED."""
    spans = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            c, s, e = f[0], int(f[1]), int(f[2])
            if c in spans:
                spans[c] = (min(spans[c][0], s), max(spans[c][1], e))
            else:
                spans[c] = (s, e)
    return spans


def norm_contig(c):
    return c[3:] if c.lower().startswith("chr") else c


def vcf_contig_map(path):
    """norm_name -> actual VCF contig name, from the tabix index (fast)."""
    try:
        out = subprocess.run(["tabix", "-l", path], check=True, capture_output=True,
                             text=True, timeout=120).stdout
    except (subprocess.SubprocessError, OSError):
        return None
    return {norm_contig(c.strip()): c.strip() for c in out.splitlines() if c.strip()}


def vcf_scope(path, spans, scope_min):
    """Does the VCF's data span the target BED? Returns a scope dict.

    A region-subset DB (e.g. a single-locus gnomAD) has records on only a few of
    the target's contigs -> low fraction -> flagged unfit. chr/no-chr naming is
    reconciled. Probes one BED span per contig (tabix stops at the first record).
    """
    if not HAVE_TABIX:
        return {"scope_checked": False, "note": "tabix unavailable; scope not gated"}
    cmap = vcf_contig_map(path)
    if cmap is None:
        return {"scope_checked": False, "note": "no usable tabix index for scope check"}
    total = len(spans)
    covered, missing = 0, []
    for c, (s, e) in spans.items():
        actual = cmap.get(norm_contig(c))
        hit = False
        if actual:
            # Stream and stop at the first record — a wide span over a genome-wide
            # DB (e.g. dbSNP) can emit huge output; we only need existence.
            proc = None
            try:
                proc = subprocess.Popen(["tabix", path, f"{actual}:{s+1}-{e}"],
                                        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
                first = proc.stdout.readline()
                hit = bool(first.strip())
            except (subprocess.SubprocessError, OSError):
                pass
            finally:
                if proc:
                    proc.stdout.close()
                    proc.terminate()
                    try:
                        proc.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        proc.kill()
        if hit:
            covered += 1
        else:
            missing.append(c)
    frac = covered / total if total else 0.0
    return {"scope_checked": True, "contigs_total": total, "contigs_covered": covered,
            "scope_fraction": round(frac, 4), "uncovered_contigs": missing[:20],
            "fit": frac >= scope_min}


def fingerprint(path, full=False):
    h = hashlib.md5()
    size = os.path.getsize(path)
    with open(path, "rb") as fh:
        if full:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        else:
            h.update(fh.read(1 << 20))
    tag = "md5" if full else "md5:first1MiB"
    return {"size_bytes": size, tag: h.hexdigest()}


def vcf_header(path, max_lines=4000):
    """Read header lines from a (b)gzipped or plain VCF."""
    opener = gzip.open if path.endswith(".gz") else open
    out = []
    try:
        with opener(path, "rt", errors="replace") as fh:
            for _ in range(max_lines):
                line = fh.readline()
                if not line or not line.startswith("#"):
                    break
                out.append(line)
    except OSError:
        pass
    return out


def detect_build(header_lines, path):
    blob = "\n".join(header_lines) + "\n" + path
    if re.search(r"GRCh38|hg38|\b38\b", blob):
        return "GRCh38"
    if re.search(r"GRCh37|hg19|\b37\b", blob):
        return "GRCh37"
    return "unknown"


def has_afr(header_lines):
    blob = "\n".join(header_lines)
    return bool(re.search(r"\bAF_afr\b|\bAF_AFR\b|gnomAD_AF_afr|\bafr\b", blob, re.I))


def indexed(path):
    for ext in (".tbi", ".csi"):
        if os.path.isfile(path + ext):
            return True
    return False


def validate_vcf(path, want_build, full_md5, need_afr=False, spans=None, scope_min=0.9):
    rec = {"path": path, "exists": os.path.isfile(path)}
    if not rec["exists"]:
        rec["status"] = "missing"
        return rec
    rec.update(fingerprint(path, full_md5))
    hdr = vcf_header(path)
    rec["build"] = detect_build(hdr, path)
    rec["indexed"] = indexed(path)
    rec["non_empty"] = rec["size_bytes"] > 0
    notes = []
    status = "ok"
    if not rec["indexed"]:
        status = "unindexed"; notes.append("no .tbi/.csi — run: tabix -p vcf")
    if rec["build"] not in (want_build, "unknown"):
        status = "build_mismatch"; notes.append(f"header build {rec['build']} != design {want_build}")
    if need_afr:
        rec["has_afr"] = has_afr(hdr)
        if not rec["has_afr"]:
            notes.append("no AFR-population AF field detected in header")
    # Genomic-scope gate: a found-but-region-limited DB silently under-annotates
    # outside its region (a missing gnomAD-AFR AF can be misread as rare-novel),
    # so a DB whose data does not span the target BED is flagged UNFIT, not ok.
    if spans and rec["indexed"]:
        sc = vcf_scope(path, spans, scope_min)
        rec["scope"] = sc
        if sc.get("scope_checked") and not sc.get("fit"):
            status = "unfit_scope"
            notes.append(f"data spans only {sc['contigs_covered']}/{sc['contigs_total']} "
                         f"target contigs (frac={sc['scope_fraction']} < {scope_min}); "
                         f"region-subset DB — replace with a genome-wide resource")
    rec["status"] = status
    rec["notes"] = "; ".join(notes)
    return rec


def find_files(dirs, patterns, maxdepth=6):
    hits = []
    for root_dir in dirs:
        root_dir = os.path.expanduser(root_dir.strip())
        if not root_dir or not os.path.isdir(root_dir):
            continue
        base_depth = root_dir.rstrip("/").count("/")
        for dirpath, dirnames, filenames in os.walk(root_dir):
            if dirpath.count("/") - base_depth > maxdepth:
                dirnames[:] = []
                continue
            dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS and not d.startswith(".Trash")]
            for fn in filenames:
                low = fn.lower()
                for pat in patterns:
                    if pat(low):
                        hits.append(os.path.join(dirpath, fn))
    # de-dup, stable
    seen, out = set(), []
    for h in hits:
        if h not in seen:
            seen.add(h); out.append(h)
    return out


def find_vep_cache(dirs, species, build, maxdepth=5):
    """A VEP cache is a dir <...>/<species>/<release>_<assembly>[/info.txt]."""
    found = []
    rx = re.compile(rf"{re.escape(species)}/(\d+)_({re.escape(build)})\b")
    for root_dir in dirs:
        root_dir = os.path.expanduser(root_dir.strip())
        if not os.path.isdir(root_dir):
            continue
        base_depth = root_dir.rstrip("/").count("/")
        for dirpath, dirnames, filenames in os.walk(root_dir):
            if dirpath.count("/") - base_depth > maxdepth:
                dirnames[:] = []
                continue
            dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
            m = rx.search(dirpath.replace(os.sep, "/"))
            if m:
                found.append({"path": dirpath, "release": int(m.group(1)), "assembly": m.group(2)})
                dirnames[:] = []  # don't descend further into a cache
    return found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dirs", required=True, help="comma-separated directories to scan")
    ap.add_argument("--species", default="homo_sapiens")
    ap.add_argument("--genome-build", default="GRCh38")
    ap.add_argument("--vep-release", type=int, default=None, help="preferred VEP cache release")
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--full-md5", action="store_true")
    # explicit overrides (validated, take precedence)
    for k in ("vep-cache", "gnomad-vcf", "dbsnp-vcf", "clinvar-vcf", "phylop-bw"):
        ap.add_argument(f"--{k}", default=None)
    ap.add_argument("--known-sites", default=None, help="comma list of BQSR known-sites VCFs")
    ap.add_argument("--target-bed", default=None,
                    help="target BED; enables the genomic-scope gate on annotation VCFs")
    ap.add_argument("--scope-min", type=float, default=0.9,
                    help="min fraction of target contigs a DB must span to be 'fit'")
    ap.add_argument("--fail-on-unfit", action="store_true",
                    help="exit 2 if any FOUND annotation DB is unfit_scope (real-run gate)")
    a = ap.parse_args()

    dirs = [d for d in a.dirs.split(",") if d.strip()]
    build = a.genome_build
    human = a.species == "homo_sapiens"
    spans = bed_spans(a.target_bed) if (a.target_bed and os.path.isfile(a.target_bed)) else None
    manifest = {"species": a.species, "genome_build": build, "scanned_dirs": dirs,
                "scope_gate": {"enabled": bool(spans), "scope_min": a.scope_min,
                               "fail_on_unfit": a.fail_on_unfit,
                               "tabix_available": HAVE_TABIX},
                "resources": {}}
    R = manifest["resources"]

    # ---- VEP cache ----
    if a.vep_cache:
        R["vep_cache"] = {"path": a.vep_cache, "status": "ok" if os.path.isdir(a.vep_cache) else "missing",
                          "source": "override"}
    else:
        caches = find_vep_cache(dirs, a.species, build)
        if caches:
            if a.vep_release:
                caches.sort(key=lambda c: (c["release"] != a.vep_release, -c["release"]))
            else:
                caches.sort(key=lambda c: -c["release"])
            best = caches[0]
            R["vep_cache"] = {"path": best["path"], "release": best["release"],
                              "assembly": best["assembly"], "status": "ok", "source": "discovered",
                              "all_found": caches}
        else:
            R["vep_cache"] = {"path": None, "status": "missing",
                              "notes": f"no {a.species}/<release>_{build} cache under scanned dirs"}

    # ---- VCF resources ----
    def pat_any(*subs):
        return lambda low: low.endswith((".vcf.gz", ".vcf.bgz")) and any(s in low for s in subs)

    sm = a.scope_min
    fit_first = lambda r: r.get("status") not in ("ok",)  # ok sorts before unfit/other

    # gnomAD (prefer AFR + fit scope)
    if a.gnomad_vcf:
        R["gnomad_vcf"] = validate_vcf(a.gnomad_vcf, build, a.full_md5, need_afr=True, spans=spans, scope_min=sm)
        R["gnomad_vcf"]["source"] = "override"
    elif human:
        cands = find_files(dirs, [pat_any("gnomad")])
        recs = [validate_vcf(c, build, a.full_md5, need_afr=True, spans=spans, scope_min=sm) for c in cands]
        recs.sort(key=lambda r: (r.get("status") != "ok", not r.get("has_afr", False)))
        R["gnomad_vcf"] = recs[0] if recs else {"path": None, "status": "missing"}
        if len(recs) > 1:
            R["gnomad_vcf"]["other_candidates"] = [{"path": r["path"], "status": r["status"]} for r in recs[1:]]
    else:
        R["gnomad_vcf"] = {"path": None, "status": "n/a", "notes": "non-human; gnomAD not applicable"}

    # dbSNP
    if a.dbsnp_vcf:
        R["dbsnp_vcf"] = validate_vcf(a.dbsnp_vcf, build, a.full_md5, spans=spans, scope_min=sm)
        R["dbsnp_vcf"]["source"] = "override"
    elif human:
        cands = find_files(dirs, [pat_any("dbsnp", "00-all", "common_all")])
        recs = [validate_vcf(c, build, a.full_md5, spans=spans, scope_min=sm) for c in cands]
        recs.sort(key=fit_first)
        R["dbsnp_vcf"] = recs[0] if recs else {"path": None, "status": "missing"}
    else:
        R["dbsnp_vcf"] = {"path": None, "status": "n/a"}

    # ClinVar
    if a.clinvar_vcf:
        R["clinvar_vcf"] = validate_vcf(a.clinvar_vcf, build, a.full_md5, spans=spans, scope_min=sm)
        R["clinvar_vcf"]["source"] = "override"
    elif human:
        cands = find_files(dirs, [pat_any("clinvar")])
        recs = [validate_vcf(c, build, a.full_md5, spans=spans, scope_min=sm) for c in cands]
        recs.sort(key=fit_first)
        R["clinvar_vcf"] = recs[0] if recs else {"path": None, "status": "missing"}
    else:
        R["clinvar_vcf"] = {"path": None, "status": "n/a"}

    # PhyloP bigWig
    if a.phylop_bw:
        R["phylop_bw"] = {"path": a.phylop_bw, "status": "ok" if os.path.isfile(a.phylop_bw) else "missing",
                          "source": "override"}
    else:
        cands = find_files(dirs, [lambda low: low.endswith((".bw", ".bigwig")) and "phylop" in low])
        R["phylop_bw"] = ({"path": cands[0], "status": "ok", **fingerprint(cands[0], a.full_md5)}
                          if cands else {"path": None, "status": "missing"})

    # BQSR known-sites
    if a.known_sites:
        ks = [validate_vcf(p.strip(), build, a.full_md5) for p in a.known_sites.split(",") if p.strip()]
        R["known_sites"] = {"status": "ok" if all(k["status"] in ("ok", "unindexed") for k in ks) else "incomplete",
                            "files": ks, "source": "override"}
    elif human:
        cands = find_files(dirs, [pat_any("mills", "1000g", "1000G".lower(), "known_indels", "dbsnp")])
        recs = [validate_vcf(c, build, a.full_md5) for c in cands]
        R["known_sites"] = {"status": "ok" if recs else "missing",
                            "files": recs} if recs else {"status": "missing", "files": []}
    else:
        R["known_sites"] = {"status": "n/a", "files": []}

    # Collect found-but-unfit annotation DBs (the silent-corruption case).
    unfit = {name: rec for name, rec in R.items()
             if isinstance(rec, dict) and rec.get("status") == "unfit_scope"}
    manifest["unfit_resources"] = sorted(unfit.keys())

    with open(a.out_json, "w") as fh:
        json.dump(manifest, fh, indent=2)

    # human-readable summary to stdout
    print("[discover_resources] resource manifest:")
    for name, rec in R.items():
        st = rec.get("status", "?")
        path = rec.get("path") if isinstance(rec, dict) else None
        extra = ""
        if name == "gnomad_vcf" and rec.get("has_afr") is not None:
            extra = f" afr={rec.get('has_afr')}"
        if name == "vep_cache" and rec.get("release"):
            extra = f" release={rec.get('release')}"
        sc = rec.get("scope") if isinstance(rec, dict) else None
        if sc and sc.get("scope_checked"):
            extra += f" scope={sc['contigs_covered']}/{sc['contigs_total']}"
        if name == "known_sites":
            path = f"{len(rec.get('files', []))} file(s)"
        mark = {"ok": "FOUND", "missing": "MISSING", "n/a": "n/a",
                "unfit_scope": "UNFIT-SCOPE"}.get(st, st.upper())
        print(f"  - {name:<12} [{mark}]{extra}  {path or ''}")
    print(f"[discover_resources] wrote {a.out_json}")

    if unfit:
        msg = ("; ".join(f"{n}: {R[n].get('notes','region-subset')}" for n in unfit))
        sys.stderr.write(f"[discover_resources] UNFIT (found but region-limited): {msg}\n")
        if a.fail_on_unfit:
            sys.stderr.write("[discover_resources] FAIL: a found annotation DB does not span "
                             "the target panel. Replace it with a genome-wide resource, or pass "
                             "--ignore-unfit / set allow_download. Refusing to under-annotate.\n")
            sys.exit(2)


if __name__ == "__main__":
    main()

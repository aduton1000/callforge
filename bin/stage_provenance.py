#!/usr/bin/env python3
"""stage_provenance.py — per-stage provenance for a CallForge stage subcommand.

Stdlib-only. Writes a `<stage>_provenance.json` describing a single standalone stage
run so it is self-describing and reproducible: inputs (path + sha256 for cheap-to-hash
files), tool versions, the params that drove the stage, the CallForge version and git
commit, and a timestamp.

  stage_provenance.py --stage anno --out anno_provenance.json \
      --input vcf=paralog.annotated.vcf.gz --input reference=test_genome.fa \
      --param vep_mode=gtf --param gnomad_af_field=AF_afr \
      --tool bcftools --tool samtools \
      --version 0.1.0 --git-commit abc1234 --profile test
"""
import argparse, hashlib, json, os, shutil, subprocess, sys, datetime

MAX_HASH_BYTES = 200 * 1024 * 1024   # skip sha256 for very large files (cheap-only)

def sha256(path):
    try:
        if os.path.getsize(path) > MAX_HASH_BYTES:
            return "skipped_large"
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()
    except OSError as e:
        return f"unreadable:{e}"

def tool_version(tool):
    exe = shutil.which(tool)
    if not exe:
        return "not_found"
    for flag in ("--version", "version", "-version"):
        try:
            r = subprocess.run([tool, flag], capture_output=True, text=True, timeout=30)
            txt = (r.stdout or r.stderr).strip().splitlines()
            if txt:
                return txt[0].strip()
        except (subprocess.SubprocessError, OSError):
            continue
    return "unknown"

def kv(pairs):
    out = {}
    for p in pairs:
        if "=" in p:
            k, v = p.split("=", 1)
            out[k] = v
    return out

def main(argv=None):
    ap = argparse.ArgumentParser(prog="stage_provenance.py")
    ap.add_argument("--stage", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--input", action="append", default=[], help="name=path (repeatable)")
    ap.add_argument("--param", action="append", default=[], help="key=value (repeatable)")
    ap.add_argument("--tool", action="append", default=[], help="tool to version (repeatable)")
    ap.add_argument("--version", default="")
    ap.add_argument("--git-commit", default="")
    ap.add_argument("--profile", default="")
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args(argv)

    inputs = {}
    for name, path in kv(a.input).items():
        base = os.path.basename(path)
        inputs[name] = {"path": path, "basename": base,
                        "sha256": sha256(path) if os.path.isfile(path) else "missing"}

    prov = {
        "stage": a.stage,
        "callforge_version": a.version,
        "git_commit": a.git_commit,
        "profile": a.profile,
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "invocation": f"callforge {a.stage}",
        "inputs": inputs,
        "params": kv(a.param),
        "tool_versions": {t: tool_version(t) for t in a.tool},
    }
    out_path = os.path.join(a.outdir, a.out)
    with open(out_path, "w") as fh:
        json.dump(prov, fh, indent=2)
    print(f"[stage_provenance] wrote {out_path} ({len(inputs)} input(s), {len(a.tool)} tool(s))")
    return 0

if __name__ == "__main__":
    sys.exit(main())

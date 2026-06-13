#!/usr/bin/env bash
# ── test/make_test_data.sh ───────────────────────────────────────────────────
# Generate the CallForge `test` profile fixture (gitignored; ~seconds).
# Needs the local GRCh38 no-alt FASTA + a CaptureForge target BED + samtools/
# bgzip/tabix on PATH. Override the two source paths via env vars if they move.
#
#   bash test/make_test_data.sh
#   nextflow run main.nf -profile test            # then run the pipeline
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"

GENOME="${CALLFORGE_TEST_GENOME:-/path/to/GRCh38_noalt.fa}"
BED="${CALLFORGE_TEST_BED:-/path/to/final_covered_targets.bed}"

[ -s "$GENOME" ]      || { echo "ERROR: GRCh38 FASTA not found: $GENOME" >&2; exit 1; }
[ -s "$GENOME.fai" ]  || samtools faidx "$GENOME"
[ -s "$BED" ]         || { echo "ERROR: target BED not found: $BED" >&2; exit 1; }

python3 "$ROOT/bin/make_test_reference.py" \
    --genome "$GENOME" --bed "$BED" \
    --ingest "$ROOT/bin/ingest_captureforge.py" \
    --outdir "$ROOT/test/data" --depth 30

echo "[make_test_data] done -> $ROOT/test/data"

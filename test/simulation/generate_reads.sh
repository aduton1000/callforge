#!/usr/bin/env bash
# generate_reads.sh — build the CallForge simulation fixture on the cluster (in-container).
#
# build_fixture.py needs samtools/bgzip/tabix (htslib). On an HPC node these live in the
# CallForge .sif, so this wrapper runs the generator INSIDE the container. The cohort plan
# (simulate_cohort.py) is pure-stdlib and is run on the host for the reviewable artifacts.
#
# Usage:
#   CALLFORGE_SIF=/path/to/callforge.sif ./generate_reads.sh <outdir> [engine]
#     <outdir>  where the fixture (ref/, fastq/, resources/, truth_manifest.json) is written.
#               Use scratch/project space — NOT inside the read-only install.
#     engine    singularity (default) | apptainer
#
# After this, run the pipeline against <outdir>/sim.params.yaml, then check_against_truth.py
# (see README.md "How to run on the cluster").
set -euo pipefail

OUTDIR="${1:?usage: generate_reads.sh <outdir> [singularity|apptainer]}"
ENGINE="${2:-singularity}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "${HERE}/../.." && pwd)"
OUTDIR="$(mkdir -p "$OUTDIR" && cd "$OUTDIR" && pwd)"

# 1) reviewable plan + sample sheet (host, no deps)
python3 "${HERE}/simulate_cohort.py" --outdir "${OUTDIR}" --fixture-dir "${OUTDIR}"

# 2) reference + reads + truth manifest (in-container — needs htslib)
if [[ -z "${CALLFORGE_SIF:-}" ]]; then
  echo "[generate_reads] CALLFORGE_SIF not set; trying host python (needs samtools/bgzip/tabix on PATH)" >&2
  python3 "${HERE}/build_fixture.py" --outdir "${OUTDIR}"
else
  case "$ENGINE" in
    singularity|apptainer) ;;
    *) echo "[generate_reads] ERROR: engine must be singularity|apptainer (got '$ENGINE')" >&2; exit 2 ;;
  esac
  echo "[generate_reads] running build_fixture.py inside ${CALLFORGE_SIF} via ${ENGINE}" >&2
  "$ENGINE" exec --bind "${REPO}" --bind "${OUTDIR}" "${CALLFORGE_SIF}" \
    python3 "${HERE}/build_fixture.py" --outdir "${OUTDIR}"
fi

echo "[generate_reads] fixture ready: ${OUTDIR}"
echo "  next: callforge-run -params-file ${OUTDIR}/sim.params.yaml -profile hpc_slurm,${ENGINE} (fill the FILL_ME params first)"
echo "  then: python3 ${HERE}/check_against_truth.py --manifest ${OUTDIR}/truth_manifest.json --results ${OUTDIR}/results"

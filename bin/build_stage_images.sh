#!/usr/bin/env bash
# =============================================================================
# build_stage_images.sh — build the per-stage container images for the tools
# that are NOT in the core callforge.sif:
#     annotate  vcfanno (+ensembl-vep)       env/annotate.yml
#     cnv       CNVkit + DNAcopy             env/cnv.yml
#     str       ExpansionHunter              env/str.yml
#     cohortqc  somalier, peddy, verifyBamID2 env/cohortqc.yml
#     burden    regenie, plink2, R SKAT/STAAR env/burden.yml
#     glnexus   GLnexus                      env/glnexus.yml
# Each becomes <out>/callforge-<stage>.sif, which the pipeline resolves via
# --stage_images_dir (callforge-run sets it from CALLFORGE_IMAGES).
#
#   bash bin/build_stage_images.sh --out /hpc/opt/callforge/images            # docker -> SIF
#   bash bin/build_stage_images.sh --out DIR --engine apptainer               # apptainer --fakeroot
#   bash bin/build_stage_images.sh --out DIR --only cnv,str --rebuild
#   bash bin/build_stage_images.sh --out DIR --docker-only                    # tags only (Mac dev)
#
# Options:
#   --out DIR         where the .sif files go (required unless --docker-only)
#   --engine E        docker (default; build natively, convert with apptainer/singularity
#                     build docker-daemon://) | apptainer (build from a .def with --fakeroot)
#   --only a,b        subset of stages (default: all six)
#   --rebuild         rebuild even if the .sif exists
#   --docker-only     build the callforge-<stage>:<version> docker tags and stop
#   --dry-run         print the plan
# Same base as the core image (mambaorg/micromamba); each env is installed into the
# image's base env; `ps` (procps-ng) is in every env because Nextflow needs it.
# =============================================================================
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VERSION="$(grep -oE "version *= *'[^']+'" "$REPO/nextflow.config" | head -1 | sed -E "s/.*'([^']+)'/\1/")"
OUT=""; ENGINE=docker; ONLY=""; REBUILD=0; DOCKER_ONLY=0; DRY=0
while [ $# -gt 0 ]; do case "$1" in
  --out) OUT="$2"; shift;; --engine) ENGINE="$2"; shift;; --only) ONLY="$2"; shift;;
  --rebuild) REBUILD=1;; --docker-only) DOCKER_ONLY=1;; --dry-run) DRY=1;;
  -h|--help) sed -n '2,30p' "$0"; exit 0;; *) echo "unknown option $1" >&2; exit 2;;
esac; shift; done
STAGES=(annotate cnv str cohortqc burden glnexus)
[ -n "$ONLY" ] && IFS=',' read -r -a STAGES <<<"$ONLY"
[ "$DOCKER_ONLY" = 1 ] || [ -n "$OUT" ] || { echo "--out DIR is required" >&2; exit 2; }
log(){ printf '\n\033[1;34m[stage-images]\033[0m %s\n' "$*"; }
run(){ echo "  + $*"; [ "$DRY" = 1 ] || "$@"; }
has(){ command -v "$1" >/dev/null 2>&1; }

RT=""; has apptainer && RT=apptainer; [ -z "$RT" ] && has singularity && RT=singularity
case "$ENGINE" in
  docker)    has docker && docker info >/dev/null 2>&1 || { echo "docker not usable — try --engine apptainer" >&2; exit 1; }
             [ "$DOCKER_ONLY" = 1 ] || [ -n "$RT" ] || { echo "apptainer/singularity needed to convert docker images to .sif" >&2; exit 1; } ;;
  apptainer) [ -n "$RT" ] || { echo "apptainer/singularity not on PATH" >&2; exit 1; } ;;
  *) echo "--engine must be docker or apptainer" >&2; exit 2;;
esac
[ -n "$OUT" ] && run mkdir -p "$OUT"
TMP="$(mktemp -d "${TMPDIR:-/tmp}/cf_stage_images.XXXXXX")"; trap 'rm -rf "$TMP"' EXIT
export APPTAINER_TMPDIR="${APPTAINER_TMPDIR:-${OUT:-$TMP}/.tmp}"; [ "$DRY" = 1 ] || mkdir -p "$APPTAINER_TMPDIR"

# tool each image must be able to run (build-time smoke test)
smoke(){ case "$1" in
  annotate) echo 'vcfanno 2>&1 | head -1; python3 --version; ps --version | head -1';;
  cnv)      echo 'cnvkit.py version; ps --version | head -1';;
  str)      echo 'ExpansionHunter --version 2>&1 | head -1; ps --version | head -1';;
  cohortqc) echo 'somalier --version 2>&1 | head -1; peddy --version 2>&1 | head -1; ps --version | head -1';;
  burden)   echo 'regenie --version 2>&1 | head -1; plink2 --version | head -1; Rscript -e "library(SKAT); library(jsonlite)"; ps --version | head -1';;
  glnexus)  echo 'glnexus_cli --version 2>&1 | head -1; bcftools --version | head -1; ps --version | head -1';;
esac; }

for st in "${STAGES[@]}"; do
  yml="$REPO/env/$st.yml"; [ -f "$yml" ] || { echo "no env/$st.yml" >&2; exit 1; }
  tag="callforge-$st:$VERSION"; sif="${OUT:+$OUT/callforge-$st.sif}"
  if [ -n "$sif" ] && [ -s "$sif" ] && [ "$REBUILD" = 0 ]; then log "$st: have $sif — skip (--rebuild to force)"; continue; fi
  log "$st  ($yml -> ${sif:-$tag})"
  if [ "$ENGINE" = docker ]; then
    cat > "$TMP/Dockerfile.$st" <<EOF
FROM mambaorg/micromamba:1.5.8
LABEL org.opencontainers.image.title="CallForge stage image: $st" org.opencontainers.image.version="$VERSION"
USER root
ENV LC_ALL=C.UTF-8 LANG=C.UTF-8 MPLCONFIGDIR=/tmp/matplotlib
COPY $st.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml && micromamba clean --all --yes
ENV PATH=/opt/conda/bin:\$PATH
RUN $(smoke "$st")
WORKDIR /work
EOF
    [ "$DRY" = 1 ] || cp "$yml" "$TMP/$st.yml"
    run docker build --platform linux/amd64 -t "$tag" -f "$TMP/Dockerfile.$st" "$TMP"
    if [ "$DOCKER_ONLY" = 0 ]; then
      [ "$DRY" = 1 ] || rm -f "$sif"
      run "$RT" build "$sif" "docker-daemon://$tag"
    fi
  else
    cat > "$TMP/$st.def" <<EOF
Bootstrap: docker
From: mambaorg/micromamba:1.5.8
%files
    $yml /tmp/env.yml
%post
    export LC_ALL=C.UTF-8 LANG=C.UTF-8
    micromamba install -y -n base -f /tmp/env.yml && micromamba clean --all --yes
%environment
    export PATH=/opt/conda/bin:\$PATH LC_ALL=C.UTF-8 LANG=C.UTF-8 MPLCONFIGDIR=/tmp/matplotlib
%test
    export PATH=/opt/conda/bin:\$PATH
    $(smoke "$st")
%labels
    Name CallForge-$st
    Version $VERSION
EOF
    [ "$DRY" = 1 ] || rm -f "$sif"
    run "$RT" build --fakeroot "$sif" "$TMP/$st.def"
  fi
  if [ -n "$sif" ] && [ "$DOCKER_ONLY" = 0 ] && [ "$DRY" = 0 ]; then
    log "$st: verify inside $sif"
    "$RT" exec "$sif" bash -c "$(smoke "$st")" | sed 's/^/    /'
    chmod 644 "$sif" 2>/dev/null || true
  fi
done
log "done. images:"; [ -n "$OUT" ] && ls -la "$OUT"/callforge-*.sif 2>/dev/null || true
[ -n "$OUT" ] && echo "  export CALLFORGE_IMAGES=$OUT   (or stage_images_dir: $OUT in the params file)"

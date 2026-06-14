#!/usr/bin/env bash
# =============================================================================
# fetch_annotation_dbs.sh
# Resumable downloader for genome-wide gnomAD v4.1 + ClinVar (GRCh38).
# Suggested location: ~/callforge/bin/fetch_annotation_dbs.sh
#
# WHY: the locally discovered gnomAD/ClinVar were PAH-region subsets — unfit for a
# multi-gene panel. This fetches genome-wide replacements.
#
# RESUMABLE: each file is fetched with `wget --continue`; after every attempt the
# BGZF/gzip stream is integrity-checked, and a truncated file (e.g. from a dropped
# connection) RESUMES from where it left off instead of restarting. A real 404 is
# detected and skipped rather than retried forever.
#
# MODES:
#   full   (default) – download the entire genome-wide DB (HUGE: ~several hundred GB
#                      to ~1 TB for genomes+exomes). This is the literal "full database".
#   panel            – tabix-slice ONLY your panel's target regions from the remote
#                      files (a few hundred MB; same data at your loci). Recommended
#                      for a targeted panel. Requires PANEL_BED.
#
# USAGE:
#   bash fetch_annotation_dbs.sh                 # full mode, ./annotation_db
#   MODE=panel PANEL_BED=/path/final_covered_targets.bed bash fetch_annotation_dbs.sh
#   DATASETS="genomes" MODE=full bash fetch_annotation_dbs.sh   # genomes only (smaller)
# =============================================================================
set -uo pipefail

# ---------------- configuration (override via env) ----------------
OUTDIR="${OUTDIR:-$(pwd)/annotation_db}"
MODE="${MODE:-full}"                       # full | panel
DATASETS="${DATASETS:-genomes exomes}"     # which gnomAD sets; "genomes" alone is fine + smaller
GNOMAD_VER="${GNOMAD_VER:-4.1}"
CHROMS="${CHROMS:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y}"
PANEL_BED="${PANEL_BED:-}"                 # required for MODE=panel (Ensembl-named, e.g. CaptureForge final_covered_targets.bed)
MAX_ATTEMPTS="${MAX_ATTEMPTS:-30}"
GNOMAD_BASE="https://storage.googleapis.com/gcp-public-data--gnomad/release/${GNOMAD_VER}/vcf"
CLINVAR_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"

mkdir -p "$OUTDIR"
LOG="$OUTDIR/fetch_$(date +%Y%m%d_%H%M%S).log"
FAILS=()

# ---------------- helpers ----------------
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found in PATH" >&2; exit 127; }; }
need wget; need bgzip; need tabix
[ "$MODE" = panel ] && { need bcftools; [ -s "$PANEL_BED" ] || { echo "MODE=panel needs a valid PANEL_BED" >&2; exit 2; }; }

log() { echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG"; }

integrity_ok() {  # validate a finished download
  local f="$1"
  [ -s "$f" ] || return 1
  case "$f" in
    *.tbi|*.csi) return 0 ;;                                  # tiny index; presence is enough
    *.bgz|*.gz)  bgzip -t "$f" 2>/dev/null || gzip -t "$f" 2>/dev/null ;;
    *)           return 0 ;;
  esac
}

# resumable fetch with integrity loop; returns 0 ok, 8 on 404/server-error, 1 on give-up
fetch() {
  local url="$1" out="$2" n=0
  while :; do
    n=$((n+1))
    wget --continue --tries=20 --retry-connrefused --waitretry=15 \
         --timeout=60 --read-timeout=300 --no-verbose \
         -O "$out" "$url" >>"$LOG" 2>&1
    local rc=$?
    if [ $rc -eq 8 ]; then log "  server error (404?) -> skip: $url"; return 8; fi
    if integrity_ok "$out"; then log "  ok: $(basename "$out")"; return 0; fi
    if [ $n -ge "$MAX_ATTEMPTS" ]; then log "  GAVE UP after $n attempts: $out"; return 1; fi
    log "  incomplete/corrupt; resuming ($n/$MAX_ATTEMPTS): $(basename "$out")"; sleep 10
  done
}

# ---------------- ClinVar (small; always full) ----------------
log "=== ClinVar (GRCh38, latest) -> $OUTDIR ==="
fetch "$CLINVAR_URL"     "$OUTDIR/clinvar.vcf.gz"     || FAILS+=("clinvar.vcf.gz")
fetch "$CLINVAR_URL.tbi" "$OUTDIR/clinvar.vcf.gz.tbi" || FAILS+=("clinvar.vcf.gz.tbi")

# ---------------- gnomAD ----------------
if [ "$MODE" = full ]; then
  log "=== gnomAD v${GNOMAD_VER} FULL genome-wide [$DATASETS] — this is very large (hundreds of GB+). Ctrl-C to abort. ==="
  for ds in $DATASETS; do
    for c in $CHROMS; do
      base="gnomad.${ds}.v${GNOMAD_VER}.sites.chr${c}.vcf.bgz"
      url="${GNOMAD_BASE}/${ds}/${base}"
      log "-- ${ds} chr${c}"
      fetch "$url"     "$OUTDIR/$base"     || FAILS+=("$base")
      fetch "$url.tbi" "$OUTDIR/$base.tbi" || FAILS+=("$base.tbi")
    done
  done

else  # ---- panel mode: remote tabix-slice only the target regions ----
  log "=== gnomAD v${GNOMAD_VER} PANEL subset from $PANEL_BED [$DATASETS] ==="
  # gnomAD is chr-prefixed; panel BED is Ensembl-named -> build a chr-prefixed BED for slicing.
  chrbed="$OUTDIR/.panel.chr.bed"
  awk 'BEGIN{OFS="\t"} {c=$1; if(c !~ /^chr/) c="chr"c; print c,$2,$3}' "$PANEL_BED" | sort -k1,1 -k2,2n > "$chrbed"
  for ds in $DATASETS; do
    out="$OUTDIR/gnomad.${ds}.v${GNOMAD_VER}.panel.vcf.gz"
    : > "$OUTDIR/.${ds}.parts"
    hdr_done=0
    for c in $CHROMS; do
      url="${GNOMAD_BASE}/${ds}/gnomad.${ds}.v${GNOMAD_VER}.sites.chr${c}.vcf.bgz"
      grep -qP "^chr${c}\t" "$chrbed" || continue
      log "-- slicing ${ds} chr${c} (remote tabix)"
      part="$OUTDIR/.${ds}.chr${c}.vcf"
      if [ $hdr_done -eq 0 ]; then tabix -h "$url" -R "$chrbed" > "$part" 2>>"$LOG" && hdr_done=1; \
      else tabix "$url" -R "$chrbed" > "$part" 2>>"$LOG"; fi || { log "  slice failed chr${c}"; FAILS+=("${ds}.chr${c}.slice"); continue; }
      echo "$part" >> "$OUTDIR/.${ds}.parts"
    done
    cat $(cat "$OUTDIR/.${ds}.parts") | bgzip > "$out" && tabix -p vcf "$out" \
      && log "  wrote $out" || FAILS+=("$out")
    rm -f $(cat "$OUTDIR/.${ds}.parts") "$OUTDIR/.${ds}.parts"
  done
  rm -f "$chrbed"
fi

# ---------------- manifest + summary ----------------
log "=== writing manifest ==="
{ echo -e "file\tbytes\tmd5(optional)"; for f in "$OUTDIR"/*.vcf.gz "$OUTDIR"/*.bgz; do
    [ -e "$f" ] || continue; printf "%s\t%s\t%s\n" "$(basename "$f")" "$(wc -c <"$f")" "${VERIFY_MD5:+$(md5sum "$f" 2>/dev/null | cut -d' ' -f1)}"; done; } > "$OUTDIR/MANIFEST.tsv"

if [ ${#FAILS[@]} -eq 0 ]; then
  log "DONE — all downloads complete and integrity-checked. See $OUTDIR/MANIFEST.tsv"
else
  log "DONE WITH ISSUES — re-run to resume these (will continue, not restart):"; printf '  %s\n' "${FAILS[@]}" | tee -a "$LOG"
  exit 1
fi

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
MERGE_DIST="${MERGE_DIST:-1000}"           # panel mode: merge intervals within this gap before slicing (fewer HTTP requests)
MAX_ATTEMPTS="${MAX_ATTEMPTS:-30}"
GNOMAD_BASE="https://storage.googleapis.com/gcp-public-data--gnomad/release/${GNOMAD_VER}/vcf"
CLINVAR_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"

mkdir -p "$OUTDIR"
LOG="$OUTDIR/fetch_$(date +%Y%m%d_%H%M%S).log"
FAILS=()

# ---------------- helpers ----------------
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found in PATH" >&2; exit 127; }; }
need bgzip; need tabix
# downloader: wget preferred, curl accepted (the callforge.sif image ships curl, not wget)
DL=""; command -v wget >/dev/null 2>&1 && DL=wget; [ -z "$DL" ] && command -v curl >/dev/null 2>&1 && DL=curl
[ -n "$DL" ] || { echo "ERROR: neither 'wget' nor 'curl' found in PATH" >&2; exit 127; }
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
  # already complete (e.g. copied in from elsewhere)? don't re-download
  if integrity_ok "$out"; then log "  present, skip: $(basename "$out")"; return 0; fi
  while :; do
    n=$((n+1))
    local rc
    if [ "$DL" = wget ]; then
      wget --continue --tries=20 --retry-connrefused --waitretry=15 \
           --timeout=60 --read-timeout=300 --no-verbose \
           -O "$out" "$url" >>"$LOG" 2>&1; rc=$?
    else
      # curl: -C - resumes, -f makes HTTP errors a non-zero rc (22 -> treated like wget's 8)
      curl -fsSL -C - --retry 20 --retry-delay 15 --retry-all-errors \
           --connect-timeout 60 -o "$out" "$url" >>"$LOG" 2>&1; rc=$?
      [ $rc -eq 22 ] && rc=8
    fi
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
  # MERGE nearby intervals so remote tabix issues a few large range-requests instead of
  # one per interval. The panel's CNV bins sit ~50 bp apart, so an unmerged BED forces
  # thousands of HTTP round-trips (request-latency bound, hours). Merging -> minutes.
  # (Extra intronic gnomAD pulled in is harmless: vcfanno keeps only matching sites.)
  if command -v bedtools >/dev/null 2>&1; then
    bedtools merge -d "$MERGE_DIST" -i "$chrbed" > "${chrbed}.m" && mv "${chrbed}.m" "$chrbed"
    log "panel BED merged to $(wc -l < "$chrbed" | tr -d ' ') ranges (-d $MERGE_DIST) for efficient remote slicing"
  else
    warn "bedtools not found — slicing unmerged intervals (may be very slow over the network)"
  fi
  for ds in $DATASETS; do
    out="$OUTDIR/gnomad.${ds}.v${GNOMAD_VER}.panel.vcf.gz"
    : > "$OUTDIR/.${ds}.parts"
    hdr_done=0
    for c in $CHROMS; do
      # Per-chromosome range subset (one small BED per chrom) for remote tabix.
      cbed="$OUTDIR/.${ds}.chr${c}.bed"
      awk -v c="chr${c}" '$1==c' "$chrbed" > "$cbed"
      [ -s "$cbed" ] || { rm -f "$cbed"; continue; }
      url="${GNOMAD_BASE}/${ds}/gnomad.${ds}.v${GNOMAD_VER}.sites.chr${c}.vcf.bgz"
      part="$OUTDIR/.${ds}.chr${c}.vcf"; done_marker="$OUTDIR/.${ds}.chr${c}.done"
      # Resume: skip a chromosome already sliced in a prior run.
      if [ -f "$done_marker" ] && [ -s "$part" ]; then
        log "-- ${ds} chr${c} already sliced (resume)"; echo "$part" >> "$OUTDIR/.${ds}.parts"; hdr_done=1; rm -f "$cbed"; continue
      fi
      log "-- slicing ${ds} chr${c} ($(wc -l < "$cbed" | tr -d ' ') ranges, remote tabix)"
      ok=0
      for attempt in 1 2 3 4 5; do
        if [ $hdr_done -eq 0 ]; then tabix -h "$url" -R "$cbed" > "$part" 2>>"$LOG"; else tabix "$url" -R "$cbed" > "$part" 2>>"$LOG"; fi
        if [ $? -eq 0 ]; then ok=1; break; fi
        log "  chr${c} slice attempt ${attempt} failed (transient?); retry in 15s"; sleep 15
      done
      rm -f "$cbed"
      if [ $ok -eq 1 ]; then touch "$done_marker"; hdr_done=1; echo "$part" >> "$OUTDIR/.${ds}.parts"
      else log "  slice FAILED chr${c} after retries"; FAILS+=("${ds}.chr${c}.slice"); fi
    done
    cat $(cat "$OUTDIR/.${ds}.parts") | bgzip > "$out" && tabix -p vcf "$out" \
      && log "  wrote $out" || FAILS+=("$out")
    # Keep parts + .done markers if any failure (so a re-run resumes); else clean up.
    if [ ${#FAILS[@]} -eq 0 ]; then rm -f $(cat "$OUTDIR/.${ds}.parts") "$OUTDIR/.${ds}.parts" "$OUTDIR"/.${ds}.chr*.done; fi
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

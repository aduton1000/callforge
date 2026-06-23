#!/usr/bin/env bash
# =============================================================================
# fetch_annotation_refs.sh — stage the REAL annotation/QC references CallForge needs,
# with RESUMABLE downloads + per-file COMPLETENESS VERIFICATION. Runs on the CLUSTER
# (needs egress). Idempotent: a file already present + verified is skipped.
#
# WHY: the real-reference simulation harness (test/simulation/) and the real cohort run
# both validate against REAL databases — nothing faked. The pipeline SPINE
# (align -> QC -> GATK call -> CNV -> STR) needs ONLY the genome (+ .dict/+bwa-mem2 index
# built on the cluster) + the STR catalog the harness provides, so you can validate the
# spine BEFORE these DBs finish downloading. Later stages need the DBs below.
#
# REFERENCE -> STAGE(S) THAT NEED IT
#   vep_cache        annotation (11)         ~15-20 GB (indexed cache, one release)
#   phylop_bigwig    annotation (11)         ~9-10 GB (100-way conservation)
#   gnomad_afr       burden (14), annotation hundreds of GB full; use MODE=panel to slice (~100s MB)
#   dbsnp            annotation (11), BQSR    ~25 GB (or Ensembl per-chrom)
#   clinvar          annotation (11)          ~200 MB
#   somalier_sites   cohort QC (12)           ~10 MB
#   giab_hg002       (real cohort GIAB only)  ~50 MB (VCF+BED) — NOT used by the sim
#   gatk_known       BQSR (stage 4)           Mills ~80 MB (+ dbSNP above)
#   ensembl_gtf      read-gen placement + VEP(gtf) ~50 MB
#
# CONTIG-NAMING CAVEAT: CallForge uses Ensembl no-alt naming (1,2,…,X — no `chr`). Every
# DB here MUST match. Sources marked [chr] are chr-prefixed (UCSC/GATK-hg38) — rename or
# pick an Ensembl-named equivalent, or the pipeline's reference-invariant check fails.
# Override any URL via the env vars below; `RELEASE` sets the Ensembl/VEP release.
#
# USAGE:
#   OUTDIR=/refs/callforge ./bin/fetch_annotation_refs.sh                 # all refs
#   OUTDIR=/refs/callforge ONLY="vep_cache clinvar somalier_sites" ./bin/fetch_annotation_refs.sh
#   MODE=panel PANEL_BED=/path/targets.bed ./bin/fetch_annotation_refs.sh # slice gnomAD to a panel
# =============================================================================
set -uo pipefail

OUTDIR="${OUTDIR:-$(pwd)/callforge_refs}"
RELEASE="${RELEASE:-110}"                 # Ensembl / VEP cache release
SPECIES="${SPECIES:-homo_sapiens}"
ASSEMBLY="${ASSEMBLY:-GRCh38}"
MODE="${MODE:-full}"                      # full | panel  (panel slices gnomAD to PANEL_BED)
PANEL_BED="${PANEL_BED:-}"
ONLY="${ONLY:-vep_cache phylop_bigwig gnomad_afr dbsnp clinvar somalier_sites giab_hg002 gatk_known ensembl_gtf}"
MAX_ATTEMPTS="${MAX_ATTEMPTS:-40}"

# ---- source URLs (override via env; marked [chr] = chr-prefixed, needs rename) ----
VEP_URL="${VEP_URL:-https://ftp.ensembl.org/pub/release-${RELEASE}/variation/indexed_vep_cache/${SPECIES}_vep_${RELEASE}_${ASSEMBLY}.tar.gz}"
GTF_URL="${GTF_URL:-https://ftp.ensembl.org/pub/release-${RELEASE}/gtf/${SPECIES}/Homo_sapiens.${ASSEMBLY}.${RELEASE}.gtf.gz}"
CLINVAR_URL="${CLINVAR_URL:-https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_${ASSEMBLY}/clinvar.vcf.gz}"
PHYLOP_URL="${PHYLOP_URL:-https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw}"   # [chr]
SOMALIER_SITES_URL="${SOMALIER_SITES_URL:-https://github.com/brentp/somalier/files/3412456/sites.hg38.nochr.vcf.gz}"  # Ensembl-named
GIAB_VCF_URL="${GIAB_VCF_URL:-https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz}"
GIAB_BED_URL="${GIAB_BED_URL:-https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed}"
MILLS_URL="${MILLS_URL:-https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz}"   # [chr]
DBSNP_URL="${DBSNP_URL:-https://ftp.ensembl.org/pub/release-${RELEASE}/variation/vcf/${SPECIES}/${SPECIES}.vcf.gz}"
GNOMAD_BASE="${GNOMAD_BASE:-https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes}"   # [chr]
CHROMS="${CHROMS:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y}"

mkdir -p "$OUTDIR"
LOG="$OUTDIR/fetch_refs_$(date +%Y%m%d_%H%M%S).log"
FAILS=()
log() { echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG"; }
need() { command -v "$1" >/dev/null 2>&1 || { log "ERROR: '$1' not in PATH"; exit 127; }; }
have() { command -v "$1" >/dev/null 2>&1; }
need wget

# ---- structural / checksum verification: a file is "done" only if this passes ----
verify() {
  local f="$1" method="structural"
  [ -s "$f" ] || return 1
  # MD5 sidecar published next to the source? (compare if we fetched one)
  if [ -s "$f.md5" ] && have md5sum; then
    local want got; want=$(awk '{print $1}' "$f.md5"); got=$(md5sum "$f" | awk '{print $1}')
    [ "$want" = "$got" ] || { log "  MD5 MISMATCH: $(basename "$f")"; return 1; }
    method="md5"
  fi
  case "$f" in
    *.vcf.gz|*.bcf)
      if have bcftools; then bcftools view -h "$f" >/dev/null 2>&1 || return 1
      else gzip -t "$f" 2>/dev/null || return 1; fi
      [ "$method" = md5 ] || method="bcftools-header" ;;
    *.tar.gz|*.tgz) tar tzf "$f" >/dev/null 2>&1 || return 1; [ "$method" = md5 ] || method="tar-listing" ;;
    *.bw|*.bigWig)  if have bigWigInfo; then bigWigInfo "$f" >/dev/null 2>&1 || return 1; method="bigWigInfo"; fi ;;
    *.gz)           gzip -t "$f" 2>/dev/null || return 1; [ "$method" = md5 ] || method="gzip-t" ;;
    *.bed)          [ -s "$f" ] || return 1 ;;
  esac
  echo "$method"; return 0
}

# resumable fetch + verify loop. rc: 0 ok, 8 server-error(404), 1 give-up
fetch() {
  local url="$1" out="$2" n=0 m
  if m=$(verify "$out" 2>/dev/null); then log "  skip (already verified: ${m:-ok}): $(basename "$out")"; return 0; fi
  while :; do
    n=$((n+1))
    wget --continue --tries=20 --retry-connrefused --waitretry=15 \
         --timeout=60 --read-timeout=600 --no-verbose -O "$out" "$url" >>"$LOG" 2>&1
    local rc=$?
    if [ $rc -eq 8 ]; then log "  server error (404?) -> skip: $url"; return 8; fi
    if m=$(verify "$out"); then
      local sz; sz=$(du -h "$out" 2>/dev/null | cut -f1)
      log "  downloaded + verified (${m}): $(basename "$out") [$sz]"; return 0
    fi
    if [ $n -ge "$MAX_ATTEMPTS" ]; then log "  GAVE UP after $n attempts: $out"; return 1; fi
    log "  incomplete/corrupt; RESUMING ($n/$MAX_ATTEMPTS): $(basename "$out")"; sleep 10
  done
}

# fetch an optional .md5 sidecar first (so verify() can use it), then the file + index
fetch_md5() { wget -q -O "$2.md5" "$1" 2>/dev/null || rm -f "$2.md5"; }
index_vcf() { have tabix && { [ -s "$1.tbi" ] || tabix -p vcf "$1" >>"$LOG" 2>&1 || tabix -C -p vcf "$1" >>"$LOG" 2>&1; }; }

wants() { case " $ONLY " in *" $1 "*) return 0;; *) return 1;; esac; }

# ----------------------------- per-reference staging -----------------------------
if wants vep_cache; then
  log "=== VEP cache ${SPECIES} ${RELEASE} ${ASSEMBLY} (~15-20 GB) ==="
  out="$OUTDIR/${SPECIES}_vep_${RELEASE}_${ASSEMBLY}.tar.gz"
  fetch_md5 "$VEP_URL.md5" "$out"
  if fetch "$VEP_URL" "$out"; then
    [ -d "$OUTDIR/vep/${SPECIES}" ] || { log "  extracting VEP cache…"; mkdir -p "$OUTDIR/vep"; tar -xzf "$out" -C "$OUTDIR/vep" >>"$LOG" 2>&1 && log "  vep_cache ready -> $OUTDIR/vep"; }
  else FAILS+=("vep_cache"); fi
fi

if wants ensembl_gtf; then
  log "=== Ensembl GTF (gene model for read-gen placement + VEP gtf-mode) (~50 MB) ==="
  fetch "$GTF_URL" "$OUTDIR/$(basename "$GTF_URL")" || FAILS+=("ensembl_gtf")
fi

if wants clinvar; then
  log "=== ClinVar ${ASSEMBLY} (~200 MB) ==="
  if fetch "$CLINVAR_URL" "$OUTDIR/clinvar.vcf.gz"; then
    fetch "$CLINVAR_URL.tbi" "$OUTDIR/clinvar.vcf.gz.tbi" || index_vcf "$OUTDIR/clinvar.vcf.gz"
  else FAILS+=("clinvar"); fi
fi

if wants phylop_bigwig; then
  log "=== PhyloP 100-way bigWig (~9-10 GB) [chr — rename chroms if your ref is Ensembl-named] ==="
  fetch "$PHYLOP_URL" "$OUTDIR/$(basename "$PHYLOP_URL")" || FAILS+=("phylop_bigwig")
fi

if wants somalier_sites; then
  log "=== somalier sites ${ASSEMBLY} (Ensembl-named) (~10 MB) ==="
  if fetch "$SOMALIER_SITES_URL" "$OUTDIR/somalier.sites.${ASSEMBLY}.vcf.gz"; then
    index_vcf "$OUTDIR/somalier.sites.${ASSEMBLY}.vcf.gz"
  else FAILS+=("somalier_sites"); fi
fi

if wants giab_hg002; then
  log "=== GIAB HG002 truth (VCF+BED) — for the REAL-cohort GIAB benchmark, NOT the sim (~50 MB) ==="
  fetch "$GIAB_VCF_URL" "$OUTDIR/HG002_${ASSEMBLY}.benchmark.vcf.gz" || FAILS+=("giab_vcf")
  fetch "$GIAB_VCF_URL.tbi" "$OUTDIR/HG002_${ASSEMBLY}.benchmark.vcf.gz.tbi" || index_vcf "$OUTDIR/HG002_${ASSEMBLY}.benchmark.vcf.gz"
  fetch "$GIAB_BED_URL" "$OUTDIR/HG002_${ASSEMBLY}.benchmark.bed" || FAILS+=("giab_bed")
fi

if wants gatk_known; then
  log "=== GATK known-sites: Mills indels (~80 MB) [chr — GATK hg38 bundle; rename if Ensembl] ==="
  if fetch "$MILLS_URL" "$OUTDIR/$(basename "$MILLS_URL")"; then
    fetch "$MILLS_URL.tbi" "$OUTDIR/$(basename "$MILLS_URL").tbi" || index_vcf "$OUTDIR/$(basename "$MILLS_URL")"
  else FAILS+=("gatk_known_mills"); fi
fi

if wants dbsnp; then
  log "=== dbSNP (Ensembl ${SPECIES} variation VCF; ~25 GB) ==="
  if fetch "$DBSNP_URL" "$OUTDIR/dbsnp.${ASSEMBLY}.vcf.gz"; then
    fetch "$DBSNP_URL.tbi" "$OUTDIR/dbsnp.${ASSEMBLY}.vcf.gz.tbi" || index_vcf "$OUTDIR/dbsnp.${ASSEMBLY}.vcf.gz"
  else FAILS+=("dbsnp"); fi
fi

if wants gnomad_afr; then
  if [ "$MODE" = panel ]; then
    [ -s "$PANEL_BED" ] || { log "MODE=panel needs PANEL_BED"; FAILS+=("gnomad_afr"); }
    if [ -s "$PANEL_BED" ]; then
      need tabix; need bgzip
      log "=== gnomAD v4.1 genomes — PANEL slice from $PANEL_BED (~100s MB) [chr source] ==="
      chrbed="$OUTDIR/.gnomad.panel.chr.bed"
      awk 'BEGIN{OFS="\t"}{c=$1; if(c!~/^chr/)c="chr"c; print c,$2,$3}' "$PANEL_BED" | sort -k1,1 -k2,2n > "$chrbed"
      have bedtools && { bedtools merge -d 1000 -i "$chrbed" > "$chrbed.m" && mv "$chrbed.m" "$chrbed"; }
      out="$OUTDIR/gnomad.genomes.v4.1.AFR.panel.vcf.gz"; : > "$OUTDIR/.gnomad.parts"; hdr=0
      for c in $CHROMS; do
        cb="$OUTDIR/.gnomad.chr${c}.bed"; awk -v k="chr${c}" '$1==k' "$chrbed" > "$cb"
        [ -s "$cb" ] || { rm -f "$cb"; continue; }
        url="${GNOMAD_BASE}/gnomad.genomes.v4.1.sites.chr${c}.vcf.bgz"
        part="$OUTDIR/.gnomad.chr${c}.vcf"; dm="$OUTDIR/.gnomad.chr${c}.done"
        if [ -f "$dm" ] && [ -s "$part" ]; then log "-- chr${c} sliced (resume)"; echo "$part" >> "$OUTDIR/.gnomad.parts"; hdr=1; rm -f "$cb"; continue; fi
        ok=0; for at in 1 2 3 4 5; do
          if [ $hdr -eq 0 ]; then tabix -h "$url" -R "$cb" > "$part" 2>>"$LOG"; else tabix "$url" -R "$cb" > "$part" 2>>"$LOG"; fi
          [ $? -eq 0 ] && { ok=1; break; }; log "  chr${c} slice retry ${at}"; sleep 15; done
        rm -f "$cb"
        [ $ok -eq 1 ] && { touch "$dm"; hdr=1; echo "$part" >> "$OUTDIR/.gnomad.parts"; } || FAILS+=("gnomad.chr${c}")
      done
      cat $(cat "$OUTDIR/.gnomad.parts") | bgzip > "$out" && index_vcf "$out" \
        && log "  downloaded + verified (tabix-slice): $(basename "$out")" || FAILS+=("gnomad_afr")
      rm -f "$chrbed"
    fi
  else
    log "=== gnomAD v4.1 genomes — FULL genome-wide is HUGE (hundreds of GB). Use MODE=panel ==="
    log "    PANEL_BED=<targets.bed> MODE=panel for a ~100s-MB slice with the AFR AF field, OR"
    log "    reuse bin/fetch_annotation_dbs.sh (genome-wide). Skipping full gnomAD here."
  fi
fi

# ----------------------------- summary -----------------------------
{ echo -e "file\tbytes"; for f in "$OUTDIR"/*.vcf.gz "$OUTDIR"/*.bw "$OUTDIR"/*.tar.gz "$OUTDIR"/*.bed; do
    [ -e "$f" ] || continue; printf "%s\t%s\n" "$(basename "$f")" "$(wc -c <"$f")"; done; } > "$OUTDIR/REFS_MANIFEST.tsv" 2>/dev/null

if [ ${#FAILS[@]} -eq 0 ]; then
  log "DONE — refs staged + verified. See $OUTDIR/REFS_MANIFEST.tsv"
  log "SPINE needs only the genome; annotation/burden/cohortqc use the DBs above (see header map)."
else
  log "DONE WITH ISSUES — re-run to RESUME these (continues, never restarts):"; printf '  %s\n' "${FAILS[@]}" | tee -a "$LOG"
  exit 1
fi

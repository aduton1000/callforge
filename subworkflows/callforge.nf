// ── subworkflows/callforge.nf ────────────────────────────────────────────────
// Top-level CallForge orchestration.
//   Phase 1: Stage 0  (inputs, resource discovery, reference invariant, CF ingest)
//   Phase 2: Stages 1-5 (raw QC -> trim -> align -> dedup/BQSR -> HsMetrics/coverage
//                        -> per-sample QC gate & quarantine)
// Later phases (calling -> CNV/STR/paralog -> annotation -> cohort QC -> GIAB ->
// burden -> reporting) attach at the marked extension point.

include { VALIDATE_SAMPLESHEET; PREPARE_REFERENCE; REFERENCE_INVARIANT;
          DISCOVER_RESOURCES;  INGEST_CAPTUREFORGE; COUNT_READS;
          PLOT_STAGE0 } from '../modules/stage0_inputs.nf'

include { FASTQC_RAW; FASTP; BWAMEM2_INDEX; BWAMEM2_ALIGN; MARKDUP; BQSR;
          BED_TO_INTERVALS; HSMETRICS; SAMTOOLS_STATS; MOSDEPTH;
          SUMMARIZE_SAMPLE_QC; QC_GATE; PLOT_ALIGN_QC; PLOT_COVERAGE;
          PLOT_QC_GATE } from '../modules/stage1_5_align_qc.nf'

include { HAPLOTYPECALLER; GENOMICSDB_IMPORT; GENOTYPE_GVCFS_DB; COMBINE_GENOTYPE;
          HARD_FILTER; CALLSET_STATS; FILTER_SUMMARY; PLOT_CALLING_QC; PLOT_FILTER_QC;
          DEEPVARIANT; GLNEXUS } from '../modules/stage6_7_calling.nf'

include { CNVKIT_BATCH; CNV_ANNOTATE; PLOT_CNV; BUILD_STR_CATALOG; EXPANSIONHUNTER;
          STR_SUMMARIZE; PLOT_STR; PARALOG_FLAG; PLOT_PARALOG } from '../modules/stage8_10_cnv_str_paralog.nf'

include { VEP; ANNOTATE_DBS; VERIFY_ANNOTATION; FLATTEN_TSV; PLOT_ANNOTATION } from '../modules/stage11_annotate.nf'

include { SOMALIER_EXTRACT; SOMALIER_RELATE; COHORT_QC; EXTRACT_CONTROL; GIAB_HAPPY;
          PLOT_GIAB } from '../modules/stage12_13_cohortqc_giab.nf'

include { BURDEN_MATRIX; BURDEN_COLLAPSE; BURDEN_REGENIE; BURDEN_SKAT;
          PLOT_BURDEN } from '../modules/stage14_burden.nf'

include { MULTIQC; DASHBOARD; PROVENANCE } from '../modules/stage15_report.nf'

// ── Pre-built reference reuse helpers ────────────────────────────────────────
// The bwa-mem2 index, .fai and .dict are properties of the REFERENCE, not the
// cohort — building them every run wastes 40-60 min on the full genome. Detect a
// co-located pre-built copy and stage it in so it is reused instead of rebuilt.

// Expected bwa-mem2 index files for the configured genome (next to the FASTA, or in
// params.genome_index_dir using the FASTA's basename so `bwa-mem2 mem <fasta>` finds them).
def bwaIdxFiles() {
    def base = file(params.genome_fasta).name
    def dir  = params.genome_index_dir ? file(params.genome_index_dir) : file(params.genome_fasta).parent
    return ['0123', 'amb', 'ann', 'bwt.2bit.64', 'pac'].collect { file("${dir}/${base}.${it}") }
}

// Co-located pre-built .fai / .dict (or the matching NO_* placeholder) to STAGE into
// PREPARE_REFERENCE. Named canonically so the in-process guard reuses them as-is.
def coRef(ext) {
    def g    = file(params.genome_fasta)
    def cand = (ext == 'dict') ? file("${g.parent}/${g.baseName}.dict")
                               : file("${params.genome_fasta}.${ext}")
    return cand.exists() ? cand : file("${projectDir}/assets/NO_${ext.toUpperCase()}")
}

workflow CALLFORGE {
    take:
    ch_reads        // tuple(sample_id, fastq_1, fastq_2)
    ch_sheet        // path : raw sample sheet
    ch_genome       // path : reference FASTA
    ch_target_bed   // path : final_covered_targets.bed
    ch_cf_metrics   // path : CaptureForge metrics.json | assets/NO_METRICS
    ch_cf_baits     // path : CaptureForge baits.csv    | assets/NO_BAITS

    main:
    // ── Stage 0 ──────────────────────────────────────────────────────────────
    VALIDATE_SAMPLESHEET( ch_sheet )
    PREPARE_REFERENCE( ch_genome, coRef('fai'), coRef('dict') )
    REFERENCE_INVARIANT( PREPARE_REFERENCE.out.fai, ch_target_bed )
    DISCOVER_RESOURCES( Channel.value(params.resource_dirs), ch_target_bed )
    INGEST_CAPTUREFORGE( ch_target_bed, ch_cf_metrics, ch_cf_baits )

    COUNT_READS( ch_reads )
    ch_readcounts = COUNT_READS.out.tsv
        .collectFile(name: 'input_readcounts.tsv', newLine: false,
                     seed: "sample_id\tread_pairs\n", sort: true)
    PLOT_STAGE0( DISCOVER_RESOURCES.out.manifest, ch_readcounts )

    // ── Stages 1-3 : raw QC -> trim -> align ─────────────────────────────────
    FASTQC_RAW( ch_reads )
    FASTP( ch_reads )

    // REUSE a complete co-located bwa-mem2 index; build it once only if missing.
    def idxFiles = bwaIdxFiles()
    def ch_bwa_idx
    if (idxFiles.every { it.exists() }) {
        log.info "[CallForge] Reusing pre-built bwa-mem2 index for ${params.genome_fasta} (BWAMEM2_INDEX skipped)."
        ch_bwa_idx = Channel.value(idxFiles)
    } else {
        log.info "[CallForge] No complete pre-built bwa-mem2 index found — building it once (BWAMEM2_INDEX)."
        BWAMEM2_INDEX( PREPARE_REFERENCE.out.fasta )
        ch_bwa_idx = BWAMEM2_INDEX.out.idx
    }
    BWAMEM2_ALIGN( FASTP.out.reads, PREPARE_REFERENCE.out.fasta, ch_bwa_idx )

    // ── Stage 4 : dedup -> BQSR (graceful skip) -> HsMetrics/coverage ────────
    MARKDUP( BWAMEM2_ALIGN.out.bam )
    BQSR( MARKDUP.out.bam, PREPARE_REFERENCE.out.fasta, PREPARE_REFERENCE.out.fai,
          PREPARE_REFERENCE.out.dict, Channel.value(params.known_sites ?: '') )
    ch_analysis_bam = BQSR.out.bam

    BED_TO_INTERVALS( ch_target_bed, PREPARE_REFERENCE.out.dict )
    HSMETRICS( ch_analysis_bam, PREPARE_REFERENCE.out.fasta, PREPARE_REFERENCE.out.fai,
               BED_TO_INTERVALS.out.il )
    SAMTOOLS_STATS( ch_analysis_bam )
    MOSDEPTH( ch_analysis_bam, ch_target_bed )

    // ── Stage 5 : per-sample QC summary -> cohort gate (flag & quarantine) ───
    ch_sex = Channel.fromPath(params.input).splitCsv(header: true)
        .map { row -> tuple(row.sample_id, (row.sex ?: 'U')) }

    ch_metrics = SAMTOOLS_STATS.out.flagstat
        .join(SAMTOOLS_STATS.out.stats)
        .join(MARKDUP.out.metrics)
        .join(HSMETRICS.out.metrics)
        .join(MOSDEPTH.out.summary)
        .join(MOSDEPTH.out.regions)
    ch_qc_in = ch_sex.join(ch_metrics)   // (sid, sex, flagstat, stats, markdup, hs, mdsum, mdreg)

    SUMMARIZE_SAMPLE_QC( ch_qc_in )
    ch_qc_json = SUMMARIZE_SAMPLE_QC.out.json.collect()
    QC_GATE( ch_qc_json )

    PLOT_ALIGN_QC( ch_qc_json, SUMMARIZE_SAMPLE_QC.out.ishist.collect().ifEmpty([]) )
    PLOT_COVERAGE( ch_qc_json, MOSDEPTH.out.regiondist.map { it[1] }.collect() )
    PLOT_QC_GATE( QC_GATE.out.summary, ch_qc_json )

    // QC-passing sample IDs (quarantined samples are excluded from joint calling
    // but kept + reported). Joined back to analysis BAMs for Phase 3.
    ch_pass_ids = QC_GATE.out.pass.splitText().map { it.trim() }.filter { it }
    ch_pass_bams = ch_analysis_bam
        .map { sid, bam, bai -> tuple(sid, bam, bai) }
        .join( ch_pass_ids.map { sid -> tuple(sid, true) } )
        .map { sid, bam, bai, ok -> tuple(sid, bam, bai) }

    // ── Stages 6-7 : SNV/indel joint calling -> hard-filter ──────────────────
    // Only QC-PASS BAMs enter here; quarantined samples are excluded from the
    // GenomicsDB / joint genotyping step (kept + reported by the gate).
    def ch_joint
    if (params.caller == 'deepvariant') {
        DEEPVARIANT( ch_pass_bams, PREPARE_REFERENCE.out.fasta, PREPARE_REFERENCE.out.fai, ch_target_bed )
        ch_dv = DEEPVARIANT.out.gvcf
        GLNEXUS( ch_dv.map { s, g, t -> g }.collect(), ch_dv.map { s, g, t -> t }.collect(), ch_target_bed )
        ch_joint = GLNEXUS.out.vcf
    } else {
        HAPLOTYPECALLER( ch_pass_bams, PREPARE_REFERENCE.out.fasta, PREPARE_REFERENCE.out.fai,
                         PREPARE_REFERENCE.out.dict, ch_target_bed )
        ch_g = HAPLOTYPECALLER.out.gvcf.map { s, g, t -> g }.collect()
        ch_t = HAPLOTYPECALLER.out.gvcf.map { s, g, t -> t }.collect()
        if (params.joint_method == 'combinegvcfs') {
            COMBINE_GENOTYPE( ch_g, ch_t, PREPARE_REFERENCE.out.fasta, PREPARE_REFERENCE.out.fai,
                              PREPARE_REFERENCE.out.dict, ch_target_bed )
            ch_joint = COMBINE_GENOTYPE.out.vcf
        } else {
            GENOMICSDB_IMPORT( ch_g, ch_t, ch_target_bed )
            GENOTYPE_GVCFS_DB( GENOMICSDB_IMPORT.out.gdb, PREPARE_REFERENCE.out.fasta,
                               PREPARE_REFERENCE.out.fai, PREPARE_REFERENCE.out.dict, ch_target_bed )
            ch_joint = GENOTYPE_GVCFS_DB.out.vcf
        }
    }

    CALLSET_STATS( ch_joint )
    HARD_FILTER( ch_joint, PREPARE_REFERENCE.out.fasta, PREPARE_REFERENCE.out.fai, PREPARE_REFERENCE.out.dict )
    FILTER_SUMMARY( HARD_FILTER.out.vcf )
    PLOT_CALLING_QC( CALLSET_STATS.out.json, FILTER_SUMMARY.out.qualdp )
    PLOT_FILTER_QC( FILTER_SUMMARY.out.counts )

    // ── Stages 8-10 : CNV + STR + paralog-aware (CaptureForge-metadata driven) ──
    ch_gene_meta = INGEST_CAPTUREFORGE.out.tsv
    ch_pass_bam_files = ch_pass_bams.map { s, b, i -> b }.collect()
    ch_pass_bai_files = ch_pass_bams.map { s, b, i -> i }.collect()

    // Stage 8 : CNV (CNVkit) with CaptureForge callability labels propagated
    CNVKIT_BATCH( ch_pass_bam_files, ch_pass_bai_files, ch_target_bed,
                  PREPARE_REFERENCE.out.fasta, PREPARE_REFERENCE.out.fai )
    CNV_ANNOTATE( CNVKIT_BATCH.out.cns, CNVKIT_BATCH.out.cnr, ch_target_bed, ch_gene_meta )
    PLOT_CNV( CNV_ANNOTATE.out.copyratio, CNV_ANNOTATE.out.summary, ch_gene_meta )

    // Stage 9 : STR (ExpansionHunter) — curated catalog if supplied, else built
    def ch_str_catalog
    if (params.str_catalog) {
        ch_str_catalog = Channel.value(file(params.str_catalog))
    } else {
        BUILD_STR_CATALOG( ch_gene_meta )
        ch_str_catalog = BUILD_STR_CATALOG.out.catalog
    }
    EXPANSIONHUNTER( ch_pass_bams, PREPARE_REFERENCE.out.fasta, PREPARE_REFERENCE.out.fai, ch_str_catalog )
    STR_SUMMARIZE( EXPANSIONHUNTER.out.json.collect() )
    PLOT_STR( STR_SUMMARIZE.out.calls, STR_SUMMARIZE.out.summary )

    // Stage 10 : paralog-aware flagging on the filtered callset
    PARALOG_FLAG( HARD_FILTER.out.vcf, ch_target_bed, ch_gene_meta )
    PLOT_PARALOG( PARALOG_FLAG.out.summary )

    // ── Stage 11 : annotation (VEP + vcfanno, manifest-driven, chr-reconciled) ──
    def ch_vep_aux = null
    def ch_vep_aux_idx = null
    if (params.vep_mode == 'gtf') {
        ch_vep_aux     = Channel.value(file(params.gtf))
        ch_vep_aux_idx = Channel.value(file(params.gtf + '.tbi'))
    } else {
        ch_vep_aux     = Channel.value(file(params.vep_cache ?: "${projectDir}/assets/NO_CACHE"))
        ch_vep_aux_idx = Channel.value(file("${projectDir}/assets/NO_GTF"))
    }
    VEP( PARALOG_FLAG.out.vcf, PREPARE_REFERENCE.out.fasta, PREPARE_REFERENCE.out.fai,
         ch_vep_aux, ch_vep_aux_idx )
    ANNOTATE_DBS( VEP.out.vcf, PREPARE_REFERENCE.out.fai, DISCOVER_RESOURCES.out.manifest )
    VERIFY_ANNOTATION( ANNOTATE_DBS.out.vcf, ANNOTATE_DBS.out.sources )
    FLATTEN_TSV( ANNOTATE_DBS.out.vcf )
    PLOT_ANNOTATION( FLATTEN_TSV.out.tsv, VERIFY_ANNOTATION.out.landing )

    // ── Stage 12 : cohort QC (somalier relatedness/sex + ancestry PCA + missingness) ──
    if (params.somalier_sites) {
        SOMALIER_EXTRACT( ch_pass_bams, PREPARE_REFERENCE.out.fasta, PREPARE_REFERENCE.out.fai,
                          Channel.value(file(params.somalier_sites)),
                          Channel.value(file(params.somalier_sites + '.tbi')) )
        SOMALIER_RELATE( SOMALIER_EXTRACT.out.somalier.collect() )
        COHORT_QC( HARD_FILTER.out.vcf, SOMALIER_RELATE.out.samples, SOMALIER_RELATE.out.pairs,
                   VALIDATE_SAMPLESHEET.out.csv )
    } else {
        log.warn "No --somalier_sites: Stage 12 cohort QC (relatedness/ancestry/sex) skipped."
    }

    // ── Stage 13 : GIAB benchmarking (hap.py, restricted to the panel BED) ──
    if (params.giab_control_id && params.giab_truth_vcf && params.giab_truth_bed) {
        EXTRACT_CONTROL( HARD_FILTER.out.vcf )
        GIAB_HAPPY( EXTRACT_CONTROL.out.vcf,
                    Channel.value(file(params.giab_truth_vcf)),
                    Channel.value(file(params.giab_truth_vcf + '.tbi')),
                    Channel.value(file(params.giab_truth_bed)),
                    ch_target_bed, PREPARE_REFERENCE.out.fasta, PREPARE_REFERENCE.out.fai )
        PLOT_GIAB( GIAB_HAPPY.out.summary )
    } else {
        log.warn "No GIAB control/truth: Stage 13 benchmarking skipped (set giab_control_id + giab_truth_vcf/bed)."
    }

    // ── Stage 14 : rare-variant burden / association (gated on phenotype) ──
    def ch_bres = null
    def ch_burden_meta = Channel.value(file("${projectDir}/assets/NO_CACHE"))
    if (params.run_burden) {
        def ch_cohortqc = params.somalier_sites ? COHORT_QC.out.json
                                                : Channel.value(file("${projectDir}/assets/NO_CACHE"))
        BURDEN_MATRIX( ANNOTATE_DBS.out.vcf, INGEST_CAPTUREFORGE.out.tsv,
                       VALIDATE_SAMPLESHEET.out.csv, ch_cohortqc )
        ch_burden_meta = BURDEN_MATRIX.out.meta
        def ch_bsum
        if (params.burden_engine == 'collapse') {
            BURDEN_COLLAPSE( BURDEN_MATRIX.out.matrix, BURDEN_MATRIX.out.pheno )
            ch_bres = BURDEN_COLLAPSE.out.results; ch_bsum = BURDEN_COLLAPSE.out.summary
        } else if (params.burden_engine == 'skat') {
            BURDEN_SKAT( BURDEN_MATRIX.out.genotypes, BURDEN_MATRIX.out.pheno, BURDEN_MATRIX.out.covar )
            ch_bres = BURDEN_SKAT.out.results; ch_bsum = BURDEN_SKAT.out.summary
        } else {
            BURDEN_REGENIE( BURDEN_MATRIX.out.genotypes, BURDEN_MATRIX.out.pheno, BURDEN_MATRIX.out.covar )
            ch_bres = BURDEN_REGENIE.out.results; ch_bsum = BURDEN_REGENIE.out.summary
        }
        PLOT_BURDEN( ch_bres, ch_bsum )
    } else {
        log.warn "run_burden=false: Stage 14 burden/association skipped."
    }

    // ── Stage 15 : reporting (MultiQC + self-contained cohort QC dashboard + provenance) ──
    // Collect every captioned plot + key summary across stages (guarded for the
    // conditional stages) into the dashboard; per-tool outputs into MultiQC.
    def ch_png = PLOT_STAGE0.out.png.mix(
        PLOT_ALIGN_QC.out.png, PLOT_COVERAGE.out.png, PLOT_QC_GATE.out.png,
        PLOT_CALLING_QC.out.png, PLOT_FILTER_QC.out.png, PLOT_CNV.out.png,
        PLOT_STR.out.png, PLOT_PARALOG.out.png, PLOT_ANNOTATION.out.png)
    def ch_cap = PLOT_STAGE0.out.captions.mix(
        PLOT_ALIGN_QC.out.captions, PLOT_COVERAGE.out.captions, PLOT_QC_GATE.out.captions,
        PLOT_CALLING_QC.out.captions, PLOT_FILTER_QC.out.captions, PLOT_CNV.out.captions,
        PLOT_STR.out.captions, PLOT_PARALOG.out.captions, PLOT_ANNOTATION.out.captions)
    def ch_sum = DISCOVER_RESOURCES.out.manifest.mix(
        QC_GATE.out.scorecard, VERIFY_ANNOTATION.out.landing, CNV_ANNOTATE.out.summary,
        STR_SUMMARIZE.out.summary, PARALOG_FLAG.out.summary)
    if (params.somalier_sites) {
        ch_png = ch_png.mix(COHORT_QC.out.png); ch_cap = ch_cap.mix(COHORT_QC.out.captions)
        ch_sum = ch_sum.mix(COHORT_QC.out.json)
    }
    if (params.giab_control_id && params.giab_truth_vcf && params.giab_truth_bed) {
        ch_png = ch_png.mix(PLOT_GIAB.out.png); ch_cap = ch_cap.mix(PLOT_GIAB.out.captions)
        ch_sum = ch_sum.mix(PLOT_GIAB.out.metrics)
    }
    if (params.run_burden && ch_bres) {
        ch_png = ch_png.mix(PLOT_BURDEN.out.png); ch_cap = ch_cap.mix(PLOT_BURDEN.out.captions)
        ch_sum = ch_sum.mix(ch_bres)
    }
    DASHBOARD( ch_png.collect(), ch_cap.collect(), ch_sum.collect() )

    def ch_mqc = FASTQC_RAW.out.zip.mix(
        FASTP.out.json, MARKDUP.out.metrics.map { it[1] }, HSMETRICS.out.metrics.map { it[1] },
        SAMTOOLS_STATS.out.stats.map { it[1] }, MOSDEPTH.out.summary.map { it[1] },
        CALLSET_STATS.out.stats)
    MULTIQC( ch_mqc.collect() )

    PROVENANCE( DISCOVER_RESOURCES.out.manifest, VALIDATE_SAMPLESHEET.out.summary,
                QC_GATE.out.summary, QC_GATE.out.quarantine, ch_burden_meta )

    emit:
    fasta         = PREPARE_REFERENCE.out.fasta
    fai           = PREPARE_REFERENCE.out.fai
    dict          = PREPARE_REFERENCE.out.dict
    manifest      = DISCOVER_RESOURCES.out.manifest
    gene_metadata = INGEST_CAPTUREFORGE.out.tsv
    analysis_bams = ch_analysis_bam
    pass_bams     = ch_pass_bams
    qc_scorecard  = QC_GATE.out.scorecard
    qc_summary    = QC_GATE.out.summary
    sheet_summary = VALIDATE_SAMPLESHEET.out.summary
    joint_vcf     = ch_joint
    filtered_vcf  = HARD_FILTER.out.vcf
    calling_stats = CALLSET_STATS.out.json
    filter_counts = FILTER_SUMMARY.out.counts
    cnv_calls     = CNV_ANNOTATE.out.calls
    cnv_summary   = CNV_ANNOTATE.out.summary
    str_calls     = STR_SUMMARIZE.out.calls
    paralog_vcf   = PARALOG_FLAG.out.vcf
    paralog_summary = PARALOG_FLAG.out.summary
    annotated_vcf   = ANNOTATE_DBS.out.vcf
    variants_tsv    = FLATTEN_TSV.out.tsv
    annotation_landing = VERIFY_ANNOTATION.out.landing
    dashboard       = DASHBOARD.out.html
    multiqc         = MULTIQC.out.report
    provenance      = PROVENANCE.out.json
}

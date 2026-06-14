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
    PREPARE_REFERENCE( ch_genome )
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
    BWAMEM2_INDEX( PREPARE_REFERENCE.out.fasta )
    BWAMEM2_ALIGN( FASTP.out.reads, PREPARE_REFERENCE.out.fasta, BWAMEM2_INDEX.out.idx )

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

    // ══════════════════════ PHASE 5+ EXTENSION POINT ═════════════════════════
    // ANNOTATE( PARALOG_FLAG.out.vcf, DISCOVER_RESOURCES.out.manifest ) ->
    // COHORT_QC -> GIAB -> BURDEN( sheet_summary gates on phenotype ) -> REPORT.

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
}

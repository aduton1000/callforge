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

    // ══════════════════════ PHASE 3+ EXTENSION POINT ═════════════════════════
    // JOINT_CALL( ch_pass_bams, reference, target_bed ) -> FILTER -> CNV/STR/
    // PARALOG -> ANNOTATE( DISCOVER_RESOURCES.out.manifest ) -> COHORT_QC ->
    // GIAB -> BURDEN( VALIDATE_SAMPLESHEET.out.summary gates on phenotype ) -> REPORT.

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
}

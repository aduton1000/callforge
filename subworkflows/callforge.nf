// ── subworkflows/callforge.nf ────────────────────────────────────────────────
// Top-level CallForge orchestration. Phase 1 wires Stage 0 (inputs & resource
// discovery + reference invariant + CaptureForge ingestion + Stage-0 QC plots)
// and exposes its outputs as channels. Later phases (align→QC gate→calling→
// CNV/STR/paralog→annotation→cohort QC→GIAB→burden→reporting) attach below the
// marked extension point, consuming these same channels.

include { VALIDATE_SAMPLESHEET; PREPARE_REFERENCE; REFERENCE_INVARIANT;
          DISCOVER_RESOURCES;  INGEST_CAPTUREFORGE; COUNT_READS;
          PLOT_STAGE0 } from '../modules/stage0_inputs.nf'

workflow CALLFORGE {
    take:
    ch_reads        // tuple(sample_id, fastq_1, fastq_2)
    ch_sheet        // path : raw sample sheet
    ch_genome       // path : reference FASTA
    ch_target_bed   // path : final_covered_targets.bed
    ch_cf_metrics   // path : CaptureForge metrics.json | assets/NO_FILE
    ch_cf_baits     // path : CaptureForge baits.csv    | assets/NO_FILE

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

    // ══════════════════════ PHASE 2+ EXTENSION POINT ═════════════════════════
    // ALIGN( ch_reads, PREPARE_REFERENCE.out ) -> DEDUP/BQSR -> HSMETRICS ->
    // QC_GATE (quarantine) -> JOINT_CALL -> FILTER -> CNV/STR/PARALOG ->
    // ANNOTATE (DISCOVER_RESOURCES.out.manifest) -> COHORT_QC -> GIAB ->
    // BURDEN (VALIDATE_SAMPLESHEET.out.summary gates on phenotype) -> REPORT.

    emit:
    samplesheet   = VALIDATE_SAMPLESHEET.out.csv
    sheet_summary = VALIDATE_SAMPLESHEET.out.summary
    fasta         = PREPARE_REFERENCE.out.fasta
    fai           = PREPARE_REFERENCE.out.fai
    dict          = PREPARE_REFERENCE.out.dict
    invariant     = REFERENCE_INVARIANT.out.report
    manifest      = DISCOVER_RESOURCES.out.manifest
    gene_metadata = INGEST_CAPTUREFORGE.out.tsv
    readcounts    = ch_readcounts
}

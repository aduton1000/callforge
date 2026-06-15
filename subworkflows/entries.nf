// ── subworkflows/entries.nf  (CallForge stage-subcommand entry workflows) ────────
// Named entry workflows so a user can run ONE analytical stage standalone via
//   nextflow run main.nf -entry <stage> --<inputs…>   (or `callforge <stage> …`)
// on externally-provided inputs (the files that normally arrive from an upstream
// stage) instead of re-running the whole pipeline.
//
// NON-DIVERGENCE RULE: every entry invokes the SAME stage modules the full pipeline
// uses (subworkflows/callforge.nf) — it only swaps upstream channel outputs for
// user-supplied inputs, adds fail-loud input validation, and writes a focused stage
// report + provenance. No stage analysis logic is reimplemented here.
//
// Built so far (Gate 1): ANNO, BURDEN. Others attach the same way in later gates.

include { PREPARE_REFERENCE; REFERENCE_INVARIANT; DISCOVER_RESOURCES;
          INGEST_CAPTUREFORGE; VALIDATE_SAMPLESHEET } from '../modules/stage0_inputs.nf'
include { FASTQC_RAW; FASTP; BWAMEM2_INDEX; BWAMEM2_ALIGN; MARKDUP; BQSR;
          BED_TO_INTERVALS; HSMETRICS; SAMTOOLS_STATS; MOSDEPTH; SUMMARIZE_SAMPLE_QC;
          QC_GATE; PLOT_ALIGN_QC; PLOT_COVERAGE; PLOT_QC_GATE } from '../modules/stage1_5_align_qc.nf'
include { HAPLOTYPECALLER; GENOMICSDB_IMPORT; GENOTYPE_GVCFS_DB; COMBINE_GENOTYPE;
          HARD_FILTER; CALLSET_STATS; FILTER_SUMMARY; PLOT_CALLING_QC; PLOT_FILTER_QC;
          DEEPVARIANT; GLNEXUS } from '../modules/stage6_7_calling.nf'
include { VEP; ANNOTATE_DBS; VERIFY_ANNOTATION; FLATTEN_TSV;
          PLOT_ANNOTATION } from '../modules/stage11_annotate.nf'
include { BURDEN_MATRIX; BURDEN_COLLAPSE; BURDEN_REGENIE; BURDEN_SKAT;
          PLOT_BURDEN } from '../modules/stage14_burden.nf'
include { VALIDATE_STAGE_INPUTS; VALIDATE_BAM; STAGE_PUBLISH;
          STAGE_REPORT; STAGE_PROVENANCE } from '../modules/stage_entry.nf'

// helper: require a param, else fail loud with a clear message
def need(cond, msg) { if (!cond) { exit 1, "ERROR: ${msg}" } }

// helper: build a (sample_id, bam, bai) channel from a glob or comma list of BAM paths.
// sample_id = filename up to the first dot (e.g. S1.analysis.bam -> S1); .bai must exist.
def bamChannel(spec) {
    def ch = spec.contains(',') ? Channel.fromPath(spec.split(',') as List, checkIfExists: true)
                                : Channel.fromPath(spec, checkIfExists: true)
    return ch.map { bam -> tuple(bam.name.split('\\.')[0], bam, file("${bam}.bai", checkIfExists: true)) }
}

// ── anno : re-annotate a VCF (VEP + vcfanno DBs + PhyloP), e.g. after a DB update ──
workflow ANNO {
    main:
    need(params.input_vcf,    "callforge anno needs --input_vcf <VCF.gz> (a paralog-flagged or filtered VCF)")
    need(params.genome_fasta, "callforge anno needs --genome_fasta <no-alt primary assembly>")
    need(params.target_bed,   "callforge anno needs --target_bed (reference invariant + DB scope gate)")
    need(file(params.input_vcf).exists(),           "--input_vcf not found: ${params.input_vcf}")
    need(file(params.input_vcf + '.tbi').exists(),  "VCF index missing: ${params.input_vcf}.tbi (run `tabix -p vcf`)")

    ch_vcf    = Channel.value( tuple(file(params.input_vcf), file(params.input_vcf + '.tbi')) )
    ch_genome = file(params.genome_fasta, checkIfExists: true)
    ch_bed    = file(params.target_bed,   checkIfExists: true)
    def noref = file("${projectDir}/assets/NO_CACHE")
    def synth = (workflow.profile?.contains('test')) ? '--synthetic' : ''

    // ── validation: reference invariant (no-alt + build/contig) + DB scope gate ──
    PREPARE_REFERENCE( ch_genome )
    REFERENCE_INVARIANT( PREPARE_REFERENCE.out.fai, ch_bed )
    DISCOVER_RESOURCES( Channel.value(params.resource_dirs), ch_bed )
    VALIDATE_STAGE_INPUTS( Channel.value('anno'), ch_vcf, PREPARE_REFERENCE.out.fai,
                           Channel.value(noref), Channel.value('') )

    // ── same Stage-11 module chain as the pipeline (VEP cache|gtf identical) ──
    def ch_vep_aux; def ch_vep_aux_idx
    if (params.vep_mode == 'gtf') {
        ch_vep_aux     = Channel.value(file(params.gtf))
        ch_vep_aux_idx = Channel.value(file(params.gtf + '.tbi'))
    } else {
        ch_vep_aux     = Channel.value(file(params.vep_cache ?: "${projectDir}/assets/NO_CACHE"))
        ch_vep_aux_idx = Channel.value(file("${projectDir}/assets/NO_GTF"))
    }
    VEP( VALIDATE_STAGE_INPUTS.out.vcf, PREPARE_REFERENCE.out.fasta, PREPARE_REFERENCE.out.fai,
         ch_vep_aux, ch_vep_aux_idx )
    ANNOTATE_DBS( VEP.out.vcf, PREPARE_REFERENCE.out.fai, DISCOVER_RESOURCES.out.manifest )
    VERIFY_ANNOTATION( ANNOTATE_DBS.out.vcf, ANNOTATE_DBS.out.sources )
    FLATTEN_TSV( ANNOTATE_DBS.out.vcf )
    PLOT_ANNOTATION( FLATTEN_TSV.out.tsv, VERIFY_ANNOTATION.out.landing )

    // ── focused report + provenance ──
    def ch_sum = VERIFY_ANNOTATION.out.landing
                    .mix(FLATTEN_TSV.out.tsv, ANNOTATE_DBS.out.sources).collect()
    STAGE_REPORT( Channel.value('anno'), Channel.value('Annotation (Stage 11)'),
                  Channel.value('stage11_annotation'),
                  PLOT_ANNOTATION.out.png.collect(), PLOT_ANNOTATION.out.captions.collect(),
                  ch_sum, Channel.value(synth) )
    STAGE_PROVENANCE( Channel.value('anno'), Channel.value('stage11_annotation'),
        Channel.value("--input vcf=${params.input_vcf} --input reference=${params.genome_fasta} " +
                      "--param vep_mode=${params.vep_mode} --param species=${params.species} " +
                      "--param gnomad_af_field=${params.gnomad_af_field} --param scope_min=${params.scope_min} " +
                      "--tool bcftools --tool samtools") )
}

// ── burden : re-run rare-variant burden with new thresholds/engine on an annotated VCF ──
workflow BURDEN {
    main:
    need(params.input_vcf,  "callforge burden needs --input_vcf <annotated VCF.gz>")
    need(params.input,      "callforge burden needs --input <sample sheet with phenotype + covariate_* columns>")
    need(params.target_bed, "callforge burden needs --target_bed (CaptureForge burden groups via metadata)")
    need(file(params.input_vcf).exists(),          "--input_vcf not found: ${params.input_vcf}")
    need(file(params.input_vcf + '.tbi').exists(), "VCF index missing: ${params.input_vcf}.tbi (run `tabix -p vcf`)")
    need(params.run_burden, "run_burden=false — set --run_burden true to run the burden stage standalone")

    ch_vcf   = Channel.value( tuple(file(params.input_vcf), file(params.input_vcf + '.tbi')) )
    ch_sheet = file(params.input,      checkIfExists: true)
    ch_bed   = file(params.target_bed, checkIfExists: true)
    def cf_metrics = params.cf_metrics_json ? file(params.cf_metrics_json) : file("${projectDir}/assets/NO_METRICS")
    def cf_baits   = params.cf_baits_csv    ? file(params.cf_baits_csv)    : file("${projectDir}/assets/NO_BAITS")
    def noref = file("${projectDir}/assets/NO_CACHE")
    def synth = (workflow.profile?.contains('test')) ? '--synthetic' : ''

    // ── validation: validate sheet, then fail-loud on missing phenotype ──
    VALIDATE_SAMPLESHEET( ch_sheet )
    INGEST_CAPTUREFORGE( ch_bed, cf_metrics, cf_baits )
    VALIDATE_STAGE_INPUTS( Channel.value('burden'), ch_vcf, Channel.value(noref),
                           VALIDATE_SAMPLESHEET.out.csv, Channel.value('--require-phenotype') )

    // ── same Stage-14 module chain as the pipeline (matrix -> engine -> plots) ──
    // ancestry PCs: from a supplied cohort_qc.json, else none (the PC-stability gate
    // inside build_burden_matrix.py applies identically).
    def ch_cohortqc = params.cohort_qc_json ? Channel.value(file(params.cohort_qc_json)) : Channel.value(noref)
    BURDEN_MATRIX( VALIDATE_STAGE_INPUTS.out.vcf, INGEST_CAPTUREFORGE.out.tsv,
                   VALIDATE_SAMPLESHEET.out.csv, ch_cohortqc )
    def ch_bres; def ch_bsum
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

    def ch_sum = ch_bsum.mix(ch_bres, BURDEN_MATRIX.out.meta).collect()
    STAGE_REPORT( Channel.value('burden'), Channel.value('Rare-variant burden (Stage 14)'),
                  Channel.value('stage14_burden'),
                  PLOT_BURDEN.out.png.collect(), PLOT_BURDEN.out.captions.collect(),
                  ch_sum, Channel.value(synth) )
    STAGE_PROVENANCE( Channel.value('burden'), Channel.value('stage14_burden'),
        Channel.value("--input vcf=${params.input_vcf} --input samplesheet=${params.input} " +
                      "--param burden_engine=${params.burden_engine} --param burden_af_max=${params.burden_af_max} " +
                      "--param n_ancestry_pcs=${params.n_ancestry_pcs} --param covariates='${params.covariates}' " +
                      "--tool bcftools") )
}

// ── align : FASTQ -> analysis-ready BAM(s) (+ alignment QC plots) ──────────────
// per-sample. Requires --target_bed because the alignment QC plots are produced via
// the SAME SUMMARIZE_SAMPLE_QC the pipeline uses (which needs HsMetrics/mosdepth).
workflow ALIGN {
    main:
    need(params.genome_fasta, "callforge align needs --genome_fasta <no-alt primary assembly>")
    need(params.target_bed,   "callforge align needs --target_bed (per-sample QC metrics feed the alignment plots)")
    def ch_reads
    if (params.input) {
        ch_reads = Channel.fromPath(params.input, checkIfExists: true).splitCsv(header: true)
            .map { row ->
                need(row.sample_id && row.fastq_1 && row.fastq_2,
                     "sample sheet row missing sample_id/fastq_1/fastq_2: ${row}")
                tuple(row.sample_id, file(row.fastq_1, checkIfExists: true),
                                     file(row.fastq_2, checkIfExists: true)) }
    } else {
        need(params.fastq_1 && params.fastq_2,
             "callforge align needs --input <sheet> OR --fastq_1 + --fastq_2 (+ --sample_id)")
        ch_reads = Channel.of( tuple((params.sample_id ?: 'sample1'),
                       file(params.fastq_1, checkIfExists: true), file(params.fastq_2, checkIfExists: true)) )
    }
    ch_genome = file(params.genome_fasta, checkIfExists: true)
    ch_bed    = file(params.target_bed,   checkIfExists: true)
    def synth = (workflow.profile?.contains('test')) ? '--synthetic' : ''

    PREPARE_REFERENCE( ch_genome )
    REFERENCE_INVARIANT( PREPARE_REFERENCE.out.fai, ch_bed )      // no-alt + build/contig gate
    BWAMEM2_INDEX( PREPARE_REFERENCE.out.fasta )
    FASTQC_RAW( ch_reads )
    FASTP( ch_reads )
    BWAMEM2_ALIGN( FASTP.out.reads, PREPARE_REFERENCE.out.fasta, BWAMEM2_INDEX.out.idx )
    MARKDUP( BWAMEM2_ALIGN.out.bam )
    BQSR( MARKDUP.out.bam, PREPARE_REFERENCE.out.fasta, PREPARE_REFERENCE.out.fai,
          PREPARE_REFERENCE.out.dict, Channel.value(params.known_sites ?: '') )
    ch_analysis = BQSR.out.bam

    BED_TO_INTERVALS( ch_bed, PREPARE_REFERENCE.out.dict )
    HSMETRICS( ch_analysis, PREPARE_REFERENCE.out.fasta, PREPARE_REFERENCE.out.fai, BED_TO_INTERVALS.out.il )
    SAMTOOLS_STATS( ch_analysis )
    MOSDEPTH( ch_analysis, ch_bed )
    ch_sex = params.input ? Channel.fromPath(params.input).splitCsv(header: true)
                                   .map { row -> tuple(row.sample_id, (row.sex ?: 'U')) }
                          : ch_reads.map { sid, a, b -> tuple(sid, 'U') }
    ch_metrics = SAMTOOLS_STATS.out.flagstat.join(SAMTOOLS_STATS.out.stats)
        .join(MARKDUP.out.metrics).join(HSMETRICS.out.metrics)
        .join(MOSDEPTH.out.summary).join(MOSDEPTH.out.regions)
    SUMMARIZE_SAMPLE_QC( ch_sex.join(ch_metrics) )
    ch_qc_json = SUMMARIZE_SAMPLE_QC.out.json.collect()
    PLOT_ALIGN_QC( ch_qc_json, SUMMARIZE_SAMPLE_QC.out.ishist.collect().ifEmpty([]) )

    // the analysis BAM is not publishDir'd by BQSR — publish it as this stage's output
    STAGE_PUBLISH( Channel.value('stage4_postalign'),
                   ch_analysis.map { sid, bam, bai -> [bam, bai] }.flatten().collect() )
    STAGE_REPORT( Channel.value('align'), Channel.value('Alignment (Stages 1-4)'),
                  Channel.value('stage4_postalign'),
                  PLOT_ALIGN_QC.out.png.collect(), PLOT_ALIGN_QC.out.captions.collect(),
                  ch_qc_json, Channel.value(synth) )
    STAGE_PROVENANCE( Channel.value('align'), Channel.value('stage4_postalign'),
        Channel.value("--input reference=${params.genome_fasta} --input target_bed=${params.target_bed} " +
                      "--param known_sites='${params.known_sites ?: ''}' " +
                      "--tool bwa-mem2 --tool samtools --tool gatk --tool fastp") )
}

// ── coverage : BAM(s) -> coverage/enrichment QC + per-sample gate (pass/quarantine) ──
// per-sample metrics -> cohort gate. Treats supplied BAMs as analysis-ready; re-runs
// MarkDuplicates ONLY to derive the duplicate-rate metric the QC summary needs.
workflow COVERAGE {
    main:
    need(params.input_bams,   "callforge coverage needs --input_bams <glob/comma list of BAM(s)>")
    need(params.genome_fasta, "callforge coverage needs --genome_fasta <reference>")
    need(params.target_bed,   "callforge coverage needs --target_bed <panel BED>")
    ch_bams   = bamChannel(params.input_bams)
    ch_genome = file(params.genome_fasta, checkIfExists: true)
    ch_bed    = file(params.target_bed,   checkIfExists: true)
    def synth = (workflow.profile?.contains('test')) ? '--synthetic' : ''

    PREPARE_REFERENCE( ch_genome )
    REFERENCE_INVARIANT( PREPARE_REFERENCE.out.fai, ch_bed )
    ch_v = VALIDATE_BAM( Channel.value('coverage'), ch_bams, PREPARE_REFERENCE.out.fai ).bam

    BED_TO_INTERVALS( ch_bed, PREPARE_REFERENCE.out.dict )
    HSMETRICS( ch_v, PREPARE_REFERENCE.out.fasta, PREPARE_REFERENCE.out.fai, BED_TO_INTERVALS.out.il )
    SAMTOOLS_STATS( ch_v )
    MOSDEPTH( ch_v, ch_bed )
    MARKDUP( ch_v )    // duplicate-rate metric for the QC summary (md.bam itself unused)
    ch_sex = params.input ? Channel.fromPath(params.input).splitCsv(header: true)
                                   .map { row -> tuple(row.sample_id, (row.sex ?: 'U')) }
                          : ch_v.map { sid, b, i -> tuple(sid, 'U') }
    ch_metrics = SAMTOOLS_STATS.out.flagstat.join(SAMTOOLS_STATS.out.stats)
        .join(MARKDUP.out.metrics).join(HSMETRICS.out.metrics)
        .join(MOSDEPTH.out.summary).join(MOSDEPTH.out.regions)
    SUMMARIZE_SAMPLE_QC( ch_sex.join(ch_metrics) )
    ch_qc_json = SUMMARIZE_SAMPLE_QC.out.json.collect()
    QC_GATE( ch_qc_json )
    PLOT_COVERAGE( ch_qc_json, MOSDEPTH.out.regiondist.map { it[1] }.collect() )
    PLOT_QC_GATE( QC_GATE.out.summary, ch_qc_json )

    def ch_png = PLOT_COVERAGE.out.png.mix(PLOT_QC_GATE.out.png).collect()
    def ch_cap = PLOT_COVERAGE.out.captions.mix(PLOT_QC_GATE.out.captions).collect()
    def ch_sum = QC_GATE.out.scorecard.mix(QC_GATE.out.summary).collect()
    STAGE_REPORT( Channel.value('coverage'), Channel.value('Coverage & QC gate (Stages 4-5)'),
                  Channel.value('stage5_qc_gate'), ch_png, ch_cap, ch_sum, Channel.value(synth) )
    STAGE_PROVENANCE( Channel.value('coverage'), Channel.value('stage5_qc_gate'),
        Channel.value("--input bams=${params.input_bams} --input reference=${params.genome_fasta} " +
                      "--input target_bed=${params.target_bed} " +
                      "--param min_mean_target_depth=${params.min_mean_target_depth} " +
                      "--param min_on_target=${params.min_on_target} --param max_dup_rate=${params.max_dup_rate} " +
                      "--tool picard --tool mosdepth --tool samtools") )
}

// ── call : QC-PASS BAM(s) -> joint hard-filtered VCF (+ calling/filter QC plots) ──
// cohort/joint. ASSUMES the supplied BAMs are QC-pass — it does NOT re-run the QC gate
// (run `callforge coverage` first and pass only the pass BAMs). Validates index+contig.
workflow CALL {
    main:
    need(params.input_bams,   "callforge call needs --input_bams <glob/comma list of QC-PASS BAM(s)>")
    need(params.genome_fasta, "callforge call needs --genome_fasta <reference>")
    need(params.target_bed,   "callforge call needs --target_bed <panel BED>")
    ch_bams   = bamChannel(params.input_bams)
    ch_genome = file(params.genome_fasta, checkIfExists: true)
    ch_bed    = file(params.target_bed,   checkIfExists: true)
    def synth = (workflow.profile?.contains('test')) ? '--synthetic' : ''

    PREPARE_REFERENCE( ch_genome )
    REFERENCE_INVARIANT( PREPARE_REFERENCE.out.fai, ch_bed )
    // pass-BAM semantics: assume QC-pass; validate index + contig agreement only.
    ch_v = VALIDATE_BAM( Channel.value('call'), ch_bams, PREPARE_REFERENCE.out.fai ).bam

    def ch_joint
    if (params.caller == 'deepvariant') {
        DEEPVARIANT( ch_v, PREPARE_REFERENCE.out.fasta, PREPARE_REFERENCE.out.fai, ch_bed )
        ch_dv = DEEPVARIANT.out.gvcf
        GLNEXUS( ch_dv.map { s, g, t -> g }.collect(), ch_dv.map { s, g, t -> t }.collect(), ch_bed )
        ch_joint = GLNEXUS.out.vcf
    } else {
        HAPLOTYPECALLER( ch_v, PREPARE_REFERENCE.out.fasta, PREPARE_REFERENCE.out.fai,
                         PREPARE_REFERENCE.out.dict, ch_bed )
        ch_g = HAPLOTYPECALLER.out.gvcf.map { s, g, t -> g }.collect()
        ch_t = HAPLOTYPECALLER.out.gvcf.map { s, g, t -> t }.collect()
        if (params.joint_method == 'combinegvcfs') {
            COMBINE_GENOTYPE( ch_g, ch_t, PREPARE_REFERENCE.out.fasta, PREPARE_REFERENCE.out.fai,
                              PREPARE_REFERENCE.out.dict, ch_bed )
            ch_joint = COMBINE_GENOTYPE.out.vcf
        } else {
            GENOMICSDB_IMPORT( ch_g, ch_t, ch_bed )
            GENOTYPE_GVCFS_DB( GENOMICSDB_IMPORT.out.gdb, PREPARE_REFERENCE.out.fasta,
                               PREPARE_REFERENCE.out.fai, PREPARE_REFERENCE.out.dict, ch_bed )
            ch_joint = GENOTYPE_GVCFS_DB.out.vcf
        }
    }
    CALLSET_STATS( ch_joint )
    HARD_FILTER( ch_joint, PREPARE_REFERENCE.out.fasta, PREPARE_REFERENCE.out.fai, PREPARE_REFERENCE.out.dict )
    FILTER_SUMMARY( HARD_FILTER.out.vcf )
    PLOT_CALLING_QC( CALLSET_STATS.out.json, FILTER_SUMMARY.out.qualdp )
    PLOT_FILTER_QC( FILTER_SUMMARY.out.counts )

    def ch_png = PLOT_CALLING_QC.out.png.mix(PLOT_FILTER_QC.out.png).collect()
    def ch_cap = PLOT_CALLING_QC.out.captions.mix(PLOT_FILTER_QC.out.captions).collect()
    def ch_sum = CALLSET_STATS.out.json.mix(FILTER_SUMMARY.out.counts).collect()
    STAGE_REPORT( Channel.value('call'), Channel.value('Joint calling + hard-filter (Stages 6-7)'),
                  Channel.value('stage7_filter'), ch_png, ch_cap, ch_sum, Channel.value(synth) )
    STAGE_PROVENANCE( Channel.value('call'), Channel.value('stage7_filter'),
        Channel.value("--input bams=${params.input_bams} --input reference=${params.genome_fasta} " +
                      "--input target_bed=${params.target_bed} --param caller=${params.caller} " +
                      "--param joint_method=${params.joint_method} --tool gatk --tool bcftools") )
}

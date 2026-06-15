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
include { VEP; ANNOTATE_DBS; VERIFY_ANNOTATION; FLATTEN_TSV;
          PLOT_ANNOTATION } from '../modules/stage11_annotate.nf'
include { BURDEN_MATRIX; BURDEN_COLLAPSE; BURDEN_REGENIE; BURDEN_SKAT;
          PLOT_BURDEN } from '../modules/stage14_burden.nf'
include { VALIDATE_STAGE_INPUTS; STAGE_REPORT; STAGE_PROVENANCE } from '../modules/stage_entry.nf'

// helper: require a param, else fail loud with a clear message
def need(cond, msg) { if (!cond) { exit 1, "ERROR: ${msg}" } }

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

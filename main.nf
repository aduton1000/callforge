#!/usr/bin/env nextflow
// ═══════════════════════════════════════════════════════════════════════════
//  CallForge — targeted-capture germline analysis pipeline
//  main.nf  (entrypoint → subworkflows/callforge.nf)
//  FASTQ → annotated joint callset → burden. Sibling of CaptureForge.
// ═══════════════════════════════════════════════════════════════════════════
nextflow.enable.dsl = 2

include { CALLFORGE } from './subworkflows/callforge.nf'

// Stage-subcommand entry workflows (run ONE stage standalone via `--stage <name>` /
// `callforge <stage>`). They reuse the same stage modules as the full pipeline.
include { ANNO; BURDEN; ALIGN; COVERAGE; CALL; CNV; STR; PARALOG;
          COHORTQC; GIAB } from './subworkflows/entries.nf'

def helpMessage() {
    log.info """
    ╔══════════════════════════════════════════════════════════════════════╗
    ║  CallForge ${workflow.manifest.version}  —  targeted-capture germline analysis  ║
    ╚══════════════════════════════════════════════════════════════════════╝
    Usage:
      nextflow run main.nf -profile <test|mac_local|hpc_slurm>[,conda|,docker|,apptainer] \\
        --input samplesheet.csv --genome_fasta noalt.fa \\
        --target_bed final_covered_targets.bed [--captureforge_dir <run>]

    Required: --input  --genome_fasta  --target_bed
    See params.example.yaml / nextflow.config for all parameters.
    """.stripIndent()
}

// Resolve an optional CaptureForge artifact (explicit path > under captureforge_dir > NO_FILE).
def resolveCF(explicit, candidates, placeholder) {
    if (explicit) {
        def f = file(explicit)
        if (f.exists()) return f
        log.warn "CaptureForge artifact not found at ${explicit}; continuing without it."
    }
    if (params.captureforge_dir) {
        def hit = candidates.collect { file("${params.captureforge_dir}/${it}") }.find { it.exists() }
        if (hit) return hit
    }
    return file("${projectDir}/assets/${placeholder}")
}

workflow {
    if (params.help) { helpMessage(); return }

    // ── Stage subcommands ────────────────────────────────────────────────────
    // `nextflow run main.nf --stage <name> --<inputs…>` runs ONE stage standalone on
    // user-supplied inputs, reusing the SAME modules as the pipeline. Exposed as
    // `callforge <stage>` by the CLI. (Param dispatch — portable across Nextflow
    // versions; the strict parser disallows `-entry` for named workflows.)
    if (params.stage) {
        // Non-destructive guard: a standalone stage must not silently overwrite the
        // pipeline's results/ (the re-annotate / re-burden use cases are COMPARISONS).
        def inplace = (params.in_place == true || "${params.in_place}".toLowerCase() == 'true')
        if (!inplace && params.outdir == 'results') {
            exit 1, "ERROR: a standalone stage run must not overwrite the pipeline's results/. " +
                    "Pass --outdir <dir> (the `callforge <stage>` CLI auto-creates " +
                    "results/standalone/<stage>_<timestamp>/), or --in_place true to update results/ in place."
        }
        if      (params.stage == 'anno')     { ANNO()     }
        else if (params.stage == 'burden')   { BURDEN()   }
        else if (params.stage == 'align')    { ALIGN()    }
        else if (params.stage == 'coverage') { COVERAGE() }
        else if (params.stage == 'call')     { CALL()     }
        else if (params.stage == 'cnv')      { CNV()      }
        else if (params.stage == 'str')      { STR()      }
        else if (params.stage == 'paralog')  { PARALOG()  }
        else if (params.stage == 'cohortqc') { COHORTQC() }
        else if (params.stage == 'giab')     { GIAB()     }
        else { exit 1, "ERROR: unknown stage '${params.stage}'. Available: align, coverage, call, cnv, str, paralog, anno, cohortqc, giab, burden" }
        return
    }

    if (!params.input)        { exit 1, "ERROR: --input sample sheet is required" }
    if (!params.genome_fasta) { exit 1, "ERROR: --genome_fasta (no-alt primary assembly) is required" }
    if (!params.target_bed)   { exit 1, "ERROR: --target_bed (final_covered_targets.bed) is required" }

    def cf_metrics = resolveCF(params.cf_metrics_json, ['stage9_qc/metrics.json', 'metrics.json'], 'NO_METRICS')
    def cf_baits   = resolveCF(params.cf_baits_csv,    ['order_package/baits.csv', 'baits.csv'], 'NO_BAITS')

    log.info """
    ─────────────────────────────────────────────────────────────────────────
    CallForge ${workflow.manifest.version}   profile=${workflow.profile}
    panel=${params.panel_name}   build=${params.genome_build}   species=${params.species}
    input=${params.input}
    genome=${params.genome_fasta}
    target_bed=${params.target_bed}
    captureforge_dir=${params.captureforge_dir}  (metrics=${cf_metrics.name}, baits=${cf_baits.name})
    caller=${params.caller}  cnv=${params.cnv_caller}  burden=${params.burden_engine}
    outdir=${params.outdir}
    ─────────────────────────────────────────────────────────────────────────
    """.stripIndent()

    ch_sheet      = file(params.input,        checkIfExists: true)
    ch_genome     = file(params.genome_fasta, checkIfExists: true)
    ch_target_bed = file(params.target_bed,   checkIfExists: true)

    // Reads channel built directly from the sheet (the validation process runs
    // in parallel as a fail-loud gate).
    ch_reads = Channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample_id || !row.fastq_1 || !row.fastq_2)
                exit 1, "ERROR: sample sheet row missing sample_id/fastq_1/fastq_2: ${row}"
            tuple(row.sample_id, file(row.fastq_1, checkIfExists: true),
                                 file(row.fastq_2, checkIfExists: true))
        }

    CALLFORGE( ch_reads, ch_sheet, ch_genome, ch_target_bed, cf_metrics, cf_baits )

    workflow.onComplete = {
        log.info ( workflow.success
            ? "\n[CallForge] Completed OK. Results: ${params.outdir}\n"
            : "\n[CallForge] FAILED (exit ${workflow.exitStatus}). See .nextflow.log\n" )
    }
}

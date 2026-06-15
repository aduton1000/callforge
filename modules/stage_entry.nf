// ── modules/stage_entry.nf  (stage-subcommand support: validate + report + provenance) ──
// THIN WRAPPERS for the standalone stage entry workflows (subworkflows/entries.nf).
// These do NOT contain any stage analysis logic — they only (a) fail-loud validate the
// externally-provided inputs, (b) render a focused per-stage report embedding that
// stage's existing QC plots/captions, and (c) write a per-stage provenance.json.
// The analytical work is done by the SAME stage modules the full pipeline uses.

process VALIDATE_STAGE_INPUTS {
    tag   "${stage}"
    label 'small'
    conda "${projectDir}/env/callforge.yml"

    input:
    val  stage
    tuple path(vcf), path(tbi)
    path fai            // reference .fai | assets/NO_CACHE (skip contig check)
    path sheet          // sample sheet  | assets/NO_CACHE (skip phenotype check)
    val  opts           // extra flags, e.g. '--require-phenotype'

    output:
    tuple path(vcf), path(tbi), emit: vcf   // pass-through gates the downstream module

    script:
    def faiarg = fai.name   != 'NO_CACHE' ? "--fai ${fai}"           : ""
    def shtarg = sheet.name != 'NO_CACHE' ? "--samplesheet ${sheet}" : ""
    """
    validate_stage_inputs.py --stage ${stage} --vcf ${vcf} ${faiarg} ${shtarg} ${opts}
    """
}

process STAGE_REPORT {
    tag   "${stage}"
    label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir { "${params.outdir}/${pubdir}" }, mode: params.publish_mode

    input:
    val  stage
    val  title
    val  pubdir
    path plots          // collected *.png from this stage's PLOT_* process
    path captions       // captions_*.tsv from this stage's PLOT_* process
    path summaries      // collected summary JSON/TSV for this stage
    val  opts           // e.g. '--synthetic'

    output:
    path "${stage}_report.md", emit: report

    script:
    def capl = (captions instanceof List ? captions : [captions]).collect { "${it}" }.join(' ')
    def suml = (summaries instanceof List ? summaries : [summaries]).collect { "--summary ${it}" }.join(' ')
    """
    stage_report.py --stage ${stage} --title "${title}" \\
        --captions ${capl} ${suml} --outdir . --out ${stage}_report.md ${opts}
    """
}

process STAGE_PROVENANCE {
    tag   "${stage}"
    label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir { "${params.outdir}/${pubdir}" }, mode: params.publish_mode

    input:
    val stage
    val pubdir
    val args            // '--input name=path ... --param k=v ... --tool bcftools ...'

    output:
    path "${stage}_provenance.json", emit: json

    script:
    // inputs are recorded/checksummed by absolute path (resolvable under -profile conda).
    """
    stage_provenance.py --stage ${stage} --out ${stage}_provenance.json --outdir . \\
        --version ${workflow.manifest.version} \\
        --git-commit \$(git -C ${projectDir} rev-parse --short HEAD 2>/dev/null || echo NA) \\
        --profile "${workflow.profile}" ${args}
    """
}

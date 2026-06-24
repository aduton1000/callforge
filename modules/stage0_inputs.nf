// ── modules/stage0_inputs.nf  (Stage 0: inputs & resource discovery) ─────────
// Validate sample sheet, prepare reference (.fai/.dict), enforce the reference
// invariant, auto-discover annotation resources, ingest CaptureForge metadata,
// count input reads, and render the Stage-0 QC plots.
// Every process pins the core conda env; the report flags emulated stages.

process VALIDATE_SAMPLESHEET {
    tag   "${params.panel_name}"
    label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage0_inputs", mode: params.publish_mode

    input:
    path sheet

    output:
    path "samplesheet.valid.csv",     emit: csv
    path "samplesheet_summary.json",  emit: summary

    script:
    """
    validate_samplesheet.py --sheet ${sheet} --check-files \\
        --out-csv samplesheet.valid.csv --out-json samplesheet_summary.json
    """
}

process PREPARE_REFERENCE {
    tag   "${params.genome_build}"
    label 'medium'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage0_inputs/reference", mode: params.publish_mode, enabled: false

    input:
    path fasta
    path fai_in      // pre-built ${fasta}.fai co-located with the FASTA, or assets/NO_FAI
    path dict_in     // pre-built ${fasta.baseName}.dict co-located with it, or assets/NO_DICT

    output:
    path "${fasta}",        emit: fasta
    path "${fasta}.fai",    emit: fai
    path "${fasta.baseName}.dict", emit: dict

    script:
    """
    # REUSE a co-located .fai / .dict if one was staged in (the workflow detects the
    # pre-built file next to the FASTA and stages it here under its canonical name, so
    # the guard below sees it). Build only when absent. Idempotent on the full genome.
    [ -s ${fasta}.fai ] || samtools faidx ${fasta}
    [ -s ${fasta.baseName}.dict ] || samtools dict ${fasta} -o ${fasta.baseName}.dict
    """
}

process REFERENCE_INVARIANT {
    tag   "${params.genome_build}"
    label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage0_inputs", mode: params.publish_mode

    input:
    path fai
    path bed

    output:
    path "reference_invariant.json", emit: report

    script:
    """
    # Fail-loud: no-alt assertion + BED/reference build match (exit 2 stops the run).
    check_reference_invariant.py --fai ${fai} --bed ${bed} \\
        --genome-build ${params.genome_build} --out-json reference_invariant.json
    """
}

process DISCOVER_RESOURCES {
    tag   "${params.species}"
    label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage0_inputs", mode: params.publish_mode

    input:
    val resource_dirs
    path bed

    output:
    path "resource_manifest.json", emit: manifest

    script:
    def vep      = params.vep_cache  ? "--vep-cache ${params.vep_cache}"   : ""
    def gnomad   = params.gnomad_vcf ? "--gnomad-vcf ${params.gnomad_vcf}" : ""
    def dbsnp    = params.dbsnp_vcf  ? "--dbsnp-vcf ${params.dbsnp_vcf}"    : ""
    def clinvar  = params.clinvar_vcf? "--clinvar-vcf ${params.clinvar_vcf}": ""
    def phylop   = params.phylop_bw  ? "--phylop-bw ${params.phylop_bw}"    : ""
    def known    = params.known_sites? "--known-sites ${params.known_sites}": ""
    def vrel     = params.vep_release? "--vep-release ${params.vep_release}": ""
    // Fail loud on a found-but-region-limited DB unless explicitly waived.
    def gate     = (params.ignore_unfit_resources || params.allow_download) ? "" : "--fail-on-unfit"
    """
    discover_resources.py --dirs "${resource_dirs}" \\
        --species ${params.species} --genome-build ${params.genome_build} ${vrel} \\
        --target-bed ${bed} --scope-min ${params.scope_min} ${gate} \\
        ${vep} ${gnomad} ${dbsnp} ${clinvar} ${phylop} ${known} \\
        --out-json resource_manifest.json
    """
}

process INGEST_CAPTUREFORGE {
    tag   "${params.panel_name}"
    label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage0_inputs", mode: params.publish_mode

    input:
    path bed
    path metrics      // assets/NO_FILE when absent
    path baits        // assets/NO_FILE when absent

    output:
    path "gene_metadata.tsv",  emit: tsv
    path "gene_metadata.json", emit: summary

    script:
    def m = metrics.name != 'NO_METRICS' ? "--metrics ${metrics}" : ""
    def b = baits.name   != 'NO_BAITS'   ? "--baits ${baits}"     : ""
    def p = params.paralog_genes ? "--paralog-genes '${params.paralog_genes}'" : ""
    """
    ingest_captureforge.py --bed ${bed} ${m} ${b} ${p} \\
        --out-tsv gene_metadata.tsv --out-json gene_metadata.json
    """
}

process COUNT_READS {
    tag   "${sample_id}"
    label 'small'
    conda "${projectDir}/env/callforge.yml"

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    path "${sample_id}.readcount.tsv", emit: tsv

    script:
    """
    n=\$(( \$(zcat -f ${r1} | wc -l) / 4 ))
    printf '%s\\t%s\\n' "${sample_id}" "\$n" > ${sample_id}.readcount.tsv
    """
}

process PLOT_STAGE0 {
    tag   "${params.panel_name}"
    label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage0_inputs/plots", mode: params.publish_mode

    input:
    path manifest
    path readcounts

    output:
    path "*.png", emit: png
    path "*.svg", emit: svg
    path "captions.tsv", emit: captions

    script:
    """
    plot_stage0.py --manifest ${manifest} --readcounts ${readcounts} --outdir .
    """
}

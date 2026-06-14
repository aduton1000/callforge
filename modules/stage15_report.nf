// ── modules/stage15_report.nf  (Stage 15: reporting) ─────────────────────────
// MultiQC over the per-tool outputs + a self-contained cohort QC dashboard (all
// captioned plots, base64-embedded) + provenance.json.

process MULTIQC {
    tag "${params.panel_name}"; label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage15_report", mode: params.publish_mode
    input:
    path qc_inputs
    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data",        optional: true, emit: data
    script:
    """
    multiqc . --filename multiqc_report.html --no-ansi -q || multiqc . --filename multiqc_report.html
    """
}

process DASHBOARD {
    tag "${params.panel_name}"; label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage15_report", mode: params.publish_mode
    input:
    path plots          // all *.png (collected)
    path captions       // all captions[_*].tsv (collected)
    path summaries      // key summary json/tsv (collected)
    output:
    path "cohort_qc_dashboard.html", emit: html
    script:
    """
    make_dashboard.py --panel ${params.panel_name} --out cohort_qc_dashboard.html
    """
}

process PROVENANCE {
    tag "${params.panel_name}"; label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage15_report", mode: params.publish_mode
    input:
    path manifest
    path sheet_summary
    path qc_summary
    path quarantine
    path burden_meta
    output:
    path "provenance.json", emit: json
    script:
    def bm = burden_meta.name != 'NO_CACHE' ? "--burden-meta ${burden_meta}" : ""
    """
    make_provenance.py --version '${workflow.manifest.version}' --profile '${workflow.profile}' \\
        --nextflow-version '${nextflow.version}' \\
        --genome-build '${params.genome_build}' --genome-fasta '${params.genome_fasta}' \\
        --manifest ${manifest} --sheet-summary ${sheet_summary} \\
        --qc-summary ${qc_summary} --quarantine ${quarantine} ${bm} \\
        --out provenance.json
    """
}

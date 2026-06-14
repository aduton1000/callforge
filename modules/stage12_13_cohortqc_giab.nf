// ── modules/stage12_13_cohortqc_giab.nf  (Stages 12-13) ──────────────────────
// Cohort QC (somalier relatedness/sex + ancestry PCA + missingness) and GIAB
// benchmarking (hap.py, restricted to the panel BED = on-target only). somalier
// uses off-target X/Y signal for sex — the autosomal-panel backstop. hap.py runs
// in the official container; scoring is BED-restricted so off-target truth never
// distorts the panel's measured sensitivity.

process SOMALIER_EXTRACT {
    tag "${sample_id}"; label 'small'
    conda "${projectDir}/env/cohortqc.yml"
    input:
    tuple val(sample_id), path(bam), path(bai)
    path fasta
    path fai
    path sites
    path sites_tbi
    output:
    path "*.somalier", emit: somalier
    script:
    """
    somalier extract --sites ${sites} -f ${fasta} -d . ${bam}
    """
}

process SOMALIER_RELATE {
    tag "${params.panel_name}"; label 'small'
    conda "${projectDir}/env/cohortqc.yml"
    publishDir "${params.outdir}/stage12_cohortqc/somalier", mode: params.publish_mode
    input:
    path somaliers
    output:
    path "cohort.samples.tsv", emit: samples
    path "cohort.pairs.tsv",   emit: pairs
    path "cohort.html",        optional: true, emit: html
    script:
    """
    somalier relate -o cohort ${somaliers}
    """
}

process COHORT_QC {
    tag "${params.panel_name}"; label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage12_cohortqc", mode: params.publish_mode
    input:
    tuple path(vcf), path(tbi)
    path somalier_samples
    path somalier_pairs
    path samplesheet
    output:
    path "cohort_qc.json", emit: json
    path "*.png",          emit: png
    path "*.svg",          emit: svg
    path "captions_cohortqc.tsv", emit: captions
    script:
    """
    export CALLFORGE_CAPTIONS=captions_cohortqc.tsv
    cohort_qc.py --vcf ${vcf} --somalier-samples ${somalier_samples} \\
        --somalier-pairs ${somalier_pairs} --samplesheet ${samplesheet} --outdir .
    """
}

process EXTRACT_CONTROL {
    tag "${params.giab_control_id}"; label 'small'
    conda "${projectDir}/env/callforge.yml"
    input:
    tuple path(vcf), path(tbi)
    output:
    tuple path("control.vcf.gz"), path("control.vcf.gz.tbi"), emit: vcf
    script:
    """
    bcftools view -s ${params.giab_control_id} -f PASS,. ${vcf} -Oz -o control.vcf.gz
    tabix -p vcf control.vcf.gz
    """
}

process GIAB_HAPPY {
    tag "${params.giab_control_id}"; label 'medium'
    container params.happy_image
    publishDir "${params.outdir}/stage13_giab", mode: params.publish_mode
    input:
    tuple path(query), path(qtbi)
    path truth_vcf
    path truth_tbi
    path truth_bed
    path panel_bed
    path fasta
    path fai
    output:
    path "happy.summary.csv", emit: summary
    path "happy.*",           emit: all
    script:
    // -T restricts to the panel BED (ON-TARGET only); -f is the truth confident BED.
    """
    /opt/hap.py/bin/hap.py ${truth_vcf} ${query} \\
        -f ${truth_bed} -T ${panel_bed} -r ${fasta} -o happy \\
        --no-roc --threads ${task.cpus} || hap.py ${truth_vcf} ${query} \\
        -f ${truth_bed} -T ${panel_bed} -r ${fasta} -o happy --no-roc --threads ${task.cpus}
    """
}

process PLOT_GIAB {
    tag "${params.panel_name}"; label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage13_giab/plots", mode: params.publish_mode
    input:
    path summary
    output:
    path "*.png", emit: png
    path "*.svg", emit: svg
    path "giab_metrics.json", emit: metrics
    path "captions_giab.tsv", emit: captions
    script:
    """
    export CALLFORGE_CAPTIONS=captions_giab.tsv
    plot_giab.py --summary ${summary} --restricted-to "panel BED (on-target)" --outdir .
    """
}

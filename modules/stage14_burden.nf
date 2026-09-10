// ── modules/stage14_burden.nf  (Stage 14: rare-variant burden / association) ──
// Build the rare+functional burden matrix collapsed by gene and CaptureForge
// gene-set (burden_group), with phenotype + covariates (covariate_* + ancestry
// PCs), then a gene-burden test. Engines: collapse (chi-square, dependency-free,
// small-N/default) | regenie (Firth, production) | skat (SKAT-O/STAAR, R, alt).
// Gated on phenotype: build_burden_matrix writes a skip marker if none.

process BURDEN_MATRIX {
    tag "${params.panel_name}"; label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage14_burden", mode: params.publish_mode
    input:
    tuple path(vcf), path(tbi)
    path metadata
    path samplesheet
    path cohort_qc        // cohort_qc.json | assets/NO_CACHE
    output:
    path "burden_matrix.tsv",        emit: matrix
    path "phenotype.tsv",            emit: pheno
    path "covariates.tsv",           emit: covar
    path "burden_meta.json",         emit: meta
    path "qualifying_variants.tsv",  emit: qual
    path "genotypes.tsv",            emit: genotypes
    script:
    def cq = cohort_qc.name != 'NO_CACHE' ? "--cohort-qc ${cohort_qc}" : ""
    def cov = params.covariates ? "--covariates ${params.covariates}" : ""
    """
    build_burden_matrix.py --vcf ${vcf} --metadata ${metadata} --samplesheet ${samplesheet} \\
        ${cq} --af-field ${params.gnomad_af_field} --af-max ${params.burden_af_max} \\
        --csq '${params.burden_csq}' --n-pcs ${params.n_ancestry_pcs} ${cov} --outdir .
    """
}

process BURDEN_COLLAPSE {
    tag "${params.panel_name}"; label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage14_burden", mode: params.publish_mode
    input:
    path matrix
    path pheno
    output:
    path "burden_results.tsv",  emit: results
    path "burden_summary.json", emit: summary
    script:
    """
    burden_collapse_test.py --matrix ${matrix} --phenotype ${pheno} --outdir .
    """
}

process BURDEN_REGENIE {
    // Production engine (Firth-corrected, covariate-adjusted). Needs a real cohort.
    tag "${params.panel_name}"; label 'medium'
    // image: --burden_image, else <stage_images_dir>/callforge-burden.sif, else the core image (tool absent!)
    container { params.burden_image ?: (params.stage_images_dir ? "${params.stage_images_dir}/callforge-burden.sif" : params.container_image) }
    conda "${projectDir}/env/burden.yml"
    publishDir "${params.outdir}/stage14_burden", mode: params.publish_mode
    input:
    path genotypes
    path pheno
    path covar
    output:
    path "burden_results.tsv",  emit: results
    path "burden_summary.json", emit: summary
    script:
    """
    run_regenie_burden.sh ${genotypes} ${pheno} ${covar} burden_results.tsv burden_summary.json
    """
}

process BURDEN_SKAT {
    // Alternative engine: SKAT-O in R (STAAR packaged in the same env/image, not yet wired). Needs a real cohort.
    tag "${params.panel_name}"; label 'medium'
    // image: --skat_image, else <stage_images_dir>/callforge-skat.sif, else the core image (tool absent!)
    container { params.skat_image ?: (params.stage_images_dir ? "${params.stage_images_dir}/callforge-skat.sif" : params.container_image) }
    conda "${projectDir}/env/skat.yml"
    publishDir "${params.outdir}/stage14_burden", mode: params.publish_mode
    input:
    path genotypes
    path pheno
    path covar
    output:
    path "burden_results.tsv",  emit: results
    path "burden_summary.json", emit: summary
    script:
    """
    run_skat_burden.R ${genotypes} ${pheno} ${covar} burden_results.tsv burden_summary.json
    """
}

process PLOT_BURDEN {
    tag "${params.panel_name}"; label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage14_burden/plots", mode: params.publish_mode
    input:
    path results
    path summary
    output:
    path "*.png", emit: png
    path "*.svg", emit: svg
    path "captions_burden.tsv", emit: captions
    script:
    """
    export CALLFORGE_CAPTIONS=captions_burden.tsv
    plot_burden.py --results ${results} --summary ${summary} --outdir .
    """
}

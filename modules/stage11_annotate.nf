// ── modules/stage11_annotate.nf  (Stage 11: annotation) ──────────────────────
// VEP (consequence/transcript/SIFT-PolyPhen/canonical/MANE) via the official
// ensembl-vep container, then manifest-driven vcfanno (gnomAD-AFR AF, dbSNP IDs,
// ClinVar) + PhyloP with chr-reconciliation, then annotation-landing verification.
// Organism-generic: VEP species/cache + vcfanno BYO databases; absent human
// resources skip gracefully. Adapted from the nf-core/raredisease VEP+vcfanno pattern.

process VEP {
    tag "${params.panel_name}"; label 'medium'
    container params.vep_image
    conda "${projectDir}/env/annotate.yml"   // local (Apple-Silicon) conda fallback: under
                                             // -profile local,conda VEP runs from this env
                                             // (ensembl-vep matched to the local cache release);
                                             // the container is used under a container engine.
    publishDir "${params.outdir}/stage11_annotation/vep", mode: params.publish_mode
    input:
    tuple path(vcf), path(tbi)
    path fasta
    path fai
    path aux          // gtf.gz (gtf mode) | VEP cache dir (cache mode)
    path aux_idx      // gtf .tbi (gtf mode) | assets/NO_FILE (cache mode)
    output:
    tuple path("vep.vcf.gz"), path("vep.vcf.gz.tbi"), emit: vcf
    script:
    def common = "--vcf --format vcf --force_overwrite --no_stats --fasta ${fasta} --symbol --canonical"
    // gtf mode uses the gene-model GTF directly (no cache, so NO --offline);
    // cache mode is fully offline against the discovered VEP cache.
    def modeargs = params.vep_mode == 'gtf' ?
        "--gtf ${aux}" :
        "--offline --cache --dir_cache ${aux} --species ${params.species} --assembly ${params.genome_build} --sift b --polyphen b --mane --hgvs --af_gnomade --af_gnomadg --max_af"
    """
    vep ${common} ${modeargs} --fork ${task.cpus} -i ${vcf} -o vep.vcf.gz --compress_output bgzip
    tabix -p vcf vep.vcf.gz
    """
}

process ANNOTATE_DBS {
    tag "${params.panel_name}"; label 'medium'
    // image: --annotate_image, else <stage_images_dir>/callforge-annotate.sif, else the core image (tool absent!)
    container { params.annotate_image ?: (params.stage_images_dir ? "${params.stage_images_dir}/callforge-annotate.sif" : params.container_image) }
    conda "${projectDir}/env/annotate.yml"
    publishDir "${params.outdir}/stage11_annotation", mode: params.publish_mode
    input:
    tuple path(vcf), path(tbi)
    path fai
    path manifest
    output:
    tuple path("annotated.vcf.gz"), path("annotated.vcf.gz.tbi"), emit: vcf
    path "annotation_sources.json", emit: sources
    script:
    """
    annotate_dbs.py --vcf ${vcf} --fai ${fai} --manifest ${manifest} \\
        --gnomad-af-field ${params.gnomad_af_field} --outdir . --out annotated.vcf.gz
    """
}

process VERIFY_ANNOTATION {
    tag "${params.panel_name}"; label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage11_annotation", mode: params.publish_mode
    input:
    tuple path(vcf), path(tbi)
    path sources
    output:
    path "annotation_landing.json",      emit: landing
    path "annotation_landing_table.tsv", emit: table
    script:
    """
    verify_annotation.py --vcf ${vcf} --sources ${sources} \\
        --gnomad-af-field ${params.gnomad_af_field} --outdir .
    """
}

process FLATTEN_TSV {
    tag "${params.panel_name}"; label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage11_annotation", mode: params.publish_mode
    input:
    tuple path(vcf), path(tbi)
    output:
    path "variants.flat.tsv", emit: tsv
    script:
    """
    flatten_vcf.py --vcf ${vcf} --gnomad-af-field ${params.gnomad_af_field} --out variants.flat.tsv
    """
}

process PLOT_ANNOTATION {
    tag "${params.panel_name}"; label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage11_annotation/plots", mode: params.publish_mode
    input:
    path tsv
    path landing
    output:
    path "*.png", emit: png
    path "*.svg", emit: svg
    path "captions_annotation.tsv", emit: captions
    script:
    """
    export CALLFORGE_CAPTIONS=captions_annotation.tsv
    plot_annotation.py --tsv ${tsv} --landing ${landing} --af-field ${params.gnomad_af_field} --outdir .
    """
}

// ── modules/stage8_10_cnv_str_paralog.nf  (Stages 8-10) ──────────────────────
// CNV (CNVkit, default) + STR (ExpansionHunter) + paralog-aware flagging, all
// driven by CaptureForge per-gene metadata. Confidence labels are PROPAGATED to
// outputs: CNV calls carry their callability (CR1/CFH = breakpoint-blind),
// paralog-region variants carry PARALOG_GENE + PARALOG_CONF.

process CNVKIT_BATCH {
    tag "${params.panel_name}"; label 'large'
    // image: --cnv_image, else <stage_images_dir>/callforge-cnv.sif, else the core image (tool absent!)
    container { params.cnv_image ?: (params.stage_images_dir ? "${params.stage_images_dir}/callforge-cnv.sif" : params.container_image) }
    conda "${projectDir}/env/cnv.yml"
    publishDir "${params.outdir}/stage8_cnv/cnvkit", mode: params.publish_mode
    input:
    path bams
    path bais
    path bed
    path fasta
    path fai
    output:
    path "cnvout/*.cnr",      emit: cnr
    path "cnvout/*.call.cns", emit: cns
    script:
    """
    # Germline cohort, two-step pooled panel-of-normals:
    #  1) build a PoN reference from ALL QC-pass samples
    cnvkit.py batch --normal *.bam --method ${params.cnv_method} \\
        --targets ${bed} --fasta ${fasta} \\
        --output-reference pon.cnn --output-dir ref_build -p ${task.cpus}
    #  2) compute copy ratios + segments + integer calls for each sample vs the PoN
    #     (segment-method haar is pure-python; cbs needs R/DNAcopy — see env/cnv.yml)
    cnvkit.py batch *.bam --reference pon.cnn --method ${params.cnv_method} \\
        --segment-method ${params.cnv_segment_method} --output-dir cnvout -p ${task.cpus}
    # ensure integer-copy-number .call.cns exists for each sample (nullglob-guarded)
    shopt -s nullglob
    for cns in cnvout/*.cns; do
        case "\$cns" in *.call.cns|*.bintest.cns) continue;; esac
        base=\$(basename "\$cns" .cns)
        [ -s "cnvout/\${base}.call.cns" ] || cnvkit.py call "\$cns" -o "cnvout/\${base}.call.cns"
    done
    """
}

process CNV_ANNOTATE {
    tag "${params.panel_name}"; label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage8_cnv", mode: params.publish_mode
    input:
    path cns
    path cnr
    path bed
    path metadata
    output:
    path "cnv_calls.tsv",     emit: calls
    path "cnv_copyratio.tsv", emit: copyratio
    path "cnv_summary.json",  emit: summary
    script:
    """
    cnv_annotate.py --cns ${cns} --cnr ${cnr} --bed ${bed} --metadata ${metadata} --outdir .
    """
}

process PLOT_CNV {
    tag "${params.panel_name}"; label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage8_cnv/plots", mode: params.publish_mode
    input:
    path copyratio
    path summary
    path metadata
    output:
    path "*.png", emit: png
    path "*.svg", emit: svg
    path "captions_cnv.tsv", emit: captions
    script:
    """
    export CALLFORGE_CAPTIONS=captions_cnv.tsv
    plot_cnv.py --copyratio ${copyratio} --summary ${summary} --metadata ${metadata} --outdir .
    """
}

process BUILD_STR_CATALOG {
    tag "${params.panel_name}"; label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage9_str", mode: params.publish_mode
    input:  path metadata
    output: path "str_catalog.json", emit: catalog
    script:
    """
    build_str_catalog.py --metadata ${metadata} --motif ${params.str_motif} --out str_catalog.json
    """
}

process EXPANSIONHUNTER {
    tag "${sample_id}"; label 'medium'
    // image: --str_image, else <stage_images_dir>/callforge-str.sif, else the core image (tool absent!)
    container { params.str_image ?: (params.stage_images_dir ? "${params.stage_images_dir}/callforge-str.sif" : params.container_image) }
    conda "${projectDir}/env/str.yml"
    publishDir "${params.outdir}/stage9_str/per_sample", mode: params.publish_mode
    input:
    tuple val(sample_id), path(bam), path(bai)
    path fasta
    path fai
    path catalog
    output:
    path "${sample_id}.json", emit: json
    path "${sample_id}.vcf",  emit: vcf
    script:
    """
    ExpansionHunter --reads ${bam} --reference ${fasta} \\
        --variant-catalog ${catalog} --output-prefix ${sample_id} --threads ${task.cpus}
    """
}

process STR_SUMMARIZE {
    tag "${params.panel_name}"; label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage9_str", mode: params.publish_mode
    input:  path eh_json
    output:
    path "str_calls.tsv",   emit: calls
    path "str_summary.json", emit: summary
    script:
    """
    str_summarize.py --eh-json ${eh_json} --outdir .
    """
}

process PLOT_STR {
    tag "${params.panel_name}"; label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage9_str/plots", mode: params.publish_mode
    input:
    path calls
    path summary
    output:
    path "*.png", emit: png
    path "*.svg", emit: svg
    path "captions_str.tsv", emit: captions
    script:
    """
    export CALLFORGE_CAPTIONS=captions_str.tsv
    plot_str.py --calls ${calls} --summary ${summary} --outdir .
    """
}

process PARALOG_FLAG {
    tag "${params.panel_name}"; label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage10_paralog", mode: params.publish_mode
    input:
    tuple path(vcf), path(tbi)
    path bed
    path metadata
    output:
    tuple path("paralog.annotated.vcf.gz"), path("paralog.annotated.vcf.gz.tbi"), emit: vcf
    path "paralog_variants.tsv", emit: variants
    path "paralog_summary.json", emit: summary
    script:
    """
    paralog_flag.py --vcf ${vcf} --bed ${bed} --metadata ${metadata} \\
        --min-mq ${params.paralog_min_mq} --outdir .
    """
}

process PLOT_PARALOG {
    tag "${params.panel_name}"; label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage10_paralog/plots", mode: params.publish_mode
    input:  path summary
    output:
    path "*.png", emit: png
    path "*.svg", emit: svg
    path "captions_paralog.tsv", emit: captions
    script:
    """
    export CALLFORGE_CAPTIONS=captions_paralog.tsv
    plot_paralog.py --summary ${summary} --outdir .
    """
}

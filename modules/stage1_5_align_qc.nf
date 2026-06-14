// ── modules/stage1_5_align_qc.nf  (Stages 1-5) ───────────────────────────────
// Raw QC -> trimming -> alignment -> dedup/BQSR -> HsMetrics/coverage ->
// per-sample QC summary -> cohort QC gate (flag & quarantine). Adapted from the
// validated nf-core/sarek alignment+dedup+BQSR pattern; per-gene depth tracks and
// the flag/quarantine scorecard are bespoke. Every process pins the core env.

process FASTQC_RAW {
    tag "${sample_id}"; label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage1_rawqc", mode: params.publish_mode
    input:  tuple val(sample_id), path(r1), path(r2)
    output: path "*_fastqc.zip", emit: zip
    script:
    """
    fastqc -q -t ${task.cpus} ${r1} ${r2}
    """
}

process FASTP {
    tag "${sample_id}"; label 'medium'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage2_trim", mode: params.publish_mode, pattern: "*.json"
    input:  tuple val(sample_id), path(r1), path(r2)
    output:
    tuple val(sample_id), path("${sample_id}.trim.R1.fastq.gz"), path("${sample_id}.trim.R2.fastq.gz"), emit: reads
    path "${sample_id}.fastp.json", emit: json
    script:
    def a1 = params.adapter_r1 ? "--adapter_sequence ${params.adapter_r1}" : ""
    def a2 = params.adapter_r2 ? "--adapter_sequence_r2 ${params.adapter_r2}" : ""
    """
    fastp -i ${r1} -I ${r2} -o ${sample_id}.trim.R1.fastq.gz -O ${sample_id}.trim.R2.fastq.gz \\
        ${a1} ${a2} --qualified_quality_phred ${params.fastp_qual} --length_required ${params.fastp_min_len} \\
        --detect_adapter_for_pe --thread ${task.cpus} \\
        --json ${sample_id}.fastp.json --html ${sample_id}.fastp.html 2> ${sample_id}.fastp.log
    """
}

process BWAMEM2_INDEX {
    tag "${fasta}"; label 'index'
    conda "${projectDir}/env/callforge.yml"
    input:  path fasta
    output: path "${fasta}.*", emit: idx
    script:
    """
    bwa-mem2 index ${fasta}
    """
}

process BWAMEM2_ALIGN {
    tag "${sample_id}"; label 'align'
    conda "${projectDir}/env/callforge.yml"
    input:
    tuple val(sample_id), path(r1), path(r2)
    path fasta
    path idx
    output: tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: bam
    script:
    def rg = "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA\\tLB:${sample_id}"
    """
    bwa-mem2 mem -t ${task.cpus} -R '${rg}' ${fasta} ${r1} ${r2} \\
        | samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam -
    samtools index ${sample_id}.sorted.bam
    """
}

process MARKDUP {
    tag "${sample_id}"; label 'medium'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage4_postalign/markdup", mode: params.publish_mode, pattern: "*.metrics.txt"
    input:  tuple val(sample_id), path(bam), path(bai)
    output:
    tuple val(sample_id), path("${sample_id}.md.bam"), path("${sample_id}.md.bam.bai"), emit: bam
    tuple val(sample_id), path("${sample_id}.markdup.metrics.txt"), emit: metrics
    script:
    """
    picard MarkDuplicates I=${bam} O=${sample_id}.md.bam M=${sample_id}.markdup.metrics.txt \\
        CREATE_INDEX=false VALIDATION_STRINGENCY=LENIENT
    samtools index ${sample_id}.md.bam
    """
}

process BQSR {
    tag "${sample_id}"; label 'medium'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage4_postalign/bqsr", mode: params.publish_mode, pattern: "*.recal.table"
    input:
    tuple val(sample_id), path(bam), path(bai)
    path fasta
    path fai
    path dict
    val known_sites      // comma-separated absolute paths, or '' to skip
    output: tuple val(sample_id), path("${sample_id}.analysis.bam"), path("${sample_id}.analysis.bam.bai"), emit: bam
    script:
    if (known_sites?.trim()) {
        def ks = known_sites.split(',').collect { "--known-sites ${it}" }.join(' ')
        """
        gatk BaseRecalibrator -I ${bam} -R ${fasta} ${ks} -O ${sample_id}.recal.table
        gatk ApplyBQSR -I ${bam} -R ${fasta} --bqsr-recal-file ${sample_id}.recal.table -O ${sample_id}.analysis.bam
        samtools index ${sample_id}.analysis.bam
        """
    } else {
        // Graceful skip: no known-sites (non-human / unavailable). Pass dedup BAM through.
        """
        echo "BQSR skipped: no known-sites supplied for ${sample_id}" > ${sample_id}.bqsr.skipped
        ln -s ${bam} ${sample_id}.analysis.bam
        samtools index -f ${sample_id}.analysis.bam || cp ${bai} ${sample_id}.analysis.bam.bai
        """
    }
}

process BED_TO_INTERVALS {
    tag "${params.panel_name}"; label 'small'
    conda "${projectDir}/env/callforge.yml"
    input:  path bed
            path dict
    output: path "targets.interval_list", emit: il
    script:
    """
    picard BedToIntervalList I=${bed} O=targets.interval_list SD=${dict}
    """
}

process HSMETRICS {
    tag "${sample_id}"; label 'medium'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage4_postalign/hsmetrics", mode: params.publish_mode
    input:
    tuple val(sample_id), path(bam), path(bai)
    path fasta
    path fai
    path il
    output: tuple val(sample_id), path("${sample_id}.hs_metrics.txt"), emit: metrics
    script:
    """
    picard CollectHsMetrics I=${bam} O=${sample_id}.hs_metrics.txt R=${fasta} \\
        BAIT_INTERVALS=${il} TARGET_INTERVALS=${il} VALIDATION_STRINGENCY=LENIENT
    """
}

process SAMTOOLS_STATS {
    tag "${sample_id}"; label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage3_align", mode: params.publish_mode
    input:  tuple val(sample_id), path(bam), path(bai)
    output:
    tuple val(sample_id), path("${sample_id}.flagstat"), emit: flagstat
    tuple val(sample_id), path("${sample_id}.stats"),    emit: stats
    script:
    """
    samtools flagstat ${bam} > ${sample_id}.flagstat
    samtools stats ${bam} > ${sample_id}.stats
    """
}

process MOSDEPTH {
    tag "${sample_id}"; label 'medium'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage4_postalign/mosdepth", mode: params.publish_mode
    input:
    tuple val(sample_id), path(bam), path(bai)
    path bed
    output:
    tuple val(sample_id), path("${sample_id}.mosdepth.summary.txt"), emit: summary
    tuple val(sample_id), path("${sample_id}.regions.bed.gz"),       emit: regions
    tuple val(sample_id), path("${sample_id}.mosdepth.region.dist.txt"), emit: regiondist
    script:
    """
    mosdepth -t ${task.cpus} --by ${bed} --no-per-base --mapq 20 ${sample_id} ${bam}
    """
}

process SUMMARIZE_SAMPLE_QC {
    tag "${sample_id}"; label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage5_qc_gate/per_sample", mode: params.publish_mode
    input:
    tuple val(sample_id), val(sex), path(flagstat), path(stats), path(markdup),
          path(hsmetrics), path(mdsum), path(mdreg)
    output:
    path "${sample_id}.qc.json", emit: json
    path "${sample_id}.qc.ishist.tsv", optional: true, emit: ishist
    script:
    """
    summarize_sample_qc.py --sample-id ${sample_id} --sex ${sex} \\
        --flagstat ${flagstat} --stats ${stats} --markdup ${markdup} \\
        --hsmetrics ${hsmetrics} --mosdepth-summary ${mdsum} --mosdepth-regions ${mdreg} \\
        --out-json ${sample_id}.qc.json
    """
}

process QC_GATE {
    tag "${params.panel_name}"; label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage5_qc_gate", mode: params.publish_mode
    input:  path qc_json
    output:
    path "qc_scorecard.tsv",        emit: scorecard
    path "qc_pass.txt",             emit: pass
    path "qc_quarantine.txt",       emit: quarantine
    path "cohort_qc_summary.json",  emit: summary
    script:
    def sex = params.enforce_sex_check ? "--enforce-sex-check" : ""
    """
    qc_gate.py --qc-json ${qc_json} \\
        --min-mean-target-depth ${params.min_mean_target_depth} \\
        --max-dup-rate ${params.max_dup_rate} --min-on-target ${params.min_on_target} \\
        --max-contamination ${params.max_contamination} ${sex} --outdir .
    """
}

process PLOT_ALIGN_QC {
    tag "${params.panel_name}"; label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage5_qc_gate/plots", mode: params.publish_mode
    input:  path qc_json
            path ishist
    output: path "*.png", emit: png
            path "*.svg", emit: svg
            path "captions_align.tsv", emit: captions
    script:
    def ish = ishist ? "--ishist ${ishist}" : ""
    """
    export CALLFORGE_CAPTIONS=captions_align.tsv
    plot_align_qc.py --qc-json ${qc_json} ${ish} --outdir .
    """
}

process PLOT_COVERAGE {
    tag "${params.panel_name}"; label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage5_qc_gate/plots", mode: params.publish_mode
    input:  path qc_json
            path region_dist
    output: path "*.png", emit: png
            path "*.svg", emit: svg
            path "captions_coverage.tsv", emit: captions
    script:
    """
    export CALLFORGE_CAPTIONS=captions_coverage.tsv
    plot_coverage.py --qc-json ${qc_json} --region-dist ${region_dist} --outdir .
    """
}

process PLOT_QC_GATE {
    tag "${params.panel_name}"; label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage5_qc_gate/plots", mode: params.publish_mode
    input:  path summary
            path qc_json
    output: path "*.png", emit: png
            path "*.svg", emit: svg
            path "captions_qc.tsv", emit: captions
    script:
    """
    export CALLFORGE_CAPTIONS=captions_qc.tsv
    plot_qc_gate.py --summary ${summary} --qc-json ${qc_json} --outdir .
    """
}

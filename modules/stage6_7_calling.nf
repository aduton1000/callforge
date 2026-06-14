// ── modules/stage6_7_calling.nf  (Stages 6-7) ────────────────────────────────
// SNV/indel joint calling + hard-filtering. Adapted from the validated
// nf-core/sarek germline path: HaplotypeCaller(GVCF) -> GenomicsDBImport /
// CombineGVCFs -> GenotypeGVCFs -> GATK hard-filters. DeepVariant + GLnexus are
// the documented container-only alternative (--caller deepvariant).
// IMPORTANT: only QC-passing BAMs are fed here; quarantined samples never enter
// the GenomicsDB / joint step (the subworkflow joins on qc_pass.txt).

process HAPLOTYPECALLER {
    tag "${sample_id}"; label 'large'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage6_calling/gvcf", mode: params.publish_mode
    input:
    tuple val(sample_id), path(bam), path(bai)
    path fasta
    path fai
    path dict
    path bed
    output:
    tuple val(sample_id), path("${sample_id}.g.vcf.gz"), path("${sample_id}.g.vcf.gz.tbi"), emit: gvcf
    script:
    """
    gatk HaplotypeCaller -R ${fasta} -I ${bam} -O ${sample_id}.g.vcf.gz \\
        -ERC GVCF -L ${bed} --interval-padding 100 -ploidy ${params.ploidy}
    """
}

process GENOMICSDB_IMPORT {
    tag "${params.panel_name}"; label 'large'
    conda "${projectDir}/env/callforge.yml"
    input:
    path gvcfs
    path tbis
    path bed
    output:
    path "genomicsdb", emit: gdb
    script:
    """
    # sample-name map from the gVCFs (sample name == RG SM == sample_id)
    for f in *.g.vcf.gz; do printf '%s\\t%s\\n' "\$(bcftools query -l \$f)" "\$f"; done > map.txt
    echo "[GenomicsDBImport] importing \$(wc -l < map.txt) samples"; cat map.txt
    gatk GenomicsDBImport --sample-name-map map.txt -L ${bed} --merge-input-intervals \\
        --genomicsdb-workspace-path genomicsdb --batch-size 50 --reader-threads ${task.cpus}
    """
}

process GENOTYPE_GVCFS_DB {
    tag "${params.panel_name}"; label 'large'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage6_calling", mode: params.publish_mode
    input:
    path gdb
    path fasta
    path fai
    path dict
    path bed
    output:
    tuple path("joint.vcf.gz"), path("joint.vcf.gz.tbi"), emit: vcf
    script:
    """
    gatk GenotypeGVCFs -R ${fasta} -V gendb://${gdb} -O joint.vcf.gz -L ${bed}
    """
}

process COMBINE_GENOTYPE {
    // Alternative joint method (--joint_method combinegvcfs): CombineGVCFs -> GenotypeGVCFs.
    tag "${params.panel_name}"; label 'large'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage6_calling", mode: params.publish_mode
    input:
    path gvcfs
    path tbis
    path fasta
    path fai
    path dict
    path bed
    output:
    tuple path("joint.vcf.gz"), path("joint.vcf.gz.tbi"), emit: vcf
    script:
    """
    args=\$(for f in *.g.vcf.gz; do echo -n " -V \$f"; done)
    gatk CombineGVCFs -R ${fasta} \$args -O combined.g.vcf.gz
    gatk GenotypeGVCFs -R ${fasta} -V combined.g.vcf.gz -O joint.vcf.gz -L ${bed}
    """
}

process HARD_FILTER {
    tag "${params.panel_name}"; label 'medium'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage7_filter", mode: params.publish_mode
    input:
    tuple path(vcf), path(tbi)
    path fasta
    path fai
    path dict
    output:
    tuple path("joint.filtered.vcf.gz"), path("joint.filtered.vcf.gz.tbi"), emit: vcf
    script:
    """
    # SNPs
    gatk SelectVariants -R ${fasta} -V ${vcf} --select-type-to-include SNP -O snps.vcf.gz
    gatk VariantFiltration -R ${fasta} -V snps.vcf.gz \\
        --filter-expression "${params.snp_filter_expr}" --filter-name CALLFORGE_SNP \\
        -O snps.filt.vcf.gz
    # indels + mixed/MNP
    gatk SelectVariants -R ${fasta} -V ${vcf} \\
        --select-type-to-include INDEL --select-type-to-include MIXED --select-type-to-include MNP \\
        -O indels.vcf.gz
    gatk VariantFiltration -R ${fasta} -V indels.vcf.gz \\
        --filter-expression "${params.indel_filter_expr}" --filter-name CALLFORGE_INDEL \\
        -O indels.filt.vcf.gz
    # recombine + sort
    gatk MergeVcfs -I snps.filt.vcf.gz -I indels.filt.vcf.gz -O joint.filtered.vcf.gz
    """
}

process CALLSET_STATS {
    tag "${params.panel_name}"; label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage6_calling", mode: params.publish_mode
    input:
    tuple path(vcf), path(tbi)
    output:
    path "joint.bcftools_stats.txt", emit: stats
    path "calling_stats.json",       emit: json
    script:
    """
    bcftools stats -s - ${vcf} > joint.bcftools_stats.txt
    parse_bcftools_stats.py --stats joint.bcftools_stats.txt --label joint_prefilter \\
        --out-json calling_stats.json
    """
}

process FILTER_SUMMARY {
    tag "${params.panel_name}"; label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage7_filter", mode: params.publish_mode
    input:
    tuple path(vcf), path(tbi)
    output:
    path "filter_counts.json", emit: counts
    path "qual_dp.tsv",        emit: qualdp
    script:
    """
    filter_counts.py --vcf ${vcf} \\
        --snp-expr "${params.snp_filter_expr}" --indel-expr "${params.indel_filter_expr}" \\
        --out-json filter_counts.json
    echo -e "#QUAL\\tDP\\tTYPE\\tFILTER" > qual_dp.tsv
    bcftools query -f '%QUAL\\t%INFO/DP\\t%TYPE\\t%FILTER\\n' ${vcf} >> qual_dp.tsv
    """
}

process PLOT_CALLING_QC {
    tag "${params.panel_name}"; label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage6_calling/plots", mode: params.publish_mode
    input:
    path stats_json
    path qual_dp
    output:
    path "*.png", emit: png
    path "*.svg", emit: svg
    path "captions_calling.tsv", emit: captions
    script:
    """
    export CALLFORGE_CAPTIONS=captions_calling.tsv
    plot_calling_qc.py --stats-json ${stats_json} --qual-dp-tsv ${qual_dp} --outdir .
    """
}

process PLOT_FILTER_QC {
    tag "${params.panel_name}"; label 'small'
    conda "${projectDir}/env/callforge.yml"
    publishDir "${params.outdir}/stage7_filter/plots", mode: params.publish_mode
    input:
    path counts_json
    output:
    path "*.png", emit: png
    path "*.svg", emit: svg
    path "captions_filter.tsv", emit: captions
    script:
    """
    export CALLFORGE_CAPTIONS=captions_filter.tsv
    plot_filter_qc.py --counts-json ${counts_json} --outdir .
    """
}

// ── DeepVariant + GLnexus (optional alternative; container-only) ──────────────
// Selected with --caller deepvariant; run under -profile ...,docker|apptainer.
// Not exercised on the `test` profile (DeepVariant ships only as a container).
process DEEPVARIANT {
    tag "${sample_id}"; label 'large'
    container 'google/deepvariant:1.6.1'
    publishDir "${params.outdir}/stage6_calling/gvcf", mode: params.publish_mode
    input:
    tuple val(sample_id), path(bam), path(bai)
    path fasta
    path fai
    path bed
    output:
    tuple val(sample_id), path("${sample_id}.g.vcf.gz"), path("${sample_id}.g.vcf.gz.tbi"), emit: gvcf
    script:
    """
    /opt/deepvariant/bin/run_deepvariant --model_type=WES \\
        --ref=${fasta} --reads=${bam} --regions=${bed} \\
        --output_gvcf=${sample_id}.g.vcf.gz --output_vcf=${sample_id}.vcf.gz \\
        --num_shards=${task.cpus}
    """
}

process GLNEXUS {
    tag "${params.panel_name}"; label 'large'
    conda "${projectDir}/env/deepvariant.yml"
    publishDir "${params.outdir}/stage6_calling", mode: params.publish_mode
    input:
    path gvcfs
    path tbis
    path bed
    output:
    tuple path("joint.vcf.gz"), path("joint.vcf.gz.tbi"), emit: vcf
    script:
    """
    glnexus_cli --config DeepVariantWES --bed ${bed} *.g.vcf.gz > joint.bcf
    bcftools view joint.bcf | bgzip > joint.vcf.gz
    tabix -p vcf joint.vcf.gz
    """
}

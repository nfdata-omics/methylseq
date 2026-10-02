process METHYLKIT_ANNOTATION {

    tag "${caller}:${comparison_id}:${result_type}"

    container 'docker.io/yussab/methylkit:1.0'
    publishDir "${params.outdir}/${caller}/annotation/${comparison_id}", mode: 'copy'
    label 'process_low'

    input:
    tuple val(caller), val(comparison_id), val(result_type), path(dmr_tsv)
    path refseq_bed 
    path cpg_bed 

    output:
    tuple val(caller), val(comparison_id), val(result_type), path("*.tsv"), emit: tsv
    tuple val(caller), val(comparison_id), val(result_type), path("*.pdf"), emit: pdf

    script:
    """
    methylkit_annotation.R \\
        --dmr_tsv    ${dmr_tsv} \\
        --refseq_bed ${refseq_bed} \\
        --cpg_bed    ${cpg_bed} \\
        --out_prefix ${comparison_id}_${caller}_${result_type}
    """
}

process METHYLKIT_ANNOTATION {

    tag "${caller}:${comparison_id}:${result_type}"

    container 'docker.io/yussab/methylkit:1.0'
    publishDir "${params.outdir}/${caller}/annotation/${comparison_id}", mode: 'copy'
    label 'process_low'

    input:
    tuple val(caller), val(comparison_id), val(result_type), path(dmr_tsv)
    path refseq_bed 
    path cpg_bed 
    path gtf_file
    val prom_up
    val prom_down
    val diff_cutoff
    val qvalue_cutoff

    output:
    tuple val(caller), val(comparison_id), val(result_type), path("*.tsv"), emit: tsv
    tuple val(caller), val(comparison_id), val(result_type), path("*.pdf"), emit: pdf
    tuple val(caller), val(comparison_id), val(result_type), path("*.txt"), emit: txt

    script:
    """
    methylkit_annotation.R \\
        --dmr_tsv    ${dmr_tsv} \\
        --caller     '${caller}' \\
        --result_type '${result_type}' \\
        --refseq_bed ${refseq_bed} \\
        --cpg_bed    ${cpg_bed} \\
        --gtf_file   ${gtf_file} \\
        --prom_up ${prom_up} \\
        --prom_down ${prom_down} \\
        --diff_cutoff ${diff_cutoff} \\
        --qvalue_cutoff ${qvalue_cutoff} \\
        --out_prefix ${comparison_id}_${caller}_${result_type}
    """
}

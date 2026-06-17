process METHYLKIT_DMA {

    tag "${comparison_id}"

    container 'docker.io/yussab/methylkit:1.0'
    publishDir "${params.outdir}/methylkit/dma", mode: 'copy'
    label 'process_low'

    input:
    path (meth_rda)
    path (methylDB_dir) , stageAs: "methylDB_dir/*"
    tuple val(comparison_id), val(case_samples), val(control_samples)
    val  diff_cutoff
    val  qvalue_cutoff
    val  overdispersion
    val  adjust
    //val  cores

    output:
    tuple val(comparison_id), path ("*_diffmeth_raw.tsv") , emit: raw_data
    tuple val(comparison_id), path ("*_diffmeth_hyper.tsv"), emit: diff_meth_hyper
    tuple val(comparison_id), path ("*_diffmeth_hypo.tsv") , emit: diff_meth_hypo
    tuple val(comparison_id), path ("*_diffmeth_all.tsv")  , emit: diff_meth_all
    tuple val(comparison_id), path ("*.pdf"), emit: pdf

    script:
    """
    methylkit_dma.R \\
    --meth_rda ${meth_rda} \\
    --comparison_id '${comparison_id}' \\
    --case_samples '${case_samples.join(",")}' \\
    --control_samples '${control_samples.join(",")}' \\
    --diff_cutoff  ${diff_cutoff} \\
    --qvalue_cutoff  ${qvalue_cutoff} \\
    --overdispersion ${overdispersion} \\
    --adjust ${adjust} \\
    --cores ${task.cpus}
    """
}

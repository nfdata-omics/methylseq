process METHYLKIT_DMA {

    //tag "${meta.id}"
    container 'docker.io/yussab/methylkit:1.0'
    //publishDir "${params.outdir}/methylkit/${meta.id}", mode: 'copy'
    publishDir "${params.outdir}/methylkit/dma", mode: 'copy'

    cpus { cores }
    memory '8 GB'
    time '24h'

    input:
    //tuple val(meta) , path (meth_rda)
    //tuple val(meta2) , path (methylDB_dir) , stageAs: "methylDB_dir/*"
    path (meth_rda)
    path (methylDB_dir) , stageAs: "methylDB_dir/*"
    val  diff_cutoff
    val  qvalue_cutoff
    val  overdispersion
    val  adjust
    val  cores

    output:
    //tuple val(meta) , path ("diffMeth_hyper.tsv"), emit: diff_meth_hyper
    //tuple val(meta) , path ("diffMeth_hypo.tsv") , emit: diff_meth_hypo
    //tuple val(meta) , path ("diffMeth_all.tsv")  , emit: diff_meth_all
    //tuple val(meta) , path ("*.pdf")
    path ("myDiff_df.tsv") , emit: raw_data
    path ("diffMeth_hyper.tsv"), emit: diff_meth_hyper
    path ("diffMeth_hypo.tsv") , emit: diff_meth_hypo
    path ("diffMeth_all.tsv")  , emit: diff_meth_all
    path ("*.pdf")

    script:
    """
    methylkit_dma.R \\
    --meth_rda ${meth_rda} \\
    --diff_cutoff  ${diff_cutoff} \\
    --qvalue_cutoff  ${qvalue_cutoff} \\
    --overdispersion ${overdispersion} \\
    --adjust ${adjust}
    """
}

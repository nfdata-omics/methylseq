process METHYLKIT_DMA {

    tag "${meta_id}"
    container 'docker.io/yussab/methylkit:1.0'
    publishDir "${params.outdir}/methylkit", mode: 'copy'

    cpus { cores }
    memory '8 GB'
    time '24h'

    input:
    path meth_rda
    path methylDB_dir , stageAs: "methylDB_dir/*"
    val  diff_cutoff
    val  qvalue_cutoff
    val  overdispersion
    val  adjust

    output:
    path "*.tsv"
    path "*.pdf"

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

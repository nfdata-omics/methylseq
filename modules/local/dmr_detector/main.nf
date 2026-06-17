process DMR_DETECTOR {

    //tag "${meta_id}"
    container 'docker.io/yussab/...'
    publishDir "${params.outdir}/dmr_detector", mode: 'copy'

    //cpus { cores }
    memory '8 GB'
    time '24h'

    input:


    output:


    script:
    """

    """
}

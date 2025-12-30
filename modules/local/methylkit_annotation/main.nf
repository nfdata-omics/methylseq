process METHYLKIT_ANNOTATION {

    tag "${meta_id}"
    container 'docker.io/yussab/methylkit:1.0'
    publishDir "${params.outdir}/methylkit", mode: 'copy'

    cpus { cores }
    memory '8 GB'
    time '24h'

    input:
    path dmr_tsv
    path refseq_bed 
    path cpg_bed 

    output:
    path "*.tsv"
    path "*.pdf"

    script:
    """
    methylkit_annotation.R \\
        --dmr_tsv    ${dmr_tsv} \\
        --refseq_bed ${refseq_bed} \\
        --cpg_bed    ${cpg_bed} \\
        --out_prefix Sallese
    """
}

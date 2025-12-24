process METHYLKIT_QC {

    tag "${meta_id}"
    container 'docker.io/yussab/methylkit:1.0'
    publishDir "${params.outdir}/methylkit", mode: 'copy'

    cpus { cores }
    memory '8 GB'
    time '24h'

    input:
    path covfiles, stageAs: "cov_dir/*"
    path metadata
    val  group_column
    val  group_case
    val  assembly
    val  cores
    val  sample_suffix

    output:
    path "meth_merged_data.rda"
    path "*.pdf"

    script:
    """
    methylkit_qc.R \\
        --cov_dir cov_dir \\
        --metadata ${metadata} \\
        --sample_suffix ${sample_suffix} \\
        --group_column ${group_column} \\
        --group_case ${group_case} \\
        --assembly ${assembly} \\
        --cores ${cores}
    """
}

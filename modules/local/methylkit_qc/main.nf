process METHYLKIT_QC {

     tag "${group_column}_${group_case}"
    container 'docker.io/yussab/methylkit:1.0'
    publishDir "${params.outdir}/methylkit/${group_column}_${group_case}", mode: 'copy'

    cpus { cores }
    memory '8 GB'
    time '24h'

    input:
    path covfiles, stageAs: "cov_dir/*"
    path metadata
    val group_column
    val group_case
    val assembly
    val cores
    val sample_suffix
    val lo_count
    val lo_perc
    val hi_count
    val hi_perc
    val destrand
    val min_per_group

    output:
    tuple val([ id: "${group_column}_${group_case}" ]), path("meth_merged_data.rda"), emit: meth_norm_rda
    tuple val([ id: "${group_column}_${group_case}" ]), path("methylDB"),             emit: methylDB_dir
    tuple val([ id: "${group_column}_${group_case}" ]), path("*.pdf"),                emit: pdf

    script:
    """
    methylkit_qc.R \\
        --cov_dir cov_dir \\
        --metadata ${metadata} \\
        --sample_suffix ${sample_suffix} \\
        --group_column ${group_column} \\
        --group_case ${group_case} \\
        --assembly ${assembly} \\
        --cores ${cores} \\
        --lo_count ${lo_count} \\
        --lo_perc ${lo_perc} \\
        --hi_count ${hi_count} \\
        --hi_perc ${hi_perc} \\
        --destrand ${destrand} \\
        --min_per_group ${min_per_group}
    """
}

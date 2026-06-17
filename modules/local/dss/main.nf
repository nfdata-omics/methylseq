process DSS_DML_DMR {

    container 'docker.io/yussab/dss:1.0-amd64'
    publishDir "${params.outdir}/dss", mode: 'copy'

    //cpus { cores }
    memory '8 GB'
    time '6h'

    input:
    path covfiles, stageAs: "cov_dir/*"
    path metadata
    val group_column
    val group_case
    //val group1
    //val group2
    //*val group_case*
    //val cores
    val sample_suffix

    //val pattern
    /*val sep
    val min_coverage
    val smoothing_single
    val smoothing_replicates
    val smoothing_span
    val p_threshold
    val delta
    val minlen
    val minCG
    val dis_merge*/

    output:
    path("*_DML.tsv"),     emit: dml
    path("*_DMR.tsv"),     emit: dmr
    path("*_summary.tsv"), emit: summary
    path("*_objects.rda"), emit: rda

    script:
    //def smoothing_span_arg = smoothing_span == null || smoothing_span == 'NA'
    //    ? ''
    //    : "--smoothing_span ${smoothing_span}"

    """
    dss.R \\
        --cov_dir cov_dir \\
        --metadata ${metadata} \\
        --sample_suffix ${sample_suffix} \\
        --group_column ${group_column} \\
        --group_case ${group_case}
    """
}

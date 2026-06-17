process DSS_DML_DMR {

    tag "${comparison_id}"

    container 'docker.io/yussab/dss:1.0-amd64'
    publishDir "${params.outdir}/dss", mode: 'copy'
    label 'process_medium'

    input:
    path covfiles, stageAs: "cov_dir/*"
    path metadata

    tuple val(comparison_id), val(case_samples), val(control_samples)

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
    tuple val(comparison_id), path("*_DML.tsv"),     emit: dml
    tuple val(comparison_id), path("*_DMR.tsv"),     emit: dmr
    tuple val(comparison_id), path("*_summary.tsv"), emit: summary
    tuple val(comparison_id), path("*_objects.rda"), emit: rda

    script:
    //def smoothing_span_arg = smoothing_span == null || smoothing_span == 'NA'
    //    ? ''
    //    : "--smoothing_span ${smoothing_span}"

    """
    dss.R \\
        --cov_dir cov_dir \\
        --metadata ${metadata} \\
        --sample_suffix ${comparison_id} \\
        --comparison_id '${comparison_id}' \\
        --case_samples '${case_samples.join(",")}' \\
        --control_samples '${control_samples.join(",")}'
    """
}

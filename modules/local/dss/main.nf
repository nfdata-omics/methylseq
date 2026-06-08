process DSS_DML_DMR {

    //container 'docker.io/yussab/dss:1.0'
    publishDir "${params.outdir}/dss", mode: 'copy'

    cpus { cores }
    memory '8 GB'
    time '24h'

    input:
    path covfiles, stageAs: "cov_dir/*"
    path metadata
    val group_column
    val group1
    val group2
    val cores
    val sample_suffix
    val pattern
    val sep
    val min_coverage
    val smoothing_single
    val smoothing_replicates
    val smoothing_span
    val p_threshold
    val delta
    val minlen
    val minCG
    val dis_merge

    output:
    path("*_DML.tsv"),     emit: dml
    path("*_DMR.tsv"),     emit: dmr
    path("*_summary.tsv"), emit: summary
    path("*_objects.rda"), emit: rda

    script:
    def smoothing_span_arg = smoothing_span == null || smoothing_span == 'NA'
        ? ''
        : "--smoothing_span ${smoothing_span}"

    """
    dss_dml_dmr.R \\
        --cov_dir cov_dir \\
        --metadata ${metadata} \\
        --sample_suffix '${sample_suffix}' \\
        --pattern '${pattern}' \\
        --group_column '${group_column}' \\
        --group1 '${group1}' \\
        --group2 '${group2}' \\
        --sep '${sep}' \\
        --min_coverage ${min_coverage} \\
        --smoothing_single ${smoothing_single} \\
        --smoothing_replicates ${smoothing_replicates} \\
        ${smoothing_span_arg} \\
        --p_threshold ${p_threshold} \\
        --delta ${delta} \\
        --minlen ${minlen} \\
        --minCG ${minCG} \\
        --dis_merge ${dis_merge} \\
        --out_prefix ${group1}_vs_${group2}
    """
}
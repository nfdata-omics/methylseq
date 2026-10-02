process DMR_DETECTOR {

    tag "${comparison_id}"

    container 'python:3.12-slim'
    publishDir "${params.outdir}/dmr_detector/${comparison_id}", mode: 'copy'
    label 'process_medium'

    input:
    path detector_script
    path methylation_matrix
    path detector_config
    tuple val(comparison_id), val(case_samples), val(control_samples)

    output:
    tuple val(comparison_id), path("hyper.txt"),     emit: hyper
    tuple val(comparison_id), path("hypo.txt"),      emit: hypo
    tuple val(comparison_id), path("hyper-win.txt"), emit: hyper_window
    tuple val(comparison_id), path("hypo-win.txt"),  emit: hypo_window

    script:
    def case_sample = case_samples instanceof List ? case_samples[0] : case_samples
    def control_sample = control_samples instanceof List ? control_samples[0] : control_samples
    if ((case_samples instanceof List && case_samples.size() != 1) || (control_samples instanceof List && control_samples.size() != 1)) {
        error "DMR_DETECTOR currently supports one case sample and one control sample per comparison: ${comparison_id}"
    }

    """
    python3 ${detector_script} \\
        ${detector_config} \\
        ${methylation_matrix} \\
        --test_sample '${case_sample}' \\
        --ctrl_sample '${control_sample}'
    """
}

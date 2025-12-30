/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { softwareVersionsToYAML } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { METHYLKIT_QC         } from '../../modules/local/methylkit_qc/main.nf'
include { METHYLKIT_DMA        } from '../../modules/local/methylkit_differential_meth_analysis/main.nf'
include { METHYLKIT_ANNOTATION } from '../../modules/local/methylkit_annotation/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DIFFERENTIAL_METHYLATION_ANALYSIS {

    take:
    ch_versions
    methylation_coverage
    metadata_ch

    main:

    // Get methylation coverage files into a single channel
    methylation_coverage
        .map { meta, cov_file -> cov_file }       // Extract the file paths
        .collect()  // Collect all files to a single directory
        .set { ch_cov_dir }

    /*ch_cov_dir
        .map { it ->
            "My values are:\n$it\n"
        }
        .collectFile(
            name: 'ch_cov_dir.txt',
            storeDir: '.',
            keepHeader: true,
            skip: 1
        )*/

    ////////////////DEV

    METHYLKIT_QC(
        ch_cov_dir,
        metadata_ch,
        params.group_column,
        params.group_case,
        params.assembly,
        params.cores,
        params.sample_suffix
        )

    //////////////

    METHYLKIT_DMA(
        METHYLKIT_QC.out.meth_norm_rda ,
        METHYLKIT_QC.out.methylDB_dir , 
        params.diff_cutoff,
        params.qvalue_cutoff,
        params.overdispersion,  
        params.adjust 
        )

    //////////////

    METHYLKIT_ANNOTATION (
        METHYLKIT_DMA.out.diff_meth_all ,
        params.refseq_bed,
        params.cpg_bed
    )


    /*
     * Collate and save software versions
     */
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_methylseq_software_mqc_versions.yml',
            sort: true,
            newLine: true
        )
        .set { ch_collated_versions }

    emit:
    versions = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

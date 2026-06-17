/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { softwareVersionsToYAML } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { METHYLKIT_QC         } from '../../modules/local/methylkit_qc/main.nf'
include { METHYLKIT_DMA        } from '../../modules/local/methylkit_differential_meth_analysis/main.nf'
include { METHYLKIT_ANNOTATION } from '../../modules/local/methylkit_annotation/main.nf'
include { DSS_DML_DMR          } from '../../modules/local/dss/main.nf'
//include { DMR_DETECTOR         } from '../../modules/local/dmr_detector/main.nf'

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

    def legacy_methods = []
    if (params.methylkit) {
        legacy_methods << 'methylkit'
    }
    if (params.dss) {
        legacy_methods << 'dss'
    }
    if (params.dmr_detector) {
        legacy_methods << 'dmr_detector'
    }

    def selected_methods = legacy_methods ?: (
        params.differential_analysis_methods instanceof Collection
            ? params.differential_analysis_methods
            : params.differential_analysis_methods.toString().tokenize(',')
    )

    selected_methods = selected_methods
        .collect { it.toString().trim().toLowerCase().replace('-', '_') }
        .findAll { it }
        .unique()

    if (selected_methods.contains('all')) {
        selected_methods = ['methylkit', 'dss']
    }

    def supported_methods = ['methylkit', 'dss', 'dmr_detector']
    def implemented_methods = ['methylkit', 'dss']
    def unknown_methods = selected_methods.findAll { !supported_methods.contains(it) }
    def unavailable_methods = selected_methods.findAll { !implemented_methods.contains(it) }

    if (!selected_methods) {
        error "No differential methylation method selected. Set --differential_analysis_methods methylkit,dss"
    }

    if (unknown_methods) {
        error "Unsupported differential methylation method(s): ${unknown_methods.join(', ')}. Supported values: ${supported_methods.join(', ')}"
    }

    if (unavailable_methods) {
        error "Differential methylation method(s) not implemented yet: ${unavailable_methods.join(', ')}"
    }

    log.info "Differential methylation method(s): ${selected_methods.join(', ')}"

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

    if (selected_methods.contains('methylkit')) {

        //Data preprocessing 
        METHYLKIT_QC(
            ch_cov_dir,
            metadata_ch,
            params.group_column,
            params.group_case,
            params.assembly,
            params.cores,
            params.sample_suffix,
            params.lo_count,
            params.lo_perc,
            params.hi_count,
            params.hi_perc,
            params.destrand,
            params.min_per_group
        )

        //Differential Methylation Analysis
        METHYLKIT_DMA(
            METHYLKIT_QC.out.meth_norm_rda ,
            METHYLKIT_QC.out.methylDB_dir , 
            params.diff_cutoff,
            params.qvalue_cutoff,
            params.overdispersion,  
            params.adjust,
            params.cores
        )

        //Annotation
        METHYLKIT_ANNOTATION (
            METHYLKIT_DMA.out.diff_meth_all ,
            params.refseq_bed,
            params.cpg_bed
        )

    }

    if (selected_methods.contains('dss')) {
    
        /////////////DEV
        DSS_DML_DMR (
            ch_cov_dir,
            metadata_ch,
            params.group_column,
            params.group_case,
            params.sample_suffix
        )
    
    }

    if (selected_methods.contains('dmr_detector')) {
    
        //Aggiungere modulo script Alberto
        //DMR_DETECTOR ()
    
    }

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
    versions = ch_collated_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

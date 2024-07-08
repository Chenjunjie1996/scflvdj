/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                 } from '../modules/nf-core/fastqc/main'

include { EXTRACT                } from '../modules/local/extract'
include { MULTIQC                } from '../modules/local/multiqc_sgr'
include { CONVERT                } from '../modules/local/convert'
include { MKVDJREF               } from '../modules/local/mkvdjref'
include { CELLRANGER             } from '../modules/local/cellranger'
include { SUMMARIZE              } from '../modules/local/summarize/summarize'
include { MATCH                  } from '../modules/local/match/match'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_scflvdj_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow scflvdj {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_samplesheet.multiMap { meta, fastq, match_barcode ->
        reads: [meta, fastq]
        match_barcode: [meta, match_barcode]
    }.set {ch_multi}


    // MODULE: Run FastQC
    if (params.run_fastqc) {
        FASTQC (
            ch_multi.reads
        )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    }

    // extract reads
    EXTRACT (
        ch_multi.reads,
        "${projectDir}/assets/",
        params.protocol,
    )
    ch_versions = ch_versions.mix(EXTRACT.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(EXTRACT.out.json.collect{it[1]})

    // convert
    CONVERT (
        EXTRACT.out.out_reads
    )

    // ref
    def imgt_name = params.imgt_name
    if (imgt_name == 'Homo_sapiens' || imgt_name == 'Mus_musculus') {
        ch_ref = "${projectDir}/assets/ref/${imgt_name}"
    } else {
        MKVDJREF ( 
            params.imgt_name,
        )
        ch_versions = ch_versions.mix(MKVDJREF.out.versions.first())
        ch_ref = MKVDJREF.out.imgt_name
    }

    // CELLRANGER VDJ
    CELLRANGER (
        CONVERT.out.convert_fq,
        ch_ref,
    )
    ch_versions = ch_versions.mix(CELLRANGER.out.versions.first())

    // SUMMARIZE
    if (params.seqtype == 'BCR' ) {
        barcode_report = MERGE.out.barcode_report_b
    } else {
        barcode_report = MERGE.out.barcode_report_t
    }
    SUMMARIZE (
        EXTRACT.out.out_reads,
        params.seqtype,
        params.coef,
        params.expected_target_cell_num,
        MERGE.out.assembled_reads,
        MERGE.out.filter_report_tsv,
        MERGE.out.annot_fa,
        barcode_report
    )
    ch_multiqc_files = ch_multiqc_files.mix(SUMMARIZE.out.json.collect{it[1]})

    // Match
    MATCH (
        ch_multi.match_barcode,
        params.seqtype,
        SUMMARIZE.out.filter_contig_csv,
        SUMMARIZE.out.filter_contig_fa
    )
    ch_multiqc_files = ch_multiqc_files.mix(SUMMARIZE.out.json.collect{it[1]})

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        "${projectDir}/multiqc_sgr/singleron_logo.png",
        "${projectDir}/multiqc_sgr/",
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

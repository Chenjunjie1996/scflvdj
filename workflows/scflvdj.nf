/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                 } from '../modules/nf-core/fastqc/main'

include { EXTRACT                } from '../modules/local/extract'
include { MULTIQC                } from '../modules/local/multiqc_sgr'
include { IMGT_DOWNLOAD          } from '../modules/local/imgt_download'
include { TRUST4                 } from '../modules/local/trust4'
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

    // ref
    def imgt_name = params.imgt_name
    if (imgt_name == 'Homo_sapiens' || imgt_name == 'Mus_musculus') {
        ch_ref = "${projectDir}/assets/imgt_ref/${imgt_name}_IMGT+C.fa"
    } else {
        IMGT_DOWNLOAD ( 
            params.imgt_name,
        )
        ch_versions = ch_versions.mix(IMGT_DOWNLOAD.out.versions.first())
        ch_ref = IMGT_DOWNLOAD.out.ref
    }

    // ch_input = Channel.fromFilePairs("/workspace/gitpod/nf-training/scflvdj_test_data/mouse_TCR/prefix*_R{1,2}.fq.gz")
    // ch_input = Channel.fromFilePairs("${EXTRACT.out.out_temp_dir}/*_R{1,2}.fq.gz")
    // ch_input | view
    // EXTRACT.out.out_temp_dir.view()
    // ch_input =  Channel.fromFilePairs("./*temp*_R{1,2}.fq.gz")
    // EXTRACT.out.out_temp_reads.multiMap { meta, fastq ->
    //     reads: [meta, fastq]
    // }.set {ch_input}
    // ch_input = ["temp0", ["$EXTRACT.out.out_temp_dir/temp0_R1.fq.gz", "$EXTRACT.out.out_temp_dir/temp0_R2.fq.gz"]]
    // EXTRACT.out.out_temp_reads.view()
    // ch_input = channel.fromFilePairs("${EXTRACT.out.out_temp_dir}".toString())
    Channel.of(
        file("${EXTRACT.out.out_temp_dir}/temp0_R1.fq.gz"),
        file("${EXTRACT.out.out_temp_dir}/temp0_R2.fq.gz"),
        file("${EXTRACT.out.out_temp_dir}/temp1_R1.fq.gz"),
        file("${EXTRACT.out.out_temp_dir}/temp1_R2.fq.gz"),
        file("${EXTRACT.out.out_temp_dir}/temp2_R1.fq.gz"),
        file("${EXTRACT.out.out_temp_dir}/temp2_R2.fq.gz"),
        file("${EXTRACT.out.out_temp_dir}/temp3_R1.fq.gz"),
        file("${EXTRACT.out.out_temp_dir}/temp3_R2.fq.gz")
    )
    .map { it -> [it.name.split('_')[0], it] }
    .groupTuple()
    .set{ ch_input }

    // trust4
    TRUST4 (
        ch_input,
        ch_ref,
    )
    ch_versions = ch_versions.mix(TRUST4.out.versions.first())

    // SUMMARIZE
    if (params.seqtype == 'BCR' ) {
        barcode_report = TRUST4.out.barcode_report_b
    } else {
        barcode_report = TRUST4.out.barcode_report_t
    }
    SUMMARIZE (
        EXTRACT.out.out_reads,
        params.seqtype,
        params.coef,
        params.expected_target_cell_num,
        TRUST4.out.assembled_reads,
        TRUST4.out.filter_report_tsv,
        TRUST4.out.annot_fa,
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

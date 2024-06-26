process TRUST4 {
    tag "$meta"
    label 'process_high'

    conda "bioconda::trust4=1.1.1"
    container 'biocontainers/trust4:1.1.1--h43eeafb_0'
    
    input:
    tuple val(meta), path(reads)
    path(ref)

    output:
    path("*_toassemble*")         , emit: candidate_reads
    path("*_report.tsv")          , emit: report_tsv
    path("*_filter_report.tsv")   , emit: filter_report_tsv
    path("*_raw.out")             , emit: raw_out
    path("*_final.out")           , emit: final_out
    path("*_cdr3.out")            , emit: cdr3_out
    path("*_assign.out")          , emit: assign_out
    path("*_barcode_report.tsv")  , emit: barcode_report
    path("*_barcode_airr.tsv")    , emit: barcode_airr
    path("*_b.csv")               , emit: barcode_report_b
    path("*_t.csv")               , emit: barcode_report_t
    path("*_assembled_reads.fa")  , emit: assembled_reads
    path("*_annot.fa")            , emit: annot_fa
    path("*_airr.tsv")            , emit: airr_tsv
    path("*_airr_align.tsv")      , emit: airr_alin
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def (r1,r2) = reads.collate(2).transpose()
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    def readformat = task.ext.readformat ?: "bc:0:25,um:26:-1"
    def run_trust4_cmd = "-u ${r2[0]} --barcode ${r1[0]} --UMI ${r1[0]} --outputReadAssignment"
    """

    run-trust4 \\
        ${run_trust4_cmd} \\
        -t $task.cpus \\
        -f ${ref} \\
        -o ${prefix} \\
        --ref ${ref} \\
        --readFormat ${readformat} \\
        $args

    awk '(\$4!~"_") && (\$4!~"?")' ${prefix}_report.tsv > ${prefix}_filter_report.tsv
    
    perl ${projectDir}/bin/trust-barcoderep-to-10X.pl ${prefix}_barcode_report.tsv ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trust4: \$(run-trust4 2>&1 | grep -o 'v[0-9.]*-r[0-9]*' | sed 's/^TRUST4//' )
    END_VERSIONS
    """
}
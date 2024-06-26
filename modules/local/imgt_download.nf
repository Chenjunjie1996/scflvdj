process IMGT_DOWNLOAD {
    tag "$imgt_name"
    label 'process_single'

    conda "bioconda::trust4=1.1.1"
    container 'biocontainers/trust4:1.1.1--h43eeafb_0'

    input:
    val imgt_name

    output:
    path "IMGT+C.fa",      emit: ref
    path  "versions.yml" , emit: versions

    script:
    """
    BuildImgtAnnot.pl $imgt_name > IMGT+C.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyfastx: \$(pyfastx --version | sed -e "s/pyfastx version //g")
    END_VERSIONS
    """
}
process CONVERT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "quay.io/singleron-rd/sccore:v0.0.0"

    input:
    tuple val(meta), path(reads)

    output:
    path "convert_fq",  emit: convert_fq
    tuple val(meta), path("*.json"),  emit: json
    path  "versions.yml" , emit: versions

    script:
    // separate forward from reverse pairs
    def (r1,r2) = reads.collate(2).transpose()
    """
    mkdir convert_fq
    convert.py \\
        --sample ${meta.id} \\
        --fq2 ${r2.join( "," )} \\
        --assets_dir ${assets_dir} \\
   

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyfastx: \$(pyfastx --version | sed -e "s/pyfastx version //g")
    END_VERSIONS
    """
}
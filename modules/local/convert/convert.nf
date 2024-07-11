process CONVERT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "quay.io/singleron-rd/sccore:v0.0.0"

    input:
    tuple val(meta), path(r1), path(r2)
    path whitelist_tenx

    output:
    tuple val(meta), path("${meta.id}_convert_fq"),  emit: convert_fq
    tuple val(meta), path("${meta.id}.barcode_convert.json"),  emit: json
    path  "versions.yml" , emit: versions

    """
    mkdir ${meta.id}_convert_fq

    convert.py \\
        --sample ${meta.id} \\
        --fq2 ${r2} \\
        --whitelist_tenx ${whitelist_tenx}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyfastx: \$(pyfastx --version | sed -e "s/pyfastx version //g")
    END_VERSIONS
    """
}
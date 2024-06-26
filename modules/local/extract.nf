process EXTRACT {
    tag "$meta.id"
    label 'process_single'

    conda 'bioconda::pyfastx=2.1.0'
    container "biocontainers/pyfastx:2.1.0--py39h3d4b85c_0"

    input:
    tuple val(meta), path(reads)
    path assets_dir
    val protocol

    output:
    tuple val(meta), path("${meta.id}_R1.fq.gz"), path("${meta.id}_R2.fq.gz"), emit: out_reads
    tuple path("temp/temp*_R1.fq.gz"), path("temp/temp*_R2.fq.gz"), emit: out_temp_reads
    path("temp/"),  emit: out_temp_dir
    tuple val(meta), path("*.json"),  emit: json
    path  "versions.yml" , emit: versions

    script:
    // separate forward from reverse pairs
    def (r1,r2) = reads.collate(2).transpose()
    """
    mkdir temp

    extract.py \\
        --sample ${meta.id} \\
        --fq1 ${r1.join( "," )} \\
        --fq2 ${r2.join( "," )} \\
        --assets_dir ${assets_dir} \\
        --protocol ${protocol} 
   

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyfastx: \$(pyfastx --version | sed -e "s/pyfastx version //g")
    END_VERSIONS
    """
}
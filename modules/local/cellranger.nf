process CELLRANGER {
    tag "${meta.id}"
    label 'process_high'

    container "nf-core/cellranger:8.0.0"

    input:
    tuple val(meta), path(fq_dir)
    path  reference

    output:
    // tuple val(meta), path("**/outs/**"), emit: outs
    tuple val(meta), path("${meta.id}/outs/metrics_summary.csv"), emit: summary
    tuple val(meta), path("${meta.id}/outs/clonotypes.csv"), emit: clonotype
    tuple val(meta), path("${meta.id}/outs/filtered_contig_annotations.csv"), emit: annotation
    tuple val(meta), path("${meta.id}/outs/filtered_contig.fasta"), emit: fasta
    tuple val(meta), path("${meta.id}/outs/all_contig.bam"), emit: bam
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CELLRANGER_VDJ module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cellranger \\
        vdj \\
        --id='${prefix}' \\
        --fastqs=${fq_dir} \\
        --reference=${reference} \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CELLRANGER_VDJ module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "${meta.id}/outs/"
    touch ${meta.id}/outs/fake_file.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(echo \$( cellranger --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}
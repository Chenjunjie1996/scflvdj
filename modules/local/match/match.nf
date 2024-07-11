process MATCH {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "quay.io/singleron-rd/sccore:v0.0.0"

    input:
    val seqtype
    tuple val(meta), path(clonotype), path(annotation), path(fasta), path(match_barcode)

    output:
    tuple val(meta), path("*.json"), emit: json
    path "${meta.id}_matched_clonotypes.csv", emit: clonotype
    path "${meta.id}_matched_contig.csv", emit: annotation
    path "${meta.id}_matched_contig.fasta", emit: fasta

    script:
    """
    match.py \\
        --sample ${meta.id} \\
        --seqtype ${seqtype} \\
        --clonotype_csv ${clonotype} \\
        --annot_csv ${annotation} \\
        --contig_fasta ${fasta} \\
        --match_barcode_file ${match_barcode}
    """
}
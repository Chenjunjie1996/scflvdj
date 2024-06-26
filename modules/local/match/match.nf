process MATCH {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "quay.io/singleron-rd/sccore:v0.0.0"

    input:
    tuple val(meta), path(match_barcode)
    val seqtype
    path filter_contig_csv
    path filter_contig_fa

    output:
    tuple val(meta), path("*.json"), emit: json
    path "${meta.id}_matched_contig.csv", emit: match_contig_csv
    path "${meta.id}_matched_contig.fasta", emit: match_contig_fa
    // path "${meta.id}_matched_clonotypes.csv", emit: match_clonotype

    script:
    """
    match.py \\
        --sample ${meta.id} \\
        --contig_csv $filter_contig_csv \\
        --contig_fasta $filter_contig_fa \\
        --seqtype ${seqtype} \\
        --match_barcode_file $match_barcode
    """
}
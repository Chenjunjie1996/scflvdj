process SUMMARIZE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "quay.io/singleron-rd/sccore:v0.0.0"

    input:
    tuple val(meta), path(r1), path(r2)
    val seqtype
    val coef
    val expected_target_cell_num
    path assembled_reads
    path filter_report_tsv
    path annot_fa
    path barcode_report

    output:
    tuple val(meta), path("*.json"), emit: json
    path "${meta.id}.count.txt", emit: umi_count_txt
    path "${meta.id}_clonotypes.csv", emit: clonotype
    path "${meta.id}_all_contig.csv", emit: all_contig_csv
    path "${meta.id}_filtered_contig.csv", emit: filter_contig_csv
    path "${meta.id}_all_contig.fasta", emit: all_contig_fa
    path "${meta.id}_filtered_contig.fasta", emit: filter_contig_fa
    
    script:

    """
    summarize.py \\
        --sample ${meta.id} \\
        --seqtype ${seqtype} \\
        --coef ${coef} \\
        --expected_target_cell_num ${expected_target_cell_num} \\
        --fq2 ${r2} \\
        --assembled_reads ${assembled_reads} \\
        --filter_report_tsv ${filter_report_tsv} \\
        --annot_fa ${annot_fa} \\
        --barcode_report ${barcode_report}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pandas: \$( python -c "import pandas;print(pandas.__version__)")
    END_VERSIONS
    """
}
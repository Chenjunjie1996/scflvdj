process SUMMARIZE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "quay.io/singleron-rd/sccore:v0.0.0"

    input:
    val seqtype
    tuple val(meta), path(summary), path(clonotype), path(annotation), path(fasta), path(json)

    output:
    tuple val(meta), path("*.json"), emit: json
    path "${meta.id}.count.txt", emit: umi_count_txt
    path "${meta.id}_clonotypes.csv", emit: clonotype
    path "${meta.id}_filtered_contig.csv", emit: annotation
    path "${meta.id}_filtered_contig.fasta", emit: fasta
    
    script:

    """
    summarize.py \\
        --sample ${meta.id} \\
        --seqtype ${seqtype} \\
        --barcode_convert_json ${json} \\
        --clonotype_csv ${clonotype} \\
        --annot_csv ${annotation} \\
        --contig_fasta ${fasta} \\
        --metrics_csv ${summary}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pandas: \$( python -c "import pandas;print(pandas.__version__)")
        pysam: \$( python -c "import pysam;print(pysam.__version__)")
    END_VERSIONS
    """
}
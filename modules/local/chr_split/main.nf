process CHR_SPLIT {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(bed)
    path chrom_sizes

    output:
    tuple val(meta), path("*.bed"), emit: beds

    script:
    """
    awk '
        NR==FNR { chroms[\$1]; next }
        \$1 in chroms { print > \$1".bed" }
    ' ${chrom_sizes} ${bed}
    """
}

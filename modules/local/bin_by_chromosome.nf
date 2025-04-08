process BIN_BY_CHROMOSOME {
    tag "$meta.id"

    input:
    tuple val(meta), path(bed)
    path(chrom_sizes)

    output:
    //tuple val(meta), path("*.bed"), emit: chrom_beds
    tuple val(meta.id), path("*.binned.bed"), emit: chrom_beds

    script:
    """
    # Count the number of chromosomes in chrom_sizes
    chrom_count=\$(wc -l < ${chrom_sizes})

    if [ \$chrom_count -eq 1 ]; then
        # If there's only one chromosome, create a new file with the chromosome name
        chrom_name=\$(cut -f1 ${chrom_sizes})
        cp ${bed} \${chrom_name}.binned.bed
    else
        # If there are multiple chromosomes, use the original awk command
        awk 'NR==FNR{list[\$1]=1; next} \$1 in list{print > \$1".binned.bed"}' ${chrom_sizes} ${bed}
    fi
    """
}

process GENERATE_CONSENSUS {
    tag "$sample_group"

    input:
    tuple val(sample_group), path(bed_files)

    output:
    path("${sample_group}_compartments_consensus.bed"), emit: consensus

    // the \ before $ ( \$ )prevents groovy interpolation

    script:
    """
    echo "track name=\\"${sample_group}\\" description=\\"${sample_group}\\" visibility=1 itemRgb=\\"On\\"" > ${sample_group}_compartments_consensus.bed

    cat ${bed_files.join(' ')} | grep -v '^track' | sort -k1,1V -k2,2n | \\
    awk '
    {
        bin = \$1 "\\t" \$2 "\\t" \$3                                       ## bin columns (chrom, start, end)
        igv_fields[bin] = \$5 "\\t" \$6 "\\t" \$7 "\\t" \$8                 ## igv columns constant
        if (\$4 == "A") a[bin]++                                            ## 4 column is the compartment column in the combined bed
        else if (\$4 == "B") b[bin]++
    }
    END {
        for (bin in igv_fields) {
            if (a[bin] > b[bin]) comp = "A"
            else if (b[bin] > a[bin]) comp = "B"
            else comp = "NA"

            color = (comp == "A") ? "90,149,143,255" : (comp == "B") ? "224,170,88,255" : "255,255,255"

            print bin "\\t" comp "\\t" igv_fields[bin] "\\t" color          ## bin (chr,start,end) | comp (A / B / NA) | igv columns + colors
        }
    }' | sort -k1,1V -k2,2n >> ${sample_group}_compartments_consensus.bed
    """
}

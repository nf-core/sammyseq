process GENERATE_CONSENSUS {
    tag "$sample_group"

    input:
    tuple val(sample_group), path(bed_files)

    output:
    path("${sample_group}_consensus_majority.bed"), emit: consensus_majority
    path("${sample_group}_consensus_strict.bed"), emit: consensus_strict

    script:
    def num_samples = bed_files.size()
    """
    # MAJORITY: more A or more B
    echo "track name=\\"${sample_group}_majority\\"" > ${sample_group}_consensus_majority.bed
    cat ${bed_files.join(' ')} | grep -v '^track' | sort -k1,1V -k2,2n | \\
    awk '{
        bin = \$1 "\\t" \$2 "\\t" \$3
        fields[bin] = \$5 "\\t" \$6 "\\t" \$7 "\\t" \$8
        if (\$4 == "A") a[bin]++
        else if (\$4 == "B") b[bin]++
    }
    END {
        for (bin in fields) {
            if (a[bin] > b[bin]) comp = "A"
            else if (b[bin] > a[bin]) comp = "B"
            else comp = "NA"
            color = (comp == "A") ? "207,207,207" : (comp == "B") ? "69,117,180" : "255,255,255"
            print bin "\\t" comp "\\t" fields[bin] "\\t" color
        }
    }' | sort -k1,1V -k2,2n >> ${sample_group}_consensus_majority.bed

    # STRICT: all A or all B
    echo "track name=\\"${sample_group}_strict\\"" > ${sample_group}_consensus_strict.bed
    cat ${bed_files.join(' ')} | grep -v '^track' | sort -k1,1V -k2,2n | \\
    awk -v n=${num_samples} '{
        bin = \$1 "\\t" \$2 "\\t" \$3
        fields[bin] = \$5 "\\t" \$6 "\\t" \$7 "\\t" \$8
        if (\$4 == "A") a[bin]++
        else if (\$4 == "B") b[bin]++
    }
    END {
        for (bin in fields) {
            if (a[bin] == n) comp = "A"
            else if (b[bin] == n) comp = "B"
            else comp = "NA"
            color = (comp == "A") ? "207,207,207" : (comp == "B") ? "69,117,180" : "255,255,255"
            print bin "\\t" comp "\\t" fields[bin] "\\t" color
        }
    }' | sort -k1,1V -k2,2n >> ${sample_group}_consensus_strict.bed
    """
}
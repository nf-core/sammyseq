process CALL_SUBCOMPARTMENTS {
    container 'docker.io/ciuki97/sammy_subcompartments_env:latest'
    label 'process_medium'

    input:
    val tsv_content
    val binsize
    path gene_bed
    path chrom_beds

    output:
    path "*_compartment.Rdata", emit: rdata
    path "*___*.Rdata", emit: aux_rdata

    script:
    """
    echo '${tsv_content.join("\n")}' > compartments_input.tsv

    call_subcompartments.R \\
        --input_file compartments_input.tsv \\
        --binsize ${binsize} \\
        --gene_bed ${gene_bed} \\
        --chrom_beds ${chrom_beds.join(",")}
    """
}

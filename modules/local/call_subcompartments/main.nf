process CALL_SUBCOMPARTMENTS {
    container 'docker.io/ciuki97/sammy_subcompartments_env:latest'
    publishDir "${params.outdir}/subcompartments", mode: 'copy'

    input:
    val tsv_content
    val binsize
    path gene_bed
    path chrom_beds

    output:
    path "Data/*_compartment.Rdata", optional: true, emit: rdata
    path "Output/**", optional: true
    path "AUX/**", optional: true

    script:
    """
    echo '${tsv_content.join("\n")}' > compartments_input.tsv

    call_subcompartments.R \\
        --input_file compartments_input.tsv \\
        --cores ${task.cpus} \\
        --binsize ${binsize} \\
        --gene_bed ${gene_bed} \\
        --chrom_beds ${chrom_beds.join(",")}
    """
}

process CALL_SUBCOMPARTMENTS {
    container 'docker.io/ciuki97/sammy_subcompartments_env:latest'
    label 'process_medium'
    errorStrategy 'terminate'
    maxRetries 0

    input:
    val tsv_content
    val binsize
    path gene_gtf
    path chrom_bed
    val patient // replica unica che si piglia da ch_chromosomes_patients.map { it[1] }

    output:
    path "*_compartment___*_*.Rdata", emit: rdata
    path "*___*.Rdata", emit: aux_rdata
    path "*_compartments.bed", emit: bed_files
    path "*_comp_eigenvector.bedgraph", emit: bedgraph_files

    script:
    """
    echo '${tsv_content.join("\n")}' > compartments_input.tsv

    call_subcompartments.R \\
        --input_file compartments_input.tsv \\
        --binsize ${binsize} \\
        --gene_gtf ${gene_gtf} \\
        --chrom_bed ${chrom_bed} \\
        --patient ${patient}
    """
}

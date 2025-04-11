process CALL_SUBCOMPARTMENTS {
    container 'docker.io/ciuki97/sammy_subcompartments_env:latest'
    label 'process_medium'
    errorStrategy 'terminate'
    maxRetries 0

    input:
    val tsv_content
    val binsize
    path gene_gtf
    tuple val(meta), val(bedLines), val(patient)

    output:
    //path "*_compartment___*_*.Rdata", emit: rdata
    //path "*___*.Rdata", emit: aux_rdata
    path "*_compartments.bed", emit: bed_files
    path "*_comp_eigenvector.bedgraph", emit: bedgraph_files
    //path "*.bed", emit: binned_bed

    script:
    def chrom_bed = "${meta.chromosome}_binned.bed"
    """
    echo '${tsv_content.join("\n")}' > compartments_input.tsv

    # Create the binned BED file
    echo '${bedLines.join("\n")}' > ${chrom_bed}

    call_subcompartments.R \\
        --input_file compartments_input.tsv \\
        --binsize ${binsize} \\
        --gene_gtf ${gene_gtf} \\
        --chrom_bed ${chrom_bed} \\
        --patient ${patient} \\
        --chromosome ${meta.chromosome}
    """
}

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

    script:
    def chrom_bed_name = "${meta.chromosome}_binned.bed"

    def echo_tsv = tsv_content.collect { it.replace('$', '\\$') }.join('\\n')
    def echo_bed = bedLines.collect { it.replace('$', '\\$') }.join('\\n')

    """
    echo -e "${echo_tsv}" > compartments_input.tsv
    echo -e "${echo_bed}" > ${chrom_bed_name}

    call_subcompartments.R \\
        --input_file compartments_input.tsv \\
        --binsize ${binsize} \\
        --gene_gtf ${gene_gtf} \\
        --chrom_bed ${chrom_bed_name} \\
        --patient ${patient} \\
        --chromosome ${meta.chromosome}
    """
}

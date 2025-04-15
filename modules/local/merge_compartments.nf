process MERGE_COMPARTMENTS {

    input:
    tuple val(sample_id), path(bed_files), path(bedgraph_files)

    output:
    path("${sample_id}_merged_compartments.bed"), emit: merged_beds
    path("${sample_id}_merged_compartments_eigen.bedgraph"), emit: merged_bedgraphs

    script:
    """
    sort -k1,1V -k2,2n ${bed_files.join(' ')} > ${sample_id}_merged_compartments.bed
    sort -k1,1V -k2,2n ${bedgraph_files.join(' ')} > ${sample_id}_merged_compartments_eigen.bedgraph
    """
}

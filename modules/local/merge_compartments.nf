process MERGE_COMPARTMENTS {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bed_files), path(bedgraph_files)

    output:
    path("${sample_id}_merged_compartments.bed"), emit: merged_beds
    path("${sample_id}_merged_compartments_eigen.bedgraph"), emit: merged_bedgraphs

    script:
    """
    head -n 1 ${bed_files[0]} > ${sample_id}_merged_compartments.bed
    cat ${bed_files.join(' ')} | grep -v '^track' | sort -k1,1V -k2,2n >> ${sample_id}_merged_compartments.bed

    head -n 1 ${bedgraph_files[0]} > ${sample_id}_merged_compartments_eigen.bedgraph
    cat ${bedgraph_files.join(' ')} | grep -v '^track' | sort -k1,1V -k2,2n >> ${sample_id}_merged_compartments_eigen.bedgraph
    """
}

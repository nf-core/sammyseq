process CHR_COMBINE_COMPARTMENTS {
    tag "$sample_id"
    label 'process_medium'

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(sample_id), path(bed_files), path(bedgraph_files)

    output:
    tuple val(sample_id), path("${sample_id}_combined_compartments.bed"), emit: combined_beds
    tuple val(sample_id), path("${sample_id}_combined_compartments_eigen.bedgraph"), emit: combined_bedgraphs
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Combine BED: pull together all the bed files and sort them by chromosome and start position. The grep commands remove any existing track or comment lines from the input files before sorting.
    { echo 'track name="${sample_id}" description="${sample_id} (Emission ordered)" visibility=1 itemRgb="On"'; cat ${bed_files.join(' ')} | grep -v '^track' | grep -v '^#' | sort -k1,1V -k2,2n; } > ${sample_id}_combined_compartments.bed

    # Combine BedGraph:pull together all the bedgraph files and sort them by chromosome and start position. The grep commands remove any existing track, browser, or comment lines from the input files before sorting.
    { echo 'track type=bedGraph name="${sample_id}_eigenvector" description="${sample_id} eigenvector" visibility=full color=90,149,143 altColor=224,170,88 priority=20'; cat ${bedgraph_files.join(' ')} | grep -v '^track' | grep -v '^browser' | grep -v '^#' | sort -k1,1V -k2,2n; } > ${sample_id}_combined_compartments_eigen.bedgraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: \$(sort --version | head -n1 | sed 's/^.* //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}_combined_compartments.bed
    touch ${sample_id}_combined_compartments_eigen.bedgraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: \$(sort --version | head -n1 | sed 's/^.* //g')
    END_VERSIONS
    """
}

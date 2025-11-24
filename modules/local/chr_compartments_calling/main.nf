process CHR_COMPARTMENTS_CALLING {
    tag "${patient}_${meta.chromosome}"
    label 'process_medium'
    errorStrategy 'terminate'
    maxRetries 0

    container 'docker.io/ciuki97/sammy_subcompartments_env:latest'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(chr_bed), val(patient), path(csv)
    val binsize
    path gtf

    output:
    tuple val(patient), path("*.bed")      , emit: bed_files
    tuple val(patient), path("*.bedgraph") , emit: bedgraph_files
    path "versions.yml"                    , emit: versions

    script:
    def args = task.ext.args ?: ''
    template 'chr_compartments_calling.R'
}
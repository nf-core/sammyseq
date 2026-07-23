process CHR_COMPARTMENTS_CALLING {
    tag "${patient}_${meta.chromosome}"
    label 'process_medium'

    container 'ghcr.io/daisymut/sammyr:0.0.0.9000'
    // conda "${moduleDir}/environment.yml"

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

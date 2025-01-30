process CUT_SIZES_GENOME {
    tag "${sample_id}"
    publishDir "${params.outdir}/genome", mode: 'copy'

    input:
    tuple val(meta), path(fai)

    output:
    path("sizes.genome") , emit: ch_sizes_genome

    script:
    """
    cut -f1,2 ${fai} > sizes.genome
    """
}


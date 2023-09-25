// Import SAMTOOLS_FAIDX module
include { SAMTOOLS_FAIDX } from '../modules/nf-core/samtools/faidx'

// Definition of the CUT_SIZES_GENOME process
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

// Definition of the subworkflow
workflow FAIDX_SUBWORKFLOW {
    take:
    path fasta_file

    main:
    SAMTOOLS_FAIDX(fasta_file)
    CUT_SIZES_GENOME(SAMTOOLS_FAIDX.out.fai)

    emit:
    path("sizes.genome"), from: CUT_SIZES_GENOME.out.ch_sizes_genome
}
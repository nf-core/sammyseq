
process RTWOSAMPLESMLE {
    //tag "$meta.id"
    tag "$meta.experimentalID $meta.ratio"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'nf-core/r_two_samples_mle:0.0.1'

    input:
    tuple val(meta), path(bam1), path(bam2), val(output_mle_name)
    path chromsizes_file

    output:
    tuple val(meta), path("*_mle.bigWig"), emit: results
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'two_samples_mle.R'

}

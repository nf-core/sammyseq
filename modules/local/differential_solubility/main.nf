process DIFFERENTIAL_SOLUBILITY {
    container 'docker.io/ciuki97/sammy_subcompartments_env:latest'
    conda "${moduleDir}/environment.yml"

    label 'process_medium'
    errorStrategy 'terminate'
    maxRetries 0

    input:
    tuple val(meta), path(samplesheet)
    path(genome_bins)
    val(binsize)
    val(comparison)

    output:
    tuple val(meta), path("*_all_bins_complete.csv"), emit: all_bins
    tuple val(meta), path("*_selected_bins_filtered.csv"), emit: selected_bins

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'differential_solubility.R'
}

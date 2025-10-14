process DIFFERENTIAL_SOLUBILITY {
    container 'docker.io/ciuki97/differential-solubility-analysis:v0.0.1'
    conda "${moduleDir}/environment.yml"

    label 'process_medium'
    errorStrategy 'terminate'
    maxRetries 0

    input:
    input:
    tuple val(meta), path(samplesheet)
    path(genome_bins)
    val(binsize)
    val(comparison)
    val(compare_groups)
    val(solubility_threshold)
    val(validation_passed)
    path(gtf)

    output:
    tuple val(meta), path("*_all_bins_complete.csv"), emit: all_bins
    tuple val(meta), path("*_selected_bins_filtered.csv"), emit: selected_bins
    tuple val(meta), path("analysis_summary.txt"), emit: report
    tuple val(meta), path("*_genes.txt"), optional: true, emit: gene_lists
    tuple val(meta), path("*.rds"), emit: analysis_results

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'differential_solubility.R'
}
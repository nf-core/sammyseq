process DIFFERENTIAL_ENRICHMENT {
    tag "${contrast_data[0]}vs${contrast_data[1]}"
    container 'docker.io/ciuki97/differential-solubility-analysis:v0.0.1'
    conda "${moduleDir}/environment.yml"

    label 'process_high'
    errorStrategy 'terminate'
    maxRetries 0

input:
tuple val(meta), path(samplesheet)
each contrast_data              // ← [test_group, ref_group]
path(genome_bins)
path(gtf)
val(binsize)
val(comparison)
val(solubility_threshold)
val(validation_passed)



    output:
    tuple val(meta), path("*_all_bins_complete.csv"), emit: all_bins
    tuple val(meta), path("*_selected_bins_filtered.csv"), emit: selected_bins
    tuple val(meta), path("rdata/*.rds"), emit: selected_bins_rds
    tuple val(meta), path("regions/*.bed"), emit: bed_regions
    tuple val(meta), path("analysis_summary.txt"), emit: report
    tuple val(meta), path("genes/*_genes.txt"), optional: true, emit: gene_lists

    when:
    task.ext.when == null || task.ext.when

script:
test_group = contrast_data[0]
ref_group = contrast_data[1]
contrast_name = "${contrast_data[0]}vs${contrast_data[1]}"

template 'differential_enrichment.R'
}
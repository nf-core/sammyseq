process DIFFERENTIAL_ENRICHMENT {
    tag "${contrast_data[0]}vs${contrast_data[1]}"
    container 'ghcr.io/daisymut/sammyr:0.0.0.9001'
    // conda "${moduleDir}/environment.yml"
    label 'process_high'

    input:
    tuple val(meta), path(samplesheet)
    each contrast_data
    path(genome_bins)
    path(gtf)
    val(binsize)
    val(comparison)
    val(solubility_threshold)

    output:
    tuple val(meta), path("*_all_bins_complete.csv")          , emit: all_bins
    tuple val(meta), path("*_selected_bins_filtered.csv")     , emit: selected_bins
    tuple val(meta), path("rdata/*.rds")                      , emit: selected_bins_rds
    tuple val(meta), path("regions/*.bed")                    , emit: bed_regions
    tuple val(meta), path("*_analysis_summary.txt")           , emit: report
    tuple val(meta), path("genes/*_genes.txt"), optional: true, emit: gene_lists
    path "versions.yml"                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    test_group = contrast_data[0]
    ref_group = contrast_data[1]
    contrast_name = "${contrast_data[0]}vs${contrast_data[1]}"

    template 'differential_enrichment.R'
}

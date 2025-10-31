//
// Identify differentially enriched solubility regions
//

include { DIFFERENTIAL_ENRICHMENT } from '../../../modules/local/differential_enrichment/main'

workflow DIFFERENTIAL_SOLUBILITY_ANALYSIS {

    take:
    mle_results_channel
    outdir
    genome_bins
    gtf
    binsize
    comparison
    solubility_threshold
    compare_groups

    main:

    //
    // Generate correct CSV samplesheet for differential analysis
    //
    ch_samplesheet = mle_results_channel
        .filter { csv_meta, mle_file -> csv_meta.csv_expid_filter }
        .collectFile(keepHeader: true, sort: true, storeDir: "${outdir}/csv") { csv_meta, mle_file ->
            ["mle_comparisons.csv", "experimental_id,sample_group,ratio,file\n${csv_meta.experimentalID},${csv_meta.sample_group},${csv_meta.ratio},${mle_file}\n"]
        }
        .map { file -> [[ id:'differential_analysis' ], file] }

    //
    // Parse contrasts for parallel processing
    //
    ch_contrasts = Channel
        .from(compare_groups.split(','))
        .map { contrast ->
            def parts = contrast.split('vs')
            def test_group = parts[0].trim()
            def ref_group = parts[1].trim()
            [test_group, ref_group]
        }

    //
    // Run differential enrichment analysis
    //

    DIFFERENTIAL_ENRICHMENT (
        ch_samplesheet,
        ch_contrasts,
        genome_bins,
        gtf,
        binsize,
        comparison,
        solubility_threshold
    )

    emit:
    all_bins          = DIFFERENTIAL_ENRICHMENT.out.all_bins
    selected_bins     = DIFFERENTIAL_ENRICHMENT.out.selected_bins
    selected_bins_rds = DIFFERENTIAL_ENRICHMENT.out.selected_bins_rds
    bed_regions       = DIFFERENTIAL_ENRICHMENT.out.bed_regions
    gene_lists        = DIFFERENTIAL_ENRICHMENT.out.gene_lists
    report            = DIFFERENTIAL_ENRICHMENT.out.report
    versions          = DIFFERENTIAL_ENRICHMENT.out.versions
}

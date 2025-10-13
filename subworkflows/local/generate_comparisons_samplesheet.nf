//
// GENERATE_MLE_RATIO_CSV
//

workflow GENERATE_COMPARISONS_SAMPLESHEET {
    take:
        mle_results     // channel: [mandatory] csv_meta, mle_file.bigWig
        outdir          // string: output directory path

    main:
        // Generate csv file from two_samples_mle.R step
        mle_results
            .filter { csv_meta, mle_file ->
                csv_meta.csv_expid_filter // filter out comparisons where experimentalID is not the same
            }
            .collectFile(keepHeader: true, sort: true, storeDir: "${outdir}/csv") { csv_meta, mle_file ->
                def experimental_id = csv_meta.experimentalID
                def sample_group = csv_meta.sample_group
                def ratio = csv_meta.ratio

                ["mle_comparisons.csv", "experimental_id,sample_group,ratio,file\n${experimental_id},${sample_group},${ratio},${mle_file}\n"]
            }
}

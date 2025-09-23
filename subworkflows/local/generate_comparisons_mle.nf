//
// Subworkflow to generate pairwise comparisons and run RTWOSAMPLESMLE analysis
//

include { RTWOSAMPLESMLE } from '../../modules/local/rtwosamplesmle/main'

workflow GENERATE_COMPARISONS_MLE {
    
    take:
    ch_bam_input        // channel: [meta, bam]
    ch_samplesheet      // channel: samplesheet data for comparison string approach
    chrom_sizes         // path: chromosome sizes file
    
    main:
    
    ch_versions = Channel.empty()
    
    // Initialize comparison channels
    comparisons_ch_s1 = Channel.empty()
    comparisons_ch_s2 = Channel.empty()
    
    if (params.comparisonFile) {    
        // comparisonFile CSV based approach
        Channel
            .fromPath(params.comparisonFile)
            .splitCsv(header: true)
            .map { row ->
                [row.sample1, row.sample1 + "_VS_" + row.sample2]
            }
            .set { comparisons_ch_s1 }
            
        Channel
            .fromPath(params.comparisonFile)
            .splitCsv(header: true)
            .map { row ->
                [row.sample2, row.sample1 + "_VS_" + row.sample2]
            }
            .set { comparisons_ch_s2 }

    } else if (params.comparison) { 
        // comparison string-based approach
        def comparison_list = params.comparison.split(',').collect { it.trim() }

        // Create comparison channels (one for sample1 and one for sample2 in each comparison)
        ch_samplesheet
            .map { meta, fastqs -> meta }
            .collect()
            .flatMap { meta_list ->
                comparison_list.collectMany { comp ->
                    def (frac1, frac2) = comp.split('vs')   // Split comparison string into two fractions
                    def samples_by_expID = meta_list.groupBy { it.experimentalID } // Group samples by experimental ID
                    // For each experimental ID, find samples for the two fractions and create a list of comparisons
                    samples_by_expID.collectMany { exp_id, samples ->
                        def s1 = samples.find { it.fraction == frac1 }?.id
                        def s2 = samples.find { it.fraction == frac2 }?.id

                        if (s1 && s2) {
                            return [[sample1: s1, sample2: s2]]
                        } else {
                            // If one of the fractions is missing, log a warning and skip this comparison
                            def missing = [!s1 ? frac1 : null, !s2 ? frac2 : null].findAll()
                            log.warn "Skipping ${comp} for ${exp_id}: missing ${missing.join(' and ')} fraction(s)"
                            return []
                        }
                    }
                }
            }
            .multiMap { row ->      // Split into separate channels for sample1 and sample2
                comparisons_ch_s1: [row.sample1, "${row.sample1}_VS_${row.sample2}"]
                comparisons_ch_s2: [row.sample2, "${row.sample1}_VS_${row.sample2}"]
            }
            .set { comparisons_ch }

        comparisons_ch_s1 = comparisons_ch.comparisons_ch_s1
        comparisons_ch_s2 = comparisons_ch.comparisons_ch_s2
    }

    // Convert bam file to input format
    // [[id:ggg, paired:true],path.bam] -> [id, bam, meta]
    ch_bam_input
        .map { meta, bam ->
            [meta.id, bam, meta]
        }
        .set { ch_bam_reformat }

    // Combine comparison channel with bam list channel
    comparisons_ch_s1
        .combine(ch_bam_reformat, by: 0)
        .map { sample1, comparison, bam, meta ->
            [comparison, bam, meta]
        }
        .set { bam1_comparison }

    comparisons_ch_s2
        .combine(ch_bam_reformat, by: 0)
        .map { sample2, comparison, bam, meta ->
            [comparison, bam, meta]
        }
        .set { bam2_comparison }

    // Join the two comparison channels and prepare for RTWOSAMPLESMLE
    bam1_comparison
        .join(bam2_comparison, remainder: false, by: 0)
        .map { comparison, bam1, meta_bam1, bam2, meta_bam2 ->

            def expid1 = meta_bam1.experimentalID
            def fraction1 = meta_bam1.fraction
            def expid2 = meta_bam2.experimentalID
            def fraction2 = meta_bam2.fraction

            def output_mle_name
            if (expid1 == expid2) {
                output_mle_name = "${expid1}_${fraction1}vs${fraction2}"
            } else {
                output_mle_name = "${expid1}_${fraction1}_vs_${expid2}_${fraction2}"
            }

            def meta_csv = [
                experimentalID: expid1,
                sample_group: meta_bam1.sample_group,
                ratio: "${fraction1}vs${fraction2}",
                csv_expid_filter: expid1 == expid2 // in GENERATE_MLE_RATIO_CSV filter out comparisons with different experimentalID
            ]

            [meta_csv, bam1, bam2, output_mle_name]
        }
        .set { comparisons_merge_ch }

    // Run RTWOSAMPLESMLE analysis
    RTWOSAMPLESMLE(
        comparisons_merge_ch,
        chrom_sizes
    )

    emit:
    mle_results = RTWOSAMPLESMLE.out.results           // channel: [meta, results]
}
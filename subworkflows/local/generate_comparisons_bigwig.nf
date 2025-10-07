//
// Subworkflow to generate pairwise comparisons and run DEEPTOOLS_BIGWIGCOMPARE analysis
//
 
include { DEEPTOOLS_BIGWIGCOMPARE } from '../../modules/nf-core/deeptools/bigwigcompare'

workflow GENERATE_COMPARISONS_BIGWIG {

    take:
    ch_bigwig_input        // channel: [meta, bigwig]
    //dovrei avere due bigwig come input, non capisco se qui è solo uno
    ch_samplesheet      // channel: samplesheet data for comparison string approach
    // chromsize??
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

    // ......
    // il blocco precedente è uguale per i due casi, si potrebbe modificare l'if

    ch_bigwig_input
        .map { meta, bigwig ->
            [meta.id, bigwig, meta]
        }
        .set { ch_bigwig_reformat }

    // Combine comparison channel with bigwig list channel
    comparisons_ch_s1
        .combine(ch_bigwig_reformat, by: 0)
        .map { sample1, comparison, bigwig, meta ->
            [comparison, bigwig, meta]
        }
        .set { bigwig1_comparison }

    comparisons_ch_s2
        .combine(ch_bigwig_reformat, by: 0)
        .map { sample2, comparison, bigwig, meta ->
            [comparison, bigwig, meta]
        }
        .set { bigwig2_comparison }

    // Join the two comparison channels and prepare for DEEPTOOLS_BIGWIGCOMPARE
    bigwig1_comparison
        .join(bigwig2_comparison, remainder: false, by: 0)
        .map { comparison, bigwig1, meta_bigwig1, bigwig2, meta_bigwig2 ->

            def expid1 = meta_bigwig1.experimentalID
            def fraction1 = meta_bigwig1.fraction
            def expid2 = meta_bigwig2.experimentalID
            def fraction2 = meta_bigwig2.fraction

            def output_bigwig_name
            if (expid1 == expid2) {
                output_bigwig_name = "${expid1}_${fraction1}vs${fraction2}"
            } else {
                output_bigwig_name = "${expid1}_${fraction1}_vs_${expid2}_${fraction2}"
            }

            def meta_csv = [
                experimentalID: expid1,
                sample_group: meta_bigwig1.sample_group,
                ratio: "${fraction1}vs${fraction2}",
                csv_expid_filter: expid1 == expid2 
            ]

            def meta_for_module = meta_csv + [ id: output_bigwig_name ]
            [ meta_for_module, bigwig1, bigwig2 ]

        }
        .set { comparisons_merge_ch }

    // Run DEEPTOOLS_BIGWIGCOMPARE analysis
    DEEPTOOLS_BIGWIGCOMPARE(
        comparisons_merge_ch,
        params.blacklist ? Channel.value([[:], file(params.blacklist)]) : Channel.value([[:], []])
    )

    emit:
    comparison_results = DEEPTOOLS_BIGWIGCOMPARE.out.output
}

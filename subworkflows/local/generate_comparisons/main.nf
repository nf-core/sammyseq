//
// Subworkflow to generate pairwise comparisons between fractions
//

include { RTWOSAMPLESMLE          } from '../../../modules/local/rtwosamplesmle/main'
include { DEEPTOOLS_BIGWIGCOMPARE } from '../../../modules/nf-core/deeptools/bigwigcompare/main'

workflow GENERATE_COMPARISONS {

    take:
    ch_input            // channel: [meta, bam || bigwig]
    ch_samplesheet      // channel: samplesheet data for comparison string approach
    chrom_sizes         // path:    chromosome sizes file (only for MLE)
    module_name         // string:  'spp' or 'bigwigcompare'

    main:
    ch_comparison_results = Channel.empty()
    ch_versions = Channel.empty()

    // Initialize comparison channels
    comparisons_ch_s1 = Channel.empty()
    comparisons_ch_s2 = Channel.empty()

    if (params.comparison_file) {
        // comparison_file CSV based approach
        Channel
            .fromPath(params.comparison_file)
            .splitCsv(header: true)
            .multiMap { row ->
                comparisons_ch_s1: [row.sample1, "${row.sample1}_VS_${row.sample2}"]
                comparisons_ch_s2: [row.sample2, "${row.sample1}_VS_${row.sample2}"]
            }
            .set { comparisons_ch }

        comparisons_ch_s1 = comparisons_ch.comparisons_ch_s1
        comparisons_ch_s2 = comparisons_ch.comparisons_ch_s2

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

    ch_input
        .map { meta, f ->
            [meta.id, f, meta]
        }
        .set { ch_input_reformat }

    // Combine comparison channel with input list channel
    comparisons_ch_s1
        .combine(ch_input_reformat, by: 0)
        .map { sample1, comparison, file, meta ->
            [comparison, file, meta]
        }
        .set { f1_comparison }

    comparisons_ch_s2
        .combine(ch_input_reformat, by: 0)
        .map { sample2, comparison, file, meta ->
            [comparison, file, meta]
        }
        .set { f2_comparison }

    // Join the two comparison channels and prepare for module
    f1_comparison
        .join(f2_comparison)
        .map { comparison, f1, meta1, f2, meta2 ->

            def expid1 = meta1.experimentalID
            def fraction1 = meta1.fraction
            def expid2 = meta2.experimentalID
            def fraction2 = meta2.fraction

            def output_name
            if (expid1 == expid2) {
                output_name = "${expid1}_${fraction1}vs${fraction2}"
            } else {
                output_name = "${expid1}_${fraction1}_vs_${expid2}_${fraction2}"
            }

            def meta_csv = [
                experimentalID: expid1,
                sample_group: meta1.sample_group,
                ratio: "${fraction1}vs${fraction2}",
                csv_expid_filter: expid1 == expid2,
                id: output_name
            ]

            [ meta_csv, f1, f2 ]

        }
        .set { comparisons_merge_ch }

    // parameter check to run comparison analysis
    if (module_name == 'spp') {
        RTWOSAMPLESMLE(
            comparisons_merge_ch.map { meta, f1, f2 -> [meta - [id: meta.id], f1, f2, meta.id] },
            chrom_sizes
        )
        emit:
        ch_comparison_results = RTWOSAMPLESMLE.out.results
    } else if (module_name == 'bigwigcompare') {
        def blacklist_ch = params.blacklist
            ? Channel.value([[:], file(params.blacklist)])
            : Channel.value([[:], []])

        DEEPTOOLS_BIGWIGCOMPARE(comparisons_merge_ch, blacklist_ch)
        emit:
        ch_comparison_results = DEEPTOOLS_BIGWIGCOMPARE.out.output
    }

    emit:
    results = ch_comparison_results
    versions = ch_versions
}

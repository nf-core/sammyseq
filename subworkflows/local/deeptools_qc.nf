/*
 * Perform full suite of deep tools analysis on bam/bigwig files
 */

include { DEEPTOOLS_MULTIBAMSUMMARY } from '../../modules/nf-core/deeptools/multibamsummary/main'
include { DEEPTOOLS_MULTIBIGWIGSUMMARY } from '../../modules/nf-core/deeptools/multibigwigsummary/main'
include { DEEPTOOLS_PLOTCORRELATION } from '../../modules/nf-core/deeptools/plotcorrelation/main'
include { DEEPTOOLS_PLOTPCA         } from '../../modules/nf-core/deeptools/plotpca/main'
include { DEEPTOOLS_PLOTFINGERPRINT as DEEPTOOLS_PLOTFINGERPRINT_GLOBAL } from '../../modules/nf-core/deeptools/plotfingerprint/main'
include { DEEPTOOLS_PLOTFINGERPRINT as DEEPTOOLS_PLOTFINGERPRINT_REGION } from '../../modules/nf-core/deeptools/plotfingerprint/main'

workflow DEEPTOOLS_QC {
    take:
    bam
    bai
    bigwig
    corr_method
    ch_blacklist

    main:
    ch_versions = Channel.empty()

    bam
        .join(bai)
        .set { ch_bam_bai }

    bigwig
        .map { row -> [row[0].id] }
        .collect()
        .map { row -> [row] }
        .set { ch_bigwig_ids }

    bigwig
        .map { row -> [row[1]] }
        .collect()
        .map { row -> [row] }
        .combine(ch_bigwig_ids)
        .map { row -> [[id: 'all_target_bigwigs'], row[0], row[1], row[0].size()] }
        .filter { row -> row[3] > 1 }
        .map { row -> [row[0], row[1], row[2]] }
        .set { ch_bigwig_all }

    DEEPTOOLS_MULTIBIGWIGSUMMARY(
        ch_bigwig_all,
        ch_blacklist
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_MULTIBIGWIGSUMMARY.out.versions)

    DEEPTOOLS_PLOTCORRELATION(
        DEEPTOOLS_MULTIBIGWIGSUMMARY.out.matrix,
        corr_method,
        "heatmap"
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_PLOTCORRELATION.out.versions)

    DEEPTOOLS_PLOTPCA(DEEPTOOLS_MULTIBIGWIGSUMMARY.out.matrix)
    ch_versions = ch_versions.mix(DEEPTOOLS_PLOTPCA.out.versions)

    ch_fingerprint_matrix_global = Channel.empty()
    ch_fingerprint_metrics_global = Channel.empty()
    ch_fingerprint_region_matrix = Channel.empty()
    ch_fingerprint_region_metrics = Channel.empty()

    if (params.plotfingerprint) {
        ch_bam_bai
            .map { row -> [row[0].id] }
            .collect()
            .map { row -> [row] }
            .set { ch_bam_ids }

        ch_bam_bai
            .map { row -> [row[1]] }
            .collect()
            .map { row -> [row] }
            .combine(ch_bam_bai.map { row -> [row[2]] }.collect().map { row -> [row] })
            .combine(ch_bam_ids)
            .map { row -> [[id: 'all_target_bams'], row[0], row[1], row[2], row[1].size()] }
            .filter { row -> row[4] > 1 }
            .map { row -> [row[0], row[1], row[2], row[3]] }
            .set { ch_bam_bai_all }

        DEEPTOOLS_MULTIBAMSUMMARY(
            ch_bam_bai_all,
            ch_blacklist
        )
        ch_versions = ch_versions.mix(DEEPTOOLS_MULTIBAMSUMMARY.out.versions)

        ch_bam_bai
            .map { meta, bam, bai ->
                def new_meta = [id: meta.id.split('_')[0]]
                [new_meta, bam, bai]
            }
            .groupTuple(by: [0])
            .map { meta, bams, bais -> [meta, bams.flatten(), bais.flatten()] }
            .set { ch_grouped_bam_bai }

        DEEPTOOLS_PLOTFINGERPRINT_GLOBAL(ch_grouped_bam_bai)
        ch_versions = ch_versions.mix(DEEPTOOLS_PLOTFINGERPRINT_GLOBAL.out.versions)
        
        ch_fingerprint_matrix_global = DEEPTOOLS_PLOTFINGERPRINT_GLOBAL.out.matrix
        ch_fingerprint_metrics_global = DEEPTOOLS_PLOTFINGERPRINT_GLOBAL.out.metrics

        if (params.region) {
            DEEPTOOLS_PLOTFINGERPRINT_REGION(ch_grouped_bam_bai)
            ch_versions = ch_versions.mix(DEEPTOOLS_PLOTFINGERPRINT_REGION.out.versions)
            ch_fingerprint_region_matrix = DEEPTOOLS_PLOTFINGERPRINT_REGION.out.matrix
            ch_fingerprint_region_metrics = DEEPTOOLS_PLOTFINGERPRINT_REGION.out.metrics
        }
    }

    emit:
    correlation_matrix         = DEEPTOOLS_PLOTCORRELATION.out.matrix
    pca_data                   = DEEPTOOLS_PLOTPCA.out.tab
    fingerprint_matrix_global  = ch_fingerprint_matrix_global
    fingerprint_metrics_global = ch_fingerprint_metrics_global
    fingerprint_matrix_region  = ch_fingerprint_region_matrix
    fingerprint_metrics_region = ch_fingerprint_region_metrics
    versions = ch_versions
}
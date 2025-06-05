/*
 * Perform full suite of deep tools analysis on bam files
*/

include { DEEPTOOLS_MULTIBAMSUMMARY } from '../../modules/nf-core/deeptools/multibamsummary/main'
include { DEEPTOOLS_PLOTCORRELATION } from '../../modules/nf-core/deeptools/plotcorrelation/main'
include { DEEPTOOLS_PLOTPCA         } from '../../modules/nf-core/deeptools/plotpca/main'
include { DEEPTOOLS_PLOTFINGERPRINT as DEEPTOOLS_PLOTFINGERPRINT_GLOBAL } from '../../modules/nf-core/deeptools/plotfingerprint/main'
include { DEEPTOOLS_PLOTFINGERPRINT as DEEPTOOLS_PLOTFINGERPRINT_REGION } from '../../modules/nf-core/deeptools/plotfingerprint/main'

workflow DEEPTOOLS_QC {
    take:
    bam         // channel: [ val(meta), [ bam ] ]
    bai         // channel: [ val(meta), [ bai ] ]
    corr_method // val
    ch_blacklist   // channel

    main:
    ch_versions = Channel.empty()

    /*
    * CHANNEL: Combine bam and bai files on id
    */
    bam
        .join(bai)
        .set { ch_bam_bai }

    /*
    * CHANNEL: Get list of sample ids
    */
    ch_bam_bai
        .map { row -> [row[0].id] }
        .collect()
        .map { row -> [row] }
        .set { ch_ids }

    /*
    * CHANNEL: Combine bam and bai files into one list
    * if we only have one file then cancel correlation and PCA
    */
    ch_bam_bai
        .map { row -> [row[1]] }
        .collect()
        .map { row -> [row] }
        .combine(ch_bam_bai.map { row -> [row[2]] }.collect().map { row -> [row] })
        .combine(ch_ids)
        .map { row -> [[id: 'all_target_bams'], row[0], row[1], row[2], row[1].size()] }
        .filter { row -> row[4] > 1 }
        .map { row -> [row[0], row[1], row[2], row[3]] }
        .set { ch_bam_bai_all }

    /*
    * MODULE: Summarise bams into bins
    */
    DEEPTOOLS_MULTIBAMSUMMARY(
        ch_bam_bai_all,
        ch_blacklist
)
    ch_versions = ch_versions.mix(DEEPTOOLS_MULTIBAMSUMMARY.out.versions)

    /*
    * MODULE: Plot correlation matrix
    */
    DEEPTOOLS_PLOTCORRELATION(
        DEEPTOOLS_MULTIBAMSUMMARY.out.matrix,
        corr_method,
        "heatmap"
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_PLOTCORRELATION.out.versions)

    /*
    * MODULE: Plot PCA's
    */
    DEEPTOOLS_PLOTPCA(DEEPTOOLS_MULTIBAMSUMMARY.out.matrix)
    ch_versions = ch_versions.mix(DEEPTOOLS_PLOTPCA.out.versions)

    /*
    * CHANNEL: Group BAM and BAI files by meta.id for PLOTFINGERPRINT
    */
    ch_bam_bai
        .map { meta, bam, bai ->
            def new_meta = [id: meta.id.split('_')[0]]  // Assuming meta.id format is like "CTRL003_S2"
            [new_meta, bam, bai]
        }
        .groupTuple(by: [0])
        .map { meta, bams, bais -> [meta, bams.flatten(), bais.flatten()] }
        .set { ch_grouped_bam_bai }

    /*
    * MODULE: Plot Fingerprint (Global)
    */
    DEEPTOOLS_PLOTFINGERPRINT_GLOBAL(ch_grouped_bam_bai)
    ch_versions = ch_versions.mix(DEEPTOOLS_PLOTFINGERPRINT_GLOBAL.out.versions)

    /*
    * MODULE: Plot Fingerprint (Region-specific, if params.region is defined)
    */
    ch_fingerprint_region_matrix = Channel.empty()
    ch_fingerprint_region_metrics = Channel.empty()
    if (params.region) {
        DEEPTOOLS_PLOTFINGERPRINT_REGION(ch_grouped_bam_bai)
        ch_versions = ch_versions.mix(DEEPTOOLS_PLOTFINGERPRINT_REGION.out.versions)
        ch_fingerprint_region_matrix = DEEPTOOLS_PLOTFINGERPRINT_REGION.out.matrix
        ch_fingerprint_region_metrics = DEEPTOOLS_PLOTFINGERPRINT_REGION.out.metrics
    }

    emit:
    correlation_matrix         = DEEPTOOLS_PLOTCORRELATION.out.matrix
    pca_data                   = DEEPTOOLS_PLOTPCA.out.tab
    fingerprint_matrix_global  = DEEPTOOLS_PLOTFINGERPRINT_GLOBAL.out.matrix
    fingerprint_metrics_global = DEEPTOOLS_PLOTFINGERPRINT_GLOBAL.out.metrics
    fingerprint_matrix_region  = ch_fingerprint_region_matrix
    fingerprint_metrics_region = ch_fingerprint_region_metrics

    versions = ch_versions                 // channel: [ versions.yml ]
}

//
// Genome binning subworkflow with keep_regions handling
//

include { GUNZIP as GUNZIP_KEEP_REGIONS_BED } from '../../modules/nf-core/gunzip/main'
include { BEDTOOLS_MAKEWINDOWS              } from '../../modules/nf-core/bedtools/makewindows/main'

workflow GENOME_BINNING {

    take:
    genome_filtered_bed    // channel: [ bed ] - genome bed file already filtered for blacklist
    keep_regions_bed_param // val: params.keep_regions_bed parameter

    main:

    ch_versions = Channel.empty()

    //
    // Load keep_regions_bed as a channel
    //
    ch_keep_regions_bed = Channel.empty()

    if (keep_regions_bed_param) {
        if (keep_regions_bed_param.endsWith('.gz')) {
            ch_keep_regions_bed = GUNZIP_KEEP_REGIONS_BED(
                [ [:], keep_regions_bed_param ]
            ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_KEEP_REGIONS_BED.out.versions)
        } else {
            ch_keep_regions_bed = Channel.fromPath(keep_regions_bed_param, checkIfExists: true)
        }
    }

    //
    // Create bins for genome
    //
    BEDTOOLS_MAKEWINDOWS(
        genome_filtered_bed.map { bed -> tuple([id: bed.simpleName], bed) }
    )
    ch_versions = ch_versions.mix(BEDTOOLS_MAKEWINDOWS.out.versions)

    ch_binned_genome = BEDTOOLS_MAKEWINDOWS.out.bed

    emit:
    binned_genome    = ch_binned_genome      // channel: [ meta, bed ] - binned genome
    keep_regions_bed = ch_keep_regions_bed   // channel: [ bed ] - keep regions bed file
    versions         = ch_versions           // channel: [ versions.yml ]
}

//
// Run genome compartmentalization analysis
//

include { CHR_SPLIT                    } from '../../../modules/local/chr_split/main'
include { CHR_COMPARTMENTS_CALLING     } from '../../../modules/local/chr_compartments_calling/main'
include { CHR_COMBINE_COMPARTMENTS     } from '../../../modules/local/chr_combine_compartments/main'
include { BUILD_CONSENSUS              } from '../../../modules/local/build_consensus/main'

workflow COMPARTMENTALIZATION_ANALYSIS {

    take:
    bigwig_tracks    // channel: [ meta, bigwig ]
    binned_genome    // channel: path BED file
    chrom_sizes      // path: chromosome sizes file
    outdir           // val: output directory
    binsize          // val: bin size
    gtf              // path: GTF file

    main:

    ch_versions = Channel.empty()

    //
    // Extract track information from BigWig files
    //
    ch_compartmentTracks = bigwig_tracks
        .map { meta, bigwig ->
            [meta.experimentalID, meta.fraction, meta.sample_group, bigwig]
        }

    //
    // Create CSV file with all tracks
    //
    ch_compartments_csv = ch_compartmentTracks
        .map { experimentalID, fraction, sample_group, bigwig ->
            "${experimentalID},${fraction},${sample_group},${bigwig}\n"
        }
        .collectFile(
            name: 'compartments_tracks.csv',
            seed: "Patient_name,Fraction,Status,File\n",
            storeDir: "${outdir}/csv",
            sort: true
        )

    //
    // Get unique samples (experimentalID)
    //
    ch_uniqueSamples = ch_compartmentTracks
        .map { experimentalID, fraction, sample_group, bigwig -> experimentalID }
        .unique()

    //
    // Create meta for BED file
    //
    ch_bed_with_meta = binned_genome
        .map { bed ->
            [ [id: bed.baseName], bed ]
        }

    //
    // Split binned genome by chromosome
    //
    CHR_SPLIT(
        ch_bed_with_meta,
        chrom_sizes
    )

    //
    // Transpose to get one item per chromosome file
    //
    ch_chromBeds = CHR_SPLIT.out.beds
        .transpose()
        .map { meta, chr_bed ->
            def chr_name = chr_bed.baseName
            def new_meta = meta + [chromosome: chr_name]
            [ new_meta, chr_bed ]
        }

    //
    // Combine each chromosome with each unique sample and CSV
    //
    ch_chromSampleTuples = ch_chromBeds
        .combine(ch_uniqueSamples)
        .combine(ch_compartments_csv)
        .map { meta, chr_bed, patient, csv ->
            [ meta, chr_bed, patient, csv ]
        }

    //
    // Call compartments for each chromosome/sample combination
    //
    CHR_COMPARTMENTS_CALLING(
        ch_chromSampleTuples,
        binsize,
        gtf
    )
    ch_versions = ch_versions.mix(CHR_COMPARTMENTS_CALLING.out.versions.first())

    //
    // Group BED files by sample
    //
    ch_beds_by_sample = CHR_COMPARTMENTS_CALLING.out.bed_files
        .groupTuple(by: 0)

    //
    // Group BedGraph files by sample
    //
    
    ch_bedgraphs_by_sample = CHR_COMPARTMENTS_CALLING.out.bedgraph_files
        .groupTuple(by: 0)

    //
    // Combine BED and BedGraph channels for same sample
    //
    ch_combine_input = ch_beds_by_sample
        .join(ch_bedgraphs_by_sample)

    //
    // Merge all chromosomes per sample
    //
    CHR_COMBINE_COMPARTMENTS(
        ch_combine_input
    )
    ch_versions = ch_versions.mix(CHR_COMBINE_COMPARTMENTS.out.versions.first())

    //
    // Extract sample_group from original tracks
    //
    ch_sample_groups = ch_compartmentTracks
        .map { experimentalID, fraction, sample_group, bigwig ->
            [experimentalID, sample_group]
        }
        .unique()

    //
    // Add sample_group to combined beds
    //
    ch_beds_with_group = CHR_COMBINE_COMPARTMENTS.out.combined_beds
        .combine(ch_sample_groups)
        .filter { patient_bed, bed, patient_group, group ->
            patient_bed == patient_group
        }
        .map { patient_bed, bed, patient_group, group ->
            [group, patient_bed, bed]
        }

    //
    // Group combined BEDs by sample_group
    //
    ch_consensus_input = ch_beds_with_group
        .map { group, patient, bed -> [group, bed] }
        .groupTuple()

    //
    // Generate majority and strict consensus
    //
    BUILD_CONSENSUS(
        ch_consensus_input
    )

    emit:
    bed_files              = CHR_COMPARTMENTS_CALLING.out.bed_files          // channel: [ patient, bed ]
    bedgraph_files         = CHR_COMPARTMENTS_CALLING.out.bedgraph_files     // channel: [ patient, bedgraph ]
    combined_beds          = CHR_COMBINE_COMPARTMENTS.out.combined_beds      // channel: [ patient, bed ]
    combined_bedgraphs     = CHR_COMBINE_COMPARTMENTS.out.combined_bedgraphs // channel: [ patient, bedgraph ]
    consensus_majority     = BUILD_CONSENSUS.out.consensus_majority          // channel: [ path(bed) ]
    consensus_strict       = BUILD_CONSENSUS.out.consensus_strict            // channel: [ path(bed) ]
    versions               = ch_versions                                     // channel: [ versions.yml ]
}

//
// Run genome compartmentalization analysis
//

include { CHROMOSOME_SPLIT     } from '../../../modules/local/chromosome_split/main'
include { COMPARTMENTS_CALLING } from '../../../modules/local/compartments_calling/main'
include { COMBINE_COMPARTMENTS } from '../../../modules/local/combine_compartments/main'

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
    // MODULE: Split binned genome by chromosome
    //
    CHROMOSOME_SPLIT(
        ch_bed_with_meta,
        chrom_sizes
    )

    //
    // Transpose to get one item per chromosome file
    //
    ch_chromBeds = CHROMOSOME_SPLIT.out.beds
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
    COMPARTMENTS_CALLING(
        ch_chromSampleTuples,
        binsize,
        gtf
    )
    ch_versions = ch_versions.mix(COMPARTMENTS_CALLING.out.versions.first())

    //
    // Group BED files by sample
    //
    ch_beds_by_sample = COMPARTMENTS_CALLING.out.bed_files
        .map { patient, bed -> [patient, bed] }
        .groupTuple()

    //
    // Group BedGraph files by sample
    //
    ch_bedgraphs_by_sample = COMPARTMENTS_CALLING.out.bedgraph_files
        .map { patient, bedgraph -> [patient, bedgraph] }
        .groupTuple()

    //
    // Combine BED and BedGraph channels for same sample
    //
    ch_combine_input = ch_beds_by_sample
        .join(ch_bedgraphs_by_sample)

    //
    // Merge all chromosomes per sample
    //
    COMBINE_COMPARTMENTS(
        ch_combine_input
    )
    ch_versions = ch_versions.mix(COMBINE_COMPARTMENTS.out.versions.first())

    emit:
    bed_files            = COMPARTMENTS_CALLING.out.bed_files       // channel: [ patient, bed ]
    bedgraph_files       = COMPARTMENTS_CALLING.out.bedgraph_files  // channel: [ patient, bedgraph ]
    combined_beds        = COMBINE_COMPARTMENTS.out.combined_beds   // channel: [ patient, bed ]
    combined_bedgraphs   = COMBINE_COMPARTMENTS.out.combined_bedgraphs // channel: [ patient, bedgraph ]
    versions             = ch_versions                              // channel: [ versions.yml ]
}
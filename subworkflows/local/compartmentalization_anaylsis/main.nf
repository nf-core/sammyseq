//
// Run genome compartmentalization analysis
//

workflow COMPARTMENTALIZATION_ANALYSIS {

    take:
    bigwig_tracks    // channel: [ meta, bigwig ] from DEEPTOOLS_BAMCOVERAGE
    binned_genome    // channel: [ meta, bed ] from GENOME_BINNING
    outdir           // val: output directory

    main:
    
    ch_compartmentTracks = bigwig_tracks
        .map { meta, bigwig -> 
            [meta.experimentalID, meta.fraction, meta.sample_group, bigwig] 
        }

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
        .map { file -> [[ id:'compartment_analysis' ], file] }

    //
    // Get unique samples (experimentalID)
    //
    ch_uniqueSamples = ch_compartmentTracks
        .map { experimentalID, fraction, sample_group, bigwig -> experimentalID }
        .unique()

    //
    // Split binned genome by chromosome
    //
    ch_genomeBins = binned_genome
        .flatMap { meta, bedFile ->
            // Read BED file and group lines by chromosome
            def lines = bedFile.getText().split('\n').findAll { it.trim() }
            def grouped = lines.groupBy { it.split('\t')[0] }
            
            // Create a separate entry for each chromosome
            grouped.collect { chrom, chromLines ->
                tuple([id: meta.id, chromosome: chrom], chromLines)
            }
        }

    //
    // Combine each chromosome with each unique sample
    //
    ch_chromSampleTuples = ch_genomeBins
        .combine(ch_uniqueSamples)
        .map { meta, bedLines, patient ->
            tuple(meta, bedLines, patient)
        }

    emit:
    tracks_csv           = ch_compartments_csv      // channel: [ meta, csv ]
    chrom_sample_tuples  = ch_chromSampleTuples     // channel: [ meta, bedLines, patient ]
}
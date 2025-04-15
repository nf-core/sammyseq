/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap            } from 'plugin/nf-schema'
include { paramsSummaryMultiqc        } from '../subworkflows/nf-core/utils_nfcore_pipeline/main'
include { softwareVersionsToYAML      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText      } from '../subworkflows/local/utils_nfcore_sammyseq_pipeline'

include { UTILS_NEXTFLOW_PIPELINE     } from '../subworkflows/nf-core/utils_nextflow_pipeline/main'
include { UTILS_NFSCHEMA_PLUGIN       } from '../subworkflows/nf-core/utils_nfschema_plugin/main'
include { UTILS_NFCORE_PIPELINE       } from '../subworkflows/nf-core/utils_nfcore_pipeline/main'

include { CAT_FASTQ                   } from '../modules/nf-core/cat/fastq'
include { TRIMMOMATIC                 } from '../modules/nf-core/trimmomatic'
include { TRIMGALORE                  } from '../modules/nf-core/trimgalore/main'
include { DEEPTOOLS_BAMCOVERAGE       } from '../modules/nf-core/deeptools/bamcoverage'
include { BEDTOOLS_MAKEWINDOWS        } from '../modules/nf-core/bedtools/makewindows/main'

include { FASTQ_ALIGN_BWAALN          } from '../subworkflows/nf-core/fastq_align_bwaaln/main.nf'
include { FASTQ_ALIGN_DNA             } from '../subworkflows/nf-core/fastq_align_dna/main'
include { FASTQ_ALIGN_BOWTIE2         } from '../subworkflows/nf-core/fastq_align_bowtie2/main'
include { BAM_MARKDUPLICATES_PICARD   } from '../subworkflows/nf-core/bam_markduplicates_picard'

include { DEEPTOOLS_MULTIBAMSUMMARY   } from '../modules/nf-core/deeptools/multibamsummary/main'
include { DEEPTOOLS_PLOTCORRELATION   } from '../modules/nf-core/deeptools/plotcorrelation/main'
include { DEEPTOOLS_PLOTPCA           } from '../modules/nf-core/deeptools/plotpca/main'
include { DEEPTOOLS_PLOTFINGERPRINT as DEEPTOOLS_PLOTFINGERPRINT_GLOBAL } from '../modules/nf-core/deeptools/plotfingerprint/main'
include { DEEPTOOLS_PLOTFINGERPRINT as DEEPTOOLS_PLOTFINGERPRINT_REGION } from '../modules/nf-core/deeptools/plotfingerprint/main'

include { DEEPTOOLS_COMPUTEMATRIX     } from '../modules/nf-core/deeptools/computematrix/main'
include { DEEPTOOLS_PLOTPROFILE       } from '../modules/nf-core/deeptools/plotprofile/main'
include { DEEPTOOLS_PLOTHEATMAP       } from '../modules/nf-core/deeptools/plotheatmap/main'

// include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_FILTER     }   from '../modules/nf-core/samtools/view/main'
// include { SAMTOOLS_SORT as SAMTOOLS_SORT_FILTERED   }   from '../modules/nf-core/samtools/sort/main'
// include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FILTERED }   from '../modules/nf-core/samtools/index/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPARE_GENOME            } from '../subworkflows/local/prepare_genome'
include { CAT_FRACTIONS             } from '../subworkflows/local/cat_fractions'
include { RTWOSAMPLESMLE            } from '../modules/local/rtwosamplesmle'
include { CALL_COMPARTMENTS         } from '../modules/local/call_compartments'
include { FILTER_BAM_SAMTOOLS       } from '../subworkflows/local/filter_bam_samtools'
include { BIGWIG_PLOT_DEEPTOOLS     } from '../subworkflows/local/bigwig_plot_deeptools'
include { DEEPTOOLS_QC              } from '../subworkflows/local/deeptools_qc'
include { MERGE_COMPARTMENTS        } from '../modules/local/merge_compartments'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.fasta) { ch_fasta =  Channel.fromPath(params.fasta) } else { exit 1, 'Fasta reference genome not specified!' }

// Modify fasta channel to include meta data
ch_fasta_meta = ch_fasta.map{ it -> [[id:it[0].baseName], it] }.collect()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SAMMYSEQ {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //

    PREPARE_GENOME (params.fasta,
                    params.aligner,
                    params.bwa_index,
                    params.bowtie2_index,
                    params.blacklist,
                    params.chrom_sizes,
                    params.fai,
                    params.binsize,
                    params.gtf,
                    params.gene_bed,
                    params.keep_regions_bed
                    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    if (params.stopAt == 'PREPARE_GENOME') {
        return
    }

    //
    // Branch channels from input samplesheet channel
    //
    ch_samplesheet
        .branch { meta, fastqs ->
            single  : fastqs.size() == 1
                return [ meta, fastqs.flatten() ]
            multiple: fastqs.size() > 1
                return [ meta, fastqs.flatten() ]
        }
        .set { ch_merge_lane }

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (ch_merge_lane.multiple)
        .reads
        .mix(ch_merge_lane.single)
        .set {ch_starter}
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())

    if (params.stopAt == 'CAT_FASTQ_lane') {
        return
    }

    //
    // Combine fractions by expID
    //

    if(params.combine_fractions){
        merged_reads = CAT_FRACTIONS(//INPUT_CHECK.out.reads_to_merge,
                                    //INPUT_CHECK.out.reads
                                    ch_starter
                                    )//.out.merged_reads
    } else {
        //merged_reads = INPUT_CHECK.out.reads
        merged_reads = ch_starter
    }

    if (params.stopAt == 'CAT_FRACTIONS') {
        return
    }

    ///
    //  TRIMMING!
    //

    ch_trimmed= Channel.empty()
    if (params.trimmer == 'trimmomatic') {
        TRIMMOMATIC(merged_reads)
        ch_trimmed=TRIMMOMATIC.out.trimmed_reads
        ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions)
    } else if (params.trimmer == 'trimgalore') {
        TRIMGALORE(merged_reads)
        ch_trimmed=TRIMGALORE.out.reads
        ch_versions = ch_versions.mix(TRIMGALORE.out.versions)
    }

    ch_fastqc_trim = ch_trimmed
                    .map{ meta, path ->
                    def id=meta.subMap('id')
                    newid=id.id + "_trim"
                    sng=meta.subMap('single_end').single_end
                    newmeta=[id: newid, single_end: sng]
                    [ newmeta ,path]
                }

    //
    // MODULE: Run FastQC
    //
    // a channel is created for the trimmed files and the id is renamed to meta, so that when passed to fastqc it does not overwrite the output files with non-trimmed ones, trimmed and untrimmed fastq channels are merged and the resulting channel is passed to FASTQC
    ch_fastqc_in = ch_fastqc_trim.mix(merged_reads)
    //ch_fastqc_in.view()
    FASTQC (
        ch_fastqc_in
        //merged_reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    if (params.stopAt == 'TRIMMING') {
        return
    }

    ch_aligned_bam = Channel.empty()
    if (params.aligner == 'bwaaln') {
        FASTQ_ALIGN_BWAALN(
            ch_trimmed,
            PREPARE_GENOME.out.bwa_index
        )
        ch_aligned_bam = FASTQ_ALIGN_BWAALN.out.bam
        ch_versions = ch_versions.mix(FASTQ_ALIGN_BWAALN.out.versions)
    } else if (params.aligner == 'bwamem') {
        FASTQ_ALIGN_DNA(
            ch_trimmed,
            PREPARE_GENOME.out.bwa_index,
            ch_fasta_meta,
            params.aligner,
            true
        )
        ch_aligned_bam = FASTQ_ALIGN_DNA.out.bam
        ch_versions = ch_versions.mix(FASTQ_ALIGN_DNA.out.versions)
    } else if (params.aligner == 'bowtie2') {
        FASTQ_ALIGN_BOWTIE2(
            ch_trimmed,
            PREPARE_GENOME.out.bowtie2_index,
            params.save_unaligned,
            false,
            PREPARE_GENOME.out.fasta.map{ fasta_ -> [ [ id:'fasta' ], fasta_ ] }
        )
        ch_aligned_bam = FASTQ_ALIGN_BOWTIE2.out.bam
        ch_versions = ch_versions.mix(FASTQ_ALIGN_BOWTIE2.out.versions)
    }

if (params.stopAt == 'ALIGNMENT') {
    return
}

    ///
    //  DUPLICATE READS REMOVAL
    //

    // MARK DUPLICATES IN BAM FILE
    BAM_MARKDUPLICATES_PICARD (
        ch_aligned_bam,
        ch_fasta_meta,
        PREPARE_GENOME.out.fai.collect { [ [id: 'fasta'], it ] }
        )
    ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_PICARD.out.versions)

    ch_mle_in = BAM_MARKDUPLICATES_PICARD.out.bam

    if (params.stopAt == 'BAM_MARKDUPLICATES_PICARD') {
        return
    }

    ch_bam_bai_combined =  BAM_MARKDUPLICATES_PICARD.out.bam
        .join(BAM_MARKDUPLICATES_PICARD.out.bai, by: [0], remainder: true)
        .map {
            meta, bam, bai  ->
                    [ meta, bam, bai ]

        }

    FILTER_BAM_SAMTOOLS(
        ch_bam_bai_combined,
        ch_fasta_meta
    )

//ch_versions = ch_versions.mix(FILTER_BAM_SAMTOOLS.out.versions)

    ch_bam_bai_filtered = FILTER_BAM_SAMTOOLS.out.bam
        .join(FILTER_BAM_SAMTOOLS.out.bai, by: [0], remainder: true)
        .map {
            meta, bam, bai  ->
                    [ meta, bam, bai ]

        }

    ch_fai_path = PREPARE_GENOME.out.fai.map { it[1] }
    //ch_fai_path.view()
    ch_fasta_path = ch_fasta_meta.map { it[1] }
    //ch_fasta_path.view()

    DEEPTOOLS_BAMCOVERAGE (
        ch_bam_bai_filtered,
        ch_fasta_path,
        ch_fai_path,
        params.blacklist ? PREPARE_GENOME.out.blacklist : []
    )

    ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE.out.versions)

    if (params.stopAt == 'DEEPTOOLS_BAMCOVERAGE') {
        return
    }

    DEEPTOOLS_QC (
        FILTER_BAM_SAMTOOLS.out.bam,
        FILTER_BAM_SAMTOOLS.out.bai,
        params.corr_method,
        params.blacklist ? PREPARE_GENOME.out.blacklist : []
    )
    ch_dt_corrmatrix     = DEEPTOOLS_QC.out.correlation_matrix
    ch_dt_pcadata        = DEEPTOOLS_QC.out.pca_data
    ch_dt_fpmatrix_global = DEEPTOOLS_QC.out.fingerprint_matrix_global
    ch_dt_fpmetrics_global = DEEPTOOLS_QC.out.fingerprint_metrics_global
    if (params.region) {
        ch_dt_fpmatrix_region = DEEPTOOLS_QC.out.fingerprint_matrix_region
        ch_dt_fpmetrics_region = DEEPTOOLS_QC.out.fingerprint_metrics_region
    }
    ch_versions = ch_versions.mix(DEEPTOOLS_QC.out.versions)

    if (params.comparisonFile) {
        // Add the suffix "_T1" to each sample ID in the comparison file

        ch_bam_input=FILTER_BAM_SAMTOOLS.out.bam

        ch_bam_input.view()

        // 1. Create a Comparisons Channels (one for sample 1 in comparison and another for sample 2 in comparison)

        Channel
            .fromPath(params.comparisonFile)
            .splitCsv(header : true)
            .map{ row ->
                //[ row.sample1 + "_T1", row.sample2 + "_T1",row.sample1 + "_T1_VS_" + row.sample2 + "_T1"]
                //[ row.sample1 + "_T1", row.sample1 + "_T1_VS_" + row.sample2 + "_T1"]
                [ row.sample1 , row.sample1 + "_VS_" + row.sample2 ]
                }
                .set { comparisons_ch_s1 }


        Channel.fromPath(params.comparisonFile)
                .splitCsv(header : true)
                .map{ row ->
                    // [ row.sample2 + "_T1", row.sample1 + "_T1",row.sample1 + "_T1_VS_" + row.sample2 + "_T1"]
                    //[ row.sample2 + "_T1", row.sample1 + "_T1_VS_" + row.sample2 + "_T1"]
                    [ row.sample2 , row.sample1 + "_VS_" + row.sample2 ]
                }
                .set { comparisons_ch_s2 }

        //2. convert bam file to input
        // [[id:ggg, paired:true],path.bam]
        ch_bam_input
                .map {meta, bam ->
                    id=meta.subMap('id')
                    [id.id, bam]

                    }
                .set { ch_bam_reformat }

        //3. combine comparison channel with bam list channel
        comparisons_ch_s1
                .combine(ch_bam_reformat , by:0)
                .map { sample1, comparison, bam ->
                    //[ comparison:comparison, sample1:sample1, sample2:sample2, bam1:bam1]
                [ comparison, bam]
                }
                .set { bam1_comparison }

                        //.view{"join= ${it}"}
        comparisons_ch_s2
                .combine(ch_bam_reformat, by:0)
                .map{ sample2, comparison, bam ->
                    //[ comparison:comparison , sample1:sample1, sample2:sample2 ,bam2:bam2 ]
                    [ comparison, bam]
                    }
                .set{ bam2_comparison }

        bam1_comparison
                .join(bam2_comparison, remainder: false, by: 0 )
                .set{comparisons_merge_ch}

        //comparisons_merge_ch
        //        .view{ "comparisons_merge_ch: ${it}" }

        //4.run rscript

        RTWOSAMPLESMLE (comparisons_merge_ch,
                        PREPARE_GENOME.out.chrom_sizes
                        )

        if (params.stopAt == 'RTWOSAMPLESMLE') {
        return
        }
    }

    ///
    /// Compartments calling
    ///

    if (params.compartmentsAnalysis) {

        ch_compartmentTracks = DEEPTOOLS_BAMCOVERAGE.out[0]
            .map { meta, bigwig -> [meta.experimentalID, meta.fraction, meta.sample_group, bigwig] }

        ch_compartmentsTSV = Channel.of(["Patient_name", "Fraction", "Status", "File"])
            .concat(ch_compartmentTracks)
            .map { it.join("\t") }
            .collect()

        ch_binsize = Channel.value(params.binsize)
        ch_uniqueSamples = ch_compartmentTracks.map { it[0] }.unique()

        def validChroms = []
        if (params.keep_regions_bed) {
             validChroms = file(params.keep_regions_bed)
                .readLines()
                .findAll { it }
                .collect { it.tokenize()[0].trim() }
                .unique()
        }

        ch_genomeBins = PREPARE_GENOME.out.binned_genome
            .flatMap { meta, bedFile ->
                def lines = bedFile.text.split('\n').findAll { it }
                def filtered = validChroms ? lines.findAll { validChroms.contains(it.split('\t')[0]) } : lines
                def grouped = filtered.groupBy { it.split('\t')[0] }

                grouped.collect { chrom, chromLines ->
                    tuple([id: meta.id, chromosome: chrom], chromLines)
                }
            }

        ch_chromSampleTuples = ch_genomeBins
            .combine(ch_uniqueSamples)
            .map { meta, bedLines, patient ->
                tuple(meta, bedLines, patient)
            }

        CALL_COMPARTMENTS(
            ch_compartmentsTSV,
            ch_binsize,
            PREPARE_GENOME.out.gtf,
            ch_chromSampleTuples
        )

        ch_compartments_beds = CALL_COMPARTMENTS.out.bed_files
            .map { file -> tuple(file.getBaseName().tokenize('_')[0..1].join('_'), file) }
            .groupTuple()

        ch_compartments_eigen_bedgraphs = CALL_COMPARTMENTS.out.bedgraph_files
            .map { file -> tuple(file.getBaseName().tokenize('_')[0..1].join('_'), file) }
            .groupTuple()

        ch_merged_compartments = ch_compartments_beds
            .join(ch_compartments_eigen_bedgraphs)
            .map { sample_id, bed_files, bedgraph_files ->
                tuple(sample_id, bed_files, bedgraph_files)
            }

        MERGE_COMPARTMENTS(ch_merged_compartments)
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'sammyseq_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    //ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.metrics.collect{it[1]}.ifEmpty([]))
    //ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.flagstat.collect{it[1]}.ifEmpty([]))
    //ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.idxstats.collect{it[1]}.ifEmpty([]))

    ch_multiqc_files = ch_multiqc_files.mix(FILTER_BAM_SAMTOOLS.out.stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FILTER_BAM_SAMTOOLS.out.flagstat.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FILTER_BAM_SAMTOOLS.out.idxstats.collect{it[1]}.ifEmpty([]))

    ch_multiqc_files = ch_multiqc_files.mix(DEEPTOOLS_QC.out.correlation_matrix.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(DEEPTOOLS_QC.out.pca_data.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(DEEPTOOLS_QC.out.correlation_matrix.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(DEEPTOOLS_QC.out.pca_data.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(DEEPTOOLS_QC.out.fingerprint_matrix_global.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(DEEPTOOLS_QC.out.fingerprint_metrics_global.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

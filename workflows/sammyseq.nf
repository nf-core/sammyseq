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
include { BAM_MARKDUPLICATES_PICARD   } from '../subworkflows/nf-core/bam_markduplicates_picard'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPARE_GENOME                } from '../subworkflows/local/prepare_genome'
include { GENOME_BINNING                } from '../subworkflows/local/genome_binning'
include { CAT_FRACTIONS                 } from '../subworkflows/local/cat_fractions'
include { FILTER_BAM_SAMTOOLS           } from '../subworkflows/local/filter_bam_samtools'
include { BIGWIG_PLOT_DEEPTOOLS         } from '../subworkflows/local/bigwig_plot_deeptools'
include { DEEPTOOLS_QC                  } from '../subworkflows/local/deeptools_qc'
include { GENERATE_COMPARISONS_MLE      } from '../subworkflows/local/generate_comparisons_mle'
include { GENERATE_COMPARISONS_SAMPLESHEET    } from '../subworkflows/local/generate_comparisons_samplesheet'
include { GENERATE_COMPARISONS_BIGWIG   } from '../subworkflows/local/generate_comparisons_bigwig'

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
                    params.tss_bed)

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
        FASTQ_ALIGN_DNA(
            ch_trimmed,
            PREPARE_GENOME.out.bowtie2_index,
            ch_fasta_meta,
            params.aligner,
            true
        )
        ch_aligned_bam = FASTQ_ALIGN_DNA.out.bam
        ch_versions = ch_versions.mix(FASTQ_ALIGN_DNA.out.versions)
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
    DEEPTOOLS_BAMCOVERAGE.out.bigwig,
    params.corr_method,
    params.blacklist ? PREPARE_GENOME.out.blacklist : []
    )
    ch_dt_corrmatrix     = DEEPTOOLS_QC.out.correlation_matrix
    ch_dt_pcadata        = DEEPTOOLS_QC.out.pca_data
    if (params.plotfingerprint) {
        ch_dt_fpmatrix_global = DEEPTOOLS_QC.out.fingerprint_matrix_global
        ch_dt_fpmetrics_global = DEEPTOOLS_QC.out.fingerprint_metrics_global
        if (params.region) {
            ch_dt_fpmatrix_region = DEEPTOOLS_QC.out.fingerprint_matrix_region
            ch_dt_fpmetrics_region = DEEPTOOLS_QC.out.fingerprint_metrics_region
        }
    }
    ch_versions = ch_versions.mix(DEEPTOOLS_QC.out.versions)

        if (params.tss_bed) {
        ch_bw_by_sample = DEEPTOOLS_BAMCOVERAGE.out.bigwig
            .map { meta, bw -> tuple(meta.experimentalID, [meta, bw]) }
            .groupTuple()
            .map { id, values ->
                def sorted = values.sort { it[1].getBaseName() }
                def sorted_bw = sorted.collect { it[1] }
                def sorted_labels = sorted.collect { it[0].id }
                tuple([id: id, labels: sorted_labels], sorted_bw)
    }

    BIGWIG_PLOT_DEEPTOOLS(
        ch_bw_by_sample,
        PREPARE_GENOME.out.tss_bed
    )

        ch_versions = ch_versions.mix(BIGWIG_PLOT_DEEPTOOLS.out.versions)
    }

    //
    // GENOME BINNING: Run only if comparisonFile or comparison is provided
    //

    if (params.comparisonFile && params.comparison) {
        error "Cannot specify both --comparisonFile and --comparison parameters. Please use only one method."
    }

    if (params.comparisonFile || params.comparison) {
        GENOME_BINNING(
            PREPARE_GENOME.out.filtered_bed,
            params.keep_regions_bed
        )
        ch_versions = ch_versions.mix(GENOME_BINNING.out.versions)
    }


    //
    // Generate comparisons
    //
    if (params.comparisonFile || params.comparison) {

        ch_comparison_results = Channel.empty()

        if (params.comparison_maker == 'spp') {
            GENERATE_COMPARISONS_MLE(
                FILTER_BAM_SAMTOOLS.out.bam,
                ch_samplesheet,
                PREPARE_GENOME.out.chrom_sizes
            )
            ch_comparison_results = GENERATE_COMPARISONS_MLE.out.mle_results

        } else if (params.comparison_maker == 'bigwigcompare') {
            GENERATE_COMPARISONS_BIGWIG(
                DEEPTOOLS_BAMCOVERAGE.out.bigwig,
                ch_samplesheet,
            )
            ch_comparison_results = GENERATE_COMPARISONS_BIGWIG.out.comparison_results

        }

        //
        // Generate ratio samplesheet CSV - common for both tools
        //
        GENERATE_COMPARISONS_SAMPLESHEET(
            ch_comparison_results,
            params.outdir
        )
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
    if (!params.skip_multiqc) {
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

    if (params.plotfingerprint) {
        ch_multiqc_files = ch_multiqc_files.mix(DEEPTOOLS_QC.out.fingerprint_matrix_global.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(DEEPTOOLS_QC.out.fingerprint_metrics_global.collect{it[1]}.ifEmpty([]))
    }

    if (params.tss_bed) {
        ch_multiqc_files  = ch_multiqc_files.mix(BIGWIG_PLOT_DEEPTOOLS.out.plotprofile_pdf.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files  = ch_multiqc_files.mix(BIGWIG_PLOT_DEEPTOOLS.out.plotprofile_table.collect{it[1]}.ifEmpty([]))
    }

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )
    }
    emit:
    multiqc_report = !params.skip_multiqc ? MULTIQC.out.report.toList() : Channel.empty() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

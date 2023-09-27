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
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowSammyseq.initialise(params, log)

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
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { PREPARE_GENOME      } from '../subworkflows/local/prepare_genome'
include { FASTQ_ALIGN_BWAALN  } from '../subworkflows/nf-core/fastq_align_bwaaln/main.nf'
include { BAM_MARKDUPLICATES_PICARD } from '../subworkflows/nf-core/bam_markduplicates_picard'
include { CUT_SIZES_GENOME } from "../modules/local/chromosomes_size"
include { RTWOSAMPLESMLE } from '../modules/local/rtwosamplesmle'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { CAT_FASTQ                   } from '../modules/nf-core/cat/fastq'
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { TRIMMOMATIC                 } from '../modules/nf-core/trimmomatic'
include { SAMTOOLS_FAIDX              } from '../modules/nf-core/samtools/faidx'
include { DEEPTOOLS_BAMCOVERAGE       } from '../modules/nf-core/deeptools/bamcoverage'
 
// include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_FILTER     }   from '../modules/nf-core/samtools/view/main'
// include { SAMTOOLS_SORT as SAMTOOLS_SORT_FILTERED   }   from '../modules/nf-core/samtools/sort/main'
// include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FILTERED }   from '../modules/nf-core/samtools/index/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow SAMMYSEQ {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //

    if (params.stopAt == 'BEGIN') {
        return
    }


    PREPARE_GENOME (params.fasta,
                    params.bwa_index)

    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    if (params.stopAt == 'PREPARE_GENOME') {
        return
    }

    INPUT_CHECK (
        file(params.input)
    )

    //INPUT_CHECK.out.reads.view()

    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    // extract fastq to merge by expID
    INPUT_CHECK.out.reads_to_merge
            .filter { meta, fastq -> fastq.size() > 1 }    
            .set { ch_reads_to_process_in_CAT_FASTQ }

    ch_reads_to_process_in_CAT_FASTQ    
            .flatMap { meta, fastq  -> 
                    meta.collect { m -> [m,fastq]}
                }
            .set { ch_to_CAT }

    CAT_FASTQ (
        ch_to_CAT
    ).reads.set { cat_fastq_output }        

    ch_cat_adjusted = CAT_FASTQ.out.reads.map { meta, fastq ->
        return [meta, [fastq]]
    } // make fastq of CAT_FASTQ a list of paths
    merged_reads = INPUT_CHECK.out.reads.mix(ch_cat_adjusted)

    //
    // MODULE: Run FastQC
    //
    //merged_reads.view()


    FASTQC (
        merged_reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    TRIMMOMATIC (
       merged_reads 
    )

    if (params.stopAt == 'TRIMMOMATIC') {
        return
    }


    //TRIMMOMATIC.out.trimmed_reads.view()
    
    FASTQ_ALIGN_BWAALN (
        TRIMMOMATIC.out.trimmed_reads,
        PREPARE_GENOME.out.bwa_index
    )

    if (params.stopAt == 'FASTQ_ALIGN_BWAALN') {
        return
    }


    // PICARD MARK_DUPLICATES
    // Index Fasta File for Markduplicates
    SAMTOOLS_FAIDX (
            ch_fasta_meta,
             [[], []]
        )

    ch_fai_for_cut = SAMTOOLS_FAIDX.out.fai.collect()

    CUT_SIZES_GENOME(ch_fai_for_cut)

    //FASTQ_ALIGN_BWAALN.out.bam.view()

    // MARK DUPLICATES IN BAM FILE
    BAM_MARKDUPLICATES_PICARD (
        FASTQ_ALIGN_BWAALN.out.bam,
        ch_fasta_meta,
        SAMTOOLS_FAIDX.out.fai.collect()
        )

    if (params.stopAt == 'BAM_MARKDUPLICATES_PICARD') {
        return
    }
    // SAMTOOLS_VIEW_FILTER (
    //                 ch_bam_sorted.join(ch_bam_sorted_bai),
    //                 ch_fasta_meta,
    //                 []
    //             )
    // ch_versions = ch_versions.mix(SAMTOOLS_VIEW_FILTER.out.versions)

    //ch_bam_from_markduplicates = BAM_MARKDUPLICATES_PICARD.bam

    //BAM_MARKDUPLICATES_PICARD.out.bam.view()
    //BAM_MARKDUPLICATES_PICARD.out.bai.view()

    //ch_bam_bai_combined = BAM_MARKDUPLICATES_PICARD.out.bam.join(BAM_MARKDUPLICATES_PICARD.out.bai, by: [0])

   ch_bam_bai_combined = BAM_MARKDUPLICATES_PICARD.out.bam
        .join(BAM_MARKDUPLICATES_PICARD.out.bai, by: [0], remainder: true)
        .map {
            meta, bam, bai  ->     
                    [ meta, bam, bai ]
             
        }

    //ch_bam_bai_combined.view()
    //ch_fasta_meta.map { it[2] }.view
    //ch_fasta_meta.collect().view()
    //SAMTOOLS_FAIDX.out.fai.collect().view()
    //SAMTOOLS_FAIDX.out.fai.map { it['path'] }.collect().view()
    ch_fai_path = SAMTOOLS_FAIDX.out.fai.map { it[1] }
    //ch_fai_path.view()
    ch_fasta_path = ch_fasta_meta.map { it[1] }
    //ch_fasta_path.view()

    DEEPTOOLS_BAMCOVERAGE (
        ch_bam_bai_combined,
        ch_fasta_path,
        ch_fai_path
    )

    if (params.stopAt == 'DEEPTOOLS_BAMCOVERAGE') {
        return
    }

    if (params.comparisonFile) {

        //BAM_MARKDUPLICATES_PICARD.out.bam.view()

        // Declare the dictionary that will store the comparisons
        def comparisons = [:]

        // Read the CSV file
        csvFile = file(params.comparisonFile)
        csvFile.eachLine { line ->
        // Skip the first line (header)
            if (line.startsWith("sample1")) {
                return
            }
    
        // Extract the samples and add the "_T1" suffix
                def (sample1, sample2) = line.split(',')
                comparisons[sample1.trim() + "_T1"] = sample2.trim() + "_T1"
                println("Comparisons Map: $comparisons")
        }

        //println("Comparisons Map: $comparisons")  // Print the comparison map

        // Filter BAMs using the updated comparisons map
        

        ch_filtered_bams = BAM_MARKDUPLICATES_PICARD.out.bam
           .filter { meta, bam -> 
                boolean isPresent = comparisons.keySet().contains(meta.id) || comparisons.values().contains(meta.id)
                if (!isPresent) {
                    //println("BAM not in comparisons list: $bam")  // Print BAMs not in the comparison list
                }
                return isPresent
            }

        // Print filtered BAMs
        //ch_filtered_bams.view()

        // Create BAM pairs for comparison
        
        ch_filtered_bams
            .groupTuple(by: [0]) // Group by the metadata
            .map { meta, bams ->
                def comparisonValue = comparisons[meta.id]
                println("Comparison Value for ${meta.id}: $comparisonValue")
            }
            .filter { it != null } // Remove any null values from the channel
            .set { ch_bam_pairs_for_comparison }
        println("Before view")

        ch_bam_pairs_for_comparison.count().view()

        ch_bam_pairs_for_comparison.view()

        println("after view")

        if (params.stopAt == 'RTWOSAMPLESMLE') {
        return
        }
    }



    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowSammyseq.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowSammyseq.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

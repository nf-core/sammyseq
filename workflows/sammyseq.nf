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
include { CAT_FRACTIONS } from '../subworkflows/local/cat_fractions'
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
                    params.bwa)

    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    if (params.stopAt == 'PREPARE_GENOME') {
        return
    }

    INPUT_CHECK (
        file(params.input)
    )

    //INPUT_CHECK.out.reads.view()

    //TODO use branch

    ch_notmerge_lane = INPUT_CHECK.out.reads
                     .map{ meta, path -> 
                        id=meta.subMap('id')
                        meta=meta
                        path=path
                        [id.id, meta, path]
                      }
                     .groupTuple()
                     .filter{ it[1].size() == 1 }
                     .map{ id, meta, path -> 
                        meta_notmerge=meta[0]
                        path_notmerge=path[0]
                        [meta_notmerge, path_notmerge]
                     }

    //INPUT_CHECK.out.reads.view{"INPUT_CHECK.out.reads : ${it}"}
    //ch_notmerge_lane.view{"ch_notmerge_lane: ${it}"}

    ch_merge_lane = INPUT_CHECK.out.reads
                     .map{ meta, path -> 
                        id=meta.subMap('id')
                        meta=meta
                        path=path
                        [id.id, meta, path]
                      }
                     .groupTuple()
                     .filter{ it[1].size() >= 2 } //filtra per numero di meta presenti dopo il tupla se hai due meta vuol dire che devi unire due campioni 
                     .map{ id, meta, path -> 
                        single = meta[0].subMap('single_end')
                        meta = meta[0]
                        def flatPath = path.flatten()
                        [ meta , flatPath ]
                      }

    //ch_merge_lane.view{"ch_merge_lane ${it}"}

    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    if (params.stopAt == 'INPUT_CHECK') {
        return
    }

    //ch_merge_lane.view{"ch_merge_lane : ${it}"}

    CAT_FASTQ (
           ch_merge_lane
    ).reads.set { cat_lane_output }

    //cat_lane_output.view()     
    ch_starter = cat_lane_output.mix(ch_notmerge_lane)

    //ch_starter.view{"ch_starter : ${it}"}

    if (params.stopAt == 'CAT_FASTQ_lane') {
        return
    }


    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    // extract fastq to merge by expID


    if(params.combine_fractions){

        merged_reads = CAT_FRACTIONS(//INPUT_CHECK.out.reads_to_merge,
                                    //INPUT_CHECK.out.reads
                                    ch_starter
                                    )//.out.merged_reads

    } else {
        //merged_reads = INPUT_CHECK.out.reads
        merged_reads = ch_starter
    }

    //merged_reads.view{"merged_reads: ${it}"}

    if (params.stopAt == 'CAT_FRACTIONS') {
        return
    }

    //
    // MODULE: Run FastQC
    //

    TRIMMOMATIC (
       merged_reads 
    )

    //TRIMMOMATIC.out.trimmed_reads.view()


    //a channel is created for the trimmed files and the id is renamed to meta, so that when passed to fastqc it does not overwrite the output files with non-trimmed ones
    ch_fastqc_trim=TRIMMOMATIC.out.trimmed_reads
                   .map{ meta, path ->
                        def id=meta.subMap('id')
                        newid=id.id + "_trim"
                        sng=meta.subMap('single_end').single_end
                        newmeta=[id: newid, single_end: sng]
                        [ newmeta ,path]
                    } 

    //trimmed and untrimmed fastq channels are merged and the resulting channel is passed to FASTQC
    ch_fastqc_in = ch_fastqc_trim.mix(merged_reads)
    //ch_fastqc_in.view()
    FASTQC (
        ch_fastqc_in
        //merged_reads
    )

    ch_versions = ch_versions.mix(FASTQC.out.versions.first())



    if (params.stopAt == 'TRIMMOMATIC') {
        return
    }

    
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

    if (params.stopAt == 'SAMTOOLS_FAIDX') {
        return
    }

    CUT_SIZES_GENOME(ch_fai_for_cut)
    //CUT_SIZES_GENOME.out.ch_sizes_genome.view()

    //FASTQ_ALIGN_BWAALN.out.bam.view()

    // MARK DUPLICATES IN BAM FILE
    BAM_MARKDUPLICATES_PICARD (
        FASTQ_ALIGN_BWAALN.out.bam,
        ch_fasta_meta,
        SAMTOOLS_FAIDX.out.fai.collect()
        )


    //BAM_MARKDUPLICATES_PICARD.out.bam.view()
       
    ch_mle_in = BAM_MARKDUPLICATES_PICARD.out.bam
    //ch_mle_in.view()

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

   ch_bam_bai_combined =  BAM_MARKDUPLICATES_PICARD.out.bam
        .join(BAM_MARKDUPLICATES_PICARD.out.bai, by: [0], remainder: true)
        .map {
            meta, bam, bai  ->     
                    [ meta, bam, bai ]
             
        }


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
        // Add the suffix "_T1" to each sample ID in the comparison file

        ch_bam_input=BAM_MARKDUPLICATES_PICARD.out.bam

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
                       CUT_SIZES_GENOME.out.ch_sizes_genome
                       )

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

    ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.metrics.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.flagstat.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.idxstats.collect{it[1]}.ifEmpty([]))

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
    NfcoreTemplate.dump_parameters(workflow, params)
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

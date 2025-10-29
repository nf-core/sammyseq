// Import SAMTOOLS_FAIDX module
include { CAT_FASTQ                   } from '../../../modules/nf-core/cat/fastq/main'

// Definition of the subworkflow
workflow CAT_FRACTIONS {
    take:
    //reads_to_merge
    //all_reads
    starting_files

    main:

        all_files = starting_files.map { meta, fastq ->
            def expID = meta.expID  // take expID from meta
            meta.remove('fraction') // remove 'fraction' field
            [expID, meta, fastq]
        }

        all_files.groupTuple( by:0 )
            .map{ expID, meta , fastq ->
                meta2 = meta[0].clone()  // take first map only and create a copy
                meta2.id = expID         // substitute id value with expID
                [ [meta2], fastq.flatten().toList() ]
            }
            .set { reads_to_merge }

            reads_to_merge
            .filter { meta, fastq -> fastq.size() > 1 }
            .set { ch_reads_to_process_in_CAT_FASTQ }

            ch_reads_to_process_in_CAT_FASTQ
            .flatMap { meta, fastq  ->
                    meta.collect { m -> [m,fastq]}
                }
            .set { ch_to_CAT }

            //ch_to_CAT.view{"ch_to_CAT : ${it}"}

            CAT_FASTQ (
                    ch_to_CAT
            ).reads.set { cat_fastq_output }

            ch_cat_adjusted = CAT_FASTQ.out.reads.map { meta, fastq ->
                    return [meta, fastq]
            } // make fastq of CAT_FASTQ a list of paths
            merged_reads = starting_files.mix(ch_cat_adjusted)

    emit:
    merged_reads = merged_reads
}
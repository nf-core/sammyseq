//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }
        .set { reads }

    reads
        
        .map { meta, fastq -> 
            def expID = meta.expID  // Prendi l'expID da meta
            meta.remove('fraction') // Rimuovi il campo 'fraction'

            [expID, meta, fastq]
        }
        .groupTuple( by:0 )
        //.view()
        .map{ expID, meta , fastq ->
            meta = meta[0].clone()  // Prendi solo la prima mappa e crea una copia
            meta.id = expID  // Sostituisci il valore nel campo id con expID 
            [ [meta], fastq.flatten() ]
            }
        .set { reads_to_merge }


    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    reads_to_merge
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.single_end = row.single_end.toBoolean()
    meta.expID = row.experimentalID
    meta.fraction = row.fraction

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return fastq_meta
}

//
// Uncompress and prepare reference genome files
//

include {
    GUNZIP as GUNZIP_FASTA
    GUNZIP as GUNZIP_GTF
    GUNZIP as GUNZIP_GFF
    GUNZIP as GUNZIP_GENE_BED
    GUNZIP as GUNZIP_BLACKLIST } from '../../modules/nf-core/gunzip/main'

include {
    UNTAR as UNTAR_BWA_INDEX
    UNTAR as UNTAR_BOWTIE2_INDEX
    UNTAR as UNTAR_CHROMAP_INDEX
    UNTAR as UNTAR_STAR_INDEX    } from '../../modules/nf-core/untar/main'

//include { GFFREAD              } from '../../modules/nf-core/gffread/main'
include { BWA_INDEX            } from '../../modules/nf-core/bwa/index/main'
//include { SAMTOOLS_FAIDX } from '../../modules/nf-core/samtools/faidx'

workflow PREPARE_GENOME {
    take:
    input_fasta
    input_bwaindex

    main:

    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    ch_fasta = Channel.empty()
    if (input_fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA ( [ [:], params.fasta ] ).gunzip.map{ it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta =  [ [:], file(input_fasta) ] 
    }

    //println(ch_fasta)
    // Make fasta file available if reference saved or IGV is run
    //if (params.save_reference || !params.skip_igv) {

    // if (params.save_reference) {    
    //     file("${params.outdir}/genome/").mkdirs()
    //     ch_fasta.copyTo("${params.outdir}/genome/")
    // }


    //
    // Prepare genome intervals for filtering by removing regions in blacklist file
    //
    ch_genome_filtered_bed = Channel.empty()

    //
    // Uncompress BWA index or generate from scratch if required
    //
    ch_bwa_index = Channel.empty()
    
        if (params.bwa_index) {
            if (params.bwa_index.endsWith('.tar.gz')) {
                ch_bwa_index = UNTAR_BWA_INDEX ( [ [:], params.bwa_index ] ).untar
                ch_versions  = ch_versions.mix(UNTAR_BWA_INDEX.out.versions)
            } else {
                ch_bwa_index = [ [:], file(params.bwa_index) ]
            }
        } else {
            ch_bwa_index = BWA_INDEX ( ch_fasta ).index
            ch_versions  = ch_versions.mix(BWA_INDEX.out.versions)
        }
    
    

    //make chromosome size index



    //
    // Uncompress Bowtie2 index or generate from scratch if required
    //


    //
    // Uncompress CHROMAP index or generate from scratch if required
    //
    //
    // Uncompress STAR index or generate from scratch if required
    //

    emit:
    fasta         = ch_fasta                  //    path: genome.fasta
    bwa_index     = ch_bwa_index              //    path: bwa/index/
 
    versions    = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
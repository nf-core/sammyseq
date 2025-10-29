include {
    GUNZIP as GUNZIP_FASTA
    GUNZIP as GUNZIP_GTF
    GUNZIP as GUNZIP_GFF
    GUNZIP as GUNZIP_TSS_BED
    GUNZIP as GUNZIP_GENE_BED
    GUNZIP as GUNZIP_CHROM_SIZES
    GUNZIP as GUNZIP_FAI
    GUNZIP as GUNZIP_BLACKLIST } from '../../../modules/nf-core/gunzip/main'

include {
    UNTAR as UNTAR_BWA_INDEX
    UNTAR as UNTAR_BOWTIE2_INDEX
    } from '../../../modules/nf-core/untar/main'


include { UNTARFILES               } from '../../../modules/nf-core/untarfiles/main'
include { GFFREAD                  } from '../../../modules/nf-core/gffread/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_CHROM_SIZES } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_FAI         } from '../../../modules/nf-core/samtools/faidx/main'
include { BWA_INDEX                } from '../../../modules/nf-core/bwa/index/main'
include { BOWTIE2_BUILD            } from '../../../modules/nf-core/bowtie2/build/main'
include { BEDTOOLS_MAKEWINDOWS     } from '../../../modules/nf-core/bedtools/makewindows/main'
include { GENOME_BLACKLIST_REGIONS } from '../../../modules/local/genome_blacklist_regions'

workflow PREPARE_GENOME {

    take:
//    genome             //  string: genome name
//    genomes            //     map: genome attributes
    fasta              //    path: path to genome fasta file
    aligner            //    string: aligner name
    gtf                //    file: /path/to/genome.gtf
    tss_bed            //    file: /path/to/tss.bed
//    gff                //    file: /path/to/genome.gff
    blacklist          //    file: /path/to/blacklist.bed
    gene_bed           //    file: /path/to/gene.bed
    bwa_index          //    file: /path/to/bwa/index/
    bowtie2_index      //    file: /path/to/bowtie2/index/
    chrom_sizes        //    file: /path/to/genome.sizes
    fai                //    file: /path/to/genome.fai
    binsize            //    binsize: genome binning

    main:

    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    ch_fasta = Channel.empty()
    if (params.fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA ( [ [:], params.fasta ] ).gunzip.map{ it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = Channel.value(file(params.fasta))
    }

    //println(ch_fasta)
    // Make fasta file available if reference saved or IGV is run
    //if (params.save_reference || !params.skip_igv) {

    // if (params.save_reference) {
    //     file("${params.outdir}/genome/").mkdirs()
    //     ch_fasta.copyTo("${params.outdir}/genome/")
    // }

    //
    // Uncompress GTF annotation file
    //
    ch_gtf = Channel.empty()
    if (params.gtf) {
        if (params.gtf.endsWith('.gz')) {
            ch_gtf      = GUNZIP_GTF ( [ [:], params.gtf ] ).gunzip.map{ it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
        } else {
            ch_gtf = Channel.value(file(params.gtf))
        }
    }

    //
    //  Uncompress bed file
    //
    ch_tss_bed = Channel.empty()
    if (params.tss_bed) {
        if (params.tss_bed.endsWith('.gz')) {
            ch_tss_bed = GUNZIP_TSS_BED ( [ [:], params.tss_bed ] ).gunzip.map{ it[1] }
            ch_versions = ch_versions.mix(GUNZIP_TSS_BED.out.versions)
        } else {
            ch_tss_bed = Channel.value(file(params.tss_bed))
        }
    }

    //
    // Uncompress blacklist file if required
    //
    ch_blacklist = Channel.empty()
    if (params.blacklist) {
        if (params.blacklist.endsWith('.gz')) {
            ch_blacklist = GUNZIP_BLACKLIST ( [ [:], params.blacklist ] ).gunzip.map{ meta, file -> tuple([id: 'blacklist'], file) }
            ch_versions  = ch_versions.mix(GUNZIP_BLACKLIST.out.versions)
        } else {
            ch_blacklist = Channel.value(tuple([id: 'blacklist'], file(params.blacklist)))
        }
    }

    //
    // Create FASTA index
    //
    ch_fai = Channel.empty()
    if (params.fai) {
        if (params.fai.endsWith('.gz')) {
            ch_fai = GUNZIP_FAI([ [:], params.fai ]).gunzip.map{ it[1] }
            ch_versions = ch_versions.mix(GUNZIP_FAI.out.versions)
        } else {
            ch_fai = Channel.value(file(params.fai))
        }
    } else {
        SAMTOOLS_FAIDX_FAI(
            ch_fasta.map { tuple([:], it) },
            [[:], []],
            false
        )
        ch_fai = SAMTOOLS_FAIDX_FAI.out.fai.map { it[1] }
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_FAI.out.versions)
    }

    //
    // Create chromosome sizes file
    //
    ch_chrom_sizes = Channel.empty()
    if (params.chrom_sizes) {
        if (params.chrom_sizes.endsWith('.gz')) {
            ch_chrom_sizes = GUNZIP_CHROM_SIZES([ [:], params.chrom_sizes ]).gunzip.map{ it[1] }
            ch_versions = ch_versions.mix(GUNZIP_CHROM_SIZES.out.versions)
        } else {
            ch_chrom_sizes = Channel.value(file(params.chrom_sizes))
        }
    } else {
        SAMTOOLS_FAIDX_CHROM_SIZES(
            ch_fasta.map { tuple([:], it) },
            [[:], []],
            true
        )
        ch_chrom_sizes = SAMTOOLS_FAIDX_CHROM_SIZES.out.sizes.map { it[1] }
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_CHROM_SIZES.out.versions)
    }


    //
    // Prepare genome intervals for filtering by removing regions in blacklist file
    //
    ch_genome_filtered_bed = Channel.empty()

    GENOME_BLACKLIST_REGIONS (
        ch_chrom_sizes,
        ch_blacklist.ifEmpty(tuple([id: 'no_blacklist'], []))
    )
    ch_genome_filtered_bed = GENOME_BLACKLIST_REGIONS.out.bed
    ch_versions = ch_versions.mix(GENOME_BLACKLIST_REGIONS.out.versions)


    //
    // Prepare BWA index
    //
    ch_bwa_index = Channel.empty()
    if (params.aligner == 'bwaaln' || params.aligner == 'bwamem') {
        if (params.bwa_index) {
            if (params.bwa_index.endsWith('.tar.gz')) {
                ch_bwa_index = UNTAR_BWA_INDEX ( [ [:], params.bwa_index ] ).untar
                ch_versions  = ch_versions.mix(UNTAR_BWA_INDEX.out.versions)
            } else {
                ch_bwa_index = [ [:], file(params.bwa_index) ]
            }
        } else {
            ch_bwa_index = BWA_INDEX ( ch_fasta.map { [ [:], it ] } ).index
            ch_versions  = ch_versions.mix(BWA_INDEX.out.versions)
        }
    }

    //
    // Uncompress Bowtie2 index or generate from scratch if required
    //
    ch_bowtie2_index = Channel.empty()
    if (params.aligner == 'bowtie2') {
        if (params.bowtie2_index) {
            if (params.bowtie2_index.endsWith('.tar.gz')) {
                ch_bowtie2_index = UNTAR_BOWTIE2_INDEX ( [ [:], params.bowtie2_index ] ).untar
                ch_versions  = ch_versions.mix(UNTAR_BOWTIE2_INDEX.out.versions)
            } else {
                ch_bowtie2_index = [ [:], file(params.bowtie2_index) ]
            }
        } else {
            ch_bowtie2_index = BOWTIE2_BUILD ( ch_fasta.map { [ [:], it ] } ).index
            ch_versions      = ch_versions.mix(BOWTIE2_BUILD.out.versions)
        }
    }

    emit:
    fasta         = ch_fasta                  //    path: genome.fasta
    fai           = ch_fai                    //    path: genome.fai
    gtf           = ch_gtf                    //    path: genome.gtf
    tss_bed       = ch_tss_bed                //    path: tss.bed
    // gene_bed      = ch_gene_bed               //    path: gene.bed
    chrom_sizes   = ch_chrom_sizes            //    path: genome.sizes
    filtered_bed  = ch_genome_filtered_bed    //    path: *.include_regions.bed
    bwa_index     = ch_bwa_index              //    path: bwa/index/
    bowtie2_index = ch_bowtie2_index          //    path: bowtie2/index/
    blacklist     = ch_blacklist
    versions      = ch_versions               // channel: [ versions.yml ]
}

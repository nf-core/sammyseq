tiling_correctly <- function(df, df_row, tile_width) {
    # Function to divide a genomic region into bins of a given size
    vect_starts <- seq(df[df_row,]$start, to = df[df_row,]$end, by = binsize)
    vect_ends <- c()

    for (start in vect_starts) {
        if (tail(vect_starts, n = 1) != start) {
            vect_ends <- c(vect_ends, start + binsize - 1)
        } else {
            vect_ends <- c(vect_ends, df[df_row,]$end)
        }
    }

    df_tiledregions <- data.frame(
        "seqnames" = rep(df[1,]$seqnames, times = length(vect_starts)),
        "start" = vect_starts,
        "end" = vect_ends
    )

    tiledregions_gr <- makeGRangesFromDataFrame(df_tiledregions, keep.extra.columns = TRUE)
    return(tiledregions_gr)
}

bins_trimmer <- function(genome_bedfile, bin_and_bl, genome, dropping_chrs = c("")) {
    # Function to trim bins by removing unwanted chromosomes
    hg38_df <- fread(genome_bedfile, data.table = FALSE)
    names(hg38_df) <- c("seqnames", "start", "end")
    hg38_df <- as.data.frame(hg38_df)

    hg38_gr <- makeGRangesFromDataFrame(df = hg38_df, start.field = "start", end.field = "end", seqnames.field = "seqnames")
    seqlevels(hg38_gr) <- hg38_df$seqnames
    seqlengths(hg38_gr) <- hg38_df$end
    genome(hg38_gr) <- genome

    hg38_gr <- keepStandardChromosomes(hg38_gr, pruning.mode = "coarse")
    hg38_gr <- dropSeqlevels(hg38_gr, dropping_chrs, pruning.mode = "coarse")

    chr_seqlengths <- seqlengths(hg38_gr)
    bin_list <- keepSeqlevels(bin_and_bl, chr, pruning.mode = "coarse")

    seqlengths(bin_list) <- seqlengths(hg38_gr)
    genome(bin_list) <- genome

    return(bin_list)
}


# Function to generate bins while excluding blacklisted regions
new_make_bins <- function(chr, binsize, genome = "mm10", bed_file = "path",
                          blacklist_path = 'path',
                          gapw = 50000,
                          chrs = c("")) {
    dropping <- dropping_chrs[-(grep(paste0("\\b", chr, "\\b"), chrs))]
    bins_gr <- bins_calculator(genome_bedfile = bed_file, binsize = binsize, genome = genome, dropping_chrs = dropping)

    binsize <- binsize

    bl_gr <- rtracklayer::import(blacklist_path)
    bl_gr <- reduce(bl_gr, min.gapwidth = gapw)

    df_2 <- as.data.frame(setdiff(bins_gr, bl_gr))

# If only one region remains after removing blacklisted regions
    if (nrow(df_2) == 1) {
        gr_list_tiled_nobl <- tiling_correctly(df_2, 1, binsize)
        gr_tofill <- gr_list_tiled_nobl
    } else { # If multiple regions remain, tile them separately
        gr_list_tiled_nobl <- lapply(seq(nrow(df_2)), FUN = function(x) {
            gr_computed <- tiling_correctly(df_2, x, binsize)
            return(gr_computed)
        })

        gr_tofill <- gr_list_tiled_nobl[[1]]

        for (gr_to_add in seq(2, length(gr_list_tiled_nobl))) {
            gr_tofill <- c(gr_tofill, gr_list_tiled_nobl[[gr_to_add]])
        }
    }

    # Merge bins and blacklist regions
    notrim_merge <- c(gr_tofill, bl_gr)
    bins_trimmed <- bins_trimmer(genome_bedfile = bed_file, bin_and_bl = notrim_merge, genome = genome, dropping_chrs = dropping)

    print("Bin list REcalculated")
    bins_trimmed <- sort(bins_trimmed)

    return(bins_trimmed)
}

txdb <- GenomicFeatures::makeTxDbFromGFF(txdb_name)
genes <- GenomicFeatures::genes(txdb)

summarizeProteinCodingGenes <- function(txdb) {
    stopifnot(is(txdb, "TxDb"))
    protein_coding_tx <- names(GenomicFeatures::cdsBy(txdb, use.names=TRUE))
    all_tx <- mcols(GenomicFeatures::transcripts(txdb, columns=c("gene_id", "tx_name")))
    all_tx$gene_id <- as.character(all_tx$gene_id)
    all_tx$is_coding <- all_tx$tx_name %in% protein_coding_tx
    tmp <- splitAsList(all_tx$is_coding, all_tx$gene_id)
    gene <- names(tmp)
    n_tx <- unname(sum(tmp, na.rm = TRUE))
    n_coding <- unname(sum(tmp))
    n_non_coding <- n_tx - n_coding
    data.frame(gene, n_tx, n_coding, n_non_coding, stringsAsFactors=FALSE)
}

# Precompute protein-coding genes list
geneid_codingdf <- summarizeProteinCodingGenes(txdb)
codingenes_grobj <- genes[genes$gene_id %in% c(geneid_codingdf[geneid_codingdf$n_coding > 0,]$gene)]
genes_gr <- codingenes_grobj

# Convert HEX colors to RGB format
rgb_str <- function(hex) {
    paste(as.vector(col2rgb(hex)), collapse = ",")
}




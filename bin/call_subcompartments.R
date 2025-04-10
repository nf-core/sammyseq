#!/usr/bin/env Rscript

suppressMessages({
  require(optparse)
  require(data.table)
  require(doParallel)
  require(GenomicRanges)
  require(factoextra)
  require(parallel)
  require(rtracklayer)
  require(ggplot2)
  require(R.utils)
  require(R.oo)
  require(R.methodsS3)
  require(maptools)
  require(patchwork)
  require(sp)
  require(dendextend)
  require(Gviz)
  require(CALDER)
  require(ape)
  require(fitdistrplus)
  require(igraph)
  require(Matrix)
  require(rARPACK)
  require(fields)
  require(strawr)
})

## https://community.seqera.io/t/source-another-r-script-in-the-bin-directory/1059
path <- Sys.getenv("PATH") |> strsplit(":")
bin_path <- tail(path[[1]], n=1)
source(file.path(bin_path, "functions_subcompartments.R"))

# CLI options
option_list <- list(
  make_option(c("-i", "--input_file"), type = "character", help = "Input TSV file"),
  make_option(c("-b", "--binsize"), type = "integer", help = "Bin size"),
  make_option(c("-g", "--gene_gtf"), type = "character", help = "Gene GTF file"),
  make_option(c("-m", "--chrom_beds"), type = "character", help = "Single chromosome BED")
)
opt <- parse_args(OptionParser(option_list = option_list))

# Set option to allow arbitrary chromosome identifiers
options(ucscChromosomeNames=FALSE)

# Output dirs

sub2_colors <- c("B" = "#4575b4", "A" = "#d73027")

# Input TSV
comp_df <- fread(opt$input_file, data.table = FALSE)

# Gene GTF
genes_gr <- import(opt$gene_gtf, format = "GTF")

# Chromosome BED BINNATO
bed_file <- opt$chrom_beds
bins_gr <- import(bed_file, format = "BED")

chroms_in_bed <- unique(seqnames(bins_gr))
chr <- as.character(chroms_in_bed[1])
subs_file <- paste0(chr, "_compartment.Rdata")

# Run SAMMY
sub_objs <- call_subcompartments_sammy(
  patients = unique(comp_df$Patient_name),
  tracks_db = comp_df,
  bins_gr = bins_gr,
  subs_file = subs_file,
  binsize = opt$binsize,
  chr = chr,
  genes_gr = genes_gr,
  keeping_bins1 = "all",
  sublevel = "sub.2",
  sub_colors = sub2_colors
)

generate_files(sub_objs, chr)

#!/usr/bin/env Rscript

suppressMessages({
    library(parallel)
    library(data.table)
    library(GenomicRanges)
    library(rtracklayer)
    library(patchwork)
    library(Gviz)
    library(CALDER)
    library(sammyR)
})

options(ucscChromosomeNames=FALSE)

## Get parameters from Nextflow template
patient    <- "${patient}"
chromosome <- "${meta.chromosome}"
binsize    <-  as.integer(${binsize})
input_file <- "${csv}"
chrom_bed  <- "${chr_bed}"
gene_gtf   <- "${gtf}"

## Load data
comp_df <- fread(input_file, data.table = FALSE)
comp_df_replica <- comp_df[comp_df\$Patient_name == patient, ]

bins_gr <- import(chrom_bed, format = "BED")
genes_gr <- import(gene_gtf, format = "GTF")

sub2_colors <- c("B" = "#4575b4", "A" = "#d73027")
subs_file <- paste0(patient, "_compartment___", chromosome, '_', binsize, ".Rdata")

# Run SAMMY
sub_objs <- call_subcompartments_sammy(
    patients = patient,
    tracks_db = comp_df_replica,
    bins_gr = bins_gr,
    subs_file = subs_file,
    binsize = binsize,
    chr = chromosome,
    genes_gr = genes_gr,
    keeping_bins1 = "all",
    sublevel = "sub.2",
    sub_colors = sub2_colors
)

generate_files(sub_objs, chromosome)

# Write versions
pkgs <- c("sammyR", "CALDER", "GenomicRanges", "rtracklayer", "Gviz", "data.table")
writeLines(c(
  '"${task.process}":',
  paste0('    r-base: "', paste(R.Version()[c("major","minor")], collapse="."), '"'),
  vapply(pkgs, function(p) paste0('    ', p, ': ', as.character(packageVersion(p))), character(1))
), "versions.yml")